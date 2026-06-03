#
# Trust-region for the discrete ISLP update.
#
function Slack_Penalty(c, g_move; safety_factor=2.0, offset=10.0)
    max_c = isempty(c) ? 0.0 : maximum(abs.(c))
    return safety_factor * max(1, g_move) * max_c + offset
end

function Trust_Region_Loop(c, A_glob, b_glob, gamma_curr, elems_design, cv,
                           obj_curr, Past, neighedge, fem_data,
                           solution_cache, MP,
                           alg::Algorithm_Params=Algorithm_Params())

    nn, ne, coord, connect, f_rho, f_kappa, mu, freqs, livres, vels, press,
    nodes_target, vA = fem_data

    nvp = length(elems_design)
    gamma_design = gamma_curr[elems_design]
    is_air = gamma_design .< 0.5

    step_accepted = false
    n_reject = 0
    last_rej_dgamma = Float64[]
    gamma_final = copy(gamma_curr)

    one_air_ext = zeros(2 * nvp)
    one_air_ext[1:nvp] = Float64.(is_air)

    one_solid_ext = zeros(2 * nvp)
    one_solid_ext[1:nvp] = Float64.(.!is_air)

    cv_min = 1.0001 / nvp

    while !step_accepted
        limit_val = ceil(Int, cv * nvp)
        g_move = isodd(limit_val) ? limit_val - 1 : limit_val
        g_move = max(0, g_move)

        # Slack must be expensive relative to the whole move budget, not just
        # to a single design coefficient. Otherwise one slack can buy many
        # favorable flips in the same MILP step.
        beta = Slack_Penalty(c, g_move;
                             safety_factor=alg.slack_safety_factor,
                             offset=alg.slack_offset)
        c_aug = vcat(c, fill(beta, nvp))

        b_loc = vcat(b_glob, g_move, g_move)
        A_loc = vcat(A_glob, transpose(one_air_ext), -transpose(one_solid_ext))

        sol_aug = LP(c_aug, A_loc, b_loc, gamma_design)
        delta_gamma = sol_aug[1:nvp]
        slack_used = sum(sol_aug[nvp+1:end])
        if slack_used > 0
            println("    -> Folgas usadas no MILP: $(slack_used) (beta=$(round(beta, digits=4))).")
        end

        if !isempty(last_rej_dgamma) && (delta_gamma == last_rej_dgamma)
            println("    -> Deadlock no Trust Region (ILP repetiu passo rejeitado).")
            println("    -> Desistindo desta iteracao (passo nulo) para evitar piora.")
            step_accepted = true
            break
        end

        if sum(abs.(delta_gamma)) == 0
            println("    -> Passo nulo do ILP. Proxima iteracao.")
            step_accepted = true
            break
        end

        pred_red = -sum(c .* delta_gamma)
        if abs(pred_red) < alg.tr_pred_tol
            println("    -> Reducao desprezivel. Aceitando.")
            step_accepted = true
            break
        end

        gamma_trial = copy(gamma_curr)
        gamma_trial[elems_design] .+= delta_gamma
        gamma_trial[elems_design] .= round.(gamma_trial[elems_design])

        P_trial = Perimiter(gamma_trial, neighedge, elems_design)
        if (Past > 0) && (P_trial > (Past + alg.tr_perimeter_tol))
            println("    -> REJEITADO (Geometria). P=$P_trial. Reduzindo cv.")
            last_rej_dgamma = copy(delta_gamma)
            cv = max(cv_min, cv * alg.tr_shrink_factor)
            n_reject += 1
            if n_reject >= alg.tr_max_reject
                println("TR exausta.")
                break
            end
            continue
        end

        cfg_hash = hash(gamma_trial)
        obj_trial = 0.0

        if haskey(solution_cache, cfg_hash)
            obj_trial = solution_cache[cfg_hash]
            println("    -> Usando o cache!")
        else
            Sweep!(nn, ne, coord, connect, gamma_trial, f_rho, f_kappa, mu,
                   freqs, livres, vels, press, MP)
            obj_trial = Objetivo(MP, nodes_target, vA)
            solution_cache[cfg_hash] = obj_trial
        end

        act_red = obj_curr - obj_trial

        if act_red < -alg.tr_physics_tol
            println("    -> REJEITADO (Fisica). SPL piorou (Red=$(round(act_red, digits=4))).")
            R = -1.0
        else
            R = abs(pred_red) < alg.tr_ratio_pred_tol ? 0.0 : act_red / pred_red
        end

        println("--- TR: Pred=$(round(pred_red, digits=4)) Act=$(round(act_red, digits=4)) R=$(round(R, digits=2)) cv=$(round(cv, digits=4))")

        if R <= alg.tr_accept_ratio_min
            last_rej_dgamma = copy(delta_gamma)
            cv = max(cv_min, cv * alg.tr_shrink_factor)
            n_reject += 1
        else
            step_accepted = true
            gamma_final .= gamma_trial

            if R > alg.tr_expand_ratio
                cv = min(alg.tr_cv_max, cv * alg.tr_expand_factor)
            elseif R < alg.tr_shrink_ratio
                cv = max(cv_min, cv * alg.tr_shrink_factor)
            end
        end

        if n_reject >= alg.tr_max_reject
            println("Trust Region exausta ($(alg.tr_max_reject) rejeicoes). Desistindo e mantendo atual.")
            step_accepted = true
            break
        end
    end

    return gamma_final, cv, step_accepted
end
