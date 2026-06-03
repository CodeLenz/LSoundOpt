#
# Penalização das variáveis de folga usadas no MILP.
#
# A penalização cresce com o número de trocas permitidas em uma iteração.
# Assim, uma folga não fica barata quando a trust-region é grande.
#
function Slack_Penalty(c, g_move; safety_factor=2.0, offset=10.0)

    # Maior ganho linear associado a uma troca de material
    max_c = isempty(c) ? 0.0 : maximum(abs.(c))

    # Penalização da folga
    return safety_factor * max(1, g_move) * max_c + offset

end

#
# Loop de trust-region para o passo discreto do ISLP.
#
# Primeiro resolve um MILP linearizado. Depois testa a geometria e recalcula
# o problema acústico para decidir se o passo deve ser aceito.
#
function Trust_Region_Loop(c, A_glob, b_glob, gamma_curr, elems_design, cv,
                           obj_curr, Past, neighedge, fem_data,
                           solution_cache, MP,
                           alg::Algorithm_Params=Algorithm_Params())

    # Desempacota os dados do problema acústico
    nn, ne, coord, connect, f_rho, f_kappa, mu, freqs, livres, vels, press,
    nodes_target, vA = fem_data

    # Variáveis de projeto
    nvp = length(elems_design)
    gamma_design = gamma_curr[elems_design]
    is_air = gamma_design .< 0.5

    # Estado inicial do loop
    step_accepted = false
    n_reject = 0
    last_rej_dgamma = Float64[]
    gamma_final = copy(gamma_curr)

    # Restrições de movimento para elementos que hoje são ar
    one_air_ext = zeros(2 * nvp)
    one_air_ext[1:nvp] = Float64.(is_air)

    # Restrições de movimento para elementos que hoje são sólido
    one_solid_ext = zeros(2 * nvp)
    one_solid_ext[1:nvp] = Float64.(.!is_air)

    # Menor cv que ainda permite pelo menos um movimento discreto
    cv_min = 1.0001 / nvp

    while !step_accepted

        # Número máximo de trocas ar/sólido nesta tentativa
        limit_val = ceil(Int, cv * nvp)
        g_move = isodd(limit_val) ? limit_val - 1 : limit_val
        g_move = max(0, g_move)

        # Custo aumentado com as folgas
        beta = Slack_Penalty(c, g_move;
                             safety_factor=alg.slack_safety_factor,
                             offset=alg.slack_offset)
        c_aug = vcat(c, fill(beta, nvp))

        # Acrescenta as restrições de movimento à linearização global
        b_loc = vcat(b_glob, g_move, g_move)
        A_loc = vcat(A_glob, transpose(one_air_ext), -transpose(one_solid_ext))

        # Resolve o MILP
        sol_aug = LP(c_aug, A_loc, b_loc, gamma_design)
        delta_gamma = sol_aug[1:nvp]

        # Monitora se as restrições topológicas precisaram de folga
        slack_used = sum(sol_aug[nvp+1:end])
        if slack_used > 0
            println("    -> Folgas usadas no MILP: $(slack_used) (beta=$(round(beta, digits=4))).")
        end

        # Evita repetir um passo que já foi rejeitado
        if !isempty(last_rej_dgamma) && (delta_gamma == last_rej_dgamma)
            println("    -> Deadlock no Trust Region (ILP repetiu passo rejeitado).")
            println("    -> Desistindo desta iteracao (passo nulo) para evitar piora.")
            step_accepted = true
            break
        end

        # Se o MILP não moveu nada, a iteração externa decide o que fazer
        if sum(abs.(delta_gamma)) == 0
            println("    -> Passo nulo do ILP. Proxima iteracao.")
            step_accepted = true
            break
        end

        # Redução prevista pelo modelo linear
        pred_red = -sum(c .* delta_gamma)
        if abs(pred_red) < alg.tr_pred_tol
            println("    -> Reducao desprezivel. Aceitando.")
            step_accepted = true
            break
        end

        # Aplica o passo candidato
        gamma_trial = copy(gamma_curr)
        gamma_trial[elems_design] .+= delta_gamma
        gamma_trial[elems_design] .= round.(gamma_trial[elems_design])

        # Teste geométrico: o perímetro não pode ultrapassar o limite imposto
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

        # Hash da topologia candidata
        cfg_hash = hash(gamma_trial)
        obj_trial = 0.0

        # Reaproveita avaliações acústicas já feitas
        if haskey(solution_cache, cfg_hash)
            obj_trial = solution_cache[cfg_hash]
            println("    -> Usando o cache!")
        else
            # Avalia o problema acústico para a topologia candidata
            Sweep!(nn, ne, coord, connect, gamma_trial, f_rho, f_kappa, mu,
                   freqs, livres, vels, press, MP)
            obj_trial = Objetivo(MP, nodes_target, vA)
            solution_cache[cfg_hash] = obj_trial
        end

        # Redução atual
        act_red = obj_curr - obj_trial

        # Razão entre redução real e redução prevista
        if act_red < -alg.tr_physics_tol
            println("    -> REJEITADO (Fisica). SPL piorou (Red=$(round(act_red, digits=4))).")
            R = -1.0
        else
            R = abs(pred_red) < alg.tr_ratio_pred_tol ? 0.0 : act_red / pred_red
        end

        println("--- TR: Pred=$(round(pred_red, digits=4)) Act=$(round(act_red, digits=4)) R=$(round(R, digits=2)) cv=$(round(cv, digits=4))")

        # Atualiza o raio da trust-region conforme a qualidade do passo
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

        # Limite de tentativas dentro da mesma iteração externa
        if n_reject >= alg.tr_max_reject
            println("Trust Region exausta ($(alg.tr_max_reject) rejeicoes). Desistindo e mantendo atual.")
            step_accepted = true
            break
        end
    end

    # Retorna a topologia aceita e o cv atualizado
    return gamma_final, cv, step_accepted

end
