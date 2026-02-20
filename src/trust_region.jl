#
# Trust-Region (Com Cache de Soluções)
#
function Trust_Region_Loop(c, A_glob, b_glob, γ_curr, elems_design, cv, obj_curr, Past, neighedge, fem_data, solution_cache)
    
    # Desempacota os dados
    nn, ne, coord, connect, fρ, fκ, μ, freqs, livres, vels, press, nodes_target, vA = fem_data
    nvp = length(elems_design)
    
    γ_design = γ_curr[elems_design]
    is_air = γ_design .< 0.5
    
    # Custo aumentado com penalidade de folga
    beta = 100.0 * maximum(abs.(c)) + 10.0
    c_aug = vcat(c, fill(beta, nvp))

    step_accepted = false
    n_reject = 0
    last_rej_dgamma = Float64[]
    
    # Inicializa γ_final como a configuração ATUAL.
    # Se o loop falhar/desistir, retornamos isso (passo nulo).
    γ_final = copy(γ_curr)
    
    one_air_ext = zeros(2*nvp)
    one_air_ext[1:nvp] = Float64.(is_air)

    one_solid_ext = zeros(2*nvp)
    one_solid_ext[1:nvp] = Float64.(.!is_air)

    cv_min = 1.0001 / nvp

    while !step_accepted
        
        limit_val = ceil(Int, cv * nvp)
        g_move = isodd(limit_val) ? limit_val - 1 : limit_val
        g_move = max(0, g_move)

        b_loc = vcat(b_glob, g_move, g_move)
        A_loc = vcat(A_glob, transpose(one_air_ext), -transpose(one_solid_ext))

        # Resolve ILP
        sol_aug = LP(c_aug, A_loc, b_loc, γ_design)
        Δγ = sol_aug[1:nvp]

        # --- CORREÇÃO DE DEADLOCK ---
        if !isempty(last_rej_dgamma) && (Δγ == last_rej_dgamma)
            println("    -> Deadlock no Trust Region (ILP repetiu passo rejeitado).")
            println("    -> Desistindo desta iteração (Passo Nulo) para evitar piora.")
            
            # Aceitamos "sair do loop", mas NÃO aplicamos o Δγ. Retorna o γ_final original.
            step_accepted = true
            break 
        end

        if sum(abs.(Δγ)) == 0
            println("    -> Passo nulo do ILP. Próxima iteração."); step_accepted = true; break
        end

        pred_red = -sum(c .* Δγ)
        if abs(pred_red) < 1e-12
            println("    -> Redução desprezível. Aceitando."); step_accepted = true; break
        end

        # Trial Step
        γ_trial = copy(γ_curr)
        γ_trial[elems_design] .+= Δγ
        γ_trial[elems_design] .= round.(γ_trial[elems_design]) 

        # Check Geométrico
        P_trial = Perimiter(γ_trial, neighedge, elems_design)
        if (Past > 0) && (P_trial > (Past + 1e-4))
            println("    -> REJEITADO (Geometria). P=$P_trial. Reduzindo cv.")
            last_rej_dgamma = copy(Δγ)
            cv = max(cv_min, cv * 0.5)
            n_reject += 1
            if n_reject >= 10
                println("TR exausta.")
                break
            end
            continue
        end

        # Check Físico (COM CACHE)
        # Gera Hash da topologia candidata
        cfg_hash = hash(γ_trial)
        obj_trial = 0.0

        if haskey(solution_cache, cfg_hash)
            # Se já calculamos essa topologia antes, recuperamos o valor
            obj_trial = solution_cache[cfg_hash]
            println("    -> Usando o cache!") 
        else
            # Se não, calculamos o FEM (custoso)
            MP_t, _, _, _ = Sweep(nn, ne, coord, connect, γ_trial, fρ, fκ, μ, freqs, livres, vels, press)
            obj_trial = Objetivo(MP_t, nodes_target, vA)
            
            # Armazena no cache para o futuro
            solution_cache[cfg_hash] = obj_trial
        end

        act_red = obj_curr - obj_trial

        # Correção de Monotonicidade
        if act_red < -1e-9
             println("    -> REJEITADO (Física). SPL piorou (Red=$(round(act_red,digits=4))).")
             R = -1.0
        else
             R = abs(pred_red) < 1e-10 ? 0.0 : act_red / pred_red
        end

        println("--- TR: Pred=$(round(pred_red,digits=4)) Act=$(round(act_red,digits=4)) R=$(round(R,digits=2)) cv=$(round(cv,digits=4))")

        if R <= 0.0
            last_rej_dgamma = copy(Δγ)
            cv = max(cv_min, cv * 0.5)
            n_reject += 1
        else
            step_accepted = true
            γ_final .= γ_trial 
            
            if R > 0.7
                cv = min(0.5, cv * 1.2)
            elseif R < 0.3
                cv = max(cv_min, cv * 0.5)
            end
        end

        if n_reject >= 10
            println("Trust Region exausta (10 rejeições). Desistindo e mantendo atual.")
            step_accepted = true # Sai do loop sem aplicar mudança (γ_final = γ_curr)
            break
        end
    end

    return γ_final, cv, step_accepted
end