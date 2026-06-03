#
# Determines the air/solid update factor.
# cv_ref is the initial value read from the YAML file.
#
function Update_Heuristics(iter, stag_count, cv, nvp, cv_ref,
                           alg::Algorithm_Params=Algorithm_Params())

    if stag_count == 1
        val = alg.heur_stag_l1_factor * cv_ref
        println("--- HEURISTICA: Recovery L1 (Reset para Nominal: $(round(val, digits=2)))")
        return min(val, alg.heur_cv_max)

    elseif stag_count == 2
        val = alg.heur_stag_l2_factor * cv_ref
        println("--- HEURISTICA: Recovery L2 (Boost: $(round(val, digits=2)))")
        return min(val, alg.heur_cv_max)

    elseif stag_count >= 3
        val = alg.heur_stag_l3_factor * cv_ref
        println("--- HEURISTICA: Recovery L3 (Kick Maximo: $(round(val, digits=2)))")
        return min(val, alg.heur_cv_max)
    end

    limit_min = alg.heur_min_moves / nvp
    if cv < limit_min
        println("--- HEURISTICA: Limite Minimo Atingido. Resetando para Nominal.")
        return cv_ref
    end

    return cv
end

#
# Previous heuristic kept for comparison with older runs.
#
function Update_Heuristics_anterior(iter, stag_count, cv, nvp)
    cv_new = cv

    if stag_count == 1
       println("--- TRAVADO L1: Reseta cv=5%")
       cv_new = 0.05
    elseif stag_count == 2
        println("--- TRAVADO L2: Aumenta cv=15%")
        cv_new = 0.15
    elseif stag_count >= 3
        println("--- TRAVADO L3: Apela cv=30%")
        cv_new = 0.30
    elseif cv < (2.0 / nvp)
        println("--- TRAVADO: Reseta cv=5%")
        cv_new = 0.05
    end

    if iter <= 5 && stag_count == 0 && cv_new > 0.05
        cv_new = 0.05
    end

    return cv_new
end
