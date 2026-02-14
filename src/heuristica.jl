#
# Determina o fator de atualização ar/sólido
#
function Update_Heuristics(iter, stag_count, cv, nvp)
    cv_new = cv

    # Recovery Logic
    if stag_count == 1;
       println("--- TRAVADO L1: Reseta cv=5%")
       cv_new = 0.05
    elseif stag_count == 2
        println("--- TRAVADO L2: Aumenta cv=15%")
        cv_new = 0.15
    elseif stag_count >= 3
        println("--- TRAVADO L3: Apela cv=30%")
        cv_new = 0.30
    elseif cv < (2.0/nvp)
        println("--- TRAVADO: Reseta cv=5%")
        cv_new = 0.05
    end

    # Warmup Logic
    if iter <= 5 && stag_count == 0 && cv_new > 0.05
        cv_new = 0.05
    end

    # Retorna o fator de atualização ar/sólido
    return cv_new
    
end