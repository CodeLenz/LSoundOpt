#
# Determina o fator de atualização ar/sólido
# cv_ref = valor inicial lido do YAML (ex: 0.05 ou 0.30)
#
function Update_Heuristics(iter, stag_count, cv, nvp, cv_ref)
    
    # Teto máximo de segurança (50% do domínio)
    # Acima disso vira "chute aleatório"
    CV_MAX = 0.50

    # 1. Recuperação de Estagnação (Baseada no cv_ref)
    if stag_count == 1
       val = cv_ref
       println("--- HEURÍSTICA: Recovery L1 (Reset para Nominal: $(round(val,digits=2)))")
       return min(val, CV_MAX)
       
    elseif stag_count == 2
        # Tenta ser mais agressivo que o inicial (1.5x a 2.0x)
        val = 1.5 * cv_ref
        println("--- HEURÍSTICA: Recovery L2 (Boost: $(round(val,digits=2)))")
        return min(val, CV_MAX)
        
    elseif stag_count >= 3
        # Agressividade máxima para sair do buraco (3.0x)
        val = 3.0 * cv_ref
        println("--- HEURÍSTICA: Recovery L3 (Kick Máximo: $(round(val,digits=2)))")
        return min(val, CV_MAX)
    end

    # 2. Proteção contra congelamento (Piso Mínimo)
    # Se o raio atual for menor que 2 elementos, reseta para
    # um valor pequeno mas operável (sugestão: 10% do valor de referência ou 5% absoluto)
    limit_min = 2.0/nvp
    if cv < limit_min
        # Reseta para o valor nominal para tentar novamente
        println("--- HEURÍSTICA: Limite Mínimo Atingido. Resetando para Nominal.")
        return cv_ref 
    end

    # 3. Comportamento Padrão
    return cv
    
end

#
# Determina o fator de atualização ar/sólido
#
function Update_Heuristics_anterior(iter, stag_count, cv, nvp)
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