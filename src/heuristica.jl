#
# Determina o fator de atualização ar/sólido.
#
# A ideia aqui é recuperar o raio de busca quando a otimização fica parada.
# O valor de referência cv_ref vem do YAML.
#
function Update_Heuristics(iter, stag_count, cv, nvp, cv_ref,
                           alg::Algorithm_Params=Algorithm_Params())

    # Primeiro nível de recuperação
    if stag_count == 1

        val = alg.heur_stag_l1_factor * cv_ref
        println("--- HEURISTICA: Recovery L1 (Reset para Nominal: $(round(val, digits=2)))")

        return min(val, alg.heur_cv_max)

    # Segundo nível de recuperação
    elseif stag_count == 2

        val = alg.heur_stag_l2_factor * cv_ref
        println("--- HEURISTICA: Recovery L2 (Boost: $(round(val, digits=2)))")

        return min(val, alg.heur_cv_max)

    # Terceiro nível de recuperação
    elseif stag_count >= 3

        val = alg.heur_stag_l3_factor * cv_ref
        println("--- HEURISTICA: Recovery L3 (Kick Maximo: $(round(val, digits=2)))")

        return min(val, alg.heur_cv_max)

    end

    # Se o cv ficou menor que o necessário para mover poucos elementos,
    # volta para o valor nominal
    limit_min = alg.heur_min_moves / nvp
    if cv < limit_min

        println("--- HEURISTICA: Limite Minimo Atingido. Resetando para Nominal.")
        return cv_ref

    end

    # Caso padrão: mantém o cv atual
    return cv

end

#
# Versão anterior da heurística.
#
# Mantida apenas para comparação com rodadas antigas.
#
function Update_Heuristics_anterior(iter, stag_count, cv, nvp)

    # Valor inicial
    cv_new = cv

    # Recuperação por estagnação
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

    # Aquecimento inicial
    if iter <= 5 && stag_count == 0 && cv_new > 0.05
        cv_new = 0.05
    end

    # Retorna o novo fator de atualização
    return cv_new

end
