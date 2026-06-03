#
# Parâmetros algorítmicos usados pela otimização discreta.
#
# Os valores abaixo são os padrões usados quando a seção "algoritmo"
# não aparece no arquivo YAML.
#
Base.@kwdef struct Algorithm_Params

    # Trust-region
    tr_accept_ratio_min::Float64 = 0.0
    tr_expand_ratio::Float64 = 0.7
    tr_shrink_ratio::Float64 = 0.3
    tr_expand_factor::Float64 = 1.2
    tr_shrink_factor::Float64 = 0.5
    tr_cv_max::Float64 = 0.5
    tr_max_reject::Int = 10

    # Tolerâncias numéricas
    tr_perimeter_tol::Float64 = 1.0e-4
    tr_physics_tol::Float64 = 1.0e-9
    tr_pred_tol::Float64 = 1.0e-12
    tr_ratio_pred_tol::Float64 = 1.0e-10

    # Penalização das folgas do MILP
    slack_safety_factor::Float64 = 2.0
    slack_offset::Float64 = 10.0

    # Heurísticas de atualização do cv
    heur_cv_max::Float64 = 0.5
    heur_stag_l1_factor::Float64 = 1.0
    heur_stag_l2_factor::Float64 = 1.5
    heur_stag_l3_factor::Float64 = 3.0
    heur_min_moves::Float64 = 2.0

    # Critério de parada por estagnação
    hard_stop_stagnation::Int = 4

end

#
# Lê um valor do dicionário YAML.
#
# Quando a chave não existe, retorna o valor padrão.
#
function _yaml_value(dados, key::String, default)

    # Chave ausente
    if !haskey(dados, key)
        return default
    end

    # Valor encontrado no YAML
    value = dados[key]

    # Inteiros
    if isa(default, Int)
        return isa(value, String) ? parse(Int, value) : Int(value)
    end

    # Reais
    return isa(value, String) ? parse(Float64, value) : Float64(value)

end

#
# Lê a seção "algoritmo" do YAML.
#
function Le_YAML_Algoritmo(arquivo::AbstractString)

    # Arquivo YAML completo
    dados = YAML.load_file(arquivo)

    # Seção de parâmetros algorítmicos
    alg = haskey(dados, "algoritmo") ? dados["algoritmo"] :
          haskey(dados, "algorithm") ? dados["algorithm"] :
          Dict{Any, Any}()

    # Monta a estrutura com valores lidos ou padrões
    params = Algorithm_Params(
        tr_accept_ratio_min = _yaml_value(alg, "tr_accept_ratio_min", 0.0),
        tr_expand_ratio = _yaml_value(alg, "tr_expand_ratio", 0.7),
        tr_shrink_ratio = _yaml_value(alg, "tr_shrink_ratio", 0.3),
        tr_expand_factor = _yaml_value(alg, "tr_expand_factor", 1.2),
        tr_shrink_factor = _yaml_value(alg, "tr_shrink_factor", 0.5),
        tr_cv_max = _yaml_value(alg, "tr_cv_max", 0.5),
        tr_max_reject = _yaml_value(alg, "tr_max_reject", 10),
        tr_perimeter_tol = _yaml_value(alg, "tr_perimeter_tol", 1.0e-4),
        tr_physics_tol = _yaml_value(alg, "tr_physics_tol", 1.0e-9),
        tr_pred_tol = _yaml_value(alg, "tr_pred_tol", 1.0e-12),
        tr_ratio_pred_tol = _yaml_value(alg, "tr_ratio_pred_tol", 1.0e-10),
        slack_safety_factor = _yaml_value(alg, "slack_safety_factor", 2.0),
        slack_offset = _yaml_value(alg, "slack_offset", 10.0),
        heur_cv_max = _yaml_value(alg, "heur_cv_max", 0.5),
        heur_stag_l1_factor = _yaml_value(alg, "heur_stag_l1_factor", 1.0),
        heur_stag_l2_factor = _yaml_value(alg, "heur_stag_l2_factor", 1.5),
        heur_stag_l3_factor = _yaml_value(alg, "heur_stag_l3_factor", 3.0),
        heur_min_moves = _yaml_value(alg, "heur_min_moves", 2.0),
        hard_stop_stagnation = _yaml_value(alg, "hard_stop_stagnation", 4),
    )

    # Verifica consistência dos parâmetros
    Validate_Algorithm_Params(params)

    # Retorna os parâmetros algorítmicos
    return params

end

#
# Testa consistência dos parâmetros algorítmicos.
#
function Validate_Algorithm_Params(p::Algorithm_Params)

    # Trust-region
    p.tr_accept_ratio_min >= 0 || throw("Algorithm_Params::tr_accept_ratio_min deve ser >= 0")
    p.tr_expand_ratio >= p.tr_shrink_ratio || throw("Algorithm_Params::tr_expand_ratio deve ser >= tr_shrink_ratio")
    p.tr_expand_factor >= 1 || throw("Algorithm_Params::tr_expand_factor deve ser >= 1")
    p.tr_shrink_factor > 0 || throw("Algorithm_Params::tr_shrink_factor deve ser > 0")
    p.tr_shrink_factor < 1 || throw("Algorithm_Params::tr_shrink_factor deve ser < 1")
    p.tr_cv_max > 0 || throw("Algorithm_Params::tr_cv_max deve ser > 0")
    p.tr_cv_max < 1 || throw("Algorithm_Params::tr_cv_max deve ser < 1")
    p.tr_max_reject >= 1 || throw("Algorithm_Params::tr_max_reject deve ser >= 1")

    # Tolerâncias
    p.tr_perimeter_tol >= 0 || throw("Algorithm_Params::tr_perimeter_tol deve ser >= 0")
    p.tr_physics_tol >= 0 || throw("Algorithm_Params::tr_physics_tol deve ser >= 0")
    p.tr_pred_tol >= 0 || throw("Algorithm_Params::tr_pred_tol deve ser >= 0")
    p.tr_ratio_pred_tol >= 0 || throw("Algorithm_Params::tr_ratio_pred_tol deve ser >= 0")

    # Folgas
    p.slack_safety_factor >= 0 || throw("Algorithm_Params::slack_safety_factor deve ser >= 0")
    p.slack_offset >= 0 || throw("Algorithm_Params::slack_offset deve ser >= 0")

    # Heurísticas
    p.heur_cv_max > 0 || throw("Algorithm_Params::heur_cv_max deve ser > 0")
    p.heur_cv_max < 1 || throw("Algorithm_Params::heur_cv_max deve ser < 1")
    p.heur_stag_l1_factor >= 0 || throw("Algorithm_Params::heur_stag_l1_factor deve ser >= 0")
    p.heur_stag_l2_factor >= 0 || throw("Algorithm_Params::heur_stag_l2_factor deve ser >= 0")
    p.heur_stag_l3_factor >= 0 || throw("Algorithm_Params::heur_stag_l3_factor deve ser >= 0")
    p.heur_min_moves > 0 || throw("Algorithm_Params::heur_min_moves deve ser > 0")

    # Parada
    p.hard_stop_stagnation >= 1 || throw("Algorithm_Params::hard_stop_stagnation deve ser >= 1")

    return nothing

end
