#
# Resolve o subproblema linear inteiro.
#
# As primeiras variáveis são os incrementos discretos Δγ.
# As variáveis restantes, quando existem, são folgas das restrições locais.
#
function LP(c, A, b, gamma_design)

    # Número de variáveis físicas de projeto
    n_design = length(gamma_design)

    # Número total de variáveis no MILP
    n_total = length(c)
    n_slack = n_total - n_design

    # Modelo inteiro resolvido com HiGHS
    model = Model(optimizer_with_attributes(HiGHS.Optimizer,
                  "output_flag" => false,
                  "presolve" => "on",
                  "time_limit" => 60.0))

    @variable(model, x[1:n_total], Int)

    # Limites das variáveis físicas:
    # γ = 0  ->  Δγ ∈ {0, 1}
    # γ = 1  ->  Δγ ∈ {-1, 0}
    for i = 1:n_design

        lb = -round(Int, gamma_design[i])
        ub =  round(Int, 1 - gamma_design[i])

        set_lower_bound(x[i], lb)
        set_upper_bound(x[i], ub)

    end

    # As folgas são inteiras e não-negativas
    if n_slack > 0
        for i = (n_design + 1):n_total
            set_lower_bound(x[i], 0.0)
        end
    end

    # Restrições linearizadas
    @constraint(model, A * x .<= b)

    # Função objetivo linear
    @objective(model, Min, c' * x)

    # Resolve o MILP
    JuMP.optimize!(model)
    status = termination_status(model)

    # Retorna o passo encontrado
    if status == MOI.OPTIMAL
        return value.(x)

    elseif status == MOI.INFEASIBLE
        println("LP: Problema infactivel!")
        return zeros(n_total)

    else
        println("LP: Status de termino nao ideal: $status")
        return has_values(model) ? value.(x) : zeros(n_total)

    end

end
