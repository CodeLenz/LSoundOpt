#
# Resolve o subproblema linear inteiro para as variaveis discretas.
#
function LP(c, A, b, gamma_design)

    # Variaveis fisicas de projeto.
    n_design = length(gamma_design)

    # Tamanho total do problema, incluindo folgas.
    n_total = length(c)
    n_slack = n_total - n_design

    model = Model(optimizer_with_attributes(HiGHS.Optimizer,
                  "output_flag" => false,
                  "presolve" => "on",
                  "time_limit" => 60.0))

    @variable(model, x[1:n_total], Int)

    # Delta gamma discreto:
    # gamma=0 -> delta gamma in {0, 1}
    # gamma=1 -> delta gamma in {-1, 0}
    for i = 1:n_design
        lb = -round(Int, gamma_design[i])
        ub =  round(Int, 1 - gamma_design[i])
        set_lower_bound(x[i], lb)
        set_upper_bound(x[i], ub)
    end

    # Folgas inteiras nao-negativas.
    if n_slack > 0
        for i = (n_design + 1):n_total
            set_lower_bound(x[i], 0.0)
        end
    end

    @constraint(model, A * x .<= b)
    @objective(model, Min, c' * x)

    JuMP.optimize!(model)
    status = termination_status(model)

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
