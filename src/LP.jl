#
# n => número de variáveis de projeto (design vars)
# O vetor 'c' e as colunas de 'A' têm tamanho n_total = n_design + n_slack
#
function LP(c, A, b, γ_design)

    # Verifica quantas variáveis de projeto temos (apenas as físicas)
    n_design = length(γ_design)
    
    # Verifica o tamanho total do problema (incluindo slacks)
    n_total = length(c)
    
    # Número de slacks
    n_slack = n_total - n_design

    # -----------------------------------------------------------
    # SOLVER SETUP
    # O problema é um MILP (Mixed-Integer Linear Programming).
    # Não é necessário usar Alpine. O HiGHS é excelente para isso.
    # -----------------------------------------------------------
    model = Model(optimizer_with_attributes(HiGHS.Optimizer, 
                  "output_flag" => false, 
                  "presolve" => "on",
                  "time_limit" => 60.0)) # Timeout de segurança

    # -----------------------------------------------------------
    # DEFINIÇÃO DE VARIÁVEIS E LIMITES
    # -----------------------------------------------------------
    
    # Cria as variáveis de decisão do modelo (x[1...n_total])
    # Vamos definir os limites individualmente
    @variable(model, x[1:n_total], Int)

    # 1. Variáveis de Projeto (Delta Gamma): Indices 1 até n_design
    # Devem respeitar a lógica binária: 
    # Se γ=0 -> Δγ ∈ {0, 1}   (Limites: 0, 1)
    # Se γ=1 -> Δγ ∈ {-1, 0}  (Limites: -1, 0)
    for i = 1:n_design
        lb = -round(Int, γ_design[i])      # 0 ou -1
        ub =  round(Int, 1 - γ_design[i])  # 1 ou 0
        set_lower_bound(x[i], lb)
        set_upper_bound(x[i], ub)
    end

    # 2. Variáveis de Folga (Xi): Indices n_design+1 até n_total
    # Devem ser apenas não-negativas: ξ ∈ {0, 1, 2, ...}
    # (Como penalizamos com Beta alto, elas tenderão a zero)
    if n_slack > 0
        for i = (n_design + 1):n_total
            set_lower_bound(x[i], 0.0)
            # Não definimos upper bound (ou definimos um valor alto seguro se necessário)
            # set_upper_bound(x[i], 10.0) 
        end
    end

    # -----------------------------------------------------------
    # RESTRIÇÕES E OBJETIVO
    # -----------------------------------------------------------

    # Monta todas as restrições lineares (Globais + Topológicas)
    # A matriz A já contém as colunas para Delta Gamma e para Xi
    @constraint(model, A * x .<= b)

    # Monta a função objetivo (Custo Físico + Penalidade Beta)
    @objective(model, Min, c' * x)

    # -----------------------------------------------------------
    # SOLUÇÃO
    # -----------------------------------------------------------
    # optimize!(model)
    # Redireciona stdout para evitar poluição no terminal se o solver for verboso
    # redirect_stdout((()->optimize!(model)), open(devnull, "w"))
    optimize!(model)

    status = termination_status(model)

    if status == MOI.OPTIMAL
        Δxopt = value.(x)
        return Δxopt
    elseif status == MOI.INFEASIBLE
        println("LP: Problema Infactível!")
        return zeros(n_total) # Retorna nulo para não quebrar, o loop principal deve tratar
    else
        println("LP: Status de término não ideal: $status")
        # Tenta recuperar o melhor valor viável encontrado
        if has_values(model)
            return value.(x)
        else
            return zeros(n_total)
        end
    end

end

#
# n => número de variáveis de projeto
#
function LP_anterior(c, A, b, γ)

   #
   # Cria o modelo vazio e associa a um otimizador
   # O Alpine utiliza mais de um pacote de otimização, dependendo 
   # da etapa que ele está realizando. O mip_solver é bem crítico para 
   # a primeira etapa, em que é contínua. O Gurobi é uma boa opção, mas 
   # precisamos instalar a licença no computador. Uma outra opção é utilizar
   # o Ipopt ou  HigHS.
   # O Cbc é o otimizador para a etapa discreta (branch and bound), que é 
   # realizada depois da etapa contínua.
   ipopt  = optimizer_with_attributes(Ipopt.Optimizer, "print_level" => 0)
   gurobi = optimizer_with_attributes(Gurobi.Optimizer, "output_flag" => false)
   highs  = optimizer_with_attributes(HiGHS.Optimizer, "output_flag" => false)
   cbc    = optimizer_with_attributes(Cbc.Optimizer, "logLevel" => 0)
    
   # minlp_solver é o solver para o problema binário (0/1)
   # nlp_solver é o solver para o problema contínuo do LP
   # mip_solver é o utilizado na primeira etapa 
   model = Model(
      optimizer_with_attributes(
         Alpine.Optimizer,
         "minlp_solver" => highs, #<- para binario. Aqui também dá para usar o cbc
         #"nlp_solver" => ipopt,  #<- para contínuo
         #"mip_solver" => gurobi,
         "mip_solver" => highs,
      ),
   )

   # Vamos descobrir o número de variáveis de projeto e de restrições 
   nc = length(c)
   mb = length(b)
   ma,na = size(A)

   # E testar a consistência dos dados
   nc==na || error("LP: dimensões inconsistentes")
   mb==ma || error("LP: dimensões inconsistentes")

   # Cria o vetor com as restrições laterais para cada variável
   Δxi = zeros(Int,na)
   Δxs = zeros(Int,na)

   for i in LinearIndices(Δxi)
       Δxi[i] = -round(Int,γ[i])
       Δxs[i] =  round(Int,1-γ[i])
   end

   # Cria um vetor de variáveis de projeto
   @variable(model, Δxi[i] <= Δx[i=1:na] <= Δxs[i], Int)

   # Monta todas as restrições ao mesmo tempo 
   @constraint(model, A * Δx .<= b)

   # Monta a função objetivo
   @objective(model, Min, c' * Δx)

   # Resolve o problema 
   #optimize!(model)
   redirect_stdout((()->optimize!(model)),open("nul", "w"))

   # Valor do objetivo
   # objective_value(model)

   # Vetor de variáveis de projeto no ponto de ótimo
   Δxopt = value(Δx)

   # Vamos testar as restrições
   @show [A*Δxopt b]
   @show termination_status(model)

   # @show c
   # @show A
   # @show xopt

   # Retorna o vetor de variáveis de projeto no ponto de ótimo
   return Δxopt

end


