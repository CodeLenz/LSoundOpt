#
# Rotina principal
#
"""
 Otim(meshfile::String,freqs=[])

 Basic input:

 meshfile -> arquivo de entrada (.msh)

 freqs -> vetor com as frquências (em Hz)

 vA -> Vetor com as sensibilidades para cada frequência

 verifica_derivada -> bool 

 Inputs -> freqs, a vector with the frequencies to sweep

"""
function Otim_ISLP(arquivo::String,freqs::Vector, vA::Vector;verifica_derivada=false)

   # Se o arquivo for um .geo, geramos um .msh utilizando a biblioteca
   # do gmsh
   if occursin(".geo",arquivo)
       
      # Gera a malha
      gmsh.initialize()
      gmsh.open(arquivo)
      gmsh.model.mesh.generate(2)
       
      # Cria o mesmo nome, mas com .msh
      mshfile = replace(arquivo,".geo"=>".msh")

      # Cria o .msh
      gmsh.write(mshfile)
      
   else 

      # Assumimos que já passaram o .msh (seria bom testar...)
      mshfile = arquivo

   end

    # Define os nomes dos arquivos de entrada (yaml) 
    arquivo_yaml = replace(mshfile,".msh"=>".yaml")

    # Nomes de arquivos que serão gravados no direorio raiz
    nomebase = basename(mshfile)

    # Arquivo .pos
    arquivo_pos  = replace(nomebase,".msh"=>".pos")

    # Define o nome dos arquivos contendo a distribuição inicial e final (otimizada)
    # das variáveis de projeto.
    # Esses arquivos serão utilizados posteriormente em 
    # Processa_FRF
    arquivo_γ_ini = replace(nomebase,".msh"=>"_γ_ini.dat")
    arquivo_γ_fin = replace(nomebase,".msh"=>"_γ_opt.dat")

    # Verificamos se existem frequências sendo informadas
    isempty(freqs) && error("Analise Harmonica:: freqs deve ser um vetor não vazio")

    # Evita passar as frequências como algo diferente de um Vetor de floats
    isa(freqs,Vector{Float64}) || error("freqs deve ser um vetor de floats")
    
    # Verifica se os arquivos de entrada existem
    isfile(mshfile) || error("Otim:: arquivo de entrada $mshfile não existe")

    # Número de frequências
    nf = length(freqs)

    # Verifica se A foi definido
    if isempty(vA)
       vA = ones(nf)
    else 
       length(vA)==nf || error("Otim:: dimensão de vA deve ser nf")
    end

    # Arquivo .yaml
    isfile(arquivo_yaml) || error("Otim:: arquivo de entrada $(arquivo_yaml) não existe")

    # Le dados da malha
    nn, coord, ne, connect, materials, nodes_open, velocities, nodes_pressure, pressures, damping, nodes_probe, nodes_target, elements_fixed, values_fixed, centroides = LSound.Parsemsh_Daniele(mshfile)

    # Lista com os elementos que são de projeto
    elements_design = setdiff(1:ne,sort!(elements_fixed))

    # Número de variáveis de projeto 
    nvp = length(elements_design)

    # Le os dados do arquivo yaml
    raio_filtro, niter, ϵ1, ϵ2,  vf, Past, fatorcv, μ = Le_YAML(arquivo_yaml)
     
   
    # Agora que queremos otimizar o SPL, vamos precisar OBRIGATÓRIAMENTE de nodes_target,
    # que vai funcionar como nodes_probe aqui
    isempty(nodes_target) && error("Otim:: nodes_target deve ter ao menos um nó informado")

    # Vamos colocar nodes_target em ordem crescente
    sort!(nodes_target)

    # Precisamos de um material
    isempty(materials) && error("Analise:: at least one material is necessary")

    # Não precisamos do centróide e de vizinhança se for para verificar a derivada
    if !verifica_derivada
    
    
      # TODO 
         # Ver cálculo automático de raio se raio_filtro for nulo
         #
         
         # Obtém os vizinhos de cada elemento da malha
         println("Determinando a vizinhança para um raio de $(raio_filtro)")
         vizinhos, pesos = Vizinhanca(ne,centroides,raio_filtro,elements_design)

    end # !verifica_derivada

    # Em teste...vizinhos de arestas para restrição de perímetro
    neighedge = NeighborEdges(ne,connect,elements_design)
        
    # Vamos inicializar o vetor de variáveis de projeto.
    # γ = 0 --> ar
    # γ = 1 --> sólido
    #
    println("Inicializando o vetor de variáveis de projeto")
    γ = zeros(ne)
   
    # Fixa os valores prescritos de densidade relativa
    Fix_γ!(γ,elements_fixed,values_fixed)
    
    # Grava o arquivo com a distribuição inicial de densidades
    writedlm(arquivo_γ_ini,γ)

    # Concatena nodes_open e nodes_pressure
    nodes_mask = sort(vcat(nodes_open,nodes_pressure))
    
    # Posições que não precisam ser calculadas no sistema de equações 
    livres = setdiff(collect(1:nn),nodes_mask)
   
    # Inicializa um arquivo de pós-processamento do gmsh
    etype = connect[:,1]
    Lgmsh_export_init(arquivo_pos,nn,ne,coord,etype,connect[:,3:end])
    
    # Grava a topologia inicial 
    Lgmsh_export_element_scalar(arquivo_pos,γ,"Iter 0")
   
    # Sweep na topologia inicial, para comparação com a otimizada
    MP,_ =  Sweep(nn,ne,coord,connect,γ,fρ,fκ,μ, freqs,livres,velocities,pressures)

    # Exporta por frequência
    for i=1:nf

      # frequência
      f = freqs[i]

      # Exporta
      Lgmsh_export_nodal_scalar(arquivo_pos,abs.(MP[:,i]),"Pressure in $f Hz [abs] - initial topology")

    end

    # Verifica derivadas
    if verifica_derivada

      # Vamos inicializar com algo aleatório
      γ = rand(ne)

      # Derivada utilizando o procedimento analítico
      MP,K,M,C = Sweep(nn,ne,coord,connect,γ,fρ,fκ,μ, freqs,livres,velocities,pressures)

      # Calcula a derivada da função objetivo em relação ao vetor γ
      dΦ = Derivada(ne,nn,γ,connect,coord,K,M,C,livres,freqs,pressures,dfρ,dfκ,μ, nodes_target,MP,elements_design,vA) 

      println("Verificando as derivadas utilizando diferenças finitas centrais...")
      println("O número efetivo de variáveis de projeto é ", length(elements_design))

      # Derivada numérica
      dnum = Verifica_derivada(γ,nn,ne,coord,connect,fρ,fκ,μ,freqs,livres,velocities,pressures,nodes_target,elements_design,vA)

      # Relativo, evitando divisão por zero
      rel = (dΦ.-dnum)./(dnum.+1E-12)
      
      # Exporta para a visualização no Gmsh
      Lgmsh_export_element_scalar(arquivo_pos,γ,"γ")
      Lgmsh_export_element_scalar(arquivo_pos,dΦ,"Analitica")
 
      Lgmsh_export_element_scalar(arquivo_pos,dnum,"Numerica")
      Lgmsh_export_element_scalar(arquivo_pos,rel,"relativa")

      #
      # Vamos validar a derivada do perímetro também 
      #
      FPerimiter(γ) =  Perimiter(γ, neighedge, elements_design)
      dPnum = df(γ,FPerimiter,elements_design)
      dP = dPerimiter(ne, γ, neighedge, elements_design)

      # Retorna as derivadas
      return dΦ, dnum, dP, dPnum 

    end

    ########################################################################################################
    ################################ Começo do loop principal de otimização topológica #####################
    ########################################################################################################

    # Volume of each element
    V = Volumes(ne,connect,coord)

    # Total volume sem a parametrização 
    # mas só dos elementos de projeto
    volume_full_projeto = sum(V[elements_design])

    # Target volume
    Vast = vf*volume_full_projeto

    # Monitora o histórico de volume, do perímetro e 
    # do SPL ao longo das  iterações 
    historico_V   = Float64[] #zeros(niter)
    historico_SLP = Float64[] #zeros(niter)
    historico_P   = Float64[] #zeros(niter)

    # INICIALIZA O FATOR DE BUSCA ADAPTATIVO (TRUST REGION)
    cv_current = fatorcv 

    #############################  Main loop ###########################
    for iter = 1:niter

        # Lista de elementos de projeto com ar e com sólido
        elements_air = Int64[]
        elements_solid = Int64[]

        # Vetores com a dimensão do número de variáveis de projeto 
        one_air = Float64[]
        one_solid = Float64[]

        # Elementos de projeto, separando por ar ou sólido
        for ele in elements_design
            if γ[ele]<0.5 
               push!(elements_air,ele)
               push!(one_air,1.0)
               push!(one_solid,0)
            else
               push!(elements_solid,ele)
               push!(one_air,0)
               push!(one_solid,1.0)
            end
        end

        # Volume atual da estrutura
        volume_atual = sum(γ[elements_design].*V[elements_design])

        # Armazena o volume no histório de volumes
        push!(historico_V , volume_atual)

        # Faz o sweep. A matriz MP tem dimensão nn × nf, ou seja, 
        # cada coluna é o vetor P para uma frequência de excitação
        MP,K,M,C =  Sweep(nn,ne,coord,connect,γ,fρ,fκ,μ,freqs,livres,velocities,pressures)

        # Calcula o SPL para esta iteração 
        objetivo = Objetivo(MP,nodes_target,vA)

        # Armazena no historico de SPL
        push!(historico_SLP, objetivo)

        # Perímetro 
        perimetro = Perimiter(γ, neighedge, elements_design)

        # Exporta por frequência
        for i=1:nf

         # frequência
         f = freqs[i]

         # Exporta
         Lgmsh_export_nodal_scalar(arquivo_pos,abs.(MP[:,i]),"Pressure in $f Hz [abs]")

        end

        # Guarda para o histórico
        push!(historico_P, perimetro)

        println("Iteração        ", iter)
        println("Objetivo        ", objetivo)
        println("Perimetro       ", perimetro)
        println("Past            ", Past)
        println("Volume atual    ", volume_atual)
        println("Volume target   ", Vast)
        println()

        # Calcula a derivada da função objetivo em relação ao vetor γ
        # somente nas posições de projeto
        dΦ = Derivada(ne,nn,γ,connect,coord,K,M,C,livres,freqs,pressures,dfρ,dfκ,μ,nodes_target,MP,elements_design,vA) 
  
        # Visualiza as ESEDS...
        Lgmsh_export_element_scalar(arquivo_pos,dΦ,"dΦ")  

        # Filtra o gradiente do objetivo (MOVIDO PARA ANTES DO TRUST REGION LOOP)
        dΦf =  Filtro(vizinhos,pesos,dΦ,elements_design)
        Lgmsh_export_element_scalar(arquivo_pos,dΦf,"Filtrada")

        # Vetor de coeficientes da função objetivo linearizada
        c = dΦf[elements_design]

        # Lado direito da restrição de volume linearizada em Δγ
        ΔV = Vast - volume_atual
       
        # Primeiro limite de restrição...volume
        b = [ΔV]

        # Mostra a linearização 
        println("Volume linearizado   ", b[1])

        # Derivada do volume
        A = vcat(V[elements_design]')

        #
        # Lógica para relaxar a restrição de perímetro
        #
         
        # Parâmetros para comparação 
        if  perimetro > 0 && Past>0

            # Variação sem a relaxação
            ΔP = Past - perimetro

            # Limite da restrição de perímetro linearizada e relaxada
            b = vcat(b, ΔP)
            
            # Mostra a linearização 
            println("Perimetro linearizado   ", ΔP)

            # Derivada do perímetro
            dP = dPerimiter(ne, γ, neighedge, elements_design)
             
            # Visualiza a derivada do perímetro
            Lgmsh_export_element_scalar(arquivo_pos,dP,"dP")  

            # Adiciona a linha em A
            A = vcat(A,transpose(dP[elements_design]))
            
        end
       
        # -----------------------------------------------------------------
        # STRATEGY B: TRUST REGION (ADAPTIVE MOVE LIMITS) SUB-LOOP
        # -----------------------------------------------------------------
        step_accepted = false
        
        # Número de rejeições 
        n_reject = 0

        flag_podrao = false

        while !step_accepted
            
            # Copias temporarias de A e b para as restricoes variaveis
            A_temp = copy(A)
            b_temp = copy(b)

            # Restrição de variação de elementos com ar (usando cv_current)
            if !isempty(elements_air)
                g_air = ceil(cv_current * length(elements_design))
                b_temp = vcat(b_temp, g_air)
                A_temp = vcat(A_temp, transpose(one_air))
            end

            # Restrição de variação de elementos sólidos (usando cv_current)
            if !isempty(elements_solid)
                g_solid  = ceil(cv_current * length(elements_design)) 
                b_temp = vcat(b_temp, g_solid)
                A_temp = vcat(A_temp, -transpose(one_solid))
            end

            # Resolve o ILP com as restrições adaptativas
            Δγ = LP(c, A_temp, b_temp, γ[elements_design])

            # Calcula a redução predita (Linearizada).
            # Como é minimização, a redução esperada é negativa do produto interno.
            pred_reduction = -sum(c .* Δγ)

            # Se a redução predita for muito proxima de 0, o algoritmo convergiu
            if abs(pred_reduction) < 1e-8
                println("Convergência local atingida (nenhuma melhora predita no LP).")
                step_accepted = true
                break
            end

            # Cria um passo de teste (Trial Step)
            γ_trial = copy(γ)
            γ_trial[elements_design] .+= Δγ
            
            # Garante que os valores trial sejam puramente binários 0 ou 1
            for ele in elements_design
                γ_trial[ele] = round(γ_trial[ele])
            end

            # Avalia a física não-linear real do passo de teste
            MP_trial, _, _, _ = Sweep(nn,ne,coord,connect,γ_trial,fρ,fκ,μ,freqs,livres,velocities,pressures)
            obj_trial = Objetivo(MP_trial,nodes_target,vA)

            # Calcula a redução real do SPL
            actual_reduction = objetivo - obj_trial
            
            # Razão de Trust-Region (Real vs Predito)
            R = actual_reduction / pred_reduction

            println("--- Trust Region Eval: Pred = $(round(pred_reduction, digits=4)), Act = $(round(actual_reduction, digits=4)), R = $(round(R, digits=2)), cv = $(round(cv_current, digits=4))")

            if R < 0.0 
                # REJEITADO: O SPL piorou! A aproximação linear falhou.
                # Encolhe o fator de busca e repete o LP na mesma iteração.
                println("    -> Passo REJEITADO. Reduzindo limites de busca (cv).")
                cv_current = max(0.002, cv_current * 0.5) 

                # Atuaiza o número de rejeições
                n_reject +=1
                
                # Se o número de rejeições for maior do que 10, termina
                if n_reject>= 10
                    flag_podrao = true
                    break
                end

            else
                # ACEITO: O SPL melhorou.
                step_accepted = true
                γ .= γ_trial  # Atualiza o design principal

                if R > 0.7 
                    # Previsão excelente: Aumenta o passo para a próxima iteração
                    println("    -> Passo ACEITO. Previsão excelente, expandindo cv para proxima iteracao.")
                    cv_current = min(0.5, cv_current * 1.2)
                elseif R > 0.3
                    println("    -> Passo ACEITO. Previsão boa, mantendo cv.")
                else
                    # Previsão ruim mas aceitável: Encolhe o passo para a próxima iteração
                    println("    -> Passo ACEITO. Previsão pobre, reduzindo cv para proxima iteracao.")
                    cv_current = max(0.002, cv_current * 0.5)
                end
            end
        end # Fim do sub-loop do Trust Region

        # Se não foi possível melhorar, terminamos
        if flag_podrao
           break
        end  

        # Grava a topologia para visualização 
        Lgmsh_export_element_scalar(arquivo_pos,γ,"Iter $iter")

    end # iterações externas

    println("Final da otimização, executando a análise SWEEP na topologia otimizada")

   # Roda o sweep na topologia otimizada e exporta para visualização 
   MP,_ =  Sweep(nn,ne,coord,connect,γ,fρ,fκ,μ,freqs,livres,velocities,pressures)

   # Exporta por frequência
   for i=1:nf

      # frequência
      f = freqs[i]

      # Exporta
      Lgmsh_export_nodal_scalar(arquivo_pos,abs.(MP[:,i]),"Pressure in $f Hz [abs]")

   end

   # Grava um arquivo com os γ finais
   writedlm(arquivo_γ_fin,γ)

   # Retorna o histórico de volume e também o da função objetivo 
   println("Retornando históricos de V e de SLP")
   return historico_V, historico_SLP, historico_P

end # main_otim