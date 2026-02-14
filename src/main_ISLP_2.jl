
"""
 Otim_ISLP(arquivo::String, freqs::Vector, vA::Vector; verifica_derivada=false)
"""
function Otim_ISLP(arquivo::String, freqs::Vector, vA::Vector; verifica_derivada=false)

    # --- 1. SETUP INICIAL ---
    mshfile, arquivos_saida = Setup_Arquivos(arquivo)
    arquivo_pos, arquivo_pos_freq, arquivo_data_opt, arquivo_γ_ini, arquivo_γ_fin = arquivos_saida

    # Inicializa vetor de pesos para o SPL (sensibilidade em relação à frequência)
    if isempty(vA)
        vA = ones(length(freqs))
    end

    # Número de frequências a considerar no SPL
    nf = length(freqs)
    
    # Leitura da malha de elementos finitos
    nn, coord, ne, connect, materials, nodes_open, velocities, nodes_pressure, pressures, damping, nodes_probe, nodes_target, elements_fixed, values_fixed, centroides = LSound.Parsemsh_Daniele(mshfile)
    
    # Gera a lista de elementos de projeto
    elements_design = setdiff(1:ne, sort!(elements_fixed))

    # Número de elementos de projeto
    nvp = length(elements_design) 

    # Leitura do YAML
    raio_filtro, niter, vf, Past, μ, fatorcv = Le_YAML(replace(mshfile,".msh"=>".yaml"))
    
    # Mapeamento Global -> Local (Necessário para montar A_topo corretamente na linearização das restrições)
    map_global_local = Dict{Int, Int}()
    for (i, ele) in enumerate(elements_design)
        map_global_local[ele] = i
    end

    # Pré-processamento de Vizinhança e Arestas
    neighedge = NeighborEdges(ne, connect, elements_design)
    if !verifica_derivada && raio_filtro > 0
         println("Determinando a vizinhança para um raio de $(raio_filtro)")
         vizinhos, pesos, _, _ = Vizinhanca(ne, centroides, raio_filtro, elements_design)
    end 

    # Variáveis de projeto iniciais (ar)
    println("Inicializando topologia...")
    γ = zeros(ne)
    Fix_γ!(γ, elements_fixed, values_fixed)

    # Grava deistribuição inicial para um arquivo
    writedlm(arquivo_γ_ini, γ)

    # Identifica graus de liberdade livres
    nodes_mask = sort(vcat(nodes_open, nodes_pressure))
    livres = setdiff(collect(1:nn), nodes_mask)
   
    # Export Inicial
    Lgmsh_export_init(arquivo_pos, nn, ne, coord, connect[:,1], connect[:,3:end])
    Lgmsh_export_init(arquivo_pos_freq, nn, ne, coord, connect[:,1], connect[:,3:end])
    Lgmsh_export_element_scalar(arquivo_pos, γ, "Iter 0")
   
    # Sweep Inicial
    MP, _, _, _ = Sweep(nn, ne, coord, connect, γ, fρ, fκ, μ, freqs, livres, velocities, pressures)
    for i=1:nf
        Lgmsh_export_nodal_scalar(arquivo_pos_freq, abs.(MP[:,i]), "Pressure $(freqs[i]) Hz - Ini")
    end

    # Rotina de verificação de derivada 
    if verifica_derivada
       println("--- Modo de Verificação de Derivada ---")

       # Inicializa as variáveis de projeto com valores aleatórias
       γ = rand(ne) 

       # Sweep 
       MP,K,M,C = Sweep(nn,ne,coord,connect,γ,fρ,fκ,μ, freqs,livres,velocities,pressures)

       # Calcula as derivadas analíticas
       dΦ = Derivada(ne,nn,γ,connect,coord,K,M,C,livres,freqs,pressures,fρ, fκ,dfρ,dfκ,μ, nodes_target,MP,elements_design,vA) 

       # Calcula as derivadas numéricas
       dnum = Verifica_derivada(γ,nn,ne,coord,connect,fρ,fκ,μ,freqs,livres,velocities,pressures,nodes_target,elements_design,vA)

       # Erro relativo
       rel = (dΦ.-dnum)./(dnum.+1E-12)

       # Exporta para o gmsh
       Lgmsh_export_element_scalar(arquivo_pos,γ,"γ_debug")
       Lgmsh_export_element_scalar(arquivo_pos,dΦ,"Analitica")
       Lgmsh_export_element_scalar(arquivo_pos,dnum,"Numerica")
       Lgmsh_export_element_scalar(arquivo_pos,rel,"relativa")
       
       # Verificação derivada perímetro

       # Função para usar nas DF
       FPerimiter(γ) =  Perimiter(γ, neighedge, elements_design)

       # Derivadas numéricas
       dPnum = df(γ,FPerimiter,elements_design)
       dP = dPerimiter(ne, γ, neighedge, elements_design)

       # Retorna os valores para verificação no terminal
       return dΦ, dnum, dP, dPnum 

    end

    ################################################################################
    #                         OTIMIZAÇÃO
    ################################################################################

    # Monta o vetor com os volumes de cada elemento da malha
    V = Volumes(ne, connect, coord)

    # Volume limite (restrição de volume)
    Vast = vf * sum(V[elements_design])

    # Históricos
    hist = (V = Float64[], SLP = Float64[], P = Float64[])
    
    # Parâmetros de Controle
    cv_current = fatorcv 
    stagnation_counter = 0

    # Loop externo 
    for iter = 1:niter
        
        # Copia as variáves de projeto atuais 
        γ_start = copy(γ)

        # Atualiza Heurísticas (Warmup / Stagnation)
        cv_current = Update_Heuristics(iter, stagnation_counter, cv_current, nvp)

        # Volume atual do projeto
        volume_atual = sum(γ[elements_design] .* V[elements_design])

        # Sweep 
        MP, K, M, C = Sweep(nn, ne, coord, connect, γ, fρ, fκ, μ, freqs, livres, velocities, pressures)

        # Objetivo e perímetro 
        objetivo = Objetivo(MP, nodes_target, vA)
        perimetro = Perimiter(γ, neighedge, elements_design)

        # Armazena histórico
        push!(hist.V, volume_atual);
        push!(hist.SLP, objetivo);
        push!(hist.P, perimetro)
        
        # Exporta o campo de pressões (apenas 1a freq para economizar disco)
        Lgmsh_export_nodal_scalar(arquivo_pos_freq, abs.(MP[:,1]), "Press $(freqs[1]) Hz - It $iter")

        # Mostra um resumo na tela, para acompanhamento
        @printf("Iter %3d | SPL: %.4f | Perim: %.4f | Vol: %.4e (T: %.4e) | cv: %.4f\n", 
                iter, objetivo, perimetro, volume_atual, Vast, cv_current)

        # Calcula a derivada do objetivo
        dΦ = Derivada(ne, nn, γ, connect, coord, K, M, C, livres, freqs, pressures, fρ, fκ, dfρ, dfκ, μ, nodes_target, MP, elements_design, vA) 

        # Exporta para o gmsh
        Lgmsh_export_element_scalar(arquivo_pos, dΦ, "dΦ")

        # Se o filtro estiver ativado, filtramos e usamos a sensibilidade filtrada
        if raio_filtro > 0 
            dΦ_filt = Filtro(vizinhos, pesos, dΦ, elements_design)
            c = dΦ_filt[elements_design] 
        else
            c = dΦ[elements_design] 
        end

        # Linearização das restrições
        # Linearização das restrições (Volume, Perímetro, Topologia + Slacks)
        A_global, b_global = Lineariza_Restricoes(V, elements_design, Vast, volume_atual, perimetro, Past, ne, γ, neighedge, map_global_local)

       
        # Empacotando dados FEM necessários para o check físico
        fem_data = (nn, ne, coord, connect, fρ, fκ, μ, freqs, livres, velocities, pressures, nodes_target, vA)
        

        # Trust Region Loop (Otimização do Passo)
        γ_new, cv_current, step_accepted = Trust_Region_Loop(c, A_global, b_global, γ, elements_design, cv_current, objetivo, Past, neighedge, fem_data)

        # Verificação do passo
        if !step_accepted
            println("Trust Region falhou. Terminando.")
            break
        end
        
        # Atualiza variáveis
        γ .= γ_new

        # Exporta para o gmsh
        Lgmsh_export_element_scalar(arquivo_pos, γ, "Iter $iter")

        # Verifica se estagnou
        if sum(abs.(γ - γ_start)) == 0

            # Incrementa o contador de stagnação
            stagnation_counter += 1

            println("--- AVISO: Estagnação detectada ($stagnation_counter/4).")

            if stagnation_counter >= 4
               println("HARD STOP.")
               break
            end
        else
            stagnation_counter = 0
        end

    end # Fim Loop

    # FINALIZAÇÃO 
    println("Finalizando...")

    # Grava a configuração otimizda
    writedlm(arquivo_γ_fin, γ)

    # Salva histórico da otimização 
    Salva_Historico(arquivo_data_opt, hist, raio_filtro)
    
    # Retorna o histórico para o terminal
    return hist.V, hist.SLP, hist.P
end

