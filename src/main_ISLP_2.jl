

# Supõe-se que LSound e outras dependências já estejam carregadas
# e que a função LP suporte o aumento de variáveis.

"""
 Otim_ISLP(arquivo::String, freqs::Vector, vA::Vector; verifica_derivada=false)
 
 Entradas:
 - arquivo: Caminho para o arquivo .geo ou .msh
 - freqs: Vetor de frequências para análise (Recomendado step=1Hz para ressonadores)
 - vA: Vetor de pesos para as frequências (geralmente ones(nf))
"""
function Otim_ISLP(arquivo::String, freqs::Vector, vA::Vector; verifica_derivada=false)

    # ---------------------------------------------------------------------------
    # 1. PREPARAÇÃO DA MALHA E ARQUIVOS
    # ---------------------------------------------------------------------------
    if occursin(".geo",arquivo)
       gmsh.initialize()
       gmsh.open(arquivo)
       gmsh.model.mesh.generate(2)
       mshfile = replace(arquivo,".geo"=>".msh")
       gmsh.write(mshfile)
    else 
       mshfile = arquivo
    end

    # Define os nomes dos arquivos de saída
    arquivo_yaml = replace(mshfile,".msh"=>".yaml")
    nomebase = basename(mshfile)
    arquivo_pos       = replace(nomebase,".msh"=>".pos")
    arquivo_pos_freq  = replace(nomebase,".msh"=>"_freq.pos")
    arquivo_data_opt  = replace(nomebase,".msh"=>".data")

    arquivo_γ_ini = replace(nomebase,".msh"=>"_γ_ini.dat")
    arquivo_γ_fin = replace(nomebase,".msh"=>"_γ_opt.dat")

    # Verificações de entrada básicas
    isempty(freqs) && error("Analise Harmonica:: freqs deve ser um vetor não vazio")
    isa(freqs,Vector{Float64}) || error("freqs deve ser um vetor de floats")
    isfile(mshfile) || error("Otim:: arquivo de entrada $mshfile não existe")
    nf = length(freqs)

    if isempty(vA)
       vA = ones(nf)
    else 
       length(vA)==nf || error("Otim:: dimensão de vA deve ser nf")
    end

    isfile(arquivo_yaml) || error("Otim:: arquivo de entrada $(arquivo_yaml) não existe")

    # Leitura da malha e parâmetros do YAML
    nn, coord, ne, connect, materials, nodes_open, velocities, nodes_pressure, pressures, damping, nodes_probe, nodes_target, elements_fixed, values_fixed, centroides = LSound.Parsemsh_Daniele(mshfile)
    elements_design = setdiff(1:ne,sort!(elements_fixed))
    nvp = length(elements_design) 
    
    # [NOVO] Mapeamento Global -> Local para montagem eficiente das restrições
    map_global_local = Dict{Int, Int}()
    for (i, ele) in enumerate(elements_design)
        map_global_local[ele] = i
    end

    raio_filtro, niter, ϵ1, ϵ2,  vf, Past, fatorcv, μ = Le_YAML(arquivo_yaml)
     
    isempty(nodes_target) && error("Otim:: nodes_target deve ter ao menos um nó informado")
    sort!(nodes_target)
    isempty(materials) && error("Analise:: at least one material is necessary")

    # Vizinhos de aresta, para cálculo do Perímetro e Restrições Topológicas
    neighedge = NeighborEdges(ne,connect,elements_design)
        
    # ---------------------------------------------------------------------------
    # 2. INICIALIZAÇÃO
    # ---------------------------------------------------------------------------
    println("Inicializando o vetor de variáveis de projeto com AR (Zeros).")
    γ = zeros(ne)
    Fix_γ!(γ,elements_fixed,values_fixed)
    writedlm(arquivo_γ_ini,γ)

    
    nodes_mask = sort(vcat(nodes_open,nodes_pressure))
    livres = setdiff(collect(1:nn),nodes_mask)
   
    # Inicializa os arquivos de visualização do gmsh (.pos)
    etype = connect[:,1]
    Lgmsh_export_init(arquivo_pos,nn,ne,coord,etype,connect[:,3:end])
    Lgmsh_export_init(arquivo_pos_freq,nn,ne,coord,etype,connect[:,3:end])
    Lgmsh_export_element_scalar(arquivo_pos,γ,"Iter 0")
   
    # Sweep Inicial
    MP,_ =  Sweep(nn,ne,coord,connect,γ,fρ,fκ,μ, freqs,livres,velocities,pressures)
    for i=1:nf
      f = freqs[i]
      Lgmsh_export_nodal_scalar(arquivo_pos_freq,abs.(MP[:,i]),"Pressure in $f Hz [abs] - initial topology")
    end

    # (Bloco de verificação de derivada omitido para brevidade, manter igual ao original se necessário)

    # --------------------------------------------------------------------------
    # 3. LOOP PRINCIPAL DE OTIMIZAÇÃO
    # --------------------------------------------------------------------------
    V = Volumes(ne,connect,coord)
    volume_full_projeto = sum(V[elements_design])
    Vast = vf*volume_full_projeto

    historico_V   = Float64[] 
    historico_SLP = Float64[] 
    historico_P   = Float64[] 

    # --- SETUP DOS PARÂMETROS DE CONTROLE ---
    cv_min = 1.0001 / nvp  
    n_warmup = 5 
    cv_warmup_cap = 0.05 
    
    cv_current = fatorcv 
    stagnation_counter = 0

    for iter = 1:niter
        
        γ_start = copy(γ)

        # --- LÓGICA DO PROGRESSIVE SHAKE DOWN ---
        if stagnation_counter > 0
            if stagnation_counter == 1
                println("--- STAGNATION RECOVERY (Nível 1): Resetando cv para 5%.")
                cv_current = 0.05
            elseif stagnation_counter == 2
                println("--- STAGNATION RECOVERY (Nível 2): Aumentando a aposta. cv = 15%.")
                cv_current = 0.15
            elseif stagnation_counter >= 3
                println("--- STAGNATION RECOVERY (Nível 3): O Martelo. cv = 30%.")
                cv_current = 0.30
            end
        elseif cv_current < (2.0 / nvp) 
            println("--- RECUPERAÇÃO PREVENTIVA: cv muito baixo. Resetando para 5%.")
            cv_current = 0.05
        end

        # --- WARM-UP LOGIC ---
        is_warmup = iter <= n_warmup
        if is_warmup && stagnation_counter == 0
            if cv_current > cv_warmup_cap
                cv_current = cv_warmup_cap
            end
        end

        # Identificação de elementos Ar/Sólido (para restrições lineares)
        elements_air = Int64[]
        elements_solid = Int64[]
        # Vetores de coeficientes estendidos (size = 2*nvp) para acomodar variáveis de folga
        # As restrições de movimento agem apenas sobre Delta_Gamma, coef 0 para Xi
        one_air_ext = zeros(2*nvp) 
        one_solid_ext = zeros(2*nvp)

        for (i, ele) in enumerate(elements_design)
            if γ[ele]<0.5 
               push!(elements_air,ele)
               one_air_ext[i] = 1.0   # Apenas na parte Delta Gamma
               # one_solid_ext[i] = 0.0 (já é zero)
            else
               push!(elements_solid,ele)
               # one_air_ext[i] = 0.0
               one_solid_ext[i] = 1.0 # Apenas na parte Delta Gamma
            end
        end

        volume_atual = sum(γ[elements_design].*V[elements_design])
        push!(historico_V , volume_atual)

        # --- ANÁLISE FEM ---
        MP,K,M,C =  Sweep(nn,ne,coord,connect,γ,fρ,fκ,μ,freqs,livres,velocities,pressures)
        objetivo = ObjetivoP(MP,nodes_target,vA)
        push!(historico_SLP, objetivo)

        perimetro = Perimiter(γ, neighedge, elements_design)
        push!(historico_P, perimetro)
        
        # Exporta resultados parciais
        for i=1:nf
            f = freqs[i]
            # Otimização visual: exportar apenas a primeira freq para economizar disco se nf for grande
            if i==1; Lgmsh_export_nodal_scalar(arquivo_pos_freq,abs.(MP[:,i]),"Pressure $f Hz - Iter $iter"); end
        end

        println("Iteração          ", iter)
        println("Objetivo (SPL)    ", objetivo)
        println("Perimetro         ", perimetro)
        println("Past (Limit)      ", Past)
        println("Volume atual      ", volume_atual)
        println()

        # --- SENSIBILIDADE ---
        dΦ = DerivadaP(ne,nn,γ,connect,coord,K,M,C,livres,freqs,pressures,fρ,fκ,dfρ,dfκ,μ,nodes_target,MP,elements_design,vA)
        Lgmsh_export_element_scalar(arquivo_pos,dΦ,"dΦ")  

        # Sensibilidade "crua" para o ILP (apenas parte Delta Gamma)
        c_phys = dΦ[elements_design] 

        # [NOVO] Montagem do Vetor de Custo Estendido (Com Penalidade para Folgas)
        # c_aug = [c_phys; c_slack]
        # Penalidade Beta alta o suficiente para ser evitada, mas finita
        beta_penalty = 100.0 * maximum(abs.(c_phys)) + 10.0
        c_aug = vcat(c_phys, fill(beta_penalty, nvp))

        # --- RESTRIÇÕES GLOBAIS (VOL & PERIM) ---
        # 1. Volume
        ΔV = Vast - volume_atual
        b_global = [ΔV]
        # A_vol deve ter zeros para as variáveis de folga
        A_vol = vcat(V[elements_design]')
        A_vol_aug = hcat(A_vol, zeros(1, nvp)) # Padding zeros

        A_global = A_vol_aug

        # 2. Perímetro (se ativo)
        if  perimetro > 0 && Past>0
            ΔP = Past - perimetro
            b_global = vcat(b_global, ΔP)
            dP = dPerimiter(ne, γ, neighedge, elements_design)
            A_per = transpose(dP[elements_design])
            A_per_aug = hcat(A_per, zeros(1, nvp)) # Padding zeros
            A_global = vcat(A_global, A_per_aug)
        end

        # --- [NOVO] RESTRIÇÕES TOPOLÓGICAS LOCAIS ---
        # Construção da matriz esparsa para No-Island e No-Hole
        # Cada elemento gera 2 restrições
        I_topo = Int[]
        J_topo = Int[]
        V_topo = Float64[]
        b_topo = Float64[]
        
        row_idx = 1
        
        for (idx_e, e_global) in enumerate(elements_design)
            vizinhos = neighedge[e_global]
            γ_e = γ[e_global]
            
            # Soma dos vizinhos (inclui fixos para cálculo do RHS)
            sum_γ_viz = 0.0
            for v in vizinhos
                sum_γ_viz += γ[v]
            end
            num_viz = length(vizinhos)

            # --- RESTRIÇÃO 1: NO ISLAND ---
            # Δγ_e - Σ Δγ_j - ξ_e <= Σ γ_j^k - γ_e^k
            
            # Coeficiente Delta_e (+1)
            push!(I_topo, row_idx); push!(J_topo, idx_e); push!(V_topo, 1.0)
            
            # Coeficiente Slack Xi_e (-1) -> Índice é nvp + idx_e
            push!(I_topo, row_idx); push!(J_topo, nvp + idx_e); push!(V_topo, -1.0)
            
            # Coeficientes Vizinhos (-1)
            for v in vizinhos
                if haskey(map_global_local, v)
                    idx_viz = map_global_local[v]
                    push!(I_topo, row_idx); push!(J_topo, idx_viz); push!(V_topo, -1.0)
                end
            end
            
            # RHS
            push!(b_topo, sum_γ_viz - γ_e)
            row_idx += 1
            
            # --- RESTRIÇÃO 2: NO HOLE ---
            # Σ Δγ_j - Δγ_e - ξ_e <= (N_viz - 1) - (Σ γ_j^k - γ_e^k)
            
            # Coeficiente Delta_e (-1)
            push!(I_topo, row_idx); push!(J_topo, idx_e); push!(V_topo, -1.0)
            
            # Coeficiente Slack Xi_e (-1)
            push!(I_topo, row_idx); push!(J_topo, nvp + idx_e); push!(V_topo, -1.0)
            
            # Coeficientes Vizinhos (+1)
            for v in vizinhos
                if haskey(map_global_local, v)
                    idx_viz = map_global_local[v]
                    push!(I_topo, row_idx); push!(J_topo, idx_viz); push!(V_topo, 1.0)
                end
            end
            
            # RHS
            push!(b_topo, (num_viz - 1) - (sum_γ_viz - γ_e))
            row_idx += 1
        end
        
        # Monta a matriz esparsa topológica (2*nvp linhas x 2*nvp colunas)
        A_topo = sparse(I_topo, J_topo, V_topo, 2*nvp, 2*nvp)

        # -----------------------------------------------------------------
        # 4. TRUST REGION LOOP (ADAPTATIVO)
        # -----------------------------------------------------------------
        step_accepted = false
        n_reject = 0
        last_rejected_dgamma = Float64[]
        
        while !step_accepted
            
            # Concatena as matrizes para o solver
            # Estrutura Final A:
            # [ Global (Vol, Per) (padded) ]
            # [ Move Limits (padded)       ]
            # [ Topological (sparse)       ]
            
            # Preparando Move Limits Temporários
            b_move = Float64[]
            A_move_rows = Vector{Float64}[] # Usando vetor de vetores para construir A_move denso depois
            
            if !isempty(elements_air)
                raw_limit = cv_current * nvp
                g_air = ceil(Int, raw_limit)
                if isodd(g_air); g_air -= 1; end; g_air = max(0, g_air)
                
                push!(b_move, g_air)
                push!(A_move_rows, one_air_ext)
            end

            if !isempty(elements_solid)
                raw_limit = cv_current * nvp
                g_solid = ceil(Int, raw_limit) 
                if isodd(g_solid); g_solid -= 1; end; g_solid = max(0, g_solid)

                push!(b_move, g_solid)
                push!(A_move_rows, -one_solid_ext) # Note o sinal negativo para <=
            end
            
            # Converte A_move para matriz
            if !isempty(A_move_rows)
                A_move = permutedims(hcat(A_move_rows...))
            else
                A_move = zeros(0, 2*nvp)
            end

            # Montagem Final do Sistema Linear
            # A_global é denso, A_move é denso, A_topo é esparso.
            # Convertendo tudo para Sparse para o solver (se suportado) ou Denso
            # Assumindo que LP aceita matriz genérica.
            
            # Concatenação Vertical
            # Atenção: vcat com esparsos e densos pode ser tricky no Julia.
            # Vamos garantir que tudo tenha 2*nvp colunas.
            
            # b_final
            b_final = vcat(b_global, b_move, b_topo)
            
            # A_final
            # A_global já tem 2*nvp cols
            # A_move já tem 2*nvp cols
            # A_topo já tem 2*nvp cols
            
            # É mais seguro converter tudo para sparse se o problema for grande
            A_final = vcat(sparse(A_global), sparse(A_move), A_topo)

            # Note que passamos 'γ[elements_design]' (apenas a parte física), 
            # mas 'c_aug', 'A_final' e 'b_final' já incluem as colunas/linhas das slacks.
            x_sol = LP(c_aug, A_final, b_final, γ[elements_design])

            # --- RESOLVE ILP (Integer Linear Programming) ---
            # Chamada ao solver. Supõe-se que ele retorne o vetor estendido x = [Delta_Gamma; Xi]
            # O solver deve tratar Xi como inteiros positivos ou contínuos >= 0.
            # Se LP for puramente binário, Xi precisaria ser discretizado, 
            # mas assumimos aqui que o backend (Alpine/HiGHS) suporta MIP.
            # x_sol = LP(c_aug, A_final, b_final, [γ[elements_design]; zeros(nvp)])

            # Extrai apenas a parte Delta Gamma
            Δγ = x_sol[1:nvp]
            
            # (Opcional) Debug das Folgas
            slack_vals = x_sol[nvp+1:end]
            max_slack = maximum(abs.(slack_vals))
            if max_slack > 1e-3
                 println("    -> Info: Restrições topológicas relaxadas (Max Slack = $max_slack).")
            end

            # [CHECK DE DEADLOCK]
            if !isempty(last_rejected_dgamma) && (Δγ == last_rejected_dgamma)
                println("    -> AVISO: Deadlock detectado. Abortando TR.")
                step_accepted = true; break 
            end

            # Predição Linear (Usando apenas custo físico)
            pred_reduction = -sum(c_phys .* Δγ)

            # Verifica mudança nula
            if sum(abs.(Δγ)) == 0
                 println("Aviso: Passo nulo. Próxima iter."); step_accepted = true; break
            elseif abs(pred_reduction) < 1e-12 
                 println("Aviso: Redução desprezível. Aceitando."); step_accepted = true; break 
            end

            # --- TRIAL STEP ---
            γ_trial = copy(γ)
            γ_trial[elements_design] .+= Δγ
            # Garante binário
            for ele in elements_design; γ_trial[ele] = round(γ_trial[ele]); end

            # --- GEOMETRIC CHECK ---
            P_trial = Perimiter(γ_trial, neighedge, elements_design)
            geom_violation = (Past > 0) && (P_trial > (Past + 1e-4))

            if geom_violation
                println("    -> REJEITADO (Geometria). P=$P_trial. Reduzindo cv.")
                last_rejected_dgamma = copy(Δγ)
                cv_current = max(cv_min, cv_current * 0.5)
                n_reject += 1
            else
                # --- PHYSICAL CHECK ---
                MP_trial, _, _, _ = Sweep(nn,ne,coord,connect,γ_trial,fρ,fκ,μ,freqs,livres,velocities,pressures)
                obj_trial = Objetivo(MP_trial,nodes_target,vA)
                actual_reduction = objetivo - obj_trial
                
                R = (abs(pred_reduction) < 1e-10) ? 0.0 : actual_reduction / pred_reduction
                println("--- Trust Region: Pred=$(round(pred_reduction,digits=4)) Act=$(round(actual_reduction,digits=4)) R=$(round(R,digits=2)) cv=$(round(cv_current,digits=4))")

                if R < 0.0 
                    println("    -> REJEITADO (Física). SPL piorou.")
                    last_rejected_dgamma = copy(Δγ)
                    cv_current = max(cv_min, cv_current * 0.5) 
                    n_reject +=1
                else
                    step_accepted = true
                    γ .= γ_trial 
                    if R > 0.7 
                         println("    -> ACEITO. Excelente."); if !is_warmup; cv_current = min(0.5, cv_current * 1.2); end
                    elseif R > 0.3
                        println("    -> ACEITO. Bom.")
                    else
                        println("    -> ACEITO. Pobre."); cv_current = max(cv_min, cv_current * 0.5)
                    end
                end
            end

            if n_reject >= 10
                println("Trust Region exausta. Desistindo.")
                step_accepted = true; break
            end
        end # Fim While Trust Region

        Lgmsh_export_element_scalar(arquivo_pos,γ,"Iter $iter")

        # --- VERIFICAÇÃO DE ESTAGNAÇÃO ---
        total_changes = sum(abs.(γ - γ_start))
        if total_changes == 0
            stagnation_counter += 1
            println("--- STAGNATION $stagnation_counter/4")
            if stagnation_counter >= 4; println("HARD STOP"); break; end
        else
            stagnation_counter = 0
        end

    end # Fim Loop iter

    # 5. FINALIZAÇÃO
    println("Finalizando e Exportando...")
    MP,_ =  Sweep(nn,ne,coord,connect,γ,fρ,fκ,μ,freqs,livres,velocities,pressures)
    writedlm(arquivo_γ_fin,γ)
    
    # Salva dados (omitido detalhes de escrita para brevidade)
    return historico_V, historico_SLP, historico_P

end # Otim_ISLP

#
# Rotina principal de Otimização Topológica (ISLP)
#
"""
 Otim_ISLP(arquivo::String, freqs::Vector, vA::Vector; verifica_derivada=false)
 
 Entradas:
 - arquivo: Caminho para o arquivo .geo ou .msh
 - freqs: Vetor de frequências para análise (Recomendado step=1Hz para ressonadores)
 - vA: Vetor de pesos para as frequências (geralmente ones(nf))
"""
function Otim_ISLP_artigo(arquivo::String, freqs::Vector, vA::Vector; verifica_derivada=false)

   # ---------------------------------------------------------------------------
   # 1. PREPARAÇÃO DA MALHA E ARQUIVOS
   # ---------------------------------------------------------------------------
   # Se o arquivo for um .geo, geramos um .msh utilizando a biblioteca do gmsh
   if occursin(".geo",arquivo)
      gmsh.initialize()
      gmsh.open(arquivo)
      gmsh.model.mesh.generate(2)
      mshfile = replace(arquivo,".geo"=>".msh")
      gmsh.write(mshfile)
   else 
      mshfile = arquivo
   end

    # Define os nomes dos arquivos de saída
    arquivo_yaml = replace(mshfile,".msh"=>".yaml")
    nomebase = basename(mshfile)
    arquivo_pos       = replace(nomebase,".msh"=>".pos")
    arquivo_pos_freq  = replace(nomebase,".msh"=>"_freq.pos")
    arquivo_data_opt  = replace(nomebase,".msh"=>".data")

    arquivo_γ_ini = replace(nomebase,".msh"=>"_γ_ini.dat")
    arquivo_γ_fin = replace(nomebase,".msh"=>"_γ_opt.dat")

    # Verificações de entrada básicas
    isempty(freqs) && error("Analise Harmonica:: freqs deve ser um vetor não vazio")
    isa(freqs,Vector{Float64}) || error("freqs deve ser um vetor de floats")
    isfile(mshfile) || error("Otim:: arquivo de entrada $mshfile não existe")
    nf = length(freqs)

    if isempty(vA)
       vA = ones(nf)
    else 
       length(vA)==nf || error("Otim:: dimensão de vA deve ser nf")
    end

    isfile(arquivo_yaml) || error("Otim:: arquivo de entrada $(arquivo_yaml) não existe")

    # Leitura da malha e parâmetros do YAML
    nn, coord, ne, connect, materials, nodes_open, velocities, nodes_pressure, pressures, damping, nodes_probe, nodes_target, elements_fixed, values_fixed, centroides = LSound.Parsemsh_Daniele(mshfile)
    elements_design = setdiff(1:ne,sort!(elements_fixed))
    nvp = length(elements_design) 

    raio_filtro, niter, ϵ1, ϵ2,  vf, Past, fatorcv, μ = Le_YAML(arquivo_yaml)
     
    isempty(nodes_target) && error("Otim:: nodes_target deve ter ao menos um nó informado")
    sort!(nodes_target)
    isempty(materials) && error("Analise:: at least one material is necessary")

    # Vizinhos de aresta, para cálculo do Perímetro
    neighedge = NeighborEdges(ne,connect,elements_design)
        
    # ---------------------------------------------------------------------------
    # 2. INICIALIZAÇÃO
    # ---------------------------------------------------------------------------
    println("Inicializando o vetor de variáveis de projeto com AR (Zeros).")
    γ = zeros(ne)
    Fix_γ!(γ,elements_fixed,values_fixed)
    writedlm(arquivo_γ_ini,γ)

    nodes_mask = sort(vcat(nodes_open,nodes_pressure))
    livres = setdiff(collect(1:nn),nodes_mask)
   
    # Inicializa os arquivos de visualização do gmsh (.pos)
    etype = connect[:,1]
    Lgmsh_export_init(arquivo_pos,nn,ne,coord,etype,connect[:,3:end])
    Lgmsh_export_init(arquivo_pos_freq,nn,ne,coord,etype,connect[:,3:end])
    Lgmsh_export_element_scalar(arquivo_pos,γ,"Iter 0")
   
    # Sweep Inicial (Análise Acústica da configuração inicial)
    MP,_ =  Sweep(nn,ne,coord,connect,γ,fρ,fκ,μ, freqs,livres,velocities,pressures)
    for i=1:nf
      f = freqs[i]
      Lgmsh_export_nodal_scalar(arquivo_pos_freq,abs.(MP[:,i]),"Pressure in $f Hz [abs] - initial topology")
    end

    # Rotina de verificação de derivada (debug)
    if verifica_derivada
       println("--- Modo de Verificação de Derivada ---")
       γ = rand(ne) # Perturbação aleatória
       MP,K,M,C = Sweep(nn,ne,coord,connect,γ,fρ,fκ,μ, freqs,livres,velocities,pressures)
       dΦ = Derivada(ne,nn,γ,connect,coord,K,M,C,livres,freqs,pressures,fρ, fκ,dfρ,dfκ,μ, nodes_target,MP,elements_design,vA) 
       dnum = Verifica_derivada(γ,nn,ne,coord,connect,fρ,fκ,μ,freqs,livres,velocities,pressures,nodes_target,elements_design,vA)
       rel = (dΦ.-dnum)./(dnum.+1E-12)
       Lgmsh_export_element_scalar(arquivo_pos,γ,"γ_debug")
       Lgmsh_export_element_scalar(arquivo_pos,dΦ,"Analitica")
       Lgmsh_export_element_scalar(arquivo_pos,dnum,"Numerica")
       Lgmsh_export_element_scalar(arquivo_pos,rel,"relativa")
       
       # Verificação derivada perímetro
       FPerimiter(γ) =  Perimiter(γ, neighedge, elements_design)
       dPnum = df(γ,FPerimiter,elements_design)
       dP = dPerimiter(ne, γ, neighedge, elements_design)
       return dΦ, dnum, dP, dPnum 
    end

    # --------------------------------------------------------------------------
    # 3. LOOP PRINCIPAL DE OTIMIZAÇÃO
    # --------------------------------------------------------------------------
    V = Volumes(ne,connect,coord)
    volume_full_projeto = sum(V[elements_design])
    Vast = vf*volume_full_projeto

    historico_V   = Float64[] 
    historico_SLP = Float64[] 
    historico_P   = Float64[] 

    # --- SETUP DOS PARÂMETROS DE CONTROLE ---
    cv_min = 1.0001 / nvp  # Move limit mínimo (pelo menos 1 elemento)
    n_warmup = 5 
    cv_warmup_cap = 0.05 
    
    cv_current = fatorcv 
    
    # Contador de Estagnação Global (Progressive Shake Down)
    stagnation_counter = 0

    for iter = 1:niter
        
        # Captura topologia no início para verificar mudanças depois
        γ_start = copy(γ)

        # --- LÓGICA DO PROGRESSIVE SHAKE DOWN ---
        if stagnation_counter > 0
            if stagnation_counter == 1
                println("--- STAGNATION RECOVERY (Nível 1): Resetando cv para 5%.")
                cv_current = 0.05
            elseif stagnation_counter == 2
                println("--- STAGNATION RECOVERY (Nível 2): Aumentando a aposta. cv = 15%.")
                cv_current = 0.15
            elseif stagnation_counter >= 3
                println("--- STAGNATION RECOVERY (Nível 3): O Martelo. cv = 30%.")
                cv_current = 0.30
            end
        
        # Recuperação Preventiva (se cv ficou muito pequeno sem estagnação)
        elseif cv_current < (2.0 / nvp) 
            println("--- RECUPERAÇÃO PREVENTIVA: cv muito baixo ($(round(cv_current,digits=5))). Resetando para 5%.")
            cv_current = 0.05
        end

        # --- WARM-UP LOGIC ---
        # Só aplicamos o teto do warm-up se NÃO estivermos tentando desestagnar.
        is_warmup = iter <= n_warmup
        if is_warmup && stagnation_counter == 0
            if cv_current > cv_warmup_cap
                cv_current = cv_warmup_cap
                println("--- WARM-UP (Iter $iter): Clamping cv to $(cv_warmup_cap)")
            end
        end

        # Identificação de elementos Ar/Sólido (para restrições lineares)
        elements_air = Int64[]
        elements_solid = Int64[]
        one_air = Float64[]
        one_solid = Float64[]

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

        volume_atual = sum(γ[elements_design].*V[elements_design])
        push!(historico_V , volume_atual)

        # --- ANÁLISE FEM ---
        MP,K,M,C =  Sweep(nn,ne,coord,connect,γ,fρ,fκ,μ,freqs,livres,velocities,pressures)
        objetivo = Objetivo(MP,nodes_target,vA)
        push!(historico_SLP, objetivo)

        perimetro = Perimiter(γ, neighedge, elements_design)
        push!(historico_P, perimetro)
        
        # Exporta resultados parciais
        for i=1:nf
         f = freqs[i]
         Lgmsh_export_nodal_scalar(arquivo_pos_freq,abs.(MP[:,i]),"Pressure in $f Hz [abs] $iter")
        end

        println("Iteração        ", iter)
        println("Objetivo (SPL)  ", objetivo)
        println("Perimetro       ", perimetro)
        println("Past (Limit)    ", Past)
        println("Volume atual    ", volume_atual)
        println("Volume target   ", Vast)
        println()

        # --- SENSIBILIDADE ---
        dΦ = Derivada(ne,nn,γ,connect,coord,K,M,C,livres,freqs,pressures,fρ,fκ,dfρ,dfκ,μ,nodes_target,MP,elements_design,vA) 
        Lgmsh_export_element_scalar(arquivo_pos,dΦ,"dΦ")  

        # Sensibilidade "crua" para o ILP
        c = dΦ[elements_design] 

        # --- RESTRIÇÕES LINEARES ---
        # 1. Volume
        ΔV = Vast - volume_atual
        b = [ΔV]
        println("Volume linearizado   ", b[1])
        A = vcat(V[elements_design]')

        # 2. Perímetro (se ativo)
        if  perimetro > 0 && Past>0
            ΔP = Past - perimetro
            b = vcat(b, ΔP)
            println("Perimetro linearizado   ", ΔP)
            dP = dPerimiter(ne, γ, neighedge, elements_design)
            Lgmsh_export_element_scalar(arquivo_pos,dP,"dP")  
            A = vcat(A,transpose(dP[elements_design]))
        end

        # -----------------------------------------------------------------
        # 4. TRUST REGION LOOP (ADAPTATIVO)
        # -----------------------------------------------------------------
        step_accepted = false
        n_reject = 0
        
        # Armazena a última solução rejeitada para evitar deadlock discreto
        last_rejected_dgamma = Float64[]
        
        while !step_accepted
            
            A_temp = copy(A)
            b_temp = copy(b)

            # --- Definição dos Limites de Movimento (Move Limits) ---
            if !isempty(elements_air)
                raw_limit = cv_current * length(elements_design)
                g_air = ceil(Int, raw_limit)
                
                # [MITIGAÇÃO DE ASSIMETRIA: ORÇAMENTO PAR]
                if isodd(g_air); g_air -= 1; end
                g_air = max(0, g_air)

                b_temp = vcat(b_temp, g_air)
                A_temp = vcat(A_temp, transpose(one_air))
            end

            if !isempty(elements_solid)
                raw_limit = cv_current * length(elements_design)
                g_solid = ceil(Int, raw_limit) 

                # [MITIGAÇÃO DE ASSIMETRIA: ORÇAMENTO PAR]
                if isodd(g_solid); g_solid -= 1; end
                g_solid = max(0, g_solid)

                b_temp = vcat(b_temp, g_solid)
                A_temp = vcat(A_temp, -transpose(one_solid))
            end

            # --- RESOLVE ILP (Integer Linear Programming) ---
            Δγ = LP(c, A_temp, b_temp, γ[elements_design])

            # [CHECK DE DEADLOCK]
            # Se a nova solução é idêntica à rejeitada anteriormente, o loop está travado.
            if !isempty(last_rejected_dgamma) && (Δγ == last_rejected_dgamma)
                println("    -> AVISO: ILP gerou passo idêntico ao rejeitado anteriormente.")
                println("    -> Redução de cv não alterou o conjunto viável inteiro. Abortando Trust Region.")
                # Sai do While, aceita o passo anterior (que é nulo se travou de primeira) e força próxima iter
                step_accepted = true 
                break 
            end

            # Predição Linear
            pred_reduction = -sum(c .* Δγ)

            # Verifica se houve mudança
            elements_changed = sum(abs.(Δγ))
            
            if elements_changed == 0
                 println("Aviso: ILP retornou passo nulo. Aceitando para próxima iter.")
                 step_accepted = true
                 break
            elseif abs(pred_reduction) < 1e-12 
                 println("Aviso: Redução predita desprezível. Aceitando.")
                 step_accepted = true
                 break 
            end

            # --- TRIAL STEP ---
            γ_trial = copy(γ)
            γ_trial[elements_design] .+= Δγ
            for ele in elements_design
                γ_trial[ele] = round(γ_trial[ele]) # Garante 0 ou 1
            end

            # --- GEOMETRIC FEASIBILITY CHECK (PERÍMETRO REAL) ---
            P_trial = Perimiter(γ_trial, neighedge, elements_design)
            geom_violation = (Past > 0) && (P_trial > (Past + 1e-4))

            if geom_violation
                println("    -> Passo REJEITADO (Geometria). P_real ($P_trial) > Past ($Past). Reduzindo cv.")
                last_rejected_dgamma = copy(Δγ) # Memoriza a falha
                cv_current = max(cv_min, cv_current * 0.5)
                n_reject += 1
            else
                # --- PHYSICAL FEASIBILITY CHECK (SPL REAL) ---
                MP_trial, _, _, _ = Sweep(nn,ne,coord,connect,γ_trial,fρ,fκ,μ,freqs,livres,velocities,pressures)
                obj_trial = Objetivo(MP_trial,nodes_target,vA)

                actual_reduction = objetivo - obj_trial
                
                if abs(pred_reduction) < 1e-10
                    R = 0.0
                else
                    R = actual_reduction / pred_reduction
                end

                println("--- Trust Region Eval: Pred = $(round(pred_reduction, digits=4)), Act = $(round(actual_reduction, digits=4)), R = $(round(R, digits=2)), cv = $(round(cv_current, digits=4))")

                if R < 0.0 
                    println("    -> Passo REJEITADO (Física). SPL piorou.")
                    last_rejected_dgamma = copy(Δγ) # Memoriza a falha
                    cv_current = max(cv_min, cv_current * 0.5) 
                    n_reject +=1
                else
                    # --- PASSO ACEITO ---
                    step_accepted = true
                    γ .= γ_trial 

                    if R > 0.7 
                        if is_warmup
                             println("    -> Passo ACEITO. Excelente.")
                        else
                             println("    -> Passo ACEITO. Excelente, expandindo cv.")
                             cv_current = min(0.5, cv_current * 1.2)
                        end
                    elseif R > 0.3
                        println("    -> Passo ACEITO. Boa, mantendo cv.")
                    else
                        println("    -> Passo ACEITO. Pobre, reduzindo cv.")
                        cv_current = max(cv_min, cv_current * 0.5)
                    end
                end
            end

            if n_reject >= 10
                println("Trust Region exausta (10 rejeições). Desistindo desta iteração.")
                step_accepted = true 
                break
            end
        end # Fim While Trust Region

        # Exporta topologia da iteração
        Lgmsh_export_element_scalar(arquivo_pos,γ,"Iter $iter")

        # --- VERIFICAÇÃO DE ESTAGNAÇÃO GLOBAL ---
        total_changes = sum(abs.(γ - γ_start))
        
        if total_changes == 0
            stagnation_counter += 1
            println("--- AVISO: Nenhuma alteração na topologia nesta iteração (Stagnation $stagnation_counter/4).")
            
            if stagnation_counter >= 4
                println("--- HARD STOP: Otimização estagnada mesmo após 'Kick' de 30%. Terminando.")
                break
            end
        else
            # Se mudou, reseta o contador
            stagnation_counter = 0
        end

    end # Fim Loop iter

    # --------------------------------------------------------------------------
    # 5. FINALIZAÇÃO
    # --------------------------------------------------------------------------
    println("Final da otimização, executando a análise SWEEP na topologia otimizada")
    MP,_ =  Sweep(nn,ne,coord,connect,γ,fρ,fκ,μ,freqs,livres,velocities,pressures)
    for i=1:nf
      f = freqs[i]
      Lgmsh_export_nodal_scalar(arquivo_pos_freq,abs.(MP[:,i]),"Pressure in $f Hz [abs]")
    end
    writedlm(arquivo_γ_fin,γ)
    
    fd = open(arquivo_data_opt,"w")
    println(fd,"Historico V")
    println(fd,historico_V)
    println(fd,"Historico SPL")
    println(fd,historico_SLP)
    println(fd,"Historico P")
    println(fd,historico_P)
    close(fd)
    

    return historico_V, historico_SLP, historico_P

end # Otim_ISLP



#
# Rotina principal
#
"""
 Otim(meshfile::String,freqs=[])
"""
function Otim_ISLP_MBMM(arquivo::String,freqs::Vector, vA::Vector;verifica_derivada=false)

   # Se o arquivo for um .geo, geramos um .msh utilizando a biblioteca do gmsh
   if occursin(".geo",arquivo)
      gmsh.initialize()
      gmsh.open(arquivo)
      gmsh.model.mesh.generate(2)
      mshfile = replace(arquivo,".geo"=>".msh")
      gmsh.write(mshfile)
   else 
      mshfile = arquivo
   end

    # Define os nomes dos arquivos
    arquivo_yaml = replace(mshfile,".msh"=>".yaml")
    nomebase = basename(mshfile)
    arquivo_pos       = replace(nomebase,".msh"=>".pos")
    arquivo_pos_freq  = replace(nomebase,".msh"=>"_freq.pos")

    arquivo_γ_ini = replace(nomebase,".msh"=>"_γ_ini.dat")
    arquivo_γ_fin = replace(nomebase,".msh"=>"_γ_opt.dat")

    # Verificações de entrada
    isempty(freqs) && error("Analise Harmonica:: freqs deve ser um vetor não vazio")
    isa(freqs,Vector{Float64}) || error("freqs deve ser um vetor de floats")
    isfile(mshfile) || error("Otim:: arquivo de entrada $mshfile não existe")
    nf = length(freqs)

    if isempty(vA)
       vA = ones(nf)
    else 
       length(vA)==nf || error("Otim:: dimensão de vA deve ser nf")
    end

    isfile(arquivo_yaml) || error("Otim:: arquivo de entrada $(arquivo_yaml) não existe")

    # Leitura da malha e YAML
    nn, coord, ne, connect, materials, nodes_open, velocities, nodes_pressure, pressures, damping, nodes_probe, nodes_target, elements_fixed, values_fixed, centroides = LSound.Parsemsh_Daniele(mshfile)
    elements_design = setdiff(1:ne,sort!(elements_fixed))
    nvp = length(elements_design) 

    raio_filtro, niter, ϵ1, ϵ2,  vf, Past, fatorcv, μ = Le_YAML(arquivo_yaml)
     
    isempty(nodes_target) && error("Otim:: nodes_target deve ter ao menos um nó informado")
    sort!(nodes_target)
    isempty(materials) && error("Analise:: at least one material is necessary")

    # Vizinhos de aresta, para cálculo do Perímetro
    neighedge = NeighborEdges(ne,connect,elements_design)
        
    # Inicialização
    println("Inicializando o vetor de variáveis de projeto")
    γ = zeros(ne)
    Fix_γ!(γ,elements_fixed,values_fixed)
    writedlm(arquivo_γ_ini,γ)

    nodes_mask = sort(vcat(nodes_open,nodes_pressure))
    livres = setdiff(collect(1:nn),nodes_mask)
   
    # Inicializa os arquivos de saída do gmsh
    etype = connect[:,1]
    Lgmsh_export_init(arquivo_pos,nn,ne,coord,etype,connect[:,3:end])
    Lgmsh_export_init(arquivo_pos_freq,nn,ne,coord,etype,connect[:,3:end])
    Lgmsh_export_element_scalar(arquivo_pos,γ,"Iter 0")
   
    # Sweep Inicial
    MP,_ =  Sweep(nn,ne,coord,connect,γ,fρ,fκ,μ, freqs,livres,velocities,pressures)
    for i=1:nf
      f = freqs[i]
      Lgmsh_export_nodal_scalar(arquivo_pos_freq,abs.(MP[:,i]),"Pressure in $f Hz [abs] - initial topology")
    end

    if verifica_derivada
       # ... (Mantido igual) ...
       γ = rand(ne)
       MP,K,M,C = Sweep(nn,ne,coord,connect,γ,fρ,fκ,μ, freqs,livres,velocities,pressures)
       dΦ = Derivada(ne,nn,γ,connect,coord,K,M,C,livres,freqs,pressures,dfρ,dfκ,μ, nodes_target,MP,elements_design,vA) 
       dnum = Verifica_derivada(γ,nn,ne,coord,connect,fρ,fκ,μ,freqs,livres,velocities,pressures,nodes_target,elements_design,vA)
       rel = (dΦ.-dnum)./(dnum.+1E-12)
       Lgmsh_export_element_scalar(arquivo_pos,γ,"γ")
       Lgmsh_export_element_scalar(arquivo_pos,dΦ,"Analitica")
       Lgmsh_export_element_scalar(arquivo_pos,dnum,"Numerica")
       Lgmsh_export_element_scalar(arquivo_pos,rel,"relativa")
       FPerimiter(γ) =  Perimiter(γ, neighedge, elements_design)
       dPnum = df(γ,FPerimiter,elements_design)
       dP = dPerimiter(ne, γ, neighedge, elements_design)
       return dΦ, dnum, dP, dPnum 
    end

    # --------------------------------------------------------------------------
    # Loop Principal de Otimização
    # --------------------------------------------------------------------------
    V = Volumes(ne,connect,coord)
    volume_full_projeto = sum(V[elements_design])
    Vast = vf*volume_full_projeto

    historico_V   = Float64[] 
    historico_SLP = Float64[] 
    historico_P   = Float64[] 

    # --- SETUP ---
    cv_min = 1.0001 / nvp  
    n_warmup = 5 
    cv_warmup_cap = 0.05 
    
    cv_current = fatorcv 
    
    # [NEW] Global Stagnation Counter
    stagnation_counter = 0

    for iter = 1:niter
        
        # Capture topology at start of iteration to check for changes later
        γ_start = copy(γ)

        # --- LÓGICA DO PROGRESSIVE SHAKE DOWN ---
        if stagnation_counter > 0
            if stagnation_counter == 1
                println("--- STAGNATION RECOVERY (Nível 1): Resetando cv para 5%.")
                cv_current = 0.05
            elseif stagnation_counter == 2
                println("--- STAGNATION RECOVERY (Nível 2): Aumentando a aposta. cv = 15%.")
                cv_current = 0.15
            elseif stagnation_counter >= 3
                println("--- STAGNATION RECOVERY (Nível 3): O Martelo. cv = 30%.")
                cv_current = 0.30
            end
        
        # Recuperação Preventiva (se cv ficou muito pequeno sem estagnação)
        elseif cv_current < (2.0 / nvp) 
            println("--- RECUPERAÇÃO PREVENTIVA: cv muito baixo ($(round(cv_current,digits=5))). Resetando para 5%.")
            cv_current = 0.05
        end

        # --- WARM-UP LOGIC ---
        is_warmup = iter <= n_warmup
        # Só aplicamos o CAP do warm-up se NÃO estivermos tentando desestagnar.
        # Se estivermos estagnados, o "Martelo" tem prioridade sobre o warm-up.
        if is_warmup && stagnation_counter == 0
            if cv_current > cv_warmup_cap
                cv_current = cv_warmup_cap
                println("--- WARM-UP (Iter $iter): Clamping cv to $(cv_warmup_cap)")
            end
        end

        # Identificação de elementos Ar/Sólido
        elements_air = Int64[]
        elements_solid = Int64[]
        one_air = Float64[]
        one_solid = Float64[]

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

        volume_atual = sum(γ[elements_design].*V[elements_design])
        push!(historico_V , volume_atual)

        # Análise FEM
        MP,K,M,C =  Sweep(nn,ne,coord,connect,γ,fρ,fκ,μ,freqs,livres,velocities,pressures)
        objetivo = Objetivo(MP,nodes_target,vA)
        push!(historico_SLP, objetivo)

        perimetro = Perimiter(γ, neighedge, elements_design)
        push!(historico_P, perimetro)
        
        for i=1:nf
         f = freqs[i]
         Lgmsh_export_nodal_scalar(arquivo_pos_freq,abs.(MP[:,i]),"Pressure in $f Hz [abs] $iter")
        end

        println("Iteração        ", iter)
        println("Objetivo        ", objetivo)
        println("Perimetro       ", perimetro)
        println("Past            ", Past)
        println("Volume atual    ", volume_atual)
        println("Volume target   ", Vast)
        println()

        # Sensibilidade
        dΦ = Derivada(ne,nn,γ,connect,coord,K,M,C,livres,freqs,pressures,dfρ,dfκ,μ,nodes_target,MP,elements_design,vA) 
        Lgmsh_export_element_scalar(arquivo_pos,dΦ,"dΦ")  

        # --- RAW SENSITIVITIES (NO FILTER) ---
        c = dΦ[elements_design] 

        # Restrições Lineares
        ΔV = Vast - volume_atual
        b = [ΔV]
        println("Volume linearizado   ", b[1])
        A = vcat(V[elements_design]')

        if  perimetro > 0 && Past>0
            ΔP = Past - perimetro
            b = vcat(b, ΔP)
            println("Perimetro linearizado   ", ΔP)
            dP = dPerimiter(ne, γ, neighedge, elements_design)
            Lgmsh_export_element_scalar(arquivo_pos,dP,"dP")  
            A = vcat(A,transpose(dP[elements_design]))
        end

        # -----------------------------------------------------------------
        # TRUST REGION LOOP
        # -----------------------------------------------------------------
        step_accepted = false
        n_reject = 0
        
        while !step_accepted
            
            A_temp = copy(A)
            b_temp = copy(b)

            if !isempty(elements_air)
                # Calcula limite baseado no cv
                raw_limit = cv_current * length(elements_design)
                g_air = ceil(Int, raw_limit) # Converte para Inteiro
                
                # --- [MITIGAÇÃO DE ASSIMETRIA: ORÇAMENTO PAR] ---
                # Se o limite for ímpar, reduz em 1 para garantir paridade.
                # Isso evita que o solver tenha que "escolher um lado" num empate.
                if isodd(g_air)
                    g_air -= 1
                end
                # Garante que não ficou negativo ou zero (se cv for muito pequeno)
                # Se ficar zero, o passo será nulo e o Shake Down resolverá na próxima.
                g_air = max(0, g_air)

                b_temp = vcat(b_temp, g_air)
                A_temp = vcat(A_temp, transpose(one_air))
            end

            if !isempty(elements_solid)
                raw_limit = cv_current * length(elements_design)
                g_solid = ceil(Int, raw_limit)

                # --- [MITIGAÇÃO DE ASSIMETRIA: ORÇAMENTO PAR] ---
                if isodd(g_solid)
                    g_solid -= 1
                end
                g_solid = max(0, g_solid)

                b_temp = vcat(b_temp, g_solid)
                A_temp = vcat(A_temp, -transpose(one_solid))
            end

            # Solve ILP
            # ... (o resto continua igual)
            # Solve ILP
            Δγ = LP(c, A_temp, b_temp, γ[elements_design])

            # Predição
            pred_reduction = -sum(c .* Δγ)

            # Convergence Checks
            elements_changed = sum(abs.(Δγ))
            
            if elements_changed == 0
                 println("Aviso: ILP retornou passo nulo. Aceitando para próxima iter.")
                 step_accepted = true
                 break
            elseif abs(pred_reduction) < 1e-12 
                 println("Aviso: Redução predita desprezível. Aceitando.")
                 step_accepted = true
                 break 
            end

            # Trial
            γ_trial = copy(γ)
            γ_trial[elements_design] .+= Δγ
            for ele in elements_design
                γ_trial[ele] = round(γ_trial[ele])
            end

            # --- GEOMETRIC FEASIBILITY CHECK ---
            P_trial = Perimiter(γ_trial, neighedge, elements_design)
            geom_violation = (Past > 0) && (P_trial > (Past + 1e-4))

            if geom_violation
                println("    -> Passo REJEITADO (Geometria). P_real ($P_trial) > Past ($Past). Reduzindo cv.")
                cv_current = max(cv_min, cv_current * 0.5)
                n_reject += 1
            else
                # Se a geometria passou, checamos a FÍSICA
                MP_trial, _, _, _ = Sweep(nn,ne,coord,connect,γ_trial,fρ,fκ,μ,freqs,livres,velocities,pressures)
                obj_trial = Objetivo(MP_trial,nodes_target,vA)

                actual_reduction = objetivo - obj_trial
                R = actual_reduction / pred_reduction

                println("--- Trust Region Eval: Pred = $(round(pred_reduction, digits=4)), Act = $(round(actual_reduction, digits=4)), R = $(round(R, digits=2)), cv = $(round(cv_current, digits=4))")

                if R < 0.0 
                    println("    -> Passo REJEITADO (Física). SPL piorou.")
                    cv_current = max(cv_min, cv_current * 0.5) 
                    n_reject +=1
                else
                    # ACEITO
                    step_accepted = true
                    γ .= γ_trial 

                    if R > 0.7 
                        if is_warmup
                             println("    -> Passo ACEITO. Excelente.")
                        else
                             println("    -> Passo ACEITO. Excelente, expandindo cv.")
                             cv_current = min(0.5, cv_current * 1.2)
                        end
                    elseif R > 0.3
                        println("    -> Passo ACEITO. Boa, mantendo cv.")
                    else
                        println("    -> Passo ACEITO. Pobre, reduzindo cv.")
                        cv_current = max(cv_min, cv_current * 0.5)
                    end
                end
            end

            if n_reject >= 10
                println("Trust Region exausta (10 rejeições). Desistindo desta iteração.")
                step_accepted = true 
                break
            end
        end # While Trust Region

        Lgmsh_export_element_scalar(arquivo_pos,γ,"Iter $iter")

        # --- CHECK FOR GLOBAL STAGNATION ---
        total_changes = sum(abs.(γ - γ_start))
        
        if total_changes == 0
            stagnation_counter += 1
            println("--- AVISO: Nenhuma alteração na topologia nesta iteração (Stagnation $stagnation_counter/4).")
            
            if stagnation_counter >= 4
                println("--- HARD STOP: Otimização estagnada mesmo após 'Kick' de 30%. Terminando.")
                break
            end
        else
            # Reset counter if we managed to move
            stagnation_counter = 0
        end

    end # Loop iter

    println("Final da otimização, executando a análise SWEEP na topologia otimizada")
    MP,_ =  Sweep(nn,ne,coord,connect,γ,fρ,fκ,μ,freqs,livres,velocities,pressures)
    for i=1:nf
      f = freqs[i]
      Lgmsh_export_nodal_scalar(arquivo_pos_freq,abs.(MP[:,i]),"Pressure in $f Hz [abs]")
    end
    writedlm(arquivo_γ_fin,γ)
    return historico_V, historico_SLP, historico_P

end # main_otim


#
# Rotina principal
#
"""
 Otim(meshfile::String,freqs=[])
"""
function Otim_ISLP_FUNFANDO(arquivo::String,freqs::Vector, vA::Vector;verifica_derivada=false)

   # Se o arquivo for um .geo, geramos um .msh utilizando a biblioteca do gmsh
   if occursin(".geo",arquivo)
      gmsh.initialize()
      gmsh.open(arquivo)
      gmsh.model.mesh.generate(2)
      mshfile = replace(arquivo,".geo"=>".msh")
      gmsh.write(mshfile)
   else 
      mshfile = arquivo
   end

    # Define os nomes dos arquivos
    arquivo_yaml = replace(mshfile,".msh"=>".yaml")
    nomebase = basename(mshfile)
    arquivo_pos       = replace(nomebase,".msh"=>".pos")
    arquivo_pos_freq  = replace(nomebase,".msh"=>"_freq.pos")

    arquivo_γ_ini = replace(nomebase,".msh"=>"_γ_ini.dat")
    arquivo_γ_fin = replace(nomebase,".msh"=>"_γ_opt.dat")

    # Verificações de entrada
    isempty(freqs) && error("Analise Harmonica:: freqs deve ser um vetor não vazio")
    isa(freqs,Vector{Float64}) || error("freqs deve ser um vetor de floats")
    isfile(mshfile) || error("Otim:: arquivo de entrada $mshfile não existe")
    nf = length(freqs)

    if isempty(vA)
       vA = ones(nf)
    else 
       length(vA)==nf || error("Otim:: dimensão de vA deve ser nf")
    end

    isfile(arquivo_yaml) || error("Otim:: arquivo de entrada $(arquivo_yaml) não existe")

    # Leitura da malha e YAML
    nn, coord, ne, connect, materials, nodes_open, velocities, nodes_pressure, pressures, damping, nodes_probe, nodes_target, elements_fixed, values_fixed, centroides = LSound.Parsemsh_Daniele(mshfile)
    elements_design = setdiff(1:ne,sort!(elements_fixed))
    nvp = length(elements_design) 

    raio_filtro, niter, ϵ1, ϵ2,  vf, Past, fatorcv, μ = Le_YAML(arquivo_yaml)
     
    isempty(nodes_target) && error("Otim:: nodes_target deve ter ao menos um nó informado")
    sort!(nodes_target)
    isempty(materials) && error("Analise:: at least one material is necessary")

    # Pré-processamento
    if !verifica_derivada
         println("Determinando a vizinhança para um raio de $(raio_filtro)")
         vizinhos, pesos = Vizinhanca(ne,centroides,raio_filtro,elements_design)
    end 
    neighedge = NeighborEdges(ne,connect,elements_design)
        
    # Inicialização
    println("Inicializando o vetor de variáveis de projeto")
    γ = zeros(ne)
    Fix_γ!(γ,elements_fixed,values_fixed)
    writedlm(arquivo_γ_ini,γ)

    nodes_mask = sort(vcat(nodes_open,nodes_pressure))
    livres = setdiff(collect(1:nn),nodes_mask)
   
    # Inicializa os arquivos de saída do gmsh
    etype = connect[:,1]
    Lgmsh_export_init(arquivo_pos,nn,ne,coord,etype,connect[:,3:end])
    Lgmsh_export_init(arquivo_pos_freq,nn,ne,coord,etype,connect[:,3:end])
    Lgmsh_export_element_scalar(arquivo_pos,γ,"Iter 0")
   
    # Sweep Inicial
    MP,_ =  Sweep(nn,ne,coord,connect,γ,fρ,fκ,μ, freqs,livres,velocities,pressures)
    for i=1:nf
      f = freqs[i]
      Lgmsh_export_nodal_scalar(arquivo_pos_freq,abs.(MP[:,i]),"Pressure in $f Hz [abs] - initial topology")
    end

    if verifica_derivada
       # ... (Mantido igual) ...
       γ = rand(ne)
       MP,K,M,C = Sweep(nn,ne,coord,connect,γ,fρ,fκ,μ, freqs,livres,velocities,pressures)
       dΦ = Derivada(ne,nn,γ,connect,coord,K,M,C,livres,freqs,pressures,dfρ,dfκ,μ, nodes_target,MP,elements_design,vA) 
       dnum = Verifica_derivada(γ,nn,ne,coord,connect,fρ,fκ,μ,freqs,livres,velocities,pressures,nodes_target,elements_design,vA)
       rel = (dΦ.-dnum)./(dnum.+1E-12)
       Lgmsh_export_element_scalar(arquivo_pos,γ,"γ")
       Lgmsh_export_element_scalar(arquivo_pos,dΦ,"Analitica")
       Lgmsh_export_element_scalar(arquivo_pos,dnum,"Numerica")
       Lgmsh_export_element_scalar(arquivo_pos,rel,"relativa")
       FPerimiter(γ) =  Perimiter(γ, neighedge, elements_design)
       dPnum = df(γ,FPerimiter,elements_design)
       dP = dPerimiter(ne, γ, neighedge, elements_design)
       return dΦ, dnum, dP, dPnum 
    end

    # --------------------------------------------------------------------------
    # Loop Principal de Otimização
    # --------------------------------------------------------------------------
    V = Volumes(ne,connect,coord)
    volume_full_projeto = sum(V[elements_design])
    Vast = vf*volume_full_projeto

    historico_V   = Float64[] 
    historico_SLP = Float64[] 
    historico_P   = Float64[] 

    # --- SETUP ---
    cv_min = 1.0001 / nvp  
    n_warmup = 5 
    cv_warmup_cap = 0.05 
    
    # Initialization
    cv_current = fatorcv 

    for iter = 1:niter

        # --- [NEW STRATEGY] STAGNATION RECOVERY ("SHAKE DOWN") ---
        # If cv is too small (e.g., less than 2 elements), boost it back to 5%
        # This applies to ALL iterations (Warm-up or Normal)
        if cv_current < (2.0 / nvp) 
            println("--- STAGNATION RECOVERY: cv too low ($(round(cv_current,digits=5))). Boosting to 5%.")
            cv_current = 0.05
        end

        # --- WARM-UP LOGIC ---
        # Only enforces the CAP (Upper limit). The Floor is handled by the Recovery above.
        is_warmup = iter <= n_warmup
        if is_warmup
            if cv_current > cv_warmup_cap
                cv_current = cv_warmup_cap
                println("--- WARM-UP (Iter $iter): Clamping cv to $(cv_warmup_cap)")
            end
        end

        # Identificação de elementos Ar/Sólido
        elements_air = Int64[]
        elements_solid = Int64[]
        one_air = Float64[]
        one_solid = Float64[]

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

        volume_atual = sum(γ[elements_design].*V[elements_design])
        push!(historico_V , volume_atual)

        # Análise FEM
        MP,K,M,C =  Sweep(nn,ne,coord,connect,γ,fρ,fκ,μ,freqs,livres,velocities,pressures)
        objetivo = Objetivo(MP,nodes_target,vA)
        push!(historico_SLP, objetivo)

        perimetro = Perimiter(γ, neighedge, elements_design)
        push!(historico_P, perimetro)
        
        for i=1:nf
         f = freqs[i]
         Lgmsh_export_nodal_scalar(arquivo_pos_freq,abs.(MP[:,i]),"Pressure in $f Hz [abs] $iter")
        end

        println("Iteração        ", iter)
        println("Objetivo        ", objetivo)
        println("Perimetro       ", perimetro)
        println("Past            ", Past)
        println("Volume atual    ", volume_atual)
        println("Volume target   ", Vast)
        println()

        # Sensibilidade
        dΦ = Derivada(ne,nn,γ,connect,coord,K,M,C,livres,freqs,pressures,dfρ,dfκ,μ,nodes_target,MP,elements_design,vA) 
        Lgmsh_export_element_scalar(arquivo_pos,dΦ,"dΦ")  

        # --- RAW SENSITIVITIES (NO FILTER) ---
        c = dΦ[elements_design] 

        # Restrições Lineares
        ΔV = Vast - volume_atual
        b = [ΔV]
        println("Volume linearizado   ", b[1])
        A = vcat(V[elements_design]')

        if  perimetro > 0 && Past>0
            ΔP = Past - perimetro
            b = vcat(b, ΔP)
            println("Perimetro linearizado   ", ΔP)
            dP = dPerimiter(ne, γ, neighedge, elements_design)
            Lgmsh_export_element_scalar(arquivo_pos,dP,"dP")  
            A = vcat(A,transpose(dP[elements_design]))
        end
       
        # -----------------------------------------------------------------
        # TRUST REGION LOOP
        # -----------------------------------------------------------------
        step_accepted = false
        n_reject = 0
        flag_podrao = false

        while !step_accepted
            
            A_temp = copy(A)
            b_temp = copy(b)

            if !isempty(elements_air)
                g_air = ceil(cv_current * length(elements_design))
                b_temp = vcat(b_temp, g_air)
                A_temp = vcat(A_temp, transpose(one_air))
            end

            if !isempty(elements_solid)
                g_solid  = ceil(cv_current * length(elements_design)) 
                b_temp = vcat(b_temp, g_solid)
                A_temp = vcat(A_temp, -transpose(one_solid))
            end

            # Solve ILP
            Δγ = LP(c, A_temp, b_temp, γ[elements_design])

            # Predição
            pred_reduction = -sum(c .* Δγ)

            # Convergence Check: Relaxed
            elements_changed = sum(abs.(Δγ))
            
            if elements_changed == 0
                 # Se não mudou nada, aceitamos o passo nulo para que o loop externo
                 # possa rodar novamente e triggar o "Shake Down" se necessário.
                 println("Aviso: ILP retornou passo nulo. Aceitando para próxima iter.")
                 step_accepted = true
                 break
            elseif abs(pred_reduction) < 1e-12 
                 # Mesmo para redução pequena, aceita e continua.
                 println("Aviso: Redução predita desprezível. Aceitando.")
                 step_accepted = true
                 break 
            end

            # Trial
            γ_trial = copy(γ)
            γ_trial[elements_design] .+= Δγ
            for ele in elements_design
                γ_trial[ele] = round(γ_trial[ele])
            end

            # --- GEOMETRIC FEASIBILITY CHECK ---
            P_trial = Perimiter(γ_trial, neighedge, elements_design)
            geom_violation = (Past > 0) && (P_trial > (Past + 1e-4))

            if geom_violation
                println("    -> Passo REJEITADO (Geometria). P_real ($P_trial) > Past ($Past). Reduzindo cv.")
                cv_current = max(cv_min, cv_current * 0.5)
                n_reject += 1
            else
                # Se a geometria passou, checamos a FÍSICA
                MP_trial, _, _, _ = Sweep(nn,ne,coord,connect,γ_trial,fρ,fκ,μ,freqs,livres,velocities,pressures)
                obj_trial = Objetivo(MP_trial,nodes_target,vA)

                actual_reduction = objetivo - obj_trial
                R = actual_reduction / pred_reduction

                println("--- Trust Region Eval: Pred = $(round(pred_reduction, digits=4)), Act = $(round(actual_reduction, digits=4)), R = $(round(R, digits=2)), cv = $(round(cv_current, digits=4))")

                if R < 0.0 
                    println("    -> Passo REJEITADO (Física). SPL piorou.")
                    cv_current = max(cv_min, cv_current * 0.5) 
                    n_reject +=1
                else
                    # ACEITO
                    step_accepted = true
                    γ .= γ_trial 

                    if R > 0.7 
                        if is_warmup
                             # No warmup, não expande acima do teto, mas ok.
                             println("    -> Passo ACEITO. Excelente.")
                        else
                             println("    -> Passo ACEITO. Excelente, expandindo cv.")
                             cv_current = min(0.5, cv_current * 1.2)
                        end
                    elseif R > 0.3
                        println("    -> Passo ACEITO. Boa, mantendo cv.")
                    else
                        println("    -> Passo ACEITO. Pobre, reduzindo cv.")
                        cv_current = max(cv_min, cv_current * 0.5)
                    end
                end
            end

            if n_reject >= 10
                println("Trust Region falhou 10 vezes. Forçando saída para reset.")
                # Em vez de break total, aceitamos passo nulo para resetar cv na proxima
                step_accepted = true 
                break
            end
        end # While Trust Region

        Lgmsh_export_element_scalar(arquivo_pos,γ,"Iter $iter")

    end # Loop iter

    println("Final da otimização, executando a análise SWEEP na topologia otimizada")
    MP,_ =  Sweep(nn,ne,coord,connect,γ,fρ,fκ,μ,freqs,livres,velocities,pressures)
    for i=1:nf
      f = freqs[i]
      Lgmsh_export_nodal_scalar(arquivo_pos_freq,abs.(MP[:,i]),"Pressure in $f Hz [abs]")
    end
    writedlm(arquivo_γ_fin,γ)
    return historico_V, historico_SLP, historico_P

end # main_otim
#
# Rotina principal
#
"""
 Otim(meshfile::String,freqs=[])
"""
function Otim_ISLP_GOOD(arquivo::String,freqs::Vector, vA::Vector;verifica_derivada=false)

   # Se o arquivo for um .geo, geramos um .msh utilizando a biblioteca do gmsh
   if occursin(".geo",arquivo)
      gmsh.initialize()
      gmsh.open(arquivo)
      gmsh.model.mesh.generate(2)
      mshfile = replace(arquivo,".geo"=>".msh")
      gmsh.write(mshfile)
   else 
      mshfile = arquivo
   end

    # Define os nomes dos arquivos
    arquivo_yaml = replace(mshfile,".msh"=>".yaml")
    nomebase = basename(mshfile)
    arquivo_pos       = replace(nomebase,".msh"=>".pos")
    arquivo_pos_freq  = replace(nomebase,".msh"=>"_freq.pos")

    arquivo_γ_ini = replace(nomebase,".msh"=>"_γ_ini.dat")
    arquivo_γ_fin = replace(nomebase,".msh"=>"_γ_opt.dat")

    # Verificações de entrada
    isempty(freqs) && error("Analise Harmonica:: freqs deve ser um vetor não vazio")
    isa(freqs,Vector{Float64}) || error("freqs deve ser um vetor de floats")
    isfile(mshfile) || error("Otim:: arquivo de entrada $mshfile não existe")
    nf = length(freqs)

    if isempty(vA)
       vA = ones(nf)
    else 
       length(vA)==nf || error("Otim:: dimensão de vA deve ser nf")
    end

    isfile(arquivo_yaml) || error("Otim:: arquivo de entrada $(arquivo_yaml) não existe")

    # Leitura da malha e YAML
    nn, coord, ne, connect, materials, nodes_open, velocities, nodes_pressure, pressures, damping, nodes_probe, nodes_target, elements_fixed, values_fixed, centroides = LSound.Parsemsh_Daniele(mshfile)
    elements_design = setdiff(1:ne,sort!(elements_fixed))
    nvp = length(elements_design) 

    raio_filtro, niter, ϵ1, ϵ2,  vf, Past, fatorcv, μ = Le_YAML(arquivo_yaml)
     
    isempty(nodes_target) && error("Otim:: nodes_target deve ter ao menos um nó informado")
    sort!(nodes_target)
    isempty(materials) && error("Analise:: at least one material is necessary")

    # Pré-processamento
    if !verifica_derivada
         println("Determinando a vizinhança para um raio de $(raio_filtro)")
         vizinhos, pesos = Vizinhanca(ne,centroides,raio_filtro,elements_design)
    end 
    neighedge = NeighborEdges(ne,connect,elements_design)
        
    # Inicialização
    println("Inicializando o vetor de variáveis de projeto")
    γ = zeros(ne)
    Fix_γ!(γ,elements_fixed,values_fixed)
    writedlm(arquivo_γ_ini,γ)

    nodes_mask = sort(vcat(nodes_open,nodes_pressure))
    livres = setdiff(collect(1:nn),nodes_mask)
   
    # Inicializa os arquivos de saída do gmsh
    etype = connect[:,1]
    Lgmsh_export_init(arquivo_pos,nn,ne,coord,etype,connect[:,3:end])
    Lgmsh_export_init(arquivo_pos_freq,nn,ne,coord,etype,connect[:,3:end])
    Lgmsh_export_element_scalar(arquivo_pos,γ,"Iter 0")
   
    # Sweep Inicial
    MP,_ =  Sweep(nn,ne,coord,connect,γ,fρ,fκ,μ, freqs,livres,velocities,pressures)
    for i=1:nf
      f = freqs[i]
      Lgmsh_export_nodal_scalar(arquivo_pos_freq,abs.(MP[:,i]),"Pressure in $f Hz [abs] - initial topology")
    end

    if verifica_derivada
       # ... (Mantido igual ao original) ...
       γ = rand(ne)
       MP,K,M,C = Sweep(nn,ne,coord,connect,γ,fρ,fκ,μ, freqs,livres,velocities,pressures)
       dΦ = Derivada(ne,nn,γ,connect,coord,K,M,C,livres,freqs,pressures,dfρ,dfκ,μ, nodes_target,MP,elements_design,vA) 
       dnum = Verifica_derivada(γ,nn,ne,coord,connect,fρ,fκ,μ,freqs,livres,velocities,pressures,nodes_target,elements_design,vA)
       rel = (dΦ.-dnum)./(dnum.+1E-12)
       Lgmsh_export_element_scalar(arquivo_pos,γ,"γ")
       Lgmsh_export_element_scalar(arquivo_pos,dΦ,"Analitica")
       Lgmsh_export_element_scalar(arquivo_pos,dnum,"Numerica")
       Lgmsh_export_element_scalar(arquivo_pos,rel,"relativa")
       FPerimiter(γ) =  Perimiter(γ, neighedge, elements_design)
       dPnum = df(γ,FPerimiter,elements_design)
       dP = dPerimiter(ne, γ, neighedge, elements_design)
       return dΦ, dnum, dP, dPnum 
    end

    # --------------------------------------------------------------------------
    # Loop Principal de Otimização
    # --------------------------------------------------------------------------
    V = Volumes(ne,connect,coord)
    volume_full_projeto = sum(V[elements_design])
    Vast = vf*volume_full_projeto

    historico_V   = Float64[] 
    historico_SLP = Float64[] 
    historico_P   = Float64[] 

    # --- SETUP DA ESTABILIDADE ---
    cv_min = 1.0001 / nvp  
    n_warmup = 1
    cv_warmup_cap = 0.05 

    cv_current = fatorcv 

    for iter = 1:niter

        # --- LÓGICA DO WARM-UP (ROBUST) ---
        is_warmup = iter <= n_warmup
        
        if is_warmup
            # FIX 1: FORCE RESET. 
            # Não deixe as rejeições anteriores matarem o warm-up. 
            # Começa toda iteração de aquecimento com 5% de liberdade.
            cv_current = cv_warmup_cap
            println("--- WARM-UP (Iter $iter): Forçando cv = $(cv_warmup_cap) para evitar estagnação.")
        
        elseif iter == (n_warmup + 1)
             # FIX 2: POST-WARMUP BOOST
             # Acabou o aquecimento? Chuta o balde para 10% para sair de mínimos locais.
             println("--- WARM-UP FINISHED. Resetting search radius to 10%.")
             cv_current = 0.1
        end

        # Identificação de elementos Ar/Sólido
        elements_air = Int64[]
        elements_solid = Int64[]
        one_air = Float64[]
        one_solid = Float64[]

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

        volume_atual = sum(γ[elements_design].*V[elements_design])
        push!(historico_V , volume_atual)

        # Análise FEM
        MP,K,M,C =  Sweep(nn,ne,coord,connect,γ,fρ,fκ,μ,freqs,livres,velocities,pressures)
        objetivo = Objetivo(MP,nodes_target,vA)
        push!(historico_SLP, objetivo)

        perimetro = Perimiter(γ, neighedge, elements_design)
        push!(historico_P, perimetro)
        
        for i=1:nf
         f = freqs[i]
         Lgmsh_export_nodal_scalar(arquivo_pos_freq,abs.(MP[:,i]),"Pressure in $f Hz [abs] $iter")
        end

        println("Iteração        ", iter)
        println("Objetivo        ", objetivo)
        println("Perimetro       ", perimetro)
        println("Past            ", Past)
        println("Volume atual    ", volume_atual)
        println("Volume target   ", Vast)
        println()

        # Sensibilidade
        dΦ = Derivada(ne,nn,γ,connect,coord,K,M,C,livres,freqs,pressures,dfρ,dfκ,μ,nodes_target,MP,elements_design,vA) 
        Lgmsh_export_element_scalar(arquivo_pos,dΦ,"dΦ")  

        # --- RAW SENSITIVITIES (NO FILTER) ---
        c = dΦ[elements_design] 

        # Restrições Lineares
        ΔV = Vast - volume_atual
        b = [ΔV]
        println("Volume linearizado   ", b[1])
        A = vcat(V[elements_design]')

        if  perimetro > 0 && Past>0
            ΔP = Past - perimetro
            b = vcat(b, ΔP)
            println("Perimetro linearizado   ", ΔP)
            dP = dPerimiter(ne, γ, neighedge, elements_design)
            Lgmsh_export_element_scalar(arquivo_pos,dP,"dP")  
            A = vcat(A,transpose(dP[elements_design]))
        end
       
        # -----------------------------------------------------------------
        # TRUST REGION LOOP
        # -----------------------------------------------------------------
        step_accepted = false
        n_reject = 0
        flag_podrao = false

        while !step_accepted
            
            A_temp = copy(A)
            b_temp = copy(b)

            if !isempty(elements_air)
                g_air = ceil(cv_current * length(elements_design))
                b_temp = vcat(b_temp, g_air)
                A_temp = vcat(A_temp, transpose(one_air))
            end

            if !isempty(elements_solid)
                g_solid  = ceil(cv_current * length(elements_design)) 
                b_temp = vcat(b_temp, g_solid)
                A_temp = vcat(A_temp, -transpose(one_solid))
            end

            # Solve ILP
            Δγ = LP(c, A_temp, b_temp, γ[elements_design])

            # Predição
            pred_reduction = -sum(c .* Δγ)

            elements_changed = sum(abs.(Δγ))
            
            # FIX 3: CONVERGENCE GUARD
            # Se não mudou nada, mas estamos no warm-up, NÃO PARE.
            # Force o loop a continuar (mesmo que seja um passo nulo, o reset de cv na próxima iter resolve).
            if elements_changed == 0
                 if is_warmup
                     println("Aviso: ILP não encontrou movimentos (Warm-up). Aceitando passo nulo e continuando.")
                     step_accepted = true
                     break # Sai do Trust Region, vai para a próxima iter (que resetará o cv)
                 else
                     println("Convergência: Nenhum elemento alterado no passo ILP.")
                     step_accepted = true
                     break # Sai do Trust Region. O loop externo vai continuar se não tiver flag de parada.
                 end
            elseif abs(pred_reduction) < 1e-12 
                 if is_warmup
                     println("Aviso: Redução desprezível (Warm-up). Continuando.")
                     step_accepted = true
                     break
                 else
                     println("Convergência local atingida (Redução < 1e-12).")
                     step_accepted = true
                     # Só paramos o OTIMIZADOR inteiro se não for warm-up
                     # Mas aqui é break do while. Para parar o loop externo, precisamos de uma flag.
                     # Vamos deixar rodar. Frequentemente é um platô temporário.
                     break 
                 end
            end

            # Trial
            γ_trial = copy(γ)
            γ_trial[elements_design] .+= Δγ
            for ele in elements_design
                γ_trial[ele] = round(γ_trial[ele])
            end

            # --- GEOMETRIC FEASIBILITY CHECK ---
            P_trial = Perimiter(γ_trial, neighedge, elements_design)
            
            geom_violation = (Past > 0) && (P_trial > (Past + 1e-4))

            if geom_violation
                println("    -> Passo REJEITADO (Geometria). P_real ($P_trial) > Past ($Past). Reduzindo cv.")
                cv_current = max(cv_min, cv_current * 0.5)
                n_reject += 1
            else
                # Se a geometria passou, checamos a FÍSICA
                MP_trial, _, _, _ = Sweep(nn,ne,coord,connect,γ_trial,fρ,fκ,μ,freqs,livres,velocities,pressures)
                obj_trial = Objetivo(MP_trial,nodes_target,vA)

                actual_reduction = objetivo - obj_trial
                R = actual_reduction / pred_reduction

                println("--- Trust Region Eval: Pred = $(round(pred_reduction, digits=4)), Act = $(round(actual_reduction, digits=4)), R = $(round(R, digits=2)), cv = $(round(cv_current, digits=4))")

                if R < 0.0 
                    println("    -> Passo REJEITADO (Física). SPL piorou.")
                    cv_current = max(cv_min, cv_current * 0.5) 
                    n_reject +=1
                else
                    # ACEITO
                    step_accepted = true
                    γ .= γ_trial 

                    if R > 0.7 
                        if is_warmup
                            println("    -> Passo ACEITO. Excelente, mantendo cv (Warm-up).")
                        else
                            println("    -> Passo ACEITO. Excelente, expandindo cv.")
                            cv_current = min(0.5, cv_current * 1.2)
                        end
                    elseif R > 0.3
                        println("    -> Passo ACEITO. Boa, mantendo cv.")
                    else
                        println("    -> Passo ACEITO. Pobre, reduzindo cv.")
                        cv_current = max(cv_min, cv_current * 0.5)
                    end
                end
            end

            if n_reject >= 10
                flag_podrao = true
                println("Trust Region falhou 10 vezes. Terminando.")
                break
            end
        end # While Trust Region

        if flag_podrao
            break
        end  

        Lgmsh_export_element_scalar(arquivo_pos,γ,"Iter $iter")

    end # Loop iter

    println("Final da otimização, executando a análise SWEEP na topologia otimizada")
    MP,_ =  Sweep(nn,ne,coord,connect,γ,fρ,fκ,μ,freqs,livres,velocities,pressures)
    for i=1:nf
      f = freqs[i]
      Lgmsh_export_nodal_scalar(arquivo_pos_freq,abs.(MP[:,i]),"Pressure in $f Hz [abs]")
    end
    writedlm(arquivo_γ_fin,γ)
    return historico_V, historico_SLP, historico_P

end # main_otim

#
# Rotina principal
#
"""
 Otim(meshfile::String,freqs=[])
"""
function Otim_ISLP_quase(arquivo::String,freqs::Vector, vA::Vector;verifica_derivada=false)

   # Se o arquivo for um .geo, geramos um .msh utilizando a biblioteca do gmsh
   if occursin(".geo",arquivo)
      gmsh.initialize()
      gmsh.open(arquivo)
      gmsh.model.mesh.generate(2)
      mshfile = replace(arquivo,".geo"=>".msh")
      gmsh.write(mshfile)
   else 
      mshfile = arquivo
   end

    # Define os nomes dos arquivos
    arquivo_yaml = replace(mshfile,".msh"=>".yaml")
    nomebase = basename(mshfile)
    arquivo_pos       = replace(nomebase,".msh"=>".pos")
    arquivo_pos_freq  = replace(nomebase,".msh"=>"_freq.pos")

    arquivo_γ_ini = replace(nomebase,".msh"=>"_γ_ini.dat")
    arquivo_γ_fin = replace(nomebase,".msh"=>"_γ_opt.dat")

    # Verificações de entrada
    isempty(freqs) && error("Analise Harmonica:: freqs deve ser um vetor não vazio")
    isa(freqs,Vector{Float64}) || error("freqs deve ser um vetor de floats")
    isfile(mshfile) || error("Otim:: arquivo de entrada $mshfile não existe")
    nf = length(freqs)

    if isempty(vA)
       vA = ones(nf)
    else 
       length(vA)==nf || error("Otim:: dimensão de vA deve ser nf")
    end

    isfile(arquivo_yaml) || error("Otim:: arquivo de entrada $(arquivo_yaml) não existe")

    # Leitura da malha e YAML
    nn, coord, ne, connect, materials, nodes_open, velocities, nodes_pressure, pressures, damping, nodes_probe, nodes_target, elements_fixed, values_fixed, centroides = LSound.Parsemsh_Daniele(mshfile)
    elements_design = setdiff(1:ne,sort!(elements_fixed))
    nvp = length(elements_design) 

    raio_filtro, niter, ϵ1, ϵ2,  vf, Past, fatorcv, μ = Le_YAML(arquivo_yaml)
     
    isempty(nodes_target) && error("Otim:: nodes_target deve ter ao menos um nó informado")
    sort!(nodes_target)
    isempty(materials) && error("Analise:: at least one material is necessary")

    # Pré-processamento
    if !verifica_derivada
         println("Determinando a vizinhança para um raio de $(raio_filtro)")
         vizinhos, pesos = Vizinhanca(ne,centroides,raio_filtro,elements_design)
    end 
    neighedge = NeighborEdges(ne,connect,elements_design)
        
    # Inicialização
    println("Inicializando o vetor de variáveis de projeto")
    γ = zeros(ne)
    Fix_γ!(γ,elements_fixed,values_fixed)
    writedlm(arquivo_γ_ini,γ)

    nodes_mask = sort(vcat(nodes_open,nodes_pressure))
    livres = setdiff(collect(1:nn),nodes_mask)
   
    # Inicializa os arquivos de saída do gmsh
    etype = connect[:,1]
    Lgmsh_export_init(arquivo_pos,nn,ne,coord,etype,connect[:,3:end])
    Lgmsh_export_init(arquivo_pos_freq,nn,ne,coord,etype,connect[:,3:end])
    Lgmsh_export_element_scalar(arquivo_pos,γ,"Iter 0")
   
    # Sweep Inicial
    MP,_ =  Sweep(nn,ne,coord,connect,γ,fρ,fκ,μ, freqs,livres,velocities,pressures)
    for i=1:nf
      f = freqs[i]
      Lgmsh_export_nodal_scalar(arquivo_pos_freq,abs.(MP[:,i]),"Pressure in $f Hz [abs] - initial topology")
    end

    if verifica_derivada
       # ... (Mantido igual ao original) ...
       γ = rand(ne)
       MP,K,M,C = Sweep(nn,ne,coord,connect,γ,fρ,fκ,μ, freqs,livres,velocities,pressures)
       dΦ = Derivada(ne,nn,γ,connect,coord,K,M,C,livres,freqs,pressures,dfρ,dfκ,μ, nodes_target,MP,elements_design,vA) 
       dnum = Verifica_derivada(γ,nn,ne,coord,connect,fρ,fκ,μ,freqs,livres,velocities,pressures,nodes_target,elements_design,vA)
       rel = (dΦ.-dnum)./(dnum.+1E-12)
       Lgmsh_export_element_scalar(arquivo_pos,γ,"γ")
       Lgmsh_export_element_scalar(arquivo_pos,dΦ,"Analitica")
       Lgmsh_export_element_scalar(arquivo_pos,dnum,"Numerica")
       Lgmsh_export_element_scalar(arquivo_pos,rel,"relativa")
       FPerimiter(γ) =  Perimiter(γ, neighedge, elements_design)
       dPnum = df(γ,FPerimiter,elements_design)
       dP = dPerimiter(ne, γ, neighedge, elements_design)
       return dΦ, dnum, dP, dPnum 
    end

    # --------------------------------------------------------------------------
    # Loop Principal de Otimização
    # --------------------------------------------------------------------------
    V = Volumes(ne,connect,coord)
    volume_full_projeto = sum(V[elements_design])
    Vast = vf*volume_full_projeto

    historico_V   = Float64[] 
    historico_SLP = Float64[] 
    historico_P   = Float64[] 

    # --- SETUP DA ESTABILIDADE ---
    cv_min = 1.0001 / nvp  
    n_warmup = 0 # Recomendo manter 5 para estabilizar
    cv_warmup_cap = 0.05 

    cv_current = fatorcv 

    for iter = 1:niter

        # --- LÓGICA DO WARM-UP ---
        is_warmup = iter <= n_warmup
        if is_warmup
            if cv_current > cv_warmup_cap
                cv_current = cv_warmup_cap
                println("--- WARM-UP (Iter $iter): Limitando cv a $(cv_warmup_cap)")
            end
        end

        # Identificação de elementos Ar/Sólido
        elements_air = Int64[]
        elements_solid = Int64[]
        one_air = Float64[]
        one_solid = Float64[]

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

        volume_atual = sum(γ[elements_design].*V[elements_design])
        push!(historico_V , volume_atual)

        # Análise FEM
        MP,K,M,C =  Sweep(nn,ne,coord,connect,γ,fρ,fκ,μ,freqs,livres,velocities,pressures)
        objetivo = Objetivo(MP,nodes_target,vA)
        push!(historico_SLP, objetivo)

        perimetro = Perimiter(γ, neighedge, elements_design)
        push!(historico_P, perimetro)
        
        for i=1:nf
         f = freqs[i]
         Lgmsh_export_nodal_scalar(arquivo_pos_freq,abs.(MP[:,i]),"Pressure in $f Hz [abs] $iter")
        end

        println("Iteração        ", iter)
        println("Objetivo        ", objetivo)
        println("Perimetro       ", perimetro)
        println("Past            ", Past)
        println("Volume atual    ", volume_atual)
        println("Volume target   ", Vast)
        println()

        # Sensibilidade
        dΦ = Derivada(ne,nn,γ,connect,coord,K,M,C,livres,freqs,pressures,dfρ,dfκ,μ,nodes_target,MP,elements_design,vA) 
        Lgmsh_export_element_scalar(arquivo_pos,dΦ,"dΦ")  

        # --- CRITICAL CHANGE 1: DISABLE FILTER FOR OPTIMIZATION ---
        # The filter encourages super-nucleation. We want RAW physics.
        # dΦf =  Filtro(vizinhos,pesos,dΦ,elements_design)
        # c = dΦf[elements_design]
        
        c = dΦ[elements_design] # Use Raw Sensitivities

        # Restrições Lineares
        ΔV = Vast - volume_atual
        b = [ΔV]
        println("Volume linearizado   ", b[1])
        A = vcat(V[elements_design]')

        if  perimetro > 0 && Past>0
            ΔP = Past - perimetro
            b = vcat(b, ΔP)
            println("Perimetro linearizado   ", ΔP)
            dP = dPerimiter(ne, γ, neighedge, elements_design)
            Lgmsh_export_element_scalar(arquivo_pos,dP,"dP")  
            A = vcat(A,transpose(dP[elements_design]))
        end
       
        # -----------------------------------------------------------------
        # TRUST REGION LOOP
        # -----------------------------------------------------------------
        step_accepted = false
        n_reject = 0
        flag_podrao = false

        while !step_accepted
            
            A_temp = copy(A)
            b_temp = copy(b)

            if !isempty(elements_air)
                g_air = ceil(cv_current * length(elements_design))
                b_temp = vcat(b_temp, g_air)
                A_temp = vcat(A_temp, transpose(one_air))
            end

            if !isempty(elements_solid)
                g_solid  = ceil(cv_current * length(elements_design)) 
                b_temp = vcat(b_temp, g_solid)
                A_temp = vcat(A_temp, -transpose(one_solid))
            end

            # Solve ILP
            Δγ = LP(c, A_temp, b_temp, γ[elements_design])

            # Predição
            pred_reduction = -sum(c .* Δγ)

            if abs(pred_reduction) < 1e-8
                println("Convergência local atingida.")
                step_accepted = true
                break
            end

            # Trial
            γ_trial = copy(γ)
            γ_trial[elements_design] .+= Δγ
            for ele in elements_design
                γ_trial[ele] = round(γ_trial[ele])
            end

            # --- CRITICAL CHANGE 2: GEOMETRIC FEASIBILITY CHECK ---
            # Antes de gastar tempo com a física (Sweep), verificamos a geometria REAL.
            # O ILP acha que nucleação é de graça (grad=0), mas aqui vemos o custo real.
            P_trial = Perimiter(γ_trial, neighedge, elements_design)
            
            # Tolerância pequena para erro numérico (1e-4)
            geom_violation = (Past > 0) && (P_trial > (Past + 1e-4))

            if geom_violation
                # REJEITADO POR GEOMETRIA
                println("    -> Passo REJEITADO (Geometria). P_real ($P_trial) > Past ($Past). Reduzindo cv.")
                cv_current = max(cv_min, cv_current * 0.5)
                n_reject += 1
            else
                # Se a geometria passou, checamos a FÍSICA
                MP_trial, _, _, _ = Sweep(nn,ne,coord,connect,γ_trial,fρ,fκ,μ,freqs,livres,velocities,pressures)
                obj_trial = Objetivo(MP_trial,nodes_target,vA)

                actual_reduction = objetivo - obj_trial
                R = actual_reduction / pred_reduction

                println("--- Trust Region Eval: Pred = $(round(pred_reduction, digits=4)), Act = $(round(actual_reduction, digits=4)), R = $(round(R, digits=2)), cv = $(round(cv_current, digits=4))")

                if R < 0.0 
                    # REJEITADO POR FÍSICA
                    println("    -> Passo REJEITADO (Física). SPL piorou.")
                    cv_current = max(cv_min, cv_current * 0.5) 
                    n_reject +=1
                else
                    # ACEITO
                    step_accepted = true
                    γ .= γ_trial 

                    if R > 0.7 
                        if is_warmup
                            println("    -> Passo ACEITO. Excelente, mantendo cv (Warm-up).")
                        else
                            println("    -> Passo ACEITO. Excelente, expandindo cv.")
                            cv_current = min(0.5, cv_current * 1.2)
                        end
                    elseif R > 0.3
                        println("    -> Passo ACEITO. Boa, mantendo cv.")
                    else
                        println("    -> Passo ACEITO. Pobre, reduzindo cv.")
                        cv_current = max(cv_min, cv_current * 0.5)
                    end
                end
            end

            if n_reject >= 10
                flag_podrao = true
                println("Trust Region falhou 10 vezes. Terminando.")
                break
            end
        end # While Trust Region

        if flag_podrao
            break
        end  

        Lgmsh_export_element_scalar(arquivo_pos,γ,"Iter $iter")

    end # Loop iter

    println("Final da otimização, executando a análise SWEEP na topologia otimizada")
    MP,_ =  Sweep(nn,ne,coord,connect,γ,fρ,fκ,μ,freqs,livres,velocities,pressures)
    for i=1:nf
      f = freqs[i]
      Lgmsh_export_nodal_scalar(arquivo_pos_freq,abs.(MP[:,i]),"Pressure in $f Hz [abs]")
    end
    writedlm(arquivo_γ_fin,γ)
    return historico_V, historico_SLP, historico_P

end # main_otim


#
# Rotina principal
#
"""
 Otim(meshfile::String,freqs=[])
"""
function Otim_ISLP_OLD(arquivo::String,freqs::Vector, vA::Vector;verifica_derivada=false)

   # Se o arquivo for um .geo, geramos um .msh utilizando a biblioteca do gmsh
   if occursin(".geo",arquivo)
      gmsh.initialize()
      gmsh.open(arquivo)
      gmsh.model.mesh.generate(2)
      mshfile = replace(arquivo,".geo"=>".msh")
      gmsh.write(mshfile)
   else 
      mshfile = arquivo
   end

    # Define os nomes dos arquivos
    arquivo_yaml = replace(mshfile,".msh"=>".yaml")
    nomebase = basename(mshfile)
    arquivo_pos       = replace(nomebase,".msh"=>".pos")
    arquivo_pos_freq  = replace(nomebase,".msh"=>"_freq.pos")

    arquivo_γ_ini = replace(nomebase,".msh"=>"_γ_ini.dat")
    arquivo_γ_fin = replace(nomebase,".msh"=>"_γ_opt.dat")

    # Verificações de entrada
    isempty(freqs) && error("Analise Harmonica:: freqs deve ser um vetor não vazio")
    isa(freqs,Vector{Float64}) || error("freqs deve ser um vetor de floats")
    isfile(mshfile) || error("Otim:: arquivo de entrada $mshfile não existe")
    nf = length(freqs)

    if isempty(vA)
       vA = ones(nf)
    else 
       length(vA)==nf || error("Otim:: dimensão de vA deve ser nf")
    end

    isfile(arquivo_yaml) || error("Otim:: arquivo de entrada $(arquivo_yaml) não existe")

    # Leitura da malha e YAML
    nn, coord, ne, connect, materials, nodes_open, velocities, nodes_pressure, pressures, damping, nodes_probe, nodes_target, elements_fixed, values_fixed, centroides = LSound.Parsemsh_Daniele(mshfile)
    elements_design = setdiff(1:ne,sort!(elements_fixed))
    nvp = length(elements_design) 

    raio_filtro, niter, ϵ1, ϵ2,  vf, Past, fatorcv, μ = Le_YAML(arquivo_yaml)
     
    isempty(nodes_target) && error("Otim:: nodes_target deve ter ao menos um nó informado")
    sort!(nodes_target)
    isempty(materials) && error("Analise:: at least one material is necessary")

    # Pré-processamento
    if !verifica_derivada
         println("Determinando a vizinhança para um raio de $(raio_filtro)")
         vizinhos, pesos = Vizinhanca(ne,centroides,raio_filtro,elements_design)
    end 
    neighedge = NeighborEdges(ne,connect,elements_design)
        
    # Inicialização
    println("Inicializando o vetor de variáveis de projeto")
    γ = zeros(ne)
    Fix_γ!(γ,elements_fixed,values_fixed)
    writedlm(arquivo_γ_ini,γ)

    nodes_mask = sort(vcat(nodes_open,nodes_pressure))
    livres = setdiff(collect(1:nn),nodes_mask)
   
    # Inicializa os arquivos de saída do gmsh
    etype = connect[:,1]
    # Para topologias
    Lgmsh_export_init(arquivo_pos,nn,ne,coord,etype,connect[:,3:end])
    # Para os modos e frequências
    Lgmsh_export_init(arquivo_pos_freq,nn,ne,coord,etype,connect[:,3:end])
    Lgmsh_export_element_scalar(arquivo_pos,γ,"Iter 0")
   
    # Sweep Inicial
    MP,_ =  Sweep(nn,ne,coord,connect,γ,fρ,fκ,μ, freqs,livres,velocities,pressures)
    for i=1:nf
      f = freqs[i]
      Lgmsh_export_nodal_scalar(arquivo_pos_freq,abs.(MP[:,i]),"Pressure in $f Hz [abs] - initial topology")
    end

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

    # --------------------------------------------------------------------------
    # Loop Principal de Otimização
    # --------------------------------------------------------------------------
    V = Volumes(ne,connect,coord)
    volume_full_projeto = sum(V[elements_design])
    Vast = vf*volume_full_projeto

    historico_V   = Float64[] 
    historico_SLP = Float64[] 
    historico_P   = Float64[] 

    # --- SETUP DA ESTABILIDADE ---
    # 1. Limite Inferior Dinâmico (1 elemento)
    cv_min = 1.0001 / nvp  
    
    # 2. Configuração do Warm-up
    n_warmup = 0         # Número de iterações de aquecimento
    cv_warmup_cap = 0.05 # Teto de 5% durante o aquecimento

    # Inicializa fator atual
    cv_current = fatorcv 

    for iter = 1:niter

        # --- LÓGICA DO WARM-UP ---
        is_warmup = iter <= n_warmup
        if is_warmup
            # Durante o warm-up, forçamos o passo a ser pequeno para estabilizar a topologia
            if cv_current > cv_warmup_cap
                cv_current = cv_warmup_cap
                println("--- WARM-UP (Iter $iter): Limitando cv a $(cv_warmup_cap)")
            end
        end

        # Identificação de elementos Ar/Sólido
        elements_air = Int64[]
        elements_solid = Int64[]
        one_air = Float64[]
        one_solid = Float64[]

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

        volume_atual = sum(γ[elements_design].*V[elements_design])
        push!(historico_V , volume_atual)

        # Análise FEM
        MP,K,M,C =  Sweep(nn,ne,coord,connect,γ,fρ,fκ,μ,freqs,livres,velocities,pressures)
        objetivo = Objetivo(MP,nodes_target,vA)
        push!(historico_SLP, objetivo)

        perimetro = Perimiter(γ, neighedge, elements_design)
        push!(historico_P, perimetro)
        
        # Exporta
        for i=1:nf
         f = freqs[i]
         Lgmsh_export_nodal_scalar(arquivo_pos_freq,abs.(MP[:,i]),"Pressure in $f Hz [abs] $iter")
        end

        println("Iteração        ", iter)
        println("Objetivo        ", objetivo)
        println("Perimetro       ", perimetro)
        println("Past            ", Past)
        println("Volume atual    ", volume_atual)
        println("Volume target   ", Vast)
        println()

        # Sensibilidade e Filtro
        dΦ = Derivada(ne,nn,γ,connect,coord,K,M,C,livres,freqs,pressures,dfρ,dfκ,μ,nodes_target,MP,elements_design,vA) 
        Lgmsh_export_element_scalar(arquivo_pos,dΦ,"dΦ")  

        dΦf =  Filtro(vizinhos,pesos,dΦ,elements_design)
        Lgmsh_export_element_scalar(arquivo_pos,dΦf,"Filtrada")
        c = dΦf[elements_design]

        # Restrições Lineares
        ΔV = Vast - volume_atual
        b = [ΔV]
        println("Volume linearizado   ", b[1])
        A = vcat(V[elements_design]')

        if  perimetro > 0 && Past>0
            ΔP = Past - perimetro
            b = vcat(b, ΔP)
            println("Perimetro linearizado   ", ΔP)
            dP = dPerimiter(ne, γ, neighedge, elements_design)
            Lgmsh_export_element_scalar(arquivo_pos,dP,"dP")  
            A = vcat(A,transpose(dP[elements_design]))
        end
       
        # -----------------------------------------------------------------
        # TRUST REGION LOOP
        # -----------------------------------------------------------------
        step_accepted = false
        n_reject = 0
        flag_podrao = false

        while !step_accepted
            
            A_temp = copy(A)
            b_temp = copy(b)

            if !isempty(elements_air)
                g_air = ceil(cv_current * length(elements_design))
                b_temp = vcat(b_temp, g_air)
                A_temp = vcat(A_temp, transpose(one_air))
            end

            if !isempty(elements_solid)
                g_solid  = ceil(cv_current * length(elements_design)) 
                b_temp = vcat(b_temp, g_solid)
                A_temp = vcat(A_temp, -transpose(one_solid))
            end

            # Solve ILP
            Δγ = LP(c, A_temp, b_temp, γ[elements_design])

            # Predição
            pred_reduction = -sum(c .* Δγ)

            if abs(pred_reduction) < 1e-8
                println("Convergência local atingida.")
                step_accepted = true
                break
            end

            # Trial
            γ_trial = copy(γ)
            γ_trial[elements_design] .+= Δγ
            for ele in elements_design
                γ_trial[ele] = round(γ_trial[ele])
            end

            # Check Physics
            MP_trial, _, _, _ = Sweep(nn,ne,coord,connect,γ_trial,fρ,fκ,μ,freqs,livres,velocities,pressures)
            obj_trial = Objetivo(MP_trial,nodes_target,vA)

            actual_reduction = objetivo - obj_trial
            R = actual_reduction / pred_reduction

            println("--- Trust Region Eval: Pred = $(round(pred_reduction, digits=4)), Act = $(round(actual_reduction, digits=4)), R = $(round(R, digits=2)), cv = $(round(cv_current, digits=4))")

            if R < 0.0 
                # REJEITADO
                println("    -> Passo REJEITADO. Reduzindo limites.")
                cv_current = max(cv_min, cv_current * 0.5) 
                n_reject +=1
                if n_reject >= 10
                    flag_podrao = true
                    break
                end
            else
                # ACEITO
                step_accepted = true
                γ .= γ_trial 

                if R > 0.7 
                    if is_warmup
                        println("    -> Passo ACEITO. Previsão excelente, mas MANTENDO cv (Warm-up).")
                        # Não expande durante o warm-up para evitar overshoot
                    else
                        println("    -> Passo ACEITO. Previsão excelente, expandindo cv.")
                        cv_current = min(0.5, cv_current * 1.2)
                    end
                elseif R > 0.3
                    println("    -> Passo ACEITO. Previsão boa, mantendo cv.")
                else
                    println("    -> Passo ACEITO. Previsão pobre, reduzindo cv.")
                    cv_current = max(cv_min, cv_current * 0.5)
                end
            end
        end # While Trust Region

        if flag_podrao
            println("Trust Region falhou 10 vezes. Terminando.")
            break
        end  

        Lgmsh_export_element_scalar(arquivo_pos,γ,"Iter $iter")

    end # Loop iter

    println("Final da otimização, executando a análise SWEEP na topologia otimizada")
    MP,_ =  Sweep(nn,ne,coord,connect,γ,fρ,fκ,μ,freqs,livres,velocities,pressures)
    for i=1:nf
      f = freqs[i]
      Lgmsh_export_nodal_scalar(arquivo_pos_freq,abs.(MP[:,i]),"Pressure in $f Hz [abs]")
    end
    writedlm(arquivo_γ_fin,γ)
    return historico_V, historico_SLP, historico_P

end # main_otim