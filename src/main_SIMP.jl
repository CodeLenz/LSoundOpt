"""
    Otim_SIMP(arquivo::String, freqs::Vector, vA::Vector;
              p_inicial=1.0, p_final=4.0, p_passo=0.5,
              verifica_derivada=false,
              restricao_volume=true,
              escala_mma=true,
              perturbacao_inicial=0.05,
              alvo_grad_mma=1.0,
              escala_mma_max=1.0e8,
              offset_obj_mma=50.0,
              algoritmo=:SLP)

Versão SIMP (contínua, via MMA/NLopt) do otimizador acústico reativo.

Reusa integralmente as rotinas da formulação binária:
  - Sweep!            (resposta harmônica multifrequência)
  - Objetivo / SPLn   (objetivo SPL médio)
  - Derivada          (sensibilidade adjunta — agnóstica à parametrização)
  - Filtro            (filtro espacial de sensibilidade)
  - Volumes           (volume/área elementar)

A ÚNICA diferença em relação a Otim_ISLP é o bloco de atualização das
variáveis de projeto: no lugar do SILP + Trust-Region, usa-se MMA
(NLopt, LD_MMA) com γ contínuo em [0,1] e um *continuation* no expoente
de penalização p (p_inicial → p_final em passos p_passo).

A parametrização SIMP é injetada via Cria_Param_SIMP (param_simp.jl),
de modo que p=1 recupera exatamente a interpolação linear (Dühring).

Objetivo da comparação (resposta ao Revisor #2): demonstrar que o SIMP,
nas MESMAS quantidades interpoladas, produz elementos cinza sem
interpretação física como mistura ar/sólido rígido — em contraste com
o resultado preto-e-branco da formulação binária.
"""
function Otim_SIMP(arquivo::String, freqs::Vector, vA::Vector;
                   p_inicial=1.0, p_final=4.0, p_passo=0.5,
                   verifica_derivada=false,
                   restricao_volume=true,
                   escala_mma=true,
                   perturbacao_inicial=0.05,
                   alvo_grad_mma=1.0,
                   escala_mma_max=1.0e8,
                   offset_obj_mma=50.0,
                   algoritmo=:SLP)

    # --- SETUP INICIAL (idêntico a Otim_ISLP) ---
    mshfile, arquivos_saida = Setup_Arquivos(arquivo)
    arquivo_pos, arquivo_pos_freq, arquivo_data_opt, arquivo_γ_ini, arquivo_γ_fin = arquivos_saida

    # Pesos do SPL por frequência
    if isempty(vA)
        vA = ones(length(freqs))
    end
    nf = length(freqs)

    # Leitura da malha
    nn, coord, ne, connect, materials, nodes_open, velocities, nodes_pressure, pressures, damping, nodes_probe, nodes_target, elements_fixed, values_fixed, centroides = LSound.Parsemsh_Daniele(mshfile)

    # Elementos de projeto
    elements_design = setdiff(1:ne, sort!(elements_fixed))
    nvp = length(elements_design)
    map_global_local = Dict{Int, Int}()
    for (i, ele) in enumerate(elements_design)
        map_global_local[ele] = i
    end

    # Leitura do YAML (raio, niter, vf, perímetro, μ, fatorcv)
    # OBS: fatorcv e o perímetro-limite não são usados no SIMP; mantidos
    #      apenas para compatibilidade de assinatura.
    raio_filtro, niter, vf, Past, μ, fatorcv = Le_YAML(replace(mshfile, ".msh" => ".yaml"))

    # Vizinhança/arestas e vizinhança do filtro
    neighedge = NeighborEdges(ne, connect, elements_design)
    if !verifica_derivada && raio_filtro > 0
        println("Determinando a vizinhança para um raio de $(raio_filtro)")
        vizinhos, pesos, _, _ = Vizinhanca(ne, centroides, raio_filtro, elements_design)
    end

    # ---- Parametrização SIMP inicial (continuation começa em p_inicial) ----
    par = Cria_Param_SIMP(p_inicial)
    fρ, fκ, dfρ, dfκ = par.fρ, par.fκ, par.dfρ, par.dfκ

    # Variáveis de projeto iniciais.
    # SIMP precisa de ponto interior viável. Começamos uniforme em vf,
    # mas com uma pequena perturbação aleatória: com p alto e γ0 uniforme,
    # a simetria deixa o gradiente quase degenerado no arranque e o MMA
    # propõe passos minúsculos. A perturbação quebra a simetria e dá ao
    # MMA uma direção de descida não-trivial já na primeira avaliação.
    # (Mantida pequena para não enviesar a topologia final.)
    println("Inicializando topologia (SIMP, γ0 = 0")#vf + perturbação)...")
    γ = zeros(ne)
    #γ[elements_design] .= vf .+ perturbacao_inicial .* (rand(nvp) .- 0.5)
    clamp!(γ, 0.0, 1.0)
    Fix_γ!(γ, elements_fixed, values_fixed)
    writedlm(arquivo_γ_ini, γ)

    # Graus de liberdade livres
    nodes_mask = sort(vcat(nodes_open, nodes_pressure))
    livres = setdiff(collect(1:nn), nodes_mask)

    # Export inicial
    Lgmsh_export_init(arquivo_pos, nn, ne, coord, connect[:, 1], connect[:, 3:end])
    Lgmsh_export_init(arquivo_pos_freq, nn, ne, coord, connect[:, 1], connect[:, 3:end])
    Lgmsh_export_element_scalar(arquivo_pos, γ, "Iter 0")

    # Matriz de resposta (reuso entre sweeps)
    MP = zeros(ComplexF64, nn, nf)
    Kd_factors = Cria_Cache_Kd(nf)

    # Sweep inicial
    K_ini, M_ini, C_ini = Sweep!(nn, ne, coord, connect, γ, fρ, fκ, μ,
                                 freqs, livres, velocities, pressures, MP,
                                 Kd_factors)
    Φ_ini = Objetivo(MP, nodes_target, vA)
    for i = 1:nf
        Lgmsh_export_nodal_scalar(arquivo_pos_freq, abs.(MP[:, i]), "Pressure $(freqs[i]) Hz - Ini")
    end

    # --- Verificação de derivada (agora com a parametrização SIMP) ---
    if verifica_derivada
        println("--- Modo de Verificação de Derivada (SIMP, p=$(p_inicial)) ---")
        γ = rand(ne)
        K, M, C = Sweep!(nn, ne, coord, connect, γ, fρ, fκ, μ, freqs, livres, velocities, pressures, MP, Kd_factors)
        dΦ = Derivada(ne, nn, γ, connect, coord, K, M, C, livres, freqs, pressures, fρ, fκ, dfρ, dfκ, μ, nodes_target, MP, elements_design, vA; Kd_factors=Kd_factors)
        dnum = Verifica_derivada(γ, nn, ne, coord, connect, fρ, fκ, μ, freqs, livres, velocities, pressures, nodes_target, elements_design, vA)
        rel = (dΦ .- dnum) ./ (dnum .+ 1e-12)
        Lgmsh_export_element_scalar(arquivo_pos, γ, "γ_debug")
        Lgmsh_export_element_scalar(arquivo_pos, dΦ, "Analitica")
        Lgmsh_export_element_scalar(arquivo_pos, dnum, "Numerica")
        Lgmsh_export_element_scalar(arquivo_pos, rel, "relativa")
        return dΦ, dnum
    end

    ############################################################################
    #                         OTIMIZAÇÃO (MMA + continuation em p)
    ############################################################################

    # Volumes elementares e volume-limite
    V = Volumes(ne, connect, coord)
    Vast = vf * sum(V[elements_design])
    restricao_volume || println("Rodando SIMP sem restrição de volume; apenas os bounds 0 <= γ <= 1 serão impostos.")

    # Históricos
    hist = (V = Float64[], SLP = Float64[], P = Float64[])

    if algoritmo == :SLP
        println("\n========== Otimização SIMP via SLP contínuo ==========")
        cv_current = fatorcv
        iter_global = 0

        p_atual = p_inicial
        while p_atual <= p_final + 1e-9

            println("\n========== Estágio SIMP-SLP: p = $(p_atual) ==========")

            par = Cria_Param_SIMP(p_atual)
            fρ, fκ, dfρ, dfκ = par.fρ, par.fκ, par.dfρ, par.dfκ

            for iter = 1:niter
                iter_global += 1

                volume_atual = sum(γ[elements_design] .* V[elements_design])

                K, M, C = Sweep!(nn, ne, coord, connect, γ, fρ, fκ, μ,
                                 freqs, livres, velocities, pressures, MP,
                                 Kd_factors)

                Φ = Objetivo(MP, nodes_target, vA)
                perim = Perimiter(γ, neighedge, elements_design)

                push!(hist.V, volume_atual)
                push!(hist.SLP, Φ)
                push!(hist.P, perim)
                Lgmsh_export_element_scalar(arquivo_pos, γ, "Iter $(iter_global)")

                if restricao_volume
                    @printf("Iter %3d | SPL: %.4f | Perim: %.4f | Vol: %.4e (T: %.4e) | cv: %.4f\n",
                            iter_global, Φ, perim, volume_atual, Vast, cv_current)
                else
                    @printf("Iter %3d | SPL: %.4f | Perim: %.4f | Vol: %.4e (sem restr.) | cv: %.4f\n",
                            iter_global, Φ, perim, volume_atual, cv_current)
                end

                dΦ = Derivada(ne, nn, γ, connect, coord, K, M, C, livres,
                              freqs, pressures, fρ, fκ, dfρ, dfκ, μ,
                              nodes_target, MP, elements_design, vA;
                              Kd_factors=Kd_factors)
                if raio_filtro > 0
                    dΦ = Filtro(vizinhos, pesos, dΦ, elements_design)
                end
                c = dΦ[elements_design]

                A_global, b_global = Lineariza_Restricoes(V, elements_design, Vast,
                                                           volume_atual, perim, Past,
                                                           ne, γ, neighedge,
                                                           map_global_local;
                                                           restricao_volume=restricao_volume)

                fem_data = (nn, ne, coord, connect, fρ, fκ, μ, freqs, livres,
                            velocities, pressures, nodes_target, vA)

                γ_new, cv_current, step_accepted = Trust_Region_Loop_SIMP(
                    c, A_global, b_global, γ, elements_design, cv_current, Φ,
                    Past, neighedge, fem_data, MP)

                if !step_accepted
                    println("Trust Region SLP falhou. Terminando.")
                    break
                end

                if maximum(abs.(γ_new[elements_design] .- γ[elements_design])) < 1e-8
                    println("    -> SLP sem mudança efetiva. Encerrando estágio p=$(p_atual).")
                    break
                end

                γ .= γ_new
            end

            Lgmsh_export_element_scalar(arquivo_pos, γ, "p = $(p_atual)")
            p_atual += p_passo
        end

        println("\nFinalizando...")
        writedlm(arquivo_γ_fin, γ)
        Salva_Historico(arquivo_data_opt, hist, raio_filtro)

        Sweep!(nn, ne, coord, connect, γ, fρ, fκ, μ, freqs, livres, velocities, pressures, MP)
        for i = 1:nf
            Lgmsh_export_nodal_scalar(arquivo_pos_freq, abs.(MP[:, i]), "Pressure $(freqs[i]) Hz - Opt")
        end

        return hist.V, hist.SLP, hist.P
    elseif algoritmo != :MMA
        error("Otim_SIMP:: algoritmo deve ser :SLP ou :MMA")
    end

    # Vetor de design reduzido (só elementos de projeto) — espaço do MMA
    x = γ[elements_design]
    Vd = V[elements_design]               # volumes dos elementos de projeto
    escala_V_mma = escala_mma ? max(abs(Vast), eps(Float64)) : 1.0
    escala_Φ_mma = 1.0
    offset_Φ_mma = escala_mma ? offset_obj_mma : 0.0
    if escala_mma
        dΦ_ini = Derivada(ne, nn, γ, connect, coord, K_ini, M_ini, C_ini, livres,
                          freqs, pressures, fρ, fκ, dfρ, dfκ, μ,
                          nodes_target, MP, elements_design, vA;
                          Kd_factors=Kd_factors)
        if raio_filtro > 0
            dΦ_ini = Filtro(vizinhos, pesos, dΦ_ini, elements_design)
        end
        grad_ref_mma = max(maximum(abs.(dΦ_ini[elements_design])), eps(Float64))
        escala_por_objetivo = 100.0 / max(abs(Φ_ini), 1.0)
        escala_por_gradiente = alvo_grad_mma / grad_ref_mma
        escala_Φ_mma = min(max(escala_por_objetivo, escala_por_gradiente), escala_mma_max)
    end
    iter_global = Ref(0)                  # contador global de iterações
    escala_mma && @printf("Escalas MMA | objetivo: %.4e | volume: %.4e | Φ0_MMA: %.4e\n",
                          escala_Φ_mma, escala_V_mma, offset_Φ_mma)

    # =========================================================================
    # Callback do OBJETIVO para o NLopt.
    # Reconstrói γ a partir de x, faz o sweep, calcula Φ e dΦ (adjunto),
    # filtra a sensibilidade e devolve valor + grad (no espaço reduzido).
    # As closures capturam fρ/fκ/dfρ/dfκ atuais (atualizados a cada estágio p).
    # =========================================================================
    function faz_callbacks(fρ, fκ, dfρ, dfκ)

        function objetivo_mma!(x::Vector, grad::Vector)
            # Reconstrói γ global
            γ[elements_design] .= x

            # Resposta harmônica
            K, M, C = Sweep!(nn, ne, coord, connect, γ, fρ, fκ, μ,
                             freqs, livres, velocities, pressures, MP,
                             Kd_factors)

            # Objetivo
            Φ = Objetivo(MP, nodes_target, vA)

            # Gradiente (adjunto) — só se o NLopt pediu
            if length(grad) > 0
                dΦ = Derivada(ne, nn, γ, connect, coord, K, M, C, livres,
                              freqs, pressures, fρ, fκ, dfρ, dfκ, μ,
                              nodes_target, MP, elements_design, vA;
                              Kd_factors=Kd_factors)
                if raio_filtro > 0
                    dΦ = Filtro(vizinhos, pesos, dΦ, elements_design)
                end
                grad .= escala_Φ_mma .* dΦ[elements_design]
            end

            # Log + histórico
            iter_global[] += 1
            perim = Perimiter(γ, neighedge, elements_design)
            volat = sum(x .* Vd)
            push!(hist.V, volat); push!(hist.SLP, Φ); push!(hist.P, perim)
            Lgmsh_export_element_scalar(arquivo_pos, γ, "Iter $(iter_global[])")
            if restricao_volume
                @printf("Iter %3d | SPL: %.4f | Perim: %.4f | Vol: %.4e (T: %.4e)\n",
                        iter_global[], Φ, perim, volat, Vast)
            else
                @printf("Iter %3d | SPL: %.4f | Perim: %.4f | Vol: %.4e (sem restr.)\n",
                        iter_global[], Φ, perim, volat)
            end

            return offset_Φ_mma + escala_Φ_mma * (Φ - Φ_ini)
        end

        # Restrição de volume escalada:  g(x) = (Σ xᵢ Vᵢ − Vast) / Vast ≤ 0.
        function restricao_volume!(x::Vector, grad::Vector)
            if length(grad) > 0
                grad .= Vd ./ escala_V_mma
            end
            return (sum(x .* Vd) - Vast) / escala_V_mma
        end

        return objetivo_mma!, restricao_volume!
    end

    # =========================================================================
    # CONTINUATION no expoente de penalização p
    # =========================================================================
    p_atual = p_inicial
    while p_atual <= p_final + 1e-9

        println("\n========== Estágio SIMP: p = $(p_atual) ==========")

        # Atualiza a parametrização para o p deste estágio
        par = Cria_Param_SIMP(p_atual)
        fρ, fκ, dfρ, dfκ = par.fρ, par.fκ, par.dfρ, par.dfκ

        # (Re)cria os callbacks com a parametrização atual
        objetivo_mma!, restricao_volume! = faz_callbacks(fρ, fκ, dfρ, dfκ)

        # ---- Configura o problema MMA no NLopt ----
        opt = NLopt.Opt(:LD_MMA, nvp)
        opt.lower_bounds = zeros(nvp)
        opt.upper_bounds = ones(nvp)
        opt.maxeval = niter
        # Convergência pelo objetivo (SPL), não pelo passo em x.
        # O xtol_rel default (1e-4) parava o MMA logo na 1ª iteração
        # externa, quando o passo ainda é pequeno mas o objetivo não
        # convergiu. ftol_rel garante parada só quando o SPL estabiliza.
        opt.ftol_rel = 1e-6
        opt.xtol_rel = 0.0          # desliga o critério de passo
        opt.min_objective = objetivo_mma!
        if restricao_volume
            NLopt.inequality_constraint!(opt, restricao_volume!, 1e-8)
        else
            println("  -> Restrição de volume desligada neste estágio.")
        end

        # Otimiza este estágio (warm-start a partir do x corrente)
        (minf, xopt, ret) = NLopt.optimize(opt, x)
        x = xopt
        Φ_estagio = escala_mma ? Φ_ini + (minf - offset_Φ_mma) / escala_Φ_mma : minf
        println("  -> p=$(p_atual): Φ = $(Φ_estagio)  (status: $(ret), MMA=$(minf))")

        # Exporta a topologia ao fim de cada estágio
        γ[elements_design] .= x
        Lgmsh_export_element_scalar(arquivo_pos, γ, "p = $(p_atual)")

        # Próximo estágio
        p_atual += p_passo
    end

    # FINALIZAÇÃO
    println("\nFinalizando...")
    γ[elements_design] .= x
    writedlm(arquivo_γ_fin, γ)
    Salva_Historico(arquivo_data_opt, hist, raio_filtro)

    # Sweep final + export dos campos de pressão otimizados
    Sweep!(nn, ne, coord, connect, γ, fρ, fκ, μ, freqs, livres, velocities, pressures, MP)
    for i = 1:nf
        Lgmsh_export_nodal_scalar(arquivo_pos_freq, abs.(MP[:, i]), "Pressure $(freqs[i]) Hz - Opt")
    end

    return hist.V, hist.SLP, hist.P
end
