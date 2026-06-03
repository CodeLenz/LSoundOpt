#
# Lineariza as restrições usadas no MILP.
#
# A matriz final tem duas partes:
#   1) restrições globais de volume/perímetro
#   2) restrições locais para evitar ilhas e buracos
#
function Lineariza_Restricoes(V, elements_design, Vast, volume_atual, perimetro_atual,
                              Past, ne, γ, neighedge, map_global_local;
                              restricao_volume=true)

    # Número de variáveis de projeto
    nvp = length(elements_design)

    # Restrições globais opcionais
    A_list = Matrix{Float64}(undef, 0, 2*nvp)
    b_list = Float64[]

    # Restrição de volume
    if restricao_volume

        ΔV = Vast - volume_atual
        push!(b_list, ΔV)

        # A matriz de volume precisa de zeros nas colunas de folga
        A_vol = vcat(V[elements_design]')
        A_vol_ext = hcat(A_vol, zeros(1, nvp))
        A_list = vcat(A_list, A_vol_ext)

    end

    # Restrição de perímetro
    if perimetro_atual > 0 && Past > 0

        ΔP = Past - perimetro_atual
        push!(b_list, ΔP)

        # Derivada do perímetro nos elementos de projeto
        dP = dPerimiter(ne, γ, neighedge, elements_design)
        A_per = transpose(dP[elements_design])
        A_per_ext = hcat(A_per, zeros(1, nvp))

        A_list = vcat(A_list, A_per_ext)

    end

    # Restrições topológicas locais
    I_topo = Int[]
    J_topo = Int[]
    V_topo = Float64[]
    b_topo = Float64[]

    row_idx = 1

    # Loop nos elementos de projeto
    for (idx_e, e_global) in enumerate(elements_design)

        # Vizinhos por aresta
        vizinhos = neighedge[e_global]
        γ_e = γ[e_global]
        sum_γ_viz = sum(γ[v] for v in vizinhos)
        num_viz = length(vizinhos)

        # No-Island: impede sólido isolado
        push!(I_topo, row_idx); push!(J_topo, idx_e); push!(V_topo, 1.0)
        push!(I_topo, row_idx); push!(J_topo, nvp + idx_e); push!(V_topo, -1.0)

        for v in vizinhos
            if haskey(map_global_local, v)
                push!(I_topo, row_idx)
                push!(J_topo, map_global_local[v])
                push!(V_topo, -1.0)
            end
        end

        push!(b_topo, sum_γ_viz - γ_e)
        row_idx += 1

        # No-Hole: impede vazio isolado
        push!(I_topo, row_idx); push!(J_topo, idx_e); push!(V_topo, -1.0)
        push!(I_topo, row_idx); push!(J_topo, nvp + idx_e); push!(V_topo, -1.0)

        for v in vizinhos
            if haskey(map_global_local, v)
                push!(I_topo, row_idx)
                push!(J_topo, map_global_local[v])
                push!(V_topo, 1.0)
            end
        end

        push!(b_topo, (num_viz - 1) - (sum_γ_viz - γ_e))
        row_idx += 1

    end

    # Monta a matriz esparsa das restrições locais
    A_topo = sparse(I_topo, J_topo, V_topo, 2*nvp, 2*nvp)

    # Concatena restrições globais e locais
    A_final = vcat(sparse(A_list), A_topo)
    b_final = vcat(b_list, b_topo)

    # Retorna o sistema A*delta <= b
    return A_final, b_final

end
