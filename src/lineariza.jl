#
# Lineariza restrições
#
function Lineariza_Restricoes(V, elements_design, Vast, volume_atual, perimetro_atual, Past, ne, γ, neighedge, map_global_local)
    
    nvp = length(elements_design)
    
    # Restrições Globais (Volume)
    ΔV = Vast - volume_atual
    b_list = [ΔV]
    
    # A matriz de volume precisa de zeros para as variáveis de folga
    # Estrutura: [Grad_Vol (1xNVP) | Zeros (1xNVP)]
    A_vol = vcat(V[elements_design]')
    A_vol_ext = hcat(A_vol, zeros(1, nvp)) 
    A_list = A_vol_ext

    # Restrição de Perímetro (se ativo)
    if perimetro_atual > 0 && Past > 0
        ΔP = Past - perimetro_atual
        b_list = vcat(b_list, ΔP)
        
        dP = dPerimiter(ne, γ, neighedge, elements_design)
        A_per = transpose(dP[elements_design])
        A_per_ext = hcat(A_per, zeros(1, nvp)) # Padding zeros para slacks
        
        A_list = vcat(A_list, A_per_ext)
    end

    # 3. Restrições Topológicas Locais (No-Island / No-Hole)
    I_topo = Int[]
    J_topo = Int[]
    V_topo = Float64[]
    b_topo = Float64[]
    
    row_idx = 1
    
    for (idx_e, e_global) in enumerate(elements_design)
        vizinhos = neighedge[e_global]
        γ_e = γ[e_global]
        sum_γ_viz = sum(γ[v] for v in vizinhos)
        num_viz = length(vizinhos)

        # No-Island: Δγ_e - Σ Δγ_j - ξ_e <= RHS
        push!(I_topo, row_idx); push!(J_topo, idx_e); push!(V_topo, 1.0)
        push!(I_topo, row_idx); push!(J_topo, nvp + idx_e); push!(V_topo, -1.0) # Slack
        for v in vizinhos
            if haskey(map_global_local, v)
                push!(I_topo, row_idx); push!(J_topo, map_global_local[v]); push!(V_topo, -1.0)
            end
        end
        push!(b_topo, sum_γ_viz - γ_e)
        row_idx += 1
        
        # No-Hole: Σ Δγ_j - Δγ_e - ξ_e <= RHS
        push!(I_topo, row_idx); push!(J_topo, idx_e); push!(V_topo, -1.0)
        push!(I_topo, row_idx); push!(J_topo, nvp + idx_e); push!(V_topo, -1.0) # Slack
        for v in vizinhos
            if haskey(map_global_local, v)
                push!(I_topo, row_idx); push!(J_topo, map_global_local[v]); push!(V_topo, 1.0)
            end
        end
        push!(b_topo, (num_viz - 1) - (sum_γ_viz - γ_e))
        row_idx += 1
    end
    
    A_topo = sparse(I_topo, J_topo, V_topo, 2*nvp, 2*nvp)
    
    # Concatena Global (denso convertido para sparse) com Topológico (sparse)
    A_final = vcat(sparse(A_list), A_topo)
    b_final = vcat(b_list, b_topo)
    
    return A_final, b_final
end