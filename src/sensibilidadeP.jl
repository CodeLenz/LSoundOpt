# ===================================================================================
# Calcula a derivada da função objetivo P-Norm Normalizada
#
function DerivadaP(ne,nn,γ::Vector{T0},connect::Matrix{T1},coord::Matrix{T0},
                  K::AbstractMatrix{T0},M::AbstractMatrix{T0},C::AbstractMatrix{T0},
                  livres::Vector{T1},freqs::Vector{T0},
                  pressures::Vector, 
                  fρ::Function, fκ::Function,
                  dfρ::Function, dfκ::Function,μ::T0,
                  nodes_target::Vector{T1},MP::Matrix{T2},
                  elements_design::Vector,A::Vector,
                  P_ref::Vector; 
                  p_norm=8.0) where {T0,T1,T2}

    d = zeros(ne)
    Nf = length(freqs)
    nt = length(nodes_target)

    # --- PRÉ-CÁLCULO DOS PESOS (η_n) ---
    Rn_values = zeros(Float64, Nf)
    P2avg_values = zeros(Float64, Nf)

    for i = 1:Nf
        P_col = MP[:, i]
        P2_avg = sum(abs2.(P_col[nodes_target])) / nt
        P2avg_values[i] = P2_avg
        
        # Razão Normalizada: P_atual / P_ref
        Rn_values[i] = A[i] * (P2_avg / P_ref[i])
    end

    # Termo da Média Generalizada (M_p)
    soma_potencias = sum(Rn_values .^ p_norm)
    Mp = soma_potencias / Nf
    
    # Constante global da derivada
    const_global = 10.0 / (log(10.0) * Nf * Mp)

    # Alocações
    λn = zeros(T2,nn)
    P = zeros(T2,nn)
    Fn = zeros(T2,nn)

    # Loop pelas frequências
    coluna = 1
    for f in freqs
        
        ωn = 2*pi*f
        
        Rn = Rn_values[coluna]
        P2avg = P2avg_values[coluna]
        
        if P2avg < 1e-20; coluna += 1; continue; end

        # --- CÁLCULO DO PESO ESPECTRAL η_n ---
        # Formula do texto Eq (46): eta = (Rn^p) / (Nf * Mp * P_atual)
        # Nota: P_ref já está "embutido" dentro de Rn^p
        
        # Simplificação numérica:
        # Fator de escala = const_global * (Rn)^p / P_atual
        # Isso equivale a derivada da log-sum-exp ponderada
        
        fator_espectral = const_global * (Rn^p_norm) / P2avg
        
        # --- SOLUÇÃO ADJUNTA ---
        P .= MP[:,coluna]
        Kd = K[livres,livres] .+ im*ωn*C[livres,livres] .- (ωn^2)*M[livres,livres]
        
        # Vetor de força adjunto base
        Fn .= F_adj(nodes_target, P)

        # Aplica a escala: fator_espectral * (1/nt)
        scale_factor = fator_espectral * (1.0/nt) 
        Fn .= Fn * scale_factor

        λn[livres] .= Kd \ Fn[livres]

        # --- SENSIBILIDADE DO ELEMENTO ---
        for ele in elements_design
            etype = connect[ele,1]
            nos, X = LSound.Nos_Coordenadas(ele,etype,coord,connect)
            pe = P[nos]
            λe = λn[nos]
            γe = γ[ele]

            dKe, dMe, dCe = Derivada_KMC(etype,γe,fρ,fκ,dfρ,dfκ,μ,X)
            dKde = dKe .+ im*ωn*dCe .- dMe*ωn^2  

            # Acumula: d += 2 * Re(λ^T * dK * p)
            d[ele] += 2.0 * real(transpose(λe) * dKde * pe) 
        end 

        coluna += 1
    end 

    return d
end