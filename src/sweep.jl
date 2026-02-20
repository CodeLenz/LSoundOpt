function Sweep(nn,ne,coord,connect,γ,fρ,fκ,μ,freqs,livres,velocities,pressures::Vector)

    # 1. Monta globais (Custo único)
    K,M,C = Monta_KMC_param(ne,coord,connect,γ,fρ,fκ,μ)
    
    # 2. OTIMIZAÇÃO: Recorta as matrizes para os DOFs livres APENAS UMA VEZ
    # Isso evita alocar memória para nós fixos dentro do loop
    K_red = K[livres, livres]
    M_red = M[livres, livres]
    C_red = C[livres, livres]

    nf = length(freqs)
    U = zeros(ComplexF64,nn)
    MP = zeros(ComplexF64,nn,nf)
    P = zeros(ComplexF64,nn)    
    P_red = zeros(ComplexF64, length(livres)) 

    contador = 1
    for f in freqs
        ω = 2*pi*f

        # 3. Monta Kd diretamente reduzida (Muito mais rápido e leve)
        # Como K, M e C têm a mesma estrutura de esparsidade, o Julia é eficiente aqui
        Kd_red = K_red .+ (im*ω).*C_red .- (ω^2).*M_red

        # Atualiza vetor de forças global (dependente da frequência)
        fill!(P, 0.0) # Zera para não acumular lixo
        LSound.Vetor_P!(0.0,velocities,coord,connect,P,ω=ω)
        
        # Recorta vetor de forças
        P_red .= P[livres]
    
        # Soluciona o sistema reduzido
        U_red = Kd_red \ P_red

        # Reconstrói o vetor global U
        U[livres] .= U_red

        # Armazena na matriz de resposta
        MP[:,contador] .= U
        contador += 1
        
    end

    return MP, K, M, C
end