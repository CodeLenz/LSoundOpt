#
# Monta o vetor adjunto para uma frequência específica
#
# nodes_target -> vector of Target nodes   (vetor com os 'nós monitorados' )
# O vetor P foi obtido via Sweep.
#
function F_adj(nodes_target::Vector{T1},P::Vector{T2}) where {T1,T2}

    # Inicializa o vetor de carregamento adjunto
    F = similar(P) 
    
    # Zera o vetor
    fill!(F,zero(T2))

    # Somatório nas posições de nodes_target
    for p in nodes_target

        F[p] += -conj(P[p])

    end

    # Retorna F adjunto, sem os termos constantes que multiplicam
    return F 

end

# ===================================================================================
# Calcula a derivada da matriz de rigidez dinâmica de  um elemento
#
# et  ->  tipo de elemento
# γe  ->  variável de projeto do elemento
# fρ ->  função que parametriza a (inversa da) densidade 
# fκ ->  função que parametriza a  (inversa do) módulo de compressibilidade
# dfρ ->  função que parametriza a derivada da (inversa da) densidade 
# dfκ ->  função que parametriza a derivada da (inversa do)  módulo de compressibilidade
# X   ->  matriz com as coordenadas do elemento 
#
function Derivada_KMC(et,γe,fρ::Function,fκ::Function,dfρ::Function,dfκ::Function,μ,X::Array)


    # Monta as matrizes dos elementos sem parametrização 
    # ou seja K^0 e M^0
    if et==3
        Ke0, Me0 = LSound.KMe_bi4(1.0,1.0,X)
    elseif et==2
        Ke0, Me0 = LSound.KMe_tri3(1.0,1.0,X)
    elseif et==4
        Ke0, Me0 = LSound.KMe_tet4(1.0,1.0,X)
    elseif et==5
        Ke0, Me0 = LSound.KMe_hex8(1.0,1.0,X)
    elseif et==7
        Ke0, Me0 = LSound.KMe_pyr5(1.0,1.0,X)
    else
        error("Derivada_KM::Elemento não definido")
    end

    # Calcula a inversa de ρ e da inversa
    # de κ 
    iρ = fρ(γe)
    iκ = fκ(γe)

    # Calcula as derivadas da inversa de ρ e da inversa
    # de κ em relação à γe
    diρ = dfρ(γe)
    diκ = dfκ(γe)

    # Calcula as derivadas
    dKe = diρ*Ke0
    dMe = diκ*Me0
    dCe = (4/3)*μ*(iρ*diκ + iκ*diρ)*Ke0

    # Devolve as derivadas
    return dKe, dMe, dCe

end



# ===================================================================================
# Pré-calcula as matrizes de sensibilidade elementar (INDEPENDENTES DA FREQUÊNCIA)
# Retorna um vetor de NamedTuples com os dados prontos para uso rápido
#
function Precalcula_Dados_Elementos(γ, connect, coord, elements_design,
                                    fρ, fκ, dfρ, dfκ, μ)
    
    # Pré-aloca o vetor de dados
    dados_elem = Vector{NamedTuple{(:nos, :dKe, :dMe, :dCe), 
                        Tuple{Vector{Int}, Matrix{Float64}, Matrix{Float64}, Matrix{Float64}}}}(undef, length(elements_design))
    
    # Usa threads aqui também para acelerar o setup
    Threads.@threads for i in eachindex(elements_design)

        # Recupera od dados do elemento
        ele = elements_design[i]
        etype = connect[ele, 1]
        
        # Recupera nós e coordenadas (operação de I/O na malha)
        nos, X = LSound.Nos_Coordenadas(ele, etype, coord, connect)
        
        # Variável de projeto
        γe = γ[ele]
        
        # Calcula as matrizes (operação pesada de geometria)
        dKe, dMe, dCe = Derivada_KMC(etype, γe, fρ, fκ, dfρ, dfκ, μ, X)
        
        # Armazena
        dados_elem[i] = (nos=nos, dKe=dKe, dMe=dMe, dCe=dCe)

    end
    
    return dados_elem
end


function Derivada(ne, nn, γ::Vector{T0}, connect::Matrix{T1}, coord::Matrix{T0},
                  K::AbstractMatrix{T0}, M::AbstractMatrix{T0}, C::AbstractMatrix{T0},
                  livres::Vector{T1}, freqs::Vector{T0},
                  pressures::Vector, 
                  fρ::Function, fκ::Function,
                  dfρ::Function, dfκ::Function, μ::T0,
                  nodes_target::Vector{T1}, MP::Matrix{T2},
                  elements_design::Vector, A::Vector, p0=20E-6) where {T0, T1, T2}

    # 1. SETUP PRÉ-LOOP (Memória e Geometria)
    dados_elementos = Precalcula_Dados_Elementos(γ, connect, coord, elements_design, fρ, fκ, dfρ, dfκ, μ)

    d = zeros(ne)
    Nf = length(freqs)
    nt = length(nodes_target)
    const_log = 10.0 / log(10.0)

    # -----------------------------------------------------------
    # OTIMIZAÇÃO DE MATRIZES: Fatiamento fora do loop (HOISTING)
    # -----------------------------------------------------------
    # Isso evita alocar memória e buscar índices milhares de vezes
    K_livre = K[livres, livres]
    M_livre = M[livres, livres]
    C_livre = C[livres, livres]
    
    # Alocações auxiliares
    λn = zeros(ComplexF64, nn)
    Fn = zeros(ComplexF64, nn)

    # Loop pelas frequências
    for (coluna, f) in enumerate(freqs)
        
        ωn = 2 * pi * f
        An = A[coluna]
        P = @view MP[:, coluna] 

        # -----------------------------------------------------------
        # MONTAGEM RÁPIDA
        # -----------------------------------------------------------
        # Agora usamos as submatrizes já cortadas.
        # A operação abaixo é rápida pois K_livre, etc, já são CSC compactos.
        Kd = K_livre .+ (im * ωn) .* C_livre .- (ωn^2) .* M_livre

        # Lado direito adjunto
        fill!(Fn, 0.0) 
        for p in nodes_target
            Fn[p] += -conj(P[p])
        end

        # Escalonamento
        P2 = sum(abs2, P[nodes_target]) 
        P2avg = P2 / nt
        scale_factor = An * const_log / (nt * P2avg)
        Fn[nodes_target] .*= scale_factor
        
        # Solver Linear
        λn[livres] .= Kd \ Fn[livres]

        # LOOP DE SENSIBILIDADE (Paralelizado)
        d_temp = zeros(Float64, length(elements_design))
        
        Threads.@threads for i in eachindex(elements_design)


            data = dados_elementos[i]
            pe = @view P[data.nos]
            λe = @view λn[data.nos]
            
            # Produto interno otimizado (dot)
            # Sens = 2 * Real( λ^T * dKde * P )
            
            # Rigidez
            val = dot(conj(λe), data.dKe, pe) 

            # Amortecimento
            val += (im * ωn) * dot(conj(λe), data.dCe, pe)

            # Massa
            val -= (ωn^2) * dot(conj(λe), data.dMe, pe)
            
            d_temp[i] += 2 * real(val)
        end
        
        # Acumula
        for i in 1:length(elements_design)
            ele = elements_design[i]
            d[ele] += d_temp[i]
        end

    end # Fim Loop Frequência

    return d ./ Nf
end

# ===================================================================================
# Calcula a derivada da função objetivo
#
# Média simples do SPL em cada frequência ---ver o NF aqui
#
function Derivada_OLD(ne,nn,γ::Vector{T0},connect::Matrix{T1},coord::Matrix{T0},
                  K::AbstractMatrix{T0},M::AbstractMatrix{T0},C::AbstractMatrix{T0},
                  livres::Vector{T1},freqs::Vector{T0},
                  pressures::Vector, 
                  fρ::Function, fκ::Function,
                  dfρ::Function, dfκ::Function,μ::T0,
                  nodes_target::Vector{T1},MP::Matrix{T2},
                  elements_design::Vector,A::Vector,p0=20E-6) where {T0,T1,T2}


    # Define o vetor de derivadas
    d = zeros(ne)

    # Número de frequências
    Nf = length(freqs)

    # Número de nós para monitorar a pressão
    nt = length(nodes_target)

    # Calcula a constante 10/(ln(10)*nt)
    cte = 10/(log(10)*nt)
 
    # Define λ fora do loop, para reaproveitar
    λn = zeros(T2,nn)

    # Aloca antes do loop
    P = MP[:,1]

    # Aloca Fn 
    Fn = similar(P)

    # Loop pelas frequências
    coluna = 1
    for f in freqs
        
        # Converte a freq para rad/s
        ωn = 2*pi*f

        # Sensibilidade para essa frequência 
        An = A[coluna]

        # Recupera as pressões para essa frequência (coluna de target)
        P .= MP[:,coluna]

        # Monta a matriz de rigidez dinâmica
        Kd = K[livres,livres]  .+im*ωn*C[livres,livres] .- (ωn^2)*M[livres,livres]

        # Monta o vetor adjunto para essa frequência
        Fn .= F_adj(nodes_target,P)

        # Calcula o Pn2
        P2 = sum((abs2.(P[nodes_target])))

        # Média (pelo número de pontos em nodes_target)
        P2avg = P2 / nt

        # Escalona F pela cte e pelo Nf
        Fn .= An*Fn*cte/(P2avg)
        
        # Soluciona o problema adjunto, obtendo λ^n
        λn[livres] .= Kd\Fn[livres]

        # Loop pelos elementos
        for ele in elements_design

            # Tipo de elemento
            etype = connect[ele,1]

            # Localizações 
            nos, X = LSound.Nos_Coordenadas(ele,etype,coord,connect)
  
            # Pressão nos nós do elementos
            pe = P[nos]

            # Vetor adjunto nos nós do elemento
            λe = λn[nos]

            # Variável de projeto do elemento
            γe = γ[ele]

            # Calcula a derivada da rigidez dinâmica do elemento
            dKe, dMe, dCe = Derivada_KMC(etype,γe,fρ,fκ,dfρ,dfκ,μ,X)

            # Derivada da matriz dinâmica do elemento
            dKde = dKe .+im*ωn*dCe .- dMe*ωn^2  

            # Derivada de Fn [dFn/dγm]
            # dFp =  Derivada_forca_pressao(nos,pressures,dKde)
            # -2*real(transpose(λe)*dFp)

            # Calcula a derivada e sobrepõe na posição do elemento
            d[ele] += 2*real(transpose(λe)*dKde*pe) 

        end # Elemento

        # Atualiza a coluna em target
        coluna += 1

    end # Frequência

    #
    # TESTE
    #
    # γ[γ.==0]    .= 0
    


    # Retorna a derivada, lembrando de dividir pelo número de 
    # frequências
    return d./Nf
         
end
