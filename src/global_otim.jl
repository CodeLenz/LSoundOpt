# ===================================================================================
# Monta_KM, versão para otimização
#
# γ   ->  vetor de variáveis de projeto [0,1] 
# fρ  ->  função que parametriza ρ
# fκ  ->  função que parametriza κ
# ===================================================================================

function Monta_KMC_param(ne,coord,connect,γ::Vector,fρ::Function,fκ::Function,μ::Float64)  
    
    # Aloca vetores para a montagem eficiente 
    # das matrizes esparsas 
    I = Int64[]
    J = Int64[]
    VK = Float64[]
    VM = Float64[]
    VC = Float64[]
 
    # Loop pelos elementos
    for ele=1:ne

        # Variável de projeto do elemento
        γe = γ[ele]

        # Calcula as inversas 
        iρ = fρ(γe)
        iκ = fκ(γe)

        # Cte for damping 
        cte = (4/3)*μ

        # Tipo de elemento
        et = connect[ele,1]

        # Descobre nos, X e Y para este elemento
        nos, X = LSound.Nos_Coordenadas(ele,et,coord,connect) 

        # Monta as matrizes dos elementos
        if et==3
           Ke, Me = LSound.KMe_bi4(iρ,iκ,X)
        elseif et==2
           Ke, Me = LSound.KMe_tri3(iρ,iκ,X)
        elseif et==4
            Ke, Me = LSound.KMe_tet4(iρ,iκ,X)    
        elseif et==5
           Ke, Me = LSound.KMe_hex8(iρ,iκ,X) 
        elseif et==7
            Ke, Me = LSound.KMe_pyr5(iρ,iκ,X)    
         else
            error("Elemento não definido")
        end
 
        # Sobreposição das locais nas globais
        for i in LinearIndices(nos)
            ni = nos[i]
            for j in LinearIndices(nos)
                push!(I,ni)
                push!(J,nos[j])
                push!(VK,Ke[i,j])
                push!(VM,Me[i,j])
                push!(VC,cte*Ke[i,j])
            end
        end

    end

    # Retorna as matrizes globais
    return sparse(I,J,VK), sparse(I,J,VM), sparse(I,J,VC)

end