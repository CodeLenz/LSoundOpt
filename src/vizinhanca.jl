#
# Retorna um vetor de vetores, com os vizinhos de cada elemento 
# considerando SOMENTE os elementos de projeto
#
function Vizinhanca(ne,centroides,raio_filtro,elements_design::Vector)

    # Número de variáveis de projeto
    np = length(elements_design)

    # Aloca a lista de saída para os vizinhos
    vizinhos = Vector{Vector{Int64}}(undef,np)

    # Aloca a lista de pesos
    pesos = Vector{Vector{Float64}}(undef,np)

    # Número mínimo e máximo de vizinhos de 
    # um elemento da malha
    n_min_viz = np
    n_max_viz = 0

    # Loop pelos elementos da malha
    contador = 1
    for ele in elements_design 

        # Centroide deste elemento
        cele = centroides[ele,:]

        # Aloca um vetor para armazenarmos os vizinhos deste elemento
        vele = Int64[]

        # Aloca um vetor para armazenarmos os pesos
        pele = Float64[]

        # Loop por todos os elementos da malha
        for viz in elements_design

            # Centroide do vizinho
            cviz = centroides[viz,:]

            # Distância entre os centróides
            dist = norm(cviz.-cele)

            # Se essa distância for menor do que o raio_filtro,
            # armazena como vizinho
            if dist < raio_filtro
               push!(vele, viz)
               push!(pele, 1 - dist/raio_filtro)
            end

        end # viz

        # Número de vizinhos deste elemento 
        nviz = length(vele)

        # Verifica número de vizinhos
        n_min_viz = min(n_min_viz,nviz) 
        n_max_viz = max(n_max_viz,nviz)

        # Armazena vele na linha ele de vizinhos
        vizinhos[contador] = copy(vele)

        # Armazena pele na linha ele de pesos
        pesos[contador] = copy(pele)

        # Acumula o contador
        contador += 1

    end # ele

    # Mostra o número mínimo e máximo de vizinhos
    println("Número mínimo de vizinhos na malha: ", n_min_viz)
    println("Número máximo de vizinhos na malha: ", n_max_viz)

    # Retorna o vetor de vetores com os vizinhos do elemento
    # e também os pesos. Esses vetores de vetores 
    # contém somente os elementos de projeto
    return vizinhos, pesos, n_min_viz, n_max_viz

end

#
# Find all elements that share edges with the current element
#
function NeighborEdges(ne, connect, elements_design)

    # 1. Dicionário para mapear Face -> Elementos
    # Chave: Vetor ordenado de nós da face (para garantir unicidade)
    # Valor: Vetor de IDs dos elementos que compartilham essa face
    face_map = Dict{Vector{Int64}, Vector{Int64}}()
    
    # Vetor de saída (lista de adjacência)
    # Usamos Dict{Int, Vector} para lidar com IDs esparsos se houver
    # Se seus IDs forem contíguos de 1 a ne, pode usar Vector{Vector}
    vizinhos = Dict{Int64, Vector{Int64}}()
    for ele in elements_design
        vizinhos[ele] = Int64[]
    end

    # Construção do Mapa de Faces ---
    for ele in elements_design
        etype = connect[ele, 1]
        nodes = connect[ele, 3:end]
        
        # Obtém as faces baseadas no tipo do elemento
        faces = Get_Element_Faces(etype, nodes)
        
        for f in faces

            # Ordena para que [1,2,3] seja igual a [3,1,2]
            sort!(f) 
            
            if !haskey(face_map, f)
                face_map[f] = Int64[]
            end
            
            push!(face_map[f], ele)
        end
    end

    # Conexão dos Vizinhos ---
    # Se uma face tem 2 elementos associados, eles são vizinhos.
    for (face, elems) in face_map
        if length(elems) == 2
            e1, e2 = elems[1], elems[2]
            
            # Adiciona bidirecionalmente (evita duplicatas se já houver)
            if !(e2 in vizinhos[e1]); push!(vizinhos[e1], e2); end
            if !(e1 in vizinhos[e2]); push!(vizinhos[e2], e1); end
        elseif length(elems) > 2
            println("Aviso: Aresta/Face não-manifold detectada (mais de 2 elementos compartilhados).")
        end
    end

    # Converte de volta para Vector{Vector} se seus IDs forem 1..ne contíguos
    # Caso contrário, retorne o Dict vizinhos.
    vizinhos_vec = Vector{Vector{Int64}}(undef, ne)
    for i in 1:ne
        if haskey(vizinhos, i)
            vizinhos_vec[i] = vizinhos[i]
        else
            vizinhos_vec[i] = Int64[]
        end
    end

    return vizinhos_vec
end

#
# Função auxiliar para extrair faces/arestas baseada no padrão Gmsh
#
function Get_Element_Faces(etype, nodes)
    faces = Vector{Vector{Int64}}()
    
    if etype == 2 # Triângulo (3 nós) -> Faces são Arestas
        push!(faces, [nodes[1], nodes[2]])
        push!(faces, [nodes[2], nodes[3]])
        push!(faces, [nodes[3], nodes[1]])
        
    elseif etype == 3 # Quadrilátero (4 nós) -> Faces são Arestas
        push!(faces, [nodes[1], nodes[2]])
        push!(faces, [nodes[2], nodes[3]])
        push!(faces, [nodes[3], nodes[4]])
        push!(faces, [nodes[4], nodes[1]])
        
    elseif etype == 4 # Tetraedro (4 nós) -> Faces são Triângulos
        push!(faces, [nodes[1], nodes[2], nodes[3]])
        push!(faces, [nodes[1], nodes[2], nodes[4]])
        push!(faces, [nodes[1], nodes[3], nodes[4]])
        push!(faces, [nodes[2], nodes[3], nodes[4]])
        
    elseif etype == 5 # Hexaedro (8 nós) -> Faces são Quads
        push!(faces, [nodes[1], nodes[2], nodes[6], nodes[5]])
        push!(faces, [nodes[2], nodes[3], nodes[7], nodes[6]])
        push!(faces, [nodes[3], nodes[4], nodes[8], nodes[7]])
        push!(faces, [nodes[4], nodes[1], nodes[5], nodes[8]])
        push!(faces, [nodes[1], nodes[2], nodes[3], nodes[4]]) # Base
        push!(faces, [nodes[5], nodes[6], nodes[7], nodes[8]]) # Topo

    elseif etype == 6 # Prisma (6 nós) -> 2 Tri, 3 Quads
        push!(faces, [nodes[1], nodes[2], nodes[3]]) # Base Tri
        push!(faces, [nodes[4], nodes[5], nodes[6]]) # Topo Tri
        push!(faces, [nodes[1], nodes[2], nodes[5], nodes[4]]) # Quad
        push!(faces, [nodes[2], nodes[3], nodes[6], nodes[5]]) # Quad
        push!(faces, [nodes[3], nodes[1], nodes[4], nodes[6]]) # Quad

    elseif etype == 7 # Pirâmide (5 nós) -> 1 Quad (Base), 4 Tri (Lados)
        push!(faces, [nodes[1], nodes[2], nodes[3], nodes[4]]) # Base Quad
        push!(faces, [nodes[1], nodes[2], nodes[5]]) # Tri
        push!(faces, [nodes[2], nodes[3], nodes[5]]) # Tri
        push!(faces, [nodes[3], nodes[4], nodes[5]]) # Tri
        push!(faces, [nodes[4], nodes[1], nodes[5]]) # Tri
    
    else
        error("Get_Element_Faces:: Elemento não suportado: $etype")
    end
    
    return faces
end

function NeighborEdgesOLD(ne,connect,elements_design)

    # Vector of vectors
    vizinhos = Vector{Vector{Int64}}(undef,ne)

    # Local vector 
    local_vector = Int64[]

    # Loop pelos elementos de projeto
    for ele in elements_design

        # Limpa o vetor local 
        empty!(local_vector)

        # Nós deste elemento 
        nos_ele = connect[ele,3:end]

        # Tipo de elemento 
        etype = connect[ele,1]

        # Dependendo do tipo de elemento, temos diferentes requisitos para 
        # comparação do número de nós que são necessários para definir uma 
        # face. Não é tão simples se considerarmos a pirâmide de 5 nós, mas 
        # para os outros podemos fazer direto 
        #
        # Válido para triângulo e para quadrado (etypes 2 e 3)
        n_lado = 2
        if etype==4
            n_lado  = 3
        elseif  etype==5
            n_lado = 4
        elseif etype>5
            error("fix this code...") 
        end

        # Número de nós do elemento -2 (para comparação)
        ncompara = length(nos_ele)-n_lado

        # Loop pelos outros elementos de projeto 
        for viz in elements_design

            # pulamos o próprio elemento 
            if ele!=viz

               # Nós deste elemento 
               nos_viz = connect[viz,3:end]

               # Verifica se temos ao menos 2 nós em comum. Vou fazer isso com o
               # comando setfiff(a,b), que retorna todos os elementos que estão em a 
               # mas não estão em b. Se 
               compara = setdiff(nos_ele,nos_viz)

               # Se tivemos (ao menos) menos 2 nós em compara
               # então eles são em comum com o candidato a vizinho
               # Guarda para vizinho de ele 
               if length(compara)<=ncompara
                  push!(local_vector,viz)
               end

            end


        end #viz

        # Armazena no vetor de vetores com os vizinhos
        vizinhos[ele] = copy(local_vector)
        
    end #ele

    # Retorna 
    return vizinhos

end