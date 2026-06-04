#
# Avalia uma topologia carregada de um arquivo .dat
#
# O arquivo .dat pode conter um vetor por elemento da malha (ne valores)
# ou apenas as variáveis dos elementos de projeto.
#
function Avalia_Solucao(meshfile::String, freqs::Vector, datfile::String; vA=Float64[])

    # Evita chamar um .geo
    occursin(".geo", meshfile) && error("Avalia_Solucao:: chamar com arquivo .msh")

    # Verificamos se existem frequências sendo informadas
    isempty(freqs) && error("Avalia_Solucao:: freqs deve ser um vetor não vazio")

    # Garante que as frequências e pesos sejam vetores de Float64
    freqs = Float64.(freqs)
    weights = isempty(vA) ? ones(length(freqs)) : Float64.(vA)
    length(weights) == length(freqs) ||
        error("Avalia_Solucao:: vA deve ter o mesmo tamanho de freqs")

    # Verifica se os arquivos de entrada existem
    isfile(meshfile) || error("Avalia_Solucao:: arquivo de entrada $meshfile não existe")
    isfile(datfile) || error("Avalia_Solucao:: arquivo de entrada $datfile não existe")

    # Arquivo .yaml associado à malha
    yamlfile = replace(meshfile, ".msh" => ".yaml")
    isfile(yamlfile) || error("Avalia_Solucao:: arquivo de entrada $yamlfile não existe")

    # Lê os dados da malha
    nn, coord, ne, connect, materials, nodes_open, velocities, nodes_pressure,
    pressures, damping, nodes_probe, nodes_target, elements_fixed, values_fixed,
    centroides = LSound.Parsemsh_Daniele(meshfile)

    # Precisamos de pelo menos um nó alvo para calcular o SPL
    isempty(nodes_target) &&
        error("Avalia_Solucao:: nodes_target deve ter ao menos um nó informado")

    # Precisamos de um material
    isempty(materials) &&
        error("Avalia_Solucao:: at least one material is necessary")

    # Vamos colocar nodes_target em ordem crescente
    sort!(nodes_target)

    # Lista com os elementos que são de projeto
    elements_design = setdiff(1:ne, sort!(elements_fixed))

    # Lê os dados do arquivo yaml
    _, _, _, _, mu, _ = Le_YAML(yamlfile)

    # Lê as variáveis de projeto carregadas do arquivo .dat
    gamma = _Carrega_Gamma_Solucao(datfile, ne, elements_design, elements_fixed, values_fixed)

    # Posições que não precisam ser calculadas no sistema de equações
    nodes_mask = sort(vcat(nodes_open, nodes_pressure))
    livres = setdiff(collect(1:nn), nodes_mask)

    # Aloca MP fora do sweep
    nf = length(freqs)
    MP = zeros(ComplexF64, nn, nf)

    # Roda o sweep na topologia carregada
    Sweep!(nn, ne, coord, connect, gamma, fρ, fκ, mu, freqs, livres, velocities, pressures, MP)

    # Calcula o SPL em cada frequência e o objetivo médio
    spl_por_frequencia = [SPLn(MP[nodes_target, i], 20E-6) for i in 1:nf]
    slp = Objetivo(MP, nodes_target, weights)

    # Monta o vetor com os volumes de cada elemento da malha
    V = Volumes(ne, connect, coord)

    # Volume atual do projeto
    volume = sum(gamma[elements_design] .* V[elements_design])

    # Perímetro atual do projeto
    neighedge = NeighborEdges(ne, connect, elements_design)
    perimetro = Perimiter(gamma, neighedge, elements_design)

    # Retorna os valores avaliados
    return (
        freqs = freqs,
        SLP = slp,
        SPL_por_frequencia = spl_por_frequencia,
        volume = volume,
        perimetro = perimetro,
        gamma = gamma,
    )
end

#
# Lê o vetor de variáveis de projeto de um arquivo .dat
#
function _Carrega_Gamma_Solucao(datfile, ne, elements_design, elements_fixed, values_fixed)

    # Leitura do arquivo .dat
    dados = vec(readdlm(datfile))
    gamma_lido = Float64.(dados)

    # Caso o arquivo contenha uma variável por elemento da malha
    if length(gamma_lido) == ne
        return gamma_lido
    end

    # Caso o arquivo contenha apenas os elementos de projeto
    if length(gamma_lido) == length(elements_design)
        gamma = zeros(ne)
        Fix_γ!(gamma, elements_fixed, values_fixed)
        gamma[elements_design] .= gamma_lido
        return gamma
    end

    error(
        "Avalia_Solucao:: $datfile contem $(length(gamma_lido)) valores; " *
        "esperado $ne ou $(length(elements_design))",
    )
end
