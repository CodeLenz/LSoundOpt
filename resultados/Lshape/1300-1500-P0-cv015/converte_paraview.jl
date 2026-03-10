using Gmsh
using WriteVTK
using WriteVTK: MeshCell

const gmsh = Gmsh.gmsh

arquivo_pos = ARGS[1]

gmsh.initialize()
gmsh.open(arquivo_pos)

# ---------------------------------------------------
# 1. Ler nós
# ---------------------------------------------------
nodeTags, nodeCoords, _ = gmsh.model.mesh.getNodes()

n_nodes = length(nodeTags)

points = reshape(nodeCoords, 3, :)

tag_to_index = Dict(nodeTags[i] => i for i in eachindex(nodeTags))

# ---------------------------------------------------
# 2. Ler elementos
# ---------------------------------------------------
types, elemTags, elemNodeTags = gmsh.model.mesh.getElements()

cells = Vector{Tuple{Symbol,Matrix{Int}}}()
cell_blocks_elem_tags = Vector{Vector{Int}}()

vtk_types = Dict(
    (2,3) => :triangle,
    (2,4) => :quad,
    (3,4) => :tetra,
    (3,8) => :hexahedron,
    (3,6) => :wedge,
    (3,5) => :pyramid
)

for i in eachindex(types)

    name, dim, order, numNodes, _, _ =
        gmsh.model.mesh.getElementProperties(types[i])

    numNodes = Int(numNodes)

    key = (dim, numNodes)

    if !haskey(vtk_types, key)
        continue
    end

    vtk_type = vtk_types[key]

    conn_vec = elemNodeTags[i]

    n_elem = length(conn_vec) ÷ numNodes

    conn = reshape(conn_vec, numNodes, n_elem)'

    for j in axes(conn,1), k in axes(conn,2)
        conn[j,k] = tag_to_index[conn[j,k]]
    end

    push!(cells, (vtk_type, conn))
    push!(cell_blocks_elem_tags, Int.(elemTags[i]))

end

all_elem_tags = vcat(cell_blocks_elem_tags...)

elem_tag_to_index = Dict(all_elem_tags[i] => i for i in eachindex(all_elem_tags))

# ---------------------------------------------------
# 3. Extrair dados das views
# ---------------------------------------------------
point_data = Dict{String,Any}()
cell_data_raw = Dict{String,Matrix{Float64}}()

for tag in gmsh.view.getTags()

    dataType, entityTags, data, time, numComp =
        gmsh.view.getModelData(tag, 0)

    if data isa Vector{Float64}
        values = reshape(data, :, 1)
    else
        values = reduce(hcat, data)'
    end

    name = "view_$tag"

    if dataType == "NodeData"

        arr = zeros(Float64, n_nodes, size(values,2))

        for i in eachindex(entityTags)

            ntag = entityTags[i]

            if haskey(tag_to_index, ntag)
                idx = tag_to_index[ntag]
                arr[idx,:] = values[i,:]
            end

        end

        point_data[name] = size(values,2) == 1 ? vec(arr) : arr

    elseif dataType == "ElementData"

        arr = zeros(Float64, length(all_elem_tags), size(values,2))

        for i in eachindex(entityTags)

            etag = entityTags[i]

            if haskey(elem_tag_to_index, etag)
                idx = elem_tag_to_index[etag]
                arr[idx,:] = values[i,:]
            end

        end

        cell_data_raw[name] = arr

    end

end

gmsh.finalize()

# ---------------------------------------------------
# 4. Distribuir CellData por bloco
# ---------------------------------------------------
cell_data = Dict{String,Vector{Any}}()

for (name, arr) in cell_data_raw

    split = Vector{Any}()

    offset = 1

    for block in cell_blocks_elem_tags

        n = length(block)

        sub = arr[offset:offset+n-1, :]

        push!(split, size(sub,2)==1 ? vec(sub) : sub)

        offset += n

    end

    cell_data[name] = split

end

# ---------------------------------------------------
# 5. Criar células VTK
# ---------------------------------------------------
vtk_cells = MeshCell[]

vtk_map = Dict(
    :triangle   => WriteVTK.VTKCellTypes.VTK_TRIANGLE,
    :quad       => WriteVTK.VTKCellTypes.VTK_QUAD,
    :tetra      => WriteVTK.VTKCellTypes.VTK_TETRA,
    :hexahedron => WriteVTK.VTKCellTypes.VTK_HEXAHEDRON,
    :wedge      => WriteVTK.VTKCellTypes.VTK_WEDGE,
    :pyramid    => WriteVTK.VTKCellTypes.VTK_PYRAMID
)

for (vtk_type, conn) in cells

    celltype = vtk_map[vtk_type]

    for i in 1:size(conn,1)
        push!(vtk_cells, MeshCell(celltype, conn[i,:]))
    end

end

# ---------------------------------------------------
# 6. Escrever VTU
# ---------------------------------------------------
saida = replace(arquivo_pos, ".pos" => ".vtu")

vtk = vtk_grid(saida, points, vtk_cells)

# NodeData
for (name, arr) in point_data
    vtk[name] = arr
end

# CellData
for (name, blocks) in cell_data
    vtk[name] = vcat(blocks...)
end

vtk_save(vtk)


println("=================================")
println("Conversão concluída")
println("Arquivo: ", saida)
println("Nós: ", n_nodes)
println("Blocos de células: ", length(cells))
println("Campos nodais: ", collect(keys(point_data)))
println("Campos elementares: ", collect(keys(cell_data)))
println("=================================")
