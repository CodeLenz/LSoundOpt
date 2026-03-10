using Gmsh
using WriteVTK
using WriteVTK: MeshCell

const gmsh = Gmsh.gmsh

arquivo = ARGS[1]

gmsh.initialize()
gmsh.open(arquivo)

# --------------------------------------------------
# nós
# --------------------------------------------------
nodeTags, nodeCoords, _ = gmsh.model.mesh.getNodes()

points = reshape(nodeCoords, 3, :)

tag_to_idx = Dict(tag => i for (i,tag) in enumerate(nodeTags))

# --------------------------------------------------
# elementos
# --------------------------------------------------
types, elemTags, elemNodes = gmsh.model.mesh.getElements()

vtk_cells = MeshCell[]

vtk_map = Dict(
    (2,3)=>WriteVTK.VTKCellTypes.VTK_TRIANGLE,
    (2,4)=>WriteVTK.VTKCellTypes.VTK_QUAD,
    (3,4)=>WriteVTK.VTKCellTypes.VTK_TETRA,
    (3,8)=>WriteVTK.VTKCellTypes.VTK_HEXAHEDRON,
    (3,6)=>WriteVTK.VTKCellTypes.VTK_WEDGE,
    (3,5)=>WriteVTK.VTKCellTypes.VTK_PYRAMID
)

for i in eachindex(types)

    _, dim, _, nnode, _, _ =
        gmsh.model.mesh.getElementProperties(types[i])

    key = (dim, Int(nnode))

    haskey(vtk_map, key) || continue

    celltype = vtk_map[key]

    conn = elemNodes[i]

    ne = length(conn) ÷ nnode

    conn = reshape(conn, nnode, ne)'

    for j in axes(conn,1), k in axes(conn,2)
        conn[j,k] = tag_to_idx[conn[j,k]]
    end

    for j in 1:size(conn,1)
        push!(vtk_cells, MeshCell(celltype, conn[j,:]))
    end

end

# --------------------------------------------------
# dados
# --------------------------------------------------
point_data = Dict()
cell_data = Dict()

for tag in gmsh.view.getTags()

    type, entityTags, data, _, _ =
        gmsh.view.getModelData(tag,0)

    values = data isa Vector{Float64} ? reshape(data,:,1) : reduce(hcat,data)'

    name = "view_$tag"

    if type == "NodeData"

        arr = zeros(length(nodeTags), size(values,2))

        for i in eachindex(entityTags)
            arr[tag_to_idx[entityTags[i]],:] = values[i,:]
        end

        point_data[name] = size(arr,2)==1 ? vec(arr) : arr

    elseif type == "ElementData"

        cell_data[name] = size(values,2)==1 ? vec(values) : values

    end

end

gmsh.finalize()

# --------------------------------------------------
# escrever VTU
# --------------------------------------------------
out = replace(arquivo,".pos"=>".vtu")

vtk = vtk_grid(out, points, vtk_cells)

for (k,v) in point_data
    vtk[k] = v
end

for (k,v) in cell_data
    vtk[k] = v
end

vtk_save(vtk)

println("Arquivo gerado: ", out)
