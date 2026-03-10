import gmsh
import meshio
import numpy as np
import sys

arquivo_pos = sys.argv[1]

gmsh.initialize()
gmsh.open(arquivo_pos)

# ---------------------------------------------------
# 1. Ler nós
# ---------------------------------------------------
nodeTags, nodeCoords, _ = gmsh.model.mesh.getNodes()

points = np.array(nodeCoords).reshape(-1, 3)

# mapear tag -> índice sequencial (VTK exige isso)
tag_to_index = {tag: i for i, tag in enumerate(nodeTags)}

# ---------------------------------------------------
# 2. Ler elementos
# ---------------------------------------------------
types, elemTags, elemNodeTags = gmsh.model.mesh.getElements()

cells = []
cell_blocks_elem_tags = []

# mapeamento gmsh -> meshio
cell_map = {
    (2,3): "triangle",
    (2,4): "quad",
    (3,4): "tetra",
    (3,8): "hexahedron",
    (3,6): "wedge",
    (3,5): "pyramid"
}

for i, etype in enumerate(types):

    name, dim, order, numNodes, _, _ = gmsh.model.mesh.getElementProperties(etype)

    key = (dim, numNodes)

    if key not in cell_map:
        continue

    vtk_type = cell_map[key]

    conn = np.array(elemNodeTags[i]).reshape(-1, numNodes)

    # converter tags para índices
    conn = np.vectorize(tag_to_index.get)(conn)

    cells.append((vtk_type, conn))
    cell_blocks_elem_tags.append(np.array(elemTags[i]))

# total de elementos
all_elem_tags = np.concatenate(cell_blocks_elem_tags)

elem_tag_to_index = {tag: i for i, tag in enumerate(all_elem_tags)}

# ---------------------------------------------------
# 3. Extrair dados das views
# ---------------------------------------------------
point_data = {}
cell_data_raw = {}

for tag in gmsh.view.getTags():

    dataType, entityTags, data, time, numComp = gmsh.view.getModelData(tag, 0)

    values = np.array(data).reshape(-1, numComp)

    name = f"view_{tag}"

    # ----------------------------
    # NodeData
    # ----------------------------
    if dataType == "NodeData":

        arr = np.zeros((len(nodeTags), numComp))

        for i, ntag in enumerate(entityTags):
            idx = tag_to_index.get(ntag)
            if idx is not None:
                arr[idx] = values[i]

        point_data[name] = arr

    # ----------------------------
    # ElementData
    # ----------------------------
    elif dataType == "ElementData":

        arr = np.zeros((len(all_elem_tags), numComp))

        for i, etag in enumerate(entityTags):
            idx = elem_tag_to_index.get(etag)
            if idx is not None:
                arr[idx] = values[i]

        cell_data_raw[name] = arr

gmsh.finalize()

# ---------------------------------------------------
# 4. Distribuir ElementData por bloco
# (meshio exige isso)
# ---------------------------------------------------
cell_data = {}

for name, arr in cell_data_raw.items():

    split = []
    offset = 0

    for block in cell_blocks_elem_tags:

        n = len(block)

        split.append(arr[offset:offset+n])

        offset += n

    cell_data[name] = split

# ---------------------------------------------------
# 5. Criar malha meshio
# ---------------------------------------------------
mesh = meshio.Mesh(
    points=points,
    cells=cells,
    point_data=point_data,
    cell_data=cell_data
)

# ---------------------------------------------------
# 6. Salvar
# ---------------------------------------------------
saida = arquivo_pos.replace(".pos", ".vtu")

mesh.write(saida)

print("=================================")
print("Conversão concluída")
print("Arquivo:", saida)
print("Nós:", len(points))
print("Blocos de células:", len(cells))
print("Campos nodais:", list(point_data.keys()))
print("Campos elementares:", list(cell_data.keys()))
print("=================================")
