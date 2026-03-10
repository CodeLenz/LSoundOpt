import gmsh
import meshio
import numpy as np
import sys

arquivo_pos = sys.argv[1]

gmsh.initialize()
gmsh.open(arquivo_pos)

# --- Nós ---
nodeTags, nodeCoords, _ = gmsh.model.mesh.getNodes()
points = np.array(nodeCoords).reshape(-1, 3)
tag_to_index = {tag: i for i, tag in enumerate(nodeTags)}

# --- Elementos 2D ---
types, elemTags, elemNodeTags = gmsh.model.mesh.getElements()

cells = []
elem_tag_list = []

for i, etype in enumerate(types):
    name, dim, order, numNodes, _, _ = gmsh.model.mesh.getElementProperties(etype)
    if dim == 2:
        conn = np.array(elemNodeTags[i]).reshape(-1, numNodes)
        conn = np.vectorize(tag_to_index.get)(conn)
        if numNodes == 3:
            cells.append(("triangle", conn))
        elif numNodes == 4:
            cells.append(("quad", conn))
        elem_tag_list.extend(elemTags[i])

# --- Dados ---
point_data = {}
cell_data = {}

for tag in gmsh.view.getTags():
    dataType, entityTags, data, time, numComp = gmsh.view.getModelData(tag, 0)
    values = np.array(data).reshape(-1, numComp)

    if dataType == "NodeData":
        arr = np.zeros((len(nodeTags), numComp))
        for i, ntag in enumerate(entityTags):
            idx = tag_to_index[ntag]
            arr[idx] = values[i]
        point_data[f"view_{tag}"] = arr

    elif dataType == "ElementData":
        cell_data[f"view_{tag}"] = [values]

gmsh.finalize()

mesh = meshio.Mesh(
    points=points,
    cells=cells,
    point_data=point_data,
    cell_data=cell_data
)

saida = arquivo_pos.replace(".pos", ".vtu")
mesh.write(saida)

print(f"Arquivo {saida} gerado.")