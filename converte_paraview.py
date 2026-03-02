import meshio
import numpy as np

# 1. Lê o arquivo MSH gerado pelo Gmsh
mesh = meshio.read("resultado_completo.msh")

# 2. Força todos os dados nodais (Point Data) a serem Float64
for key, data in mesh.point_data.items():
    mesh.point_data[key] = data.astype(np.float64)

# 3. Força todos os dados elementares (Cell Data) a serem Float64
# O meshio armazena cell_data como uma lista de arrays (um por bloco de elementos)
for key, data_list in mesh.cell_data.items():
    mesh.cell_data[key] = [d.astype(np.float64) for d in data_list]

# 4. Salva no formato VTU (Ideal para o ParaView)
meshio.write("resultado_final.vtu", mesh)

print("Conversão para VTU (Float64) concluída!")