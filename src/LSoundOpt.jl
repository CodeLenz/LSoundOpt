module LSoundOpt

   # Carrega o LSound
   using LSound

   # Carrega os pacotes para leitura de dados
   using YAML
   using Lgmsh
   using Gmsh
   using DelimitedFiles

   # Pacotes de otimização
   using JuMP
   using Alpine
   using Gurobi
   using Ipopt
   using HiGHS
   using Cbc

   # Inclui os arquivos do pacote

   # Lê os dados do .yaml
   include("le_yaml.jl")

   # Parametrização do material 
   include("param_duhring.jl")

   # Impõe valores fixos 
   include("fixos.jl") 

   # Monta as matrizes globais para otimização 
   include("global_otim.jl")

   # Monta vetor com áreas/volumes dos elementos da malha
   include("volumes.jl")

   # Vizinhança e filtro 
   include("vizinhanca.jl")
   include("filtro_espacial.jl")

   # Diferenças finitas (para verificação)
   include("df.jl")

   # Script para verificação dos gradientes
   include("verifica_derivada.jl")

   # Calcula a função objetivo SPL
   include("objetivo.jl")

   # Calcula a derivada da função objetivo
   include("sensibilidade.jl")

   # Calcula o perímetro e seu gradiente
   include("perimiter.jl")

   # Soluciona o problema de otimização ISLP
   include("LP.jl")

   # Realiza um sweep 
   include("sweep.jl")

   # Processa a FRF da solução 
   include("processa_FRF.jl")

   # Arquivo principal 
   include("main_ISLP.jl")


end
