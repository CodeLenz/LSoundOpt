 # Adiciona os pacotes
  using  Pkg
  Pkg.add(url="https://github.com/CodeLenz/Lgmsh")
  Pkg.add(url="https://github.com/CodeLenz/LSound")
  Pkg.add(url="https://github.com/CodeLenz/LSoundOpt")

  # Roda o programa 
  using LSoundOpt

  # Roda o programa 
  x = LSoundOpt.Otim_ISLP("muffler3D.msh",collect(380:5.0:420),[])