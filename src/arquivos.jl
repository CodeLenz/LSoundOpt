
#
# Nome dos arquivos de saída, baseados no nome da entrada
#
function Setup_Arquivos(arquivo)

    # Processa .geo ou .msh
    if occursin(".geo", arquivo)

       # Chama o gmsh para gerar o .msh
       gmsh.initialize(); 
       gmsh.open(arquivo); 
       gmsh.model.mesh.generate(2)
       mshfile = replace(arquivo, ".geo"=>".msh"); 
       gmsh.write(mshfile)
    else 
       mshfile = arquivo
    end
    
    # Gera os arquivos de saída com os seus nomes
    nomebase = basename(mshfile)
    saidas = (
        replace(nomebase, ".msh"=>".pos"),
        replace(nomebase, ".msh"=>"_freq.pos"),
        replace(nomebase, ".msh"=>".data"),
        replace(nomebase, ".msh"=>"_γ_ini.dat"),
        replace(nomebase, ".msh"=>"_γ_opt.dat")
    )
    return mshfile, saidas
end


#
# Salva histórico e dados do filtro
#
function Salva_Historico(arquivo, hist, raio)

    fd = open(arquivo,"w")
    println(fd, "Historico V")
    println(fd, hist.V)
    println(fd, "Historico SPL")
    println(fd, hist.SLP)
    println(fd, "Historico P")
    println(fd, hist.P)
    println(fd, "Raio do filtro ", raio)
    close(fd)
    
end