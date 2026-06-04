
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
    #nomebase = basename(mshfile)
    nomebase = mshfile
    saidas = (
        replace(nomebase, ".msh"=>".pos"),
        replace(nomebase, ".msh"=>"_freq.pos"),
        replace(nomebase, ".msh"=>".data"),
        replace(nomebase, ".msh"=>"_γ_ini.dat"),
        replace(nomebase, ".msh"=>"_γ_opt.dat"),
        replace(nomebase, ".msh"=>"_runinfo.txt")
    )
    return mshfile, saidas
end


#
# Salva histórico e dados do filtro
#
function Salva_Historico(arquivo, hist, raio; final=nothing)

    fd = open(arquivo,"w")
    println(fd, "Historico V")
    println(fd, hist.V)
    println(fd, "Historico SPL")
    println(fd, hist.SLP)
    println(fd, "Historico P")
    println(fd, hist.P)
    println(fd, "Raio do filtro ", raio)

    if hasproperty(hist, :cv)
        println(fd, "Historico cv")
        println(fd, hist.cv)
    end

    if hasproperty(hist, :moves)
        println(fd, "Historico trocas")
        println(fd, hist.moves)
    end

    if hasproperty(hist, :stagnation)
        println(fd, "Historico estagnacao")
        println(fd, hist.stagnation)
    end

    if final !== nothing
        println(fd, "Resultado final avaliado")
        println(fd, "Volume ", final.volume)
        println(fd, "SPL ", final.objetivo)
        println(fd, "Perimetro ", final.perimetro)
        println(fd, "cv ", final.cv)
    end

    close(fd)
    
end

#
# Executa um comando externo opcional e retorna "" em caso de falha.
#
function _Read_Cmd_Or_Empty(cmd::Cmd)

    try
        return strip(read(pipeline(cmd, stderr=devnull), String))
    catch
        return ""
    end

end

#
# Salva metadados suficientes para reproduzir uma rodada de otimização.
#
function Salva_RunInfo(arquivo, dados)

    fd = open(arquivo, "w")

    println(fd, "LSoundOpt run info")
    println(fd, "timestamp: ", Libc.strftime("%Y-%m-%dT%H:%M:%S", time()))
    println(fd)

    println(fd, "[arquivos]")
    println(fd, "entrada: ", dados.arquivo_entrada)
    println(fd, "malha: ", dados.mshfile)
    println(fd, "yaml: ", dados.arquivo_yaml)
    println(fd, "topologias_gmsh: ", dados.arquivo_pos)
    println(fd, "frequencias_gmsh: ", dados.arquivo_pos_freq)
    println(fd, "historico: ", dados.arquivo_data)
    println(fd, "gamma_inicial: ", dados.arquivo_gamma_ini)
    println(fd, "gamma_final: ", dados.arquivo_gamma_fin)
    println(fd)

    println(fd, "[problema]")
    println(fd, "frequencias: ", dados.freqs)
    println(fd, "pesos_vA: ", dados.vA)
    println(fd, "restricao_volume: ", dados.restricao_volume)
    println(fd, "nn: ", dados.nn)
    println(fd, "ne: ", dados.ne)
    println(fd, "nvp: ", dados.nvp)
    println(fd, "n_elementos_fixos: ", dados.n_fixed)
    println(fd, "n_nodes_open: ", dados.n_nodes_open)
    println(fd, "n_nodes_pressure: ", dados.n_nodes_pressure)
    println(fd, "n_nodes_target: ", dados.n_nodes_target)
    println(fd)

    println(fd, "[parametros_yaml]")
    println(fd, "raio_filtro: ", dados.raio_filtro)
    println(fd, "niter: ", dados.niter)
    println(fd, "volfrac: ", dados.vf)
    println(fd, "perimetro_limite: ", dados.Past)
    println(fd, "mu: ", dados.mu)
    println(fd, "fatorcv: ", dados.fatorcv)
    println(fd, "volume_limite: ", dados.Vast)
    println(fd)

    println(fd, "[parametros_algoritmo]")
    for name in fieldnames(typeof(dados.alg))
        println(fd, string(name), ": ", getfield(dados.alg, name))
    end
    println(fd)

    println(fd, "[resultado]")
    println(fd, "motivo_parada: ", dados.stop_reason)
    println(fd, "iteracoes_executadas: ", dados.iterations_done)
    println(fd, "objetivo_final: ", dados.final.objetivo)
    println(fd, "volume_final: ", dados.final.volume)
    println(fd, "perimetro_final: ", dados.final.perimetro)
    println(fd, "cv_final: ", dados.final.cv)
    println(fd)

    println(fd, "[git]")
    repo = _Read_Cmd_Or_Empty(`git -C $(dirname(abspath(dados.mshfile))) rev-parse --show-toplevel`)
    if isempty(repo)
        println(fd, "repo: ")
        println(fd, "commit: ")
        println(fd, "branch: ")
        println(fd, "status_short: ")
    else
        println(fd, "repo: ", repo)
        println(fd, "commit: ", _Read_Cmd_Or_Empty(`git -C $repo rev-parse HEAD`))
        println(fd, "branch: ", _Read_Cmd_Or_Empty(`git -C $repo rev-parse --abbrev-ref HEAD`))
        println(fd, "status_short:")
        status = _Read_Cmd_Or_Empty(`git -C $repo status --short`)
        if !isempty(status)
            println(fd, status)
        end
    end
    println(fd)

    println(fd, "[yaml_original]")
    if isfile(dados.arquivo_yaml)
        print(fd, read(dados.arquivo_yaml, String))
    end

    close(fd)

end
