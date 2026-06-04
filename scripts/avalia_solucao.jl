using LSoundOpt

function parse_freqs(spec::AbstractString)
    if occursin(":", spec)
        parts = parse.(Float64, split(spec, ":"))
        if length(parts) == 2
            return collect(parts[1]:parts[2])
        elseif length(parts) == 3
            return collect(parts[1]:parts[2]:parts[3])
        end
    elseif occursin(",", spec)
        return parse.(Float64, split(spec, ","))
    else
        return [parse(Float64, spec)]
    end

    error("Faixa de frequencias invalida: $spec")
end

function usage()
    println("Uso:")
    println("  julia --project=. scripts/avalia_solucao.jl arquivo.msh inicio:passo:fim solucao.dat")
    println()
    println("Exemplos:")
    println("  julia --project=. scripts/avalia_solucao.jl msh/LShape/lshape.msh 1300:5:1500 msh/LShape/lshape_γ_opt.dat")
    println("  julia --project=. scripts/avalia_solucao.jl msh/LShape/lshape.msh 1300,1400,1500 msh/LShape/lshape_γ_opt.dat")
end

if length(ARGS) != 3
    usage()
    exit(1)
end

meshfile = ARGS[1]
freqs = parse_freqs(ARGS[2])
datfile = ARGS[3]

resultado = Avalia_Solucao(meshfile, freqs, datfile)

println("Arquivo .msh: ", meshfile)
println("Arquivo .dat: ", datfile)
println("SLP medio: ", resultado.SLP)
println("Volume: ", resultado.volume)
println("Volume total de projeto: ", resultado.volume_total_projeto)
println("Fracao de volume: ", resultado.fracao_volume)
println("Perimetro: ", resultado.perimetro)
println()
println("freq,SLP")
for (freq, spl) in zip(resultado.freqs, resultado.SPL_por_frequencia)
    println(freq, ",", spl)
end
