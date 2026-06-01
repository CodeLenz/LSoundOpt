#
# Validação das derivadas por DFC
#
function Verifica_derivada(γ,nn,ne,coord,connect,fρ,fκ,μ,freqs,livres,velocities,pressures,nodes_target,elements_design,A)
    
    # Vamos validar a derivada usando diferenças finitas
    function f_(γ,nn,ne,coord,connect,fρ,fκ,μ, freqs,livres,velocities,pressures,nodes_target,A)

        # Aloca a matriz MP
        MP = zeros(ComplexF64, nn, length(freqs))

        Sweep!(nn,ne,coord,connect,γ,fρ,fκ,μ,freqs,livres,velocities,pressures,MP) 

        # Calcula a função objetivo SPL_w
        objetivo = Objetivo(MP,nodes_target,A)

        return objetivo

    end
    f(γ) = f_(γ,nn,ne,coord,connect,fρ,fκ,μ, freqs,livres,velocities,pressures,nodes_target,A)

    # Calcula a derivada por DFC
    d_numerica = df(γ,f,elements_design,1E-6)

end