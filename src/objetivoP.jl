# ===================================================================================
# Função Objetivo P-Norm ABSOLUTA (Sem normalização por espectro de referência)
#
function ObjetivoP(MP::Matrix, nodes_target::Vector, A::Vector, p0=20E-6; p_norm=8.0)

    Nf = size(MP, 2)
    Rn_values = zeros(Float64, Nf)
    
    coluna = 1
    for P in eachcol(MP)
        # Pressão quadrática média absoluta
        P2_avg = sum(abs2.(P[nodes_target])) / length(nodes_target)
        
        # Razão adimensional absoluta (P_rms / p0)^2
        Rn = P2_avg / (p0^2)
        
        # Armazena
        Rn_values[coluna] = A[coluna] * Rn
        coluna += 1
    end

    # Média Generalizada das Potências Adimensionais
    # Se p_norm for alto (ex: 8), isso foca nos picos de pressão
    soma_potencias = sum(Rn_values .^ p_norm)
    Mp = soma_potencias / Nf

    # Retorna em dB
    return (10.0 / p_norm) * log10(Mp)
end