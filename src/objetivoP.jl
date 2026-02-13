# ===================================================================================
# Função Objetivo P-Norm Normalizada
#
# Minimiza a média generalizada da TRANSMISSÃO RELATIVA (P_curr / P_ref).
# Isso aproxima a maximização do Insertion Loss mínimo na banda.
#
function ObjetivoP(MP::Matrix, nodes_target::Vector, A::Vector, P_ref::Vector; p_norm=8.0)

    Nf = size(MP, 2)
    Rn_values = zeros(Float64, Nf)
    
    coluna = 1
    for P in eachcol(MP)
        # Pressão quadrática média atual
        P2_avg = sum(abs2.(P[nodes_target])) / length(nodes_target)
        
        # Razão de Potência Normalizada (R_hat do texto)
        # R_n = P_atual / P_referencia
        # Nota: p0 cancela, então usamos direto as médias quadráticas
        Rn = P2_avg / P_ref[coluna]
        
        # Aplica peso de frequência se houver
        Rn_values[coluna] = A[coluna] * Rn
        coluna += 1
    end

    # Agregação P-Norm
    soma_potencias = sum(Rn_values .^ p_norm)
    Mp = soma_potencias / Nf

    # Retorna em dB: 10 * log10( (Mp)^(1/p) )
    Φ = (10.0 / p_norm) * log10(Mp)

    return Φ
end