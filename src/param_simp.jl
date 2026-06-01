#=
Parametrização SIMP (Solid Isotropic Material with Penalization)
para o problema acústico reativo.

IMPORTANTE — comparação justa com a formulação binária:
Mantemos a MESMA escolha física da formulação Dühring (param_duhring.jl):
interpolam-se os INVERSOS de ρ e κ entre ar (γ=0) e sólido (γ=1).
A única diferença é a penalização por lei de potência γ^p, que
substitui a interpolação linear (p=1 recupera exatamente o Dühring).

    f□(γ) = 1/□_ar + γ^p * (1/□_s - 1/□_ar)
   df□(γ) = p * γ^(p-1) * (1/□_s - 1/□_ar)

Note que, ao contrário do Dühring, as derivadas NÃO são constantes:
passam a depender de γ. Como as funções são chamadas elemento a
elemento dentro de Derivada_KMC / Monta_KMC_param, isso não exige
nenhuma alteração estrutural no restante do código.

O expoente p é capturado por closure (ver Cria_Param_SIMP), de modo
que a assinatura fρ(γe), dfρ(γe), etc. permanece idêntica à de
param_duhring.jl e os arquivos Sweep!/Derivada/Monta_KMC_param não
precisam saber que estão usando SIMP.
=#

# -----------------------------------------------------------------------------
# Fábrica de funções de parametrização SIMP.
#
# Recebe o expoente de penalização p e as propriedades de base, e devolve
# um NamedTuple com as quatro funções (fρ, fκ, dfρ, dfκ) já com p fixado,
# prontas para serem passadas a Sweep!, Derivada, etc.
#
# p        -> expoente de penalização (p >= 1; p=1 == Dühring linear)
# ρ_ar,ρ_s -> densidades de ar e sólido
# κ_ar,κ_s -> módulos de compressibilidade de ar e sólido
#             (mesmos valores default de param_duhring.jl)
# -----------------------------------------------------------------------------
function Cria_Param_SIMP(p::Real; ρ_ar=1.21, ρ_s=2700.0,
                                   κ_ar=142355.29, κ_s=675e8)

    # Consistência do expoente
    p >= 1 || error("Cria_Param_SIMP:: expoente p deve ser >= 1")

    # Diferenças (inversos) pré-calculadas — constantes do problema
    Δiρ = (1/ρ_s - 1/ρ_ar)
    Δiκ = (1/κ_s - 1/κ_ar)
    iρ_ar = 1/ρ_ar
    iκ_ar = 1/κ_ar

    # --- Inverso de ρ ---
    function fρ_simp(γe)
        0 <= γe <= 1 || error("fρ_simp:: γe inválido ($γe)")
        return iρ_ar + (γe^p)*Δiρ
    end

    function dfρ_simp(γe)
        0 <= γe <= 1 || error("dfρ_simp:: γe inválido ($γe)")
        # Cuidado numérico: para p>1 e γe=0, γe^(p-1)=0 (ok);
        # para p=1, expoente 0 dá 1 (recupera o caso linear).
        return p*(γe^(p-1))*Δiρ
    end

    # --- Inverso de κ ---
    function fκ_simp(γe)
        0 <= γe <= 1 || error("fκ_simp:: γe inválido ($γe)")
        return iκ_ar + (γe^p)*Δiκ
    end

    function dfκ_simp(γe)
        0 <= γe <= 1 || error("dfκ_simp:: γe inválido ($γe)")
        return p*(γe^(p-1))*Δiκ
    end

    return (fρ=fρ_simp, fκ=fκ_simp, dfρ=dfρ_simp, dfκ=dfκ_simp)
end
