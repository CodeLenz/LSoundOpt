const KdFactor = typeof(lu(sparse([1], [1], ComplexF64[1.0])))

Cria_Cache_Kd(nf::Integer) = Vector{KdFactor}(undef, nf)
