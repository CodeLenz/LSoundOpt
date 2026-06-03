using LSoundOpt
using LinearAlgebra
using SparseArrays
using TOML
using Test
using YAML

const ROOT = dirname(@__DIR__)

@testset "Project metadata" begin
    project = TOML.parsefile(joinpath(ROOT, "Project.toml"))

    @test haskey(project, "sources")
    @test project["sources"]["LSound"] == Dict(
        "url" => "https://github.com/CodeLenz/LSound",
        "rev" => "main",
    )
    @test project["sources"]["Lgmsh"] == Dict(
        "url" => "https://github.com/CodeLenz/Lgmsh",
        "rev" => "main",
    )
end

@testset "Public API is discrete-only" begin
    @test isdefined(LSoundOpt, :Otim_ISLP)
    @test !isdefined(LSoundOpt, :Otim_SIMP)
    @test !isdefined(LSoundOpt, :LP_Continuo)
    @test !isdefined(LSoundOpt, :Cria_Param_SIMP)
end

@testset "YAML input parsing" begin
    mktempdir() do dir
        yaml = joinpath(dir, "case.yaml")
        write(yaml, """
        versao: 1.0
        volfrac: "0.75"
        raio: "0.0125"
        fatorcv: "0.15"
        perimetro: "42"
        niter: "7"
        """)

        raio, niter, vf, perimetro, mu, fatorcv = LSoundOpt.Le_YAML(yaml)

        @test raio == 0.0125
        @test niter == 7
        @test vf == 0.75
        @test perimetro == 42.0
        @test mu == 0.0
        @test fatorcv == 0.15
    end
end

@testset "YAML algorithm parameters" begin
    mktempdir() do dir
        yaml = joinpath(dir, "case.yaml")
        write(yaml, """
        versao: 1.0
        volfrac: "0.75"
        raio: "0.0125"
        fatorcv: "0.15"
        algoritmo:
          tr_expand_ratio: "0.8"
          tr_shrink_factor: "0.25"
          tr_max_reject: "3"
          slack_safety_factor: "5.0"
          heur_stag_l2_factor: "2.0"
          hard_stop_stagnation: "6"
        """)

        alg = LSoundOpt.Le_YAML_Algoritmo(yaml)

        @test alg.tr_expand_ratio == 0.8
        @test alg.tr_shrink_ratio == 0.3
        @test alg.tr_shrink_factor == 0.25
        @test alg.tr_max_reject == 3
        @test alg.slack_safety_factor == 5.0
        @test alg.heur_stag_l2_factor == 2.0
        @test alg.hard_stop_stagnation == 6

        default_yaml = joinpath(dir, "default.yaml")
        write(default_yaml, """
        versao: 1.0
        volfrac: 0.75
        raio: 0.0125
        fatorcv: 0.15
        """)

        defaults = LSoundOpt.Le_YAML_Algoritmo(default_yaml)
        @test defaults.tr_expand_ratio == 0.7
        @test defaults.tr_max_reject == 10
        @test defaults.hard_stop_stagnation == 4
    end
end

@testset "LShape YAML template is complete" begin
    yaml = joinpath(ROOT, "msh", "LShape", "lshape.yaml")
    dados = YAML.load_file(yaml)
    alg = dados["algoritmo"]

    @test LSoundOpt.Le_YAML(yaml) == (0.008, 100, 0.99, 0.0, 1.83e-5, 0.15)
    @test LSoundOpt.Le_YAML_Algoritmo(yaml) isa LSoundOpt.Algorithm_Params

    expected_algorithm_keys = [
        "tr_accept_ratio_min",
        "tr_expand_ratio",
        "tr_shrink_ratio",
        "tr_expand_factor",
        "tr_shrink_factor",
        "tr_cv_max",
        "tr_max_reject",
        "tr_perimeter_tol",
        "tr_physics_tol",
        "tr_pred_tol",
        "tr_ratio_pred_tol",
        "slack_safety_factor",
        "slack_offset",
        "heur_cv_max",
        "heur_stag_l1_factor",
        "heur_stag_l2_factor",
        "heur_stag_l3_factor",
        "heur_min_moves",
        "hard_stop_stagnation",
    ]

    @test all(key -> haskey(alg, key), expected_algorithm_keys)
end

@testset "SPL objective" begin
    p0 = 1.0
    @test LSoundOpt.SPLn([1.0 + 0im, -1.0 + 0im], p0) ≈ 0.0 atol=1e-12
    @test LSoundOpt.SPLn([10.0 + 0im, -10.0 + 0im], p0) ≈ 20.0 atol=1e-12

    MP = ComplexF64[
        1.0 10.0
        2.0 20.0
        -1.0 -10.0
    ]
    nodes_target = [1, 3]

    @test LSoundOpt.Objetivo(MP, nodes_target, [1.0, 1.0], p0) ≈ 10.0 atol=1e-12
    @test LSoundOpt.Objetivo(MP, nodes_target, [2.0, 1.0], p0) ≈ 10.0 atol=1e-12
    @test LSoundOpt.Objetivo(MP, nodes_target, [1.0, 2.0], p0) ≈ 20.0 atol=1e-12
end

@testset "Discrete LP step" begin
    delta = LSoundOpt.LP(
        [-1.0, 1.0],
        spzeros(0, 2),
        Float64[],
        [0.0, 1.0],
    )

    @test delta == [1.0, -1.0]

    constrained_delta = LSoundOpt.LP(
        [-1.0, 1.0],
        sparse([1.0 0.0]),
        [0.0],
        [0.0, 1.0],
    )

    @test constrained_delta == [0.0, -1.0]
end

@testset "Slack penalty scaling" begin
    c = [-3.0, 2.0, -1.0]
    beta_small = LSoundOpt.Slack_Penalty(c, 1)
    beta_large = LSoundOpt.Slack_Penalty(c, 20)

    @test beta_large > beta_small
    @test beta_large > 2 * 20 * maximum(abs.(c))
    @test LSoundOpt.Slack_Penalty(Float64[], 20) == 10.0
end

@testset "Algorithmic heuristic parameters" begin
    alg = LSoundOpt.Algorithm_Params(
        heur_cv_max=0.4,
        heur_stag_l1_factor=0.5,
        heur_stag_l2_factor=2.0,
        heur_stag_l3_factor=4.0,
        heur_min_moves=3.0,
    )

    @test LSoundOpt.Update_Heuristics(1, 1, 0.1, 100, 0.1, alg) == 0.05
    @test LSoundOpt.Update_Heuristics(1, 2, 0.1, 100, 0.1, alg) == 0.2
    @test LSoundOpt.Update_Heuristics(1, 3, 0.1, 100, 0.2, alg) == 0.4
    @test LSoundOpt.Update_Heuristics(1, 0, 0.02, 100, 0.15, alg) == 0.15
end

@testset "Constraint linearization dimensions" begin
    elements_design = [1, 2, 3]
    ne = 3
    V = [1.0, 2.0, 3.0]
    gamma = [0.0, 1.0, 0.0]
    neighedge = [[2], [1, 3], [2]]
    map_global_local = Dict(1 => 1, 2 => 2, 3 => 3)

    A, b = LSoundOpt.Lineariza_Restricoes(
        V, elements_design, 3.0, 2.0, 0.0, 0.0, ne, gamma, neighedge,
        map_global_local;
        restricao_volume=true,
    )

    @test size(A) == (1 + 2 * length(elements_design), 2 * length(elements_design))
    @test length(b) == size(A, 1)
    @test Matrix(A[1:1, :]) == [1.0 2.0 3.0 0.0 0.0 0.0]
    @test b[1] == 1.0

    A_no_volume, b_no_volume = LSoundOpt.Lineariza_Restricoes(
        V, elements_design, 3.0, 2.0, 0.0, 0.0, ne, gamma, neighedge,
        map_global_local;
        restricao_volume=false,
    )

    @test size(A_no_volume) == (2 * length(elements_design), 2 * length(elements_design))
    @test length(b_no_volume) == size(A_no_volume, 1)
end
