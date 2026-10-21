#=
Tests of the pick-and-freeze design (pf_common.jl).

    julia +release --project=Code/sa/shapley Code/sa_pf/test/runtests.jl
=#

using Test, Statistics
using LinearAlgebra: I
const MODEL = "coralpop"
include(joinpath(@__DIR__, "..", "pf_common.jl"))

@testset "pick-and-freeze design" begin
    D = pf_design(0.5, 64)
    d = ncol(D.A)

    @testset "walks equal SAShE's pick-and-freeze table" begin
        S = CallablePickAndFreezeSample(D.A, D.B, D.perms)
        for n in (1, 17, 64)
            W = walk_rows(D, n)
            block = S.samples[((n - 1) * (d + 1) + 1):(n * (d + 1)), :]
            @test Matrix(select(W, names(D.A))) == Matrix{Float64}(block)
        end
        W = walk_rows(D, 5)
        @test all(W.dhw_scenario .== PF_DHW_SLOT)
        @test all(k -> all(W[!, k] .== D.constants[k]), keys(D.constants))
        @test Vector(W[1, names(D.A)]) == Vector{Float64}(D.A[5, :]) && Vector(W[end, names(D.A)]) == Vector{Float64}(D.B[5, :])
    end

    @testset "permutations: uniform, shared order common to every model" begin
        P = pf_permutations(4000, 21, 176)
        @test all(sort(P[n, :]) == 1:176 for n in 1:4000)
        Q = pf_permutations(4000, 21, 60)                     # another model: other d
        rel(P, n) = filter(<=(21), P[n, :])
        @test all(rel(P, n) == rel(Q, n) for n in 1:4000)     # same walk order of the shared inputs
        pos1 = [findfirst(==(1), P[n, :]) for n in 1:4000]    # where input 1 falls: uniform on 1..176
        @test abs(mean(pos1) - 88.5) < 3 && abs(std(pos1) - std(1:176)) < 3
    end

    @testset "zero-DHW design" begin
        Z = pf_design(0.5, 64; set="zero")
        dz = ncol(Z.A)
        @test dz == d - 8 && !any(startswith("dhw_"), names(Z.A)) && Z.n_shared == 13
        @test Z.A == select(D.A, Not(r"^dhw_")) && Z.B == select(D.B, Not(r"^dhw_"))   # same parameter sets
        S = CallablePickAndFreezeSample(Z.A, Z.B, Z.perms)
        W = walk_rows(Z, 9)
        @test Matrix(select(W, names(Z.A))) == Matrix{Float64}(S.samples[(8 * (dz + 1) + 1):(9 * (dz + 1)), :])
        @test all(W.dhw_scenario .== 1.0)                     # the all-zero member
        @test all(walk_rows(D, 9).dhw_scenario .== PF_DHW_SLOT)
        ## the 13 shared (domain) inputs walk in the same order as in any other model's zero design
        rel(P, n) = filter(<=(13), P[n, :])
        Q = pf_permutations(64, 13, 40)
        @test all(rel(Z.perms, n) == rel(Q, n) for n in 1:64)
    end

    @testset "A and B independent, inputs independent" begin
        Dc = pf_design(0.5, 4096)
        c = [cor(Dc.A[!, k], Dc.B[!, k]) for k in names(Dc.A)]
        @test maximum(abs.(c)) < 4 / sqrt(4096)              # A[n] and B[n] unrelated
        C = cor(Matrix(Dc.A))
        @test maximum(abs.(C - I)) < 5 / sqrt(4096)         # no input dependence within A
    end

    @testset "estimator on a toy model with known effects" begin
        ## Y = linear in 4 inputs of the real design (independent uniforms): Sh_i = beta_i^2 Var(x_i)
        Dt = pf_design(0.5, 2048)
        cols = ["dhw_trend_per_decade", "m_cover_total", "dhw_weather", "growth_rate_Other"]
        β = [1.0, 4.0, 0.005, 0.5]
        f(r) = sum(β[i] * r[cols[i]] for i in 1:4)
        S = CallablePickAndFreezeSample(Dt.A, Dt.B, Dt.perms)
        Y = [f(r) for r in eachrow(S.samples)]
        Φ, Φ², _ = SAShE.analyze(S, Y)
        sh = SAShE.shapley_effects(Φ)
        moe = SAShE.margin_of_error(Φ, Φ²)
        truth = [β[i]^2 * var(vcat(Dt.A[!, cols[i]], Dt.B[!, cols[i]])) for i in 1:4]
        idx = [findfirst(==(c), names(Dt.A)) for c in cols]
        ## noise is about the same for every input, so test against each effect's own 95 % bound (x 1.5)
        @test all(abs.(sh[idx] .- truth) .< 1.5 .* moe[idx])
        others = setdiff(1:ncol(Dt.A), idx)
        @test count(abs.(sh[others]) .> moe[others]) <= 0.1 * length(others)   # no effect: inside bounds
        @test abs(sum(sh) - sum(truth)) < 0.1 * sum(truth)
    end

    @testset "outputs from the state summary (sa_outputs.jl)" begin
        @test native_group_class([[3], [1, 2], [4, 5]]) == [2, 1, 3]                # CoralPop
        @test native_group_class([[2], [1], [3], [5], [4], Int[]]) == [1, 1, 2, 3, 3, 0]   # CoCoNet, tt -> none
        @test_throws ErrorException native_group_class([[2, 4]])                     # spans two classes
        ## 2 years x 2 reefs, reef 2 has three times the habitable area
        taxa = zeros(2, 3, 2)
        taxa[:, 1, :] .= 0.1; taxa[:, 2, :] .= 0.2; taxa[:, 3, :] .= [0.1 0.3; 0.1 0.3]   # Branching, Foliose, Other
        cc = class_covers(taxa, [2, 1, 3])
        @test cc[:, 1, :] == taxa[:, 2, :] && cc[:, 2, :] == taxa[:, 1, :] && cc[:, 3, :] == taxa[:, 3, :]
        S = Dict("cover" => dropdims(sum(taxa; dims=2); dims=2), "class_cover" => cc,
            "juveniles" => fill(0.05, 2, 2), "diversity" => [0.4 0.6; 0.4 0.6], "weights" => [1.0, 3.0],
            "thermal_mean" => [5.0 6.0; 7.0 8.0], "thermal_weight" => [1.0 1.0; 1.0 0.0])
        Y = derive_outputs(S; thermal_groups=["a", "b"])
        @test Y["relative_cover"] ≈ [0.25 * 0.4 + 0.75 * 0.6, 0.25 * 0.4 + 0.75 * 0.6]
        tot = 0.25 * 0.4 + 0.75 * 0.6
        @test Y["acropora_share"] ≈ fill(0.2 / tot, 2) && Y["massive_share"] ≈ fill((0.25 * 0.1 + 0.75 * 0.3) / tot, 2)
        @test Y["coral_diversity"] ≈ fill(0.25 * 0.4 + 0.75 * 0.6, 2) && Y["relative_juveniles"] ≈ fill(0.05, 2)
        @test Y["thermal_a"] == [5.0, 7.0] && isequal(Y["thermal_b"], [6.0, NaN])      # no weight -> missing
        @test Y["spatial_cv"] ≈ fill(std([0.4, 0.6]) / 0.5, 2)
        @test !haskey(derive_outputs(Dict(k => (k == "weights" ? v[1:1] : v[:, ntuple(_ -> :, ndims(v) - 2)..., 1:1])
            for (k, v) in S if k in ("cover", "class_cover", "weights"))), "spatial_cv")   # one reef: no CV
        @test output_order(keys(Y); thermal_groups=["b", "a"]) ==
              ["relative_cover", "acropora_share", "massive_share", "relative_juveniles", "coral_diversity",
               "thermal_b", "thermal_a", "spatial_cv"]
        ## store -> chunk -> merge -> derive: the same outputs back
        store = state_store(S, 3, 2)
        store_state!(store, S, 2, 1)
        @test all(isnan, store["cover"][1, :, :, 1]) && store["cover"][2, :, :, 1] == Float32.(S["cover"])
        state = merge_states([Dict("S" => store, "base" => [4, 1])], 4, 3)
        @test size(state["class_cover"]) == (4, 3, 2, 3, 2)
        Ys = derive_stored(state; thermal_groups=["a", "b"])
        @test isapprox(Ys["relative_cover"][4, 2, :], Y["relative_cover"]; rtol=1e-6)   # stored as Float32
        @test all(isnan, Ys["relative_cover"][4, 1, :])
        @test all(isnan, Ys["relative_cover"][2, :, :])                               # base sample never stored
        ## years with a missing value are left out; thermal years together
        Yv = Dict("relative_cover" => ones(2, 3, 2), "thermal_a" => ones(2, 3, 2), "thermal_b" => ones(2, 3, 2))
        Yv["thermal_b"][1, 2, 2] = NaN
        v = valid_years(Yv, collect(keys(Yv)), [1, 2])
        @test v["relative_cover"] == [1, 2] && v["thermal_a"] == [1] && v["thermal_b"] == [1]
        @test valid_years(Yv, collect(keys(Yv)), [2])["thermal_a"] == [1, 2]          # walk 1 not analysed
    end

    @testset "CoralPop state summary on real runs" begin
        for domain in SA_DOMAINS
            dom = load_sa_domain(domain)
            draws = load_dhw_draws(joinpath(PF_INPUTS, domain, "dhw_draws.jld2"))
            r = pf_run(dom, walk_rows(D, 1)[1, :], draws)
            St = r.state
            L = size(St["cover"], 2)
            @test size(St["class_cover"]) == (26, 3, L) && length(St["weights"]) == L
            @test maximum(abs.(dropdims(sum(St["class_cover"]; dims=2); dims=2) .- St["cover"])) < 1e-12
            Yr = derive_outputs(St; thermal_groups=PF_THERMAL_GROUPS)
            @test all(0 .<= Yr["acropora_share"] .<= 1) && all(0 .<= Yr["massive_share"] .<= 1)
            @test haskey(Yr, "spatial_cv") == (L > 1) && !any(is_thermal, keys(Yr))   # CoralPop: no tolerance
            z = pf_run_zero(dom, zero_rows(D, 1)[1, :])
            @test derive_outputs(z.state)["relative_cover"][1] ≈ Yr["relative_cover"][1]   # same initial state
        end
    end
end
