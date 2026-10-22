#=
Tests of the Shapley pipeline (shapley_core.jl, groups.jl, and indicator_series in sa_common.jl).

    julia +release -t auto --project=Code/sa/shapley Code/sa/shapley/test/runtests.jl

The given-data estimator is noisy, so the analytic checks average a few replicates and use loose
tolerances; the exact check is the bit-for-bit equality with SAShE.analyze.
=#

using Test, SAShE, DataFrames, Statistics
using Random: Xoshiro
include(joinpath(@__DIR__, "..", "shapley_core.jl"))
include(joinpath(@__DIR__, "..", "groups.jl"))

"Independent normal inputs with SDs `sd`, `N` rows, as a DataFrame x1, x2, ..."
normal_inputs(rng, N, sd) = DataFrame(randn(rng, N, length(sd)) .* sd', ["x$i" for i in eachindex(sd)])

"Mean effects over `reps` datasets: `f(rng) -> ShapleyResult` (first output column)."
average(f, reps) = mean(f(Xoshiro(100 + r)).effect[:, 1] for r in 1:reps)

@testset "Shapley pipeline" begin
    @testset "plan reproduces SAShE.analyze exactly" begin
        rng = Xoshiro(3)
        X = DataFrame(rand(rng, 300, 4), [:a, :b, :c, :d])
        Y = hcat(X.a .+ 2 .* X.b .* X.c, sin.(3 .* X.d) .+ X.a, X.c .^ 2)
        names_, players = singleton_players(names(X))
        plan = neighbour_plan(X, players; names=names_, n_perm=200, seed=11, threaded=false, search=:sashe)
        r, _ = shapley_from_plan(plan, Y)
        for t in 1:3
            Φₙ, Φ²ₙ = SAShE.analyze(DataModel(X, Y[:, t]), 200, PickAndFreeze(); rng=Xoshiro(11))
            Φ, lb, ub = SAShE.shapley_effects(Φₙ, Φ²ₙ)
            @test r.effect[:, t] == Φ && r.lb[:, t] == lb && r.ub[:, t] == ub
            @test sum(r.effect[:, t]) ≈ var(Y[:, t])
        end
        # threaded: other random streams, same estimator
        rt, _ = shapley_from_plan(neighbour_plan(X, players; names=names_, n_perm=200, seed=11,
            search=:sashe), Y)
        @test rt.total ≈ r.total
        # the incremental search finds the same neighbours (continuous inputs: no ties)
        inc = neighbour_plan(X, players; names=names_, n_perm=200, seed=11)
        Xm = Matrix(X)
        @test all(inc.nn[i, m] == SAShE._nearest_neighbour_indices(Xm, inc.s[m], inc.σ[1:i, m], 1)[1]
                  for m in 1:200 for i in 1:3)
        ## with ties (repeated rows, as the reused DHW members): uniform among them, never the row itself
        Xt = DataFrame(a=repeat(1.0:10.0, 30), b=rand(rng, 300))
        pt = neighbour_plan(Xt, [[1], [2]]; n_perm=2000, seed=2, terminal=:estimate)
        ms = [m for m in 1:2000 if pt.σ[1, m] == 1]          # column a entered first
        first_nn = pt.nn[1, ms]
        @test all(first_nn .!= pt.s[ms]) && all(Xt.a[first_nn] .== Xt.a[pt.s[ms]])
        @test length(unique(first_nn)) > 200
    end

    ## Y = x1 + 2 x2 + 0.5 x3 + x4 with SDs (1, 1, 2, 1): Shapley effects = (beta sd)^2 = (1, 4, 1, 1)
    β, sd = [1.0, 2.0, 0.5, 1.0], [1.0, 1.0, 2.0, 1.0]
    truth = (β .* sd) .^ 2
    lin(X) = Matrix(X) * β
    N, n_perm, reps = 3000, 1500, 3

    @testset "individual, grouped sum, group-as-player" begin
        groups = (["g12", "g34"], [[1, 2], [3, 4]])
        res = [begin
            X = normal_inputs(Xoshiro(100 + r), N, sd)
            nm, pl = singleton_players(names(X))
            ind, sums = shapley_from_plan(neighbour_plan(X, pl; names=nm, n_perm, seed=r),
                lin(X); sums=Dict("s" => groups))
            gp, _ = shapley_from_plan(neighbour_plan(X, groups[2]; names=groups[1], n_perm, seed=r), lin(X))
            (ind.effect[:, 1], sums["s"].effect[:, 1], gp.effect[:, 1], sum(ind.effect), sums["s"])
        end for r in 1:reps]
        ind, gsum, gplay = mean(first.(res)), mean(getindex.(res, 2)), mean(getindex.(res, 3))
        @test all(abs.(ind .- truth) .< 0.2 .* truth)
        @test res[1][2] ≈ [sum(res[1][1][1:2]), sum(res[1][1][3:4])]   # sum of the members' effects
        @test all(abs.(gsum .- [5, 2]) .< 0.2 .* [5, 2])
        @test all(abs.(gplay .- [5, 2]) .< 0.15 .* [5, 2])
        s = res[1][5]
        @test sum(s.effect) ≈ res[1][4]                           # grouped sums add up to Var(Y)
        @test all(s.lb .<= s.effect .<= s.ub)
    end

    @testset "restricted game" begin
        ## players x1, x2 only; their effects sum to Var(E[Y | x1, x2]) = 1 + 4
        e = average(reps) do rng
            X = normal_inputs(rng, N, sd)
            first(shapley_from_plan(neighbour_plan(X, [[1], [2]]; n_perm, terminal=:estimate), lin(X)))
        end
        @test all(abs.(e .- [1, 4]) .< 0.2 .* [1, 4])
    end

    @testset "chain" begin
        ## M = x1 + 2 x2 + e (upstream x1, x2; e ~ N(0,1) is not an input, like a weather realisation),
        ## Y = M + 0.5 x3 + x4: Sh_Y(M) = Var(M) = 6, paths x1 -> 1, x2 -> 4, unexplained part of M -> 1
        out = [begin
            rng = Xoshiro(200 + r)
            X = normal_inputs(rng, N, sd)
            M = reshape(X.x1 .+ 2 .* X.x2 .+ randn(rng, N), :, 1)
            Y = M .+ reshape(0.5 .* X.x3 .+ X.x4, :, 1)
            chain_analysis(X, Y, M, [1, 2], [[3], [4]]; upstream_names=["x1", "x2"],
                other_names=["x3", "x4"], n_perm, seed=r)
        end for r in 1:reps]
        shY = mean(o.B.effect[o.m_index, 1] for o in out)
        path = mean(o.path[:, 1] for o in out)
        @test abs(shY - 6) < 0.15 * 6
        @test all(abs.(path .- [1, 4]) .< 0.25 .* [1, 4])
        @test abs(mean(o.path_unexplained[1] for o in out) - 1) < 0.4
        @test all(abs.(mean(o.B.effect[1:2, 1] for o in out) .- [1, 1]) .< 0.25)
        lagged = chain_analysis(normal_inputs(Xoshiro(1), 500, sd), rand(500, 3), rand(500, 3), [1, 2],
            [[3], [4]]; upstream_names=["x1", "x2"], other_names=["x3", "x4"], n_perm=20, lag=1)
        @test all(isnan, lagged.B.effect[:, 1]) && !any(isnan, lagged.B.effect[:, 2])
    end

    @testset "groups file" begin
        dir = mktempdir()
        write(joinpath(dir, "base.toml"), """
            [[scheme]]
            name = "s"
              [[scheme.group]]
              name = "A"
              match = ["^a_"]
            """)
        write(joinpath(dir, "ok.toml"), """
            extends = "base.toml"
            [[scheme]]
            name = "s"
              [[scheme.group]]
              name = "B"
              match = ["^b_", "^c\$"]
              [[scheme.group]]
              name = "Empty"
              match = ["^zzz"]
            """)
        G = load_groupings(joinpath(dir, "ok.toml"), ["a_1", "b_1", "c", "a_2"])
        @test G["s"].names == ["A", "B"] && G["s"].members == [[1, 4], [2, 3]]
        @test group_of(G["s"], 3) == "B"
        @test_throws ErrorException load_groupings(joinpath(dir, "ok.toml"), ["a_1", "d"])      # d in no group
        write(joinpath(dir, "dup.toml"), """
            extends = "base.toml"
            [[scheme]]
            name = "s"
              [[scheme.group]]
              name = "B"
              match = ["_1\$"]
            """)
        @test_throws ErrorException load_groupings(joinpath(dir, "dup.toml"), ["a_1", "b_1"])  # a_1 twice
        ## the real CoralPop file covers the 175 inputs
        cols = split(readline(joinpath(@__DIR__, "..", "..", "..", "..", "CoralPOP", "model_run_local",
            "SingleReef_SobolDHW", "results.csv")), ",")
        X_names = cols[3:findfirst(==("rc_final"), cols) - 1]
        @test length(X_names) == 175
        GC = load_groupings(joinpath(@__DIR__, "..", "groups", "coralpop.toml"), X_names)
        @test length(GC["block"].names) == 3 && length(GC["process"].names) == 12
    end
end

## indicator_series vs indicator_output on the local CoralPop runs (needs ADRIA + CoralPopInterface)
const LOCAL = joinpath(@__DIR__, "..", "..", "..", "..", "CoralPOP", "model_run_local", "SingleReef_SobolDHW", "result_set.jld2")
if isfile(LOCAL)
    ENV["CMR_SA_RUN_SUBDIR"] = "model_run_local"
    include(joinpath(@__DIR__, "..", "..", "sa_common.jl"))
    include(joinpath(@__DIR__, "..", "..", "models", "coralpop.jl"))
    @testset "indicator_series" begin
        R = load_sa_results("coralpop", "SingleReef_SobolDHW")
        rs = load_sa_result_set("coralpop", "SingleReef_SobolDHW")
        S = indicator_series(R, rs, :relative_cover)
        @test size(S) == (length(R["years"]), nrow(R["X"]))
        @test all(S[t, :] == indicator_output(R, rs, :relative_cover; year=y) for (t, y) in enumerate(R["years"]))
        @test vec(mean(S; dims=1)) == indicator_output(R, rs, :relative_cover; year=:mean)
    end
end
