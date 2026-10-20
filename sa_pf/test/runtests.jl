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
end
