#=
shapley_core.jl -- given-data Shapley effects (SAShE.jl's nearest-neighbour pick-and-freeze estimator)
for many outputs at once: every year of a time series, from one neighbour search.

    include(joinpath(@__DIR__, "shapley_core.jl"))

SAShE's `analyze(DataModel(X, y), n_perm, PickAndFreeze())` draws, per permutation, an order σ of the
players and a reference row s, then for each coalition u = σ[1:i] takes the nearest neighbour nn of s in
X restricted to u and estimates W_u = Var(E[Y | X_u]) by y[s] * y[nn] - mean(y)^2. The neighbours
depend on X, σ, s and u only, never on y. `neighbour_plan` does that (expensive) search once and
`shapley_from_plan` evaluates any number of outputs (columns of a matrix, e.g. the 26 years) from it.
Run unthreaded with the same rng, this reproduces `SAShE.analyze` exactly (test/runtests.jl checks it).

Players can be single inputs (the usual analysis), groups of inputs (group-as-player: a coalition is
the union of its groups' columns) or any subset of the inputs (restricted game, `terminal=:estimate`).
=#

using SAShE, DataFrames, Statistics
using Random: AbstractRNG, Xoshiro, randperm

"""
    NeighbourPlan

The neighbour search behind one Shapley analysis (see `neighbour_plan`): `names[k]` is player `k`, made
of the X columns `players[k]`; for permutation `m`, `σ[:, m]` is the player order, `s[m]` the reference
row and `nn[i, m]` the neighbour of `s[m]` for the coalition `σ[1:i, m]` (0 where the last step uses
Var(Y) instead, `terminal == :var_y`).
"""
struct NeighbourPlan
    names::Vector{String}
    players::Vector{Vector{Int}}
    σ::Matrix{Int}
    s::Vector{Int}
    nn::Matrix{Int}
    terminal::Symbol
    n_rows::Int
end

n_players(p::NeighbourPlan) = length(p.players)
n_permutations(p::NeighbourPlan) = length(p.s)

"Players as single columns: player `j` is column `j` of X, named after it."
singleton_players(names::AbstractVector) = (string.(names), [[j] for j in eachindex(names)])

"""
    neighbour_plan(X, players; names, n_perm, seed=1, terminal=:var_y, threaded=true) -> NeighbourPlan

The nearest-neighbour search of SAShE's `DataModel` pick-and-freeze estimator, over `n_perm` random
permutations of `players` (vectors of column indices of `X`; a `DataFrame` or matrix, one row per run).

- `terminal = :var_y`: the full coalition is every input, so its W is Var(Y) (as in SAShE) and the
  effects sum to Var(Y). Use it when `players` cover all columns of X.
- `terminal = :estimate`: the full coalition is only the players' columns (a restricted game, e.g. the
  parameters of one group); its W = Var(E[Y | X_players]) is estimated by the same neighbour step and the
  effects sum to that estimate.
- `search = :incremental` (default): a coalition only grows, so the squared standardized distances
  from the reference row are accumulated one column at a time and the neighbour is their argmin, with
  ties broken uniformly at random (as SAShE does). Same neighbours as SAShE's search (standardized
  Euclidean, reference row excluded), O(runs x inputs) per permutation instead of a BallTree per step.
  `search = :sashe` calls SAShE's own `_nearest_neighbour_indices` at every step.
- `threaded = true`: permutations run on all threads, permutation `m` with its own rng seeded from `(seed, m)`;
  `false` uses one `Xoshiro(seed)` in sequence. With `search = :sashe` that reproduces
  `SAShE.analyze(...; rng=Xoshiro(seed))` exactly.
"""
function neighbour_plan(
    X, players::Vector{Vector{Int}}; names::Vector{String}=["p$k" for k in eachindex(players)],
    n_perm::Integer, seed::Integer=1, terminal::Symbol=:var_y, search::Symbol=:incremental,
    threaded::Bool=true, verbose::Bool=false,
)::NeighbourPlan
    terminal in (:var_y, :estimate) || throw(ArgumentError("terminal must be :var_y or :estimate"))
    Xm = Matrix{Float64}(X isa DataFrame ? Matrix(X) : X)
    n_rows, n_cols = size(Xm)
    p = length(players)
    p >= 1 || throw(ArgumentError("need at least one player"))
    all(c -> all(in(1:n_cols), c), players) || throw(ArgumentError("player columns outside 1:$n_cols"))
    terminal == :var_y && p < 2 && throw(ArgumentError("terminal=:var_y needs at least 2 players"))
    search in (:incremental, :sashe) || throw(ArgumentError("search must be :incremental or :sashe"))
    bad = findfirst(j -> any(isnan, @view Xm[:, j]), 1:n_cols)
    bad === nothing || throw(ArgumentError("column $bad of X contains NaN"))
    Xs = SAShE._standardize_columns(Xm)     # the scaling SAShE applies to every coalition

    σ = zeros(Int, p, n_perm)
    s = zeros(Int, n_perm)
    nn = zeros(Int, p, n_perm)
    last_step = terminal == :var_y ? p - 1 : p
    done = Threads.Atomic{Int}(0)
    report = max(1, n_perm ÷ 10)

    function one_permutation!(m, rng)
        σ[:, m] = randperm(rng, p)
        s[m] = rand(rng, 1:n_rows)
        if search == :sashe
            for i in 1:last_step
                cols = reduce(vcat, players[@view σ[1:i, m]])   # == σ[1:i] for single-column players
                nn[i, m] = SAShE._nearest_neighbour_indices(Xm, s[m], cols, 1; rng=rng)[1]
            end
        else
            _incremental_neighbours!(view(nn, :, m), Xs, players, view(σ, :, m), s[m], last_step, rng)
        end
        k = Threads.atomic_add!(done, 1) + 1
        verbose && k % report == 0 && println("    $(k)/$(n_perm) permutations")
    end

    if threaded
        Threads.@threads :dynamic for m in 1:n_perm
            one_permutation!(m, Xoshiro(hash((seed, m))))
        end
    else
        rng = Xoshiro(seed)
        for m in 1:n_perm
            one_permutation!(m, rng)
        end
    end
    return NeighbourPlan(names, players, σ, s, nn, terminal, n_rows)
end

"""
Neighbours of row `s` for the growing coalitions of one permutation `σ`: `d2` accumulates the squared
standardized distance of every row to `s` over the columns added so far; the neighbour after step `i` is
the row (not `s`) with the smallest `d2`, ties broken uniformly at random (reservoir sampling).
"""
function _incremental_neighbours!(nn_m, Xs, players, σ_m, s, last_step, rng)
    n_rows = size(Xs, 1)
    d2 = zeros(n_rows)
    for i in 1:last_step
        for c in players[σ_m[i]]
            xs = Xs[s, c]
            @inbounds @simd for r in 1:n_rows
                d2[r] += (Xs[r, c] - xs)^2
            end
        end
        best, n_ties, pick = Inf, 0, 0
        @inbounds for r in 1:n_rows
            r == s && continue
            v = d2[r]
            if v < best
                best, n_ties, pick = v, 1, r
            elseif v == best
                n_ties += 1
                rand(rng, 1:n_ties) == 1 && (pick = r)
            end
        end
        nn_m[i] = pick
    end
    return nn_m
end

"""
    ShapleyResult

Shapley effects of `names` (players) on each output column `t`: `effect[k, t]` with its 95 % bounds
`lb`/`ub` (SAShE's `shapley_effects`; anti-conservative for given data, see SAShE's docs), `var_y[t]`
the output variance and `total[t]` the sum of the effects (= `var_y` for a full analysis, the estimated
Var(E[Y | X_players]) for a restricted game).
"""
struct ShapleyResult
    names::Vector{String}
    effect::Matrix{Float64}
    lb::Matrix{Float64}
    ub::Matrix{Float64}
    var_y::Vector{Float64}
    total::Vector{Float64}
end

"Effects as shares of the output variance (NaN where the variance is 0)."
share_of_variance(r::ShapleyResult) = r.effect ./ _nonzero(r.var_y)'
"Effects as shares of their own total (within a restricted game; NaN where the total is 0)."
share_of_total(r::ShapleyResult) = r.effect ./ _nonzero(r.total)'
_nonzero(v) = [x == 0 ? NaN : x for x in v]

"""
    shapley_from_plan(plan, Y; sums=Dict()) -> (ShapleyResult, Dict{String,ShapleyResult})

Shapley effects of the players of `plan` on each column of `Y` (`[runs, outputs]`, rows in the order of
the X the plan was built on; a vector is one output). `sums` maps a name to a grouping of the players
(`(group_names, [player indices per group])`); for each, the members' per-permutation increments are
added before the confidence bounds, giving the grouped-sum effects with valid bounds.
"""
function shapley_from_plan(plan::NeighbourPlan, Y::AbstractVecOrMat{<:Real};
    sums::AbstractDict=Dict{String,Tuple{Vector{String},Vector{Vector{Int}}}}())
    Ym = Matrix{Float64}(Y isa AbstractVector ? reshape(Y, :, 1) : Y)
    size(Ym, 1) == plan.n_rows || throw(ArgumentError("Y has $(size(Ym, 1)) rows, the plan $(plan.n_rows)"))
    p, n_perm, T = n_players(plan), n_permutations(plan), size(Ym, 2)
    n_perm > 1 || throw(ArgumentError("need at least 2 permutations for confidence bounds"))

    alloc() = (zeros(p, T), zeros(p, T), zeros(p, T))
    E, L, U = alloc()
    S = Dict(k => (zeros(length(g[1]), T), zeros(length(g[1]), T), zeros(length(g[1]), T)) for (k, g) in sums)
    var_y, total = zeros(T), zeros(T)
    Φ = zeros(p, n_perm)
    Φ² = zeros(p, n_perm)
    for t in 1:T
        y = Ym[:, t]
        Ȳ = mean(y)
        vy = var(y)
        fill!(Φ, 0.0)
        fill!(Φ², 0.0)
        for m in 1:n_perm
            s = plan.s[m]
            prevW = 0.0
            for i in 1:p
                W = plan.nn[i, m] == 0 ? vy : y[s] * y[plan.nn[i, m]] - Ȳ^2
                Δ = W - prevW
                k = plan.σ[i, m]
                Φ[k, m] = Δ / n_perm
                Φ²[k, m] = Δ^2 / n_perm
                prevW = W
            end
        end
        E[:, t], L[:, t], U[:, t] = SAShE.shapley_effects(Φ, Φ²)
        var_y[t] = vy
        total[t] = sum(@view E[:, t])
        for (key, (_, members)) in sums
            Φg = reduce(vcat, [sum(@view(Φ[g, :]); dims=1) for g in members])
            Φg² = (Φg .* n_perm) .^ 2 ./ n_perm
            S[key][1][:, t], S[key][2][:, t], S[key][3][:, t] = SAShE.shapley_effects(Φg, Φg²)
        end
    end
    main = ShapleyResult(plan.names, E, L, U, var_y, total)
    grouped = Dict(k => ShapleyResult(sums[k][1], v[1], v[2], v[3], var_y, total) for (k, v) in S)
    return main, grouped
end

"""
    chain_analysis(X, Y, M, upstream, others; upstream_names, other_names, m_name="M", lag=0,
                   n_perm, seed=1, threaded=true) -> (; A, B, path, path_unexplained, m_index)

Follow-on analysis inputs -> intermediate `M` -> output `Y` (both `[runs, years]`, same years).

- **A** (inputs -> M): restricted game of the `upstream` columns of `X` on `M`: each upstream input's
  effect on M, summing to Var(E[M | X_upstream]) (`share_of_total(A)` = its share).
- **B** (M -> Y): for each year t, Shapley effects on `Y[:, t]` of the players `others` (column sets of
  `X`, which together must be every column not in `upstream`) plus `M[:, t - lag]` in place of the
  upstream columns. `B.effect[m_index, t]` is Sh_Y(M).
- **path**: `path[i, t]` = Sh_Y(M)(t) x upstream input i's effect on M / Var(M) (year t - lag): Y's
  variance attributed to input i through M. `path_unexplained[t]` = Sh_Y(M)(t) x (1 - Var(E[M | upstream])
  / Var(M)): the part through M that the upstream inputs do not describe (e.g. the weather realisation of
  a DHW trajectory). Dividing by Var(M) rather than by the noisy estimate of Var(E[M | upstream]) keeps
  this stable.

Years with t - lag < 1 are NaN.
"""
function chain_analysis(
    X, Y::AbstractMatrix, M::AbstractMatrix, upstream::Vector{Int}, others::Vector{Vector{Int}};
    upstream_names::Vector{String}, other_names::Vector{String}, m_name::String="M", lag::Integer=0,
    n_perm::Integer, seed::Integer=1, threaded::Bool=true, verbose::Bool=false,
)
    Xm = Matrix{Float64}(X isa DataFrame ? Matrix(X) : X)
    N, n_cols = size(Xm)
    T = size(Y, 2)
    size(M) == size(Y) || throw(ArgumentError("M is $(size(M)), Y $(size(Y)): they need the same runs and years"))
    keep = setdiff(1:n_cols, upstream)
    sort(reduce(vcat, others; init=Int[])) == keep ||
        throw(ArgumentError("`others` must cover every column not in `upstream`, each once"))

    planA = neighbour_plan(Xm, [[c] for c in upstream]; names=upstream_names, n_perm, seed,
        terminal=:estimate, threaded)
    A, _ = shapley_from_plan(planA, M)
    shareA = share_of_variance(A)
    unexplA = 1 .- A.total ./ _nonzero(A.var_y)

    pos = Dict(c => k for (k, c) in enumerate(keep))
    playersB = [[[pos[c] for c in g] for g in others]; [[length(keep) + 1]]]
    namesB = [other_names; m_name]
    p = length(playersB)
    E, L, U = fill(NaN, p, T), fill(NaN, p, T), fill(NaN, p, T)
    var_y, total = fill(NaN, T), fill(NaN, T)
    path = fill(NaN, length(upstream), T)
    path_unexplained = fill(NaN, T)
    for t in 1:T
        tM = t - lag
        tM >= 1 || continue
        Xt = hcat(Xm[:, keep], M[:, tM])
        planB = neighbour_plan(Xt, playersB; names=namesB, n_perm, seed, terminal=:var_y, threaded)
        r, _ = shapley_from_plan(planB, Y[:, t])
        E[:, t], L[:, t], U[:, t] = r.effect, r.lb, r.ub
        var_y[t], total[t] = r.var_y[1], r.total[1]
        path[:, t] = r.effect[p] .* shareA[:, tM]
        path_unexplained[t] = r.effect[p] * unexplA[tM]
        verbose && println("    chain year $t/$T done")
    end
    B = ShapleyResult(namesB, E, L, U, var_y, total)
    return (; A, B, path, path_unexplained, m_index=p)
end
