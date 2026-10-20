#=
sa_outputs.jl -- the outputs of the sensitivity analysis, defined the same way for every coral model.

    include(joinpath(@__DIR__, "sa_outputs.jl"))      # no package dependencies beyond Statistics

Every run keeps only a small **state summary** (a model's adapter, adapters/<model>.jl, builds it from
the run's in-memory output; no raw files are kept), and the outputs are derived from it afterwards:

    "cover"          [T, L]     relative cover per reef: coral cover / habitable area (k x area)
    "class_cover"    [T, 3, L]  relative cover per harmonised class (HARMONISED_CLASSES)
    "juveniles"      [T, L]     ADRIAIndicators relative_juveniles, the model's own juvenile classes
                                (absent for models without colony sizes)
    "diversity"      [T, L]     ADRIAIndicators coral_diversity over the model's own groups (absent
                                where not provided)
    "weights"        [L]        habitable area (m2) of each reef in this run (after the multipliers)
    "thermal_mean"   [T, G]     mean thermal tolerance per native group, on the model's own scale and
    "thermal_weight" [T, G]     with its own weighting (absent for models without heritable tolerance)

Outputs (one value per year; `derive_outputs`):

    relative_cover      habitable-area weighted mean of "cover"
    acropora_share      Acropora cover / total cover (domain level)
    massive_share       massive cover / total cover (domain level)
    relative_juveniles  habitable-area weighted mean of "juveniles"
    coral_diversity     habitable-area weighted mean of "diversity"
    thermal_<group>     "thermal_mean" of each native group (NaN where its weight is 0)
    spatial_cv          SD / mean of "cover" across the reefs (domains with more than one reef)
=#

using Statistics

"Harmonised coral classes every model with morphology can provide (CoralPop cannot split them further)."
const HARMONISED_CLASSES = ["Acropora", "corymbose_non_Acropora", "massive"]
"Class of each of ADRIA's 5 functional groups: tabular + corymbose Acropora, corymbose non-Acropora, small + large massives."
const ADRIA_GROUP_CLASS = [1, 1, 2, 3, 3]

"""
    native_group_class(adria_groups) -> Vector{Int}

The harmonised class of each native group, from the ADRIA groups it maps onto (one vector per native
group, as the interfaces' `*_ADRIA_GROUPS` tables). 0 for a group that maps onto none (e.g. CoCoNet's
intervention-only `tt`); an error if a group spans two classes.
"""
function native_group_class(adria_groups)::Vector{Int}
    return map(collect(adria_groups)) do gs
        gs = collect(gs)
        isempty(gs) && return 0
        cs = unique(ADRIA_GROUP_CLASS[gs])
        length(cs) == 1 || error("native group maps onto ADRIA groups $gs, which span classes $cs")
        cs[1]
    end
end

"""
    class_covers(taxa, class_of_group) -> Array [T, 3, L]

Sum the native groups' relative cover `taxa` (`[T, G, L]`) into the harmonised classes.
"""
function class_covers(taxa::AbstractArray{<:Real,3}, class_of_group::AbstractVector{<:Integer})
    T, G, L = size(taxa)
    G == length(class_of_group) || error("$G groups in the cover, $(length(class_of_group)) in the class map")
    cc = zeros(T, length(HARMONISED_CLASSES), L)
    for g in 1:G
        c = class_of_group[g]
        c > 0 && (cc[:, c, :] .+= @view taxa[:, g, :])
    end
    return cc
end

"""
    output_order(names; thermal_groups=String[]) -> Vector{String}

Output names in display order: the fixed outputs, then the thermal ones in the model's group order
(`thermal_groups`), then spatial_cv, then anything else alphabetically.
"""
function output_order(names; thermal_groups=String[])
    names = String.(collect(names))
    fixed = ["relative_cover", "acropora_share", "massive_share", "relative_juveniles", "coral_diversity"]
    out = [n for n in fixed if n in names]
    append!(out, [thermal_name(g) for g in thermal_groups if thermal_name(g) in names])
    append!(out, sort([n for n in names if is_thermal(n) && !(n in out)]))
    "spatial_cv" in names && push!(out, "spatial_cv")
    append!(out, sort([n for n in names if !(n in out)]))
    return out
end

is_thermal(name) = startswith(String(name), "thermal_")
thermal_name(group) = "thermal_" * String(group)

"Readable label of an output."
function output_label(name)::String
    n = String(name)
    is_thermal(n) && return "Thermal tolerance: " * n[9:end]
    return get(Dict(
        "relative_cover" => "Relative cover", "acropora_share" => "Acropora share",
        "massive_share" => "Massive share", "relative_juveniles" => "Relative juveniles",
        "coral_diversity" => "Coral diversity", "spatial_cv" => "Between-reef CV of cover",
    ), n, n)
end

"""
    derive_outputs(S; thermal_groups=String[]) -> Dict{String,Vector{Float64}}

The outputs of one run (one value per year) from its state summary `S` (a Dict with the keys of this
file's header). `thermal_groups` names the columns of "thermal_mean".
"""
function derive_outputs(S::AbstractDict; thermal_groups=String[])::Dict{String,Vector{Float64}}
    w = Float64.(S["weights"])
    wn = w ./ sum(w)
    dmean(A) = vec(sum(Float64.(A) .* reshape(wn, 1, :); dims=2))           # [T, L] -> [T]
    out = Dict{String,Vector{Float64}}()
    cover = S["cover"]
    out["relative_cover"] = dmean(cover)
    cls = [dmean(@view S["class_cover"][:, c, :]) for c in eachindex(HARMONISED_CLASSES)]
    tot = reduce(+, cls)
    out["acropora_share"] = cls[1] ./ tot                                    # NaN where there is no cover
    out["massive_share"] = cls[3] ./ tot
    haskey(S, "juveniles") && (out["relative_juveniles"] = dmean(S["juveniles"]))
    haskey(S, "diversity") && (out["coral_diversity"] = dmean(S["diversity"]))
    if haskey(S, "thermal_mean")
        M, W = S["thermal_mean"], S["thermal_weight"]
        size(M, 2) == length(thermal_groups) || error("thermal_mean has $(size(M, 2)) groups, $(length(thermal_groups)) names")
        for (g, name) in enumerate(thermal_groups)
            out[thermal_name(name)] = [W[t, g] > 0 ? Float64(M[t, g]) : NaN for t in axes(M, 1)]
        end
    end
    if size(cover, 2) > 1
        out["spatial_cv"] = [std(@view cover[t, :]) / mean(@view cover[t, :]) for t in axes(cover, 1)]
    end
    return out
end

## ---------------------------------------------------------------------------------------------------
## Storage of many runs' state summaries
## ---------------------------------------------------------------------------------------------------

"""
    state_store(S, n_runs, n_walks) -> Dict key => Array{Float32} [n_runs, size(S[key])..., n_walks]

Empty (NaN) storage for the state summaries of `n_runs` runs per walk and `n_walks` walks, shaped like
the summary `S` of one run.
"""
state_store(S::AbstractDict, n_runs::Integer, n_walks::Integer) =
    Dict(k => fill(NaN32, n_runs, size(v)..., n_walks) for (k, v) in S)

"Store run `i` of walk `j`."
function store_state!(store::AbstractDict, S::AbstractDict, i::Integer, j::Integer)
    for (k, v) in S
        store[k][i, ntuple(_ -> Colon(), ndims(v))..., j] .= v
    end
    return store
end

"""
    derive_stored(state; thermal_groups) -> Dict output => Array{Float64} [n_walks, n_runs, T]

The outputs of every stored run (`state`: key => `[n_walks, n_runs, dims...]`, as 04_collect.jl saves
it). A walk whose summary is missing (NaN cover) gives NaN outputs.
"""
function derive_stored(state::AbstractDict; thermal_groups=String[])
    N, R = size(state["cover"], 1), size(state["cover"], 2)
    T = size(state["cover"], 3)
    Y = Dict{String,Array{Float64,3}}()
    for n in 1:N, i in 1:R
        S = Dict(k => v[n, i, ntuple(_ -> Colon(), ndims(v) - 2)...] for (k, v) in state)
        any(isnan, S["cover"]) && continue
        for (k, y) in derive_outputs(S; thermal_groups)
            haskey(Y, k) || (Y[k] = fill(NaN, N, R, T))
            Y[k][n, i, :] = y
        end
    end
    return Y
end

"""
    merge_states(parts, N, n_runs) -> Dict key => Array{Float32} [N, n_runs, dims...]

The chunks' state stores (`S`: key => [n_runs, dims..., n_here]) in base-sample order.
"""
function merge_states(parts, N::Integer, n_runs::Integer)
    i = findfirst(p -> p["S"] !== nothing, parts)
    i === nothing && error("no chunk has any run's state summary")
    T0 = parts[i]["S"]
    state = Dict(k => fill(NaN32, N, n_runs, size(v)[2:(end - 1)]...) for (k, v) in T0)
    for p in parts
        p["S"] === nothing && continue
        for (j, n) in enumerate(p["base"]), (k, v) in p["S"]
            state[k][n, ntuple(_ -> Colon(), ndims(v) - 1)...] = selectdim(v, ndims(v), j)
        end
    end
    return state
end

"""
    valid_years(Y, outs, ok) -> Dict output => year indices to analyse

A year with a missing value (NaN) in any analysed walk is left out of that output. For the thermal
outputs a year is left out of all of them together (a group without weight in that year drops the year).
"""
function valid_years(Y, outs, ok)
    good(o) = [t for t in axes(Y[o], 3) if !any(isnan, @view Y[o][ok, :, t])]
    v = Dict(o => good(o) for o in outs)
    th = filter(is_thermal, outs)
    if !isempty(th)
        common = reduce(intersect, [v[o] for o in th])
        foreach(o -> v[o] = common, th)
    end
    return v
end
