#=
output_select.jl -- the output of interest of a Shapley analysis, as a [runs, years] matrix, from a
model's existing run outputs (results.jld2, result_set.jld2, indicator_outputs.csv, the domain's DHW).

    include(joinpath(@__DIR__, "output_select.jl"))     # after sa_common.jl and models/<model>.jl
    data = RunData(model, domain)
    out = output_matrix(Dict("indicator" => "relative_cover"), data)    # out.Y, out.years, out.tag, out.label

Output spec (the `[output]` table of a config, or a chain's `intermediate`):

    source    = "result_set"          # default: an ADRIA indicator per year (result_set.jld2)
    indicator = "relative_cover"      # any keys(rs.outcomes): reef_condition_index, relative_shelter_volume,
                                      #   juvenile_indicator, relative_loc_taxa_cover (needs `group`), ...
    reef      = "domain"              # habitable-area weighted mean over reefs, or a reef id ("reef_03")
    group     = 1                     # morphology of the *_taxa_* indicators (CoralPop: 1-3 Branching, Foliose, Other)

    source    = "series"              # results.jld2's domain relative cover per year (no result set needed)
    source    = "indicator_outputs"   # one scalar: `column` of indicator_outputs.csv (e.g. "coral_diversity_2050")
    source    = "dhw"                 # realised DHW of each run's member, `stat` = "annual" | "cumulative" |
                                      #   "max_to_date", `reef` = "domain" (mean over reefs) or a reef id

    years     = "all"                 # or a list, e.g. [2030, 2040, 2050]
    transform = "none"                # "change_from_start" (minus the 2025 value) | "relative_to_start" (divided by it)
    tag       = "..."                 # optional: folder name of the results (default built from the above)
=#

"Lazy access to one model/domain run set (results, result set, DHW are loaded once, when first needed)."
mutable struct RunData
    model::String
    domain::String
    set::String
    dir::String
    R::Dict{String,Any}
    n_rows::Int
    rs::Any
    dhw::Union{Nothing,Array{Float64,3}}
end

"""
    RunData(model, domain; set="main", n_rows=0)

`set` is the run set: `"main"` (results.jld2, result_set.jld2) or `"control"` (the zero-DHW runs,
control_results.jld2, control_result_set.jld2). `n_rows > 0` keeps only the first `n_rows` runs (a nested
Sobol' prefix).
"""
function RunData(model::AbstractString, domain::AbstractString; set::AbstractString="main", n_rows::Integer=0)
    R = load_sa_results(model, domain; set)
    N = nrow(R["X"])
    n = n_rows > 0 ? min(n_rows, N) : N
    return RunData(model, domain, set, model_run_dir(model, domain), R, n, nothing, nothing)
end

rows_used(d::RunData) = 1:d.n_rows
inputs(d::RunData) = d.R["X"][rows_used(d), :]
years(d::RunData) = Int.(d.R["years"])

function result_set(d::RunData)
    d.rs === nothing && (d.rs = load_sa_result_set(d.model, d.domain; set=d.set))
    return d.rs
end

function dhw_array(d::RunData)
    d.dhw === nothing && (d.dhw = Float64.(Array(load_sa_domain(d.domain).dhw_scens)))
    return d.dhw
end

_reef(spec) = (r = get(spec, "reef", "domain"); r == "domain" ? :domain : String(r))
_group(spec) = (g = get(spec, "group", nothing); g === nothing ? nothing : Int(g))

"`[years, runs]` series of a result-set, `series` or `dhw` spec, over all simulated years."
function _series(spec, d::RunData)::Matrix{Float64}
    src = get(spec, "source", "result_set")
    if src == "result_set"
        haskey(spec, "indicator") || error("output spec needs `indicator` (source = \"result_set\")")
        return indicator_series(d.R, result_set(d), Symbol(spec["indicator"]); reef=_reef(spec), group=_group(spec))
    elseif src == "series"
        return Float64.(d.R["series"])
    elseif src == "dhw"
        D = dhw_array(d)                                     # [timesteps, locations, members]
        ## control runs used the all-zero member 1 (their results keep the design's member index)
        members = d.set == "control" ? fill(1, length(d.R["dhw_scenario"])) : Int.(d.R["dhw_scenario"])
        reef = _reef(spec)
        l = reef === :domain ? (1:size(D, 2)) : [findfirst(==(reef), d.R["reef_ids"])]
        l == [nothing] && error("unknown reef $reef; reefs: $(join(d.R["reef_ids"], ", "))")
        A = dropdims(mean(D[:, l, members]; dims=2); dims=2)     # [timesteps, runs]
        stat = get(spec, "stat", "annual")
        stat == "annual" && return A
        stat == "cumulative" && return cumsum(A; dims=1)
        stat == "max_to_date" && return accumulate(max, A; dims=1)
        error("unknown DHW stat '$stat' (annual, cumulative, max_to_date)")
    end
    error("unknown output source '$src' (result_set, series, indicator_outputs, dhw)")
end

"Default folder name of an output spec."
function output_tag(spec)::String
    haskey(spec, "tag") && return String(spec["tag"])
    src = get(spec, "source", "result_set")
    base = src == "result_set" ? String(spec["indicator"]) :
           src == "series" ? "series_relative_cover" :
           src == "indicator_outputs" ? String(spec["column"]) :
           "dhw_" * get(spec, "stat", "annual")
    parts = [base]
    src == "indicator_outputs" || push!(parts, string(get(spec, "reef", "domain")))
    _group(spec) === nothing || push!(parts, "g$(_group(spec))")
    tr = get(spec, "transform", "none")
    tr == "none" || push!(parts, tr)
    return join(parts, "__")
end

"Readable label of an output spec, for figure titles."
function output_label(spec)::String
    tag = output_tag(spec)
    return replace(tag, "__" => ", ", "_" => " ")
end

"""
    output_matrix(spec, data::RunData) -> (; Y, years, tag, label)

`Y` is `[runs, outputs]` (rows in the order of `inputs(data)`, one column per year in `years`; a single
column with year 0 for `indicator_outputs`).
"""
function output_matrix(spec::AbstractDict, d::RunData)
    if get(spec, "source", "result_set") == "indicator_outputs"
        haskey(spec, "column") || error("output spec needs `column` (source = \"indicator_outputs\")")
        d.set == "main" || error("indicator_outputs.csv exists for the main runs only (set '$(d.set)')")
        T = CSV.read(joinpath(d.dir, "indicator_outputs.csv"), DataFrame)
        Int.(T.row) == Int.(d.R["rows"]) || error("indicator_outputs.csv rows differ from results.jld2")
        Y = reshape(Float64.(T[rows_used(d), spec["column"]]), :, 1)
        return (; Y, years=[0], tag=output_tag(spec), label=output_label(spec))
    end
    S = _series(spec, d)[:, rows_used(d)]
    tr = get(spec, "transform", "none")
    S = tr == "none" ? S :
        tr == "change_from_start" ? S .- S[1:1, :] :
        tr == "relative_to_start" ? S ./ S[1:1, :] :
        error("unknown transform '$tr' (none, change_from_start, relative_to_start)")
    yrs = years(d)
    sel = get(spec, "years", "all")
    idx = sel == "all" ? eachindex(yrs) : [something(findfirst(==(y), yrs), 0) for y in sel]
    any(==(0), idx) && error("years $(sel) not all simulated ($(first(yrs))-$(last(yrs)))")
    return (; Y=permutedims(S[idx, :]), years=yrs[idx], tag=output_tag(spec), label=output_label(spec))
end
