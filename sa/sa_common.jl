#=
sa_common.jl -- shared pieces of the Sobol' sensitivity runs, used by every script in Code/sa and by
every coral model.

    include(joinpath(@__DIR__, "sa_common.jl"))

Three factor blocks make up one model's design (see Code/sa/README.md):

  1. DHW       one categorical factor `dhw_scenario` picking a member of the shared DHW ensemble
               (built from a Sobol' sample of the 7 generator parameters in `dhw_spec()`);
               the member's 7 parameters are attached to each row after the runs
  2. domain    13 multipliers on the loaded domain (`multiplier_spec()`, `apply_domain_multipliers`),
               bands from multiplier_bands.toml (`multiplier_bands`)
  3. model     the model's own native parameters (e.g. `coralpop_native_spec()`)

Every spec is a DataFrame in the `fieldname`/`val`/`dist`/`dist_params`/`is_constant` shape that
`ADRIA.sample(spec, n)` takes (scrambled Sobol', `n` a power of 2).

Needs ADRIA, DataFrames, Distributions and Statistics loaded by the caller.
=#

using ADRIA, DataFrames, Distributions, Statistics
using DataFrames: select, select!, transform  # ADRIA also exports a `select`
using LinearAlgebra: diagind
import Random, JLD2, TOML

# ---------------------------------------------------------------------------------------------------
# Locations
# ---------------------------------------------------------------------------------------------------

"CoralModelReview root (the folder holding Code/, domains/, CoralPOP/, ...)."
const CMR_DIR = normpath(joinpath(@__DIR__, "..", ".."))

"Shared synthetic domains. `CMR_SA_DOMAINS_DIR` overrides it (the HPC jobs point it at the uploaded copy)."
const SA_DOMAINS_DIR = normpath(get(ENV, "CMR_SA_DOMAINS_DIR", joinpath(CMR_DIR, "domains")))

"Root under which each model's `<ModelFolder>/model_run/` lives. `CMR_SA_RUNS_DIR` overrides it."
const SA_RUNS_ROOT = normpath(get(ENV, "CMR_SA_RUNS_DIR", CMR_DIR))

"The two synthetic domains every model is run on."
const SA_DOMAINS = ["SingleReef_SobolDHW", "NineReefs_SobolDHW"]

"Model name (as passed to the scripts) => its folder in CoralModelReview."
const MODEL_FOLDERS = Dict(
    "coralpop" => "CoralPOP", "coralblox" => "CoralBlox", "cscape" => "CScape", "kora" => "KORA",
    "corset" => "CORSET", "reefmod" => "ReefMod", "coconet" => "CoCoNet"
)

"Climate label of the synthetic domains' DHW file (`DHWs/dhwRCP45.nc`)."
const SA_RCP = "45"

"Simulated years."
const SA_YEARS = (2025, 2050)

domain_path(domain::AbstractString) = joinpath(SA_DOMAINS_DIR, domain)
sa_inputs_dir(domain::AbstractString) = mkpath(joinpath(domain_path(domain), "sa_inputs"))

"Run folder name inside each model folder: `model_run`, or e.g. `model_run_local` for local test runs (`CMR_SA_RUN_SUBDIR`)."
const SA_RUN_SUBDIR = get(ENV, "CMR_SA_RUN_SUBDIR", "model_run")

"`<CMR>/<ModelFolder>/model_run/<domain>` (or `CMR_SA_RUN_SUBDIR` in place of `model_run`), created if needed."
function model_run_dir(model::AbstractString, domain::AbstractString)::String
    haskey(MODEL_FOLDERS, model) || error("Unknown model '$model'; known: $(join(sort(collect(keys(MODEL_FOLDERS))), ", ")).")
    return mkpath(joinpath(SA_RUNS_ROOT, MODEL_FOLDERS[model], SA_RUN_SUBDIR, domain))
end

# ---------------------------------------------------------------------------------------------------
# Factor specs
# ---------------------------------------------------------------------------------------------------

function _spec(rows)::DataFrame
    return DataFrame(
        fieldname=Symbol[r[1] for r in rows], val=Float64[r[2] for r in rows],
        dist=[r[3] for r in rows], dist_params=[r[4] for r in rows], is_constant=Bool[r[5] for r in rows]
    )
end
_uniform(name, lo, hi) = (Symbol(name), (lo + hi) / 2, Uniform, (Float64(lo), Float64(hi)), false)

"""
    DHW_FACTORS

The 7 sampled DHW-generator parameters (SyntheticReefDomain's `DHWConfig`), their ranges and meaning.
For member `s` and year `t`:

    R[t]      = baseline + trend_per_decade * (t - 1) / 10 + anomaly[t] + event[t]
    DHW[t, r] = max(0, R[t] + local_sd * e[t, r])

`anomaly` is AR(1) with stationary SD `interannual_sd` and lag-1 correlation `ar1`; with probability
`event_prob` per year a heatwave adds U(5, 5 + event_mag_width) DHW (the lower bound is fixed at 5);
`e` is standard normal reef-level noise.
"""
const DHW_FACTORS = [
    (:dhw_baseline, 1.75, 5.25, "regional mean DHW in 2025 (DHW)"),
    (:dhw_trend_per_decade, 0.5, 4.0, "linear warming rate of the regional mean (DHW/decade)"),
    (:dhw_interannual_sd, 1.0, 5.0, "stationary SD of the regional year-to-year anomaly (DHW)"),
    (:dhw_ar1, 0.15, 0.45, "lag-1 autocorrelation of that anomaly"),
    (:dhw_event_prob, 0.075, 0.225, "annual probability of a marine-heatwave pulse"),
    (:dhw_event_mag_width, 1.0, 5.0, "spread of heatwave size: a heatwave adds U(5, 5 + width) DHW"),
    (:dhw_local_sd, 1.0, 5.0, "SD of independent reef-level noise (DHW)"),
]
const DHW_EVENT_MAG_MIN = 5.0
const DHW_PARAM_NAMES = [f[1] for f in DHW_FACTORS]

"The DHW-generator parameters as an `ADRIA.sample` spec (7 Uniform factors)."
dhw_spec()::DataFrame = _spec([_uniform(n, lo, hi) for (n, lo, hi, _) in DHW_FACTORS])

"Number of stochastic DHW members in the shared ensemble (the DHW file holds one more: member 1 = all zero)."
const N_DHW_MEMBERS = 1024

"""
    dhw_sample(; n=N_DHW_MEMBERS, seed=1234) -> DataFrame

The Sobol' sample of the 7 DHW parameters behind the shared ensemble: row `i` is DHW member `i + 1`.
"""
function dhw_sample(; n::Int=N_DHW_MEMBERS, seed::Int=1234)::DataFrame
    Random.seed!(seed)
    return ADRIA.sample(dhw_spec(), n)
end

"""
    dhw_selector_spec(members) -> DataFrame

The DHW block of a model design: one factor `dhw_scenario`, uniform over the DHW-file member indices
`members` (e.g. `2:1025`, the stochastic members; member 1 is the zero-DHW control).
"""
function dhw_selector_spec(members::AbstractUnitRange{<:Integer})::DataFrame
    lo, hi = Float64(first(members)), Float64(last(members))
    return _spec([(:dhw_scenario, lo, DiscreteUniform, (lo, hi), false)])
end

"""
    MULTIPLIER_FACTORS

The 13 domain multipliers, applied by [`apply_domain_multipliers`](@ref). Their bands come from
`multiplier_bands.toml` ([`multiplier_bands`](@ref)); the built-in default is `MULTIPLIER_BOUNDS`.
"""
const MULTIPLIER_FACTORS = [
    (:m_depth, "median depth"),
    (:m_cover_total, "total initial coral cover (relative to k)"),
    (:m_cover_fg1, "initial-cover weight of tabular Acropora"),
    (:m_cover_fg2, "initial-cover weight of corymbose Acropora"),
    (:m_cover_fg3, "initial-cover weight of corymbose non-Acropora"),
    (:m_cover_fg4, "initial-cover weight of small massives"),
    (:m_cover_fg5, "initial-cover weight of large massives"),
    (:m_ub, "wave exposure (orbital velocity Ub)"),
    (:m_k, "carrying capacity k (habitable fraction)"),
    (:m_area, "site area"),
    (:m_self_recruit, "connectivity diagonal (self-recruitment)"),
    (:m_conn_strength, "connectivity off-diagonal (between-reef exchange)"),
    (:m_retention, "total larvae retained in the system (mean connectivity row sum)"),
]
const MULTIPLIER_NAMES = [f[1] for f in MULTIPLIER_FACTORS]
const MULTIPLIER_BOUNDS = (0.5, 1.5)

"The multiplier band file read by 03/04; `CMR_SA_MULTIPLIER_BANDS` points at another one."
const MULTIPLIER_BANDS_FILE = normpath(get(ENV, "CMR_SA_MULTIPLIER_BANDS", joinpath(@__DIR__, "multiplier_bands.toml")))

"""
    multiplier_bands(path=MULTIPLIER_BANDS_FILE) -> Dict{Symbol,Tuple{Float64,Float64}}

The (lo, hi) band of every domain multiplier, from a TOML file with an optional `default = [lo, hi]` and
a `[bands]` table of `m_<name> = [lo, hi]` (see Code/sa/multiplier_bands.toml). Multipliers not listed
get `default` (else `MULTIPLIER_BOUNDS`). `lo == hi` fixes the multiplier. Unknown names, a missing
explicitly named file, and bands with lo <= 0 or lo > hi are errors.

Only 03_multiplier_sample.jl and 04_design.jl read it. The bands then travel with the design
(`design.jld2` `spec`), so the runs and 06/07 never depend on this file.
"""
function multiplier_bands(path::AbstractString=MULTIPLIER_BANDS_FILE)::Dict{Symbol,Tuple{Float64,Float64}}
    if !isfile(path)
        haskey(ENV, "CMR_SA_MULTIPLIER_BANDS") && error("CMR_SA_MULTIPLIER_BANDS: $path not found")
        @warn "$path not found: every multiplier gets the built-in band $MULTIPLIER_BOUNDS"
        return Dict(n => MULTIPLIER_BOUNDS for n in MULTIPLIER_NAMES)
    end
    t = TOML.parsefile(path)
    extra = setdiff(keys(t), ["default", "bands"])
    isempty(extra) || error("$path: unknown keys $(join(sort(collect(extra)), ", ")) (expected `default` and [bands])")
    function band(x, what)
        (x isa AbstractVector && length(x) == 2 && all(v -> v isa Real, x)) ||
            error("$path: $what must be [lo, hi], got $(repr(x))")
        lo, hi = Float64(x[1]), Float64(x[2])
        (isfinite(lo) && isfinite(hi) && 0 < lo <= hi) || error("$path: $what needs 0 < lo <= hi, got [$lo, $hi]")
        return (lo, hi)
    end
    default = haskey(t, "default") ? band(t["default"], "default") : MULTIPLIER_BOUNDS
    listed = get(t, "bands", Dict{String,Any}())
    unknown = setdiff(Symbol.(collect(keys(listed))), MULTIPLIER_NAMES)
    isempty(unknown) || error("$path: unknown multipliers $(join(unknown, ", ")); known: $(join(MULTIPLIER_NAMES, ", "))")
    return Dict(n => (haskey(listed, String(n)) ? band(listed[String(n)], String(n)) : default)
                for n in MULTIPLIER_NAMES)
end

"""
    multiplier_spec(bands=multiplier_bands()) -> DataFrame

The domain multipliers as an `ADRIA.sample` spec: `Uniform(lo, hi)` per multiplier, or a constant where
lo == hi. `bands` is a Dict from [`multiplier_bands`](@ref) or one `(lo, hi)` for all 13.
"""
function multiplier_spec(bands::AbstractDict=multiplier_bands())::DataFrame
    function row(n)
        lo, hi = bands[n]
        return lo == hi ? (n, Float64(lo), Uniform, (Float64(lo), Float64(hi)), true) : _uniform(n, lo, hi)
    end
    return _spec([row(n) for n in MULTIPLIER_NAMES])
end
multiplier_spec(band::Tuple{<:Real,<:Real})::DataFrame =
    multiplier_spec(Dict(n => (Float64(band[1]), Float64(band[2])) for n in MULTIPLIER_NAMES))

"Factor-block label of a design column, for grouped summaries."
function factor_block(name)::String
    s = string(name)
    (s == "dhw_scenario" || startswith(s, "dhw_")) && return "DHW"
    startswith(s, "m_") && return "domain"
    return "model"
end

"A spec with its `dist` column as text, so the JLD2/CSV copy does not depend on Distributions' types."
spec_for_disk(spec::DataFrame) = transform(spec, :dist => ByRow(string) => :dist)

"Inverse of [`spec_for_disk`](@ref): the `dist` column back to Distributions types (e.g. a design's saved spec)."
function spec_from_disk(spec::DataFrame)::DataFrame
    to_type(d) = d isa AbstractString ? getfield(Distributions, Symbol(last(split(d, '.')))) : d
    return transform(spec, :dist => ByRow(to_type) => :dist)
end

# ---------------------------------------------------------------------------------------------------
# Domain multipliers
# ---------------------------------------------------------------------------------------------------

"Number of ADRIA size classes per functional group in `init_coral_cover`."
const N_SIZE_CLASSES = 7

"""
    scale_connectivity(C, m_self, m_strength, m_retention) -> (C′, n_capped)

Scale the diagonal of the `[source, sink]` matrix `C` by `m_self` and the off-diagonal by `m_strength`,
then rescale the whole matrix so its mean row sum (the share of larvae that stays anywhere in the
system) is `m_retention` times that of `C`. Rows that would then exceed 1 have their exports cut so the
row sums to 1 (a diagonal above 1 is capped at 1); `n_capped` counts them.
"""
function scale_connectivity(C::AbstractMatrix{<:Real}, m_self::Real, m_strength::Real, m_retention::Real)
    base_total = mean(sum(C; dims=2))
    D = Matrix{Float64}(C) .* m_strength
    D[diagind(D)] .= C[diagind(C)] .* m_self
    cur = mean(sum(D; dims=2))
    cur > 0 && (D .*= base_total * m_retention / cur)
    n_capped = 0
    for i in axes(D, 1)
        rs = sum(@view D[i, :])
        rs <= 1 && continue
        n_capped += 1
        D[i, i] = min(D[i, i], 1.0)
        others = [j for j in axes(D, 2) if j != i]
        exports = sum(@view D[i, others])
        exports > 0 && (D[i, others] .*= (1 - D[i, i]) / exports)
    end
    return D, n_capped
end

"""
    scale_initial_cover(cover, m_total, m_groups) -> cover′

`cover` is ADRIA's `init_coral_cover` data, `[groups*sizes, locations]` relative to k, with each
group's `N_SIZE_CLASSES` size classes contiguous. Per location: the total C becomes min(C * m_total, 1);
the group shares w_g become w_g * m_groups[g], renormalised; every size class of group g is scaled by the
same factor, which keeps the size distribution within the group.
"""
function scale_initial_cover(cover::AbstractMatrix{<:Real}, m_total::Real, m_groups::AbstractVector{<:Real})
    n_rows, n_locs = size(cover)
    n_groups = length(m_groups)
    n_rows == n_groups * N_SIZE_CLASSES || error(
        "init_coral_cover has $n_rows rows; expected $n_groups groups x $N_SIZE_CLASSES sizes."
    )
    out = Matrix{Float64}(cover)
    rows(g) = ((g - 1) * N_SIZE_CLASSES + 1):(g * N_SIZE_CLASSES)
    for l in 1:n_locs
        c = [sum(@view cover[rows(g), l]) for g in 1:n_groups]
        C = sum(c)
        C > 0 || continue
        w = (c ./ C) .* m_groups
        w ./= sum(w)
        C′ = min(C * m_total, 1.0)
        for g in 1:n_groups
            c[g] > 0 && (out[rows(g), l] .*= C′ * w[g] / c[g])
        end
    end
    return out
end

"""
    apply_domain_multipliers(dom, row) -> (dom′, info)

A copy of the ADRIA domain `dom` with the 13 multipliers of `row` (a design row, `NamedTuple` or
`Dict`) applied (see `MULTIPLIER_FACTORS`). `info` holds `n_conn_capped`, the rows of the connectivity
matrix whose sum had to be capped at 1.
"""
function apply_domain_multipliers(dom::ADRIA.Domain, row)
    m(k) = Float64(row[k])
    d = deepcopy(dom)
    ld = d.loc_data
    ld.depth_med = Float64.(ld.depth_med) .* m(:m_depth)
    ld.k = min.(Float64.(ld.k) .* m(:m_k), 1.0)
    ld.area = Float64.(ld.area) .* m(:m_area)
    :ub in propertynames(ld) && (ld.ub = Float64.(ld.ub) .* m(:m_ub))
    :ub_med in propertynames(ld) && (ld.ub_med = Float64.(ld.ub_med) .* m(:m_ub))   # read by C~scape
    d.wave_scens.data .*= m(:m_ub)
    groups = [m(Symbol("m_cover_fg$g")) for g in 1:5]
    d.init_coral_cover.data .= scale_initial_cover(d.init_coral_cover.data, m(:m_cover_total), groups)
    C, n_capped = scale_connectivity(d.conn.data, m(:m_self_recruit), m(:m_conn_strength), m(:m_retention))
    d.conn.data .= C
    return d, (; n_conn_capped=n_capped)
end

# ---------------------------------------------------------------------------------------------------
# DHW ensemble summaries
# ---------------------------------------------------------------------------------------------------

"""
    dhw_member_summaries(dhw) -> DataFrame

Per member of a `[timesteps, locations, members]` DHW array: mean, max, mean number of years per reef
with DHW >= 4 / 8 / 12, and the reef-mean DHW in the final year.
"""
function dhw_member_summaries(dhw::AbstractArray{<:Real,3})::DataFrame
    S = size(dhw, 3)
    years_at(x, thr) = mean(sum(x .>= thr; dims=1))
    return DataFrame(
        member=1:S,
        dhw_mean=[mean(@view dhw[:, :, s]) for s in 1:S],
        dhw_max=[maximum(@view dhw[:, :, s]) for s in 1:S],
        years_dhw_ge4=[years_at(@view(dhw[:, :, s]), 4) for s in 1:S],
        years_dhw_ge8=[years_at(@view(dhw[:, :, s]), 8) for s in 1:S],
        years_dhw_ge12=[years_at(@view(dhw[:, :, s]), 12) for s in 1:S],
        dhw_final_year=[mean(@view dhw[end, :, s]) for s in 1:S],
    )
end

"Load one of the shared synthetic domains in ADRIA."
load_sa_domain(domain::AbstractString) = ADRIA.load_domain(domain_path(domain), SA_RCP)

# ---------------------------------------------------------------------------------------------------
# Reading results (for the Shapley analysis)
# ---------------------------------------------------------------------------------------------------

_set_suffix(set) = set == "main" ? "" : "_$set"

"""
    load_sa_results(model, domain; set="main") -> Dict

`results.jld2` of a model's runs (06_collect.jl): `X` (Shapley inputs), `Y` (scalar outputs), `rows`,
`series`, `reef_final`, `reef_area`, `reef_k`, `years`, `reef_ids`, ... `set` = "main", "b2", ...,
or "control".
"""
function load_sa_results(model, domain; set="main")
    name = set == "control" ? "control_results" : "results" * _set_suffix(set)
    return JLD2.load(joinpath(model_run_dir(model, domain), name * ".jld2"))
end

"""
    load_sa_result_set(model, domain; set="main")

The model's ADRIA result set over the same runs, in the same order as `load_sa_results` (06_collect.jl):
`rs.outcomes[:<indicator>]` is `[timesteps, locations, scenarios]` (or with a `groups` axis for the
`*_taxa_*` ones). Needs the model's package loaded first (`include("Code/sa/models/<model>.jl")`).
"""
function load_sa_result_set(model, domain; set="main")
    name = set == "control" ? "control_result_set" : "result_set" * _set_suffix(set)
    return JLD2.load(joinpath(model_run_dir(model, domain), name * ".jld2"), "result_set")
end

"""
    domain_mean(x, weights) -> Matrix [timesteps, scenarios]

Weighted mean over locations of a `[timesteps, locations, scenarios]` indicator, with per-run weights
`[locations, scenarios]` -- e.g. the habitable area `R["reef_k"] .* R["reef_area"]` from
`load_sa_results`, which differs between runs (domain multipliers).
"""
function domain_mean(x::AbstractArray{<:Real,3}, weights::AbstractMatrix{<:Real})
    T, L, N = size(x)
    size(weights) == (L, N) || error("weights are $(size(weights)); expected ($L, $N)")
    w = weights ./ sum(weights; dims=1)
    return dropdims(sum(Float64.(x) .* reshape(w, 1, L, N); dims=2); dims=2)
end

"Axis names of a result-set outcome (`:timesteps`, `:groups`, `:locations`, `:scenarios`)."
outcome_axes(A) = [p for p in propertynames(A) if p in (:timesteps, :groups, :locations, :scenarios)]

"""
    indicator_output(R, rs, indicator; year=2050, reef=:domain, group=nothing) -> Vector (one value per run)

One Shapley output from the result set `rs` (`load_sa_result_set`) and results `R` (`load_sa_results`),
in the order of `R["X"]`:

- `indicator`: any outcome in `keys(rs.outcomes)`, e.g. `:reef_condition_index`, `:relative_cover`,
  `:relative_shelter_volume`, `:juvenile_indicator`, `:coral_diversity`. The per-morphology ones
  (`:relative_loc_taxa_cover`, `:relative_taxa_cover`, ... with a `groups` axis) need `group`
  (1-3 = CoralPop's Branching, Foliose, Other).
- `year`: a year (2025-2050), `:mean` (mean over all years) or `:min`
- `reef`: `:domain` (habitable-area weighted mean over reefs, with each run's own k x area) or a reef id
  such as `"reef_03"`. `relative_taxa_*` are already domain-level (`reef` must be `:domain`).

    y = indicator_output(R, rs, :reef_condition_index; year=2040)
    y = indicator_output(R, rs, :relative_loc_taxa_cover; year=:mean, group=1)   # Branching
    # SAShE: DataModel(R["X"], y)
"""
function indicator_output(R, rs, indicator::Symbol; year=last(R["years"]), reef=:domain, group=nothing)
    S = indicator_series(R, rs, indicator; reef, group)
    year === :mean && return vec(mean(S; dims=1))
    year === :min && return vec(minimum(S; dims=1))
    t = findfirst(==(year), R["years"])
    t === nothing && error("year $year not simulated ($(first(R["years"]))-$(last(R["years"])))")
    return S[t, :]
end

"""
    indicator_series(R, rs, indicator; reef=:domain, group=nothing) -> Matrix [timesteps, runs]

The whole time series behind [`indicator_output`](@ref) (same `indicator`, `reef` and `group`), one
column per run in the order of `R["X"]`, one row per year of `R["years"]`.
"""
function indicator_series(R, rs, indicator::Symbol; reef=:domain, group=nothing)
    O = rs.outcomes[indicator]
    ax = outcome_axes(O)
    A = Array(O)
    if :groups in ax
        group === nothing && error("$indicator has a groups axis: pass group=1, 2 or 3 (Branching, Foliose, Other)")
        A = selectdim(A, findfirst(==(:groups), ax), group)
        ax = filter(!=(:groups), ax)
    else
        group === nothing || error("$indicator has no groups axis; leave out `group`")
    end
    size(A, ndims(A)) == length(R["rows"]) || error("result set has $(size(A, ndims(A))) runs, results $(length(R["rows"]))")
    if !(:locations in ax)                              # already domain-level: [timesteps, scenarios]
        reef === :domain || error("$indicator is domain-level only; use reef=:domain")
        return Float64.(A)
    elseif reef === :domain
        return domain_mean(A, R["reef_k"] .* R["reef_area"])
    end
    l = findfirst(==(String(reef)), R["reef_ids"])
    l === nothing && error("unknown reef $reef; reefs: $(join(R["reef_ids"], ", "))")
    return Float64.(A[:, l, :])
end
