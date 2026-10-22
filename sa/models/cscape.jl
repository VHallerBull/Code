#=
models/cscape.jl -- C~scape's part of the sensitivity runs: its factor block and one-row runner.

Included by Code/sa_pf (pf_common.jl) after sa_common.jl, in the Code/sa/models/cscape environment.
Every model file defines:

    MODEL_LABEL                     display name
    model_spec_block() -> DataFrame the model's native factors (ADRIA.sample spec)
    run_row(dom, row) -> NamedTuple one design row on the (unmodified) loaded domain `dom`:
                                    (; scalars::NamedTuple, domain_series::Vector, reef_final::Vector,
                                       reef_area::Vector, reef_k::Vector, out)
    save_row_outputs(fpath, id, r)  write run `id`'s full outputs (`r` from run_row) under `fpath`
    build_result_set(fpath, scenarios, dom, spec) -> result set from those files (06_collect.jl)
    save_result_set(rs, path)

C~scape is only run through the pick-and-freeze design (Code/sa_pf), which keeps no run files, so
`save_row_outputs`, `build_result_set` and `save_result_set` are not provided (they error).

Factors (CScape/CScape_sampled_parameters.md, agreed 2026-10-22): 62 native factors, none fixed --
- 18 posterior draw indices `draw_{growth,surv,fecundity}_<ft>`, DiscreteUniform(1, n_draws), the draws
  used as fitted (no band);
- 36 per-FT constants `{fert_prob,cont_prob,size1yr_mean,size1yr_sd,juvi_cap,bleaching_suscept}_<ft>`,
  defaults from the demography files, ±band (`fert_prob` capped at 1);
- 8 global factors: `recruitment_multiplier`, `brooder_retention`, `heritability`, `heat_tolerance_sd` and
  the four bleaching coefficients, ±band; `bleaching_slope` at most ±20 % (±50 % would multiply the
  mortality term at 8 DHW above tolerance by 4-64 instead of ~16).

Held fixed: region Offshore_Northern, plasticity 0, heritability SD 0, initial tolerance mean 0, tolerance
classes -5..8 DHW (14), seed 3381, no post-IPM scales, temperature-growth off, no cyclones/COTS/external
eggs. With heritability SD 0 and one connectivity matrix the seed has no effect: the model is
deterministic.

Each run: the domain multipliers, `domain_to_cscape_inputs` (the row's DHW member, the row's brooder
retention), one demography per FT from the row's draws and constants (`CScape.build_coral_demography`
on the demography files loaded once per process), a one-row scenario table for the global factors, then
`CScape.initialise_simulation` + `run_all_years!` and `CScape.build_output` -- in memory, no output file.
=#

## ---- demography files: container path first, this machine's as the fallback ----
function _first_existing(paths::AbstractString...)
    i = findfirst(ispath, paths)
    return i === nothing ? last(paths) : paths[i]
end
haskey(ENV, "CSCAPE_STRICT_DEMOGRAPHY") || (ENV["CSCAPE_STRICT_DEMOGRAPHY"] = "true")

import CScape
using CscapeInterface
import ADRIAIndicators

const MODEL_LABEL = "C~scape"

"Folder of the converted demography files (`CSCAPE_DEMOG_DIR`, as Code/engines.jl)."
const CSCAPE_DEMOG_DIR = get(ENV, "CSCAPE_DEMOG_DIR", _first_existing(
    "/opt/coral/data/cscape_demog",
    raw"C:\Users\vhaller\Documents\GitHub\CscapeInterface.jl\sandbox\ADRIAIntegration\projects\AdjustedRecruitment10_Disturbance\data"))
const CSCAPE_DEMOG_STAMP = "ALLDAT_260903"
"Functional types, in the interface's (and every per-FT factor's) order."
const CSCAPE_SA_FTS = String.(CscapeInterface.CSCAPE_FUNCTIONAL_TYPES)
const CSCAPE_DEMOG_NAMES = Dict("acro_table" => "acro_table", "acro_corym" => "acro_corym",
    "corym_non_acro" => "corym_non_acro", "small_massive" => "small_massive",
    "large_massive" => "large_massive", "brooder" => "corym_non_acro_BR")
cscape_demog_path(ft) = joinpath(CSCAPE_DEMOG_DIR, "demog_julia_$(CSCAPE_DEMOG_NAMES[ft])$(CSCAPE_DEMOG_STAMP).rds")

const CSCAPE_SA_REGION = "Offshore_Northern"
const CSCAPE_SA_SCENARIO_ID = 1010001
const CSCAPE_SA_SEED = 3381
"Tolerance classes (`HeatToleranceGroups`, DHW offsets from the population's mean response)."
const CSCAPE_SA_HEAT_TOLERANCE_GROUPS = "-5_8_14"
const CSCAPE_SA_HEAT_CLASSES = collect(range(-5.0, 8.0; length=14))
"Largest relative band of `bleaching_slope` (CScape_sampled_parameters.md, section 2.3)."
const BLEACHING_SLOPE_MAX_BAND = 0.2

## ---- demography files, loaded once per process ----
const _CSCAPE_RAW = Ref{Union{Nothing,Vector{Any}}}(nothing)
"The six demography files (RData dicts, `CSCAPE_SA_FTS` order), read on first use."
function cscape_raw_demography()
    if _CSCAPE_RAW[] === nothing
        _CSCAPE_RAW[] = Any[CScape.safe_rdata_load(cscape_demog_path(ft)) for ft in CSCAPE_SA_FTS]
    end
    return _CSCAPE_RAW[]
end

"Number of posterior draws of each regression of FT `i`: (growth, surv, fecundity)."
function cscape_n_draws(i)
    raw = cscape_raw_demography()[i]
    ng = size(raw["growth_mu"]["draws"], 1)
    size(raw["growth_sigma"]["draws"], 1) == ng ||
        error("$(CSCAPE_SA_FTS[i]): growth mu and sigma have different numbers of draws (one distributional model expected)")
    return (growth=ng, surv=size(raw["surv_mu"]["draws"], 1), fecundity=size(raw["fecundity_eggs"], 1))
end

const _CSCAPE_IS_JUVENILE = Ref{Union{Nothing,Matrix{Bool}}}(nothing)
"Size classes with zero fecundity in every draw, `[ft x sizes]` (as `CScape.load_output`)."
function cscape_is_juvenile()
    if _CSCAPE_IS_JUVENILE[] === nothing
        _CSCAPE_IS_JUVENILE[] = permutedims(reduce(hcat,
            [vec(sum(Matrix{Float64}(r["fecundity_eggs"]); dims=1) .== 0) for r in cscape_raw_demography()]))
    end
    return _CSCAPE_IS_JUVENILE[]
end

"The file defaults of the per-FT constants (`CScape.DEMOGRAPHY_CONSTANTS`) of FT `i`."
function cscape_constant_defaults(i)
    raw = cscape_raw_demography()[i]
    distri = Float64.(raw["Distri_Setter_Cont"])
    return (fert_prob=Float64(raw["fert_prob"]), cont_prob=Float64(raw["cont_prob"]),
        size1yr_mean=distri[1], size1yr_sd=distri[2], juvi_cap=Float64(raw["juvi_cap"]),
        bleaching_suscept=Float64(raw["bleaching_suscept"]))
end

## ---- factors ----
"Global factors: name, default, which band rule (`:direct`, or `:slope` = at most ±20 %)."
const CSCAPE_GLOBAL_FACTORS = [
    (:recruitment_multiplier, 10.0, :direct), (:brooder_retention, 0.03, :direct),
    (:heritability, 0.3, :direct), (:heat_tolerance_sd, 1.91, :direct),
    (:bleaching_intercept, CScape.BLEACHING_COEFS_DEFAULT[1], :direct),
    (:bleaching_slope, CScape.BLEACHING_COEFS_DEFAULT[2], :slope),
    (:bleaching_depth_intercept, CScape.BLEACHING_COEFS_DEFAULT[3], :direct),
    (:bleaching_depth_slope, CScape.BLEACHING_COEFS_DEFAULT[4], :direct),
]
"The bleaching factors in `RunConfig.bleaching_coefs` order."
const CSCAPE_BLEACHING_FACTORS = (:bleaching_intercept, :bleaching_slope, :bleaching_depth_intercept,
    :bleaching_depth_slope)

"""
    cscape_native_spec(; band=0.5) -> DataFrame

C~scape's 62 native factors (module header): posterior draws over all draws whatever the band; the
constants and global factors U(v(1 - b), v(1 + b)), `fert_prob` capped at 1, `bleaching_slope` with
b = min(band, 0.2). Defaults and draw counts are read from the demography files.
"""
function cscape_native_spec(; band::Real=0.5)::DataFrame
    b = Float64(band)
    rows = Any[]
    for kind in (:growth, :surv, :fecundity), (i, ft) in enumerate(CSCAPE_SA_FTS)
        n = getproperty(cscape_n_draws(i), kind)
        push!(rows, (Symbol("draw_$(kind)_$ft"), 1.0, DiscreteUniform, (1.0, Float64(n)), false))
    end
    for c in CScape.DEMOGRAPHY_CONSTANTS, (i, ft) in enumerate(CSCAPE_SA_FTS)
        v = getproperty(cscape_constant_defaults(i), c)
        hi = c === :fert_prob ? min(v * (1 + b), 1.0) : v * (1 + b)
        push!(rows, (Symbol("$(c)_$ft"), v, Uniform, (v * (1 - b), hi), false))
    end
    for (name, v, rule) in CSCAPE_GLOBAL_FACTORS
        bb = rule === :slope ? min(b, BLEACHING_SLOPE_MAX_BAND) : b
        push!(rows, (name, v, Uniform, (v * (1 - bb), v * (1 + bb)), false))
    end
    return _spec(rows)
end

"C~scape's native factors at ±50 % (`cscape_native_spec`)."
model_spec_block()::DataFrame = cscape_native_spec(; band=0.5)

## ---- one run ----
"""
The one-row scenario table of a design row: the fixed settings above plus the row's recruitment
multiplier, heritability (SD 0), initial tolerance SD (mean 0) and bleaching coefficients. The file
names are placeholders: the domain translation and the demography override replace their content.
"""
function cscape_scenario_table(row)::DataFrame
    f(k) = Float64(row[k])
    tbl = DataFrame(
        ID=[CSCAPE_SA_SCENARIO_ID], num_repeats=[1], Region=[CSCAPE_SA_REGION],
        Year_Start=[SA_YEARS[1]], Year_End=[SA_YEARS[2]],
        init_cover=[join(fill("0.1", length(CSCAPE_SA_FTS)), "_")], fts=[join(CSCAPE_SA_FTS, "/")],
        HeatToleranceGroups=[CSCAPE_SA_HEAT_TOLERANCE_GROUPS], HeatToleranceInit=["0_$(f(:heat_tolerance_sd))"],
        Heritability=["$(f(:heritability))_0"], Plasticity=[0.0], temp_growth=[0],
        Disturbance_file=["placeholder.RData"], Connectivity_file=["placeholder.RData"],
        Spatial_file=["placeholder.RData"],
        Growth_Surv_file=[join([basename(cscape_demog_path(ft)) for ft in CSCAPE_SA_FTS], "/")],
        recruitment_multiplier=[f(:recruitment_multiplier)],
    )
    for (col, k) in zip(CScape.BLEACHING_COEF_COLUMNS, CSCAPE_BLEACHING_FACTORS)
        tbl[!, col] = [f(k)]
    end
    return tbl
end

const _CSCAPE_PROJECT = Ref{Union{Nothing,String}}(nothing)
"""
This process's C~scape project folder (a temporary folder, made on first use): only the empty
intervention schedule `RunConfig` setup reads. Every run's inputs are passed in memory.
"""
function cscape_project_dir()
    if _CSCAPE_PROJECT[] === nothing
        dir = mktempdir(; prefix="cscape_pf_")
        mkpath(joinpath(dir, "data"))
        tbl = cscape_scenario_table(Dict(k => v for (k, v, _) in CSCAPE_GLOBAL_FACTORS))
        CScape.write_scenario_interventions(tbl, CSCAPE_SA_SCENARIO_ID, dir)
        _CSCAPE_PROJECT[] = dir
    end
    return _CSCAPE_PROJECT[]
end

"The juvenile-cap log of every run goes nowhere (it would grow by lines per run and year)."
const CSCAPE_NULL_LOG = Sys.iswindows() ? "NUL" : "/dev/null"

"One `CoralDemography` per FT from the row's draws and constants, at the run's site covariates."
function cscape_row_demography(row, spatial::DataFrame)
    raws = cscape_raw_demography()
    return [CScape.build_coral_demography(cscape_demog_path(ft), CSCAPE_SA_REGION, spatial;
                draw_growth=Int(row[Symbol("draw_growth_$ft")]), draw_surv=Int(row[Symbol("draw_surv_$ft")]),
                draw_fecundity=Int(row[Symbol("draw_fecundity_$ft")]),
                constants=NamedTuple(c => Float64(row[Symbol("$(c)_$ft")]) for c in CScape.DEMOGRAPHY_CONSTANTS),
                raw=raws[i])
            for (i, ft) in enumerate(CSCAPE_SA_FTS)]
end

"""
    run_row(dom, row) -> (; scalars, domain_series, reef_final, reef_area, reef_k, out)

Apply the row's domain multipliers to a copy of `dom`, translate it for C~scape with the row's DHW
member and brooder retention, run C~scape with the row's demography and global factors (module header)
and summarise relative cover (coral cover / habitable area kappa x area, 2025-2050):

- `scalars`: `rc_final` (2050), `rc_mean` (2025-2050 mean), `rc_min`, `rc_change` (2050 minus 2025),
  and `n_conn_capped`
- `domain_series`: domain relative cover per year (habitable-area weighted)
- `reef_final`: each reef's relative cover in 2050
- `reef_area`, `reef_k`: each reef's area (m2) and kappa (fraction) in this run, after the multipliers
- `out`: the run's `CscapeOutput` (in memory, with its 2024 initial state)
"""
function run_row(dom, row)
    d, info = apply_domain_multipliers(dom, row)
    inp = domain_to_cscape_inputs(d; dhw_idx=Int(row[:dhw_scenario]), wave_idx=0, cyclone_idx=0,
        brood_site_value=Float64(row[:brooder_retention]))
    demog = cscape_row_demography(row, inp.spatial_override)
    config = CScape.build_scenario_input(cscape_scenario_table(row), CSCAPE_SA_SCENARIO_ID,
        cscape_project_dir(), CSCAPE_NULL_LOG, 1; use_cached_ipm=false, write_interventions=false,
        seed=CSCAPE_SA_SEED)
    state = CScape.initialise_simulation(config; verbose=false, inp..., demography_override=demog)
    CScape.run_all_years!(state)
    out = CscapeOutput(CScape.build_output(state; is_juvenile=cscape_is_juvenile()))

    rows = [findfirst(==(y), out.years) for y in SA_YEARS[1]:SA_YEARS[2]]
    any(isnothing, rows) && error("C~scape output covers $(first(out.years))-$(last(out.years)), not SA_YEARS")
    rc = ADRIAIndicators.relative_cover(CscapeInterface.to_adria_format(out).cover[Int.(rows), :, :, :])  # [T, reefs]
    reef_k, reef_area = copy(out.kappa), out.area[:, 3]
    w = reef_k .* reef_area
    series = vec(rc * w) ./ sum(w)
    scalars = (; rc_final=series[end], rc_mean=mean(series), rc_min=minimum(series),
        rc_change=series[end] - series[1], n_conn_capped=info.n_conn_capped)
    return (; scalars, domain_series=series, reef_final=rc[end, :], reef_area, reef_k, out)
end

save_row_outputs(fpath::AbstractString, id::Integer, r) =
    error("models/cscape.jl: C~scape runs keep no files (pick-and-freeze only, Code/sa_pf)")
build_result_set(fpath, scenarios, dom, spec) =
    error("models/cscape.jl: no result set: C~scape is run through Code/sa_pf only")
save_result_set(rs, path) = error("models/cscape.jl: no result set: C~scape is run through Code/sa_pf only")
