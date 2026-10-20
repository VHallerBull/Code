#=
models/coralpop.jl -- CoralPop's part of the Sobol' runs: its factor block and one-row runner.

Included by Code/sa/04-07 after sa_common.jl (run in the Code/sa/models/coralpop environment).
Every model file defines:

    MODEL_LABEL                     display name
    model_spec_block() -> DataFrame the model's native factors (ADRIA.sample spec)
    run_row(dom, row) -> NamedTuple one design row on the (unmodified) loaded domain `dom`:
                                    (; scalars::NamedTuple, domain_series::Vector, reef_final::Vector,
                                       reef_area::Vector, reef_k::Vector, out)
    save_row_outputs(fpath, id, r)  write run `id`'s full outputs (`r` from run_row) under `fpath`
    build_result_set(fpath, scenarios, dom, spec) -> result set from those files (06_collect.jl)
    save_result_set(rs, path)

Full outputs (05_run_chunk.jl, per run, under <run dir>/runs_<set>/):
- adria_exports/CoralPop_indicators_scenario_<row>_seed_<seed>.jld2  every ADRIA indicator per year and
  reef (`CMR_SA_KEEP_INDICATORS`, default true); 06 stacks them into result_set.jld2
- model_outputs/CoralPop_scenario_<row>_seed_<seed>.jld2  CoralPop's raw output: per-bin colony area,
  cover, benthos, larvae (`CMR_SA_KEEP_RAW`, default true)

CoralPop runs with every spawner larva (a reef's own included) travelling through the domain's
connectivity (`larvae_via_connectivity=true`), so the connectivity multipliers act on all spawner
recruitment. The 17 local-larvae factors are therefore held at their defaults by the spec (155 varying).
=#

## ---- the Python engine: container paths first, this machine's as the fallback ----
function _first_existing(paths::AbstractString...)
    i = findfirst(ispath, paths)
    return i === nothing ? last(paths) : paths[i]
end
ENV["JULIA_CONDAPKG_BACKEND"] = "Null"
haskey(ENV, "JULIA_PYTHONCALL_EXE") || (ENV["JULIA_PYTHONCALL_EXE"] = _first_existing(
    "/opt/pyenv/bin/python", raw"C:\Users\vhaller\.julia\conda\3\x86_64\envs\coral\python.exe"))
haskey(ENV, "CORALPOP_ENGINE_DIR") || (ENV["CORALPOP_ENGINE_DIR"] = _first_existing(
    "/opt/sa/repos/CoralPop-Engine", raw"C:\Users\vhaller\Documents\GitHub\CoralPop-Engine"))

using CoralPopInterface

const MODEL_LABEL = "CoralPop"
const LARVAE_VIA_CONNECTIVITY = true
const KEEP_INDICATORS = get(ENV, "CMR_SA_KEEP_INDICATORS", "true") == "true"
const KEEP_RAW = get(ENV, "CMR_SA_KEEP_RAW", "true") == "true"

"CoralPop's native factors (`coralpop_native_spec`, +-50 %, connectivity larvae, 2025-2050)."
model_spec_block()::DataFrame = coralpop_native_spec(;
    band=0.5, larvae_via_connectivity=LARVAE_VIA_CONNECTIVITY, rcp=parse(Int, SA_RCP),
    year_start=SA_YEARS[1], year_end=SA_YEARS[2]
)

"""
    run_row(dom, row) -> (; scalars, domain_series, reef_final)

Apply the row's domain multipliers to a copy of `dom`, translate it for CoralPop with the row's DHW
member, run CoralPop with the row's native parameters and summarise relative cover (coral cover as a
fraction of the habitable area k * area, after the multipliers):

- `scalars`: `rc_final` (2050), `rc_mean` (2025-2050 mean), `rc_min`, `rc_change` (2050 minus 2025),
  final-year % cover of each CoralPop morphology (area-weighted), and `n_conn_capped`
- `domain_series`: domain relative cover per year
- `reef_final`: each reef's relative cover in 2050
- `reef_area`, `reef_k`: each reef's area (m2) and k (fraction) in this run, after the multipliers
- `out`: the full `CoralPopOutput`, for [`save_row_outputs`](@ref)
"""
function run_row(dom, row)
    d, info = apply_domain_multipliers(dom, row)
    inp = domain_to_coralpop_inputs(d; dhw_idx=Int(row[:dhw_scenario]), wave_idx=0, cyclone_idx=0)
    out = run_coralpop_scenario(inp, coralpop_scenario_inputs(row);
        larvae_via_connectivity=LARVAE_VIA_CONNECTIVITY)

    cover = out.cover_pct                       # [T1, 3 morphologies, L], % of reef area
    area = out.reef_area_m2
    habitable = sum(out.k_pct ./ 100 .* area)
    coral_area = dropdims(sum(cover ./ 100 .* reshape(area, 1, 1, :); dims=(2, 3)); dims=(2, 3))
    series = coral_area ./ habitable
    reef_final = vec(sum(cover[end, :, :]; dims=1)) ./ out.k_pct
    morph_final = [sum(cover[end, m, :] .* area) / sum(area) for m in axes(cover, 2)]

    scalars = (;
        rc_final=series[end], rc_mean=mean(series), rc_min=minimum(series),
        rc_change=series[end] - series[1],
        (Symbol("cover_final_", lowercase(CORALPOP_MORPHS[m])) => morph_final[m] for m in eachindex(morph_final))...,
        n_conn_capped=info.n_conn_capped,
    )
    return (; scalars, domain_series=series, reef_final, reef_area=copy(area), reef_k=out.k_pct ./ 100, out)
end

"""
    save_row_outputs(fpath, id, r) -> Nothing

Write run `id`'s ADRIA indicators (`adria_exports/`) and raw CoralPop output (`model_outputs/`) under
`fpath`, in the file layout `load_coralpop_results` reads. `CMR_SA_KEEP_INDICATORS` / `CMR_SA_KEEP_RAW`
(default true) switch each off.
"""
function save_row_outputs(fpath::AbstractString, id::Integer, r)
    seed = Int(r.out.metadata["seed"])
    KEEP_INDICATORS && save_coralpop_indicators(coralpop_indicators_path(fpath, id, seed),
        coralpop_indicators(r.out), r.out)
    KEEP_RAW && save_coralpop_output(coralpop_output_path(fpath, id, seed), r.out)
    return nothing
end

"""
    build_result_set(fpath, scenarios, dom, spec) -> CoralPopResultSet

Stack the indicator files under `fpath/adria_exports/` into one `CoralPopResultSet` (an ADRIA
`ResultSet`: `rs.outcomes[:reef_condition_index]` etc. are `[timesteps, locations, scenarios]`, the
`ADRIA.metrics` functions work on it). `scenarios` must have an `ID` column (the design row) and is
stored as `rs.inputs`; scenarios come out sorted by ID.

The result set's `loc_area` / `loc_max_coral_cover` are those of the first run only: every run has
its own (domain multipliers), stored per run in results.jld2 (`reef_area`, `reef_k`).
"""
build_result_set(fpath, scenarios, dom, spec) =
    load_coralpop_results(fpath; scenarios, dom, spec, scenario_ids=Int.(scenarios.ID),
        name="CoralPop SA $(basename(dirname(fpath)))", show_progress=false)

save_result_set(rs, path) = save_coralpop_results(rs, path)
