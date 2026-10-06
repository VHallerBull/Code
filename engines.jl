#=
engines.jl -- the per-engine code shared by `08_RUN_ENGINE.jl` (run ONE engine, e.g. as one HPC job)
and `09_ANALYSE_SIX_ENGINES.jl` (load every engine's saved output and compare).

    include(joinpath(@__DIR__, "engines.jl"))

The run/load logic is lifted from `06_RUNNING_A_NEW_PROJECT.jl` (CoralBlox, C~scape) and
`07_COMPARE_SIX_ENGINES.jl` (CoralPop, Kora, CORSET, ReefMod) rather than rewritten, so an engine
here is driven exactly the way those scripts drive it. What is new: every engine writes to its OWN
folder, so the six can run as six independent jobs without touching each other's files:

    ModelRuns/<project>/<engine>/
        result.jld2 (or result_path.txt for CoralBlox's ADRIA result store)
        samples.jld2         scenarios, varying factors, factor spec, and the `settings` used
        benchmark.csv        load / warm-up / main-run time, allocations, peak memory, disk
        size_structure.jld2  domain-level size distribution over time, where the engine has one
        run/                 the engine's own output folder

Configuration is read from `CMR_<NAME>` environment variables (see `engine_config`), with defaults
that work both inside the HPC container (`/opt/coral/...`) and on this Windows machine.

This file sets the Python (CoralPop) and MATLAB (ReefMod) environment variables BEFORE loading any
package -- PythonCall reads them at load time and Distributed workers inherit them -- so include it
before any `using` of the engine packages.
=#

@isdefined(MODEL_RUNS_DIR) || include(joinpath(@__DIR__, "paths.jl"))

_cfg(name, default) = get(ENV, "CMR_" * name, default)

"The first of `paths` that exists, else the last one (so an error names a sensible location)."
function _first_existing(paths::AbstractString...)
    i = findfirst(ispath, paths)
    return i === nothing ? last(paths) : paths[i]
end

## ---- external interpreters: container paths first, this machine's as the fallback ----
ENV["JULIA_CONDAPKG_BACKEND"] = "Null"
haskey(ENV, "JULIA_PYTHONCALL_EXE") || (ENV["JULIA_PYTHONCALL_EXE"] = _first_existing(
    "/opt/pyenv/bin/python", raw"C:\Users\vhaller\.julia\conda\3\x86_64\envs\coral\python.exe"))
haskey(ENV, "CORALPOP_ENGINE_DIR") || (ENV["CORALPOP_ENGINE_DIR"] = _first_existing(
    "/opt/coral/repos/CoralPop-Engine", raw"C:\Users\vhaller\Documents\GitHub\CoralPop-Engine"))
haskey(ENV, "REEFMOD_MATLAB_EXE") || (ENV["REEFMOD_MATLAB_EXE"] = "matlab")
haskey(ENV, "REEFMOD_GBR_DIR") || (ENV["REEFMOD_GBR_DIR"] = _first_existing(
    "/opt/coral/repos/ReefMod-GBR", raw"C:\Users\vhaller\Documents\GitHub\ReefMod-GBR"))

using ADRIA, CscapeInterface, CoralPopInterface, KoraInterface, CORSETInterface, ReefModInterface,
      DataFrames, Statistics, Printf, Dates, Random
import CScape, JLD2, XLSX, MAT

# =============================================================================
# Configuration
# =============================================================================

"""
    engine_config() -> NamedTuple

Every setting an engine run depends on, from `CMR_<NAME>` environment variables. Defaults: the ADRIA
Test domain at its median-severity climate member (DHW 40 / wave 16 of 50, computed as the
k-area-weighted mean over years and sites), 8 Sobol samples per engine, ADRIA's default CoralBlox
parameters (no calibration file), serial engines (one scenario at a time -- the fairest basis for
the per-scenario benchmarks).
"""
function engine_config()
    calib = _cfg("CALIBRATION_PATH", "")
    ## CMR_NO_DISTURBANCE=true: no heat stress and no cyclones for any engine (ADRIA's member 0 =
    ## "none"), i.e. intrinsic dynamics only. The wave member stays: no engine here treats waves as a
    ## disturbance (CoralBlox's wave damage is disabled), but C~scape derives each site's wave exposure
    ## (ub_med) from it, which must not differ between a disturbed and an undisturbed run.
    no_disturbance = _cfg("NO_DISTURBANCE", "false") == "true"
    return (
        project=_cfg("PROJECT_NAME", "TestDomain_N8"),
        domain_path=_cfg("DOMAIN_PATH", _first_existing(
            "/opt/coral/data/Test_domain",
            raw"C:\Users\vhaller\Documents\GitHub\ADRIA.jl\ADRIA\test\data\Test_domain")),
        rcp=_cfg("RCP", "45"),
        dhw_idx=no_disturbance ? 0 : parse(Int, _cfg("DHW_IDX", "40")),
        wave_idx=parse(Int, _cfg("WAVE_IDX", "16")),
        cyclone_idx=no_disturbance ? 0 : parse(Int, _cfg("CYCLONE_IDX", "1")),
        n_samples=parse(Int, _cfg("N_SAMPLES", "8")),
        ## ADRIA's Sobol' sampler is Owen-scrambled, i.e. random: seeding it makes a run's scenarios reproducible.
        seed=parse(Int, _cfg("SEED", "1234")),
        calibration_path=calib,
        ## 0 = run in this process / one scenario at a time. C~scape always needs >= 1 Distributed worker.
        n_workers=parse(Int, _cfg("N_WORKERS", "0")),
        ## Threads each C~scape worker gets (CScape.jl's per-site `Threads.@threads` loop).
        cscape_threads=parse(Int, _cfg("CSCAPE_THREADS", "1")),
        ## A 1-scenario run before the timed run: separates time-to-first-result (compilation, worker
        ## and MATLAB start-up) from steady-state cost in the benchmarks.
        warmup=_cfg("WARMUP", "true") == "true",
        ## The six converted `demog_julia_*` files C~scape reads (R/brms is never run here).
        cscape_demog_dir=_cfg("CSCAPE_DEMOG_DIR", _first_existing(
            "/opt/coral/data/cscape_demog",
            raw"C:\Users\vhaller\Documents\GitHub\CscapeInterface.jl\sandbox\ADRIAIntegration\projects\AdjustedRecruitment10_Disturbance\data")),
        ## CORSET's consumer biomass: PLACEHOLDERS (kg/km^2), no ADRIA domain carries these -- see 07.
        corset_herbivore=parse(Float64, _cfg("CORSET_INITIAL_HERBIVORE_BIOMASS", "4.0e4")),
        corset_small_piscivore=parse(Float64, _cfg("CORSET_INITIAL_SMALL_PISCIVORE_BIOMASS", "4.0e3")),
        corset_large_piscivore=parse(Float64, _cfg("CORSET_INITIAL_LARGE_PISCIVORE_BIOMASS", "1.0e3")),
        corset_urchin=parse(Float64, _cfg("CORSET_INITIAL_URCHIN_BIOMASS", "1.0e4")),
        corset_use_cyclones=_cfg("CORSET_USE_CYCLONES_AS_HURRICANES", "false") == "true",
        corset_use_dhw=_cfg("CORSET_USE_DHW_AS_BLEACHING", "false") == "true",
        ## true: CoralPop and Kora take all larvae from the domain's connectivity matrix (a site's own
        ## through the diagonal), as CoralBlox and C~scape do, instead of their own retention/self-seeding.
        larvae_via_connectivity=_cfg("LARVAE_VIA_CONNECTIVITY", "false") == "true",
    )
end

"Settings recorded with every engine's output; `09` refuses to compare engines whose `FAIRNESS_KEYS` differ."
function settings_dict(cfg, years)::Dict{String,Any}
    return Dict{String,Any}(
        "DOMAIN_NAME" => basename(normpath(cfg.domain_path)), "DOMAIN_PATH" => cfg.domain_path,
        "RCP" => cfg.rcp, "DHW_IDX" => cfg.dhw_idx, "WAVE_IDX" => cfg.wave_idx,
        "CYCLONE_IDX" => cfg.cyclone_idx, "YEAR_START" => first(years), "YEAR_END" => last(years),
        "CALIBRATION_PATH" => cfg.calibration_path, "N_SAMPLES" => cfg.n_samples, "SEED" => cfg.seed,
        "N_WORKERS" => cfg.n_workers, "THREADS" => Threads.nthreads(),
        "CORSET_USE_DHW_AS_BLEACHING" => cfg.corset_use_dhw,
        "CORSET_USE_CYCLONES_AS_HURRICANES" => cfg.corset_use_cyclones,
        "LARVAE_VIA_CONNECTIVITY" => cfg.larvae_via_connectivity,
        "HOST" => gethostname(), "DATE" => string(now()),
    )
end
## The domain is compared by folder name, not full path, so results produced on the HPC can be
## re-analysed after copying them to another machine.
const FAIRNESS_KEYS = ("DOMAIN_NAME", "RCP", "DHW_IDX", "WAVE_IDX", "CYCLONE_IDX", "YEAR_START", "YEAR_END",
    "CALIBRATION_PATH")

# =============================================================================
# The six engines
# =============================================================================

## key => (label, colour, functional-group labels). Colours match 07 (Okabe-Ito).
const ENGINE_INFO = Dict(
    :coralblox => (label="CoralBlox", color="#0072B2", groups=string.(ADRIA.functional_group_names())),
    :cscape => (label="C~scape", color="#D55E00", groups=CscapeInterface.CSCAPE_FUNCTIONAL_TYPES),
    :coralpop => (label="CoralPop", color="#009E73", groups=CORALPOP_MORPHS),
    :kora => (label="Kora", color="#CC79A7", groups=KORA_GROUPS),
    :corset => (label="CORSET", color="#E69F00", groups=CORSET_GROUPS),
    :reefmod => (label="ReefMod", color="#56B4E9", groups=REEFMOD_GROUPS),
)
const ENGINE_ORDER = [:coralblox, :cscape, :coralpop, :kora, :corset, :reefmod]

"""
Which of ADRIA's 5 functional groups each engine-native group maps onto (indices into
`ADRIA.functional_group_names()`), reusing each interface's own translation table -- the same mapping
its domain translator uses for initial cover. CORSET's brooding/spawning split has no morphological
counterpart, so it is left out of the harmonised taxa figure.
"""
const ADRIA_GROUP_MAP = Dict(
    :coralblox => [[g] for g in 1:5],
    :kora => [[g] for g in KORA_ADRIA_GROUP],
    :cscape => [[g] for g in CscapeInterface.CSCAPE_FT_TO_ADRIA_GROUP],
    :coralpop => CORALPOP_MORPH_ADRIA_GROUPS,
    :reefmod => [[g] for g in REEFMOD_ADRIA_GROUPS],
)

engine_dir(project::AbstractString, key::Symbol) = mkpath(joinpath(model_run_dir(project), String(key)))

# =============================================================================
# Domain
# =============================================================================

"""
    load_shared_domain(cfg) -> (dom, dom_fixed, years)

The ADRIA domain every engine shares. `dom_fixed` has the DHW/wave/cyclone members fixed as factors
(how CoralBlox receives them); the other engines take the same indices through their translators.
"""
function load_shared_domain(cfg)
    dom = ADRIA.load_domain(cfg.domain_path, cfg.rcp; calib_params_fn=cfg.calibration_path)
    n_dhw = size(dom.dhw_scens, 3)
    0 <= cfg.dhw_idx <= n_dhw ||
        error("DHW_IDX=$(cfg.dhw_idx) is outside this domain's 1:$n_dhw DHW members (0 = no heat stress).")
    dom_fixed = deepcopy(dom)
    ADRIA.fix_factor!(dom_fixed, :dhw_scenario, cfg.dhw_idx)
    ADRIA.fix_factor!(dom_fixed, :wave_scenario, cfg.wave_idx)
    ADRIA.fix_factor!(dom_fixed, :cyclone_mortality_scenario, cfg.cyclone_idx)
    return dom, dom_fixed, collect(ADRIA.timesteps(dom))
end

# =============================================================================
# C~scape project folder (06 sections 1.1-1.3, converted demography only)
# =============================================================================

const CSCAPE_BASE_SCENARIO_ID = 1010001
const CSCAPE_DEMOG_STAMP = "ALLDAT_260903"
const CSCAPE_DEMOG_NAMES = Dict(
    "acro_table" => "acro_table", "acro_corym" => "acro_corym", "corym_non_acro" => "corym_non_acro",
    "brooder" => "corym_non_acro_BR", "small_massive" => "small_massive", "large_massive" => "large_massive",
)
const CSCAPE_FTS = ["acro_table", "acro_corym", "corym_non_acro", "small_massive", "large_massive", "brooder"]

"""
    setup_cscape_project(dir, cfg, years)

Make `dir` a C~scape project: the six converted `demog_julia_*` files in `dir/data` and a one-row
`ScenarioID.xlsx` -- the same row `06_RUNNING_A_NEW_PROJECT.jl` writes (Offshore_Northern,
`recruitment_multiplier = 10`). Idempotent.
"""
function setup_cscape_project(dir::AbstractString, cfg, years)
    data_dir = mkpath(joinpath(dir, "data"))
    input_names = Dict(ft => "demog_inputs_$(CSCAPE_DEMOG_NAMES[ft])$(CSCAPE_DEMOG_STAMP).rds" for ft in CSCAPE_FTS)
    for ft in CSCAPE_FTS
        f = replace(input_names[ft], "demog_inputs_" => "demog_julia_")
        src = joinpath(cfg.cscape_demog_dir, f)
        isfile(src) || error("Converted C~scape demography '$f' not found in CSCAPE_DEMOG_DIR=$(cfg.cscape_demog_dir).")
        dst = joinpath(data_dir, f)
        isfile(dst) || cp(src, dst)
    end
    xlsx = joinpath(dir, "ScenarioID.xlsx")
    isfile(xlsx) && return dir
    row = DataFrame(
        ID=[CSCAPE_BASE_SCENARIO_ID], num_repeats=[1], folder=["scenarios"], Purpose=[""],
        Intervention=["No"], date=[Dates.format(today(), "yyyy/mm/dd")], Region=["Offshore_Northern"],
        Year_Start=[first(years)], Year_End=[last(years)],
        init_cover=["0.1_0.1_0.1_0.1_0.1_0.1"],                     # placeholder: domain translation supplies it
        fts=[join(CSCAPE_FTS, "/")], HeatToleranceGroups=["-5_8_14"], HeatToleranceInit=["0_1.91"],
        Heritability=["0.5_0.01"], Plasticity=[0.0], counterfactual=[missing],
        Disturbance_file=["placeholder.RData"], Connectivity_file=["placeholder.RData"],
        Spatial_file=["placeholder.RData"], Geometry_file=["placeholder.RData"],
        Growth_Surv_file=[join([input_names[ft] for ft in CSCAPE_FTS], "/")],
        Fogging_reducer=[missing], Fogging_sites=[missing], Fogging_start=[missing], TradeOff=[missing],
        var"Deployment area"=[missing], TotalCorals=[missing], species=[missing], Enhancement=[missing],
        InterventionYears_start=[missing], duration=[missing], frequency=[missing], Reef_siteids=[missing],
        temp_growth=[1], species_proportions=[missing], recruitment_multiplier=[10.0],
    )
    XLSX.writetable(xlsx, row; overwrite=true, sheetname="ScenarioID")
    return dir
end

# =============================================================================
# Factor specs and runs
# =============================================================================

"""
    make_spec(key, dom_fixed, cfg, years, cscape_dir) -> DataFrame

Each engine's own ADRIA-sampleable factor spec, with the same bounds 06/07 use (C~scape narrowed to
+-10% to match CoralBlox's calibrated band).
"""
function make_spec(key::Symbol, dom_fixed, cfg, years, cscape_dir)::DataFrame
    rcp, ys, ye = parse(Int, cfg.rcp), first(years), last(years)
    key === :coralblox && return coralblox_spec(dom_fixed)
    key === :cscape && return cscape_spec(setup_cscape_project(cscape_dir, cfg, years), CSCAPE_BASE_SCENARIO_ID;
        mortality_bounds=(-0.1, 0.1), growth_bounds=(0.9, 1.1), recruitment_bounds=(0.9, 1.1))
    key === :coralpop && return coralpop_spec(; rcp, year_start=ys, year_end=ye)
    key === :kora && return kora_spec(; rcp, year_start=ys, year_end=ye)
    key === :corset && return corset_spec(; rcp, year_start=ys, year_end=ye)
    key === :reefmod && return reefmod_spec(; rcp, year_start=ys, year_end=ye)
    error("Unknown engine :$key -- expected one of $(ENGINE_ORDER).")
end

"""
    run_engine(key, dom, dom_fixed, scens, spec, out, cfg, years) -> result set

Run `scens` through engine `key`, writing the engine's own files under `out`. Raw outputs are kept
(`keep_raw=true`) where an engine supports it, because the size-structure extraction reads them.
"""
function run_engine(key::Symbol, dom, dom_fixed, scens::DataFrame, spec, out::AbstractString, cfg, years)
    mkpath(out)
    if key === :coralblox
        ## 06's order: setup() first (it returns early if ADRIA_OUTPUT_DIR is already set and resets
        ## ADRIA_LOG_COVER without a config.toml), THEN the overrides. The cover log is the per
        ## group/size cover the size-structure figures read.
        ADRIA.setup()
        ENV["ADRIA_OUTPUT_DIR"] = out
        ENV["ADRIA_LOG_COVER"] = "true"
        return ADRIA.run_scenarios(dom_fixed, scens, cfg.rcp)
    elseif key === :cscape
        setup_cscape_project(out, cfg, years)
        return ADRIA.run_scenarios(dom, scens, cfg.rcp, :cscape;
            cscape_fpath=out, cscape_base_scenario_id=CSCAPE_BASE_SCENARIO_ID,
            use_domain_translation=true, dhw_idx=cfg.dhw_idx, wave_idx=cfg.wave_idx, cyclone_idx=cfg.cyclone_idx,
            n_workers=max(1, cfg.n_workers), exeflags=["-t", string(cfg.cscape_threads)],
            export_adria=false, calc_indicators=true)
    elseif key === :coralpop
        return ADRIA.run_scenarios(dom, scens, cfg.rcp, :coralpop;
            coralpop_fpath=out, dhw_idx=cfg.dhw_idx, wave_idx=cfg.wave_idx, cyclone_idx=cfg.cyclone_idx,
            n_workers=cfg.n_workers, spec, keep_raw=true, larvae_via_connectivity=cfg.larvae_via_connectivity)
    elseif key === :kora
        return ADRIA.run_scenarios(dom, scens, cfg.rcp, :kora;
            kora_fpath=out, dhw_idx=cfg.dhw_idx, wave_idx=cfg.wave_idx, cyclone_idx=cfg.cyclone_idx,
            n_workers=cfg.n_workers, spec, keep_raw=true, larvae_via_connectivity=cfg.larvae_via_connectivity)
    elseif key === :corset
        return ADRIA.run_scenarios(dom, scens, cfg.rcp, :corset;
            corset_fpath=out, dhw_idx=cfg.dhw_idx, wave_idx=cfg.wave_idx, cyclone_idx=cfg.cyclone_idx,
            n_workers=cfg.n_workers, spec,
            translator_kwargs=(;
                initial_herbivore_biomass=cfg.corset_herbivore,
                initial_small_piscivore_biomass=cfg.corset_small_piscivore,
                initial_large_piscivore_biomass=cfg.corset_large_piscivore,
                initial_urchin_biomass=cfg.corset_urchin,
                use_cyclones_as_hurricanes=cfg.corset_use_cyclones,
                use_dhw_as_bleaching=cfg.corset_use_dhw))
    elseif key === :reefmod
        return ADRIA.run_scenarios(dom, scens, cfg.rcp, :reefmod;
            reefmod_fpath=out, dhw_idx=cfg.dhw_idx, n_workers=cfg.n_workers, spec, keep_raw=true)
    end
    error("Unknown engine :$key")
end

# =============================================================================
# Saving and loading
# =============================================================================

"Store the engine's result set in `dir` (CoralBlox: the path of ADRIA's own result store, relative to `dir`)."
function save_engine_result(key::Symbol, rs, dir::AbstractString)
    f = joinpath(dir, "result.jld2")
    if key === :coralblox
        write(joinpath(dir, "result_path.txt"), relpath(ADRIA.result_location(rs), dir))
    elseif key === :cscape
        CscapeInterface.save_results(rs, f)
    elseif key === :coralpop
        save_coralpop_results(rs, f)
    elseif key === :kora
        save_kora_results(rs, f)
    elseif key === :corset
        save_corset_results(rs, f)
    elseif key === :reefmod
        save_reefmod_results(rs, f)
    end
    return nothing
end

"Has `dir` a saved result for engine `key`?"
has_engine_result(key::Symbol, dir::AbstractString) =
    isfile(joinpath(dir, key === :coralblox ? "result_path.txt" : "result.jld2")) &&
    isfile(joinpath(dir, "samples.jld2"))

"""
    load_engine_result(key, dir, years) -> result set

Reload what `save_engine_result` stored. C~scape's result is re-aligned onto `years` (see
`align_cscape`) so the same ADRIA metric functions run on it as on the other five.
"""
function load_engine_result(key::Symbol, dir::AbstractString, years)
    f = joinpath(dir, "result.jld2")
    key === :coralblox && return ADRIA.load_results(normpath(joinpath(dir, strip(read(joinpath(dir, "result_path.txt"), String)))))
    key === :cscape && return align_cscape(CscapeInterface.load_results(CscapeInterface.CScapeResultSet, f), years)
    key === :coralpop && return load_coralpop_results(f)
    key === :kora && return load_kora_results(f)
    key === :corset && return load_corset_results(f)
    key === :reefmod && return load_reefmod_results(f)
    error("Unknown engine :$key")
end

"""
    align_cscape(rs, years)

C~scape stores its initial state under the year BEFORE it starts (2024) and simulates 2025..2099; the
other engines label their initial state 2025. Trim every outcome to the common `years` and rebuild it
as an ADRIA `DataCube` (Symbol-keyed properties, which ADRIA's Metric wrapper needs) -- lifted
unchanged from 07, where the full reasoning is.
"""
function align_cscape(rs, years)
    t_axis(x) = findfirst(==(:timesteps), ADRIA.axes_names(x))
    function trim_cube(x)
        nm = ADRIA.axes_names(x)
        labels = Any[ADRIA.axis_labels(x, n) for n in nm]
        ti = t_axis(x)
        ti === nothing && return x
        keep = findall(in(years), labels[ti])
        length(keep) == length(years) || error(
            "C~scape outcome covers years $(first(labels[ti]))-$(last(labels[ti])), " *
            "which does not include all of $(first(years))-$(last(years)).")
        idx = ntuple(d -> d == ti ? keep : Colon(), ndims(x))
        labels[ti] = labels[ti][keep]
        return ADRIA.DataCube(Array(x.data)[idx...]; NamedTuple{nm}(Tuple(labels))...)
    end
    el = rs.env_layer_md
    new_env = ADRIA.EnvLayer(String(el.dpkg_path), String(el.loc_data_fn), String(el.loc_id_col),
        String(el.cluster_id_col), String(el.init_coral_cov_fn), String(el.connectivity_fn),
        String(el.DHW_fn), String(el.wave_fn), first(years):last(years))
    outcomes = Dict{Symbol,Any}(k => trim_cube(v) for (k, v) in rs.outcomes)
    return CscapeInterface.CScapeResultSet(
        rs.name, rs.RCP, rs.loc_ids, rs.loc_area, rs.loc_max_coral_cover, rs.loc_centroids, new_env,
        rs.connectivity_data, rs.loc_data, rs.inputs, rs.sim_constants, rs.model_spec, outcomes,
        rs.coral_size_diameter, rs.raw_output_paths)
end

"A factor spec with `dist` stored by name (Distributions types don't round-trip through JLD2 reliably)."
spec_for_disk(spec::DataFrame) = DataFrame(
    fieldname=spec.fieldname, val=spec.val, dist=string.(spec.dist),
    dist_params=spec.dist_params, is_constant=spec.is_constant)

# =============================================================================
# Size structure: one common format for every engine that has one
# =============================================================================
#
# size_structure.jld2 holds
#   "years"  [T]
#   "groups" [G]           the engine's own group labels
#   "d_cm"   [G, B]        bin-centre colony diameter (cm)
#   "value"  [T, G, B, S]  domain-level amount per group/bin/scenario
#   "units"                what `value` measures -- NOT the same for every engine:
#                          CoralBlox/Kora/CoralPop: fraction of the domain's habitable area (cover);
#                          C~scape: IPM colony density (summed over sites); ReefMod: colony counts
#                          (summed over reefs). `09` converts counts to cover (count x bin-centre
#                          colony area) and compares size-distribution SHAPES (each class's share
#                          of the engine's cover), not absolute amounts, across engines.

"Domain total of `x[:, :, :, loc]`-style cover in m^2 as a fraction of the total habitable area."
_habitable_area(rs) = sum(ADRIA.loc_k_area(rs))

"""
    extract_size_structure(key, rs, run_dir, dom, cfg, scens) -> Dict or nothing

Read each scenario's raw output (or CoralBlox's cover log) into the common format above. `nothing`
for CORSET (no size structure). Engine-specific array layouts are documented inline; any failure is
reported and returns `nothing` rather than failing the run.
"""
function extract_size_structure(key::Symbol, rs, run_dir::AbstractString, dom, cfg, scens::DataFrame)
    try
        if key === :coralblox
            ## coral_cover_log [T, coral_id, sites, scenarios]: fraction of each site's habitable area,
            ## coral_id in ADRIA's coral-spec order (each group's 7 size classes contiguous -- checked:
            ## the domain's initial cover reshaped this way reproduces relative_taxa_cover exactly).
            ## ADRIA writes the log inside its time loop AFTER growth and BEFORE that year's bleaching/
            ## cyclone mortality, from the 2nd timestep on; the 1st (never written, all zero) is the
            ## initial state, filled here from the domain's initial cover (same units, same order).
            cl = rs.coral_cover_log
            (cl === nothing || size(cl, 3) == 1) && return nothing   # logging was off
            a = Float64.(Array(cl))
            T, C, L, S = size(a)
            a[1, :, :, :] .= reshape(Float64.(Array(dom.init_coral_cover)), C, L, 1)    # same for every scenario
            n_sizes = size(ADRIA.bin_edges(), 2) - 1
            G = C ÷ n_sizes
            w = ADRIA.loc_k_area(rs)
            dom_cover = dropdims(sum(a .* reshape(w, 1, 1, L, 1); dims=3); dims=3) ./ sum(w)   # [T, C, S]
            value = permutedims(reshape(dom_cover, T, n_sizes, G, S), (1, 3, 2, 4))         # [T, G, B, S]
            edges = ADRIA.bin_edges(; unit=:cm)
            d_cm = (edges[:, 1:(end - 1)] .+ edges[:, 2:end]) ./ 2
            return Dict("years" => collect(ADRIA.timesteps(rs)), "groups" => ENGINE_INFO[key].groups,
                "d_cm" => d_cm, "value" => value,
                "units" => "fraction of domain habitable area (cover, after growth / before each year's disturbance)")
        elseif key === :kora
            ## KoraOutput.cover_m2 [T, G, B, L] (m^2 of colony surface area).
            outs = [load_kora_output(p) for p in rs.raw_output_paths]
            o1 = outs[1]
            hab = sum(o1.reef_area_m2 .* o1.k)
            value = cat([dropdims(sum(o.cover_m2; dims=4); dims=4) ./ hab for o in outs]...; dims=4)
            return Dict("years" => o1.years, "groups" => ENGINE_INFO[key].groups,
                "d_cm" => kora_bin_diameters_cm(), "value" => value,
                "units" => "fraction of domain habitable area (cover)")
        elseif key === :coralpop
            ## CoralPopOutput.area_m2 [T1, n_bins, 3 morphs, L] (m^2), T1 includes the initial year.
            outs = [load_coralpop_output(p) for p in rs.raw_output_paths]
            o1 = outs[1]
            hab = sum(o1.reef_area_m2 .* o1.k_pct ./ 100)
            value = cat([permutedims(dropdims(sum(o.area_m2; dims=4); dims=4), (1, 3, 2)) ./ hab for o in outs]...; dims=4)
            years = length(o1.years) == size(value, 1) ? o1.years : collect(first(o1.years) .+ (0:(size(value, 1) - 1)))
            return Dict("years" => years, "groups" => ENGINE_INFO[key].groups, "d_cm" => Float64.(rs.coral_size_diameter),
                "value" => value, "units" => "fraction of domain habitable area (cover)")
        elseif key === :cscape
            return _cscape_size_structure(rs, run_dir, dom, cfg, scens)
        elseif key === :reefmod
            return _reefmod_size_structure(rs)
        end
    catch err
        @warn "Could not extract $(ENGINE_INFO[key].label)'s size structure; the size figures will skip it" exception = (err, catch_backtrace())
    end
    return nothing
end

## C~scape: out_array (years, sites, intervention, ft, heat_bin, metric); metric VAL_SIZE = the 100 IPM
## size classes (colony density), intervention slot 3 = combined. Loaded exactly as 06 section 5.4 does.
function _cscape_size_structure(rs, run_dir, dom, cfg, scens)
    files = sort(filter(f -> occursin(r"^Array_scenario_.*\.jld2$", f), readdir(joinpath(run_dir, "model_outputs"))))
    isempty(files) && return nothing
    cs_scen, _, _ = CscapeInterface.cscape_sample_to_scenarios(run_dir, CSCAPE_BASE_SCENARIO_ID, scens)
    cs_scen_jl = CscapeInterface._demog_julia_scenarios(run_dir; scenarios=cs_scen)
    inp = CscapeInterface.domain_to_cscape_inputs(dom; dhw_idx=cfg.dhw_idx, wave_idx=cfg.wave_idx, cyclone_idx=cfg.cyclone_idx)
    vals, d_cm, years = Array{Float64,3}[], nothing, nothing
    for f in files
        m = match(r"Array_scenario_(\d+)_rep_(\d+)_draw_(.+)\.jld2$", f)
        out = CScape.load_output(run_dir, parse(Int, m[1]); rep_index=parse(Int, m[2]), draw=parse(Int, m[3]),
            scenarios=cs_scen_jl, spatial_override=inp.spatial_override)
        oa = out.out_array
        push!(vals, dropdims(sum(oa[:, :, 3, :, :, CScape.Layout.VAL_SIZE]; dims=(2, 4)); dims=(2, 4)))   # [T, ft, 100]
        d_cm = out.meshpoints
        years = collect(ADRIA.timesteps(rs))[1:min(end, size(oa, 1))]
    end
    value = cat(vals...; dims=4)                                                                     # [T, G, B, S]
    ## C~scape's array includes its extra initial year (see align_cscape): keep the last length(years).
    value = value[(end - length(years) + 1):end, :, :, :]
    d = size(d_cm, 1) == size(value, 2) ? d_cm : repeat(reshape(vec(d_cm[1, :]), 1, :), size(value, 2), 1)
    return Dict("years" => years, "groups" => ENGINE_INFO[:cscape].groups, "d_cm" => Float64.(d),
        "value" => value, "units" => "IPM colony density (summed over sites)")
end

## ReefMod: RESULT.coral_{juv,adol,adult}_count [R, T6+1, 6, bins] (colony counts per reef grid). Bin edges
## from ReefMod-GBR/settings/PARAMETERS.m ("UPDATE March 2022"): juveniles 1-3-5 cm, adolescents
## 5-9-13-17 cm, adults 17-31-45-59-73-87-101 cm. Annual values at the same 6-monthly indices
## ReefModInterface's own indicators use.
function _reefmod_size_structure(rs)
    edges = [1, 3, 5, 9, 13, 17, 31, 45, 59, 73, 87, 101]
    mids = (edges[1:(end - 1)] .+ edges[2:end]) ./ 2
    vals = Array{Float64,3}[]
    for p in rs.raw_output_paths
        R = MAT.matread(p)["RESULT"]
        haskey(R, "coral_adol_count") || return nothing
        counts = cat(Float64.(R["coral_juv_count"]), Float64.(R["coral_adol_count"]), Float64.(R["coral_adult_count"]); dims=4)
        idx = reefmod_annual_indices(size(counts, 2) - 1)[1:(end - 1)]
        push!(vals, dropdims(sum(counts[:, idx, :, :]; dims=1); dims=1))                       # [T, 6, 11]
    end
    value = cat(vals...; dims=4)
    return Dict("years" => collect(ADRIA.timesteps(rs))[1:size(value, 1)], "groups" => ENGINE_INFO[:reefmod].groups,
        "d_cm" => repeat(reshape(mids, 1, :), 6, 1), "value" => value, "units" => "colony counts (summed over reefs)")
end

# =============================================================================
# Benchmarking
# =============================================================================

"""
    timed_call(f) -> (value, NamedTuple)

`f()` under `@timed`: wall time, GC time, allocations, and the share of time spent compiling (Julia
1.12's `compile_time`). Peak resident memory of THIS process is `Sys.maxrss()`; worker processes
(C~scape/CoralPop Distributed workers, ReefMod's MATLAB subprocesses) are not included in it -- the
Slurm accounting (`sacct MaxRSS`) collected by `hpc/jobs/collect_benchmarks.sh` is.
"""
function timed_call(f)
    t = @timed f()
    return t.value, (time_s=t.time, gc_s=t.gctime, alloc_gb=t.bytes / 1e9,
        compile_s=hasproperty(t, :compile_time) ? Float64(t.compile_time) : NaN)
end

"Total size of every file under `dir`, GB."
function dir_size_gb(dir::AbstractString)::Float64
    isdir(dir) || return 0.0
    total = 0
    for (root, _, files) in walkdir(dir), f in files
        total += filesize(joinpath(root, f))
    end
    return total / 1e9
end

"Write `rows` (a vector of NamedTuples with the same fields) as CSV."
function write_csv(path::AbstractString, rows::AbstractVector)
    isempty(rows) && return nothing
    cols = keys(first(rows))
    open(path, "w") do io
        println(io, join(cols, ","))
        for r in rows
            println(io, join((v isa AbstractFloat ? @sprintf("%.6g", v) : replace(string(v), "," => ";") for v in values(r)), ","))
        end
    end
    return nothing
end

"Read a CSV written by `write_csv` (or any simple comma-separated file) into a DataFrame of strings/numbers."
function read_csv(path::AbstractString)::DataFrame
    lines = filter(!isempty, readlines(path))
    hdr = split(lines[1], ",")
    rows = [split(l, ","; limit=length(hdr)) for l in lines[2:end]]
    parsev(s) = (x = tryparse(Float64, s); x === nothing ? String(s) : x)
    return DataFrame([Symbol(h) => [parsev(r[i]) for r in rows] for (i, h) in enumerate(hdr)])
end
