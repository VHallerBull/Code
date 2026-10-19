#=
07_COMPARE_SIX_ENGINES.jl

CoralBlox vs C~scape vs CoralPop vs Kora vs ReefMod vs CORSET on ONE ADRIA domain and ONE shared
climate/cyclone member, each sampled with ADRIA's own sampler from its own factor spec, each returning
an ADRIA `ResultSet`, and each analysed with the SAME ADRIA metric functions and reef indicators.
Writes figures and a summary table to `../ModelRuns/SixEngineComparison/`.

    julia +1.12 --project=. 07_COMPARE_SIX_ENGINES.jl        # from this folder -- Kora needs Julia 1.12+

How the six engines are obtained
- CoralBlox: run here (fast) with `ADRIA.run_scenarios(dom, scens, RCP)`, or reloaded (`CORALBLOX_MODE=:load`).
- CoralPop:  run here through ADRIA's `:coralpop` engine (needs CoralPopInterface + the Python
             `coralpop` engine, see the Python setup below), or reloaded (`CORALPOP_MODE=:load`).
- Kora:      run here through ADRIA's `:kora` engine (needs KoraInterface; pure Julia, no external
             interpreter), or reloaded (`KORA_MODE=:load`).
- CORSET:    run here through ADRIA's `:corset` engine (needs CORSETInterface; pure Julia). CORSET is a
             from-scratch reimplementation -- its own source was never published, see
             `CORSET.jl/CAVEATS.md` -- with no coral size structure and no fish/urchin/benthic state of
             its own, so an ADRIA domain supplies far less of what it needs than for the other three.
             The `CORSET_INITIAL_*_BIOMASS` constants below are PLACEHOLDER consumer biomass, not survey
             data for this domain -- see `CORSETInterface.jl/src/domain_translator.jl` before trusting
             CORSET's absolute trajectories here, as opposed to its relative response to disturbance.
             Reloaded with `CORSET_MODE=:load`.
- ReefMod:   run here through ADRIA's `:reefmod` engine (needs ReefModInterface AND a real MATLAB on
             `PATH` -- see the setup below), or reloaded (`REEFMOD_MODE=:load`). Unlike
             CoralPop/Kora/CORSET, each scenario is a `matlab -batch ...` SUBPROCESS, not an
             in-process call -- see `ReefModInterface.jl/src/reefmod_bridge.jl`. No local MATLAB
             license? `ReefModInterface.jl/jobs/submit_reefmod_pipeline.sh` runs the same thing on an
             HPC's Slurm-scheduled license instead; copy its `reefmod_result.jld2` into
             `ModelRuns/<PROJECT_NAME>/` and use `REEFMOD_MODE=:load` here either way. The eventual
             target is all six engines (this one included) inside one Apptainer image on the HPC, at
             which point `:run` here and the HPC pipeline converge to the same MATLAB.
- C~scape:   NOT run here -- it needs the R-fitted demography, ~4 GB per scenario and a long run. It is
             loaded from the result set that `06_RUNNING_A_NEW_PROJECT.jl` saves (`my_run.jld2` +
             `samples.jld2`) in its own `ModelRuns/<CSCAPE_PROJECT_NAME>/` folder. If that is missing, the
             comparison continues with the other engines and says so.

Fairness checks made below (the script stops rather than compare unlike things): all engines must use
the same locations, the same years, and -- for C~scape, whose settings were recorded when its result
was produced -- the same domain, RCP, DHW/wave/cyclone member and disturbance switch.

Time axis: C~scape stores its initial state under the year BEFORE it starts (2024) and simulates 2025-2099;
CoralBlox, CoralPop, Kora, CORSET and ReefMod label their initial state 2025. The comparison uses the
common years (2025-2099), so C~scape's 2025 state has already had one year of dynamics while the others'
has not -- a one-year offset in where each engine starts, not in which year's climate it sees.

What is compared: domain-mean (k-area weighted) relative cover, juveniles, shelter volume, juvenile
indicator, Reef Condition Index and Reef Tourism Index per engine, taxa cover per engine (CoralBlox, Kora
and C~scape share ADRIA's 5 functional groups; CoralPop and ReefMod keep their own 3/6 groups, CORSET its
own 2 -- brooding/spawning), and PAWN sensitivity of final-year cover to each engine's own sampled
factors. The engines' factor spaces are different by design, so PAWN rankings are compared per engine,
not factor by factor.

Kora has no cyclone-mortality pathway (a documented limitation of the published model, not a bug here --
see `KoraInterface.jl/src/domain_translator.jl`), so its trajectories reflect DHW-only disturbance even
when `CYCLONE_IDX` is set for the other engines. CORSET by default ALSO ignores `DHW_IDX`/`CYCLONE_IDX`
(`CORSET_USE_CYCLONES_AS_HURRICANES`/`CORSET_USE_DHW_AS_BLEACHING` default `false`) -- turning them on
substitutes ADRIA's cyclone/DHW layers for CORSET's own hurricane/bleaching pulses, a modelling judgement
call documented in `domain_translator.jl`'s header, not a validated equivalence between the two models'
disturbance definitions. ReefMod has a real cyclone (and CoTS) pathway but it is switched OFF for this
comparison (not wired to ADRIA's own cyclone forcing in this first pass -- see
`ReefModInterface.jl/src/domain_translator.jl`), so its trajectories are also DHW-only disturbance.
ReefMod's `:relative_shelter_volume` (and the Reef Condition/Tourism indices derived from it) are not yet
computed (`:absolute_shelter_volume` only, see `ReefModInterface.jl/src/indicators.jl`); CORSET computes
NONE of relative_juveniles/relative_shelter_volume/juvenile_indicator/reef_condition_index/
reef_tourism_index at all, because none are defined without coral size structure (see
`CORSETInterface.jl/src/indicators.jl`) -- the shelter-volume validity check in section 8 below omits
both from those panels the same way it already handles any engine outside ADRIA's `[0, 1]` convention.

CAUTION: PAWN needs N/S > 80 (see `ADRIAanalysis.sensitivity.pawn`), i.e. N_SAMPLES > 800 at S=10. At the
default N_SAMPLES the PAWN panels are illustrative only.

Environment variables `CMR_<NAME>` override the CONFIG constants below (handy for a quick test run).
=#

include(joinpath(@__DIR__, "paths.jl"))       # CODE_DIR, MODEL_RUNS_DIR, model_run_dir(name)

# =============================================================================
# CONFIG -- edit these
# =============================================================================

_cfg(name, default) = get(ENV, "CMR_" * name, default)

## Python for CoralPop. PythonCall reads these when it loads, so they MUST be set before the `using`
## block below (and Distributed workers inherit them). A conda environment with numpy + pandas:
ENV["JULIA_CONDAPKG_BACKEND"] = "Null"
haskey(ENV, "JULIA_PYTHONCALL_EXE") ||
    (ENV["JULIA_PYTHONCALL_EXE"] = raw"C:\Users\vhaller\.julia\conda\3\x86_64\envs\coral\python.exe")
haskey(ENV, "CORALPOP_ENGINE_DIR") ||
    (ENV["CORALPOP_ENGINE_DIR"] = raw"C:\Users\vhaller\Documents\GitHub\CoralPop-Engine")

## MATLAB for ReefMod. Read by ReefModInterface.run_reefmod_scenarios (not at `using` time -- each
## scenario is a fresh subprocess -- but set here anyway, alongside the Python setup above, so every
## engine's external-tool config lives in one place). Defaults: `matlab` on PATH; the sibling
## `ReefMod-GBR` clone.
haskey(ENV, "REEFMOD_MATLAB_EXE") || (ENV["REEFMOD_MATLAB_EXE"] = "matlab")
haskey(ENV, "REEFMOD_GBR_DIR") ||
    (ENV["REEFMOD_GBR_DIR"] = raw"C:\Users\vhaller\Documents\GitHub\ReefMod-GBR")

const PROJECT_NAME = _cfg("PROJECT_NAME", "SixEngineComparison")           # output folder under ModelRuns/
const CSCAPE_PROJECT_NAME = _cfg("CSCAPE_PROJECT_NAME", "AdjustedRecruitment10_Disturbance")   # 06's project

## An ADRIA domain package on disk, and the RCP -- keep identical to 06 (checked below for C~scape).
const DOMAIN_PATH = _cfg("DOMAIN_PATH", raw"C:\Users\vhaller\Downloads\Moore_2025-08-22_v080")
const RCP = _cfg("RCP", "45")

## CoralBlox calibration (`calibrated_params.nc`) or `nothing` for ADRIA's defaults -- keep identical to 06.
const CORALBLOX_CALIBRATION_OVERRIDES_PATH = _cfg("CALIBRATION_PATH", "C:\\Users\\vhaller\\Australian Institute of Marine Science\\AIMS Decision Support - General\\01_PROJECTS\\RRAP M&DS\\ADRIA\\ADRIA Domain data\\CoralBlox-ADRIA GBR Calibrated Parameters\\2026-08-11\\calibrated_params.nc")

## Turn off DHW/wave/cyclone forcing (intrinsic dynamics only) -- keep identical to 06.
const NO_DISTURBANCE = _cfg("NO_DISTURBANCE", "false") == "true"
## The single climate member all engines share (see 06 for why 64: draw 1 is a near-worst-case DHW draw).
## WAVE_IDX must stay 1 for the Moore domain (its wave_scens has one all-zero member); CoralPop ignores it.
## `CMR_DHW_IDX`/`CMR_WAVE_IDX`/`CMR_CYCLONE_IDX` override the members (e.g. 40/16/1, the ADRIA Test
## domain's median-severity DHW and wave members -- see engines.jl).
const DHW_IDX, WAVE_IDX, CYCLONE_IDX = NO_DISTURBANCE ? (0, 0, 0) :
    parse.(Int, (_cfg("DHW_IDX", "64"), _cfg("WAVE_IDX", "1"), _cfg("CYCLONE_IDX", "1")))

## ADRIA.sample's default Sobol' sampler requires a power of 2.
const N_SAMPLES = parse(Int, _cfg("N_SAMPLES", "8"))

## :run = run it here; :load = reload what an earlier run of THIS script saved. C~scape: :load (from
## 06) or :skip -- it can never be :run from this script, see the header.
const CORALBLOX_MODE = Symbol(_cfg("CORALBLOX_MODE", "run"))
const CORALPOP_MODE = Symbol(_cfg("CORALPOP_MODE", "run"))
const KORA_MODE = Symbol(_cfg("KORA_MODE", "run"))
const CORSET_MODE = Symbol(_cfg("CORSET_MODE", "run"))
const CSCAPE_MODE = Symbol(_cfg("CSCAPE_MODE", "load"))
## ReefMod defaults to :run (needs a real `matlab` on PATH -- workstation license, Apptainer with the
## host's MATLAB bind-mounted, or an HPC `module load matlab`). Set CMR_REEFMOD_MODE=load to reuse an
## earlier run's saved result instead (e.g. iterating on the other five engines or the analysis step
## without re-paying ReefMod's MATLAB cost -- see jobs.txt/REEFMOD_MODE notes in ReefModInterface.jl's
## README for why that's sometimes deliberate), or =skip to leave ReefMod out entirely.
const REEFMOD_MODE = Symbol(_cfg("REEFMOD_MODE", "run"))

## CoralPop parallelism: Distributed workers, each with its own Python interpreter (0 = this process).
const CORALPOP_N_WORKERS = parse(Int, _cfg("CORALPOP_N_WORKERS", "2"))
## Kora/CORSET parallelism: a `Threads.@threads` cap, not separate processes (pure Julia) -- start this
## script with `julia -t N` for a value above 1 to do anything.
const KORA_N_WORKERS = parse(Int, _cfg("KORA_N_WORKERS", "1"))
const CORSET_N_WORKERS = parse(Int, _cfg("CORSET_N_WORKERS", "1"))
## ReefMod parallelism: concurrent `matlab` subprocesses (each its own OS process and MATLAB license
## checkout -- like CoralPop's Distributed workers, not Kora/CORSET's thread cap). Keep this within
## however many MATLAB licenses/cores are actually available wherever this runs.
const REEFMOD_N_WORKERS = parse(Int, _cfg("REEFMOD_N_WORKERS", "1"))

## CORSET needs consumer (fish/urchin) biomass that no ADRIA domain carries at all -- see
## `CORSETInterface.jl/src/domain_translator.jl`. These are PLACEHOLDER values (kg/km^2) for a
## structural comparison, NOT survey data for the Moore domain or any other -- replace them with real
## site data before drawing any conclusion from CORSET's absolute cover/biomass here.
const CORSET_INITIAL_HERBIVORE_BIOMASS = parse(Float64, _cfg("CORSET_INITIAL_HERBIVORE_BIOMASS", "4.0e4"))
const CORSET_INITIAL_SMALL_PISCIVORE_BIOMASS = parse(Float64, _cfg("CORSET_INITIAL_SMALL_PISCIVORE_BIOMASS", "4.0e3"))
const CORSET_INITIAL_LARGE_PISCIVORE_BIOMASS = parse(Float64, _cfg("CORSET_INITIAL_LARGE_PISCIVORE_BIOMASS", "1.0e3"))
const CORSET_INITIAL_URCHIN_BIOMASS = parse(Float64, _cfg("CORSET_INITIAL_URCHIN_BIOMASS", "1.0e4"))
## Off by default -- see the header's note on substituting ADRIA's cyclone/DHW layers for CORSET's own
## hurricane/bleaching pulses.
const CORSET_USE_CYCLONES_AS_HURRICANES = _cfg("CORSET_USE_CYCLONES_AS_HURRICANES", "false") == "true"
const CORSET_USE_DHW_AS_BLEACHING = _cfg("CORSET_USE_DHW_AS_BLEACHING", "false") == "true"

for (n, m, ok) in (("CORALBLOX_MODE", CORALBLOX_MODE, (:run, :load)), ("CORALPOP_MODE", CORALPOP_MODE, (:run, :load)),
    ("KORA_MODE", KORA_MODE, (:run, :load)), ("CORSET_MODE", CORSET_MODE, (:run, :load)),
    ("CSCAPE_MODE", CSCAPE_MODE, (:load, :skip)), ("REEFMOD_MODE", REEFMOD_MODE, (:run, :load, :skip)))
    m in ok || error("$n must be one of $ok, got :$m")
end

# =============================================================================
# Packages
# =============================================================================

using ADRIA, CscapeInterface, CoralPopInterface, KoraInterface, CORSETInterface, ReefModInterface,
      DataFrames, Statistics, Printf
using ADRIAanalysis
## ADRIAviz's Makie extension needs ADRIAviz + Makie/GeoMakie/GraphMakie loaded; CairoMakie (a pure
## software rasteriser) rather than WGLMakie because this is a headless script writing PNGs. `import`, not
## `using`, so Makie's exports don't clash with other names.
import ADRIAviz, CairoMakie, GeoMakie, GraphMakie
CairoMakie.activate!()
import JLD2, YAXArrays

const M = ADRIA.metrics

@info "Python used for CoralPop" coralpop_python_info()

# =============================================================================
# 0. Project folders
# =============================================================================

const FPATH = model_run_dir(PROJECT_NAME)
const IMAGES_DIR = joinpath(FPATH, "images")
const SAMPLES_FILE = joinpath(FPATH, "samples.jld2")
mkpath(IMAGES_DIR)
@info "Project folder ready" FPATH

# =============================================================================
# 1. The shared ADRIA domain
# =============================================================================

dom = ADRIA.load_domain(
    DOMAIN_PATH, RCP; calib_params_fn=something(CORALBLOX_CALIBRATION_OVERRIDES_PATH, "")
)
YEARS = collect(ADRIA.timesteps(dom))

## CoralBlox reads DHW/wave/cyclone through fixed factors; CoralPop and C~scape take the same member
## indices through their domain translators. One trajectory for everyone, not each drawing its own.
dom_nodist = deepcopy(dom)
ADRIA.fix_factor!(dom_nodist, :dhw_scenario, DHW_IDX)
ADRIA.fix_factor!(dom_nodist, :wave_scenario, WAVE_IDX)
ADRIA.fix_factor!(dom_nodist, :cyclone_mortality_scenario, CYCLONE_IDX)

## Everything each engine contributes to the comparison
struct EngineRun
    key::Symbol
    label::String
    color::String
    rs::Any                     # an ADRIA.ResultSet
    factors::DataFrame          # the sampled factor table, one row per scenario, VARYING factors only
    groups::Vector{String}      # the engine's own functional-group labels
end
engines = EngineRun[]

# =============================================================================
# 2. CoralBlox
# =============================================================================

const CB_RESULT_PATH_FILE = joinpath(FPATH, "coralblox_result_path.txt")

ADRIA.setup()       # first: it returns early (leaving other defaults unset) if ADRIA_OUTPUT_DIR is already set
ENV["ADRIA_OUTPUT_DIR"] = joinpath(FPATH, "coralblox")
mkpath(ENV["ADRIA_OUTPUT_DIR"])

if CORALBLOX_MODE === :run
    cb_spec = coralblox_spec(dom_nodist)
    scens_cb = ADRIA.sample(cb_spec, N_SAMPLES)
    rs_cb = ADRIA.run_scenarios(dom_nodist, scens_cb, RCP)
    cb_varying = cb_spec.fieldname[.!cb_spec.is_constant]
    write(CB_RESULT_PATH_FILE, ADRIA.result_location(rs_cb))
    JLD2.save(SAMPLES_FILE, merge(isfile(SAMPLES_FILE) ? JLD2.load(SAMPLES_FILE) : Dict{String,Any}(),
        Dict("scens_cb" => scens_cb, "cb_varying" => String.(cb_varying))))
else
    isfile(CB_RESULT_PATH_FILE) || error("CORALBLOX_MODE=:load but no earlier run recorded in $CB_RESULT_PATH_FILE")
    rs_cb = ADRIA.load_results(String(strip(read(CB_RESULT_PATH_FILE, String))))
    s = JLD2.load(SAMPLES_FILE)
    scens_cb, cb_varying = s["scens_cb"], Symbol.(s["cb_varying"])
end
push!(engines, EngineRun(:coralblox, "CoralBlox", "#0072B2", rs_cb,
    DataFrames.select(scens_cb, cb_varying), string.(ADRIA.functional_group_names())))

# =============================================================================
# 3. CoralPop -- through ADRIA's own entry point
# =============================================================================
# `:coralpop` is provided by the ADRIACoralPopExt extension (active because CoralPopInterface is loaded).

const CP_RESULT_FILE = joinpath(FPATH, "coralpop_result.jld2")

if CORALPOP_MODE === :run
    cp_spec = coralpop_spec(; rcp=parse(Int, RCP), year_start=first(YEARS), year_end=last(YEARS))
    scens_cp = ADRIA.sample(cp_spec, N_SAMPLES)
    cp_dir = joinpath(FPATH, "coralpop")
    ## A fresh directory per run is tidier, but not needed for correctness: only this run's scenario IDs are
    ## loaded back and same-ID files are overwritten. OneDrive can hold freshly written files locked while it
    ## syncs them, which makes the delete fail, so don't stop on that.
    try
        rm(cp_dir; force=true, recursive=true)
    catch err
        @warn "Could not clear the previous CoralPop output folder (OneDrive may be syncing it); continuing" cp_dir exception = err
    end
    rs_cp = ADRIA.run_scenarios(
        dom, scens_cp, RCP, :coralpop;
        coralpop_fpath=cp_dir, dhw_idx=DHW_IDX, wave_idx=WAVE_IDX, cyclone_idx=CYCLONE_IDX,
        n_workers=CORALPOP_N_WORKERS, spec=cp_spec
    )
    save_coralpop_results(rs_cp, CP_RESULT_FILE)
    cp_varying = cp_spec.fieldname[.!cp_spec.is_constant]
    JLD2.save(SAMPLES_FILE, merge(JLD2.load(SAMPLES_FILE),
        Dict("scens_cp" => scens_cp, "cp_varying" => String.(cp_varying))))
else
    isfile(CP_RESULT_FILE) || error("CORALPOP_MODE=:load but no earlier run saved at $CP_RESULT_FILE")
    rs_cp = load_coralpop_results(CP_RESULT_FILE)
    cp_varying = Symbol.(JLD2.load(SAMPLES_FILE)["cp_varying"])
end
push!(engines, EngineRun(:coralpop, "CoralPop", "#009E73", rs_cp,
    DataFrames.select(rs_cp.inputs, cp_varying), CORALPOP_MORPHS))

# =============================================================================
# 4. Kora -- through ADRIA's own entry point
# =============================================================================
# `:kora` is provided by the ADRIAKoraExt extension (active because KoraInterface is loaded).
# Kora is pure Julia (no Python bridge) and its 5 groups/7 size classes ARE ADRIA's own, so this
# section needs no group-aggregation or size-rebinning caveats, unlike CoralPop's.

const KORA_RESULT_FILE = joinpath(FPATH, "kora_result.jld2")

if KORA_MODE === :run
    kora_spec_ = kora_spec(; rcp=parse(Int, RCP), year_start=first(YEARS), year_end=last(YEARS))
    scens_kora = ADRIA.sample(kora_spec_, N_SAMPLES)
    kora_dir = joinpath(FPATH, "kora")
    try
        rm(kora_dir; force=true, recursive=true)
    catch err
        @warn "Could not clear the previous Kora output folder (OneDrive may be syncing it); continuing" kora_dir exception = err
    end
    rs_kora = ADRIA.run_scenarios(
        dom, scens_kora, RCP, :kora;
        kora_fpath=kora_dir, dhw_idx=DHW_IDX, wave_idx=WAVE_IDX, cyclone_idx=CYCLONE_IDX,
        n_workers=KORA_N_WORKERS, spec=kora_spec_
    )
    save_kora_results(rs_kora, KORA_RESULT_FILE)
    kora_varying = kora_spec_.fieldname[.!kora_spec_.is_constant]
    JLD2.save(SAMPLES_FILE, merge(JLD2.load(SAMPLES_FILE),
        Dict("scens_kora" => scens_kora, "kora_varying" => String.(kora_varying))))
else
    isfile(KORA_RESULT_FILE) || error("KORA_MODE=:load but no earlier run saved at $KORA_RESULT_FILE")
    rs_kora = load_kora_results(KORA_RESULT_FILE)
    kora_varying = Symbol.(JLD2.load(SAMPLES_FILE)["kora_varying"])
end
push!(engines, EngineRun(:kora, "Kora", "#CC79A7", rs_kora,
    DataFrames.select(rs_kora.inputs, kora_varying), KORA_GROUPS))

# =============================================================================
# 5. C~scape -- loaded from 06's outputs
# =============================================================================

const CS_DIR = joinpath(MODEL_RUNS_DIR, CSCAPE_PROJECT_NAME)
const CS_RESULT_FILE = joinpath(CS_DIR, "my_run.jld2")
const CS_SAMPLES_FILE = joinpath(CS_DIR, "samples.jld2")

if CSCAPE_MODE === :load
    if isfile(CS_RESULT_FILE) && isfile(CS_SAMPLES_FILE)
        cs_samples = JLD2.load(CS_SAMPLES_FILE)

        ## The comparison is only meaningful if 06 used the same domain and climate member as this script.
        cfg = Dict{String,Any}(
            "DOMAIN_PATH" => DOMAIN_PATH, "RCP" => RCP, "DHW_IDX" => DHW_IDX, "WAVE_IDX" => WAVE_IDX,
            "CYCLONE_IDX" => CYCLONE_IDX, "NO_DISTURBANCE" => NO_DISTURBANCE,
            "YEAR_START" => first(YEARS), "YEAR_END" => last(YEARS)
        )
        mismatch = [k for (k, v) in cfg if get(cs_samples["settings"], k, missing) != v]
        isempty(mismatch) || error(
            "06's C~scape run (", CS_DIR, ") used different settings than this script for: ",
            join(["$k (06: $(get(cs_samples["settings"], k, missing)), 07: $(cfg[k]))" for k in mismatch], "; "),
            ". Re-run 06 with matching settings or change this script's CONFIG."
        )

        rs_cs_raw = CscapeInterface.load_results(CscapeInterface.CScapeResultSet, CS_RESULT_FILE)

        ## Year alignment. C~scape stores its INITIAL state under the year before it starts (2024) and
        ## simulates 2025..2099; CoralBlox and CoralPop label their initial state 2025 and simulate 2026..2099.
        ## All three have then applied each year's forcing to the state carrying that year's label, so the
        ## engines are compared on the common years (dropping C~scape's extra initial state). One consequence
        ## to keep in mind: C~scape's 2025 state has already had one year of dynamics, the others' has not.
        ##
        ## ADRIA's Metric wrapper also annotates results with Symbol-keyed properties, which the plain
        ## `YAXArray`s in a CScapeResultSet do not support (Dict{String,Any}). Each outcome is therefore
        ## rebuilt as an ADRIA `DataCube` (same axes and data, trimmed to the common years) so the SAME
        ## wrapped metric functions run on all three engines.
        function align_cscape(rs, years)
            t_axis(x) = findfirst(==(:timesteps), ADRIA.axes_names(x))
            function trim_cube(x)
                nm = ADRIA.axes_names(x)
                labels = Any[ADRIA.axis_labels(x, n) for n in nm]
                ti = t_axis(x)
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
        rs_cs = align_cscape(rs_cs_raw, YEARS)

        cs_spec_disk = cs_samples["cs_spec"]
        cs_varying = cs_spec_disk.fieldname[.!cs_spec_disk.is_constant]
        push!(engines, EngineRun(:cscape, "C~scape", "#D55E00", rs_cs,
            DataFrames.select(cs_samples["sampled_cs"], cs_varying), CscapeInterface.CSCAPE_FUNCTIONAL_TYPES))
    else
        @warn "C~scape results not found -- continuing without it. " *
              "Run 06_RUNNING_A_NEW_PROJECT.jl first (it writes my_run.jld2 and samples.jld2)." CS_DIR
    end
end

# =============================================================================
# 6. CORSET -- through ADRIA's own entry point
# =============================================================================
# `:corset` is provided by the ADRIACORSETExt extension (active because CORSETInterface is loaded).
# CORSET is pure Julia (no Python bridge) like Kora, but unlike Kora/CoralPop/C~scape it has no coral
# size structure and no fish/urchin/benthic state an ADRIA domain carries at all -- the
# `CORSET_INITIAL_*_BIOMASS` CONFIG constants above are REQUIRED placeholders, not domain-derived.
# See `CORSETInterface.jl/src/domain_translator.jl`'s file header before trusting this section's
# absolute output.

const CORSET_RESULT_FILE = joinpath(FPATH, "corset_result.jld2")

if CORSET_MODE === :run
    corset_spec_ = corset_spec(; rcp=parse(Int, RCP), year_start=first(YEARS), year_end=last(YEARS))
    scens_corset = ADRIA.sample(corset_spec_, N_SAMPLES)
    corset_dir = joinpath(FPATH, "corset")
    try
        rm(corset_dir; force=true, recursive=true)
    catch err
        @warn "Could not clear the previous CORSET output folder (OneDrive may be syncing it); continuing" corset_dir exception = err
    end
    rs_corset = ADRIA.run_scenarios(
        dom, scens_corset, RCP, :corset;
        corset_fpath=corset_dir, dhw_idx=DHW_IDX, wave_idx=WAVE_IDX, cyclone_idx=CYCLONE_IDX,
        n_workers=CORSET_N_WORKERS, spec=corset_spec_,
        translator_kwargs=(;
            initial_herbivore_biomass=CORSET_INITIAL_HERBIVORE_BIOMASS,
            initial_small_piscivore_biomass=CORSET_INITIAL_SMALL_PISCIVORE_BIOMASS,
            initial_large_piscivore_biomass=CORSET_INITIAL_LARGE_PISCIVORE_BIOMASS,
            initial_urchin_biomass=CORSET_INITIAL_URCHIN_BIOMASS,
            use_cyclones_as_hurricanes=CORSET_USE_CYCLONES_AS_HURRICANES,
            use_dhw_as_bleaching=CORSET_USE_DHW_AS_BLEACHING
        )
    )
    save_corset_results(rs_corset, CORSET_RESULT_FILE)
    corset_varying = corset_spec_.fieldname[.!corset_spec_.is_constant]
    JLD2.save(SAMPLES_FILE, merge(JLD2.load(SAMPLES_FILE),
        Dict("scens_corset" => scens_corset, "corset_varying" => String.(corset_varying))))
else
    isfile(CORSET_RESULT_FILE) || error("CORSET_MODE=:load but no earlier run saved at $CORSET_RESULT_FILE")
    rs_corset = load_corset_results(CORSET_RESULT_FILE)
    corset_varying = Symbol.(JLD2.load(SAMPLES_FILE)["corset_varying"])
end
push!(engines, EngineRun(:corset, "CORSET", "#E69F00", rs_corset,
    DataFrames.select(rs_corset.inputs, corset_varying), CORSET_GROUPS))

# =============================================================================
# 7. ReefMod -- through ADRIA's own entry point
# =============================================================================
# `:reefmod` is provided by the ADRIAReefModExt extension (active because ReefModInterface is
# loaded). Unlike CoralPop/Kora/CORSET, ReefMod is real MATLAB: each scenario is run as a `matlab
# -batch ...` SUBPROCESS, not an in-process call -- see `ReefModInterface.jl/src/reefmod_bridge.jl`.
# Needs a real `matlab` on PATH (ENV["REEFMOD_MATLAB_EXE"], set above) and the local `ReefMod-GBR`
# clone with ReefModInterface.jl's override hooks applied (ENV["REEFMOD_GBR_DIR"], also set above;
# see that clone's LOCAL_CHANGES.md). No MATLAB here? Run
# `ReefModInterface.jl/jobs/submit_reefmod_pipeline.sh` on an HPC with a Slurm-scheduled license
# instead, copy its `reefmod_result.jld2` to `REEFMOD_RESULT_FILE` below, and use
# `REEFMOD_MODE=:load` -- `load_reefmod_results` reads either source's output the same way.

const REEFMOD_RESULT_FILE = joinpath(FPATH, "reefmod_result.jld2")

if REEFMOD_MODE === :run
    reefmod_spec_ = reefmod_spec(; rcp=parse(Int, RCP), year_start=first(YEARS), year_end=last(YEARS))
    scens_reefmod = ADRIA.sample(reefmod_spec_, N_SAMPLES)
    reefmod_dir = joinpath(FPATH, "reefmod")
    try
        rm(reefmod_dir; force=true, recursive=true)
    catch err
        @warn "Could not clear the previous ReefMod output folder (OneDrive may be syncing it); continuing" reefmod_dir exception = err
    end
    rs_reefmod = ADRIA.run_scenarios(
        dom, scens_reefmod, RCP, :reefmod;
        reefmod_fpath=reefmod_dir, dhw_idx=DHW_IDX, n_workers=REEFMOD_N_WORKERS, spec=reefmod_spec_
    )
    save_reefmod_results(rs_reefmod, REEFMOD_RESULT_FILE)
elseif REEFMOD_MODE === :load
    if isfile(REEFMOD_RESULT_FILE)
        rs_reefmod = load_reefmod_results(REEFMOD_RESULT_FILE)
    else
        @warn "ReefMod results not found -- continuing without it. " *
              "Run with REEFMOD_MODE=:run (needs matlab on PATH), or copy a " *
              "ReefModInterface.jl/jobs/submit_reefmod_pipeline.sh result to REEFMOD_RESULT_FILE." REEFMOD_RESULT_FILE
    end
    ## CAUTION: unlike C~scape's dedicated settings check above, a :load here only gets the generic
    ## loc_ids/years check below -- if REEFMOD_RESULT_FILE came from a separate HPC pipeline run with a
    ## different DHW_IDX than this script's, that silently goes undetected (loc_ids/years don't depend
    ## on it). Safe when :load is reloading THIS project's own earlier :run (the common case); double
    ## check DHW_IDX by hand when copying in a result produced elsewhere.
end
if @isdefined(rs_reefmod)
    ## `reefmod_varying` comes from the result set's OWN `model_spec` (set from `spec=` whichever way
    ## `rs_reefmod` was produced -- :run above, or `load_reefmod_results(...; spec=...)` in
    ## `sandbox/03_load_hpc_results.jl`), not a separate side-channel -- works the same for either source.
    reefmod_varying = rs_reefmod.model_spec.fieldname[.!rs_reefmod.model_spec.is_constant]
    push!(engines, EngineRun(:reefmod, "ReefMod", "#56B4E9", rs_reefmod,
        DataFrames.select(rs_reefmod.inputs, reefmod_varying), REEFMOD_GROUPS))
end

## --- every engine must describe the same locations and years ---
for e in engines[2:end]
    collect(e.rs.loc_ids) == collect(engines[1].rs.loc_ids) ||
        error("$(e.label) and $(engines[1].label) cover different locations -- not comparable.")
    collect(ADRIA.timesteps(e.rs)) == YEARS ||
        error("$(e.label) covers years $(first(ADRIA.timesteps(e.rs)))-$(last(ADRIA.timesteps(e.rs))), expected $(first(YEARS))-$(last(YEARS)).")
end
@info "Engines in the comparison" labels = [e.label for e in engines] n_scenarios = [ADRIA.n_scenarios(e.rs) for e in engines]

# =============================================================================
# 8. The same ADRIA metrics for every engine
# =============================================================================

"k-area weighted domain mean of a `[timesteps, locations, scenarios]` metric -> `[timesteps, scenarios]`."
function domain_mean(x::AbstractArray{<:Real,3}, rs)::Matrix{Float64}
    w = ADRIA.loc_k_area(rs)
    w = w ./ sum(w)
    return dropdims(sum(Float64.(x) .* reshape(w, 1, :, 1); dims=2); dims=2)
end

## metric name => (title, ADRIA metric function). Every engine goes through ADRIA's own function, so the
## definitions cannot drift between engines.
const METRICS = [
    (:relative_cover, "Relative cover", M.relative_cover),
    (:relative_juveniles, "Relative juveniles", M.relative_juveniles),
    (:relative_shelter_volume, "Relative shelter volume", M.relative_shelter_volume),
    (:juvenile_indicator, "Juvenile indicator", M.juvenile_indicator),
    (:reef_condition_index, "Reef Condition Index", M.reef_condition_index),
    (:reef_tourism_index, "Reef Tourism Index", M.reef_tourism_index),
]

## domain[engine][metric] = [timesteps, scenarios], or `missing` when the engine's result set doesn't
## have that outcome at all -- ReefMod's indicators are a first pass (see ReefModInterface.jl/README.md's
## Status table): it has :relative_cover/:relative_taxa_cover/:absolute_shelter_volume and (sometimes)
## :relative_juveniles, but not yet :relative_shelter_volume/:juvenile_indicator/:reef_condition_index/
## :reef_tourism_index. Every other engine computes all six through ADRIAIndicators and is unaffected;
## this `try`/`catch` only ever fires for an engine with a genuinely incomplete `outcomes` dict.
domain = Dict(
    e.key => Dict(name => try
        domain_mean(Array(f(e.rs)), e.rs)
    catch err
        err isa Union{KeyError,MethodError} || rethrow()
        missing
    end for (name, _, f) in METRICS) for e in engines
)

## --- validity check: relative shelter volume is defined on [0, 1] (and must exist at all) ---
## Wrapping a result set in ADRIA's metric function only re-labels what the engine stored; each engine computed its
## own shelter volume. ADRIA's current definition (used for CoralBlox, and followed by CoralPop) normalises against a
## 95 cm tabular Acropora and clamps to [0, 1]. An engine whose stored values fall outside [0, 1] (or doesn't have
## them at all, e.g. ReefMod so far) used a different normalisation, and everything built on shelter volume (the
## Reef Condition and Tourism indices) inherits that. Plotting such an engine next to the others would be misleading,
## so it is left out of those panels -- with a warning -- unless you set INCLUDE_ENGINES_WITH_INVALID_SHELTER_VOLUME.
## Nothing is changed in the engines' own results.
const SHELTER_DEPENDENT = (:relative_shelter_volume, :reef_condition_index, :reef_tourism_index)
const INCLUDE_ENGINES_WITH_INVALID_SHELTER_VOLUME = _cfg("INCLUDE_INVALID_SHELTER", "false") == "true"

invalid_shelter = Set{Symbol}()
for e in engines
    rsv = domain[e.key][:relative_shelter_volume]
    if rsv === missing
        push!(invalid_shelter, e.key)
        @warn "$(e.label): no :relative_shelter_volume outcome at all. Its shelter volume and the Reef " *
              "Condition / Reef Tourism indices derived from it are omitted from those panels."
    elseif maximum(rsv) > 1 + 1e-9
        push!(invalid_shelter, e.key)
        @warn "$(e.label): relative shelter volume reaches $(round(maximum(rsv); sigdigits=3)) " *
              "($(round(100 * count(>(1 + 1e-9), rsv) / length(rsv); digits=1))% of values > 1), outside its defined range [0, 1]. " *
              "It was normalised differently from ADRIA's current definition, so its shelter volume and the Reef Condition / " *
              "Reef Tourism indices derived from it are not comparable with the other engines" *
              (INCLUDE_ENGINES_WITH_INVALID_SHELTER_VOLUME ? " (plotted anyway: INCLUDE_ENGINES_WITH_INVALID_SHELTER_VOLUME)." : " and are omitted from those panels.")
    end
end
## --- stale-result check: the C~scape juvenile-indicator units bug ---
## Before the fix in CscapeInterface (commit 0786416) the juvenile indicator was exactly 1e4 x too small (~1e-5),
## so a result saved by an older 06 shows juveniles but a juvenile indicator that is essentially zero. Skipped for
## any engine missing :relative_juveniles/:juvenile_indicator entirely (the `domain[...] === missing` case already
## keeps it out of those panels via `plottable`, no separate warning needed here).
for e in engines
    rel_juv, juv_ind = domain[e.key][:relative_juveniles], domain[e.key][:juvenile_indicator]
    (rel_juv === missing || juv_ind === missing) && continue
    if maximum(rel_juv) > 1e-2 && maximum(juv_ind) < 1e-3
        @warn "$(e.label): juvenile indicator is ~0 (max $(round(maximum(juv_ind); sigdigits=2))) although relative juveniles reach " *
              "$(round(maximum(rel_juv); sigdigits=2)). This looks like a result saved before the juvenile-indicator units fix. " *
              "Re-run 06, or repair its indicator files with CscapeInterface.fix_juvenile_indicators(\"<project>/adria_exports\") and " *
              "reload/save the result set. The juvenile indicator and stored Reef Condition Index for this engine are unreliable."
    end
end

"Is `e`'s `metric` fit to plot beside the other engines?"
function plottable(e, metric)
    domain[e.key][metric] === missing && return false
    return INCLUDE_ENGINES_WITH_INVALID_SHELTER_VOLUME || !(metric in SHELTER_DEPENDENT && e.key in invalid_shelter)
end

# =============================================================================
# 9. Figures
# =============================================================================

## --- 9.1 all engines, one panel per metric: mean over scenarios with a min-max band ---
fig_cmp = CairoMakie.Figure(; size=(1300, 750))
for (i, (name, title, _)) in enumerate(METRICS)
    omitted = [e.label for e in engines if !plottable(e, name)]
    ax = CairoMakie.Axis(fig_cmp[fld1(i, 3), mod1(i, 3)];
        title=isempty(omitted) ? title : "$title\n($(join(omitted, ", ")) omitted: not available or outside [0, 1])",
        titlesize=isempty(omitted) ? 16 : 13, xlabel="Year")
    for e in engines
        plottable(e, name) || continue
        a = domain[e.key][name]
        c = CairoMakie.Makie.to_color(e.color)
        CairoMakie.band!(ax, YEARS, vec(minimum(a; dims=2)), vec(maximum(a; dims=2)); color=(c, 0.25))
        CairoMakie.lines!(ax, YEARS, vec(mean(a; dims=2)); color=c, linewidth=2.5, label=e.label)
    end
    i == 1 && CairoMakie.axislegend(ax; position=:lt)
end
CairoMakie.Label(fig_cmp[0, :], "Domain-mean metrics, mean and min-max over $(N_SAMPLES) sampled scenarios per engine"; fontsize=18)
CairoMakie.save(joinpath(IMAGES_DIR, "compare_metrics.png"), fig_cmp; px_per_unit=3)

## --- 9.2 each engine through ADRIA's own scenario plot ---
for e in engines
    dc = ADRIA.DataCube(domain[e.key][:relative_cover]; timesteps=YEARS, scenarios=1:ADRIA.n_scenarios(e.rs))
    fig = if "guided" in names(e.rs.inputs)
        ## CoralBlox and CoralPop inputs carry the intervention columns ADRIA groups scenarios by
        ADRIA.viz.scenarios(e.rs.inputs, dc; axis_opts=Dict{Symbol,Any}(:title => "$(e.label) relative cover"))
    else
        ## C~scape's inputs are its ScenarioID table: supply one "all scenarios" group directly
        f = CairoMakie.Figure(; size=(800, 300))
        g = f[1, 1] = CairoMakie.GridLayout()
        ax = CairoMakie.Axis(g[1, 1]; title="$(e.label) relative cover")
        ADRIA.viz.scenarios!(g, ax, dc, Dict(:all => trues(ADRIA.n_scenarios(e.rs))))
        f
    end
    CairoMakie.save(joinpath(IMAGES_DIR, "trajectories_$(e.key).png"), fig; px_per_unit=3)
end

## --- 9.3 taxa cover, each engine keeping its own functional groups ---
## CoralBlox and Kora use ADRIA's own 5 groups and need no special-casing; CoralPop, ReefMod and
## CORSET keep their own 3/6/2 groups (series labels only); C~scape's `:species` axis needs the
## workaround below.
for e in engines
    fig = if e.key === :cscape
        ## `relative_taxa_cover`'s Metric wrapper fails on CScapeResultSet's `:species` axis; use the plain
        ## function and rename the axis, and pass C~scape's own labels (taxonomy! hard-codes ADRIA's 5 names)
        taxa = YAXArrays.renameaxis!(M._relative_taxa_cover(e.rs), :species => :groups)
        f = CairoMakie.Figure(; size=(800, 600))
        g = f[1, 1] = CairoMakie.GridLayout()
        ADRIA.viz.taxonomy!(g, taxa, Dict(:CScape => trues(ADRIA.n_scenarios(e.rs))); series_opts=Dict(:labels => e.groups))
        f
    elseif e.key === :coralpop || e.key === :reefmod || e.key === :corset
        ADRIA.viz.taxonomy(e.rs.inputs, M.relative_taxa_cover(e.rs); series_opts=Dict(:labels => e.groups))
    else
        ADRIA.viz.taxonomy(e.rs.inputs, M.relative_taxa_cover(e.rs))
    end
    CairoMakie.save(joinpath(IMAGES_DIR, "taxa_$(e.key).png"), fig; px_per_unit=3)
end

# =============================================================================
# 10. PAWN sensitivity of final-year cover to each engine's own factors
# =============================================================================
## `S` scaled down from N_SAMPLES so slicing doesn't degenerate; this does NOT fix the small-sample bias
## (see the CAUTION in the header).
const PAWN_S = clamp(N_SAMPLES ÷ 2, 1, 10)

pawn_results = Dict{Symbol,Any}()
for e in engines
    y = domain[e.key][:relative_cover][end, :]
    pawn_results[e.key] = ADRIAanalysis.sensitivity.pawn(e.factors, y; S=PAWN_S)
    fig = ADRIA.viz.pawn(
        pawn_results[e.key];
        axis_opts=Dict{Symbol,Any}(:title => "$(e.label) ecology-driver sensitivity (PAWN, S=$PAWN_S)")
    )
    CairoMakie.save(joinpath(IMAGES_DIR, "pawn_$(e.key).png"), fig; px_per_unit=3)
end

# =============================================================================
# 11. Summary table
# =============================================================================

rows = NamedTuple[]
for e in engines, (name, title, _) in METRICS
    plottable(e, name) || continue                                       # see the shelter-volume validity check above
    a = domain[e.key][name]
    for t in unique((1, cld(length(YEARS), 2), length(YEARS)))          # first, middle and last year
        push!(rows, (engine=e.label, metric=String(name), year=YEARS[t],
            mean=mean(a[t, :]), min=minimum(a[t, :]), max=maximum(a[t, :])))
    end
end
summary = DataFrame(rows)
open(joinpath(FPATH, "summary.csv"), "w") do io
    println(io, join(names(summary), ","))
    for r in eachrow(summary)
        println(io, join((r[c] isa AbstractFloat ? @sprintf("%.6g", r[c]) : string(r[c]) for c in names(summary)), ","))
    end
end
println(DataFrames.subset(summary, :metric => ByRow(==("relative_cover"))))

@info "Comparison written" FPATH IMAGES_DIR
