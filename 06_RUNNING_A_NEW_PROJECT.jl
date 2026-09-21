#=
06_RUNNING_A_NEW_PROJECT.jl

Runnable counterpart to `06_RUNNING_A_NEW_PROJECT.md` -- the same data-folder-to-plots
walkthrough, as an executable script instead of copy-pasteable snippets. Creates a fresh
project under `../ModelRuns/` (see `paths.jl`), populates it, runs both engines, and plots the
comparison. EVERYTHING it writes goes to `ModelRuns/<FPATH's folder name>/`: the C~scape project,
the CoralBlox results (`coralblox/`), the saved result set (`my_run.jld2`), the sampled factor
tables (`samples.jld2`, read by `07_COMPARE_THREE_ENGINES.jl`) and the figures (`images/`).

HOW TO USE: edit every value in the "CONFIG -- edit these" section below for your own domain/
demography files, then run this file top to bottom from this folder (`julia --project=.
06_RUNNING_A_NEW_PROJECT.jl`, or step through it in a REPL started with `--project=.` -- it's a
linear script, not a function, precisely so you can `include()` it partway or re-run a later
section without re-running the (slow) simulation steps above it). See
`06_RUNNING_A_NEW_PROJECT.md` for the narrative explanation of every step and, in the
CscapeInterface.jl repo, `sandbox/ADRIAIntegration/04_USER_GUIDE.md` for the reference docs on
every function called here.

The environment (`Project.toml`/`Manifest.toml` next to this file) uses the local clones of ADRIA,
CScape, CscapeInterface and CoralPopInterface -- see `README.md`/`SOURCES.md`.
=#

include(joinpath(@__DIR__, "paths.jl"))       # CODE_DIR, MODEL_RUNS_DIR, model_run_dir(name)

using ADRIA, CscapeInterface, DataFrames, XLSX, Dates, Statistics, Plots
using ADRIAanalysis
## ADRIAviz's Makie extension (ADRIAvizMakieExt) needs ADRIAviz + all three of Makie/GeoMakie/
## GraphMakie loaded to activate -- CairoMakie brings in Makie itself. Only needed for section
## 6's `ADRIA.viz.pawn` figure and sections 5.1/5.2's `ADRIA.viz.scenarios`/`taxonomy` figures.
## `import`, not `using` -- Makie/CairoMakie export `plot`/`plot!`/`scatter`/etc., the same
## names Plots.jl exports and that other sections below already call bare; `using` both would
## make every earlier `plot(...)`/`plot!(...)` call ambiguous. `import` still loads the packages
## (activating the extension) without bringing their exports into this namespace.
##
## CairoMakie, not WGLMakie (PAWNSetup.jl's choice elsewhere in this repo) -- WGLMakie needs a
## live browser/websocket session to rasterise a figure for `save`, which times out in this
## kind of headless batch script (confirmed: `WGLMakie.save` hung for 100s and errored when
## tested here). CairoMakie is a pure software rasteriser with no display dependency, the
## standard choice for exactly this "run a script, write static PNGs" use case.
import ADRIAviz, CairoMakie, GeoMakie, GraphMakie
CairoMakie.activate!()
import CScape
import JLD2
## `YAXArrays.renameaxis!` -- used in 5.2 to rename `rs_cs`'s `:species` dim to `:groups`
## (`ADRIA.viz.taxonomy!` hard-requires that exact name; see 5.2's note).
import YAXArrays

# =============================================================================
# CONFIG -- edit these
# =============================================================================

## Where the six real, brms-fitted demography inputs actually live on disk today (NOT yet
## inside the project folder -- this script copies them in below). Map each C~scape functional
## type to its own source file; adjust names/date stamp to your own IPM-fitting workflow output.
const DEMOG_SOURCE_DIR = raw"C:\Users\vhaller\Documents\GitHub\CscapeInterface.jl\sandbox\ADRIAIntegration\projects\NewTest\data"
const DEMOG_INPUT_NAMES = Dict(
    "acro_table"     => "demog_inputs_acro_tableALLDAT_260903.rds",
    "acro_corym"     => "demog_inputs_acro_corymALLDAT_260903.rds",
    "corym_non_acro" => "demog_inputs_corym_non_acroALLDAT_260903.rds",
    "brooder"        => "demog_inputs_corym_non_acro_BRALLDAT_260903.rds",
    "small_massive"  => "demog_inputs_small_massiveALLDAT_260903.rds",
    "large_massive"  => "demog_inputs_large_massiveALLDAT_260903.rds",
)

## OPTIONAL: a folder that already holds the six converted `demog_julia_*` files for THESE SAME fitted
## inputs (e.g. the `data/` folder of an earlier project). If set, they are copied into the new project
## so the slow, one-time R/brms conversion in step 1.2 is skipped. Only valid while DEMOG_INPUT_NAMES
## still name the same fitted files (same date stamp) -- the converted file names are derived from them.
## Set to `nothing` to always convert from scratch.
const CONVERTED_DEMOG_DIR = raw"C:\Users\vhaller\Documents\GitHub\CscapeInterface.jl\sandbox\ADRIAIntegration\projects\AdjustedRecruitment10_Disturbance\data"

## An ADRIA domain package on disk (a folder or a `.zip` ADRIA can load directly).
const DOMAIN_PATH = raw"C:\Users\vhaller\Downloads\Moore_2025-08-22_v080"
const RCP = "45"

## Base scenario metadata -- must overlap the domain's own timeframe (checked below).
const REGION = "Offshore_Northern"
const YEAR_START = 2025
const YEAR_END = 2099
const BASE_SCENARIO_ID = 1010001

## Set to `nothing` to leave Coral/GrowthAcceleration/DepthAttenuation on ADRIA's defaults, or
## a path to a CoralBlox calibration NetCDF (`calibrated_params.nc`) -- passed straight through
## to `ADRIA.load_domain`'s `calib_params_fn` keyword, which builds the calibrated component
## instances baked into the domain's model spec (see `ADRIA.load_calib_params` /
## `ADRIA/src/io/calib_params.jl`). NOT a `Serialization.deserialize`-able Dict -- the file is
## a real NetCDF/HDF5 dataset.
const CORALBLOX_CALIBRATION_OVERRIDES_PATH = "C:\\Users\\vhaller\\Australian Institute of Marine Science\\AIMS Decision Support - General\\01_PROJECTS\\RRAP M&DS\\ADRIA\\ADRIA Domain data\\CoralBlox-ADRIA GBR Calibrated Parameters\\2026-08-11\\calibrated_params.nc"

## Turn off DHW/wave/cyclone forcing to isolate intrinsic growth/survival/recruitment dynamics
## (matches `NoDisturbanceComparison.md`). Set to `false` to run with disturbances on instead.
const NO_DISTURBANCE = false

## Pin C~scape's mortality/growth/recruitment scale factors to their neutral (no-op) values
## instead of sampling them, for an apples-to-apples "intrinsic dynamics" comparison against
## CoralBlox. `mortality_scale` is additive (no-op = 0.0); the other two are multiplicative
## (no-op = 1.0).
const NEUTRALISE_CSCAPE_SCALE_FACTORS = false

## ADRIA.sample's default Sobol' sampler requires a power of 2. (`CMR_N_SAMPLES` overrides it for a quick test run.)
const N_SAMPLES = parse(Int, get(ENV, "CMR_N_SAMPLES", "8"))

## C~scape parallelism: one Distributed worker with N threads, so CScape.jl's own per-site
## `Threads.@threads` loop (annual_loop.jl) actually runs in parallel -- addprocs alone leaves
## workers on Julia's default of 1 thread regardless of this driver's own thread count.
## Each worker holds a full C~scape scenario in memory, so too many workers -> OutOfMemoryError
## (all 8 scenarios failed at 8 workers). Lower this further if it still OOMs.
const CSCAPE_N_WORKERS = 2
const CSCAPE_EXEFLAGS = ["-t", "4"]

# =============================================================================
# 0. Project folder
# =============================================================================

## Everything this script writes lives under ModelRuns/<project name> (see paths.jl).
## (`CMR_PROJECT_NAME` overrides the name, e.g. to keep a test run apart from a real one.)
const PROJECT_NAME = get(ENV, "CMR_PROJECT_NAME", "AdjustedRecruitment10_Disturbance")
const FPATH = model_run_dir(PROJECT_NAME)
const DATA_DIR = joinpath(FPATH, "data")
const IMAGES_DIR = joinpath(FPATH, "images")
mkpath(DATA_DIR)
@info "Project folder ready" FPATH

# =============================================================================
# 1.1 Place the real demography inputs
# =============================================================================

const CSCAPE_FTS = ["acro_table", "acro_corym", "corym_non_acro", "small_massive", "large_massive", "brooder"]

for ft in CSCAPE_FTS
    src = joinpath(DEMOG_SOURCE_DIR, DEMOG_INPUT_NAMES[ft])
    isfile(src) || error(
        "Demography input for '$ft' not found at $src -- edit DEMOG_SOURCE_DIR/" *
        "DEMOG_INPUT_NAMES at the top of this script before running."
    )
    dst = joinpath(DATA_DIR, DEMOG_INPUT_NAMES[ft])
    isfile(dst) || cp(src, dst)
end
@info "Real demography inputs staged in $DATA_DIR"

# =============================================================================
# 1.2 Convert demography to Julia-readable form
# =============================================================================
# One-time per file (skipped if already converted from a previous run of this script) --
# spins up an embedded R session via RCall to run brms::as_draws_matrix() on each file's two
# model objects. See src/demography_export.jl for exactly what gets converted and why.

for ft in CSCAPE_FTS
    ## Same naming convention `cscape_spec`/`_demog_julia_scenarios` expect elsewhere
    ## (naive "demog_inputs_" -> "demog_julia_" substring swap, suffix preserved) -- see
    ## `sensitivity_specs.jl`'s `_n_fitted_draws` and `cscape_jl_bridge.jl`'s
    ## `_demog_julia_scenarios`.
    julia_path = joinpath(DATA_DIR, replace(DEMOG_INPUT_NAMES[ft], "demog_inputs_" => "demog_julia_"))
    if isfile(julia_path)
        @info "Already converted, skipping" ft
        continue
    end
    ## Reuse a conversion done for an earlier project instead of repeating the slow R step.
    if CONVERTED_DEMOG_DIR !== nothing && isfile(joinpath(CONVERTED_DEMOG_DIR, basename(julia_path)))
        cp(joinpath(CONVERTED_DEMOG_DIR, basename(julia_path)), julia_path)
        @info "Reused an existing conversion (skipping R/brms)" ft from = CONVERTED_DEMOG_DIR
        continue
    end
    @info "Converting demography (R/brms, one-time)..." ft
    export_demography_for_julia(joinpath(DATA_DIR, DEMOG_INPUT_NAMES[ft]), julia_path)
end

# =============================================================================
# 1.3 Build a minimal ScenarioID.xlsx
# =============================================================================
# Spatial_file/Connectivity_file/Disturbance_file/init_cover are placeholders -- their content
# is never read once use_domain_translation=true (step 4) supplies the real thing from `dom`.
# fts/Growth_Surv_file/Region/Year_Start/Year_End are what actually get used.

const SCENARIO_XLSX = joinpath(FPATH, "ScenarioID.xlsx")

if isfile(SCENARIO_XLSX)
    @info "ScenarioID.xlsx already exists, leaving it as-is" SCENARIO_XLSX
else
    row = DataFrame(
        ID                      = [BASE_SCENARIO_ID],
        num_repeats             = [1],
        folder                  = ["scenarios"],
        Purpose                 = [""],
        Intervention            = ["No"],
        date                    = [Dates.format(Dates.today(), "yyyy/mm/dd")],
        Region                  = [REGION],
        Year_Start              = [YEAR_START],
        Year_End                = [YEAR_END],
        init_cover              = ["0.1_0.1_0.1_0.1_0.1_0.1"],  # placeholder, ignored (step 4)
        fts                     = [join(CSCAPE_FTS, "/")],
        HeatToleranceGroups     = ["-5_8_14"],
        HeatToleranceInit       = ["0_1.91"],
        Heritability            = ["0.5_0.01"],
        Plasticity              = [0.0],
        counterfactual          = [missing],
        Disturbance_file        = ["placeholder.RData"],
        Connectivity_file       = ["placeholder.RData"],
        Spatial_file            = ["placeholder.RData"],
        Geometry_file           = ["placeholder.RData"],
        Growth_Surv_file        = [join([DEMOG_INPUT_NAMES[ft] for ft in CSCAPE_FTS], "/")],
        Fogging_reducer         = [missing],
        Fogging_sites           = [missing],
        Fogging_start           = [missing],
        TradeOff                = [missing],
        var"Deployment area"    = [missing],
        TotalCorals             = [missing],
        species                 = [missing],
        Enhancement             = [missing],
        InterventionYears_start = [missing],
        duration                = [missing],
        frequency               = [missing],
        Reef_siteids            = [missing],
        temp_growth             = [1],
        species_proportions     = [missing],
        ## Independent of the CoralBlox-aligned recruitment_scale sensitivity factor below --
        ## see CScape.jl's RunConfig docstring. A single scalar (not per-ft), always applied,
        ## carried through unchanged into every sampled scenario by cscape_sample_to_scenarios.
        recruitment_multiplier  = [10.0],
    )
    XLSX.writetable(SCENARIO_XLSX, row; overwrite=true, sheetname="ScenarioID")
    @info "Wrote base scenario row" SCENARIO_XLSX
end

# =============================================================================
# 2. Load the ADRIA domain
# =============================================================================

dom = ADRIA.load_domain(
    DOMAIN_PATH, RCP;
    calib_params_fn=something(CORALBLOX_CALIBRATION_OVERRIDES_PATH, "")
)

## Sanity-check the scenario metadata actually overlaps this domain's own timeframe.
dom_years = collect(ADRIA.timesteps(dom))
(YEAR_START in dom_years && YEAR_END in dom_years) || @warn(
    "YEAR_START/YEAR_END ($YEAR_START-$YEAR_END) fall outside dom's own timeframe " *
    "($(first(dom_years))-$(last(dom_years))) -- update the CONFIG section."
)

## Both engines are meant to compare under one shared disturbance trajectory, not each draw
## independently (see InitialComparison.md/SENSITIVITY_COMPARISON.md/05_PHASE_D.md) --
## `coralblox_spec` itself enforces this for CoralBlox by fixing dhw_scenario/wave_scenario/
## cyclone_mortality_scenario to a constant (it does NOT expose them as ecology-driver factors
## to sample), and CScape's own run_scenarios(:cscape) only accepts one dhw_idx/wave_idx/
## cyclone_idx for its whole sampled batch anyway (see ADRIACscapeExt._run_scenarios_cscape).
##
## DHW_IDX=1 (this domain's own default) turned out to be a near-worst-case draw -- checked
## directly against all 550 available dhw_scens members, draw 1 has a domain/time-mean of
## ~23.9 DHW-weeks/yr, the 548th most severe of 550 (99.6th percentile) against a domain
## median of ~14.6. Using it made CScape's population collapse to zero, which is a property of
## that specific extreme draw, not a modelling bug. DHW_IDX=64 is that domain's median-severity
## draw instead (~14.6 DHW-weeks/yr). Cyclone draw 1 is already unremarkable (60th percentile,
## 0.044 vs domain median 0.042) so CYCLONE_IDX is left at 1. WAVE_IDX must stay 1 regardless --
## this domain has no wave data file, so wave_scens only has one (all-zero) member.
const DHW_IDX, WAVE_IDX, CYCLONE_IDX = NO_DISTURBANCE ? (0, 0, 0) : (64, 1, 1)

dom_nodist = deepcopy(dom)
ADRIA.fix_factor!(dom_nodist, :dhw_scenario, DHW_IDX)
ADRIA.fix_factor!(dom_nodist, :wave_scenario, WAVE_IDX)
ADRIA.fix_factor!(dom_nodist, :cyclone_mortality_scenario, CYCLONE_IDX)

# =============================================================================
# 3. Build a sampled scenario set per engine
# =============================================================================

cb_spec  = coralblox_spec(dom_nodist)   # fixes everything except ecology-driver factors
scens_cb = ADRIA.sample(cb_spec, N_SAMPLES)

## `cscape_spec`'s own defaults (mortality (-1,1), growth (0.5,1.5), recruitment (0,2)) are
## 5-10x wider than CoralBlox's calibrated ecology-driver bounds (ADRIA's `new_bounds`, a +-10%
## band around each factor's calibrated value from CORALBLOX_CALIBRATION_OVERRIDES_PATH -- see
## coral_factors.jl/Corals.jl). Narrowed here to a comparable +-10% so neither engine's spread
## is an artifact of a wider/narrower sampling band than the other's. `growth_scale`/
## `recruitment_scale` are multiplicative (no-op=1.0), so +-10% is literally (0.9,1.1); there's
## no exact unit correspondence for `mortality_scale` (additive, no-op=0.0) since CoralBlox's
## side is multiplicative there too, so (-0.1,0.1) is a same-order-of-magnitude proxy, not a
## derived equivalence.
cs_spec = cscape_spec(
    FPATH, BASE_SCENARIO_ID;
    mortality_bounds=(-0.1, 0.1), growth_bounds=(0.9, 1.1), recruitment_bounds=(0.9, 1.1)
)
if NEUTRALISE_CSCAPE_SCALE_FACTORS
    ## One factor per functional type now (see `cscape_spec`'s per-`ft` expansion) --
    ## neutralise each `<prefix>_<ft>` column individually.
    for (prefix, neutral_val) in (("mortality_scale", 0.0), ("growth_scale", 1.0), ("recruitment_scale", 1.0)),
        ft in CSCAPE_FTS

        idx = findfirst(==(Symbol("$(prefix)_$(ft)")), cs_spec.fieldname)
        cs_spec[idx, [:val, :dist_params, :is_constant]] = (neutral_val, (neutral_val, neutral_val), true)
    end
end
sampled = ADRIA.sample(cs_spec, N_SAMPLES)

## Persist the sampled factor tables and the forcing choices, so a later script (07) can reload
## this run's results and analyse them without re-sampling. (`dist` holds Distributions types, so
## it is stored by name.)
spec_for_disk(spec) = DataFrame(
    fieldname=spec.fieldname, val=spec.val, dist=string.(spec.dist),
    dist_params=spec.dist_params, is_constant=spec.is_constant
)
JLD2.save(
    joinpath(FPATH, "samples.jld2"),
    Dict(
        "scens_cb" => scens_cb, "cb_spec" => spec_for_disk(cb_spec),
        "sampled_cs" => sampled, "cs_spec" => spec_for_disk(cs_spec),
        "settings" => Dict{String,Any}(
            "DOMAIN_PATH" => DOMAIN_PATH, "RCP" => RCP, "DHW_IDX" => DHW_IDX, "WAVE_IDX" => WAVE_IDX,
            "CYCLONE_IDX" => CYCLONE_IDX, "YEAR_START" => YEAR_START, "YEAR_END" => YEAR_END,
            "NO_DISTURBANCE" => NO_DISTURBANCE, "N_SAMPLES" => N_SAMPLES,
            "CORALBLOX_CALIBRATION_OVERRIDES_PATH" => something(CORALBLOX_CALIBRATION_OVERRIDES_PATH, "")
        )
    )
)

# =============================================================================
# 4. Run both engines
# =============================================================================
# CoralBlox writes its results to ADRIA's output directory. ADRIA's default is `./Outputs` relative
# to the current working directory; it is redirected to this project's `coralblox/` folder below so
# nothing lands next to the code.

## ADRIA.setup() returns early (leaving its other defaults unset) if ENV["ADRIA_OUTPUT_DIR"] is
## already set, and resets ENV["ADRIA_LOG_COVER"] to false when no config.toml is found -- so call
## it explicitly first, THEN override the env vars, or the overrides are silently lost / the other
## defaults missing.
ADRIA.setup()
ENV["ADRIA_LOG_COVER"] = "true"      # per-species raw cover (rs_cb.coral_cover_log)
ENV["ADRIA_OUTPUT_DIR"] = joinpath(FPATH, "coralblox")
mkpath(ENV["ADRIA_OUTPUT_DIR"])
rs_cb = ADRIA.run_scenarios(dom_nodist, scens_cb, RCP)

rs_cs = ADRIA.run_scenarios(
    dom, sampled, RCP, :cscape;
    cscape_fpath=FPATH, cscape_base_scenario_id=BASE_SCENARIO_ID,
    use_domain_translation=true, dhw_idx=DHW_IDX, wave_idx=WAVE_IDX, cyclone_idx=CYCLONE_IDX,
    n_workers=CSCAPE_N_WORKERS, exeflags=CSCAPE_EXEFLAGS,
    export_adria=false, calc_indicators=true
)

save_results(rs_cs, joinpath(FPATH, "my_run.jld2"))
@info "Both engines run" fpath=FPATH

# =============================================================================
# 5. Plotting
# =============================================================================


mkpath(IMAGES_DIR)

## --- 5.1 Domain-mean trajectories ---
rel_cover_cb = Array(ADRIA.metrics.relative_cover(rs_cb))          # [timesteps, locations, scenarios]
rel_cover_cs = Array(rs_cs.outcomes[:relative_cover])              # [timesteps, locations, scenarios]

years_cb = collect(ADRIA.timesteps(rs_cb))
years_cs = collect(timesteps(rs_cs))

## `ADRIA.viz.scenarios!`'s confidence-band renderer (`_confints`) indexes its `outcomes`
## argument as `[timesteps, scenario-boolmask]` -- it needs a domain-already-summarised 2D
## `[timesteps, scenarios]` YAXArray (confirmed by testing against a saved run: passing the raw
## per-location 3D `relative_cover` throws a BoundsError there), not the raw per-location metric.
## Reduce over `:locations` by name (YAXArrays/DimensionalData dispatch on `dims=:locations`,
## same pattern as the positional `dims=2` reduction already used elsewhere in this script for
## the plain `Array` versions) rather than reaching for an ADRIA "scenario_"-summary metric,
## since those go through the same Metric-wrapper machinery that's already broken for `rs_cs`
## (see 5.2's note below).
domain_cover_cb = dropdims(mean(ADRIA.metrics.relative_cover(rs_cb); dims=:locations); dims=:locations)
domain_cover_cs = dropdims(mean(rs_cs.outcomes[:relative_cover]; dims=:locations); dims=:locations)

## `ADRIA.viz.scenarios(rs, outcomes)` groups scenarios by intervention type (counterfactual/
## unguided/guided -- derived from `rs.inputs.guided`/`.fogging`/`.SRM`/`.N_mc_settlers`/
## `.mcb_duration`, see ADRIAviz's `_scenario_types`). `rs_cb.inputs` is a genuine ADRIA scenario
## table (from `coralblox_spec`/`ADRIA.sample`) with all those columns (fixed constants, but
## present), so the plain top-level call works directly.
fig_traj_cb = ADRIA.viz.scenarios(
    rs_cb, domain_cover_cb; axis_opts=Dict{Symbol,Any}(:title => "CoralBlox relative cover")
)
CairoMakie.save(joinpath(IMAGES_DIR, "newproject_trajectories_coralblox.png"), fig_traj_cb; px_per_unit=3)

## `rs_cs.inputs` is CScape's own ScenarioID.xlsx-derived table (ID/Region/fts/... columns) --
## none of the ADRIA intervention columns `_scenario_types` needs exist on it, so the top-level
## call would throw a column-not-found error. Drop to the lower-level mutating form instead and
## supply a single "all scenarios" group directly, rather than deriving one from columns that
## don't exist on this side. Both the domain-reduction above and this workaround were verified
## against a saved run's CScapeResultSet before being written here.
fig_traj_cs = let
    n = size(domain_cover_cs, :scenarios)
    scen_groups = Dict(:CScape => trues(n))
    fig = CairoMakie.Figure(; size=(800, 300))
    g = fig[1, 1] = CairoMakie.GridLayout()
    ax = CairoMakie.Axis(g[1, 1]; title="C~scape relative cover")
    ADRIA.viz.scenarios!(g, ax, domain_cover_cs, scen_groups)
    fig
end
CairoMakie.save(joinpath(IMAGES_DIR, "newproject_trajectories_cscape.png"), fig_traj_cs; px_per_unit=3)

## --- 5.2 Species (functional-group) cover ---
## Note: call the underlying `_relative_taxa_cover`, not the Metric-wrapped `relative_taxa_cover`,
## for `rs_cs` -- the wrapped version throws inside `fill_axes_metadata!` on CScapeResultSet's
## renamed axes (a pre-existing ADRIA issue).
cb_groups = ADRIA.functional_group_names()
cs_fts = CscapeInterface.CSCAPE_FUNCTIONAL_TYPES

## Same CoralBlox-works/CScape-needs-a-workaround split as 5.1, for the same reason
## (`ADRIA.viz.taxonomy`'s scenario grouping hits the same missing-columns issue on `rs_cs`).
fig_taxa_cb = ADRIA.viz.taxonomy(rs_cb.inputs, ADRIA.metrics.relative_taxa_cover(rs_cb))
CairoMakie.save(joinpath(IMAGES_DIR, "newproject_species_cover_coralblox.png"), fig_taxa_cb; px_per_unit=3)

## Two more CScape-specific fixups on top of the scen_groups workaround, both confirmed against
## a saved run before being written here:
## 1. `taxonomy!` hard-requires a dimension literally named `:groups` (`relative_taxa_cover.
##    groups`) -- `_relative_taxa_cover(rs_cs)`'s functional-type dimension is named `:species`
##    (CScape's own convention), so it errors with "type YAXArray has no field groups" unless
##    renamed first.
## 2. `taxonomy!` also hard-codes `ADRIA.functional_group_names()` (CoralBlox's 5 groups) as the
##    legend labels regardless of the data passed in -- with CScape's 6 functional types this
##    throws a BoundsError indexing the 6th name into a 5-element list. `series_opts[:labels]`
##    overrides that default, so pass CScape's own names explicitly instead.
fig_taxa_cs = let
    taxa_cs = YAXArrays.renameaxis!(ADRIA.metrics._relative_taxa_cover(rs_cs), :species => :groups)
    n = size(taxa_cs, :scenarios)
    scen_groups = Dict(:CScape => trues(n))
    fig = CairoMakie.Figure(; size=(800, 600))
    g = fig[1, 1] = CairoMakie.GridLayout()
    ADRIA.viz.taxonomy!(g, taxa_cs, scen_groups; series_opts=Dict(:labels => cs_fts))
    fig
end
CairoMakie.save(joinpath(IMAGES_DIR, "newproject_species_cover_cscape.png"), fig_taxa_cs; px_per_unit=3)

## --- 5.3 Juveniles ---
juv_cb = Array(ADRIA.metrics.relative_juveniles(rs_cb))
juv_cs = Array(rs_cs.outcomes[:relative_taxa_juveniles])   # already domain-summed over locations

domain_juv_cb = vec(mean(dropdims(sum(juv_cb; dims=2); dims=2); dims=2))
domain_juv_cs = vec(mean(dropdims(sum(juv_cs; dims=2); dims=2); dims=2))

p_juv = plot(years_cb, domain_juv_cb; label="CoralBlox", lw=2, color=:coral)
plot!(p_juv, years_cs, domain_juv_cs; label="C~scape", lw=2, color=:steelblue)
xlabel!(p_juv, "Year"); ylabel!(p_juv, "Domain-total juvenile relative cover")
savefig(p_juv, joinpath(IMAGES_DIR, "newproject_juveniles.png"))

@info "Plots written" IMAGES_DIR

## --- 5.4 Advanced: size structure and thermal tolerance (raw out_array) ---
## Not behind a stable public function yet -- see 06_RUNNING_A_NEW_PROJECT.md §5.4 for the
## full explanation. Left as an opt-in block since it reaches past CScapeResultSet's normal
## interface and needs `size(out_array)` sanity-checked against CScape.Layout for your own run.
const RUN_ADVANCED_PLOTS = true

if RUN_ADVANCED_PLOTS
    ## Load one raw C~scape output. `CScape.load_output` expands the raw 2-slot on-disk array
    ## (non_deployed, deployed) into the 3-slot form with slot 3 = "combined" -- the same slot
    ## `rs_cs`'s indicators (and so plots 5.1-5.3) are computed from. One Array_scenario_*.jld2
    ## per sampled scenario (~3.6 GB each, ~5+ GB once expanded) -- change `ARRAY_SCENARIO_IDX`
    ## to inspect a different scenario.
    ARRAY_SCENARIO_IDX = 1
    COMBINED_IDX = 3   # {non_deployed, deployed, combined}
    cs_scenarios, _, _ = CscapeInterface.cscape_sample_to_scenarios(FPATH, BASE_SCENARIO_ID, sampled)
    cs_scenarios_jl = CscapeInterface._demog_julia_scenarios(FPATH; scenarios=cs_scenarios)
    cs_domain_inputs = CscapeInterface.domain_to_cscape_inputs(dom; dhw_idx=DHW_IDX, wave_idx=WAVE_IDX, cyclone_idx=CYCLONE_IDX)
    cs_array_file = sort(filter(f -> endswith(f, ".jld2"), readdir(joinpath(FPATH, "model_outputs"); join=true)))[ARRAY_SCENARIO_IDX]
    cs_m = match(r"Array_scenario_(\d+)_rep_(\d+)_draw_(.+)\.jld2$", basename(cs_array_file))
    cs_out = CScape.load_output(
        FPATH, parse(Int, cs_m[1]); rep_index=parse(Int, cs_m[2]), draw=parse(Int, cs_m[3]),
        scenarios=cs_scenarios_jl, spatial_override=cs_domain_inputs.spatial_override
    )
    out_array = cs_out.out_array   # (years, sites, intervention, ft, heat_bin, metric)
    size_by_year = dropdims(
        sum(out_array[:, :, COMBINED_IDX, :, :, CScape.Layout.VAL_SIZE]; dims=(2, 3, 4));
        dims=(2, 3, 4)
    )   # [years, 100 size classes]

    diam_cm = cs_out.meshpoints[1, :]

    p_size = plot(diam_cm, size_by_year[1, :]; label="$(years_cs[1])", xlabel="Colony diameter (cm)", ylabel="Cover")
    plot!(p_size, diam_cm, size_by_year[end, :]; label="$(years_cs[end])")
    savefig(p_size, joinpath(IMAGES_DIR, "newproject_size_structure.png"))

    n_bins = size(out_array, 5)
    bin_centres = range(-1, 1; length=n_bins)
    cover_by_bin = dropdims(sum(out_array[:, :, COMBINED_IDX, :, :, 1]; dims=(2, 3)); dims=(2, 3))
    weighted_tol = vec(sum(cover_by_bin .* bin_centres'; dims=2) ./ sum(cover_by_bin; dims=2))

    p_tol = plot(years_cs, weighted_tol; label="C~scape mean thermal tolerance", lw=2)
    savefig(p_tol, joinpath(IMAGES_DIR, "newproject_thermal_tolerance.png"))

    @info "Advanced plots written"
end

# =============================================================================
# 6. PAWN sensitivity analysis
# =============================================================================
## Global sensitivity of each engine's own ecology-driver factors (the same `is_constant=false`
## columns `ADRIA.sample` actually varied -- see `cb_spec`/`cs_spec` above) against a
## per-scenario outcome. Uses `ADRIAanalysis.sensitivity.pawn` (already a CscapeInterface.jl
## dependency) for the computation and the real `ADRIA.viz.pawn` (from `ADRIAviz`, added as a
## local dev dependency for this) for the figure -- a heatmap of every Si summary statistic
## (min/lb/mean/median/ub/max/std/cv) per factor, column-normalised to [0,1], sorted by
## `:median` descending. Saved with Makie's own `save` (via `CairoMakie.save`), not Plots'
## `savefig` -- the two plotting packages are loaded independently (see the `import` above) and
## don't share a figure type.
##
## CAUTION: PAWN needs N/S > 80 (N=N_SAMPLES, S=number of conditioning slices) to avoid biased
## results (Puy et al. 2020 -- see `ADRIAanalysis.sensitivity.pawn`'s docstring). At the
## default S=10 that means N_SAMPLES > 800. `PAWN_S` below is scaled down from N_SAMPLES so the
## slicing doesn't degenerate at small sample counts, but this does NOT fix the underlying
## bias -- at low N_SAMPLES (as configured above) these results are illustrative only, not a
## statistically valid sensitivity ranking. Increase N_SAMPLES before drawing real conclusions.
const PAWN_S = clamp(N_SAMPLES ÷ 2, 1, 10)

## Per-scenario outcome (NOT yet averaged across scenarios, unlike `mean_cb`/`mean_cs` in
## section 5.1): final-year domain-mean relative cover, one value per sampled scenario.
y_cover_cb = vec(mean(rel_cover_cb[end, :, :]; dims=1))
y_cover_cs = vec(mean(rel_cover_cs[end, :, :]; dims=1))

cb_varying = cb_spec.fieldname[.!cb_spec.is_constant]
cs_varying = cs_spec.fieldname[.!cs_spec.is_constant]

pawn_x_cb = DataFrames.select(scens_cb, cb_varying)
pawn_x_cs = DataFrames.select(sampled, cs_varying)

pawn_cb = ADRIAanalysis.sensitivity.pawn(pawn_x_cb, y_cover_cb; S=PAWN_S)
pawn_cs = ADRIAanalysis.sensitivity.pawn(pawn_x_cs, y_cover_cs; S=PAWN_S)

fig_pawn_cb = ADRIA.viz.pawn(
    pawn_cb; axis_opts=Dict{Symbol,Any}(:title => "CoralBlox ecology-driver sensitivity (PAWN, S=$PAWN_S)")
)
CairoMakie.save(joinpath(IMAGES_DIR, "newproject_pawn_coralblox.png"), fig_pawn_cb; px_per_unit=3)

fig_pawn_cs = ADRIA.viz.pawn(
    pawn_cs; axis_opts=Dict{Symbol,Any}(:title => "C~scape ecology-driver sensitivity (PAWN, S=$PAWN_S)")
)
CairoMakie.save(joinpath(IMAGES_DIR, "newproject_pawn_cscape.png"), fig_pawn_cs; px_per_unit=3)

@info "PAWN sensitivity analysis complete" images=IMAGES_DIR n_samples=N_SAMPLES S=PAWN_S
