# Running a new project end-to-end, by hand

A single, linear walkthrough — data folder to plots — for running C~scape directly off a
real ADRIA `Domain` (Phase D) and comparing it against CoralBlox, without any of the
higher-level driver scripts. This is exactly the manual process behind
`NoDisturbanceComparison.md`; read that document for the *results* of running this, and
`04_USER_GUIDE.md` for the reference documentation of every function used below (this guide
just sequences them into one story and fills in the "how do I get a data folder in the first
place" step that guide assumes is already done).

Use this when you have a real ADRIA domain package and a set of real, brms-fitted demography
`.rds` files, and want to run C~scape (and optionally CoralBlox, for comparison) against them
from a bare Julia session.

> **Prefer to run this instead of copy-pasting it?** [`06_RUNNING_A_NEW_PROJECT.jl`](06_RUNNING_A_NEW_PROJECT.jl)
> is the same walkthrough as one executable script. Edit its `CONFIG` section (domain path,
> demography source files, sample size, ...) and run it top to bottom — it creates and
> populates `../ModelRuns/<project name>/` (a sibling of this `Code/` folder, see `paths.jl`),
> runs both engines, and writes plots to that project's `images/`. Run it with
> `julia --project=. 06_RUNNING_A_NEW_PROJECT.jl` from this folder.

---

## 0. Prerequisites

- The local sibling repos on the right branches (see `SOURCES.md` in this folder) and this
  folder's own Julia environment instantiated: `julia --project=. setup_environment.jl`. Its
  `Project.toml` lists the local clones under `[sources]`, so no manual `Pkg.develop` is needed
  (the older `CscapeInterface.jl`-project workflow is described in that repo's `04_USER_GUIDE.md`
  §1). Everything below assumes `julia --project=.` from this folder.
- An ADRIA domain package on disk (a folder or the `.zip` ADRIA can load directly), e.g. the
  RRAP M&DS `Moore_2025-08-22_v080` package.
- A set of real per-functional-type demography files,
  `demog_inputs_<ft>ALLDAT_<date>.rds` (brms `growth_regression`/`surv_regression` fits, one
  per functional type), for the 6 C~scape functional types:
  `acro_table`, `acro_corym`, `corym_non_acro`, `small_massive`, `large_massive`, `brooder`.
  These come from a separate IPM-fitting workflow (see `README.md`'s "Demographic files"
  section) — this guide starts from having them, not from fitting them.
- R with `brms` installed and reachable via RCall (only needed once, for the demography
  conversion step below — the simulation itself never touches R).

---

## 1. Set up the project's data folder

Unlike the classic C~scape workflow (`README.md`'s Workflow 1), running off a `Domain`
(Phase D) means you do **not** need to hand-build spatial/connectivity/disturbance `.RData`
files — those are translated from `dom` itself at run time. The project folder you still need
is much smaller: just real demography, and a `ScenarioID.xlsx` with one base row to carry
metadata the translator doesn't supply (functional-type list, demography filenames, region
name, year range).

```
MyProject/
├── ScenarioID.xlsx          # one base scenario row — built in step 1.3 below
└── data/
    ├── demog_inputs_acro_tableALLDAT_<date>.rds        # real brms fits, one per ft (step 1.1)
    ├── demog_inputs_acro_corymALLDAT_<date>.rds
    ├── demog_inputs_corym_non_acroALLDAT_<date>.rds
    ├── demog_inputs_corym_non_acro_BRALLDAT_<date>.rds # this ft's own internal `FT` field
    │                                                    # reads "corym_non_acro_BR", not
    │                                                    # "brooder" -- confirm which real
    │                                                    # file is meant to back "brooder"
    │                                                    # before relying on the filename alone
    ├── demog_inputs_small_massiveALLDAT_<date>.rds
    ├── demog_inputs_large_massiveALLDAT_<date>.rds
    ├── demog_julia_acro_table.rds        # converted, Julia-readable — step 1.2 writes these
    ├── demog_julia_acro_corym.rds
    ├── demog_julia_corym_non_acro.rds
    ├── demog_julia_brooder.rds
    ├── demog_julia_small_massive.rds
    └── demog_julia_large_massive.rds
```

### 1.1 Place the real demography inputs

Copy the six `demog_inputs_<ft>*.rds` files into `MyProject/data/`. Each is a `brmsfit`-based
R object Julia cannot load directly — don't try to `RData.load` these yet, that's what the
next step is for. If a filename's functional type isn't obvious (e.g. a `_BR` suffix), open it
in R and check its own internal `FT` field rather than guessing from the filename — confirmed
directly against real production files, the on-disk name and the `FT` field don't always agree
(see `NoDisturbanceComparison.md`'s "Real demography" section for a worked example: a file
named `..._BR...` internally reports `FT == "corym_non_acro_BR"`, which is C~scape's
`brooder`).

### 1.2 Convert demography to Julia-readable form

This is a one-time step per demography file, using the function moved into
`CscapeInterface.jl` for exactly this purpose (previously a standalone script in
`CScape.jl/data-raw/`):

```julia
using CscapeInterface

data_dir = "MyProject/data"
fts = ["acro_table", "acro_corym", "corym_non_acro", "small_massive", "large_massive", "brooder"]
# Map each ft to its real input filename -- adjust to your own file-naming convention/date stamp.
input_names = Dict(
    "acro_table"     => "demog_inputs_acro_tableALLDAT_260903.rds",
    "acro_corym"     => "demog_inputs_acro_corymALLDAT_260903.rds",
    "corym_non_acro" => "demog_inputs_corym_non_acroALLDAT_260903.rds",
    "brooder"        => "demog_inputs_corym_non_acro_BRALLDAT_260903.rds",
    "small_massive"  => "demog_inputs_small_massiveALLDAT_260903.rds",
    "large_massive"  => "demog_inputs_large_massiveALLDAT_260903.rds",
)

for ft in fts
    export_demography_for_julia(
        joinpath(data_dir, input_names[ft]),
        joinpath(data_dir, "demog_julia_$(ft).rds"),
    )
end
```

Each call spins up an embedded R session (via RCall) to run `brms::as_draws_matrix()` on that
file's two model objects, then writes a plain-list `.rds` Julia reads directly via RData.jl.
This is slow-ish (R/brms startup cost) but one-time — you never need R again for this project
after `data/demog_julia_<ft>.rds` exist. See
[`src/demography_export.jl`](../../src/demography_export.jl) for exactly what gets converted
and why R is unavoidable for this one step.

### 1.3 Build a minimal `ScenarioID.xlsx`

The domain translator (§3) supplies spatial/connectivity/disturbance/init-cover data straight
from `dom`, so the scenario table's own `Spatial_file`/`Connectivity_file`/`Disturbance_file`/
`init_cover` values are placeholders — their *content* is never read, they just need to be
non-empty strings so `CScape.jl`'s scenario-table validation doesn't reject the row. What
*does* matter: `fts`, `Growth_Surv_file` (the **original** `demog_inputs_...` names, not the
converted `demog_julia_...` ones — the translator does that filename substitution itself),
`Region`, and a `Year_Start`/`Year_End` range that actually overlaps your domain's own
timeframe (`ADRIA.timesteps(dom)`).

```julia
using DataFrames, XLSX, Dates

fpath = "MyProject"

row = DataFrame(
    ID                      = [1010001],
    RCP                     = [1],
    folder                  = ["scenarios"],
    Purpose                 = [""],
    Intervention            = ["No"],
    date                    = [Dates.format(Dates.today(), "yyyy/mm/dd")],
    Region                  = ["Moore"],                    # must match your build_demography region
    Year_Start              = [2025],
    Year_End                = [2099],
    init_cover              = ["0.1_0.1_0.1_0.1_0.1_0.1"],  # placeholder -- ignored when
                                                             # use_domain_translation=true (§4)
    fts                     = ["acro_table/acro_corym/corym_non_acro/small_massive/large_massive/brooder"],
    HeatToleranceGroups     = ["-5_8_14"],
    HeatToleranceInit       = ["0_1.91"],
    Heritability            = ["0.5_0.01"],
    Plasticity              = [0.0],
    counterfactual          = [missing],
    Disturbance_file        = ["placeholder.RData"],        # content never read (§3 translates this)
    Connectivity_file       = ["placeholder.RData"],         # content never read
    Spatial_file            = ["placeholder.RData"],         # content never read
    Geometry_file           = ["placeholder.RData"],
    Growth_Surv_file        = ["demog_inputs_acro_tableALLDAT_260903.rds/demog_inputs_acro_corymALLDAT_260903.rds/demog_inputs_corym_non_acroALLDAT_260903.rds/demog_inputs_small_massiveALLDAT_260903.rds/demog_inputs_large_massiveALLDAT_260903.rds/demog_inputs_corym_non_acro_BRALLDAT_260903.rds"],
    Fogging_reducer         = [missing],
    Fogging_sites           = [missing],
    Fogging_start           = [missing],
    TradeOff                = [missing],
    "Deployment area"       = [missing],
    TotalCorals             = [missing],
    species                 = [missing],
    Enhancement             = [missing],
    InterventionYears_start = [missing],
    duration                = [missing],
    frequency               = [missing],
    Reef_siteids            = [missing],
    temp_growth             = [1],
    species_proportions     = [missing],
)

XLSX.writetable(joinpath(fpath, "ScenarioID.xlsx"), row; overwrite=true, sheetname="ScenarioID")
```

> **`Growth_Surv_file` order must match `fts` order** — both are `/`-separated and read
> positionally. Double-check this if a functional type ends up with another one's demography.

You now have everything `cscape_spec`/`ADRIA.run_scenarios(..., :cscape)` need to look up
(§4): a base scenario row at `ID=1010001` to clone per sampled draw, and real
`demog_julia_<ft>.rds` files for `CScape.build_demography` to load.

---

## 2. Load the ADRIA domain

```julia
using ADRIA

dom = ADRIA.load_domain("path/to/MyProject_domain", "45")   # RCP as a string
```

Optionally fix disturbance factors off, to isolate intrinsic growth/survival/recruitment
dynamics from climate forcing (as `NoDisturbanceComparison.md` does):

```julia
dom_nodist = deepcopy(dom)
ADRIA.fix_factor!(dom_nodist, :dhw_scenario, 0)
ADRIA.fix_factor!(dom_nodist, :wave_scenario, 0)
ADRIA.fix_factor!(dom_nodist, :cyclone_mortality_scenario, 0)
```

---

## 3. Build a sampled scenario set per engine

If you have calibrated CoralBlox parameters (a `calibrated_params.nc`-style calibration
exercise), apply them to `dom_nodist` with `ADRIA.fix_factor!` per factor **before** calling
`coralblox_spec` — see `NoDisturbanceComparison.md`'s "Setup" table and `ManualCode.jl` for a
worked example loading them from a `Serialization.deserialize`d `Dict`:

```julia
using Serialization

overrides = Serialization.deserialize("coralblox_calibration_overrides.jls")
for (k, v) in overrides
    ADRIA.fix_factor!(dom_nodist, Symbol(k), v)
end
```

```julia
using CscapeInterface

fpath = "MyProject"
base_scenario_id = 1010001

# CoralBlox
cb_spec  = coralblox_spec(dom_nodist)   # fixes everything except ecology-driver factors
scens_cb = ADRIA.sample(cb_spec, 8)     # N must be a power of 2 for the default Sobol' sampler

# C~scape
cs_spec = cscape_spec(fpath, base_scenario_id)
sampled = ADRIA.sample(cs_spec, 8)
```

To run C~scape's own fitted demography **unperturbed** (rather than perturbed by
`mortality_scale`/`growth_scale`/`recruitment_scale`'s default sampled ranges — useful for an
apples-to-apples "intrinsic dynamics" comparison against CoralBlox's calibrated parameters),
pin those three factors to their neutral values before sampling
(`mortality_scale` is *additive*, so its no-op value is `0.0`; the other two are
*multiplicative*, so theirs is `1.0`):

```julia
for (fname, neutral_val) in (("mortality_scale", 0.0), ("growth_scale", 1.0), ("recruitment_scale", 1.0))
    idx = findfirst(==(Symbol(fname)), cs_spec.fieldname)
    cs_spec[idx, [:val, :dist_params, :is_constant]] = (neutral_val, (neutral_val, neutral_val), true)
end
sampled = ADRIA.sample(cs_spec, 8)   # re-sample after pinning
```

---

## 4. Run both engines

CoralBlox writes results to `./Outputs` relative to your current working directory (not
confined to a scratch dir) — `cd` somewhere disposable first if that matters to you.

```julia
# ADRIA.setup() resets ENV["ADRIA_LOG_COVER"] to false when no config.toml is found -- call it
# explicitly first if you want per-species raw cover (rs_cb.coral_cover_log) later, then
# override the env var, or the override is silently lost.
ADRIA.setup(); ENV["ADRIA_LOG_COVER"] = "true"
rs_cb = ADRIA.run_scenarios(dom_nodist, scens_cb, "45")

rs_cs = ADRIA.run_scenarios(
    dom, sampled, "45", :cscape;
    cscape_fpath=fpath, cscape_base_scenario_id=base_scenario_id,
    use_domain_translation=true, dhw_idx=0, wave_idx=0, cyclone_idx=0,
    # n_workers=1 + exeflags=["-t","8"]: one Distributed worker, 8 threads, so CScape.jl's own
    # per-site Threads.@threads loop (annual_loop.jl) actually runs in parallel -- addprocs
    # alone leaves workers on Julia's default of 1 thread regardless of the driver's -t flag.
    n_workers=1, exeflags=["-t", "8"], export_adria=false, calc_indicators=true
)
```

`dhw_idx=wave_idx=cyclone_idx=0` is C~scape's "no disturbance" convention, matching
`dhw_scenario=wave_scenario=cyclone_mortality_scenario=0` on the CoralBlox side (§2). Use
`>0` indices to select a real disturbance scenario member from `dom`'s own DHW/wave/cyclone
data instead. `rs_cs` is a `CscapeInterface.CScapeResultSet`; `rs_cb` is ADRIA's own
`ResultSet` — different types, same accessor interface (`.inputs`, `ADRIA.metrics.*`, etc.).

Save either result set so later analysis doesn't require re-running the simulation:

```julia
save_results(rs_cs, joinpath(fpath, "my_run.jld2"))
rs_cs = load_results(CscapeInterface.CScapeResultSet, joinpath(fpath, "my_run.jld2"))  # reload
```

---

## 5. Plotting

All examples below use `Plots.jl` for simplicity; swap in `ADRIA.viz.*`/`WGLMakie` (see
`README.md`'s Workflow 3) if you want ADRIA's own scenario/taxonomy/PAWN figure styles instead.

### 5.1 Domain-mean trajectories

```julia
using ADRIA, Statistics, Plots

rel_cover_cb = Array(ADRIA.metrics.relative_cover(rs_cb))          # [timesteps, locations, scenarios]
rel_cover_cs = Array(rs_cs.outcomes[:relative_cover])              # [timesteps, locations, scenarios]

years_cb = collect(ADRIA.timesteps(rs_cb))
years_cs = collect(timesteps(rs_cs))   # CscapeInterface re-exports timesteps for CScapeResultSet

mean_cb = vec(mean(dropdims(mean(rel_cover_cb; dims=2); dims=2); dims=2))   # mean over locs, then scenarios
mean_cs = vec(mean(dropdims(mean(rel_cover_cs; dims=2); dims=2); dims=2))

plot(years_cb, mean_cb; label="CoralBlox", lw=2, color=:coral)
plot!(years_cs, mean_cs; label="C~scape", lw=2, color=:steelblue)
xlabel!("Year"); ylabel!("Mean relative cover")
```

### 5.2 Species (functional-group) cover

C~scape's per-functional-type equivalent of CoralBlox's `relative_taxa_cover` is the
**underlying** function, not the `Metric`-wrapped one — the wrapped version throws inside
`fill_axes_metadata!` on `CScapeResultSet`'s renamed axes (a pre-existing ADRIA issue, not
specific to this workflow):

```julia
taxa_cb = Array(ADRIA.metrics.relative_taxa_cover(rs_cb))        # [timesteps, groups, scenarios]
taxa_cs = Array(ADRIA.metrics._relative_taxa_cover(rs_cs))       # call the underlying function directly

mean_taxa_cb = dropdims(mean(taxa_cb; dims=3); dims=3)   # [timesteps, groups]
mean_taxa_cs = dropdims(mean(taxa_cs; dims=3); dims=3)

cb_groups = ADRIA.functional_group_names(dom)            # 5 CoralBlox groups, in order
cs_fts    = CscapeInterface.CSCAPE_FUNCTIONAL_TYPES       # 6 C~scape functional types, in order

p1 = plot(title="CoralBlox")
for (i, g) in enumerate(cb_groups)
    plot!(p1, years_cb, mean_taxa_cb[:, i]; label=string(g), lw=2)
end
p2 = plot(title="C~scape")
for (i, ft) in enumerate(cs_fts)
    plot!(p2, years_cs, mean_taxa_cs[:, i]; label=ft, lw=2)
end
plot(p1, p2; layout=(1, 2), size=(1100, 400))
```

Remember the group correspondence isn't 1:1: C~scape's `brooder` and `corym_non_acro` both
draw their cyclone-mortality/init-cover data from CoralBlox's single `corymbose_non_Acropora`
group (`CscapeInterface.CSCAPE_FT_TO_ADRIA_GROUP`) — expect to sum those two when comparing
against that one CoralBlox group.

### 5.3 Juveniles

C~scape's juvenile indicator is precomputed — no extra extraction needed:

```julia
juv_cb = Array(ADRIA.metrics.relative_juveniles(rs_cb))            # [timesteps, locations, scenarios]
juv_cs = Array(rs_cs.outcomes[:relative_taxa_juveniles])           # [timesteps, groups, scenarios] (already domain-summed)

domain_juv_cb = vec(mean(dropdims(sum(juv_cb; dims=2); dims=2); dims=2))
domain_juv_cs = vec(mean(dropdims(sum(juv_cs; dims=2); dims=2); dims=2))

plot(years_cb, domain_juv_cb; label="CoralBlox", lw=2, color=:coral)
plot!(years_cs, domain_juv_cs; label="C~scape", lw=2, color=:steelblue)
xlabel!("Year"); ylabel!("Domain-total juvenile relative cover")
```

### 5.4 Advanced: size structure and thermal tolerance (raw `out_array`)

These two aren't behind a stable public function yet — `NoDisturbanceComparison.md` built
them by reading C~scape's raw `model_outputs/Array_scenario_*_draw_*.jld2` output directly.
Treat this section as a worked recipe to adapt, not a guaranteed-stable API: **check
`size(out_array)` against `CScape.Layout`'s constants for your own run before trusting the
numbers**, since this reaches past `CScapeResultSet`'s normal interface.

`out_array`'s shape is `(years, sites, intervention, functional_types, heat_tolerance_bins,
metrics)` — `intervention` is a combined dimension where index `3` is "combined" (see
`CScape.jl/src/io/load_output.jl`), and `metrics` is C~scape's fixed 107-slot layout
(`CScape.Layout`): index `1` is `cover`, indices `7:106` (`CScape.Layout.VAL_SIZE`) are the
100 IPM mesh size classes.

```julia
import CScape
using JLD2, Statistics

# One scenario/draw's raw array -- adjust the glob/filename to whatever run_cscape_native_parallel
# actually wrote under fpath/model_outputs/.
array_file = only(filter(f -> endswith(f, ".jld2"), readdir(joinpath(fpath, "model_outputs"); join=true)))
out_array = JLD2.load(array_file, "out_array")   # (years, sites, intervention, ft, heat_bin, metric)

combined_idx = 3   # CScape.jl's "combined" intervention slot

# --- Size structure: sum cover-in-each-size-class across sites/ft/heat_bin, per year ---
size_by_year = dropdims(
    sum(out_array[:, :, combined_idx, :, :, CScape.Layout.VAL_SIZE]; dims=(2, 3, 4));
    dims=(2, 3, 4)
)   # [years, 100 size classes]

# Convert mesh-point index -> colony diameter (cm) using the same demography build used for the run
demography = CScape.build_demography(
    cs_fts, replace.(split(row.Growth_Surv_file[1], "/"), "demog_inputs_" => "demog_julia_"),
    fpath, row.Region[1], domain_inputs.spatial_override
)
diam_cm = demography[1].meshpoints_diam   # same 100-point mesh shared across functional types

plot(diam_cm, size_by_year[1, :]; label="2025", xlabel="Colony diameter (cm)", ylabel="Cover")
plot!(diam_cm, size_by_year[end, :]; label="2099")

# --- Thermal tolerance: cover-weighted mean over the heat-tolerance bins ---
# Bins are evenly spaced and centred at 0 -- for the common 3-bin case that's [-1, 0, 1];
# confirm your own run's bin count/centres against config.heat_tolerance before reusing this.
n_bins = size(out_array, 5)
bin_centres = range(-1, 1; length=n_bins)

cover_by_bin = dropdims(sum(out_array[:, :, combined_idx, :, :, 1]; dims=(2, 3)); dims=(2, 3))  # [years, bins]
weighted_tol = vec(sum(cover_by_bin .* bin_centres'; dims=2) ./ sum(cover_by_bin; dims=2))

plot(years_cs, weighted_tol; label="C~scape mean thermal tolerance", lw=2)
```

For the CoralBlox side of a thermal-tolerance comparison, `rs_cb.coral_dhw_tol_log` is already
in DHW units (not the same scale as the bin-centred value above — compare shapes, not raw
numbers, unless you convert one to the other).

---

## 6. Troubleshooting

See `04_USER_GUIDE.md` §6 for the accumulated list of concrete errors hit building this
pipeline (power-of-2 sample sizes, `CScapeResultSet` name collisions, stale `Outputs/`
folders, mismatched `Year_Start`/`Year_End`, etc.) — everything there applies unchanged when
following this guide.

---

## 7. See also

- `04_USER_GUIDE.md` — the reference documentation for every function used above, organised
  by the phase that introduced it
- `NoDisturbanceComparison.md` — the actual results and findings from running this exact
  pipeline against the real `Moore` domain, including the full plotting/extraction
  methodology this guide's §5.4 is adapted from
- `05_PHASE_D.md` — the domain-translator's field-by-field design (what gets translated from
  `dom`, what doesn't, and why)
- [`src/demography_export.jl`](../../src/demography_export.jl) — the demography conversion
  function used in step 1.2
