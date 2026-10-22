# Code/sa: Sobol' designs for the Shapley-effect sensitivity analysis

One plain ADRIA Sobol' sample per coral model, run on two synthetic domains. The Shapley analysis
itself uses SAShE.jl's given-data estimator (`DataModel(X, Y)`) on each model's outputs. It runs in
`shapley/` (its own environment; see `shapley/README.md`): every year, ungrouped, grouped, within
groups and along follow-on chains, with the output chosen in a config. The PAWN figures here are only
a check.

## Factor blocks

| Block | Factors | Shared? | Defined in |
|---|---|---|---|
| DHW | `dhw_scenario`: one member (2-1025) of a 1024-member DHW ensemble whose 7 generator parameters are a Sobol' sample. The 7 parameters replace the member index in `X` after the runs. | all models | `sa_common.jl` (`DHW_FACTORS`), `01_build_domains.jl` |
| domain | 13 multipliers, U(lo, hi) per multiplier from `multiplier_bands.toml` (default U(0.5, 1.5); lo = hi fixes one): depth, total initial cover, 5 functional-group cover weights, wave exposure, k, area, self-recruitment, between-reef connectivity strength, total larval retention | all models | `sa_common.jl` (`MULTIPLIER_FACTORS`, `multiplier_bands`, `apply_domain_multipliers`), `multiplier_bands.toml` |
| model | the model's native parameters (CoralPop: `coralpop_native_spec`, 155 incl. seed) | per model | `models/<model>.jl` |

## Changing the ranges

- **Domain multiplier bands:** edit `multiplier_bands.toml` (or point `CMR_SA_MULTIPLIER_BANDS` at
  another file), then re-run `03_multiplier_sample.jl` to check the effect and `04_design.jl <model>`,
  and upload the new designs. The bands are saved with the design (`design.jld2` `spec` and
  `settings["multiplier_bands"]`). The runs, `06_collect.jl` and `07_analyse.jl` take them from there,
  so **no image rebuild** is needed. A fixed multiplier (lo = hi) is applied to every run but is not a
  design factor, so it is absent from `X`.
- **Base reef values** (area, k, depth, wave exposure, initial cover, connectivity of each synthetic
  domain): edit `domain_configs` / `COMMON` in `01_build_domains.jl`, then run it with `--force`, then
  `02`, `03` and `04`, and upload the domains and designs. Domains are bound into the container at run
  time, so no rebuild is needed for this either.

Member 1 of the DHW file is all zero. The **control** runs the first 1024 design rows on it; it is not
part of the Shapley table.

## Sample sizes

**d** is the number of sampled factors; **N** is the number of rows, i.e. model runs.
- **DHW ensemble:** N = 1024 members for d = 7. These are generated trajectories, not model runs.
- **Design:** N = 8192 = 2¹³ per domain for d = 169 (175 inputs once the member index is swapped for
  its 7 parameters).
- **Control:** 1024 extra runs.
- **Convergence:** check on the nested prefixes, and add independent 8192-run blocks if needed
  (`04_design.jl <model> --block 2`).

## Scripts

| Script | Env | Does |
|---|---|---|
| `01_build_domains.jl` | `Code/sa` | `domains/SingleReef_SobolDHW` (one median reef), `domains/NineReefs_SobolDHW` (3x3 grid, seeded per-reef IQR draws): distance-decay connectivity, 2025-2050, the shared DHW ensemble |
| `02_dhw_ensemble.jl` | `Code/sa` | `sa_inputs/dhw_design.csv` (member -> 7 parameters + summaries), DHW QC figures |
| `03_multiplier_sample.jl` | `Code/sa` | `sa_inputs/multiplier_sample.csv`, what the multipliers (bands from `multiplier_bands.toml`) do to each domain |
| `04_design.jl <model>` | `models/<model>` | `<ModelFolder>/model_run/<domain>/design.jld2/.csv` |
| `05_run_chunk.jl <model> <domain> <set> <chunk> <size>` | `models/<model>` | one chunk of runs -> `chunks/` (one HPC array task) |
| `06_collect.jl <model>` | `models/<model>` | `results.jld2/.csv` (`X`, `Y`, series), `result_set.jld2` (every ADRIA indicator per year, reef and run), `control_*` |
| `07_analyse.jl <model>` | `models/<model>` | `figures/`, `pawn.csv`, `FIGURES.md` |
| `08_indicator_outputs.jl <model>` | `models/<model>` | `indicator_outputs.csv`: every indicator's domain value in 2030 / 2040 / 2050 / mean, one row per run (same order as `X`) |
| `test/runtests.jl` | `Code/sa` | helper tests |
| `shapley/run_shapley.jl <config>` | `shapley` | Shapley effects per year -> `<run dir>/shapley/<output>/` (tables, figures, `SHAPLEY.md`) |

## Shapley on any indicator and any year

Every run keeps its full outputs under `<run dir>/runs_<set>/`:
- `adria_exports/`: the ADRIA indicators, one file per run;
- `model_outputs/`: CoralPop's raw output (size bins, benthos, larvae).

`CMR_SA_KEEP_INDICATORS=false` / `CMR_SA_KEEP_RAW=false` switch either off. `06_collect.jl` stacks the
indicators into `result_set.jld2`, whose scenarios are in the same order as `X` in `results.jld2`:

```julia
include("Code/sa/sa_common.jl"); include("Code/sa/models/coralpop.jl")
R  = load_sa_results("coralpop", "NineReefs_SobolDHW")       # X, Y, reef_area, reef_k, years, ...
rs = load_sa_result_set("coralpop", "NineReefs_SobolDHW")    # ADRIA result set, same runs/order
keys(rs.outcomes)                                             # relative_cover, reef_condition_index, ...
rci = Array(rs.outcomes[:reef_condition_index])               # [years, reefs, runs]
dom_rci = domain_mean(rci, R["reef_k"] .* R["reef_area"])     # [years, runs], per-run habitable-area weights
y_2040 = dom_rci[findfirst(==(2040), R["years"]), :]          # one output per run, any year
# the same in one call, plus other choices:
y_2040 = indicator_output(R, rs, :reef_condition_index; year=2040)
y_reef = indicator_output(R, rs, :relative_cover; year=:mean, reef="reef_03")
y_morph = indicator_output(R, rs, :relative_loc_taxa_cover; year=2050, group=1)   # Branching
# SAShE: DataModel(R["X"], y_2040)
```

Local test runs live in `CoralPOP/model_run_local/` (set `ENV["CMR_SA_RUN_SUBDIR"] = "model_run_local"`
before including sa_common.jl); see its README.

Use `R["reef_k"] .* R["reef_area"]` for domain means: the multipliers change area and k in every
run, while `rs.loc_area` / `rs.loc_max_coral_cover` are only those of the first run.

`Code/sa` is the shared environment, local only (it includes SyntheticReefDomain). `models/<model>` is
the model's run environment; it is what goes into the model's own Apptainer image.

HPC: `hpc/RUNNING_SA_ON_HPC.md` (stage, build, upload, `jobs/sa_submit.sh`, retrieve).

## Adding the next model

1. Write `models/<model>.jl` with `MODEL_LABEL`, `model_spec_block()` and `run_row(dom, row)`; see
   `models/coralpop.jl`.
2. Create the `models/<model>/` folder:
   - `Project.toml` / `Manifest.toml`;
   - `hpc/<model>_sa.def`;
   - `hpc/repos.txt`.
3. Run `04_design.jl <model>`, stage, build, upload, then `sa_submit.sh --model <model>`.
