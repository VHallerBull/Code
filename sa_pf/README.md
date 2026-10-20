# Code/sa_pf: pick-and-freeze Shapley design (option B)

The callable pick-and-freeze design (Goda 2021; SAShE `CallablePickAndFreezeSample`). Each base sample is
a walk of d + 1 runs from a Monte Carlo point A[n] to an independent point B[n], swapping one input at a
time in a random order. The resulting estimates have no bias at any number of inputs and come with valid
confidence intervals. Code/sa (code, designs, outputs) is **not** changed; this folder only reads
`sa_common.jl` and `models/<model>.jl`.

**Option B:** the DHW block is the 7 generator parameters plus `dhw_weather` (1..1024, the column of
fixed random draws behind ensemble member w + 1). Every run generates its own trajectory
(`dhw_generator.jl`). With a member's parameters and weather index it reproduces the domain's DHW file
exactly, which `01_dhw_draws.jl` checks for all 1024 members.

**Shared across models:**
- The DHW and domain inputs get their A and B values from model-independent Monte Carlo streams.
- Their relative order within each walk is identical for every model.

So every model sees the same 9 trajectories and domain settings per base sample.

**Sampling:** plain Monte Carlo. Two `ADRIA.sample` (scrambled Sobol') draws with different seeds are
correlated row by row, up to |r| = 0.98, which pick-and-freeze cannot tolerate (`pf_common.jl`,
`sample_varying`).

## Outputs: the same for every model

Runs keep no files of their own: at 10^5-10^6 runs per domain the raw outputs cannot be stored or read
back. Each run instead returns a small **state summary**, built by the model's adapter
(`adapters/<model>.jl`) from its in-memory output. The outputs are derived from it after the runs
(`sa_outputs.jl`), so a changed definition needs no rerun.

| Output | Definition | Models |
|---|---|---|
| `relative_cover` | ADRIA relative cover (cover / habitable area k x area), habitable-area weighted over reefs | all |
| `acropora_share`, `massive_share` | Acropora (ADRIA groups 1+2) and massive (4+5) cover / total cover, on 3 harmonised classes (the third: corymbose non-Acropora) | all with morphology (not CORSET) |
| `relative_juveniles` | ADRIAIndicators `relative_juveniles` with the model's own juvenile classes | models with colony sizes (not CoCoNet, CORSET) |
| `coral_diversity` | ADRIAIndicators `coral_diversity` over the model's own groups | not CORSET |
| `thermal_<group>` | mean thermal tolerance per native group on the model's own scale and weighting (as Code/09's `thermal_mean_native.png`) | models with heritable tolerance (not CoralPop, CORSET) |
| `spatial_cv` | SD / mean of per-reef relative cover | domains with more than one reef |

State summary per run (Float32, about 6 KB on NineReefs): per reef and year the relative cover and the
3 class covers, `relative_juveniles` and `coral_diversity`; each reef's habitable area; per group and
year the mean tolerance and its weight (where the model has them). A year with a missing value in any
analysed walk is left out of that output, and of all thermal outputs together
(`analysis/dropped_years.csv`).

**A new model** needs only `adapters/<model>.jl`: `PF_CLASS_OF_GROUP` (from its `*_ADRIA_GROUPS` table via
`native_group_class`), `PF_THERMAL_GROUPS`, and `state_summary(r)` from its `run_row` result.

| Script | Env | Does |
|---|---|---|
| `01_dhw_draws.jl` | `Code/sa` | `pf_inputs/<domain>/dhw_draws.jld2` + exact-reproduction check |
| `02_design.jl <model> <n_base> <band>...` | `Code/sa/shapley` | `<Model>/model_run_pf_pilot/band<NN>/design_pf.jld2` (A, B, permutations) |
| `03_run_chunk.jl <model> <proc> <n_procs> <band:domain>...` | `Code/sa/shapley` | one process's walks -> `chunks/` (state summaries) |
| `03b_run_zero.jl <model> <proc> <n_procs> <band:domain>...` | `Code/sa/shapley` | zero-DHW reference: each walk's A and B sets without heat stress (2 runs per base sample) -> `zero_chunks/` |
| `04_collect.jl <model> <band>...` | `Code/sa/shapley` | `pf_results.jld2` (+ `pf_zero.jld2`) per band and domain: the state summaries and the outputs derived from them |
| `05_analyse.jl <model> <band>...` | `Code/sa/shapley` | `analysis/`: effects, group effects, noise, base samples needed, convergence figures, `PILOT.md` |
| `06_figures.jl <model> [year...]` | `Code/sa/shapley` | `analysis/final_figures/`: output trajectories (with the zero-DHW reference), per-input and per-group figures |
| `run_pilot.sh [n_base] [n_procs]` | | 02 -> 03 -> 03b (parallel processes) -> 04 for bands 0.1, 0.2, 0.5 |
| `sa_outputs.jl` | (any) | output definitions, state-summary storage and derivation (no package dependencies) |
| `adapters/<model>.jl` | | the model's state summary |
| `test/runtests.jl` | `Code/sa/shapley` | walks = SAShE's table, permutations, A/B independence, toy model, output derivations, CoralPop adapter on real runs |

**Bands:** ±band applies to the 13 domain multipliers and the model's native factors. The DHW generator
ranges stay at SamplingSetup.docx Table 2.

**Base samples needed:** N = (1.96 σ / (δ Var(Y)))², where σ is an input's increment SD and δ the target
95 % half-width as a share of Var(Y). Runs = N × (d + 1).

**Zero-DHW reference:** the A and B parameter sets of every walk rerun on the all-zero DHW member (as Code/sa's
control). They are drawn under the DHW trajectories in `outputs_time_<domain>.png` and are **not** part of the
Shapley analysis, whose DHW inputs never reach zero. A sensitivity analysis under no heat stress would need its
own design without the DHW inputs (155 + 13 inputs, 169 runs per base sample).
