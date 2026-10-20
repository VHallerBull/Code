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

| Script | Env | Does |
|---|---|---|
| `01_dhw_draws.jl` | `Code/sa` | `pf_inputs/<domain>/dhw_draws.jld2` + exact-reproduction check |
| `02_design.jl <model> <n_base> <band>...` | `Code/sa/shapley` | `<Model>/model_run_pf_pilot/band<NN>/design_pf.jld2` (A, B, permutations) |
| `03_run_chunk.jl <model> <proc> <n_procs> <band:domain>...` | `Code/sa/shapley` | one process's walks -> `chunks/` |
| `03b_run_zero.jl <model> <proc> <n_procs> <band:domain>...` | `Code/sa/shapley` | zero-DHW reference: each walk's A and B sets without heat stress (2 runs per base sample) -> `zero_chunks/` |
| `04_collect.jl <model> <band>...` | `Code/sa/shapley` | `pf_results.jld2` (+ `pf_zero.jld2`) per band and domain |
| `05_analyse.jl <model> <band>...` | `Code/sa/shapley` | `analysis/`: effects, group effects, noise, base samples needed, convergence figures, `PILOT.md` |
| `06_figures.jl <model> [year...]` | `Code/sa/shapley` | `analysis/final_figures/`: output trajectories (with the zero-DHW reference), per-input and per-group figures |
| `run_pilot.sh [n_base] [n_procs]` | | 02 -> 03 -> 03b (parallel processes) -> 04 for bands 0.1, 0.2, 0.5 |
| `test/runtests.jl` | `Code/sa/shapley` | walks = SAShE's table, permutations, A/B independence, toy model |

**Bands:** ±band applies to the 13 domain multipliers and the model's native factors. The DHW generator
ranges stay at SamplingSetup.docx Table 2.

**Base samples needed:** N = (1.96 σ / (δ Var(Y)))², where σ is an input's increment SD and δ the target
95 % half-width as a share of Var(Y). Runs = N × (d + 1).

**Zero-DHW reference:** the A and B parameter sets of every walk rerun on the all-zero DHW member (as Code/sa's
control). They are drawn under the DHW trajectories in `outputs_time_<domain>.png` and are **not** part of the
Shapley analysis, whose DHW inputs never reach zero. A sensitivity analysis under no heat stress would need its
own design without the DHW inputs (155 + 13 inputs, 169 runs per base sample).
