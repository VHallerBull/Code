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
| `02_design.jl <model> <n_base> <band>... [--sets=main,zero]` | `Code/sa/shapley` | `<Model>/model_run_pf_pilot/band<NN>/design_pf.jld2` and `design_pf_zero.jld2` (A, B, permutations) |
| `03_run_chunk.jl <model> <proc> <n_procs> <band:domain[:zero]>...` | `Code/sa/shapley` | one process's walks -> `chunks/` (`:zero`: the zero-DHW design -> `pf_zero_dhw/chunks/`) |
| `03b_run_zero.jl <model> <proc> <n_procs> <band:domain>...` | `Code/sa/shapley` | zero-DHW reference: each walk's A and B sets without heat stress (2 runs per base sample) -> `zero_chunks/` |
| `04_collect.jl <model> <band>...` | `Code/sa/shapley` | `pf_results.jld2` (+ `pf_zero.jld2`) per band and domain |
| `05_analyse.jl <model> <band>...` | `Code/sa/shapley` | `analysis/`: effects, group effects, noise, base samples needed, convergence figures, `PILOT.md` |
| `06_figures.jl <model> [year...]` | `Code/sa/shapley` | `analysis/final_figures/`: output trajectories (with the zero-DHW reference), per-input and per-group figures |
| `run_pilot.sh [n_base] [n_procs] [sets]` | | 02 -> 03 (parallel processes; sets `main zero` by default) -> 04 for bands 0.1, 0.2, 0.5 |
| `test/runtests.jl` | `Code/sa/shapley` | walks = SAShE's table, permutations, A/B independence, toy model |

**Bands:** ±band applies to the 13 domain multipliers and the model's native factors. The DHW generator
ranges stay at SamplingSetup.docx Table 2.

**Base samples needed:** N = (1.96 σ / (δ Var(Y)))², where σ is an input's increment SD and δ the target
95 % half-width as a share of Var(Y). Runs = N × (d + 1).

**Zero-DHW design (set `zero`):** the sensitivity analysis without heat stress, the counterpart of the CoralBlox
paper's "DHW = 0" column. `02_design.jl` writes it next to the main design (`design_pf_zero.jld2`): the DHW
inputs are dropped (d = 13 + 155 = 168 for CoralPop, 169 runs per base sample), every run uses the domain's
all-zero DHW member, and the other inputs keep exactly the main design's A and B values. Its walks' end points
are therefore the main design's parameter sets without heat stress, i.e. the zero-DHW reference. Run it with
`band:domain:zero` jobs in `03_run_chunk.jl` (results under `<band>/<domain>/pf_zero_dhw/`); 04-06 handle both
sets (`set` column in the tables, `_zero` figure names).

**Zero-DHW reference:** the A and B parameter sets of every walk without heat stress, drawn behind the outputs in
`outputs_time_<domain>.png`; not part of the main Shapley analysis, whose DHW inputs never reach zero. It comes
from `03b_run_zero.jl` (2 runs per base sample) or, when the zero design was run, from its walks' end points
(then 03b is not needed).
