# Code/sa/shapley: Shapley effects per year, ungrouped, grouped, within groups and along chains

Given-data Shapley effects (SamplingSetup.docx §2.8) of any output of a model's Sobol' runs, for every
simulated year. The estimator is SAShE.jl's nearest-neighbour pick-and-freeze estimator (`DataModel`). It
reads the files the run pipeline already writes: `results.jld2`, `result_set.jld2`,
`indicator_outputs.csv` and the domain's DHW.

```
julia +release -t auto --project=Code/sa/shapley Code/sa/shapley/run_shapley.jl Code/sa/shapley/configs/coralpop_relative_cover.toml
```

Run from the CoralModelReview folder. `-t auto` uses every core.

## Changing the output of interest

Edit `[output]` in the config, or override it on the command line:

```
... run_shapley.jl configs/coralpop_relative_cover.toml --output.indicator=reef_condition_index
... --output.indicator=relative_loc_taxa_cover --output.group=1          # Branching cover
... --output.reef=reef_03 --domains='["NineReefs_SobolDHW"]'              # one reef
... --output.transform=change_from_start
... --output.source=indicator_outputs --output.column=coral_diversity_2050   # a ready-made scalar
```

`output_select.jl` lists every option. Each output gets its own folder,
`<ModelFolder>/<run_subdir>/<domain>/shapley/<tag>/` (e.g. `relative_cover__domain`), so outputs never
overwrite each other. Copy the config to keep a set of choices.

## Run sets: main and zero-DHW control

`sets = ["main", "control"]` runs everything twice. The second run uses the zero-DHW control runs
(`control_results.jld2`, `control_result_set.jld2`) and writes to `shapley/<tag>__control/`; figure titles
say "(zero-DHW control)". The control `X` from `06_collect.jl` holds the DHW parameters of each row's
design member, but those runs used the all-zero member. They are therefore left out
(`[control] drop_inputs = ["^dhw_"]`), and chains through DHW are skipped.

## Analyses

| Analysis | Players | Effects sum to | File |
|---|---|---|---|
| ungrouped | every input (175 for CoralPop) | Var(Y) | `ungrouped.csv` |
| grouped, sum | members' effects added per permutation | Var(Y) | `grouped_sum_<scheme>.csv` |
| grouped, as player | each group is one player | Var(Y) | `grouped_player_<scheme>.csv` |
| convergence | ungrouped + grouped on nested prefixes of the runs | Var(Y) | `convergence.csv` |
| within group | only the group's inputs (restricted game) | Var(E[Y \| group]) | `within_<scheme>.csv` |
| chain | see below | | `chain_<name>.csv` |

"Sum" and "as player" differ only through interactions across groups. Shapley splits an interaction
equally among its players: one between two DHW inputs and one model input gives DHW 2/3 in the sum, but
1/2 as a player. With 175 noisy individual effects, the group-as-player estimate is also much less noisy.

Within-group effects rank a group's inputs by what they explain together, with every other input treated
as noise. `share_of_total` is each input's share of the group's explained variance.

**Chains** (inputs → intermediate M → output Y) are set in the config with `[[chain]]`. The upstream
group (`replaces`) is replaced by M in the inputs, the same substitution the design uses for DHW member →
7 parameters:
- **chain_A**: the restricted game of the upstream group on M.
- **chain_B**: the Shapley effects on Y_t of M_{t−lag} plus the other groups (or inputs).
- **chain_path**: Sh_Y(M) × (input's effect on M / Var(M)). This is Y's variance attributed to each
  upstream input through M, plus a row "unexplained part of M": what flows through M but is not described
  by the upstream inputs (for realised DHW: the weather realisation).
- **chain_upstream_group**: the upstream group as one player on Y, for comparison. Sh_Y(M) above it
  means M carries variance its parameters do not describe (for realised DHW: the weather realisation).
- **intermediate_\***: M analysed like any output.

For CoralPop the configured chain is DHW parameters → cumulative realised DHW → cover. For a model with
a thermal-tolerance indicator, use `intermediate = { indicator = "thermal_tolerance" }`.

`terminal = "var_y"` (default, as SAShE) makes the effects sum to Var(Y). Variance no input describes,
such as the weather realisation or seed noise, is then spread over all inputs. The inert inputs
(`inert` in the groups file, drawn dotted in `top_inputs.png`) show that floor. With
`terminal = "estimate"`, the full-analysis total is the estimated Var(E[Y | X]) instead, and
`variance.png` shows the unexplained remainder. That estimate is rough with 175 inputs.

## Cost

One neighbour search (a BallTree per coalition step) serves all years and outputs of an analysis. This
gives the same result as calling SAShE once per year with the same seed, at 1/26 of the cost.
`run_shapley.jl` prints the cost per permutation; set `n_permutations` from it. Within-group analyses
with large groups (60 mortality bins) and chains with `players = "individual"` are the expensive parts.
**Convergence** (`[convergence]` in the config, on by default) repeats the ungrouped and grouped analyses
on nested Sobol' prefixes of the runs (default N/8, N/4, N/2, N) for a few years. It writes
`convergence.csv` (column `n_rows`), `convergence_inputs.png` and `convergence_<scheme>.png`. Effects still
moving between the last two sizes are not converged (SamplingSetup.docx §2.6). `--estimator.n_rows=2048`
runs the whole analysis on one prefix.

## Files

| File | Does |
|---|---|
| `run_shapley.jl` | the runner: config → tables, `shapley.jld2`, figures, `SHAPLEY.md` |
| `shapley_core.jl` | `neighbour_plan`, `shapley_from_plan` (+ grouped sums), `chain_analysis` |
| `output_select.jl` | `RunData`, `output_matrix`: the output spec → `[runs, years]` |
| `groups.jl`, `groups/*.toml` | grouping schemes (regexes on input names); `shared.toml` holds the DHW and domain groups of every model |
| `plots.jl` | figures |
| `configs/*.toml` | analysis settings |
| `test/runtests.jl` | analytic checks, bit-for-bit equality with `SAShE.analyze`, groups-file checks |

## Adding a model

1. Write `groups/<model>.toml` (`extends = "shared.toml"`, the model's own groups, `inert`).
2. Add the model's interface package to `Project.toml` (`[deps]` + `[sources]`), so its result set loads.
3. Copy a config and set `model`.
