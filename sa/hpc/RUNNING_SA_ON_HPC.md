# Sobol' sensitivity runs on the AIMS HPC (one image per model)

These runs are separate from the six/seven-engine comparison (`Code/hpc`, `coral_engines.sif`). Each
coral model gets its own slim image, `<model>_sa.sif`. All models share the same synthetic domains, DHW
ensemble and domain-multiplier block. Everything below is for **CoralPop**; to do another model,
replace `coralpop` once its `Code/sa/models/<model>.jl` and `models/<model>/` exist.

| Where | What |
|---|---|
| Windows | `CMR` = `C:\Users\vhaller\DecisionSupportProjects\CoralModelReview`. Domains, designs and (retrieved) results live here. |
| WSL | `/mnt/c/Users/vhaller/DecisionSupportProjects/CoralModelReview`. The image is built in `~/sa-build/<model>/`. |
| HPC | `ssh vhaller@hpc-l001.aims.gov.au`, then `~/CoralModelReview/sa/{images,jobs,domains,logs,CoralPOP/model_run}`. |

## 1. Inputs on Windows (shared steps once, then once per model)

```bash
cd /c/Users/vhaller/DecisionSupportProjects/CoralModelReview
julia +release --project=Code/sa Code/sa/01_build_domains.jl        # domains/SingleReef_SobolDHW, NineReefs_SobolDHW
julia +release --project=Code/sa Code/sa/02_dhw_ensemble.jl         # domains/*/sa_inputs/dhw_design.csv + QC figures
julia +release --project=Code/sa Code/sa/03_multiplier_sample.jl    # domains/*/sa_inputs/multiplier_sample.csv
julia +release --project=Code/sa/models/coralpop Code/sa/04_design.jl coralpop   # CoralPOP/model_run/*/design.jld2
```

Smoke test (16 rows, about a minute):

```bash
julia +release --project=Code/sa/models/coralpop Code/sa/05_run_chunk.jl coralpop SingleReef_SobolDHW main 1 16
```

## 2. Stage and build the image (WSL)

```bash
bash /mnt/c/Users/vhaller/DecisionSupportProjects/CoralModelReview/Code/sa/hpc/stage_model.sh coralpop
bash ~/sa-build/coralpop/build_model_image.sh coralpop
```

`stage_model.sh` copies the clones in `models/coralpop/hpc/repos.txt` and `Code/sa`. It records each
repo's commit, and any uncommitted patch, in `STAGED_SOURCES.md`. The build precompiles everything and
ends with the image's `%test`.

## 3. Upload (WSL)

```bash
HPC=vhaller@hpc-l001.aims.gov.au
CMR=/mnt/c/Users/vhaller/DecisionSupportProjects/CoralModelReview
ssh $HPC "mkdir -p ~/CoralModelReview/sa/{images,jobs,domains,logs}"
rsync -avP ~/sa-build/coralpop/coralpop_sa.sif ~/sa-build/coralpop/coralpop_sa.sif.sha256 $HPC:~/CoralModelReview/sa/images/
rsync -avP --chmod=go-w ~/sa-build/coralpop/upload/jobs/ $HPC:~/CoralModelReview/sa/jobs/
rsync -avP --chmod=go-w "$CMR/domains/SingleReef_SobolDHW" "$CMR/domains/NineReefs_SobolDHW" $HPC:~/CoralModelReview/sa/domains/
rsync -avP --chmod=go-w --exclude='chunks/' "$CMR/CoralPOP/model_run/" $HPC:~/CoralModelReview/sa/CoralPOP/model_run/
```

WSL shows `/mnt/c` files as 777, hence `--chmod=go-w`. The domains only need uploading once for all
models.

**New multiplier bands or base reef values need no new image.** The bands are read only on Windows by
`03`/`04` from `Code/sa/multiplier_bands.toml`, and saved inside the design. The base values come from the
domains (`01_build_domains.jl --force`). Re-run 01-04 as needed (see `Code/sa/README.md`, "Changing the
ranges"), then repeat only the last two `rsync` lines. Images built before 2026-10-20 still take the bands
for `06_collect.jl`'s result set from their own copy of the code (U(0.5, 1.5)), so rebuild once with the
current `Code/sa`.

## 4. Run (HPC)

```bash
ssh vhaller@hpc-l001.aims.gov.au
cd ~/CoralModelReview/sa/images && sha256sum -c coralpop_sa.sif.sha256
# test: 16 rows + 16 control rows per domain, in 2 tasks each, then collect + figures
bash ~/CoralModelReview/sa/jobs/sa_submit.sh --model coralpop --chunk 8 --n 16
# full run: 8192 rows (8 tasks of 1024) + 1024 control rows per domain
bash ~/CoralModelReview/sa/jobs/sa_submit.sh --model coralpop
squeue -u vhaller
```

Per domain, `sa_submit.sh` submits three jobs:
- the main array (`sa_chunk.slurm`);
- the control array;
- `sa_post.slurm`: `06_collect.jl`, `07_analyse.jl` and `08_indicator_outputs.jl`, after both arrays succeed (`afterok`).

The job ids go to `CoralPOP/model_run/<domain>/jobs.txt`. Options: `--domains`, `--chunk`, `--n`,
`--n-control`, `--set b2` (extension block) and `--dry-run`.

Expected cost (1 CPU per task): CoralPop takes about 0.04 s per run on SingleReef and about 0.2 s on
NineReefs (26 years). One 1024-run task is therefore a few minutes, plus about a minute to start Julia.

**A failed task** leaves its chunk file missing. `06_collect.jl` then lists the missing rows and stops
(so `sa_post` fails). Rerun only those tasks, then collect by hand:

```bash
cd ~/CoralModelReview/sa/logs
sbatch --array=3,7 --export=ALL,MODEL=coralpop,DOMAIN=NineReefs_SobolDHW,SET=main,CHUNK_SIZE=1024 ../jobs/sa_chunk.slurm
sbatch --export=ALL,MODEL=coralpop,DOMAIN=NineReefs_SobolDHW ../jobs/sa_post.slurm
```

**More runs** (if the Shapley estimates have not converged):
1. On Windows, run `04_design.jl coralpop --block 2` (another independent 8192 rows, `design_b2.*`).
2. Upload it.
3. Run `sa_submit.sh --model coralpop --set b2`.

`06_collect.jl` writes `results_b2.*`; the Shapley analysis simply uses both tables.

## 5. Retrieve (WSL)

```bash
# everything except the chunk files: results, result sets, figures AND every run's full outputs (~3 GB)
rsync -avP --exclude='chunks/' $HPC:~/CoralModelReview/sa/CoralPOP/model_run/ "$CMR/CoralPOP/model_run/"
# or without the per-run files (runs_*/): about 0.6 GB, enough for Shapley on any indicator and year
rsync -avP --exclude='chunks/' --exclude='runs_*/' $HPC:~/CoralModelReview/sa/CoralPOP/model_run/ "$CMR/CoralPOP/model_run/"
```

Every run keeps its full outputs on the HPC under `runs_<set>/`:
- `adria_exports/`: one ADRIA-indicator file per run;
- `model_outputs/`: one raw CoralPop file per run.

`06_collect.jl` stacks the indicators into `result_set.jld2`, so you only need `runs_*` locally to
recompute indicators or to read the raw size structure. Rough sizes for CoralPop, both domains (from
the smoke test):

| File | SingleReef | NineReefs |
|---|---|---|
| `runs_*/adria_exports` | 0.23 GB | 0.52 GB |
| `runs_*/model_outputs` (raw) | 0.30 GB | 1.6 GB |
| `result_set.jld2` + control | 0.08 GB | 0.36 GB |
| everything else | 0.1 GB | 0.1 GB |

To skip writing the per-run files, set `CMR_SA_KEEP_RAW=false` and/or `CMR_SA_KEEP_INDICATORS=false`
before `sa_submit.sh` (sbatch passes the environment on). Without indicators there is no
`result_set.jld2`.

Per domain, `CoralPOP/model_run/<domain>/` then holds:

| File | Contents |
|---|---|
| `results.jld2` / `results.csv` | `X` = 175 inputs (the DHW member index replaced by its 7 DHW generator parameters), `Y` = scalar outputs, plus series and per-run reef area / k. The input for SAShE, e.g. `DataModel(R["X"], R["Y"].rc_final)`. |
| `result_set.jld2` | The CoralPop ADRIA result set over the same runs in the same order: all 16 indicators (relative cover, RCI, shelter volume, juveniles, diversity, ...) per year × reef × run. Shapley on any of them at any year: see `Code/sa/README.md`. |
| `control_results.*` | The same parameter sets with zero DHW. Not part of the Shapley table. |
| `figures/`, `FIGURES.md`, `pawn.csv` | Ensemble fans and the PAWN check. |
