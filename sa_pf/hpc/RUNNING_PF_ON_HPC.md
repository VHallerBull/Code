# Pick-and-freeze runs on the AIMS HPC

How to run the pick-and-freeze walks (`Code/sa_pf`) on the HPC. Same machines, folders and tools as the
Sobol' runs (`Code/sa/hpc/RUNNING_SA_ON_HPC.md`), with three differences:
- **The image** is `<model>_pf.sif`. It has the same definition and packages as `<model>_sa.sif`, plus
  `Code/sa_pf` in `/opt/sa/Code/sa_pf`. The Sobol' image has no `sa_pf`, so it cannot run these.
- **Outputs** go to `<ModelFolder>/model_run_pf/`. It never touches `model_run/` or the local pilot's
  `model_run_pf_pilot/`.
- **Collecting** happens on the HPC (`04_collect.jl`). The analysis and figures (`05`, `06`) run on Windows.

The example below is CoralPop, input band ±50 %, 2,500 base samples, SingleReef only, both sets:
- **Varied DHW:** 177 runs per base sample, 442,500 runs.
- **Zero DHW:** 169 runs per base sample, 422,500 runs.

In total that is 865,000 runs, about 15 CPU-h at the HPC's about 0.066 s per run.

| Where | What |
|---|---|
| Windows (Git Bash) | `cd /c/Users/vhaller/DecisionSupportProjects/CoralModelReview`. The design is written here and the results analysed here. |
| WSL | `CMR=/mnt/c/Users/vhaller/DecisionSupportProjects/CoralModelReview`. The image is built in `~/sa-build/coralpop_pf/`. |
| HPC | `ssh vhaller@hpc-l001.aims.gov.au`, then `~/CoralModelReview/sa/{images,jobs,domains,pf_inputs,logs,CoralPOP/model_run_pf}`. |

## 1. Design (Windows)

```bash
cd /c/Users/vhaller/DecisionSupportProjects/CoralModelReview
CMR_PF_SUBDIR=model_run_pf julia +release --project=Code/sa/shapley Code/sa_pf/02_design.jl coralpop 2500 0.5
```

This writes `CoralPOP/model_run_pf/band50/design_pf.{jld2,toml}` (varied DHW) and
`design_pf_zero.{jld2,toml}` (zero DHW), about 20 MB in all. The design does not depend on the domain.

The DHW draws `pf_inputs/<domain>/dhw_draws.jld2` already exist from the pilot. If they are missing, run
`julia +release --project=Code/sa Code/sa_pf/01_dhw_draws.jl`.

## 2. Stage and build the image (WSL, once per code change)

```bash
bash /mnt/c/Users/vhaller/DecisionSupportProjects/CoralModelReview/Code/sa_pf/hpc/stage_pf.sh coralpop
BUILD_DIR=~/sa-build/coralpop_pf bash ~/sa-build/coralpop_pf/build_model_image.sh coralpop
```

`stage_pf.sh` runs `Code/sa/hpc/stage_model.sh` into `~/sa-build/coralpop_pf`. That copies the repos and
`Code/sa`, records their commits in `STAGED_SOURCES.md`, and copies the Sobol' job scripts. It then adds
`Code/sa_pf` and the `pf_*` job scripts. The Sobol' build folder `~/sa-build/coralpop` is not touched.

The build is the same as for the Sobol' image: mostly Julia precompilation, ending with the image's
`%test`. Its closing upload hints name `coralpop_sa.sif`; use the commands in step 3 instead.

Commit `Code/sa_pf` first, so that `STAGED_SOURCES.md` records a commit rather than a patch.

**When to rebuild:**

| Change | Rebuild? |
|---|---|
| Any change to `Code/sa_pf/*.jl` or `adapters/` | Yes |
| New designs: another N, another band | No. Re-run step 1 and upload the designs |
| Job scripts | No. Re-upload `jobs/` |

## 3. Upload (WSL)

```bash
HPC=vhaller@hpc-l001.aims.gov.au
CMR=/mnt/c/Users/vhaller/DecisionSupportProjects/CoralModelReview
B=~/sa-build/coralpop_pf
ssh $HPC "mkdir -p ~/CoralModelReview/sa/{images,jobs,domains,pf_inputs,logs,CoralPOP/model_run_pf/band50}"
# the image, renamed to coralpop_pf.sif, with a checksum for that name
(cd $B && sha256sum coralpop_sa.sif | sed 's/coralpop_sa.sif/coralpop_pf.sif/' > coralpop_pf.sif.sha256)
rsync -avP $B/coralpop_sa.sif $HPC:~/CoralModelReview/sa/images/coralpop_pf.sif
rsync -avP $B/coralpop_pf.sif.sha256 $HPC:~/CoralModelReview/sa/images/
rsync -avP $B/upload/STAGED_SOURCES_coralpop.md $HPC:~/CoralModelReview/sa/images/STAGED_SOURCES_coralpop_pf.md
# job scripts (sa_common.sh is shared with the Sobol' jobs), domain, DHW draws, designs
rsync -avP --chmod=go-w $B/upload/jobs/ $HPC:~/CoralModelReview/sa/jobs/
rsync -avP --chmod=go-w "$CMR/domains/SingleReef_SobolDHW" $HPC:~/CoralModelReview/sa/domains/
rsync -avP --chmod=go-w "$CMR/pf_inputs/SingleReef_SobolDHW" $HPC:~/CoralModelReview/sa/pf_inputs/
rsync -avP --chmod=go-w "$CMR"/CoralPOP/model_run_pf/band50/design_pf* $HPC:~/CoralModelReview/sa/CoralPOP/model_run_pf/band50/
```

WSL shows `/mnt/c` files as 777, hence `--chmod=go-w`. The domain may already be there from the Sobol'
runs; rsync then copies nothing.

## 4. Smoke test (HPC, about 2 minutes)

The test runs one walk of each set into a separate folder, `model_run_pf_test`, so the real `chunks/`
folder stays empty. `04_collect.jl` refuses chunks that overlap.

```bash
ssh vhaller@hpc-l001.aims.gov.au
cd ~/CoralModelReview/sa/images && sha256sum -c coralpop_pf.sif.sha256
cd ~/CoralModelReview/sa
mkdir -p CoralPOP/model_run_pf_test/band50 && cp CoralPOP/model_run_pf/band50/design_pf* CoralPOP/model_run_pf_test/band50/
module load singularity/singularity.module
singularity exec --bind $PWD images/coralpop_pf.sif env CMR_SA_DOMAINS_DIR=$PWD/domains CMR_SA_RUNS_DIR=$PWD \
    CMR_PF_INPUTS_DIR=$PWD/pf_inputs CMR_PF_SUBDIR=model_run_pf_test \
    julia --startup-file=no /opt/sa/Code/sa_pf/03_run_chunk.jl coralpop --bases=1:1 \
    0.5:SingleReef_SobolDHW 0.5:SingleReef_SobolDHW:zero
rm -rf CoralPOP/model_run_pf_test
```

Expect two lines like `wrote .../walks_00001-00001.jld2 (1 walks, 0 failed, ...)`. The progress lines
give the time per run (about 0.07 s on SingleReef).

## 5. Run (HPC)

```bash
bash ~/CoralModelReview/sa/jobs/pf_submit.sh --model coralpop --band 0.5 --domains SingleReef_SobolDHW --dry-run   # check
bash ~/CoralModelReview/sa/jobs/pf_submit.sh --model coralpop --band 0.5 --domains SingleReef_SobolDHW
squeue -u vhaller
```

`pf_submit.sh` reads the base samples (n_base = 2500) from `design_pf*.toml` and submits:
- **`pf_coralpop_band50`:** an array of 50 tasks. Each task runs 50 base samples of both sets (8,650
  runs, about 20 minutes), into
  `CoralPOP/model_run_pf/band50/SingleReef_SobolDHW/chunks/walks_<first>-<last>.jld2` and
  `.../pf_zero_dhw/chunks/`.
- **`pf_coralpop_band50_collect`:** `04_collect.jl`, after every task has succeeded (`afterok`). It
  writes `pf_results.jld2` in both folders and prints the walks found, failed and missing.

The job ids go to `CoralPOP/model_run_pf/band50/jobs.txt`. The Slurm logs go to `~/CoralModelReview/sa/logs`,
and each task's time and memory to `CoralPOP/model_run_pf/band50/logs/`.

**Options:**

| Option | Use |
|---|---|
| `--sets main` or `--sets zero` | Run only one set |
| `--walks K` | Base samples per task (default 50) |
| `--time` | Wall time per task (default 1:30 h) |
| `--first F --last L` | Run only part of the design |
| `--no-collect` | Submit the walks without the collect job |

**A failed task** leaves its chunk files missing, and the collect job does not start. Find the failed tasks
with `sacct -j <job id> --format=JobID,State,Elapsed`. Task t covers base samples 50 (t − 1) + 1 to 50 t.
Rerun only those tasks, then collect:

```bash
cd ~/CoralModelReview/sa/logs
sbatch --array=7,12 --export=ALL,MODEL=coralpop,WALKS=50,FIRST=1,JOBS="0.5:SingleReef_SobolDHW 0.5:SingleReef_SobolDHW:zero" ../jobs/pf_chunk.slurm
sbatch --export=ALL,MODEL=coralpop,BAND=0.5 ../jobs/pf_collect.slurm
```

A finished chunk file is never rerun (the script skips it), so rerunning a whole task is safe.

## 6. Retrieve (WSL), then analyse (Windows)

```bash
rsync -avP --exclude='chunks/' $HPC:~/CoralModelReview/sa/CoralPOP/model_run_pf/ "$CMR/CoralPOP/model_run_pf/"
```

This retrieves about 1.6 GB: the two `pf_results.jld2` files, plus the job and time logs. The chunk
files, about 0.7 GB, stay on the HPC; add them only if you want to collect locally.

```bash
cd /c/Users/vhaller/DecisionSupportProjects/CoralModelReview
CMR_PF_SUBDIR=model_run_pf julia +release --project=Code/sa/shapley Code/sa_pf/05_analyse.jl coralpop 0.5
CMR_PF_SUBDIR=model_run_pf julia +release --project=Code/sa/shapley Code/sa_pf/06_figures.jl coralpop
```

The results go to `CoralPOP/model_run_pf/analysis/`:
- `effects.csv`, `group_effects.csv`, `noise.csv`;
- `required_base_samples.csv`, which shows whether 2,500 base samples reached ±0.05;
- figures and `PILOT.md`;
- `final_figures/` with `FIGURES.md`.

## More base samples later

The design is prefix-stable: a design with more base samples starts with exactly the same 2,500 walks
(`pf_common.jl`, tested in `test/runtests.jl`). To extend to, say, 5,000:
1. **Windows:** re-run step 1 with `5000` in place of `2500`. This overwrites the design files.
2. **Upload:** upload the designs again (the last `rsync` line of step 3).
3. **HPC:** submit only the new range:
   ```bash
   bash ~/CoralModelReview/sa/jobs/pf_submit.sh --model coralpop --band 0.5 --domains SingleReef_SobolDHW --first 2501
   ```
   `04_collect.jl` then merges the old and new chunk files.

Do not change `--walks` between runs of the same design if that would make chunk ranges overlap: the
collect step stops on any base sample found in two files.

## NineReefs

The same steps with `NineReefs_SobolDHW` in the domain and `pf_inputs` uploads and in `--domains`. Runs
are about 3.5 times slower, so use `--walks 15` and keep about 1:30 h per task.
