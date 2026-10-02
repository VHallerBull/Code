# Running the six-engine comparison on the HPC

How to run CoralBlox, C~scape, CoralPop, Kora, CORSET and ReefMod on the ADRIA Test domain on the AIMS
HPC (`hpc-h001`, Slurm), from a fresh build of the container to the comparison figures. This covers what
to build, what to upload, and which scripts to run in which order.

## Contents

1. [Overview](#1-overview)
2. [One-time setup](#2-one-time-setup)
3. [Before every build](#3-before-every-build)
4. [Build the image (WSL)](#4-build-the-image-wsl)
5. [Upload](#5-upload)
6. [Run, in this order](#6-run-in-this-order)
7. [Retrieve the results](#7-retrieve-the-results)
8. [Configuration reference](#8-configuration-reference)
9. [Scaling up to the full run](#9-scaling-up-to-the-full-run)
10. [Troubleshooting](#10-troubleshooting)
11. [Reproducibility](#11-reproducibility)
12. [File reference](#12-file-reference)

---

## 1. Overview

```
 Windows laptop                      WSL (Ubuntu, same laptop)               AIMS HPC (hpc-h001, Slurm)
 --------------                      -------------------------               --------------------------
 Documents/GitHub/<12 clones>  --.
 CoralModelReview/Code/        --+-> stage_sources.sh -> ~/coral-build/
 ADRIA Test domain, C~scape    --'   build_image.sh   -> coral_engines.sif --rsync--> ~/CoralModelReview/
   demography                                          upload/jobs/        --rsync-->   coral_engines.sif
                                                                                        jobs/  runs/  logs/
                                                                                                 |
                                                                     preflight.sh, reefmod_check.sh (once)
                                                                     submit_all.sh <project>
                                                                                                 |
                       ModelRuns/<project>/images, FIGURES.md, *.csv  <--rsync--  runs/<project>/
```

**What runs where.** Everything runs inside one Apptainer/Singularity image (`coral_engines.sif`):
Julia 1.12.5 with the whole `Code/` environment precompiled, R, Python 3.12 for CoralPop, the ADRIA Test
domain and C~scape's converted demography. **MATLAB is the exception.** It cannot be put in the image
(licensing), so the jobs bind the HPC's own MATLAB installation into the container. ReefMod then runs
through ADRIA's `:reefmod` engine like the other five, using the HPC's MATLAB licence. If that MATLAB
cannot start inside the container, a fallback runs ReefMod's MATLAB step on the host instead
([6e](#6e-fallback-reefmod-with-matlab-outside-the-container)). Building needs no MATLAB at all.

**The job chain** (`submit_all.sh`). Six independent engine jobs, then one analysis job that starts once
all of them have finished, whether they succeeded or not:

```
cmr_coralblox ─┐
cmr_cscape    ─┤
cmr_coralpop  ─┤   each: 08_RUN_ENGINE.jl <engine>      afterany
cmr_kora      ─┼──────────────────────────────────────────────────>  cmr_analyse
cmr_corset    ─┤   -> runs/<project>/<engine>/                       collect_benchmarks.sh (sacct)
cmr_reefmod   ─┘                                                     09_ANALYSE_SIX_ENGINES.jl
                                                                     -> runs/<project>/images, FIGURES.md
```

**The test run.** The ADRIA Test domain (10 sites, 2025-2099, RCP 4.5), with **8 scenarios per engine**
(ADRIA's Sobol' sampler needs a power of 2; 8 is the nearest to "about 10"). The shared climate member is
the domain's median-severity one: DHW member 40, wave member 16 (cyclone data in this domain is all
zero). Each engine gets 1 CPU and runs serially, so the per-scenario benchmarks are comparable.

---

## 2. One-time setup

### WSL (laptop)

Apptainer 1.5.3, `rsync`, `git` and `python3` are already installed in the Ubuntu WSL distribution,
and your user has an `/etc/subuid` entry. That is everything an unprivileged `apptainer build --fakeroot`
needs, so no `sudo` is required. On a fresh WSL:

```bash
sudo add-apt-repository -y ppa:apptainer/ppa && sudo apt-get update && sudo apt-get install -y apptainer rsync
```

### HPC

```bash
ssh hpc-h001
mkdir -p ~/CoralModelReview/{runs,logs,jobs}
```

After the first upload ([5](#5-upload)), run the checks once from the login node:

```bash
cd ~/CoralModelReview && bash jobs/preflight.sh
```

`preflight.sh` changes nothing. It reports:

- **Layout:** the folders, the image and its checksum.
- **Singularity:** the module (`singularity/singularity.module`, since this cluster has no `apptainer`
  command; the `.sif` runs under either tool) and whether the image starts.
- **Slurm:** `sbatch`, `sacct` (needed for the memory/CPU benchmarks), the `defq` partition and
  `/usr/bin/time`.
- **MATLAB:** the available `matlab` modules, the licence variables, and MATLAB's version plus whether
  the Statistics toolbox (needed by ReefMod) and the Compiler are licensed.

ReefMod-GBR's README asks for **MATLAB R2023b or earlier**. If several modules exist, pick one with
`export MATLAB_MODULE=matlab/R2023b` (the jobs read `MATLAB_MODULE`). Note the number of MATLAB licence
seats too: it caps how many ReefMod scenarios can run at once.

---

## 3. Before every build

The image contains a **copy** of the local clones, so build only from code you're happy to run.

1. **Commit** what you want in the run, in every clone. `stage_sources.sh` records each clone's branch
   and commit in the image (`STAGED_SOURCES.md`). Any uncommitted changes are saved as
   `patches/<repo>.patch` so the run can still be reconstructed, but committed is cleaner. As of this
   writing, CScape.jl and CscapeInterface.jl have uncommitted changes.
2. **Resolve the environment** if any clone's `Project.toml` changed (Julia 1.12: `julia +release`;
   the juliaup default on this laptop is 1.11):
   ```bash
   cd Code
   julia +release --project=. -e "using Pkg; Pkg.resolve(); Pkg.instantiate()"
   ```
3. **Record the sources** and commit them with the code (`record_sources.jl` updates `SOURCES.md`):
   ```bash
   julia +release --project=. record_sources.jl
   git add -A && git commit -m "..."
   ```
4. **Optional local smoke test** (no MATLAB needed, about 30-60 min, mostly Kora). Do this before
   spending an hour on the build:
   ```powershell
   $env:CMR_MODEL_RUNS_DIR="C:\Users\vhaller\cmrsmoke"; $env:CMR_PROJECT_NAME="S"; $env:CMR_N_SAMPLES="2"
   foreach ($e in "coralblox","corset","coralpop","cscape","kora") { julia +release --project=. 08_RUN_ENGINE.jl $e }
   julia +release --project=. 09_ANALYSE_SIX_ENGINES.jl
   ```
   Then look at `C:\Users\vhaller\cmrsmoke\S\FIGURES.md`. Keep the output folder **short and outside
   OneDrive**: Kora's and C~scape's output file names are long, so a deep folder breaks Windows' 260-character
   path limit (the run then fails while saving), and OneDrive locks freshly written files. A smoke test from
   2026-10-02 is kept in `ModelRuns/LocalSmoke_TestDomain_N2/` as a preview of the figures.

---

## 4. Build the image (WSL)

Open an Ubuntu (WSL) terminal. The paths contain spaces, so quote them:

```bash
CODE="/mnt/c/Users/vhaller/OneDrive - Australian Institute of Marine Science/DecisionSupport/CoralModelReview/Code"
bash "$CODE/hpc/stage_sources.sh"     # ~2 min:  assembles ~/coral-build (prints STAGED_SOURCES.md at the end)
bash "$CODE/hpc/build_image.sh"       # ~50 min: builds + tests ~/coral-build/coral_engines.sif
```

**`stage_sources.sh`** copies the 12 engine clones (code only: no `.git`, `sandbox/`, `Outputs/`,
`docs/`; CscapeInterface's sandbox alone is ~97 GB of earlier runs), `Code/`, the Test domain and
C~scape's six `demog_julia_*.rds` files into `~/coral-build/coral`. It works on WSL's own disk because
building from `/mnt/c` (OneDrive) is slow and unreliable. It also:

- rewrites `Code/Project.toml` and `Manifest.toml`'s `[sources]` from `../../../../Documents/GitHub/<repo>`
  to `../repos/<repo>`;
- strips Windows line endings from the shell scripts;
- prepares `~/coral-build/upload/` (the job scripts, plus the ReefMod copy the fallback route needs).

Settings you can override (prefix the command):

| Variable | Default | |
|---|---|---|
| `GITHUB_SRC` | `/mnt/c/Users/vhaller/Documents/GitHub` | where the clones are |
| `DOMAIN_SRC` | `$GITHUB_SRC/ADRIA.jl/ADRIA/test/data/Test_domain` | the domain to bake in |
| `DEMOG_SRC` | CscapeInterface's `AdjustedRecruitment10_Disturbance/data` | converted C~scape demography |
| `CALIB_SRC` | *(none)* | a CoralBlox `calibrated_params.nc` to bake in (then set `CMR_CALIBRATION_PATH=/opt/coral/data/calibrated_params.nc` when running) |
| `BUILD_DIR` | `~/coral-build` | |

**`build_image.sh`** runs `apptainer build --fakeroot` (log: `~/coral-build/build.log`). The image is
based on `julia:1.12.5-bookworm` plus R, the MATLAB system libraries, micromamba Python and a
precompiled Julia depot. The build downloads packages, so it needs internet. It then runs the image's
self-test (all six engine packages load; Python is found), writes `coral_engines.sif.sha256` and
prints the upload commands. The image is about 1.8 GB (first build: 47 min on this laptop).

**Optional:** test the image in WSL before uploading. This is the same smoke test as in 3.4, inside
the image:

```bash
mkdir -p ~/coral-test
for e in coralblox corset coralpop cscape kora; do
  apptainer exec --bind ~/coral-test ~/coral-build/coral_engines.sif \
    env CMR_MODEL_RUNS_DIR=$HOME/coral-test CMR_PROJECT_NAME=Smoke CMR_N_SAMPLES=2 \
    julia /opt/coral/Code/08_RUN_ENGINE.jl $e
done
apptainer exec --bind ~/coral-test ~/coral-build/coral_engines.sif \
  env CMR_MODEL_RUNS_DIR=$HOME/coral-test CMR_PROJECT_NAME=Smoke julia /opt/coral/Code/09_ANALYSE_SIX_ENGINES.jl
```

---

## 5. Upload

Everything else (the code, the engines, the domain, the demography) is **inside the image**. Only
these files go to the HPC:

| What | From (WSL) | To (HPC) | When |
|---|---|---|---|
| The image | `~/coral-build/coral_engines.sif` (+ `.sha256`) | `~/CoralModelReview/` | every build |
| Job scripts | `~/coral-build/upload/jobs/` | `~/CoralModelReview/jobs/` | every build (they change rarely) |
| Provenance | `~/coral-build/upload/STAGED_SOURCES.md` | `~/CoralModelReview/` | every build |
| ReefMod for host MATLAB | `~/coral-build/upload/reefmod_host/` | `~/CoralModelReview/reefmod_host/` | **only** for the fallback route (6e) |

From WSL (`build_image.sh` prints these with your paths):

```bash
ssh hpc-h001 "mkdir -p ~/CoralModelReview/{runs,logs,jobs}"
rsync -avP ~/coral-build/coral_engines.sif ~/coral-build/coral_engines.sif.sha256 hpc-h001:~/CoralModelReview/
rsync -avP ~/coral-build/upload/jobs/ hpc-h001:~/CoralModelReview/jobs/
rsync -avP ~/coral-build/upload/STAGED_SOURCES.md hpc-h001:~/CoralModelReview/
# fallback route only:
rsync -avP ~/coral-build/upload/reefmod_host/ hpc-h001:~/CoralModelReview/reefmod_host/
```

`rsync -P` resumes an interrupted transfer of the multi-GB image. The layout on the HPC is then:

```
~/CoralModelReview/
    coral_engines.sif  coral_engines.sif.sha256  STAGED_SOURCES.md
    jobs/       common.sh run_engine.slurm analyse.slurm collect_benchmarks.sh submit_all.sh
                preflight.sh reefmod_check.sh reefmod_A/B/C.slurm
    runs/       one folder per project (results)
    logs/       Slurm's <job name>_<job id>.out files
    reefmod_host/   (fallback route only)
```

---

## 6. Run, in this order

### 6a. Preflight (login node, once per image)

```bash
cd ~/CoralModelReview && bash jobs/preflight.sh
```

Fix anything marked `[FAIL]` before continuing ([2](#2-one-time-setup) and [10](#10-troubleshooting)).

### 6b. ReefMod check (compute node, once; about 15-30 min)

This is the first time ReefMod runs against a real MATLAB, so check it interactively before the batch
run:

```bash
srun -p defq -c 1 --mem=16G -t 01:00:00 --pty bash ~/CoralModelReview/jobs/reefmod_check.sh
```

It does three things:

1. Starts the host MATLAB **inside** the container. If that fails, use the fallback (6e) and skip the
   rest of this check.
2. Runs 2 ReefMod scenarios on the Test domain through ADRIA's `:reefmod` engine, into
   `runs/ReefModCheck/`.
3. Runs `hpc/reefmod_inspect.jl` on the raw output. This checks the items ReefMod-GBR's
   `LOCAL_CHANGES.md` lists as "still to verify":
   - `coral_pct2D`'s shape and that its units are percent;
   - the DHW shape and step count;
   - the 6-monthly to annual indexing;
   - juvenile and size-frequency counts.

Every line should read `[ok]`. A `[CHECK]` line points at an assumption in ReefModInterface that a real
run contradicts. Fix that before trusting ReefMod's numbers; the other five engines are unaffected.

> Already fixed (ReefMod-GBR commit `9be6b5c`): upstream's June 2026 MCB update made every run with
> restoration switched off (the way ReefModInterface runs it) stop with "Unrecognized field name
> 'doing_cooling'". `LOCAL_CHANGES.md` item 5 has the details.

### 6c. Submit the test run (login node)

```bash
cd ~/CoralModelReview
bash jobs/submit_all.sh TestDomain_N8
```

This submits the six engine jobs (`cmr_<engine>`) and `cmr_analyse`, which waits for them
(`--dependency=afterany`). It writes the job ids to `runs/TestDomain_N8/jobs.txt` and prints the
monitoring commands. Resources per engine for this test run (all 1 CPU):

| Engine | Memory | Time limit | Notes |
|---|---|---|---|
| coralblox | 8 GB | 1 h | seconds per scenario |
| cscape | 24 GB | 6 h | one Distributed worker holding a scenario, plus the driver |
| coralpop | 8 GB | 2 h | Python inside the Julia process |
| kora | 16 GB | 6 h | individual-based: the slowest engine here (~5-7 GB and ~10-15 min per scenario locally) |
| corset | 8 GB | 2 h | |
| reefmod | 16 GB | 12 h | one MATLAB at a time, including a MATLAB start-up per scenario |

These are generous first guesses. The run's own benchmarks (figure section K) show what each engine
actually used, so the next run can request less.

Variations (every `CMR_*` variable in your shell is passed on to the jobs):

```bash
bash jobs/submit_all.sh TestDomain_N8 --engines "kora corset"      # a subset
CMR_N_SAMPLES=16 bash jobs/submit_all.sh TestDomain_N16             # more scenarios
CMR_DHW_IDX=19 bash jobs/submit_all.sh TestDomain_N8_dhw19          # another climate member
bash jobs/submit_all.sh TestDomain_N8 --reefmod-fallback            # ReefMod via host MATLAB (6e)
```

Use a new project name for every distinct configuration. Re-running a project replaces each engine's
previous output.

### 6d. Monitor

```bash
squeue -u $USER                                                     # queued / running
sacct -j $(awk '{print $2}' runs/TestDomain_N8/jobs.txt | paste -sd,) --format=JobName,State,Elapsed,MaxRSS
tail -f logs/cmr_kora_<jobid>.out                                   # one engine's log
ls runs/TestDomain_N8/*/benchmark.csv                               # engines that have finished
```

Each engine writes `runs/<project>/<engine>/` containing:

- `result.jld2`, or `result_path.txt` for CoralBlox;
- `samples.jld2`;
- `benchmark.csv`;
- `size_structure.jld2`;
- `run/`, the engine's own files.

Each also writes `runs/<project>/time_<engine>.log` (GNU time: wall time and peak memory).

**Re-run one engine** (for example after a failure), then the analysis:

```bash
cd ~/CoralModelReview/logs
sbatch --job-name=cmr_kora --mem=16G --time=04:00:00 --export=ALL,ENGINE=kora,CMR_PROJECT_NAME=TestDomain_N8 ../jobs/run_engine.slurm
# when it has finished, add its id to runs/TestDomain_N8/jobs.txt as "kora <jobid>" (so its accounting is collected), then:
sbatch --export=ALL,CMR_PROJECT_NAME=TestDomain_N8 ../jobs/analyse.slurm
```

The analysis can be re-run any time. It compares whichever engines have a result and lists the missing
ones in `FIGURES.md`.

### 6e. Fallback: ReefMod with MATLAB outside the container

Use this route only if 6b step 1 failed. First upload `reefmod_host/` ([5](#5-upload)), then:

```bash
REEFMOD_ARRAY_CONCURRENCY=4 bash jobs/submit_all.sh TestDomain_N8 --reefmod-fallback
```

ReefMod then runs in three stages:

- **A** (`reefmod_A.slurm`, container): samples ReefMod's factors and writes one MATLAB input per
  scenario.
- **B** (`reefmod_B.slurm`, **host** MATLAB, Slurm array): runs one scenario per task, at most
  `REEFMOD_ARRAY_CONCURRENCY` at once (one licence each).
- **C** (`reefmod_C.slurm`, container): builds the result in the same `runs/<project>/reefmod/` layout.

The analysis treats the result like any other engine's. Only ReefMod's benchmark differs: it sums the
array tasks' wall times.

---

## 7. Retrieve the results

The analysis job writes into `runs/<project>/`:

| File | |
|---|---|
| `FIGURES.md` | every figure with a caption, the settings, the caveats, and a table of what each engine represents |
| `images/A_sanity` … `images/K_benchmarks` | the figures (PNG) |
| `summary.csv` | per engine × metric × year (first, 2050, last): median, IQR, min, max over scenarios |
| `benchmark_summary.csv` | per engine: load time, time to first result, time per scenario, allocations, peak memory (Julia and Slurm), CPU time, disk |
| `slurm_accounting.csv` | raw `sacct` output for the project's jobs |

Copy them back to `ModelRuns/` (the raw `*/run/` folders can stay on the HPC). From WSL:

```bash
DEST="/mnt/c/Users/vhaller/OneDrive - Australian Institute of Marine Science/DecisionSupport/CoralModelReview/ModelRuns/TestDomain_N8"
mkdir -p "$DEST"
rsync -avP --exclude='*/run/' --exclude='*/warmup/' hpc-h001:~/CoralModelReview/runs/TestDomain_N8/ "$DEST/"
```

Leaving out `run/` keeps the copy small. Every `<engine>/result.jld2` and `samples.jld2` comes back,
so `09_ANALYSE_SIX_ENGINES.jl` can be re-run on the laptop with `CMR_PROJECT_NAME=TestDomain_N8`. The
one exception is CoralBlox, whose ADRIA result store sits in `coralblox/run/`: add
`--include='coralblox/run/***'` before the excludes if you want to re-analyse locally.

---

## 8. Configuration reference

All scripts read `CMR_<NAME>` environment variables; defaults are in `engine_config()` in `engines.jl`.

| Variable | Default | Meaning |
|---|---|---|
| `CMR_PROJECT_NAME` | `TestDomain_N8` | output folder under `runs/` (set by `submit_all.sh`) |
| `CMR_MODEL_RUNS_DIR` | `~/CoralModelReview/runs` (jobs) / `../ModelRuns` (laptop) | where projects are written |
| `CMR_DOMAIN_PATH` | `/opt/coral/data/Test_domain` (image) | ADRIA domain; the folder must keep its own name |
| `CMR_RCP` | `45` | |
| `CMR_DHW_IDX` / `CMR_WAVE_IDX` / `CMR_CYCLONE_IDX` | `40` / `16` / `1` | the climate member all engines share |
| `CMR_N_SAMPLES` | `8` | scenarios per engine (power of 2) |
| `CMR_SEED` | `1234` | seed for ADRIA's scrambled Sobol' sampler |
| `CMR_N_WORKERS` | `0` | parallelism; its meaning differs per engine (CoralPop/C~scape: Distributed workers; Kora/CORSET: thread cap, needs more CPUs; ReefMod: concurrent MATLABs) |
| `CMR_CSCAPE_THREADS` | `1` | threads per C~scape worker |
| `CMR_WARMUP` | `true` | run 1 scenario first so time to first result and steady-state cost are measured separately |
| `CMR_CALIBRATION_PATH` | *(empty: ADRIA defaults)* | CoralBlox calibration NetCDF |
| `CMR_CORSET_INITIAL_*_BIOMASS`, `CMR_CORSET_USE_*` | see `engines.jl` | CORSET's placeholder consumer biomass and disturbance switches |
| `MATLAB_MODULE` | `matlab` | the module the ReefMod jobs load |
| `REEFMOD_ARRAY_CONCURRENCY` | `4` | fallback route: simultaneous MATLAB tasks |
| `CMR_BASE`, `SIF`, `RUNS` | `~/CoralModelReview`, `$CMR_BASE/coral_engines.sif`, `$CMR_BASE/runs` | locations (`jobs/common.sh`) |

---

## 9. Scaling up to the full run

1. Read `benchmark_summary.csv` and figure section K from the test run: seconds per scenario, peak
   memory and disk per scenario, per engine.
2. Choose `N_SAMPLES`. PAWN needs N/S > 80, i.e. **N ≥ 1024** at S = 10. Below that, the sensitivity
   figures are labelled illustrative.
3. Size each engine's job as **time ≈ N × seconds-per-scenario ÷ parallelism**, with memory from peak
   memory × workers. Raise `--cpus-per-task` and `CMR_N_WORKERS` together (edit `resources()` in
   `submit_all.sh`). ReefMod's parallelism is capped by MATLAB licences.
4. For a real domain, re-stage with `DOMAIN_SRC=...`, or bind a domain folder into the container and
   set `CMR_DOMAIN_PATH`. Choose the climate member the same way (median severity: see the computation
   in `engines.jl`'s `engine_config` docstring), and remember the domain's own `WAVE_IDX`/`ub_med`
   caveats (06's notes).

---

## 10. Troubleshooting

| Symptom | Cause / fix |
|---|---|
| `sbatch: error: ... account` | Add `#SBATCH --account=...` to the job scripts. None was needed on `hpc-h001` so far. |
| `ERROR: image not found` | Upload the `.sif`, or set `SIF=`. |
| ReefMod: `matlab exited with code ...` mentioning missing `lib*.so` | A MATLAB system library is missing from the image. Add it to `coral_engines.def`'s apt list and rebuild, or use the fallback (6e). `container ldd $MATLAB_ROOT/bin/glnxa64/MATLAB \| grep "not found"` (after `source jobs/common.sh; setup_matlab`) lists what is missing. |
| ReefMod: licence error / "Cannot find a valid license" | Is the licence server reachable from compute nodes? Check `MLM_LICENSE_FILE`/`LM_LICENSE_FILE` after `module load matlab` (preflight prints them). They pass into the container unchanged; a licence *file* outside the MATLAB folder is bound automatically by `setup_matlab`. |
| ReefMod: "Unrecognized field name 'doing_cooling'" | The image was built from a ReefMod-GBR without commit `9be6b5c`. Re-stage and rebuild. |
| `Read-only file system` writing under `/opt/...` | A script wrote outside `CMR_MODEL_RUNS_DIR`. Check that `CMR_MODEL_RUNS_DIR` points into a bound folder (the jobs set it). |
| Julia recompiles at job start, or fails writing `~/.julia_coral` | The HPC CPU isn't covered by the image's `JULIA_CPU_TARGET`, or the home quota is full. A one-off recompile into `~/.julia_coral` is harmless. |
| `OutOfMemory` / job `OUT_OF_MEMORY` in `sacct` | Raise `--mem` for that engine in `submit_all.sh`. Kora needs ~6-7 GB per scenario on the Test domain, C~scape more on bigger domains. |
| `$'\r': command not found` | A script has Windows line endings. Re-stage (it strips them), or run `sed -i 's/\r$//' jobs/*`. |
| 09: "were run with different settings" | Two engines in one project were run with different `CMR_*` settings. Re-run the odd one out, or use separate projects. |
| A figure missing from `images/` | It failed. `FIGURES.md` lists failed figures with the error, and the other figures are still written. |
| `DHW_IDX=... is outside this domain's 1:50 DHW members` | The Test domain has 50 members (the Moore domain's 64 does not exist here). |
| Laptop only: `SystemError: opening file ...` while an engine saves | Windows' 260-character path limit: use a short `CMR_MODEL_RUNS_DIR` (see 3.4). |

---

## 11. Reproducibility

- `STAGED_SOURCES.md` inside the image (`/opt/coral/STAGED_SOURCES.md`; also uploaded next to the `.sif`)
  records every clone's branch and commit. `/opt/coral/patches/` holds any uncommitted changes that were
  baked in.
- `coral_engines.sif.sha256` identifies the image. Keep it with the results (copy it into
  `runs/<project>/`).
- Every `<engine>/samples.jld2` stores the sampled scenarios, the factor spec and the full settings
  (domain, member, years, N, host, date). `benchmark.csv` records the host and Slurm job id.
- Re-running with the same image and settings gives the same scenarios: ADRIA's Sobol' sampler is
  Owen-scrambled (random), so `08` seeds it with `CMR_SEED` (default 1234, recorded in `samples.jld2`),
  and each engine seeds its own stochasticity from the scenario's `seed` factor. The scenarios actually
  run are stored in `samples.jld2` either way.

---

## 12. File reference

| File | Where it runs | What it does |
|---|---|---|
| `Code/engines.jl` | image / laptop | per-engine config, factor specs, run calls, save/load, size-structure extraction, benchmark helpers |
| `Code/08_RUN_ENGINE.jl <engine>` | image / laptop | runs one engine (1-scenario warm-up plus N scenarios) and writes `<project>/<engine>/` |
| `Code/09_ANALYSE_SIX_ENGINES.jl` | image / laptop | loads all engines, checks fairness, writes figures A-K, `summary.csv`, `benchmark_summary.csv`, `FIGURES.md` |
| `hpc/stage_sources.sh` | WSL | assembles the build folder from the clones |
| `hpc/coral_engines.def` | WSL (build) | the image definition |
| `hpc/build_image.sh` | WSL | builds and tests the image, prints upload commands |
| `hpc/reefmod_inspect.jl` | image | checks a real ReefMod output against ReefModInterface's assumptions |
| `hpc/reefmod_fallback.jl prepare\|finish` | image | the Julia stages of the fallback ReefMod route |
| `jobs/common.sh` | HPC | shared locations, `container`, `setup_matlab`, `timed` |
| `jobs/preflight.sh` | HPC login | one-time checks |
| `jobs/reefmod_check.sh` | HPC compute (srun) | first real-MATLAB ReefMod run and its inspection |
| `jobs/submit_all.sh` | HPC login | submits the engine jobs plus the analysis |
| `jobs/run_engine.slurm` | HPC compute | one engine (`ENGINE=...`) |
| `jobs/analyse.slurm` | HPC compute | `collect_benchmarks.sh` then `09` |
| `jobs/collect_benchmarks.sh` | HPC | `sacct` for the project's jobs into `slurm_accounting.csv` |
| `jobs/reefmod_A/B/C.slurm` | HPC compute | the fallback ReefMod route |
