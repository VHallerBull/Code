# Running the six-engine comparison on the HPC

How to run CoralBlox, C~scape, CoralPop, Kora, CORSET and ReefMod on the ADRIA Test domain on the AIMS
HPC (`hpc-l001.aims.gov.au`, Slurm), from a fresh build of the container to the comparison figures. This covers what
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
9. [Other domains and the full run](#9-other-domains-and-the-full-run)
10. [Troubleshooting](#10-troubleshooting)
11. [Reproducibility](#11-reproducibility)
12. [File reference](#12-file-reference)

---

## 1. Overview

```
 Windows laptop                      WSL (Ubuntu, same laptop)               AIMS HPC (hpc-l001, Slurm)
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
ssh vhaller@hpc-l001.aims.gov.au
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
4. **Optional local smoke test** (no MATLAB needed, about 15-30 min, mostly C~scape; add 30+ min the first time after a dependency changed, while Julia precompiles). Do this before
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
| `CALIB_SRC` | the shared drive's `CoralBlox-ADRIA GBR Calibrated Parameters/2026-08-11/calibrated_params.nc` | the CoralBlox calibration to bake in (required; `engines.jl` uses the baked-in copy by default) |
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
HPC=vhaller@hpc-l001.aims.gov.au
ssh $HPC "mkdir -p ~/CoralModelReview/{runs,logs,jobs}"
rsync -avP ~/coral-build/coral_engines.sif ~/coral-build/coral_engines.sif.sha256 $HPC:~/CoralModelReview/
rsync -avP --chmod=go-w ~/coral-build/upload/jobs/ $HPC:~/CoralModelReview/jobs/
rsync -avP ~/coral-build/upload/STAGED_SOURCES.md $HPC:~/CoralModelReview/
```

```bash
# fallback route only (HPC set as above):
rsync -avP --chmod=go-w ~/coral-build/upload/reefmod_host/ $HPC:~/CoralModelReview/reefmod_host/
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
| kora | 8 GB | 2 h | individual-based, but simulates one 250 m² section per site and extrapolates (~3 GB and ~1 s per scenario locally) |
| corset | 8 GB | 2 h | |
| reefmod | 16 GB | 12 h | one MATLAB at a time, including a MATLAB start-up per scenario |

These are generous first guesses. The run's own benchmarks (figure section K) show what each engine
actually used, so the next run can request less.

Variations (every `CMR_*` variable in your shell is passed on to the jobs):

```bash
bash jobs/submit_all.sh TestDomain_N8 --engines "kora corset"      # a subset
CMR_N_SAMPLES=16 bash jobs/submit_all.sh TestDomain_N16             # more scenarios
CMR_DHW_IDX=19 bash jobs/submit_all.sh TestDomain_N8_dhw19          # another climate member
CMR_NO_DISTURBANCE=true bash jobs/submit_all.sh TestDomain_N8_nodist  # no heat stress/cyclones, all six engines
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

**Re-run one engine** (for example after a failure), then the analysis, with the same `CMR_*`
settings the project was submitted with (domain, members, ...):

```bash
cd ~/CoralModelReview
bash jobs/submit_all.sh TestDomain_N8 --engines kora
```

`jobs.txt` keeps the other engines' job ids and replaces kora's, so the analysis (which waits for the
re-run) still collects everyone's Slurm accounting.

The analysis can be re-run on its own any time. It compares whichever engines have a result and lists
the missing ones in `FIGURES.md`:

```bash
cd ~/CoralModelReview/logs
sbatch --export=ALL,CMR_PROJECT_NAME=TestDomain_N8 ../jobs/analyse.slurm
```

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
HPC=vhaller@hpc-l001.aims.gov.au
DEST="/mnt/c/Users/vhaller/OneDrive - Australian Institute of Marine Science/DecisionSupport/CoralModelReview/ModelRuns/TestDomain_N8"
mkdir -p "$DEST"
rsync -avP --exclude='*/run/' --exclude='*/warmup/' $HPC:~/CoralModelReview/runs/TestDomain_N8/ "$DEST/"
```

Leaving out `run/` keeps the copy small. Every `<engine>/result.jld2` and `samples.jld2` comes back,
so `09_ANALYSE_SIX_ENGINES.jl` can be re-run on the laptop with `CMR_PROJECT_NAME=TestDomain_N8` (for
another domain, also set `CMR_DOMAIN_PATH` to a local copy of it: [9a](#9a-running-a-different-domain)). The
one exception is CoralBlox, whose ADRIA result store sits in `coralblox/run/`: add
`--include='coralblox/run/***'` before the excludes if you want to re-analyse locally.

---

## 8. Configuration reference

All scripts read `CMR_<NAME>` environment variables; defaults are in `engine_config()` in `engines.jl`.

| Variable | Default | Meaning |
|---|---|---|
| `CMR_PROJECT_NAME` | `TestDomain_N8` | output folder under `runs/` (set by `submit_all.sh`) |
| `CMR_MODEL_RUNS_DIR` | `~/CoralModelReview/runs` (jobs) / `../ModelRuns` (laptop) | where projects are written |
| `CMR_DOMAIN_PATH` | `/opt/coral/data/Test_domain` (image) | ADRIA domain; the folder must keep its own name. Another domain needs no rebuild ([9a](#9a-running-a-different-domain)) |
| `CMR_RCP` | `45` | |
| `CMR_DHW_IDX` / `CMR_WAVE_IDX` / `CMR_CYCLONE_IDX` | `40` / `16` / `1` | the climate member all engines share (Test domain values; each domain has its own, see 9a) |
| `CMR_NO_DISTURBANCE` | `false` | `true` = no heat stress and no cyclones for every engine (DHW and cyclone member 0); the wave member stays, as it only sets C~scape's site wave exposure. Use its own project name. ReefMod gets an all-zero DHW (it only bleaches at DHW >= 3) and CORSET no bleaching years; images built before 2026-10-08 stop ReefMod with "Needs dhw_idx > 0", so there add `--engines "coralblox cscape coralpop kora corset"` |
| `CMR_LARVAE_VIA_CONNECTIVITY` | `false` | `true` = CoralPop and Kora take all larvae from the domain's connectivity matrix (a site's own via the diagonal), like CoralBlox and C~scape, instead of CoralPop's own retention / Kora's self-seeding. Use its own project name |
| `CMR_N_SAMPLES` | `8` | scenarios per engine (power of 2) |
| `CMR_SEED` | `1234` | seed for ADRIA's scrambled Sobol' sampler |
| `CMR_N_WORKERS` | `0` | parallelism; its meaning differs per engine (CoralPop/C~scape: Distributed workers; Kora/CORSET: thread cap, needs more CPUs; ReefMod: concurrent MATLABs) |
| `CMR_CSCAPE_THREADS` | `1` | threads per C~scape worker |
| `CMR_WARMUP` | `true` | run 1 scenario first so time to first result and steady-state cost are measured separately |
| `CMR_CALIBRATION_PATH` | `/opt/coral/data/calibrated_params.nc` (baked in) | CoralBlox calibration NetCDF. Set it to the empty string to run on ADRIA's built-in parameters instead (every run before 2026-10-12 did, by mistake); a path that does not exist stops the CoralBlox job (the other engines don't use the file) |
| `CMR_CORSET_USE_DHW_AS_BLEACHING` | `true` | CORSET's heat stress: each year whose across-site mean DHW > 4 is one bleaching event at every site, a random 0-50 % coral loss per site (CORSET.jl's placeholder, not scaled by DHW). DHW member 0 or an all-zero member gives none. `false` = CORSET without heat stress (the default before 2026-10-08) |
| `CMR_CORSET_USE_CYCLONES_AS_HURRICANES` | `false` | `true` = years whose mean cyclone mortality > 0.1 become CORSET hurricane events |
| `CMR_CORSET_INITIAL_*_BIOMASS` | see `engines.jl` | CORSET's placeholder consumer biomass |
| `MATLAB_MODULE` | `matlab` | the module the ReefMod jobs load |
| `REEFMOD_ARRAY_CONCURRENCY` | `4` | fallback route: simultaneous MATLAB tasks |
| `CMR_BASE`, `SIF`, `RUNS` | `~/CoralModelReview`, `$CMR_BASE/coral_engines.sif`, `$CMR_BASE/runs` | locations (`jobs/common.sh`) |

---

## 9. Other domains and the full run

### 9a. Running a different domain

The domain is a setting, not part of the build, so a new domain needs **no rebuild**. `CMR_DOMAIN_PATH`
points the jobs at any ADRIA domain folder (the Test domain inside the image is only the default).
Upload the folder, check that the container can see it, and submit with that domain's own climate
members. The example uses the Moore domain.

1. **Upload** (WSL). Leave off the trailing `/` after the folder name so it keeps its name on the HPC:
   ```bash
   HPC=vhaller@hpc-l001.aims.gov.au
   ssh $HPC "mkdir -p ~/CoralModelReview/domains"
   rsync -avP --chmod=go-w "/mnt/c/Users/vhaller/OneDrive - Australian Institute of Marine Science/DecisionSupport/CoralModelReview/domains/SingleReef"
   ```
   Don't rename the folder: every engine records the domain by its folder name, and `09` refuses to
   compare engines whose domain names differ. `--chmod=go-w` is there because WSL shows files under
   `/mnt/c` as writable by everyone, and `rsync -a` would carry that over to the shared cluster.
2. **Check the container can see it** (login node):
   ```bash
   module load singularity/singularity.module
   singularity exec ~/CoralModelReview/coral_engines.sif ls $HOME/CoralModelReview/domains/SingleReef
   ```
   It should list `connectivity`, `DHWs`, `datapackage.json` and so on. Singularity normally mounts
   your home folder inside the container, but the job scripts only bind `runs/` explicitly. If you get
   "No such file", run `export SINGULARITY_BIND=$HOME/CoralModelReview/domains` before submitting; the
   jobs inherit it.
3. **Submit** with the domain and no disturbance run:
   ```bash
   cd ~/CoralModelReview
   CMR_DOMAIN_PATH=$HOME/CoralModelReview/domains/SingleReef \
     CMR_N_SAMPLES=8 CMR_DHW_IDX=1 CMR_WAVE_IDX=1 CMR_CYCLONE_IDX=1 bash jobs/submit_all.sh SingleReef_noDist
   ```

   **Submit** with the domain and its climate members:
   ```bash
   cd ~/CoralModelReview
   CMR_DOMAIN_PATH=$HOME/CoralModelReview/domains/SingleReef \
     CMR_N_SAMPLES=8 CMR_DHW_IDX=2 CMR_WAVE_IDX=1 CMR_CYCLONE_IDX=1 bash jobs/submit_all.sh SingleReef_1
   ```

Points to watch:

- **The climate members belong to the domain.** The defaults (DHW 40 / wave 16) are the Test domain's
  median-severity members. For Moore, `06_RUNNING_A_NEW_PROJECT.jl`'s notes use DHW 64 (its median
  severity), wave 1 (its only wave member) and cyclone 1. For another domain, choose them the same way
  (median severity: see the computation in `engines.jl`'s `engine_config` docstring). A DHW member the
  domain doesn't have stops the run with an error.
- **Use `$HOME`, not `~`,** in `CMR_DOMAIN_PATH`, so it is a full path everywhere.
- **The job sizes are set for the 10-site Test domain.** The memory and time limits in `resources()`
  (`jobs/submit_all.sh`) will probably be too small for Moore, especially for C~scape and ReefMod.
  Raise them, or first do a sizing run (`CMR_N_SAMPLES=2`, project `Moore_N2`) and size the real run
  from its benchmarks ([9b](#9b-scaling-up-to-the-full-run)).
- **Re-running the analysis by hand** (6d) needs nothing extra. `09` reads the domain, RCP and climate
  members from the saved results, and finds the domain at the path the engines recorded. To
  re-analyse on the laptop after copying the results back ([7](#7-retrieve-the-results)), set
  `CMR_DOMAIN_PATH` to a local copy of the same domain folder (for example
  `C:\Users\vhaller\Downloads\Moore_2025-08-22_v080`).

Rebuilding only makes sense to bake a domain into the image so it is self-contained: re-stage with
`DOMAIN_SRC=...` and rebuild (section 4, ~50 min). Even then you set
`CMR_DOMAIN_PATH=/opt/coral/data/<folder>`, because only `Test_domain` is the default.

### 9b. Scaling up to the full run

1. Read `benchmark_summary.csv` and figure section K from the test run: seconds per scenario, peak
   memory and disk per scenario, per engine.
2. Choose `N_SAMPLES`. PAWN needs N/S > 80, i.e. **N ≥ 1024** at S = 10. Below that, the sensitivity
   figures are labelled illustrative.
3. Size each engine's job as **time ≈ N × seconds-per-scenario ÷ parallelism**, with memory from peak
   memory × workers. Raise `--cpus-per-task` and `CMR_N_WORKERS` together (edit `resources()` in
   `submit_all.sh`). ReefMod's parallelism is capped by MATLAB licences.
4. For a real domain, set it up as in [9a](#9a-running-a-different-domain) (its own climate members
   included), and remember the domain's own `WAVE_IDX`/`ub_med` caveats (06's notes).

### 9c. Several domains and climate members at once

`jobs/submit_batch.sh` runs `submit_all.sh` once per domain × DHW member. Each pair is its own project,
named `<domain>/DHW<member>`, so the results are grouped by domain: `runs/SingleReef/DHW1/`,
`runs/SingleReef/DHW2/`, and so on. Upload the domains first (9a step 1), check the list with
`--dry-run`, then submit the same command without it:

```bash
cd ~/CoralModelReview
CMR_N_SAMPLES=8 CMR_WAVE_IDX=1 CMR_CYCLONE_IDX=1 bash jobs/submit_batch.sh \
  --domains "SingleReef ThreeReefsLine ThreeReefsRing NineReefsLine NineReefsRing" --dhw "1 2 17" --dry-run
```

Check the completion of the jobs

```bash
cd ~/CoralModelReview && bash jobs/batch_status.sh
```

Download the folders

```bash
rsync -avP --include='coralblox/run/***' --exclude='*/run/' --exclude='*/warmup/' $HPC:~/CoralModelReview/runs/{SingleReef,ThreeReefsLine,ThreeReefsRing,NineReefsLine,NineReefsRing} "/mnt/c/Users/vhaller/OneDrive - Australian Institute of Marine Science/DecisionSupport/CoralModelReview/ModelRuns/"
```

- Every listed domain gets every listed DHW member. For different members per domain, run it once per
  group of domains.
- Every `CMR_*` variable in your shell applies to all projects; `--engines "..."` and the other
  `submit_all.sh` options are passed on.
- The synthetic domains (SyntheticReefDomain.jl) have 51 DHW members: 1 is all zero (no heat stress),
  2-51 are random draws. Their median-severity members: SingleReef 17, ThreeReefsLine/Ring 46,
  NineReefsLine/Ring 31. Their cyclone data is all zero.
- Each project is 7 jobs. Every ReefMod job holds one MATLAB licence while it runs, so a large batch can
  start more MATLABs than there are licence seats (the ReefMod runs then fail with a licence error).
  If that happens, leave ReefMod out with `--engines "coralblox cscape coralpop kora corset"` and
  submit it in smaller batches with `--engines reefmod --no-analyse`, then re-run the analysis (6d).
- Retrieve a whole domain at once (7's command with the domain folder as the source):
  `rsync -avP --exclude='*/run/' --exclude='*/warmup/' $HPC:~/CoralModelReview/runs/SingleReef/ "$DEST/SingleReef/"`.
  To re-analyse one on the laptop, set `CMR_PROJECT_NAME=SingleReef/DHW2`.

---

## 10. Troubleshooting

| Symptom | Cause / fix |
|---|---|
| `sbatch: error: ... account` | Add `#SBATCH --account=...` to the job scripts. None was needed on `hpc-l001` so far. |
| `ERROR: image not found` | Upload the `.sif`, or set `SIF=`. |
| ReefMod: `matlab exited with code ...` mentioning missing `lib*.so` | A MATLAB system library is missing from the image. Add it to `coral_engines.def`'s apt list and rebuild, or use the fallback (6e). `container ldd $MATLAB_ROOT/bin/glnxa64/MATLAB \| grep "not found"` (after `source jobs/common.sh; setup_matlab`) lists what is missing. |
| ReefMod: licence error / "Cannot find a valid license" | Is the licence server reachable from compute nodes? Check `MLM_LICENSE_FILE`/`LM_LICENSE_FILE` after `module load matlab` (preflight prints them). They pass into the container unchanged; a licence *file* outside the MATLAB folder is bound automatically by `setup_matlab`. |
| Trying a ReefMod-GBR fix without rebuilding | Re-stage (or rsync the clone into `~/coral-build/upload/reefmod_host/ReefMod-GBR/`), upload `reefmod_host/` (5), and prefix the command with `REEFMOD_GBR_DIR=$HOME/CoralModelReview/reefmod_host/ReefMod-GBR SINGULARITY_BIND=$HOME/CoralModelReview/reefmod_host`. Only the MATLAB side changes; rebuild before the real run so the image records what ran. |
| ReefMod: "Index in position 3 exceeds array bounds. Index must not exceed 100" in `settings_GBR_NEW` | The image's ReefMod-GBR lacks `LOCAL_CHANGES.md` item 7 (ADRIA's seed used as the cyclone-series index). Re-stage and rebuild, or use the row above. |
| ReefMod: "Unrecognized field name 'doing_cooling'" | The image was built from a ReefMod-GBR without commit `9be6b5c`. Re-stage and rebuild. |
| `Read-only file system` writing under `/opt/...` | A script wrote outside `CMR_MODEL_RUNS_DIR`. Check that `CMR_MODEL_RUNS_DIR` points into a bound folder (the jobs set it). |
| Julia recompiles at job start, or fails writing `~/.julia_coral` | The HPC CPU isn't covered by the image's `JULIA_CPU_TARGET`, or the home quota is full. A one-off recompile into `~/.julia_coral` is harmless. |
| `OutOfMemory` / job `OUT_OF_MEMORY` in `sacct` | Raise `--mem` for that engine in `submit_all.sh`. C~scape needs the most on bigger domains; Kora needs ~3 GB on the Test domain (more if its section size is raised). |
| `IOError: rm(".../warmup/...): directory not empty (ENOTEMPTY)`, mostly CoralBlox | The warm-up's result files were still open when `08` deleted them; on the HPC's NFS home they survive as `.nfs*` files. Fixed in `08` from 2026-10-08 (it collects them first and never fails over this). With an older image, re-run that engine with the fixed `08_RUN_ENGINE.jl` bound over the image's: upload it to `~/CoralModelReview/patches/`, `export SINGULARITY_BIND=$HOME/CoralModelReview/patches/08_RUN_ENGINE.jl:/opt/coral/Code/08_RUN_ENGINE.jl`, then `bash jobs/submit_all.sh <project> --engines coralblox` with the project's settings. |
| `$'\r': command not found` | A script has Windows line endings. Re-stage (it strips them), or run `sed -i 's/\r$//' jobs/*`. |
| 09: "The engines ran on the domain `...`, but ... is a different folder" | 09 can't find the domain the engines ran on, usually after copying results to another machine. Set `CMR_DOMAIN_PATH` to a copy of that domain's folder (same folder name). |
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
| `jobs/submit_batch.sh` | HPC login | `submit_all.sh` for every domain × DHW member, as projects `<domain>/DHW<member>` (9c) |
| `jobs/run_engine.slurm` | HPC compute | one engine (`ENGINE=...`) |
| `jobs/analyse.slurm` | HPC compute | `collect_benchmarks.sh` then `09` |
| `jobs/collect_benchmarks.sh` | HPC | `sacct` for the project's jobs into `slurm_accounting.csv` |
| `jobs/reefmod_A/B/C.slurm` | HPC compute | the fallback ReefMod route |
