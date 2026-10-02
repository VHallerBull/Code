# CoralModelReview / Code

Scripts and the Julia environment for running and comparing the coral models under review
(CoralBlox, C~scape, CoralPop, Kora, ReefMod, CORSET) on a shared ADRIA domain.

CORSET is a from-scratch reimplementation, not a port of the original model (CORSET's own source was
never published) -- see `CORSET.jl/CAVEATS.md` before treating its output as more than structurally
representative, and `CORSETInterface.jl/src/domain_translator.jl` for what it needs that an ADRIA
domain cannot supply (consumer biomass above all).

```
CoralModelReview/
    Code/          <- this folder (under git): scripts + Julia environment
    ModelRuns/     <- everything the scripts write, one sub-folder per project (NOT under git)
```

| File | What it is |
|---|---|
| `06_RUNNING_A_NEW_PROJECT.jl` / `.md` | CoralBlox vs C~scape from a data folder to plots. Produces the C~scape results 07 reuses. |
| `07_COMPARE_SIX_ENGINES.jl` | CoralBlox vs C~scape vs CoralPop vs Kora vs ReefMod vs CORSET: same domain and climate member, ADRIA's sampler and metrics for all six, PAWN per engine. |
| `engines.jl` | The per-engine code shared by 08 and 09 (factor specs, run calls, save/load, size structure, benchmark helpers), lifted from 06/07. |
| `08_RUN_ENGINE.jl <engine>` | Runs ONE engine into `ModelRuns/<project>/<engine>/` with benchmarks -- one HPC job per engine. Defaults: ADRIA Test domain, DHW/wave member 40/16, 8 scenarios. |
| `09_ANALYSE_SIX_ENGINES.jl` | The full comparison of whatever 08 produced: figures A-K (sanity, domain metrics, taxa, size structure, disturbance response, spatial, agreement, sensitivity, CORSET extras, benchmarks), `summary.csv`, `FIGURES.md`. |
| `hpc/` | The Apptainer image and Slurm jobs for running all six engines on the AIMS HPC -- start with **`hpc/RUNNING_ON_HPC.md`**. |
| `paths.jl` | Shared locations (`MODEL_RUNS_DIR`, `model_run_dir(name)`); `include` it from every script. |
| `Project.toml`, `Manifest.toml` | The Julia environment. |
| `setup_environment.jl` | Checks the local package clones exist, instantiates the environment. |
| `record_sources.jl`, `SOURCES.md` | Which branch/commit of each local clone the environment was verified against. |

## Running things

```bash
cd Code
julia --project=. setup_environment.jl            # once, or after a fresh clone
julia --project=. 06_RUNNING_A_NEW_PROJECT.jl     # slow: runs C~scape (needs R/brms only if demography isn't already converted)
julia --project=. 07_COMPARE_SIX_ENGINES.jl       # loads 06's C~scape result and the ReefMod HPC pipeline's
                                                   # result; runs CoralBlox, CoralPop, Kora and CORSET here
```

**On the HPC (all six engines, including ReefMod):** `08`/`09` inside one Apptainer image, one Slurm job per
engine plus an analysis job -- see `hpc/RUNNING_ON_HPC.md` for building, uploading and running, in order.
Locally, `08`/`09` run the same way for every engine except ReefMod (no MATLAB here):

```bash
julia +release --project=. 08_RUN_ENGINE.jl kora          # one engine (CMR_N_SAMPLES=2 for a quick test)
julia +release --project=. 09_ANALYSE_SIX_ENGINES.jl      # compare whatever has been run
```

ReefMod has no local MATLAB license, so it isn't run by this environment at all: its whole pipeline
(ADRIA sampling through to a saved result set) runs on the HPC via
`../../../Documents/GitHub/ReefModInterface.jl/jobs/submit_reefmod_pipeline.sh` (see that package's
README), and `07` just loads what it saved -- the same shape C~scape already uses.

Each script has a `CONFIG` section at the top; edit it rather than the code below. Outputs go to
`../ModelRuns/<project name>/` (`images/`, saved result sets, sampled factor tables, `summary.csv`). Any
`CONFIG` constant that reads `_cfg("NAME", ...)`/`get(ENV, "CMR_NAME", ...)` can be overridden with an environment
variable `CMR_NAME` for a quick test run, e.g. `CMR_N_SAMPLES=2 CMR_PROJECT_NAME=scratch`.

**Order:** run 06 (for C~scape) and the ReefMod HPC pipeline (for ReefMod) before 07 -- 07 reads
`ModelRuns/<06's project>/{my_run,samples}.jld2` and `ModelRuns/<REEFMOD_PROJECT_NAME>/{reefmod_result,job_meta}.jld2`.
If either is missing, 07 still runs and compares whichever engines it has, and says so. 07 refuses to
compare an engine that used a different domain, RCP, climate member or years.

## The environment and its sources

`Project.toml` lists the ADRIA, CScape, CscapeInterface, CoralPopInterface, Kora, KoraInterface,
CORSET, CORSETInterface and ReefModInterface packages as **local clones** under `[sources]` (paths relative to this folder,
currently `../../../../Documents/GitHub/...`). The versions of everything else were seeded from
`CscapeInterface.jl`'s working manifest, so this environment resolves to the same packages the earlier
scripts ran with. ReefMod's own MATLAB source (`ReefModInterface.jl`'s sibling `ReefMod-GBR` clone) is
not a Julia dependency and plays no part in this environment -- it only needs to exist, locally and on
the HPC, for `ReefModInterface.jl`'s own pipeline.

Because the model packages are development clones rather than releases, `Manifest.toml` records *where* they are, not
*which code* is there. `SOURCES.md` records the branch and commit of each clone, and flags any with uncommitted
changes (results from those cannot be reproduced from git). Re-run `julia --project=. record_sources.jl` after checking
that the scripts work with the clones' current state, and commit the result with the run you care about.

CoralPop additionally needs a Python interpreter with numpy and pandas and the `CoralPop-Engine` folder.
`07` sets `JULIA_PYTHONCALL_EXE`, `JULIA_CONDAPKG_BACKEND=Null` and `CORALPOP_ENGINE_DIR` **before** loading any
packages (PythonCall reads them at load time and Distributed workers inherit them); edit those lines at the top
of `07` if your paths differ.

## Notes

- This folder sits in OneDrive. A git repository's `.git` folder is synchronised like any other, and OneDrive can
  interfere with git's internal files if two machines work on the folder at once. Work on one machine at a time, or
  keep the repository's history somewhere else (`git remote add origin ...; git push`).
- `ModelRuns/` is outside the repository on purpose: outputs are large (hundreds of MB per C~scape project).
  It is also inside OneDrive, which will try to sync every file the scripts write. In testing this made deleting a
  freshly written output folder fail with "permission denied" while OneDrive held the files (07 now warns and
  carries on instead of stopping), and syncing gigabytes of model output is slow and uses quota. If that bothers you,
  pause OneDrive while a big run is writing, or set `ModelRuns` to "Free up space" / exclude it from sync.
- `Manifest.toml` is committed for reproducibility of the registered packages; if it becomes noisy, regenerate it with
  `julia --project=. -e 'using Pkg; Pkg.update()'` and commit the result deliberately.
