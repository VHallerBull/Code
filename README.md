# CoralModelReview / Code

Scripts and the Julia environment for running and comparing the coral models under review
(CoralBlox, C~scape, CoralPop) on a shared ADRIA domain.

```
CoralModelReview/
    Code/          <- this folder (under git): scripts + Julia environment
    ModelRuns/     <- everything the scripts write, one sub-folder per project (NOT under git)
```

| File | What it is |
|---|---|
| `06_RUNNING_A_NEW_PROJECT.jl` / `.md` | CoralBlox vs C~scape from a data folder to plots. Produces the C~scape results 07 reuses. |
| `07_COMPARE_THREE_ENGINES.jl` | CoralBlox vs C~scape vs CoralPop: same domain and climate member, ADRIA's sampler and metrics for all three, PAWN per engine. |
| `paths.jl` | Shared locations (`MODEL_RUNS_DIR`, `model_run_dir(name)`); `include` it from every script. |
| `Project.toml`, `Manifest.toml` | The Julia environment. |
| `setup_environment.jl` | Checks the local package clones exist, instantiates the environment. |
| `record_sources.jl`, `SOURCES.md` | Which branch/commit of each local clone the environment was verified against. |

## Running things

```bash
cd Code
julia --project=. setup_environment.jl            # once, or after a fresh clone
julia --project=. 06_RUNNING_A_NEW_PROJECT.jl     # slow: runs C~scape (needs R/brms only if demography isn't already converted)
julia --project=. 07_COMPARE_THREE_ENGINES.jl     # loads 06's C~scape result; runs CoralBlox and CoralPop
```

Each script has a `CONFIG` section at the top; edit it rather than the code below. Outputs go to
`../ModelRuns/<project name>/` (`images/`, saved result sets, sampled factor tables, `summary.csv`). Any
`CONFIG` constant that reads `_cfg("NAME", ...)`/`get(ENV, "CMR_NAME", ...)` can be overridden with an environment
variable `CMR_NAME` for a quick test run, e.g. `CMR_N_SAMPLES=2 CMR_PROJECT_NAME=scratch`.

**Order:** run 06 before 07 (07 reads `ModelRuns/<06's project>/my_run.jld2` and `samples.jld2`). If those are
missing, 07 still runs and compares CoralBlox with CoralPop, and says so. 07 refuses to compare if 06 used a
different domain, RCP, climate member or years.

## The environment and its sources

`Project.toml` lists the ADRIA, CScape, CscapeInterface and CoralPopInterface packages as **local clones** under
`[sources]` (paths relative to this folder, currently `../../../../Documents/GitHub/...`). The versions of everything
else were seeded from `CscapeInterface.jl`'s working manifest, so this environment resolves to the same packages the
earlier scripts ran with.

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
- `Manifest.toml` is committed for reproducibility of the registered packages; if it becomes noisy, regenerate it with
  `julia --project=. -e 'using Pkg; Pkg.update()'` and commit the result deliberately.
