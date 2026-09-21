#=
paths.jl -- where things live, shared by every script in this folder.

    include(joinpath(@__DIR__, "paths.jl"))

CoralModelReview/
    Code/         <- this folder: scripts + the Julia environment; under git
    ModelRuns/    <- everything the scripts WRITE (one sub-folder per project); NOT under git
=#

"The folder holding these scripts and the Julia environment (`Project.toml`/`Manifest.toml`)."
const CODE_DIR = @__DIR__

"""
Root output folder for model runs: `CoralModelReview/ModelRuns`, next to `Code/`. Each script writes to its
own sub-folder (`model_run_dir("SomeProject")`). Kept outside `Code/` so large outputs never end up in git.
"""
const MODEL_RUNS_DIR = normpath(joinpath(CODE_DIR, "..", "ModelRuns"))

"""
    model_run_dir(name) -> String

`MODEL_RUNS_DIR/name`, created if it does not exist.
"""
function model_run_dir(name::AbstractString)::String
    dir = joinpath(MODEL_RUNS_DIR, name)
    mkpath(dir)
    return dir
end
