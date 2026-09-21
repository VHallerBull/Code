#=
setup_environment.jl -- one-time (or after a fresh clone) setup of this folder's Julia environment.

    julia --project=. setup_environment.jl

1. Checks that every local package clone listed under [sources] in Project.toml exists (they live
   outside this folder, under Documents/GitHub -- see README.md).
2. Instantiates the environment from Manifest.toml and precompiles it.
3. Prints the branch/commit each clone is on, next to what SOURCES.md recorded, so drift is visible.

Nothing is modified except Julia's own package depot.
=#

using Pkg

## PythonCall (a dependency of CoralPopInterface) must not try to build its own conda environment.
ENV["JULIA_CONDAPKG_BACKEND"] = "Null"

const CODE = @__DIR__
project = Pkg.TOML.parsefile(joinpath(CODE, "Project.toml"))

## --- 1. the local clones must be present ------------------------------------------------------
missing_clones = String[]
clones = Dict{String,String}()
for (name, spec) in project["sources"]
    dir = normpath(joinpath(CODE, spec["path"]))
    clones[name] = dir
    isfile(joinpath(dir, "Project.toml")) || push!(missing_clones, "  $name  ->  $dir")
end
isempty(missing_clones) || error(
    "These local package clones (Project.toml [sources]) were not found:\n" * join(missing_clones, "\n") *
    "\nClone them (see SOURCES.md for the branches) or edit the [sources] paths in Project.toml."
)
println("All ", length(clones), " local package clones found.")

## --- 2. instantiate ---------------------------------------------------------------------------
Pkg.activate(CODE)
Pkg.instantiate()
Pkg.precompile()

## --- 3. show where each clone is ---------------------------------------------------------------
git(dir, args...) = try
    strip(read(Cmd(`git -C $dir $(collect(args))`), String))
catch
    "?"
end
println("\nclone                       branch                commit    working tree")
for name in sort(collect(keys(clones)))
    repo = git(clones[name], "rev-parse", "--show-toplevel")
    repo == "?" && (println(rpad(name, 28), "(not a git checkout)"); continue)
    branch = git(repo, "branch", "--show-current")
    commit = git(repo, "rev-parse", "--short", "HEAD")
    dirty = isempty(git(repo, "status", "--porcelain")) ? "clean" : "uncommitted changes"
    println(rpad(name, 28), rpad(branch, 22), rpad(commit, 10), dirty)
end
println("\nCompare with SOURCES.md; run `julia --project=. record_sources.jl` to refresh it.")
