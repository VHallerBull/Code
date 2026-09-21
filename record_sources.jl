#=
record_sources.jl -- (re)write SOURCES.md: the branch and commit of every local package clone this
environment was last verified against.

    julia --project=. record_sources.jl

The packages are local development clones, not registered releases, so Manifest.toml alone cannot
reproduce a run: it records paths, not the code at those paths. Run this after you have checked
that the scripts work with the clones' current state, and commit the result together with the run.
Clones with uncommitted changes are flagged -- results produced from them are not reproducible from git.
=#

using Pkg, Dates

const CODE = @__DIR__
project = Pkg.TOML.parsefile(joinpath(CODE, "Project.toml"))

git(dir, args...) = try
    strip(read(Cmd(`git -C $dir $(collect(args))`), String))
catch
    ""
end

## several packages can live in one repo (ADRIA, ADRIAanalysis, ADRIAviz): group by repository
repos = Dict{String,Vector{String}}()
for (name, spec) in project["sources"]
    dir = normpath(joinpath(CODE, spec["path"]))
    top = git(dir, "rev-parse", "--show-toplevel")
    push!(get!(repos, isempty(top) ? dir : normpath(top), String[]), name)
end

io = IOBuffer()
println(io, "# Source clones")
println(io)
println(io, "The Julia packages this environment uses are local development clones (see `Project.toml` `[sources]`).")
println(io, "Recorded by `record_sources.jl` on ", Dates.format(now(), "yyyy-mm-dd HH:MM"), ".")
println(io)
println(io, "| Repository (local path) | Packages used | Branch | Commit | Uncommitted changes |")
println(io, "|---|---|---|---|---|")
for repo in sort(collect(keys(repos)))
    branch = git(repo, "branch", "--show-current")
    commit = git(repo, "log", "-1", "--format=%h %s")
    dirty = isempty(git(repo, "status", "--porcelain", "--untracked-files=no")) ? "no" : "**yes**"
    println(io, "| `", repo, "` | ", join(sort(repos[repo]), ", "), " | ", branch, " | ", commit, " | ", dirty, " |")
end
println(io)
println(io, "Other local inputs the scripts read (not packages):")
println(io, "- the ADRIA domain package and CoralBlox calibration file named in each script's CONFIG section")
println(io, "- the CoralPop Python engine, `CoralPop-Engine` (branch `engine-api`), used through CoralPopInterface")
engine = normpath(joinpath(CODE, "..", "..", "..", "..", "Documents", "GitHub", "CoralPop-Engine"))
if isdir(engine)
    println(io, "  - ", engine, ": ", git(engine, "branch", "--show-current"), " @ ", git(engine, "log", "-1", "--format=%h %s"),
        isempty(git(engine, "status", "--porcelain", "--untracked-files=no")) ? "" : " **(uncommitted changes)**")
end
write(joinpath(CODE, "SOURCES.md"), String(take!(io)))
println("Wrote SOURCES.md")
