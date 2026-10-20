#=
04_collect.jl -- merge the chunk files of each band and domain into one results file.

    julia +release --project=Code/sa/shapley Code/sa_pf/04_collect.jl <model> <band> [band ...]

Writes <band>/<domain>/pf_results.jld2: `Y` (Dict output => [n_base, d + 1, years], walks in base-sample
order), `dhw` (domain-mean DHW, same shape), `base`, `failed`, `years`, `seconds` (CPU seconds of the
runs), `n_runs`, `settings`. A walk with a failed run is NaN and left out of the analysis.
=#

length(ARGS) >= 2 || error("usage: 04_collect.jl <model> <band> [band ...]")
const MODEL = ARGS[1]
const BANDS = parse.(Float64, ARGS[2:end])
include(joinpath(@__DIR__, "pf_common.jl"))

for band in BANDS, domain in SA_DOMAINS, set in PF_SETS
    dir = pf_set_dir(MODEL, band, domain, set)
    cdir = joinpath(dir, "chunks")
    files = isdir(cdir) ? filter(f -> endswith(f, ".jld2"), readdir(cdir; join=true)) : String[]
    isempty(files) && (set == "main" && println("$(band_label(band))/$domain: no chunks"); continue)
    parts = [JLD2.load(f) for f in files]
    design = JLD2.load(pf_design_path(MODEL, band, set))
    N, d = design["n_base"], ncol(design["A"])
    years = parts[1]["years"]
    T = length(years)
    keys_ = collect(keys(parts[1]["Y"]))
    Y = Dict(k => fill(NaN, N, d + 1, T) for k in keys_)
    dhw = fill(NaN, N, d + 1, T)
    have = falses(N)
    for p in parts, (j, n) in enumerate(p["base"])
        for k in keys_
            Y[k][n, :, :] = p["Y"][k][:, :, j]
        end
        dhw[n, :, :] = p["dhw"][:, :, j]
        have[n] = true
    end
    failed = sort(unique(vcat([p["failed"] for p in parts]...)))
    missing_ = findall(.!have)
    seconds = sum(p["seconds"] for p in parts)
    n_runs = sum(p["n_runs"] for p in parts)
    JLD2.jldsave(joinpath(dir, "pf_results.jld2"); Y, dhw, base=collect(1:N), failed, missing=missing_, years,
        seconds, n_runs, settings=Dict("model" => MODEL, "band" => band, "domain" => domain, "set" => set, "n_base" => N,
            "d" => d, "names" => design["names"], "n_shared" => design["n_shared"]))
    println("$(band_label(band))$(set_suffix(set))/$domain: $(N - length(missing_)) of $N walks, $(length(failed)) failed, ",
        "$n_runs runs, $(round(seconds / 3600; digits=2)) CPU-h ($(round(seconds / n_runs; digits=3)) s/run)")

    ## zero-DHW reference (03b_run_zero.jl), if run: pf_zero.jld2 with Y = Dict output => [n_base, 2, years]
    zdir = joinpath(dir, "zero_chunks")
    zfiles = isdir(zdir) ? filter(f -> endswith(f, ".jld2"), readdir(zdir; join=true)) : String[]
    isempty(zfiles) && continue
    zparts = [JLD2.load(f) for f in zfiles]
    Z = Dict(k => fill(NaN, N, 2, T) for k in keys(zparts[1]["Y"]))
    for p in zparts, (j, n) in enumerate(p["base"]), k in keys(Z)
        Z[k][n, :, :] = p["Y"][k][:, :, j]
    end
    zfailed = sort(unique(vcat([p["failed"] for p in zparts]...)))
    JLD2.jldsave(joinpath(dir, "pf_zero.jld2"); Y=Z, failed=zfailed, years,
        seconds=sum(p["seconds"] for p in zparts), n_runs=sum(p["n_runs"] for p in zparts))
    println("  zero-DHW reference: $(sum(p["n_runs"] for p in zparts)) runs, $(length(zfailed)) failed")
end
