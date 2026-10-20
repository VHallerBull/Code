#=
04_collect.jl -- merge the chunk files of each band and domain into one results file.

    julia +release --project=Code/sa/shapley Code/sa_pf/04_collect.jl <model> <band> [band ...]

Writes <band>/<domain>/pf_results.jld2:
    `state`  the runs' state summaries (sa_outputs.jl): Dict key => [n_base, d + 1, dims...]
    `Y`      the outputs derived from them: Dict output => [n_base, d + 1, years] (sa_outputs.jl)
    `dhw`    domain-mean DHW, [n_base, d + 1, years]
    `base`, `failed`, `missing`, `years`, `seconds` (CPU seconds of the runs), `n_runs`,
    `thermal_groups`, `settings`.
Walks in base-sample order; a walk with a failed run is NaN and left out of the analysis. The outputs can
be re-derived from `state` with a changed definition (`derive_stored`) without rerunning anything.

If the zero-DHW reference was run (03b_run_zero.jl), also pf_zero.jld2 with `state` and `Y`
([n_base, 2, ...], the A and B parameter sets without heat stress).
=#

length(ARGS) >= 2 || error("usage: 04_collect.jl <model> <band> [band ...]")
const MODEL = ARGS[1]
const BANDS = parse.(Float64, ARGS[2:end])
include(joinpath(@__DIR__, "pf_common.jl"))

chunk_files(dir) = isdir(dir) ? filter(f -> endswith(f, ".jld2"), readdir(dir; join=true)) : String[]

for band in BANDS, domain in SA_DOMAINS
    dir = pf_run_dir(MODEL, band, domain)
    files = chunk_files(joinpath(dir, "chunks"))
    isempty(files) && (println("$(band_label(band))/$domain: no chunks"); continue)
    parts = [JLD2.load(f) for f in files]
    design = JLD2.load(joinpath(pf_band_dir(MODEL, band), "design_pf.jld2"))
    N, d = design["n_base"], ncol(design["A"])
    years = parts[1]["years"]
    T = length(years)
    thermal_groups = parts[1]["thermal_groups"]
    state = merge_states(parts, N, d + 1)
    dhw = fill(NaN, N, d + 1, T)
    have = falses(N)
    for p in parts, (j, n) in enumerate(p["base"])
        dhw[n, :, :] = p["dhw"][:, :, j]
        have[n] = true
    end
    failed = sort(unique(vcat([p["failed"] for p in parts]...)))
    missing_ = findall(.!have)
    seconds = sum(p["seconds"] for p in parts)
    n_runs = sum(p["n_runs"] for p in parts)
    Y = derive_stored(state; thermal_groups)
    JLD2.jldsave(joinpath(dir, "pf_results.jld2"); state, Y, dhw, base=collect(1:N), failed, missing=missing_,
        years, seconds, n_runs, thermal_groups, settings=Dict("model" => MODEL, "band" => band, "domain" => domain,
            "n_base" => N, "d" => d, "names" => design["names"], "n_shared" => design["n_shared"]))
    println("$(band_label(band))/$domain: $(N - length(missing_)) of $N walks, $(length(failed)) failed, ",
        "$n_runs runs, $(round(seconds / 3600; digits=2)) CPU-h ($(round(seconds / n_runs; digits=3)) s/run); ",
        "outputs: ", join(output_order(keys(Y); thermal_groups), ", "))

    ## zero-DHW reference (03b_run_zero.jl), if run
    zfiles = chunk_files(joinpath(dir, "zero_chunks"))
    isempty(zfiles) && continue
    zparts = [JLD2.load(f) for f in zfiles]
    zstate = merge_states(zparts, N, 2)
    zfailed = sort(unique(vcat([p["failed"] for p in zparts]...)))
    JLD2.jldsave(joinpath(dir, "pf_zero.jld2"); state=zstate, Y=derive_stored(zstate; thermal_groups), failed=zfailed,
        years, thermal_groups, seconds=sum(p["seconds"] for p in zparts), n_runs=sum(p["n_runs"] for p in zparts))
    println("  zero-DHW reference: $(sum(p["n_runs"] for p in zparts)) runs, $(length(zfailed)) failed")
end
