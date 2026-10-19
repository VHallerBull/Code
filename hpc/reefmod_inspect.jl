#=
reefmod_inspect.jl <engine dir> -- check a REAL ReefMod run against what ReefModInterface.jl assumes.

    julia reefmod_inspect.jl ~/CoralModelReview/runs/ReefModCheck/reefmod

Run by jobs/reefmod_check.sh after the first real MATLAB run. Covers the items ReefMod-GBR's
LOCAL_CHANGES.md lists under "Still to verify once a real run exists", and prints a verdict per item:

1. RESULT.coral_pct2D's shape is [reefs, T6+1, 6] and its units are PERCENT (ReefModInterface divides
   by 100): the first time step should reproduce the initial cover handed to ReefMod.
2. The external DHW reached ReefMod with the right shape (the driver errors otherwise; reaching the
   output at all means the shape check passed) and the run covers the expected number of 6-monthly steps.
3. 6-monthly -> annual indexing: the result set's first year equals the raw initial state.
4. Juvenile counts and size-frequency counts are present (needed for juveniles / size structure).
=#

include(joinpath(@__DIR__, "..", "engines.jl"))

dir = isempty(ARGS) ? error("usage: julia reefmod_inspect.jl <ModelRuns/<project>/reefmod>") : ARGS[1]
cfg = engine_config()
dom, _, YEARS = load_shared_domain(cfg)
inp = domain_to_reefmod_inputs(dom; dhw_idx=cfg.dhw_idx)
outs = reefmod_output_paths(joinpath(dir, "run"))
isempty(outs) && error("No ReefMod_scenario_*.mat under $(joinpath(dir, "run", "outputs"))")
verdict(ok, msg) = println(ok ? "  [ok]    " : "  [CHECK] ", msg)

for p in outs
    println("\n== ", basename(p))
    raw = MAT.matread(p)
    R, META = raw["RESULT"], get(raw, "META", Dict{String,Any}())
    println("  RESULT fields: ", join(sort(collect(keys(R))), ", "))
    cov = Float64.(R["coral_pct2D"])
    nts = Int(get(META, "nb_time_steps", size(cov, 2) - 1))
    println("  coral_pct2D size ", size(cov), "; nb_time_steps ", nts, "; range ", extrema(cov))

    ## 1. shape + units
    verdict(size(cov, 1) == length(inp.site_ids) && size(cov, 3) == 6,
        "shape [reefs=$(size(cov, 1)) (domain $(length(inp.site_ids))), T+1=$(size(cov, 2)), groups=$(size(cov, 3))]")
    init_given = vec(sum(inp.init_cover; dims=1))           # proportion of reef area handed to ReefMod
    init_out = vec(sum(cov[:, 1, :]; dims=2))                # ReefMod's own first step
    ratio = init_out ./ max.(init_given, 1e-12)
    println("  initial total cover per reef, given (proportion): ", round.(init_given; digits=3))
    println("  initial total cover per reef, ReefMod output:      ", round.(init_out; digits=3))
    verdict(all(r -> 80 < r < 120, ratio[init_given .> 0.01]),
        "coral_pct2D is in PERCENT (output/given ~ 100: $(round.(ratio; digits=1))) -- ReefModInterface divides by 100")

    ## 2. DHW + steps
    verdict(nts == 2 * length(YEARS), "6-monthly steps $nts == 2 x $(length(YEARS)) ADRIA years (external_DHW shape check passed: the run completed)")

    ## 3. annual indexing
    idx = reefmod_annual_indices(nts)
    println("  annual indices into the T+1 axis (first 5, last): ", idx[1:min(5, end)], " ... ", idx[end])
    verdict(idx[1] == 1, "the first ADRIA year is ReefMod's initial state (index 1)")

    ## 4. juveniles + size frequency
    for k in ("coral_juv_count", "coral_adol_count", "coral_adult_count")
        haskey(R, k) ? println("  $k size ", size(R[k]), ", total ", sum(R[k])) : println("  $k MISSING")
    end
    verdict(haskey(R, "coral_juv_count"), "juvenile counts present (relative_juveniles can be computed)")
    println("  grid area (m^2): ", get(META, "total_area_cm2", NaN) / 1e4)
end

## the result set 08 saved, through the same metric function the other engines use
if has_engine_result(:reefmod, dir)
    rs = load_engine_result(:reefmod, dir, YEARS)
    rc = Array(ADRIA.metrics.relative_cover(rs))
    println("\n== saved ReefModResultSet: relative_cover size ", size(rc), ", first-year domain mean ",
        round(mean(rc[1, :, :]); digits=4), ", final ", round(mean(rc[end, :, :]); digits=4))
    verdict(all(isfinite, rc) && maximum(rc) <= 1.0 + 1e-6, "relative cover finite and <= 1 (fraction of k)")
end
println("\nIf every line is [ok], record that in ReefMod-GBR/LOCAL_CHANGES.md ('Still to verify') and run the batch.")
