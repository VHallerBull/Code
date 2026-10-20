#=
01_dhw_draws.jl -- save the random draws behind each domain's DHW ensemble, so any run can generate its
own trajectory from 7 parameters + a weather index (option B of the pick-and-freeze design).

    julia +release --project=Code/sa Code/sa_pf/01_dhw_draws.jl

Reads domains/<domain>/generator_config.toml (read only), recreates the draws exactly as
SyntheticReefDomain's `build_dhw` does, and writes pf_inputs/<domain>/dhw_draws.jld2 (shared by every
model). Checks that `generate_dhw` with each member's parameters reproduces `build_dhw` exactly and the
domain's DHW file (as ADRIA loads it) to file precision.
=#

using SyntheticReefDomain, DataFrames, Statistics
import JLD2, ADRIA
include(joinpath(@__DIR__, "dhw_generator.jl"))

const CMR = normpath(joinpath(@__DIR__, "..", ".."))
const DOMAINS = ["SingleReef_SobolDHW", "NineReefs_SobolDHW"]
const RCP = 45

for domain in DOMAINS
    cfg = load_config(joinpath(CMR, "domains", domain, "generator_config.toml"))
    reefs = build_reefs(cfg)
    yrs = years(cfg)
    T, n, S = length(yrs), nrow(reefs), cfg.n_dhw_scenarios

    ## the draws, in build_dhw's order (common random numbers: one regional stream, one per reef)
    rng = SyntheticReefDomain.component_rng(cfg.seed, "dhw_regional")
    z = randn(rng, T, S)
    u_event = rand(rng, T, S)
    u_mag = rand(rng, T, S)
    local_noise = Array{Float64}(undef, T, n, S)
    for r in 1:n
        local_noise[:, r, :] .= randn(SyntheticReefDomain.component_rng(cfg.seed, "dhw_local_$(r)"), T, S)
    end
    p = scenario_params(cfg, RCP)
    prescribed = collect(in.(yrs, Ref(Set(p.event_years))))
    out = mkpath(joinpath(CMR, "pf_inputs", domain))
    path = joinpath(out, "dhw_draws.jld2")
    JLD2.jldsave(path; z, u_event, u_mag, local_noise, exposure=Float64.(reefs.dhw_exposure),
        offset=Float64.(reefs.dhw_offset), prescribed, years=collect(yrs), seed=cfg.seed)

    ## check 1: every member against build_dhw (same arithmetic: must be identical)
    D = load_dhw_draws(path)
    ref, _ = build_dhw(cfg, reefs)
    ref = ref[RCP]                                    # [T, n, S + 1]; member 1 = zero
    params(s) = (; dhw_baseline=p.baseline[s], dhw_trend_per_decade=p.trend_per_decade[s],
        dhw_interannual_sd=p.interannual_sd[s], dhw_ar1=p.ar1[s], dhw_event_prob=p.event_prob[s],
        dhw_event_mag_width=p.event_magnitude_max[s] - PF_DHW_EVENT_MAG_MIN, dhw_local_sd=p.local_sd[s])
    maxdiff = maximum(maximum(abs.(generate_dhw(params(s), s, D) .- ref[:, :, s + 1])) for s in 1:S)
    ## check 2: the domain's DHW file as ADRIA loads it
    file = Array(ADRIA.load_domain(joinpath(CMR, "domains", domain), string(RCP)).dhw_scens)
    filediff = maximum(maximum(abs.(generate_dhw(params(s), s, D) .- file[:, :, s + 1])) for s in 1:S)
    println("$domain: $T years x $n reefs x $S members -> $path")
    println("  max |generated - build_dhw| = $maxdiff;  max |generated - DHW file| = $filediff (file eltype $(eltype(file)))")
    maxdiff <= 1e-12 || error("$domain: generated trajectories differ from build_dhw")
    filediff <= 1e-4 * maximum(file) || error("$domain: generated trajectories differ from the DHW file")
end
