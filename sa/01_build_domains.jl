#=
01_build_domains.jl -- build the two shared synthetic domains, each with the shared DHW ensemble.

    julia +release --project=Code/sa Code/sa/01_build_domains.jl [--force]

Writes CoralModelReview/domains/SingleReef_SobolDHW and NineReefs_SobolDHW (ADRIA data packages, made
with SyntheticReefDomain.jl). Each domain's generator_config.toml records exactly how it was built.

- SingleReef_SobolDHW: one reef at median values (Moore IQR midpoints; area 100,000 m2; ub 0.22 m/s,
  the ECORRAP median).
- NineReefs_SobolDHW: 3 x 3 grid, 5 km apart, per-reef seeded uniform draws within the Moore IQRs (area
  and ub re-centred on 100,000 m2 and 0.22 m/s), so the reefs differ but are identical for every model.
- Both: distance-decay connectivity (self 0.2, strength 0.05, decay 10 km), cyclones off, 2025-2050,
  DHW reef exposure 1 on every reef, and 1024 stochastic DHW members whose 7 generator parameters are
  one Sobol' sample (`dhw_sample()`, the same for both domains) plus member 1 = all-zero DHW.

`--force` replaces existing domain folders.
=#

using SyntheticReefDomain
using CairoMakie  # domain README figures
include(joinpath(@__DIR__, "sa_common.jl"))

const FORCE = "--force" in ARGS

"The shared DHW settings: one Sobol' row per stochastic member."
function sa_dhw_config(X::DataFrame)::DHWConfig
    col(name) = Float64.(X[!, name])
    return DHWConfig(;
        baseline=col(:dhw_baseline),
        trend_per_decade=col(:dhw_trend_per_decade),
        interannual_sd=col(:dhw_interannual_sd),
        ar1=col(:dhw_ar1),
        event_prob=col(:dhw_event_prob),
        local_sd=col(:dhw_local_sd),
        event_magnitude=[DHW_EVENT_MAG_MIN, DHW_EVENT_MAG_MIN + 3.0],  # shared default; bounds below
        event_magnitude_min=fill(DHW_EVENT_MAG_MIN, nrow(X)),
        event_magnitude_max=DHW_EVENT_MAG_MIN .+ col(:dhw_event_mag_width),
        reef_exposure=1.0,
        reef_offset=0.0,
        climate=Dict(parse(Int, SA_RCP) => Dict{String,Any}())
    )
end

const COMMON = (;
    timeframe=collect(SA_YEARS),
    seed=1234,
    n_dhw_scenarios=N_DHW_MEMBERS,
    connectivity=ConnectivityConfig(; topology="distance_decay", self_recruitment=0.2, strength=0.05,
        decay_km=10.0),
    cyclones=CycloneConfig(; enabled=false),
)
const COMPOSITION = [0.26, 0.25, 0.07, 0.19, 0.23]  # ADRIA functional-group split of initial cover

function domain_configs(dhw::DHWConfig)
    single = DomainConfig(;
        name="SingleReef_SobolDHW", COMMON...,
        reefs=ReefConfig(; n=1, area_m2=400_000.0, k_percent=54.0, depth_med=8.5, ub_med=0.22,
            cb_calib_group=6),
        cover=CoverConfig(; composition=COMPOSITION, jitter=0.0, total_percent=25.0),
        dhw
    )
    nine = DomainConfig(;
        name="NineReefs_SobolDHW", COMMON...,
        reefs=ReefConfig(; n=9, layout="grid", spacing_km=5.0,
            area_m2=(min=100_000.0, max=400_000.0), k_percent=(min=40.0, max=80.0),
            depth_med=(min=2.0, max=15.0), ub_med=(min=0.06, max=0.38), cb_calib_group=6),
        cover=CoverConfig(; composition=COMPOSITION, jitter=0.0, total_percent=(min=10.0, max=40.0)),
        dhw
    )
    return [single, nine]
end

X = dhw_sample()
println("DHW ensemble: $(nrow(X)) Sobol' members of $(length(DHW_PARAM_NAMES)) parameters (+ member 1 = zero DHW)")
for cfg in domain_configs(sa_dhw_config(X))
    validate(cfg)
    out = domain_path(cfg.name)
    if isdir(out) && !FORCE
        println("  $(cfg.name): exists, skipped (use --force to rebuild)")
        continue
    end
    path = generate_domain(cfg; outdir=out, force=FORCE)
    println("  $(cfg.name): written to $path")
end
