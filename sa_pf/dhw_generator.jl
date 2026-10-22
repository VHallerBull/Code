#=
dhw_generator.jl -- one DHW trajectory from the 7 generator parameters and a weather index, exactly as
SyntheticReefDomain's `build_dhw` makes ensemble member `w + 1` (no package dependencies).

    include(joinpath(@__DIR__, "dhw_generator.jl"))
    D = load_dhw_draws(path)                          # written by 01_dhw_draws.jl
    traj = generate_dhw(p, w, D)                      # [years, reefs]

`p` has the fields (or keys) of `DHW_PARAM_NAMES` in sa_common.jl: dhw_baseline, dhw_trend_per_decade,
dhw_interannual_sd, dhw_ar1, dhw_event_prob, dhw_event_mag_width, dhw_local_sd. `w` (1..1024) selects
the column of fixed random draws (AR(1) innovations, heatwave occurrence and size, reef-level noise):
the weather realisation of stochastic ensemble member `w + 1`. With that member's parameters the result
equals the domain's DHW file (checked by 01_dhw_draws.jl).
=#

import JLD2

"Lower bound of the heatwave size: a heatwave adds U(5, 5 + width) DHW (SamplingSetup.docx 2.3)."
const PF_DHW_EVENT_MAG_MIN = 5.0

"""
    load_dhw_draws(path) -> NamedTuple

`z`, `u_event`, `u_mag` ([years, members]), `local_noise` ([years, reefs, members]), `exposure`, `offset`
([reefs]) and `prescribed` ([years], prescribed heatwave years) of one domain.
"""
load_dhw_draws(path::AbstractString) = (; (Symbol(k) => v for (k, v) in JLD2.load(path))...)

_get(p, k::Symbol) = p isa AbstractDict ? p[k] : p[k]

"""
    generate_dhw(p, w, D) -> Matrix{Float64} [years, reefs]

The formula of SyntheticReefDomain `build_dhw` for one scenario, with draws column `w`:

    R[t]      = baseline + trend * (t - 1) / 10 + anomaly[t] + event[t]
    DHW[t, r] = max(0, exposure[r] * R[t] + offset[r] + local_sd * noise[t, r])
"""
function generate_dhw(p, w::Integer, D)::Matrix{Float64}
    T, n = size(D.local_noise, 1), size(D.local_noise, 2)
    1 <= w <= size(D.z, 2) || throw(ArgumentError("weather index $w outside 1:$(size(D.z, 2))"))
    sd, ar1 = Float64(_get(p, :dhw_interannual_sd)), Float64(_get(p, :dhw_ar1))
    lo = PF_DHW_EVENT_MAG_MIN
    hi = PF_DHW_EVENT_MAG_MIN + Float64(_get(p, :dhw_event_mag_width))
    baseline, trend = Float64(_get(p, :dhw_baseline)), Float64(_get(p, :dhw_trend_per_decade))
    prob, local_sd = Float64(_get(p, :dhw_event_prob)), Float64(_get(p, :dhw_local_sd))
    innovation_sd = sd * sqrt(1 - ar1^2)
    anomaly = sd * D.z[1, w]
    R = Vector{Float64}(undef, T)
    for t in 1:T
        t > 1 && (anomaly = ar1 * anomaly + innovation_sd * D.z[t, w])
        is_event = D.prescribed[t] || D.u_event[t, w] < prob
        event = is_event ? lo + (hi - lo) * D.u_mag[t, w] : 0.0
        R[t] = baseline + trend * (t - 1) / 10 + anomaly + event
    end
    return [max(0.0, D.exposure[r] * R[t] + D.offset[r] + local_sd * D.local_noise[t, r, w]) for t in 1:T, r in 1:n]
end
