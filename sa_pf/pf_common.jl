#=
pf_common.jl -- the pick-and-freeze (Goda 2021; SAShE `CallablePickAndFreezeSample`) design, option B:
the DHW block is the 7 generator parameters plus a weather index, and every run generates its own
trajectory. Separate from Code/sa: it only includes (reads) sa_common.jl and models/<model>.jl.

    MODEL = "coralpop"
    include(joinpath(@__DIR__, "pf_common.jl"))

Inputs (d = 176 for CoralPop): 7 DHW parameters (fixed ranges, SamplingSetup.docx Table 2) +
`dhw_weather` (1..1024) + 13 domain multipliers U(1 - band, 1 + band) + the model's native factors at
+-band. The first 21 (DHW + domain) are shared by every model: their A and B values come from one
model-independent sample, and their relative order in each walk is the same for every model, so every
model sees exactly the same DHW trajectories and domain settings in each base sample.

Run with the Code/sa/shapley environment (ADRIA, CoralPopInterface, SAShE).
=#

ENV["CMR_SA_RUN_SUBDIR"] = get(ENV, "CMR_SA_RUN_SUBDIR", "model_run_pf_pilot")  # never model_run
include(joinpath(@__DIR__, "..", "sa", "sa_common.jl"))
include(joinpath(@__DIR__, "..", "sa", "models", "$MODEL.jl"))
include(joinpath(@__DIR__, "dhw_generator.jl"))
using SAShE
using Random: Xoshiro, randperm
import JLD2

"Folder of the shared pick-and-freeze inputs (DHW draws), one per domain."
const PF_INPUTS = joinpath(CMR_DIR, "pf_inputs")
"Output folder of the pilot runs (never Code/sa's model_run / model_run_local)."
const PF_SUBDIR = "model_run_pf_pilot"

band_label(band::Real) = "band$(lpad(round(Int, 100 * band), 2, '0'))"
pf_band_dir(model, band) = mkpath(joinpath(SA_RUNS_ROOT, MODEL_FOLDERS[model], PF_SUBDIR, band_label(band)))
pf_run_dir(model, band, domain) = mkpath(joinpath(pf_band_dir(model, band), domain))

"""
Design sets: `"main"` (the DHW block as inputs, every run with its own trajectory) and `"zero"` (the
sensitivity analysis without heat stress: no DHW inputs, every run on the domain's all-zero member).
"""
const PF_SETS = ("main", "zero")
set_suffix(set) = set == "main" ? "" : "_$set"
pf_design_path(model, band, set="main") = joinpath(pf_band_dir(model, band), "design_pf$(set_suffix(set)).jld2")
"Run folder of a set: the domain folder for `main`, its `pf_zero_dhw/` subfolder for `zero`."
pf_set_dir(model, band, domain, set="main") =
    set == "main" ? pf_run_dir(model, band, domain) : mkpath(joinpath(pf_run_dir(model, band, domain), "pf_zero_dhw"))

"Slot of the domain's DHW array that every run overwrites with its own trajectory (in memory only)."
const PF_DHW_SLOT = 2

## ---------------------------------------------------------------------------------------------------
## Specs
## ---------------------------------------------------------------------------------------------------

"The weather index: which column of fixed random draws (= weather of ensemble member w + 1)."
weather_spec()::DataFrame = _spec([(:dhw_weather, 512.0, DiscreteUniform, (1.0, Float64(N_DHW_MEMBERS)), false)])

"Shared block: 7 DHW parameters (fixed ranges), the weather index, 13 multipliers at +-band."
shared_spec(band::Real)::DataFrame = vcat(dhw_spec(), weather_spec(), multiplier_spec((1 - band, 1 + band)))

"The model's native factors at +-band (CoralPop: the same call as models/coralpop.jl, other band)."
function model_spec_band(band::Real)::DataFrame
    MODEL == "coralpop" || error("model_spec_band: add $MODEL")
    return coralpop_native_spec(; band=Float64(band), larvae_via_connectivity=LARVAE_VIA_CONNECTIVITY,
        rcp=parse(Int, SA_RCP), year_start=SA_YEARS[1], year_end=SA_YEARS[2])
end

## ---------------------------------------------------------------------------------------------------
## Design
## ---------------------------------------------------------------------------------------------------

"Seeds (Monte Carlo streams): the shared block's are model independent (same A/B DHW and domain values for every model)."
const PF_SEEDS = (; shared_A=101, shared_B=102, model_A=201, model_B=202, order_shared=301, order_model=401)

"""
    sample_varying(spec, n, seed) -> DataFrame

Plain Monte Carlo draws of the varying factors of `spec` (each through its distribution's quantile
function), from `Xoshiro(seed)`. Not ADRIA's Sobol' sampler: two of its samples with different seeds
are strongly correlated row by row (A[n] vs B[n] up to |r| = 0.98, checked on 2048 rows), which breaks
pick-and-freeze, whose A[n] and B[n] must be independent. Goda (2021) and SAShE's documented workflow use
Monte Carlo too.
"""
function sample_varying(spec::DataFrame, n::Integer, seed::Integer)::DataFrame
    v = spec[.!spec.is_constant, :]
    U = rand(Xoshiro(seed), n, nrow(v))
    cols = map(enumerate(eachrow(v))) do (j, r)
        dist = r.dist == DiscreteUniform ? DiscreteUniform(Int.(r.dist_params)...) : r.dist(r.dist_params...)
        r.fieldname => Float64.(quantile.(dist, U[:, j]))
    end
    return DataFrame(cols)
end

"""
    pf_permutations(n_base, n_shared, d) -> Matrix{Int} [n_base, d]

One uniform random order of the `d` inputs per base sample, in which the relative order of the first
`n_shared` (the shared block) depends only on `n_base` and `n_shared`: drawn from its own stream, then
written into the positions a uniform permutation of all `d` gives them. The result is still a uniform
permutation, and every model walks the shared inputs in the same order.
"""
function pf_permutations(n_base::Integer, n_shared::Integer, d::Integer;
    seed_shared=PF_SEEDS.order_shared, seed_model=PF_SEEDS.order_model)::Matrix{Int}
    rs, rm = Xoshiro(seed_shared), Xoshiro(seed_model)
    P = zeros(Int, n_base, d)
    for n in 1:n_base
        shared_order = randperm(rs, n_shared)
        p = randperm(rm, d)
        pos = findall(<=(n_shared), p)
        p[pos] = shared_order
        P[n, :] = p
    end
    return P
end

"""
    pf_design(band, n_base; set="main") -> NamedTuple

`A`, `B` (`[n_base, d]`, varying inputs only), `perms` (`[n_base, d]`), `constants` (fixed factors and
their values), `spec`, `n_shared`, `names`, `set`.

`set = "zero"`: the design without heat stress. The DHW inputs (7 parameters + weather) are dropped and
every run uses the all-zero DHW member; the other inputs keep exactly the main design's A and B values
(so the walks' end points are the main design's parameter sets without DHW), and new walk orders are
drawn for the remaining d = 13 + model inputs (the 13 shared ones in the same order for every model).
"""
function pf_design(band::Real, n_base::Integer; set::AbstractString="main")
    set in PF_SETS || error("unknown set '$set' ($(join(PF_SETS, ", ")))")
    ss, ms = shared_spec(band), model_spec_band(band)
    A = hcat(sample_varying(ss, n_base, PF_SEEDS.shared_A), sample_varying(ms, n_base, PF_SEEDS.model_A))
    B = hcat(sample_varying(ss, n_base, PF_SEEDS.shared_B), sample_varying(ms, n_base, PF_SEEDS.model_B))
    spec = vcat(ss, ms)
    n_shared = count(.!ss.is_constant)
    if set == "zero"
        dhw = filter(startswith("dhw_"), names(A))
        select!(A, Not(dhw))
        select!(B, Not(dhw))
        spec = spec[.!startswith.(String.(spec.fieldname), "dhw_"), :]
        n_shared -= length(dhw)
    end
    d = ncol(A)
    constants = Dict(String(r.fieldname) => r.val for r in eachrow(spec) if r.is_constant)
    perms = pf_permutations(n_base, n_shared, d)
    return (; A, B, perms, constants, spec=spec_for_disk(spec), n_shared, names=names(A), band, n_base, set)
end

"""
    walk_rows(design, n) -> DataFrame (d + 1 rows)

The runs of base sample `n`: row 0 = A[n], then each input in the order `perms[n, :]` swapped to its
B[n] value (SAShE's `_build_pick_freeze_block` layout), plus the constants and `dhw_scenario`:
`PF_DHW_SLOT` (the run's own trajectory goes there) when the design has the DHW inputs, else 1 (the
all-zero member: the `zero` set).
"""
function walk_rows(design, n::Integer)::DataFrame
    A, B, p = design.A, design.B, design.perms
    d = ncol(A)
    rows = DataFrame([Symbol(c) => fill(Float64(A[n, c]), d + 1) for c in names(A)])
    for k in 1:d
        c = p[n, k]
        rows[(k + 1):end, c] .= Float64(B[n, c])
    end
    for (k, v) in design.constants
        rows[!, k] .= v
    end
    rows.dhw_scenario .= "dhw_weather" in names(A) ? Float64(PF_DHW_SLOT) : 1.0
    return rows
end

## ---------------------------------------------------------------------------------------------------
## One run
## ---------------------------------------------------------------------------------------------------

"Outputs kept per run (domain level, habitable-area weighted with the run's own k x area)."
const PF_OUTPUTS = [:relative_cover, :relative_shelter_volume, :coral_evenness, :juvenile_indicator]

"""
    pf_run(dom, row, draws) -> (; series::Dict{Symbol,Vector}, dhw::Vector)

Write the row's own DHW trajectory into `dom`'s slot `PF_DHW_SLOT` (in memory; `run_row` copies the
domain before applying the multipliers), run the model, and return the domain series of `PF_OUTPUTS`
plus the domain-mean DHW per year. A row without `dhw_weather` (the `zero` set, `dhw_scenario` = 1) runs
on the all-zero member and returns DHW = 0.
"""
function pf_run(dom, row, draws)
    if !(:dhw_weather in propertynames(row))                     # zero set: the all-zero member
        r = run_row(dom, row)
        return (; series=pf_series(r), dhw=zeros(length(r.domain_series)))
    end
    traj = generate_dhw(row, Int(row[:dhw_weather]), draws)
    dom.dhw_scens.data[:, :, PF_DHW_SLOT] .= traj
    return (; series=pf_series(run_row(dom, row)), dhw=vec(mean(traj; dims=2)))
end

"Domain series of `PF_OUTPUTS` from a `run_row` result (habitable-area weights of that run)."
function pf_series(r)
    ind = coralpop_indicators(r.out)
    w = r.reef_k .* r.reef_area
    w = w ./ sum(w)
    series = Dict{Symbol,Vector{Float64}}()
    for k in PF_OUTPUTS
        A = getfield(ind, k)
        ndims(A) == 2 || error("$k is not [years, reefs]")
        series[k] = vec(sum(Float64.(A) .* reshape(w, 1, :); dims=2))
    end
    series[:relative_cover_engine] = Float64.(r.domain_series)
    return series
end

## ---------------------------------------------------------------------------------------------------
## Zero-DHW reference
## ---------------------------------------------------------------------------------------------------

"""
    zero_rows(design, n) -> DataFrame (2 rows)

The A[n] and B[n] parameter sets (start and end of walk n) with the DHW switched off: `dhw_scenario` = 1,
the domain's all-zero member (as the zero-DHW control of Code/sa). A reference without heat stress for
the output trajectories; not part of the Shapley analysis.
"""
function zero_rows(design, n::Integer)::DataFrame
    W = walk_rows(design, n)
    Z = W[[1, nrow(W)], :]
    Z.dhw_scenario .= 1.0
    return Z
end

"One zero-DHW run: `run_row` on the member-1 (all-zero) DHW, no trajectory generated."
pf_run_zero(dom, row) = (; series=pf_series(run_row(dom, row)))
