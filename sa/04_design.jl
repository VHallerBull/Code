#=
04_design.jl -- the Sobol' design of one coral model (the same design for both synthetic domains).

    julia +release --project=Code/sa/models/<model> Code/sa/04_design.jl <model> [--block k]

Factors (d varying), one row = one model run:

    dhw_scenario       DHW member 2..1025 of the shared ensemble (member 1, all-zero DHW, is the
                       separate control); its 7 generator parameters are attached in 06_collect.jl
    m_*                the 13 domain multipliers (sa_common.jl), bands from multiplier_bands.toml
                       (`CMR_SA_MULTIPLIER_BANDS`); fixed ones (lo == hi) are constants
    <model factors>    `model_spec_block()` from models/<model>.jl

`ADRIA.sample(spec, N)` (Owen-scrambled Sobol'), N = 8192 by default (`CMR_SA_N`, a power of 2), seed
1234. Writes <CMR>/<ModelFolder>/model_run/<domain>/design.jld2 (`design`, `spec`, `varying`,
`settings`, incl. the multiplier bands) and design.csv. Everything downstream (05 runs, 06/07) takes the
bands from the design's `spec`, so new bands need a new design, not a new image. The zero-DHW control runs the first `CMR_SA_N_CONTROL` (1024) rows with
`dhw_scenario = 1` (05_run_chunk.jl, set `control`).

`--block k` (k >= 2) writes an independent extension block design_b<k>.jld2 (seed 1234 + k - 1), for
when the Shapley estimates need more runs: given-data estimators simply take the extra rows.
=#

const MODEL = get(ARGS, 1, "")
isempty(MODEL) && error("usage: 04_design.jl <model> [--block k]")
const BLOCK = let i = findfirst(==("--block"), ARGS)
    i === nothing ? 1 : parse(Int, ARGS[i + 1])
end

include(joinpath(@__DIR__, "sa_common.jl"))
include(joinpath(@__DIR__, "models", "$MODEL.jl"))
using CSV, Dates
import JLD2

const N = parse(Int, get(ENV, "CMR_SA_N", "8192"))
ispow2(N) || error("CMR_SA_N must be a power of 2 (ADRIA's Sobol' sampler), got $N")
const SEED = 1234 + BLOCK - 1

bands = multiplier_bands()
spec = vcat(dhw_selector_spec(2:(N_DHW_MEMBERS + 1)), multiplier_spec(bands), model_spec_block())
dups = [n for n in unique(spec.fieldname) if count(==(n), spec.fieldname) > 1]
isempty(dups) || error("Duplicate factor names across blocks: $(join(dups, ", "))")
varying = spec.fieldname[.!spec.is_constant]

Random.seed!(SEED)
design = ADRIA.sample(spec, N)
insertcols!(design, 1, :row => 1:N)

settings = Dict{String,Any}(
    "model" => MODEL, "model_label" => MODEL_LABEL, "n" => N, "block" => BLOCK, "seed" => SEED,
    "n_varying" => length(varying),
    "n_per_block" => Dict(b => count(f -> factor_block(f) == b, varying) for b in ("DHW", "domain", "model")),
    "dhw_members" => "2:$(N_DHW_MEMBERS + 1) (member 1 = zero-DHW control)",
    "years" => collect(SA_YEARS), "rcp" => SA_RCP,
    "multiplier_bands" => Dict(String(n) => collect(bands[n]) for n in MULTIPLIER_NAMES),
    "multiplier_bands_file" => MULTIPLIER_BANDS_FILE,
    "created" => string(now()), "julia" => string(VERSION),
)
println("domain multipliers (", MULTIPLIER_BANDS_FILE, "):")
for n in MULTIPLIER_NAMES
    lo, hi = bands[n]
    println("  ", rpad(n, 18), lo == hi ? "fixed at $lo" : "U($lo, $hi)")
end
tag = BLOCK == 1 ? "" : "_b$BLOCK"
for domain in SA_DOMAINS
    dir = model_run_dir(MODEL, domain)
    JLD2.jldsave(joinpath(dir, "design$tag.jld2"); design, spec=spec_for_disk(spec),
        varying=String.(varying), settings=merge(settings, Dict("domain" => domain)))
    CSV.write(joinpath(dir, "design$tag.csv"), design)
    println("$domain: $(joinpath(dir, "design$tag.jld2")) -- $N rows x $(length(varying)) varying factors ",
        "(DHW $(settings["n_per_block"]["DHW"]), domain $(settings["n_per_block"]["domain"]), ",
        "$MODEL_LABEL $(settings["n_per_block"]["model"]))")
end
