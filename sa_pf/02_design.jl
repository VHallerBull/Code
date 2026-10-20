#=
02_design.jl -- the pick-and-freeze design of one model for one or more input bands.

    julia +release --project=Code/sa/shapley Code/sa_pf/02_design.jl <model> <n_base> <band> [band ...]

e.g. `02_design.jl coralpop 256 0.1 0.2 0.5`. Writes
<ModelFolder>/model_run_pf_pilot/band<NN>/design_pf.jld2 (`A`, `B`, `perms`, `constants`, `spec`,
`n_shared`, `names`, `band`, `n_base`), the same for both domains. Runs = n_base x (d + 1).
=#

length(ARGS) >= 3 || error("usage: 02_design.jl <model> <n_base> <band> [band ...]")
const MODEL = ARGS[1]
const N_BASE = parse(Int, ARGS[2])
const BANDS = parse.(Float64, ARGS[3:end])
include(joinpath(@__DIR__, "pf_common.jl"))

for band in BANDS
    D = pf_design(band, N_BASE)
    d = ncol(D.A)
    path = joinpath(pf_band_dir(MODEL, band), "design_pf.jld2")
    JLD2.jldsave(path; D...)
    println("$(band_label(band)): $N_BASE base samples x $(d + 1) runs = $(N_BASE * (d + 1)) runs; ",
        "$d inputs ($(D.n_shared) shared: DHW + weather + domain), $(length(D.constants)) constants -> $path")
end
