#=
02_design.jl -- the pick-and-freeze design of one model for one or more input bands.

    julia +release --project=Code/sa/shapley Code/sa_pf/02_design.jl <model> <n_base> <band> [band ...] [--sets=main,zero]

e.g. `02_design.jl coralpop 256 0.1 0.2 0.5`. Writes, per band and the same for both domains,
<ModelFolder>/model_run_pf_pilot/band<NN>/design_pf.jld2 (set `main`: with the DHW inputs, d = 176 for
CoralPop, 83 for C~scape) and design_pf_zero.jld2 (set `zero`: no heat stress, no DHW inputs, d = 168 /
75), each with `A`,
`B`, `perms`, `constants`, `spec`, `n_shared`, `names`, `band`, `n_base`, `set`, plus its sizes in
design_pf[_zero].toml. Runs = n_base x (d + 1). The folder is `CMR_PF_SUBDIR` (default model_run_pf_pilot).
A larger n_base gives the same first walks (pf_common.jl), so a design can be extended later.
=#

length(ARGS) >= 3 || error("usage: 02_design.jl <model> <n_base> <band> [band ...]")
const MODEL = ARGS[1]
const N_BASE = parse(Int, ARGS[2])
const SETS = let i = findfirst(startswith("--sets="), ARGS)
    i === nothing ? ["main", "zero"] : String.(split(ARGS[i][8:end], ","))
end
const BANDS = parse.(Float64, filter(a -> !startswith(a, "--"), ARGS[3:end]))
include(joinpath(@__DIR__, "pf_common.jl"))

for band in BANDS, set in SETS
    D = pf_design(band, N_BASE; set)
    d = ncol(D.A)
    path = pf_design_path(MODEL, band, set)
    JLD2.jldsave(path; D...)
    ## the sizes as plain text, for the HPC submit script (pf_submit.sh reads n_base without Julia)
    open(replace(path, r"\.jld2$" => ".toml"), "w") do io
        println(io, "n_base = $N_BASE\nd = $d\nruns_per_base = $(d + 1)\nband = $band\nset = \"$set\"")
    end
    println("$(band_label(band)) $set: $N_BASE base samples x $(d + 1) runs = $(N_BASE * (d + 1)) runs; ",
        "$d inputs ($(D.n_shared) shared), $(length(D.constants)) constants -> $path")
end
