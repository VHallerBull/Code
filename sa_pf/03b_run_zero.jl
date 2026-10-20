#=
03b_run_zero.jl -- the zero-DHW reference: every walk's A and B parameter sets rerun with no heat stress.

    julia +release --project=Code/sa/shapley Code/sa_pf/03b_run_zero.jl <model> <proc> <n_procs> <band:domain> [band:domain ...]

Same split as 03_run_chunk.jl (process `proc` takes base samples proc, proc + n_procs, ...), 2 runs per
base sample. Writes <band>/<domain>/zero_chunks/proc_<pp>_of_<nn>.jld2 (`base`, `Y` = Dict output =>
[2, years, n_base_here] for the A and B sets, `failed`, `seconds`, `n_runs`, `years`); skipped if it
exists. These runs give the reference trajectories without heat stress (06_figures.jl); they are not part
of the Shapley analysis.
=#

length(ARGS) >= 4 || error("usage: 03b_run_zero.jl <model> <proc> <n_procs> <band:domain> ...")
const MODEL = ARGS[1]
const PROC, N_PROCS = parse(Int, ARGS[2]), parse(Int, ARGS[3])
const JOBS = [(parse(Float64, split(j, ":")[1]), String(split(j, ":")[2])) for j in ARGS[4:end]]
include(joinpath(@__DIR__, "pf_common.jl"))

const DOMS = Dict{String,Any}()

for (band, domain) in JOBS
    out_dir = mkpath(joinpath(pf_run_dir(MODEL, band, domain), "zero_chunks"))
    out = joinpath(out_dir, "proc_$(lpad(PROC, 2, '0'))_of_$(lpad(N_PROCS, 2, '0')).jld2")
    isfile(out) && (println("$(band_label(band))/$domain: $out exists, skipped"); continue)
    design = JLD2.load(joinpath(pf_band_dir(MODEL, band), "design_pf.jld2"))
    design = (; (Symbol(k) => v for (k, v) in design)...)
    dom = get!(() -> load_sa_domain(domain), DOMS, domain)
    years = collect(ADRIA.timesteps(dom))
    bases = collect(PROC:N_PROCS:design.n_base)
    Y = Dict(string(k) => fill(NaN, 2, length(years), length(bases)) for k in [PF_OUTPUTS; :relative_cover_engine])
    failed = Int[]
    t0 = time()
    for (j, n) in enumerate(bases)
        rows = zero_rows(design, n)
        try
            for i in 1:2
                r = pf_run_zero(dom, rows[i, :])
                for (k, v) in r.series
                    Y[string(k)][i, :, j] = v
                end
            end
        catch err
            push!(failed, n)
            @warn "$(band_label(band))/$domain zero-DHW base sample $n failed" exception = (err, catch_backtrace())
        end
    end
    seconds = time() - t0
    JLD2.jldsave(out * ".tmp"; base=bases, Y, failed, seconds, n_runs=2 * length(bases), years)
    mv(out * ".tmp", out; force=true)
    println("wrote $out ($(2 * length(bases)) zero-DHW runs, $(length(failed)) failed, $(round(seconds; digits=1)) s)")
end
