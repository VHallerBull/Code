#=
03_run_chunk.jl -- run one process's share of the pick-and-freeze walks, for one or more jobs.

    julia +release --project=Code/sa/shapley Code/sa_pf/03_run_chunk.jl <model> <proc> <n_procs> <band:domain> [band:domain ...]

Process `proc` (1..n_procs) takes base samples proc, proc + n_procs, ... of each job (a whole walk of
d + 1 runs each), e.g. `03_run_chunk.jl coralpop 3 12 0.1:SingleReef_SobolDHW 0.5:NineReefs_SobolDHW:zero`
(`:zero` = the design without heat stress, written under <band>/<domain>/pf_zero_dhw/).
Each job writes <band>/<domain>/chunks/proc_<pp>_of_<nn>.jld2 (skipped if it exists, so a rerun resumes):
`base` (base-sample indices), `S` (the runs' state summaries, sa_outputs.jl: Dict key => [d + 1, dims...,
n_base_here]), `dhw` (domain-mean DHW, [d + 1, years, n_base_here]), `failed` (base samples with a failed
run), `seconds`, `n_runs`, `years`, `thermal_groups`. No other run output is kept.
CoralPop runs through Python, so parallelism is one process per share (no threads).
=#

length(ARGS) >= 4 || error("usage: 03_run_chunk.jl <model> <proc> <n_procs> <band:domain> ...")
const MODEL = ARGS[1]
const PROC, N_PROCS = parse(Int, ARGS[2]), parse(Int, ARGS[3])
"Jobs `band:domain` (set main) or `band:domain:zero` (the zero-DHW design)."
const JOBS = [(parse(Float64, p[1]), String(p[2]), length(p) > 2 ? String(p[3]) : "main")
              for p in (split(j, ":") for j in ARGS[4:end])]
include(joinpath(@__DIR__, "pf_common.jl"))

const DOMS = Dict{String,Any}()
const DRAWS = Dict{String,Any}()

for (band, domain, set) in JOBS
    out_dir = mkpath(joinpath(pf_set_dir(MODEL, band, domain, set), "chunks"))
    out = joinpath(out_dir, "proc_$(lpad(PROC, 2, '0'))_of_$(lpad(N_PROCS, 2, '0')).jld2")
    isfile(out) && (println("$(band_label(band))$(set_suffix(set))/$domain: $out exists, skipped"); continue)
    design = JLD2.load(pf_design_path(MODEL, band, set))
    design = (; (Symbol(k) => v for (k, v) in design)...)
    dom = get!(() -> load_sa_domain(domain), DOMS, domain)
    draws = get!(() -> load_dhw_draws(joinpath(PF_INPUTS, domain, "dhw_draws.jld2")), DRAWS, domain)
    years = collect(ADRIA.timesteps(dom))
    T, d = length(years), ncol(design.A)
    bases = collect(PROC:N_PROCS:design.n_base)
    S = nothing                      # state summaries: key => [d + 1, dims..., walks], allocated at the first run
    dhw = fill(NaN, d + 1, T, length(bases))
    failed = Int[]
    t0 = time()
    for (j, n) in enumerate(bases)
        rows = walk_rows(design, n)
        try
            for i in 1:nrow(rows)
                r = pf_run(dom, rows[i, :], draws)
                S === nothing && (S = state_store(r.state, d + 1, length(bases)))
                store_state!(S, r.state, i, j)
                dhw[i, :, j] = r.dhw
            end
        catch err
            push!(failed, n)
            @warn "$(band_label(band))$(set_suffix(set))/$domain base sample $n failed" exception = (err, catch_backtrace())
        end
        el = time() - t0
        println("$(band_label(band))$(set_suffix(set))/$domain proc $PROC: $j/$(length(bases)) walks, ",
            "$(round(el / (j * (d + 1)); digits=3)) s/run, $(round(el / 60; digits=1)) min")
        flush(stdout)
    end
    seconds = time() - t0
    JLD2.jldsave(out * ".tmp"; base=bases, S, dhw, failed, seconds, n_runs=length(bases) * (d + 1), years,
        thermal_groups=PF_THERMAL_GROUPS)
    mv(out * ".tmp", out; force=true)
    println("wrote $out ($(length(bases)) walks, $(length(failed)) failed, $(round(seconds / 60; digits=1)) min)")
end
