#=
05_run_chunk.jl -- run one chunk of a model's Sobol' design on one synthetic domain.

    julia --project=Code/sa/models/<model> Code/sa/05_run_chunk.jl <model> <domain> <set> <chunk> <chunk_size>

`<set>`: `main` (design.jld2), `b<k>` (extension block design_b<k>.jld2) or `control` (the first
`CMR_SA_N_CONTROL` = 1024 rows of design.jld2 with `dhw_scenario = 1`, the all-zero DHW member).
Chunk `c` holds rows `(c-1)*chunk_size+1 : c*chunk_size` of the set. On the HPC this is one Slurm array
task (`chunk = SLURM_ARRAY_TASK_ID`).

Writes
- <run dir>/chunks/<set>_<cccc>.jld2: `rows`, `scalars` (DataFrame, one row per run), `series`
  ([years, runs] domain relative cover), `reef_final`, `reef_area`, `reef_k` ([reefs, runs]),
  `failures`, `seconds`;
- <run dir>/runs_<set>/: every run's full outputs (`save_row_outputs` in models/<model>.jl; for
  CoralPop the ADRIA indicator file and the raw output, one file each per run).
A failing row is recorded and the chunk continues.
=#

length(ARGS) == 5 || error("usage: 05_run_chunk.jl <model> <domain> <set> <chunk> <chunk_size>")
const MODEL, DOMAIN, SET = ARGS[1], ARGS[2], ARGS[3]
const CHUNK, CHUNK_SIZE = parse(Int, ARGS[4]), parse(Int, ARGS[5])

include(joinpath(@__DIR__, "sa_common.jl"))
include(joinpath(@__DIR__, "models", "$MODEL.jl"))
import JLD2

const N_CONTROL = parse(Int, get(ENV, "CMR_SA_N_CONTROL", "1024"))

function set_design(dir, set)::DataFrame
    file = set in ("main", "control") ? "design.jld2" : "design_$set.jld2"
    path = joinpath(dir, file)
    isfile(path) || error("$path not found: run 04_design.jl $MODEL first.")
    d = JLD2.load(path, "design")
    set == "control" || return d
    c = d[1:min(N_CONTROL, nrow(d)), :]
    c.dhw_scenario .= 1.0
    return c
end

dir = model_run_dir(MODEL, DOMAIN)
design = set_design(dir, SET)
rows = ((CHUNK - 1) * CHUNK_SIZE + 1):min(CHUNK * CHUNK_SIZE, nrow(design))
isempty(rows) && (println("chunk $CHUNK is past the end of the $SET set ($(nrow(design)) rows); nothing to do"); exit(0))

println("$MODEL_LABEL / $DOMAIN / $SET: rows $(first(rows))-$(last(rows)) of $(nrow(design))")
t_load = @elapsed dom = load_sa_domain(DOMAIN)
println("  domain loaded in $(round(t_load; digits=1)) s")

done = Int[]
scalars = NamedTuple[]
series = Vector{Float64}[]
reef_final = Vector{Float64}[]
reef_area = Vector{Float64}[]
reef_k = Vector{Float64}[]
failures = DataFrame(row=Int[], error=String[])
full_dir = joinpath(dir, "runs_$SET")      # every run's full outputs
t0 = time()
for (k, i) in enumerate(rows)
    row = design[i, :]
    try
        r = run_row(dom, row)
        save_row_outputs(full_dir, row.row, r)
        push!(done, row.row); push!(scalars, r.scalars)
        push!(series, r.domain_series); push!(reef_final, r.reef_final)
        push!(reef_area, r.reef_area); push!(reef_k, r.reef_k)
    catch err
        push!(failures, (; row=Int(row.row), error=sprint(showerror, err, catch_backtrace())))
        @warn "row $(row.row) failed" exception = err
    end
    if k == 1 || k % 100 == 0 || k == length(rows)
        el = time() - t0
        println("  $k/$(length(rows)) rows, $(round(el; digits=1)) s ($(round(el / k; digits=3)) s/row), $(nrow(failures)) failed")
        flush(stdout)
    end
end
seconds = time() - t0

out_dir = mkpath(joinpath(dir, "chunks"))
out = joinpath(out_dir, "$(SET)_$(lpad(CHUNK, 4, '0')).jld2")
tmp = out * ".tmp"
JLD2.jldsave(tmp; rows=done, scalars=DataFrame(scalars),
    series=isempty(series) ? zeros(0, 0) : reduce(hcat, series),
    reef_final=isempty(reef_final) ? zeros(0, 0) : reduce(hcat, reef_final),
    reef_area=isempty(reef_area) ? zeros(0, 0) : reduce(hcat, reef_area),
    reef_k=isempty(reef_k) ? zeros(0, 0) : reduce(hcat, reef_k),
    failures, seconds, years=collect(ADRIA.timesteps(dom)),
    reef_ids=String.(dom.loc_ids))
mv(tmp, out; force=true)
println("wrote $out ($(length(done)) runs, $(nrow(failures)) failed, $(round(seconds; digits=1)) s)")
