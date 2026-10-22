#=
06_collect.jl -- merge a model's chunk files into one result table per set (main, b<k>, control).

    julia --project=Code/sa/models/<model> Code/sa/06_collect.jl <model> [domain ...]

For every set found in <run dir>/chunks/ it writes, next to design.jld2:

    results.jld2 / results.csv                (main)
    results_b<k>.jld2 / .csv                  (extension blocks)
    control_results.jld2 / .csv               (zero-DHW control)

The JLD2 holds `X` (the inputs for the Shapley analysis: every varying design factor except the
`dhw_scenario` member index, plus that member's 7 DHW generator parameters from
domains/<domain>/sa_inputs/dhw_design.csv), `Y` (the scalar outputs), `rows`, `dhw_scenario`, `series`
([years, runs]), `reef_final`, `reef_area`, `reef_k` ([reefs, runs]; area and k differ between runs
through the multipliers), `years`, `reef_ids`, `failures`, `settings`. The CSV is `row`,
`dhw_scenario`, `X`, `Y` side by side. SAShE: `DataModel(res["X"], res["Y"].rc_final)`.

If the runs kept their full outputs (runs_<set>/adria_exports/, see models/<model>.jl) it also writes

    result_set.jld2 / result_set_b<k>.jld2 / control_result_set.jld2

the model's ADRIA result set over the same runs in the same order: every indicator per year and reef,
`rs.outcomes[:reef_condition_index]` etc. `[timesteps, locations, scenarios]`, with `rs.inputs` = `row`
(`ID`), `dhw_scenario` and `X`. Load with `load_sa_result_set` (sa_common.jl).

Missing or failed rows are listed and make the script exit with an error (after writing what it has).
=#

const MODEL = get(ARGS, 1, "")
isempty(MODEL) && error("usage: 06_collect.jl <model> [domain ...]")
const DOMAINS = length(ARGS) > 1 ? ARGS[2:end] : nothing

include(joinpath(@__DIR__, "sa_common.jl"))
include(joinpath(@__DIR__, "models", "$MODEL.jl"))
using CSV
import JLD2

const N_CONTROL = parse(Int, get(ENV, "CMR_SA_N_CONTROL", "1024"))

result_name(set) = set == "main" ? "results" : set == "control" ? "control_results" : "results_$set"
result_set_name(set) = set == "main" ? "result_set" : set == "control" ? "control_result_set" : "result_set_$set"
design_file(set) = set in ("main", "control") ? "design.jld2" : "design_$set.jld2"

"""
Spec of the Shapley inputs `X` for the result set: the 7 DHW parameters in place of the member index,
then the design's own multiplier and model rows (its saved `spec`), so the recorded bounds are the ones
the design was drawn from -- not whatever this copy of the code (e.g. the image's) defaults to.
"""
analysis_spec(design_spec::DataFrame) =
    vcat(dhw_spec(), filter(r -> r.fieldname != :dhw_scenario, spec_from_disk(design_spec)))

"Stack the runs' indicator files into the model's result set, in the order of `rows`."
function collect_result_set(dir, domain, set, rows, dhw_scenario, X, spec)
    fpath = joinpath(dir, "runs_$set")
    ind_dir = joinpath(fpath, "adria_exports")
    isdir(ind_dir) || (println("  (no full outputs in runs_$set/ -- result set skipped)"); return)
    n_files = count(f -> endswith(f, ".jld2"), readdir(ind_dir))
    n_files == length(rows) || @warn "$domain/$set: $n_files indicator files for $(length(rows)) runs; using the runs' rows"
    scenarios = hcat(DataFrame(ID=rows, dhw_scenario=dhw_scenario, rcp=fill(parse(Int, SA_RCP), length(rows))), X)
    t = @elapsed rs = build_result_set(fpath, scenarios, load_sa_domain(domain), spec)
    Int.(rs.inputs.ID) == rows || error("$domain/$set: result-set scenario order differs from results.jld2")
    path = save_result_set(rs, joinpath(dir, result_set_name(set) * ".jld2"))
    println("  -> $(basename(path)): $(length(rows)) scenarios, outcomes ",
        join(sort(string.(collect(keys(rs.outcomes)))), ", "), " ($(round(t; digits=1)) s)")
    return nothing
end

function collect_set(dir, domain, set)
    files = sort(filter(f -> startswith(f, set * "_") && endswith(f, ".jld2"), readdir(joinpath(dir, "chunks"))))
    D = JLD2.load(joinpath(dir, design_file(set)))
    design, varying, settings = D["design"], Symbol.(D["varying"]), D["settings"]
    set == "control" && (design = design[1:min(N_CONTROL, nrow(design)), :])

    parts = [JLD2.load(joinpath(dir, "chunks", f)) for f in files]
    rows = reduce(vcat, [p["rows"] for p in parts])
    order = sortperm(rows)
    rows = rows[order]
    allunique(rows) || error("$domain/$set: duplicated rows across chunk files")
    Y = reduce(vcat, [p["scalars"] for p in parts if nrow(p["scalars"]) > 0])[order, :]
    series = reduce(hcat, [p["series"] for p in parts if !isempty(p["series"])])[:, order]
    reef_final = reduce(hcat, [p["reef_final"] for p in parts if !isempty(p["reef_final"])])[:, order]
    stack_(k) = all(p -> haskey(p, k), parts) ?
                reduce(hcat, [p[k] for p in parts if !isempty(p[k])])[:, order] : nothing
    reef_area, reef_k = stack_("reef_area"), stack_("reef_k")
    failures = reduce(vcat, [p["failures"] for p in parts])
    missing_rows = setdiff(design.row, rows)

    ## inputs: the varying factors, with the DHW member index replaced by its 7 generator parameters
    d = design[indexin(rows, design.row), :]
    X = select(d, [v for v in varying if v != :dhw_scenario])
    dhw = CSV.read(joinpath(sa_inputs_dir(domain), "dhw_design.csv"), DataFrame)
    idx = indexin(Int.(d.dhw_scenario), dhw.member)
    any(isnothing, idx) && error("$domain/$set: DHW members missing from dhw_design.csv")
    for p in DHW_PARAM_NAMES
        X[!, p] = Float64.(coalesce.(dhw[idx, p], 0.0))     # the control member has no parameters: 0
    end
    select!(X, DHW_PARAM_NAMES, Not(DHW_PARAM_NAMES))

    res = Dict{String,Any}(
        "X" => X, "Y" => Y, "rows" => rows, "dhw_scenario" => Int.(d.dhw_scenario),
        "series" => series, "reef_final" => reef_final, "reef_area" => reef_area, "reef_k" => reef_k,
        "years" => parts[1]["years"],
        "reef_ids" => parts[1]["reef_ids"], "failures" => failures,
        "settings" => merge(settings, Dict("set" => set, "n_runs" => length(rows),
            "seconds" => sum(p["seconds"] for p in parts)))
    )
    name = result_name(set)
    JLD2.save(joinpath(dir, "$name.jld2"), res)
    CSV.write(joinpath(dir, "$name.csv"), hcat(DataFrame(row=rows, dhw_scenario=Int.(d.dhw_scenario)), X, Y))
    println("$domain/$set: $(length(rows)) of $(nrow(design)) runs -> $name.jld2 ",
        "($(size(X, 2)) inputs, $(size(Y, 2)) outputs; $(length(files)) chunks, ",
        "$(round(res["settings"]["seconds"] / 3600; digits=2)) CPU-h)")
    if !isempty(missing_rows)
        println("  MISSING $(length(missing_rows)) rows (failed: $(nrow(failures))), first: $(first(missing_rows, 10))")
    end
    collect_result_set(dir, domain, set, rows, Int.(d.dhw_scenario), X, analysis_spec(D["spec"]))
    return isempty(missing_rows)
end

ok = true
for domain in something(DOMAINS, SA_DOMAINS)
    dir = model_run_dir(MODEL, domain)
    isdir(joinpath(dir, "chunks")) || (println("$domain: no chunks/ yet"); continue)
    sets = unique([first(split(f, '_')) for f in readdir(joinpath(dir, "chunks")) if endswith(f, ".jld2")])
    for set in sets
        global ok &= collect_set(dir, domain, set)
    end
end
if !ok
    get(ENV, "CMR_SA_ALLOW_PARTIAL", "false") == "true" ?
        @warn("some rows are missing (CMR_SA_ALLOW_PARTIAL=true: a test run of part of the design)") :
        error("some rows are missing; rerun their chunks (05_run_chunk.jl) and collect again")
end
