#=
08_indicator_outputs.jl -- a ready-made table of Shapley outputs: many indicators x years, one row per run.

    julia --project=Code/sa/models/<model> Code/sa/08_indicator_outputs.jl <model> [domain ...]

Reads results.jld2 + result_set.jld2 (06_collect.jl) and writes <run dir>/indicator_outputs.csv:
`row`, then one column per indicator and year, `<indicator>_<year>` for 2030, 2040, 2050 and
`<indicator>_mean` (mean over 2025-2050), each the habitable-area weighted domain mean of that run.
Rows are in the order of `X` in results.jld2, so for SAShE:

    R = load_sa_results(model, domain); Y = CSV.read(".../indicator_outputs.csv", DataFrame)
    DataModel(R["X"], Y.reef_condition_index_2040)

Other choices (one reef, another year, the minimum): `indicator_output` in sa_common.jl.
=#

const MODEL = get(ARGS, 1, "")
isempty(MODEL) && error("usage: 08_indicator_outputs.jl <model> [domain ...]")
const DOMAINS = length(ARGS) > 1 ? ARGS[2:end] : nothing

include(joinpath(@__DIR__, "sa_common.jl"))
include(joinpath(@__DIR__, "models", "$MODEL.jl"))
using CSV

const YEARS = [2030, 2040, 2050, :mean]

for domain in something(DOMAINS, SA_DOMAINS)
    dir = model_run_dir(MODEL, domain)
    isfile(joinpath(dir, "result_set.jld2")) || (println("$domain: no result_set.jld2 (run 06_collect.jl)"); continue)
    R = load_sa_results(MODEL, domain)
    rs = load_sa_result_set(MODEL, domain)
    out = DataFrame(row=R["rows"])
    for k in sort(collect(keys(rs.outcomes)))
        groups = :groups in outcome_axes(rs.outcomes[k]) ? (1:3) : (nothing,)
        for g in groups, y in YEARS
            name = g === nothing ? "$(k)_$y" : "$(k)_$(lowercase(CORALPOP_MORPHS[g]))_$y"
            out[!, name] = indicator_output(R, rs, k; year=y, group=g)
        end
    end
    CSV.write(joinpath(dir, "indicator_outputs.csv"), out)
    println("$domain: indicator_outputs.csv -- $(nrow(out)) runs x $(ncol(out) - 1) outputs ",
        "($(length(rs.outcomes)) indicators, per-morphology ones split, x $(join(string.(YEARS), ", ")))")
end
