#=
run_shapley.jl -- Shapley effects of one output of a model's Sobol' runs, for every year, ungrouped,
grouped, within groups and along follow-on chains (see Code/sa/shapley/README.md).

    julia +release -t auto --project=Code/sa/shapley Code/sa/shapley/run_shapley.jl <config.toml> [--key=value ...]

`--key=value` overrides a config entry; nested keys use dots and values are TOML, e.g.
`--output.indicator=reef_condition_index`, `--domains='["NineReefs_SobolDHW"]'`, `--estimator.n_rows=128`.

Runs every run set in `sets`: "main" (result_set.jld2) and "control" (control_result_set.jld2, the
zero-DHW runs). Writes, per domain and set, to <ModelFolder>/<run_subdir>/<domain>/shapley/<output tag>/
(main) or .../<output tag>__control/ (control):
    ungrouped.csv, grouped_sum_<scheme>.csv, grouped_player_<scheme>.csv, within_<scheme>.csv,
    chain_<name>.csv      long tables: analysis, scheme, group, player, block, year, effect, lb, ub,
                          var_y, total, share (= effect / var_y), share_of_total (= effect / total)
    convergence.csv       ungrouped and grouped effects on nested prefixes of the runs (column n_rows)
    shapley.jld2          the same tables, the config and run metadata
    figures/, SHAPLEY.md
=#

using TOML

"Config with `--a.b=value` overrides applied (values parsed as TOML, else kept as text)."
function read_config(args)
    isempty(args) && error("usage: run_shapley.jl <config.toml> [--key=value ...]")
    cfg = TOML.parsefile(args[1])
    for a in args[2:end]
        m = match(r"^--([^=]+)=(.*)$", a)
        m === nothing && error("cannot read argument '$a' (expected --key=value)")
        val = try
            TOML.parse("v = " * m[2])["v"]
        catch
            String(m[2])
        end
        keys_ = split(m[1], ".")
        d = cfg
        for k in keys_[1:end-1]
            d = get!(d, k, Dict{String,Any}())
        end
        d[keys_[end]] = val
    end
    return cfg
end

const CFG = read_config(ARGS)
const MODEL = CFG["model"]
ENV["CMR_SA_RUN_SUBDIR"] = get(CFG, "run_subdir", "model_run")   # before sa_common.jl reads it

include(joinpath(@__DIR__, "..", "sa_common.jl"))
include(joinpath(@__DIR__, "..", "models", "$MODEL.jl"))
include(joinpath(@__DIR__, "shapley_core.jl"))
include(joinpath(@__DIR__, "groups.jl"))
include(joinpath(@__DIR__, "output_select.jl"))
include(joinpath(@__DIR__, "plots.jl"))
using CSV
import JLD2

const EST = CFG["estimator"]
const AN = get(CFG, "analyses", Dict{String,Any}())
const N_PERM = Int(EST["n_permutations"])
const SEED = Int(get(EST, "seed", 1234))
const N_ROWS = Int(get(EST, "n_rows", 0))
const TERMINAL = Symbol(get(EST, "terminal", "var_y"))
const GROUPS_FILE = joinpath(@__DIR__, get(CFG, "groups", "groups/$MODEL.toml"))
n_perm_for(analysis) = Int(get(AN, "n_permutations_$analysis", N_PERM))

"Block (DHW / domain / model) of an input, or of a set of inputs when they share one."
block_of(input_names) = (b = unique(factor_block.(input_names)); length(b) == 1 ? b[1] : "mixed")

"""
    long_table(r, years; analysis, scheme="", group="", blocks) -> DataFrame

One row per player and year of the `ShapleyResult` `r`. `blocks[k]` is player k's block.
"""
function long_table(r::ShapleyResult, yrs; analysis, scheme="", group="", blocks)
    sv, st = share_of_variance(r), share_of_total(r)
    rows = [(; analysis, scheme, group, player=r.names[k], block=blocks[k], year=yrs[t],
        effect=r.effect[k, t], lb=r.lb[k, t], ub=r.ub[k, t], var_y=r.var_y[t], total=r.total[t],
        share=sv[k, t], share_of_total=st[k, t]) for t in eachindex(yrs) for k in eachindex(r.names)]
    return DataFrame(rows)
end

"Expand `within` entries (a scheme or \"scheme:group\") to (scheme, group index) pairs."
function within_targets(groupings)
    out = Tuple{String,Int}[]
    for ref in get(AN, "within", String[])
        if occursin(":", ref)
            has_group(groupings, ref) || (println("within '$ref': group not among this set's inputs, skipped"); continue)
            G, g = resolve_group(groupings, ref)
            push!(out, (G.scheme, g))
        else
            haskey(groupings, ref) || error("unknown scheme '$ref' in `within`")
            append!(out, [(ref, g) for g in eachindex(groupings[ref].names)])
        end
    end
    return out
end

"Whether \"scheme:group\" names a group with inputs in `groupings`."
function has_group(groupings, ref)
    scheme, name = split(ref, ":"; limit=2)
    return haskey(groupings, scheme) && name in groupings[scheme].names
end

"""
Run counts of the convergence check: `convergence.n_rows` if given, else the nested Sobol' prefixes
N/8, N/4, N/2, N (at least 32 runs).
"""
function convergence_sizes(N, conv)
    haskey(conv, "n_rows") && return sort([Int(n) for n in conv["n_rows"] if n <= N])
    return unique(filter(>=(32), [N ÷ 2^k for k in 3:-1:0]))
end

"Years of the convergence check: `convergence.years` if given, else first, middle and last."
function convergence_years(yrs, conv)
    haskey(conv, "years") && return [Int(y) for y in conv["years"] if y in yrs]
    return unique([yrs[1], yrs[cld(length(yrs), 2)], yrs[end]])
end

"`M` (`[runs, all years]`) re-indexed to the output years minus `lag`; missing years are NaN columns."
function align_years(M, m_years, y_years, lag)
    out = fill(NaN, size(M, 1), length(y_years))
    for (t, y) in enumerate(y_years)
        k = findfirst(==(y - lag), m_years)
        k === nothing || (out[:, t] = M[:, k])
    end
    return out
end

"Run-set label for titles."
set_label(set) = set == "main" ? "" : set == "control" ? " (zero-DHW control)" : " ($set)"

function run_domain(domain, set="main")
    println("\n== $MODEL / $domain / $set ==")
    t_load = @elapsed begin
        d = RunData(MODEL, domain; set, n_rows=N_ROWS)
        X = inputs(d)
        pats = Regex.(get(get(CFG, set, Dict{String,Any}()), "drop_inputs", String[]))
        drop = [n for n in names(X) if any(r -> occursin(r, n), pats)]
        isempty(drop) || (X = select(X, Not(drop)); println("set '$set': inputs left out: $(join(drop, ", "))"))
        xn = names(X)
        out = output_matrix(CFG["output"], d)
    end
    println("loaded runs and output in $(round(t_load; digits=1)) s")
    tag = out.tag * (d.n_rows < nrow(d.R["X"]) ? "__n$(d.n_rows)" : "") * (set == "main" ? "" : "__$set")
    dir = mkpath(joinpath(d.dir, "shapley", tag))
    G = load_groupings(GROUPS_FILE, xn)
    schemes = String.(get(AN, "grouped", String[]))
    tables = Dict{String,DataFrame}()
    timing = Dict{String,Float64}()
    println("$(d.n_rows) runs x $(length(xn)) inputs; output $(out.label), $(length(out.years)) years -> $dir")

    ## 1. ungrouped (+ grouped sums from the same increments)
    plan = nothing
    if get(AN, "ungrouped", true) || !isempty(schemes) || !isempty(get(CFG, "chain", []))
        nm, pl = singleton_players(xn)
        println("ungrouped: $(length(pl)) players, $N_PERM permutations, $(Threads.nthreads()) threads")
        timing["ungrouped_plan"] = @elapsed plan = neighbour_plan(X, pl; names=nm, n_perm=N_PERM,
            seed=SEED, terminal=TERMINAL, verbose=true)
        println("  neighbour plan: $(round(timing["ungrouped_plan"]; digits=1)) s " *
                "($(round(1000 * timing["ungrouped_plan"] / N_PERM; digits=1)) ms per permutation)")
        sums = Dict(s => (G[s].names, G[s].members) for s in schemes)
        ind, gsum = shapley_from_plan(plan, out.Y; sums)
        tables["ungrouped"] = long_table(ind, out.years; analysis="ungrouped", blocks=factor_block.(xn))
        for s in schemes
            tables["grouped_sum_$s"] = long_table(gsum[s], out.years; analysis="grouped_sum", scheme=s,
                blocks=[block_of(xn[m]) for m in G[s].members])
        end
    end

    ## 2. group-as-player
    for s in schemes
        println("group-as-player '$s': $(length(G[s].names)) players")
        timing["player_$s"] = @elapsed planG = neighbour_plan(X, G[s].members; names=G[s].names,
            n_perm=n_perm_for("players"), seed=SEED, terminal=TERMINAL)
        r, _ = shapley_from_plan(planG, out.Y)
        tables["grouped_player_$s"] = long_table(r, out.years; analysis="grouped_player", scheme=s,
            blocks=[block_of(xn[m]) for m in G[s].members])
    end

    ## 3. within groups (restricted game)
    for (s, g) in within_targets(G)
        cols = G[s].members[g]
        gname = G[s].names[g]
        println("within '$s:$gname': $(length(cols)) players")
        timing["within_$s:$gname"] = @elapsed planW = neighbour_plan(X, [[c] for c in cols];
            names=xn[cols], n_perm=n_perm_for("within"), seed=SEED, terminal=:estimate)
        r, _ = shapley_from_plan(planW, out.Y)
        t = long_table(r, out.years; analysis="within", scheme=s, group=gname, blocks=factor_block.(xn[cols]))
        tables["within_$s"] = vcat(get(tables, "within_$s", DataFrame()), t)
    end

    ## 4. follow-on chains
    for ch in get(CFG, "chain", [])
        name = ch["name"]
        lag = Int(get(ch, "lag", 0))
        if !has_group(G, ch["replaces"])
            println("chain '$name': upstream group $(ch["replaces"]) not among this set's inputs, skipped")
            continue
        end
        mspec = merge(Dict{String,Any}(ch["intermediate"]), Dict{String,Any}("years" => "all"))
        t_m = @elapsed mo = output_matrix(mspec, d)
        println("  intermediate loaded in $(round(t_m; digits=1)) s")
        if all(t -> var(@view mo.Y[:, t]) == 0, axes(mo.Y, 2))
            println("chain '$name': the intermediate does not vary in this set, skipped")
            continue
        end
        Gs, g = resolve_group(G, ch["replaces"])
        upstream = Gs.members[g]
        if get(ch, "players", "groups") == "groups"
            ks = [k for k in eachindex(Gs.names) if k != g]
            others, other_names = Gs.members[ks], Gs.names[ks]
        else
            keep = setdiff(eachindex(xn), upstream)
            others, other_names = [[c] for c in keep], xn[keep]
        end
        m_name = "M: " * mo.label
        println("chain '$name': $(xn[upstream]) -> $(mo.label) (lag $lag) -> $(out.label); " *
                "$(length(others) + 1) stage-B players")
        M = align_years(mo.Y, mo.years, out.years, lag)
        valid = findall(t -> !any(isnan, @view M[:, t]), eachindex(out.years))
        timing["chain_$name"] = @elapsed res = chain_analysis(X, out.Y[:, valid], M[:, valid], upstream,
            others; upstream_names=xn[upstream], other_names, m_name, n_perm=n_perm_for("chain"), seed=SEED)
        yv = out.years[valid]
        tA = long_table(res.A, yv; analysis="chain_A", scheme=Gs.scheme, group=Gs.names[g],
            blocks=factor_block.(xn[upstream]))     # indexed by the output year (M taken at year - lag)
        tB = long_table(res.B, yv; analysis="chain_B", scheme=Gs.scheme,
            blocks=[[block_of(xn[m]) for m in others]; "intermediate"])
        paths = vcat(res.path, permutedims(res.path_unexplained))
        P = ShapleyResult([xn[upstream]; "unexplained part of M"], paths, fill(NaN, size(paths)),
            fill(NaN, size(paths)), res.B.var_y, res.B.total)
        tP = long_table(P, yv; analysis="chain_path", scheme=Gs.scheme, group=Gs.names[g],
            blocks=[factor_block.(xn[upstream]); "intermediate"])
        ## the upstream group as one player on Y, from the main analysis (computed here if not there)
        gp = get(tables, "grouped_player_$(Gs.scheme)", nothing)
        if gp === nothing
            planG = neighbour_plan(X, Gs.members; names=Gs.names, n_perm=n_perm_for("players"), seed=SEED,
                terminal=TERMINAL)
            gp = long_table(first(shapley_from_plan(planG, out.Y)), out.years; analysis="grouped_player",
                scheme=Gs.scheme, blocks=fill("", length(Gs.names)))
        end
        tG = filter(r -> r.player == Gs.names[g] && r.year in yv, gp)
        tG.analysis .= "chain_upstream_group"
        ## the intermediate itself, analysed like any output (reuses the ungrouped neighbour plan)
        Mi, Mg = shapley_from_plan(plan, mo.Y; sums=Dict(s => (G[s].names, G[s].members) for s in schemes))
        tMi = long_table(Mi, mo.years; analysis="intermediate_ungrouped", blocks=factor_block.(xn))
        tMg = [long_table(Mg[s], mo.years; analysis="intermediate_grouped_sum", scheme=s,
            blocks=[block_of(xn[m]) for m in G[s].members]) for s in schemes]
        tables["chain_$name"] = vcat(tA, tB, tP, tG, tMi, tMg...; cols=:union)
        println("  $(round(timing["chain_$name"]; digits=1)) s")
    end

    ## 5. convergence over nested Sobol' prefixes of the runs
    conv = get(CFG, "convergence", Dict{String,Any}())
    if get(conv, "enabled", true)
        ns = convergence_sizes(d.n_rows, conv)
        cy = convergence_years(out.years, conv)
        ti = [findfirst(==(y), out.years) for y in cy]
        println("convergence: runs $(join(ns, ", ")); years $(join(cy, ", "))")
        parts = DataFrame[]
        timing["convergence"] = @elapsed for n in ns
            Xn, Yn = X[1:n, :], out.Y[1:n, ti]
            nm, pl = singleton_players(xn)
            pn = neighbour_plan(Xn, pl; names=nm, n_perm=N_PERM, seed=SEED, terminal=TERMINAL)
            ind, gsum = shapley_from_plan(pn, Yn; sums=Dict(s => (G[s].names, G[s].members) for s in schemes))
            push!(parts, insertcols!(long_table(ind, cy; analysis="ungrouped", blocks=factor_block.(xn)), 1, :n_rows => n))
            for s in schemes
                bl = [block_of(xn[m]) for m in G[s].members]
                push!(parts, insertcols!(long_table(gsum[s], cy; analysis="grouped_sum", scheme=s, blocks=bl), 1, :n_rows => n))
                pg = neighbour_plan(Xn, G[s].members; names=G[s].names, n_perm=n_perm_for("players"), seed=SEED,
                    terminal=TERMINAL)
                push!(parts, insertcols!(long_table(first(shapley_from_plan(pg, Yn)), cy; analysis="grouped_player",
                    scheme=s, blocks=bl), 1, :n_rows => n))
            end
        end
        tables["convergence"] = vcat(parts...)
        println("  $(round(timing["convergence"]; digits=1)) s")
    end

    ## write
    for (k, t) in tables
        CSV.write(joinpath(dir, "$k.csv"), t)
    end
    meta = Dict("model" => MODEL, "domain" => domain, "label" => out.label * set_label(set), "tag" => tag, "set" => set,
        "n_rows" => d.n_rows, "n_inputs" => length(xn), "years" => out.years, "timing" => timing,
        "threads" => Threads.nthreads(), "model_label" => MODEL_LABEL,
        "inert" => inert_inputs(GROUPS_FILE))
    JLD2.save(joinpath(dir, "shapley.jld2"), "tables", tables, "config", CFG, "meta", meta)
    t_fig = @elapsed figs = make_figures(dir, tables, meta)
    println("figures in $(round(t_fig; digits=1)) s")
    write_index(dir, figs, meta)
    println("$domain: $(length(tables)) tables, $(length(figs)) figures in $dir")
end

for domain in CFG["domains"], set in get(CFG, "sets", ["main"])
    run_domain(domain, set)
end
