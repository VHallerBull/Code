#=
07_analyse.jl -- figures and a PAWN check for one model's Sobol' runs (not the Shapley analysis,
which is done separately with SAShE.jl on results.jld2).

    julia --project=Code/sa/models/<model> Code/sa/07_analyse.jl <model> [domain ...]

Reads results.jld2 (and control_results.jld2 if present) from 06_collect.jl and writes
<run dir>/figures/*.png, pawn.csv (PAWN median per input and output) and FIGURES.md.
=#

const MODEL = get(ARGS, 1, "")
isempty(MODEL) && error("usage: 07_analyse.jl <model> [domain ...]")
const DOMAINS = length(ARGS) > 1 ? ARGS[2:end] : nothing

include(joinpath(@__DIR__, "sa_common.jl"))
using CSV, CairoMakie
using ADRIAanalysis: sensitivity
import JLD2

const BLOCK_COLOURS = Dict("DHW" => :firebrick, "domain" => :darkorange, "model" => :steelblue, "dummy" => :grey60)
const PAWN_S = 10
## CoralPop reads neither depth nor wave exposure: these multipliers are inert (a noise floor for PAWN)
const DUMMIES = MODEL == "coralpop" ? [:m_depth, :m_ub] : Symbol[]

"PAWN median per input (pawn's Si columns: min, lb, mean, median, ub, max, std, cv)."
pawn_median(X::DataFrame, y) = Array(sensitivity.pawn(Matrix{Float64}(X), Float64.(y), names(X); S=PAWN_S))[:, 4]

function fan!(ax, yrs, S; colour=:steelblue, label="")
    q(p) = [quantile(view(S, t, :), p) for t in axes(S, 1)]
    band!(ax, yrs, q(0.025), q(0.975); color=(colour, 0.18))
    band!(ax, yrs, q(0.25), q(0.75); color=(colour, 0.4))
    lines!(ax, yrs, q(0.5); color=colour, linewidth=2, label=label)
end

function analyse(domain)
    dir = model_run_dir(MODEL, domain)
    path = joinpath(dir, "results.jld2")
    isfile(path) || (println("$domain: no results.jld2 yet (run 06_collect.jl)"); return)
    R = JLD2.load(path)
    X, Y, S, yrs = R["X"], R["Y"], R["series"], R["years"]
    ctrl_path = joinpath(dir, "control_results.jld2")
    C = isfile(ctrl_path) ? JLD2.load(ctrl_path) : nothing
    fig_dir = mkpath(joinpath(dir, "figures"))
    figs = Tuple{String,String}[]
    savefig(name, fig, caption) = (save(joinpath(fig_dir, name), fig; px_per_unit=2); push!(figs, (name, caption)))
    label = "$(R["settings"]["model_label"]) / $domain (N = $(nrow(X)))"

    ## 1. relative-cover ensemble, with the zero-DHW control
    fig = Figure(; size=(850, 450))
    ax = Axis(fig[1, 1]; xlabel="Year", ylabel="Relative cover (domain)", title=label)
    fan!(ax, yrs, S; label="all runs")
    C === nothing || fan!(ax, yrs, C["series"]; colour=:seagreen, label="zero-DHW control")
    axislegend(ax; position=:lt)
    savefig("cover_fan.png", fig, "Domain relative cover: median, 50 % and 95 % bands; green = the same parameter sets with zero DHW")

    ## 2. by DHW warming-trend tercile
    tr = X.dhw_trend_per_decade
    cuts = quantile(tr, [1 / 3, 2 / 3])
    fig = Figure(; size=(850, 450))
    ax = Axis(fig[1, 1]; xlabel="Year", ylabel="Relative cover (domain)", title="$label: by DHW trend tercile")
    for (i, (sel, col, nm)) in enumerate([(tr .<= cuts[1], :seagreen, "low trend"),
        (cuts[1] .< tr .<= cuts[2], :darkorange, "mid trend"), (tr .> cuts[2], :firebrick, "high trend")])
        fan!(ax, yrs, S[:, sel]; colour=col, label=nm)
    end
    axislegend(ax; position=:lt)
    savefig("cover_fan_by_trend.png", fig, "Domain relative cover split by tercile of the DHW warming trend")

    ## 3. output distributions
    outs = [c for c in names(Y) if c != "n_conn_capped"]
    fig = Figure(; size=(300 * min(length(outs), 4), 260 * cld(length(outs), 4)))
    for (i, c) in enumerate(outs)
        ax = Axis(fig[fld(i - 1, 4) + 1, mod(i - 1, 4) + 1]; title=c)
        hist!(ax, Y[!, c]; bins=40, color=(:steelblue, 0.7))
    end
    savefig("outputs_hist.png", fig, "Distribution of each output over all runs")

    ## 4. PAWN check: every input against every output. ADRIA's PAWN is a scaled KS statistic (not bounded
    ## by 1), so a random dummy input is added to show the noise floor.
    Random.seed!(99)
    Xd = hcat(X, DataFrame(random_dummy=rand(nrow(X))))
    P = DataFrame(input=names(Xd), block=[n == "random_dummy" ? "dummy" : factor_block(n) for n in names(Xd)])
    for c in outs
        P[!, c] = pawn_median(Xd, Y[!, c])
    end
    CSV.write(joinpath(dir, "pawn.csv"), P)
    floor_ = P.rc_final[P.input .== "random_dummy"][1]
    dummies = filter(d -> String(d) in names(X), DUMMIES)       # a fixed multiplier is not an input
    inert = isempty(dummies) ? "" : " Inert inputs ($(join(dummies, ", "))) should sit near it."
    P = P[P.input .!= "random_dummy", :]
    top = sort(P, :rc_final; rev=true)[1:min(30, nrow(P)), :]
    fig = Figure(; size=(800, 750))
    ax = Axis(fig[1, 1]; yticks=(1:nrow(top), top.input), yreversed=true, xlabel="PAWN (median), rc_final",
        title="$label: top $(nrow(top)) inputs (check only)")
    barplot!(ax, 1:nrow(top), top.rc_final; direction=:x, color=[BLOCK_COLOURS[b] for b in top.block])
    vlines!(ax, [floor_]; color=:black, linestyle=:dash)
    Legend(fig[1, 2], [PolyElement(; color=BLOCK_COLOURS[b]) for b in ("DHW", "domain", "model")],
        ["DHW", "domain", "model"])
    savefig("pawn_top30.png", fig, "PAWN of final-year relative cover, top 30 inputs; dashed = a random " *
        "dummy input, the noise floor.$inert")

    ## 5. PAWN of every output, top 25 inputs by their maximum
    P.maxval = [maximum(P[i, outs]) for i in 1:nrow(P)]
    T = sort(P, :maxval; rev=true)[1:min(25, nrow(P)), :]
    fig = Figure(; size=(950, 750))
    ax = Axis(fig[1, 1]; xticks=(1:length(outs), outs), yticks=(1:nrow(T), T.input), yreversed=true,
        xticklabelrotation=pi / 4, title="$label: PAWN per output")
    hm = heatmap!(ax, 1:length(outs), 1:nrow(T), Matrix(T[:, outs])'; colormap=:viridis)
    Colorbar(fig[1, 2], hm; label="PAWN (median; random-dummy floor ≈ $(round(floor_; digits=2)))")
    savefig("pawn_heatmap.png", fig, "PAWN per output for the 25 inputs with the largest value on any output")

    ## 6. convergence over nested Sobol' prefixes
    ns = filter(n -> n >= 256, [nrow(X) ÷ 2^k for k in 3:-1:0])
    length(ns) >= 2 || @goto done   # too few runs (a test run)
    top10 = Symbol.(top.input[1:min(10, nrow(top))])
    conv = [pawn_median(X[1:n, :], Y.rc_final[1:n]) for n in ns]
    fig = Figure(; size=(800, 450))
    ax = Axis(fig[1, 1]; xlabel="runs (nested Sobol' prefix)", ylabel="PAWN, rc_final", xscale=log2,
        xticks=(ns, string.(ns)), title="$label: convergence of the top 10")
    idx = indexin(top10, Symbol.(names(X)))
    for (j, f) in enumerate(top10)
        lines!(ax, ns, [c[idx[j]] for c in conv]; label=string(f))
    end
    Legend(fig[1, 2], ax)
    savefig("pawn_convergence.png", fig, "PAWN of rc_final for the top 10 inputs on nested prefixes of the design")

    @label done
    open(joinpath(dir, "FIGURES.md"), "w") do io
        println(io, "# $label\n")
        println(io, "$(nrow(X)) runs, $(ncol(X)) inputs (the DHW member replaced by its 7 generator parameters). ",
            "Inputs/outputs for the Shapley analysis: `results.jld2` (`X`, `Y`). PAWN here is a check only.\n")
        for (f, c) in figs
            println(io, "## $f\n\n$c\n\n![](figures/$f)\n")
        end
    end
    println("$domain: $(length(figs)) figures, pawn.csv and FIGURES.md in $dir")
end

foreach(analyse, something(DOMAINS, SA_DOMAINS))
