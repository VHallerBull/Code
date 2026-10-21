#=
06_figures.jl -- the final figures of a pick-and-freeze run, per input and per group (the Shapley
counterparts of Code/sa/07_analyse.jl's PAWN figures, plus time and input-band comparisons).

    julia +release --project=Code/sa/shapley Code/sa_pf/06_figures.jl <model> [year ...]

Reads <ModelFolder>/<CMR_PF_SUBDIR, default model_run_pf_pilot>/analysis/effects.csv and group_effects.csv (05_analyse.jl), and the
runs' outputs in <band>/<domain>/pf_results.jld2 (04_collect.jl), and
writes analysis/final_figures/*.png and FIGURES.md. `year`s are the snapshot years of the bar charts
(default 2030 and 2050). Shares are effect / Var(Y) of the year; whiskers are 95 % intervals. In the
heatmaps, cells whose interval includes 0 are blank: not resolved at this number of base samples.

Per domain, and for each design set (`main`; `zero` = no heat stress, file names ending in _zero):
  outputs_time                   the outputs themselves: median, 50 % and 95 % of runs per year, per band, with
                                 the zero-DHW reference behind them (main set only)
  inputs_top30_<output>_<year>   top 30 inputs, one column per input band (like pawn_top30)
  inputs_outputs_<year>          top 25 inputs x outputs, one panel per band (like pawn_heatmap)
  inputs_time_<output>           top 30 inputs x years, one panel per band
  groups_time_<output>           group shares over time, stacked, one panel per band
  groups_bands_<year>            every group, all bands side by side, one panel per output
  within_groups_<output>_<year>  the leading members of each group, one figure per band
=#

using CSV, DataFrames, Statistics, CairoMakie
import JLD2
include(joinpath(@__DIR__, "sa_outputs.jl"))                # output names, order and labels

const MODEL = get(ARGS, 1, "coralpop")
const SNAP_YEARS = length(ARGS) > 1 ? parse.(Int, ARGS[2:end]) : [2030, 2050]
const CMR = normpath(joinpath(@__DIR__, "..", ".."))
const FOLDER = Dict("coralpop" => "CoralPOP")[MODEL]
"Root of the model folders (as sa_common.jl's SA_RUNS_ROOT: `CMR_SA_RUNS_DIR` overrides it)."
const RUNS_ROOT = normpath(get(ENV, "CMR_SA_RUNS_DIR", CMR))
"Run folder inside the model folder (as pf_common.jl's PF_SUBDIR: `CMR_PF_SUBDIR`)."
const PF_SUBDIR = get(ENV, "CMR_PF_SUBDIR", "model_run_pf_pilot")
const ANALYSIS = joinpath(RUNS_ROOT, FOLDER, PF_SUBDIR, "analysis")
const OUT = mkpath(joinpath(ANALYSIS, "final_figures"))
const LABEL = Dict("coralpop" => "CoralPop")[MODEL]

## ---- palette (dataviz reference instance, light mode) ----
const CATS = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300", "#6250d6", "#e34948"]
const OTHER, INK, INK2, GRID, AXIS, SURF = "#898781", "#0b0b0b", "#52514e", "#e1e0d9", "#c3c2b7", "#fcfcfb"
"Heatmap cells whose interval includes 0 are left blank (surface); any blue means resolved."
const UNRESOLVED = SURF
"Factor blocks: the first three categorical slots (validated as a set)."
const BLOCK_COL = Dict("DHW" => CATS[2], "domain" => CATS[3], "model" => CATS[1])
"Input bands are ordered: one hue, light to dark (ordinal ramp, no lighter than step 250)."
const BAND_COL = ["#86b6ef", "#3987e5", "#184f95"]
const SEQ = cgrad(["#cde2fb", "#86b6ef", "#3987e5", "#1c5cab", "#0d366b"])
const MAX_SERIES = 7

set_theme!(Theme(; fontsize=13, backgroundcolor=SURF, textcolor=INK,
    Axis=(; backgroundcolor=SURF, xgridcolor=GRID, ygridcolor=GRID, topspinevisible=false, rightspinevisible=false,
        leftspinecolor=AXIS, bottomspinecolor=AXIS, xtickcolor=AXIS, ytickcolor=AXIS, xticklabelcolor=INK2,
        yticklabelcolor=INK2, xlabelcolor=INK2, ylabelcolor=INK2, titlealign=:left, titlefont=:regular, titlesize=14),
    Legend=(; framevisible=false, labelcolor=INK2), Lines=(; linewidth=2)))

block(name) = startswith(name, "dhw_") ? "DHW" : startswith(name, "m_") ? "domain" : "model"
bandlab(b) = "±$(round(Int, 100b)) %"

E = CSV.read(joinpath(ANALYSIS, "effects.csv"), DataFrame)
Gr = CSV.read(joinpath(ANALYSIS, "group_effects.csv"), DataFrame)
for D in (E, Gr)
    "set" in names(D) || (D.set = fill("main", nrow(D)))          # tables written before the zero set existed
    D.set = String.(D.set)
end
"Design sets present: `main` (with the DHW inputs) and, if run, `zero` (no heat stress)."
const SETS = [s for s in ("main", "zero") if s in E.set]
"File-name suffix and title addition of the set being drawn (set in the loop below)."
const SFX, STT = Ref(""), Ref("")
for D in (E, Gr)
    D.slb = D.lb ./ D.var_y
    D.sub = D.ub ./ D.var_y
    D.resolved = D.slb .> 0                                 # 95 % interval above 0
end
const BANDS = sort(unique(E.band))
const DOMAINS = sort(unique(E.domain); rev=true)            # SingleReef first
"Every output analysed (sa_outputs.jl order); each domain draws the ones it has (spatial_cv: NineReefs)."
const OUTPUTS = output_order(unique(E.output))
const NBASE = Dict(b => 256 for b in BANDS)                 # pilot; read from noise.csv below
noise = CSV.read(joinpath(ANALYSIS, "noise.csv"), DataFrame)
for g in groupby(noise, :band); NBASE[g.band[1]] = g.n_base[1]; end

figs = Tuple{String,String}[]
save_(name, fig, cap) = (save(joinpath(OUT, name), fig; px_per_unit=2); push!(figs, (name, STT[] == "" ? cap : "Zero DHW (no heat stress, no DHW inputs): " * cap)))
title_(dom, txt) = "$LABEL / $dom$(STT[]): $txt"

"Top-k rows of `D` (one row per player) by the largest share in any panel, ranked."
function top_players(D, k; col=:share)
    pk = combine(groupby(D, :player), col => (v -> maximum(v)) => :pk)
    return first(sort(pk, :pk; rev=true), k).player
end

E.player = E.input
Gr.player = Gr.group

## period aggregates (05_analyse.jl): variance-weighted shares per period and the time-mean outputs
const HAS_PERIODS = isfile(joinpath(ANALYSIS, "period_effects.csv"))
P, Pg = if HAS_PERIODS
    p, pg = (CSV.read(joinpath(ANALYSIS, f), DataFrame) for f in ("period_effects.csv", "period_group_effects.csv"))
    p.player, pg.player = p.input, pg.group
    for D in (p, pg)
        D.slb, D.sub = D.lb ./ D.var_y, D.ub ./ D.var_y
    end
    p, pg
else
    DataFrame(), DataFrame()
end
"Period series of the period figures: the two sub-periods (light, dark blue), the whole period (ink), the time mean (orange)."
const PERIOD_COL = Dict("2025-2035" => BAND_COL[1], "2036-2050" => BAND_COL[3], "2025-2050" => INK2)
const MEAN_COL = CATS[2]
period_label(p) = replace(p, "-" => "–")

"Zero-DHW reference colour (categorical slot 6, as the green control band of Code/sa/07_analyse.jl)."
const ZERO_COL = CATS[6]

"Median line with 50 % and 95 % bands of `S` ([runs, years]) over `yrs`."
function fan!(ax, yrs, S, col)
    ## a year left out of an output (missing values, e.g. a thermal group without weight) is a gap
    q(p) = [(v = filter(!isnan, view(S, :, t)); isempty(v) ? NaN : quantile(v, p)) for t in axes(S, 2)]
    band!(ax, yrs, q(0.025), q(0.975); color=(col, 0.18))
    band!(ax, yrs, q(0.25), q(0.75); color=(col, 0.38))
    lines!(ax, yrs, q(0.5); color=col)
end

"""
Output trajectories of one domain: median, 50 % and 95 % bands per year of every output and of the
domain-mean DHW, one column per input band (rows share a y-axis so the spread compares across bands).
Only the A and B runs of each walk are used: 2N independent draws from the input distribution.
"""
function output_fans(dom)
    R, Z = Dict{Float64,Any}(), Dict{Float64,Any}()
    for b in BANDS
        dir = joinpath(RUNS_ROOT, FOLDER, PF_SUBDIR, "band$(lpad(round(Int, 100b), 2, '0'))", dom)
        isfile(joinpath(dir, "pf_results.jld2")) && (R[b] = JLD2.load(joinpath(dir, "pf_results.jld2")))
        if isfile(joinpath(dir, "pf_zero.jld2"))                    # 03b_run_zero.jl: A and B sets, no DHW
            Z[b] = JLD2.load(joinpath(dir, "pf_zero.jld2"))
        elseif isfile(joinpath(dir, "pf_zero_dhw", "pf_results.jld2"))  # zero design: its walks' end points
            zr = JLD2.load(joinpath(dir, "pf_zero_dhw", "pf_results.jld2"))
            Z[b] = Dict("Y" => Dict(k => v[:, [1, size(v, 2)], :] for (k, v) in zr["Y"]))
        end
    end
    isempty(R) && return
    Ys = first(values(R))["Y"]
    rows = [[(o, output_label(o), CATS[1]) for o in OUTPUTS if haskey(Ys, o)]; [("dhw", "DHW (domain mean)", BLOCK_COL["DHW"])]]
    fig = Figure(; size=(max(330 * length(BANDS) + 60, 900), 210 * length(rows) + 90))
    for (i, (key, lab, col)) in enumerate(rows)
        axs = Axis[]
        for (j, b) in enumerate(BANDS)
            haskey(R, b) || continue
            r = R[b]
            A = key == "dhw" ? r["dhw"] : r["Y"][key]                       # [N, d + 1, T]
            ok = findall(n -> !any(isnan, @view A[n, :, :]), axes(A, 1))
            S = vcat(A[ok, 1, :], A[ok, end, :])                            # A and B runs: [2N, T]
            yrs = r["years"]
            ax = Axis(fig[i, j]; ylabel=j == 1 ? lab : "", xlabel=i == length(rows) ? "Year" : "",
                title=i == 1 ? bandlab(b) : "")
            if key != "dhw" && haskey(Z, b) && haskey(Z[b]["Y"], key)         # zero-DHW reference, behind
                ZA = Z[b]["Y"][key]
                zok = findall(n -> !any(isnan, @view ZA[n, :, :]), axes(ZA, 1))
                fan!(ax, yrs, vcat(ZA[zok, 1, :], ZA[zok, 2, :]), ZERO_COL)
            end
            fan!(ax, yrs, S, col)
            push!(axs, ax)
        end
        length(axs) > 1 && linkyaxes!(axs...)
    end
    entries = [LineElement(; color=INK2, linewidth=2), PolyElement(; color=(INK2, 0.38)), PolyElement(; color=(INK2, 0.18)),
        PolyElement(; color=CATS[1]), PolyElement(; color=BLOCK_COL["DHW"])]
    labels = ["median", "50 % of runs", "95 % of runs", "runs", "DHW"]
    if !isempty(Z)
        push!(entries, PolyElement(; color=ZERO_COL))
        push!(labels, "same parameter sets, zero DHW")
    end
    Legend(fig[length(rows) + 1, :], entries, labels; orientation=:horizontal)
    nb = 2 * NBASE[BANDS[1]]
    Label(fig[0, :], title_(dom, "output trajectories by input band ($nb independent runs per band)"); halign=:left, fontsize=14)
    save_("outputs_time_$dom.png", fig, "$dom: every output and the domain-mean DHW per year: median, 50 % and 95 % " *
        "of runs, one column per input band (rows share a y-axis). From the A and B runs of the walks, " *
        "$nb independent draws per band. Green: the same parameter sets with zero DHW (reference; not part of " *
        "the sensitivity analysis).")
end

for dom in DOMAINS, set in SETS
    SFX[], STT[] = set == "main" ? ("", "") : ("_$set", " (zero DHW)")
    set == "main" && output_fans(dom)
    Ed = E[(E.domain .== dom) .& (E.set .== set), :]
    Gd = Gr[(Gr.domain .== dom) .& (Gr.set .== set), :]
    nrow(Ed) == 0 && continue
    OUTS = [o for o in OUTPUTS if o in Ed.output]          # this domain's outputs

    ## 1. top 30 inputs, one column per band (bars + 95 % whiskers, coloured by block)
    for out in OUTS, yr in SNAP_YEARS
        S = Ed[(Ed.output .== out) .& (Ed.year .== yr), :]
        nrow(S) == 0 && continue
        players = top_players(S, 30)
        fig = Figure(; size=(330 * length(BANDS) + 330, 760))
        xmax = maximum(S.sub[in.(S.player, Ref(players))]) * 1.05
        xmin = min(0.0, minimum(S.slb[in.(S.player, Ref(players))]))
        for (j, b) in enumerate(BANDS)
            T = S[S.band .== b, :]
            T = T[indexin(players, T.player), :]
            ax = Axis(fig[1, j]; yticks=(1:30, j == 1 ? players : fill("", 30)), yreversed=true,
                xlabel="Share of Var(Y)", title=bandlab(b), limits=((xmin, xmax), nothing), ygridvisible=false)
            barplot!(ax, 1:30, T.share; direction=:x, color=[BLOCK_COL[block(p)] for p in T.player], gap=0.25)
            errorbars!(ax, T.share, 1:30, T.share .- T.slb, T.sub .- T.share; direction=:x, color=INK2, whiskerwidth=5)
            vlines!(ax, [0.0]; color=AXIS)
        end
        Legend(fig[1, length(BANDS) + 1], [PolyElement(; color=BLOCK_COL[k]) for k in ("DHW", "domain", "model")],
            ["DHW", "domain", "model"])
        Label(fig[0, :], title_(dom, "$(output_label(out)) in $yr, top 30 inputs (N = $(NBASE[BANDS[1]]) base samples)");
            halign=:left, fontsize=14)
        save_("inputs_top30_$(out)_$(yr)_$dom$(SFX[]).png", fig, "$dom, $(output_label(out)) in $yr: Shapley share of " *
            "Var(Y) of the 30 inputs largest in any band, with 95 % intervals; one column per input band, same order.")
    end

    ## 2. top 25 inputs x outputs, one panel per band (unresolved cells grey)
    for yr in SNAP_YEARS
        S = Ed[Ed.year .== yr, :]
        players = top_players(S, 25)
        cmax = maximum(S.share[in.(S.player, Ref(players))])
        fig = Figure(; size=(260 * length(BANDS) + 330, 760))
        local hm
        for (j, b) in enumerate(BANDS)
            M = fill(NaN, length(OUTS), 25)
            for (i, o) in enumerate(OUTS), (k, p) in enumerate(players)
                r = S[(S.band .== b) .& (S.output .== o) .& (S.player .== p), :]
                nrow(r) == 1 && r.resolved[1] && (M[i, k] = r.share[1])
            end
            ax = Axis(fig[1, j]; xticks=(1:length(OUTS), [output_label(o) for o in OUTS]), xticklabelrotation=pi / 4,
                yticks=(1:25, j == 1 ? players : fill("", 25)), yreversed=true, title=bandlab(b))
            hm = heatmap!(ax, 1:length(OUTS), 1:25, M; colormap=SEQ, colorrange=(0, cmax), nan_color=UNRESOLVED)
        end
        Colorbar(fig[1, length(BANDS) + 1], hm; label="Share of Var(Y) (blank: interval includes 0)")
        Label(fig[0, :], title_(dom, "top 25 inputs per output, $yr"); halign=:left, fontsize=14)
        save_("inputs_outputs_$(yr)_$dom$(SFX[]).png", fig, "$dom, $yr: share of Var(Y) of the 25 inputs largest on any " *
            "output and band; blank cells are not resolved (95 % interval includes 0).")
    end

    ## 3. top 30 inputs x years, one panel per band
    for out in OUTS
        S = Ed[Ed.output .== out, :]
        S26 = S[S.year .> minimum(S.year), :]                   # 2025 is the initial state: rank/scale without it
        players = top_players(S26, 30)
        yrs = sort(unique(S.year))
        cmax = maximum(S26.share[in.(S26.player, Ref(players))])
        fig = Figure(; size=(330 * length(BANDS) + 330, 760))
        local hm
        for (j, b) in enumerate(BANDS)
            M = fill(NaN, length(yrs), 30)
            T = S[S.band .== b, :]
            for r in eachrow(T)
                k = findfirst(==(r.player), players)
                (k === nothing || !r.resolved) && continue
                M[findfirst(==(r.year), yrs), k] = r.share
            end
            ax = Axis(fig[1, j]; yticks=(1:30, j == 1 ? players : fill("", 30)), yreversed=true, title=bandlab(b),
                xlabel="Year")
            hm = heatmap!(ax, yrs, 1:30, M; colormap=SEQ, colorrange=(0, cmax), nan_color=UNRESOLVED, highclip=SEQ[end])
        end
        Colorbar(fig[1, length(BANDS) + 1], hm; label="Share of Var(Y) (blank: interval includes 0)")
        Label(fig[0, :], title_(dom, "$(output_label(out)), top 30 inputs over time"); halign=:left, fontsize=14)
        save_("inputs_time_$(out)_$dom$(SFX[]).png", fig, "$dom, $(output_label(out)): share of Var(Y) per year of the 30 " *
            "inputs largest in any year and band (ranked and scaled on 2026-2050; 2025, the initial state, saturates); blank = not resolved.")
    end

    ## 4. group shares over time, stacked (top 7 groups + Other), one panel per band
    for out in OUTS
        S = Gd[Gd.output .== out, :]
        groups = top_players(S, MAX_SERIES)
        colors = Dict(g => CATS[i] for (i, g) in enumerate(groups))
        fig = Figure(; size=(380 * length(BANDS) + 280, 420))
        for (j, b) in enumerate(BANDS)
            T = S[S.band .== b, :]
            yrs = sort(unique(T.year))
            val(g, y) = max(0.0, sum(T.share[(T.player .== g) .& (T.year .== y)]))
            layers = [[val(g, y) for y in yrs] for g in groups]
            other = [max(0.0, sum(T.share[(T.year .== y) .& .!in.(T.player, Ref(groups))])) for y in yrs]
            push!(layers, other)
            ax = Axis(fig[1, j]; xlabel="Year", ylabel=j == 1 ? "Share of Var(Y)" : "", title=bandlab(b),
                limits=(nothing, (0, nothing)))
            base = zeros(length(yrs))
            for (k, l) in enumerate(layers)
                c = k <= length(groups) ? colors[groups[k]] : OTHER
                band!(ax, yrs, base, base .+ l; color=c)
                lines!(ax, yrs, base .+ l; color=SURF, linewidth=1)     # 1 px surface gap between layers
                base = base .+ l
            end
            tot = [sum(T.share[T.year .== y]) for y in yrs]
            lines!(ax, yrs, tot; color=INK, linestyle=:dash, linewidth=1.5)
        end
        Legend(fig[1, length(BANDS) + 1], [[PolyElement(; color=colors[g]) for g in groups]; PolyElement(; color=OTHER);
            LineElement(; color=INK, linestyle=:dash)], [groups; "Other ($(length(unique(S.player)) - length(groups)))";
            "sum of all effects"])
        Label(fig[0, :], title_(dom, "$(output_label(out)), group shares over time ('process' groups)"); halign=:left, fontsize=14)
        save_("groups_time_$(out)_$dom$(SFX[]).png", fig, "$dom, $(output_label(out)): share of Var(Y) per group (sum of its " *
            "members' effects), stacked per year; dashed = sum of all effects (≈ 1 if converged). Small negative " *
            "estimates are drawn as 0.")
    end

    ## 5. every group, all bands side by side, one panel per output (dots + 95 % whiskers)
    for yr in SNAP_YEARS
        S = Gd[Gd.year .== yr, :]
        order = combine(groupby(S[S.output .== OUTS[1], :], :player), :share => maximum => :pk)
        groups = sort(order, :pk; rev=true).player
        ng = length(groups)
        fig = Figure(; size=(300 * length(OUTS) + 400, 560))
        for (i, out) in enumerate(OUTS)
            ax = Axis(fig[1, i]; yticks=(1:ng, i == 1 ? groups : fill("", ng)), yreversed=true, xlabel="Share of Var(Y)",
                title=output_label(out), ygridvisible=false)
            vlines!(ax, [0.0]; color=AXIS)
            for (j, b) in enumerate(BANDS)
                T = S[(S.output .== out) .& (S.band .== b), :]
                nrow(T) == 0 && continue                            # year left out of this output
                T = T[indexin(groups, T.player), :]
                y = (1:ng) .+ (j - (length(BANDS) + 1) / 2) * 0.22
                errorbars!(ax, T.share, y, T.share .- T.slb, T.sub .- T.share; direction=:x, color=BAND_COL[j], whiskerwidth=4)
                scatter!(ax, T.share, y; color=BAND_COL[j], markersize=9, strokecolor=SURF, strokewidth=1.5)
            end
        end
        Legend(fig[1, length(OUTS) + 1], [MarkerElement(; marker=:circle, color=c, markersize=10) for c in BAND_COL[1:length(BANDS)]],
            bandlab.(BANDS); title="Input band", titlefont=:regular)
        Label(fig[0, :], title_(dom, "group shares in $yr by input band ('process' groups)"); halign=:left, fontsize=14)
        save_("groups_bands_$(yr)_$dom$(SFX[]).png", fig, "$dom, $yr: share of Var(Y) per group with 95 % interval, for each " *
            "input band; groups ordered by their largest share on $(output_label(OUTS[1])).")
    end

    ## 6. within groups: the leading members of each group (one figure per band)
    for out in OUTS[1:1], yr in SNAP_YEARS[end:end], (j, b) in enumerate(BANDS)
        S = Ed[(Ed.output .== out) .& (Ed.year .== yr) .& (Ed.band .== b), :]
        gs = combine(groupby(S, :group), :share => sum => :tot, nrow => :n)
        gs = sort(gs[gs.n .> 1, :], :tot; rev=true)
        nc = 4
        nr = cld(nrow(gs), nc)
        fig = Figure(; size=(300 * nc, 230 * nr + 60))
        for (k, g) in enumerate(gs.group)
            T = sort(S[S.group .== g, :], :share; rev=true)
            T = first(T, 6)
            ax = Axis(fig[fld(k - 1, nc) + 1, mod(k - 1, nc) + 1]; yticks=(1:nrow(T), T.player), yreversed=true,
                title="$g ($(round(gs.tot[k]; digits=2)))", titlesize=12, yticklabelsize=10,
                xticks=WilkinsonTicks(3; k_min=2), xticklabelsize=10,
                limits=((min(0, minimum(T.slb)), max(1e-6, maximum(T.sub)) * 1.08), (0.4, 6.6)), ygridvisible=false)
            barplot!(ax, 1:nrow(T), T.share; direction=:x, color=BLOCK_COL[block(T.player[1])], gap=0.3)
            errorbars!(ax, T.share, 1:nrow(T), T.share .- T.slb, T.sub .- T.share; direction=:x, color=INK2, whiskerwidth=4)
            vlines!(ax, [0.0]; color=AXIS)
        end
        Label(fig[0, :], title_(dom, "$(output_label(out)) in $yr, leading members of each group, $(bandlab(b)) " *
            "(group share in brackets)"); halign=:left, fontsize=14)
        save_("within_groups_$(out)_$(yr)_band$(round(Int, 100b))_$dom$(SFX[]).png", fig, "$dom, $(output_label(out)) in $yr, " *
            "$(bandlab(b)): the 6 largest members of each group (share of Var(Y), 95 % intervals; each panel on its own " *
            "scale, to show the ranking within the group); the group's total share in brackets. Shares from the main design: interactions with inputs outside the group are split.")
    end

    ## 7. period aggregates: groups and top inputs, one series per period (+ the time-mean output), per band
    HAS_PERIODS || continue
    Pd = P[(P.domain .== dom) .& (P.set .== set), :]
    Pgd = Pg[(Pg.domain .== dom) .& (Pg.set .== set), :]
    for out in OUTS, b in BANDS
        sel(D) = D[(D.band .== b) .& (((D.output .== out) .& (D.kind .== "variance-weighted")) .| (D.output .== "mean_$out")), :]
        Si, Sg = sel(Pd), sel(Pgd)
        nrow(Sg) == 0 && continue
        ## series: (period, legend label, colour, is the time-mean output)
        series = [(p, "$(period_label(p)), variance-weighted", PERIOD_COL[p], false)
                  for p in ("2025-2035", "2036-2050", "2025-2050") if p in Sg.period[Sg.output .== out]]
        mean_p = unique(Sg.period[Sg.output .== "mean_$out"])
        isempty(mean_p) || push!(series, (only(mean_p), "mean $(lowercase(output_label(out))) $(period_label(only(mean_p)))",
            MEAN_COL, true))
        pick(D, s) = s[4] ? D[D.output .== "mean_$out", :] : D[(D.output .== out) .& (D.period .== s[1]), :]
        whole = Sg[(Sg.output .== out) .& (Sg.period .== "2025-2050"), :]
        groups = sort(whole, :share; rev=true).player
        inputs = let pk = combine(groupby(Si, :player), :share => maximum => :pk)
            first(sort(pk, :pk; rev=true), 20).player
        end
        ns = length(series)
        fig = Figure(; size=(1250, 40 * max(length(groups), 20) * 0.8 + 180))
        for (k, (players, D, ttl)) in enumerate(((groups, Sg, "Groups ('process')"), (inputs, Si, "Top 20 inputs")))
            np_ = length(players)
            ax = Axis(fig[1, k]; yticks=(1:np_, players), yreversed=true, xlabel="Share of Var(Y)", title=ttl,
                ygridvisible=false)
            vlines!(ax, [0.0]; color=AXIS)
            for (j, s) in enumerate(series)
                T = pick(D, s)
                T = T[indexin(players, T.player), :]
                y = (1:np_) .+ (j - (ns + 1) / 2) * (0.8 / ns)
                errorbars!(ax, T.share, y, T.share .- T.slb, T.sub .- T.share; direction=:x, color=s[3], whiskerwidth=3)
                scatter!(ax, T.share, y; color=s[3], markersize=8, strokecolor=SURF, strokewidth=1.2)
            end
        end
        Legend(fig[2, :], [MarkerElement(; marker=:circle, color=s[3], markersize=10) for s in series], [s[2] for s in series];
            orientation=:horizontal, tellwidth=false)
        bsfx = length(BANDS) > 1 ? "_band$(round(Int, 100b))" : ""
        Label(fig[0, :], title_(dom, "$(output_label(out)), Shapley shares over periods" *
            (length(BANDS) > 1 ? ", $(bandlab(b))" : "")); halign=:left, fontsize=14)
        save_("periods_$(out)$(bsfx)_$dom$(SFX[]).png", fig, "$dom, $(output_label(out)): variance-weighted Shapley " *
            "shares per period (sum of the effects over the period's years / sum of Var(Y); years with more variance " *
            "weigh more), with 95 % intervals" * (isempty(mean_p) ? "" : ", and the Shapley shares of the output's " *
            "mean over $(period_label(only(mean_p))) as an output of its own (its variance includes the covariances " *
            "between years)") * ". Groups ordered by the 2025–2050 share; inputs: the 20 largest in any series.")
    end
end

open(joinpath(OUT, "FIGURES.md"), "w") do io
    println(io, "# $LABEL pick-and-freeze: final figures (draft, pilot data)\n")
    println(io, "Shares are Shapley effect / Var(Y); 95 % intervals; blank heatmap cells are not resolved. ",
        "Pilot: $(NBASE[BANDS[1]]) base samples per band and domain.\n")
    for (f, c) in figs
        println(io, "## $f\n\n$c\n\n![]($f)\n")
    end
end
println("$(length(figs)) figures -> $OUT")
