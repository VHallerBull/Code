#=
plots.jl -- figures of one Shapley run (run_shapley.jl): `make_figures(dir, tables, meta)` writes
figures/*.png and returns (file, caption) pairs; `write_index` writes SHAPLEY.md.

Shares are effects / Var(Y) of the year (or / the group total for within-group panels). Effects from
the given-data estimator can be slightly negative; lines show that honestly instead of hiding it in a
stack. Colours: fixed categorical order (dataviz reference palette), at most 7 series + "Other".
=#

using CairoMakie

const CATS = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300", "#6250d6", "#e34948"]
const OTHER = "#898781"
const INK, INK2, GRID, AXIS = "#0b0b0b", "#52514e", "#e1e0d9", "#c3c2b7"
const MAX_SERIES = 7

set_theme!(Theme(;
    fontsize=13, backgroundcolor="#fcfcfb", textcolor=INK,
    Axis=(; backgroundcolor="#fcfcfb", xgridcolor=GRID, ygridcolor=GRID, topspinevisible=false,
        rightspinevisible=false, leftspinecolor=AXIS, bottomspinecolor=AXIS, xtickcolor=AXIS, ytickcolor=AXIS,
        xticklabelcolor=INK2, yticklabelcolor=INK2, xlabelcolor=INK2, ylabelcolor=INK2,
        titlealign=:left, titlefont=:regular, titlesize=14),
    Legend=(; framevisible=false, labelcolor=INK2, patchsize=(18, 10)),
    Lines=(; linewidth=2),
))

"""
    wide(t, col) -> (players, years, M)

`t` (one analysis' long table) as a `[players, years]` matrix of column `col`, players in table order.
"""
function wide(t::DataFrame, col::Symbol)
    players = unique(t.player)
    yrs = sort(unique(t.year))
    M = fill(NaN, length(players), length(yrs))
    pi_ = Dict(p => i for (i, p) in enumerate(players))
    yi = Dict(y => i for (i, y) in enumerate(yrs))
    for r in eachrow(t)
        M[pi_[r.player], yi[r.year]] = r[col]
    end
    return players, yrs, M
end

_peak(v) = (w = filter(!isnan, v); isempty(w) ? -Inf : maximum(w))

"Indices of the `n` rows of `M` with the largest peak, in their original (table) order."
top_rows(M, n) = sort(partialsortperm([_peak(M[i, :]) for i in axes(M, 1)], 1:min(n, size(M, 1)); rev=true))

"""
    draw_series!(ax, t; col=:share, lo=:lb, hi=:ub, scale=:var_y, n=MAX_SERIES, colors=nothing)

Lines of the top `n` players of long table `t` over the years (bars if there is only one year), with
95 % bands, the rest summed as "Other". `scale` is the column the effects and bounds are divided by.
Returns (labels, colours) for a legend.
"""
function draw_series!(ax, t::DataFrame; scale=:var_y, n=MAX_SERIES, fixed=nothing, other=true)
    s = t[!, scale]
    t = transform(t, [:effect, :lb, :ub] .=> (c -> c ./ s) .=> [:v, :vlo, :vhi])
    players, yrs, V = wide(t, :v)
    _, _, LO = wide(t, :vlo)
    _, _, HI = wide(t, :vhi)
    keep = fixed === nothing ? top_rows(V, n) : [findfirst(==(p), players) for p in fixed if p in players]
    rest = setdiff(axes(V, 1), keep)
    labels, cols = String[], Any[]
    single = length(yrs) == 1
    vals = [V[k, :] for k in keep]
    if other && !isempty(rest)
        push!(vals, vec(sum(V[rest, :]; dims=1)))
    end
    for (j, v) in enumerate(vals)
        is_other = j > length(keep)
        c = is_other ? OTHER : CATS[mod1(j, length(CATS))]
        lab = is_other ? "Other ($(length(rest)))" : players[keep[j]]
        if single
            barplot!(ax, [j], v; direction=:x, color=c, gap=0.25)
            is_other || errorbars!(ax, v, [j], v .- LO[keep[j], :], HI[keep[j], :] .- v;
                direction=:x, color=INK2, whiskerwidth=6)
        else
            is_other || band!(ax, yrs, LO[keep[j], :], HI[keep[j], :]; color=(c, 0.15))
            lines!(ax, yrs, v; color=c, linestyle=is_other ? :dash : :solid)
        end
        push!(labels, lab)
        push!(cols, c)
    end
    if single
        ax.yticks = (1:length(labels), labels)
        ax.yreversed = true
    end
    return labels, cols, single
end

legend_entries(cols, labels) = [LineElement(; color=c, linewidth=3, linestyle=startswith(l, "Other") ? :dash : :solid)
                                for (c, l) in zip(cols, labels)]

function make_figures(dir::AbstractString, tables::Dict{String,DataFrame}, meta)
    fig_dir = mkpath(joinpath(dir, "figures"))
    figs = Tuple{String,String}[]
    save_(name, fig, caption) = (save(joinpath(fig_dir, name), fig; px_per_unit=2); push!(figs, (name, caption)))
    title = "$(meta["model_label"]) / $(meta["domain"]): $(meta["label"]) (N = $(meta["n_rows"]))"
    ylab = "Share of Var(Y)"

    ## 1. variance of the output over time (shares hide when the variance itself is tiny)
    base = something(get(tables, "ungrouped", nothing), first(values(tables)))
    v = combine(groupby(base, :year), :var_y => first => :var_y, :total => first => :total)
    if nrow(v) > 1
        fig = Figure(; size=(800, 320))
        ax = Axis(fig[1, 1]; xlabel="Year", ylabel="Var(Y)", title="$title: output variance")
        lines!(ax, v.year, v.var_y; color=CATS[1])
        if any(abs.(v.total .- v.var_y) .> 1e-12 .* abs.(v.var_y))
            lines!(ax, v.year, v.total; color=CATS[2])
            Legend(fig[1, 2], legend_entries(CATS[1:2], ["Var(Y)", "explained by inputs"]), ["Var(Y)", "explained by inputs"])
        end
        save_("variance.png", fig, "Variance of the output per year: the scale the shares below are relative to")
    end

    ## 2. grouped: sum of members vs group as one player, same colours per group
    for k in sort(collect(keys(tables)))
        startswith(k, "grouped_sum_") || continue
        s = replace(k, "grouped_sum_" => "")
        ts, tp = tables[k], get(tables, "grouped_player_$s", nothing)
        _, _, V = wide(ts, :share)
        fixed = unique(ts.player)[top_rows(V, MAX_SERIES)]
        fig = Figure(; size=(tp === nothing ? 750 : 1250, 430))
        ax1 = Axis(fig[1, 1]; xlabel="Year", ylabel=ylab, title="sum of member effects")
        labels, cols, single = draw_series!(ax1, ts; fixed)
        if tp !== nothing
            ax2 = Axis(fig[1, 2]; xlabel="Year", ylabel=ylab, title="group as one player")
            draw_series!(ax2, tp; fixed)
            single || linkyaxes!(ax1, ax2)
        end
        single || Legend(fig[1, tp === nothing ? 2 : 3], legend_entries(cols, labels), labels)
        Label(fig[0, :], "$title: groups '$s'"; halign=:left, fontsize=14, color=INK)
        save_("grouped_$s.png", fig, "Groups of scheme '$s': share of the output variance per year; left = sum of " *
              "the members' individual effects, right = the group as one player (bands: 95 % CI, anti-conservative)")
    end

    ## 3. individual inputs: top 7 over time, inert inputs as the noise floor; heatmap of the top 30
    if haskey(tables, "ungrouped")
        t = tables["ungrouped"]
        inert = [String(x) for x in get(meta, "inert", String[]) if x in t.player]
        fig = Figure(; size=(900, 450))
        ax = Axis(fig[1, 1]; xlabel="Year", ylabel=ylab, title="$title: top inputs")
        labels, cols, single = draw_series!(ax, t; other=false)
        if !single
            entries = legend_entries(cols, labels)
            for x in inert
                r = sort(filter(r -> r.player == x, t), :year)
                lines!(ax, r.year, r.share; color=OTHER, linestyle=:dot, linewidth=1.5)
            end
            if !isempty(inert)
                push!(entries, LineElement(; color=OTHER, linestyle=:dot))
                push!(labels, "inert: $(join(inert, ", "))")
            end
            Legend(fig[1, 2], entries, labels)
        end
        save_("top_inputs.png", fig, "The $(MAX_SERIES) inputs with the largest share in any year (bands: 95 % CI); " *
              "dotted = inputs inert for this model, an empirical noise floor")

        players, yrs, V = wide(t, :share)
        if length(yrs) > 1
            idx = top_rows(V, 30)
            idx = idx[sortperm([_peak(V[i, :]) for i in idx]; rev=true)]
            fig = Figure(; size=(900, 750))
            ax = Axis(fig[1, 1]; xlabel="Year", yticks=(1:length(idx), players[idx]), yreversed=true,
                title="$title: top 30 inputs")
            hm = heatmap!(ax, yrs, 1:length(idx), permutedims(V[idx, :]);
                colormap=cgrad(["#fcfcfb", "#86b6ef", "#2a78d6", "#104281"]), colorrange=(0, max(1e-12, maximum(filter(!isnan, V[idx, :])))))
            Colorbar(fig[1, 2], hm; label=ylab)
            save_("heatmap_inputs.png", fig, "Share of the output variance per year for the 30 inputs with the largest peak")
        end
    end

    ## 4. within groups (restricted game): one panel per group with more than one member
    for k in sort(collect(keys(tables)))
        startswith(k, "within_") || continue
        s = replace(k, "within_" => "")
        t = tables[k]
        gs = [g for g in unique(t.group) if length(unique(t.player[t.group .== g])) > 1]
        isempty(gs) && continue
        nc = min(3, length(gs))
        nr = cld(length(gs), nc)
        fig = Figure(; size=(430 * nc, 330 * nr))
        for (i, g) in enumerate(gs)
            tg = t[t.group .== g, :]
            closed = combine(groupby(tg, :year), [:total, :var_y] => ((a, b) -> first(a) / first(b)) => :c).c
            sub = length(closed) > 1 ? "explains $(round(Int, 100 * minimum(closed)))–$(round(Int, 100 * maximum(closed))) % of Var(Y)" :
                  "explains $(round(Int, 100 * closed[1])) % of Var(Y)"
            ax = Axis(fig[fld(i - 1, nc) + 1, mod(i - 1, nc) + 1]; xlabel="Year", ylabel=ylab,
                title="$g ($sub)", titlesize=12)
            labels, cols, single = draw_series!(ax, tg; n=4)
            single || axislegend(ax, legend_entries(cols, labels), labels; position=:lt, labelsize=10,
                backgroundcolor=(:white, 0.7))
        end
        Label(fig[0, :], "$title: within groups of '$s' (only the group's inputs as players)"; halign=:left, fontsize=14)
        save_("within_$s.png", fig, "Within each group of '$s', only the group's inputs as players (restricted " *
              "game): each member's effect as a share of Var(Y); they add up to the group's Var(E[Y | group]) / " *
              "Var(Y), whose range over the years is in the title. Top 4 members + Other. (Shares of the group " *
              "total are in the table, `share_of_total`; they are unstable where that total is near 0.)")
    end

    ## 5. chains
    for k in sort(collect(keys(tables)))
        startswith(k, "chain_") || continue
        name = replace(k, "chain_" => "")
        t = tables[k]
        tB = t[t.analysis .== "chain_B", :]
        mname = tB.player[findfirst(==("intermediate"), tB.block)]
        tM = sort(tB[tB.player .== mname, :], :year)
        tG = sort(t[t.analysis .== "chain_upstream_group", :], :year)
        tP, tA = t[t.analysis .== "chain_path", :], t[t.analysis .== "chain_A", :]
        nrow(tM) > 1 || continue
        fig = Figure(; size=(1350, 420))
        ax1 = Axis(fig[1, 1]; xlabel="Year", ylabel=ylab, title="intermediate vs upstream group on Y")
        lines!(ax1, tM.year, tM.share; color=CATS[1], label="M as player (stage B)")
        nrow(tG) > 0 && lines!(ax1, tG.year, tG.share; color=CATS[2], label="upstream group as player")
        axislegend(ax1; position=:lt, labelsize=11)
        ax2 = Axis(fig[1, 2]; xlabel="Year", ylabel=ylab, title="path: upstream input → M → Y")
        l2, c2, _ = draw_series!(ax2, tP; other=false, n=8)
        ax3 = Axis(fig[1, 3]; xlabel="Year", ylabel="Share of Var(M)", title="stage A: upstream inputs on M")
        draw_series!(ax3, tA; other=false, fixed=l2)
        Legend(fig[1, 4], legend_entries(c2, l2), l2)
        Label(fig[0, :], "$title: chain '$name' ($(replace(mname, "M: " => "M = ")))"; halign=:left, fontsize=14)
        save_("chain_$name.png", fig, "Chain '$name'. Left: Y's variance share of the intermediate M (with the upstream " *
              "group replaced by M) vs the upstream group as one player; M above the group = M carries variance " *
              "its parameters don't describe (e.g. weather realisation). Middle: share of Var(Y) through M per " *
              "upstream input. Right: each upstream input's share of Var(M) (restricted game)")
    end
    ## 6. convergence over nested prefixes of the runs
    haskey(tables, "convergence") && append!(figs, convergence_figures(fig_dir, tables["convergence"], meta, title))
    return figs
end

"""
One panel of a convergence figure: the share of Var(Y) of `players` in long table `t` (one year, several
`n_rows`) against the number of runs, with 95 % bands. `colors[k]` and `styles[k]` per player.
"""
function convergence_panel!(ax, t, players, colors; styles=fill(:solid, length(players)))
    for (k, p) in enumerate(players)
        r = sort(t[t.player .== p, :], :n_rows)
        nrow(r) == 0 && continue
        styles[k] == :solid && band!(ax, r.n_rows, r.lb ./ r.var_y, r.ub ./ r.var_y; color=(colors[k], 0.12))
        lines!(ax, r.n_rows, r.share; color=colors[k], linestyle=styles[k], linewidth=styles[k] == :solid ? 2 : 1.5)
        scatter!(ax, r.n_rows, r.share; color=colors[k], markersize=8, strokecolor="#fcfcfb", strokewidth=1.5)
    end
end

"The `n` players with the largest peak share at the largest run count, in rank order."
function top_at_full(t, n)
    full = t[t.n_rows .== maximum(t.n_rows), :]
    pk = combine(groupby(full, :player), :share => (v -> maximum(filter(!isnan, v); init=-Inf)) => :pk)
    return first(sort(pk, :pk; rev=true), n).player
end

function convergence_figures(fig_dir, c::DataFrame, meta, title)
    figs = Tuple{String,String}[]
    ns = sort(unique(c.n_rows))
    length(ns) >= 2 || return figs
    cy = sort(unique(c.year))
    axis_(f, i, j, ttl) = Axis(f[i, j]; xscale=log2, xticks=(ns, string.(ns)), xlabel="Runs (nested Sobol' prefix)",
        ylabel="Share of Var(Y)", title=ttl)

    ## individual inputs: the top 7 at the largest N, plus the inert ones (noise floor)
    cu = c[c.analysis .== "ungrouped", :]
    top = top_at_full(cu, MAX_SERIES)
    inert = [String(x) for x in get(meta, "inert", String[]) if x in cu.player && !(x in top)]
    players = [top; inert]
    cols = [CATS[1:length(top)]; fill(OTHER, length(inert))]
    styles = [fill(:solid, length(top)); fill(:dot, length(inert))]
    fig = Figure(; size=(380 * length(cy) + 260, 420))
    for (j, y) in enumerate(cy)
        convergence_panel!(axis_(fig, 1, j, string(y)), cu[cu.year .== y, :], players, cols; styles)
    end
    labels = [top; isempty(inert) ? String[] : ["inert: $(join(inert, ", "))"]]
    entries = [[LineElement(; color=c_, linewidth=3) for c_ in cols[1:length(top)]];
               isempty(inert) ? LineElement[] : [LineElement(; color=OTHER, linestyle=:dot)]]
    Legend(fig[1, length(cy) + 1], entries, labels)
    Label(fig[0, :], "$title: convergence, top inputs"; halign=:left, fontsize=14)
    save(joinpath(fig_dir, "convergence_inputs.png"), fig; px_per_unit=2)
    push!(figs, ("convergence_inputs.png", "Convergence over the runs: share of Var(Y) of the $(MAX_SERIES) " *
        "inputs largest at the full sample, on nested Sobol' prefixes, for selected years (bands: 95 % CI; dotted: " *
        "inert inputs, the noise floor). Estimates still moving between the last two sizes are not converged."))

    ## groups: sum of members (top row) and as one player (bottom row), same colours per group
    for s in unique(c.scheme[c.analysis .== "grouped_player"])
        cp = c[(c.analysis .== "grouped_player") .& (c.scheme .== s), :]
        cs = c[(c.analysis .== "grouped_sum") .& (c.scheme .== s), :]
        groups = top_at_full(cp, MAX_SERIES)
        gcols = CATS[1:length(groups)]
        fig = Figure(; size=(380 * length(cy) + 260, 760))
        for (j, y) in enumerate(cy)
            convergence_panel!(axis_(fig, 1, j, "$y: sum of member effects"), cs[cs.year .== y, :], groups, gcols)
            convergence_panel!(axis_(fig, 2, j, "$y: group as one player"), cp[cp.year .== y, :], groups, gcols)
        end
        Legend(fig[1:2, length(cy) + 1], [LineElement(; color=c_, linewidth=3) for c_ in gcols], groups)
        Label(fig[0, :], "$title: convergence, groups '$s'"; halign=:left, fontsize=14)
        save(joinpath(fig_dir, "convergence_$s.png"), fig; px_per_unit=2)
        push!(figs, ("convergence_$s.png", "Convergence over the runs for the groups of '$s' (the $(length(groups)) " *
            "largest as players at the full sample): top = sum of member effects, bottom = group as one player."))
    end
    return figs
end

function write_index(dir, figs, meta)
    open(joinpath(dir, "SHAPLEY.md"), "w") do io
        println(io, "# $(meta["model_label"]) / $(meta["domain"]): Shapley effects of $(meta["label"])\n")
        println(io, "$(meta["n_rows"]) runs, $(meta["n_inputs"]) inputs, years $(first(meta["years"]))–$(last(meta["years"])). ",
            "Given-data nearest-neighbour estimator (SAShE.jl, pick-and-freeze). Tables: the `*.csv` files ",
            "next to this one (long format; `share` = effect / Var(Y), `share_of_total` = effect / sum of the ",
            "analysis' effects). Confidence bands are SAShE's and anti-conservative for given data.\n")
        println(io, "Chains: the intermediate M is a function of inputs that stay in the analysis (e.g. thermal ",
            "tolerance of model parameters); the estimator allows such dependence, and credit is shared between M ",
            "and those inputs.\n")
        for (f, c) in figs
            println(io, "## $f\n\n$c\n\n![](figures/$f)\n")
        end
    end
end
