#=
05_analyse.jl -- Shapley effects of the pick-and-freeze pilot, their noise, and the base samples a full
run needs.

    julia +release --project=Code/sa/shapley Code/sa_pf/05_analyse.jl <model> <band> [band ...]

For each band, domain, output and year: SAShE's estimator on the walks (`analyze(S, Y)`), each input's
increment SD sigma_i over the base samples, and Var(Y) from the A and B runs. With 95 % half-width
MOE = 1.96 sigma / sqrt(N), a target half-width delta (as a share of Var(Y)) needs

    N = (1.96 sigma / (delta Var(Y)))^2   base samples,   runs = N (d + 1).

Writes <ModelFolder>/model_run_pf_pilot/analysis/: effects.csv, group_effects.csv, noise.csv,
required_base_samples.csv,
figures/*.png and PILOT.md.
=#

length(ARGS) >= 2 || error("usage: 05_analyse.jl <model> <band> [band ...]")
const MODEL = ARGS[1]
const BANDS = parse.(Float64, ARGS[2:end])
include(joinpath(@__DIR__, "pf_common.jl"))
include(joinpath(@__DIR__, "..", "sa", "shapley", "groups.jl"))
include(joinpath(@__DIR__, "..", "sa", "shapley", "plots.jl"))   # theme, colours (read only)
using CSV

const OUT = mkpath(joinpath(SA_RUNS_ROOT, MODEL_FOLDERS[MODEL], PF_SUBDIR, "analysis"))
const FIG = mkpath(joinpath(OUT, "figures"))
## The outputs are whatever 04_collect.jl derived (sa_outputs.jl): they differ between models (thermal
## tolerance, juveniles) and domains (spatial_cv only with more than one reef).
const DELTAS = [0.05, 0.02, 0.01]
const TOP_K = 10
const SCHEME = "process"
"CPU seconds per run on the AIMS HPC (the 8192-run given-data runs), for the cost table."
const HPC_S_PER_RUN = Dict("SingleReef_SobolDHW" => 0.066, "NineReefs_SobolDHW" => 0.236)

"Base samples for a 95 % half-width of δ Var(Y); 0 where the output does not vary (nothing to resolve)."
n_required(σ, vy, δ) = vy > 0 && isfinite(σ) ? ceil(Int, (1.96 * σ / (δ * vy))^2) : 0

"Per-base-sample increments f [d, N] and SAShE's effects / bounds for one output column `y` (walk order)."
function increments(S, y, N)
    Φ, Φ², _ = SAShE.analyze(S, y)
    eff, lb, ub = SAShE.shapley_effects(Φ, Φ²)
    return Φ .* N, eff, lb, ub
end

effects_rows, noise_rows, req_rows, group_rows = NamedTuple[], NamedTuple[], NamedTuple[], NamedTuple[]
running = Dict{Tuple{Float64,String,String},Any}()       # (band, domain, output) -> cumulative curves at 2050
running_top = Dict{Tuple{Float64,String},Any}()          # (band, domain) -> running shares, relative cover
seconds_per_run = Dict{Tuple{Float64,String},Float64}()
dropped_rows = NamedTuple[]                              # years left out of an output (missing values)

for band in BANDS
    design = JLD2.load(joinpath(pf_band_dir(MODEL, band), "design_pf.jld2"))
    names_ = design["names"]
    d = length(names_)
    G = load_groupings(joinpath(@__DIR__, "..", "sa", "shapley", "groups", "$MODEL.toml"), names_)[SCHEME]
    for domain in SA_DOMAINS
        path = joinpath(pf_run_dir(MODEL, band, domain), "pf_results.jld2")
        isfile(path) || (println("$(band_label(band))/$domain: no results"); continue)
        R = JLD2.load(path)
        years = R["years"]
        seconds_per_run[(band, domain)] = R["seconds"] / R["n_runs"]
        ok = findall(n -> !any(isnan, R["Y"]["relative_cover"][n, :, :]), 1:size(R["Y"]["relative_cover"], 1))
        N = length(ok)
        S = CallablePickAndFreezeSample(design["A"][ok, :], design["B"][ok, :], design["perms"][ok, :])
        println("$(band_label(band))/$domain: $N walks")
        outs = output_order(keys(R["Y"]); thermal_groups=get(R, "thermal_groups", String[]))
        vy_idx = valid_years(R["Y"], outs, ok)
        for out in outs
            Yall = R["Y"][out][ok, :, :]                                   # [N, d + 1, T]
            for t in setdiff(eachindex(years), vy_idx[out])
                push!(dropped_rows, (; band, domain, output=out, year=years[t]))
            end
            for t in vy_idx[out]
                yr = years[t]
                Ym = Yall[:, :, t]
                y = vec(permutedims(Ym))                                   # walk by walk, SAShE's order
                vy = var(vcat(Ym[:, 1], Ym[:, end]))                       # A and B runs: 2N independent
                f, eff, lb, ub = increments(S, y, N)
                σ = vec(std(f; dims=2))
                fg = reduce(vcat, [sum(f[m, :]; dims=1) for m in G.members])
                σg = vec(std(fg; dims=2))
                effg = vec(mean(fg; dims=2))
                σtot = std(vec(sum(f; dims=1)))
                top = partialsortperm(abs.(eff), 1:TOP_K; rev=true)
                for i in 1:d
                    push!(effects_rows, (; band, domain, output=out, year=yr, input=names_[i],
                        group=group_of(G, i), effect=eff[i], lb=lb[i], ub=ub[i], var_y=vy,
                        share=eff[i] / vy, sigma=σ[i]))
                end
                moe_g = 1.96 .* σg ./ sqrt(N - 1)                      # as SAShE's margin_of_error
                for (g, gname) in enumerate(G.names)
                    push!(group_rows, (; band, domain, output=out, year=yr, scheme=SCHEME, group=gname,
                        n_members=length(G.members[g]), effect=effg[g], lb=effg[g] - moe_g[g],
                        ub=effg[g] + moe_g[g], var_y=vy, share=effg[g] / vy, sigma=σg[g]))
                end
                push!(noise_rows, (; band, domain, output=out, year=yr, var_y=vy, n_base=N,
                    sigma_max_input=maximum(σ), sigma_median_input=median(σ), sigma_max_top=maximum(σ[top]),
                    sigma_max_group=maximum(σg), sigma_total=σtot, sum_effects=sum(eff),
                    max_share=maximum(eff) / vy, cover_total_share=eff[findfirst(==("m_cover_total"), names_)] / vy))
                if yr == last(years)
                    cs = cumsum(vec(sum(f; dims=1))) ./ (1:N)
                    cv = [k < 2 ? NaN : var(vcat(Ym[1:k, 1], Ym[1:k, end])) for k in 1:N]
                    running[(band, domain, out)] = (; cs, cv)
                end
                if out == "relative_cover"
                    rt = get!(running_top, (band, domain), Dict{Int,Any}())
                    rt[yr] = (; f, vy, names=names_)
                end
            end
        end
    end
end

effects = DataFrame(effects_rows)
noise = DataFrame(noise_rows)
dropped = DataFrame(dropped_rows)
isempty(dropped_rows) || println("years left out (missing values): ", nrow(dropped), " output-years; see dropped_years.csv")
CSV.write(joinpath(OUT, "dropped_years.csv"), isempty(dropped_rows) ?
    DataFrame(band=Float64[], domain=String[], output=String[], year=Int[]) : dropped)
const OUTPUTS = output_order(unique(noise.output))
CSV.write(joinpath(OUT, "effects.csv"), effects)
CSV.write(joinpath(OUT, "noise.csv"), noise)
CSV.write(joinpath(OUT, "group_effects.csv"), DataFrame(group_rows))

## required base samples: worst year (all 26 years must meet the target) and the median year
for g in groupby(noise, [:band, :domain, :output]), δ in DELTAS
    band, domain, out = g.band[1], g.domain[1], g.output[1]
    d1 = length(unique(effects.input)) + 1
    for (crit, col) in [("every input", :sigma_max_input), ("top $TOP_K inputs", :sigma_max_top),
        ("groups ($SCHEME)", :sigma_max_group), ("sum of effects", :sigma_total)]
        ns = [n_required(r[col], r.var_y, δ) for r in eachrow(g)]
        nmax, nmed = maximum(ns), round(Int, median(ns))
        worst = g.year[argmax(ns)]
        spr = get(seconds_per_run, (band, domain), NaN)
        push!(req_rows, (; band, domain, output=out, delta=δ, criterion=crit, n_base_worst_year=nmax,
            worst_year=worst, n_base_median_year=nmed, runs_worst_year=nmax * d1,
            cpu_h_local=nmax * d1 * spr / 3600, cpu_h_hpc=nmax * d1 * HPC_S_PER_RUN[domain] / 3600))
    end
end
req = DataFrame(req_rows)
CSV.write(joinpath(OUT, "required_base_samples.csv"), req)

## ---------------------------------------------------------------------------------------------------
## Figures
## ---------------------------------------------------------------------------------------------------
figs = Tuple{String,String}[]
savefig_(name, fig, cap) = (save(joinpath(FIG, name), fig; px_per_unit=2); push!(figs, (name, cap)))
bands_present = sort(unique(noise.band))

## 1. paper-style convergence: sum of effects vs variance, cumulative over base samples (final year)
for domain in SA_DOMAINS
    any(k -> k[2] == domain, keys(running)) || continue
    douts = [o for o in OUTPUTS if any(k -> k[2] == domain && k[3] == o, keys(running))]
    fig = Figure(; size=(330 * length(bands_present) + 40, 230 * length(douts) + 60))
    for (i, out) in enumerate(douts), (j, band) in enumerate(bands_present)
        haskey(running, (band, domain, out)) || continue
        r = running[(band, domain, out)]
        n = 1:length(r.cs)
        ax = Axis(fig[i, j]; xscale=log10, xlabel=i == length(douts) ? "Base samples" : "",
            ylabel=j == 1 ? output_label(out) : "", title=i == 1 ? "±$(round(Int, 100band)) %" : "")
        lines!(ax, n, r.cs; color=CATS[1], label="Sum of Shapley effects")
        lines!(ax, n, r.cv; color=CATS[8], label="Variance of Y")
    end
    Legend(fig[length(douts) + 1, :], [LineElement(; color=CATS[1], linewidth=3), LineElement(; color=CATS[8], linewidth=3)],
        ["Sum of Shapley effects (cumulative)", "Variance of Y (cumulative)"]; orientation=:horizontal)
    Label(fig[0, :], "$MODEL_LABEL / $domain: pick-and-freeze convergence, year 2050"; halign=:left, fontsize=14)
    savefig_("convergence_sum_$domain.png", fig, "$domain, year 2050: running estimate of the sum of all " *
        "Shapley effects (mean of (Y_A - Y_B)^2 / 2) and of Var(Y) over the A and B runs, per input band. " *
        "The two estimate the same quantity independently (as in the CoralBlox paper's Fig. S4); this checks " *
        "the total only, not the individual effects.")
end

## 2. running estimates of the top inputs (relative cover), rows = bands, columns = years
for domain in SA_DOMAINS
    any(k -> k[2] == domain, keys(running_top)) || continue
    yrs_show = [2025, 2037, 2050]
    fig = Figure(; size=(380 * length(yrs_show) + 250, 280 * length(bands_present) + 60))
    for (i, band) in enumerate(bands_present), (j, yr) in enumerate(yrs_show)
        haskey(running_top, (band, domain)) || continue
        r = running_top[(band, domain)][yr]
        N = size(r.f, 2)
        top = partialsortperm(abs.(vec(mean(r.f; dims=2))), 1:5; rev=true)
        ax = Axis(fig[i, j]; xscale=log10, title="±$(round(Int, 100band)) %, $yr", titlesize=12,
            xlabel=i == length(bands_present) ? "Base samples" : "", ylabel=j == 1 ? "Share of Var(Y)" : "")
        n = 10:N
        for (c, k) in enumerate(top)
            m = cumsum(r.f[k, :]) ./ (1:N)
            s = [std(r.f[k, 1:q]) for q in n]
            band!(ax, n, (m[n] .- 1.96 .* s ./ sqrt.(n)) ./ r.vy, (m[n] .+ 1.96 .* s ./ sqrt.(n)) ./ r.vy; color=(CATS[c], 0.15))
            lines!(ax, n, m[n] ./ r.vy; color=CATS[c], label=r.names[k])
        end
        axislegend(ax; position=:rt, labelsize=9, backgroundcolor=(:white, 0.7))
    end
    Label(fig[0, :], "$MODEL_LABEL / $domain: running Shapley shares, relative cover (top 5 at the pilot size)";
        halign=:left, fontsize=14)
    savefig_("running_top_$domain.png", fig, "$domain, relative cover: running estimate of the 5 largest " *
        "shares against the number of base samples, with 95 % bands, per input band and year.")
end

## 3. base samples needed per year (relative cover and the other outputs), target ±0.02 of Var(Y)
for domain in SA_DOMAINS
    sub = noise[noise.domain .== domain, :]
    nrow(sub) == 0 && continue
    douts = [o for o in OUTPUTS if o in sub.output]
    nrows_ = cld(length(douts), 2)
    fig = Figure(; size=(1250, 350 * nrows_))
    for (i, out) in enumerate(douts)
        ax = Axis(fig[fld(i - 1, 2) + 1, mod(i - 1, 2) + 1]; yscale=log10, xlabel="Year",
            ylabel="Base samples for ±0.02", title=output_label(out))
        for (j, band) in enumerate(bands_present)
            g = sort(sub[(sub.output .== out) .& (sub.band .== band), :], :year)
            for (col, ls) in ((:sigma_max_input, :solid), (:sigma_max_group, :dash))
                n = n_required.(g[!, col], g.var_y, 0.02)
                keep = n .> 0                                      # log axis: years without variance left out
                any(keep) && lines!(ax, g.year[keep], n[keep]; color=CATS[j], linestyle=ls)
            end
        end
    end
    entries = [[LineElement(; color=CATS[j], linewidth=3) for j in eachindex(bands_present)];
               [LineElement(; color=INK2), LineElement(; color=INK2, linestyle=:dash)]]
    labels = [["±$(round(Int, 100b)) %" for b in bands_present]; ["every input", "every group ($SCHEME)"]]
    Legend(fig[1:nrows_, 3], entries, labels)
    Label(fig[0, :], "$MODEL_LABEL / $domain: base samples for a 95 % half-width of 0.02 Var(Y)"; halign=:left, fontsize=14)
    savefig_("required_base_samples_$domain.png", fig, "$domain: base samples needed per year so that every " *
        "input's (solid) or every group's (dashed) 95 % half-width is 0.02 of Var(Y); runs = base samples x 177.")
end

## 4. group shares over time at the pilot size (relative cover), one panel per band
for domain in SA_DOMAINS
    sub = effects[(effects.domain .== domain) .& (effects.output .== "relative_cover"), :]
    nrow(sub) == 0 && continue
    gs = combine(groupby(sub, [:band, :year, :group]), :effect => sum => :effect, :var_y => first => :var_y)
    gs.share = gs.effect ./ gs.var_y
    order = combine(groupby(gs, :group), :share => (v -> maximum(abs.(v))) => :pk)
    groups = first(sort(order, :pk; rev=true), MAX_SERIES).group
    fig = Figure(; size=(380 * length(bands_present) + 250, 400))
    for (j, band) in enumerate(bands_present)
        ax = Axis(fig[1, j]; xlabel="Year", ylabel=j == 1 ? "Share of Var(Y)" : "", title="±$(round(Int, 100band)) %")
        for (c, gname) in enumerate(groups)
            r = sort(gs[(gs.band .== band) .& (gs.group .== gname), :], :year)
            lines!(ax, r.year, r.share; color=CATS[c])
        end
    end
    Legend(fig[1, length(bands_present) + 1], [LineElement(; color=CATS[c], linewidth=3) for c in eachindex(groups)], groups)
    Label(fig[0, :], "$MODEL_LABEL / $domain: group shares of relative cover (pilot, sum of member effects)";
        halign=:left, fontsize=14)
    savefig_("group_shares_$domain.png", fig, "$domain, relative cover: Shapley share of each '$SCHEME' group " *
        "(sum of member effects) per year and input band, at the pilot size.")
end

## ---------------------------------------------------------------------------------------------------
## Summary
## ---------------------------------------------------------------------------------------------------
open(joinpath(OUT, "PILOT.md"), "w") do io
    println(io, "# $MODEL_LABEL pick-and-freeze pilot (option B)\n")
    println(io, "Bands: ", join(["±$(round(Int, 100b)) %" for b in bands_present], ", "),
        ". Base samples per band and domain: ", join(unique(noise.n_base), ", "), "; runs per base sample: ",
        length(unique(effects.input)) + 1, ". Tables: effects.csv, noise.csv, required_base_samples.csv.\n")
    println(io, "## Base samples needed (relative cover, worst year)\n")
    println(io, "| Band | Domain | Criterion | ±0.05 | ±0.02 | ±0.01 |\n|---|---|---|---|---|---|")
    for g in groupby(req[req.output .== "relative_cover", :], [:band, :domain, :criterion])
        v = Dict(r.delta => r.n_base_worst_year for r in eachrow(g))
        println(io, "| ±$(round(Int, 100g.band[1])) % | $(g.domain[1]) | $(g.criterion[1]) | ",
            join([string(v[δ]) for δ in DELTAS], " | "), " |")
    end
    println(io)
    for (f, c) in figs
        println(io, "## $f\n\n$c\n\n![](figures/$f)\n")
    end
end
println("analysis -> $OUT")
