#=
05_analyse.jl -- Shapley effects of the pick-and-freeze runs, their noise, and the base samples a full
run needs.

    julia +release --project=Code/sa/shapley Code/sa_pf/05_analyse.jl <model> <band> [band ...]

For each design set (`main`: with the DHW inputs; `zero`: no heat stress, if it was run), band, domain,
output and year: SAShE's estimator on the walks (`analyze(S, Y)`), each input's increment SD sigma_i over
the base samples, and Var(Y) from the A and B runs. With 95 % half-width MOE = 1.96 sigma / sqrt(N), a
target half-width delta (as a share of Var(Y)) needs

    N = (1.96 sigma / (delta Var(Y)))^2   base samples,   runs = N (d + 1).

Writes <ModelFolder>/model_run_pf_pilot/analysis/: effects.csv, group_effects.csv, noise.csv,
required_base_samples.csv (all with a `set` column), figures/*.png (the zero set's with suffix `_zero`)
and PILOT.md.
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
const OUTPUTS = ["relative_cover", "relative_shelter_volume", "coral_evenness", "juvenile_indicator"]
const OUTPUT_LABEL = Dict("relative_cover" => "Relative cover", "relative_shelter_volume" => "Shelter volume",
    "coral_evenness" => "Evenness", "juvenile_indicator" => "Juveniles")
const DELTAS = [0.05, 0.02, 0.01]
const TOP_K = 10
const SCHEME = "process"
"CPU seconds per run on the AIMS HPC (the 8192-run given-data runs), for the cost table."
const HPC_S_PER_RUN = Dict("SingleReef_SobolDHW" => 0.066, "NineReefs_SobolDHW" => 0.236)
"Readable set label for titles."
set_title(set) = set == "main" ? "" : " (zero DHW)"

n_required(σ, vy, δ) = vy > 0 ? ceil(Int, (1.96 * σ / (δ * vy))^2) : 0

"Per-base-sample increments f [d, N] and SAShE's effects / bounds for one output column `y` (walk order)."
function increments(S, y, N)
    Φ, Φ², _ = SAShE.analyze(S, y)
    eff, lb, ub = SAShE.shapley_effects(Φ, Φ²)
    return Φ .* N, eff, lb, ub
end

effects_rows, noise_rows, req_rows, group_rows = NamedTuple[], NamedTuple[], NamedTuple[], NamedTuple[]
running = Dict{Tuple{String,Float64,String,String},Any}()   # (set, band, domain, output) -> cumulative curves, 2050
running_top = Dict{Tuple{String,Float64,String},Any}()      # (set, band, domain) -> running shares, relative cover
seconds_per_run = Dict{Tuple{String,Float64,String},Float64}()
n_inputs = Dict{String,Int}()

for set in PF_SETS, band in BANDS
    dpath = pf_design_path(MODEL, band, set)
    isfile(dpath) || continue
    design = JLD2.load(dpath)
    names_ = design["names"]
    d = length(names_)
    n_inputs[set] = d
    G = load_groupings(joinpath(@__DIR__, "..", "sa", "shapley", "groups", "$MODEL.toml"), names_)[SCHEME]
    for domain in SA_DOMAINS
        path = joinpath(pf_set_dir(MODEL, band, domain, set), "pf_results.jld2")
        isfile(path) || (println("$(band_label(band))$(set_suffix(set))/$domain: no results"); continue)
        R = JLD2.load(path)
        years = R["years"]
        seconds_per_run[(set, band, domain)] = R["seconds"] / R["n_runs"]
        ok = findall(n -> !any(isnan, R["Y"]["relative_cover"][n, :, :]), 1:size(R["Y"]["relative_cover"], 1))
        N = length(ok)
        S = CallablePickAndFreezeSample(design["A"][ok, :], design["B"][ok, :], design["perms"][ok, :])
        println("$(band_label(band))$(set_suffix(set))/$domain: $N walks, $d inputs")
        for out in OUTPUTS
            Yall = R["Y"][out][ok, :, :]                                   # [N, d + 1, T]
            for (t, yr) in enumerate(years)
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
                    push!(effects_rows, (; set, band, domain, output=out, year=yr, input=names_[i],
                        group=group_of(G, i), effect=eff[i], lb=lb[i], ub=ub[i], var_y=vy,
                        share=eff[i] / vy, sigma=σ[i]))
                end
                moe_g = 1.96 .* σg ./ sqrt(N - 1)                      # as SAShE's margin_of_error
                for (g, gname) in enumerate(G.names)
                    push!(group_rows, (; set, band, domain, output=out, year=yr, scheme=SCHEME, group=gname,
                        n_members=length(G.members[g]), effect=effg[g], lb=effg[g] - moe_g[g],
                        ub=effg[g] + moe_g[g], var_y=vy, share=effg[g] / vy, sigma=σg[g]))
                end
                push!(noise_rows, (; set, band, domain, output=out, year=yr, var_y=vy, n_base=N,
                    sigma_max_input=maximum(σ), sigma_median_input=median(σ), sigma_max_top=maximum(σ[top]),
                    sigma_max_group=maximum(σg), sigma_total=σtot, sum_effects=sum(eff),
                    max_share=maximum(eff) / vy, cover_total_share=eff[findfirst(==("m_cover_total"), names_)] / vy))
                if yr == last(years)
                    cs = cumsum(vec(sum(f; dims=1))) ./ (1:N)
                    cv = [k < 2 ? NaN : var(vcat(Ym[1:k, 1], Ym[1:k, end])) for k in 1:N]
                    running[(set, band, domain, out)] = (; cs, cv)
                end
                if out == "relative_cover"
                    rt = get!(running_top, (set, band, domain), Dict{Int,Any}())
                    rt[yr] = (; f, vy, names=names_)
                end
            end
        end
    end
end

effects = DataFrame(effects_rows)
noise = DataFrame(noise_rows)
CSV.write(joinpath(OUT, "effects.csv"), effects)
CSV.write(joinpath(OUT, "noise.csv"), noise)
CSV.write(joinpath(OUT, "group_effects.csv"), DataFrame(group_rows))

## required base samples: worst year (all 26 years must meet the target) and the median year
for g in groupby(noise, [:set, :band, :domain, :output]), δ in DELTAS
    set, band, domain, out = g.set[1], g.band[1], g.domain[1], g.output[1]
    d1 = n_inputs[set] + 1
    for (crit, col) in [("every input", :sigma_max_input), ("top $TOP_K inputs", :sigma_max_top),
        ("groups ($SCHEME)", :sigma_max_group), ("sum of effects", :sigma_total)]
        ns = [n_required(r[col], r.var_y, δ) for r in eachrow(g)]
        nmax, nmed = maximum(ns), round(Int, median(ns))
        worst = g.year[argmax(ns)]
        spr = get(seconds_per_run, (set, band, domain), NaN)
        push!(req_rows, (; set, band, domain, output=out, delta=δ, criterion=crit, n_base_worst_year=nmax,
            worst_year=worst, n_base_median_year=nmed, runs_per_base=d1, runs_worst_year=nmax * d1,
            cpu_h_local=nmax * d1 * spr / 3600, cpu_h_hpc=nmax * d1 * HPC_S_PER_RUN[domain] / 3600))
    end
end
req = DataFrame(req_rows)
CSV.write(joinpath(OUT, "required_base_samples.csv"), req)

## ---------------------------------------------------------------------------------------------------
## Figures (one set of figures per design set; the zero set's file names end in _zero)
## ---------------------------------------------------------------------------------------------------
figs = Tuple{String,String}[]
savefig_(name, fig, cap) = (save(joinpath(FIG, name), fig; px_per_unit=2); push!(figs, (name, cap)))
sets_present = [s for s in PF_SETS if s in noise.set]

for set in sets_present
    sfx, stt = set_suffix(set), set_title(set)
    bands_present = sort(unique(noise.band[noise.set .== set]))

    ## 1. paper-style convergence: sum of effects vs variance, cumulative over base samples (final year)
    for domain in SA_DOMAINS
        any(k -> k[1] == set && k[3] == domain, keys(running)) || continue
        fig = Figure(; size=(330 * length(bands_present) + 40, 230 * length(OUTPUTS) + 60))
        for (i, out) in enumerate(OUTPUTS), (j, band) in enumerate(bands_present)
            haskey(running, (set, band, domain, out)) || continue
            r = running[(set, band, domain, out)]
            n = 1:length(r.cs)
            ax = Axis(fig[i, j]; xscale=log10, xlabel=i == length(OUTPUTS) ? "Base samples" : "",
                ylabel=j == 1 ? OUTPUT_LABEL[out] : "", title=i == 1 ? "±$(round(Int, 100band)) %" : "")
            lines!(ax, n, r.cs; color=CATS[1])
            lines!(ax, n, r.cv; color=CATS[8])
        end
        Legend(fig[length(OUTPUTS) + 1, :], [LineElement(; color=CATS[1], linewidth=3), LineElement(; color=CATS[8], linewidth=3)],
            ["Sum of Shapley effects (cumulative)", "Variance of Y (cumulative)"]; orientation=:horizontal)
        Label(fig[0, :], "$MODEL_LABEL / $domain$stt: pick-and-freeze convergence, year 2050"; halign=:left, fontsize=14)
        savefig_("convergence_sum_$domain$sfx.png", fig, "$domain$stt, year 2050: running estimate of the sum of " *
            "all Shapley effects (mean of (Y_A - Y_B)^2 / 2) and of Var(Y) over the A and B runs, per input band. " *
            "The two estimate the same quantity independently (as in the CoralBlox paper's Fig. S4); this checks " *
            "the total only, not the individual effects.")
    end

    ## 2. running estimates of the top inputs (relative cover), rows = bands, columns = years
    for domain in SA_DOMAINS
        any(k -> k[1] == set && k[3] == domain, keys(running_top)) || continue
        yrs_show = [2025, 2037, 2050]
        fig = Figure(; size=(380 * length(yrs_show) + 250, 280 * length(bands_present) + 60))
        for (i, band) in enumerate(bands_present), (j, yr) in enumerate(yrs_show)
            haskey(running_top, (set, band, domain)) || continue
            r = running_top[(set, band, domain)][yr]
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
        Label(fig[0, :], "$MODEL_LABEL / $domain$stt: running Shapley shares, relative cover (top 5 at the pilot size)";
            halign=:left, fontsize=14)
        savefig_("running_top_$domain$sfx.png", fig, "$domain$stt, relative cover: running estimate of the 5 " *
            "largest shares against the number of base samples, with 95 % bands, per input band and year.")
    end

    ## 3. base samples needed per year, target ±0.02 of Var(Y)
    for domain in SA_DOMAINS
        sub = noise[(noise.set .== set) .& (noise.domain .== domain), :]
        nrow(sub) == 0 && continue
        fig = Figure(; size=(1250, 700))
        for (i, out) in enumerate(OUTPUTS)
            ax = Axis(fig[fld(i - 1, 2) + 1, mod(i - 1, 2) + 1]; yscale=log10, xlabel="Year",
                ylabel="Base samples for ±0.02", title=OUTPUT_LABEL[out])
            for (j, band) in enumerate(bands_present)
                g = sort(sub[(sub.output .== out) .& (sub.band .== band) .& (sub.var_y .> 0), :], :year)
                nrow(g) == 0 && continue
                lines!(ax, g.year, max.(1, n_required.(g.sigma_max_input, g.var_y, 0.02)); color=CATS[j])
                lines!(ax, g.year, max.(1, n_required.(g.sigma_max_group, g.var_y, 0.02)); color=CATS[j], linestyle=:dash)
            end
        end
        entries = [[LineElement(; color=CATS[j], linewidth=3) for j in eachindex(bands_present)];
                   [LineElement(; color=INK2), LineElement(; color=INK2, linestyle=:dash)]]
        labels = [["±$(round(Int, 100b)) %" for b in bands_present]; ["every input", "every group ($SCHEME)"]]
        Legend(fig[1:2, 3], entries, labels)
        Label(fig[0, :], "$MODEL_LABEL / $domain$stt: base samples for a 95 % half-width of 0.02 Var(Y)"; halign=:left, fontsize=14)
        savefig_("required_base_samples_$domain$sfx.png", fig, "$domain$stt: base samples needed per year so that " *
            "every input's (solid) or every group's (dashed) 95 % half-width is 0.02 of Var(Y); runs = base " *
            "samples x $(n_inputs[set] + 1).")
    end

    ## 4. group shares over time at the pilot size (relative cover), one panel per band
    for domain in SA_DOMAINS
        sub = effects[(effects.set .== set) .& (effects.domain .== domain) .& (effects.output .== "relative_cover"), :]
        nrow(sub) == 0 && continue
        gs = combine(groupby(sub, [:band, :year, :group]), :effect => sum => :effect, :var_y => first => :var_y)
        gs.share = gs.effect ./ gs.var_y
        order = combine(groupby(gs, :group), :share => (v -> maximum(abs.(filter(isfinite, v)); init=0.0)) => :pk)
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
        Label(fig[0, :], "$MODEL_LABEL / $domain$stt: group shares of relative cover (pilot, sum of member effects)";
            halign=:left, fontsize=14)
        savefig_("group_shares_$domain$sfx.png", fig, "$domain$stt, relative cover: Shapley share of each " *
            "'$SCHEME' group (sum of member effects) per year and input band, at the pilot size.")
    end
end

## ---------------------------------------------------------------------------------------------------
## Summary
## ---------------------------------------------------------------------------------------------------
open(joinpath(OUT, "PILOT.md"), "w") do io
    println(io, "# $MODEL_LABEL pick-and-freeze pilot (option B)\n")
    println(io, "Sets: ", join(["$s ($(n_inputs[s]) inputs, $(n_inputs[s] + 1) runs per base sample)" for s in sets_present], ", "),
        ". Bands: ", join(["±$(round(Int, 100b)) %" for b in sort(unique(noise.band))], ", "),
        ". Base samples per band and domain: ", join(unique(noise.n_base), ", "),
        ". Tables: effects.csv, group_effects.csv, noise.csv, required_base_samples.csv.\n")
    println(io, "## Base samples needed (relative cover, worst year)\n")
    println(io, "| Set | Band | Domain | Criterion | ±0.05 | ±0.02 | ±0.01 |\n|---|---|---|---|---|---|---|")
    for g in groupby(req[req.output .== "relative_cover", :], [:set, :band, :domain, :criterion])
        v = Dict(r.delta => r.n_base_worst_year for r in eachrow(g))
        println(io, "| $(g.set[1]) | ±$(round(Int, 100g.band[1])) % | $(g.domain[1]) | $(g.criterion[1]) | ",
            join([string(v[δ]) for δ in DELTAS], " | "), " |")
    end
    println(io)
    for (f, c) in figs
        println(io, "## $f\n\n$c\n\n![](figures/$f)\n")
    end
end
println("analysis -> $OUT")
