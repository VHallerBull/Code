#=
02_dhw_ensemble.jl -- describe and check the shared DHW ensemble of each synthetic domain.

    julia +release --project=Code/sa Code/sa/02_dhw_ensemble.jl

Run after 01. Per domain, writes to domains/<domain>/sa_inputs/:

- dhw_design.csv   one row per DHW-file member: `member`, `control` (member 1 = all-zero DHW), the 7
                   Sobol'-sampled generator parameters (`dhw_*`, see `DHW_FACTORS` in sa_common.jl;
                   empty for the control) and summaries of the resulting trajectories. 06_collect.jl
                   attaches the 7 parameters to every model run from this file.
- figures/dhw_*.png  the QC figures listed in sa_inputs/DHW_QC.md

and checks that the domain's DHW file was generated from this sample (dhw_scenario_params.csv).
=#

using CSV, CairoMakie
using ADRIAanalysis: sensitivity
include(joinpath(@__DIR__, "sa_common.jl"))

const X = dhw_sample()

function dhw_design(dom)::DataFrame
    D = Array(dom.dhw_scens)                    # [timesteps, locations, members]
    size(D, 3) == nrow(X) + 1 || error("DHW file has $(size(D, 3)) members; expected $(nrow(X) + 1).")
    all(iszero, @view D[:, :, 1]) || error("DHW member 1 is not the all-zero control.")
    design = DataFrame(member=1:size(D, 3), control=[s == 1 for s in 1:size(D, 3)])
    for p in DHW_PARAM_NAMES
        design[!, p] = [missing; X[!, p]]
    end
    return hcat(design, select(dhw_member_summaries(D), Not(:member)))
end

"The domain's own record of its DHW parameters must equal the Sobol' sample."
function check_against_domain(domain, design)
    srd = CSV.read(joinpath(domain_path(domain), "dhw_scenario_params.csv"), DataFrame)
    stoch = srd[.!srd.control, :]
    pairs = [:baseline => :dhw_baseline, :trend_per_decade => :dhw_trend_per_decade,
        :interannual_sd => :dhw_interannual_sd, :ar1 => :dhw_ar1, :event_prob => :dhw_event_prob,
        :local_sd => :dhw_local_sd]
    for (a, b) in pairs
        isapprox(stoch[!, a], X[!, b]; rtol=1e-9) || error("$domain: $a in dhw_scenario_params.csv differs from the sample.")
    end
    isapprox(stoch.event_magnitude_max .- stoch.event_magnitude_min, X.dhw_event_mag_width; rtol=1e-9) ||
        error("$domain: heatwave size range differs from the sample.")
    all(stoch.event_magnitude_min .== DHW_EVENT_MAG_MIN) || error("$domain: heatwave minimum is not $DHW_EVENT_MAG_MIN.")
    return nothing
end

function qc_figures(domain, dom, design, fig_dir)
    D = Array(dom.dhw_scens)[:, :, 2:end]
    yrs = collect(ADRIA.timesteps(dom))
    reef_mean = dropdims(mean(D; dims=2); dims=2)       # [timesteps, members]
    q(p) = [quantile(reef_mean[t, :], p) for t in axes(reef_mean, 1)]
    files = String[]

    fig = Figure(; size=(800, 450))
    ax = Axis(fig[1, 1]; xlabel="Year", ylabel="DHW (mean over reefs)",
        title="$domain: $(size(D, 3)) DHW members")
    band!(ax, yrs, q(0.025), q(0.975); color=(:steelblue, 0.2), label="95 %")
    band!(ax, yrs, q(0.25), q(0.75); color=(:steelblue, 0.45), label="50 %")
    lines!(ax, yrs, q(0.5); color=:navy, label="median")
    for s in 1:8
        lines!(ax, yrs, reef_mean[:, s]; color=(:grey30, 0.5), linewidth=0.7)
    end
    axislegend(ax; position=:lt)
    push!(files, "dhw_fan.png"); save(joinpath(fig_dir, files[end]), fig; px_per_unit=2)

    stoch = design[.!design.control, :]
    summaries = [:dhw_mean, :dhw_max, :years_dhw_ge4, :years_dhw_ge8, :dhw_final_year]
    ## a random dummy input gives PAWN's noise floor (ADRIA's PAWN is a scaled KS statistic, not bounded by 1)
    Random.seed!(99)
    names_ = [string.(DHW_PARAM_NAMES); "random dummy"]
    Xm = hcat(Matrix{Float64}(stoch[:, DHW_PARAM_NAMES]), rand(nrow(stoch)))
    fig = Figure(; size=(1100, 420))
    for (i, s) in enumerate(summaries)
        ## pawn's Si columns: min, lb, mean, median, ub, max, std, cv
        med = Array(sensitivity.pawn(Xm, Float64.(stoch[!, s]), names_))[:, 4]
        ax = Axis(fig[1, i]; title=string(s), yticks=(1:length(names_), names_),
            xlabel="PAWN (median)", yticklabelsvisible=i == 1)
        barplot!(ax, 1:length(med), med; direction=:x,
            color=[j == length(med) ? :grey60 : :steelblue for j in eachindex(med)])
        vlines!(ax, [med[end]]; color=:black, linestyle=:dash)
    end
    Label(fig[0, :], "$domain: what drives each DHW summary (PAWN, S = 10; dashed = random dummy, the noise floor)")
    push!(files, "dhw_pawn.png"); save(joinpath(fig_dir, files[end]), fig; px_per_unit=2)

    fig = Figure(; size=(1000, 600))
    for (i, p) in enumerate(DHW_PARAM_NAMES)
        ax = Axis(fig[fld(i - 1, 4) + 1, mod(i - 1, 4) + 1]; xlabel=string(p), ylabel=i in (1, 5) ? "mean DHW" : "")
        scatter!(ax, stoch[!, p], stoch.dhw_mean; markersize=3, color=(:steelblue, 0.5))
    end
    Label(fig[0, :], "$domain: mean DHW (2025-2050, all reefs) per member against each parameter")
    push!(files, "dhw_scatter.png"); save(joinpath(fig_dir, files[end]), fig; px_per_unit=2)
    return files
end

for domain in SA_DOMAINS
    isdir(domain_path(domain)) || error("$(domain_path(domain)) not found: run 01_build_domains.jl first.")
    dom = load_sa_domain(domain)
    design = dhw_design(dom)
    check_against_domain(domain, design)
    out = sa_inputs_dir(domain)
    CSV.write(joinpath(out, "dhw_design.csv"), design)
    fig_dir = mkpath(joinpath(out, "figures"))
    files = qc_figures(domain, dom, design, fig_dir)
    open(joinpath(out, "DHW_QC.md"), "w") do io
        println(io, "# $domain: DHW ensemble\n")
        println(io, "$(nrow(X)) stochastic members (2-$(nrow(X) + 1)) from one Sobol' sample of the 7 generator ",
            "parameters (`dhw_design.csv`), plus member 1 = all-zero DHW (the control). Heatwaves add ",
            "U($(DHW_EVENT_MAG_MIN), $(DHW_EVENT_MAG_MIN) + width) DHW.\n")
        println(io, "| Parameter | Range | Meaning |\n|---|---|---|")
        for (n, lo, hi, what) in DHW_FACTORS
            println(io, "| `$n` | $lo-$hi | $what |")
        end
        println(io, "\n| Figure | Shows |\n|---|---|")
        captions = ["reef-mean DHW per year: median, 50 % and 95 % bands, 8 example members",
            "PAWN of DHW summaries against the 7 parameters and a random dummy (check that each one is above the dummy)",
            "mean DHW per member against each parameter"]
        for (f, c) in zip(files, captions)
            println(io, "| ![](figures/$f) | $c |")
        end
    end
    s = design[.!design.control, :]
    println("$domain: dhw_design.csv written; mean DHW $(round(minimum(s.dhw_mean); digits=1))-",
        "$(round(maximum(s.dhw_mean); digits=1)), years with DHW>=8 per reef $(round(minimum(s.years_dhw_ge8); digits=1))-",
        "$(round(maximum(s.years_dhw_ge8); digits=1))")
end
