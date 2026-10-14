#=
09_ANALYSE_SIX_ENGINES.jl -- the full six-engine comparison, from what `08_RUN_ENGINE.jl` saved.

    julia +1.12 --project=. 09_ANALYSE_SIX_ENGINES.jl

Loads every `ModelRuns/<project>/<engine>/` folder that has a result (missing engines are reported
and skipped), checks that all of them used the same domain, climate member and years, and writes:

    images/<section>/*.png   the figures, one folder per section (A-K below)
    summary.csv              per engine x metric x year: median, IQR, min, max over scenarios
    benchmark_summary.csv    per engine: time, memory, disk (from benchmark.csv + Slurm accounting)
    FIGURES.md               every figure with a one-line caption, plus the caveats behind them

Sections
    A  sanity        initial state per site vs the domain, the forcing every engine received
    B  domain        domain-mean metrics (cover, juveniles, shelter, evenness, diversity, indices)
    C  engines       each engine's scenario spread; final-year distributions
    D  taxa          native functional groups; harmonised onto ADRIA's 5 groups; composition
    E  size          size-distribution shape over time; mean colony diameter
    F  disturbance   response to the DHW trajectory: loss vs DHW, recovery, years above thresholds
    G  spatial       per-site trajectories, maps, between-engine spread, site-ranking agreement
    H  agreement     engine x engine correlation; within- vs between-engine variance
    I  sensitivity   PAWN / TSA / outcome maps per engine (illustrative at small N)
    J  corset        CORSET's algae and consumer state
    K  benchmarks    load / first-result / per-scenario time, CPU time, peak memory, disk

Each figure is made independently: one that fails is reported in FIGURES.md and the console, and the
rest are still written. Engine identity is colour + line style + marker throughout (the palette's
closest pair under deuteranopia is CoralPop vs Kora, Delta E 7.6, so colour alone is not enough).

Configuration: the same `CMR_<NAME>` variables as 08 (`engine_config` in engines.jl);
`CMR_PROJECT_NAME` picks the project to analyse. `CMR_ONLY_FIGURES` (comma-separated figure file names,
e.g. `composition_initial`) builds only those figures and merges their lines into the existing FIGURES.md. The domain, RCP and climate member are read from the
saved results instead. `CMR_DOMAIN_PATH`, if set, only says where that domain's folder is; otherwise the
path the engines recorded is used, or the default when that doesn't exist here (results copied from
the HPC).
=#

include(joinpath(@__DIR__, "engines.jl"))
using ADRIAanalysis
## CairoMakie: a pure software rasteriser for a headless script writing PNGs (see 06 for why not
## WGLMakie). `import`, so Makie's exports don't clash with ADRIA's.
import ADRIAviz, CairoMakie, GeoMakie, GraphMakie
CairoMakie.activate!()
const Mk = CairoMakie
const M = ADRIA.metrics

cfg = engine_config()
const PROJECT_DIR = model_run_dir(cfg.project)
const IMG = mkpath(joinpath(PROJECT_DIR, "images"))
@info "Analysing" PROJECT_DIR

# =============================================================================
# 0. Load every engine that has a result
# =============================================================================

## The domain, RCP and climate member come from the saved results, not from this shell's CMR_* variables,
## so the forcing figures and FIGURES.md always describe what the engines actually ran. CMR_DOMAIN_PATH
## only says where the domain folder is (needed after copying results to another machine); its folder
## name must still match the one the engines ran on.
function run_settings(project_dir)
    for key in ENGINE_ORDER
        d = joinpath(project_dir, String(key))
        has_engine_result(key, d) && return JLD2.load(joinpath(d, "samples.jld2"), "settings")
    end
    error("No engine results under $project_dir -- run 08_RUN_ENGINE.jl first.")
end
ran = run_settings(PROJECT_DIR)
ran_path = get(ran, "DOMAIN_PATH", "")
ran_domain = (haskey(ENV, "CMR_DOMAIN_PATH") || !ispath(ran_path)) ? cfg.domain_path : ran_path
ran_name = get(ran, "DOMAIN_NAME", basename(normpath(ran_domain)))
basename(normpath(ran_domain)) == ran_name || error("The engines ran on the domain `$ran_name`, but " *
    "$ran_domain is a different folder. Set CMR_DOMAIN_PATH to that domain's folder.")
ran_cfg = (domain_path=ran_domain, rcp=get(ran, "RCP", cfg.rcp), dhw_idx=get(ran, "DHW_IDX", cfg.dhw_idx),
    wave_idx=get(ran, "WAVE_IDX", cfg.wave_idx), cyclone_idx=get(ran, "CYCLONE_IDX", cfg.cyclone_idx))
differs = [k for k in (:rcp, :dhw_idx, :wave_idx, :cyclone_idx) if ran_cfg[k] != cfg[k]]
isempty(differs) || @warn "This shell's CMR_* settings differ from the saved results ($(join(differs, ", "))); using the results' values." ran_cfg
cfg = merge(cfg, ran_cfg)

dom, _, YEARS = load_shared_domain(cfg)
const NY = length(YEARS)

struct Engine
    key::Symbol
    label::String
    color::Any
    linestyle::Symbol
    marker::Symbol
    rs::Any
    factors::DataFrame
    settings::Dict{String,Any}
    dir::String
end

const STYLE = Dict(:coralblox => (:solid, :circle), :cscape => (:dash, :rect), :coralpop => (:dot, :utriangle),
    :kora => (:dashdot, :diamond), :corset => (:dashdotdot, :dtriangle), :reefmod => (:solid, :star5))

engines = Engine[]
missing_engines = String[]
for key in ENGINE_ORDER
    d = joinpath(PROJECT_DIR, String(key))
    if !has_engine_result(key, d)
        push!(missing_engines, ENGINE_INFO[key].label)
        continue
    end
    s = JLD2.load(joinpath(d, "samples.jld2"))
    rs = load_engine_result(key, d, YEARS)
    varying = Symbol.(s["varying"])
    X = DataFrames.select(s["scens"], intersect(varying, Symbol.(names(s["scens"]))))
    push!(engines, Engine(key, ENGINE_INFO[key].label, Mk.Makie.to_color(ENGINE_INFO[key].color),
        STYLE[key]..., rs, X, s["settings"], d))
end
isempty(engines) && error("No engine results under $PROJECT_DIR -- run 08_RUN_ENGINE.jl first.")
isempty(missing_engines) || @warn "Not in this comparison (no saved result): $(join(missing_engines, ", "))"

## --- fairness: same domain, member, years for everyone; same locations and timesteps ---
ref = engines[1]
for e in engines[2:end]
    bad = [k for k in FAIRNESS_KEYS if get(e.settings, k, missing) != get(ref.settings, k, missing)]
    isempty(bad) || error("$(e.label) and $(ref.label) were run with different settings: " *
        join(["$k ($(get(ref.settings, k, missing)) vs $(get(e.settings, k, missing)))" for k in bad], "; "))
    collect(e.rs.loc_ids) == collect(ref.rs.loc_ids) || error("$(e.label) and $(ref.label) cover different locations.")
    collect(ADRIA.timesteps(e.rs)) == YEARS || error("$(e.label) covers different years than the domain.")
end
n_scen = Dict(e.key => ADRIA.n_scenarios(e.rs) for e in engines)
@info "Engines in the comparison" labels = [e.label for e in engines] n_scen

const LOC_IDS = collect(ref.rs.loc_ids)
const NL = length(LOC_IDS)
const KW = ADRIA.loc_k_area(ref.rs) ./ sum(ADRIA.loc_k_area(ref.rs))   # k-area weights, shared by all engines
const MEMBER = cfg.dhw_idx == 0 ? "no heat stress (DHW member 0), wave member $(cfg.wave_idx)" :
    "DHW member $(cfg.dhw_idx), wave member $(cfg.wave_idx)"
const NSTR = join(["$(e.label) $(n_scen[e.key])" for e in engines], ", ")

# =============================================================================
# Helpers: figure registry, metric extraction
# =============================================================================

figure_index = NamedTuple[]
failed_figures = String[]

## CMR_ONLY_FIGURES="composition_initial,taxa_harmonised": build just those figures (by file name, no .png)
## and leave summary.csv / benchmark_summary.csv alone; their lines are merged into the existing FIGURES.md.
const ONLY_FIGURES = filter(!isempty, strip.(split(get(ENV, "CMR_ONLY_FIGURES", ""), ',')))

"Build and save one figure; failures are recorded, not fatal."
function figure(f::Function, section::AbstractString, name::AbstractString, caption::AbstractString)
    isempty(ONLY_FIGURES) || name in ONLY_FIGURES || return nothing
    try
        fig = f()
        fig === nothing && return nothing
        dir = mkpath(joinpath(IMG, section))
        Mk.save(joinpath(dir, name * ".png"), fig; px_per_unit=2)
        push!(figure_index, (section=section, file="images/$section/$name.png", caption=caption))
    catch err
        push!(failed_figures, "$section/$name: $(sprint(showerror, err))")
        @warn "Figure $section/$name failed" exception = (err, catch_backtrace())
    end
    return nothing
end

"`[timesteps, locations, scenarios]` -> k-area weighted domain mean `[timesteps, scenarios]`."
domain_mean(x::AbstractArray{<:Real,3}) = dropdims(sum(Float64.(x) .* reshape(KW, 1, :, 1); dims=2); dims=2)

"Engine `e`'s metric `f` (named `m`) as a plain `[timesteps, locations, scenarios]` array, or `missing`."
function loc_metric(e::Engine, f, m)
    try
        a = Array(f(e.rs))
        ndims(a) == 3 && size(a, 1) == NY && size(a, 2) == NL || return missing
        return Float64.(a)
    catch err
        ## an engine without this outcome (KeyError/MethodError) is expected; anything else is reported
        err isa Union{KeyError,MethodError} || @warn "$(e.label): could not compute $m" exception = err
        return missing
    end
end

const METRICS = [
    (:relative_cover, "Relative cover", M.relative_cover),
    (:relative_juveniles, "Relative juveniles", M.relative_juveniles),
    (:juvenile_indicator, "Juvenile indicator", M.juvenile_indicator),
    (:relative_shelter_volume, "Relative shelter volume", M.relative_shelter_volume),
    (:absolute_shelter_volume, "Absolute shelter volume", M.absolute_shelter_volume),
    (:coral_evenness, "Coral evenness", M.coral_evenness),
    (:coral_diversity, "Coral diversity", M.coral_diversity),
    (:reef_condition_index, "Reef Condition Index", M.reef_condition_index),
    (:reef_tourism_index, "Reef Tourism Index", M.reef_tourism_index),
    (:reef_fish_index, "Reef Fish Index", M.reef_fish_index),
    (:reef_biodiversity_condition_index, "Reef Biodiversity Condition Index", M.reef_biodiversity_condition_index),
]
const METRIC_TITLE = Dict(m => t for (m, t, _) in METRICS)

## The taxa-based figures don't use the domain metrics (slow to compute); a filtered run of only those skips them.
const TAXA_ONLY_FIGURES = ("taxa_native", "taxa_harmonised", "composition_initial", "composition_final")
const NEED_METRICS = isempty(ONLY_FIGURES) || any(!in(TAXA_ONLY_FIGURES), ONLY_FIGURES)
loc = Dict(e.key => Dict(m => (NEED_METRICS ? loc_metric(e, f, m) : missing) for (m, _, f) in METRICS) for e in engines)
dommean = Dict(e.key => Dict(m => (a === missing ? missing : domain_mean(a)) for (m, a) in loc[e.key]) for e in engines)

## Relative shelter volume (and the indices built on it) is defined on [0, 1]; an engine outside that
## used a different normalisation and is left out of those panels -- 07's rule, unchanged.
const SHELTER_DEPENDENT = (:relative_shelter_volume, :reef_condition_index, :reef_tourism_index,
    :reef_fish_index, :reef_biodiversity_condition_index)
invalid_shelter = Set{Symbol}()
for e in engines
    rsv = dommean[e.key][:relative_shelter_volume]
    if rsv === missing || maximum(rsv) > 1 + 1e-9
        push!(invalid_shelter, e.key)
    end
end
plottable(e, m) = dommean[e.key][m] !== missing && !(m in SHELTER_DEPENDENT && e.key in invalid_shelter)

q(x, p) = quantile(x, p)
"Median and IQR over scenarios of `[timesteps, scenarios]`."
med_iqr(a) = (vec(median(a; dims=2)), [q(r, 0.25) for r in eachrow(a)], [q(r, 0.75) for r in eachrow(a)])

"One engine's median line + IQR band + thin per-scenario lines on `ax`."
function spread!(ax, e::Engine, a; scen_lines=true, label=e.label)
    m, lo, hi = med_iqr(a)
    scen_lines && for s in axes(a, 2)
        Mk.lines!(ax, YEARS, a[:, s]; color=(e.color, 0.18), linewidth=0.6)
    end
    Mk.band!(ax, YEARS, lo, hi; color=(e.color, 0.22))
    Mk.lines!(ax, YEARS, m; color=e.color, linewidth=2, linestyle=e.linestyle, label=label)
end

"A legend entry per engine: colour + line style + marker."
function engine_legend!(pos, es=engines; kw...)
    elems = [[Mk.LineElement(color=e.color, linestyle=e.linestyle, linewidth=2),
        Mk.MarkerElement(color=e.color, marker=e.marker, markersize=9)] for e in es]
    Mk.Legend(pos, elems, [e.label for e in es]; framevisible=false, kw...)
end

## One climate member as [years, locs]; member 0 is ADRIA's "none" (all zero).
member_slice(scens, idx) = idx == 0 ? zeros(size(scens, 1), size(scens, 2)) : Array(scens[:, :, idx])
loc_dhw = member_slice(dom.dhw_scens, cfg.dhw_idx)                                # [years, locs]
domain_dhw = vec(sum(loc_dhw .* KW'; dims=2))                                     # [years]

# =============================================================================
# A. Sanity: initial state and forcing
# =============================================================================

figure("A_sanity", "initial_cover_by_site",
    "Relative cover in the first year at each site, per engine (median over scenarios); black bars = the domain's own initial cover. Engines that translate the domain faithfully sit on the bar.") do
    fig = Mk.Figure(; size=(1100, 450))
    ax = Mk.Axis(fig[1, 1]; xticks=(1:NL, LOC_IDS), xticklabelrotation=π / 4, ylabel="Relative cover (fraction of k)",
        title="Initial relative cover by site")
    init = vec(sum(Array(dom.init_coral_cover); dims=1))                        # fraction of k per site
    Mk.scatter!(ax, 1:NL, init; marker=:hline, markersize=28, color=:black, label="Domain")
    off = range(-0.3, 0.3; length=length(engines))
    for (i, e) in enumerate(engines)
        a = loc[e.key][:relative_cover]
        a === missing && continue
        Mk.scatter!(ax, (1:NL) .+ off[i], vec(median(a[1, :, :]; dims=2)); color=e.color, marker=e.marker, markersize=10)
    end
    engine_legend!(fig[1, 2])
    fig
end

figure("A_sanity", "forcing", "The climate forcing all engines share ($MEMBER): DHW per site (grey) and k-area-weighted domain mean (black); wave orbital velocity likewise. Cyclone mortality in the Test domain is all zero.") do
    fig = Mk.Figure(; size=(1000, 600))
    ax1 = Mk.Axis(fig[1, 1]; ylabel="DHW (°C-weeks)", title="Degree heating weeks, $MEMBER")
    for l in 1:NL
        Mk.lines!(ax1, YEARS, loc_dhw[:, l]; color=(:grey50, 0.5), linewidth=0.8)
    end
    Mk.lines!(ax1, YEARS, domain_dhw; color=:black, linewidth=2)
    Mk.hlines!(ax1, [4, 8]; color=:grey30, linestyle=:dash, linewidth=1)
    ax2 = Mk.Axis(fig[2, 1]; ylabel="Ub (m/s)", xlabel="Year", title="Wave orbital velocity")
    w = member_slice(dom.wave_scens, cfg.wave_idx)
    for l in 1:NL
        Mk.lines!(ax2, YEARS, w[:, l]; color=(:grey50, 0.5), linewidth=0.8)
    end
    Mk.lines!(ax2, YEARS, vec(sum(w .* KW'; dims=2)); color=:black, linewidth=2)
    Mk.linkxaxes!(ax1, ax2)
    fig
end

# =============================================================================
# B. Domain metrics
# =============================================================================

avail_metrics = [(m, t) for (m, t, _) in METRICS if any(plottable(e, m) for e in engines)]

figure("B_domain", "domain_metrics", "Domain-mean (k-area weighted) metrics: median (line), interquartile range (band) and every scenario (thin lines) per engine. Engines without a metric, or whose shelter volume is outside [0, 1], are omitted from that panel (named in the panel title).") do
    nc = 3
    nr = cld(length(avail_metrics), nc)
    fig = Mk.Figure(; size=(1400, 380 * nr))
    for (i, (m, t)) in enumerate(avail_metrics)
        omitted = [e.label for e in engines if !plottable(e, m)]
        ax = Mk.Axis(fig[fld1(i, nc), mod1(i, nc)]; title=isempty(omitted) ? t : "$t\n(not shown: $(join(omitted, ", ")))",
            titlesize=isempty(omitted) ? 15 : 12, xlabel=fld1(i, nc) == nr ? "Year" : "")
        for e in engines
            plottable(e, m) && spread!(ax, e, dommean[e.key][m])
        end
    end
    engine_legend!(fig[nr + 1, :]; orientation=:horizontal)
    Mk.Label(fig[0, :], "Domain-mean metrics ($NSTR scenarios; $MEMBER)"; fontsize=18)
    fig
end

figure("B_domain", "cover_indexed", "Relative cover indexed to each engine's own first year (1 = initial cover): compares the shape of the trajectories independently of differences in starting cover.") do
    fig = Mk.Figure(; size=(1000, 500))
    ax = Mk.Axis(fig[1, 1]; xlabel="Year", ylabel="Cover / initial cover", title="Relative cover, indexed to the first year")
    for e in engines
        a = dommean[e.key][:relative_cover]
        a === missing && continue
        spread!(ax, e, a ./ max.(a[1:1, :], 1e-12))
    end
    Mk.hlines!(ax, [1.0]; color=:grey40, linestyle=:dash)
    engine_legend!(fig[1, 2])
    fig
end

# =============================================================================
# C. Per engine
# =============================================================================

figure("C_engines", "cover_small_multiples", "Domain relative cover per engine, every scenario (thin) with median (bold), on a shared y-axis.") do
    n = length(engines)
    fig = Mk.Figure(; size=(1300, 330 * cld(n, 3)))
    axs = []
    for (i, e) in enumerate(engines)
        a = dommean[e.key][:relative_cover]
        ax = Mk.Axis(fig[fld1(i, 3), mod1(i, 3)]; title="$(e.label) ($(n_scen[e.key]) scenarios)")
        push!(axs, ax)
        a === missing && continue
        spread!(ax, e, a)
    end
    Mk.linkyaxes!(axs...)
    fig
end

for (yr_label, t) in (("final", NY), ("2050", findfirst(==(2050), YEARS)))
    t === nothing && continue
    figure("C_engines", "cover_distribution_$yr_label", "Distribution over scenarios of domain relative cover in $(YEARS[t]), per engine (box = IQR, whiskers = 1.5 IQR, dots = scenarios).") do
        fig = Mk.Figure(; size=(900, 500))
        ax = Mk.Axis(fig[1, 1]; xticks=(1:length(engines), [e.label for e in engines]), ylabel="Relative cover",
            title="Domain relative cover in $(YEARS[t])")
        for (i, e) in enumerate(engines)
            a = dommean[e.key][:relative_cover]
            a === missing && continue
            v = a[t, :]
            Mk.boxplot!(ax, fill(i, length(v)), v; color=(e.color, 0.35), strokecolor=e.color, strokewidth=1.5, width=0.5)
            Mk.scatter!(ax, i .+ 0.12 .* (rand(length(v)) .- 0.5), v; color=e.color, marker=e.marker, markersize=8)
        end
        fig
    end
end

# =============================================================================
# D. Taxa
# =============================================================================

"Engine `e`'s relative taxa cover as `[timesteps, groups, scenarios]` (domain level), or `missing`."
function taxa_cover(e::Engine)
    try
        a = if e.key === :cscape
            Array(M._relative_taxa_cover(e.rs))
        elseif e.key === :corset
            ## ADRIA's Metric wrapper expects [t, groups, scenarios]; CORSET stores [t, groups, locations, scenarios]
            Array(e.rs.outcomes[:relative_taxa_cover])
        else
            Array(M.relative_taxa_cover(e.rs))
        end
        if ndims(a) == 4                                                   # [t, groups, locations, scenarios]
            a = dropdims(sum(a .* reshape(KW, 1, 1, :, 1); dims=3); dims=3)
        end
        size(a, 1) == NY + 1 && (a = a[2:end, :, :])                       # an extra initial year
        return Float64.(a)
    catch err
        @warn "No taxa cover for $(e.label)" exception = err
        return missing
    end
end
taxa = Dict(e.key => taxa_cover(e) for e in engines)
const GROUP_COLORS = Mk.Makie.to_color.(["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300"])

figure("D_taxa", "taxa_native", "Each engine's cover by its OWN functional groups (median over scenarios, IQR band). Group sets differ between engines, so colours mean different groups in different panels; see each panel's legend.") do
    es = [e for e in engines if taxa[e.key] !== missing]
    fig = Mk.Figure(; size=(1400, 360 * cld(length(es), 3)))
    for (i, e) in enumerate(es)
        a = taxa[e.key]
        r, c = fld1(i, 3), mod1(i, 3)
        ax = Mk.Axis(fig[r, c][1, 1]; title=e.label, xlabel="Year", ylabel="Relative cover")
        gl = ENGINE_INFO[e.key].groups
        for g in axes(a, 2)
            m, lo, hi = med_iqr(a[:, g, :])
            col = GROUP_COLORS[mod1(g, length(GROUP_COLORS))]
            Mk.band!(ax, YEARS, lo, hi; color=(col, 0.2))
            Mk.lines!(ax, YEARS, m; color=col, linewidth=2, label=g <= length(gl) ? replace(gl[g], "_" => " ") : "group $g")
        end
        Mk.axislegend(ax; position=:rt, labelsize=9, framevisible=false)
    end
    fig
end

## CoralPop's morphology -> ADRIA group mapping as its run recorded it; results from before it was
## recorded (2026-10-07) used [[2, 3], [1], [4, 5]].
parse_groups(s) = [parse.(Int, split(m.captures[1], r"[,\s]+"; keepempty=false)) for m in eachmatch(r"\[([0-9,\s]*)\]", s)]
function group_map(e::Engine)
    e.key === :coralpop || return ADRIA_GROUP_MAP[e.key]
    s = get(e.settings, "CORALPOP_MORPH_GROUPS", nothing)
    return s === nothing ? [[2, 3], [1], [4, 5]] : parse_groups(s)
end
const CORALPOP_MAP = (i = findfirst(e -> e.key === :coralpop, engines); i === nothing ? nothing : group_map(engines[i]))

"Engine taxa cover summed onto ADRIA's 5 groups `[t, 5, scenarios]` (CoralPop's mixed morphs split evenly)."
function harmonised(e::Engine)
    (taxa[e.key] === missing || !haskey(ADRIA_GROUP_MAP, e.key)) && return missing
    a = taxa[e.key]
    out = zeros(NY, 5, size(a, 3))
    for (g, targets) in enumerate(group_map(e))
        g > size(a, 2) && break
        for t in targets
            out[:, t, :] .+= a[:, g, :] ./ length(targets)
        end
    end
    return out
end
harm = Dict(e.key => harmonised(e) for e in engines)
const ADRIA_GROUPS = replace.(string.(ADRIA.functional_group_names()), "_" => " ")

figure("D_taxa", "taxa_harmonised", "Cover per ADRIA functional group, every engine mapped onto ADRIA's 5 groups with its interface's own translation table (CoralPop's mixed morphologies split evenly across their target groups; C~scape's brooders added to corymbose non-Acropora; CORSET has no morphological groups and is not shown). An approximation, not an equivalence.") do
    es = [e for e in engines if harm[e.key] !== missing]
    fig = Mk.Figure(; size=(1400, 700))
    for g in 1:5
        ax = Mk.Axis(fig[fld1(g, 3), mod1(g, 3)]; title=ADRIA_GROUPS[g], xlabel="Year", ylabel="Relative cover")
        for e in es
            spread!(ax, e, harm[e.key][:, g, :]; scen_lines=false)
        end
    end
    engine_legend!(fig[3, :], es; orientation=:horizontal)
    fig
end

for (name, t, caption) in (
    ("composition_initial", 1, "Share of first-year cover per harmonised ADRIA group, per engine (median scenario). Shows what each engine starts from: how faithfully it translates the domain's initial group composition. C~scape's first compared year has already had one year of dynamics."),
    ("composition_final", NY, "Share of final-year cover per harmonised ADRIA group, per engine (median scenario). Shows which groups each engine ends up dominated by."))
    figure("D_taxa", name, caption) do
        es = [e for e in engines if harm[e.key] !== missing]
        fig = Mk.Figure(; size=(1000, 500))
        ax = Mk.Axis(fig[1, 1]; xticks=(1:length(es), [e.label for e in es]), ylabel="Share of cover", title="Composition in $(YEARS[t]) (harmonised groups)")
        for (i, e) in enumerate(es)
            v = vec(median(harm[e.key][t, :, :]; dims=2))
            v = sum(v) > 0 ? v ./ sum(v) : v
            Mk.barplot!(ax, fill(i, 5), v; stack=1:5, color=GROUP_COLORS[1:5], gap=0.3, strokecolor=:white, strokewidth=1)
        end
        Mk.Legend(fig[1, 2], [Mk.PolyElement(color=c) for c in GROUP_COLORS[1:5]], ADRIA_GROUPS; framevisible=false)
        fig
    end
end

# =============================================================================
# E. Size structure
# =============================================================================

sizes = Dict{Symbol,Any}()
for e in engines
    f = joinpath(e.dir, "size_structure.jld2")
    isfile(f) && (sizes[e.key] = JLD2.load(f))
end
size_engines = [e for e in engines if haskey(sizes, e.key)]

## C~scape and ReefMod store colony COUNTS per size class, the other engines cover. Counts are dominated
## by small colonies (a 2 cm colony has 1/400 of a 40 cm colony's area), so count shares and cover shares
## are not comparable. Weight counts by the area of a colony at the bin-centre diameter, so every engine's
## size structure below is cover-based (approximate within each bin).
function as_cover!(ss)
    occursin(r"count|density"i, ss["units"]) || return ss
    area_m2 = π .* (ss["d_cm"] ./ 2) .^ 2 ./ 1e4                       # [G, B]
    ss["value"] = ss["value"] .* reshape(area_m2, 1, size(area_m2)..., 1)
    ss["units"] = "colony area, m^2 (from $(ss["units"]) x bin-centre colony area)"
    return ss
end
foreach(as_cover!, values(sizes))

## Every engine bins colonies differently (per group, per functional type, 7 to 100 classes), so pooling
## raw bins gives a saw-tooth. Each (group, bin) amount is assigned to the common diameter class
## containing its bin centre instead; all engines then share one set of classes.
const SIZE_EDGES = [0.0, 2.5, 5.0, 10.0, 20.0, 40.0, 80.0, 160.0, Inf]
const SIZE_LABELS = ["<2.5", "2.5-5", "5-10", "10-20", "20-40", "40-80", "80-160", ">160"]

"Share of the engine's total in each common size class at year index `t`, all groups pooled, median scenario."
function common_size_share(ss, t)
    v = dropdims(median(ss["value"][t, :, :, :]; dims=3); dims=3)    # [G, B]
    d = ss["d_cm"]                                                    # [G, B]
    out = zeros(length(SIZE_LABELS))
    for i in eachindex(v)
        out[clamp(searchsortedlast(SIZE_EDGES, d[i]), 1, length(SIZE_LABELS))] += v[i]
    end
    tot = sum(out)
    return tot > 0 ? out ./ tot : out
end
year_index(ss, yr) = argmin(abs.(ss["years"] .- yr))

figure("E_size", "size_distribution", "Size distribution of every engine on common diameter classes: share of the engine's coral cover in each class (all groups pooled, median scenario) in the first, middle and last year. C~scape's and ReefMod's colony counts are converted to cover (count x area of a colony at the bin-centre diameter), so every engine is on the same basis. CoralBlox's is logged after growth, before that year's disturbance. CORSET has no size structure.") do
    isempty(size_engines) && return nothing
    fig = Mk.Figure(; size=(1500, 480))
    nk = length(SIZE_LABELS)
    for (j, yr) in enumerate(unique((first(YEARS), YEARS[cld(NY, 2)], last(YEARS))))
        ax = Mk.Axis(fig[1, j]; title=string(yr), xticks=(1:nk, SIZE_LABELS), xticklabelrotation=π / 5,
            xlabel="Colony diameter (cm)", ylabel=j == 1 ? "Share of cover" : "")
        for e in size_engines
            ss = sizes[e.key]
            Mk.scatterlines!(ax, 1:nk, common_size_share(ss, year_index(ss, yr)); color=e.color,
                linestyle=e.linestyle, marker=e.marker, markersize=9, linewidth=2)
        end
    end
    engine_legend!(fig[1, 4], size_engines)
    fig
end

figure("E_size", "size_structure_over_time", "How each engine's size structure shifts over the run: share of the engine's coral cover in each common diameter class (rows) per year (columns), median scenario; C~scape's and ReefMod's counts converted to cover. Darker = larger share.") do
    isempty(size_engines) && return nothing
    n = length(size_engines)
    fig = Mk.Figure(; size=(1400, 300 * cld(n, 2)))
    nk = length(SIZE_LABELS)
    local hm
    for (i, e) in enumerate(size_engines)
        ss = sizes[e.key]
        T = size(ss["value"], 1)
        M_ = reduce(hcat, [common_size_share(ss, t) for t in 1:T])     # [classes, T]
        ax = Mk.Axis(fig[fld1(i, 2), mod1(i, 2)]; title=e.label, yticks=(1:nk, SIZE_LABELS), ylabel="Diameter (cm)")
        hm = Mk.heatmap!(ax, ss["years"][1:T], 1:nk, permutedims(M_); colormap=:Blues, colorrange=(0, 1))
    end
    Mk.Colorbar(fig[:, 3], hm; label="Share of cover")
    fig
end

figure("E_size", "mean_colony_diameter", "Cover-weighted mean colony diameter over time (all groups pooled; median over scenarios; C~scape's and ReefMod's counts converted to cover). A rising line = populations shifting to larger colonies (recruitment failure or growth), falling = recovery from recruitment or loss of large colonies.") do
    isempty(size_engines) && return nothing
    fig = Mk.Figure(; size=(1000, 500))
    ax = Mk.Axis(fig[1, 1]; xlabel="Year", ylabel="Mean diameter (cm)", yscale=log10, title="Mean colony diameter")
    for e in size_engines
        ss = sizes[e.key]
        v = ss["value"]                                  # [T, G, B, S]
        d = ss["d_cm"]
        T = min(size(v, 1), NY)
        md = [let w = v[t, :, :, s]; sum(w) > 0 ? sum(w .* d) / sum(w) : NaN end for t in 1:T, s in axes(v, 4)]
        m = [(v = filter(!isnan, md[t, :]); isempty(v) ? NaN : median(v)) for t in 1:T]
        Mk.lines!(ax, ss["years"][1:T], m; color=e.color, linestyle=e.linestyle, linewidth=2, label=e.label)
    end
    engine_legend!(fig[1, 2], size_engines)
    fig
end

# =============================================================================
# F. Disturbance response
# =============================================================================

"Annual proportional change in domain cover `(c[t] - c[t-1]) / c[t-1]` `[T-1, S]`."
annual_change(a) = (a[2:end, :] .- a[1:(end - 1), :]) ./ max.(a[1:(end - 1), :], 1e-9)

figure("F_disturbance", "dhw_and_cover", "Domain DHW (top) and relative cover (bottom, median per engine) on a shared time axis: where each engine's losses and recoveries line up with heat stress.") do
    fig = Mk.Figure(; size=(1100, 650))
    ax1 = Mk.Axis(fig[1, 1]; ylabel="Domain DHW", title=MEMBER)
    Mk.barplot!(ax1, YEARS, domain_dhw; color=:grey60, gap=0.1)
    Mk.hlines!(ax1, [4, 8]; color=:grey30, linestyle=:dash, linewidth=1)
    ax2 = Mk.Axis(fig[2, 1]; ylabel="Relative cover", xlabel="Year")
    for e in engines
        a = dommean[e.key][:relative_cover]
        a === missing || spread!(ax2, e, a; scen_lines=false)
    end
    Mk.linkxaxes!(ax1, ax2)
    Mk.rowsize!(fig.layout, 1, Mk.Relative(0.35))
    engine_legend!(fig[2, 2])
    fig
end

figure("F_disturbance", "cover_change_vs_dhw", "Annual proportional change in domain cover against that year's domain DHW, every scenario-year per engine, with the median per DHW bin (line). Shows each engine's effective heat-stress response curve.") do
    fig = Mk.Figure(; size=(1100, 550))
    ax = Mk.Axis(fig[1, 1]; xlabel="Domain DHW in the year (°C-weeks)", ylabel="Proportional change in cover",
        title="Cover change vs heat stress")
    bins = 0:2:ceil(maximum(domain_dhw) / 2) * 2
    for e in engines
        a = dommean[e.key][:relative_cover]
        a === missing && continue
        ch = annual_change(a)
        x = repeat(domain_dhw[2:end], 1, size(ch, 2))
        Mk.scatter!(ax, vec(x), vec(ch); color=(e.color, 0.25), marker=e.marker, markersize=5)
        mids, meds = Float64[], Float64[]
        for (lo, hi) in zip(bins[1:(end - 1)], bins[2:end])
            sel = (lo .<= vec(x) .< hi)
            count(sel) >= 3 || continue
            push!(mids, (lo + hi) / 2)
            push!(meds, median(vec(ch)[sel]))
        end
        Mk.scatterlines!(ax, mids, meds; color=e.color, linestyle=e.linestyle, marker=e.marker, markersize=9, linewidth=2)
    end
    Mk.hlines!(ax, [0.0]; color=:grey40, linestyle=:dash)
    engine_legend!(fig[1, 2])
    fig
end

## Recovery: median annual proportional growth in years whose domain DHW is below 4 (no bleaching
## expected), per scenario.
low_dhw = findall(<(4), domain_dhw[2:end])
recovery = Dict(e.key => (dommean[e.key][:relative_cover] === missing || isempty(low_dhw) ? missing :
    vec(median(annual_change(dommean[e.key][:relative_cover])[low_dhw, :]; dims=1))) for e in engines)

figure("F_disturbance", "recovery_rate", "Recovery rate: median annual proportional growth of domain cover in years with domain DHW < 4 ($(length(low_dhw)) of $(NY - 1) years in this member), one dot per scenario.") do
    isempty(low_dhw) && return nothing
    es = [e for e in engines if recovery[e.key] !== missing]
    fig = Mk.Figure(; size=(900, 480))
    ax = Mk.Axis(fig[1, 1]; xticks=(1:length(es), [e.label for e in es]), ylabel="Annual proportional growth",
        title="Recovery between heat-stress years")
    for (i, e) in enumerate(es)
        v = recovery[e.key]
        Mk.boxplot!(ax, fill(i, length(v)), v; color=(e.color, 0.35), strokecolor=e.color, width=0.5)
        Mk.scatter!(ax, fill(i, length(v)), v; color=e.color, marker=e.marker, markersize=8)
    end
    Mk.hlines!(ax, [0.0]; color=:grey40, linestyle=:dash)
    fig
end

const COVER_THRESHOLDS = (0.1, 0.2, 0.3)
figure("F_disturbance", "years_above_threshold", "Number of years (of $NY) with domain relative cover above 10%, 20% and 30% of carrying capacity, per scenario.") do
    fig = Mk.Figure(; size=(1300, 450))
    for (j, th) in enumerate(COVER_THRESHOLDS)
        ax = Mk.Axis(fig[1, j]; xticks=(1:length(engines), [e.label for e in engines]), xticklabelrotation=π / 6,
            title="Cover > $(Int(100th))%", ylabel=j == 1 ? "Years" : "")
        for (i, e) in enumerate(engines)
            a = dommean[e.key][:relative_cover]
            a === missing && continue
            v = vec(sum(a .> th; dims=1))
            Mk.scatter!(ax, i .+ 0.15 .* (rand(length(v)) .- 0.5), v; color=e.color, marker=e.marker, markersize=9)
            Mk.scatter!(ax, [i], [median(v)]; marker=:hline, markersize=30, color=:black)
        end
        Mk.ylims!(ax, -1, NY + 1)
    end
    fig
end

# =============================================================================
# G. Spatial
# =============================================================================

loc_median(e, m=:relative_cover) = loc[e.key][m] === missing ? missing : dropdims(median(loc[e.key][m]; dims=3); dims=3)   # [T, L]

figure("G_spatial", "site_trajectories", "Relative cover at each site, median over scenarios per engine. Connectivity in the Test domain is sparse (only reef_05 and reef_06 export larvae; Kora and ReefMod do not use connectivity at all), which shapes site-to-site differences.") do
    fig = Mk.Figure(; size=(1500, 300 * cld(NL, 5)))
    axs = []
    for l in 1:NL
        ax = Mk.Axis(fig[fld1(l, 5), mod1(l, 5)]; title=LOC_IDS[l], titlesize=12)
        push!(axs, ax)
        for e in engines
            a = loc_median(e)
            a === missing || Mk.lines!(ax, YEARS, a[:, l]; color=e.color, linestyle=e.linestyle, linewidth=1.8)
        end
    end
    Mk.linkyaxes!(axs...)
    engine_legend!(fig[cld(NL, 5) + 1, :]; orientation=:horizontal)
    fig
end

"Site centroids (lon, lat) from the domain geometry."
function site_xy()
    c = ADRIA.centroids(dom)
    return [p[1] for p in c], [p[2] for p in c]
end

figure("G_spatial", "map_final_cover", "Final-year relative cover per site (median scenario) for each engine, on a common colour scale. Dot positions are site centroids.") do
    x, y = site_xy()
    fin = Dict(e.key => loc_median(e)[end, :] for e in engines if loc_median(e) !== missing)
    hi = maximum(maximum(v) for v in values(fin))
    n = length(fin)
    fig = Mk.Figure(; size=(1300, 360 * cld(n, 3)))
    k = 0
    for e in engines
        haskey(fin, e.key) || continue
        k += 1
        ax = Mk.Axis(fig[fld1(k, 3), mod1(k, 3)]; title=e.label, aspect=Mk.DataAspect(), xticklabelsvisible=false, yticklabelsvisible=false)
        Mk.scatter!(ax, x, y; color=fin[e.key], colorrange=(0, hi), colormap=:Blues, markersize=22, strokecolor=:grey30, strokewidth=0.8)
    end
    Mk.Colorbar(fig[:, 4]; colorrange=(0, hi), colormap=:Blues, label="Relative cover, $(last(YEARS))")
    fig
end

figure("G_spatial", "map_between_engine_spread", "Between-engine disagreement per site: standard deviation across engines of final-year median relative cover (left) and the coefficient of variation (right).") do
    x, y = site_xy()
    M_ = reduce(hcat, [loc_median(e)[end, :] for e in engines if loc_median(e) !== missing])   # [L, engines]
    sd = vec(std(M_; dims=2))
    cv = sd ./ max.(vec(mean(M_; dims=2)), 1e-9)
    fig = Mk.Figure(; size=(1100, 450))
    for (j, (v, lab)) in enumerate(((sd, "SD across engines"), (cv, "CV across engines")))
        ax = Mk.Axis(fig[1, 2j - 1]; title=lab, aspect=Mk.DataAspect(), xticklabelsvisible=false, yticklabelsvisible=false)
        Mk.scatter!(ax, x, y; color=v, colormap=:Oranges, markersize=22, strokecolor=:grey30, strokewidth=0.8)
        Mk.Colorbar(fig[1, 2j]; colormap=:Oranges, colorrange=(0, max(maximum(v), 1e-9)))
    end
    fig
end

"Kendall's tau-b between two vectors."
function kendall_tau(a, b)
    n = length(a)
    c = d = ta = tb = 0
    for i in 1:(n - 1), j in (i + 1):n
        s1, s2 = sign(a[i] - a[j]), sign(b[i] - b[j])
        if s1 == 0 && s2 == 0
        elseif s1 == 0
            ta += 1
        elseif s2 == 0
            tb += 1
        elseif s1 == s2
            c += 1
        else
            d += 1
        end
    end
    den = sqrt((c + d + ta) * (c + d + tb))
    return den == 0 ? NaN : (c - d) / den
end

"Engine x engine matrix of `f(vector_i, vector_j)`."
pairwise_matrix(vs, f) = [f(vs[i], vs[j]) for i in eachindex(vs), j in eachindex(vs)]

function matrix_heatmap(Mx, labels, title, cbar_label)
    fig = Mk.Figure(; size=(720, 620))
    n = length(labels)
    ax = Mk.Axis(fig[1, 1]; xticks=(1:n, labels), yticks=(1:n, labels), xticklabelrotation=π / 4, title=title,
        aspect=Mk.DataAspect(), yreversed=true)
    hm = Mk.heatmap!(ax, 1:n, 1:n, Mx; colormap=:RdBu, colorrange=(-1, 1))
    for i in 1:n, j in 1:n
        isnan(Mx[i, j]) || Mk.text!(ax, i, j; text=@sprintf("%.2f", Mx[i, j]), align=(:center, :center), fontsize=12,
            color=abs(Mx[i, j]) > 0.6 ? :white : :black)
    end
    Mk.Colorbar(fig[1, 2], hm; label=cbar_label)
    return fig
end

figure("G_spatial", "site_rank_agreement", "Do engines agree on WHICH sites end up in best condition? Kendall's tau between engines' site rankings by final-year median relative cover (1 = identical ranking, 0 = unrelated, -1 = reversed). Relevant for site prioritisation.") do
    es = [e for e in engines if loc_median(e) !== missing]
    vs = [loc_median(e)[end, :] for e in es]
    matrix_heatmap(pairwise_matrix(vs, kendall_tau), [e.label for e in es], "Site-ranking agreement, $(last(YEARS))", "Kendall's tau")
end

figure("G_spatial", "site_final_cover_by_engine", "Final-year relative cover at each site, each engine (median over scenarios): the between-engine spread at every site.") do
    fig = Mk.Figure(; size=(1100, 450))
    ax = Mk.Axis(fig[1, 1]; xticks=(1:NL, LOC_IDS), xticklabelrotation=π / 4, ylabel="Relative cover", title="Final-year cover by site")
    off = range(-0.3, 0.3; length=length(engines))
    for (i, e) in enumerate(engines)
        a = loc_median(e)
        a === missing || Mk.scatter!(ax, (1:NL) .+ off[i], a[end, :]; color=e.color, marker=e.marker, markersize=10)
    end
    engine_legend!(fig[1, 2])
    fig
end

# =============================================================================
# H. Agreement and uncertainty
# =============================================================================

figure("H_agreement", "trajectory_correlation", "Pearson correlation between engines' median domain-cover trajectories (does cover rise and fall in the same years?).") do
    es = [e for e in engines if dommean[e.key][:relative_cover] !== missing]
    vs = [med_iqr(dommean[e.key][:relative_cover])[1] for e in es]
    matrix_heatmap(pairwise_matrix(vs, (a, b) -> cor(a, b)), [e.label for e in es], "Trajectory correlation", "Pearson r")
end

figure("H_agreement", "variance_partition", "Where does the uncertainty come from? Within-engine variance (mean over engines of the variance across each engine's sampled scenarios = parameter uncertainty) vs between-engine variance (variance of the engines' mean trajectories = structural uncertainty) of domain relative cover over time.") do
    es = [e for e in engines if dommean[e.key][:relative_cover] !== missing]
    within = vec(mean(reduce(hcat, [vec(var(dommean[e.key][:relative_cover]; dims=2)) for e in es]); dims=2))
    between = vec(var(reduce(hcat, [vec(mean(dommean[e.key][:relative_cover]; dims=2)) for e in es]); dims=2))
    fig = Mk.Figure(; size=(1000, 480))
    ax = Mk.Axis(fig[1, 1]; xlabel="Year", ylabel="Variance of relative cover", title="Parameter vs structural uncertainty")
    Mk.lines!(ax, YEARS, within; color=:grey30, linewidth=2, linestyle=:dash, label="Within engines (parameters)")
    Mk.lines!(ax, YEARS, between; color=:black, linewidth=2.5, label="Between engines (structure)")
    Mk.axislegend(ax; position=:lt, framevisible=false)
    fig
end

# =============================================================================
# I. Sensitivity
# =============================================================================

const PAWN_S = clamp(minimum(values(n_scen)) ÷ 2, 1, 10)
const SA_NOTE = minimum(values(n_scen)) < 800 ? " ILLUSTRATIVE ONLY: PAWN needs N/S > 80 (N >= 800 at S = 10); this run has N = $(minimum(values(n_scen)))." : ""

for e in engines
    a = dommean[e.key][:relative_cover]
    (a === missing || ncol(e.factors) == 0) && continue
    y = a[end, :]
    figure("I_sensitivity", "pawn_$(e.key)", "PAWN sensitivity of $(e.label)'s final-year domain cover to its own sampled factors (S = $PAWN_S).$SA_NOTE") do
        si = ADRIAanalysis.sensitivity.pawn(e.factors, y; S=PAWN_S)
        ADRIA.viz.pawn(si; axis_opts=Dict{Symbol,Any}(:title => "$(e.label): PAWN, final-year cover"))
    end
    figure("I_sensitivity", "outcome_map_$(e.key)", "Outcome map for $(e.label): final-year domain cover across the range of its four most influential factors (by PAWN median).$SA_NOTE") do
        ## pawn returns a [factors, Si] cube: factors in the factor table's column order, Si = min, lb,
        ## mean, median, ub, max, std, cv
        si = ADRIAanalysis.sensitivity.pawn(e.factors, y; S=PAWN_S)
        med = Array(si)[:, 4]
        top = Symbol.(names(e.factors))[sortperm(med; rev=true)][1:min(4, length(med))]
        ADRIA.viz.outcome_map(e.factors, y, top)
    end
    figure("I_sensitivity", "tsa_$(e.key)", "Temporal sensitivity analysis for $(e.label): PAWN sensitivity of domain cover to each factor, year by year (how the dominant factors change over the run).$SA_NOTE") do
        si = ADRIAanalysis.sensitivity.tsa(e.factors, a)
        ADRIA.viz.tsa(e.rs, si; axis_opts=Dict{Symbol,Any}(:title => "$(e.label): temporal sensitivity of domain cover"))
    end
end

# =============================================================================
# J. CORSET extras
# =============================================================================

corset = findfirst(e -> e.key === :corset, engines)
if corset !== nothing
    ec = engines[corset]
    extra = [(k, replace(String(k), "_" => " ")) for k in (:macroturf_cover, :macroalgae_cover, :grazed_eac_cover,
        :herbivore_biomass, :small_piscivore_biomass, :large_piscivore_biomass, :urchin_biomass) if haskey(ec.rs.outcomes, k)]
    figure("J_corset", "corset_state", "CORSET's non-coral state, domain mean (median, IQR, scenarios). Consumer biomass starts from PLACEHOLDER values (no ADRIA domain carries fish or urchin biomass), so absolute levels are not meaningful; relative change is.") do
        isempty(extra) && return nothing
        fig = Mk.Figure(; size=(1400, 330 * cld(length(extra), 3)))
        for (i, (k, t)) in enumerate(extra)
            a = Float64.(Array(ec.rs.outcomes[k]))
            a = ndims(a) == 3 ? domain_mean(a) : a
            ax = Mk.Axis(fig[fld1(i, 3), mod1(i, 3)]; title=t)
            spread!(ax, ec, a)
        end
        fig
    end
end

# =============================================================================
# K. Benchmarks
# =============================================================================

## Slurm accounting (written by hpc/jobs/collect_benchmarks.sh): MaxRSS / TotalCPU / Elapsed per job.
"Slurm sizes like 1234K / 5.2G -> GB."
function slurm_gb(s)
    s = strip(string(s))
    (isempty(s) || s == "missing") && return NaN
    m = match(r"^([\d.]+)([KMGT]?)", s)
    m === nothing && return NaN
    return parse(Float64, m[1]) * Dict("" => 1e-9, "K" => 1e-6, "M" => 1e-3, "G" => 1.0, "T" => 1e3)[m[2]]
end
"Slurm durations [D-]HH:MM:SS[.ms] or MM:SS.ms -> seconds."
function slurm_s(s)
    s = strip(string(s))
    isempty(s) && return NaN
    days = 0
    occursin('-', s) && ((d, s) = split(s, '-'; limit=2); days = parse(Int, d))
    parts = parse.(Float64, split(s, ':'))
    secs = foldl((acc, x) -> acc * 60 + x, parts)
    return days * 86400 + secs
end

bench_rows = DataFrame[]
for e in engines
    f = joinpath(e.dir, "benchmark.csv")
    isfile(f) && push!(bench_rows, read_csv(f))
end
bench = isempty(bench_rows) ? DataFrame() : reduce((a, b) -> vcat(a, b; cols=:union), bench_rows)

sacct_file = joinpath(PROJECT_DIR, "slurm_accounting.csv")
if !isempty(bench) && isfile(sacct_file)
    ## sacct -P output, '|' separated; the engine job steps are named cmr_<engine>.
    lines = filter(!isempty, readlines(sacct_file))
    hdr = split(lines[1], '|')
    rows = [Dict(zip(hdr, split(l, '|'))) for l in lines[2:end]]
    for col in (:sacct_maxrss_gb, :sacct_cpu_s, :sacct_elapsed_s)
        bench[!, col] = fill(NaN, nrow(bench))
    end
    for (i, eng) in enumerate(bench.engine)
        ## the job's steps (batch, extern, ...) all share the parent job id; take the max MaxRSS / CPU
        jid = string(bench.slurm_job_id[i])
        jid = endswith(jid, ".0") ? jid[1:(end - 2)] : jid
        own(id) = id == jid || startswith(id, jid * ".")              # the job and its steps (1234.batch, ...)
        mine = isempty(jid) ? Dict{String,String}[] : [r for r in rows if own(get(r, "JobID", ""))]
        if isempty(mine)
            mine = [r for r in rows if get(r, "JobName", "") == "cmr_$(eng)"]
        end
        isempty(mine) && continue
        bench.sacct_maxrss_gb[i] = maximum(slurm_gb(get(r, "MaxRSS", "")) for r in mine)
        bench.sacct_cpu_s[i] = maximum(slurm_s(get(r, "TotalCPU", "")) for r in mine)
        bench.sacct_elapsed_s[i] = maximum(slurm_s(get(r, "Elapsed", "")) for r in mine)
    end
end

bcol(k) = Dict(e.key => e for e in engines)[Symbol(k)]
if !isempty(bench)
    isempty(ONLY_FIGURES) && write_csv(joinpath(PROJECT_DIR, "benchmark_summary.csv"), [NamedTuple(r) for r in eachrow(bench)])

    "Horizontal bars of one benchmark column per engine (log axis when values span > 20x)."
    function bench_panel!(pos, col, title, unit)
        hasproperty(bench, col) || return false
        v = Float64[x isa Number ? x : NaN for x in bench[!, col]]
        ok = findall(x -> isfinite(x) && x > 0, v)
        isempty(ok) && return false
        es = [bcol(k) for k in bench.engine[ok]]
        logx = maximum(v[ok]) / minimum(v[ok]) > 20
        ax = Mk.Axis(pos; yticks=(1:length(ok), [e.label for e in es]), xlabel=unit, title=title,
            xscale=logx ? log10 : identity)
        lo = logx ? minimum(v[ok]) / 2 : 0.0
        Mk.barplot!(ax, 1:length(ok), v[ok]; direction=:x, color=[e.color for e in es], gap=0.35, fillto=lo)
        fmt(x) = x >= 100 ? @sprintf(" %.0f", x) : x >= 1 ? @sprintf(" %.1f", x) : @sprintf(" %.2g", x)
        Mk.text!(ax, v[ok], 1:length(ok); text=fmt.(v[ok]), align=(:left, :center), fontsize=11)
        Mk.xlims!(ax, lo, logx ? maximum(v[ok]) * 6 : maximum(v[ok]) * 1.18)   # room for the labels
        return true
    end

    figure("K_benchmarks", "time", "Run time per engine: package loading, the 1-scenario warm-up (time to first result incl. compilation and worker/MATLAB start-up), steady-state wall time per scenario, and per scenario-site-year. Test run: 1 CPU per engine, engines serial. ReefMod's per-scenario time includes a MATLAB start-up per scenario.") do
        fig = Mk.Figure(; size=(1300, 750))
        bench_panel!(fig[1, 1], :t_load_s, "Package load", "seconds")
        bench_panel!(fig[1, 2], :t_first_s, "First result (warm-up, 1 scenario)", "seconds")
        bench_panel!(fig[2, 1], :t_per_scen_s, "Steady state, per scenario", "seconds")
        bench_panel!(fig[2, 2], :t_per_scen_site_year_ms, "Per scenario · site · year", "milliseconds")
        fig
    end

    figure("K_benchmarks", "memory_disk", "Peak memory (Slurm MaxRSS of the whole job when available -- includes worker processes and MATLAB -- else the Julia process's own peak), CPU time per scenario (Slurm TotalCPU / N) and output size on disk per scenario.") do
        fig = Mk.Figure(; size=(1300, 750))
        hasproperty(bench, :sacct_maxrss_gb) && any(isfinite, bench.sacct_maxrss_gb) ?
            bench_panel!(fig[1, 1], :sacct_maxrss_gb, "Peak memory (Slurm MaxRSS)", "GB") :
            bench_panel!(fig[1, 1], :maxrss_gb, "Peak memory (Julia process)", "GB")
        if hasproperty(bench, :sacct_cpu_s)
            bench[!, :sacct_cpu_per_scen_s] = bench.sacct_cpu_s ./ bench.n_scenarios
            bench_panel!(fig[1, 2], :sacct_cpu_per_scen_s, "CPU time per scenario (Slurm)", "core-seconds")
        else
            bench_panel!(fig[1, 2], :t_main_alloc_gb, "Allocations, main run", "GB")
        end
        bench_panel!(fig[2, 1], :disk_per_scen_gb, "Output on disk per scenario", "GB")
        bench_panel!(fig[2, 2], :t_main_compile_s, "Compilation inside the main run", "seconds")
        fig
    end
end

# =============================================================================
# Summary table and figure index
# =============================================================================

rows = NamedTuple[]
for e in engines, (m, t, _) in METRICS
    plottable(e, m) || continue
    a = dommean[e.key][m]
    for ti in unique((1, findfirst(==(2050), YEARS), NY))
        ti === nothing && continue
        v = a[ti, :]
        push!(rows, (engine=e.label, metric=String(m), year=YEARS[ti], median=median(v),
            q25=q(v, 0.25), q75=q(v, 0.75), min=minimum(v), max=maximum(v)))
    end
end
isempty(ONLY_FIGURES) && write_csv(joinpath(PROJECT_DIR, "summary.csv"), rows)

## C~scape's demography check as that run recorded it (results from before the check existed carry none --
## and were run with survival models whose cloglog link had been exported as logit, so their cover is too low).
const CSCAPE_DEMOGRAPHY = let i = findfirst(e -> e.key === :cscape, engines)
    i === nothing ? nothing : get(engines[i].settings, "CSCAPE_DEMOGRAPHY", "NOT RECORDED -- run before the demography check existed")
end
const CSCAPE_DEMOGRAPHY_OK = CSCAPE_DEMOGRAPHY !== nothing && startswith(CSCAPE_DEMOGRAPHY, "verified")
CSCAPE_DEMOGRAPHY === nothing || CSCAPE_DEMOGRAPHY_OK ||
    @warn "C~scape demography is not verified against brms -- treat its results with caution" CSCAPE_DEMOGRAPHY

## CORSET's disturbance switches as that run recorded them (older results don't record them: off).
corset_switch(name) = corset !== nothing && get(engines[corset].settings, name, false) === true
const CORSET_DHW = corset_switch("CORSET_USE_DHW_AS_BLEACHING") ?
    "on: each year whose across-site mean DHW > 4 is one bleaching event at every site (a random 0-50 % coral loss per site, not scaled by DHW)" :
    "off (CORSET_USE_DHW_AS_BLEACHING)"
const CORSET_CYC = corset_switch("CORSET_USE_CYCLONES_AS_HURRICANES") ? "on" : "off (CORSET_USE_CYCLONES_AS_HURRICANES)"
## Where CoralPop's and Kora's larvae came from (CMR_LARVAE_VIA_CONNECTIVITY, as each run recorded it).
via_matrix(key) = (i = findfirst(e -> e.key === key, engines);
    i !== nothing && get(engines[i].settings, "LARVAE_VIA_CONNECTIVITY", false) === true)
const CORALPOP_LARVAE = via_matrix(:coralpop) ? "all spawner larvae via the matrix (own retention off)" :
    "between sites via the matrix; own larvae via CoralPop's retention"
const KORA_LARVAE = via_matrix(:kora) ? "all larvae via the matrix (self_seed off)" : "none: closed, self-seeding sites"

## Which processes each engine represents in THIS comparison (from 07's header and the interfaces).
const CAPABILITIES = [
    ("CoralBlox", "5 ADRIA groups x 7", "yes", "yes", "yes", "yes"),
    ("C~scape", "6 types, 100 IPM classes", "yes", "yes", "yes", "yes"),
    ("CoralPop", "3 morphologies x bins", CORALPOP_LARVAE, "yes", "yes (via domain)", "yes"),
    ("Kora", "5 ADRIA groups x 7", KORA_LARVAE, "yes", "no (no cyclone pathway)", "yes"),
    ("CORSET", "none (brooding/spawning)", "yes", CORSET_DHW, CORSET_CYC, "no"),
    ("ReefMod", "6 groups, 11 bins (counts)", "yes (ADRIA's matrix)", "yes", "switched off", "provisional"),
]

## A filtered run only adds (or refreshes) its own lines in the existing FIGURES.md, in the figure's section,
## after the last line of that section; everything else in the file is kept as it was.
function merge_figure_lines!(path, figs)
    isfile(path) || return false
    lines = readlines(path)
    for f in figs
        line = "- [`$(basename(f.file))`]($(f.file)) -- $(f.caption)"
        key = "- [`$(basename(f.file))`]"
        i = findfirst(l -> startswith(l, key), lines)
        if i !== nothing
            lines[i] = line
            continue
        end
        h = findfirst(==("### " * replace(f.section, "_" => " ", count=1)), lines)
        h === nothing && return false
        j = h + 1
        while j + 1 <= length(lines) && !startswith(lines[j + 1], "###") && !startswith(lines[j + 1], "## ")
            j += 1
        end
        while j > h && isempty(lines[j]) j -= 1 end                          # insert before the blank line
        insert!(lines, j + 1, line)
    end
    write(path, join(lines, "\n") * "\n")
    return true
end

if !isempty(ONLY_FIGURES)
    ok = merge_figure_lines!(joinpath(PROJECT_DIR, "FIGURES.md"), figure_index)
    ok || @warn "FIGURES.md missing or without the figures' section; not updated (run 09 without CMR_ONLY_FIGURES)"
    @info "Selected figures written" PROJECT_DIR ONLY_FIGURES n_figures = length(figure_index) failed = length(failed_figures)
    isempty(failed_figures) || (@warn "Some figures failed" failed_figures; exit(1))
    exit(0)
end

open(joinpath(PROJECT_DIR, "FIGURES.md"), "w") do io
    println(io, "# Six-engine comparison: $(cfg.project)\n")
    println(io, "Generated $(Dates.format(now(), "yyyy-mm-dd HH:MM")) by `09_ANALYSE_SIX_ENGINES.jl`.\n")
    println(io, "- Domain: `$(ref.settings["DOMAIN_NAME"])`, RCP $(cfg.rcp), years $(first(YEARS))-$(last(YEARS)), $NL sites")
    println(io, "- Climate: $MEMBER, cyclone member $(cfg.cyclone_idx) (domain mean DHW $(round(mean(domain_dhw); digits=2)); $(count(>(4), domain_dhw)) of $NY years above 4)")
    println(io, "- Scenarios: $NSTR")
    isempty(missing_engines) || println(io, "- **Not included** (no saved result): $(join(missing_engines, ", "))")
    isempty(invalid_shelter) || println(io, "- Shelter volume / indices omitted for: $(join([ENGINE_INFO[k].label for k in invalid_shelter], ", ")) (not computed, or outside [0, 1])")
    println(io, "\n## Caveats\n")
    println(io, "- With few scenarios per engine, spreads and every sensitivity figure are illustrative only.$(SA_NOTE)")
    println(io, "- C~scape stores its initial state one year earlier; its first compared year has had one year of dynamics (see 07).")
    CSCAPE_DEMOGRAPHY === nothing || println(io, CSCAPE_DEMOGRAPHY_OK ? "- C~scape demography: $CSCAPE_DEMOGRAPHY." :
        "- **WARNING -- C~scape demography NOT verified against brms** ($CSCAPE_DEMOGRAPHY). Results produced with the 260903 files before this check had survival (cloglog models exported as logit) too low, so C~scape's cover was far below carrying capacity and its massives died out: do not trust them; re-run C~scape.")
    println(io, "- Kora has no cyclone pathway; ReefMod's cyclones and CoTS are switched off. CORSET heat stress: $CORSET_DHW; cyclones: $CORSET_CYC. The Test domain's cyclone data is all zero anyway.")
    CORALPOP_MAP === nothing || println(io, "- CoralPop morphologies from ADRIA groups: ",
        join(["$m = " * join(ADRIA_GROUPS[g], " + ") for (m, g) in zip(CORALPOP_MORPHS, CORALPOP_MAP)], "; "),
        ". Its brooders are the part of Branching CoralPop's upstream split assigns (Branching - Foliose / 9.9, recomputed yearly).")
    println(io, "- Larval supply: CoralBlox and C~scape take larvae only from the domain's connectivity matrix (C~scape's brooders also retain locally). CoralPop: $CORALPOP_LARVAE. Kora: $KORA_LARVAE. Sites the matrix does not supply get no recruits from matrix-only engines.")
    println(io, "- Kora simulates at most $(round(Int, KoraInterface.KORA_SECTION_AREA_M2)) m² of habitable area per site and scales its cover up to the whole site (KoraInterface's section extrapolation), so its time and memory (K) are for those sections, not whole sites.")
    let i = findfirst(e -> e.key === :coralblox, engines)
        if i !== nothing
            cp = get(engines[i].settings, "CALIBRATION_PATH", "")
            println(io, "- CoralBlox parameters: ", isempty(cp) ? "ADRIA's built-in defaults (no calibration file)" :
                "calibration `$(basename(cp))` (the file only affects CoralBlox; it is not compared across engines)", ".")
        end
    end
    println(io, "- CORSET's consumer biomass is a placeholder; ReefMod's juveniles are provisional (see ReefModInterface indicators.jl).")
    println(io, "- The harmonised-taxa mapping is each interface's own translation table, an approximation.")
    println(io, "- Size structure (E) is compared as cover for every engine: C~scape's and ReefMod's colony counts are weighted by the area of a colony at each bin's centre diameter (approximate within a bin). Compare shapes (shares), not amounts.")
    println(io, "\n## What each engine represents here\n")
    println(io, "| Engine | Size structure | Connectivity | DHW bleaching | Cyclones | Juveniles |\n|---|---|---|---|---|---|")
    for c in CAPABILITIES
        println(io, "| ", join(c[[1, 2, 3, 4, 5, 6]], " | "), " |")
    end
    println(io, "\n## Figures\n")
    for sec in unique(f.section for f in figure_index)
        println(io, "### ", replace(sec, "_" => " ", count=1), "\n")
        for f in figure_index
            f.section == sec && println(io, "- [`$(basename(f.file))`]($(f.file)) -- $(f.caption)")
        end
        println(io)
    end
    if !isempty(failed_figures)
        println(io, "## Figures that failed\n")
        foreach(s -> println(io, "- ", s), failed_figures)
    end
end

@info "Comparison written" PROJECT_DIR n_figures = length(figure_index) failed = length(failed_figures)
isempty(failed_figures) || @warn "Some figures failed" failed_figures
