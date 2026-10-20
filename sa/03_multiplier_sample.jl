#=
03_multiplier_sample.jl -- preview what the 13 domain multipliers do to each synthetic domain.

    julia +release --project=Code/sa Code/sa/03_multiplier_sample.jl

A standalone Sobol' sample of the multipliers (N = 1024, QC only -- the values actually run come from
each model's joint design in 04_design.jl) is applied to each domain with `apply_domain_multipliers`
(sa_common.jl), and the resulting domain attributes are written to
domains/<domain>/sa_inputs/multiplier_sample.csv plus figures/multipliers_*.png.
=#

using CSV, CairoMakie
include(joinpath(@__DIR__, "sa_common.jl"))

const N_PREVIEW = 1024

Random.seed!(1234)
const BANDS = multiplier_bands()
println("multiplier bands from ", MULTIPLIER_BANDS_FILE, ": ",
    join(["$n " * (BANDS[n][1] == BANDS[n][2] ? "=$(BANDS[n][1])" : "$(BANDS[n])") for n in MULTIPLIER_NAMES], ", "))
M = ADRIA.sample(multiplier_spec(BANDS), N_PREVIEW)

"Domain attributes after applying one row of multipliers (area-weighted means over reefs)."
function attributes(dom)
    ld = dom.loc_data
    a = Float64.(ld.area)
    w = a ./ sum(a)
    cover = Array(dom.init_coral_cover)                        # [groups*sizes, locations]
    g_tot = [vec(sum(cover[((g - 1) * N_SIZE_CLASSES + 1):(g * N_SIZE_CLASSES), :]; dims=1)) for g in 1:5]
    total = sum(g_tot)
    C = Matrix(dom.conn.data)
    rows = vec(sum(C; dims=2))
    diagC = [C[i, i] for i in axes(C, 1)]
    return (;
        depth_m=sum(w .* ld.depth_med), k=sum(w .* ld.k), area_m2=mean(a),
        cover_total=sum(w .* total),                                       # relative to k
        (Symbol("share_fg$g") => sum(w .* g_tot[g]) / max(sum(w .* total), eps()) for g in 1:5)...,
        ub=mean(Array(dom.wave_scens)[1, :, 1]),
        self_recruitment=mean(diagC), exports=mean(rows .- diagC), retention=mean(rows),
        max_row_sum=maximum(rows),
    )
end

for domain in SA_DOMAINS
    dom = load_sa_domain(domain)
    base = attributes(dom)
    out = DataFrame()
    n_capped = 0
    for r in eachrow(M)
        d, info = apply_domain_multipliers(dom, r)
        n_capped += info.n_conn_capped > 0
        push!(out, merge(NamedTuple(r), attributes(d), (; n_conn_capped=info.n_conn_capped)))
    end
    dir = sa_inputs_dir(domain)
    CSV.write(joinpath(dir, "multiplier_sample.csv"), out)

    cols = [:depth_m, :k, :area_m2, :cover_total, :share_fg1, :share_fg3, :ub, :self_recruitment,
        :exports, :retention, :max_row_sum]
    fig = Figure(; size=(1200, 650))
    for (i, c) in enumerate(cols)
        ax = Axis(fig[fld(i - 1, 4) + 1, mod(i - 1, 4) + 1]; title=string(c))
        x = out[!, c]
        allequal(x) ? vlines!(ax, [first(x)]; color=:steelblue, linewidth=3) :   # a fixed attribute
                      hist!(ax, x; bins=30, color=(:steelblue, 0.7))
        vlines!(ax, [getfield(base, c)]; color=:black, linestyle=:dash)
    end
    Label(fig[0, :], "$domain: domain attributes under $N_PREVIEW multiplier draws (dashed = unmodified domain)")
    save(joinpath(mkpath(joinpath(dir, "figures")), "multipliers_attributes.png"), fig; px_per_unit=2)
    println("$domain: multiplier_sample.csv written; connectivity rows capped at 1 in $n_capped of $N_PREVIEW draws")
end
