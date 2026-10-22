#=
adapters/cscape.jl -- C~scape's state summary (sa_outputs.jl) from one run's in-memory output.

Included by pf_common.jl after models/cscape.jl. Every model adapter defines:

    PF_CLASS_OF_GROUP   harmonised class of each native group (sa_outputs.jl `native_group_class`)
    PF_THERMAL_GROUPS   names of the "thermal_mean" columns (empty: no heritable tolerance)
    state_summary(r)    the state summary of one run, from `run_row`'s result `r`

C~scape: 6 functional types (CSCAPE_FUNCTIONAL_TYPES), mapped onto ADRIA's groups by the interface's own
table (CSCAPE_FT_TO_ADRIA_GROUP = [1, 2, 3, 4, 5, 3]: the brooder counts as corymbose non-Acropora). Cover,
juveniles and diversity are the ADRIAIndicators definitions on `to_adria_format`'s cover (the combined
intervention slot, relative to kappa x area), as CscapeInterface's `calculate_indicators` computes them,
without the indicators the outputs do not use. Juveniles are C~scape's own juvenile size classes
(`is_juvenile`). Thermal tolerance: per FT, the cover-weighted mean of C~scape's tolerance classes (offset
DHW, `CSCAPE_SA_HEAT_CLASSES` from models/cscape.jl), cover in m2 summed over sites, as Code/engines.jl's
`_cscape_thermal`.

Needs from models/cscape.jl: `CScape` and `CscapeInterface` loaded, `CSCAPE_SA_HEAT_CLASSES` (the tolerance
class values of the runs' `HeatToleranceGroups`), and `run_row` returning `out`, `reef_k`, `reef_area`.

C~scape stores its initial state under the year before the start (2024); only SA_YEARS are kept, as
Code/engines.jl's `align_cscape` does, so year 1 is 2025 as for the other models.
=#

import ADRIAIndicators

const PF_CLASS_OF_GROUP = native_group_class([[g] for g in CscapeInterface.CSCAPE_FT_TO_ADRIA_GROUP])
const PF_THERMAL_GROUPS = String.(CscapeInterface.CSCAPE_FUNCTIONAL_TYPES)

"Rows of a C~scape output's year axis that are SA_YEARS (error if any is missing)."
function _cscape_sa_rows(years)
    rows = [findfirst(==(y), years) for y in SA_YEARS[1]:SA_YEARS[2]]
    any(isnothing, rows) && error("C~scape output covers years $(first(years))-$(last(years)), not all of SA_YEARS $SA_YEARS")
    return Int.(rows)
end

"""
    cscape_thermal(out, rows) -> (mean [T, ft], weight [T, ft])

Cover-weighted mean tolerance class per FT and its weight (m2 of cover over all sites), for the year rows
`rows` of the CscapeOutput `out` (combined intervention slot).
"""
function cscape_thermal(out, rows)
    cls = CSCAPE_SA_HEAT_CLASSES
    oa = out.out_array                                                # [years, sites, int, ft, class, values]
    size(oa, 5) == length(cls) || error("C~scape output has $(size(oa, 5)) heat classes, expected $(length(cls))")
    area = out.area[:, 3]
    c = oa[rows, :, 3, :, :, CScape.Layout.VAL_COVER] ./ 100 .* reshape(area, 1, :, 1, 1)   # m2 [T, site, ft, class]
    c = dropdims(sum(c; dims=2); dims=2)                              # [T, ft, class]
    W = dropdims(sum(c; dims=3); dims=3)
    M = [W[t, g] > 0 ? sum(cls .* c[t, g, :]) / W[t, g] : NaN for t in axes(c, 1), g in axes(c, 2)]
    return M, W
end

"""
    state_summary(r) -> Dict{String,Array{Float64}}

`r` is `run_row`'s result (`r.out` the CscapeOutput; `r.reef_k`, `r.reef_area` this run's k and area
after the multipliers).
"""
function state_summary(r)
    out = r.out
    rows = _cscape_sa_rows(out.years)
    cover = CscapeInterface.to_adria_format(out).cover[rows, :, :, :]                # [T, ft, sizes, reefs], relative to k x area
    taxa = ADRIAIndicators.relative_loc_taxa_cover(cover)             # [T, ft, reefs]
    M, W = cscape_thermal(out, rows)
    return Dict{String,Array{Float64}}(
        "cover" => Float64.(ADRIAIndicators.relative_cover(cover)),                    # [T, reefs]
        "class_cover" => class_covers(taxa, PF_CLASS_OF_GROUP),                        # [T, 3, reefs]
        "juveniles" => Float64.(ADRIAIndicators.relative_juveniles(cover, out.is_juvenile)),   # [T, reefs]
        "diversity" => Float64.(ADRIAIndicators.coral_diversity(taxa)),                # [T, reefs]
        "weights" => Float64.(r.reef_k .* r.reef_area),                                # [reefs]
        "thermal_mean" => M,                                                           # [T, ft]
        "thermal_weight" => W,                                                         # [T, ft]
    )
end
