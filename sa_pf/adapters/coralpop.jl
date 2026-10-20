#=
adapters/coralpop.jl -- CoralPop's state summary (sa_outputs.jl) from one run's in-memory output.

Included by pf_common.jl after models/coralpop.jl. Every model adapter defines:

    PF_CLASS_OF_GROUP   harmonised class of each native group (sa_outputs.jl `native_group_class`)
    PF_THERMAL_GROUPS   names of the "thermal_mean" columns (empty: no heritable tolerance)
    state_summary(r)    the state summary of one run, from `run_row`'s result `r`

CoralPop: 3 morphologies (Branching = corymbose non-Acropora, Foliose = tabular + corymbose Acropora,
Other = small + large massives; CORALPOP_MORPH_ADRIA_GROUPS). The indicators come from
`coralpop_indicators` (ADRIAIndicators definitions): relative_juveniles with CoralPop's juvenile bins
(CORALPOP_JUVENILE_BINS, 0-10 cm), coral_diversity over the 3 morphologies. CoralPop has only a fixed
bleaching acclimation, no heritable tolerance, so no thermal outputs.
=#

const PF_CLASS_OF_GROUP = native_group_class(CORALPOP_MORPH_ADRIA_GROUPS)
const PF_THERMAL_GROUPS = String[]

"""
    state_summary(r) -> Dict{String,Array{Float64}}

`r` is `run_row`'s result (`r.out` the CoralPopOutput; `r.reef_k`, `r.reef_area` this run's k and area
after the multipliers).
"""
function state_summary(r)
    ind = coralpop_indicators(r.out)
    taxa = Float64.(ind.relative_loc_taxa_cover)                    # [T, morphologies, reefs]
    return Dict{String,Array{Float64}}(
        "cover" => Float64.(ind.relative_cover),                    # [T, reefs]
        "class_cover" => class_covers(taxa, PF_CLASS_OF_GROUP),     # [T, 3, reefs]
        "juveniles" => Float64.(ind.relative_juveniles),            # [T, reefs]
        "diversity" => Float64.(ind.coral_diversity),               # [T, reefs]
        "weights" => Float64.(r.reef_k .* r.reef_area),             # [reefs]
    )
end
