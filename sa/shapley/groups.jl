#=
groups.jl -- grouping schemes of the Shapley inputs, read from groups/<model>.toml.

    include(joinpath(@__DIR__, "groups.jl"))
    G = load_groupings(joinpath(@__DIR__, "groups", "coralpop.toml"), names(R["X"]))
    G["process"].names, G["process"].members     # group names, X column indices per group

File format (schemes and groups keep their order):

    extends = "shared.toml"          # optional: its groups come first in schemes of the same name
    inert   = ["m_depth"]            # optional: inputs without effect in this model (noise floor in figures)

    [[scheme]]
    name = "process"
      [[scheme.group]]
      name  = "DHW"
      match = ["^dhw_"]              # regexes on the input names

Every input must fall in exactly one group of each scheme; anything else is an error that lists the
offending inputs. Groups that match no input (e.g. a shared group a model lacks) are dropped.
=#

using TOML

"One grouping scheme: group `names[g]` holds the X columns `members[g]`."
struct Grouping
    scheme::String
    names::Vector{String}
    members::Vector{Vector{Int}}
end

"Schemes of a groups file as `scheme => [(group, patterns)]`, in file order, `extends` merged in first."
function _read_schemes(path::AbstractString)::Vector{Pair{String,Vector{Tuple{String,Vector{Regex}}}}}
    d = TOML.parsefile(path)
    base = haskey(d, "extends") ? _read_schemes(joinpath(dirname(path), d["extends"])) : Pair{String,Vector{Tuple{String,Vector{Regex}}}}[]
    for sch in get(d, "scheme", [])
        groups = [(String(g["name"]), Regex.(g["match"])) for g in get(sch, "group", [])]
        i = findfirst(p -> p.first == sch["name"], base)
        i === nothing ? push!(base, sch["name"] => groups) : append!(base[i].second, groups)
    end
    return base
end

"""
    load_groupings(path, input_names) -> Dict{String,Grouping}

The grouping schemes of `path` resolved against the inputs `input_names` (the columns of `X`).
"""
function load_groupings(path::AbstractString, input_names::AbstractVector)::Dict{String,Grouping}
    inputs = string.(input_names)
    out = Dict{String,Grouping}()
    for (scheme, groups) in _read_schemes(path)
        hits = [[j for (j, n) in enumerate(inputs) if any(r -> occursin(r, n), pats)] for (_, pats) in groups]
        count_in = [count(h -> j in h, hits) for j in eachindex(inputs)]
        missing_ = inputs[count_in .== 0]
        double = inputs[count_in .> 1]
        isempty(missing_) && isempty(double) || error(
            "groups file $path, scheme '$scheme': " *
            (isempty(missing_) ? "" : "in no group: $(join(missing_, ", ")). ") *
            (isempty(double) ? "" : "in more than one group: $(join(double, ", ")).")
        )
        keep = .!isempty.(hits)
        out[scheme] = Grouping(scheme, first.(groups)[keep], hits[keep])
    end
    return out
end

"Inputs a groups file marks as `inert` (no effect in this model by construction: a noise floor), `extends` included."
function inert_inputs(path::AbstractString)::Vector{String}
    d = TOML.parsefile(path)
    base = haskey(d, "extends") ? inert_inputs(joinpath(dirname(path), d["extends"])) : String[]
    return unique([base; String.(get(d, "inert", String[]))])
end

"Scheme names of a groups file, in file order."
scheme_names(path::AbstractString) = first.(_read_schemes(path))

"Name of the group of input column `j` in grouping `G`."
group_of(G::Grouping, j::Integer) = G.names[findfirst(m -> j in m, G.members)]

"""
    resolve_group(groupings, "scheme:group") -> (Grouping, group index)
"""
function resolve_group(groupings::Dict{String,Grouping}, ref::AbstractString)
    parts = split(ref, ":"; limit=2)
    length(parts) == 2 || error("group reference '$ref' must be \"scheme:group\"")
    scheme, name = parts
    haskey(groupings, scheme) || error("unknown scheme '$scheme' in '$ref'; known: $(join(keys(groupings), ", "))")
    G = groupings[scheme]
    g = findfirst(==(name), G.names)
    g === nothing && error("no group '$name' in scheme '$scheme'; groups: $(join(G.names, ", "))")
    return G, g
end
