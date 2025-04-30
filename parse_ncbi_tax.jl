#=
Parse the NCBI dump data to a more regular format, with useless clades removed
etc.
The input should be the latest from https://ftp.ncbi.nih.gov/pub/taxonomy/

It creates a list of child/parent pairs, with the following changes from the raw data:

* Removes all descendants from known catch-all trash clades like "environmental samples"

* For any pair (child, parent) where the parent is not a canonical rank (i.e. one of
  the standard seven ranks used to e.g. GTDB), set the parent to the nearest canonical
  ancestor.
  E.g. if a genus has a subfamily parent, set its parent to the family instead.
  For this purpose, the universal ancestor of all life ("LUCA", tax number 1) is canonical.
  
* Remove all nodes which are not descendant from LUCA
=#

function main(args)
    if length(args) != 3
        println(stderr, "Usage: julia parse_ncbi_tax.jl outfile names.dmp nodes.dmp")
        exit(1)
    end
    (outpath, namespath, nodespath) = args
    ispath(outpath) && error("Outpath exists: \"$(outpath)\"")
    for f in [namespath, nodespath]
        isfile(f) || error("No such file: \"$(f)\"")
    end
    open(outpath, "w") do out
        make_ncbi(out, namespath, nodespath)
    end
end

function make_ncbi(out::IO, names_path::AbstractString, nodes_path::AbstractString)
    nodes = open(parse_nodes, nodes_path)
    parent_dict = build_parent_dict(nodes)

    # We remove these nodes, because they don't refer to real clades.
    to_remove = [
        57727, # environmental samples
    ] |> Set

    # Remove the nodes in to_remove and all their descendants
    parent_dict = remove_descendants(parent_dict, to_remove)

    # Get names, but only for those we haven't filtered away
    keep_ids = Set([c.id for c in keys(parent_dict)])
    names = open(names_path) do io
        get_names(io, keep_ids)
    end

    # Rename some important names
    # This is called "Bacteria <bacteria>" for some reason.
    names[2] = (NameTypes.scientific_name, "Bacteria")

    # Remove non-canonical ranks like infraorder.
    parent_dict = make_parent_canonical(parent_dict, RANK_TO_CANONICAL_INDEX)

    # Remove clades not directly descendant from the universal common ancestor
    parent_dict = remove_holey_descendants(parent_dict)

    fields = [
        (child.id, child.rank, parent.id, names[child.id][2])
            for (child, parent) in parent_dict
    ]

    # Make the output file compress better by sorting according to name
    sort!(fields; by = last)
    println(out, "child_id\tchild_rank\tparent_id\tname")
    for i in fields
        println(out, join(i, '\t'))
    end
    return
end

# All ranks in NCBI file
module Ranks
    names = [
        "acellular root",
        "biotype",
        "cellular root",
        "clade",
        "class",
        "cohort",
        "domain",
        "family",
        "forma",
        "forma specialis",
        "genotype",
        "genus",
        "infraclass",
        "infraorder",
        "isolate",
        "kingdom",
        "morph",
        "no rank",
        "order",
        "parvorder",
        "pathogroup",
        "phylum",
        "realm",
        "section",
        "series",
        "serogroup",
        "serotype",
        "species",
        "species group",
        "species subgroup",
        "strain",
        "subclass",
        "subcohort",
        "subfamily",
        "subgenus",
        "subkingdom",
        "suborder",
        "subphylum",
        "subsection",
        "subspecies",
        "subtribe",
        "subvariety",
        "superclass",
        "superfamily",
        "superorder",
        "superphylum",
        "tribe",
        "varietas",
    ]
    syms = map(i -> replace(i, " " => "_"), names)
    @eval @enum Rank::UInt8 $(Symbol.(syms)...)
    const STR_TO_ENUM = Dict(n => i for (n, i) in zip(names, instances(Rank)))
end

# All clade name types in NCBI
module NameTypes
    # Sorted in order from worst to best
    names = [
        "in-part",
        "acronym",
        "type material",
        "blast name",
        "includes",
        "synonym",
        "authority",
        "genbank acronym",
        "equivalent name",
        "genbank common name",
        "common name",
        "scientific name",
    ]
    syms = map(i -> replace(i, " " => "_", "-" => "_"), names)
    @eval @enum NameType::UInt8 $(Symbol.(syms)...)
    const STR_TO_ENUM = Dict(n => i for (n, i) in zip(names, instances(NameType)))
end

using .Ranks: Ranks, Rank

# These are the ranks we care to keep
const RANK_TO_CANONICAL_INDEX = Dict(
    Ranks.species => 1,
    Ranks.genus => 2,
    Ranks.family => 3,
    Ranks.order => 4,
    Ranks.class => 5,
    Ranks.phylum => 6,
    Ranks.domain => 7,
)

using .NameTypes: NameTypes, NameType

Base.tryparse(::Type{Rank}, s::AbstractString) = get(Ranks.STR_TO_ENUM, s, nothing)

function Base.parse(::Type{Rank}, s::AbstractString)
    r = tryparse(Rank, s)
    r === nothing && error(lazy"Could not parse as Rank: \"$(s)\"")
    r
end

Base.tryparse(::Type{NameType}, s::AbstractString) = get(NameTypes.STR_TO_ENUM, s, nothing)

struct Node
    id::Int
    parent::Int
    rank::Rank
end

is_top(x::Node) = x.id == 1

"Remove all descendants of the nodes in `to_remove`"
function remove_descendants(
        parents::Dict{Node, Node},
        to_remove::Set{Int},
    )::Dict{Node, Node}
    result = empty(parents)
    for (child, parent) in parents
        should_add = true
        p = child
        while !is_top(p)
            if p.id ∈ to_remove
                should_add = false
                break
            end
            p = parents[p]
        end
        if should_add
            result[child] = parent
        end
    end
    @assert length(result) <= length(parents)
    return result
end

"""
    parse_nodes(io::IO)::Vector{Node}

Parse the NCBI nodes dump file, and return a vector of `Node`s.
"""
function parse_nodes(io::IO)::Vector{Node}
    result = Node[]
    for line in eachline(io)
        line = chopsuffix(rstrip(line), "\t|")
        isempty(line) && continue
        (id_str, parent_id_str, rankname) = split(line, "\t|\t")
        node = Node(parse(Int, id_str), parse(Int, parent_id_str), parse(Rank, rankname))
        push!(result, node)
    end
    return result
end

function build_parent_dict(nodes::Vector{Node})::Dict{Node, Node}
    parent_of = Dict{Node, Node}() # only top node has no parent
    node_by_id = Dict{Int, Node}()
    for node in nodes
        if get!(node_by_id, node.id, node) !== node
            error("Duplicate id")
        end
    end
    for node in nodes
        is_top(node) && continue
        parent_of[node] = node_by_id[node.parent]
    end
    return parent_of
end

"""
Change the parent to be the closest ancestor of a canonical clade.

Lots of the entries in the tax database are partially unknown, e.g.
it may be a genus which is the child of an unlabeled clade, which is
the child of a known family.
In this function, we remove all intermediate non-canonical clades, such that
every child's parent is either the top node, or a canonical clade.

Further, some canonical children have parents of the wrong rank, e.g.
a species whose direct parent is a family (skipping the genus rank).
We remove these.
"""
function make_parent_canonical(
    parent_of::Dict{Node, Node},
    rank_to_index::Dict{Rank, <:Integer}
)
    new_parents = empty(parent_of)
    for (child, parent) in parent_of
        new_parent = parent
        parent_index = get(rank_to_index, new_parent.rank, nothing)
        while !is_top(new_parent) && isnothing(parent_index)
            new_parent = parent_of[new_parent]
            parent_index = get(rank_to_index, new_parent.rank, nothing)
        end
        # Tax id 1 (the top) is "all life"
        parent_index = if is_top(new_parent)
            8
        else
            parent_index
        end

        @assert is_top(new_parent) || child.rank != new_parent.rank

        # A lot of clades have incomplete rankings - e.g. a species
        # being the child of a "no rank" being the child of an order.
        # We skip these here.
        child_index = get(rank_to_index, child.rank, nothing)
        if child_index === nothing || child_index == parent_index - 1
            new_parents[child] = new_parent
        end
    end
    return new_parents
end


"""Remove any descendants that can't be traced directly to the top.
We do this by traversing from the top down, and only include reached nodes"""
function remove_holey_descendants(parent_of::Dict{Node, Node})
    new_parents = empty(parent_of)
    children_of = Dict{Node, Vector{Node}}()
    for (child, parent) in parent_of
        push!(get!(valtype(children_of), children_of, parent), child)
    end
    top = first(Iterators.filter(is_top, values(parent_of)))
    stack = Node[top]
    while !isempty(stack)
        parent = pop!(stack)
        for child in get(children_of, parent, Node[])
            new_parents[child] = parent
            push!(stack, child)
        end
    end
    return new_parents
end


"""
    get_names(io::IO, tax_ids::Union{Nothing, Set{Int}})::Dict{Int, Tuple{NameType, String}}

Parse the NCBI names dump fil, and return a dictionary mapping tax ids to names,
only for the taxids in `tax_ids`. If `tax_ids` is `nothing`, all tax ids are included.
"""
function get_names(
        io::IO,
        tax_ids::Union{Nothing, Set{Int}},
    )::Dict{Int, Tuple{NameType, String}}
    # Store all names for a given tax id (tax ids may have multiple names)
    # E.g. "human", "Homo sapiens" etc.
    names_of_tax = Dict{Int, Vector{Tuple{NameType, String}}}()
    for line in eachline(io)
        line = chopsuffix(rstrip(line), "\t|")
        # The NCBI file has two kinds of names: Unique names, and non-unique names.
        # For the latter, different clades may have the same name.
        # The unique name is empty if the name itself is already unique.
        (tax_str, non_unique_name, unique_name, class_str) = split(line, "\t|\t")
        name = isempty(unique_name) ? non_unique_name : unique_name
        tax_id = parse(Int, tax_str)
        # Only include tax ids the user asked for
        !isnothing(tax_ids) && tax_id ∉ tax_ids && continue
        # The name type is e.g. "common name" or "scientific name".
        # The NCBI file contains multiple rows with different names for the same clade,
        # so if there are multiple, we pick the kind of name we like the best.
        class = tryparse(NameType, class_str)
        isnothing(class) && error(
            lazy"Error: Could not parse name type $class_str for tax"
        )
        push!(
            get!(valtype(names_of_tax), names_of_tax, parse(Int, tax_str)),
            (class, String(name)),
        )
    end
    # Pick the kind of name we like the best, where multiple names are present.
    # E.g. we prefer "Opacifrons coxata" over "Opacifrons coxata (Stenhammar, 1854)".
    result = Dict(
        k => argmax(((class, name),) -> Integer(class), v) for (k, v) in names_of_tax
    )
    # Assert names are unique
    names = Set{String}()
    for (_, (_, name)) in result
        if UInt8('\t') ∈ codeunits(name)
            error(lazy"Cannot have tab in name, but it is present in $(name)")
        end
        L = length(names)
        push!(names, name)
        length(names) == L && error(
            lazy"Error: Name \"$(name)\" is not unique"
        )
    end
    return result
end

main(ARGS)
