
function category_list(A, cat)
    list = PreallocatedUnorderedList(CartesianIndex{2}; maxlength=length(A))
    for I in CartesianIndices(A)
        if A[I] == cat
            push!(list, I)
        end
    end
    return list
end

adjust(A, B, categories) = map(cat -> adjust(A, B, cat), categories)
adjust(A, B, cat::Int) = count(==(cat), A) - count(==(cat), B)

abstract type NeutralLUCModel end

"""
    RandomConstraintMatch <: NeutralLUCModel

Fill random locations in an `init` array to match
the number of pixels in each category in the `target`.

Adapted from: "Neutral models of land-use change as 
benchmarks in the assessment of model performance"
Alex Hagen-Zanker and Gilles Lajoie (2007)

This model can not be run on GPU arrays.

# Example
```julia
result = sim(RandomConstraintMatch(), init, target)
```
"""
struct RandomConstraintMatch <: NeutralLUCModel end

"""
    GrowingClusters <: NeutralLUCModel

Grow locations in an `init` array to match
the number of pixels in each category in the `target`.

Random seeds are created when there are no locations to grow.

This model can not be run on GPU arrays.

Adapted from: "Neutral models of land-use change as 
benchmarks in the assessment of model performance"
Alex Hagen-Zanker and Gilles Lajoie (2007)

# Example
```julia
result = sim(GrowingClusters(), init, target)
```
"""
struct GrowingClusters  <: NeutralLUCModel end


function sim(x::RandomConstraintMatch, init::AbstractDimArray, target::AbstractDimArray; kw...)
    newdata = sim(x, parent(init), parent(target); kw...) 
    return rebuild(init; data=newdata, name=:result)
end
function sim(::RandomConstraintMatch, init::AbstractArray, target::AbstractArray; 
    categories::AbstractVector=sort(union(init))
)
    result = copy(init)
    catlists = map(c -> category_list(init, c), categories)
    celllist = PreallocatedUnorderedList(CartesianIndex{2}; maxlength=length(init))
    adj = adjust(init, target, categories)
    for cat in categories
        adj[cat] > 0 || continue
        catlist = catlists[cat]
        for sample in 1:adj[cat]
            i = rand(1:length(catlist))
            push!(celllist, catlist[i])
            deleteat!(catlist, i)
        end
    end
    for cat in categories
        adj[cat] < 0 || continue
        for sample in 1:-adj[cat]
            write_random!(result, celllist, cat)
        end
    end
    return result
end


function write_random!(result, celllist, cat)
    i = rand(1:length(celllist))
    result[celllist[i]] = cat
    deleteat!(celllist, i)
end

function sim(x::GrowingClusters, init::AbstractDimArray, target::AbstractDimArray; kw...)
    newdata = sim(x, parent(init), parent(target); kw...) 
    return rebuild(init; data=newdata, name=:result)
end
function sim(::GrowingClusters, init::AbstractArray, target::AbstractArray; 
    categories::AbstractVector=sort(union(init)),
    neighborhood=VonNeumann{1}(),
)
    result = Neighborhoods.pad_array_out(init, neighborhood; padval=missing)
    ax = Neighborhoods.inner_axes(init, neighborhood)
    inds = CartesianIndices(ax)
    edgetype = NamedTuple{(:cell, :cat, :hoodcat),Tuple{CartesianIndex{2},Int,Int}}
    edgelist = PreallocatedUnorderedList(edgetype; maxlength=length(neighborhood) * length(inds))
    adj = map(c -> adjust(init, target, c), categories)
    while any(!=(0), adj)
        for cell in inds
            hood = unsafe_updatewindow(neighborhood, result, cell)
            cat = result[cell]
            # add_hood_edges!(edgelist, hood, adj, cell, cat)
            for hoodcat in hood
                if (adj[cat] > 0) && (adj[hoodcat] < 0)
                    edge = (; cell, cat, hoodcat)
                    push!(edgelist, edge)
                end
            end
        end
        if isempty(edgelist)
            # Randomly seed
            for cat in categories
                if adj[cat] < 0
                    local cell, randomcat
                    while true 
                        cell = rand(inds)
                        randomcat = result[cell]
                        adj[randomcat] > 0 && break
                    end
                    result[cell] = cat
                    adj[randomcat] -= 1
                    adj[cat] += 1
                end
            end
        end
        @show adj
        while !isempty(edgelist)
            i = rand(1:length(edgelist))
            edge = edgelist[i]
            deleteat!(edgelist, i)
            # write_if_matches!(result, adj, edge)
            stillover = adj[edge.cat] > 0
            stillunder = adj[edge.hoodcat] < 0
            stillneed = stillunder && stillunder
            stillsame = result[edge.cell] == edge.cat
            if stillneed && stillsame
                adj[edge.cat] -= 1
                adj[edge.hoodcat] += 1
                result[edge.cell] = edge.hoodcat
            end
        end
        @show adj
        println()
    end
    return view(result, Neighborhoods.inner_axes(result, neighborhood)...)
end

# These are mostly separated as function barriers for type stability
function write_if_matches!(result, adj, edge)
    stillover = adj[edge.cat] > 0
    stillunder = adj[edge.hoodcat] < 0
    stillsame = result[edge.cell] == edge.cat
    if stillover && stillunder && stillsame
        adj[edge.cat] -= 1
        adj[edge.hoodcat] += 1
        result[edge.cell] = edge.hoodcat
    end
end

# function add_hood_edges!(edgelist, hood, adj, cell, cat)
#     for hoodcat in hood
#         if (adj[cat] > 0) && (adj[hoodcat] < 0)
#             push!(edgelist, (; cell, cat, hoodcat))
#         end
#     end
# end
