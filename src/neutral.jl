
# See: Neutral models of land-use change as benchmarks in the assessment of model performance
# Alex Hagen-Zanker and Gilles Lajoie 

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
struct RandomConstraintMatch <: NeutralLUCModel end
struct GrowingClusters  <: NeutralLUCModel end


function sim(::RandomConstraintMatch, init, target, categories)
    result = copy(init)
    catlists = map(c -> category_list(init, c), categories)
    celllist = PreallocatedUnorderedList(CartesianIndex{2}; maxlength=length(init))
    adj = adjust(init, target, categories)
    for cat in categories
        adj[cat] > 0 || continue
        catlist = catlists[cat]
        for sample in 1:adj[cat]
            read_random!(celllist, catlist, cat)
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

function read_random!(celllist, catlist, cat)
    i = rand(1:length(catlist))
    push!(celllist, catlist[i])
    deleteat!(catlist, i)
end

function write_random!(result, celllist, cat)
    i = rand(1:length(celllist))
    result[celllist[i]] = cat
    deleteat!(celllist, i)
end

function sim(::GrowingClusters, init, target, categories; neighborhood=VonNeumann{1}())
    result = copy(init)
    ax = unpad_axes(init, neighborhood)
    inds = CartesianIndices(ax)
    edgetype = NamedTuple{(:cell, :cat, :hoodcat),Tuple{CartesianIndex{2},Int,Int}}
    edgelist = UnorderedList(edgetype, length(neighborhood) * length(inds))
    adj = [map(cat -> adjust(view(init, ax...), view(target, ax...), cat), categories)...]
    while !all(adj .== 0)
        for cell in inds
            hood = unsafe_updatewindow(neighborhood, result, cell)
            cat = result[cell] 
            add_hood_edges!(edgelist, hood, adj, cell, cat)
        end
        if isempty(edgelist)
            for cat in categories 
                if adj[cat] < 0
                    while adj[randomcat] <= 0 
                        cell = rand(inds)
                        randomcat = result[cell]
                    end
                    result[cell] = cat
                    adj[randomcat] -= 1
                    adj[cat] += 1
                end
            end 
        end
        while !isempty(edgelist)
            i = rand(1:length(edgelist))
            edge = edgelist[i]
            deleteat!(edgelist, i)
            write_if_matches!(result, adj, edge)
        end 
    end
    return result
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

function add_hood_edges!(edgelist, hood, adj, cell, cat)
    for hoodcat in hood
        if (adj[cat] > 0) && (adj[hoodcat] < 0)
            push!(edgelist, (; cell, cat, hoodcat))
        end 
    end
end
