"""
    LandCoverChangeMode

Methods for randomly selecting a new land-use category.
"""
abstract type SuitabilityScore end

"""
    MaxProbability <: LandCoverChangeMode

Select from the highest probability land-use change or
the current state.
"""
struct MaxProbability <: SuitabilityScore end

# Return the category number (1..N)
score(::MaxProbability, suitabilities) = findmax(values(suitabilities))

"""
    ChangeProbabilities <: LandCoverChangeRandimisation

Choose the new land-use from the probabilities of shifting to all
possible land uses.
"""
struct ChangeProbabilities <: SuitabilityScore end

function score(::ChangeProbabilities, suitabilities::NamedTuple{K}, lc) where K
    probabilities = map(s -> s / sum(suitabilities), suitabilities)
    cumulative_probabilities = reduce(Base.tail(probabilities); init=(first(probabilities),)) do acc, x
       (acc..., x + last(acc))
    end
    return NamedVector{K}(cumulative_probabilities)
end

"""
    RawSuitabilities <: LandCoverChangeRandimisation

Choose the new land-use from the probabilities of shifting to all
possible land uses.
"""
struct RawSuitabilities <: SuitabilityScore end

score(::RawSuitabilities, suitabilities) = suitabilities

"""
    LandCoverCount{R,W} <: NeighborhoodRule
    LandCoverCount{R,W}(neighborhood, categories)

Simple rule that counts category values in the neighborhood.

- `R` must be a single `Symbol` representing the grid that holds
    `Int` land-use values
- `W` must be a single `Symbol` representing the grid to write counts to,
    with the `eltype` of `NameTuple` matching the `categories` `NameTuple`.
"""
struct LandCoverCount{R,W,N,C} <: NeighborhoodRule{R,W}
    neighborhood::N
    categories::C
end
LandCoverCount(; neighborhood=Moore(1), categories) = LandCoverCount(neighborhood, categories)
categories(rule::LandCoverCount) = rule.categories

function DynamicGrids.applyrule(data, rule::LandCoverCount, lc, I)
    _count_categories(neighborhood(rule), categories(rule))
end

function _count_categories(itr, categories)
    # The simplest method is fastest for low numbers of categories
    map(categories) do e
        count(x -> x === e, itr)
    end
end

struct LandCoverSuitability{R,W,F,P,M} <: CellRule{R,W}
    functions::F
    params::P
end
function LandCoverSuitability{R,W}(; functions, params=nothing, mode=RawSuitabilities()) where {R<:Tuple{},W}
    LandCoverSuitability{R,W}(functions, params, mode)
end

DynamicGrids.applyrule(data, rule::LandCoverSuitability, _, I) = _suitability(data, rule, I)

function _suitability(data, rule, I)
    cellparams = map(rule.params) do p
        get(data, p, I)
    end
    map(f -> f(cellparams), rule.functions)
end

"""
    LandCoverPotential{R,W} <: NeighborhoodRule
    LandCoverPotential{R,W}(suitabilities, categories,

Calculate land-use potential from the land-use categories
of the surrounding neighborhoods, and suitability parameters.
"""
struct LandCoverPotential{R<:Tuple{<:Any,<:Any},W,N,NM,SP,SF,C,I,RM} <: NeighborhoodRule{R,W}
    categories::C
    neighborhood::N
    suitability_parameters::SP
    suitability_functions::SF
    neighbor_mode::NM
    result_mode::RM
end
function LandCoverPotential(;
    categories,
    neighborhood=Moore(1),
    suitability_paremeters=nothing,
    suitability_functions=nothing,
    neighbor_mode=NeighborCount(),
    result_mode=ChangeProbabilities(),
)
    LandCoverPotential(
        categories, neighborhood, suitability_paremeters,
        suitability_functions, neighbor_mode, result_mode,
    )
end

function DynamicGrids.applyrule(data, rule::LandCoverPotential, lc, I)
    # Get auxiliary suitability data
    cell_suitabilities = _suitability(data, rule, I)
    # Count all the categories
    potential = _category_potential(rule.neigbor_mode, )
    # Compine the neighborhood-based potential
    # with the land-use suitability and pressure
    cell_potential = map(*, suitabilities, potential)

    return score(rule.mode, cell_potential)
end

abstract type NeighborhoodInfluence end
struct NeighborCount <: NeighborhoodInfluence end
struct NeighborDistance{F} <: NeighborhoodInfluence
    f::F
end
NeighborDistance() = NeighborDistance()

# Utils
function _category_potential(nd::NeighborDistance, rule)
    # This should compile away - `distances` are known at compile time
    dist_weights = map(nd.f, distances(rule))

    return map(categories) do c
        nsum = reduce(zip(neighbors(rule), dist_weights); init=zero(first(dist_weights))) do acc, (n, d)
            n === c ? acc + d : acc
        end
        nsum / length(sum(distances(rule)))
    end
end

function _category_potential(::NeighborCount, rule)
    map(_countcategories(neighborhood(rule), categories(rule))) do x
        x / length(neighborhood(rule)) 
    end
end

struct LandCoverChange{R<:Tuple{<:Any,<:Any},W,S,P,X} <: NeighborhoodRule{R,W}
    suitability::S
    pressure::P
    weights::X
end

function DynamicGrids.modifyrule(data, rule::LandCoverChange{<:Tuple{S,P}}) where {S,P}
    _get_thresholds(data[P], rule)
end
function _get_thresholds(A, rule::LandCoverChange{R}) where R
    # steps = map(x -> ?, rule.suitability)
    thresholds = map(x -> false, rule.suitability)
    nbelow = map(x -> 0, rule.suitability)
    while true
        counts = mapreduce((acc, xs) -> map(+, acc, xs), A) do val, cat
            val > thresholds[icategory] ? setindex(nbelow, 1, cat) : nbelow
        end
        all(map(counts, targets) do c, t
            abs(c - t) > rule.tol
        end) && break
        thresholds = map(counts, targets, thresholds, steps) do th, s
            c > t ? th + abs(s) : th - abs(s)
        end
    end

    rule = @set rule.thresholds = thresholds
    return rule
end

function DynamicGrids.applyrule(
    data, rule::LandCoverChange{R}, (state, (potential, category)), I
) where R
    # TODO add noise
    potential > rule.thresholds[category] ? category : state
end

# influences = (
#     forest = (
#         forest = 1.0,
#         cleared = 0.1,
#         urban = 0.01,
#     ),
#     cleared = (
#         forest = 0.1,
#         cleared = 1.0,
#         urban = 0.2
#     ),
#     urban = (
#         forest =
#         cleared =
#         urban = 1.0
#     ),
# )

# LandCoverSuitability(
#     params=(slope=Aux{:slope}(), river=Aux{:river}()),
#     suitability = (
#         forest=(_ -> 1.0)
#         cleared=(p ->
#         urban=(p ->
#     ),
# )



# function checkMacroDemand(thisYearsMacroDemand):
#     tot = 0
#     for luNr in thisYearsMacroDemand.keys()
#         tot += thisYearsMacroDemand[luNr]
#     end
#     return tot
# end

#function applyMacroDemand(macroDemand, year, luNrsDyn, luNrsStat, potArrs, arrLU)
#function applyMacroDemand(macroDemand, year, neighborhood, luNrsDyn, luNrsStat, potArrs, arrLU)
#    """ Make the new land use map by changing the required cells for current time step,
#    for the dynamic land uses.
#    """

#    n = 8 # neighborhood size

#    thisYearsMacroDemand = macroDemand[year]
#    #lastYearsMacroDemand = macroDemand[year-1]

#    # Make a new arr for the new land use array, where all dynamic LU are replaced by number 1.
#    newArrLU = Utils.findAndReplace(arrLU.copy(), find=luNrsDyn, replaceWith=1, within=[(n,n),(-n,-n)])

#    # Cut away the surrounding "frame" of neighborhood cells which are not supposed to be calculated.
#    for nr in potArrs.keys()
#        potArrs[nr] = potArrs[nr][n:-n, n:-n]
#        # print "max", potArrs[nr].max()
#    end
#    # This raster keeps track of which land use has been changed already so that it wont change again.
#    arrLUChange = numpy.zeros(arrLU.shape, dtype=numpy.integer) # I would like int0/unit0 here instead

#    # Iterate
#    while checkMacroDemand(thisYearsMacroDemand) > 0
#        mx = 0 # int Note! Maybe this should be initialized as None?
#        luNr = None # int
#        tempMax = None # int
#        # Get max value of all potential arrays and find out which array (LuNr) holds it
#        for tLu in potArrs.keys()
#            # If the land use demand is satisfied for one land use,
#            # then don't try to find the potential for this LU. When
#            # all LU are satisfied, the loop with stop automatically (see while-statement)
#            if thisYearsMacroDemand[tLu]<=0:
#                continue
#            tempMax = potArrs[tLu].max()
#            if tempMax > mx
#                mx = tempMax # save current max val
#                luNr = tLu # save lu (which has the highest potential)
#            end
#        if luNr == None
#            # print "Breaking when macroDemandStatus is: %i" %(checkMacroDemand(thisYearsMacroDemand))
#            break
#        end

#        # Find out the xy-location of the max value
#        #print "lu", lu, "tempMax", tempMax, "mx", mx, "luNr", luNr
#        potArr = potArrs[luNr] # get potArr for this land use
#        sortedArr = potArr.argsort() # sort it according to potential

#        rowIndex = 0
#        maxValList = []
#        for column in sortedArr # (column is an entire column (1D-array) )
#            # Get the column index which has the max value, for row nr = rowIndex
#            mxColIndex = column[-1]

#            # get the max value, of the argsorted array, for comparison with max
#            val = potArr[rowIndex, mxColIndex]

#            if val == mx: # if there is more than one max-value... choose randomly
#                maxValList.append((rowIndex, mxColIndex))
#            end
#            rowIndex+=1
#        end

#        # One (or more) locations of the max potential value found,
#        # inserted into the list maxValList. In case many values - pick a random.
#        random.shuffle(maxValList)
#        row, col = maxValList[0]
#        luRow, luCol = row+n, col+n # Add the neighborhood size to the row and col.

#        #print "maxVal=", potArr[row, col]

#        potArr[row, col] = -999 # make sure it's not again selected for land use change.

#        # Don't allow LU change if LU has already been assigned once.
#        if arrLUChange[luRow, luCol] == 1
#            continue
#        end

#        nrExistingLU = arrLU[luRow, luCol]
#        if nrExistingLU not in luNrsStat
#            # Allow land use change
#            newArrLU[luRow, luCol] = luNr
#            arrLUChange[luRow, luCol] = 1 # register as "changed"

#            # If replacing dynamic land use (which is not its own), increment the replaced LU demand by one.
#            if nrExistingLU in luNrsDyn and nrExistingLU!=luNr:
#                thisYearsMacroDemand[nrExistingLU] += 1

#            # Always decrease decrement the value when satisfying LU demand.
#            thisYearsMacroDemand[luNr] -= 1

#            # ...and decrease loop counters value by 1.
#            # If all demand satisfied for this land use (empty) - take away lu potential arr
#            # This means, it won't look for LU change potential in this array anymore.
#            #if thisYearsMacroDemand[luNr] <= 0:
#                #print "del potArrs[%i]" % luNr
#                #del potArrs[luNr]
#        else
#            # Continue, without subtracting from the demand counter - since the location of
#            # existing land use was not allowed to change
#            continue
#        end
#        #print macroDemandStatus
#    # Areas not assigned any dynamic land use gets value nr 1 - i.e. same as vacant lu.
#    return newArrLU
#end

#function makeMacroDemandDict(macroDemandNodes):
#    """ output should look like this:
#        dict = {1998 : {3 : 187, 4 : 459}, 1999 : {3 : 193, 4 : 498}}"""
#    pass
#endVery interesting.
