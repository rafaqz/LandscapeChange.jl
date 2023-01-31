module LandscapeChange

using Reexport

@reexport using DynamicGrids

using DynamicGrids.Neighborhoods
using DynamicGrids.Setfield
using DynamicGrids.StaticArrays
using DynamicGrids.DimensionalData

export LandCoverCount, LandCoverPotential, LandCoverSuitability, LandCoverChange

export WhiteEngalinUljeeWeights, WhiteEngalinUljeeUpdate

export NeutralLUCModel, RandomConstraintMatch, GrowingClusters

export SuitabilityScore, MaxProbability, ChangeProbabilities, RawSuitabilities

export PreallocatedUnorderedList, NamedVector

export sim, fuzzy_scores, category_change, category_persistance

import DynamicGrids.Neighborhoods: neighbors, unsafe_neighbors, neighborhood,
    kernel, kernelproduct, offsets, positions, radius, distances, distance_zones,
    setneighbors, update_neighborhood, unsafe_update_neighborhood,
    boundary, padding, source, dest, switch, padval

include("unordered_list.jl")
include("named_vector.jl")
include("rules.jl")
include("neutral.jl")
include("scores.jl")
include("statistics.jl")
include("white_engalin.jl")

end
