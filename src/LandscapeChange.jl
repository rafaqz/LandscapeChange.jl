module LandscapeChange

using Reexport

@reexport using DynamicGrids

import Distributions
using DynamicGrids.Stencils
using DynamicGrids.Setfield
using DynamicGrids.StaticArrays
using DynamicGrids.DimensionalData
using DynamicGrids.DimensionalData.LookupArrays

using Rasters

export LandCoverCount, LandCoverPotential, LandCoverSuitability, LandCoverChange

export WhiteEngalinUljeeWeights, WhiteEngalinUljeeUpdate, BottomUp

export NeutralLUCModel, RandomConstraintMatch, GrowingClusters

export SuitabilityScore, MaxProbability, ChangeProbabilities, RawSuitabilities

export PreallocatedUnorderedList, NamedVector

export sim, fuzzy_scores, cover_change, cover_persistence, cover_fraction

export stripe_raster, namedvector_raster, compile_timeline, cross_validate_timeline

import DynamicGrids.Stencils: neighbors, unsafe_neighbors, stencil, unsafe_stencil,
    kernel, kernelproduct, offsets, positions, radius, distances, distance_zones,
    boundary, padding, source, dest, switch, padval

include("unordered_list.jl")
include("named_vector.jl")
include("models/rules.jl")
include("models/bottom_up.jl")
include("models/white_engalin.jl")
# include("models/historical.jl")
include("analysis/neutral.jl")
include("analysis/scores.jl")
include("analysis/statistics.jl")
include("input/transitions.jl")
include("input/generate.jl")

end
