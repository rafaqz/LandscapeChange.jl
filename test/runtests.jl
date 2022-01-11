using LandscapeChange
using NeutralLandscapes
using DynamicGrids
using Rasters
using Plots
using Test
using StaticArrays

@testset "" begin
    proj = EPSG(4326)
    projdims = X(Projected(20:0.01:22; crs=proj)),
               Y(Projected(9:0.01:11; crs=proj))
    random_landscape = Raster(rand(MidpointDisplacement(), projdims); missingval=nothing)
    c = NamedVector((forest=1, cleared=2, settled=3))

    landscape = Rasters.classify(random_landscape, 
         0..0.5 => c.forest, 0.5..0.8 => c.cleared, 0.8..1.0 => c.settled;
         others=typemin(Int), missingval=typemin(Int)
    ) 
    # plot(plot(random_landscape), plot(landscape; c=:viridis))

    counts = broadcast(_ -> zero(c), landscape)

    lu_count_rule = LandCoverCount{:landscape,:counts}(Window(1), c)
    extent = Extent((; landscape, counts); tspan=1:10)
    sd = DynamicGrids.SimData(extent, lu_count_rule)
    sd = step!(sd)
    sd[:counts]

    suitability_params = (
        riverdist=Aux{:distance_to_rivers}(),
        coastdist=Aux{:distance_to_coast}(),
        elevation=Aux{:elevation}(),
        steepness=Aux{:steepness}(),
    )
    suitability_params = (;)

    suitability_functions = (
        forest=lc -> (forest=0.9, cleared=0.1, settled=0.0)[lc],
        cleared=lc -> (forest=0.1, cleared=0.9, settled=0.01)[lc],
        settlement=lc -> (forest=0.0, cleared=0.05, settled=0.95)[lc],
    )

    lu_suitability_rule = LandCoverSuitability{:landscape,:suitability}(;
        functions=suitability_functions, mode=ChangeProbabilities()
    )

    using BenchmarkTools
    f(xs) = map(x -> x^2, xs)
    @btime f(distances(Window(1)))


    suitability = broadcast(_ -> zero(c) .* 0.0, landscape)
    extent = Extent((; landscape, counts, suitability); tspan=1:10)
    sd = DynamicGrids.SimData(extent, lu_suitability_rule)
    sd = step!(sd)
    sd[:suitability]
end

@testset "NamedVector" begin
    using BenchmarkTools
    NamedVector{(:a,:b,:c),3,Int}(NamedTuple{(:a,:b,:c)}((1,2,3)))
    @btime NamedVector{(:a,:b,:c),3}((1.0,2,3))
    @btime NamedVector{(:a,:b,:c)}((1,2,3))
    # This doesn't completely compile away for some reason
    @btime NamedVector{(:a,:b,:c),3}(NamedTuple{(:a,:b,:c)}((1,2,3)))
    @code_warntype NamedVector(NamedTuple{(:a,:b,:c)}((1,2,3)))
    @btime map(x -> 2x, $c)
end
