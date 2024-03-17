using LandscapeChange
# using NeutralLandscapes
# using DynamicGrids
# using Stencils
# using Rasters
# using Plots
using Test
# using StaticArrays

include("apply_transitions.jl")

# @testset "NamedVector" begin
#     c = NamedVector((forest=1, cleared=2, settled=3))
#     using BenchmarkTools
#     NamedVector{(:a,:b,:c),3,Int}((1,2,3))
#     @btime NamedVector{(:a,:b,:c),3}((1.0,2,3))
#     @btime NamedVector{(:a,:b,:c)}((1,2,3))
#     # This doesn't completely compile away for some reason
#     @btime NamedVector{(:a,:b,:c),3}((1,2,3))
#     @time NamedVector(NamedTuple{(:a,:b,:c)}((1,2,3)))
#     @btime map(x -> 2x, $c)
# end

# function shrink(A, category;
#     replacements=Tuple(filter(!=(category), sort(union(A)))),
#     neighborhood=Moore{1}(),
# )
#     Neighborhoods.broadcast_neighborhood(neighborhood, A) do h, x
#         if (x == category) && any(map(!=(category), h)) 
#             replacement_counts = map(replacements) do r
#                 count(==(r), h)
#             end
#             _, i = findmax(replacement_counts)
#             replacements[i]
#         else
#             x
#         end
#     end
# end
# function grow(A, category; neighborhood=Moore{1}())
#     Neighborhoods.broadcast_neighborhood(neighborhood, A) do h, x
#         any(map(==(category), h)) ? category : x
#     end
# end


# # Setup
# function gen_landscape(dims, name)
#     A = rand(NearestNeighborCluster(0.005), projdims)
#     return rebuild(round.(Int, NeutralLandscapes.classify(A, ones(3))); name)
# end

# proj = EPSG(4326)
# projdims = X(Projected(20:0.01:22; crs=proj)),
#            Y(Projected(9:0.01:11; crs=proj))
# init = gen_landscape(projdims, :init)
# target = rebuild(grow(init, 2; neighborhood=Moore{4}()); name=:target)
# target = rebuild(shrink(init, 2; neighborhood=Moore{3}()); name=:target)
# gen_landscape(projdims, :target)
# plot(plot.((init, target); c=:viridis)...)

# @testset "neutral models" begin
#     cats = sort(union(init))
#     LandscapeChange.adjust(init, target, cats)
#     @test length(LandscapeChange.category_list(init, 1)) == count(==(1), init)
#     @test "Category coverage matches after `RandomConstraintMatch`" begin
#         rcm = sim(RandomConstraintMatch(), init, target)
#         @test LandscapeChange.adjust(rcm, target, cats) == zero(cats)
#         change = rebuild(rcm .- init, name=:change)
#         error = rebuild(rcm .- target; name=:error)
#         @time fuzzy = fuzzy_scores(rcm, target; neighborhood=Moore(2))
#         plot(plot.((rcm, target))...)
#         plot(plot.((init, target, rcm); axis=nothing, c=:viridis)..., 
#              plot.((change, error, fuzzy); axis=nothing, c=:magma)...
#         )
#     end
#     @test "Category coverage matches after `GrowingClusters`" begin
#         LandscapeChange._weight(2, Moore{1}())
#         gc = sim(GrowingClusters(), init, target)
#         @test LandscapeChange.adjust(gc, target, cats) == zero(cats)
#         change = rebuild(gc .- init; name=:change)
#         error = rebuild(gc .- target; name=:error)
#         fuzzy_scores(gc, target);
#         plot(plot.((init, target, gc); axis=nothing, c=:viridis)..., 
#              plot.((change, error, fuzzy); axis=nothing, c=:magma)...)
#     end
# end

# @testset "" begin
#     counts = broadcast(_ -> zero(c), landscape)
#     lu_count_rule = LandCoverCount{:landscape,:counts}(Window(1), c)
#     extent = Extent((; landscape, counts); tspan=1:10)
#     sd = DynamicGrids.SimData(extent, lu_count_rule)
#     sd = step!(sd)

#     suitability_params = (
#         riverdist=Aux{:distance_to_rivers}(),
#         coastdist=Aux{:distance_to_coast}(),
#         elevation=Aux{:elevation}(),
#         steepness=Aux{:steepness}(),
#     )
#     suitability_params = (;)

#     suitability_functions = (
#         forest=lc -> (forest=0.9, cleared=0.1, settled=0.0)[lc],
#         cleared=lc -> (forest=0.1, cleared=0.9, settled=0.01)[lc],
#         settlement=lc -> (forest=0.0, cleared=0.05, settled=0.95)[lc],
#     )

#     lu_suitability_rule = LandCoverSuitability{:landscape,:suitability}(;
#         functions=suitability_functions, mode=ChangeProbabilities()
#     )

#     using BenchmarkTools
#     f(xs) = map(x -> x^2, xs)
#     @btime f(distances(Window(1)))


#     suitability = broadcast(_ -> zero(c) .* 0.0, landscape)
#     extent = Extent((; landscape, counts, suitability); tspan=1:10)
#     sd = DynamicGrids.SimData(extent, lu_suitability_rule)
#     sd = step!(sd)
#     sd[:suitability]
# end

# using CUDA
# A = LandscapeChange.NamedTuple{(:prob,:state,:inds)}.((zip(rand(Float32, 4000000), rand(1:3, 4000000), 1:4000000)))
# c = CuArray(A)
# @time partialsort!(A, 1:2000000; by=last);
# @time partialsort!(c, 1:2000000; by=last);
# @time partialsort!(c, 1:2000000; by=f);

# f(x) = x.state == 1 ? zero(x.prob) : x.prob



