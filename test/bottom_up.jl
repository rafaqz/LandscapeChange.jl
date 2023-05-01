using LandscapeChange, NeutralLandscapes, DimensionalData, Plots
using CUDAKernels, CUDA
using Revise, GLMakie
using Distributions


includet("/home/raf/.julia/dev/LandscapeChange/src/makie.jl")

# Define the states we are working with
states = NamedVector(a=1, b=2, c=3)

s = 500
# Define suitability maps
suitability_c = rand(PerlinNoise((3, 3)), X(s), Y(s); name="suitability_a") ./ 5 .+ 0.3
suitability_b = rand(PerlinNoise((4, 4)), X(s), Y(s); name="suitability_b") ./ 5 .+ 0.3
suitability_a = rand(PerlinNoise((4, 3)), X(s), Y(s); name="suitability_c") ./ 5 .+ 0.3
# plot(
#     heatmap(suitability_a; aspect_ratio=1), 
#     heatmap(suitability_b; aspect_ratio=1),
#     heatmap(suitability_c; aspect_ratio=1),
# )
suitability = map(suitability_a, suitability_b, suitability_c) do a, b, c
    NamedVector(; a, b, c)
end

# init_state = UInt8.((suitability_b .> (maximum(suitability_b) - 0.05)) .+ 1)
init_state = rand(1:3, size(suitability_a)) 

# Set an initial state where b has started in
# its absolute most suitable areas
# Plots.heatmap(init_state)
        
# Calculate change rates to drive the simulation
growth = vcat(0:0.5:10, 10:-0.5:0) # 10x growth over 10 timesteps
growth = 10:-0.2:0 # 10x growth over 10 timesteps
start_b = count(==(states.b), init_state)
start_c = count(==(states.c), init_state)
b_counts = DimArray(round.(Int, growth .* start_b), Ti)
c_counts = DimArray(round.(Int, growth .* start_c), Ti)
a_counts = length(suitability_a) .- b_counts .- c_counts
counts = DimArray(map(a_counts, b_counts, c_counts) do a, b, c
    NamedVector(; a, b, c)
end, Ti)
fixed = false # or a raster mask of fixed landcover

b = (; bounds=(0.0, 1.0))
# transitions = (;
#     a=(; a=(d1=Param(0.8; b...), d2=Param(0.1; b...), d3=Param(0.0; b...)), b=(d1=0.0, d2=0.0, d3=0.0), c=(d1=0.0, d2=0.0, d3=0.0)),
#     b=(; a=(d1=0.0, d2=0.0, d3=0.0), b=(d1=Param(0.8; b...), d2=Param(0.8; b...), d3=Param(0.8; b...)), c=(d1=0.0, d2=0.0, d3=0.0)),
#     c=(; a=(d1=0.0, d2=0.0, d3=0.0), b=(d1=0.0, d2=0.0, d3=0.0), c=(d1=Param(1.0; b...), d2=Param(0.8; b...), d3=Param(0.8; b...))),
# )

kw = (; bounds=(0.01, 1.0))
transitions = (;
    a=(; a=Exponential( Param(1.0; kw...)), b=0.0, c=0.0),
    b=(; a=0.0, b=Exponential(Param(1.0; kw...)), c=0.0),
    c=(; a=0.0, b=0.0, c=Exponential(Param(1.0; kw...))),
)
A = rand(1:3, 100, 100)
B = rand(1:3, 100, 100)
C = rand(1:3, 100, 100)
cover_change(A, B; categories=(x=1, y=2, z=3)) |> pairs
cover_persistence(A, B; categories=(z=1, x=2, y=3))

# using GLMakie, Makie
# GLMakie.activate!()

# x = Observable(Float32.(parent(init_state.state)))
# p = Makie.heatmap(x)
# fig = p.figure
# grid = GridLayout(fig[2, 1])

weight_rule = LandscapeChange.BottomUp(;
    neighborhood=Window(4),
    states,                                     
    inertia=(a=Param(0.8; b..., label="inertia a"), b=Param(0.8; b..., label="inertia b"), c=Param(0.8; b..., label="inertia c")),
    transitions,
    suitability=Aux{:suitability}(),
    pressure=(a=Param(0.0; b..., label="pressure a"), b=Param(0.0; b..., label="pressure b"), c=Param(0.0; b..., label="pressure c")),
    fixed,
    perturbation=Param(0.1; bounds=(0.0, 5.0), label="perturbation"),
)
ruleset = Ruleset(weight_rule; proc=CPUGPU());

output = MakieOutput(parent(init_state);
    aux=(; counts, suitability),
    tspan=1:.1:200,
    store=false,
    boundary=Wrap(),
    padval=1,
    ruleset,
) do fig, axis, frame
    hm = Makie.heatmap!(axis, frame; colorrange=(1, 3), colormap=:Juarez)
    hm = Makie.heatmap!(axis, frame; colorrange=(1, 3), colormap=:Juarez)
    Colorbar(fig[1, 2], hm)
end
display(output.fig)

# output.fig[3, 1] = ax

output = ResultOutput(parent(init_state);
    aux=(; suitability),
    tspan=1:10,
    store=false,
    padval=1,
    ruleset)
sim!(output, ruleset);


# using ProfileView
# @profview sim!(output, ruleset);
# using Cthulhu
# rs = StaticRuleset(ruleset)
# @descend isinferred(sd)
# sim!(output, ruleset);
# @descend 
# @time sim!(output, ruleset; proc=DynamicGrids.CPUGPU());
# using CUDAKernels
# @time sim!(output, ruleset; proc=DynamicGrids.CuGPU());
