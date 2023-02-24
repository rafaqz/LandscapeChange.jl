using LandscapeChange, NeutralLandscapes, DimensionalData, Plots
using CUDAKernels, CUDA

# Define the states we are working with
states = NamedVector(a=1, b=2)

s = 500
# Define suitability maps
suitability_a = rand(PerlinNoise((3, 3)), X(s), Y(s); name="suitability a") ./ 5 .+ 0.3
suitability_b = rand(PerlinNoise((4, 4)), X(s), Y(s); name="suitability b")
plot(
    heatmap(suitability_a; aspect_ratio=1), 
    heatmap(suitability_b; aspect_ratio=1),
)
suitability = map(suitability_a, suitability_b) do a, b
    NamedVector(; a, b)
end

init_state = (
    state = init_state = UInt8.((suitability_b .> (maximum(suitability_b) - 0.05)) .+ 1),
    weight = map(_ -> zero(LandscapeChange.StateWeight), suitability_a),
)

# Set an initial state where b has started in
# its absolute most suitable areas
heatmap(init_state.state)
        
# Calculate change rates to drive the simulation
growth = 0:0.2:10 # 10x growth over 10 timesteps
start_b = count(==(states.b), init_state.state)
b_counts = DimArray(round.(Int, growth .* start_b), Ti)
a_counts = length(suitability_a) .- b_counts
counts = DimArray(map(a_counts, b_counts) do a, b
    NamedVector(; a, b)
end, Ti)

SArray(inertia)
fixed = false # or a raster mask of fixed landcover

# transition_distributions = = (;
#     a=(; a=Poisson(0.5), b=nothing,
#     b=(; a=nothing, b=Geometric(0.2),
# )
# transitions = map(transition_distributions) do layers
#     map(layers) do layer
#         if isnothing(layer)
#         else
#         end
#     end
# end
transitions = (;
    a=(; a=(0.2, 0.1, 0.00), b=(0.0, 0.0, 0.0)),
    b=(; a=(0.0, 0.0, 0.0), b=(1.2, 0.8, 0.07)),
)

A = rand(1:3, 100, 100)
B = rand(1:3, 100, 100)
cover_change(A, B; categories=(x=1, y=2, z=3)) |> pairs
cover_persistence(A, B; categories=(z=1, x=2, y=3))

weight_rule = WhiteEngalinUljeeWeights{:state,:weight}(;
    neighborhood=Window(1),
    states,                                     
    inertia,
    transitions,
    suitability=Aux{:suitability}(),
    fixed,
    perturbation=0.001,
)
update_rule = WhiteEngalinUljeeUpdate{Tuple{},Tuple{:weight,:state}}(Aux(:counts), states)

rs = (weight_rule, update_rule)

using GLMakie
GLMakie.activate!()

x = Observable(Float32.(parent(init_state.state)))
Makie.heatmap(x)

output = TransformedOutput(init_state; 
    aux=(; counts, suitability),
    tspan=1:0.1:50,
    store=false,
    padval=(state=0x01, weight=zero(eltype(init_state.weight))), 
) do layers
    x[] .= Float32.(parent(layers.state))
    notify(x)
    sleep(0.001)
    nothing
end
# WGLMakie.activate!()
@time sim!(output, rs);

function JSServe.serialize_binary(session::JSServe.Session, msg::JSServe.SerializedMessage)
    return transcode(Noop, JSServe.MsgPack.pack(msg))
end

# anim = @animate for A in output
#     plot(
#         heatmap(A.state; aspect_ratio=:equal, title=:state), 
#         heatmap(getfield.(A.weight, :weight); aspect_ratio=:equal, title=:weight),
#     )
# end
# gif(anim, "weu.gif", fps=5)
