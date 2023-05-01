"""
    BottomUp <: NeighborhoodRule

Follows "The use of constrained cellula automata for high-resolution modelling
of urban land-use dynamics", BottomUp and Engalin (1996)

Abstract:
> "A cellular automaton is specified to give a spatially detailed representation of the evolution
of urban land-use patterns. Cell states represent land uses, and transition rules express the like-
lihood of a change from one state to another as a function both of existing land use in the 113-cell
neighbourhood of the cell and of the inherent suitability of the cell for each possible use."

Here we can use any neighborhood shape size.

The core equation is:
``P_{hj} = vs_j \\left(1 + \\sum_{k,i,d} m_{kd}I_{id} \\right) + H_j``

where:

- ``P_{hj}`` is the transition potential from state h to state j ;
- ``m_{kd}``, is the weighting parameter applied to cells with state k in distance zone d; ( 1, if the state of cell / is equal to A:,
- ``I_{id} = \\left\\{\\begin{cases} \\frac{x^2-x}{x},& \\text{if} x\\geq 1, if the state of the cell i is equal to k 0, otherwise & \\text{otherwise} \\end{cases}``
    where ``i`` is the index of cells within the distance zone
    (``I_{id}`` simply selects the proper weight, ``m_{kd}``, to apply to the cell at location ``i, d``);
- ``H_j`` is an inertia parameter, with ``H_j > 0`` if ``j = h``, and ``H_j = 0`` otherwise;
- ``s_j`` is the suitability of the cell state for ``j``, ``0 \\leq s_j \\leq 1;
- ``v`` is a stochastic disturbance term ``v = l + [-ln(rand)]^\\alpha``

## Arguments
- `transition_potential`: `P` in the equation, a `NamedVector` of `NamedVector` - making a named matrix of potentials.
- `inertia`: `H` in the equation, a `NamedVector` of inertia values for each category.
- `suitability`: `s` in the equation - an `Aux` suitability array, another `Grid`, or a constant value.
    All values must be between [0, 1]
- `perturbation`: perturbation scalar. A single value, `NamedVector`, an `Aux` value or another `Grid`.
"""
struct BottomUp{R,W,N,St<:NamedVector,T,L,Su,H,P,F,α} <: NeighborhoodRule{R,W}
    neighborhood::N
    states::St
    transitions::T
    logic::L
    suitability::Su
    inertia::H
    pressure::P
    fixed::F
    perturbation::α
end
function BottomUp{R,W}(;
    neighborhood,
    states,
    transitions,
    logic=map(ts -> map(_ -> true, ts), transitions), # all true by default
    suitability,
    inertia=0.0,
    pressure=0.0,
    fixed=false,
    perturbation=0.0,
) where {R,W}
    length(transitions) == length(states) || throw(ArgumentError("Number of transitions $(length(transitions)) does not match number of states $(length(states))"))
    # length(logic) == length(states) || throw(ArgumentError("Number of logic $(length(logic)) does not match number of states $(length(states))"))
    length(inertia) == length(states) || throw(ArgumentError("Number of inertia values $(length(inertia)) does not match number of states $(length(states))"))
    ((pressure isa Function) || length(pressure) == length(states)) || throw(ArgumentError("Number of pressure values $(length(pressure)) does not match number of states $(length(states))"))
    # length(pressure) == length(states) || throw(ArgumentError("Number of pressure values $(length(pressure)) does not match number of states $(length(states))"))
    all(t -> length(t) == length(transitions), transitions) || throw(ArgumentError("transition lengths do not match"))
    BottomUp{R,W}(neighborhood, states, transitions, logic, suitability, inertia, pressure, fixed, perturbation)
end

function DynamicGrids.applyrule(data, rule::BottomUp, current_state::T, I) where T<:Integer
    current_state == zero(current_state) && return current_state
    # Cells with fixed states dont change
    get(data, rule.fixed, I) && return current_state
    # Get the future state from aux
    future_state::Int64 = get(data, Aux{:history}(), I)
    # In case there is a masked cell just change state
    future_state == zero(current_state) && return future_state
    # Don't change states if its already the future state 
    current_state == future_state && return current_state

    pressures = get(data, rule.pressure, I)
    suitabilities = get(data, rule.suitability, I)
    perturbation = get(data, rule.perturbation, I)
    inertias = get(data, rule.inertia, I)
    α = rule.perturbation

    # Calculate transition probabilities for each active state
    neighbor_weights = map(_ -> 0.0, suitabilities)
    nbrs = neighbors(rule)
    dist_zones = distance_zones(rule)
    for i in 1:length(dist_zones)
        if nbrs[i] != zero(current_state)
            neighbor_weights = map(neighbor_weights, rule.states) do σm, state
                σm + rule.transitions[Int(state)][nbrs[i]][dist_zones[i]+1]
            end
        end
    end

    # For sum, suitability and inertia of each class
    # calculate transition potentials P_hj
    transition_potentials = map(
            values(neighbor_weights),
            values(suitabilities),
            values(inertias),
            values(rule.states),
            values(pressures)
        ) do weight, suitability, inertia, potential_state, pressure
        # If we cant logically switch to this state return the typemin
        rule.logic.indirect[future_state][potential_state] || return typemin(typeof(weight))
        rule.logic.direct[potential_state][current_state] || return typemin(typeof(weight))
        # Define a stochastic disturbance term `v`
        # @fastmath v = 1.0 + (-log(rand()))^α
        # v = (1.0 + rand())::Float64
        # v = 1 + (rand() ^ 2 * rand((-1, 1))) * α
        noise = 1 + (rand()^5 * 2 - 1) * rule.perturbation
        # Inertia only applies when state == current_state
        H = inertia * (potential_state == current_state)
        return suitability * (1 + weight) * pressure + H + noise
    end

    # The index with highest probability is the next state
    _, i = findmax(SVector(transition_potentials))
    # Just set the cell to the next state now
    return convert(T, i)
end

function DynamicGrids.modifyrule(rule::BottomUp{grid}, data) where grid
    if rule.pressure isa Function
        pressure = rule.pressure(data, rule)
        println(stdout, pressure)
        @set! rule.pressure = pressure
    end
    @set rule.transitions = transitions_to_vectors(rule)
end

function transitions_to_vectors(rule)
    dists = sort(union(distances(neighborhood(rule))))
    transitions = map(rule.transitions) do xs
        map(xs) do x
            map(1:length(dists)) do i
                if x isa Distributions.Distribution
                    # Normalised exponential
                    Distributions.pdf(x, dists[i]) ./ Distributions.pdf(x, 0)
                else
                    x
                end
            end |> SVector{length(dists)}
        end
    end
end

# function DynamicGrids.validaterule(rule::BottomUp, data)
#     get(data, rule.fixed, I) && return h
#     pressure = get(data, rule.pressure, I)
#     suitability = get(data, rule.suitability, I)
#     perturbation = get(data, rule.perturbation, I)
#     inertia = get(data, rule.inertia, I)
#     α = rule.perturbation
# end
