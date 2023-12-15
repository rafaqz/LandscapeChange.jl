"""
    BottomUp <: NeighborhoodRule

Follows "The use of constrained cellula automata for high-resolution modelling
of urban land-use dynamics", BottomUp and Engalin (1996)

Abstract:
> "A cellular automaton is specified to give a spatially detailed representation of the evolution
of urban land-use patterns. Cell states represent land uses, and transition rules express the like-
lihood of a change from one state to another as a function both of existing land use in the 113-cell
stencil of the cell and of the inherent suitability of the cell for each possible use."

Here we can use any stencil shape and size.

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
struct BottomUp{R,W,Ste,Sta<:NamedVector,T,L,Su,Hi,H,P,F,α} <: NeighborhoodRule{R,W}
    stencil::Ste
    states::Sta
    transitions::T
    logic::L
    suitability::Su
    history::Hi
    inertia::H
    pressure::P
    fixed::F
    perturbation::α
end
function BottomUp{R,W}(;
    stencil,
    states,
    transitions,
    logic,
    suitability,
    history,
    inertia=0.0,
    pressure=0.0,
    fixed=false,
    perturbation=0.0,
) where {R,W}
    length(transitions) == length(states) || throw(ArgumentError("Number of transitions $(length(transitions)) does not match number of states $(length(states))"))
    # length(logic) == length(states) || throw(ArgumentError("Number of logic $(length(logic)) does not match number of states $(length(states))"))
    ((pressure isa Function) || length(pressure) == length(states)) || throw(ArgumentError("Number of pressure values $(length(pressure)) does not match number of states $(length(states))"))
    # length(pressure) == length(states) || throw(ArgumentError("Number of pressure values $(length(pressure)) does not match number of states $(length(states))"))
    all(t -> length(t) == length(transitions), transitions) || throw(ArgumentError("transition lengths do not match"))
    BottomUp{R,W}(stencil, states, transitions, logic, suitability, history, inertia, pressure, fixed, perturbation)
end

function DynamicGrids.applyrule(data, rule::BottomUp, current_state::T, I) where T<:Integer
    current_state == zero(current_state) && return current_state::T
    # Cells with fixed states dont change
    get(data, rule.fixed, I) && return current_state::T
    # Get the future state from aux
    future_states = get(data, rule.history, I)
    # In case there is no future state
    any(future_states) || return current_state::T
    # Don't change states if its already the future state 
    future_states[current_state] && return current_state::T

    pressures = get(data, rule.pressure, I)
    suitabilities = get(data, rule.suitability, I)
    perturbation = get(data, rule.perturbation, I)
    inertias = get(data, rule.inertia, I)

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
            values(pressures),
        ) do weight, suitability, inertia, potential_state, pressure
        if potential_state != current_state 
            # If we cant directly switch return the typemin
            rule.logic.direct[potential_state][current_state] || return typemin(typeof(weight))
            # But maybe we can switch indirectly
            found = false
            for (to, s) in enumerate(future_states)
                s || continue
                from = potential_state
                if rule.logic.indirect[to][from] 
                    found = true
                end
            end 
            found || return typemin(typeof(weight))
        end
        # rule.logic.indirect[i][potential_state] || return typemin(typeof(weight))
        # any(map(rule.logic.indirect) do ind
            # any(ind[potential_state] .& future_states)
        # end) || 
        # Define a stochastic disturbance term `v`
        # @fastmath v = 1.0 + (-log(rand()))^α
        # v = (1.0 + rand())::Float64
        # v = 1 + (rand() ^ 2 * rand((-1, 1))) * α
        # v = 0.5 + (rand()^5 * 2 - 1) * rule.perturbation
        
        # noise = 1 + (rand()^5 * 2 - 1) * rule.perturbation
        # Inertia only applies when state == current_state
        # weight determines which cells will shift first
        # pressure determines if any will shift at all
        # weight mus always be less than inertia
        if isinf(pressure)
            return pressure
        else
            v = 1.0 + -log(rand()^rule.perturbation)
            return v * weight * pressure + inertia * (potential_state == current_state)
        end
        # v * pressure
    end
    # @show transition_potentials

    # The index with highest probability is the next state
    _, i = findmax(transition_potentials)
    # i != current_state && @show i transition_potentials
    # Just set the cell to the next state now
    return convert(T, i)::T
end

function DynamicGrids.modifyrule(rule::BottomUp{grid}, data) where grid
    if rule.inertia isa Number
        @set! rule.inertia = map(_ -> rule.inertia, rule.states) 
    end
    if rule.pressure isa Function
        pressure = rule.pressure(data, rule)
        @set! rule.pressure = pressure
    end
    # @show pairs(NamedTuple(rule.pressure))
    @set rule.transitions = transitions_to_vectors(rule)
end

function transitions_to_vectors(rule)
    dists = sort(union(distances(stencil(rule))))
    transitions = map(rule.transitions) do xs
        map(xs) do x
            map(dists) do d
                if x isa Distributions.Distribution
                    # Normalised exponential
                    Distributions.pdf(x, d) ./ Distributions.pdf(x, 0)
                elseif x isa Function
                    x(d)
                else
                    x
                end
            end |> SVector{length(dists)}
        end
    end
    return transitions
end

function calc_pressure(leverage, current, predicted, allowed)
    #leverage * (p - n)# / max(a, oneunit(a))
    # 10 * sign(x) * log(abs(x))
    # leverage * (p / n)
    # leverage * max(1.0, predicted) / max(predicted - current, 1.0)#,  max(1.0, predicted) / smoothing)
    needed = (predicted - current)
    # if needed > 0
        (leverage * needed / max(allowed, 1.0))
    # else
        # (leverage *  / max(allowed, 1.0))
    # end
end

leverage = 10
smoothing = 1
current = 1000
predicted = 4500
allowed = 10000
calc_pressure(leverage, current, predicted, allowed)

# function DynamicGrids.validaterule(rule::BottomUp, data)
#     get(data, rule.fixed, I) && return h
#     pressure = get(data, rule.pressure, I)
#     suitability = get(data, rule.suitability, I)
#     perturbation = get(data, rule.perturbation, I)
#     inertia = get(data, rule.inertia, I)
#     α = rule.perturbation
# end


