"""
    HistoricalWeights <: NeighborhoodRule

Follows "The use of constrained cellula automata for high-resolution modelling
of urban land-use dynamics", White and Engalin (1996)

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
struct HistoricalWeights{R,W,N,St<:NamedVector,T,Su,H,F,α} <: NeighborhoodRule{R,W}
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
function HistoricalWeights{R,W}(; 
    stencil,
    states,
    transitions,
    logic=map(ts -> map(_ -> true, ts), transitions), # all true by default
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
    HistoricalWeights{R,W}(stencil, states, transitions, logic, suitability, history, inertia, pressure, fixed, perturbation)
end

function DynamicGrids.applyrule(data, rule::HistoricalWeights, h::Integer, I)
    # Cells with fixed states dont change
    get(data, rule.fixed, I) && return active

    suitability = get(data, rule.suitability, I)
    perturbation = get(data, rule.perturbation, I)
    inertia = get(data, rule.inertia, I)
    α = rule.perturbation
    
    # Calculate transition probabilities for each active state
    Σw = map(_ -> 0.0, suitability)
    for (k, d) in zip(neighbors(rule), distance_zones(rule))
        Σw = map(Σw, rule.states) do σm, j
            σm + rule.transitions[Int(j)][Int(k)][Int(d)]
        end
    end
    # For sum, suitability and inertia of each class
    # calculate transition potentials P_hj
    probabilities = map(Σw, suitability, inertia, rule.states) do σw, s, Hmax, j
        # Define a stochastic disturbance term `v`
        @fastmath v = 1.0 + (-log(rand()))^α
        # v = (1.0 + rand())::Float64
        # v = 1 + rand() * α
        # Inertia only applies when j == h 
        H = Hmax * (j == h)
        # Equation 1
        v * s * (1 + σw) + H
    end
    # Choose the most probable next state
    return Weights(probability, CartesianIndex(I))
end

struct Weights{W}
    weights::W
    I::CartesianIndex{2}
end
Weights(; weights, I) = Weights(weights, I) 

struct HistoricalUpdate{R,W,N,S} <: SetGridRule{R,W} 
    nrequired::N
    states::S
end

Base.zero(::T) where T<:StateWeight = Base.zero(T)
Base.zero(::Type{<:StateWeight}) =
    Base.zero(StateWeight{UInt8,Float64})
Base.zero(::Type{<:StateWeight{S,W}}) where {S,W} =
    StateWeight(zero(S), zero(W), CartesianIndex(0,0))

function DynamicGrids.applyrule!(data, rule::HistoricalUpdate{Tuple{},Tuple{Weight,State}}
) where {Weight,State}
    w = deepcopy(parent(source(data[Weight])))
    raw = vec(w)
    nrequired = get(data, rule.nrequired)

    # For how many of each state are reuired for this timestep
    foreach(nrequired, rule.states) do nreq, state
        # @show state nreq
        # Sort by maximum probability for this state, for however many we require
        partialsort!(raw, 1:nreq;
            by=x -> x.state == state ? x.weight : zero(x.weight),
            rev=true,
        )
        # Set the most probable indices to the new state
        most_probable = view(raw, 1:nreq)
        foreach(most_probable) do x
            if checkbounds(Bool, data[State], x.I)
                data[State][x.I] = state
                w[x.I] = zero(x)
            end
        end
    end
    return nothing
end


