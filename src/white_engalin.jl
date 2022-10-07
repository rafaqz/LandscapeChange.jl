"""
    WhiteEngalin <: NeighborhoodRule

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
- `transition_potential`:: `P` in the equation, a `NamedVector` of `NamedVector` or vectors.
- `inertia`:: `H` in the equation, a `NamedVector` of inertia values for each category.
- `suitability`:: `s` in the equation, an `Aux` suitability matrix, andother `Grid`, or a constant value.
"""
struct WhiteEngalinWeights{R,W,P,S,H<:NamedVector,A<:NamedVector,F<:NamedVector,α<:Number} <: NeighborhoodRule{R,W}
    transition_potential::P
    suitability::S
    inertia::H
    active::A
    fixed::F
    pertubation::α
end

function DynamicGrids.applyrule(data, rule::WhiteEngalinWeights, h, I)
    (; transition_potential, inertia, suitability, active, fixed) = rule
    # Cells with fixed states dont change
    k in fixed && return k
    suit = get(data, suitability, I)
    # Calculate transition probabilities for each active state
    Σw = ntuple(_ -> 0.0, length(active))
    for (k, d) in zip(neighbors(rule), distance_zones(rule)) 
        Σw = map(Σm, active) do σm, j
            σm + transition_potential[j][k][d] 
        end
    end
    P_hj = map(Σw, suit, inertia) do σw, s, Hmax
        # Stochastic term
        v = 1 + (-ln(rand()))^α
        # Inertia only applies when j == h
        H = Hmax * j == h
        # Main equation
        v * s * (1 + σw) + H
    end
    # Choose the most propable next state
    i, probability = findmax(P_hj)
    return NamedArray(; state=active[i], probability, index=I)
end

struct WhiteEngalinUpdate{R,W} <: SetGridRule{R,W} end

function DynamicGrids.applyrule(data, rule::WhiteEngalinUpdate, (weights, dest))
    required = get(data, rule.required)
    for nrequired in required
        # Sort by maximum probability for this state, for however many we require
        partialsort!(weights, 1:nrequired; 
            by=x -> x.state == state ? x.probability : zero(x.probability)
        )
        # Set the most probable indices to the new state
        broadcast(view(weights, 1:nrequired)) do x
            dest[x.I] = x.state
        end
    end
end
