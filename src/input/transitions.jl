
function reverse_transitions(transitions::NamedVector{K}) where K
    map(NamedVector{K}(K)) do k
        map(transitions) do l
            l[k]
        end
    end
end

"""
    indirect_transitions(transitions::NamedVector{K})

Calculate all indirect transitions in the graph recursively.
"""
function indirect_transitions(transitions::NamedVector{K}) where K
    indices = NamedVector{K}(ntuple(identity, length(K)))
    return indirect_transitions(indices, transitions)
end
function indirect_transitions(indices, transitions)
    map(indices) do s1
        map(indices) do s2
            _can_change(indices, transitions, s1, s2)
        end
    end
end

all_transitions(trans::NamedVector) = (trans, next_transitions(trans)...)
function next_transitions(transitions::NamedVector{K}) where K
    indices = NamedVector{K}(ntuple(identity, length(K)))
    next = map(indices) do s1
        mapreduce(.|, indices) do s2
            if transitions[s1][s2]
                map(indices) do s3
                    transitions[s2][s3]
                end
            else
                map(_ -> false, indices)
            end
        end
    end
    if next == transitions
        ()
    else
        (next, next_transitions(next)...)
    end
end

function _can_change(indices, transitions, to, from, checked=(), path=())
    if transitions[to][from]
        return true
    else
        return map(indices) do s
            if s == to || s in checked
                false
            elseif transitions[to][s]
                _can_change(indices, transitions, s, from, (checked..., to), (path..., from))
            else
                false
            end
        end |> any
    end
end

"""
    namedvector_raster(x)

Collapse a stack into a raster of `NamedVector`.
"""
namedvector_raster(stack::RasterStack) = namedvector_raster(NamedTuple(stack))
function namedvector_raster(stack::NamedTuple{K}) where K
    map(stack...) do args...
        NamedVector{K}(args)
    end
end

"""
    stripe_raster(x, states)

Convert a raster of `NamedVector{Bool}` to a striped raster
that can be easily visualised.

This is a lossy method only usful for quick visualisation.
Never use the output of this function as data.

All `true` values will become stripes. Contiguous regions must be
as wide as the number of true categories for stripes to be
visible.
"""
function stripe_raster(ser::AbstractRasterSeries, states)
    map(ser) do raster
        stripe_raster(raster, states)
    end
end
function stripe_raster(raster::AbstractRaster, states; stripedim=X())
    A = rebuild(similar(raster, Int); missingval=0)
    for I in DimIndices(A)
        v = raster[I...]
        c = count(v)
        x = if c == 0
            0
        elseif c == 1
            states[findfirst(Tuple(v))]
        else
            # Stripe the possible values
            stripe = mod1(val(dims(I, stripedim)), c)
            findall(Tuple(v))[stripe]
        end
        A[I...] = x
    end
    return A
end

"""
    cross_validate_timeline(transitions, series::RasterSeries{<:Raster{<:NamedVector,2}}; kw...)
    cross_validate_timeline(transitions, raster::Raster{<:NamedVector,3}; kw...)

Compile a series of rasters into a coherent timeline,
following the logic defined in `transitions`.

Returns a `RasterStack` containing both the compiled time series,
and the error associated with forcing `transitions`.

# Keywords

-`simplify`: assume continuity between single valued pixels
    accross time. If any intermediate pixels include that value
    as well as others, assume there were no changes. This may
    vary between a realistic assumtion and a very unrealistic
    assumtion, depending on categories and the length of time.
-`transitions`: force all transition logic, forwards and backwards,
    recording all transition errors to a separate raster.
"""
cross_validate_timeline(transitions, ser::RasterSeries; kw...) =
    cross_validate_timeline(transitions, Rasters.combine(ser, Ti); kw...)
function cross_validate_timeline(transitions::NamedVector{K}, source::Raster{T};
    kw...
) where {K,T}
    dest = copy(source)
    N = size(timeline, Ti())
    logic = (; 
        transitions=all_transitions(transitions),
        reversed=all_transitions(reverse_transitions(transitions)),
    )

    # Function barrier
    return _cross_validate_timeline!(dest, source, logic, Val(N))
end
function _cross_validate_timeline!(dest::Raster{T}, source::Raster{T}, logic, ::Val{N}) where {N,K,T}
    for xy in DimIndices(dims(source, (X(), Y())))
        # Read a single timeline as a Tuple
        timeline = ntuple(i -> source[xy..., Ti(i)], N)
        # Update it with transition logic
        updated = apply_transitions(timeline, logic)
        # Write it to the dest array
        for (i, x) in enumerate(updated)
            dest[xy..., Ti(i)] = x
        end
    end
    return dest
    # missingval = (timeline=zero(eltype(timeline)), error=zero(eltype(error)))
    # x = RasterStack((; timeline, error); missingval)
    # return rebuild(x; missingval)
end

function apply_transitions(
    timeline::Tuple,
    logic::NamedTuple, # = (; transitions, reversed, reversed_indirect)
)
    up = down = (;
        forced = zero(first(timeline)),
        uncertain = zero(first(timeline)),
    )
    newtimeline, _, _ = _combine(timeline, logic, up, down, firstindex(timeline), lastindex(timeline))
    return newtimeline
end

struct Forced end
struct Uncertain end

Base.@assume_effects :foldable function _combine(::Tuple{}, logic, upfrom::L, downfrom::L, i, j) where L
    return (), downfrom, upfrom
end
Base.@assume_effects :foldable function _combine(timeline::Tuple{<:Any}, logic, upfrom::L, downfrom::L, i, j) where L
    t = first(timeline)
    # println("\n++++++++ Start single center +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++")
    # @show t downfrom upfrom

    # Combine timeline and carried over
    _, up = _merge_to(upfrom, t, logic.reversed)
    _, down = _merge_to(downfrom, t, logic.transitions)
    final = if any(up.forced)
        if any(down.forced)
            map(|, up.forced, down.forced)
        else
            combined = map(&, up.forced, down.forced)
            if any(combined)
                combined
            else
                # TODO this could be reduced further
                map(|, up.forced, down.forced)
            end
        end
    else
        if any(down.forced)
            combined = map(&, up.uncertain, down.forced)
            if any(combined)
                combined
            else
                map(|, up.uncertain, down.forced)
            end
        else
            combined = map(&, up.uncertain, down.uncertain)
            if any(combined)
                combined
            else
                map(|, up.uncertain, down.uncertain)
            end
        end
    end
    # @show final
    # println("++++++++ End single center +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n")
    return (final,), down, up
end
Base.@assume_effects :foldable function _combine(timeline::Tuple{T,Vararg}, logic, upfrom::L, downfrom::L, i, j) where {T,L}
    # println("\n=== ", (i, j), " ============================================================================")
    # Get the current timeline values
    tup = first(timeline)
    tdown = last(timeline)

    # Combine timeline and previous states up and down
    _, merged_up = _merge_to(upfrom, tup, logic.reversed)
    _, merged_down = _merge_to(downfrom, tdown, logic.transitions)

    # @show upfrom tup merged_up
    # println("*")
    # @show downfrom tdown merged_down
    # println("---- ", (i, j), " ----------------------------------------------------------------------")

    # Recursively move inwards in the timeline
    innertimeline, upto, downto = _combine(timeline[begin+1:end-1], logic, merged_up, merged_down, i+1, j-1)

    # println("X")
    # println("\n=== ", (i, j), " ============================================================================")
    # ntup = count(tup)
    # ntdown = count(tdown)
    # up_is_forced = ntup == 1
    # down_is_forced = ntdown == 1
    # idm = any(downto.forced) ? :forced : :uncertain
    # dm = down_is_forced ? :forced : :uncertain
    # ium = any(upto.forced) ? :forced : :uncertain
    # um = up_is_forced ? :forced : :uncertain

    # # Update the timeline with merged values from forwards/backwards passes
    # println("Up to:")
    finalup, next_up_to = _merge_to(upto, merged_up, logic.transitions)
    # println("Down to:")
    finaldown, next_down_to = _merge_to(downto, merged_down, logic.reversed)
    new_timeline = (finalup, innertimeline..., finaldown)

    # @show tup upfrom merged_up upto finalup next_up_to
    # println("*")
    # @show tdown downfrom merged_down downto finaldown next_down_to
    # println("---- ", (i, j), " ----------------------------------------------------------------------")
    # println()
    return new_timeline, next_up_to, next_down_to
end

# Get source and dest modes
Base.@assume_effects :foldable function _merge_to(source::T, dest::T, transitions) where {T<:NamedTuple}
    sm, s = any(source.forced) ? (Forced(), source.forced) : (Uncertain(), source.uncertain)
    dm, d = any(dest.forced) ? (Forced(), dest.forced) : (Uncertain(), dest.uncertain)
    return _merge_to(sm, dm, s, d, transitions)
end
Base.@assume_effects :foldable function _merge_to(source::L, dest::T, transitions) where {L<:NamedTuple,T<:NamedVector}
    sm, s = any(source.forced) ? (Forced(), source.forced) : (Uncertain(), source.uncertain)
    dm, d = (count(dest) == 1 ? Forced() : Uncertain()), dest
    return _merge_to(sm, dm, s, d, transitions)
end
# Handle assigning state to forced or uncertain
Base.@assume_effects :foldable @inline function _merge_to(
    sourcemode, destmode, source::T, dest::T, transitions
) where {T<:NamedVector}
    final = _merge(sourcemode, destmode, source, dest, transitions)
    next_dest = if sourcemode isa Uncertain
        # Move single category from uncertain to forced
        if destmode isa Forced || count(final) == 1
            (forced=final, uncertain=zero(final))
        else
            (forced=zero(final), uncertain=final)
        end
    else
        (forced=final, uncertain=zero(final))
    end
    return final, next_dest
end

Base.@assume_effects :foldable function _merge(
    sourcemode, destmode, source::T, dest::T, transitions::Tuple
)::T where {T<:NamedVector{K}} where K
    # println("********************** $sourcemode $destmode **************************")
    # Just return identical source and dest
    source == dest && return dest
    # Handle completely missing source or dest
    any(source) || return dest
    # We minimise forced values in the source throught its possible 
    # transitions, they may resolve to one or two distinct timelines. 
    any(dest) || return _minimise_states(sourcemode, source, transitions...)

    # Otherwise find transitions between uncertain states.
    # This loop handles merging previous and current state.
    # If a previous state can move to a current state directly or indirectly
    # it is replaced with the dest state, otherwise kept.
    # Where no source can transition to them, dest states are dropped.
    # For Forced mode, all dest state is added aftwerwards.
    bitmasks = generate_bitmasks(source)
    merged = zero(source)
    for i in eachindex(source)
        if source[i]
            # Loop over all transition distances until one/multiple possible transitions are found
            transitioned = zero(source)
            for t in transitions
                attempt = map(&, dest, t[i])
                if any(attempt)
                    transitioned = attempt
                    break
                end
            end
            new_states = if any(transitioned)
                transitioned
            else
                # Nothing matching was found
                if sourcemode == Forced()
                    # If source is forced return it 
                    bitmasks[i]
                else
                    # Otherwise remove this state
                    zero(source)
                end
            end
            merged = map(|, merged, new_states)
        end
    end
    # This only acts on Forced, which will add dest state back afterwards
    minimum = _minimise_states(destmode, merged, transitions...)
    # Finish the merge for specific sourcemode and destmode combinations
    return _merge1(sourcemode, destmode, minimum, source, dest)
end

@inline function _merge1(::Forced, ::Forced, out, source, dest)
    # dest is forced so keep it all
    return map(|, out, dest)
end
@inline function _merge1(::Uncertain, ::Forced, out, source, dest)
    # Don't add extra uncertain states
    if map(&, source, dest) == dest
        return dest
    else
        # Dest is forced so keep it all
        return map(|, out, dest)
    end
end
@inline function _merge1(::Forced, ::Uncertain, out, source, dest)
    # Don't add extra uncertain states
    map(&, source, out) == source ? source : out
end
@inline function _merge1(::Uncertain, ::Uncertain, out, source, dest)
    # See if source and dest share values without needing transitions
    shared = map(&, source, dest)
    any(shared) ? shared : out
end

function _minimise_states(::Forced, source::NamedVector{K}, transitions...) where K
    bitmasks = generate_bitmasks(source)
    indices = NamedVector{K}(ntuple(identity, length(K)))
    mapreduce(.|, source, indices, bitmasks) do s1, i1, bm
        if s1
            possible_transitions = map(indices, source) do i2, s2
                s2 && transitions[end][i1][i2]
            end
            if count(possible_transitions) > 1
                zero(source)
            else
                bm
            end
        else
            zero(source)
        end
    end
end
_minimise_states(::Uncertain, source::NamedVector, transitions...) = source

@generated function generate_bitmasks(::NamedVector{K}) where K
    nvs = ntuple(length(K)) do i
        vals = falses(length(K))
        vals[i] = true
        :(NamedVector{K}($(Tuple(vals))))
    end
    expr = Expr(:tuple, nvs...)
    return :(NamedVector{K}($expr))
end

