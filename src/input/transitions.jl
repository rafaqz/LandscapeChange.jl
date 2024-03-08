"""
    indirect_transitions(transitions::NamedVector{K})

Calculate all indirect transitions in the graph recursively.
"""
function indirect_transitions(transitions::NamedVector{K}) where K
    states = NamedVector{K}(ntuple(identity, length(K)))
    return indirect_transitions(states, transitions)
end
function indirect_transitions(states, transitions)
    map(states) do s1
        map(states) do s2
            _can_change(states, transitions, s1, s2)
        end
    end
end

function reverse_transitions(transitions::NamedVector{K}) where K
    map(NamedVector{K}(K)) do k
        map(transitions) do l
            l[k]
        end
    end
end

function _can_change(states, transitions, to, from, checked=(), path=())
    if transitions[to][from]
        return true
    else
        return map(states) do s
            if s == to || s in checked
                false
            elseif transitions[to][s]
                _can_change(states, transitions, s, from, (checked..., to), (path..., from))
            else
                false
            end
        end |> any
    end
end

"""
    namedvector_raster(x, states)

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
    indirect=indirect_transitions(transitions),
    reversed=reverse_transitions(transitions),
    reversed_indirect=reverse_transitions(indirect),
    kw...
) where {K,T}
    states = NamedVector{K}(ntuple(identity, length(K)))
    dest = copy(source)
    N = size(timeline, Ti())
    logic = (; transitions, indirect, reversed, reversed_indirect)

    # Function barrier
    return _cross_validate_timeline!(dest, source, logic, Val(N))
end
function _cross_validate_timeline_inner!(dest::Raster{T}, source::Raster{T}, logic, ::Val{N}) where {N,K,T}
    for xy in DimIndices(dims(raster, (X(), Y())))
        pixel_timeline = ntuple(i -> raster[xy..., Ti(t)], N)
        filled = apply_both_transitions(pixel_timeline, logic)
        # Write to dest
        for (i, x) in enumerate(pixel_timeline)
            dest[xy..., Ti(i)] = x
        end
    end
    return dest
    # missingval = (timeline=zero(eltype(timeline)), error=zero(eltype(error)))
    # x = RasterStack((; timeline, error); missingval)
    # return rebuild(x; missingval)
end

function apply_both_transitions(
    timeline::Tuple,
    logic::NamedTuple, # = (; transitions, indirect, reversed, reversed_indirect)
) where K
    lower = upper = (;
        forced = zero(first(timeline)),
        uncertain = zero(first(timeline)),
    )
    newtimeline, _, _ = _combine(timeline, logic, lower, upper)
    return newtimeline
end

Base.@assume_effects :foldable _combine(timeline::Tuple{}, logic, lower::L, upper::L) where L = (), upper, lower
Base.@assume_effects :foldable function _combine(timeline::Tuple{T,Vararg}, logic, lower::L, upper::L) where {T,L}
    # Get the current timeline values
    tl = first(timeline)
    tu = last(timeline)
    ctl = count(tl)
    ctu = count(tu)
    lower_forced = ctl == 1
    upper_forced = ctu == 1

    # @show lower upper
    # @show tl tu
    # @show ctl ctu
    # Update the forced and uncertain categories based on the current values
    if lower_forced # We have 1 certain category: update forced and remove uncertain
        forced, match= merge_forced(tl, lower.forced, logic.transitions, logic.indirect)::Tuple{T,Bool}
        uncertain = zero(lower.uncertain)::T
        lower = (; forced, uncertain)::L
    elseif ctl != 0 # We have multiple uncertain categories: update uncertain
        # forced, match = merge_forced(tu, lower.forced, logic.transitions, logic.indirect)
        forced = lower.forced::T
        uncertain, match = merge_uncertain(tl, lower.uncertain, logic.reversed, logic.reversed_indirect)::Tuple{T,Bool}
        lower = (; forced, uncertain)::L
    end
    if upper_forced # We have 1 certain category: update forced and remove uncertain
        # TODO what if uncertain cant convert to certain
        forced, match = merge_forced(tu, upper.forced, logic.transitions, logic.indirect)::Tuple{T,Bool}
        uncertain = zero(upper.uncertain)::T # No uncertainty
        upper = (; forced, uncertain)::L
    elseif ctu != 0 # We have multiple uncertain categories: update uncertain
        # Existing forced states continue
        # forced, match = merge_forced(tu, lower.forced, logic.transitions, logic.indirect)
        forced = upper.forced
        # Merge uncertain states over transitions
        uncertain, match = merge_uncertain(tu, upper.uncertain, logic.transitions, logic.indirect)::Tuple{T,Bool}
        upper = (; forced, uncertain)::L
    end
    # @show lower upper

    # Recursively move inwards in the timeline
    innertimeline, innerlower, innerupper = _combine(timeline[begin+1:end-1], logic, lower, upper)

    # Update the timeline with combined values from forwards/backwards passes

    # println()
    # @show any(innerlower.forced) lower_forced
    if any(innerlower.forced)
        # @show "innerlower forced"
        if lower_forced
            finallower, match = merge_forced(innerlower.forced, tl, logic.transitions, logic.indirect)::Tuple{T,Bool}
            returnlower = (forced=finallower, uncertain=zero(finallower))::L
        else
            finallower, match = merge_forced(innerlower.forced, tl, logic.transitions, logic.indirect)::Tuple{T,Bool}
            if count(finallower) == 1
                returnlower = (forced=finallower, uncertain=zero(finallower))::L
            else
                returnlower = (forced=innerlower.forced, uncertain=finallower)::L
            end
        end
    else
        # @show "innerlower uncertain"
        if lower_forced
            finallower, match = merge_forced(innerlower.uncertain, tl, logic.transitions, logic.indirect)::Tuple{T,Bool}
            returnlower = (forced = tl, uncertain=finallower)::L
        else
            finallower, match = merge_uncertain(innerlower.uncertain, tl, logic.transitions, logic.indirect)::Tuple{T,Bool}
            returnlower = (forced = zero(finallower), uncertain=finallower)::L
        end
    end

    # println()
    # @show any(innerupper.forced) upper_forced
    if any(innerupper.forced)
        # @show "innerupper forced"
        if upper_forced
            finalupper, match = merge_forced(innerupper.forced, tu, logic.reversed, logic.reversed_indirect)::Tuple{T,Bool}
            # @show innerupper.forced tu finalupper logic.reversed
            returnupper = (forced=finalupper, uncertain=zero(finalupper))::L
        else
            finalupper, match = merge_uncertain(innerupper.forced, tu, logic.reversed, logic.reversed_indirect)::Tuple{T,Bool}
            if count(finalupper) == 1
                returnupper = (forced=finalupper, uncertain=zero(finalupper))::L
            else
                returnupper = (forced=innerupper.forced, uncertain=finalupper)::L
            end
        end
    else
        # @show "innerupper uncertain"
        if upper_forced
            finalupper, match = merge_uncertain(innerupper.uncertain, tu, logic.reversed, logic.reversed_indirect)::Tuple{T,Bool}
            returnupper = (forced=tu, uncertain=finalupper)::L
        else
            finalupper, match = merge_uncertain(innerupper.uncertain, tu, logic.reversed, logic.reversed_indirect)::Tuple{T,Bool}
            returnupper = (forced=zero(finalupper), uncertain=finalupper)::L
        end
    end
    # println()
    # @show finallower finalupper
    # println()

    # Update timeline to the best combined
    # values from forward and backawards passes
    
    # Check if there is one or two times to write to.
    # At the mid point we will only have one if we have
    # an odd number of layers.
    if length(timeline) == 1
        # We only have one line, not two. So combine it.
        
        # Keep all the forced categories
        combined = map(|, returnlower.forced, returnupper.forced)
        # If there are no forced, keep the minimum uncertain categories 
        if !any(combined)
            combined = map(&, finallower, finalupper)
        end
        # If there are no shared uncertain categories, 
        # keep all the uncertain categories 
        if !any(combined)
            combined = map(|, finallower, finalupper)
        end
        return (combined,), returnlower, returnupper
    else
        # We have two times, just return finallower and finalupper
        return (finallower, innertimeline..., finalupper), returnlower, returnupper
    end

    # Return updated forced and uncertain values    # Return updated forced and uncertain values    # Return updated forced and uncertain values    # Return updated forced and uncertain values    # Return updated forced and uncertain values    # Return updated forced and uncertain values    # Return updated forced and uncertain values    # Return updated forced and uncertain values
end

Base.@assume_effects :total function merge_uncertain(
    source::T, dest::T, transitions::L, indirect_transitions::L,
)::Tuple{T,Bool} where {T<:NamedVector{K},L<:NamedVector{K}} where K
    # Handle completely missing values
    any(source) || return dest, true
    any(dest) || return source, true
    # First see if they share values already and assume continuity
    shared = map(&, source, dest)
    any(shared) && return shared, true
    bitmasks = generate_bitmasks(source)
    # Otherwise find all possible transitions between uncertain states
    # We need to lock down dest in a let block for type stability
    let dest=dest
        @inline function f((acc, match), (s, direct_dest, indirect_dest, bitmask))
            m = match
            x = if s
                direct = map(&, dest, direct_dest)
                indirect = map(&, dest, indirect_dest)
                xs = if any(direct)
                    direct
                else
                    if any(indirect)
                        indirect
                    else
                        # @warn "Broken logic in merge_uncertain: keep this category, unmatch"
                        # @show source dest
                        # @show bitmask
                        m = true
                        bitmask
                    end
                end
                map(|, xs, acc)
            else
                acc
            end
            return (x, m)
        end
        foldl(f, zip(source, transitions, indirect_transitions, bitmasks); init=(zero(source), true))
    end
end

Base.@assume_effects :foldable function merge_forced(
    source::T, dest::T, transitions::L, indirect_transitions::L,
)::Tuple{T,Bool} where {T<:NamedVector{K},L<:NamedVector{K}} where K
    any(source) || return dest, true
    any(dest) || return source, true
    bitmasks = generate_bitmasks(source)
    f = let dest=dest
        function ((acc, match), (s, direct_dest, indirect_dest, bitmask))
            x = if s
                direct = map(&, dest, direct_dest)
                # @show source dest direct_dest direct
                xs = if any(direct)
                    direct
                else
                    indirect = map(&, dest, indirect_dest)
                    if any(indirect)
                        indirect
                    else
                        # @warn "Broken logic in merge_forced: merging"
                        # Broken logic: keep this category, unmatch
                        match = false
                        map(|, bitmask, dest)
                    end
                end
                map(|, xs, acc)
            else
                acc
            end
            return x, match
        end
    end
    return foldl(f, zip(source, transitions, indirect_transitions, bitmasks); init=(zero(source), true))
end

const LOG = []

Base.@assume_effects :foldable function merge_all(
    source::NamedVector{K}, dest::NamedVector{K},
    transitions::NamedVector{K}, indirect::NamedVector{K},
) where K
    any(source) || return dest
    bitmasks = generate_bitmasks(source)
    reduce(zip(source, transitions, indirect, bitmasks); init=(zero(source), true)) do (acc, match), (s, direct_dest, indirect_dest, bitmask)
        x = if s
            direct = map(&, dest, direct_dest)
            xs = if any(direct)
                direct
            else
                indirect = map(&, dest, indirect_dest)
                if any(indirect)
                    indirect
                else
                    # Broken logic: keep this category, unmatch
                    match = false
                    bitmask
                end
            end
            map(|, xs, acc)
        else
            acc
        end
        return x, match
    end
end

@generated function generate_bitmasks(::NamedVector{K}) where K
    nvs = ntuple(length(K)) do i
        vals = falses(length(K))
        vals[i] = true
        :(NamedVector{K}($(Tuple(vals))))
    end
    expr = Expr(:tuple, nvs...)
    return :(NamedVector{K}($expr))
end
