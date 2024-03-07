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
function cross_validate_timeline(transitions::NamedVector{K}, raster::Raster{T}; 
    kw...
) where {K,T}
    states = NamedVector{K}(ntuple(identity, length(K)))
    # Get time on the first dimension for performance
    timeline = copy(raster)#permutedims(raster, (Ti(), X(), Y()))
    error = fill!(similar(timeline), zero(T))
    forward = fill!(timeline[X(1), Y(1)], zero(T))
    backward = fill!(timeline[X(1), Y(1)], zero(T))
    possible = fill!(timeline[X(1), Y(1)], zero(T))
    for xy in DimIndices(dims(raster, (X(), Y())))
        pixel_timeline = view(timeline, xy...)
        pixel_error = view(error, xy...)
        forward .= pixel_timeline
        backward .= pixel_timeline
        cross_validate_timeline!(pixel_timeline, transitions;
            pixel_error, forward, backward, possible, kw...
        )
    end
    missingval = (timeline=zero(eltype(timeline)), error=zero(eltype(error)))
    x = RasterStack((; timeline, error); missingval)
    return rebuild(x; missingval)
end

function cross_validate_timeline!(pixel_timeline, transitions;
    indirect=indirect_transitions(transitions),
    reversed=reverse_transitions(transitions),
    reversed_indirect=reverse_transitions(indirect),
    pixel_error=fill!(copy(pixel_timeline), zero(eltype(pixel_timeline))),
    forward=copy(pixel_timeline),
    backward=copy(pixel_timeline),
    possible=fill!(copy(pixel_timeline), zero(eltype(pixel_timeline))),
    force=map(_ -> false, transitions), 
    simplify=true, 
    cull=true, 
    fill_only=false,
)
    rev = lastindex(pixel_timeline):-1:firstindex(pixel_timeline)
    if fill_only
        _fill_empty_times!(pixel_timeline)
    else
        # Forward
        apply_transitions!(forward, reversed, reversed_indirect)
        # Backwards
        apply_transitions!(view(backward, rev), transitions, indirect)
        apply_both_transitions!(view(backward, rev), transitions, indirect, reversed, reversed_indirect)

        # if simplify
            # minimise_uncertainty!(forward, pixel_timeline, possible, transitions, reversed, indirect, reversed_indirect)
            # minimise_uncertainty!(view(backward, rev), view(pixel_timeline, rev), view(possible, rev), reversed, transitions, reversed_indirect, indirect)
            # minimise_uncertainty!(forward)
            # minimise_uncertainty!(view(backward, rev))
        # end

        for i in eachindex(possible)
            f, b = forward[i], backward[i]
            possible[i] = map(|, f, b)
        end
        # _fill_empty_times!(possible)

        # if cull
        #     # Forward
        #     cull_transitions!(forward, pixel_timeline, possible, reversed, reversed_indirect)
        #     # Backwards
        #     cull_transitions!(view(backward, rev), view(pixel_timeline, rev), view(possible, rev), transitions, indirect)
        # end

        # if simplify
        #     minimise_uncertainty!(forward, pixel_timeline, possible, transitions, reversed, indirect, reversed_indirect)
        #     minimise_uncertainty!(view(backward, rev), view(pixel_timeline, rev), view(possible, rev), reversed, transitions, reversed_indirect, indirect)
        # end

        # if any(force)
        #     for i in eachindex(possible)
        #         x = map(&, force, possible[i])
        #         if any(x)
        #             possible[i] = x
        #         end
        #     end
        # end

        # Finalise output
        for i in eachindex(pixel_timeline)
            pixel_timeline[i] = possible[i]
            # before = pixel_timeline[i]
            # after = possible[i]
            # Any shared state is enough to say there was no error
            # shared = map(&, before, after)
            # Write errors
            # if !any(shared)
                # pixel_error[i] = before
            # end
            # Write timeline
            # pixel_timeline[i] = after
        end
    end

    return pixel_timeline, pixel_error
end

function apply_transitions!(
    timeline::AbstractArray{<:NamedVector{K}}, 
    transitions::NamedVector{K}, 
    indirect::NamedVector{K},
) where K
    a = first(timeline)
    for i in firstindex(timeline)+1:lastindex(timeline)
        b = timeline[i]
        if any(b)
            a, match = merge_all(a, b, transitions, indirect)
            timeline[i] = a
        end
    end
    return timeline
end


function apply_both_transitions!(
    timeline::AbstractArray{<:NamedVector{K}}, 
    logic::NamedTuple,# = (; transitions, indirect, reversed, reversed_indirect)
) where K
    lower = upper = (;
        forced = zero(first(timeline)),
        uncertain = zero(first(timeline)),
    )
    i = firstindex(timeline)
    j = lastindex(timeline)
    _combine!(timeline, logic, lower, upper, i, j)
    return timeline
end

function _combine!(timeline, logic, lower, upper, i, j)
    if i > j 
        return upper, lower
    end

    # Get the current timeline values
    tl = timeline[i]
    tu = timeline[j]
    ctl = count(tl) 
    ctu = count(tu) 
    lower_forced = ctl == 1
    upper_forced = ctu == 1

    # @show lower upper 
    @show tl tu
    @show ctl ctu 
    # Update the forced and uncertain categories based on the current values
    if lower_forced # We have 1 certain category: update forced and remove uncertain
        forced, match = merge_forced(tl, lower.forced, logic.transitions, logic.indirect) 
        uncertain = zero(lower.uncertain)
        lower = (; forced, uncertain)
    elseif ctl != 0 # We have multiple uncertain categories: update uncertain
        # forced, match = merge_forced(tu, lower.forced, logic.transitions, logic.indirect) 
        forced = lower.forced
        uncertain, match = merge_uncertain(tl, lower.uncertain, logic.reversed, logic.reversed_indirect) 
        lower = (; forced, uncertain) 
    end
    if upper_forced # We have 1 certain category: update forced and remove uncertain
        # TODO what if uncertain cant convert to certain
        forced, match = merge_forced(tu, upper.forced, logic.transitions, logic.indirect) 
        uncertain = zero(upper.uncertain) # No uncertainty
        upper = (; forced, uncertain)
    elseif ctu != 0 # We have multiple uncertain categories: update uncertain
        # Existing forced states continue
        # forced, match = merge_forced(tu, lower.forced, logic.transitions, logic.indirect) 
        forced = upper.forced
        # Merge uncertain states over transitions
        uncertain, match = merge_uncertain(tu, upper.uncertain, logic.transitions, logic.indirect) 
        upper = (; forced, uncertain) 
    end
    @show lower upper

    # Recursively move inwards in the timeline
    innerlower, innerupper = _combine!(timeline, logic, lower, upper, i+1, j-1)

    # Update the timeline with combined values from forwards/backwards passes
    
    println()
    @show any(innerlower.forced) lower_forced
    if any(innerlower.forced)
        @show "innerlower forced"
        if lower_forced
            finallower, match = merge_forced(innerlower.forced, tl, logic.transitions, logic.indirect)
        else
            finallower, match = merge_uncertain(innerlower.forced, tl, logic.transitions, logic.indirect)
        end
        returnlower = (forced=finallower, uncertain=zero(finallower))
    else
        @show "innerlower uncertain"
        if lower_forced
            finallower, match = merge_forced(innerlower.uncertain, tl, logic.transitions, logic.indirect)
        else
            finallower, match = merge_uncertain(innerlower.uncertain, tl, logic.transitions, logic.indirect)
        end
        returnlower = (forced = zero(finallower), uncertain=finallower)
    end

    println()
    @show any(innerupper.forced) upper_forced
    if any(innerupper.forced) 
        @show "innerupper forced"
        if upper_forced
            finalupper, match = merge_forced(innerupper.forced, tu, logic.reversed, logic.reversed_indirect)
        else
            finalupper, match = merge_uncertain(innerupper.forced, tu, logic.reversed, logic.reversed_indirect)
        end
        returnupper = (forced = zero(finalupper), uncertain=finalupper)
    else
        @show "innerupper uncertain"
        if upper_forced
            finalupper, match = merge_uncertain(innerupper.uncertain, tu, logic.reversed, logic.reversed_indirect)
        else
            finalupper, match = merge_uncertain(innerupper.uncertain, tu, logic.reversed, logic.reversed_indirect)
        end
        returnupper = (forced = zero(finalupper), uncertain=finalupper)
    end
    println()
    @show finallower finalupper
    println()
    # Update timeline to the best combined 
    # values from forward and backawards pass
    timeline[i] = finallower
    timeline[j] = finalupper

    # Return updated forced and uncertain values    # Return updated forced and uncertain values    # Return updated forced and uncertain values    # Return updated forced and uncertain values    # Return updated forced and uncertain values    # Return updated forced and uncertain values    # Return updated forced and uncertain values    # Return updated forced and uncertain values
    return returnlower, returnupper
end

Base.@assume_effects :foldable function merge_uncertain(
    source::NamedVector{K}, dest::NamedVector{K}, 
    transitions::NamedVector{K}, indirect::NamedVector{K},
) where K
    # Handle completely missing values
    any(source) || return dest, true
    any(dest) || return source, true
    # First see if they share values already and assume continuity
    shared = map(&, source, dest)
    any(shared) && return shared, true
    # Otherwise find all possible transitions between uncertain states
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
                    @warn "Broken logic in merge_uncertain: keep this category, unmatch"
                    @show source dest
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

Base.@assume_effects :foldable function merge_forced(
    source::NamedVector{K}, dest::NamedVector{K}, 
    transitions::NamedVector{K}, indirect_transitions::NamedVector{K},
) where K
    any(source) || return dest, true
    any(dest) || return source, true
    bitmasks = generate_bitmasks(source)
    reduce(zip(source, transitions, indirect_transitions, bitmasks); 
        init=(zero(source), true)
    ) do (acc, match), (s, direct_dest, indirect_dest, bitmask)
        x = if s
            direct = map(&, dest, direct_dest)
            xs = if any(direct) 
                direct
            else
                indirect = map(&, dest, indirect_dest)
                if any(indirect) 
                    indirect
                else
                    @warn "Broken logic in merge_forced: merging"
                    # Broken logic: keep this category, unmatch
                    match = false
                    map(|, source, dest)
                end
            end
            map(|, xs, acc)
        else
            acc
        end
        return x, match
    end
end

# Fill empty slices with previous data
# function _fill_empty_times!(timeline::AbstractVector)
#     last_non_empty_i = typemax(Int)
#     last_non_empty_categories = zero(eltype(timeline))
#     for i in eachindex(timeline)
#         present_categories = timeline[i]
#         any(present_categories) || continue # No category data for this time slice
#         if last_non_empty_i < (i - 1)
#             # There has been a gap in data, fill it with the combination of
#             # the last non empty categories and the present category
#             fill_categories = map(|, last_non_empty_categories, present_categories)
#             for n in last_non_empty_i+1:i-1
#                 timeline[n] = fill_categories
#             end
#         end
#         last_non_empty_i = i
#         last_non_empty_categories = present_categories
#     end
# end

# function minimise_uncertainty!(timeline)
#     bitmasks = generate_bitmasks(first(timeline))
#     singles = zero(first(timeline))
#     lastsingle_i = firstindex(timeline)
#     for i in eachindex(timeline)
#         x = timeline[i]
#         if count(x) == 1
#             singles = map(|, singles, x)
#             # Only keep states that where alone somewhere
#             for i in lastsingle_i+1:i-1
#                 timeline[i] = timeline[i] .& singles
#             end
#             lastsingle_i = i
#         end
#     end
#     rev = lastindex(timeline):-1:firstindex(timeline)
#     simplify_end!(timeline)
#     simplify_end!(view(timeline, rev))
#     return timeline
# end

# function minimise_uncertainty!(
#     timeline::AbstractVector{<:NamedVector}, known::AbstractVector{<:NamedVector}, possible::AbstractVector{<:NamedVector}, transitions, reversed, indirect, reversed_indirect
# )
#     past = first(known)
#     last_single_i = typemax(Int)
#     last_single_categories = zeros(eltype(timeline))

#     for i in eachindex(known)
#         known_categories = known[i]
#         # Only proceed if we have got to a single category, or the end
#         count(known_categories) == 1 || continue
#         current_categories = timeline[i]
#         if last_single_i < (i - 1)
#             # the present category as one of the possibilities
#             fillrange = last_single_i+1:i-1
#             # Replace the intermediate uncertain categories
#             for n in fillrange
#                 categories_at_n = possible[n]
#                 shared = map(&, last_single_categories, categories_at_n, current_categories)
#                 if any(shared)
#                     timeline[n] = shared
#                     if count(shared) == 1
#                         last_single_categories = shared
#                         last_single_i = n
#                     end
#                 else # Check possible direct transitions
#                     cats = merge_all(last_single_categories, categories_at_n, reversed, reversed_indirect)
#                     possible_direct = merge_all(current_categories, cats, transitions, indirect)
#                     if any(possible_direct)
#                         timeline[n] = possible_direct
#                     else # Check possible indirect transitions
#                         cats_indirect = merge_all(last_single_categories, categories_at_n, reversed_indirect, reversed)
#                         possible_indirect = merge_all(current_categories, cats_indirect, indirect, transitions)
#                         # Here we nuke impossible intermediate values
#                         timeline[n] = possible_indirect
#                     end
#                 end
#             end
#         end
#         last_single_i = i
#         last_single_categories = current_categories
#     end

#     # rev = lastindex(timeline):-1:firstindex(timeline)
#     # simplify_end!(timeline)
#     # simplify_end!(view(timeline, rev))

#     return timeline
# end

Base.@assume_effects :foldable function minimise_uncertainty!(
    timeline::AbstractVector{<:NamedVector{K}}, known::AbstractVector{<:NamedVector}, possible::AbstractVector{<:NamedVector}, transitions, reversed, indirect, reversed_indirect
) where K
    map(K) do k 
        for i in eachindex(timeline)
        end
    end
end

const LOG = []

# At the end assume nothing changes from the previous 
# step unless there are no shared categories
function simplify_end!(timeline)
    count(timeline[end]) > 1 || return timeline
    i = lastindex(timeline)
    while i > 0
        if count(timeline[i]) == 1
            # Check that every step until the end contains this category
            for j in i+1:lastindex(timeline)
                any(map(&, timeline[i], timeline[j])) || return timeline
            end
            # Simplify every step
            for j in i+1:lastindex(timeline)
                timeline[j] = timeline[i]
            end
            return timeline
        end
        i -= 1
    end
    return timeline
end

function cull_transitions!(
    timeline::AbstractArray{<:NamedVector{K}}, 
    known::AbstractArray{<:NamedVector{K}}, 
    possible::AbstractArray{<:NamedVector{K}}, 
    transitions::NamedVector{K}, 
    indirect::NamedVector{K},
) where K
    timeline[1] = a = map(&, first(known), first(possible))
    for i in firstindex(timeline)+1:lastindex(timeline)
        b = map(&, known[i], possible[i])
        a, match = merge_all(a, b, transitions, indirect)
        timeline[i] = a
    end
    return timeline
end

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
