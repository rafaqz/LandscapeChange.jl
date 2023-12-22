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
    rev = lastindex(pixel_timeline):-1:1
    if fill_only
        _fill_empty_times!(pixel_timeline)
    else
        # Forward
        _apply_transitions!(forward, reversed, reversed_indirect)
        # Backwards
        _apply_transitions!(view(backward, rev), transitions, indirect)

        if simplify
            # _remove_intermediate_uncertainty!(forward, pixel_timeline, possible, transitions, reversed, indirect, reversed_indirect)
            # _remove_intermediate_uncertainty!(view(backward, rev), view(pixel_timeline, rev), view(possible, rev), reversed, transitions, reversed_indirect, indirect)
            _remove_intermediate_uncertainty!(forward)
            _remove_intermediate_uncertainty!(view(backward, rev))
        end

        for i in eachindex(possible)
            f, b = forward[i], backward[i]
            x = map(&, f, b)
            possible[i] = any(x) ? x : map(|, f, b)
        end

        if cull
            # Forward
            _apply_transitions!(forward, pixel_timeline, possible, reversed, reversed_indirect)
            # Backwards
            _apply_transitions!(view(backward, rev), view(pixel_timeline, rev), view(possible, rev), transitions, indirect)

            for i in eachindex(possible)
                f, b = forward[i], backward[i]
                x = map(&, f, b)
                possible[i] = any(x) ? x : map(|, f, b)
            end
        end
        if any(force)
            for i in eachindex(possible)
                x = map(&, force, possible[i])
                if any(x)
                    possible[i] = x
                end
            end
        end


        # Finalise output
        for i in eachindex(pixel_timeline)
            before = pixel_timeline[i]
            after = possible[i]
            # Any shared state is enough to say there was no error
            shared = map(&, before, after)
            # Write errors
            if !any(shared)
                pixel_error[i] = before
            end
            # Write timeline
            pixel_timeline[i] = after
        end
    end

    return pixel_timeline, pixel_error
end

# Fill empty slices with previous data
function _fill_empty_times!(timeline::AbstractVector)
    last_non_empty_i = typemax(Int)
    last_non_empty_categories = zero(eltype(timeline))
    for i in eachindex(timeline)
        present_categories = timeline[i]
        any(present_categories) || continue # No category data for this time slice
        if last_non_empty_i < (i - 1)
            # There has been a gap in data, fill it with the combination of
            # the last non empty categories and the present category
            fill_categories = map(|, last_non_empty_categories, present_categories)
            for n in last_non_empty_i+1:i-1
                timeline[n] = fill_categories
            end
        end
        last_non_empty_i = i
        last_non_empty_categories = present_categories
    end
end

function _remove_intermediate_uncertainty!(timeline)
    bitmasks = _bitmasks(first(timeline))
    singles = zero(first(timeline))
    lastsingle_i = firstindex(timeline)
    for i in eachindex(timeline)
        x = timeline[i]
        if count(x) == 1
            singles = map(|, singles, x)
            # Only keep states that where alone somewhere
            for i in lastsingle_i+1:i-1
                timeline[i] = timeline[i] .& singles
            end
            lastsingle_i = i
        end
    end
    _simplify_end!(timeline)
    return timeline
end

function _remove_intermediate_uncertainty!(
    dest::AbstractVector, known, possible, transitions, reversed, indirect, reversed_indirect
)
    past = first(known)
    last_single_i = typemax(Int)
    last_single_categories = zeros(eltype(dest))
    for i in eachindex(known)
        known_categories = known[i]
        # Only proceed if we have got to a single category, or the end
        count(known_categories) == 1 || continue
        if last_single_i < (i - 1)
            # the present category as one of the possibilities
            fillrange = last_single_i+1:i-1
            for n in fillrange
                possible[n]
            end
            # Replace the intermediate uncertain categories
            for n in fillrange
                categories_at_n = possible[n]
                shared = map(&, last_single_categories, categories_at_n, current_categories)
                if any(shared)
                    timeline[n] = shared
                    if count(shared) == 1
                        last_single_categories = shared
                        last_single_i = n
                    end
                else # Check possible direct transitions
                    cats = _merge_all_possible(last_single_categories, categories_at_n, reversed)
                    possible = _merge_all_possible(current_categories, cats, transitions)
                    if any(possible)
                        timeline[n] = possible
                    else # Check possible indirect transitions
                        cats_indirect = _merge_all_possible(last_single_categories, categories_at_n, reversed_indirect)
                        possible_indirect = _merge_all_possible(current_categories, cats_indirect, indirect)
                        # Here we nuke impossible intermediate values
                        timeline[n] = possible_indirect
                    end
                end
            end
        end
        last_single_i = i
        last_single_categories = current_categories
    end

    # _simplify_end!(timeline)

    return timeline
end

# At the end assume nothing changes from the previous 
# step unless there are no shared categories
function _simplify_end!(timeline)
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

function _apply_transitions!(
    timeline::AbstractArray{<:NamedVector{K}}, 
    transitions::NamedVector{K}, indirect::NamedVector{K},
) where K
    a = first(timeline)
    for i in firstindex(timeline)+1:lastindex(timeline)
        b = timeline[i]
        a = _merge_all(a, b, transitions, indirect)
        timeline[i] = a
    end
    return timeline
end

function _apply_transitions!(
    timeline::AbstractArray{<:NamedVector{K}}, 
    known::AbstractArray{<:NamedVector{K}}, 
    possible::AbstractArray{<:NamedVector{K}}, 
    transitions::NamedVector{K}, indirect::NamedVector{K},
) where K
    timeline[1] = a = map(&, first(known), first(possible))
    for i in firstindex(timeline)+1:lastindex(timeline)
        b = map(&, known[i], possible[i])
        a = _merge_all(a, b, transitions, indirect)
        timeline[i] = a
    end
    return timeline
end

Base.@assume_effects :foldable function _merge_all(
    source::NamedVector{K}, dest::NamedVector{K}, 
    transitions::NamedVector{K}, indirect::NamedVector{K},
) where K
    any(source) || return dest
    bitmasks = _bitmasks(source)
    reduce(zip(source, transitions, indirect, bitmasks); init=zero(source)) do acc, (s, direct_dest, indirect_dest, bitmask)
        x = if s
            dir = map(&, dest, direct_dest)
            xs = if any(dir) 
                dir
            else
                ind = map(&, dest, indirect_dest)
                if any(ind) 
                    ind
                else
                    # Broken logic: keep this category
                    bitmask
                end
            end
            map(|, xs, acc)
        else
            acc
        end
        return x
    end
end

# Base.@assume_effects :foldable function _merge_all_possible(
#     source::NamedVector{K}, dest::NamedVector{K}, transitions::NamedVector{K},
# ) where K
#     any(source) || return dest
#     reduce(zip(source, transitions); init=zero(source)) do acc, (s, direct_dest)
#         x = if s
#             map(|, map(&, dest, direct_dest), acc)
#         else
#             acc
#         end
#         return x
#     end
# end

@generated function _bitmasks(::NamedVector{K}) where K
    nvs = ntuple(length(K)) do i
        vals = falses(length(K))
        vals[i] = true 
        :(NamedVector{K}($(Tuple(vals))))
    end
    expr = Expr(:tuple, nvs...)
    return :(NamedVector{K}($expr))
end
