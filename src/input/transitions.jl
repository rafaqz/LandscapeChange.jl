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

-`assume_continuity`: assume continuity between single valued pixels
    accross time. If any intermediate pixels include that value
    as well as others, assume there were no changes. This may
    vary between a realistic assumtion and a very unrealistic
    assumtion, depending on categories and the length of time.
-`transitions`: force all transition logic, forwards and backwards,
    recording all transition errors to a separate raster.
"""
cross_validate_timeline(transitions, ser::RasterSeries; kw...) =
    cross_validate_timeline(transitions, Rasters.combine(ser, Ti); kw...)
function cross_validate_timeline(transitions::NamedVector{K}, raster::Raster{T}; kw...
) where {K,T}
    states = NamedVector{K}(ntuple(identity, length(K)))
    # Get time on the first dimension for performance
    timeline = copy(raster)#permutedims(raster, (Ti(), X(), Y()))
    error = fill!(similar(timeline), zero(T))
    pixel_timeline_copy1 = fill!(timeline[X(1), Y(1)], zero(T))
    pixel_timeline_copy2 = fill!(timeline[X(1), Y(1)], zero(T))
    for xy in DimIndices(dims(raster, (X(), Y())))
        pixel_timeline = view(timeline, xy...)
        pixel_error = view(error, xy...)
        pixel_timeline_copy1 .= pixel_timeline
        pixel_timeline_copy2 .= pixel_timeline
        cross_validate_timeline!(pixel_timeline, transitions;
            pixel_error, pixel_timeline_copy1, pixel_timeline_copy2, kw...
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
    pixel_error=copy(pixel_timeline),
    pixel_timeline_copy1=copy(pixel_timeline),
    pixel_timeline_copy2=copy(pixel_timeline),
    assume_continuity=true, fill_only=false,
)
    rev = lastindex(pixel_timeline):-1:1
    if fill_only
        _fill_empty_times!(pixel_timeline)
    else
        # Forward
        _apply_transitions!(pixel_timeline_copy1, pixel_error, reversed, reversed_indirect)
        # Backwards
        _apply_transitions!(view(pixel_timeline_copy2, rev), view(pixel_error, rev), transitions, indirect)
        if assume_continuity
            _remove_intermediate_uncertainty!(pixel_timeline_copy1, transitions, reversed, indirect, reversed_indirect)
            _remove_intermediate_uncertainty!(view(pixel_timeline_copy2, rev), reversed, transitions, reversed_indirect, indirect)
        end
        broadcast!(pixel_timeline, pixel_timeline_copy1, pixel_timeline_copy2) do a, b
            ands = map(&, a, b)
            if any(ands)
                ands
            else
                map(|, a, b)
            end
        end
    end
    if assume_continuity
        _remove_intermediate_uncertainty!(pixel_timeline, transitions, reversed, indirect, reversed_indirect)
        _remove_intermediate_uncertainty!(view(pixel_timeline_copy2, rev), reversed, transitions, reversed_indirect, indirect)
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

function _remove_intermediate_uncertainty!(
    timeline::AbstractVector, transitions, reversed, indirect, reversed_indirect
)
    past = first(timeline)
    last_single_i = typemax(Int)
    last_single_categories = zeros(eltype(timeline))
    for i in eachindex(timeline)
        current_categories = timeline[i]
        # prev_categories = timeline[i - 1]
        # next_categories = timeline[i + 1]
        # @show current_categories last_single_i
        # Only proceed if we have got to a single category, or the end
        if !(count(current_categories) == 1 || (i == lastindex(timeline)))
            continue
        end
        if last_single_i < (i - 1)
            # the present category as one of the possibilities
            fillrange = last_single_i+1:i-1
            # Replace the intermediate uncertain categories
            for n in fillrange
                categories_at_n = timeline[n]
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
                        cats_indirect  = _merge_all_possible(last_single_categories, categories_at_n, reversed_indirect)
                        possible_indirect = _merge_all_possible(current_categories, cats_indirect, indirect)
                        timeline[n] = possible_indirect
                    end
                end
            end
        end
        last_single_i = i
        last_single_categories = current_categories
    end
end


function _merge_all_possible(source, dest, transitions, indirect)
    direct = _merge_all_possible(source, dest, transitions)
    if any(direct)
        return direct
    else
        return _merge_all_possible(source, dest, indirect)
    end
end
function _merge_all_possible(source, dest, transitions)
    reduce(zip(source, transitions); init=zero(source)) do acc, (s, possible_dest)
        if s
            map(|, map(&, dest, possible_dest), acc)
        else
            acc
        end
    end
end

function _apply_transitions!(timeline, errors, transitions, indirect)
    a = first(timeline)
    for i in firstindex(timeline)+1:lastindex(timeline)
        b = timeline[i]
        matching_b = _merge_all_possible(a, b, transitions, indirect)
        a, error = if any(matching_b)
            matching_b, zero(a)
        else
            if any(a)
                if any(b)
                    # Nothing matches and `a` is certain, so
                    # propagate `a` and return `b` in error
                    if count(a) == 1
                        a, b
                    else
                        # Never force-propagate uncertain states
                        map(|, a, b), map(|, a, b)
                    end
                else
                    # b is empty, just propagate `a`
                    a, zero(b)
                end
            else
                # `a` is empty, just propagate `b`
                b, zero(a)
            end
        end
        timeline[i] = a
        errors[i] = errors[i] .| error
    end
end
