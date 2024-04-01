
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

# The first transitions is the "unit" transition - not going anywhere.
all_transitions(trans::NamedVector) = (generate_bitmasks(trans), trans, next_transitions(trans)...)
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
function stripe_raster(xs, states)
    map(xs) do raster
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

function color_raster(xs, states)
    map(xs) do raster
        color_raster(raster, states)
    end
end
function color_raster(raster::AbstractRaster, states; stripedim=X())
    A = rebuild(similar(raster, Float64); missingval=0.0)
    for I in DimIndices(A)
        v = raster[I...]
        c = count(v)
        x = if c == 0
            0.0
        else
            sum(v .* states ./ c)
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
cross_validate_timeline(ser::RasterSeries, transitions; kw...) =
    cross_validate_timeline(Rasters.combine(ser, Ti), transitions; kw...)
function cross_validate_timeline(source::Raster, transitions::NamedVector{K}) where K
    dest = copy(source)
    N = size(source, Ti())
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
end

function apply_transitions(timeline::Tuple, logic::NamedTuple) # (; transitions, reversed, reversed_indirect)
    state = _init_state(first(timeline))
    newtimeline, _ = _combine(timeline, logic, state, firstindex(timeline))
    return newtimeline
end

_init_state(x) = (; forced=zero(x), latent_forced=zero(x), uncertain=one(x), latent_uncertain=zero(x))
_maybefillmissing(t) = any(t) ? t : map(x -> true, t)

struct Forced end
struct Uncertain end

Base.@assume_effects :foldable function _combine(::Tuple{}, logic, from::L, i) where L
    return (), _init_state(first(from))
end
Base.@assume_effects :foldable function _combine(timeline::Tuple{T,Vararg}, logic, up::L, i) where {T,L}
    println("\n=== ", i, " ============================================================================")
    # Get the current timeline values
    t = _maybefillmissing(first(timeline))

    # Combine timeline and previous states up and down
    _, merged_up = _merge_to(up, t, logic.reversed)
    @show t up merged_up
    println("=== ", i, " ============================================================================\n")

    # Recursively move inwards in the timeline
    completed_timeline, down = _combine(Base.tail(timeline), logic, merged_up, i+1)

    println("\n---- ", i, " ----------------------------------------------------------------------")
    # Update the timeline with merged values from forwards/backwards passes
    _, merged_down = _merge_to(down, t, logic.transitions)

    final = _finalise(merged_down, merged_up)#, down, up)

    @show t up down merged_up merged_down final
    println("---- ", i, " ----------------------------------------------------------------------\n")

    new_timeline = (final, completed_timeline...)
    merged = (; merged_down[(:forced, :latent_forced)]..., uncertain=final, latent_uncertain=t)
    
    return new_timeline, merged
end

function _finalise(up, down)
    if any(up.forced)
        if any(down.forced)
            final = map(|, up.forced, down.forced)
        else
            final = if any(map(&, up.forced, down.uncertain))
                up.forced
            else
                map(|, up.forced, down.uncertain)
            end
        end
    else
        if any(down.forced)
            final = if any(map(&, up.uncertain, down.forced))
                down.forced
            else
                map(|, up.uncertain, down.forced)
            end
        else
            combined = map(&, up.uncertain, down.uncertain)
            final = if any(combined)
                combined
            else
                map(|, up.uncertain, down.uncertain)
            end
        end
    end
    # @show up down final
    return all(final) ? map(x -> false, final) : final
end

# Get source and dest modes
Base.@assume_effects :foldable function _merge_to(source::L, dest::T, transitions) where {L<:NamedTuple,T<:NamedVector}
    sourcemode, s = any(source.forced) ? (Forced(), source.forced) : (Uncertain(), source.uncertain)
    destmode = (count(dest) == 1 ? Forced() : Uncertain())
    final = _merge(sourcemode, destmode, s, dest, transitions)
    next_dest = if sourcemode isa Forced
        if destmode isa Forced
            (; forced=final, latent_forced=zero(final), uncertain=dest, latent_uncertain=source.uncertain)
        else
            # Here we remove unnessesary forcing of uncertain parameters of the 
            # same graph distance, by filtering them by which is also in the
            # passed in uncertain categories.
            shared_forced = map(&, final, source.forced)
            # If the forced state has changed
            next_dest = if !any(shared_forced) 
                shared_source_uncertain = map(&, final, source.uncertain)
                shared_dest_uncertain = map(&, dest, final)
                next_dest = if any(shared_source_uncertain)
                    final = uncertain = shared_source_uncertain
                    (; forced=zero(dest), latent_forced=source.forced, uncertain, latent_uncertain=source.uncertain)
                elseif any(shared_dest_uncertain)
                    final = uncertain = shared_dest_uncertain
                    (; forced=zero(final), latent_forced=source.forced, uncertain, latent_uncertain=source.uncertain)
                else
                    (; forced=zero(final), latent_forced=source.forced, uncertain=dest, latent_uncertain=source.uncertain)
                end
            else
                (; forced=final, latent_forced=zero(final), uncertain=dest, latent_uncertain=source.uncertain)
            end
            next_dest
        end
    else
        if destmode isa Forced
            (forced=final, latent_forced=zero(final), uncertain=dest, latent_uncertain=source.uncertain)
        else
            (forced=zero(final), latent_forced=zero(final), uncertain=final, latent_uncertain=source.uncertain)
        end
    end

    @show source dest final next_dest
    println()
    return final, next_dest
end

Base.@assume_effects :foldable function _merge(
    sourcemode, destmode, source::T, dest::T, transitions::Tuple
)::T where {T<:NamedVector{K}} where K
    # println("********************** $sourcemode $destmode **************************")
    # Just return identical source and dest
    source == dest && return dest

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
        # @show merged
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

