
"""
    compile_timeline(masks, categories; kw...)

Generate 
"""
function compile_timeline(masks, categories;
    path="/home/raf/PhD/Mascarenes",
    files = define_map_files(; path),
    category_names=nothing,
)
    println("Generating raster slices...")
    return _compile_timeline(masks, categories, files; category_names)
end

function _compile_timeline(masks::NamedTuple, categories, file_list::NamedTuple; category_names)
    # Load all rasters and make masks layers for all land-cover classes.
    # We get the numbers form the saved ".json" file for the project.
    map(file_list, masks[keys(file_list)]) do file_list, mask
        _compile_timeline(mask, categories, file_list; category_names)
    end
end
function _compile_timeline(mask::AbstractArray, categories, file_list::NamedTuple{Keys}; category_names) where Keys
    forced = []
    files = map(file_list, NamedTuple{Keys}(Keys)) do (image_path, data), key
        raster_path = splitext(image_path)[1] * ".tif"
        raw = rebuild(mask .* _fix_order(Raster(raster_path)); missingval=0)
        if data isa Pair{<:AbstractArray}
            original_names = data[1]
            categories = data[2]
        elseif data isa Pair{<:Function}
            original_names = data[1](image_path)
            categories = data[2]
        else
            throw(ArgumentError("Error at file key $key: $(data[1]) is not a valid specification value"))
        end
        times = _get_times(categories)
        specification = _as_namedtuples(categories, category_names, times)
        # Gapfill NamedTuples with nothing
        grouped = map(specification) do (time, categories)
            _category_raster(raw, original_names, categories, mask, forced, time)
        end
        (; raw, grouped, times, specification, original_names)
    end
    alltimes = sort!(union(map(f -> f.times, files)...))
    timeline_dict = Dict{Int,Any}()
    for file in files
        for time in alltimes
            i = findfirst(==(time), file.times)
            if !isnothing(i)
                if haskey(timeline_dict, time)
                    timeline_dict[time] = map(.|, timeline_dict[time], file.grouped[i])
                else
                    timeline_dict[time] = file.grouped[i]
                end
            end
        end
    end
    timeline_pairs = sort!(collect(pairs(timeline_dict)); by=first)
    stacks = map(identity, RasterStack.(last.(timeline_pairs)))
    if eltype(stacks) == Any
        return nothing
    else
        years = first.(timeline_pairs)
        time = Ti(years; sampling=Intervals(End()), span=Irregular(1500, last(years)))
        timeline = RasterSeries(stacks, time)
        for (time, (cat_key, forced_rast)) in forced
            st = timeline[At(time)]
            # Wipe all stack layers where forced raster is true
            map(Rasters.DimensionalData.layers(st)) do layer
                broadcast!(layer, layer, forced_rast) do l, f
                    f ? false : l
                end
            end
            # Add forced values to its category
            cat_rast = getproperty(st, cat_key)
            broadcast!(cat_rast, cat_rast, forced_rast) do c, f
                f ? true : c
            end
        end

        return (; files, timeline)
    end
end

function _category_raster(raster::Raster, layer_names::Vector, layers::NamedTuple{Keys}, mask, forced, time) where Keys
    map(layers, NamedTuple{Keys}(Keys)) do layer, key
        _category_raster(raster, layer_names, layer, mask, forced, time, key)
    end
end
function _category_raster(raster::Raster, layer_names::Vector, layer_components::Vector, mask, forced, time, key)::Raster{Bool}
    layers = map(l -> _category_raster(raster, layer_names, l, mask, forced, time, key), layer_components)
    out = rebuild(Bool.(broadcast(|, layers...)); missingval=false) .& mask
    @assert missingval(out) == false
    return out
end
function _category_raster(raster::Raster, layer_names::Vector, xs::Tuple{<:Function,Vararg}, mask, forced, time, key)::Raster{Bool}
    f, args... = xs
    vals = map(args) do layer
        _category_raster(raster, layer_names, layer, mask, forced, time, key)
    end
    return map(f, vals...) .& mask
end
function _category_raster(raster::Raster, layer_names::Vector, category::Symbol, mask, forced, time, key)::Raster{Bool}
    if category === :mask
        return mask
    else
        error(":$category directive not understood")
    end
end
function _category_raster(raster::Raster, layer_names::Vector, category::Nothing, mask, forced, time, key)::Raster{Bool}
    return map(_ -> false, raster)
end
function _category_raster(raster::Raster, layer_names::Vector, category::String, mask, forced, time, key)::Raster{Bool}
    I = findall(==(category), map(String, layer_names))
    if length(I) == 0
        error("could not find $category in $(layer_names)")
    end
    # Get all values matching the first category as a mask
    out = rebuild(Bool.(raster .== first(I)); missingval=false)
    # Add pixels for any subsequent categories
    foreach(I[2:end]) do i
        out .|= raster .== first(i)
    end
    @assert eltype(out) == Bool
    @assert missingval(out) == false
    return out .& mask
end
function _category_raster(raster::Raster, layer_names::Vector, x::Pair, mask, forced, time, key)
    x[1] == :force || error("$(x[1]) not recognised")
    cr = _category_raster(raster, layer_names, x[2], mask, forced, time, key)
    push!(forced, time => key => cr)
    return cr
end
function _category_raster(raster::Raster, layer_names::Vector, x, mask, forced, time, key)
    error("slice must be a NamedTuple, String or Vector{String}, got a $(typeof(x)) - $x")
end

function _get_times(categories)
    times = Set{Int}()
    for cat in categories
        if cat isa Pair{Int}
            push!(times, cat[1])
        elseif cat isa Vector
            foreach(cat) do (time, val)
                push!(times, time)
            end
        end
    end
    return sort!(collect(times))
end

function _as_namedtuples(categories::NamedTuple, category_names, times)
    map(times) do time
        time => map(category_names) do k
            haskey(categories, k) || return nothing
            category = categories[k]
            if category isa Pair
                if last(category) isa Symbol
                    return nothing
                elseif first(category) == time
                    return last(category)
                else
                    return nothing
                end
            elseif isnothing(category)
                return nothing
            else
                i = findfirst(c -> first(c) == time, category)
                if isnothing(i)
                    return nothing
                else
                    return last(category[i])
                end
            end
        end
    end
end

_fix_order(A) = reorder(A, ForwardOrdered)
