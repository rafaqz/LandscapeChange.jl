
"""
    PreallocatedUnorderedList <: AbstractVector

    PreallocatedUnorderedList(::Type{T}, maxlen::Int)
    PreallocatedUnorderedList(values::AbstractArray, len::Int)

A fast fixed-length unordered list. The maximum length is fixed 
and the order not specified so that `push!` does not need to allocate memory,
and `pop!` and `deleteat!` do not not remove memory - they just move values
from the end of the array to the specified indices and shorten the apparent
array length.
"""
mutable struct PreallocatedUnorderedList{T,A<:AbstractArray{T}} <: AbstractVector{T}
    values::A
    len::Int
end
PreallocatedUnorderedList(values::AbstractVector) = PreallocatedUnorderedList(values, length(values))
function PreallocatedUnorderedList(::Type{T}; maxlength::Int) where T
    values = Vector{T}(undef, maxlength)
    PreallocatedUnorderedList(values, 0)
end
Base.length(list::PreallocatedUnorderedList) = list.len
Base.size(list::PreallocatedUnorderedList) = (length(list),)
Base.@propagate_inbounds function Base.getindex(list::PreallocatedUnorderedList, i) 
    @boundscheck checkbounds(list, i) 
    @inbounds list.values[i]
end
Base.@propagate_inbounds function Base.setindex!(list::PreallocatedUnorderedList, x, i) 
    @boundscheck checkbounds(list, i)
    @inbounds list.values[i] = x
end
function Base.push!(list::PreallocatedUnorderedList, x)
    newlen = list.len + 1
    list.len = newlen
    list.values[newlen] = x 
end
function Base.pop!(list::PreallocatedUnorderedList)
    val = list.values[1]
    deleteat!(list, 1)
    return val
end
function Base.deleteat!(list::PreallocatedUnorderedList, i::Int)
    val = list.values[i]
    len = list.len
    list.values[i] = list.values[len]
    list.len = len - oneunit(len)
    return list
end

