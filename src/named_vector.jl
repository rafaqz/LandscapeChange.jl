struct NamedVector{K,L,T,V} <: FieldVector{L,T}
    val::V
    NamedVector{K,L,T,V}(val::V) where {K,L,T,V} = new{K,L,T,V}(val)
    NamedVector{K,L,T,V}(val::V) where {K,L,T,V<:StaticArray} = new{K,L,T,V}(val)
    NamedVector{K,L,T,V}(val::V) where {K,L,T,V<:Tuple} = new{K,L,T,V}(val)
end
NamedVector{K,L,T,V}(x1::T, x2::T, xs::T...) where {K,L,T,V} = NamedVector{K,L,T,V}((x1, x2, xs...))

NamedVector{K,L,T}(t::V) where {K,L,T,V<:Tuple} = NamedVector{K,L,T,V}(t)
NamedVector{K,L,T}(t::V) where {K,L,T,V<:NamedTuple{K,V1}} where V1 = NamedVector{K,L,T,V1}(t)
NamedVector(val::NamedTuple{K,V}) where {K,V<:Tuple{Vararg{Any,L}}} where L = NamedVector{K,L}(val)
NamedVector(val::NamedVector) = val
NamedVector{K}(t::V) where {K,V<:Tuple{Vararg{<:Any,L}}} where L = NamedVector{K,L}(t)
# Type promotion step...
function NamedVector{K,L}(t::V) where {K,L,V<:Tuple{Vararg{<:Any,L}}}
    t_uniform_type = map(x -> convert(promote_type(map(typeof, t)...), x), t)
    T = promote_type(map(typeof, t)...)
    NamedVector{K,L,T,NTuple{L,T}}(t_uniform_type)
end
function NamedVector{K,L}(val::NT) where {K,L,NT<:NamedTuple{K,<:Tuple{Vararg{<:Any,L}}}}
    NamedVector{K,L}(values(val))
end
function NamedVector{K,L}(v::V) where {V <: StaticVector{L,T}} where {K,L,T}
    NamedVector{K,L,T,V}(v)
end
NamedVector{K}(v::StaticVector{L,T}) where {K,L,T} = NamedVector{K,L}(v)
NamedVector{K}(v::AbstractVector{T}) where {K,T} = NamedVector{K}(StaticVector{length(K),T}(v))
NamedVector(; kw...) = NamedVector{keys(kw),length(kw)}(values(kw)) 


Base.@assume_effects :foldable function Base.getproperty(nv::NamedVector{names}, k::Symbol) where names
    idx = findfirst(n -> k == n, names)
    if idx === nothing
        _noidx_error(nv, k)
    end
    return parent(nv)[idx]
end

@noinline _noidx_error(::T, k) where T = error("type $T has no field $(k)")

Base.NamedTuple(a::NamedVector{K}) where K = NamedTuple{K}(Tuple(a)) 
# Base.pairs(a::NamedVector) = Base.pairs(NamedTuple(a))
Base.zero(a::NamedVector) = map(x -> zero(x), a)
# Base.zero(a::Type{<:NamedVector{K,L,T}}) where {K,L,T} = NamedVector{K}(ntuple(zero(T), L))   
Base.one(a::NamedVector) = map(x -> one(x), a)
Base.oneunit(a::NamedVector) = map(x -> oneunit(x), a)
Base.parent(a::NamedVector) = getfield(a, :val)
Base.values(a::NamedVector) = values(parent(a))
# Base.getproperty(a::NamedVector{K}, x::Symbol) where K = getproperty(NamedTuple{K}(parent(a)), x)
Base.propertynames(a::NamedVector{K}) where K = K
Base.@propagate_inbounds Base.getindex(a::NamedVector, i::Int) = getindex(parent(a), i)
Base.@propagate_inbounds Base.setindex!(a::NamedVector, x, i::Int) = (setindex!(parent(a), i, x); a)
Base.@propagate_inbounds Base.setindex(a::NamedVector, x, i::Int) = (setindex(parent(a), i, x); a)
Base.getindex(a::NamedVector, x::Symbol) = getproperty(a, x)
Base.getindex(a::NamedVector, x::Tuple) = NamedVector(getindex(NamedTuple(a), x))
Base.convert(::Type{NamedVector}, nt::NamedTuple) = NamedVector(nt)
# This is slow and apparently not needed
# function Base.map(f, a1::NamedVector{K,L}, as::NamedVector...) where {K,L}
#     map(as) do a
#         K == keys(a) || _throw_nv_mismatch(K, keys(a))
#     end
#     NamedVector{K}(map(f, map(SVector{L}, (a1, as...))...))
# end
Base.show(io::IO, a::NamedVector) = Base.show(io, parent(a))
Base.merge(a::NamedVector) = a

@noinline _throw_nv_mismatch(K1, K2) = throw(ArgumentError("NamedVector `map`: keys $K1 not equal to keys $K2"))


function DynamicGrids.ConstructionBase.constructorof(::Type{<:NamedVector{K,L,T}}) where {K,L,T}
    NamedVector{K,L,T}
end

# Base.merge(a::NamedVector, b::NamedVector, args...) =
    # NamedVector(merge(merge(parent(a), parent(b)), args...)
Base.merge(a::NamedVector, b::NamedTuple, args...) =
    NamedVector(merge(merge(parent(a), b), args...))
Base.merge(a::NamedTuple, b::NamedVector, args...) =
    merge(merge(a, parent(b)), args...)

StaticArrays.similar_type(::Type{A}, ::Type{T}, S::Size) where {T,A<:NamedVector} =
    _namedvector_similar_size(A, T, S, Size(A))

StaticArrays.Size(::Type{<:NamedVector{K}}) where K = StaticArrays.Size(length(K))

function _namedvector_similar_size(
    A::Type{<:NamedVector{K,L}}, ::Type{T}, NewSize::S, OldSize::S
) where {K,L,T,S}
    return NamedVector{K,L,T}
end
function _namedvector_similar_size(A, T::Type, NewSize, OldSize)
    StaticArrays.default_similar_type(T, NewSize, StaticArrays.length_val(NewSize))
end

_maybeparent(a::NamedVector) = parent(a)
_maybeparent(nt::NamedTuple) = nt

Base.tail(a::NamedVector{K,L,T}) where {K,L,T} = NamedVector(Base.tail(parent(a)))
