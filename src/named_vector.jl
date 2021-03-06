

struct NamedVector{K,L,T,V} <: FieldVector{L,T}
    nt::NamedTuple{K,V}
    NamedVector{K,L,T,V}(nt::NT) where {K,L,T,V<:Tuple,NT<:NamedTuple{K,Tuple{Vararg{T,L}}}} = new{K,L,T,V}(nt)
end
NamedVector{K,L,T,V}(t::Tuple{Vararg{T,L}}) where {K,L,T,V<:Tuple} = NamedVector{K,L,T,V}(NamedTuple{K}(t))
NamedVector{K,L,T,V}(x::T, xs::T...) where {K,L,T,V} = NamedVector{K,L,T,V}((x, xs...))

NamedVector{K,L,T}(t::V) where {K,L,T,V<:Tuple} = NamedVector{K,L,T,V}(t)
NamedVector{K,L,T}(t::V) where {K,L,T,V<:NamedTuple{K,V1}} where V1 = NamedVector{K,L,T,V1}(t)
NamedVector(nt::NamedTuple{K,V}) where {K,V<:Tuple{Vararg{Any,L}}} where L = NamedVector{K,L}(nt)
NamedVector{K}(t::V) where {K,V<:Tuple{Vararg{<:Any,L}}} where L = NamedVector{K,L}(t)
# Type promotion step...
function NamedVector{K,L}(t::V) where {K,L,V<:Tuple{Vararg{<:Any,L}}}
    t_uniform_type = map(x -> convert(promote_type(map(typeof, t)...), x), t)
    T = promote_type(map(typeof, t)...)
    NamedVector{K,L,T,NTuple{L,T}}(NamedTuple{K,NTuple{L,T}}(t_uniform_type))
end
function NamedVector{K,L}(nt::NT) where {K,L,NT<:NamedTuple{K,<:Tuple{Vararg{<:Any,L}}}}
    # TODO: why doesn't this compile away...
    NamedVector{K,L}(Tuple(nt))
end
function NamedVector{K,L}(v::SVector{L,T}) where {K,L,T}
    NamedVector{K,L,T,NTuple{L,T}}(NamedTuple{K,NTuple{L,T}}(Tuple(v)))
end
NamedVector(; kw...) = NamedVector{keys(kw),length(kw)}(values(kw)) 

Base.parent(a::NamedVector) = getfield(a, :nt)
Base.values(a::NamedVector) = values(parent(a))
Base.getproperty(a::NamedVector, x::Symbol) = getproperty(parent(a), x)
Base.@propagate_inbounds Base.getindex(a::NamedVector, i::Int) = getindex(parent(a), i)
Base.@propagate_inbounds Base.setindex!(a::NamedVector, x, i::Int) = (setindex!(parent(a), i, x); a)
# Base.map(f, as::Union{NamedVector{K},NamedTuple{K}}...) where K = begin
Base.map(f, a1::NamedVector{K,L}, as::NamedVector{K}...) where {K,L} = begin
    NamedVector{K,L}(map(f, map(SVector, (a1, as...))...))
end

_maybeparent(a::NamedVector) = parent(a)
_maybeparent(nt::NamedTuple) = nt

Base.tail(a::NamedVector{K,L,T}) where {K,L,T} = NamedVector(Base.tail(parent(a)))
