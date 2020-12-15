signp(x::Real) = x<0 ? -one(x) : one(x)

checkvec(n::Number) = n:n
checkvec(n::AbstractVector) = n

blasfloat(x::Int32) = Float32(x)
blasfloat(x::Int16) = Float32(x)
blasfloat(x::Number) = Float64(x)

matsize1(n,q) = 2 + maximum(n) + ceil(Int,sqrt(abs(q)))
matsize2(m,q) = 2 + maximum(m)÷2 + ceil(Int,sqrt(abs(q)))
matsize3(n,q) = 2 + maximum(n) + ceil(Int,2*sqrt(abs(q)))
matsize4(m,q) = 2 + maximum(m)÷2 + ceil(Int,2*sqrt(abs(q)))

Base.filter(::typeof(iseven), r::UnitRange{T}) where {T<:Integer} =
    StepRange{T,T}((first(r)+T(1)) & ~T(1), T(2), last(r) & ~T(1))
Base.filter(::typeof(isodd), r::UnitRange{T}) where {T<:Integer} =
    StepRange{T,T}(first(r) & ~T(1) + T(1), T(2), (last(r)+T(1)) & ~T(1) - T(1))
