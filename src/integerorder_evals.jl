"""
    mat_C1(q,N)

Generate the matrix corresponding to even, π-periodic solutions to the
angular Mathieu equation.

* `q`: real or complex parameter of Mathieu’s equation
* `N`: truncated matrix dimension (matrix is N×N)
"""
function mat_C1(q::Float64,N::Int)
    dd = [(2k+0.0)^2 for k = 0:N-1]
    du = fill(q,N-1)
    du[1] *= sqrt(2)
    return SymTridiagonal(dd,du)
end

"""
    mat_C2(q,N)

Generate the matrix corresponding to even, 2π-periodic (π-antiperiodic)
solutions to the angular Mathieu equation.

* `q`: real or complex parameter of Mathieu’s equation
* `N`: truncated matrix dimension (matrix is N×N)
"""
function mat_C2(q::Float64,N::Int)
    dd = [(2k+1.0)^2 for k = 0:N-1]
    du = fill(q,N-1)
    dd[1] += q
    return SymTridiagonal(dd,du)
end

"""
    mat_C3(q,N)

Generate the matrix corresponding to odd, 2π-periodic (π-antiperiodic)
solutions to the angular Mathieu equation.

* `q`: real or complex parameter of Mathieu’s equation
* `N`: truncated matrix dimension (matrix is N×N)
"""
function mat_C3(q::Float64,N::Int)
    dd = [(2k+1.0)^2 for k = 0:N-1]
    du = fill(q,N-1)
    dd[1] -= q
    return SymTridiagonal(dd,du)
end

"""
    mat_C4(q,N)

Generate the matrix corresponding to odd, π-periodic solutions to the
angular Mathieu equation.

* `q`: real or complex parameter of Mathieu’s equation
* `N`: truncated matrix dimension (matrix is N×N)
"""
function mat_C4(q::Float64,N::Int)
    dd = [(2k+2.0)^2 for k = 0:N-1]
    du = fill(q,N-1)
    return SymTridiagonal(dd,du)
end

"""
    ap(n,q,N)

Compute the characteristic value a_{2n}(q) corresponding to an even, π-periodic
solution to the angular Mathieu equation.
"""
function ap(n::Union{Int,AbstractVector{Int}},q,N::Int=8+maximum(n)+ceil(Int,sqrt(abs(q))))
    C = mat_C1(q,N)
    a = eigvals(C)[n+1]
    return a
end

"""
    aa(n,q,N)

Compute the characteristic value a_{2n+1}(q) corresponding to an even,
π-antiperiodic (2π-periodic) solution to the angular Mathieu equation.
"""
function aa(n::Union{Int,AbstractVector{Int}},q,N::Int=13+maximum(n)+ceil(Int,sqrt(abs(q))))
    C = mat_C2(q,N)
    a = eigvals(C)[n+1]
    return a
end

"""
    a(m,q,N)

Compute the characteristic value a_m(q). m must be 2n or 2n+1, with n=0,1,2...
"""
function a(m::Union{Int,AbstractVector{Int}},q)
    if isempty(filter(isodd,m))
        a = ap(div(m,2),q)
    elseif isempty(filter(iseven,m))
        a = aa(div(m,2),q)
    else
        a = zeros(length(m))
        a[iseven.(m)] = ap(div(filter(iseven,m),2),q)
        a[isodd.(m)] = aa(div(filter(isodd,m),2),q)
    end
    return a
end

"""
    ba(n,q,N)

Compute the characteristic value b_{2n+1}(q) corresponding to an odd,
π-antiperiodic (2π-periodic) solution to the angular Mathieu equation.
"""
function ba(n::Union{Int,AbstractVector{Int}},q,N::Int=6+maximum(n)+ceil(Int,sqrt(abs(q))))
    C = mat_C3(q,N)
    b = eigvals(C)[n+1]
    return b
end

"""
    bp(n,q,N)

Compute the characteristic value b_{2n+2}(q) corresponding to an odd, π-periodic
solution to the angular Mathieu equation.
"""
function bp(n::Union{Int,AbstractVector{Int}},q,N::Int=7+maximum(n)+ceil(Int,sqrt(abs(q))))
    C = mat_C4(q,N)
    b = eigvals(C)[n+1]
    return b
end

"""
    b(m,q,N)

Compute the characteristic value b_m(q). m must be 2n+1 or 2n+2, with n=0,1,2...
"""
function b(m::Union{Int,AbstractVector{Int}},q)
    if isempty(filter(isodd,m))
        a = bp(div(m,2)-2,q)
    elseif isempty(filter(iseven,m))
        a = ba(div(m,2)-1,q)
    else
        a = zeros(length(m))
        a[isodd.(m)] = ba(div(filter(isodd,m)-1,2),q)
        a[iseven.(m)] = bp(div(filter(iseven,m)-2,2),q)
    end
    return a
end

function abp(m::Union{Int,AbstractVector{Int}},q)
    if isempty(filter(isodd,m))
        a = ap(div(m,2),q)
    elseif isempty(filter(iseven,m))
        a = bp(div(m,2),q)
    else
        a = zeros(length(m))
        a[iseven.(m)] = ap(div(filter(iseven,m),2),q)
        a[isodd.(m)] = bp(div(filter(isodd,m),2),q)
    end
    return a
end
