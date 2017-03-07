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
