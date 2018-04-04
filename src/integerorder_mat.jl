"""
    mat_C1(q,N)

Generate the matrix corresponding to even, π-periodic solutions to the
angular Mathieu equation.

* `q`: real parameter of Mathieu’s equation
* `N`: truncated matrix dimension (matrix is N×N)
"""
mat_C1(q::Number,N::Integer) = mat_C1(blasfloat(q),N)
function mat_C1(q::T,N::Integer) where T<:BlasReal
    dd = [T(2k)^2 for k = 0:N-1]
    du = fill(q,N-1)
    du[1] *= √2
    return SymTridiagonal(dd,du)
end

"""
    mat_C2(q,N)

Generate the matrix corresponding to even, 2π-periodic (π-antiperiodic)
solutions to the angular Mathieu equation.

* `q`: real parameter of Mathieu’s equation
* `N`: truncated matrix dimension (matrix is N×N)
"""
mat_C2(q::Number,N::Integer) = mat_C2(blasfloat(q),N)
function mat_C2(q::T,N::Integer) where T<:BlasReal
    dd = [T(2k+1)^2 for k = 0:N-1]
    dd[1] += q
    du = fill(q,N-1)
    return SymTridiagonal(dd,du)
end

"""
    mat_C3(q,N)

Generate the matrix corresponding to odd, 2π-periodic (π-antiperiodic)
solutions to the angular Mathieu equation.

* `q`: real parameter of Mathieu’s equation
* `N`: truncated matrix dimension (matrix is N×N)
"""
mat_C3(q::Number,N::Integer) = mat_C3(blasfloat(q),N)
function mat_C3(q::T,N::Integer) where T<:BlasReal
    dd = [T(2k+1)^2 for k = 0:N-1]
    dd[1] -= q
    du = fill(q,N-1)
    return SymTridiagonal(dd,du)
end

"""
    mat_C4(q,N)

Generate the matrix corresponding to odd, π-periodic solutions to the
angular Mathieu equation.

* `q`: real parameter of Mathieu’s equation
* `N`: truncated matrix dimension (matrix is N×N)
"""
mat_C4(q::Number,N::Integer) = mat_C4(blasfloat(q),N)
function mat_C4(q::T,N::Integer) where T<:BlasReal
    dd = [T(2k+2)^2 for k = 0:N-1]
    du = fill(q,N-1)
    return SymTridiagonal(dd,du)
end
