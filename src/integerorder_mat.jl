"""
    mat_C1(q,N)

Generate the matrix corresponding to even, π-periodic solutions to the
angular Mathieu equation.

* `q`: real or complex parameter of Mathieu’s equation
* `N`: truncated matrix dimension (matrix is N×N)
"""
function mat_C1(q::T,N::Integer) where T<:Number
    dd = [T(2k)^2 for k = 0:N-1]
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
function mat_C2(q::T,N::Integer) where T<:Number
    dd = [T(2k+1)^2 for k = 0:N-1]
    dd[1] += q
    du = fill(q,N-1)
    return SymTridiagonal(dd,du)
end

"""
    mat_C3(q,N)

Generate the matrix corresponding to odd, 2π-periodic (π-antiperiodic)
solutions to the angular Mathieu equation.

* `q`: real or complex parameter of Mathieu’s equation
* `N`: truncated matrix dimension (matrix is N×N)
"""
function mat_C3(q::T,N::Integer) where T<:Number
    dd = [T(2k+1)^2 for k = 0:N-1]
    dd[1] -= q
    du = fill(q,N-1)
    return SymTridiagonal(dd,du)
end

"""
    mat_C4(q,N)

Generate the matrix corresponding to odd, π-periodic solutions to the
angular Mathieu equation.

* `q`: real or complex parameter of Mathieu’s equation
* `N`: truncated matrix dimension (matrix is N×N)
"""
function mat_C4(q::T,N::Integer) where T<:Number
    dd = [T(2k+2)^2 for k = 0:N-1]
    du = fill(q,N-1)
    return SymTridiagonal(dd,du)
end

"""
    mat_per(q,N)

Generate the matrix corresponding to both even and odd, π-periodic solutions to
the angular Mathieu equation.

* `q`: real or complex parameter of Mathieu’s equation
* `N`: truncated matrix dimension (matrix is N×N)
"""
function mat_per(q::T,N::Integer) where T<:Number
    dd = [T(2k)^2 for k = -div(N,2):div(N-1,2)]
    du = fill(q,N-1)
    return SymTridiagonal(dd,du)
end

"""
    mat_aper(q,N)

Generate the matrix corresponding to both even and odd, π-antiperiodic solutions
to the angular Mathieu equation.

* `q`: real or complex parameter of Mathieu’s equation
* `N`: truncated matrix dimension (matrix is N×N)
"""
function mat_aper(q::T,N::Integer) where T<:Number
    dd = [T(2k+1)^2 for k = -div(N,2):div(N-1,2)]
    du = fill(q,N-1)
    return SymTridiagonal(dd,du)
end



#function mat_even(q::T,N::Integer) where T<:Number
#    dd = [T(k)^2 for k = 0:N-1]
#    dd[2] += q/1.3
#    du = fill(q,N-1)
#    #du[1] *= sqrt(2)
#    return SymTridiagonal(dd,du)
#end
