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
    a_p(n,q,N)

Compute the characteristic value a_{2n}(q) corresponding to an even, π-periodic
solution to the angular Mathieu equation.
"""
function a_p(n::Union{Int,AbstractVector{Int}},q::Float64,N::Int=6+maximum(n)+ceil(Int,sqrt(q)))
    C = mat_C1(q,N)
    a = eigvals(C)[n+1]
    return a
end

"""
    a_a(n,q,N)

Compute the characteristic value a_{2n+1}(q) corresponding to an even,
π-antiperiodic (2π-periodic) solution to the angular Mathieu equation.
"""
function a_a(n::Union{Int,AbstractVector{Int}},q::Float64,N::Int=6+maximum(n)+ceil(Int,sqrt(q)))
    C = mat_C2(q,N)
    a = eigvals(C)[n+1]
    return a
end

"""
    b_a(n,q,N)

Compute the characteristic value b_{2n+1}(q) corresponding to an odd,
π-antiperiodic (2π-periodic) solution to the angular Mathieu equation.
"""
function b_a(n::Union{Int,AbstractVector{Int}},q::Float64,N::Int=6+maximum(n)+ceil(Int,sqrt(q)))
    C = mat_C3(q,N)
    b = eigvals(C)[n+1]
    return b
end

"""
    b_p(n,q,N)

Compute the characteristic value b_{2n+2}(q) corresponding to an odd, π-periodic
solution to the angular Mathieu equation.
"""
function b_p(n::Union{Int,AbstractVector{Int}},q::Float64,N::Int=6+maximum(n)+ceil(Int,sqrt(q)))
    C = mat_C4(q,N)
    b = eigvals(C)[n+1]
    return b
end

function ce_p(n::Union{Int,AbstractVector{Int}},q::Float64,z::Real,N::Int=6+maximum(n)+ceil(Int,sqrt(q)))
    C1 = mat_C1(q,N) # Initialize matrix
    an = eigvals(C1)[n+1:n+1] # Calculate characteristic value
    An = eigvecs(C1,an) # Calculate eigenvector
    An[2:end] *= sqrt(2) # Renormalize eigenvector
    An /= sqrt(norm(An)^2+An[1]^2)
    # Calculate function at z
    res = zero(z)
    kk = [2k+zero(z) for k = 0:N]
    for i = 1:N+1
        res += vn[i]*cos(kk[i]*z)
    end
    return res
end

function ce_a(n::Union{Int,AbstractVector{Int}},q::Float64,z::Real,N=n+10)
    C2 = mat_C2(q,N) # Initialize matrix
    an = eigvals(C2)[n+1:n+1] # Calculate characteristic value
    An = eigvecs(C2,an) # Calculate eigenvector
    # Calculate function at z
    res = zero(z)
    kk = [2k+one(z) for k = 0:N]
    for i = 1:N+1
        res += An[i]*cos(kk[i]*z)
    end
    return res
end
