function cep_coeffs(n::Order,q::Number,N::Integer)
    iszero(q) && return q0_coeffs1(n,q,N)
    C1 = mat_C1(q,N) # Generate matrix
    Au = coeffs(C1,n)
    Au[1,:] .*= √0.5 # Renormalize eigenvectors
    return Au
end

function cea_coeffs(n::Order,q::Number,N::Integer)
    iszero(q) && return q0_coeffs2(n,q,N)
    C2 = mat_C2(q,N) # Generate matrix
    return coeffs(C2,n)
end

function sea_coeffs(n::Order,q::Number,N::Integer)
    iszero(q) && return q0_coeffs2(n,q,N)
    C3 = mat_C3(q,N) # Generate matrix
    return coeffs(C3,n)
end

function sep_coeffs(n::Order,q::Number,N::Integer)
    iszero(q) && return q0_coeffs2(n,q,N)
    C4 = mat_C4(q,N) # Generate matrix
    return coeffs(C4,n)
end

function coeffs(C::AbstractMatrix,n::Order)
    u = sort!(unique(n))
    au = eigvals(C)[u.+1] # Calculate characteristic values
    #au = LAPACK.stebz!('A','E',0.0,0.0,0,0,2*eps(C.ev[1]),C.dv,C.ev)[1][u+1]
    Au = eigvecs(C,au) # Calculate eigenvectors
    Au .*= signp.(sum(Au,1)) # Make sure we have the correct sign
    return Au[:,indexin(checkvec(n),u)]
end

function q0_coeffs1(n::Order,q::Number,N::Integer)
    A = zeros(typeof(blasfloat(q)),N,length(n))
    for (i,n) in enumerate(n); A[n+1,i] = n==0 ? 1/√(eltype(A)(2)) : 1; end
    return A
end

function q0_coeffs2(n::Order,q::Number,N::Integer)
    A = zeros(typeof(blasfloat(q)),N,length(n))
    for (i,n) in enumerate(n); A[n+1,i] = one(eltype(A)); end
    return A
end

function per_coeffs(m::Order,q::Number,N::Integer)
    oddind = findall(isodd,m)
    isempty(oddind)  && return cep_coeffs(div.(m,2),q,N)
    evenind = findall(iseven,m)
    isempty(evenind) && return sep_coeffs(div.(m,2),q,N)
    P1 = cep_coeffs(div.(m[evenind],2),q,N)
    P2 = sep_coeffs(div.(m[oddind],2),q,N)
    P = Matrix{eltype(P1)}(undef,N,length(m))
    P[:,evenind] = P1
    P[:,oddind] = P2
    return P
end
