antisign(x::Real) = ifelse(x<0, one(x), oftype(x,-1))
checkvec(n::Number) = [n]
checkvec(n::AbstractVector) = n

function q0_coeffs1(n,q,N)
    A = zeros(eltype(q),N,length(n))
    for n in n; A[n+1,n+1] = n==0?1/sqrt(2):1.0; end
    return A
end

function q0_coeffs2(n,q,N)
    A = zeros(eltype(q),N,length(n))
    for n in n; A[n+1,n+1] = 1.0; end
    return A
end
