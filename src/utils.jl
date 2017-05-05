signp(x::Real) = ifelse(x<0, oftype(x,-1), one(x))
checkvec(n::Number) = [n]
checkvec(n::AbstractVector) = n

function q0_coeffs1(n,q,N)
    A = zeros(eltype(q),N,length(n))
    for (i,n) in enumerate(n)
        A[n+1,i] = n==0 ? 1/sqrt(2) : 1
    end
    return A
end

function q0_coeffs2(n,q,N)
    A = zeros(eltype(q),N,length(n))
    for (i,n) in enumerate(n)
        A[n+1,i] = 1
    end
    return A
end
