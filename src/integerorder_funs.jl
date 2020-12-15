"""
    ce(m,q,z)

Compute the cosine-elliptic Mathieu function \$\\text{ce}_m(q,z)\$. `m` is the order of the solution and must be an integer \$m=0,1,2,3…\$.

Pass a vector, e.g. `m=0:4` and/or `z=LinRange(0,π,101)`, to compute the function for multiple orders and variables efficiently.
"""
function ce(m::Order,q::Number,z,N::Integer=matsize4(m,q))
    ie = map(iseven,m)
    all(ie)  && return cep(div.(m,2),q,z,N) # all even
    !any(ie) && return cea(div.(m,2),q,z,N) # all odd
    ce1 = cep(div.(filter(iseven,m),2),q,z,N)
    ce2 = cea(div.(filter(isodd, m),2),q,z,N)
    ce = similar(ce1,length(z),length(m))
    ce[:,ie] = ce1
    ce[:,.!ie] = ce2
    return ce
end

"""
    se(m,q,z)

Compute the sine-elliptic Mathieu function \$\\text{se}_m(q,z)\$. `m` is the order of the solution and must be an integer, excluding 0, \$m=1,2,3,4…\$.

Pass a vector, e.g. `m=1:5` and/or `z=LinRange(0,π,101)`, to compute the function for multiple orders and variables efficiently.
"""
function se(m::Order,q::Number,z,N::Integer=matsize4(m,q))
    ie = map(iseven,m)
    all(ie)  && return sep(div.(m.-2,2),q,z,N) # all even
    !any(ie) && return sea(div.(m.-1,2),q,z,N) # all odd
    se1 = sep(div.(filter(iseven,m).-2,2),q,z,N)
    se2 = sea(div.(filter(isodd, m).-1,2),q,z,N)
    se = similar(se1,length(z),length(m))
    se[:,ie] = se1
    se[:,.!ie] = se2
    return se
end

cep(n::Order,q::Number,z::Number,N::Integer=matsize3(n,q)) = cep(n,q,z:z,N)[1]
function cep(n::Order,q::Number,z::AbstractVector,N::Integer=matsize3(n,q))
    An = cep_coeffs(n,q,N) # Calculate coefficients
    # Evaluate cosine for all z values
    C = [cos(2k*z) for z in z, k in 0:N-1]
    # Calculate function at z, for all n
    return C * An
end

cea(n::Order,q::Number,z::Number,N::Integer=matsize3(n,q)) = cea(n,q,z:z,N)[1]
function cea(n::Order,q::Number,z::AbstractVector,N::Integer=matsize3(n,q))
    An = cea_coeffs(n,q,N) # Calculate coefficients
    # Evaluate cosine for all z values
    C = [cos((2k+1)*z) for z in z, k in 0:N-1]
    # Calculate function at z, for all n
    return C * An
end

sea(n::Order,q::Number,z::Number,N::Integer=matsize3(n,q)) = sea(n,q,z:z,N)[1]
function sea(n::Order,q::Number,z::AbstractVector,N::Integer=matsize3(n,q))
    Bn = sea_coeffs(n,q,N) # Calculate coefficients
    # Evaluate sine for all z values
    S = [sin((2k+1)*z) for z in z, k in 0:N-1]
    # Calculate function at z, for all n
    return S * Bn
end

sep(n::Order,q::Number,z::Number,N::Integer=matsize3(n,q)) = sep(n,q,z:z,N)[1]
function sep(n::Order,q::Number,z::AbstractVector,N::Integer=matsize3(n,q))
    Bn = sep_coeffs(n,q,N) # Calculate coefficients
    # Evaluate sine for all z values
    S = [sin((2k+2)*z) for z in z, k in 0:N-1]
    # Calculate function at z, for all n
    return S * Bn
end

fun_per(m::Order,q::Number,z::Number,N::Integer=matsize4(m,q)) = fun_per(m,q,z:z,N)[1]
function fun_per(m::Order,q::Number,z::AbstractVector,N::Integer=matsize4(m,q))
    ie = map(iseven,m)
    all(ie)  && return cep(div.(m,2),q,z,N)
    !any(ie) && return sep(div.(m,2),q,z,N)
    fp1 = cep(div.(filter(iseven,m),2),q,z,N)
    fp2 = sep(div.(filter(isodd, m),2),q,z,N)
    fp = similar(fp1,length(z),length(m))
    fp[:,ie] = fp1
    fp[:,.!ie] = fp2
    return fp
end
