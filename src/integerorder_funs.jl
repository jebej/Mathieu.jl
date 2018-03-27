function cep(n::Order,q,z,N::Int=matsize3(n,q))
    An = cep_coeffs(n,q,N) # Calculate coefficients
    # Evaluate cosine for all z values
    C = [cos(2k*z) for z in z, k in 0:N-1]
    # Calculate function at z, for all n
    return C * An
end

function cea(n::Order,q,z,N::Int=matsize3(n,q))
    An = cea_coeffs(n,q,N) # Calculate coefficients
    # Evaluate cosine for all z values
    C = [cos((2k+1)*z) for z in z, k in 0:N-1]
    # Calculate function at z, for all n
    return C * An
end

function sea(n::Order,q,z,N::Int=matsize3(n,q))
    Bn = sea_coeffs(n,q,N) # Calculate coefficients
    # Evaluate sine for all z values
    S = [sin((2k+1)*z) for z in z, k in 0:N-1]
    # Calculate function at z, for all n
    return S * Bn
end

function sep(n::Order,q,z,N::Int=matsize3(n,q))
    Bn = sep_coeffs(n,q,N) # Calculate coefficients
    # Evaluate sine for all z values
    S = [sin((2k+2)*z) for z in z, k in 0:N-1]
    # Calculate function at z, for all n
    return S * Bn
end

function fun_per(m::Order,q,z,N::Int=matsize3(n,q))
    oddind = find(isodd,m)
    isempty(oddind)  && return cep(div.(m,2),q,z,N)
    evenind = find(iseven,m)
    isempty(evenind) && return sep(div.(m,2),q,z,N)
    a = zeros(length(z),length(m))
    a[:,evenind] = cep(div.(m[evenind],2),q,z,N)
    a[:,oddind]  = sep(div.(m[oddind],2),q,z,N)
    return a
end
