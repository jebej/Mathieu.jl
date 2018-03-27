function cep_prime(n::Order,q,z,N::Int=matsize3(n,q))
    An = cep_coeffs(n,q,N) # Calculate coefficients
    # Evaluate sine for all z values
    C = [-2k*sin(2k*z) for z in z, k in 0:N-1]
    # Calculate function at z, for all n
    return C * An
end

function cea_prime(n::Order,q,z,N::Int=matsize3(n,q))
    An = cea_coeffs(n,q,N) # Calculate coefficients
    # Evaluate sine for all z values
    C = [-(2k+1)*sin((2k+1)*z) for z in z, k in 0:N-1]
    # Calculate function at z, for all n
    return C * An
end

function sea_prime(n::Order,q,z,N::Int=matsize3(n,q))
    Bn = sea_coeffs(n,q,N) # Calculate coefficients
    # Evaluate cosine for all z values
    S = [(2k+1)*cos((2k+1)*z) for z in z, k in 0:N-1]
    # Calculate function at z, for all n
    return S * Bn
end

function sep_prime(n::Order,q,z,N::Int=matsize3(n,q))
    Bn = sep_coeffs(n,q,N) # Calculate coefficients
    # Evaluate cosine for all z values
    S = [(2k+2)*cos((2k+2)*z) for z in z, k in 0:N-1]
    # Calculate function at z, for all n
    return S * Bn
end
