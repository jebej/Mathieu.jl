cep_prime(n::Order,q::Number,z::Number,N::Integer=matsize3(n,q)) = cep_prime(n,q,z:z,N)[1]
function cep_prime(n::Order,q::Number,z::AbstractVector,N::Integer=matsize3(n,q))
    An = cep_coeffs(n,q,N) # Calculate coefficients
    # Evaluate sine for all z values
    C = [-2k*sin(2k*z) for z in z, k in 0:N-1]
    # Calculate function at z, for all n
    return C * An
end

cea_prime(n::Order,q::Number,z::Number,N::Integer=matsize3(n,q)) = cea_prime(n,q,z:z,N)[1]
function cea_prime(n::Order,q::Number,z::AbstractVector,N::Integer=matsize3(n,q))
    An = cea_coeffs(n,q,N) # Calculate coefficients
    # Evaluate sine for all z values
    C = [-(2k+1)*sin((2k+1)*z) for z in z, k in 0:N-1]
    # Calculate function at z, for all n
    return C * An
end

sea_prime(n::Order,q::Number,z::Number,N::Integer=matsize3(n,q)) = sea_prime(n,q,z:z,N)[1]
function sea_prime(n::Order,q::Number,z::AbstractVector,N::Integer=matsize3(n,q))
    Bn = sea_coeffs(n,q,N) # Calculate coefficients
    # Evaluate cosine for all z values
    S = [(2k+1)*cos((2k+1)*z) for z in z, k in 0:N-1]
    # Calculate function at z, for all n
    return S * Bn
end

sep_prime(n::Order,q::Number,z::Number,N::Integer=matsize3(n,q)) = sep_prime(n,q,z:z,N)[1]
function sep_prime(n::Order,q::Number,z::AbstractVector,N::Integer=matsize3(n,q))
    Bn = sep_coeffs(n,q,N) # Calculate coefficients
    # Evaluate cosine for all z values
    S = [(2k+2)*cos((2k+2)*z) for z in z, k in 0:N-1]
    # Calculate function at z, for all n
    return S * Bn
end
