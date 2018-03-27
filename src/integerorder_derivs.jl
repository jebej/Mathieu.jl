function cep_prime(n::Order,q,z,N::Int=matsize3(n,q))
    An = cep_coeffs(n,q,N) # Calculate coefficients
    ce = zeros(length(z),length(n))
    for j = 1:length(n), k = 0:N-1 # Calculate function at z, for all n
        ce[:,j] .-= An[k+1,j].*2k.*sin.(2k.*z)
    end
    return ce
end

function cea_prime(n::Order,q,z,N::Int=matsize3(n,q))
    An = cea_coeffs(n,q,N) # Calculate coefficients
    ce = zeros(length(z),length(n))
    for j = 1:length(n), k = 0:N-1 # Calculate function at z, for all n
        ce[:,j] .-= An[k+1,j].*(2k+1).*sin.((2k+1).*z)
    end
    return ce
end

function sea_prime(n::Order,q,z,N::Int=matsize3(n,q))
    Bn = sea_coeffs(n,q,N) # Calculate coefficients
    se = zeros(length(z),length(n))
    for j = 1:length(n), k = 0:N-1 # Calculate function at z, for all n
        se[:,j] .+= Bn[k+1,j].*(2k+1).*cos.((2k+1).*z)
    end
    return se
end

function sep_prime(n::Order,q,z,N::Int=matsize3(n,q))
    Bn = sep_coeffs(n,q,N) # Calculate coefficients
    se = zeros(length(z),length(n))
    for j = 1:length(n), k = 0:N-1 # Calculate function at z, for all n
        se[:,j] .+= Bn[k+1,j].*(2k+2).*cos.((2k+2).*z)
    end
    return se
end
