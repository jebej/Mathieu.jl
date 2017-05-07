function cep_prime(n::Order,q,z,N::Int=8+maximum(n)+ceil(Int,sqrt(abs(q))))
    An = cep_coeffs(checkvec!(n),float(q),N) # Calculate coefficients
    ce = zeros(length(z),length(n))
    for j = 1:length(n), k = 0:N-1 # Calculate function at z, for all n
        ce[:,j] .-= An[k+1,j].*2k.*sin.(2k.*z)
    end
    return ce
end

function cea_prime(n::Order,q,z,N::Int=13+maximum(n)+ceil(Int,sqrt(abs(q))))
    An = cea_coeffs(checkvec!(n),float(q),N) # Calculate coefficients
    ce = zeros(length(z),length(n))
    for j = 1:length(n), k = 0:N-1 # Calculate function at z, for all n
        ce[:,j] .-= An[k+1,j].*(2k+1).*sin.((2k+1).*z)
    end
    return ce
end

function sea_prime(n::Order,q,z,N::Int=6+maximum(n)+ceil(Int,sqrt(abs(q))))
    Bn = sea_coeffs(checkvec!(n),float(q),N) # Calculate coefficients
    se = zeros(length(z),length(n))
    for j = 1:length(n), k = 0:N-1 # Calculate function at z, for all n
        se[:,j] .+= Bn[k+1,j].*(2k+1).*cos.((2k+1).*z)
    end
    return se
end

function sep_prime(n::Order,q,z,N::Int=7+maximum(n)+ceil(Int,sqrt(abs(q))))
    Bn = sep_coeffs(checkvec!(n),float(q),N) # Calculate coefficients
    se = zeros(length(z),length(n))
    for j = 1:length(n), k = 0:N-1 # Calculate function at z, for all n
        se[:,j] .+= Bn[k+1,j].*(2k+2).*cos.((2k+2).*z)
    end
    return se
end
