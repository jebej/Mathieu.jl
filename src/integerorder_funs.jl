function cep(n::Order,q,z,N::Int=8+maximum(n)+ceil(Int,sqrt(abs(q))))
    An = cep_coeffs(checkvec!(n),float(q),N) # Calculate coefficients
    ce = zeros(length(z),length(n))
    for j = 1:length(n), k = 0:N-1 # Calculate function at z, for all n
        ce[:,j] .+= An[k+1,j].*cos.(2k.*z)
    end
    return ce
end

function cea(n::Order,q,z,N::Int=13+maximum(n)+ceil(Int,sqrt(abs(q))))
    An = cea_coeffs(checkvec!(n),float(q),N) # Calculate coefficients
    ce = zeros(length(z),length(n))
    for j = 1:length(n), k = 0:N-1 # Calculate function at z, for all n
        ce[:,j] .+= An[k+1,j].*cos.((2k+1).*z)
    end
    return ce
end

function sea(n::Order,q,z,N::Int=6+maximum(n)+ceil(Int,sqrt(abs(q))))
    Bn = sea_coeffs(checkvec!(n),float(q),N) # Calculate coefficients
    se = zeros(length(z),length(n))
    for j = 1:length(n), k = 0:N-1 # Calculate function at z, for all n
        se[:,j] .+= Bn[k+1,j].*sin.((2k+1).*z)
    end
    return se
end

function sep(n::Order,q,z,N::Int=7+maximum(n)+ceil(Int,sqrt(abs(q))))
    Bn = sep_coeffs(checkvec!(n),float(q),N) # Calculate coefficients
    se = zeros(length(z),length(n))
    for j = 1:length(n), k = 0:N-1 # Calculate function at z, for all n
        se[:,j] .+= Bn[k+1,j].*sin.((2k+2).*z)
    end
    return se
end

function fun_per(m::Order,q,z,N::Int=8+maximum(m)+ceil(Int,sqrt(abs(q))))
    oddorders = filter(isodd,m)
    isempty(oddorders)  && return cep(div.(m,2),q,z,N)
    evenorders = filter(iseven,m)
    isempty(evenorders) && return sep(div.(m,2),q,z,N)
    a = zeros(length(z),length(m))
    a[:,iseven.(m)] = cep(div.(evenorders,2),q,z,N)
    a[:,isodd.(m)]  = sep(div.(oddorders,2),q,z,N)
    return a
end
