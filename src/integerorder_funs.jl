function cep(n::Order,q,z,N::Int=8+maximum(n)+ceil(Int,sqrt(abs(q))))
    An = cep_coeffs(checkvec(n),float(q),N) # Calculate coefficients
    ce = zeros(length(z),length(n))
    for j = 1:length(n), k = 0:N-1 # Calculate function at z, for all n
        ce[:,j] .+= An[k+1,j].*cos.(2k.*z)
    end
    return ce
end

function cea(n::Order,q,z,N::Int=13+maximum(n)+ceil(Int,sqrt(abs(q))))
    An = cea_coeffs(checkvec(n),float(q),N) # Calculate coefficients
    ce = zeros(length(z),length(n))
    for j = 1:length(n), k = 0:N-1 # Calculate function at z, for all n
        ce[:,j] .+= An[k+1,j].*cos.((2k+1).*z)
    end
    return ce
end

function sea(n::Order,q,z,N::Int=6+maximum(n)+ceil(Int,sqrt(abs(q))))
    Bn = sea_coeffs(checkvec(n),float(q),N) # Calculate coefficients
    se = zeros(length(z),length(n))
    for j = 1:length(n), k = 0:N-1 # Calculate function at z, for all n
        se[:,j] .+= Bn[k+1,j].*sin.((2k+1).*z)
    end
    return se
end

function sep(n::Order,q,z,N::Int=7+maximum(n)+ceil(Int,sqrt(abs(q))))
    Bn = sep_coeffs(checkvec(n),float(q),N) # Calculate coefficients
    se = zeros(length(z),length(n))
    for j = 1:length(n), k = 0:N-1 # Calculate function at z, for all n
        se[:,j] .+= Bn[k+1,j].*sin.((2k+2).*z)
    end
    return se
end

function fun_per(m::Order,q,z,N::Int=8+maximum(m)+ceil(Int,sqrt(abs(q))))
    oddind = find(isodd,m)
    isempty(oddind)  && return cep(div.(m,2),q,z,N)
    evenind = find(iseven,m)
    isempty(evenind) && return sep(div.(m,2),q,z,N)
    a = zeros(length(z),length(m))
    a[:,evenind] = cep(div.(m[evenind],2),q,z,N)
    a[:,oddind]  = sep(div.(m[oddind],2),q,z,N)
    return a
end
