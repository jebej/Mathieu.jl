function cep_coeffs(n::AbstractVector{Int},q,N::Int)
    C1 = mat_C1(q,N) # Generate matrix
    an = eigvals(C1)[n+1] # Calculate characteristic values
    An = eigvecs(C1,an) # Calculate eigenvectors
    An[2:N,:] .*= sqrt(2) # Renormalize eigenvectors
    An ./= sqrt.(sum(An.^2,1).+An[1:1,:].^2)
    return An
end

function cep_coeffs(n::Int,q,N::Int)
    C1 = mat_C1(q,N) # Generate matrix
    an = eigvals(C1)[n+1] # Calculate characteristic values
    An = eigvecs(C1,an) # Calculate eigenvectors
    An[2:N] .*= sqrt(2) # Renormalize eigenvectors
    An ./= sqrt(norm(An)^2+An[1]^2)
    return An
end

function cep(n::Int,q,z,N::Int=6+n+ceil(Int,sqrt(abs(q))))
    An = cep_coeffs(n:n,q,N) # Calculate coefficients
    ce = zero(z)
    kk = [2k+0.0 for k = 0:N-1]
    for i = 1:N # Calculate function at z
        ce += An[i].*cos.(kk[i].*z)
    end
    return ce
end

function cep(n::AbstractVector{Int},q,z,N::Int=6+maximum(n)+ceil(Int,sqrt(abs(q))))
    An = cep_coeffs(n,q,N) # Calculate coefficients
    ce = zeros(length(z),length(n))
    kk = [2k+0.0 for k = 0:N-1]
    for j = 1:length(n), i = 1:N # Calculate function at z, for all n
        ce[:,j] .+= An[i].*cos.(kk[i].*z)
    end
    return ce
end

function cep_prime(n::Int,q,z,N::Int=6+n+ceil(Int,sqrt(abs(q))))
    An = cep_coeffs(n:n,q,N) # Calculate coefficients
    ce = zero(z)
    kk = [2k+0.0 for k = 0:N-1]
    for i = 1:N # Calculate function at z
        ce -= An[i].*kk[i].*sin.(kk[i].*z)
    end
    return ce
end

function cep_prime(n::AbstractVector{Int},q,z,N::Int=6+maximum(n)+ceil(Int,sqrt(abs(q))))
    An = cep_coeffs(n,q,N) # Calculate coefficients
    ce = zeros(length(z),length(n))
    kk = [2k+0.0 for k = 0:N-1]
    for j = 1:length(n), i = 1:N # Calculate function at z, for all n
        ce[:,j] .-= An[i].*kk[i].*sin.(kk[i].*z)
    end
    return ce
end

function cea_coeffs(n::AbstractVector{Int},q,N::Int)
    C2 = mat_C2(q,N) # Generate matrix
    an = eigvals(C2)[n+1] # Calculate characteristic values
    An = eigvecs(C2,an) # Calculate eigenvectors
    return An
end

function cea(n::Int,q,z,N::Int=6+n+ceil(Int,sqrt(abs(q))))
    An = cea_coeffs(n:n,q,N) # Calculate coefficients
    ce = zero(z)
    kk = [2k+1.0 for k = 0:N-1]
    for i = 1:N # Calculate function at z
        ce += An[i].*cos.(kk[i].*z)
    end
    return ce
end

function cea(n::AbstractVector{Int},q,z,N::Int=6+maximum(n)+ceil(Int,sqrt(abs(q))))
    An = cea_coeffs(n,q,N) # Calculate coefficients
    ce = zeros(length(z),length(n))
    kk = [2k+1.0 for k = 0:N-1]
    for j = 1:length(n), i = 1:N # Calculate function at z, for all n
        ce[:,j] .+= An[i].*cos.(kk[i].*z)
    end
    return ce
end

function cea_prime(n::Int,q,z,N::Int=6+n+ceil(Int,sqrt(abs(q))))
    An = cea_coeffs(n:n,q,N) # Calculate coefficients
    ce = zero(z)
    kk = [2k+1.0 for k = 0:N-1]
    for i = 1:N # Calculate function at z
        ce -= An[i].*kk[i].*sin.(kk[i].*z)
    end
    return ce
end

function cea_prime(n::AbstractVector{Int},q,z,N::Int=6+maximum(n)+ceil(Int,sqrt(abs(q))))
    An = cea_coeffs(n,q,N) # Calculate coefficients
    ce = zeros(length(z),length(n))
    kk = [2k+1.0 for k = 0:N-1]
    for j = 1:length(n), i = 1:N # Calculate function at z, for all n
        ce[:,j] .-= An[i].*kk[i].*sin.(kk[i].*z)
    end
    return ce
end

function sea_coeffs(n::AbstractVector{Int},q,N::Int)
    C3 = mat_C3(q,N) # Generate matrix
    bn = eigvals(C3)[n+1] # Calculate characteristic values
    Bn = eigvecs(C3,bn) # Calculate eigenvectors
    return Bn
end

function sea(n::Int,q,z,N::Int=6+n+ceil(Int,sqrt(abs(q))))
    Bn = sea_coeffs(n:n,q,N) # Calculate coefficients
    se = zero(z)
    kk = [2k+1.0 for k = 0:N-1]
    for i = 1:N # Calculate function at z
        se += Bn[i].*sin.(kk[i].*z)
    end
    return se
end

function sea(n::AbstractVector{Int},q,z,N::Int=6+maximum(n)+ceil(Int,sqrt(abs(q))))
    Bn = sea_coeffs(n,q,N) # Calculate coefficients
    se = zeros(length(z),length(n))
    kk = [2k+1.0 for k = 0:N-1]
    for j = 1:length(n), i = 1:N # Calculate function at z, for all n
        se[:,j] .+= Bn[i].*sin.(kk[i].*z)
    end
    return ce
end

function sea_prime(n::Int,q,z,N::Int=6+n+ceil(Int,sqrt(abs(q))))
    Bn = sea_coeffs(n:n,q,N) # Calculate coefficients
    se = zero(z)
    kk = [2k+1.0 for k = 0:N-1]
    for i = 1:N # Calculate function at z
        se += Bn[i].*kk[i].*cos.(kk[i].*z)
    end
    return se
end

function sea_prime(n::AbstractVector{Int},q,z,N::Int=6+maximum(n)+ceil(Int,sqrt(abs(q))))
    Bn = sea_coeffs(n,q,N) # Calculate coefficients
    se = zeros(length(z),length(n))
    kk = [2k+1.0 for k = 0:N-1]
    for j = 1:length(n), i = 1:N # Calculate function at z, for all n
        se[:,j] .+= Bn[i].*kk[i].*cos.(kk[i].*z)
    end
    return se
end

function sep_coeffs(n::AbstractVector{Int},q,N::Int)
    C4 = mat_C4(q,N) # Generate matrix
    bn = eigvals(C4)[n+1] # Calculate characteristic values
    Bn = eigvecs(C4,bn) # Calculate eigenvectors
    return Bn
end

function sep(n::Int,q,z,N::Int=6+n+ceil(Int,sqrt(abs(q))))
    Bn = sep_coeffs(n:n,q,N) # Calculate coefficients
    se = zero(z)
    kk = [2k+2.0 for k = 0:N-1]
    for i = 1:N # Calculate function at z
        se += Bn[i].*sin.(kk[i].*z)
    end
    return se
end

function sep(n::AbstractVector{Int},q,z,N::Int=6+maximum(n)+ceil(Int,sqrt(abs(q))))
    Bn = sep_coeffs(n,q,N) # Calculate coefficients
    se = zeros(length(z),length(n))
    kk = [2k+2.0 for k = 0:N-1]
    for j = 1:length(n), i = 1:N # Calculate function at z, for all n
        se[:,j] .+= Bn[i].*sin.(kk[i].*z)
    end
    return se
end

function sep_prime(n::Int,q,z,N::Int=6+n+ceil(Int,sqrt(abs(q))))
    Bn = sep_coeffs(n:n,q,N) # Calculate coefficients
    se = zero(z)
    kk = [2k+2.0 for k = 0:N-1]
    for i = 1:N # Calculate function at z
        se += Bn[i].*kk[i].*cos.(kk[i].*z)
    end
    return se
end

function sep_prime(n::AbstractVector{Int},q,z,N::Int=6+maximum(n)+ceil(Int,sqrt(abs(q))))
    Bn = sep_coeffs(n,q,N) # Calculate coefficients
    se = zeros(length(z),length(n))
    kk = [2k+2.0 for k = 0:N-1]
    for j = 1:length(n), i = 1:N # Calculate function at z, for all n
        se[:,j] .+= Bn[i].*kk[i].*cos.(kk[i].*z)
    end
    return se
end
