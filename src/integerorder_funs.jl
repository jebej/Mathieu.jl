function cep_coeffs(n::AbstractVector{Int},q,N::Int)
    u = sort!(unique(n))
    C1 = mat_C1(q,N) # Generate matrix
    au = eigvals(C1)[u+1] # Calculate characteristic values
    Au = eigvecs(C1,au) # Calculate eigenvectors
    Au[2:N,:] .*= sqrt(2) # Renormalize eigenvectors
    Au ./= sqrt.(sumabs2(Au,1).+Au[1:1,:].^2)
    return Au[:,indexin(n,u)]
end

function cep(n::Int,q,z,N::Int=8+n+ceil(Int,sqrt(abs(q))))
    An = cep_coeffs(n:n,q,N) # Calculate coefficients
    ce = zero(z)
    for k = 0:N-1 # Calculate function at z
        ce += An[k+1].*cos.(2k.*z)
    end
    return ce
end

function cep(n::AbstractVector{Int},q,z,N::Int=8+maximum(n)+ceil(Int,sqrt(abs(q))))
    An = cep_coeffs(n,q,N) # Calculate coefficients
    ce = zeros(length(z),length(n))
    for j = 1:length(n), k = 0:N-1 # Calculate function at z, for all n
        ce[:,j] .+= An[k+1,j].*cos.(2k.*z)
    end
    return ce
end

function cep_prime(n::Int,q,z,N::Int=8+n+ceil(Int,sqrt(abs(q))))
    An = cep_coeffs(n:n,q,N) # Calculate coefficients
    ce = zero(z)
    for k = 0:N-1 # Calculate function at z
        ce -= An[k+1].*2k.*sin.(2k.*z)
    end
    return ce
end

function cep_prime(n::AbstractVector{Int},q,z,N::Int=8+maximum(n)+ceil(Int,sqrt(abs(q))))
    An = cep_coeffs(n,q,N) # Calculate coefficients
    ce = zeros(length(z),length(n))
    for j = 1:length(n), k = 0:N-1 # Calculate function at z, for all n
        ce[:,j] .-= An[k+1,j].*2k.*sin.(2k.*z)
    end
    return ce
end

function cea_coeffs(n::AbstractVector{Int},q,N::Int)
    u = sort!(unique(n))
    C2 = mat_C2(q,N) # Generate matrix
    au = eigvals(C2)[u+1] # Calculate characteristic values
    Au = eigvecs(C2,au) # Calculate eigenvectors
    return Au[:,indexin(n,u)]
end

function cea(n::Int,q,z,N::Int=6+n+ceil(Int,sqrt(abs(q))))
    An = cea_coeffs(n:n,q,N) # Calculate coefficients
    ce = zero(z)
    for k = 0:N-1 # Calculate function at z
        ce += An[k+1].*cos.((2k+1).*z)
    end
    return ce
end

function cea(n::AbstractVector{Int},q,z,N::Int=13+maximum(n)+ceil(Int,sqrt(abs(q))))
    An = cea_coeffs(n,q,N) # Calculate coefficients
    ce = zeros(length(z),length(n))
    for j = 1:length(n), k = 0:N-1 # Calculate function at z, for all n
        ce[:,j] .+= An[k+1,j].*cos.((2k+1).*z)
    end
    return ce
end

function cea_prime(n::Int,q,z,N::Int=13+n+ceil(Int,sqrt(abs(q))))
    An = cea_coeffs(n:n,q,N) # Calculate coefficients
    ce = zero(z)
    for k = 0:N-1 # Calculate function at z
        ce -= An[k+1].*(2k+1).*sin.((2k+1).*z)
    end
    return ce
end

function cea_prime(n::AbstractVector{Int},q,z,N::Int=13+maximum(n)+ceil(Int,sqrt(abs(q))))
    An = cea_coeffs(n,q,N) # Calculate coefficients
    ce = zeros(length(z),length(n))
    for j = 1:length(n), k = 0:N-1 # Calculate function at z, for all n
        ce[:,j] .-= An[k+1,j].*(2k+1).*sin.((2k+1).*z)
    end
    return ce
end

function sea_coeffs(n::AbstractVector{Int},q,N::Int)
    u = sort!(unique(n))
    C3 = mat_C3(q,N) # Generate matrix
    bu = eigvals(C3)[u+1] # Calculate characteristic values
    Bu = eigvecs(C3,bu) # Calculate eigenvectors
    return Bu[:,indexin(n,u)]
end

function sea(n::Int,q,z,N::Int=6+n+ceil(Int,sqrt(abs(q))))
    Bn = sea_coeffs(n:n,q,N) # Calculate coefficients
    se = zero(z)
    for k = 0:N-1 # Calculate function at z
        se += Bn[k+1].*sin.((2k+1).*z)
    end
    return se
end

function sea(n::AbstractVector{Int},q,z,N::Int=6+maximum(n)+ceil(Int,sqrt(abs(q))))
    Bn = sea_coeffs(n,q,N) # Calculate coefficients
    se = zeros(length(z),length(n))
    for j = 1:length(n), k = 0:N-1 # Calculate function at z, for all n
        se[:,j] .+= Bn[k+1,j].*sin.((2k+1).*z)
    end
    return ce
end

function sea_prime(n::Int,q,z,N::Int=6+n+ceil(Int,sqrt(abs(q))))
    Bn = sea_coeffs(n:n,q,N) # Calculate coefficients
    se = zero(z)
    for k = 0:N-1 # Calculate function at z
        se += Bn[k+1].*(2k+1).*cos.((2k+1).*z)
    end
    return se
end

function sea_prime(n::AbstractVector{Int},q,z,N::Int=6+maximum(n)+ceil(Int,sqrt(abs(q))))
    Bn = sea_coeffs(n,q,N) # Calculate coefficients
    se = zeros(length(z),length(n))
    for j = 1:length(n), k = 0:N-1 # Calculate function at z, for all n
        se[:,j] .+= Bn[k+1,j].*(2k+1).*cos.((2k+1).*z)
    end
    return se
end

function sep_coeffs(n::AbstractVector{Int},q,N::Int)
    u = sort!(unique(n))
    C4 = mat_C4(q,N) # Generate matrix
    bu = eigvals(C4)[u+1] # Calculate characteristic values
    Bu = eigvecs(C4,bu) # Calculate eigenvectors
    return Bu[:,indexin(n,u)]
end

function sep(n::Int,q,z,N::Int=6+n+ceil(Int,sqrt(abs(q))))
    Bn = sep_coeffs(n:n,q,N) # Calculate coefficients
    se = zero(z)
    for k = 0:N-1 # Calculate function at z
        se += Bn[k+1].*sin.((2k+2).*z)
    end
    return se
end

function sep(n::AbstractVector{Int},q,z,N::Int=7+maximum(n)+ceil(Int,sqrt(abs(q))))
    Bn = sep_coeffs(n,q,N) # Calculate coefficients
    se = zeros(length(z),length(n))
    for j = 1:length(n), k = 0:N-1 # Calculate function at z, for all n
        se[:,j] .+= Bn[k+1,j].*sin.((2k+2).*z)
    end
    return se
end

function sep_prime(n::Int,q,z,N::Int=7+n+ceil(Int,sqrt(abs(q))))
    Bn = sep_coeffs(n:n,q,N) # Calculate coefficients
    se = zero(z)
    for k = 0:N-1 # Calculate function at z
        se += Bn[k+1].*(2k+2).*cos.((2k+2).*z)
    end
    return se
end

function sep_prime(n::AbstractVector{Int},q,z,N::Int=7+maximum(n)+ceil(Int,sqrt(abs(q))))
    Bn = sep_coeffs(n,q,N) # Calculate coefficients
    se = zeros(length(z),length(n))
    for j = 1:length(n), k = 0:N-1 # Calculate function at z, for all n
        se[:,j] .+= Bn[k+1,j].*(2k+2).*cos.((2k+2).*z)
    end
    return se
end


function csp_coeffs(m::AbstractVector{Int},q,N::Int)
    if isempty(filter(isodd,m))
        P = cep_coeffs(div(m,2),q,N)
    elseif isempty(filter(iseven,m))
        P = sep_coeffs(div(m,2),q,N)
    else
        P = zeros(N,length(m))
        P[:,iseven.(m)] = cep_coeffs(div(filter(iseven,m),2),q,N)
        P[:,isodd.(m)] = sep_coeffs(div(filter(isodd,m),2),q,N)
    end
    return P
end


function int_cspm_cspn_prime(mm,nn,q,z,N::Int=8+maximum([mm;nn])+ceil(Int,sqrt(abs(q))))
    P = csp_coeffs([mm;nn],q,N)
    map(product(1:length(mm),1:length(nn))) do x
        m = mm[x[1]]
        n = nn[x[2]]
        Pm = P[:,x[1]]
        Pn = P[:,length(mm)+x[2]]
        res = zero(z)
        if iseven(m)&iseven(n)
            for k=0:N-1, l=0:N-1
                if k==l
                    res += Pm[k+1].*Pn[l+1].*cos.(4k.*z)./4
                else
                    res += Pm[k+1].*Pn[l+1].*0.5l.*( cos.(2(k+l).*z)./(k+l) .- cos.(2(k-l).*z)./(k-l) )
                end
            end
        elseif iseven(m)&isodd(n)
            for k=0:N-1, l=0:N-1
                if k==l+1
                    res += Pm[k+1].*Pn[l+1].*( k.*z .+ sin.(4k.*z)./4 )
                else
                    res += Pm[k+1].*Pn[l+1].*0.5(l+1).*( (k+l+1).*sin.(2(k-l-1).*z) .+ (k-l-1).*sin.(2(k+l+1).*z) )./(k^2-l^2-2l-1)
                end
            end
        elseif isodd(m)&iseven(n)
            for k=0:N-1, l=0:N-1
                if k+1==l
                    res += Pm[k+1].*Pn[l+1].*( sin.(4(k+1).*z)./4 .- (k+1).*z )
                else
                    res += Pm[k+1].*Pn[l+1].*0.5l.*( sin.(2(k+l+1).*z)./(k+l+1) .- sin.(2(k-l+1).*z)./(k-l+1) )
                end
            end
        else
            for k=0:N-1, l=0:N-1
                if k==l
                    res -= Pm[k+1].*Pn[l+1].*cos(4(k+1).*z)./4
                else
                    res -= Pm[k+1].*Pn[l+1].*0.5(l+1).*( (k+l+2).*cos(2(k-l).*z) .+ (k-l).*cos.(2(k+l+2).*z) )./(k^2-l^2+2k-2l)
                end
            end
        end
        return res[2]-res[1]
    end
end
