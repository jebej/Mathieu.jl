function cep_coeffs(n::Order,q,N::Int)
    if q==0.0; return q0_coeffs1(n,q,N); end
    u = sort!(unique(n))
    C1 = mat_C1(q,N) # Generate matrix
    au = eigvals(C1)[u+1] # Calculate characteristic values
    Au = eigvecs(C1,au) # Calculate eigenvectors
    Au[2:N,:] .*= sqrt(2) # Renormalize eigenvectors
    Au ./= sqrt.(sum(abs2,Au,1).+Au[1:1,:].^2)
    Au .*= signp.(sum(Au,1))
    return Au[:,indexin(n,u)]
end

function cea_coeffs(n::Order,q,N::Int)
    if q==0.0; return q0_coeffs2(n,q,N); end
    u = sort!(unique(n))
    C2 = mat_C2(q,N) # Generate matrix
    au = eigvals(C2)[u+1] # Calculate characteristic values
    Au = eigvecs(C2,au) # Calculate eigenvectors
    Au .*= signp.(sum(Au,1))
    return Au[:,indexin(n,u)]
end

function sea_coeffs(n::Order,q,N::Int)
    if q==0.0; return q0_coeffs2(n,q,N); end
    u = sort!(unique(n))
    C3 = mat_C3(q,N) # Generate matrix
    bu = eigvals(C3)[u+1] # Calculate characteristic values
    Bu = eigvecs(C3,bu) # Calculate eigenvectors
    Bu .*= signp.(sum(Bu,1))
    return Bu[:,indexin(n,u)]
end

function sep_coeffs(n::Order,q,N::Int)
    if q==0.0; return q0_coeffs2(n,q,N); end
    u = sort!(unique(n))
    C4 = mat_C4(q,N) # Generate matrix
    bu = eigvals(C4)[u+1] # Calculate characteristic values
    Bu = eigvecs(C4,bu) # Calculate eigenvectors
    Bu .*= signp.(sum(Bu,1))
    return Bu[:,indexin(n,u)]
end

function per_coeffs(m::AbstractVector{Int},q,N::Int)
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
