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
    oddorders = filter(isodd,m)
    isempty(oddorders)  && return cep_coeffs(div.(m,2),q,N)
    evenorders = filter(iseven,m)
    isempty(evenorders) && return sep_coeffs(div.(m,2),q,N)
    P = zeros(N,length(m))
    P[:,iseven.(m)] = cep_coeffs(div.(evenorders,2),q,N)
    P[:,isodd.(m)] = sep_coeffs(div.(oddorders,2),q,N)
    return P
end


# Not working
#function per_coeffs2(m::AbstractVector{Int},q,N::Int)
#    #if q==0.0; return q0_coeffs1(m,q,N); end
#    u = sort!(unique(m))
#    MP = mat_per(q,N) # Generate matrix
#    au = eigvals(MP)[u+1] # Calculate characteristic values
#    Au = eigvecs(MP,au) # Calculate eigenvectors
#    Au./= sqrt(2) # Renormalize eigenvectors
#    Au .*= signp.(sum(Au,1))
#    #Au[div(N,2)+2:end,:] .*= sign.(sum(Au[div(N,2)+2:end,:],1))*(-1).^(u+1)
#    Au[div(N,2)+2:end,:] .+= Au[div(N,2):-1:1,:].*((-1).^(u)).'
#    Au = Au[div(N,2)+1:end,:]
#    return Au[:,indexin(m,u)]
#end
#
#function aper_coeffs2(m::AbstractVector{Int},q,N::Int)
#    #if q==0.0; return q0_coeffs2(m,q,N); end
#    u = sort!(unique(m))
#    MA = mat_aper(q,N) # Generate matrix
#    au = eigvals(MA)[u+1] # Calculate characteristic values
#    Au = eigvecs(MA,au) # Calculate eigenvectors
#    Au .*= signp.(sum(Au,1))
#    return Au[:,indexin(m,u)]
#end
