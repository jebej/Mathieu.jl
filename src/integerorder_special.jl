function int_p_pprime(mm,nn,q,z,N::Int=8+maximum([mm;nn])+ceil(Int,sqrt(abs(q))))
    P = per_coeffs([mm;nn],q,N)
    map(product(1:length(mm),1:length(nn))) do x
        m = mm[x[1]]
        n = nn[x[2]]
        Pm = P[:,x[1]]
        Pn = P[:,length(mm)+x[2]]
        res = zero(z)
        if iseven(m)&&iseven(n)
            for k=0:N-1, l=0:N-1
                if k==l
                    res .+= Pm[k+1].*Pn[l+1].*cos.(4k.*z)./4
                else
                    res .+= Pm[k+1].*Pn[l+1].*0.5l.*( cos.(2(k+l).*z)./(k+l) .- cos.(2(k-l).*z)./(k-l) )
                end
            end
        elseif iseven(m)&&isodd(n)
            for k=0:N-1, l=0:N-1
                if k==l+1
                    res .+= Pm[k+1].*Pn[l+1].*( k.*z .+ sin.(4k.*z)./4 )
                else
                    res .+= Pm[k+1].*Pn[l+1].*0.5(l+1).*( sin.(2(k-l-1).*z)./(k-l-1) .+ sin.(2(k+l+1).*z)./(k+l+1) )
                end
            end
        elseif isodd(m)&&iseven(n)
            for k=0:N-1, l=0:N-1
                if k+1==l
                    res .+= Pm[k+1].*Pn[l+1].*( sin.(4(k+1).*z)./4 .- (k+1).*z )
                else
                    res .+= Pm[k+1].*Pn[l+1].*0.5l.*( sin.(2(k+l+1).*z)./(k+l+1) .- sin.(2(k-l+1).*z)./(k-l+1) )
                end
            end
        else
            for k=0:N-1, l=0:N-1
                if k==l
                    res .-= Pm[k+1].*Pn[l+1].*cos.(4(k+1).*z)./4
                else
                    res .-= Pm[k+1].*Pn[l+1].*0.5(l+1).*( cos.(2(k-l).*z)./(k-l) .+ cos.(2(k+l+2).*z)./(k+l+2) )
                end
            end
        end
        return (res[2]-res[1]) #(-1)^(1+m÷2+n÷2)*
    end
end

# optimized version of the above when z = [0,π]
function intpi_p_pprime(nmax::Int,q,N::Int=8+nmax+ceil(Int,sqrt(abs(q))))
    P = per_coeffs(0:nmax-1,q,N)
    M = zeros(eltype(P),nmax,nmax)
    @inbounds for j = 2:2:nmax, i = 1:2:nmax
        el = zero(eltype(P))
        for k=1:N-1
            el += P[k+1,i]*P[k,j]*k*π
        end
        el = (-1)^(i÷2+j÷2)*el # hmmmm
        M[i,j] = el; M[j,i] = -el
    end
    return M
end
