function int_cspm_cspn_prime(mm,nn,q,z,N::Int=8+maximum([mm;nn])+ceil(Int,sqrt(abs(q))))
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
                    res .-= Pm[k+1].*Pn[l+1].*cos(4(k+1).*z)./4
                else
                    res .-= Pm[k+1].*Pn[l+1].*0.5(l+1).*( cos(2(k-l).*z)./(k-l) .+ cos.(2(k+l+2).*z)./(k+l+2) )
                end
            end
        end
        return res[2]-res[1]
    end
end
