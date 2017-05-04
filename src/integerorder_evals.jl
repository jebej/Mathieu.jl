"""
    ap(n,q,N)

Compute the characteristic value a_{2n}(q) corresponding to an even, π-periodic
solution to the angular Mathieu equation.
"""
function ap(n::Order,q::BlasReal,N::Int=8+maximum(n)+ceil(Int,sqrt(abs(q))))
    C = mat_C1(q,N)
    a = eigvals(C)[n+1]
    return a
end

"""
    aa(n,q,N)

Compute the characteristic value a_{2n+1}(q) corresponding to an even,
π-antiperiodic (2π-periodic) solution to the angular Mathieu equation.
"""
function aa(n::Order,q::BlasReal,N::Int=13+maximum(n)+ceil(Int,sqrt(abs(q))))
    C = mat_C2(q,N)
    a = eigvals(C)[n+1]
    return a
end

"""
    a(m,q,N)

Compute the characteristic value a_m(q). m must be 2n or 2n+1, with n=0,1,2...
"""
function a(m::Order,q)
    if isempty(filter(isodd,m))
        a = ap(div(m,2),q)
    elseif isempty(filter(iseven,m))
        a = aa(div(m,2),q)
    else
        a = zeros(length(m))
        a[iseven.(m)] = ap(div(filter(iseven,m),2),q)
        a[isodd.(m)] = aa(div(filter(isodd,m),2),q)
    end
    return a
end

"""
    ba(n,q,N)

Compute the characteristic value b_{2n+1}(q) corresponding to an odd,
π-antiperiodic (2π-periodic) solution to the angular Mathieu equation.
"""
function ba(n::Order,q::BlasReal,N::Int=6+maximum(n)+ceil(Int,sqrt(abs(q))))
    C = mat_C3(q,N)
    b = eigvals(C)[n+1]
    return b
end

"""
    bp(n,q,N)

Compute the characteristic value b_{2n+2}(q) corresponding to an odd, π-periodic
solution to the angular Mathieu equation.
"""
function bp(n::Order,q::BlasReal,N::Int=7+maximum(n)+ceil(Int,sqrt(abs(q))))
    C = mat_C4(q,N)
    b = eigvals(C)[n+1]
    return b
end

"""
    b(m,q,N)

Compute the characteristic value b_m(q). m must be 2n+1 or 2n+2, with n=0,1,2...
"""
function b(m::Order,q)
    if isempty(filter(isodd,m))
        a = bp(div(m,2)-2,q)
    elseif isempty(filter(iseven,m))
        a = ba(div(m,2)-1,q)
    else
        a = zeros(length(m))
        a[isodd.(m)] = ba(div(filter(isodd,m)-1,2),q)
        a[iseven.(m)] = bp(div(filter(iseven,m)-2,2),q)
    end
    return a
end

"""
    char_per(m,q,N)

Compute the characteristic value a_{m}(q) for even m, or b_{2m}(q) for odd m.
These eigenvalues characteristic values to even and odd, respectively,
π-periodic solutions to the angular Mathieu equation.
"""
function char_per(n::Order,q::BlasReal,N::Int=9+maximum(n)+ceil(Int,sqrt(abs(q))))
    M = mat_per(q,N)
    v = eigvals(M)[n+1]
    return v
end

"""
    char_aper(m,q,N)

Compute the characteristic value a_{m+1}(q) for even m, or b_{2m+1}(q) for odd
m. These characteristic values correspond to even and odd, respectively,
π-antiperiodic solutions to the angular Mathieu equation.
"""
function char_aper(n::Order,q::BlasReal,N::Int=8+maximum(n)+ceil(Int,sqrt(abs(q))))
    M = mat_aper(q,N)
    v = eigvals(M)[n+1]
    return v
end

#function abp(m::Order,q)
#    if isempty(filter(isodd,m))
#        a = ap(div(m,2),q)
#    elseif isempty(filter(iseven,m))
#        a = bp(div(m,2),q)
#    else
#        a = zeros(length(m))
#        a[iseven.(m)] = ap(div(filter(iseven,m),2),q)
#        a[isodd.(m)] = bp(div(filter(isodd,m),2),q)
#    end
#    return a
#end
