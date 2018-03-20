"""
    ap(n,q,N)

Compute the characteristic value \$a_{2n}(q)\$ corresponding to an even, π-periodic solution to the angular Mathieu equation, for \$n=0,1,2,3,…\$.

Pass a vector, e.g. `n=0:4`, to obtain multiple characteristic values efficiently.
"""
function ap(n::Order,q::Number,N::Integer=8+maximum(n)+ceil(Int,sqrt(abs(q))))
    C = mat_C1(q,N)
    a = eigvals(C)[n+1]
    return a
end

"""
    aa(n,q,N)

Compute the characteristic value \$a_{2n+1}(q)\$ corresponding to an even, π-antiperiodic (2π-periodic) solution to the angular Mathieu equation, for \$n=0,1,2,3,…\$.

Pass a vector, e.g. `n=0:4`, to obtain multiple characteristic values efficiently.
"""
function aa(n::Order,q::Number,N::Integer=8+maximum(n)+ceil(Int,sqrt(abs(q))))
    C = mat_C2(q,N)
    a = eigvals(C)[n+1]
    return a
end

"""
    a(m,q,N)

Compute the characteristic value \$a_m(q)\$. `m` must be \$2n\$ or \$2n+1\$, with \$n=0,1,2…\$

Pass a vector, e.g. `m=0:4`, to obtain multiple characteristic values efficiently.
"""
function a(m::Order,q::Number,N::Integer=8+maximum(m)+ceil(Int,sqrt(abs(q))))
    oddind = find(isodd,m)
    isempty(oddind)  && return ap(div.(m,2),q,N)
    evenind = find(iseven,m)
    isempty(evenind) && return aa(div.(m,2),q,N)
    a = zeros(length(m))
    a[evenind] = ap(div.(m[evenind],2),q,N)
    a[oddind] = aa(div.(m[oddind],2),q,N)
    return a
end

"""
    ba(n,q,N)

Compute the characteristic value \$b_{2n+1}(q)\$ corresponding to an odd, π-antiperiodic (2π-periodic) solution to the angular Mathieu equation, for \$n=0,1,2,3,…\$.

Pass a vector, e.g. `n=0:4`, to obtain multiple characteristic values efficiently.
"""
function ba(n::Order,q::Number,N::Integer=8+maximum(n)+ceil(Int,sqrt(abs(q))))
    C = mat_C3(q,N)
    b = eigvals(C)[n+1]
    return b
end

"""
    bp(n,q,N)

Compute the characteristic value \$b_{2n+2}(q)\$ corresponding to an odd, π-periodic solution to the angular Mathieu equation, for \$n=0,1,2,3,…\$.

Pass a vector, e.g. `n=0:4`, to obtain multiple characteristic values efficiently.
"""
function bp(n::Order,q::Number,N::Integer=8+maximum(n)+ceil(Int,sqrt(abs(q))))
    C = mat_C4(q,N)
    b = eigvals(C)[n+1]
    return b
end

"""
    b(m,q,N)

Compute the characteristic value \$b_m(q)\$. `m` must be \$2n+1\$ or \$2n+2\$, with \$n=0,1,2…\$

Pass a vector, e.g. `m=1:5`, to obtain multiple characteristic values efficiently.
"""
function b(m::Order,q::Number,N::Integer=8+maximum(m)+ceil(Int,sqrt(abs(q))))
    oddind = find(isodd,m)
    isempty(oddind)  && return bp(div.(m,2)-2,q,N)
    evenind = find(iseven,m)
    isempty(evenind) && return ba(div.(m,2)-1,q,N)
    b = zeros(length(m))
    b[evenind] = bp(div.(m[evenind]-2,2),q,N)
    b[oddind] = ba(div.(m[oddind]-1,2),q,N)
    return b
end

"""
    char_per(m,q,N)

Compute the characteristic value \$a_{m}(q)\$ for even `m`, or \$b_{2m}(q)\$ for odd `m`. These characteristic values correspond respectively to even and odd π-periodic solutions to the angular Mathieu equation, for \$n=0,1,2,3,…\$.

Pass a vector, e.g. `m=0:4`, to obtain multiple characteristic values efficiently.
"""
function char_per(m::Order,q::Number,N::Integer=8+maximum(m)+ceil(Int,sqrt(abs(q))))
    oddind = find(isodd,m)
    isempty(oddind)  && return ap(div.(m,2),q,N)
    evenind = find(iseven,m)
    isempty(evenind) && return bp(div.(m,2),q,N)
    a = zeros(length(m))
    a[evenind] = ap(div.(m[evenind],2),q,N)
    a[oddind]  = bp(div.(m[oddind],2),q,N)
    return a
end

"""
    char_aper(m,q,N)

Compute the characteristic value \$a_{m+1}(q)\$ for even `m`, or \$b_{2m+1}(q)\$ for odd `m`. These characteristic values correspond respectively to even and odd π-antiperiodic (2π-periodic) solutions to the angular Mathieu equation, for \$n=0,1,2,3,…\$.

Pass a vector, e.g. `m=0:4`, to obtain multiple characteristic values efficiently.
"""
function char_aper(m::Order,q::Number,N::Integer=8+maximum(m)+ceil(Int,sqrt(abs(q))))
    oddind = find(isodd,m)
    isempty(oddind)  && return aa(div.(m,2),q,N)
    evenind = find(iseven,m)
    isempty(evenind) && return ba(div.(m,2),q,N)
    a = zeros(length(m))
    a[evenind] = aa(div.(m[evenind],2),q,N)
    a[oddind]  = ba(div.(m[oddind],2),q,N)
    return a
end



#function char_per2(m::Order,q::Number,N::Integer=8+maximum(m)+ceil(Int,sqrt(abs(q))))
#    M = mat_per(q,N)
#    v = eigvals(M)[m+1]
#    return v
#end
#function char_aper2(m::Order,q::Number,N::Integer=8+maximum(m)+ceil(Int,sqrt(abs(q))))
#    M = mat_aper(q,N)
#    v = eigvals(M)[m+1]
#    return v
#end
