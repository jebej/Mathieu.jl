"""
    a(m,q)

Compute the characteristic value \$a_m(q)\$ corresponding to an even (cosine-elliptic) solution to the angular Mathieu equation. `m` is the order of the solution and must be an integer \$m=0,1,2,3…\$.

Pass a vector, e.g. `m=0:4`, to compute multiple characteristic values efficiently.
"""
function a(m::Order,q::Number,N::Integer=matsize2(m,q))
    ie = map(iseven,m)
    all(ie)  && return ap(div.(m,2),q,N) # all even
    !any(ie) && return aa(div.(m,2),q,N) # all odd
    a = zeros(length(m))
    a[ie]   = ap(div.(filter(iseven,m),2),q,N)
    a[.!ie] = aa(div.(filter(isodd, m),2),q,N)
    return a
end

"""
    b(m,q)

Compute the characteristic value \$b_m(q)\$ corresponding to an odd (sine-elliptic) solution to the angular Mathieu equation. `m` is the order of the solution and must be an integer, excluding 0, \$m=1,2,3,4…\$.

Pass a vector, e.g. `m=1:5`, to compute multiple characteristic values efficiently.
"""
function b(m::Order,q::Number,N::Integer=matsize2(m,q))
    ie = map(iseven,m)
    all(ie)  && return bp(div.(m.-2,2),q,N) # all even
    !any(ie) && return ba(div.(m.-1,2),q,N) # all odd
    b = zeros(length(m))
    b[ie]   = bp(div.(filter(iseven,m).-2,2),q,N)
    b[.!ie] = ba(div.(filter(isodd, m).-1,2),q,N)
    return b
end


"""
    ap(n,q)

Compute the characteristic value \$a_{2n}(q)\$ corresponding to an even (cosine-elliptic), π-periodic solution to the angular Mathieu equation, for \$n=0,1,2,3,…\$.

Pass a vector, e.g. `n=0:4`, to compute multiple characteristic values efficiently.
"""
function ap(n::Order,q::Number,N::Integer=matsize1(n,q))
    C = mat_C1(q,N)
    a = eigvals(C)[n.+1] # stev aka STERF
    #a = LAPACK.stebz!('I','E',0.0,0.0,1,maximum(n)+1,2*eps(q),C.dv,C.ev)[1][n+1]
    return a
end

"""
    aa(n,q)

Compute the characteristic value \$a_{2n+1}(q)\$ corresponding to an even (cosine-elliptic), π-antiperiodic (2π-periodic) solution to the angular Mathieu equation, for \$n=0,1,2,3,…\$.

Pass a vector, e.g. `n=0:4`, to compute multiple characteristic values efficiently.
"""
function aa(n::Order,q::Number,N::Integer=matsize1(n,q))
    C = mat_C2(q,N)
    a = eigvals(C)[n.+1]
    return a
end

"""
    ba(n,q)

Compute the characteristic value \$b_{2n+1}(q)\$ corresponding to an odd (sine-elliptic), π-antiperiodic (2π-periodic) solution to the angular Mathieu equation, for \$n=0,1,2,3,…\$.

Pass a vector, e.g. `n=0:4`, to compute multiple characteristic values efficiently.
"""
function ba(n::Order,q::Number,N::Integer=matsize1(n,q))
    C = mat_C3(q,N)
    b = eigvals(C)[n.+1]
    return b
end

"""
    bp(n,q)

Compute the characteristic value \$b_{2n+2}(q)\$ corresponding to an odd (sine-elliptic), π-periodic solution to the angular Mathieu equation, for \$n=0,1,2,3,…\$.

Pass a vector, e.g. `n=0:4`, to compute multiple characteristic values efficiently.
"""
function bp(n::Order,q::Number,N::Integer=matsize1(n,q))
    C = mat_C4(q,N)
    b = eigvals(C)[n.+1]
    return b
end


"""
    char_per(m,q)

Compute the characteristic value \$a_{m}(q)\$ for even `m`, or \$b_{m+1}(q)\$ for odd `m`. These characteristic values correspond respectively to even and odd π-periodic solutions to the angular Mathieu equation, for \$n=0,1,2,3,…\$.

|   m |  value |
| --- |--------|
|   0 |    a_0 |
|   1 |    b_2 |
|   2 |    a_2 |
|   3 |    b_4 |
|   4 |    a_4 |
|   ⋮ |       ⋮ |

Pass a vector, e.g. `m=0:4`, to compute multiple characteristic values efficiently.
"""
function char_per(m::Order,q::Number,N::Integer=matsize2(m,q))
    ie = map(iseven,m)
    all(ie)  && return ap(div.(m,2),q,N)
    !any(ie) && return bp(div.(m,2),q,N)
    a = zeros(length(m))
    a[ie]   = ap(div.(filter(iseven,m),2),q,N)
    a[.!ie] = bp(div.(filter(isodd, m),2),q,N)
    return a
end

"""
    char_aper(m,q)

Compute the characteristic value \$a_{m+1}(q)\$ for even `m`, or \$b_{m}(q)\$ for odd `m`. These characteristic values correspond respectively to even and odd π-antiperiodic (2π-periodic) solutions to the angular Mathieu equation, for \$n=0,1,2,3,…\$.

|   m |  value |
| --- |--------|
|   0 |    a_1 |
|   1 |    b_1 |
|   2 |    a_3 |
|   3 |    b_3 |
|   4 |    a_5 |
|   ⋮ |       ⋮ |

Pass a vector, e.g. `m=0:4`, to compute multiple characteristic values efficiently.
"""
function char_aper(m::Order,q::Number,N::Integer=matsize2(m,q))
    ie = map(iseven,m)
    all(ie)  && return aa(div.(m,2),q,N)
    !any(ie) && return ba(div.(m,2),q,N)
    a = zeros(length(m))
    a[ie]   = aa(div.(filter(iseven,m),2),q,N)
    a[.!ie] = ba(div.(filter(isodd, m),2),q,N)
    return a
end
