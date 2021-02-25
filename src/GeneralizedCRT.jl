module GeneralizedCRT

export crt

"""
    crt(a, b, p, q)

Given `0 <= a < p` and `0 <= b < q` with `a == b mod(gcd(p, q)`
return unique solution `x ≡ a mod p and y ≡ b mod q` 
with `0 <= x < lcm(p, q)`.

Return `x, lcm(p, q)`.
"""
function crt(a::A, b::B, p::C, q::D) where {A,B,C,D}
    crt(promote(a, b, p, q)...)
end
function crt(a::T, b::T, p::T, q::T) where T
    c, u = gcdx(p, q)
    q1 = div(q, c)
    u = mod(u, q1)
    iszero(mod(b - a, c)) || throw(ArgumentError("remainders must have same rest when diveded by ($a, $b, gcd($p, $q))"))
    bac = div(b - a, c)
    x = mod(bac * u, q1) * p + a
    lcm = p * q1
    return x, lcm
end

"""
    crt(a::Vector, m::Vector)

Given `0 <= a[i] < p[i] for all i` with `a[i] == a[j] mod(gcd(p[i], p[j]) for i != j`,
return unique solution `x ≡ a[i] mod p[i] for all i`
with `0 <= x < lcm(p, q)`.

Return `x, lcm(p, q)`.

"""
function crtm(a::AbstractVector{T}, m::AbstractVector{S}) where {S,T}
    n = length(a)
    n == length(m) || throw(ArgumentError("vectors of same size required"))
    xI, lcmI = one(T), one(S)
    for i = 1:n
        xI, lcmI = crt(xI, a[i], lcmI, m[i])
    end
    xI, lcmI
end

const THRESHOLD = 100

function crt(a::AbstractVector{T}, m::AbstractVector{S}) where {S,T}
    n = length(a)
    n == length(m) || throw(ArgumentError("vectors of same size required"))
    if n < THRESHOLD
        crtm(a, m)
    else
        n2 = (n + 1) ÷ 2
        xI, lcmI = crt(view(a,1:n2), view(m, 1:n2))
        xJ, lcmJ = crt(view(a, n2+1:n), view(m, n2+1:n))
        crt(xI, xJ, lcmI, lcmJ)
    end
end

end
