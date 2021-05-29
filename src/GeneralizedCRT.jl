module GeneralizedCRT

export crt
using Base.Threads
import Base.Threads.@spawn

const THRESHOLD1 = 7  # split longer input arrays to use binary instead of sequential algo
const THRESHOLD2 = 15 # split longer input arrays to different cpu threads

"""
    crt(a, b, p, q)

Given `0 <= a < p` and `0 <= b < q` with `a == b mod(gcd(p, q)`
return unique solution `x ≡ a mod p` and `x ≡ b mod q` 
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
    lcm = Base.checked_mul(p, q1)
    return x, lcm
end

"""
    crt(a::Vector, m::Vector; nthreads, splitstreads, splitbinary)

Given `0 <= a[i] < p[i] for all i` with `a[i] == a[j] mod(gcd(p[i], p[j]) for i != j`,
return unique solution `x ≡ a[i] mod p[i] for all i`
with `0 <= x < lcm(p, q)`.

Optional arguments:
Use `nthreads=nthreads()` threads if `length(a) > splitthreads=($THRESHOLD2)`.
Use binary algorithm is `length(a) > splitbinary=($THRESHOLD1)`.

Return `x, lcm(p, q)`.
"""
function crt(a::AbstractVector{T}, m::AbstractVector{S}; nt=nthreads(), splitthreads=THRESHOLD2, splitbinary=THRESHOLD1) where {S,T}
    n = length(a)
    n == length(m) || throw(ArgumentError("vectors of same size required"))
    if nt <= 1 || n <= splitthreads
        crts(a, m; splitbinary)
    else
        R = promote_type(S, T)
        t0 = time_ns()
        ni = max(n ÷ nt, 2)
        nt = (n + ni - 1) ÷ ni
        tasks = Vector{Any}(undef, nt)
        n2 = 0
        nth = nt
        for i = 1:nt
            n1 = n2 + 1
            n2 = n2 + (n - n2) ÷ nth
            nth -= 1
            av = view(a, n1:n2)
            mv = view(m, n1:n2)
            tasks[i] = @spawn crts($av, $mv; splitbinary=$splitbinary)
        end
        a0 = Vector{R}(undef, nt)
        m0 = Vector{R}(undef, nt)
        for i = 1:nt
            res = try
                fetch(tasks[i])
            catch ex
                ex isa TaskFailedException ? throw(ex.task.result) : rethrow()
            end
            a0[i] = first(res)
            m0[i] = last(res)
        end
        t1 = time_ns()
        x, l = crt(a0, m0)
        t2 = time_ns()
        x, l
    end
end

function crtm(a::AbstractVector{T}, m::AbstractVector{S}) where {S,T}
    n = length(a)
    # println("crtm($n, $(length(m)))")
    n == length(m) || throw(ArgumentError("vectors of same size required"))
    R = promote_type(S, T)
    n == 0 && return zero(R), one(R)
    xI, lcmI = promote(a[1], m[1])
    xI = mod(xI, lcmI)
    for i = 2:n
        xI, lcmI = crt(xI, a[i], lcmI, m[i])
    end
    xI, lcmI
end

function crts(a::AbstractVector{T}, m::AbstractVector{S}; splitbinary) where {S,T}
    n = length(a)
    n == length(m) || throw(ArgumentError("vectors of same size required"))
    if n <= max(splitbinary, 3)
        crtm(a, m)
    else
        n2 = (n + 1) ÷ 2
        xI, lcmI = crts(view(a,1:n2), view(m, 1:n2); splitbinary)
        xJ, lcmJ = crts(view(a, n2+1:n), view(m, n2+1:n); splitbinary)
        crt(xI, xJ, lcmI, lcmJ)
    end
end

end # module
