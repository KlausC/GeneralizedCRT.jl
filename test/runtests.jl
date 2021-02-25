using GeneralizedCRT
using Test

@testset "GeneralizedCRT" begin

@testset "crt2" begin
    a, b, p, q = 1, 11, 20, 50
    @test crt(a, b, p, q) == (61, 100) 
    a = 0
    @test_throws ArgumentError crt(a, b, p, q)
end

@testset "crt_short" begin
    a = [1, 11]
    m = [20, 50]
    @test crt(a, m) == crt(a..., m...)

    a = [1, 1, 1, 1, 1, 0]
    m = [2, 3, 4, 5, 6, 7]
    @test crt(a, m) == (301, 420)

    a = [5, 3, 1]
    m = [6, 8, 10]
    x, = crt(a, m)
    @test a == mod.(x, m)

    a = [2, 5, 7]
    m = [6, 9, 15]
    @test_throws ArgumentError crt(a, m)

    a = vec([2  0  3  3  17  0  0  0  3  1])
    m = vec([5  13  17  12  25  3  21  1  8  2])
    x, lc = crt(a, m)
    @test (x, lc) == (594867, 928200)
    @test a == mod.(x, m)
end

using Random
rng = MersenneTwister(1)

@testset "crt_long" begin
    n = GeneralizedCRT.THRESHOLD * 1000 - 100
    m = rand(rng, 1:999, n)
    l = lcm(big.(m))
    x = rand(rng, 1:l)
    a = mod.(x, m)
    @test crt(a, m) == (x, l)
end

end
