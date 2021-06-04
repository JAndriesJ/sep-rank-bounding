module testMoments
using Test
using LinearAlgebra

srcDir = dirname(dirname(@__FILE__))*"\\src\\Moments\\"
include(srcDir*"Moments.jl")
using .Moments


@testset "make_mom_expo_mat" begin
    mon_expo_mat = make_mon_expo_mat(3,(1,1),true)

    @test diag(mon_expo_mat) == [ [0, 0, 0],
                                  [2, 0, 0],
                                  [0, 2, 0],
                                  [0, 0, 2]]
end

@testset "make_mom_expo_keys" begin
    n = 3
    t = (1,1)
    mom_expo_keys = make_mom_expo_keys(n,t)
    Hard_code = [[0, 0, 0],
        [1, 0, 0],
        [0, 1, 0],
        [0, 0, 1],
        [2, 0, 0],
        [1, 1, 0],
        [1, 0, 1],
        [0, 2, 0],
        [0, 1, 1],
        [0, 0, 2]]

    @test mom_expo_keys ==  Hard_code

    for i in 1:3
        n = rand(3:6,1)[1]
        t = rand(2:4,1)[1]
        mom_expo_keys = make_mom_expo_keys(n,(t,t))
        @test size(unique(mom_expo_keys))[1]  == size(mom_expo_keys)[1]
    end
end

@testset "make_xxᵀ_tens_yyᵀ" begin
    d = (3,3)
    xxᵀ_tens_yyᵀ = make_xxᵀ_tens_yyᵀ(d)
    @test size(xxᵀ_tens_yyᵀ) == (prod(d),prod(d))

    @test diag(xxᵀ_tens_yyᵀ) == [[2, 0, 0, 2, 0, 0],
    [2, 0, 0, 0, 2, 0],
    [2, 0, 0, 0, 0, 2],
    [0, 2, 0, 2, 0, 0],
    [0, 2, 0, 0, 2, 0],
    [0, 2, 0, 0, 0, 2],
    [0, 0, 2, 2, 0, 0],
    [0, 0, 2, 0, 2, 0],
    [0, 0, 2, 0, 0, 2]]

    @test sum.(xxᵀ_tens_yyᵀ) == repeat([4],prod(d),prod(d))

    for i in 1:5
    d = rand(2:5,2)
    xxᵀ_tens_yyᵀ = Moments.make_xxᵀ_tens_yyᵀ(d)
    halfsum(arr) = [sum(arr[1:d[1]]),sum(arr[d[1]+1:end])]
    @test halfsum.(xxᵀ_tens_yyᵀ) == repeat([[2,2]],prod(d),prod(d))
    end
end

@testset "xx̄ᵀ_tens_yȳᵀ" begin
    d = rand(2:5,2)
    xx̄ᵀ_tens_yȳᵀ = Moments.make_xx̄ᵀ_tens_yȳᵀ(d)
end


@testset "get_even_sub_mats" begin
    d₁,d₂ = rand(2:4,2)
    n = 2*d₁+ 2*d₂
    MM = Dict()
    MM[1] = Moments.make_mon_expo_mat(n,(0,0),true)
    MM[2] = Moments.make_mon_expo_mat(n,(1,1),true)
    MM[3] = Moments.make_mon_expo_mat(n,(2,2),true)
    MM[4] = Moments.make_mon_expo_mat(n,(3,3),true)
    MM[5] = Moments.make_mon_expo_mat(n,(4,4),true)

    SB = Dict()
    SB[1] = Moments.get_even_sub_mats(MM[1],(d₁,d₂))
    SB[2] = Moments.get_even_sub_mats(MM[2],(d₁,d₂))
    SB[3] = Moments.get_even_sub_mats(MM[3],(d₁,d₂))
    SB[4] = Moments.get_even_sub_mats(MM[4],(d₁,d₂))
    SB[5] = Moments.get_even_sub_mats(MM[5],(d₁,d₂))

    for ind in 1:5
        @test sum([ collect(size(SB[ind][key])) for key in keys(SB[ind])]) == collect(size(MM[ind]))
    end
end


end
