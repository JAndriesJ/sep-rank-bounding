module testMoments
using Test
using LinearAlgebra

srcDir = pwd()*"\\src\\"
include(srcDir*"Moments.jl")
using .Moments


@testset "make_mom_expo_mat" begin
    mon_expo_mat = make_mon_expo_mat(3,1,true)

    @test diag(mon_expo_mat) == [ [0, 0, 0],
                                  [2, 0, 0],
                                  [0, 2, 0],
                                  [0, 0, 2]]
end

@testset "make_xxᵀ_tens_yyᵀ" begin
    d = 3
    xxᵀ_tens_yyᵀ = make_xxᵀ_tens_yyᵀ(d)
    @test size(xxᵀ_tens_yyᵀ) == (d^2,d^2)

    @test diag(xxᵀ_tens_yyᵀ) == [[2, 0, 0, 2, 0, 0],
    [2, 0, 0, 0, 2, 0],
    [2, 0, 0, 0, 0, 2],
    [0, 2, 0, 2, 0, 0],
    [0, 2, 0, 0, 2, 0],
    [0, 2, 0, 0, 0, 2],
    [0, 0, 2, 2, 0, 0],
    [0, 0, 2, 0, 2, 0],
    [0, 0, 2, 0, 0, 2]]

    @test sum.(xxᵀ_tens_yyᵀ) == repeat([4],d^2,d^2)
end

@testset "make_mom_expo_keys" begin
    n = 3
    t = 1
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
        mom_expo_keys = make_mom_expo_keys(n,t)
        @test size(unique(mom_expo_keys))[1]  == size(mom_expo_keys)[1]
    end
end

@testset "zero propagation" begin
    d = 2
    n = 2*d
    t = 2
    isLeq = false

    ρ =  [0 0 0 0
          0 2 0 0
          0 0 0 0
          0 0 0 3]


    mon_expo_mat  = make_mon_expo_mat(n,t,false)
    xxᵀ_tens_yyᵀ  = make_xxᵀ_tens_yyᵀ(d)

    @test size(mon_expo_mat) == (10,10)
end

# @testset "make_mon_expo_arr" begin
#     rand_n = rand(1:10,1)
#     rand_t = rand(0:4,1)
#     for n ∈ rand_n, t ∈ rand_t
#         mon_arr_full = make_mon_expo_arr(n,t,true)
#         @test size(mon_arr_full)[2] == n
#         @test size(mon_arr_full)[1] == binomial(n +t,t)
#
#         mon_arr_part = make_mon_expo_arr(n,t,false)
#         @test size(mon_arr_part)[2] == n
#         # @test size(mon_arr_part)[1] == binomial(n +t,t)
#     end
# end
#
# @testset "make_mon_expo" begin
#     MonBase = make_mon_expo(2,3,false)
#     @test MonBase ==  [ [3, 0],[2, 1],[1, 2],[0, 3]]
#
#     MonBase = make_mon_expo(2,4,false)
#     @test MonBase ==  [[4, 0],[3, 1],[2, 2],[1, 3],[0, 4]]
# end


end
