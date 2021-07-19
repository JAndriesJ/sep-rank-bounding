module testMoments
using Test
using LinearAlgebra

srcDir = dirname(dirname(@__FILE__))*"\\src\\Moments\\"
include(srcDir*"Moments.jl")
import .Moments


const la = LinearAlgebra
const mom = Moments

@testset "eᵢ" begin
    n =7
    for k in 1:n
        e = mom.eᵢ(n,k)
        @test length(e) == n
        @test maximum(e) == 1
        @test minimum(e) == 0
        @test sum(e) == 1
    end
end

@testset "var_kron" begin
    A = repeat([[1,0,1]],2,2)
    B = repeat([[0,1,0]],2,2)
    A_tens_B =  mom.var_kron(A,B)
    @test A_tens_B == repeat([[1,1,1]],4,4)

    n = rand(2:4)
    m = rand(2:5)
    k = n*m
    A = [rand(0:2,5) for i in 1:n, j in 1:n]
    B = [rand(0:2,5) for i in 1:m, j in 1:m]
    A_tens_B =  mom.var_kron(A,B)
    @test size(A_tens_B) == (k,k)


    A =  mom.make_mon_expo(4,(1,1))
    B =  mom.make_mon_expo(4,(1,1);isle = false)
    C = mom.var_kron(A,B)
    for i in 1:4, j in 1:4
        rc(k) = (1+(k-1)*4):(4+(k-1)*4)
        @assert C[rc(i),rc(j)] == [A[i,j]] .+ B
    end

    # Need more
end

@testset "make_mom_expo_mat" begin
    mon_expo_mat =  mom.make_mon_expo(3,(1,1);isle = true)
    @test la.diag(mon_expo_mat) == [ [0, 0, 0],
                                  [2, 0, 0],
                                  [0, 2, 0],
                                  [0, 0, 2]]


    n = 4
    t = 2

    M_vec =  mom.make_mon_expo(n,t)
    M_mat =  mom.make_mon_expo(n,(t,t))

    ## Complex moments

    t = 2
    d = (2,2)
    n = sum(2 .* d)
    M_vec =  mom.make_mon_expo(d,t)
    @assert length(M_vec) == binomial(n+t,t)
    M_mat =  mom.make_mon_expo(d,(t,t))
    @assert M_vec == M_mat[1,:]

end

@testset "expo_conj" begin
    @test [mom.expo_conj((2,2), 1:8)] == [[3, 4, 1, 2, 7, 8, 5, 6]]
    @test [mom.expo_conj((3,3), 1:12)] == [[4, 5, 6, 1, 2, 3, 10, 11, 12, 7, 8, 9]]
end

@testset "make_mom_expo_keys" begin
    n = 3
    t = (1,1)
    mek = mom.make_mon_expo(n,2)
    @test mek == [[0, 0, 0],
                            [1, 0, 0],
                            [0, 1, 0],
                            [0, 0, 1],
                            [2, 0, 0],
                            [1, 1, 0],
                            [1, 0, 1],
                            [0, 2, 0],
                            [0, 1, 1],
                            [0, 0, 2]]

    for i in 1:3
        n = rand(3:6,1)[1]
        t = rand(2:4,1)[1]
        mek = mom.make_mon_expo(n,t*2)
        @test size(unique(mek))[1]  == size(mek)[1]
    end
end

## Real

@testset "make_xxᵀ_tens_yyᵀ" begin
    d = (3,3)
    xxᵀtensyyᵀ = mom.make_xxᵀ_tens_yyᵀ(d)
    @test size(xxᵀtensyyᵀ) == (prod(d),prod(d))

    @test diag(xxᵀtensyyᵀ) ==  [[2, 0, 0, 2, 0, 0],
                                [2, 0, 0, 0, 2, 0],
                                [2, 0, 0, 0, 0, 2],
                                [0, 2, 0, 2, 0, 0],
                                [0, 2, 0, 0, 2, 0],
                                [0, 2, 0, 0, 0, 2],
                                [0, 0, 2, 2, 0, 0],
                                [0, 0, 2, 0, 2, 0],
                                [0, 0, 2, 0, 0, 2]]

    @test sum.(xxᵀtensyyᵀ) == repeat([4],prod(d),prod(d))

    for i in 1:5
    d = rand(2:5,2)
    xxᵀtensyyᵀ = mom.make_xxᵀ_tens_yyᵀ(d)
    halfsum(arr) = [sum(arr[1:d[1]]),sum(arr[d[1]+1:end])]
    @test halfsum.(xxᵀtensyyᵀ) == repeat([[2,2]],prod(d),prod(d))
    end
end

@testset "get_ℝ_block_diagᵀ" begin
    d = (3,3)
    t = (4,4)
    ℝ_block_diag = mom.get_ℝ_block_diag(d,t)

    @test sum([size(ℝ_block_diag[b])[1] for b ∈ keys(ℝ_block_diag)]) == binomial(sum(d)+t[1],t[1])
end

## Complex

@testset "make_xx̄ᵀ_tens_yȳᵀ" begin
    d = rand(2:5,2)
    xx̄ᵀ_tens_yȳᵀ = mom.make_xx̄ᵀ_tens_yȳᵀ(d)
end

@testset "get_ℂ_block_diagᵀ" begin


end


end
