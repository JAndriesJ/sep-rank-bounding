module testUtils
using Test

srcDir = pwd()*"\\src\\"

include(srcDir*"Utils.jl")
using .Utils

@testset "eᵢ" begin
    n =7
    for k in 1:n
        e = eᵢ(n,k)
        @test length(e) == n
        @test maximum(e) == 1
        @test minimum(e) == 0
        @test sum(e) == 1
    end
end

@testset "var_kron" begin
    A = repeat([[1,0,1]],2,2)
    B = repeat([[0,1,0]],2,2)
    A_tens_B = var_kron(A,B)
    @test A_tens_B == repeat([[1,1,1]],4,4)

    n = rand(2:4)
    m = rand(2:5)
    k = n*m
    A = [rand(0:2,5) for i in 1:n, j in 1:n]
    B = [rand(0:2,5) for i in 1:m, j in 1:m]
    A_tens_B = var_kron(A,B)
    @test size(A_tens_B) == (k,k)
    # Need more
end

# index_to_var
end  # module testUtils
