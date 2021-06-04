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

@testset "extr_d" begin
    ex_names = ["BPd₁3d₂3",    "CD1d₁2d₂3",    "CD1d₁2d₂4",    "CD1d₁2d₂5",    "CD2d₁3d₂4",    "CDsep1d₁2d₂2","CDsep2d₁2d₂2",
    "CDsep3d₁3d₂3","CDsep4d₁3d₂3", "CDsep5d₁3d₂3","DNY1d₁3d₂3",   "DNY2d₁2d₂2",  "HHH1d₁2d₂2",   "HHH2d₁2d₂2",
    "HKd₁3d₂3",    "Randd₁4d₂4",   "wiki1d₁2d₂2",  "wiki2d₁3d₂3"]
    ex_nums = [(3, 3),(2, 3),(2, 4),(2, 5),(3, 4),(2, 2),(2, 2),(3, 3),(3, 3),
               (3, 3),(3, 3),(2, 2),(2, 2),(2, 2),(3, 3),(4, 4),(2, 2),(3, 3)]
    for ex_ind in 1:length(ex_names)
        @test  Utils.extr_d(ex_names[ex_ind]) == ex_nums[ex_ind]
    end
    # Need more
end


# index_to_var
end  # module testUtils
