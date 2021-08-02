module testSep_Model
using Test
using JuMP

srcDir = pwd()*"\\src\\"
include(srcDir*"Examples.jl")
include(srcDir*"sep_Model.jl")
using .Examples
using .sep_Model

@testset "Model test" begin
    ρ_dict = Examples.get_examples()
    ρ_sep  = ρ_dict["sep"]
    ρ = ρ_sep["sep3d3r3"]
    t = 3
    con_list = "S₁S₂S₃sG"

    sep_mod = Modelξₜˢᵉᵖ(ρ,t,con_list)

    @test typeof(sep_mod) == JuMP.Model
    @test size(sep_mod[:Lx])[1] == 924
# end

# @testset "Model export test" begin
    export_model(sep_mod,pwd()*"\\Defualt_sep_mod.dat-s")
#end

#@testset "Model export test" begin
    sep_mod_read = read_model(pwd()*"\\Defualt_sep_mod.dat-s")
    # @test sep_mod_read  == sep_mod
end


end
