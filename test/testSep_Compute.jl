module testSep_Compute
using Test
using JuMP
srcDir = pwd()*"\\src\\"
include(srcDir*"sep_Model.jl")
include(srcDir*"sep_Compute.jl")
using .sep_Model
using .sep_Compute

@testset "Compute sep.-rank test" begin
    sep_mod = read_model(pwd()*"\\test\\Defualt_sep_mod.dat-s")
    sep_mod_opt = Computeξₜˢᵉᵖ(sep_mod)

    @test result_count(sep_mod_opt) == 1
    @test string(primal_status(sep_mod_opt)) == "FEASIBLE_POINT"
    @test string(dual_status(sep_mod_opt)) == "FEASIBLE_POINT"
    @test (objective_value(sep_mod_opt) - 2.828427006418139) < 0.0000001
end

# Lx = sep_mod_opt[:Lx]
# rec_mom_mat

# using MathOptInterface
# const MOI = MathOptInterface
#
# # this only works if you have a linear objective function (the model has a ScalarAffineFunction as its objective)
# obj = MOI.get(Mod, MOI.ObjectiveFunction{MOI.ScalarAffineFunction{Float64}}())
#
# # take a look at the terms
# obj.terms
# # from this you could extract your vector c
# c = zeros(4)
# for term in obj.terms
#     c[term.variable_index.value] = term.coefficient
# end
# @show(c)

end  # module testSep_Compute
