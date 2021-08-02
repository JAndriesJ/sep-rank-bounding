module testSep_Constraints
# sep_rank_proj_path = dirname(dirname(@__FILE__))
# Pkg.activate(sep_rank_proj_path)
using Test
using JuMP

srcDir = pwd()*"\\src\\"
include(srcDir*"Utils_states.jl")
include(srcDir*"Examples.jl")
include(srcDir*"Moments.jl")
include(srcDir*"R_constraints.jl")

using .Utils_states
using .Examples
using .Moments
using .R_constraints


examples = Examples.get_examples()

@testset "idx2var" begin
    t = (2,2)
    d = (2,2)
    n = sum(d)
    model = Model()
    mom_list = R_constraints.make_mon_expo_keys(n,t[1]) # Define variables in the moment matrix.
    @variable(model, Lx[mom_list] ) #
    @test length(Lx) == 70

    M_full = Moments.make_mon_expo(4,t;isle = true)
    M      = M_full[1:4,1:4]
    A      = reshape(1:16,4,4)
    HC =    [Lx[[0, 0, 0, 0]]  Lx[[1, 0, 0, 0]]  Lx[[0, 1, 0, 0]]  Lx[[0, 0, 1, 0]]
             Lx[[1, 0, 0, 0]]  Lx[[2, 0, 0, 0]]  Lx[[1, 1, 0, 0]]  Lx[[1, 0, 1, 0]]
             Lx[[0, 1, 0, 0]]  Lx[[1, 1, 0, 0]]  Lx[[0, 2, 0, 0]]  Lx[[0, 1, 1, 0]]
             Lx[[0, 0, 1, 0]]  Lx[[1, 0, 1, 0]]  Lx[[0, 1, 1, 0]]  Lx[[0, 0, 2, 0]]]
    @test Utils_cons.idx2var(Lx, M) ==  HC
    @test Utils_cons.idx2var(Lx, M, A) == A .* HC
    M_dict = Dict(1 => M,
                  2 => M_full[1:4,5:8],
                  3 => M_full[5:8,1:4])
    A_dict = Dict(1 => A,
                  2 => 1 ./ A,
                  3 => -A )
    @test Utils_cons.idx2var(Lx, M_dict, A_dict) == Utils_cons.idx2var(Lx,M,A) +
                                              Utils_cons.idx2var(Lx,M_full[1:4,5:8],1 ./ A) +
                                              Utils_cons.idx2var(Lx,M_full[5:8,1:4],-A)
end




## Separable states
@testset "Randd₁4d₂4" begin
    ex = "Randd₁4d₂4"
    d = Utils_states.extr_d(ex)
    @test d == (4,4)
    ρ = examples["sep"][ex]

    n = sum(d)
    model = Model()
    t = (3,3)
    list_of_keys = Moments.make_mom_expo_keys(n,t) # Define variables in the moment matrix.
    @variable(model, Lx[list_of_keys] ) #
    @test length(sep_Constraints.R_cons.make_loc_cons_S₁(ρ,d,t,Lx)) == 32
    @test length(sep_Constraints.R_cons.make_loc_cons_S₂(ρ,d,t,Lx)) == 8
    @test length(sep_Constraints.R_cons.make_loc_cons_S₃(ρ,d,t,Lx)) == 2
    @test length(sep_Constraints.R_cons.make_weakG_con(ρ,d,t,Lx)) == 4
    @test length(sep_Constraints.R_cons.make_G_con(ρ,d,t,Lx)) == 4
end

@testset "HHH1d₁2d₂2" begin
    ex = "HHH1d₁2d₂2"
    d = Utils_states.extr_d(ex)
    @test d == (2,2)
    ρ = examples["sep"][ex]

    n = sum(d)
    model = Model()
    t = (2,2)
    list_of_keys = Moments.make_mom_expo_keys(n,t) # Define variables in the moment matrix.
    @variable(model, Lx[list_of_keys] ) #
    @test length(sep_Constraints.R_cons.make_loc_cons_S₁(ρ,d,t,Lx)) == 16
    @test length(sep_Constraints.R_cons.make_loc_cons_S₂(ρ,d,t,Lx)) == 8
    @test length(sep_Constraints.R_cons.make_loc_cons_S₃(ρ,d,t,Lx)) == 2
    @test length(sep_Constraints.R_cons.make_weakG_con(ρ,d,t,Lx)) == 0
    @test length(sep_Constraints.R_cons.make_G_con(ρ,d,t,Lx)) == 4
end

@testset "DNY1d₁3d₂3" begin
    ex = "DNY1d₁3d₂3"
    d = Utils_states.extr_d(ex)
    @test d == (3,3)
    ρ = examples["sep"][ex]

    n = sum(d)
    model = Model()
    t = (4,4)
    list_of_keys = Moments.make_mom_expo_keys(n,t) # Define variables in the moment matrix.
    @variable(model, Lx[list_of_keys] ) #
    @test length(sep_Constraints.R_cons.make_loc_cons_S₁(ρ,d,t,Lx)) == 24
    @test length(sep_Constraints.R_cons.make_loc_cons_S₂(ρ,d,t,Lx)) == 8
    @test length(sep_Constraints.R_cons.make_loc_cons_S₃(ρ,d,t,Lx)) == 2
    @test length(sep_Constraints.R_cons.make_weakG_con(ρ,d,t,Lx)) == 8
    @test length(sep_Constraints.R_cons.make_G_con(ρ,d,t,Lx)) == 4
end


end















































































# ρ_dict = Examples.get_examples()
# ρ_sep  = ρ_dict["sep"]
# ρ = ρ_sep["sep3d3r3"]
#
# d =  Int(sqrt(size(ρ)[1]))
# n = 2*d
# t = 3
#
# @testset "Model setup" begin
#     model = Model()
# #end
#
# #@testset "Making variables" begin
#     list_of_keys = Moments.make_mom_expo_keys(n,t) # Define variables in the moment matrix.
#     @variable(model, Lx[list_of_keys] ) # Create the variables, i.e., entries of moment matrix.
#
#     @test size(Lx) == size(list_of_keys)
#     @test supertype(typeof(Lx)) == AbstractVector{VariableRef}
# #end
#
# #@testset "Build the moment matrix and constrain it to be PSD  : L([x,y]≦ₜᵀ[x,y]≦ₜ) ⪰ 0" begin
#     mom_matₜ_expo = Moments.make_mon_expo_mat(n,t,ρ,true)
#     mom_matₜ      = Utils.index_to_var(Lx, mom_matₜ_expo)
#
#     mom_matₜ₋₁_expo = Moments.make_mon_expo_mat(n,t-1,true)
#
#     @test size(mom_matₜ) == size(mom_matₜ_expo)
#     @test supertype(typeof(mom_matₜ)) ==  DenseArray{VariableRef,2}
# #end
#
# #@testset "Fourth order Moment constraints: L(xxᵀ ⊗ yyᵀ) = ρ," begin
#     xxᵀ_tens_yyᵀ    = Moments.make_xxᵀ_tens_yyᵀ(d,ρ)
#     L_xxᵀ_tens_yyᵀ  = Utils.index_to_var(Lx,xxᵀ_tens_yyᵀ)
#
#     @test typeof(L_xxᵀ_tens_yyᵀ) == Matrix{VariableRef}
#     @test size(L_xxᵀ_tens_yyᵀ) == size(Moments.prop_zero_diags(ρ))
# #end
#
# #@testset "Localizing g constraints: S₁" begin
#     loc_con_S₁ = make_loc_cons_S₁(ρ,t,d,Lx)
#     @test length(loc_con_S₁) == n
#     rand_key = [key for key in keys(loc_con_S₁)][rand(1:n)]
#     @test typeof(loc_con_S₁[rand_key]) == Matrix{GenericAffExpr{Float64,VariableRef}}
#     @test size(loc_con_S₁[rand_key]) == size(mom_matₜ₋₁_expo)
# #end
#
# #@testset "Localizing g constraints: S₂" begin
#     loc_con_S₂ = make_loc_cons_S₂(ρ,t,d,Lx)
#     @test length(loc_con_S₂) == 2
#     rand_key = [key for key in keys(loc_con_S₂)][rand(1:2)]
#     @test typeof(loc_con_S₂[rand_key]) == Matrix{GenericAffExpr{Float64,VariableRef}}
#     @test size(loc_con_S₂[rand_key]) == size(mom_matₜ₋₁_expo)
# #end
#
# #@testset "Localizing g constraints: S₃" begin
#     loc_con_S₃, g₂ = make_loc_cons_S₃(ρ,t,d,Lx)
#     @test length(loc_con_S₃) == 1
#     @test typeof(loc_con_S₃[1]) == Matrix{GenericAffExpr{Float64,VariableRef}}
#     @test size(loc_con_S₃[1]) == size(mom_matₜ₋₁_expo)
# #end
#
# #@testset "weak G Constraints" begin
#     weakG_con = make_weakG_con(ρ,t,d,Lx)
#     for key in keys(weakG_con)
#       @test size(weakG_con[key]) == size(make_mon_expo_mat(n,key,false)).*size(ρ)[1]
#       @test typeof(weakG_con[1]) == Matrix{GenericAffExpr{Float64,VariableRef}}
#     end
# #end
#
# #@testset "G Constraints" begin
#     G_con = make_G_con(ρ,t,d,Lx)
#     @test typeof(G_con) == Matrix{GenericAffExpr{Float64,VariableRef}}
#     @test size(G_con ) == size(make_mon_expo_mat(n,t-2,true)).*size(ρ)[1]
# end
#
#
#
#
# end
