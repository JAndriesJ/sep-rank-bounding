module testSep_Constraints
# sep_rank_proj_path = dirname(dirname(@__FILE__))
# Pkg.activate(sep_rank_proj_path)
using Test
using JuMP

srcDir = pwd()*"\\src\\"
include(srcDir*"Utils.jl")
include(srcDir*"Examples.jl")
include(srcDir*"Moments.jl")
include(srcDir*"sep_Constraints.jl")

using .Utils
using .Examples
using .Moments
using .sep_Constraints

ρ_dict = Examples.get_examples()
ρ_sep  = ρ_dict["sep"]
const ρ = ρ_sep["sep3d3r3"]

const d =  Int(sqrt(size(ρ)[1]))
const n = 2*d
const t = 3

@testset "Model setup" begin
    model = Model()
#end

#@testset "Making variables" begin
    list_of_keys = Moments.make_mom_expo_keys(n,t) # Define variables in the moment matrix.
    @variable(model, Lx[list_of_keys] ) # Create the variables, i.e., entries of moment matrix.

    @test size(Lx) == size(list_of_keys)
    @test supertype(typeof(Lx)) == AbstractVector{VariableRef}
#end

#@testset "Build the moment matrix and constrain it to be PSD  : L([x,y]≦ₜᵀ[x,y]≦ₜ) ⪰ 0" begin
    mom_matₜ_expo = Moments.make_mon_expo_mat(n,t,ρ,true)
    mom_matₜ      = Utils.index_to_var(Lx, mom_matₜ_expo)

    mom_matₜ₋₁_expo = Moments.make_mon_expo_mat(n,t-1,true)

    @test size(mom_matₜ) == size(mom_matₜ_expo)
    @test supertype(typeof(mom_matₜ)) ==  DenseArray{VariableRef,2}
#end

#@testset "Fourth order Moment constraints: L(xxᵀ ⊗ yyᵀ) = ρ," begin
    xxᵀ_tens_yyᵀ    = Moments.make_xxᵀ_tens_yyᵀ(d,ρ)
    L_xxᵀ_tens_yyᵀ  = Utils.index_to_var(Lx,xxᵀ_tens_yyᵀ)

    @test typeof(L_xxᵀ_tens_yyᵀ) == Matrix{VariableRef}
    @test size(L_xxᵀ_tens_yyᵀ) == size(Moments.prop_zero_diags(ρ))
#end

#@testset "Localizing g constraints: S₁" begin
    loc_con_S₁ = make_loc_cons_S₁(ρ,t,d,Lx)
    @test length(loc_con_S₁) == n
    rand_key = [key for key in keys(loc_con_S₁)][rand(1:n)]
    @test typeof(loc_con_S₁[rand_key]) == Matrix{GenericAffExpr{Float64,VariableRef}}
    @test size(loc_con_S₁[rand_key]) == size(mom_matₜ₋₁_expo)
#end

#@testset "Localizing g constraints: S₂" begin
    loc_con_S₂ = make_loc_cons_S₂(ρ,t,d,Lx)
    @test length(loc_con_S₂) == 2
    rand_key = [key for key in keys(loc_con_S₂)][rand(1:2)]
    @test typeof(loc_con_S₂[rand_key]) == Matrix{GenericAffExpr{Float64,VariableRef}}
    @test size(loc_con_S₂[rand_key]) == size(mom_matₜ₋₁_expo)
#end

#@testset "Localizing g constraints: S₃" begin
    loc_con_S₃, g₂ = make_loc_cons_S₃(ρ,t,d,Lx)
    @test length(loc_con_S₃) == 1
    @test typeof(loc_con_S₃[1]) == Matrix{GenericAffExpr{Float64,VariableRef}}
    @test size(loc_con_S₃[1]) == size(mom_matₜ₋₁_expo)
#end

#@testset "weak G Constraints" begin
    weakG_con = make_weakG_con(ρ,t,d,Lx)
    for key in keys(weakG_con)
      @test size(weakG_con[key]) == size(make_mon_expo_mat(n,key,false)).*size(ρ)[1]
      @test typeof(weakG_con[1]) == Matrix{GenericAffExpr{Float64,VariableRef}}
    end
#end

#@testset "G Constraints" begin
    G_con = make_G_con(ρ,t,d,Lx)
    @test typeof(G_con) == Matrix{GenericAffExpr{Float64,VariableRef}}
    @test size(G_con ) == size(make_mon_expo_mat(n,t-2,true)).*size(ρ)[1]
end


end
