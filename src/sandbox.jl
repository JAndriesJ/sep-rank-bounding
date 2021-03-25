    sourceDir = pwd()*"\\src\\"
    include(sourceDir*"data.jl")
    include(sourceDir*"moments.jl")
    include(sourceDir*"constraints.jl")
    include(sourceDir*"model.jl")
    include(sourceDir*"compute.jl")

    include(sourceDir*"batch_run.jl")

    using .Examples
    using .moments
    using .constraints
    using .sep_model
    using .compute

    using .batch_run
##
## batch
# run_batch()
## states examples:
d = 3
ρ_sep = get_sep_example(d)

ρ_ent = get_ent_example(d)


# Pkg.add("CSV")
# Pkg.add("StringDistances")

# using LinearAlgebra
#
# using JuMP
# model = Model()
# list_of_keys = make_mom_expo_keys(n, t) # Define variables in the moment matrix.
# @variable(model, Lx[list_of_keys] ) # ????
#
#
# d = 2
# ρ = ρᵈʳ[2,2]
# n       = 2*d
# MB      = make_mon_expo_mat(n,t-2)
# sqrt_tr_ρ   = sqrt(tr(ρ))
# s       = size(MB)[1]
# # expos of  (∑xᵢ²)
# XsqExpos = repeat([2*vcat(ones(Int64,d),zeros(Int64,d))], s, s)
# # expos of  (∑yᵢ²)
# YsqExpos = reverse(XsqExpos, dims = 2)
# #  √(ρₘₐₓ)⋅L([x,y]ₜ₋₁[x,y]ₜ₋₁ᵀ) - L((∑xᵢ²)⋅[x,y]ₜ₋₁[x,y]ₜ₋₁ᵀ )
# g₁ = sqrt_tr_ρ*index_to_var(Lx, MB) - index_to_var(Lx, MB + XsqExpos)
# #  √(ρₘₐₓ)⋅L([x,y]ₜ₋₁[x,y]ₜ₋₁ᵀ) - L((∑yᵢ²)⋅[x,y]ₜ₋₁[x,y]ₜ₋₁ᵀ )
# g₂ = sqrt_tr_ρ*index_to_var(Lx, MB) - index_to_var(Lx, MB + YsqExpos)
# loc_con = Dict(1 => g₁, 2 => g₂)
# return loc_con
#
#
# occursin("S₂", "S₂G")
# ρᵈʳ      = generate_random_states(2:4,2:5)
# con_list_1 = ["S₁"]
# con_list_2 = ["S₂"]
# t = 3
# d = 4
#
#
# sep_mod_1 = Modelξₜˢᵉᵖ(ρᵈʳ[d,2],t,d,con_list_1)
# ov_1 = Computeξₜˢᵉᵖ(sep_mod_1)
# sep_mod_2 = Modelξₜˢᵉᵖ(ρᵈʳ[d,2],t,d,con_list_2)
# ov_2 = Computeξₜˢᵉᵖ(sep_mod_2)

# quick_export_model(sep_mod_1)
