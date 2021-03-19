module seprank
sourceDir = "C:\\Users\\andries\\all-my-codes\\ju-sep-rank\\src\\"
include(sourceDir*"data.jl")
include(sourceDir*"moments.jl")
include(sourceDir*"constraints.jl")
include(sourceDir*"model.jl")
include(sourceDir*"compute.jl")

using .Exmaples
using .moments
using .constraints
using .sep_model
using .compute

t = 2
d = 4
ρ =  gen_rand_state(d,3)
con_list = []

sep_mod_1 = Modelξₜˢᵉᵖ(ρ,t,d,con_list)
Lx,model = Computeξₜˢᵉᵖ(sep_mod_1)

end

# include("C:\\Users\\andries\\all-my-codes\\ju-sep-rank\\src\\ju-sep-rank.jl)
