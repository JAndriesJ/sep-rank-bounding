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
r = 3
ρ =  gen_rand_state(d,r)
con_list = ["S₁","G"]

sep_mod_1 = Modelξₜˢᵉᵖ(ρ,t,d,con_list)
ov = Computeξₜˢᵉᵖ(sep_mod_1)
quick_export_model(sep_mod_1)
