dirname(@__FILE__)

### Begin here
println("""Initialize project:""")
    srcDir  = dirname(@__FILE__)*"\\"
    include(srcDir*"Examples\\Examples.jl")
    include(srcDir*"Examples\\Utils_states.jl")
    include("Moments.jl")
    include(srcDir*"Constraints\\C_constraints.jl")
    include(srcDir*"Constraints\\R_constraints.jl")
    include(srcDir*"Constraints\\Utils_cons.jl")
    include(srcDir*"Model\\R_sep_Model.jl")
    include(srcDir*"Model\\C_sep_Model.jl")
    include(srcDir*"Model\\Utils_Model.jl")
    include(srcDir*"sep_Compute.jl")


    using .Examples
    using .Utils_states
    using .Moments
    using .Utils_cons
    using .C_constraints
    using .R_constraints
    using .Utils_Model
    using .R_sep_Model
    using .C_sep_Model
    using .sep_Compute

## Examples
println("""Load examples:
examples     = Examples.get_examples()
examples_all = merge(examples["ent"],examples["sep"])
""")
examples     = Examples.get_examples()
display(examples_all)

println("""To get a tabular overview of the examples use the following code:
df = Examples.get_example_overview(false)
show(df)""")
df      = Examples.get_example_overview(false)
show(df)

println("""You can choose a specific example from the above list as follows:
ρ = examples_all["Eq25HHH96"]""")
ρ = examples_all["Eq25HHH96"]
display(ρ)
println("""Alternatively you can give your own example as follows:
ρ = [1 0 0 0
     0 0 0 0
     0 0 0 0
     0 0 0 1]""")
## The real model.
println("""Specify: ρ, registers: d, constraint list: cl,level of hierarchy: t
ρ =  [1 0 0 0
      0 0 0 0
      0 0 0 0
      0 0 0 1]
d  = (2,2)
cl = "S∞ sG" or "S₂ sG" or "S₂₁ sG"
t  = (2,2)""")
ρ  =  [1. 0. 0. 0.
       0. 0. 0. 0.
       0. 0. 0. 0.
       0. 0. 0. 1.]
d  = (2,2)
cl = "S∞ sG"
t  = (3,3)
println("""Generate the model:
sep_mod = R_sep_Model.Modelξₜˢᵉᵖ(ρ,d,t;con_list = cl)
""")
sep_mod = R_sep_Model.Modelξₜˢᵉᵖ(ρ,d,t;con_list = cl)
display(sep_mod)

## Solve the OP.
println("We initialize the computation module")
sep_mod = sep_Compute.Computeξₜˢᵉᵖ(sep_mod)
Lx = sep_mod[:Lx]
## A closer look at the real constraints
println("""Let us take a deeper look into the constraints:""")
println("These are $(binomial(sum(d)+2*t[1],2*t[1])) the moments present in the hierachy at level t = (2,2)")
mon_expo_keys = R_constraints.make_mon_expo_keys(sum(d),t[1])
display(mon_expo_keys)

println("""These are the blocks of the 'L([x,y]ₜ[x,y]ₜᵀ)'.
 """)

zero_moms = R_constraints.zeroprop(ρ,d,(3,3),Lx)


PSD_con = R_constraints.make_PSD_con(d,t,Lx)
display(PSD_con)
println("""
"ee": moments of the form xᵅyᵝ ⋅ xᵞyᵟ where |α|,|β|,|γ|,|δ| are all even.
"eo": moments of the form xᵅyᵝ ⋅ xᵞyᵟ where |α|,|γ| even and |β|,|δ| odd.
"oe": moments of the form xᵅyᵝ ⋅ xᵞyᵟ where |α|,|γ| odd and |β|,|δ| even.
"oo": moments of the form xᵅyᵝ ⋅ xᵞyᵟ where |α|,|β|,|γ|,|δ| are all odd.
 """)
display(PSD_con["ee"])
display(PSD_con["eo"])
display(PSD_con["oe"])
display(PSD_con["oo"])

println("This is 'L(xxᵀ⊗yyᵀ)'.")
ord4_con = R_constraints.make_ord4_con(d,Lx)
display(ord4_con)

println("This is the 'L(( √ρₘₐₓ - xᵢ²)⋅η), L(( √ρₘₐₓ - yⱼ²)⋅η) ⪰ 0 ∀ i ∈ [d₁], j ∈ [d₂]' constraints.")
cons_S_inf = R_constraints.make_loc_cons_S_inf(ρ,d,t,Lx)
display(cons_S_inf)
display(cons_S_inf[("ee", "x²_1")])
display(cons_S_inf[("eo", "x²_1")])
display(cons_S_inf[("oe", "x²_1")])
# display(cons_S_inf[("oo", "x_1")])

println("This is the L(( (Tr(ρ)) - ∑xᵢ²)⋅η), L(( (Tr(ρ)) - ∑yⱼ²)⋅η) ⪰ 0' constraints.")
cons_S₂    = R_constraints.make_loc_cons_S₂(ρ,d,t,Lx)
display(cons_S₂)
display(cons_S₂[("ee", "x")])
display(cons_S₂[("eo", "x")])
display(cons_S₂[("oe", "x")])

println("This is the ' L(( (Tr(ρ)) - ∑xᵢ²)⋅η) ⪰ 0, L(( 1 - ∑yⱼ²)⋅η) = 0' constraints.")
cons_S₂₁x,cons_S₂₁y   = R_constraints.make_loc_cons_S₂₁(ρ,d,t,Lx)
display(cons_S₂₁x)
display(cons_S₂₁y)


println("This is the 'ρ⊗L(η) - L( (xxᵀ⊗yyᵀ) ⊗ η ) ⪰ 0 ∀ η ∈ even-degree-principle-submatrices of ([x,y]₌ₜ₋₂[x,y]₌ₜ₋₂ᵀ)' constraints.")
G_con = R_constraints.make_G_con(ρ,d,t,Lx)
display(G_con["ee"])
## Inspecting the solution
using JuMP
# println("We can recover the moments as with the following command: Dict(zip(Lx.data,value.(Lx.data)))")
# display(Dict(zip(Lx.data,value.(Lx.data))))
#
# display(value.(PSD_con["ee"]))
#
# println("The fourth order constraint should equal ρ")
# display(value.(ord4_con))
#
# println("One of the S_inf constraints is")
# display(value.(cons_S_inf[("oe", "x²_2")]))
#
# println("The G constriant blocks are")
# display(value.(G_con["ee"]))
## The Real model quick
ex = ["Eq22HHH96a" "RANDa3" "p12HK14b" "T1CD12i"][1]
ρ = Examples.get_example(7)
d = (3,3)#filter(:Example => ==(ex),df).Size[1]
t = (2,2)
cl = "S₂"#"S∞ sG" # "S₂ sG" # "S₂₁ sG"
sepr_mod = R_sep_Model.Modelξₜˢᵉᵖ(ρ,d,t;con_list = cl)
sepr_mod = sep_Compute.Computeξₜˢᵉᵖ(sepr_mod)


display(filter(:Example => ==(ex),df))

#Lrx = sepr_mod[:Lx]
## The complex model quick
include(srcDir*"Model\\C_sep_Model.jl")
using .C_sep_Model
ex = ["RANDℂb2" "Eq2-3BP00b" "RANDa3"][3]
ρ = Examples.get_example(ex)
d = filter(:Example => ==(ex),df).Size[1]
t = (2,2)
cl = "S∞sG" # "S₂ sG" # "S₂₁ sG"
sep_mod = C_sep_Model.Modelξₜˢᵉᵖ(ρ,d,t;con_list = cl)
sep_mod = sep_Compute.Computeξₜˢᵉᵖ(sep_mod)
display(filter(:Example => ==(ex),df))

Lx = sep_mod[:Lx]
## Running things in batches
include(srcDir*"batch.jl")
using .batch
t = (2,2)
batch.batch_model(t,Examples.get_examples(),df)
batch.batch_Computeξₜˢᵉᵖ("C:\\Users\\andries\\all-my-codes\\ju-sep-rank\\assets\\bounds\\")
nar = batch.unstack_constraints(df)


## A closer look at the complex constraints
println("""Let us take a deeper look into the constraints:""")
include(srcDir*"Constraints\\C_constraints.jl")
using .C_constraints
#println("These are $(binomial(sum(d)+2*t[1],2*t[1])) the moments present in the hierachy at level t = (2,2)")
mon_expo_keys = C_constraints.make_mon_expo_keys(d,t[1])
display(mon_expo_keys)


println("""These are the blocks of the 'L([x,y]ₜ[x,y]ₜᵀ)'.
 """)

#zero_moms = C_constraints.zeroprop(ρ,d,(3,3),Lx)
PSD_con = C_constraints.make_PSD_con(d,t,Lx)
display(PSD_con)

println("This is 'L(xx*⊗yy*)'.")
ord4_con = C_constraints.make_ord4_con(d,Lx)
display(ord4_con)

println("This is the 'L(( √ρₘₐₓ - xᵢ²)⋅η), L(( √ρₘₐₓ - yⱼ²)⋅η) ⪰ 0 ∀ i ∈ [d₁], j ∈ [d₂]' constraints.")
cons_S_inf = C_constraints.make_loc_cons_S_inf(ρ,d,t,Lx)
display(cons_S_inf)


println("This is the L(( (Tr(ρ)) - ∑xᵢ²)⋅η), L(( (Tr(ρ)) - ∑yⱼ²)⋅η) ⪰ 0' constraints.")
cons_S₂    = C_constraints.make_loc_cons_S₂(ρ,d,t,Lx)
display(cons_S₂)


println("This is the ' L(( (Tr(ρ)) - ∑xᵢ²)⋅η) ⪰ 0, L(( 1 - ∑yⱼ²)⋅η) = 0' constraints.")
cons_S₂₁x,cons_S₂₁y   = C_constraints.make_loc_cons_S₂₁(ρ,d,t,Lx)
display(cons_S₂₁x)
display(cons_S₂₁y)


println("This is the 'ρ⊗L(η) - L( (xxᵀ⊗yyᵀ) ⊗ η ) ⪰ 0 ∀ η ∈ even-degree-principle-submatrices of ([x,y]₌ₜ₋₂[x,y]₌ₜ₋₂ᵀ)' constraints.")
G_con = C_constraints.make_G_con(ρ,d,t,Lx)
display(G_con)



include(srcDir*"Moments\\Moments.jl")
using .Moments
const mom = Moments

expo_conj(d, ααᶥββᶥ) = vcat(ααᶥββᶥ[d[1]+1:2*d[1]],
                            ααᶥββᶥ[1:d[1]],
                            ααᶥββᶥ[1+2*d[1]+d[2]:end],
                            ααᶥββᶥ[1+2*d[1]:2*d[1]+d[2]])

d = (2,2)
t = (2,2)
ρ  =  [1. 0. 0. 0.
       0. 0. 0. 0.
       0. 0. 0. 0.
       0. 0. 0. 1.]




# include(srcDir*"Constraints\\Utils_cons.jl")
# include(srcDir*"\\Constraints\\C_constraints.jl")
# using .Utils_cons
# using .C_constraints
# using LinearAlgebra
#
# CPSD_con = C_constraints.make_PSD_con(d,t,Lx)
# eigmin(value.(CPSD_con["Default"]))
# eigvals(value.(CPSD_con["Default"]))
# Cord4_con = C_constraints.make_ord4_con(d,Lx)
# value.(Cord4_con["real"])
# value.(Cord4_con["imag"])
#
# Cloc_cons_S₂ = C_constraints.make_loc_cons_S₂(ρ,d,t,Lx)
#
# value.(Cloc_cons_S₂[("Default", "x")])
# Cloc_cons_S₂[("Default", "x")][2,2]
#
# for i in [1,2,3], j in 1:8
#     println(i,j)
# end

#
# include(srcDir*"Examples\\Examples.jl")
# include(srcDir *"Examples\\Utils_states.jl")
# using .Examples
# using .Utils_states
# exd     = Examples.get_examples()
# ex = merge(exd["ent"],exd["sep"])
# keys(ex)
# ρ = ex["RANDa3"]
# Export_example_overview = false
# df           = Examples.get_example_overview(Export_example_overview);


# using LinearAlgebra
# function makediagone(ρ)
#     mpli = 1 ./ map(x -> iszero(x) ? 1 : x,  sqrt.([ρ[i,i] for i in 1:size(ρ)[1]]))
#     return diagm(mpli)*ρ*diagm(mpli)
# end
# ρ = [i+j*im for i in 1:4, j in 1:4]
# makediagone(ρ)
