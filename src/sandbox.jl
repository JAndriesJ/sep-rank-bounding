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

## states examples:
d = 3
n = 2*d

## Example states
t = 3
# ρ_sep = get_sep_example(d)
# run_batch(t,ρ_sep,"sep t = $t")

ρ_rand  = generate_random_states(2:4,2:9,343)
run_batch(t,ρ_rand,"rand_t=$t")

# ρ_ent = get_ent_example(d)
# run_batch(t,ρ_ent,"sep t = $t")




# ρ = ρ_dict[3, 2, "s"]
# sep_mod = Modelξₜˢᵉᵖ(ρ,t,d, "S₁")
# sep_mod_opt =  Computeξₜˢᵉᵖ(sep_mod)
