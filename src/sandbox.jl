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
for norma in ["no"]#, "tr", "dg"]
    for t in 2:2
        # ρ_rand         = generate_random_states(2:4,2:9,norma)
        # run_batch(t, ρ_rand, "rand_t=$t"*norma)

        ρ_sep = get_sep_example(norma)
        run_batch(t,ρ_sep, "sep_t=$t"*norma)
    end
end

#

# ρ_sep = get_sep_example()
# run_batch(t,ρ_sep,"sep_t=$t")
#
#
# for t in  2:2
#     ρ_sep = get_sep_example()
#     run_batch(t,ρ_sep,"sep_t=$t")
#
#     ρ_rand  = generate_random_states(2:4,2:9,343)
#     run_batch(t,ρ_rand,"rand_t=$t")
#
#     # ρ_ent = get_ent_example(d)
#     # run_batch(t,ρ_ent,"sep t = $t")
# end

# ρ =  ρ_sep[(3, 5, "s1")]
#
# sep_mod = Modelξₜˢᵉᵖ(ρ,t,3, "S₁")
# sep_mod_opt =  Computeξₜˢᵉᵖ(sep_mod)

# ρ = ρ_sep[(2, 2, "s")]
# ρ_sep[(3, 3, "s")]
# ρ_sep[(3, 4, "s")]
# ρ_sep[(3, 5, "s1")]
# ρ_sep[(3, 5, "s2")]


# These are also separable...yo!!!!!





# ρ = ρ_dict[3, 2, "s"]
# sep_mod = Modelξₜˢᵉᵖ(ρ,t,d, "S₁")
# sep_mod_opt =  Computeξₜˢᵉᵖ(sep_mod)
