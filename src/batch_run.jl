module batch_run
using DataFrames
sourceDir = "C:\\Users\\andries\\all-my-codes\\ju-sep-rank\\src\\"
include(sourceDir*"data.jl")
include(sourceDir*"model.jl")
include(sourceDir*"compute.jl")

using .Examples
using .sep_model
using .compute

# using StringDistances
# using CSV


export run_batch

function run_batch()
    d_set = 2:6
    r_set = 2:9
    ρᵈʳ      = generate_random_states(d_set,r_set,343)
    con_list_list = ["S₁","S₂","S₁G","S₂G","S₁wG","S₂wG"] #
    t = 3

    # bounds_array = Dict()

    touch("C:\\Users\\andries\\all-my-codes\\ju-sep-rank\\bounds.txt")
    for con in con_list_list
        # bounds_array[con] = Dict()
        for key in keys(ρᵈʳ)
            ρ = ρᵈʳ[key]
            sep_mod = Modelξₜˢᵉᵖ(ρ,t,key[1], con)
            ov =  Computeξₜˢᵉᵖ(sep_mod)
            # pstat = primal_status(model)
            # dstat =  dual_status(model)
            # ov =  objective_value(sep_mod_opt)
            # con_ind = findall(con, con_list_list, Levenshtein())[1]
            # bounds_array[con][key...] = ov
            open("C:\\Users\\andries\\all-my-codes\\ju-sep-rank\\bounds.txt","a") do io
                write(io, con*"$key $ov \n")
            end
        end

    # bounds_mat = [ bounds_array[con][d,r] for d in d_set, r in r_set ]
    # d_names = vec(string.(["d="], d_set))
    # r_names = vec(string.(["r="], r_set))
    # df = DataFrame(bounds_mat,r_names)
    # df.d_size = d_names
    # select!(df, r"d", :)
    #
    # file_name = con*".csv"
    # CSV.write(file_name, df)
    # end
    end
end
end
