module batch_run
using DataFrames
sourceDir = pwd()*"\\src\\"
include(sourceDir*"data.jl")
include(sourceDir*"model.jl")
include(sourceDir*"compute.jl")

using .Examples
using .sep_model
using .compute
using JuMP
# using StringDistances
# using CSV

export run_batch


function run_batch(t::Int,ρ_dict,save_file)
    con_list_list = ["S₁","S₂","S₁G","S₂G","S₁wG","S₂wG"] #
    # bounds_array = Dict()
    file_loc = pwd()*"\\$save_file bounds.txt"
    touch(file_loc)
    for con in con_list_list
        # bounds_array[con] = Dict()
        for key in keys(ρ_dict)
            try
                ρ = ρ_dict[key]
                sep_mod = Modelξₜˢᵉᵖ(ρ,t,key[1], con)
                sep_mod_opt =  Computeξₜˢᵉᵖ(sep_mod)
                pstat = primal_status(sep_mod_opt)
                dstat =  dual_status(sep_mod_opt)
                ov =  objective_value(sep_mod_opt)

                open(file_loc,"a") do io
                    write(io, con*",$key,$pstat,$dstat,$ov \n")
                    # write(io, con*",$key,test \n")
                end

            catch err
                open(file_loc,"a") do io
                    write(io, con*",$key, Error,Error,Nan, \n")
                    # write(io, con*",$key,test \n")
                end
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
