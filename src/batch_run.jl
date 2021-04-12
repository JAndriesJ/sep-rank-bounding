module Batch_run

using DataFrames
using JuMP

sourceDir = pwd()*"\\"
include(sourceDir*"data.jl")
include(sourceDir*"model.jl")
include(sourceDir*"compute.jl")

using .Examples
using .sep_model
using .compute


include("C:\\Users\\andries\\all-my-codes\\ju-text-manipulation-tools\\read_n_proc_tools.jl")
using .read_n_proc_tools

export batch_read_run







# function run_batch(t::Int, ρ_dict, save_file, con_list_list = ["S₁","S₂","S₃","S₁sG","S₂sG","S₃sG","S₁wG","S₂wG","S₃wG"])
#     file_loc = pwd()*"\\$save_file bounds.txt"
#     touch(file_loc)
#     for con in con_list_list
#         for key in keys(ρ_dict)
#             ρ = ρ_dict[key]
#             sep_mod     =  Modelξₜˢᵉᵖ(ρ,t, con)
#             sep_mod_opt =  Computeξₜˢᵉᵖ(sep_mod)
#             if typeof(sep_mod_opt) == Float64
#                 open(file_loc,"a") do io
#                     write(io, con*",$key, Error,Error,Nan, \n")
#                 end
#             else
#                 pstat = primal_status(sep_mod_opt)
#                 dstat =  dual_status(sep_mod_opt)
#                 ov =  objective_value(sep_mod_opt)
#
#                 open(file_loc,"a") do io
#                     write(io, con*",$key,$pstat,$dstat,$ov \n")
#                 end
#             end
#         end
#     end
# end
#


end
