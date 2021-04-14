module sep_Compute

using MosekTools
using JuMP

srcDir = dirname(@__FILE__)*"\\"
include(srcDir*"Utils.jl")
include(srcDir*"Moments.jl")
include(srcDir*"sep_Model.jl")
using .Utils
using .Moments
using .sep_Model


export Computeξₜˢᵉᵖ,
       rec_mom_mat,
       mass_read_comp

"""
Inpute: Jump model of separability problem.
Output: MOSEK solution response.
"""
function Computeξₜˢᵉᵖ(model)
    set_optimizer(model, Mosek.Optimizer)
    # optimize
    optimize!(model)

    # output results
    println("Primal: ", primal_status(model))
    println("Dual: ", dual_status(model))
    println("Objective: ", objective_value(model))
    return model
end


function mass_read_comp(boundsDir)
    datsFiles = [file for file in readdir(boundsDir,join = true) if contains(file,".dat-s")]
    file_loc = boundsDir*"Summary.csv"
    touch(file_loc)
    open(file_loc,"a") do io
        write(io, "|ex|rank|con|Primal|Dual|obj_val|\n")
        write(io, "|---|---|---|---|---|---|\n")
    end
    extract_sep_bound_file_meta(s) = s[1:4],s[6],s[8],split(split(s,"_")[end],".")[1]

    for file in datsFiles
        mod = sep_Model.read_model(file)
        sep_mod_opt = sep_Compute.Computeξₜˢᵉᵖ(mod)

        ex,d,rank,con   =  extract_sep_bound_file_meta(basename(file))

        pstat = primal_status(sep_mod_opt)
        dstat =  dual_status(sep_mod_opt)
        ov    =  objective_value(sep_mod_opt)

        open(file_loc,"a") do io
            write(io, "|$ex|$rank|$con|$pstat|$dstat|$ov| \n")
        end
    end
end



"""
Input: The variable list Lx returened by a solver.
Output: Moment matrix associated with the solution.
"""
function rec_mom_mat(Lx_vals)
    n = length(Lx_vals.axes[1][1])
    t = Int(sum(Lx_vals.axes[1][end])/2)
    MB_exp  = Moments.make_mon_expo_mat(n,t)
    MB      = Utils.index_to_var(Lx_vals, MB_exp)
    mom_mat = value.(MB)
    return mom_mat
end

end
