module sep_Compute

using MosekTools
using JuMP

export Computeξₜˢᵉᵖ,
       rec_mom_mat

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


function batch_read_run(Dats_dir)
    t = split(split(Dats_dir,"t=")[end],"\\")[1]
    files_in_dir = read_n_proc_tools.collect_files_in_dir(Dats_dir,".dat-s")
    file         = files_in_dir[1]
    file_name    = read_n_proc_tools.extract_file_name(file);
    ex,d,r,con   =  read_n_proc_tools.extract_sep_bound_file_meta(file_name);


    file_loc = Dats_dir*ex[1:3]*"t=$t"*".csv"
    touch(file_loc)
    open(file_loc,"a") do io
        write(io, "ex,rank,con,Primal,Dual,obj_val\n")
    end

    for file in files_in_dir
        mod = sep_model.read_dat_s_model(file)
        sep_mod_opt = Computeξₜˢᵉᵖ(mod)

        file_name    = read_n_proc_tools.extract_file_name(file);
        ex,d,rank,con   =  read_n_proc_tools.extract_sep_bound_file_meta(file_name);


        pstat = primal_status(sep_mod_opt)
        dstat =  dual_status(sep_mod_opt)
        ov =  objective_value(sep_mod_opt)

        open(file_loc,"a") do io
            write(io, "$ex,$rank,$con,$pstat,$dstat,$ov \n")
        end
    end
end


"""
Input: The variable list Lx returened by a solver.
Output: Moment matrix associated with the solution.
"""
function rec_mom_mat(n::Int64,t::Int64,Lx)
    MB_exp  = make_mon_expo_mat(n,t)
    MB      = index_to_var(Lx, MB_exp)
    mom_mat = value.(MB)
    return mom_mat
end

end
