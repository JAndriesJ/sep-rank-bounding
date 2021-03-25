module compute
export Computeξₜˢᵉᵖ,
       rec_mom_mat

using MosekTools
using JuMP


function Computeξₜˢᵉᵖ(model)
    set_optimizer(model, Mosek.Optimizer)
    # optimize
    optimize!(model)

    # output results
    println("Primal: ", primal_status(model))
    println("Dual: ", dual_status(model))
    println("Objective: ", objective_value(model))
    # println("Variables: ", "x = ",value(x))
    if string(primal_status(model)) == "NO_SOLUTION"
        return NaN
    else
        return model
    end
end


function rec_mom_mat(n::Int64,t::Int64,Lx)
    MB_exp  = make_mon_expo_mat(n,t)
    MB      = index_to_var(Lx, MB_exp)
    mom_mat = value.(MB)
    return mom_mat
end


end
