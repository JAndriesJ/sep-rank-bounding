module sep_Compute

using MosekTools
using JuMP
using LinearAlgebra


srcDir = dirname(dirname(@__FILE__))*"\\src\\"
include(srcDir*"Examples\\Examples.jl")
include(srcDir*"Model\\Utils_Model.jl")
include(srcDir*"Model\\C_sep_Model.jl")
include(srcDir*"Model\\R_sep_Model.jl")
using .Utils_Model
using .Examples
using .C_sep_Model
using .R_sep_Model

export Computeξₜˢᵉᵖ,
       get_sol_vals,
       get_sol_min_eigval,
       batch_Computeξₜˢᵉᵖ,
       quick_run_spec_example


"""
Inpute: Jump model of separability problem.
Output: MOSEK solution response.
"""
function Computeξₜˢᵉᵖ(model)
    JuMP.set_optimizer(model, Mosek.Optimizer)
    # optimize
    optimize!(model)

    # output results
    println("Primal: ", primal_status(model))
    println("Dual: ", dual_status(model))
    println("Objective: ", objective_value(model))
    return model
end

# get_sol_dict(Lx) = Dict(zip(Lx.data,value.(Lx).data))
get_sol_vals(arr) = JuMP.value.(arr)
get_sol_min_eigval(arr) = LinearAlgebra.eigmin(JuMP.value.(arr))




end
