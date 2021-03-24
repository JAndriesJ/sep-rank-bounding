module sep_model

export Modelξₜˢᵉᵖ,
       export_model

using JuMP
using LinearAlgebra
using Dates


include("C:\\Users\\andries\\all-my-codes\\ju-sep-rank\\src\\moments.jl")
using .moments
include("C:\\Users\\andries\\all-my-codes\\ju-sep-rank\\src\\constraints.jl")
using .constraints


function Modelξₜˢᵉᵖ(ρ,t,d,con_list)
    n = 2*d
    model = Model()
    list_of_keys = make_mom_expo_keys(n, t) # Define variables in the moment matrix.
    @variable(model, Lx[list_of_keys] ) # ????
# Build the moment matrix and constrain it to be PSD: L([x,y]≦ₜᵀ[x,y]≦ₜ) ⪰ 0 (doubble check this!!!!!!)
    mom_matₜ_expo = make_mon_expo_mat(n,t,true)
    mom_matₜ      = index_to_var(Lx, mom_matₜ_expo)
    @constraint(model, Symmetric(mom_matₜ) in PSDCone())
# Fourth order Moment constraints: L(xxᵀ ⊗ yyᵀ) = ρ,
    xxᵀ_tens_yyᵀ       = make_xxᵀ_tens_yyᵀ(d)
    L_xxᵀ_tens_yyᵀ  = index_to_var(Lx,xxᵀ_tens_yyᵀ)
    @constraint(model, fix_con,L_xxᵀ_tens_yyᵀ .==  ρ)
# Localizing g constraint: L ≥ 0 on M₂ₜ(S)
    if occursin("S₁",con_list)
        loc_con = make_loc_cons(ρ,t,d,Lx)
    elseif occursin("S₂",con_list)
        loc_con = make_loc_cons_var_1(ρ,t,d,Lx)
    end
    for key in keys(loc_con)
        @constraint(model, Symmetric(loc_con[key]) in PSDCone())
    end
# weak G Constraints
    if occursin("wG",con_list)
        println("----------------Weak G-constraints are active")
        weakG_con = make_weakG_con(ρ,t,d,Lx)
        for key in keys(weakG_con)
           @constraint(model, Symmetric(weakG_con[key]) in PSDCone())
        end
    end
# G Constraints
    if  occursin("G",con_list)
        println("----------------G-constraints are active")
        G_con                 = make_G_con(ρ,t,d,Lx)
        @constraint(model, Symmetric(G_con) in PSDCone())
    end
#  Set objective
    @objective(model, Min, Lx[zeros(n)])
    return model
end

"""
Exmports the JuMP model to SDPA format
"""
function export_model(model,file_name::String)
    JuMP.write_to_file(model, file_name, format=MOI.FileFormats.FORMAT_SDPA)
end

function quick_export_model(model)
    timestamp = replace(string(now()),":"=>"-")[6:16]
    saveLoc = "C:\\Users\\andries\\all-my-codes\\ju-sep-rank\\DAT-s\\"
    export_model(sep_mod_1,saveLoc*timestamp*"dat-s")
end


end
