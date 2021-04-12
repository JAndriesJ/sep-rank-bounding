module sep_Model

using JuMP
using LinearAlgebra

srcDir = dirname(@__FILE__)*"\\"
include(srcDir*"Utils.jl")
include(srcDir*"Moments.jl")
include(srcDir*"sep_Constraints.jl")

using .Utils
using .Moments
using .sep_Constraints

export Modelξₜˢᵉᵖ,
       export_model,
       batch_model,
       read_model

"""
The model
"""
function Modelξₜˢᵉᵖ(ρ,t,con_list)
    d =  Int(sqrt(size(ρ)[1]))
    n = 2*d
    model = Model()
    list_of_keys = make_mom_expo_keys(n,t) # Define variables in the moment matrix.
    @variable(model, Lx[list_of_keys] ) # ????
 # Build the moment matrix and constrain it to be PSD: L([x,y]≦ₜᵀ[x,y]≦ₜ) ⪰ 0 (doubble check this!!!!!!)
    mom_matₜ_expo = make_mon_expo_mat(n,t,ρ,true)
    mom_matₜ      = Utils.index_to_var(Lx, mom_matₜ_expo)
    @constraint(model, Symmetric(mom_matₜ) in PSDCone())
 # Fourth order Moment constraints: L(xxᵀ ⊗ yyᵀ) = ρ,
    xxᵀ_tens_yyᵀ    = make_xxᵀ_tens_yyᵀ(d,ρ)
    L_xxᵀ_tens_yyᵀ  = Utils.index_to_var(Lx,xxᵀ_tens_yyᵀ)
    @constraint(model, fix_con,L_xxᵀ_tens_yyᵀ .==  prop_zero_diags(ρ))
 # Localizing g constraint: L ≥ 0 on M₂ₜ(S)
    if occursin("S₁",con_list)
        loc_con = make_loc_cons_S₁(ρ,t,d,Lx)
    elseif occursin("S₂",con_list)
        loc_con = make_loc_cons_S₂(ρ,t,d,Lx)
    elseif occursin("S₃",con_list)
        loc_con, g₂  =  make_loc_cons_S₃(ρ,t,d,Lx)
        @constraint(model, fix_S₃_con, g₂ .==  zeros(size(g₂)))
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
    if  occursin("sG",con_list)
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
function export_model(model,file_name::String = pwd()*"\\Defualt_sep_mod.dat-s")
    JuMP.write_to_file(model, file_name, format=MOI.FileFormats.FORMAT_SDPA)
end


function batch_model(t::Int,ρ_dict,save_dir, con_list_list = ["S₁","S₂","S₃","S₁sG","S₂sG","S₃sG","S₁wG","S₂wG","S₃wG"])
    for key in keys(ρ_dict)
        ρ = ρ_dict[key]
        for con in con_list_list
            sep_mod = Modelξₜˢᵉᵖ(ρ,t, con)

            save_name = key*"_"*con*".dat-s"
            save_path = save_dir*save_name
            export_model(sep_mod, save_path)
        end
    end
end

"""
Input: .dat-s filepath
Output: JuMP model
"""
function read_model(file_name::String)
    @assert splitext(file_name)[end] == ".dat-s" "Can only load models from .dat-s (SDPA) files!"
    return JuMP.read_from_file(file_name, format=MOI.FileFormats.FORMAT_SDPA)
end


end
