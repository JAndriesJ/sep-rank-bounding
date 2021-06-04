module Utils_Model
using JuMP


include(dirname(dirname(@__FILE__))*"\\Examples\\Utils_states.jl")
include(dirname(dirname(@__FILE__))*"\\Model\\R_sep_Model.jl")
include(dirname(dirname(@__FILE__))*"\\Model\\C_sep_Model.jl")
using .Utils_states
using .R_sep_Model
using .C_sep_Model

export export_model,
       batch_model,
       read_model


""" Exmports the JuMP model to SDPA format """
function export_model(model,file_name::String = pwd()*"\\Defualt_sep_mod.dat-s")
    JuMP.write_to_file(model, file_name, format=MOI.FileFormats.FORMAT_SDPA)
end

""" Reads a JuMP model from SDPA format text document """
function read_model(file_name::String)
    @assert splitext(file_name)[end] == ".dat-s" "Can only load models from .dat-s (SDPA) files!"
    return JuMP.read_from_file(file_name, format=MOI.FileFormats.FORMAT_SDPA)
end

""" Produces models for a set of constraints."""
function batch_model(t,ρ_dict,save_dir, con_list= ["S₁","S₂","S₃","S₁sG","S₂sG","S₃sG","S₁wG","S₂wG","S₃wG"],isR = true)
    for key in keys(ρ_dict)
        ρ = ρ_dict[key]
        d₁,d₂ = Utils_states.extr_d(key)
        d = (d₁,d₂)
        for con in con_list
            print(key*"__"*con)
            if isR
                sep_mod = R_sep_Model.Modelξₜˢᵉᵖ(ρ,d,t,con)
                save_name = "ℝ"*key*"_"*con*".dat-s"
            else
                sep_mod = C_sep_Model.Modelξₜˢᵉᵖ(ρ,d,t,con)
                save_name = "ℂ"*key*"_"*con*".dat-s"
            end


            save_path = save_dir*save_name
            export_model(sep_mod, save_path)
        end
    end
end


end
