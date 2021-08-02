module Utils_Model
using JuMP
using DataFrames

include(dirname(dirname(@__FILE__))*"\\Model\\R_sep_Model.jl")
include(dirname(dirname(@__FILE__))*"\\Model\\C_sep_Model.jl")
using .R_sep_Model
using .C_sep_Model

export export_model,
       batch_model,
       read_model


""" Exmports the JuMP model to SDPA format """
export_model(model,file_name::String = pwd()*"\\Defualt_sep_mod.dat-s") =
    JuMP.write_to_file(model, file_name, format=MOI.FileFormats.FORMAT_SDPA)


""" Reads a JuMP model from SDPA format text document """
function read_model(file_name::String)
    @assert splitext(file_name)[end] == ".dat-s" "Can only load models from .dat-s (SDPA) files!"
    return JuMP.read_from_file(file_name, format=MOI.FileFormats.FORMAT_SDPA)
end



end
