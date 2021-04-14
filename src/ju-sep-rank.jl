# Pkg.activate(sep_rank_proj_path)

# Packages that are required
# for pak in [] #["LinearAlgebra", "JuMP", "DataFrames", "CSV", "MosekTools"]
#     Pkg.add(pak)
# end

sep_rank_proj_path = dirname(dirname(@__FILE__))
sourceDir = sep_rank_proj_path*"\\src\\"
testDir   = sep_rank_proj_path*"\\test\\"
for mod in ["Examples","sep_Model","sep_Compute"]#, "Moments", "sep_Constraints"]
    include(sourceDir*mod*".jl")
end
using .Examples
using .sep_Model
using .sep_Compute

# include(testDir*"runTests.jl")

ρ_dict = Examples.get_examples()
ρ_sep  = ρ_dict["sep"]
ρ            = ρ_sep["sep4d3r4"]
t            = 2
con_list     = "S₁"

sep_mod      = sep_Model.Modelξₜˢᵉᵖ(ρ,t,con_list)
sep_mod_opt  = sep_Compute.Computeξₜˢᵉᵖ(sep_mod)
Lx_vals = values.(sep_mod_opt[:Lx])
mom_mat_vals = sep_Compute.rec_mom_mat(Lx_vals)

# boundsDir = pwd()*"\\bounds\\"
# sep_Model.batch_model(t,ρ_sep,boundsDir)
# mass_read_comp(boundsDir)



# kat = dirname(@__FILE__)
# #
# basename("C:\\Users\\andries\\all-my-codes\\ju-text-manipulation-tools\\read_n_proc_tools.jl")
# poes = splitdir(kat)
# cnut = splitpath("C:\\Users\\andries\\all-my-codes\\ju-text-manipulation-tools\\read_n_proc_tools.jl")
# #  splitext("C:\\Users\\andries\\all-my-codes\\ju-text-manipulation-tools\\read_n_proc_tools.jl")
# kip = splitext("C:\\Users\\andries\\all-my-codes\\ju-text-manipulation-tools\\read_n_proc_tools.jl")[end]
