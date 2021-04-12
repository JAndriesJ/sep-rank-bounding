    sep_rank_proj_path = dirname(dirname(@__FILE__))
    Pkg.activate(sep_rank_proj_path)

    # Packages that are required
    # for pak in [] #["LinearAlgebra", "JuMP", "DataFrames", "CSV", "MosekTools"]
    #     Pkg.add(pak)
    # end

    sourceDir = sep_rank_proj_path*"\\src\\"
    testDir   = sep_rank_proj_path*"\\test\\"
    for mod in ["Examples","sep_Model","sep_Compute"]#, "Moments", "sep_Constraints", , "sep_Compute"]
        include(sourceDir*mod*".jl")
    end
    using .Examples  # Moments sep_Constraints sep_Model sep_Compute
    using .sep_Model
    using .sep_Compute

    ρ_dict = get_examples()

    # ρ_rand = ρ_dict["rand"]
    ρ_sep  = ρ_dict["sep"]
    # ρ_ent  = ρ_dict["ent"]

    ρ            = ρ_sep["sep3d3r3"]
    t            = 2
    con_list     = "S₁wG"
    sep_mod      = sep_Model.Modelξₜˢᵉᵖ(ρ,t,con_list)


sep_mod_opt  = sep_Compute.Computeξₜˢᵉᵖ(sep_mod)






# kat = dirname(@__FILE__)
# #
# basename("C:\\Users\\andries\\all-my-codes\\ju-text-manipulation-tools\\read_n_proc_tools.jl")
# poes = splitdir(kat)
# cnut = splitpath("C:\\Users\\andries\\all-my-codes\\ju-text-manipulation-tools\\read_n_proc_tools.jl")
# #  splitext("C:\\Users\\andries\\all-my-codes\\ju-text-manipulation-tools\\read_n_proc_tools.jl")
# kip = splitext("C:\\Users\\andries\\all-my-codes\\ju-text-manipulation-tools\\read_n_proc_tools.jl")[end]
