module batch
srcDir  = dirname(@__FILE__)*"\\";

include(srcDir*"Model\\R_sep_Model.jl")
include(srcDir*"Model\\C_sep_Model.jl")
include(srcDir*"Model\\Utils_Model.jl")
include(srcDir*"sep_Compute.jl")
using .R_sep_Model
using .C_sep_Model
using .Utils_Model
using .sep_Compute


using JuMP

export batch_model
default_sdir = dirname(dirname(@__FILE__))*"\\assets\\bounds\\"
""" Produces models for a set of constraints."""
function batch_model(t,ρdic,df;sdir = default_sdir,con_list= ["S∞sG","S₂sG","S₂₁sG"])
    for ex in keys(ρdic)
        ρ = ρdic[ex]
        d = filter(:Example => ==(ex),df).Size[1]
        for con ∈  con_list, isR ∈ ['R', 'C']
            (!filter(:Example => ==(ex),df).isReal[1] && (isR == 'R')) ? continue : 0
            println("Currently $isR-modeling example $ex with constraints $con")
            sep_mod = (isR == 'R') ? R_sep_Model.Modelξₜˢᵉᵖ(ρ,d,t;con_list=con) : C_sep_Model.Modelξₜˢᵉᵖ(ρ,d,t;con_list=con,noBlock = true)

            Utils_Model.export_model(sep_mod, sdir*"$(ex)_$(con)_$(isR)_$t.dat-s")
        end
    end
end


# """ Reads all .dat-s files and computes bounds storing result summary in a .csv file"""
function batch_Computeξₜˢᵉᵖ(boundsDir = default_sdir)
    datsFiles = [file for file in readdir(boundsDir,join = true) if contains(file,".dat-s")]
    file_loc = boundsDir*"Summary.csv"
    touch(file_loc)
    open(file_loc,"a") do io
        write(io, "|Example|Model|Constraint|Primal|Dual|obj_val|\n")
    end
    for file in datsFiles
        name = split(basename(file),".")[1]
        ex   = split(name,"_")[1]
        con  = split(name,"_")[2]
        num  = split(name,"_")[3]
        # try
        mod  = Utils_Model.read_model(file)
        omod = sep_Compute.Computeξₜˢᵉᵖ(mod)
        # catch
        #     open(file_loc,"a") do io; write(io, "|$ex|$con|err.|err.|err.| \n"); end
        #     continue
        # end
        pstat   = JuMP.primal_status(omod)
        dstat   = JuMP.dual_status(omod)
        ov_temp = round(JuMP.objective_value(omod),digits=3)
        ov = (string(pstat) == string(dstat) == "FEASIBLE_POINT") ?
             ov_temp :
             (contains(string(pstat)*string(dstat),"CER") ? "*" : NaN )

        open(file_loc,"a") do io
            write(io, "|$ex|$num|$con|$pstat|$dstat|$ov|\n")
        end
    end
end

using CSV
using DataFrames

function unstack_constraints(df,bDir = default_sdir)
    temp_df  = CSV.read("C:\\Users\\andries\\all-my-codes\\ju-sep-rank\\assets\\bounds\\"*"Summary.csv", DataFrame)
    temp2_df = select(innerjoin(df, temp_df, on = :Example),[:Example,:isReal,:isSeparable,:Model,:Size,:Bi_rank,:Constraint,:obj_val])
    temp4_df = unstack(temp2_df, :Model, :obj_val)

    temp4R_df = select(temp4_df,[:Example,:isReal,:isSeparable,:Size,:Bi_rank,:Constraint,:R])
    temp4Rus_df = unstack(temp4R_df, :Constraint, :R)
    rename!(temp4Rus_df,Dict(:S₂sG => :RS₂sG, :S₂₁sG => :RS₂₁sG, :S∞sG=>:RS∞sG));

    temp4C_df = select(temp4_df,[:Example,:isReal,:isSeparable,:Size,:Bi_rank,:Constraint,:C])
    temp4Cus_df = unstack(temp4C_df, :Constraint, :C)
    rename!(temp4Cus_df,Dict(:S₂sG => :CS₂sG, :S₂₁sG => :CS₂₁sG, :S∞sG=>:CS∞sG));

    temp5_df = innerjoin(temp4Cus_df,select(temp4Rus_df,[:Example,:RS₂sG, :RS₂₁sG, :RS∞sG]), on = :Example)

    CSV.write(bDir*"SummaryUnstacked.csv", temp5_df, delim="&")
    return temp5_df
end



end  # module
## Selecting all real examples
# df_real = filter(:isReal => ==(true),df)
# examples_real = Dict()
# for exa in df_real.ex
#     examples_real[exa] = examples_all[exa]
# end
# ##
# save_dir = sep_rank_proj_path*"\\PreRunModels\\R-models\\t=$(t[1])\\"
# Utils_Model.batch_model(t,examples_real,df_real,save_dir,["S₁sG","S₂sG","S₃sG"],true)
# sep_Compute.batch_Computeξₜˢᵉᵖ(save_dir)
# unstack_constraints(save_dir,df_real)
# ##
# temp_df  = CSV.read("C:\\Users\\andries\\all-my-codes\\ju-sep-rank\\PreRunModels\\C-models\\t=3\\SummaryUnstacked.csv", DataFrame, delim = "&")
# temp_dfr = CSV.read("C:\\Users\\andries\\all-my-codes\\ju-sep-rank\\PreRunModels\\R-models\\t=3\\SummaryUnstacked.csv", DataFrame, delim = "&")
#
# poes_df = select(innerjoin(temp_dfr, temp_df, on = :ex,makeunique=true),[:ex,:siz,:birank,:S₁sG,:S₂sG,:S₃sG,:S₁sG_1,:S₂sG_1,:S₃sG_1,])
# rename!(poes_df ,Dict(:S₁sG_1 => :RS₁sG, :S₂sG_1 => :RS₂sG, :S₃sG_1 => :RS₃sG))
# CSV.write("C:\\Users\\andries\\all-my-codes\\ju-sep-rank\\PreRunModels\\t=2.csv", poes_df, delim="&")

# function quick_run_spec_example(ex,t = (2,2),con_list = "S₁sG",isReal = true)
#     examples     = Examples.get_examples()
#     examples_all = merge(examples["ent"],examples["sep"])
#     df = Examples.get_example_overview(false)
#
#     ρ = examples_all[ex]
#     ρ_meta_data = filter(:ex => ==(ex),df)
#     d =  ρ_meta_data.siz[1]
#     if isReal
#         Lx, mod = R_sep_Model.Modelξₜˢᵉᵖ(ρ,d,t,con_list)
#     else
#         Lx, mod = C_sep_Model.Modelξₜˢᵉᵖ(ρ,d,t,con_list)
#     end
#     mod_opt = sep_Compute.Computeξₜˢᵉᵖ(mod)
#     return JuMP.values(Lx), mod_opt
# end
