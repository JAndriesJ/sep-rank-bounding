# using Revise
sep_rank_proj_path = dirname(dirname(@__FILE__))
srcDir  = sep_rank_proj_path*"\\src\\";
# testDir = sep_rank_proj_path*"\\test\\";

using JuMP

include(srcDir*"Examples\\Examples.jl")
include(srcDir*"Moments\\Moments.jl")
include(srcDir*"Constraints\\Utils_cons.jl")
include(srcDir*"Constraints\\C_constraints.jl")
include(srcDir*"Model\\C_sep_Model.jl")
include(srcDir*"Model\\Utils_Model.jl")
include(srcDir*"sep_Compute.jl")

using .Examples
using .Moments
using .Utils_cons
using .C_constraints
using .C_sep_Model
using .Utils_Model
using .sep_Compute




t = (3,3)
examples = Examples.get_examples()
save_dir = sep_rank_proj_path*"\\PreRunModels\\"

Utils_Model.batch_model(t,examples["ent"]["real"],save_dir,["S₁sG"],false)
Utils_Model.batch_model(t,examples["sep"]["real"],save_dir,["S₁sG"],false)
Utils_Model.batch_model(t,examples["ent"]["imag"],save_dir,["S₁sG"],false)
Utils_Model.batch_model(t,examples["sep"]["imag"],save_dir,["S₁sG"],false)

boundsDir = save_dir*"t=2\\ent\\real\\"
boundsDir = save_dir*"t=2\\sep\\real\\"
boundsDir = save_dir*"t=3\\ent\\real\\"
boundsDir = save_dir*"t=3\\sep\\real\\"
boundsDir = save_dir*"t=4\\ent\\real\\"
boundsDir = save_dir*"t=4\\sep\\real\\"

boundsDir = save_dir*"t=2\\ent\\imag\\"
boundsDir = save_dir*"t=2\\sep\\imag\\"
boundsDir = save_dir*"t=3\\ent\\imag\\"
boundsDir = save_dir*"t=3\\sep\\imag\\"


boundsDir = save_dir*"C-models\\t=2\\"
if true
    datsFiles = [file for file in readdir(boundsDir,join = true) if contains(file,".dat-s")]
    file_loc = boundsDir*"Summary.csv"
    touch(file_loc)
    open(file_loc,"a") do io
        write(io, "|ex|d|con|Primal|Dual|obj_val|\n")
    end
    for file in datsFiles
        mod = Utils_Model.read_model(file)
        sep_mod_opt = sep_Compute.Computeξₜˢᵉᵖ(mod)

        # ex,d,rank,con   =  extract_sep_bound_file_meta(basename(file))
        name = split(basename(file),".")[1]
        con = split(name,"_")[2]
        name = split(name,"_")[1]
        dstr = name[(end-9):end]
        ex = name[1:end-10]

        pstat = JuMP.primal_status(sep_mod_opt)
        dstat = JuMP.dual_status(sep_mod_opt)
        ov    = round(JuMP.objective_value(sep_mod_opt),digits=2)

        open(file_loc,"a") do io
            write(io, "|$ex|$dstr|$con|$pstat|$dstat|$ov| \n")
        end
    end
end






## Block diag ℂ
CMM = Moments.make_mon_expo_mat(n,(2,0))
kwartsum(arr) = [sum(arr[1:d[1]]),sum(arr[d[1]+1:2*d[1]]),sum(arr[(2*d[1]+1):(2*d[1]+d[2])]),sum(arr[(2*d[1]+d[2]+1):end])]
halfdiff(arr) = [arr[1]- arr[2],arr[3]- arr[4]]
ks = kwartsum.(CMM)
hf = halfdiff.(ks)

bind = unique(hf)

nar = Dict()
for b in bind
    temp1  = findall([ent == b for ent in hf])
    nar[b] = [t[1]  for t in temp1 ]
end


findall(nar[[0,2]])

fieldnames(CartesianIndex)

using Combinatorics
k = (2,2,0)
multinomial(k...)
