# using Revise
    sep_rank_proj_path = dirname(dirname(@__FILE__))
    # Pkg.activate(sep_rank_proj_path)

    using LinearAlgebra
    using JuMP
    using MosekTools
    using DataFrames
    using CSV
    using Test

    srcDir  = sep_rank_proj_path*"\\src\\"
    testDir = sep_rank_proj_path*"\\test\\"

    include(srcDir*"Utils.jl")
    include(srcDir*"Examples.jl")
    include(srcDir*"Moments.jl")
    include(srcDir*"sep_Constraints.jl")
    include(srcDir*"sep_Model.jl")
    include(srcDir*"sep_Compute.jl")

    using .Utils
    using .Examples
    using .Moments
    using .sep_Constraints
    using .sep_Model
    using .sep_Compute


# include(testDir*"runTests.jl")



    ρ_dict = Examples.get_examples()
    ρ_sep  = ρ_dict["sep"]
    ρ_ent  = ρ_dict["ent"]
    # ρ_rand  = ρ_dict["rand"]
    ρ            = ρ_sep["sep4d3r4"]
    # model = Model()
    # list_of_keys = Moments.make_mom_expo_keys(n,t) # Define variables in the moment matrix.
    # @variable(model, Lx[list_of_keys] ) # ????

t = 2
sep_mod     = sep_Model.Modelξₜˢᵉᵖ(ρ,t,"S₂wGsG")
sep_mod_opt = sep_Compute.Computeξₜˢᵉᵖ(sep_mod )

#

for t in 6:5
    boundsDir = "C:\\Users\\andries\\Desktop\\Dissertation\\1.2. Projects\\1. Separable rank\\Bounds\\t=$t\\"
    boundsDirSep = boundsDir*"sep\\"
    sep_Model.batch_model(t,ρ_sep,boundsDirSep)
    sep_Compute.mass_read_comp(boundsDirSep)

    boundsDirEnt = boundsDir*"ent\\"
    sep_Model.batch_model(t,ρ_ent,boundsDirEnt)
    mass_read_comp(boundsDirEnt)
end


# ## Modeling phase
# if true
#     ρ         = Examples.get_examples()
#     t = 3
#     for state_type in ["rand"] # ["sep","ent"]#,
#         Data_dir = "C:\\Users\\andries\\Desktop\\Dissertation\\1.2. Projects\\1. Separable rank\\Bounds\\t=$t\\$state_type\\"
#         batch_run.batch_model(t,ρ[state_type],Data_dir)
#         batch_run.batch_read_run(Data_dir)
#     end
# end
#
## Post processing

# csv_path = "C:\\Users\\andries\\Desktop\\Dissertation\\1.2. Projects\\1. Separable rank\\Bounds\\t=2\\ent\\entt=2.csv"
# bounds_dir = "C:\\Users\\andries\\Desktop\\Dissertation\\1.2. Projects\\1. Separable rank\\Bounds\\"
# tiepe = "sep"
# t = 2


function proc_csv(csv_path)
    raw_str = read(csv_path , String)
    io = IOBuffer(raw_str)
    df = CSV.File(io,
                   delim='|',
                   missingstring="NA") |>
          DataFrame
    return df
end

function make_df_ex_dict(df)
    df_ex_dict = Dict()
    for ex in unique(df.ex)
        df_ex_dict[ex] = df[df.ex .== ex,:]
    end
    return df_ex_dict
end

function make_ordered_ov_list(df_dict)
    ov_list = Dict()
    for ex in keys(df_dict)
        df_temp   = df_dict[ex]
        rank      = df_dict[ex].rank[1]
        ov_list[ex]  = Any[]
        push!(ov_list[ex], ex,rank)
        df_ka = df_dict[ex]

        for con in ["S₁","S₂","S₃","S₁wG","S₂wG","S₃wG","S₁sG","S₂sG","S₃sG"]
            row = df_ka[df_ka.con .== con,:]
            P    = row.Primal[1]
            D    = row.Dual[1]
            ov   = round(row.obj_val[1], digits = 4)
            con  = row.con[1]

            if P == D == "FEASIBLE_POINT" && ov ≤ rank
                push!(ov_list[ex], ov)
            elseif P == "INFEASIBILITY_CERTIFICATE"  || D ==  "INFEASIBILITY_CERTIFICATE" || ov > rank
                if P == "INFEASIBILITY_CERTIFICATE"
                    push!(ov_list[ex],"p-INF.c")
                elseif D ==  "INFEASIBILITY_CERTIFICATE"
                    push!(ov_list[ex],"d-INF.c")
                else ov > rank
                    push!(ov_list[ex],"r-INF.c")
                end
            else
                push!(ov_list[ex],"U.R.S")
            end
        end
    end
    return ov_list
end

function proc_df_dict(df)
    df_dict = make_df_ex_dict(df)
    poes = make_ordered_ov_list(df_dict)
    df_new = DataFrame( ex = String[], rank = Int[],
                        S₁ = Any[], S₂ = Any[], S₃ = Any[],
                        S₁wG = Any[], S₂wG = Any[], S₃wG = Any[],
                        S₁sG = Any[], S₂sG = Any[], S₃sG = Any[])

    for ex in keys(poes)
        push!(df_new,poes[ex])
    end
    sort!(df_new,:ex)
    return df_new
end


# for tiepe in ["sep","ent"] # ,
#     for t = 2:5
tiepe = "ent"
t = 4
csv_path = bounds_dir*"t=$t\\$tiepe\\Summary.csv"


df = proc_csv(csv_path)
df_ess = df[:,[:ex, :rank, :con, :Primal, :Dual, :obj_val]]

df_new = proc_df_dict(df_ess)

save_path = bounds_dir*"t=$t\\$tiepe\\Summary_clean.csv"
CSV.write(save_path, df_new , delim = "|")
#     end
# end


# Data_dir = "C:\\Users\\andries\\Desktop\\Dissertation\\1.2. Projects\\1. Separable rank\\Bounds\\t=2\\ent\\"




# files_in_dir = read_n_proc_tools.collect_files_in_dir(Data_dir,".dat-s");
# file         = files_in_dir[11];
# file_name    = read_n_proc_tools.extract_file_name(file);
# t,ex,con        = read_n_proc_tools.extract_sep_bound_file_meta(file_name);
#
# # file_loc = Data_dir*"t=$t summary.csv"
#
# mod          = sep_model.read_dat_s_model(file)
# sep_mod_opt  = Computeξₜˢᵉᵖ(mod)
#
# pstat = primal_status(sep_mod_opt)
# dstat =  dual_status(sep_mod_opt)
# ov =  objective_value(sep_mod_opt)
#
# file_loc = Data_dir*"narXXX420BlazeIt.txt"
# touch(file_loc)
# open(file_loc,"a") do io
#     write(io, "$ex,$con,$pstat,$dstat,$ov \n")
#     # write(io, con*",$key,test \n")
# end
#

#
# if true
#     t = 5
#     Dats_dir = save_dir*"t=$t\\"
#     # batch_model(t,ρ_sep, Dats_dir, ["S₁","S₂","S₃","S₁wG","S₂wG","S₃wG","S₁sG","S₂sG","S₃sG"])
#     # batch_model(t,ρ_ent, Dats_dir, ["S₁","S₂","S₃","S₁wG","S₂wG","S₃wG","S₁sG","S₂sG","S₃sG"])
#
#     # computing phase
#     batch_run.batch_read_run(Dats_dir*"sep\\")
#     # batch_run.batch_read_run(Dats_dir*"ent\\")
# end
#
#
#
# # Data cleaning phase
# function clean_csv_Andries1(df)
#     df_selected = Dict()
#     df_new = sort(filter(row -> row.con == "S₁", df),:ex)
#     df_new  = select(df_new,[:ex, :obj_val])
#     rename!(df_new,:obj_val => :S₁)
#     df_new.S₁ = round.(df_new.S₁, digits = 4)
#     for con in ["S₂","S₃","S₁wG","S₂wG","S₃wG","S₁sG","S₂sG","S₃sG"]
#         df_selected_temp =  sort(filter(row -> row.con == con, df),:ex)
#         df_selected_temp = select(df_selected_temp ,[:ex, :obj_val])
#         df_new[con]      = round.(df_selected_temp.obj_val, digits = 4)
#     end
#     return df_new
# end
#
# for t in 4:4
#     for ex in ["sep"]#, "ent"]
#         Dats_dir = save_dir*"t=$t\\"
#         csv_file = Dats_dir*"$ex\\t=$t summary.csv"
#         df = DataFrame(CSV.File(csv_file))
#
#         df_new = clean_csv_Andries1(df)
#
#         new_csv_file = Dats_dir*"$ex\\t=$t clean_summary.csv"
#         CSV.write(new_csv_file, df_new, delim = "|")
#     end
# end
#
#
# t = 4
# ex = "sep"
# Dats_dir = save_dir*"t=$t\\"
# csv_file = Dats_dir*"$ex\\t=$t summary.csv"
# df = DataFrame(CSV.File(csv_file))
#
#
#
# df = DataFrame(ex = String[], con = String[], feas = String[], obj_val = String[] )
#
#
# function search_lines_for(file_path)
#     for line in readlines(file_path)
#         if contains(line,"phase.value")
#             feas = split(line," = ")[end]
#         end
#
#         if contains(line,"objValPrimal")
#             objVal = split(line," = +")[end]
#             break
#         end
#     end
#     return feas, objVal
# end
#
# pdf_files = read_n_proc_tools.collect_files_in_dir(Dats_dir*"sep\\", "pdf")
# file_path = pdf_files[1]
#
#
#
# for pdf_file in pdf_files
#
#     pdf_name = split(split(pdf_file,"=")[end],".")[1]
#     name = split(pdf_name,"_")
#     ex = name[1][2:end]
#     con = name[end]
#
#
#     push!(df, [ex con feas objVal])
# end
#
#
#
#
#
#
#
#
# for line in readlines(file_path)
#     if contains(line,"phase.value")
#         println("phase.value found")
#         feas = split(line," = ")[end]
#     end
#
#     if contains(line,"objValPrimal")
#         println("objValPrimal found")
#         objVal = split(line," = +")[end]
#         break
#     end
# end
#
#
#
# function proc_sep_something(csv_file)
#     raw_str = read(csv_file, String)
#     io      = IOBuffer(raw_str)
#     df = CSV.File(io,
#                    delim=',',
#                    header = [:ex, :con, :Primal, :Dual, :obj_val],
#                    missingstring="NA") |>
#           DataFrame
#
#     ex_list = ["sep1","sep2","sep3","sep4","sep5","sep6"]
#     con_list = ["S₁","S₂","S₃","S₁wG","S₂wG","S₃wG","S₁sG","S₂sG","S₃sG"]
#     mat = zeros(length(ex_list),length(con_list))
#
#     for row in eachrow(df)
#         row_ind = findfirst(row.ex .== ex_list)
#         col_ind = findfirst(row.con .== con_list)
#         try
#             mat[row_ind, col_ind] =  round(parse(Float64, row.obj_val),digits = 4)
#         catch
#             mat[row_ind, col_ind] = NaN
#         end
#     end
#
#     df_new = DataFrame(mat, con_list)
#     df_new.ex = ex_list
#     return df_new
# end
#
#
# """
# Returns all files in a specified directory that contains an identifier in the name.
# """
# function collect_files_in_dir(dir_loc, itentifier =".txt")
#     files = readdir(dir_loc, join=true)
#     files = [file for file in files if contains(file,itentifier)]
#     return files
# end
#
#
#
#
#
#
#
#
#
#
#
#
#
# # cd(dirname(@__FILE__))
# #     sourceDir = dirname(@__FILE__)*"\\"
# #     include(sourceDir*"data.jl")
# #     include(sourceDir*"moments.jl")
# #     using .Examples
# #     using .moments
# #     using LinearAlgebra
# #
# #
# # d = 2
# # n = 2*d
# # t = 2
# # isLeq = false
# #
# # ρ =  [0 0 0 0
# #       0 2 0 0
# #       0 0 0 0
# #       0 0 0 3]
# #
# #
# # mon_expo_mat  = make_mon_expo_mat(n,t,false)
# # xxᵀ_tens_yyᵀ  = make_xxᵀ_tens_yyᵀ(d)
# # mom_expo_keys = make_mom_expo_keys(n,t)
# #
# # mon_expo_mat_redu  = make_mon_expo_mat(n,t,ρ,false)
# # xxᵀ_tens_yyᵀ_redu  = make_xxᵀ_tens_yyᵀ(d,ρ)
# # mom_expo_keys_redu = make_mom_expo_keys(n,t,ρ)
#
#
#
#
#
#
# # nar = CSV.File(csv_file)
# # df = DataFrame(nar, header = [:ex :con :Primal :Dual :obj_val])
#
# # file_txt = readlines(csv_file)
#
#
# # mon_expo = make_mon_expo(n,t, isLeq)
# # list_of_keys = make_mom_expo_keys(n, t) # Define variables in the moment matrix.
# #
# #
# #
# # mom_matₜ_expo = make_mon_expo_mat(n,t,true)
# #
# # make_loc_cons_S₁(ρ,t,d,Lx)
# # make_loc_cons_S₂(ρ,t,d,Lx)
# # make_weakG_con(ρ,t,d,Lx)
# #  make_G_con(ρ,t,d,Lx)
# #
#
#
# #
# # using DataFrames
# # using CSV
# # using Pipe
# #
# #
# # save_dir = "C:\\Users\\andries\\Desktop\\Dissertation\\1.2. Projects\\1. Separable rank\\Bounds\\"
# # csv_file = save_dir*"t=2\\t=2 summary.csv"
# # # nar = CSV.File(csv_file)
# # # df = DataFrame(nar, header = [:ex :con :Primal :Dual :obj_val])
# #
# # # file_txt = readlines(csv_file)
# #
#
# #
# # csv_file_2 = save_dir*"t=2\\t=2 summary.csv"
# # df_2 = proc_sep_something(csv_file_2)
# #
# # csv_file_3 = save_dir*"t=3\\t=3 summary.csv"
# # df_3 = proc_sep_something(csv_file_3)
# #
# # csv_file_4 = save_dir*"t=4\\t=4 summary.csv"
# # df_4 = proc_sep_something(csv_file_4)
# #
# #
# # CSV.write("asdfasdf.txt", df_4, delim = "|")
#
# # include("C:\\Users\\andries\\all-my-codes\\ju-text-manipulation-tools\\read_n_proc_tools.jl")
# # using .read_n_proc_tools
# #
# # files_in_dir = read_n_proc_tools.collect_files_in_dir(Dats_dir, itentifier =".dat-s")
# # d = 2
# # r = 2
# # ρ = gen_rand_state(d,r,"no")
# # # ρ = gen_rand_state_int(d,r,"no")
# # # ρ = ρ/minimum(ρ)
# # ρ = round.(ρ,digits = 4)
# #
# # D = round.(diagm(rand(d^2)),digits = 4)
# # ρ_D = D*ρ*transpose(D)
# # ρ_D = round.(ρ_D,digits = 4)
# #
# # ρ_D_pt1 = partial_transpose_per(ρ_D , 1, [3,3])
# # ρ_D_pt2 = partial_transpose_per(ρ_D , 2, [3,3])
# #
# #
# #
# # ishermitian(ρ)
# # # ρ_evals = eigvals(ρ)
# # ishermitian(ρ_D)
# # issymmetric(ρ_D)
# # Symmetric(ρ_D) - ρ_D
# # # ρ_D_evals = eigvals(ρ_D)
# # ishermitian(ρ_D_pt1)
# # ishermitian(ρ_D_pt2)
#
#
# #
# # t = 4
# # save_dir  = dirname(@__FILE__)[1:end-3]*"DAT-s\\"
# # batch_model(t,ρ_sep, save_dir, ["S₁sG","S₂sG"])
#
#
# ##
# if false
#     for norma in ["tr"] #,"tr",
#         for t in 3:3
#             ρ_sep = get_sep_example(norma)
#             run_batch(t,ρ_sep, norma*"_sep_t=$t")
#
#             # ρ_rand         = generate_random_states(2:4,2:9,norma)
#             # run_batch(t, ρ_rand, norma*"_rand_t=$t")
#         end
#     end
# end
#
#
# #
# # ρ_rand_dg =  generate_random_states(2:4,2:9,"dg")
# # ρ_rand_no =  generate_random_states(2:4,2:9,"no")
# #
# # ex_keys = [key  for  key in  keys(ρ_rand_no) ]
# # ρ_dg = ρ_rand_dg[ex_keys[5]]
# # ρ_no = ρ_rand_no[ex_keys[5]]
# #
# #
# # t = 2
# # sep_mod_dg = Modelξₜˢᵉᵖ(ρ_dg,t, "S₁")
# # sep_mod_opt_dg=  Computeξₜˢᵉᵖ(sep_mod_dg)
# #
# # sep_mod_no = Modelξₜˢᵉᵖ(ρ_no,t, "S₁")
# # sep_mod_opt_no=  Computeξₜˢᵉᵖ(sep_mod_no)
#
# # ρ_sep = get_sep_example()
# # run_batch(t,ρ_sep,"sep_t=$t")
# #
# #
# # for t in  2:2
# #     ρ_sep = get_sep_example()
# #     run_batch(t,ρ_sep,"sep_t=$t")
# #
# #     ρ_rand  = generate_random_states(2:4,2:9,343)
# #     run_batch(t,ρ_rand,"rand_t=$t")
# #
# #     # ρ_ent = get_ent_example(d)
# #     # run_batch(t,ρ_ent,"sep t = $t")
# # end
#
# # ρ =  ρ_sep[(3, 5, "s1")]
# #
# # sep_mod = Modelξₜˢᵉᵖ(ρ,t,3, "S₁")
# # sep_mod_opt =  Computeξₜˢᵉᵖ(sep_mod)
#
# # ρ = ρ_sep[(2, 2, "s")]
# # ρ_sep[(3, 3, "s")]
# # ρ_sep[(3, 4, "s")]
# # ρ_sep[(3, 5, "s1")]
# # ρ_sep[(3, 5, "s2")]
#
#
# # These are also separable...yo!!!!!
#
#
#
#
#
# # ρ = ρ_dict[3, 2, "s"]
# # sep_mod = Modelξₜˢᵉᵖ(ρ,t,d, "S₁")
# # sep_mod_opt =  Computeξₜˢᵉᵖ(sep_mod)
