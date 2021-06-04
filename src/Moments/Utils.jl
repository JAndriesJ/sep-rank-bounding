module Utils

using LinearAlgebra
using CSV
using DataFrames

export eᵢ,
       var_kron,
       post_proc

"""The standard basis vector eᵢ in dimension n"""
eᵢ(n::Int,i::Int) = [j==i for j in 1:n]

"""
input: Dictionary:
    keys: i,j ∈ [n], n::Int
    values: square array A_{i,j}
output: Array A = (A_{i,j})_i,j
"""
function assemble_dict(dict_of_blocks)
    n = Int(sqrt(length(keys(dict_of_blocks))))
    if n == 1
        return dict_of_blocks[1,1]
    end
    block = []
    row_block = []
    for i in 1:n
        for j in 1:n
            if j == 1
                row_block = dict_of_blocks[i, j]
            else
                row_block = hcat(row_block, dict_of_blocks[i,j])
            end
        end

        if i == 1
            block = row_block
        elseif i > 1
            block = vcat(block, row_block)
        end

    end
    return block
end

"""
input: A,B (arrays of integer tupples)
output: exponent array of A ⊗ B
"""
function var_kron(A,B)
    n1,n2 = size(A)

    D = Dict()
    for i in 1:n1
        for j in 1:n2
            C = repeat( [A[i,j]] , inner = (1,1), outer = size(B))
            D[i,j] = C + B
        end
    end
    return assemble_dict(D)
end

end







# ## Post proccessing
# ## DO I NEED THIS SHIT ????
# function update_bounds()
#     ρ_dict = Examples.get_examples()
#     for t in 2:4
#     # t = 2
#     bounds_dir = "C:\\Users\\andries\\Desktop\\Dissertation\\1.2. Projects\\1. Separable rank\\Bounds\\"
#         for ex in ["rand"]
#     # ex = "sep"
#             save_dir = bounds_dir*"t=$t\\$ex\\"
#             sep_Model.batch_model(t,ρ_dict[ex],save_dir)
#             sep_Compute.mass_read_comp(save_dir)
#             csv_path = save_dir*"Summary.csv"
#             df       = Utils.post_proc(csv_path)
#         end
#     end
# end
#
#
# function proc_csv(csv_path)
#     raw_str = read(csv_path , String)
#     io = IOBuffer(raw_str)
#     df = CSV.File(io,
#                    delim='|',
#                    missingstring="NA") |>
#           DataFrame
#     return df
# end
#
# function make_df_ex_dict(df)
#     df_ex_dict = Dict()
#     for ex in unique(df.ex)
#         df_ex_dict[ex] = df[df.ex .== ex,:]
#     end
#     return df_ex_dict
# end
#
# function make_ordered_ov_list(df_dict)
#     ov_list = Dict()
#     for ex in keys(df_dict)
#         df_temp   = df_dict[ex]
#         rank      = df_dict[ex].rank[1]
#         ov_list[ex]  = Any[]
#         push!(ov_list[ex], ex,rank)
#         df_ka = df_dict[ex]
#
#         for con in ["S₁","S₂","S₃","S₁wG","S₂wG","S₃wG","S₁sG","S₂sG","S₃sG"]
#             row = df_ka[df_ka.con .== con,:]
#             P    = row.Primal[1]
#             D    = row.Dual[1]
#             ov   = round(row.obj_val[1], digits = 4)
#             con  = row.con[1]
#
#             if P == D == "FEASIBLE_POINT" && ov ≤ rank
#                 push!(ov_list[ex], ov)
#             elseif P == "INFEASIBILITY_CERTIFICATE"  || D ==  "INFEASIBILITY_CERTIFICATE" || ov > rank
#                 if P == "INFEASIBILITY_CERTIFICATE"
#                     push!(ov_list[ex],"p-INF.c")
#                 elseif D ==  "INFEASIBILITY_CERTIFICATE"
#                     push!(ov_list[ex],"d-INF.c")
#                 else ov > rank
#                     push!(ov_list[ex],"r-INF.c")
#                 end
#             else
#                 push!(ov_list[ex],"U.R.S.")
#             end
#         end
#     end
#     return ov_list
# end
#
# function post_proc(csv_path)
#     df      = proc_csv(csv_path)
#     df_dict = make_df_ex_dict(df)
#     poes    = make_ordered_ov_list(df_dict)
#     df_new  = DataFrame( ex = String[], rank = Int[],
#                         S₁ = Any[], S₂ = Any[], S₃ = Any[],
#                         S₁wG = Any[], S₂wG = Any[], S₃wG = Any[],
#                         S₁sG = Any[], S₂sG = Any[], S₃sG = Any[])
#
#     for ex in keys(poes)
#         push!(df_new,poes[ex])
#     end
#     sort!(df_new,:ex)
#     save_path = splitext(csv_path)[1]*"_clean.csv"
#     CSV.write(save_path, df_new , delim = "|")
#     return df_new
# end
