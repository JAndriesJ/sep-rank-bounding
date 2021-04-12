module Utils

export eᵢ,
       var_kron,
       index_to_var

"""
output: eₖ ∈ {0,1}ⁿ i.e. the standard basis vector
"""
eᵢ(n::Int,i::Int) = [j==i for j in 1:n]
# function get_std_base_vec(n::Int,k::Int)
#     @assert n >= k
#     eₖ = zeros(Int64,n)
#     eₖ[k] = 1
#     return eₖ
# end

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



"""Takes an array of exponents α's and gives array of same shape L(xᵅ)
Input: Array B of multi-indices: α
Output: L of x of array of indices: L(x.(B))
 """
function index_to_var(var, index_array)
    sub = α -> var[α]
    var_array = sub.(index_array)
    return var_array
end


end  # module Utilities
