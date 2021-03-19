module moments
using LinearAlgebra

export make_mon_expo_arr,
       make_mon_expo,
       get_std_base_vec,
       make_mon_expo_mat,
       make_xxᵀ_tens_yyᵀ,
       make_mom_expo_keys,
       assemble_dict,
       var_kron


"""
output: exponents α ∈ Nⁿₜ of [x]≦ₜ of [x]₌ₜ (array of integers)
"""
function make_mon_expo_arr(n::Int,t::Int, isLeq::Bool = true)
    if t == 0 # [x]₌₀
        return zeros(Int32,1,n)
    else # [x]₌ₜ
        temp = make_mon_expo_arr(n,t-1,isLeq)
        e₁ = hcat(1, zeros(Int32,1,n-1))
        output = e₁ .+ temp
        for i = 1:n-1
            output = vcat(output,circshift(e₁,(0,i)) .+ temp)
        end

        if isLeq # [x]≦ₜ
            output = vcat(temp, output)
        end
        return unique(output,dims=1)
    end
end

"""
output: exponents α ∈ Nⁿₜ of [x]≦ₜ of [x]₌ₜ (array of arrays of integers)
"""
function make_mon_expo(n::Int,t::Int, isLeq::Bool = true)
    mon_expo_arr = make_mon_expo_arr(n,t,isLeq)
    mon_expo     = [r  for r in  eachrow(mon_expo_arr)]
    return mon_expo
end


"""
output: eₖ ∈ {0,1}ⁿ i.e. the standard basis vector
"""
function get_std_base_vec(n::Int,k::Int)
    @assert n >= k
    eₖ = zeros(Int64,n)
    eₖ[k] = 1
    return eₖ
end

"""
output: exponents α ∈ Nⁿₜ of [x]≦ₜ[x]≦ₜᵀ or [x]₌ₜ[x]₌ₜᵀ where , x = (x₁,x₂,...,xₙ)
"""
function make_mon_expo_mat(n::Int,t::Tuple{Int64,Int64},isLeq::Bool = true)
    mon1      = make_mon_expo(n,t[1], isLeq)
    mon2      = make_mon_expo(n,t[2], isLeq)

    mon1_mon2ᵀₜ_vec = [ mi+ mj for mi in mon1 for mj in mon2]
    mon1_mon2ᵀₜ    = reshape(mon1_mon2ᵀₜ_vec, (length(mon1), length(mon2)) )
    return mon1_mon2ᵀₜ
end
function make_mon_expo_mat(n::Int,t::Int,isLeq::Bool = true)
    xxᵀₜ     = make_mon_expo_mat(n,(t,t),isLeq)
    return xxᵀₜ
end

"""
returns the exponents of xxᵀ⊗yyᵀ
"""
function make_xxᵀ_tens_yyᵀ(d::Int64)
    n = 2*d
    pre_expo = make_mon_expo_mat(n,(1,0),false)
    x = pre_expo[1:d]
    y = pre_expo[d+1:n]

    xxᵀ_expo = x .+ reshape(x,1,d)
    yyᵀ_expo = y .+ reshape(y,1,d)

    var_kron(xxᵀ_expo,yyᵀ_expo)
end


"""
output: array: unique exponents in [x]≦ₜ[x]≦ₜᵀ γ ∈ N_2t^n, values are indeces in
                    moment matrix array of (α,β) ∈ (N_2t^n)^2 such that α + β = γ
"""
function make_mom_expo_keys(n::Int,t::Int)
    mon_vec = make_mon_expo(n,t)
    return unique( [ α + β for α in mon_vec, β in mon_vec ])
end

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
