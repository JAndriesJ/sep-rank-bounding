module Moments

using LinearAlgebra

srcDir = dirname(@__FILE__)
include(srcDir*"\\Utils.jl")
using .Utils

export prop_zero_diags,
       make_mon_expo_mat,
       make_xxᵀ_tens_yyᵀ,
       make_mom_expo_keys

       # assemble_dict,
       # make_mon_expo_arr,
       # make_mon_expo,

"""
Input: ρ with possible zeros on the diagonal.
Output: ρ with rows/cols corresponding to zeros diagonals deleted.
"""
function prop_zero_diags(ρ)
    maks = findall(iszero.(diag(ρ)))
    return ρ[setdiff(1:end,maks), setdiff(1:end,maks)]
end


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

function make_mon_expo_mat(n::Int,t::Int,ρ, isLeq::Bool = true)
    d = Int(n/2)
    z_dags             = iszero.(diag(ρ))
    zero_diag          = diag(make_xxᵀ_tens_yyᵀ(d))[z_dags]

    mon_expo_mat       = make_mon_expo_mat(n,t,isLeq)
    row_col_purge_list = [ mom in zero_diag   for mom in diag(mon_expo_mat)]
    row_col_keep_list  = setdiff(1:size(mon_expo_mat)[1], findall(row_col_purge_list))

    return mon_expo_mat[row_col_keep_list, row_col_keep_list]
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

    Utils.var_kron(xxᵀ_expo,yyᵀ_expo)
end

function make_xxᵀ_tens_yyᵀ(d::Int64,ρ)
    xxᵀ_tens_yyᵀ = make_xxᵀ_tens_yyᵀ(d)
    maks = findall(iszero.(diag(ρ)))
    return  xxᵀ_tens_yyᵀ[setdiff(1:end,maks), setdiff(1:end,maks)]
end

"""
output: array: unique exponents in [x]≦ₜ[x]≦ₜᵀ γ ∈ N_2t^n, values are indeces in
                    moment matrix array of (α,β) ∈ (N_2t^n)^2 such that α + β = γ
"""
function make_mom_expo_keys(n::Int,t::Int)
    mon_vec = make_mon_expo(n,t)
    return unique( [ α + β for α in mon_vec, β in mon_vec ])
end

# make_mom_expo_keys(n::Int,t::Int,ρ) = unique(vec(make_mon_expo_mat(n,t,ρ)))

end
