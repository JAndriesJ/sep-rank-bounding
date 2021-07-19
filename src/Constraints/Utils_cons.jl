module Utils_cons

export eᵢ,
       idx2var,
       idx2varxx̄ᵀtyȳᵀ,
       get_Gᴿ


"""The standard basis vector eᵢ in dimension n"""
eᵢ(n::Int,i::Int) = [j==i for j in 1:n]


"""Takes an array of exponents α's and gives array of same shape L(xᵅ)
Input: Array B of multi-indices: α
Output: L of x of array of indices: L(x.(B))
 """

idx2var(var, index_array) = map(α -> var[α],index_array)

idx2var(var, index_array::Array{Array{Int64,1},2},coef_array) = coef_array .* map(α -> var[α],index_array)

#  idx2var(var, index_dict::Dict{Any,Any},coef_dict::Dict{Any,Any})
function idx2var(var, index_dict,coef_dict)
    kl = [keys(coef_dict)...]
    k_1 = pop!(kl)
    Res = idx2var(var, index_dict[k_1],coef_dict[k_1])
    for k in kl
        Res += idx2var(var, index_dict[k],coef_dict[k])
    end
    return Res
end


## Complex side
"""Returns Lᴿ(ℜe(xx̄ᵀ⊗yȳᵀ)), Lᴿ(ℑm(xx̄ᵀ⊗yȳᵀ)) """
function idx2varxx̄ᵀtyȳᵀ(Lx,xx̄ᵀtyȳᵀ)
    T =  Dict()
    for r in ["real", "imag"]
        sleut = [keys(xx̄ᵀtyȳᵀ[r])...]
        T[r] =   idx2var(Lx,xx̄ᵀtyȳᵀ[r][popfirst!(sleut)])
        for k in sleut
            (k[1] == '+') ? T[r] += idx2var(Lx,xx̄ᵀtyȳᵀ[r][k]) : T[r] -= idx2var(Lx,xx̄ᵀtyȳᵀ[r][k])
        end
    end
    return T
end


""" ?"""
function get_Gᴿ(xx̄ᵀ_tens_yȳᵀ) #???????
    Rxx̄ᵀ_tens_yȳᵀ = xx̄ᵀ_tens_yȳᵀ["real"]
    Ixx̄ᵀ_tens_yȳᵀ = xx̄ᵀ_tens_yȳᵀ["imag"]
    rkeys = [keys(Rxx̄ᵀ_tens_yȳᵀ)...]
    ikeys = [keys(Ixx̄ᵀ_tens_yȳᵀ)...]

    Gℝ = Dict()
    sgn_mat = Dict()
    for k in 1:8
        Re =    xx̄ᵀ_tens_yȳᵀ["real"][rkeys[k]]
        Im =    xx̄ᵀ_tens_yȳᵀ["imag"][ikeys[k]]
        Gℝ[k] = vcat(hcat(Re,Im), hcat(Im,Re))

        rsgn = rkeys[k][1] == '+' ? 1 : -1
        isgn = ikeys[k][1] == '+' ? 1 : -1
        sgn_mat[k] = [rsgn -isgn; isgn rsgn]
    end
    return Gℝ,sgn_mat
end



end
