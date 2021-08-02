module Utils_cons

export eᵢ,
       idx2var,
       idx2var_arr,
       idx2var_dic,
       idx2varxx̄ᵀtyȳᵀ,
       get_Gᴿ


"""The standard basis vector eᵢ in dimension n"""
eᵢ(n::Int,i::Int) = [j==i for j in 1:n]


"""Takes an array of exponents α's and gives array of same shape L(xᵅ)
Input: Array B of multi-indices: α
Output: L of x of array of indices: L(x.(B))
 """

idx2var(var,index) = isempty(index) ?  0 : map(α -> var[α],index)
idx2var(var,index,coef) = sum(coef .* idx2var(var, index))
idx2var(var,index,coef,g) = idx2var(var,g .+ index,coef)

make_gen_aff(arr) = reshape([arr...],size(arr)...)
idx2var_arr(var,index,coef) = make_gen_aff(map((x,y)-> Utils_cons.idx2var(var,x,y),index,coef))
idx2var_arr(var,index,coef,g) = make_gen_aff(map((x,y)-> Utils_cons.idx2var(var,x,y,g),index,coef))

function idx2var_dic(var,index,coef)
    kl = [keys(coef)...]
    f(b) = make_gen_aff(idx2var_arr(var, index[b],coef[b]))
    return Dict(zip(kl,map(b->f(b),kl)))
end
function idx2var_dic(var,index,coef,g)
    kl = [keys(coef)...]
    f(b) = make_gen_aff(idx2var_arr(var, index[b],coef[b],g))
    return Dict(zip(kl,map(b->f(b),kl)))
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
