module Utils_cons

export index_to_var,
       eᵢ,
       get_Lxx̄ᵀ_tens_yȳᵀ,
       get_Gℝ


"""The standard basis vector eᵢ in dimension n"""
eᵢ(n::Int,i::Int) = [j==i for j in 1:n]


"""Takes an array of exponents α's and gives array of same shape L(xᵅ)
Input: Array B of multi-indices: α
Output: L of x of array of indices: L(x.(B))
 """
function index_to_var(var, index_array)
    sub = α -> var[α]
    var_array = sub.(index_array)
    return var_array
end

## Complex side
"""Returns L^ℝ(ℜ(xx̄ᵀ⊗yȳᵀ)), L^ℝ(ℑ(xx̄ᵀ⊗yȳᵀ)) """
function get_Lxx̄ᵀ_tens_yȳᵀ(xx̄ᵀ_tens_yȳᵀ,Lx)
    T =  Dict()
    for k in ["real", "imag"]
        sleut = [key for key in keys(xx̄ᵀ_tens_yȳᵀ[k])]
        T[k] =   Utils_cons.index_to_var(Lx,xx̄ᵀ_tens_yȳᵀ[k][popfirst!(sleut)])
        for key in sleut
            key = sleut[1]
            if key[1] == '+'
                T[k] = T[k] + Utils_cons.index_to_var(Lx,xx̄ᵀ_tens_yȳᵀ[k][key])
            elseif key[1] == '-'
                T[k] = T[k] - Utils_cons.index_to_var(Lx,xx̄ᵀ_tens_yȳᵀ[k][key])
            end
        end
    end
    return T
end
""" ?"""
function get_Gℝ(xx̄ᵀ_tens_yȳᵀ) 
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
# """"""
# function get_L()
#
# end
