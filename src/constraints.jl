module constraints
export index_to_var,
       make_loc_cons,
       make_weakG_con,
       make_G_con

include("C:\\Users\\andries\\all-my-codes\\ju-sep-rank\\src\\moments.jl")
using .moments


"""Takes an array of exponents α's and gives array of same shape L(xᵅ)
Input: Array B of multi-indices: α
Output: L of x of array of indices: L(x.(B))
 """
function index_to_var(var, index_array)
    sub = α -> var[α]
    var_array = sub.(index_array)
    return var_array
end

"""
L ≥ 0 on M₂ₜ(S_ρ)
⟺
L(g⋅[x,y]ₜ₋₁[x,y]ₜ₋₁ᵀ) ⪰ 0 for
g = √(ρₘₐₓ) - xᵢ² , √(ρₘₐₓ) - yᵢ² and i ∈ [d]
"""
function make_loc_cons(ρ,t::Int64,d::Int64,Lx)
    n       = 2*d
    MB      = make_mon_expo_mat(n,t-1)
    sqr_ρ   = sqrt(maximum(ρ))
    s       = size(MB)[1]
    make_loc_con(eᵢ) = sqr_ρ*index_to_var(Lx, MB) - index_to_var(Lx, MB + repeat([2*eᵢ], s, s))
    loc_con = Dict()
    for k in 1:d
            eₖ = get_std_base_vec(n,k)
            eₖ₊d = get_std_base_vec(n,k+d)
            # Constraint: off diagonal of L ≧ 0 on 𝑀(Sᶜᵖ)   : L((√(ρₘₐₓ) - xᵢ²))⋅[x,y]ₜ₋₁[x,y]ₜ₋₁ᵀ) ⪰ 0
            loc_con[(k,1)] = make_loc_con(eₖ)
            # Constraint: off diagonal of L ≧ 0 on 𝑀(Sᶜᵖ)   : L((√(ρₘₐₓ) - yᵢ²))⋅[x,y]ₜ₋₁[x,y]ₜ₋₁ᵀ) ⪰ 0
            loc_con[(k,2)] = make_loc_con(eₖ₊d)
    end
    return loc_con
end

"""
Input: A(data matrix),t(Integer),Lx(JuMP variable)
Output: L((M-([x]₌₁[x]₌₁ᵀ))⊗([x]₌ₗ[x]₌ₗᵀ)))for l ∈ 0,1,...,t-1.
= M⊗L([x]₌ₗ[x]₌ₗᵀ) - L(([x]₌₁[x]₌₁ᵀ)⊗([x]₌ₗ[x]₌ₗᵀ))
"""
function make_weakG_con(ρ,t::Int64,d::Int64,Lx)
    n = 2*d
    weakG_con = Dict()
    xxᵀ_tens_yyᵀ         = make_xxᵀ_tens_yyᵀ(d)            # exponents of xxᵀ⊗yyᵀ
    for ℓ in 1:(t-2)
        LMBexp_ℓ          = make_mon_expo_mat(n,ℓ,false)   # exponents of [x,y]₌ₗ[x,y]₌ₗᵀ
        LMBexp_1ℓ         = var_kron(xxᵀ_tens_yyᵀ,LMBexp_ℓ)# exponents of(xxᵀ⊗yyᵀ)⊗([x,y]₌ₗ[x,y]₌ₗᵀ)
        LMB_ℓ             = index_to_var(Lx,LMBexp_ℓ)      # L([x,y]₌ₗ[x,y]₌ₗᵀ)
        RTerm             = index_to_var(Lx,LMBexp_1ℓ)     # L((xxᵀ⊗yyᵀ)⊗([x,y]₌ₗ[x,y]₌ₗᵀ))

        Lterm             = kron(ρ,LMB_ℓ)                  # ρ⊗L([x,y]₌ₗ[x,y]₌ₗᵀ)
        weakG_con[ℓ]      = Lterm - Rterm                  # ρ⊗L([x]₌ₗ[x]₌ₗᵀ) - L((xxᵀ⊗yyᵀ)⊗([x,y]₌ₗ[x,y]₌ₗᵀ)) ⪰ 0, ℓ ∈ 0,1,t-deg(g)/2
    end
    return weakG_con
end

"""M(Gρ ⊗ L) ⪰ 0 constraints
Where: Gρ := ρ - xxᵀ⊗yyᵀ
M(Gρ ⊗ L) = L(Gρ ⊗ [x, y]ₜ₋₂[x, y]ᵀₜ₋₂)
output: ρ⊗L([x, y]ₜ₋₂[x, y]ᵀₜ₋₂) - L( (xxᵀ⊗yyᵀ) ⊗ ([x, y]ₜ₋₂[x, y]ᵀₜ₋₂) )
"""
function make_G_con(ρ,t::Int64,d::Int64,Lx)
    n = 2*d
    xxᵀ_tens_yyᵀ       = make_xxᵀ_tens_yyᵀ(d)          # exponents of xxᵀ⊗yyᵀ
    LMBexpₜ₋₂          = make_mon_expo_mat(n,t-2,true) # exponents of [x, y]ₜ₋₂[x, y]ᵀₜ₋₂
    LMBₜ₋₂             = index_to_var(Lx,LMBexpₜ₋₂)    # L([x, y]ₜ₋₂[x, y]ᵀₜ₋₂)
    Lterm              = kron(ρ,LMBₜ₋₂)               # ρ⊗L([x, y]ₜ₋₂[x, y]ᵀₜ₋₂)

    LMBexp_1ₜ₋₂        = var_kron(xxᵀ_tens_yyᵀ,LMBexpₜ₋₂)  # exponents of (xxᵀ⊗yyᵀ) ⊗ ([x, y]ₜ₋₂[x, y]ᵀₜ₋₂)
    Rterm              = index_to_var(Lx,LMBexp_1ₜ₋₂)    # L( (xxᵀ⊗yyᵀ) ⊗ ([x, y]ₜ₋₂[x, y]ᵀₜ₋₂) )

    G_con = Lterm - Rterm                              # ρ⊗L([x, y]ₜ₋₂[x, y]ᵀₜ₋₂) - L( (xxᵀ⊗yyᵀ) ⊗ ([x, y]ₜ₋₂[x, y]ᵀₜ₋₂) )
    return G_con
end



end
