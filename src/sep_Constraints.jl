module sep_Constraints

using LinearAlgebra

srcDir = dirname(@__FILE__)*"\\"
include(srcDir*"Utils.jl")
include(srcDir*"Moments.jl")
using .Utils
using .Moments

export make_loc_cons_S₁,
       make_loc_cons_S₂,
       make_loc_cons_S₃,
       make_weakG_con,
       make_G_con

"""
L ≥ 0 on M₂ₜ(S_ρ¹)
⟺
L(g⋅[x,y]ₜ₋₁[x,y]ₜ₋₁ᵀ) ⪰ 0 for
g = √(ρₘₐₓ) - xᵢ² , √(ρₘₐₓ) - yᵢ² and i ∈ [d]
"""
function make_loc_cons_S₁(ρ,t::Int64,d::Int64,Lx)
    n       = 2*d
    MB      = make_mon_expo_mat(n,t-1)
    sqr_ρ   = sqrt(maximum(ρ))
    s       = size(MB)[1]
    make_loc_con(eᵢ) = sqr_ρ*Utils.index_to_var(Lx, MB) - Utils.index_to_var(Lx, MB + repeat([2*eᵢ], s, s))
    loc_con = Dict()
    for k in 1:d
            eₖ = Utils.eᵢ(n,k)
            eₖ₊d = Utils.eᵢ(n,k+d)
            # Constraint:  L ≧ 0 on 𝑀(Sᶜᵖ)   : L((√(ρₘₐₓ) - xᵢ²))⋅[x,y]ₜ₋₁[x,y]ₜ₋₁ᵀ) ⪰ 0
            loc_con[(k,1)] = make_loc_con(eₖ)
            # Constraint:  L ≧ 0 on 𝑀(Sᶜᵖ)   : L((√(ρₘₐₓ) - yᵢ²))⋅[x,y]ₜ₋₁[x,y]ₜ₋₁ᵀ) ⪰ 0
            loc_con[(k,2)] = make_loc_con(eₖ₊d)
    end
    return loc_con
end

"""
L ≥ 0 on M₂ₜ(S_ρ²)
⟺
L(g⋅[x,y]ₜ₋₁[x,y]ₜ₋₁ᵀ) ⪰ 0 for
g = √(Tr(ρ)) - ∑xᵢ² , √(Tr(ρ)) - ∑yᵢ²
"""
function make_loc_cons_S₂(ρ,t::Int64,d::Int64,Lx)
    n       = 2*d
    MB      = make_mon_expo_mat(n,t-1)
    sqrt_tr_ρ   = sqrt(tr(ρ))
    s       = size(MB)[1]
    # L((∑xᵢ²)⋅[x,y]ₜ₋₁[x,y]ₜ₋₁ᵀ )
    xRterm = sum([Utils.index_to_var(Lx, MB .+ [2* Utils.eᵢ(n,k)]) for k in 1:d ])
    # L((∑yᵢ²)⋅[x,y]ₜ₋₁[x,y]ₜ₋₁ᵀ )
    yRterm = sum([Utils.index_to_var(Lx, MB .+ [2* Utils.eᵢ(n,k+d)]) for k in 1:d ])
    #  √(ρₘₐₓ)⋅L([x,y]ₜ₋₁[x,y]ₜ₋₁ᵀ) - L((∑xᵢ²)⋅[x,y]ₜ₋₁[x,y]ₜ₋₁ᵀ )
    g₁ = sqrt_tr_ρ*Utils.index_to_var(Lx, MB) - xRterm
    #  √(ρₘₐₓ)⋅L([x,y]ₜ₋₁[x,y]ₜ₋₁ᵀ) - L((∑yᵢ²)⋅[x,y]ₜ₋₁[x,y]ₜ₋₁ᵀ )
    g₂ = sqrt_tr_ρ*Utils.index_to_var(Lx, MB) - yRterm
    loc_con = Dict(1 => g₁, 2 => g₂)
    return loc_con
end


"""
L ≥ 0 on M₂ₜ(S_ρ³)
⟺
L((Tr(ρ) - ∑xᵢ²)⋅[x,y]ₜ₋₁[x,y]ₜ₋₁ᵀ) ⪰ 0
L((sqrt(∑yᵢ²) - 1)⋅[x,y]ₜ₋₁[x,y]ₜ₋₁ᵀ) = 0
"""
function make_loc_cons_S₃(ρ,t::Int64,d::Int64,Lx)
    n       = 2*d
    MB      = make_mon_expo_mat(n,t-1)
    tr_ρ    = tr(ρ)
    s       = size(MB)[1]
    MBLx    = Utils.index_to_var(Lx, MB)
    # L((∑xᵢ²)⋅[x,y]ₜ₋₁[x,y]ₜ₋₁ᵀ )
    xRterm = sum([Utils.index_to_var(Lx, MB .+ [2* Utils.eᵢ(n,k)]) for k in 1:d ])
    #  tr(ρ)⋅L([x,y]ₜ₋₁[x,y]ₜ₋₁ᵀ) - L((∑xᵢ²)⋅[x,y]ₜ₋₁[x,y]ₜ₋₁ᵀ ) ⪰ 0
    g₁ = tr_ρ*MBLx - xRterm
    loc_con = Dict(1 => g₁)

    # L((∑yᵢ²)⋅[x,y]ₜ₋₁[x,y]ₜ₋₁ᵀ )
    yRterm = sum([Utils.index_to_var(Lx, MB .+ [2* Utils.eᵢ(n,k+d)]) for k in 1:d ])
    # L([x,y]ₜ₋₁[x,y]ₜ₋₁ᵀ) - L((∑yᵢ²)⋅[x,y]ₜ₋₁[x,y]ₜ₋₁ᵀ ) = 0
    g₂ = MBLx - yRterm

    return loc_con, g₂
end


"""
Input: A(data matrix),t(Integer),Lx(JuMP variable)
Output: L((ρ - ([x,y]₌₁[x,y]₌₁ᵀ)) ⊗ ([x,y]₌ₗ[x,y]₌ₗᵀ)))for l ∈ 1,...,t-2.
= ρ⊗L([x,y]₌ₗ[x,y]₌ₗᵀ) - L(([x,y]₌₁[x,y]₌₁ᵀ)⊗([x,y]₌ₗ[x,y]₌ₗᵀ))
"""
function make_weakG_con(ρ,t::Int64,d::Int64,Lx)
    n = 2*d
    weakG_con = Dict()
    xxᵀ_tens_yyᵀ         = make_xxᵀ_tens_yyᵀ(d)            # exponents of xxᵀ⊗yyᵀ
    for ℓ in 1:(t-2)
        LMBexp_ℓ          = make_mon_expo_mat(n,ℓ,false)   # exponents of [x,y]₌ₗ[x,y]₌ₗᵀ
        LMBexp_1ℓ         = var_kron(xxᵀ_tens_yyᵀ,LMBexp_ℓ)# exponents of  (xxᵀ⊗yyᵀ)⊗([x,y]₌ₗ[x,y]₌ₗᵀ)
        LMB_ℓ             = Utils.index_to_var(Lx,LMBexp_ℓ)      # L([x,y]₌ₗ[x,y]₌ₗᵀ)
        Rterm             = Utils.index_to_var(Lx,LMBexp_1ℓ)     # L((xxᵀ⊗yyᵀ)⊗([x,y]₌ₗ[x,y]₌ₗᵀ))

        Lterm             = kron(ρ,LMB_ℓ)                  # ρ⊗L([x,y]₌ₗ[x,y]₌ₗᵀ)
        weakG_con[ℓ]      = Lterm - Rterm                  # ρ⊗L([x,y]₌ₗ[x,y]₌ₗᵀ) - L((xxᵀ⊗yyᵀ)⊗([x,y]₌ₗ[x,y]₌ₗᵀ)) ⪰ 0, ℓ ∈ ,1,t-deg(G)/2
    end
    return weakG_con
end

"""M(Gρ ⊗ L) ⪰ 0 constraints
Where: Gρ := ρ - xxᵀ⊗yyᵀ
M(Gρ ⊗ L) = L(Gρ ⊗ [x, y]ₜ₋₂[x, y]ᵀₜ₋₂) ⪰ 0
output: ρ⊗L([x, y]ₜ₋₂[x, y]ᵀₜ₋₂) - L( (xxᵀ⊗yyᵀ) ⊗ ([x, y]ₜ₋₂[x, y]ᵀₜ₋₂) )
"""
function make_G_con(ρ,t::Int64,d::Int64,Lx)
    n = 2*d
    xxᵀ_tens_yyᵀ       = make_xxᵀ_tens_yyᵀ(d)          # exponents of xxᵀ⊗yyᵀ
    LMBexpₜ₋₂          = make_mon_expo_mat(n,t-2,true) # exponents of [x, y]ₜ₋₂[x, y]ᵀₜ₋₂
    LMBₜ₋₂             = Utils.index_to_var(Lx,LMBexpₜ₋₂)    # L([x, y]ₜ₋₂[x, y]ᵀₜ₋₂)
    Lterm              = kron(ρ,LMBₜ₋₂)               # ρ⊗L([x, y]ₜ₋₂[x, y]ᵀₜ₋₂)

    LMBexp_1ₜ₋₂        = var_kron(xxᵀ_tens_yyᵀ,LMBexpₜ₋₂)  # exponents of (xxᵀ⊗yyᵀ) ⊗ ([x, y]ₜ₋₂[x, y]ᵀₜ₋₂)
    Rterm              = Utils.index_to_var(Lx,LMBexp_1ₜ₋₂)    # L( (xxᵀ⊗yyᵀ) ⊗ ([x, y]ₜ₋₂[x, y]ᵀₜ₋₂) )

    G_con = Lterm - Rterm                              # ρ⊗L([x, y]ₜ₋₂[x, y]ᵀₜ₋₂) - L( (xxᵀ⊗yyᵀ) ⊗ ([x, y]ₜ₋₂[x, y]ᵀₜ₋₂) )
    return G_con
end



end
