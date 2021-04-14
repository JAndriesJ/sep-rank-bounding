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
Simple utility function
"""
function get_stuff(ρ,t::Int64)
   d       = Int(sqrt(size(ρ)[1]))
   n       = 2*d
   MB      = Moments.make_mon_expo_mat_perm(n,t-1)
   return d,n,MB
end


"""
L ≥ 0 on M₂ₜ(S_ρ¹)
⟺
L(g⋅η) ⪰ 0 for
η ∈ even-degree-principle-submatrices of [x,y]ₜ₋₁[x,y]ₜ₋₁ᵀ
g = √(ρₘₐₓ) - xᵢ² , √(ρₘₐₓ) - yᵢ² and i ∈ [d]
"""
function make_loc_cons_S₁(ρ,t::Int64,Lx)
    d,n,MB = get_stuff(ρ,t)
    sqr_ρ   = sqrt(maximum(ρ))
    make_loc_con(eᵢ,key) = sqr_ρ*Utils.index_to_var(Lx, MB[key]) - Utils.index_to_var(Lx, MB[key] + repeat([2*eᵢ], size(MB[key])[1], size(MB[key])[1]))
    loc_con = Dict()
    for key in keys(MB)
        for k in 1:d
                eₖ = Utils.eᵢ(n,k)
                eₖ₊d = Utils.eᵢ(n,k+d)
                # Constraint:  L ≧ 0 on 𝑀(Sᶜᵖ)   : L((√(ρₘₐₓ) - xᵢ²))⋅[x,y]ₜ₋₁[x,y]ₜ₋₁ᵀ) ⪰ 0
                loc_con[(key,k,1)] = make_loc_con(eₖ,key)
                # Constraint:  L ≧ 0 on 𝑀(Sᶜᵖ)   : L((√(ρₘₐₓ) - yᵢ²))⋅[x,y]ₜ₋₁[x,y]ₜ₋₁ᵀ) ⪰ 0
                loc_con[(key,k,2)] = make_loc_con(eₖ₊d,key)
        end
    end
    return loc_con
end

"""
L ≥ 0 on M₂ₜ(S_ρ²)
⟺
L(g⋅η) ⪰ 0 for
η ∈ even-degree-principle-submatrices of [x,y]ₜ₋₁[x,y]ₜ₋₁ᵀ
g = √(Tr(ρ)) - ∑xᵢ² , √(Tr(ρ)) - ∑yᵢ²
"""
function make_loc_cons_S₂(ρ,t::Int64,Lx)
    d,n,MB = get_stuff(ρ,t)
    loc_con = Dict()
    for key in keys(MB)
        sqrt_tr_ρ   = sqrt(tr(ρ))
        # L((∑xᵢ²)⋅[x,y]ₜ₋₁[x,y]ₜ₋₁ᵀ )
        xRterm = sum([Utils.index_to_var(Lx, MB[key] .+ [2* Utils.eᵢ(n,k)]) for k in 1:d ])
        # L((∑yᵢ²)⋅[x,y]ₜ₋₁[x,y]ₜ₋₁ᵀ )
        yRterm = sum([Utils.index_to_var(Lx, MB[key] .+ [2* Utils.eᵢ(n,k+d)]) for k in 1:d ])
        #  √(ρₘₐₓ)⋅L([x,y]ₜ₋₁[x,y]ₜ₋₁ᵀ) - L((∑xᵢ²)⋅[x,y]ₜ₋₁[x,y]ₜ₋₁ᵀ )
        g₁ = sqrt_tr_ρ*Utils.index_to_var(Lx, MB[key]) - xRterm
        #  √(ρₘₐₓ)⋅L([x,y]ₜ₋₁[x,y]ₜ₋₁ᵀ) - L((∑yᵢ²)⋅[x,y]ₜ₋₁[x,y]ₜ₋₁ᵀ )
        g₂ = sqrt_tr_ρ*Utils.index_to_var(Lx, MB[key]) - yRterm
        loc_con[key,1] = g₁
        loc_con[key,2] = g₂
    end
    return loc_con
end

"""
L ≥ 0 on M₂ₜ(S_ρ³)
⟺
L((Tr(ρ) - ∑xᵢ²)⋅η) ⪰ 0
η ∈ even-degree-principle-submatrices of [x,y]ₜ₋₁[x,y]ₜ₋₁ᵀ
L((sqrt(∑yᵢ²) - 1)⋅[x,y]ₜ₋₁[x,y]ₜ₋₁ᵀ) = 0
"""
function make_loc_cons_S₃(ρ,t::Int64,Lx)
    d,n,MB = get_stuff(ρ,t)
    tr_ρ    = tr(ρ)

    loc_con = Dict()
    g₂ = Dict()
    for key in keys(MB)
        # L([x,y]ₜ₋₁[x,y]ₜ₋₁ᵀ)
        MBLx    = Utils.index_to_var(Lx, MB[key])
        # L((∑xᵢ²)⋅[x,y]ₜ₋₁[x,y]ₜ₋₁ᵀ )
        xRterm = sum([Utils.index_to_var(Lx, MB[key] .+ [2* Utils.eᵢ(n,k)]) for k in 1:d ])
        #  tr(ρ)⋅L([x,y]ₜ₋₁[x,y]ₜ₋₁ᵀ) - L((∑xᵢ²)⋅[x,y]ₜ₋₁[x,y]ₜ₋₁ᵀ ) ⪰ 0
        loc_con[key,1] =  tr_ρ*MBLx - xRterm

        # L((∑yᵢ²)⋅[x,y]ₜ₋₁[x,y]ₜ₋₁ᵀ )
        yRterm = sum([Utils.index_to_var(Lx, MB[key] .+ [2* Utils.eᵢ(n,k+d)]) for k in 1:d ])
        # L([x,y]ₜ₋₁[x,y]ₜ₋₁ᵀ) - L((∑yᵢ²)⋅[x,y]ₜ₋₁[x,y]ₜ₋₁ᵀ ) = 0
        g₂[key] = MBLx - yRterm
    end

    return loc_con, g₂
end


""" SOMETHING IS WRONG HERE!!!!!
Input: A(data matrix),t(Integer),Lx(JuMP variable)
Output:
L((ρ - ([x,y]₌₁[x,y]₌₁ᵀ)) ⊗ ([x,y]₌ₗ[x,y]₌ₗᵀ)))
for l ∈ 1,...,t-2.

= ρ⊗L([x,y]₌ₗ[x,y]₌ₗᵀ) - L(([x,y]₌₁[x,y]₌₁ᵀ)⊗([x,y]₌ₗ[x,y]₌ₗᵀ)) ⪰ 0
note

"""
function make_weakG_con(ρ,t::Int64,Lx)
    d = Int(sqrt(size(ρ)[1]))
    n = 2*d

    weakG_con = Dict()
    xxᵀ_tens_yyᵀ = Moments.make_xxᵀ_tens_yyᵀ(d)            # exponents of xxᵀ⊗yyᵀ
    for ℓ in 1:(t-2)
        LMBexp_ℓ = Moments.make_mon_expo_mat_perm(n,ℓ,false)   # exponents of [x,y]₌ₗ[x,y]₌ₗᵀ
        for key in keys(LMBexp_ℓ)
            LMBexp_1ℓ         = Moments.var_kron(xxᵀ_tens_yyᵀ,LMBexp_ℓ[key])  #  # exponents of  (xxᵀ⊗yyᵀ)⊗([x,y]₌ₗ[x,y]₌ₗᵀ)
            LMB_ℓ             = Utils.index_to_var(Lx,LMBexp_ℓ[key])  # L([x,y]₌ₗ[x,y]₌ₗᵀ)
            Rterm             = Utils.index_to_var(Lx,LMBexp_1ℓ) # L((xxᵀ⊗yyᵀ)⊗([x,y]₌ₗ[x,y]₌ₗᵀ))

            Lterm             = kron(ρ,LMB_ℓ)                  # ρ⊗L([x,y]₌ₗ[x,y]₌ₗᵀ)
            weakG_con[key,ℓ]      = Lterm - Rterm                  # ρ⊗L([x,y]₌ₗ[x,y]₌ₗᵀ) - L((xxᵀ⊗yyᵀ)⊗([x,y]₌ₗ[x,y]₌ₗᵀ)) ⪰ 0, ℓ ∈ ,1,t-deg(G)/2
        end
    end
    return weakG_con
end

"""M(Gρ ⊗ L) ⪰ 0 constraints
Where: Gρ := ρ - xxᵀ⊗yyᵀ
M(Gρ ⊗ L) = L(Gρ ⊗ [x, y]ₜ₋₂[x, y]ᵀₜ₋₂) ⪰ 0
output: ρ⊗L([x, y]ₜ₋₂[x, y]ᵀₜ₋₂) - L( (xxᵀ⊗yyᵀ) ⊗ ([x, y]ₜ₋₂[x, y]ᵀₜ₋₂) )
"""
function make_G_con(ρ,t::Int64,Lx)
    d = Int(sqrt(size(ρ)[1]))
    n = 2*d
    xxᵀ_tens_yyᵀ       = Moments.make_xxᵀ_tens_yyᵀ(d)          # exponents of xxᵀ⊗yyᵀ
    LMBexpₜ₋₂          = Moments.make_mon_expo_mat_perm(n,t-2,true) # exponents of [x, y]ₜ₋₂[x, y]ᵀₜ₋₂
    G_con = Dict()
    for key in keys(LMBexpₜ₋₂ )
        LMBₜ₋₂             = Utils.index_to_var(Lx,LMBexpₜ₋₂[key]) # L([x, y]ₜ₋₂[x, y]ᵀₜ₋₂)
        Lterm              = kron(ρ,LMBₜ₋₂)                        # ρ⊗L([x, y]ₜ₋₂[x, y]ᵀₜ₋₂)

        LMBexp_1ₜ₋₂        = Moments.var_kron(xxᵀ_tens_yyᵀ,LMBexpₜ₋₂[key]) # exponents of (xxᵀ⊗yyᵀ) ⊗ ([x, y]ₜ₋₂[x, y]ᵀₜ₋₂)
        Rterm              = Utils.index_to_var(Lx,LMBexp_1ₜ₋₂)    # L( (xxᵀ⊗yyᵀ) ⊗ ([x, y]ₜ₋₂[x, y]ᵀₜ₋₂) )

        G_con[key] = Lterm - Rterm                              # ρ⊗L([x, y]ₜ₋₂[x, y]ᵀₜ₋₂) - L( (xxᵀ⊗yyᵀ) ⊗ ([x, y]ₜ₋₂[x, y]ᵀₜ₋₂) )
    end
    return G_con
end


end
