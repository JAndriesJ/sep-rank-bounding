module R_constraints
    include(dirname(@__FILE__)*"\\Utils_cons.jl")
    include(dirname(dirname(@__FILE__))*"\\Moments\\Moments.jl")
    using .Utils_cons
    using .Moments
    const uc = Utils_cons
    const mom = Moments

    export make_mon_expo_keys,
           make_PSD_con,
           make_ord4_con,
           make_loc_cons_S_inf,
           make_loc_cons_S₂,
           make_loc_cons_S₂₁,
           make_G_con

    make_mon_expo_keys(n,t::Int) = Moments.make_mon_expo(n,t*2)

    """L([x,y]ₜ[x,y]ₜᵀ) ⪰ 0"""
    function make_PSD_con(d,t,Lx)
        MB = mom.get_ℝ_block_diag(d,t)
        Bs = [keys(MB)...]
        PSD_con = Dict(zip(Bs, map(b -> uc.idx2var(Lx, MB[b]),Bs)))
        return PSD_con
    end

    # """ρ = aaᵀ⊗bbᵀ , ρᵢⱼᵢⱼ = 0 ⟹ aᵢbⱼ = 0 ⟹ L(xᵢyⱼ) = 0 """
    #
    # function zeroprop(ρ,d,t,Lx)
    #     nair((i,j)) = vcat(mom.eᵢ(d[1],i),mom.eᵢ(d[2],j))
    #     zdiags =[(i,j) for i in 1:d[1], j in 1:d[2] if ρ[i+(j-1)*d[1],i+(j-1)*d[1]] == 0 ]
    #     return uc.idx2var(Lx, [k .* m  for m ∈  map(nair, zdiags) for k in vcat(1,3:t[1])])
    # end

    """L(xxᵀ⊗yyᵀ)"""
    make_ord4_con(d,Lx) = uc.idx2var(Lx,mom.make_xxᵀ_tens_yyᵀ(d))

    """ L(g⋅η) ⪰ 0 for η ∈ even-degree-principle-submatrices of [x,y]ₜ₋₁[x,y]ₜ₋₁ᵀ
                        g = √(ρₘₐₓ) - xᵢ² , √(ρₘₐₓ) - yᵢ² and i ∈ [d] """
    function make_loc_con(eᵢ,sqr_ρ,B,Lx)
        Lxᵢ = uc.idx2var(Lx, B + repeat([2*eᵢ], size(B)...)) #  L(xᵢ²⋅η)
        Lρ = sqr_ρ*uc.idx2var(Lx, B)  # L(√ρₘₐₓ⋅η)
        return  Lρ - Lxᵢ
    end

    function make_loc_cons_S_inf(ρ,d,t,Lx)
        n     = sum(d)
        MB    = mom.get_ℝ_block_diag(d,t.- 1)
        sqr_ρ = sqrt(maximum(ρ))

        loc_con = Dict()
        for b in keys(MB)
            for k in 1:d[1]
                loc_con[(b,"x²_$k")] = make_loc_con(uc.eᵢ(n,k),sqr_ρ,MB[b],Lx) # Constraint:  L((√(ρₘₐₓ) - xᵢ²))⋅η) ⪰ 0
            end
            for k in 1:d[2]
                loc_con[(b,"y²_$k")] = make_loc_con(uc.eᵢ(n,k+d[1]),sqr_ρ,MB[b],Lx) # Constraint:  L((√(ρₘₐₓ) - yᵢ²))⋅η) ⪰ 0
            end
        end
        return loc_con
    end

    """L(g⋅η) ⪰ 0 ∀ η ∈ even-degree-principle-submatrices of [x,y]ₜ₋₁[x,y]ₜ₋₁ᵀ
                       g = √(Tr(ρ)) - ∑xᵢ² , √(Tr(ρ)) - ∑yᵢ² """
    function make_loc_cons_S₂(ρ,d,t,Lx)
        d₁,d₂ = d
        n     = sum(d)
        MB    = mom.get_ℝ_block_diag(d,t.- 1)
        loc_con = Dict() ; g₂ = Dict()
        stρ   = sqrt(sum([ρ[i,i] for i in 1:n]))
        for b in keys(MB)
            xRterm = sum([uc.idx2var(Lx, MB[b] .+ 2* [uc.eᵢ(n,k)])    for k in 1:d₁ ] ) # L((∑xᵢ²)⋅[x,y]ₜ₋₁[x,y]ₜ₋₁ᵀ )
            yRterm = sum([uc.idx2var(Lx, MB[b] .+ 2* [uc.eᵢ(n,k+d₁)]) for k in 1:d₂ ]) # L((∑yᵢ²)⋅[x,y]ₜ₋₁[x,y]ₜ₋₁ᵀ )
            loc_con[b,"x"] = stρ*uc.idx2var(Lx, MB[b]) - xRterm # √(Tr(ρ))⋅L([x,y]ₜ₋₁[x,y]ₜ₋₁ᵀ) - L((∑xᵢ²)⋅[x,y]ₜ₋₁[x,y]ₜ₋₁ᵀ )
            loc_con[b,"y"] = stρ*uc.idx2var(Lx, MB[b]) - yRterm # √(Tr(ρ))⋅L([x,y]ₜ₋₁[x,y]ₜ₋₁ᵀ) - L((∑yᵢ²)⋅[x,y]ₜ₋₁[x,y]ₜ₋₁ᵀ )
            g₂[b] = xRterm - yRterm
        end
        return loc_con, g₂
    end

    """ L((Tr(ρ) - ∑xᵢ²)⋅η) ⪰ 0 ∀ η ∈ "even-degree"-principle-submatrices of [x,y]ₜ₋₁[x,y]ₜ₋₁ᵀ
        L((∑yᵢ² - 1)⋅η) = 0 """
    function make_loc_cons_S₂₁(ρ,d,t,Lx)
        d₁,d₂ = d
        n     = sum(d)
        MB    = mom.get_ℝ_block_diag(d,t.- 1)
        trρ   = sum([ρ[i,i] for i in 1:size(ρ)[1]])
        loc_con = Dict()
        g₂ = Dict()
        for b in keys(MB)
            MBLx   = uc.idx2var(Lx, MB[b]) # L([x,y]ₜ₋₁[x,y]ₜ₋₁ᵀ)
            xRterm = sum([uc.idx2var(Lx, MB[b] .+ [2*uc.eᵢ(n,k)]) for k in 1:d₁ ]) # L((∑xᵢ²)⋅[x,y]ₜ₋₁[x,y]ₜ₋₁ᵀ )
            yRterm = sum([uc.idx2var(Lx, MB[b] .+ [2*uc.eᵢ(n,k+d₁)]) for k in 1:d₂ ]) # L((∑yᵢ²)⋅[x,y]ₜ₋₁[x,y]ₜ₋₁ᵀ )
            loc_con[b] = MBLx - xRterm #  tr(ρ)L([x,y]ₜ₋₁[x,y]ₜ₋₁ᵀ) - L((∑xᵢ²)⋅[x,y]ₜ₋₁[x,y]ₜ₋₁ᵀ ) ⪰ 0
            g₂[b]      = MBLx - yRterm # tr(ρ)L([x,y]ₜ₋₁[x,y]ₜ₋₁ᵀ) - L((∑yᵢ²)⋅[x,y]ₜ₋₁[x,y]ₜ₋₁ᵀ ) = 0
        end
        return loc_con, g₂
    end

    """ ρ⊗L(η) - L( (xxᵀ⊗yyᵀ) ⊗ η ) ⪰ 0 ∀ η ∈ even-degree-principle-submatrices of ([x,y]₌ₜ₋₂[x,y]₌ₜ₋₂ᵀ) """
    function make_G_con(ρ,d,t,Lx)
        MB          = mom.get_ℝ_block_diag(d,t.- 2)  # exponents of [x, y]ₜ₋₂[x, y]ᵀₜ₋₂
        xxᵀtensyyᵀ  = mom.make_xxᵀ_tens_yyᵀ(d)       # exponents of xxᵀ⊗yyᵀ

        G_con = Dict()
        for b in keys(MB)
            Lterm       = kron(ρ, uc.idx2var(Lx,MB[b]))  # ρ⊗L([x, y]ₜ₋₂[x, y]ᵀₜ₋₂)
            Rterm       = uc.idx2var(Lx, mom.var_kron(xxᵀtensyyᵀ,MB[b]))  # L( (xxᵀ⊗yyᵀ) ⊗ ([x, y]ₜ₋₂[x, y]ᵀₜ₋₂) )

            G_con[b]    = Lterm - Rterm # ρ⊗L([x,y]ₜ₋₂[x,y]ᵀₜ₋₂) - L( (xxᵀ⊗yyᵀ) ⊗ ([x,y]ₜ₋₂[x,y]ᵀₜ₋₂) )
        end
        return G_con
    end
end
