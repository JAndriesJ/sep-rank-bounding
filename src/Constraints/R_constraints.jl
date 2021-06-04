module R_constraints
    using LinearAlgebra

    include(dirname(@__FILE__)*"\\Utils_cons.jl")
    include(dirname(dirname(@__FILE__))*"\\Moments\\Moments.jl")

    using .Utils_cons
    using .Moments

    export make_mom_expo_keys,
           make_PSD_MM_con,
           make_ord4_con,
           make_loc_cons_S₁,
           make_loc_cons_S₂,
           make_loc_cons_S₃,
           make_weakG_con,
           make_G_con

    get_MB_r(d,t,isLeq=true) = Moments.get_ℝ_block_diag(Moments.make_mon_expo_mat(sum(d),t,isLeq),d)
    utl_r(d,t) = sum(d),get_MB_r(d,t)
    make_mom_expo_keys(n,t) = Moments.make_mom_expo_keys(n,t)

    """PSD Moment matrix blocks"""
    function make_PSD_MM_con(d,t,Lx)
        mom_matₜ_expo_blocks = get_MB_r(d,t)
        PSD_MM_con = Dict()
        for block in keys(mom_matₜ_expo_blocks)
            if isempty(mom_matₜ_expo_blocks[block])
                continue
            end
            PSD_MM_con[block] = Utils_cons.index_to_var(Lx, mom_matₜ_expo_blocks[block])
        end
        return PSD_MM_con
    end

    """L(xxᵀ⊗yyᵀ)"""
    function make_ord4_con(d,Lx)
        xxᵀ_tens_yyᵀ    = Moments.make_xxᵀ_tens_yyᵀ(d)
        L_xxᵀ_tens_yyᵀ  = Utils_cons.index_to_var(Lx,xxᵀ_tens_yyᵀ)
        return L_xxᵀ_tens_yyᵀ
    end

    """ L(g⋅η) ⪰ 0 for η ∈ even-degree-principle-submatrices of [x,y]ₜ₋₁[x,y]ₜ₋₁ᵀ
                        g = √(ρₘₐₓ) - xᵢ² , √(ρₘₐₓ) - yᵢ² and i ∈ [d] """
    function make_loc_cons_S₁(ρ,d,t,Lx)
        n,MB    = utl_r(d,t .- 1)
        sqr_ρ   = sqrt(maximum(ρ))
        function make_loc_con(eᵢ,key)
            sz = size(MB[key])
            Lxᵢ = Utils_cons.index_to_var(Lx, MB[key] + repeat([2*eᵢ], sz...)) #  L(xᵢ²⋅η)
            Lρ = sqr_ρ*Utils_cons.index_to_var(Lx, MB[key])  # L(√ρₘₐₓ⋅η)
            return  Lρ - Lxᵢ
        end
        loc_con = Dict()
        for key in keys(MB)
            for k in 1:d[1]
                eₖ = Utils_cons.eᵢ(n,k)
                loc_con[(key,k,1)] = make_loc_con(eₖ,key) # Constraint:  L ≧ 0 on 𝑀(Sᶜᵖ)   : L((√(ρₘₐₓ) - xᵢ²))⋅η) ⪰ 0
            end
            for k in 1:d[2]
                eₖ₊d₁ = Utils_cons.eᵢ(n,k+d[1])
                loc_con[(key,k,2)] = make_loc_con(eₖ₊d₁,key) # Constraint:  L ≧ 0 on 𝑀(Sᶜᵖ)   : L((√(ρₘₐₓ) - yᵢ²))⋅η) ⪰ 0
            end
        end
        return loc_con
    end

    """L(g⋅η) ⪰ 0 ∀ η ∈ even-degree-principle-submatrices of [x,y]ₜ₋₁[x,y]ₜ₋₁ᵀ
                       g = √(Tr(ρ)) - ∑xᵢ² , √(Tr(ρ)) - ∑yᵢ² """
    function make_loc_cons_S₂(ρ,d,t,Lx)
        n,MB    = utl_r(d,t .- 1)
        loc_con = Dict()
        sqrt_tr_ρ   = sqrt(tr(ρ))
        for key in keys(MB)
            xRterm = sum([Utils_cons.index_to_var(Lx, MB[key] .+ [2* Utils_cons.eᵢ(n,k)]) for k in 1:d[1] ]) # L((∑xᵢ²)⋅[x,y]ₜ₋₁[x,y]ₜ₋₁ᵀ )
            yRterm = sum([Utils_cons.index_to_var(Lx, MB[key] .+ [2* Utils_cons.eᵢ(n,k+d[1])]) for k in 1:d[2] ]) # L((∑yᵢ²)⋅[x,y]ₜ₋₁[x,y]ₜ₋₁ᵀ )
            g₁ = sqrt_tr_ρ*Utils_cons.index_to_var(Lx, MB[key]) - xRterm #  √(ρₘₐₓ)⋅L([x,y]ₜ₋₁[x,y]ₜ₋₁ᵀ) - L((∑xᵢ²)⋅[x,y]ₜ₋₁[x,y]ₜ₋₁ᵀ )
            g₂ = sqrt_tr_ρ*Utils_cons.index_to_var(Lx, MB[key]) - yRterm #  √(ρₘₐₓ)⋅L([x,y]ₜ₋₁[x,y]ₜ₋₁ᵀ) - L((∑yᵢ²)⋅[x,y]ₜ₋₁[x,y]ₜ₋₁ᵀ )
            loc_con[key,1] = g₁
            loc_con[key,2] = g₂
        end
        return loc_con
    end

    """ L((Tr(ρ) - ∑xᵢ²)⋅η) ⪰ 0 ∀ η ∈ "even-degree"-principle-submatrices of [x,y]ₜ₋₁[x,y]ₜ₋₁ᵀ
        L((∑yᵢ² - 1)⋅η) = 0 """
    function make_loc_cons_S₃(ρ,d,t,Lx)
        n,MB    = utl_r(d,t .- 1)
        tr_ρ    = tr(ρ)

        loc_con = Dict()
        g₂ = Dict()
        for key in keys(MB)
            MBLx    = Utils_cons.index_to_var(Lx, MB[key]) # L([x,y]ₜ₋₁[x,y]ₜ₋₁ᵀ)
            xRterm = sum([Utils_cons.index_to_var(Lx, MB[key] .+ [2* Utils_cons.eᵢ(n,k)]) for k in 1:d[1] ]) # L((∑xᵢ²)⋅[x,y]ₜ₋₁[x,y]ₜ₋₁ᵀ )
            loc_con[key,1] =  tr_ρ*MBLx - xRterm #  tr(ρ)⋅L([x,y]ₜ₋₁[x,y]ₜ₋₁ᵀ) - L((∑xᵢ²)⋅[x,y]ₜ₋₁[x,y]ₜ₋₁ᵀ ) ⪰ 0

            yRterm = sum([Utils_cons.index_to_var(Lx, MB[key] .+ [2* Utils_cons.eᵢ(n,k+d[1])]) for k in 1:d[2] ]) # L((∑yᵢ²)⋅[x,y]ₜ₋₁[x,y]ₜ₋₁ᵀ )
            g₂[key] = MBLx - yRterm # L([x,y]ₜ₋₁[x,y]ₜ₋₁ᵀ) - L((∑yᵢ²)⋅[x,y]ₜ₋₁[x,y]ₜ₋₁ᵀ ) = 0
        end
        return loc_con, g₂
    end

    """ρ⊗L(η) - L((xxᵀ⊗yyᵀ) ⊗ η) ⪰ 0 ∀ η ∈ even-degree-principle-submatrices of ([x,y]₌ₗ[x,y]₌ₗᵀ) """
    function make_weakG_con(ρ,d,t,Lx)
        n,MB    =  utl_r(d,t .- 1)

        weakG_con = Dict()
        xxᵀ_tens_yyᵀ = Moments.make_xxᵀ_tens_yyᵀ(d)                # exponents of xxᵀ⊗yyᵀ
        for ℓ in 1:(t[1] - 2)
            LMBexp_ℓ = get_MB_r(d,(ℓ,ℓ),false)  # exponents of [x,y]₌ₗ[x,y]₌ₗᵀ

            for key in keys(LMBexp_ℓ)
                LMBexp_1ℓ         = Moments.var_kron(xxᵀ_tens_yyᵀ,LMBexp_ℓ[key])  #  # exponents of  (xxᵀ⊗yyᵀ)⊗([x,y]₌ₗ[x,y]₌ₗᵀ)
                LMB_ℓ             = Utils_cons.index_to_var(Lx,LMBexp_ℓ[key])  # L([x,y]₌ₗ[x,y]₌ₗᵀ)
                Rterm             = Utils_cons.index_to_var(Lx,LMBexp_1ℓ) # L((xxᵀ⊗yyᵀ)⊗([x,y]₌ₗ[x,y]₌ₗᵀ))

                Lterm             = kron(ρ,LMB_ℓ)                  # ρ⊗L([x,y]₌ₗ[x,y]₌ₗᵀ)
                weakG_con[key,ℓ]  = Lterm - Rterm                  # ρ⊗L([x,y]₌ₗ[x,y]₌ₗᵀ) - L((xxᵀ⊗yyᵀ)⊗([x,y]₌ₗ[x,y]₌ₗᵀ)) ⪰ 0, ℓ ∈ ,1,t-deg(G)/2
            end
        end
        return weakG_con
    end

    """ ρ⊗L(η) - L( (xxᵀ⊗yyᵀ) ⊗ η )⪰ 0 ∀ η ∈ even-degree-principle-submatrices of ([x,y]₌ₜ₋₂[x,y]₌ₜ₋₂ᵀ) """
    function make_G_con(ρ,d,t,Lx)
        n,LMBexpₜ₋₂   =  utl_r(d,t .- 2)  # exponents of [x, y]ₜ₋₂[x, y]ᵀₜ₋₂
        xxᵀ_tens_yyᵀ  = Moments.make_xxᵀ_tens_yyᵀ(d)       # exponents of xxᵀ⊗yyᵀ

        G_con = Dict()
        for key in keys(LMBexpₜ₋₂ )
            LMBₜ₋₂             = Utils_cons.index_to_var(Lx,LMBexpₜ₋₂[key])    # L([x, y]ₜ₋₂[x, y]ᵀₜ₋₂)
            Lterm              = kron(ρ,LMBₜ₋₂)                                # ρ⊗L([x, y]ₜ₋₂[x, y]ᵀₜ₋₂)

            LMBexp_1ₜ₋₂        = Moments.var_kron(xxᵀ_tens_yyᵀ,LMBexpₜ₋₂[key]) # exponents of (xxᵀ⊗yyᵀ) ⊗ ([x, y]ₜ₋₂[x, y]ᵀₜ₋₂)
            Rterm              = Utils_cons.index_to_var(Lx,LMBexp_1ₜ₋₂)       # L( (xxᵀ⊗yyᵀ) ⊗ ([x, y]ₜ₋₂[x, y]ᵀₜ₋₂) )

            G_con[key]         = Lterm - Rterm                                 # ρ⊗L([x, y]ₜ₋₂[x, y]ᵀₜ₋₂) - L( (xxᵀ⊗yyᵀ) ⊗ ([x, y]ₜ₋₂[x, y]ᵀₜ₋₂) )
        end
        return G_con
    end
end
