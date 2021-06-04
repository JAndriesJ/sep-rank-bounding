module C_constraints

    using LinearAlgebra

    include(dirname(@__FILE__)*"\\Utils_cons.jl")
    include(dirname(dirname(@__FILE__))*"\\Moments\\Moments.jl")
    using .Utils_cons
    using .Moments

    const uc = Utils_cons
    const mom = Moments

    export make_PSD_MM_con,
           make_ord4_con,
           make_loc_cons_S₁,
           make_loc_cons_S₂,
           make_loc_cons_S₃,
           make_weakG_con,
           make_Gℝ_con


    get_MB_c(d,t,isLeq=true) = mom.get_ℂ_block_diag(mom.make_mon_expo_mat(sum(2 .*d),t,isLeq),d)
    utl_c(d,t) = sum(2 .*d),get_MB_c(d,t)
    make_mom_expo_keys(n,t) = mom.make_mom_expo_keys(n,t)

    """PSD Moment matrix blocks"""
    function make_PSD_MM_con(d,t,Lx)
        MB = get_MB_c(d,t)
        PSD_MM_con = Dict()
        for block in keys(MB)
            if isempty(MB[block])
                continue
            end
            PSD_MM_con[block] = uc.index_to_var(Lx, MB[block])
        end
        return PSD_MM_con
    end

    """L(xx*⊗yy*)"""
    function make_ord4_con(d,Lx)
        xx̄ᵀ_tens_yȳᵀ    = mom.make_xx̄ᵀ_tens_yȳᵀ(d)
        L_xx̄ᵀ_tens_yȳᵀ  = uc.get_Lxx̄ᵀ_tens_yȳᵀ(xx̄ᵀ_tens_yȳᵀ,Lx)
        return L_xx̄ᵀ_tens_yȳᵀ
    end

    """ L^ℝ(g^ℝ⋅[uₓ,vₓ,uᵥ,vᵥ]ₜ₋₁[uₓ,vₓ,uᵥ,vᵥ]ₜ₋₁ᵀ) ⪰ 0
        g^ℝ = √(ρₘₐₓ) - (uₓᵢ² + vₓᵢ²) , i ∈ [d_1],
                √(ρₘₐₓ) - (u_yⱼ² + v_yⱼ²) , j ∈ [d_2] """
    function make_loc_cons_S₁(ρ,d,t,Lx)
        n,MB    = utl_c(d,t .- 1)
        sqr_ρ   = sqrt(maximum(norm.(ρ)))
        function make_loc_con(eu,ev,key)
            sz = size(MB[key])
            Lxᵢ = uc.index_to_var(Lx, MB[key] + repeat([2*eu], sz...)) #  L(uₓᵢ²⋅ηₜ₋₁)
            Lxⱼ = uc.index_to_var(Lx, MB[key] + repeat([2*ev], sz...)) #  L(vₓᵢ²⋅ηₜ₋₁)
            Lρ = sqr_ρ*uc.index_to_var(Lx, MB[key])  # L(√ρₘₐₓ ⋅ ηₜ₋₁)
            return  Lρ - (Lxᵢ + Lxⱼ)  # L( (√ρₘₐₓ - (uₓᵢ² + vₓᵢ²))⋅ηₜ₋₁)
        end
        loc_con = Dict()
        for b in keys(MB)
            for k in 1:d[1]
                eₖ   = uc.eᵢ(n,k)
                eₖ₊d = uc.eᵢ(n,k+d[1])
                loc_con[(b,k,1)] = make_loc_con(eₖ,eₖ₊d,b) # Constraint: L( (√ρₘₐₓ - (uₓₖ² + vₓₖ²))⋅ηₜ₋₁ ) ⪰ 0 for k ∈ [d₁]
            end
            for k in 1:d[2]
                eₖ₊₂d = uc.eᵢ(n,k+2*d[1])
                eₖ₊₃d = uc.eᵢ(n,k+2*d[1]+d[2])
                loc_con[(b,k,2)] = make_loc_con(eₖ₊₂d,eₖ₊₃d,b) # Constraint: L( (√ρₘₐₓ - (u_yₖ² + v_yₖ²))⋅ηₜ₋₁ᵀ ) ⪰ 0 for k ∈ [d₂]
            end
        end
        return loc_con
    end

    """L^ℝ(g^ℝ⋅[uₓ,vₓ,uᵥ,vᵥ]ₜ₋₁[uₓ,vₓ,uᵥ,vᵥ]ₜ₋₁ᵀ) ⪰ 0
        g^{ℝ} = √(ρₘₐₓ) - ∑ᵈᵢ(uₓᵢ² + vₓᵢ²),
                √(ρₘₐₓ) - ∑ᵈᵢ(u_yⱼ² + v_yⱼ²)"""
    function make_loc_cons_S₂(ρ,d,t,Lx)
        n,MB    = utl_c(d,t.-1)
        loc_con = Dict()
        sqrt_tr_ρ   = sqrt(real(tr(ρ)))

        tmp1(MB,key,n,d) = sum([uc.index_to_var(Lx, MB[key] .+ [2* uc.eᵢ(n,k)]) for k in 1:d])     # ∑ₖ L((uₖ²⋅ηₜ₋₁)
        tmp2(MB,key,n,d) = sum([uc.index_to_var(Lx, MB[key] .+ [2* uc.eᵢ(n,k+d)]) for k in 1:d])   # ∑ₖ L((vₖ²⋅ηₜ₋₁)
        for key in keys(MB)
            xRterm = tmp1(MB,key,n,d[1]) + tmp2(MB,key,n,d[1])                 #  L((||uₓ||² + ||vₓ||²) ⋅ ηₜ₋₁ )
            yRterm = tmp1(MB,key,n,d[2]) + tmp2(MB,key,n,d[2])                 #  L((||u_y||² + ||v_y||²) ⋅ ηₜ₋₁ )
            loc_con[(key,1)] = sqrt_tr_ρ*uc.index_to_var(Lx, MB[key]) - xRterm #  √Tr(ρ) ⋅ L(ηₜ₋₁) - L( (||uₓ||² + ||vₓ||²)⋅ηₜ₋₁ )
            loc_con[(key,2)] = sqrt_tr_ρ*uc.index_to_var(Lx, MB[key]) - yRterm #  √Tr(ρ) ⋅ L(ηₜ₋₁) - L( (||u_y||² + ||v_y||²)⋅ηₜ₋₁ )
        end
        return loc_con
    end

    """L^ℝ(g^ℝ⋅[uₓ,vₓ,uᵥ,vᵥ]ₜ₋₁[uₓ,vₓ,uᵥ,vᵥ]ₜ₋₁ᵀ) ⪰ 0
        for
        g^ℝ = Tr(ρ) - ∑ᵈᵢ(uₓᵢ² + vₓᵢ²),
                ∑ᵈᵢ(u_yⱼ² + v_yⱼ²)  - 1"""
    function make_loc_cons_S₃(ρ,d,t,Lx)
        n,MB    = utl_c(d,t.-1)
        tr_ρ    = real(tr(ρ))

        loc_con = Dict()
        loc_con_eq = Dict()
        tmp1(MB,key,n,d_x) = sum([uc.index_to_var(Lx, MB[key] .+ [2* uc.eᵢ(n,k)]) for k in 1:d_x]) # L((||u_||²⋅ηₜ₋₁)
        tmp2(MB,key,n,d_x) = sum([uc.index_to_var(Lx, MB[key] .+ [2* uc.eᵢ(n,k+d_x)]) for k in 1:d_x]) # L((||v_||²⋅ηₜ₋₁)
        for key in keys(MB)
            xRterm = tmp1(MB,key,n,d[1]) + tmp2(MB,key,n,d[1]) # L((||uₓ||² + ||vₓ||²)⋅ηₜ₋₁ )
            yRterm = tmp1(MB,key,n,d[2]) + tmp2(MB,key,n,d[2]) # L((||uᵥ||² + ||vᵥ||²)⋅ηₜ₋₁ )
            loc_con[key]    = tr_ρ*uc.index_to_var(Lx, MB[key]) - xRterm #  Tr(ρ)⋅L(ηₜ₋₁) - L((||uₓ||² + ||vₓ||²)⋅ηₜ₋₁ ) ≧ 0
            loc_con_eq[key] = uc.index_to_var(Lx, MB[key]) - yRterm #     1⋅L(ηₜ₋₁) - L((||uᵥ||² + ||vᵥ||²)⋅ηₜ₋₁ ) = 0
        end
        return loc_con, loc_con_eq
    end

    """ρ⊗L(η) - L( (xx*⊗yy*) ⊗ η )⪰ 0 ∀ η ∈ blocks of ([x,̄x,y,̄y]ₜ₋₂[x,̄x,y,̄y]*ₜ₋₂)
     L( ρℝ ⊗ η)  - L( Gℝ[k] ⊗ η)⪰ 0 ∀ η ∈ blocks of [uₓ,vₓ,u_y,v_y]ₜ₋₂[uₓ,vₓ,u_y,v_y]^Tₜ₋₂  """
    function make_Gℝ_con(ρ,d,t,Lx)
        xx̄ᵀ_tens_yȳᵀ = Moments.make_xx̄ᵀ_tens_yȳᵀ(d)
        Gℝ,sgn_mat   = uc.get_Gℝ(xx̄ᵀ_tens_yȳᵀ)
        n,MB         = utl_c(d,(t .- 2)) # exponents of ηₜ₋₂ of Mₜ₋₂(L)
        ρℝ           = [real(ρ) -imag(ρ); imag(ρ) real(ρ)]
        LGℝη         = Dict()
        for B in keys(MB)
            LGℝηTEMP = uc.index_to_var(Lx, Moments.var_kron(Gℝ[1],MB[B]))  # L( Gℝ[k] ⊗ ηₜ₋₂) with wrong signs
            LGℝη[B]      =   kron(sgn_mat[1],ones(d[1]*d[2] .*size(MB[B]))) .*LGℝηTEMP           #  L( Gℝ[k] ⊗ ηₜ₋₂)
            for k in 2:8
                LGℝηTEMP = uc.index_to_var(Lx, Moments.var_kron(Gℝ[k],MB[B]))  # L( Gℝ[k]⋅ ηₜ₋₂) with wrong signs
                LGℝη[B]    = LGℝη[B] + kron(sgn_mat[k],ones(d[1]*d[2] .*size(MB[B]))) .* LGℝηTEMP  # L( Gℝ[k] ⊗ ηₜ₋₂)
            end
            LρℝηTEMP = kron(ρℝ,uc.index_to_var(Lx,MB[B]))            # L( ρℝ ⊗ ηₜ₋₂)
            LGℝη[B] = LρℝηTEMP - LGℝη[B] # L( ρℝ ⊗ ηₜ₋₂)  - ∑ₖL( Gℝ[k] ⊗ ηₜ₋₂)
        end
        return  LGℝη
    end



################################################################################

    """
        XXXXXXXXXXXXXXXXXXXXXXXx SKIPPED FOR NOW XXXXXXXXXXXXXXXXXXXXX
        Input: ρ(data matrix),t(Integer),Lx(JuMP variable)
        Output:
        L((ρ - xxᵀ⊗yyᵀ) ⊗ η)) ⪰ 0
        for l ∈ 1,...,t-2.
        η ∈ even-degree-principle-submatrices of ([x,y]₌ₗ[x,y]₌ₗᵀ)
        ⟺
        ρ⊗L(η) - L((xxᵀ⊗yyᵀ) ⊗ η) ⪰ 0
        η ∈ even-degree-principle-submatrices of ([x,y]₌ₗ[x,y]₌ₗᵀ)
    """
    function make_weakG_con(ρ,d,t,Lx)
        n,MB =  utl_c(d,t)
        xx̄ᵀ_tens_yȳᵀ      = mom.make_xx̄ᵀ_tens_yȳᵀ(d)      # exponents of α of xx̄ᵀ⊗yȳᵀ

        weakG_con = Dict()
        for ℓ in 1:(t[1] - 2)
        LMBexp_ℓ = get_MB_c(d,(ℓ,ℓ),false)  # exponents of ηₜ₋₂ of Mₜ₋₂(L)
            for key in keys(LMBexp_ℓ )
                MomMat = LMBexp_ℓ[key]
                if isempty(MomMat)
                    continue
                end
                LMBₜ₋₂             = uc.index_to_var(Lx,MomMat)         # L(ηₜ₋₂)
                rLterm             = kron(real(ρ),LMBₜ₋₂)                  # real(ρ)⊗L(ηₜ₋₂)
                iLterm             = kron(imag(ρ),LMBₜ₋₂)                  # imag(ρ)⊗L(ηₜ₋₂)
                rRterm, iRterm     = uc.get_Lxx̄ᵀ_tens_yȳᵀ(xx̄ᵀ_tens_yȳᵀ, Lx, MomMat)  # exponents of 🍌 ⊗ [uₓ,vₓ,uᵥ,vᵥ]ₜ₋₂[uₓ,vₓ,uᵥ,vᵥ]ᵀₜ₋₂

                weakG_con["r",key] = rLterm - rRterm                           # real(ρ)⊗L(ηₜ₋₂) - L( 🍌 ⊗ [uₓ,vₓ,uᵥ,vᵥ]ₜ₋₂[uₓ,vₓ,uᵥ,vᵥ]ᵀₜ₋₂ )
                weakG_con["i",key] = iLterm - iRterm                           # imag(ρ)⊗L(ηₜ₋₂) - L( 🍌 ⊗ [uₓ,vₓ,uᵥ,vᵥ]ₜ₋₂[uₓ,vₓ,uᵥ,vᵥ]ᵀₜ₋₂ )
            end
        end
        return weakG_con
    end




end
