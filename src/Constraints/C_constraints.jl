module C_constraints
    include(dirname(@__FILE__)*"\\Utils_cons.jl")
    include(dirname(dirname(@__FILE__))*"\\Moments\\Moments.jl")
    using .Utils_cons
    using .Moments
    using LinearAlgebra
    const la = LinearAlgebra
    const uc = Utils_cons
    const mom = Moments

    export make_mon_expo_keys,
           make_PSD_con,
           make_ord4_con,
           make_loc_cons_S_inf,
           make_loc_cons_S₂,
           make_loc_cons_S₂₁,
           make_Gᴿ_con

    make_mon_expo_keys(d,t::Int) = Moments.make_mon_expo(d,t*2)

    """L([xᵣₑ,xᵢₘ,yᵣₑ,yᵢₘ]ₜ[xᵣₑ,xᵢₘ,yᵣₑ,yᵢₘ]ₜᵀ) ⪰ 0"""
    function make_PSD_con(d,t,Lx;noBlock=false)
        MMexᴿ,MMCoefᴿ = Moments.get_ℂ_block_diag(d,t;noBlock=noBlock)
        println(keys(MMexᴿ))
        return  Utils_cons.idx2var_dic(Lx,MMexᴿ,MMCoefᴿ)
    end


    """L(xx*⊗yy*) ⪰ 0 ⟺ """
    make_ord4_con(d,Lx) = uc.idx2varxx̄ᵀtyȳᵀ(Lx,mom.make_xx̄ᵀ_tens_yȳᵀ(d))

    """ Lᴿ(gᴿ⋅[xᵣₑ,xᵢₘ,yᵣₑ,yᵢₘ]ₜ₋₁[xᵣₑ,xᵢₘ,yᵣₑ,yᵢₘ]ₜ₋₁ᵀ) ⪰ 0
        gᴿ = √(ρₘₐₓ) - ((xᵣₑ)ᵢ² + (xᵢₘ)ᵢ²) , i ∈ [d_1],
        or   √(ρₘₐₓ) - ((yᵣₑ)ⱼ² + (yᵢₘ)ⱼ²) , j ∈ [d_2] """
    function make_loc_con(Lx,sqρ,eᵤ,eᵥ,B,C)
        Lxᵢxⱼ = uc.idx2var_arr(Lx,B,C,[2*eᵤ]) + uc.idx2var_arr(Lx,B,C,[2*eᵥ]) #  L((xᵣₑ)ᵢ²+ (xᵢₘ)ᵢ²⋅ηₜ₋₁)
        return  sqρ*uc.idx2var_arr(Lx,B,C) - Lxᵢxⱼ  # L((√ρₘₐₓ - ((xᵣₑ)ᵢ²+(xᵢₘ)ᵢ²))⋅ηₜ₋₁)
    end
    function make_loc_cons_S_inf(ρ,d,t,Lx;noBlock=false)
        d₁,d₂ = d ; n = sum(2 .*d) ; sqρ   = sqrt(maximum(norm.(ρ)))
        MMexᴿ,MMCoefᴿ = mom.get_ℂ_block_diag(d,t .- 1;noBlock=noBlock)
        loc_con = Dict()
        for b in keys(MMCoefᴿ)
            for k in 1:d₁ # Constraint:  Lᴿ( (√ρₘₐₓ-((xᵣₑ)ᵢ²-(xᵢₘ)ᵢ²))⋅ηₜ₋₁)) ⪰ 0 for k ∈ [d₁]
                loc_con[(b,"x²_$k")] = make_loc_con(Lx,sqρ,uc.eᵢ(n,k),uc.eᵢ(n,k+d₁),MMexᴿ[b],MMCoefᴿ[b])
            end
            for k in 1:d₂ # Constraint:  Lᴿ( (√ρₘₐₓ-((yᵣₑ)ᵢ²+(yᵢₘ)ᵢ²))⋅ηₜ₋₁) ⪰ 0 for k ∈ [d₂]
                loc_con[(b,"y²_$k")] = make_loc_con(Lx,sqρ,uc.eᵢ(n,k+2*d₁),uc.eᵢ(n,k+2*d₁+d₂),MMexᴿ[b],MMCoefᴿ[b])
            end
        end
        return loc_con
    end

    """Lᴿ(gᴿ⋅[xᵣₑ,xᵢₘ,yᵣₑ,yᵢₘ]ₜ₋₁[xᵣₑ,xᵢₘ,yᵣₑ,yᵢₘ]ₜ₋₁ᵀ) ⪰ 0
        gᴿ   = √Tr(ρ) - ∑ᵈᵢ((xᵣₑ)ᵢ² + (xᵢₘ)ᵢ²),
                √Tr(ρ) - ∑ᵈᵢ((yᵣₑ)ᵢ² + (yᵢₘ)ᵢ²)"""
    tmp(Lx,B,C,n,s,l) = sum([uc.idx2var_arr(Lx,B,C,[2* uc.eᵢ(n,k)] ) for k in s:(s+l)]) #### SOmething is wrong
    function make_loc_cons_S₂(ρ,d,t,Lx;noBlock=false)
        d₁,d₂ = d ; n = sum(2 .*d) ; sqrt_tr_ρ = sqrt(real(tr(ρ)))
        MMexᴿ,MMCoefᴿ = mom.get_ℂ_block_diag(d, t.- 1,noBlock=noBlock)

        loc_con = Dict(); g₂ = Dict()
        for b in keys(MMCoefᴿ)
            xRterm = tmp(Lx,MMexᴿ[b],MMCoefᴿ[b],n,1,d₁)    + tmp(Lx,MMexᴿ[b],MMCoefᴿ[b],n,d₁,d₁)       #  Lᴿ(∑ᵈᵢ((xᵣₑ)ᵢ² + (xᵢₘ)ᵢ² ) ⋅ ηₜ₋₁ )
            yRterm = tmp(Lx,MMexᴿ[b],MMCoefᴿ[b],n,2*d₁,d₂) + tmp(Lx,MMexᴿ[b],MMCoefᴿ[b],n,2*d₁+d₂,d₂)  #  Lᴿ(∑ᵈᵢ((yᵣₑ)ᵢ² + (yᵢₘ)ᵢ² ) ⋅ ηₜ₋₁ )
            loc_con[b,"x"] = sqrt_tr_ρ*uc.idx2var_arr(Lx,MMexᴿ[b],MMCoefᴿ[b]) - xRterm #  √Tr(ρ) ⋅ Lᴿ(ηₜ₋₁) - Lᴿ(∑ᵈᵢ((xᵣₑ)ᵢ² + (xᵢₘ)ᵢ² ) ⋅ ηₜ₋₁ )
            # loc_con[b,"y"] = sqrt_tr_ρ*uc.idx2var_arr(Lx,MMexᴿ[b],MMCoefᴿ[b]) - yRterm #  √Tr(ρ) ⋅ Lᴿ(ηₜ₋₁) - Lᴿ(∑ᵈᵢ((yᵣₑ)ᵢ² + (yᵢₘ)ᵢ² ) ⋅ ηₜ₋₁ )
            g₂[b] = xRterm - yRterm
        end
        return loc_con, g₂
    end

    """Lᴿ(gᴿ⋅[xᵣₑ,xᵢₘ,yᵣₑ,yᵢₘ]ₜ₋₁[xᵣₑ,xᵢₘ,yᵣₑ,yᵢₘ]ₜ₋₁ᵀ) ⪰ 0
        g^ℝ = Tr(ρ) - ∑ᵈᵢ(uₓᵢ² + vₓᵢ²),
                ∑ᵈᵢ(u_yⱼ² + v_yⱼ²)  - 1"""
    function make_loc_cons_S₂₁(ρ,d,t,Lx;noBlock=false)
        d₁,d₂ = d
        n     = sum(2 .*d)
        # MB     = mom.get_ℂ_block_diag(d,t .- 1)
        MMexᴿ,MMCoefᴿ = mom.get_ℂ_block_diag(d, t.- 1,noBlock=noBlock)
        tr_ρ  = real(tr(ρ))

        loc_con    = Dict()
        loc_con_eq = Dict()
        for b in keys(MMCoefᴿ)
            xRterm = tmp(Lx,MMexᴿ[b],MMCoefᴿ[b],n,1,d₁)    + tmp(Lx,MMexᴿ[b],MMCoefᴿ[b],n,d₁,d₁)        #  L(∑ᵈᵢ((xᵣₑ)ᵢ² + (xᵢₘ)ᵢ² ) ⋅ ηₜ₋₁ )
            yRterm = tmp(Lx,MMexᴿ[b],MMCoefᴿ[b],n,2*d₁,d₂) + tmp(Lx,MMexᴿ[b],MMCoefᴿ[b],n,2*d₁+d₂,d₂)   #  L(∑ᵈᵢ((yᵣₑ)ᵢ² + (yᵢₘ)ᵢ² ) ⋅ ηₜ₋₁ )
            loc_con[b]    = tr_ρ*uc.idx2var_arr(Lx,MMexᴿ[b],MMCoefᴿ[b]) - xRterm       #  Tr(ρ)⋅L(ηₜ₋₁) - L(∑ᵈᵢ((xᵣₑ)ᵢ² + (xᵢₘ)ᵢ² ) ⋅ ηₜ₋₁ ) ⪰ 0
            loc_con_eq[b] = uc.idx2var_arr(Lx, MMexᴿ[b],MMCoefᴿ[b]) - yRterm            #  1⋅L(ηₜ₋₁) -  L(∑ᵈᵢ((yᵣₑ)ᵢ² + (yᵢₘ)ᵢ² ) ⋅ ηₜ₋₁ ) = 0
        end
        return loc_con, loc_con_eq
    end

    """ρ⊗L(η) - L( (xx*⊗yy*) ⊗ η )⪰ 0 ∀ η ∈ blocks of ([x,̄x,y,̄y]ₜ₋₂[x,̄x,y,̄y]*ₜ₋₂)
     L( ρℝ ⊗ η)  - L( Gℝ[k] ⊗ η)⪰ 0 ∀ η ∈ blocks of [uₓ,vₓ,u_y,v_y]ₜ₋₂[uₓ,vₓ,u_y,v_y]^Tₜ₋₂  """



    function make_Gᴿ_con(ρ,d,t,Lx;noBlock=false)
        n = sum(2 .*d) ; D = prod(d) ; ρᴿ = [real(ρ) -imag(ρ); imag(ρ) real(ρ)]
        Gᴿ,sm        = uc.get_Gᴿ(mom.make_xx̄ᵀ_tens_yȳᵀ(d)) #
        MMexᴿ,MMCoefᴿ = mom.get_ℂ_block_diag(d, t.- 2,noBlock=noBlock)

        LGℝη         = Dict()
        tmp2(B,C,k) = la.kron(sm[k],ones(D .*size(B))) .*
            uc.idx2var_arr(
                Lx,
                mom.var_kron_C(Gᴿ[k],B),
                mom.var_kron_C(fill(0.0,size(Gᴿ[k])...),C)) # L(Gᴿ[k]⋅ηₜ₋₂)
        for b in keys(MMCoefᴿ)
            TEMP1 = sum([tmp2(MMexᴿ[b],MMCoefᴿ[b],k) for k in 1:8]) # ∑ₖL(Gᴿ[k]⋅ηₜ₋₂)
            TEMP2 = la.kron(ρᴿ,uc.idx2var_arr(Lx,MMexᴿ[b],MMCoefᴿ[b]))            #
            LGℝη[b]   = TEMP2 - TEMP1 # L(ρᴿ ⊗ ηₜ₋₂)  - ∑ₖL(Gᴿ[k] ⊗ ηₜ₋₂)
        end
        return  LGℝη
    end
end
