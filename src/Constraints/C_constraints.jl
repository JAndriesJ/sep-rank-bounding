module C_constraints
    include(dirname(@__FILE__)*"\\Utils_cons.jl")
    include(dirname(dirname(@__FILE__))*"\\Moments\\Moments.jl")
    using .Utils_cons
    using .Moments
    const uc = Utils_cons
    const mom = Moments

    export make_mon_expo_keys,
           make_PSD_con,
           zeroprop,
           make_ord4_con,
           make_loc_cons_S_inf,
           make_loc_cons_S₂,
           make_loc_cons_S₂₁,
           make_Gᴿ_con

    make_mon_expo_keys(d,t::Int) = Moments.make_mon_expo(d,t*2)

    """L([xᵣₑ,xᵢₘ,yᵣₑ,yᵢₘ]ₜ[xᵣₑ,xᵢₘ,yᵣₑ,yᵢₘ]ₜᵀ) ⪰ 0"""
    function make_PSD_con(d,t,Lx)
        MB = mom.get_ℂ_block_diag(d,t)
        PSD_con = Dict(zip([keys(MB)...], map(b -> uc.idx2var(Lx, MB[b]),[keys(MB)...])))
        return PSD_con
    end

    function zeroprop(ρ,d,t,Lx)
        nar((i,j)) = [(i,0,j,0) (0,i,j,0) (i,0,0,j) (0,i,0,j)]
        nair((i,j,k,l)) = vcat(mom.eᵢ(d[1],i),mom.eᵢ(d[1],j),mom.eᵢ(d[2],k),mom.eᵢ(d[2],l))
        zdiags =[(i,j) for i in 1:d[1], j in 1:d[2] if ρ[i+(j-1)*d[1],i+(j-1)*d[1]] == 0 ]
        return uc.idx2var(Lx, [k .* m  for m ∈  map(nair, hcat(nar.(zdiags)...)) for k in vcat(1,3:t[1])])
    end

    """L(xx*⊗yy*) ⪰ 0 ⟺ """
    make_ord4_con(d,Lx) = uc.idx2varxx̄ᵀtyȳᵀ(Lx,mom.make_xx̄ᵀ_tens_yȳᵀ(d))

    """ Lᴿ(gᴿ⋅[xᵣₑ,xᵢₘ,yᵣₑ,yᵢₘ]ₜ₋₁[xᵣₑ,xᵢₘ,yᵣₑ,yᵢₘ]ₜ₋₁ᵀ) ⪰ 0
        gᴿ = √(ρₘₐₓ) - ((xᵣₑ)ᵢ² + (xᵢₘ)ᵢ²) , i ∈ [d_1],
        or   √(ρₘₐₓ) - ((yᵣₑ)ⱼ² + (yᵢₘ)ⱼ²) , j ∈ [d_2] """
    function make_loc_con(Lx,sqρ,eᵤ,eᵥ,B)
        sz  = size(B)
        Lxᵢ = uc.idx2var(Lx, B .+ [2*eᵤ]) #  L((xᵣₑ)ᵢ²⋅ηₜ₋₁)
        Lxⱼ = uc.idx2var(Lx, B .+ [2*eᵥ]) #  L((xᵢₘ)ᵢ²⋅ηₜ₋₁)
        return  sqρ*uc.idx2var(Lx, B) - (Lxᵢ + Lxⱼ)  # L( (√ρₘₐₓ - (xᵣₑ)ᵢ² -(xᵢₘ)ᵢ²)⋅ηₜ₋₁)
    end

    function make_loc_cons_S_inf(ρ,d,t,Lx)
        d₁,d₂   = d
        n       = sum(2 .*d)
        MB      = mom.get_ℂ_block_diag(d,t .- 1)
        sqρ     = sqrt(maximum(norm.(ρ)))

        loc_con = Dict()
        for b in keys(MB)
            for k in 1:d₁ # Constraint:  Lᴿ( (√ρₘₐₓ - (xᵣₑ)ᵢ² - (xᵢₘ)ᵢ²)⋅ηₜ₋₁)) ⪰ 0 for k ∈ [d₁]
                loc_con[(b,"x²_$k")] = make_loc_con(Lx,sqρ,uc.eᵢ(n,k), uc.eᵢ(n,k+d₁), MB[b])
            end
            for k in 1:d₂ # Constraint:  Lᴿ( (√ρₘₐₓ - (yᵣₑ)ᵢ² - (yᵢₘ)ᵢ²)⋅ηₜ₋₁) ⪰ 0 for k ∈ [d₂]
                loc_con[(b,"y²_$k")] = make_loc_con(Lx,sqρ,uc.eᵢ(n,k+2*d₁),uc.eᵢ(n,k+2*d₁+d₂), MB[b])
            end
        end
        return loc_con
    end

    """Lᴿ(gᴿ⋅[xᵣₑ,xᵢₘ,yᵣₑ,yᵢₘ]ₜ₋₁[xᵣₑ,xᵢₘ,yᵣₑ,yᵢₘ]ₜ₋₁ᵀ) ⪰ 0
        gᴿ   = √Tr(ρ) - ∑ᵈᵢ((xᵣₑ)ᵢ² + (xᵢₘ)ᵢ²),
                √Tr(ρ) - ∑ᵈᵢ((yᵣₑ)ᵢ² + (yᵢₘ)ᵢ²)"""
    tmp(Lx,B,n,s,l) = sum([uc.idx2var(Lx, B .+ [2* uc.eᵢ(n,k)]) for k in s:(s+l)])
    function make_loc_cons_S₂(ρ,d,t,Lx)
        d₁,d₂       = d
        n           = sum(2 .*d)
        MB          = mom.get_ℂ_block_diag(d, t.- 1)
        sqrt_tr_ρ   = sqrt(real(tr(ρ)))

        loc_con = Dict()
        for b in keys(MB)
            xRterm = tmp(Lx,MB[b],n,1,d₁)    + tmp(Lx,MB[b],n,d₁,d₁)        #  L(∑ᵈᵢ((xᵣₑ)ᵢ² + (xᵢₘ)ᵢ² ) ⋅ ηₜ₋₁ )
            yRterm = tmp(Lx,MB[b],n,2*d₁,d₂) + tmp(Lx,MB[b],n,2*d₁+d₂,d₂)   #  L(∑ᵈᵢ((yᵣₑ)ᵢ² + (yᵢₘ)ᵢ² ) ⋅ ηₜ₋₁ )
            loc_con[b,"x"] = sqrt_tr_ρ*uc.idx2var(Lx, MB[b]) - xRterm #  √Tr(ρ) ⋅ L(ηₜ₋₁) - L(∑ᵈᵢ((xᵣₑ)ᵢ² + (xᵢₘ)ᵢ² ) ⋅ ηₜ₋₁ )
            loc_con[b,"y"] = sqrt_tr_ρ*uc.idx2var(Lx, MB[b]) - yRterm #  √Tr(ρ) ⋅ L(ηₜ₋₁) - L(∑ᵈᵢ((yᵣₑ)ᵢ² + (yᵢₘ)ᵢ² ) ⋅ ηₜ₋₁ )
        end
        return loc_con
    end

    """Lᴿ(gᴿ⋅[xᵣₑ,xᵢₘ,yᵣₑ,yᵢₘ]ₜ₋₁[xᵣₑ,xᵢₘ,yᵣₑ,yᵢₘ]ₜ₋₁ᵀ) ⪰ 0
        g^ℝ = Tr(ρ) - ∑ᵈᵢ(uₓᵢ² + vₓᵢ²),
                ∑ᵈᵢ(u_yⱼ² + v_yⱼ²)  - 1"""
    function make_loc_cons_S₂₁(ρ,d,t,Lx)
        d₁,d₂ = d
        n     = sum(2 .*d)
        MB     = mom.get_ℂ_block_diag(d,t .- 1)
        tr_ρ  = real(tr(ρ))

        loc_con    = Dict()
        loc_con_eq = Dict()
        for b in keys(MB)
            xRterm = tmp(Lx,MB[b],n,1,d₁)    + tmp(Lx,MB[b],n,d₁,d₁)        #  L(∑ᵈᵢ((xᵣₑ)ᵢ² + (xᵢₘ)ᵢ² ) ⋅ ηₜ₋₁ )
            yRterm = tmp(Lx,MB[b],n,2*d₁,d₂) + tmp(Lx,MB[b],n,2*d₁+d₂,d₂)   #  L(∑ᵈᵢ((yᵣₑ)ᵢ² + (yᵢₘ)ᵢ² ) ⋅ ηₜ₋₁ )
            loc_con[b]    = tr_ρ*uc.idx2var(Lx, MB[b]) - xRterm       #  Tr(ρ)⋅L(ηₜ₋₁) - L(∑ᵈᵢ((xᵣₑ)ᵢ² + (xᵢₘ)ᵢ² ) ⋅ ηₜ₋₁ ) ⪰ 0
            loc_con_eq[b] = uc.idx2var(Lx, MB[b]) - yRterm            #  1⋅L(ηₜ₋₁) -  L(∑ᵈᵢ((yᵣₑ)ᵢ² + (yᵢₘ)ᵢ² ) ⋅ ηₜ₋₁ ) = 0
        end
        return loc_con, loc_con_eq
    end

    """ρ⊗L(η) - L( (xx*⊗yy*) ⊗ η )⪰ 0 ∀ η ∈ blocks of ([x,̄x,y,̄y]ₜ₋₂[x,̄x,y,̄y]*ₜ₋₂)
     L( ρℝ ⊗ η)  - L( Gℝ[k] ⊗ η)⪰ 0 ∀ η ∈ blocks of [uₓ,vₓ,u_y,v_y]ₜ₋₂[uₓ,vₓ,u_y,v_y]^Tₜ₋₂  """

    using LinearAlgebra
    const la = LinearAlgebra

    function make_Gᴿ_con(ρ,d,t,Lx)
        n = sum(2 .*d)
        D = prod(d)
        Gᴿ,sm        = uc.get_Gᴿ(mom.make_xx̄ᵀ_tens_yȳᵀ(d)) #
        MB           = mom.get_ℂ_block_diag(d,t .- 2)
        ρᴿ           = [real(ρ) -imag(ρ); imag(ρ) real(ρ)]

        LGℝη         = Dict()
        tmp2(B,k) = la.kron(sm[k],ones(D .*size(B))) .* uc.idx2var(Lx, mom.var_kron(Gᴿ[k],B)) # L(Gᴿ[k]⋅ηₜ₋₂)
        for b in keys(MB)
            TEMP1 = tmp2(MB[b],1) # L(Gᴿ[1]⋅ηₜ₋₂)
            for k in 2:8
                TEMP1  += tmp2(MB[b],k) # L(Gᴿ[k]⋅ηₜ₋₂)
            end
            TEMP2 = la.kron(ρᴿ,uc.idx2var(Lx,MB[b]))            #
            LGℝη[b]   = TEMP2 - TEMP1 # L(ρᴿ ⊗ ηₜ₋₂)  - ∑ₖL(Gᴿ[k] ⊗ ηₜ₋₂)
        end
        return  LGℝη
    end
end
