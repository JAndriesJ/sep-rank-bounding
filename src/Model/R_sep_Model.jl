module R_sep_Model
    using JuMP
    # using LinearAlgebra

    srcDir = dirname(dirname(@__FILE__))
    include(srcDir*"\\Constraints\\R_Constraints.jl")

    using .R_constraints
    const rcon = R_constraints

    export R_Modelξₜˢᵉᵖ
    """The model"""
    function Modelξₜˢᵉᵖ(ρ,d,t;con_list="S_infsG")
        model = Model()
        @variable(model, Lx[rcon.make_mon_expo_keys(sum(d),t[1])] ) # Create variables
        function set_con(c)
            for k in keys(c)
                 size(c[k]) == (1, 1) ?
                 @constraint(model, c[k] .>= 0) :
                 @constraint(model, c[k] in PSDCone())
            end
        end
## PSD Moment matrix
        PSD_MM_con = rcon.make_PSD_con(d,t,Lx)
        set_con(PSD_MM_con) # L(ηᵦ) ⪰ 0 ∀ ηᵦ ∈  ⨁ᵦ ηᵦ = L([x,y]ₜ[x,y]ₜᵀ)
## L(xxᵀ ⊗ yyᵀ) = ρ
        @constraint(model, rcon.make_ord4_con(d,Lx) .==  ρ)
## Zero propagation: ρ = aaᵀ⊗bbᵀ , ρᵢⱼᵢⱼ = 0 ⟹ aᵢbⱼ = 0 ⟹ L(xᵢyⱼ) = 0
        # zero_moms = rcon.zeroprop(ρ,d,t,Lx)
        # @constraint(model,zero_moms .== zeros(size(zero_moms)))
## Localizing g constraint: L ≥ 0 on M₂ₜ(S)
        if     occursin("S∞",con_list)
            println("----------S∞ constraints active----------")
            set_con(rcon.make_loc_cons_S_inf(ρ,d,t,Lx))
        elseif occursin("S₂",con_list)
            println("----------S₂ constraints active----------")
            loc_con, g₂  =  rcon.make_loc_cons_S₂(ρ,d,t,Lx)
            for k in keys(g₂)
                @constraint(model, g₂[k] .==  zeros(size(g₂[k])))
            end
            set_con(loc_con)
        elseif occursin("S₂₁",con_list)
            println("----------S₂₁ constraints active----------")
            loc_con, g₂  =  rcon.make_loc_cons_S₂₁(ρ,d,t,Lx)
            for k in keys(g₂)
                @constraint(model, g₂[k] .==  zeros(size(g₂[k])))
            end
            set_con(loc_con)
        end

## G Constraints
        if  occursin("sG",con_list)
            println("----------G-constraints are active----------")
            set_con(rcon.make_G_con(ρ,d,t,Lx)) # L((ρ-xxᵀ⊗yyᵀ)⊗([x,y]ₜ₋₂[x,y]ᵀₜ₋₂))
        end

        @objective(model, Min, Lx[zeros(sum(d))]) #  Set objective
        return model
    end
end
