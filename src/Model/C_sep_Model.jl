module C_sep_Model
    using JuMP
    # using LinearAlgebra

    srcDir = dirname(dirname(@__FILE__))
    include(srcDir*"\\Constraints\\C_constraints.jl")
    using .C_constraints
    const ccon = C_constraints
    export Modelξₜˢᵉᵖ

    """The model"""
    function Modelξₜˢᵉᵖ(ρ,d,t;con_list ="S∞ sG",noBlock = false)
        model = JuMP.Model()
        @variable(model, Lx[ccon.make_mon_expo_keys(d,t[1])] )## Create variables
        function set_con(c)
            for k in keys(c)
                 size(c[k]) == (1, 1) ?
                 @constraint(model, c[k] .>= 0) :
                 @constraint(model, c[k] in PSDCone())
            end
        end
##      PSD Moment matrix blocks
        PSD_con = ccon.make_PSD_con(d,t,Lx;noBlock=noBlock)
        set_con(PSD_con)
## Fourth order Moment constraints: L(xxᵀ ⊗ yyᵀ) = ρ,
        L_xx̄ᵀ_tens_yȳᵀ = ccon.make_ord4_con(d,Lx)
        @constraint(model, L_xx̄ᵀ_tens_yȳᵀ["real"] .== real(ρ))
        @constraint(model, L_xx̄ᵀ_tens_yȳᵀ["imag"] .== imag(ρ))
## Localizing g constraint: L ≥ 0 on M₂ₜ(S)
        if     occursin("S∞",con_list)
            set_con(ccon.make_loc_cons_S_inf(ρ,d,t,Lx;noBlock=noBlock))
        elseif occursin("S₂",con_list)
            loc_con, loc_con_eq  = ccon.make_loc_cons_S₂(ρ,d,t,Lx;noBlock=noBlock)
            for k in keys(loc_con_eq)
                @constraint(model, loc_con_eq[k] .==  zeros(size(loc_con_eq[k])))
            end
            set_con(loc_con)
        elseif occursin("S₂₁",con_list)
            loc_con, loc_con_eq  = ccon.make_loc_cons_S₂₁(ρ,d,t,Lx)
            for k in keys(loc_con_eq)
                @constraint(model, loc_con_eq[k] .==  zeros(size(loc_con_eq[k])))
            end
            set_con(loc_con)
        end
## G Constraints
        if  occursin("sG",con_list)
            #println("----------------G-constraints are active")
            set_con(ccon.make_Gᴿ_con(ρ,d,t,Lx;noBlock=noBlock))
        end
        #  Set objective
        @objective(model, Min, Lx[zeros(sum(2 .* d))])
        return model
    end
end
