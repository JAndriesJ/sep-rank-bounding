module R_sep_Model
    using JuMP
    using LinearAlgebra

    srcDir = dirname(dirname(@__FILE__))
    include(srcDir*"\\Constraints\\R_Constraints.jl")

    using .R_constraints

    export R_Modelξₜˢᵉᵖ

    """The model"""
    function Modelξₜˢᵉᵖ(ρ,d,t,con_list)
        n     = sum(d)
        model = Model()
## Create variables
        list_of_keys = R_constraints.make_mom_expo_keys(n,t) # Define variables in the moment matrix.
        @variable(model, Lx[list_of_keys] ) #
## PSD Moment matrix
        PSD_MM_con = R_constraints.make_PSD_MM_con(d,t,Lx)
        for block in keys(PSD_MM_con)
            @constraint(model, Symmetric(PSD_MM_con[block]) in PSDCone())
        end
## Fourth order Moment constraints: L(xxᵀ ⊗ yyᵀ) = ρ,
        L_xxᵀ_tens_yyᵀ  = R_constraints.make_ord4_con(d,Lx)
        @constraint(model, fix_con,L_xxᵀ_tens_yyᵀ .==  ρ)
## Localizing g constraint: L ≥ 0 on M₂ₜ(S)
        loc_con = Dict()
        if     occursin("S₁",con_list)
            loc_con = R_constraints.make_loc_cons_S₁(ρ,d,t,Lx)
        elseif occursin("S₂",con_list)
            loc_con = R_constraints.R_cons.make_loc_cons_S₂(ρ,d,t,Lx)
        elseif occursin("S₃",con_list)
            loc_con, g₂  =  R_constraints.make_loc_cons_S₃(ρ,d,t,Lx)
            for key in keys(g₂)
                @constraint(model, g₂[key] .==  zeros(size(g₂[key])))
            end
        end
        for key in keys(loc_con)
            if isempty(loc_con[key])
                continue
            end
            if size(loc_con[key]) == (1, 1)
                @constraint(model, loc_con[key] .>= 0)
            else
                @constraint(model, Symmetric(loc_con[key]) in PSDCone())
            end
        end
## weak G Constraints
        if occursin("wG",con_list)
            println("----------------Weak G-constraints are active")
            weakG_con = R_constraints.make_weakG_con(ρ,d,t,Lx)
            for key in keys(weakG_con)
                if isempty(weakG_con[key])
                    continue
                end
                if size(weakG_con[key]) == (1, 1)
                    @constraint(model, weakG_con[key] .>= 0)
                else
                    @constraint(model, Symmetric(weakG_con[key]) in PSDCone())
                end
            end
        end
## G Constraints
        if  occursin("sG",con_list)
            println("----------------G-constraints are active")
            G_con = R_constraints.make_G_con(ρ,d,t,Lx)
            for key in keys(G_con)
                if isempty(G_con[key])
                    continue
                end
                @constraint(model, Symmetric(G_con[key]) in PSDCone())
            end
        end
        #  Set objective
        @objective(model, Min, Lx[zeros(n)])
        return model
    end
end
