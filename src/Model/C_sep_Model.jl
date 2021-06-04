module C_sep_Model
    using JuMP
    using LinearAlgebra

    srcDir = dirname(dirname(@__FILE__))
    include(srcDir*"\\Constraints\\C_constraints.jl")
    using .C_constraints

    export Modelξₜˢᵉᵖ


    """The model"""
    function Modelξₜˢᵉᵖ(ρ,d,t,con_list)
        n = sum(2 .* d)
        model = JuMP.Model()
        list_of_keys = C_constraints.make_mom_expo_keys(n,t) # Define variables in the moment matrix.
        @variable(model, Lx[list_of_keys] ) #
##      PSD Moment matrix blocks
        PSD_MM_con = C_constraints.make_PSD_MM_con(d,t,Lx)
        for b in keys(PSD_MM_con)
            @constraint(model, Symmetric(PSD_MM_con[b]) in PSDCone())
        end
## Fourth order Moment constraints: L(xxᵀ ⊗ yyᵀ) = ρ,
        L_xx̄ᵀ_tens_yȳᵀ = C_constraints.make_ord4_con(d,Lx)
        print(d)
        print(size(L_xx̄ᵀ_tens_yȳᵀ["real"]))
        @constraint(model, L_xx̄ᵀ_tens_yȳᵀ["real"] .== real(ρ) )
        @constraint(model, L_xx̄ᵀ_tens_yȳᵀ["imag"] .== imag(ρ) )

## Localizing g constraint: L ≥ 0 on M₂ₜ(S)
        if con_list != ""
            if     occursin("S₁",con_list)
                loc_con      = C_constraints.make_loc_cons_S₁(ρ,d,t,Lx)
            elseif occursin("S₂",con_list)
                loc_con      = C_constraints.make_loc_cons_S₂(ρ,d,t,Lx)
            elseif occursin("S₃",con_list)
                loc_con, loc_con_eq  = C_constraints.make_loc_cons_S₃(ρ,d,t,Lx)
                for key in keys(loc_con_eq)
                    @constraint(model, loc_con_eq[key] .==  zeros(size(loc_con_eq[key])))
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
        end
## weak G Constraints
        # if occursin("wG",con_list)
        #     println("----------------Weak G-constraints are active")
        #     weakG_con = C_constraints.make_weakG_con(ρ,t,d,Lx)
        #     for key in keys(weakG_con)
        #         if isempty(weakG_con[key])
        #             continue
        #         end
        #         if size(weakG_con[key]) == (1, 1)
        #             @constraint(model, weakG_con[key] .>= 0)
        #         else
        #             @constraint(model, Symmetric(weakG_con[key]) in PSDCone())
        #         end
        #     end
        #
        # end
## G Constraints
        if  occursin("sG",con_list)
            println("----------------G-constraints are active")
            G_con                 = C_constraints.make_Gℝ_con(ρ,d,t,Lx)
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
