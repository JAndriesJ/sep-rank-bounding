module Examples
using LinearAlgebra
using Random
using DataFrames
using CSV

srcDir = dirname(@__FILE__)*"\\"
# include(srcDir *"Utils.jl")
include(srcDir *"Utils_states.jl")
include(srcDir *"Examples_sep.jl")
include(srcDir *"Examples_ent.jl")

# using .Utils
using .Utils_states
using .Examples_sep
using .Examples_ent

export get_examples,
       get_example_overview

"""Returns dict with three keys: ["rand", "sep", "ent"] each entry is also a dict."""
function get_examples()
   ρ         = Dict()
   ρ["sep"]  = Examples_sep.get_sep_example()
   ρ["ent"]  = Examples_ent.get_ent_example()
   return ρ
end


"""Returns a DataFrame object containing the"""
function get_example_overview( isWrite = false)
    ρ_dict = get_examples()
    df  = DataFrame( ex = String[], siz = Any[],isSep = String[], isReal = Bool[], rank = Int[],
                        isHerm = Bool[], isPPt = Bool[]);
    for key1 in keys(ρ_dict)
        seap = key1
        for key2 in keys(ρ_dict[key1])
            ρ      =  ρ_dict[key1][key2]
            isr    = isreal(ρ)
            ex     = key2
            vrank  = rank(ρ)
            isHerm = Hermitian(ρ) ≈ ρ

            d₁,d₂ = Utils_states.extr_d(key2)
            siz = d₁,d₂
            size(ρ)[1] == d₁*d₂
            if ~isHerm
                isPPt  = false
            else
                isPPt  =  Utils_states.isPPᵀ(ρ, [d₁,d₂] )
            end
            push!(df, [ex,siz,seap ,isr,vrank,isHerm,isPPt] )
        end
    end
    sort!(df,[:ex,:isSep])
    if isWrite
        CSV.write(pwd()*"\\ExMetaSep.md", df , delim = "|")
    end
    return df
end



end
