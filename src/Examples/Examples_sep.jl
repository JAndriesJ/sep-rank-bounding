module Examples_sep
using LinearAlgebra
using Random

srcDir = dirname(@__FILE__)*"\\"
include(srcDir *"Utils_states.jl")
using .Utils_states

export get_sep_examples,
       gen_ρ_rand,
       get_ρ_CDsep,
       get_ρ_DNY1,
       get_ρ_HHH1

## Separable States
"""
Separable state examples:
"""
function get_sep_example()
    ρ = Dict()
    ρ["real"] = Dict()
    ρ["imag"] = Dict()

    ρ["real"]["Rand25d₁4d₂4"]  = gen_ρ_rand(4, 25)
    ρ["real"]["Rand20d₁4d₂4"]  = gen_ρ_rand(4, 20)
    ρ["real"]["Rand15d₁4d₂4"]  = gen_ρ_rand(4, 15)
    ρ["real"]["CDsep3d₁3d₂3"] = get_ρ_CDsep(3)
    ρ["real"]["CDsep4d₁3d₂3"] = get_ρ_CDsep(4)
    ρ["real"]["CDsep5d₁3d₂3"] = get_ρ_CDsep(5)
    ρ["real"]["DNY1d₁3d₂3"]   = get_ρ_DNY1()


    ρ["imag"]["Randℂ20d₁4d₂4"]  = gen_ρ_rand(4, 20,true)
    ρ["imag"]["Randℂ15d₁4d₂4"]  = gen_ρ_rand(4, 15,true)
    ρ["imag"]["Randℂ10d₁4d₂4"]  = gen_ρ_rand(4, 10,true)
    ρ["imag"]["Randℂ5d₁4d₂4"]  = gen_ρ_rand(4, 5,true)

    # ρ["HHH1d₁2d₂2"]   = get_ρ_HHH1(0.5)
    # ρ["CDsep1d₁2d₂2"] = get_ρ_CDsep(1)
    # ρ["CDsep2d₁2d₂2"] = get_ρ_CDsep(2)

    ρ_new = Dict()
    for k in ["real", "imag"]
        ρ_new[k] = Dict()
        for key in keys(ρ[k])
            ρ_new[k][key] = Utils_states.maketraceone(ρ[k][key])
        end
    end

    return ρ_new
end

## Examples
"""generates a matrix of the form ∑ʳaᵀa⊗bᵀb, with a,b ∈ uniform random entry-wise.  d ∈ {2,3,4} , r ∈ [9]"""
function gen_ρ_rand(d::Integer, r::Integer, isℂ = false)
    # seed = 343
    # Random.seed!(seed)
    if isℂ
        a = [randn(1,d) + im*randn(1,d) for i in 1:r ]
        b = [randn(1,d) + im*randn(1,d) for i in 1:r ]
    else
        a = [randn(1,d) for i in 1:r ]
        b = [randn(1,d) for i in 1:r ]
    end
    ρ = zeros(d^2,d^2)
    for i in 1:r
        A = transpose(a[i])*conj(a[i])
        B = transpose(b[i])*conj(b[i])
        ρ = kron(A,B) + ρ
    end
    return ρ
end

""" http://arxiv.org/abs/1210.0111v2"""
function get_ρ_CDsep(i)
    ρ = Dict()
    # |00><00| + |11><11|
    ρ[1] =  psep(1,1,2) + psep(2,2,2)

    # |00><00| + |11><11| + (|0> + |1>)⊗(|0> + |1>)(<0| + <1|)⊗(<0| + <1|)
    ρ[2] = psep(1,1,2) + psep(2,2,2) + sq(ψ(2))

    # |00><00| + |11><11| + |12><12|
    ρ[3] = psep(1,1,3) + psep(2,2,3) + psep(2,3,3)

    # |00><00| + |01><01| + |11><11| + |12><12|
    ρ[4] = psep(1,1,3) + psep(1,2,3) + psep(2,2,3) + psep(2,3,3)

    # |00><00| + |01><01| + |02><02| + |11><11| + |12><12|
    ρ[5] = psep(1,1,3) + psep(1,2,3) +  psep(1,3,3) + psep(2,2,3) + psep(2,3,3)
    return ρ[i]
end

"""https://arxiv.org/abs/2011.08132v1"""
function get_ρ_DNY1()
    H_flat = [i₁*j₁ + i₂*j₂ for i₁ in 1:3 for i₂ in 1:3 for j₁ in 1:3   for j₂ in 1:3 ]
    return  reshape(H_flat,9,9)
end
""" http://arxiv.org/abs/quant-ph/9605038v2 """
function get_ρ_HHH1(p  = 0.5)
    # Random.seed!(343)
    a  = 0.5; b = 0.8; # arbitary possitive numbers
    ψ₁ = a*Utils_states.ψ(2,2,2) + b*Utils_states.ψ(1,1,2)
    ψ₂ = a*Utils_states.ψ(2,1,2) + b*Utils_states.ψ(1,2,2)
    ρ_temp = p*Utils_states.sq(ψ₁) + (1 - p)*Utils_states.sq(ψ₂)
    return ρ_temp
end

end
