module Examples_sep
using LinearAlgebra
using Random

srcDir = dirname(@__FILE__)*"\\"
include(srcDir *"Utils_states.jl")
using .Utils_states

export get_sep_examples

## Separable States
"""
Separable state examples:
"""
function get_sep_example()
    ρ = Dict()
##
    ρ["Ex3-8DNY20"]   = get_ρ_DNY20()

##  CD12
    ρ["T1CD12i"]    = get_ρ_T1CD12(1)
    ρ["T1CD12ii"]   = get_ρ_T1CD12(2)
    # ρ["T1CD12iii"]  = get_ρ_T1CD12(3)

    ρ["T2CD12i"]    = get_ρ_T2CD12(1)
    ρ["T2CD12ii"]   = get_ρ_T2CD12(2)
    # ρ["T2CD12iii"]  = get_ρ_T2CD12(3)
    # ρ["T2CD12iv"]  = get_ρ_T2CD12(4)
    ρ["T2CD12v"]    = get_ρ_T2CD12(5)

    # ρ["Ex26CD12"]   = get_ρ_Ex26CD12() # Wrong shape

## Random batch
 ## Real
    ρ["RANDa1"]  = gen_ρ_RAND(3, 5)
    ρ["RANDa2"]  = gen_ρ_RAND(3, 10)
    ρ["RANDa3"]  = gen_ρ_RAND(3, 15)
    ρ["RANDa4"]  = gen_ρ_RAND(3, 20)
    ρ["RANDa5"]  = gen_ρ_RAND(3, 25)

    ρ["RANDb1"]  = gen_ρ_RAND(4, 5)
    ρ["RANDb2"]  = gen_ρ_RAND(4, 10)
    ρ["RANDb3"]  = gen_ρ_RAND(4, 15)
    ρ["RANDb4"]  = gen_ρ_RAND(4, 20)
    ρ["RANDb5"]  = gen_ρ_RAND(4, 25)
 ## Complex
    ρ["RANDℂa1"]  = gen_ρ_RAND(3,  5,true)
    ρ["RANDℂa2"]  = gen_ρ_RAND(3, 10,true)
    ρ["RANDℂa3"]  = gen_ρ_RAND(3, 15,true)
    ρ["RANDℂa4"]  = gen_ρ_RAND(3, 20,true)
    ρ["RANDℂa5"]  = gen_ρ_RAND(3, 25,true)

    ρ["RANDℂb1"]  = gen_ρ_RAND(4,  5,true)
    ρ["RANDℂb2"]  = gen_ρ_RAND(4, 10,true)
    ρ["RANDℂb3"]  = gen_ρ_RAND(4, 15,true)
    ρ["RANDℂb4"]  = gen_ρ_RAND(4, 20,true)
    ρ["RANDℂb5"]  = gen_ρ_RAND(4, 25,true)


    for key in keys(ρ)
        ρ[key] = Utils_states.maketraceone(ρ[key])
    end
    return ρ
end

## Examples
"""generates a matrix of the form ∑ʳaᵀa⊗bᵀb, with a,b ∈ uniform random entry-wise.  d ∈ {2,3,4} , r ∈ [9]"""
function gen_ρ_RAND(d::Integer, r::Integer, isℂ = false)
    seed = 343
    Random.seed!(seed)
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
function get_ρ_T1CD12(i)
    ρ = Dict()
    # |00><00| + |11><11|
    ρ[1] =  psep(1,1,2) + psep(2,2,2)

    # |00><00| + |11><11| + (|0> + |1>)⊗(|0> + |1>)(<0| + <1|)⊗(<0| + <1|)
    ρ[2] = psep(1,1,2) + psep(2,2,2) + sq(ψ(2))


    return ρ[i]
end

""" http://arxiv.org/abs/1210.0111v2"""
function get_ρ_T2CD12(i)
    ρ = Dict()

    # |00><00| + |11><11| + |12><12|
    ρ[1] = psep(1,1,3) + psep(2,2,3) + psep(2,3,3)

    # |00><00| + |01><01| + |11><11| + |12><12|
    ρ[2] = psep(1,1,3) + psep(1,2,3) + psep(2,2,3) + psep(2,3,3)

    # |00><00| + |01><01| + |02><02| + |11><11| + |12><12|
    ρ[5] = psep(1,1,3) + psep(1,2,3) +  psep(1,3,3) + psep(2,2,3) + psep(2,3,3)

    return ρ[i]
end

""" http://arxiv.org/abs/1210.0111v2"""
function get_ρ_Ex26CD12()

    ρ = [2 0 0 0 0 0
         0 4 2 0 0 2
         0 2 2 1 -1 0
         0 0 1 2 1 -1
         0 0 -1 1 5 1
         0 2 0 0 1 2]

    return ρ
end

"""https://arxiv.org/abs/2011.08132v1"""
function get_ρ_DNY20()
    H_flat = [i₁*j₁ + i₂*j₂ for i₁ in 1:3 for i₂ in 1:3 for j₁ in 1:3   for j₂ in 1:3 ]
    return  reshape(H_flat,9,9)
end

# """ http://arxiv.org/abs/quant-ph/9605038v2 """
# function get_ρ_HHH1(p  = 0.5)
#     # Random.seed!(343)
#     a  = 0.5; b = 0.8; # arbitary possitive numbers
#     ψ₁ = a*Utils_states.ψ(2,2,2) + b*Utils_states.ψ(1,1,2)
#     ψ₂ = a*Utils_states.ψ(2,1,2) + b*Utils_states.ψ(1,2,2)
#     ρ_temp = p*Utils_states.sq(ψ₁) + (1 - p)*Utils_states.sq(ψ₂)
#     return ρ_temp
# end

end
