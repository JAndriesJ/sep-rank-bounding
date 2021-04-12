module Examples
using LinearAlgebra
using Random

srcDir = dirname(@__FILE__)*"\\"
include(srcDir *"Utils.jl")
using .Utils


export get_examples,
       partial_transpose_per

"""
Returns dict with three keys: ["rand", "sep", "ent"]
each entry is also a dict.
"rand": ["randd3r3"  "randd3r6"  "randd2r2"  "randd3r5"  "randd3r1"  "randd2r1"  "randd2r6"  "randd3r7"  "randd2r7"  "randd2r5"  "randd2r4"  "randd3r2"  "randd3r4"  "randd2r3"]
"sep":  ["sep4d3r4",  "sep5d3r5",  "sep3d3r3",  "sep1d2r2",  "sep2d2r3",  "sep6d3r2"]
"ent":  ["ent4d3r4",  "ent2d2r2",  "ent3d2r2",  "ent1d2r1",  "ent5d2r2"]
"""
function get_examples()
   ρ         = Dict()
   ρ["rand"] = get_rand_states(2:3,1:7)
   ρ["sep"]  = get_sep_example()
   ρ["ent"]  = get_ent_example()
   return ρ
end


function partial_transpose_per(x, sys, dims)
	  n = length(dims)
	  d = prod(dims)
	  s = n - sys + 1
	  p = collect(1:2n)
	  p[s] = n + s
	  p[n + s] = s
	  rdims = reverse(dims)
	  r = reshape(x, (rdims..., rdims...))
	  return reshape(permutedims(r,p),(d,d))
end


"""
Makes the trace equal to one via scaling
"""
maketraceone(ρ) = ρ/tr(ρ)


## Random matrices
"""
generates a matrix of the form ∑ʳaᵀa⊗bᵀb, with a,b ∈ [0,1]ᵈ uniform random entry-wise.  d ∈ {2,3,4} , r ∈ [9]
"""
function gen_rand_state(d::Integer, r::Integer)
    seed = 343
    Random.seed!(seed)
    a = [rand(1,d) for i in 1:r ]
    b = [rand(1,d) for i in 1:r ]
    ρ = zeros(d^2,d^2)
    for i in 1:r
        A = transpose(a[i])*a[i]
        B = transpose(b[i])*b[i]
        ρ = kron(A,B) + ρ
    end
    return maketraceone(ρ)
end

"""
Generates some random matrices of the separabel structure:
Note: Fixed random seed.
"""
function get_rand_states(d_range,r_range)
    ρᵈʳ = Dict()
    for  d in d_range
        for r in r_range
            ρᵈʳ_temp = gen_rand_state(d,r)
            ρᵈʳ["rand"*"d$d"*"r$r"] = ρᵈʳ_temp
        end
    end
    return  ρᵈʳ
end


## Separable States

function ψ(i₀::Int,i₁::Int,n::Int)
    e₀ = Utils.eᵢ(n,i₀)
    e₁ = Utils.eᵢ(n,i₁)
    return  kron(e₀,e₁)
end

function ψ(n::Int)
    e₀ = Utils.eᵢ(n,1)
    e₁ = Utils.eᵢ(n,2)
    e = e₀ + e₁
    return  kron(e,e)
end

sq(ϕ) = ϕ*transpose(ϕ)
psep(i₀::Int,i₁::Int,n::Int) = sq(ψ(i₀,i₁,n))


"""
Separable state examples:
"""
function get_sep_example()
    ρ = Dict()
    ##  http://arxiv.org/abs/1210.0111v2 Table I
    # sep 1
    # |00><00| + |11><11|
    ρ_temp = psep(1,1,2) + psep(2,2,2)
    ρ["sep1d2r2"] = ρ_temp

    # sep 2
    # |00><00| + |11><11| + (|0> + |1>)⊗(|0> + |1>)(<0| + <1|)⊗(<0| + <1|)
    ψₑₑ = ψ(2)
    ρ_temp = psep(1,1,2) + psep(2,2,2) + sq(ψₑₑ)
    ρ["sep2d2r3"] = ρ_temp

    # sep 3
    # |00><00| + |11><11| + |12><12|
    ρ_temp = psep(1,1,3) + psep(2,2,3) + psep(2,3,3)
    ρ["sep3d3r3"] = ρ_temp

    # sep 4
    # |00><00| + |01><01| + |11><11| + |12><12|
    ρ_temp = psep(1,1,3) + psep(1,2,3) + psep(2,2,3) + psep(2,3,3)
    ρ["sep4d3r4"] = ρ_temp

    # sep 5
    # |00><00| + |01><01| + |02><02| + |11><11| + |12><12|
    ρ_temp = psep(1,1,3) + psep(1,2,3) +  psep(1,3,3) + psep(2,2,3) + psep(2,3,3)
    ρ["sep5d3r5"] =  ρ_temp

    # sep 6
    #  Hᵢ₁ᵢ₂ⱼ₁ⱼ₂ = i₁j₁ + i₂j₂
    H_flat = [i₁*j₁ + i₂*j₂ for i₁ in 1:3 for i₂ in 1:3 for j₁ in 1:3   for j₂ in 1:3 ]
    ρ_temp = reshape(H_flat,9,9)
    ρ["sep6d3r2"] = ρ_temp

    # sep 6 fact
    #  Hᵢ₁ᵢ₂ⱼ₁ⱼ₂ = i₁j₁ + i₂j₂
    # λ₁ u¹₁ (u¹₁)ᵀ ⊗ u²₁(u²₁)ᵀ  + λ₂ u¹₂ (u¹₂)ᵀ ⊗ u²₂(u²₂)ᵀ
    # λ₁ = λ₂ = 42 ;
    # u¹₁ = u²₂ = [sqrt(14)/14, sqrt(14)/7, 3/sqrt(14)] ;
    # u²₁ = u¹₂ = [sqrt(3)/3, sqrt(3)/3, sqrt(3)/3] ;
    # H_2 = λ₁ * kron(u¹₁ *transpose(u¹₁) , u²₁*transpose(u²₁))  + λ₂ * kron(u¹₂ *transpose(u¹₂), u²₂*transpose(u²₂))

    ρ_new = Dict()
    for key in keys(ρ)
        ρ_new[key] = maketraceone(ρ[key])
    end

    return ρ_new
end

## Entangled States

"""
Entangled state examples:
"""
function get_ent_example()
    ρ = Dict()
    ## Entangled states examples
    # # Example 1: https://en.wikipedia.org/wiki/Quantum_entanglement
    # e₀ = Utils.eᵢ(2,1)
    # e₁ = Utils.eᵢ(2,2)
    # ρ[2,"ent1"] = (1/sqrt(2))*(e₀*transpose(e₁) - e₁*transpose(e₀))
    ϕ = [1, 0, 0, 1]

    ρ_temp = ϕ*transpose(ϕ)
    ρ["ent1d2r1"] = ρ_temp

    # Example 2: http://arxiv.org/abs/quant-ph/9605038v2  page 9 eq. (22)
    Random.seed!(343)
    a  = 0.5; b = 0.8; # arbitary possitive numbers
    p  = rand()
    ψ₁ = a*ψ(2,2,2) + b*ψ(1,1,2)
    ψ₂ = a*ψ(2,1,2) + b*ψ(1,2,2)
    ρ_temp = p*sq(ψ₁) + (1 - p)*sq(ψ₂)
    ρ["ent2d2r2"] = ρ_temp

    # Example 3: http://arxiv.org/abs/quant-ph/9605038v2  page 10 eq. (25) fails PPT
    a =  1/sqrt(2) ; p = rand();
    ψₐ = a*(ψ(1,2,2) - ψ(2,1,2))
    ρ_temp = p*sq(ψₐ) + (1 - p)*psep(1,1,2)
    ρ["ent3d2r2"] = ρ_temp

    # ent 4
    # = |00><00| + |02><02| + 2|11><11| + (|01> + |10>)(<01| + <10|)
    # C²⊗C³
    ψ₀₁ = ψ(1,2,3)
    ψ₁₀ = ψ(2,1,3)
    ρ_temp = psep(1,1,3) + psep(1,3,3) + 2*psep(2,2,3) + sq(ψ₀₁ + ψ₁₀)
    ρ["ent4d3r4"] = ρ_temp


    # https://arxiv.org/abs/2011.08132v1 Example 3.7
    # ent 5
    # Hᵢ₁ᵢ₂ⱼ₁ⱼ₂ = i₁ + i₂ + j₁ + j₂
    H_flat = [i₁ + i₂ + j₁ + j₂ for i₁ in 1:2 for i₂ in 1:2 for j₁ in 1:2 for j₂ in 1:2 ]
    ρ_temp = reshape(H_flat,4,4)
    ρ["ent5d2r2"] = ρ_temp


    for key in keys(ρ)
        ρ[key] =  maketraceone(ρ[key])
    end
    return ρ
end


end
