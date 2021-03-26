module Examples
using LinearAlgebra
using Random
include(pwd()*"\\src\\moments.jl")
using .moments
# include(sourceDir*"moments.jl")
# using .moments

export gen_rand_state,
       generate_random_states,
       get_sep_example,
       get_ent_example

# This does not work well when there are zeros on the diagonal
makediagone(ρ)  = diagm(1 ./ sqrt.(diag(ρ)) ) * ρ * diagm(1 ./ sqrt.(diag(ρ)) )
maketraceone(ρ) = ρ/tr(ρ)


"""
generates a matrix of the form ∑aᵀa⊗bᵀb, with a,b random vectors.
"""
function gen_rand_state(d::Integer, r::Integer,seed = 343)
    Random.seed!(seed)
    a = [rand(1,d) for i in 1:r ]
    b = [rand(1,d) for i in 1:r ]
    ρ = zeros(d^2,d^2)
    for i in 1:r
        A = transpose(a[i])*a[i]
        B = transpose(b[i])*b[i]
        ρ = kron(A,B) + ρ
    end
    # TODO make the trace equal to 1
    return  maketraceone(ρ) #makediagone(ρ)
end

"""
Generates some random matrices of the separabel structure:
Note: Fixed random seed.
"""
function generate_random_states(d_range = 2:4,r_range = 2:9,seed = 343)
    ρᵈʳ = Dict()
    for  d in d_range
        for r in r_range
            ρᵈʳ[d,r] = gen_rand_state(d,r,seed)
        end
    end
    return  ρᵈʳ
end


function ψ(i₀::Int,i₁::Int,n::Int)
    e₀ = get_std_base_vec(n,i₀)
    e₁ = get_std_base_vec(n,i₁)
    return  kron(e₀,e₁)
end

function ψ(n::Int)
    e₀ = get_std_base_vec(n,1)
    e₁ = get_std_base_vec(n,2)
    e = e₀ + e₁
    return  kron(e,e)
end
sq(ϕ) = ϕ*transpose(ϕ)
sqs(ϕ₁,ϕ₂) = sq(ϕ₁) + sq(ϕ₂)
psep(i₀::Int,i₁::Int,n::Int) = sq(ψ(i₀,i₁,n))


"""
Separable state examples:
"""

function get_sep_example(d = 3)
    ρ = Dict()
    ## Separable: http://arxiv.org/abs/1210.0111v2 Table I
    # |00><00| + |11><11|
    ρ[d,2,"s"] = psep(1,1,d) + psep(2,2,d)

    # http://arxiv.org/abs/1210.0111v2 Table II
    # |00><00| + |11><11| + |12><12|
    ρ[d,3,"s"] = psep(1,1,d) + psep(2,2,d) + psep(2,3,d)

    # http://arxiv.org/abs/1210.0111v2 Table II
    # |00><00| + |01><01| + |11><11| + |12><12|
    ρ[d,4,"s"] = psep(1,1,d) + psep(1,2,d) + psep(2,2,d) + psep(2,3,d)

    #  http://arxiv.org/abs/1210.0111v2  Example 13
    # = |00><00| + |02><02| + 2|11><11| + (|01> + |10>)(|01> + <10|)
    ψ₀₁ = ψ(1,2,d)
    ψ₁₀ = ψ(2,1,d)
    ρ[d,5,"s"] = psep(1,1,d) + psep(1,3,d) + 2*psep(2,2,d) + sq(ψ₀₁ + ψ₁₀)

    # http://arxiv.org/abs/1210.0111v2 Table II
    # |00><00| + |01><01| + |02><02| + |11><11| + |12><12|
    ρ[d,5,"s"] =  psep(1,1,d) + psep(1,2,d) +  psep(1,3,d) + psep(2,2,d) + psep(2,3,d)

    for key in keys(ρ)
        ρ[key] = maketraceone(ρ[key])
    end
    return ρ
end


"""
Entangled state examples:
"""
function get_ent_example(d = 3)
    ρ = Dict()
    ## Entangled states examples
    # Example 1: https://en.wikipedia.org/wiki/Quantum_entanglement
    # ρ = (1/sqrt(2))*(e₀*transpose(e₁) - e₁*transpose(e₀))

    # Example 1:  http://arxiv.org/abs/1210.0111v2 Table I
    # |00><00| + |11><11| + |01><10| (TODO Double Check this)
    ψₑₑ = ψ(d)
    ρ[d,"e1"] = psep(1,1,d) + psep(2,2,d) + sq(ψₑₑ)

    # Example 2: http://arxiv.org/abs/quant-ph/9605038v2  page 9 eq. (22)
    a  = 0.5; b = 0.8; # arbitary possitive numbers
    p  = rand()
    ψ₁ = a*ψ(2,2,d) + b*ψ(1,1,d)
    ψ₂ = a*ψ(2,1,d) + b*ψ(1,2,d)
    ρ[d,"e2"] = p*sq(ψ₁) + (1 - p)*sq(ψ₂)

    # Example 3: http://arxiv.org/abs/quant-ph/9605038v2  page 10 eq. (25)
    a =  1/sqrt(2) ; p = rand();
    ψₐ = a*(ψ(1,2,d) - ψ(2,1,d))
    ρ[d,"e3"] = p*sq(ψₐ) + (1 - p)*psep(1,1,d)

    for key in keys(ρ)
        ρ[key] = maketraceone(ρ[key])
    end
    return ρ
end



end
