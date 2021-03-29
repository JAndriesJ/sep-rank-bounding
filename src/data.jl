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

"""
Makes the entries on the diagonal equal to one via scaling
"""
function makediagone(ρ)
    norm_vec = 1 ./ sqrt.(diag(ρ))
    norm_vec[norm_vec.==Inf] .= 1 # Incase there are zeros on the diagonal
    ρ = diagm(norm_vec ) * ρ * diagm(norm_vec)
end
"""
Makes the trace equal to one via scaling
"""
maketraceone(ρ) = ρ/tr(ρ)

"""
 Wraps all the normalizations I use.
"""
function multi_norma(ρ, norm)
    if  norm == "tr"
        return  maketraceone(ρ)
    elseif norm == "dg"
        return  makediagone(ρ)
    else norm == "no"
        return ρ
    end
end

#
# """
# generates a matrix of the form ∑ʳaᵀa⊗bᵀb, with a,b ∈ [0,1]ᵈ uniform random entry-wise.  d ∈ {2,3,4} , r ∈ [9]
# """
function gen_rand_state(d::Integer, r::Integer,norm::String)
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
    return multi_norma(ρ,norm)
end

"""
Generates some random matrices of the separabel structure:
Note: Fixed random seed.
"""
function generate_random_states(d_range,r_range, norm)

    ρᵈʳ = Dict()
    for  d in d_range
        for r in r_range
            ρᵈʳ[d,r,"rand"] = gen_rand_state(d,r, norm)
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
psep(i₀::Int,i₁::Int,n::Int) = sq(ψ(i₀,i₁,n))


"""
Separable state examples:
"""

function get_sep_example(norma)
    ρ = Dict()
    ##  http://arxiv.org/abs/1210.0111v2 Table I
    # |00><00| + |11><11|
    ρ[2,2,"sep1"] = psep(1,1,2) + psep(2,2,2)

    # |00><00| + |11><11| + |12><12|
    ρ[3,3,"sep2"] = psep(1,1,3) + psep(2,2,3) + psep(2,3,3)

    # |00><00| + |01><01| + |11><11| + |12><12|
    ρ[3,4,"sep3"] = psep(1,1,3) + psep(1,2,3) + psep(2,2,3) + psep(2,3,3)

    # |00><00| + |01><01| + |02><02| + |11><11| + |12><12|
    ρ[3,5,"sep4"] =  psep(1,1,3) + psep(1,2,3) +  psep(1,3,3) + psep(2,2,3) + psep(2,3,3)


    # Give errors
    # |00><00| + |11><11| + |01><10|
    ψₑₑ = ψ(2)
    ρ[2,3,"sep5"] = psep(1,1,2) + psep(2,2,2) + sq(ψₑₑ)

    # = |00><00| + |02><02| + 2|11><11| + (|01> + |10>)(|01> + <10|)
    ψ₀₁ = ψ(1,2,3)
    ψ₁₀ = ψ(2,1,3)
    ρ[3,5,"sep6"] = psep(1,1,3) + psep(1,3,3) + 2*psep(2,2,3) + sq(ψ₀₁ + ψ₁₀)

    for key in keys(ρ)
        ρ[key] =  multi_norma(ρ[key],norma)
    end
    return ρ
end


"""
Entangled state examples:
"""
function get_ent_example()
    ρ = Dict()
    ## Entangled states examples
    # Example 1: https://en.wikipedia.org/wiki/Quantum_entanglement
    e₀ = get_std_base_vec(2,1)
    e₁ = get_std_base_vec(2,2)
    ρ[2,"ent1"] = (1/sqrt(2))*(e₀*transpose(e₁) - e₁*transpose(e₀))

    # Example 2: http://arxiv.org/abs/quant-ph/9605038v2  page 9 eq. (22)
    a  = 0.5; b = 0.8; # arbitary possitive numbers
    p  = rand()
    ψ₁ = a*ψ(2,2,2) + b*ψ(1,1,2)
    ψ₂ = a*ψ(2,1,2) + b*ψ(1,2,2)
    ρ[2,"e2"] = p*sq(ψ₁) + (1 - p)*sq(ψ₂)

    # Example 3: http://arxiv.org/abs/quant-ph/9605038v2  page 10 eq. (25)
    a =  1/sqrt(2) ; p = rand();
    ψₐ = a*(ψ(1,2,2) - ψ(2,1,2))
    ρ[2,"e3"] = p*sq(ψₐ) + (1 - p)*psep(1,1,2)


    for key in keys(ρ)
        # ρ[key] = maketraceone(ρ[key])
        # ρ[key] = makediagone(ρ[key])
    end
    return ρ
end



end
