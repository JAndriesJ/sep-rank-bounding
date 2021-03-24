module Examples
using LinearAlgebra
using Random


export gen_rand_state,
       generate_random_states

makediagone(ρ) = diagm(1 ./ sqrt.(diag(ρ))) * ρ * diagm( 1 ./ sqrt.(diag(ρ)))

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

    return  ρ/r #makediagone(ρ)
end

"""
Generates some random matrices of the separabel structure:
Note: Fixed random seed.
"""
function generate_random_states(d_range,r_range,seed = 343)

    ρᵈʳ = Dict()
    for  d in d_range
        for r in r_range
            ρᵈʳ[d,r] = gen_rand_state(d,r,seed)
        end
    end
    return  ρᵈʳ
end



end
