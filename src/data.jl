module Exmaples
export gen_rand_state

"""
generates a matrix of the form ∑aᵀa⊗bᵀb, with a,b random vectors.
"""
function gen_rand_state(d::Int64, r::Int64)
    a = [rand(1,d) for i in 1:r ]
    b = [rand(1,d) for i in 1:r ]
    ρ = zeros(d^2,d^2)
    for i in 1:r
        A = transpose(a[1])*a[1]
        B = transpose(b[1])*b[1]
        ρ = kron(A,B) + ρ
    end
    return ρ
end


end
