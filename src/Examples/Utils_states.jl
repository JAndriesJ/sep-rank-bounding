module Utils_states
using LinearAlgebra

export maketraceone,
       take_Pᵀ,
       isPPᵀ,
       ψ,
       sq,
       psep
       # makediagone,

""" Returns the partial transpose of the sysemt. """
function take_Pᵀ(ρ,sys,d)
    d₁,d₂ = d
    D = [(i₁,i₂,j₁,j₂)   for i₁ in 1:d₁, i₂ in 1:d₁, j₁ in 1:d₂, j₂ in 1:d₂]
    rcs =  [(i,j)   for j in 1:d₂  for i in 1:d₁]
    DB₀ = [D[ij[1],kl[1],ij[2],kl[2]] for ij in rcs , kl in rcs ]
    DB₁ = [D[kl[1],ij[1],ij[2],kl[2]] for ij in rcs , kl in rcs ]
    DB₂ = [D[ij[1],kl[1],kl[2],ij[2]] for ij in rcs , kl in rcs ]
    Zar = Dict(zip(DB₀, ρ))
    if sys == 1
        return map(a -> Zar[a],DB₁)
    elseif sys == 2
        return map(a -> Zar[a],DB₂)
    else
        return ρ
    end
end
"""
Checks if a sysetem is a PPT via checking if the smallest eigenvalues > -1.e-16
"""
function isPPᵀ(ρ,d)
    for sys ∈ 0:2
        ρᵀ = take_Pᵀ(ρ, sys, d)
        isPSDr = real(LinearAlgebra.eigvals(ρᵀ)[1]) > -1.e-15
        isPSDi = imag(LinearAlgebra.eigvals(ρᵀ)[1]) <  1.e-10
        ~(isPSDr && isPSDi) ? (return false) : 0
    end
    return true
end
"""Divides ρ by its trace"""
maketraceone(ρ) = ρ/tr(ρ)

""" diagM(1/√diag(ρ))*ρ*diagM(1/√diag(ρ))"""
# function makediagone(ρ)
#     mpli = 1 ./ map(x -> iszero(x) ? 1 : x,  sqrt.([ρ[i,i] for i in 1:size(ρ)[1]]))
#     return diagm(mpli)*ρ*diagm(mpli)
# end

eᵢ(n::Int,i::Int) = [j==i for j in 1:n] #"""The standard basis vector eᵢ in dimension n"""
ψ(i₀::Int,i₁::Int,n::Int) = kron(eᵢ(n,i₀),eᵢ(n,i₁))
ψ(n::Int) = kron(eᵢ(n,1) + eᵢ(n,2),eᵢ(n,1) + eᵢ(n,2))
sq(ϕ) = ϕ*transpose(conj(ϕ))
psep(i₀::Int,i₁::Int,n::Int) = sq(ψ(i₀,i₁,n))

end
