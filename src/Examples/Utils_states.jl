module Utils_states
using LinearAlgebra

export maketraceone,
       take_Pᵀ,
       isPPᵀ,
       ψ,
       sq,
       psep,
       extr_d

"""Divides ρ by its trace"""
maketraceone(ρ) = ρ/tr(ρ)

"""
Returns the partial transpose of the sysemt.
Usage: partial_transpose_per(x ∈ ℝᵈˣᵈ, sys ∈ 1,2,...,m , dims = [d₁,d₂,...,dₘ])
"""
function take_Pᵀ(x, sys, dims)
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
Checks if a sysetem is a PPT via checking if the smallest eigenvalues > -1.e-16
"""
function isPPᵀ(x,dims)
    isPSDr = real(LinearAlgebra.eigvals(x)[1]) > -1.e-16
    isPSDi = imag(LinearAlgebra.eigvals(x)[1]) <  1.e-10
    if ~(isPSDr && isPSDi)
        return false
    end
    for sys ∈ length(dims)
        isPSDr = real(LinearAlgebra.eigvals(x)[1]) > -1.e-16
        isPSDi = imag(LinearAlgebra.eigvals(x)[1]) <  1.e-10
        if  ~(isPSDr && isPSDi)
            return false
        end
    end
    return true
end

"""The standard basis vector eᵢ in dimension n"""
eᵢ(n::Int,i::Int) = [j==i for j in 1:n]

function ψ(i₀::Int,i₁::Int,n::Int)
    e₀ = eᵢ(n,i₀)
    e₁ = eᵢ(n,i₁)
    return  kron(e₀,e₁)
end

function ψ(n::Int)
    e₀ = eᵢ(n,1)
    e₁ = eᵢ(n,2)
    e = e₀ + e₁
    return  kron(e,e)
end
sq(ϕ) = ϕ*transpose(conj(ϕ))
psep(i₀::Int,i₁::Int,n::Int) = sq(ψ(i₀,i₁,n))


## String tricks
"""
Extracts that dimensions of the systems from the name.
"""
function extr_d(str)
    pap = split(str[(end-9):end],"d")
    d₁ = parse(Int,pap[2][end])
    d₂ = parse(Int,pap[3][end])
    return d₁, d₂
end


end
