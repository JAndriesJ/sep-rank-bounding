module Examples_ent
using LinearAlgebra
using Random

srcDir = dirname(@__FILE__)*"\\"
include(srcDir *"Utils_states.jl")
using .Utils_states
const us = Utils_states


export get_ent_examples

## Entangled States
"""
Entangled state examples:
"""
function get_ent_example()
    ρ = Dict()

    ρ["Eq22HHH96a"] = get_ρ_Eq22HHH96(0.25,0.5,0.25)
    ρ["Eq22HHH96b"] = get_ρ_Eq22HHH96(0.5,0.25,0.25)
    ρ["Eq22HHH96c"] = get_ρ_Eq22HHH96(0.75,0.5,0.85)

    ρ["Eq25HHH96"]  = get_ρ_Eq25HHH96(0.3)

    seed = 343
    Random.seed!(seed)
    ρ["Eq2-3BP00a"] = get_ρ_Eq2_3BP00(randn(8) + im *randn(8)...)
    ρ["Eq2-3BP00b"] = get_ρ_Eq2_3BP00(randn(8) + im *randn(8)...)
    ρ["Eq2-3BP00c"] = get_ρ_Eq2_3BP00(randn(8) + im *randn(8)...)
    ρ["Eq2-3BP00d"] = get_ρ_Eq2_3BP00(randn(8) + im *randn(8)...)
    ρ["Eq2-3BP00e"] = get_ρ_Eq2_3BP00(randn(8) + im *randn(8)...)

    ρ["Eq14Hor97a"] = get_ρ_Eq14Hor97(0.25)
    ρ["Eq14Hor97b"] = get_ρ_Eq14Hor97(0.5)
    ρ["Eq14Hor97c"] = get_ρ_Eq14Hor97(0.75)

    # ρ["Eq32Hor97a"] = get_ρ_Eq32Hor97(0.25) # Wrong shape
    # ρ["Eq32Hor97b"] = get_ρ_Eq32Hor97(0.5)
    # ρ["Eq32Hor97c"] = get_ρ_Eq32Hor97(0.75)

    ρ["p12HK14a"]   = get_ρ_p12HK14(0.5)
    ρ["p12HK14b"]   = get_ρ_p12HK14(2)
    ρ["p12HK14c"]   = get_ρ_p12HK14(3)

    for key in keys(ρ)
        ρ[key] = Utils_states.maketraceone(ρ[key])
    end


    return ρ
end


## Examples

"""
##    Separability of Mixed States: Necessary and Sufficient Conditions
##    Michal Horodecki, Pawel Horodecki, Ryszard Horodecki
## http://arxiv.org/abs/quant-ph/9605038v2  page 9 eq. (22)
"""
function get_ρ_Eq22HHH96(a,b,p)
    @assert 0 < a
    @assert 0 < b
    @assert a != b
    @assert p != 0.5
    ψ₁ = a*us.ψ(2,2,2) + b*us.ψ(1,1,2)
    ψ₂ = a*us.ψ(2,1,2) + b*us.ψ(1,2,2)
    ρ_temp = p*us.sq(ψ₁) + (1 - p)*us.sq(ψ₂)
    return ρ_temp
end

"""
##   page 10 eq. (25) fails PPT
"""
function get_ρ_Eq25HHH96(p)
    @assert p != 0
    a =  1/sqrt(2) ;
    ψₐ = a*(us.ψ(1,2,2) - us.ψ(2,1,2))
    ρ_temp = p*us.sq(ψₐ) + (1 - p)*us.psep(1,1,2)
    return ρ_temp
end

"""
## 1999, Bruß,Peres,Construction of quantum states with bound entanglement.pdf
## https://arxiv.org/abs/quant-ph/9911056v1
"""
function get_ρ_Eq2_3BP00(m,s,n,a,b,c,t,d)
    # randn(8) + im *randn(8))
    V₁ =   [m  0  s 0  n  0  0  0  0]
    V₂ =   [0  a  0 b  0  c 0  0  0]
    V₃ =   [conj(n)  0   0  0  -conj(m)  0    t      0     0]
    V₄ = [0    conj(b)  0 -conj(a)   0   0   0      d      0]
    ρ =  adjoint(V₁)*V₁  + adjoint(V₂)*V₂  +adjoint(V₃)*V₃  + adjoint(V₄)*V₄
    return ρ
end

function get_ρ_Eq14Hor97(a)
    @assert a > 0
    @assert a < 1
    b = (1+a)/2
    c = (sqrt(1-a^2)/2)
    ρ =     (1/(8*a+1))*[a 0 0 0 a 0 0 0 a
                         0 a 0 0 0 0 0 0 0
                         0 0 a 0 0 0 0 0 0
                         0 0 0 a 0 0 0 0 0
                         a 0 0 0 a 0 0 0 a
                         0 0 0 0 0 a 0 0 0
                         0 0 0 0 0 0 b 0 c
                         0 0 0 0 0 0 0 a 0
                         a 0 0 0 a 0 c 0 b]

    return ρ
end

function get_ρ_Eq32Hor97(b)
    @assert b > 0
    @assert b < 1
    a = (1+b)/2
    c = (sqrt(1-b^2)/2)
    ρ =     (1/(7*b+1))*[b 0 0 0 0 b 0 0
                         0 b 0 0 0 0 b 0
                         0 0 b 0 0 0 0 b
                         0 0 0 b 0 0 0 0
                         0 0 0 0 a 0 0 c
                         b 0 0 0 0 b 0 0
                         0 b 0 0 0 0 b 0
                         0 0 b 0 c 0 0 a]

    return ρ
end

"""
##Separable states with unique decompositions,
## Kil-Chan Ha, Seung-Hyeok Kye
## https://arxiv.org/abs/1210.1088v5
"""
function get_ρ_p12HK14(b)
    @assert b >  0
    @assert b != 1.0
    c = 3*(1+b+1/b)
    ρ =  (1/c)*[1 0 0 0 1 0 0 0 1
                0 b 0 1 0 0 0 0 0
                0 0 1/b 0 0 0 1 0 0
                0 1 0 1/b 0 0 0 0 0
                1 0 0 0 1 0 0 0 1
                0 0 0 0 0 b 0 1 0
                0 0 1 0 0 0 b 0 0
                0 0 0 0 0 1 0 1/b 0
                1 0 0 0 1 0 0 0 1]
    return ρ
end

function get_ρ_Eq5TAH12()

end

function get_ρ_Eq6CR11(a,λ₁,λ₂)

end


"""
## 2013,Chen,Dokovic,Dimensions, lengths and separability in finite-dimensional quantum systems
## https://arxiv.org/abs/1206.3775v4
"""
function get_ρ_CD1(N)
    n₁ = 2
    e₀  = us.eᵢ(n₁,1)
    e₁  = us.eᵢ(n₁,2)

    f = Dict()
    for i ∈ 1:N
        f[i-1] = us.eᵢ(N,i)
    end

    aₖ = Dict()
    bₖ = Dict()
    for k ∈ 1:N
        aₖ[k] = e₀ + (k - 1)*e₁
        bₖ[k] = f[0] + f[k - 1]
    end
    for k ∈ (N+1):(2*N - 2)
        aₖ[k] = e₀ + (k - 1)*e₁
        bₖ[k] = f[0] + f[2*N - k - 1] + f[2*N - k]
    end
    aₖ[2*N - 1] = e₀ + (N - 1)*im*e₁
    bₖ[2*N - 1] = im*f[0] + f[N - 1]

    aₖ[2*N] = e₀
    bₖ[2*N] = f[N-1]

    ρ = zeros(2*N,2*N)
    for key in keys(aₖ)
        c = kron(us.sq(aₖ[key]),us.sq(bₖ[key]))
        ρ = ρ + c
    end
    return ρ
end

function get_ρ_CD2()
    n₁ = 3
    n₂ = 4
    e₀  = us.eᵢ(n₁,1)
    e₁  = us.eᵢ(n₁,2)
    e₂  = us.eᵢ(n₁,3)

    f₀  = us.eᵢ(n₂,1)
    f₁  = us.eᵢ(n₂,2)
    f₂  = us.eᵢ(n₂,3)
    f₃  = us.eᵢ(n₂,4)

    a = Dict()
    b = Dict()

    a[1] = e₀
    a[2] = e₀ + e₁
    a[3] = e₀ - e₁
    a[4] = e₀ + im*e₁
    a[5] = e₀ + e₂
    a[6] = e₀ - e₂
    a[7] = e₀ + e₁ + e₂
    a[8] = e₀ - e₁ + e₂
    a[9] = e₀ + (1 + im)*e₂
    a[10] = e₀ + im*e₁ - e₂
    a[11] = e₀ + e₁ + im*e₂
    a[12] = e₀ + im*e₁ + e₂
    a[13] = e₀ + im*e₁ + im*e₂
    a[14] = e₀ - im*e₁


    b[1] = f₀
    b[2] = f₀ + f₁
    b[3] = f₀ + f₂
    b[4] = f₀ - f₃
    b[5] = f₀ + f₃
    b[6] = f₀ + f₁ + f₃
    b[7] = f₀ - f₂
    b[8] = f₀ - im*f₁
    b[9] = f₀ + im*f₁
    b[10] = f₀ + f₁ + f₂
    b[11] = f₀ + f₁ + im*f₂
    b[12] = f₀ + im*f₂ + f₃
    b[13] = f₀ + f₂ + f₃
    b[14] = f₀ - f₁

    ρ = zeros(12,12)
    for key in keys(a)
        c = kron(us.sq(a[key]),us.sq(b[key]))
        ρ = ρ + c
    end
    return ρ
end





end


# """# # Example 1: https://en.wikipedia.org/wiki/Quantum_entanglement"""
# function get_ρ_wiki(i)
#     ρ = Dict()
#     ϕ = [1, 0, 0, 1]
#     ρ[1] = ϕ*transpose(ϕ)
#
#     # = |00><00| + |02><02| + 2|11><11| + (|01> + |10>)(<01| + <10|)
#     ψ₀₁ = us.ψ(1,2,3)
#     ψ₁₀ = us.ψ(2,1,3)
#     ρ[2] = us.psep(1,1,3) + us.psep(1,3,3) + 2*us.psep(2,2,3) + us.sq(ψ₀₁ + ψ₁₀)
#     return ρ[i]
# end


# """
# ## Separability of Hermitian Tensors and PSD Decompositions
# ## Mareike Dressler, Jiawang Nie, Zi Yang
# ## https://arxiv.org/abs/2011.08132v1 Example 3.7
# # Hᵢ₁ᵢ₂ⱼ₁ⱼ₂ = i₁ + i₂ + j₁ + j₂
# """
# function get_ρ_DNY2()
#     H_flat = [i₁ + i₂ + j₁ + j₂ for i₁ in 1:2 for i₂ in 1:2 for j₁ in 1:2 for j₂ in 1:2 ]
#     return  reshape(H_flat,4,4)
# end
