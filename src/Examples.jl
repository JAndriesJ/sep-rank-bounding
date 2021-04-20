module Examples
using LinearAlgebra
using Random

srcDir = dirname(@__FILE__)*"\\"
include(srcDir *"Utils.jl")
using .Utils


export get_examples,
       partial_transpose_per,
       gen_rand_state,
       isPPT,
       get_ρ_CD1,
       get_ρ_CD2

"""
Returns dict with three keys: ["rand", "sep", "ent"]
each entry is also a dict.
"rand": ["randd3r3"  "randd3r6"  "randd2r2"  "randd3r5"  "randd3r1"  "randd2r1"  "randd2r6"  "randd3r7"  "randd2r7"  "randd2r5"  "randd2r4"  "randd3r2"  "randd3r4"  "randd2r3"]
"sep":  ["sep4d3r4",  "sep5d3r5",  "sep3d3r3",  "sep1d2r2",  "sep2d2r3",  "sep6d3r2"]
"ent":  ["ent4d3r4",  "ent2d2r2",  "ent3d2r2",  "ent1d2r1",  "ent5d2r2"]
"""
function get_examples()
   ρ         = Dict()
   ρ["rand"] = get_rand_states([4],5:20)
   ρ["sep"]  = get_sep_example()
   ρ["ent"]  = get_ent_example()

   ρ["HK"]  = get_ρ_HK()
   ρ_temp = Dict()
   for i in 1:10
       ρ_temp[i]  = get_ρ_BP(rand(1:4,8,1))
   end
   ρ["BP"]  =  ρ_temp
   ρ["CD1"] = Dict("CD13" => get_ρ_CD1(3),
                   "CD14" => get_ρ_CD1(4))
   ρ["CD2"] = get_ρ_CD2()
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

function isPPT(x)
    d = Int(sqrt(size(x)[1]))
    dims = [d,d]
    is_1_PPT = eigvals(Examples.partial_transpose_per(x,1,[d,d]))[1] > -1.e-16
    is_2_PPT = eigvals(Examples.partial_transpose_per(x,2,[d,d]))[1] > -1.e-16
    return is_1_PPT && is_2_PPT
end

"""
Makes the trace equal to one via scaling
"""
maketraceone(ρ) = ρ/tr(ρ)



## Random matrices
"""
generates a matrix of the form ∑ʳaᵀa⊗bᵀb, with a,b ∈ uniform random entry-wise.  d ∈ {2,3,4} , r ∈ [9]
"""
function gen_rand_state(d::Integer, r::Integer)
    seed = 343
    Random.seed!(seed)
    a = [randn(1,d) for i in 1:r ]
    b = [randn(1,d) for i in 1:r ]
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
            rankρ = rank(ρᵈʳ_temp)
            ρᵈʳ["rand($d,$r)d$d"*"r$rankρ"] = ρᵈʳ_temp
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
##  Qubit-qudit states with positive partial transpose
##  Lin Chen, Dragomir Z. Djokovic
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

## Separability of Hermitian Tensors and PSD Decompositions
## Mareike Dressler, Jiawang Nie, Zi Yang
##  https://arxiv.org/abs/2011.08132v1 Example 3.8
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

    # ent 4
    # = |00><00| + |02><02| + 2|11><11| + (|01> + |10>)(<01| + <10|)
    # C²⊗C³
    ψ₀₁ = ψ(1,2,3)
    ψ₁₀ = ψ(2,1,3)
    ρ_temp = psep(1,1,3) + psep(1,3,3) + 2*psep(2,2,3) + sq(ψ₀₁ + ψ₁₀)
    ρ["ent4d3r4"] = ρ_temp

##    Separability of Mixed States: Necessary and Sufficient Conditions
##    Michal Horodecki, Pawel Horodecki, Ryszard Horodecki
## http://arxiv.org/abs/quant-ph/9605038v2  page 9 eq. (22)
    Random.seed!(343)
    a  = 0.5; b = 0.8; # arbitary possitive numbers
    p  = rand()
    ψ₁ = a*ψ(2,2,2) + b*ψ(1,1,2)
    ψ₂ = a*ψ(2,1,2) + b*ψ(1,2,2)
    ρ_temp = p*sq(ψ₁) + (1 - p)*sq(ψ₂)
    ρ["ent2d2r2"] = ρ_temp


##   page 10 eq. (25) fails PPT
    a =  1/sqrt(2) ; p = rand();
    ψₐ = a*(ψ(1,2,2) - ψ(2,1,2))
    ρ_temp = p*sq(ψₐ) + (1 - p)*psep(1,1,2)
    ρ["ent3d2r2"] = ρ_temp



## Separability of Hermitian Tensors and PSD Decompositions
## Mareike Dressler, Jiawang Nie, Zi Yang
## https://arxiv.org/abs/2011.08132v1 Example 3.7
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


## 1999, Bruß,Peres,Construction of quantum states with bound entanglement.pdf
## https://arxiv.org/abs/quant-ph/9911056v1
function get_ρ_BP(V)
    m,s,n,a,b,c,t,d = V
    # m,s,n,a,b,c,t,d = [1,1,1,1,1,1,1,1]
    V₁ =  [m  0  s
           0  n  0
           0  0  0]

    V₂ =   [0  a  0
            b  0  c
            0  0  0]

    V₃ =   [conj(n)  0   0
            0  -conj(m)  0
            t      0     0]
    V₄ = [0    conj(b)  0
      -conj(a)   0      0
          0      d      0]
    return kron(V₁,V₁)  + kron(V₂,V₂)  + kron(V₃,V₃)  + kron(V₄,V₄)
end

## 2013,Chen,Dokovic,Dimensions, lengths and separability in finite-dimensional quantum systems
function get_ρ_CD1(N)
    n₁ = 2
    e₀  = Utils.eᵢ(n₁,1)
    e₁  = Utils.eᵢ(n₁,2)

    f = Dict()
    for i ∈ 1:N
        f[i-1] = Utils.eᵢ(N,i)
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
        c = kron(sq(aₖ[key]),sq(bₖ[key]))
        ρ = ρ + c
    end
    return ρ
end

function get_ρ_CD2()
    n₁ = 3
    n₂ = 4
    e₀  = Utils.eᵢ(n₁,1)
    e₁  = Utils.eᵢ(n₁,2)
    e₂  = Utils.eᵢ(n₁,3)

    f₀  = Utils.eᵢ(n₂,1)
    f₁  = Utils.eᵢ(n₂,2)
    f₂  = Utils.eᵢ(n₂,3)
    f₃  = Utils.eᵢ(n₂,4)

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
        c = kron(sq(a[key]),sq(b[key]))
        ρ = ρ + c
    end
    return ρ
end

##Separable states with unique decompositions,
## Kil-Chan Ha, Seung-Hyeok Kye
## https://arxiv.org/abs/1210.1088v5
function  get_ρ_HK(b_vec = [0.3220677161622946, 2.200060128843455, 0.029041219347510496, -1.049106738126062,  1.7134943546489485, 2, 3, 1])
    ρ = Dict()
    for i in 1:length(b_vec)
        b = b_vec[i]
        ρ_temp = gen_Ha_Kye_ex(b)
        r      = rank(ρ_temp)
        ρ["HK$i"*"d3"*"r$r"] =  ρ_temp
    end
    return ρ
end


function gen_Ha_Kye_ex(b)
    # @assert b == 1.0
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



end
