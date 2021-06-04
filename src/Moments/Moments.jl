module Moments

using LinearAlgebra

srcDir = dirname(@__FILE__)
include(srcDir*"\\Utils.jl")
using .Utils

export make_mon_expo_mat,
       make_mom_expo_keys,
       get_ℝ_block_diag,
       make_xxᵀ_tens_yyᵀ,
       get_ℂ_block_diag,
       make_xx̄ᵀ_tens_yȳᵀ


"""output: exponents α ∈ Nⁿₜ of [x]≦ₜ of [x]₌ₜ (array of integers)"""
function make_mon_expo_arr(n::Int,t::Int, isLeq::Bool = true)
    if t < 0
        error("The ...")
    elseif t == 0 # [x]₌₀
        return zeros(Int32,1,n)
    else # [x]₌ₜ
        temp = make_mon_expo_arr(n,t-1,isLeq)
        e₁ = hcat(1, zeros(Int32,1,n-1))
        output = e₁ .+ temp
        for i = 1:n-1
            output = vcat(output,circshift(e₁,(0,i)) .+ temp)
        end

        if isLeq # [x]≦ₜ
            output = vcat(temp, output)
        end
        return unique(output,dims=1)
    end
end

"""output: exponents α ∈ Nⁿₜ of [x]≦ₜ of [x]₌ₜ (array of arrays of integers)"""
function make_mon_expo_vect(n::Int,t::Int, isLeq::Bool = true)
    mon_expo_arr = make_mon_expo_arr(n,t,isLeq)
    mon_expo     = [r  for r in  eachrow(mon_expo_arr)]
    return mon_expo
end

"""output: exponents α ∈ Nⁿₜ of [x]≦ₜ₁[x]≦ₜ₂ᵀ or [x]₌ₜ₁[x]₌ₜ₂ᵀ where , x = (x₁,x₂,...,xₙ)"""
function make_mon_expo_mat(n::Int,t,isLeq::Bool = true)
    mon1      = make_mon_expo_vect(n,t[1], isLeq)
    mon2      = make_mon_expo_vect(n,t[2], isLeq)

    mon1_mon2ᵀₜ_vec = [ mi+ mj for mi in mon1 for mj in mon2]
    mon1_mon2ᵀₜ    = reshape(mon1_mon2ᵀₜ_vec, (length(mon1), length(mon2)) )
    return mon1_mon2ᵀₜ
end

""" { (a,b) | xᵃx̄ᵇ ∈ [x]≦ₜ[x]≦ₜ* }"""
function make_mom_expo_keys(n::Int,t)
    mon_mat = make_mon_expo_mat(n,t)
    return unique( [ α for α in mon_mat ])
end


## Real
""" returns the exponents of xxᵀ⊗yyᵀ """
function make_xxᵀ_tens_yyᵀ(d)
    n = sum(d)
    pre_expo = make_mon_expo_mat(n,(1,0),false)
    x = pre_expo[1:d[1]]
    y = pre_expo[d[1]+1:n]

    xxᵀ_expo = x .+ reshape(x,1,d[1])
    yyᵀ_expo = y .+ reshape(y,1,d[2])

    Utils.var_kron(xxᵀ_expo,yyᵀ_expo)
end

"""Returns diagonal block of the momoment matrix with partition:
    xᵃ⁺ᶜyᵇ⁺ᵈ where |a+c|,|b+d| ∈ 2N """
function get_ℝ_block_diag(mom_matₜ_expo,d)
    mom_vecₜ_expo = mom_matₜ_expo[1,:]
    halfsum(arr) = [sum(arr[1:2*d[1]]),sum(arr[2*d[1]+1:end])]
    isoddd(p) = isodd(p[1]),isodd(p[2])

    isevev(pair) = pair == (false,false)
    isodod(pair) = pair == (true,true)
    isevod(pair) = pair == (false,true)
    isodev(pair) = pair == (true,false)

    select_first(p) = p[1]
    evev = select_first.(findall(isevev.(isoddd.(halfsum.(mom_vecₜ_expo)))))
    odod = select_first.(findall(isodod.(isoddd.(halfsum.(mom_vecₜ_expo)))))
    evod = select_first.(findall(isevod.(isoddd.(halfsum.(mom_vecₜ_expo)))))
    odev = select_first.(findall(isodev.(isoddd.(halfsum.(mom_vecₜ_expo)))))

    ℝ_block_diag = Dict("evev" => mom_matₜ_expo[evev,evev],
                         "odod" => mom_matₜ_expo[odod,odod],
                         "evod" => mom_matₜ_expo[evod,evod],
                         "odev" => mom_matₜ_expo[odev,odev])
    return ℝ_block_diag
end

## Complex
########### THIS IS TOO MUCH????
""""""
function get_ℂ_block_diag(mom_matₜ_expo,d)
    ℂ_block_diag = Dict("Default" => mom_matₜ_expo)
end


""" returns ℜ(xx*⊗yy*) and ℑ(xx*⊗yy*)  """
function make_xx̄ᵀ_tens_yȳᵀ(d)
    n = sum(2 .* d)
    pre_expo = make_mon_expo_mat(n,(1,0),false)
    uₓ = pre_expo[1:d[1]]
    vₓ = pre_expo[d[1]+1:2*d[1]]
    uᵥ = pre_expo[2*d[1]+1:2*d[1]+d[2]]
    vᵥ = pre_expo[2*d[1]+d[2]+1:2*d[1]+2*d[2]]

    U(vec,dₐ) = vec .+ reshape(vec,1,dₐ)
    W(vec1,vec2,dₐ) = vec1 .+ reshape(vec2,1,dₐ)

    Uₓ = U(uₓ,d[1])      # uₓuₓᵀ
    Vₓ = U(vₓ,d[1])      # vₓvₓᵀ

    Uᵥ = U(uᵥ,d[2])      #
    Vᵥ = U(vᵥ,d[2])      #

    Wₓᵀ = W(vₓ,uₓ,d[1])  # uₓvₓᵀ
    Wₓ  = W(uₓ,vₓ,d[1])  # vₓuₓᵀ

    Wᵥ  = W(vᵥ,uᵥ,d[2])  #
    Wᵥᵀ = W(uᵥ,vᵥ,d[2])  #

    Real_dict =  Dict(("+ Uₓ⊗Uᵥ")   => Utils.var_kron(Uₓ,Uᵥ),     # + Uₓ⊗Uᵥ
                      ("+ Uₓ⊗Vᵥ")   => Utils.var_kron(Uₓ,Vᵥ),     # + Uₓ⊗Vᵥ
                      ("+ Vₓ⊗Uᵥ")   => Utils.var_kron(Vₓ,Uᵥ),     # + Vₓ⊗Uᵥ
                      ("+ Vₓ⊗Vᵥ")   => Utils.var_kron(Vₓ,Vᵥ),     # + Vₓ⊗Vᵥ
                      ("- Wₓ⊗Wᵥ")   => Utils.var_kron(Wₓ,Wᵥ),     # - Wₓ⊗Wᵥ
                      ("+ Wₓ⊗Wᵥᵀ")  => Utils.var_kron(Wₓ,Wᵥᵀ),    # + Wₓ⊗Wᵥᵀ
                      ("+ Wₓᵀ⊗Wᵥ")  => Utils.var_kron(Wₓᵀ,Wᵥ),    # + Wₓᵀ⊗Wᵥ
                      ("- Wₓᵀ⊗Wᵥᵀ") => Utils.var_kron(Wₓᵀ,Wᵥᵀ))   # - Wₓᵀ⊗Wᵥᵀ

    Imag_dict =  Dict(("+ Uₓ⊗Wᵥ") => Utils.var_kron(Uₓ,Wᵥ),     # + Uₓ⊗Wᵥ
                      ("- Uₓ⊗Wᵥᵀ") => Utils.var_kron(Uₓ,Wᵥᵀ),   # - Uₓ⊗Wᵥᵀ
                      ("+ Vₓ⊗Wᵥ") => Utils.var_kron(Vₓ,Wᵥ),     # + Vₓ⊗Wᵥ
                      ("- Vₓ⊗Wᵥᵀ") => Utils.var_kron(Vₓ,Wᵥᵀ),   # - Vₓ⊗Wᵥᵀ
                      ("+ Wₓ⊗Uᵥ") => Utils.var_kron(Wₓ,Uᵥ),     # + Wₓ⊗Uᵥ
                      ("+ Wₓ⊗Vᵥ") => Utils.var_kron(Wₓ,Vᵥ),     # + Wₓ⊗Vᵥ
                      ("- Wₓᵀ⊗Uᵥ") => Utils.var_kron(Wₓᵀ,Uᵥ),   # - Wₓᵀ⊗Uᵥ
                      ("- Wₓᵀ⊗Vᵥ") => Utils.var_kron(Wₓᵀ,Vᵥ))   # - Wₓᵀ⊗Vᵥ

    return Dict("real" => Real_dict,
                "imag" => Imag_dict)
end




end
