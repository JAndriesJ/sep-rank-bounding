module Moments

export make_mon_expo,
       get_ℝ_block_diag,
       make_xxᵀ_tens_yyᵀ,
       get_ℂ_block_diag,
       make_xx̄ᵀ_tens_yȳᵀ,
       var_kron,
       var_kron_C,
       eᵢ,
       get_ℜℑααᶥββᶥᴿ,
       get_xx̄yȳMM_blocks,
       get_γγᶥζζᶥ_δδᶥηηᶥ



## Utils
"""The standard basis vector eᵢ in dimension n"""
eᵢ(n::Int,i::Int) = [Int(j==i) for j in 1:n]

""" input: A,B (arrays of integer tupples) output: A ⊗ B """
function var_kron(A,B)
   C_temp = [ a + b for  b in B, a ∈ A]
   rcs =  [(i,j) for i in 1:size(A)[1] for j in 1:size(B)[1]]
   return [C_temp[ij[2],kl[2],ij[1],kl[1]] for ij in rcs , kl in rcs ]
end

function var_kron_C(A,B)
   C_temp = [ [a] + b for  b in B, a ∈ A]
   rcs =  [(i,j) for i in 1:size(A)[1] for j in 1:size(B)[1]]
   return [C_temp[ij[2],kl[2],ij[1],kl[1]] for ij in rcs , kl in rcs ]
end

## Moments
"""[x]≦ₜ or [x]₌ₜ"""
function make_mon_expo(n::Int,t::Int; isle::Bool = true)
    @assert typeof(t) == Int64
    t == 0 ? (return [eᵢ(n,0)]) : 0
    tmp = make_mon_expo(n,t-1;isle=isle)
    M_vec = reshape([m + eᵢ(n,i) for i ∈ 1:n, m ∈ tmp],:,1)
    return unique(isle ? vcat(tmp,M_vec) : M_vec)
end

"""[x]≦ₜ[x]ᵀ≦ₜ or [x]₌ₜ[x]ᵀ₌ₜ"""
function make_mon_expo(n::Int,t::Tuple{Int,Int}; isle::Bool = true)
    M_vec1      = make_mon_expo(n,t[1]; isle=isle)
    M_vec2      = make_mon_expo(n,t[2]; isle=isle)
    M_mat_flat  = [ mi+ mj for mi in M_vec1 for mj in M_vec2]
    M_mat       = reshape(M_mat_flat, (length(M_vec1), length(M_vec2)) )
    return M_mat
end

## Complex moments
"""[x,̄x]≦ₜ or [x,̄x]ᵀ₌ₜ"""
make_mon_expo(d::Tuple{Int,Int},t::Int; isle::Bool = true) = make_mon_expo(sum(2 .* d),t; isle=isle)

"""[x,̄x]≦ₜ[x,̄x]*≦ₜ or [x,̄x]₌ₜ[x,̄x]*₌ₜ
but is actually [xᵣₑ,xᵢₘ]≦ₜ[xᵣₑ,xᵢₘ]ᵀ≦ₜ or [xᵣₑ,xᵢₘ]₌ₜ[xᵣₑ,xᵢₘ]ᵀ₌ₜ
"""
function make_mon_expo(d::Tuple{Int,Int},t::Tuple{Int,Int}; isle::Bool = true)
    M_vec1      = make_mon_expo(d,t[1]; isle=isle)
    M_vec2      = make_mon_expo(d,t[1]; isle=isle)
    M_mat      = [ mi+ mj for mi in M_vec1, mj in M_vec2]
    # M_mat       = reshape(M_mat_flat, (length(M_vec1), length(M_vec2)) )
    return M_mat
end
## Real
"""xxᵀ⊗yyᵀ"""
function make_xxᵀ_tens_yyᵀ(d)
    n = sum(d)
    pre_expo = make_mon_expo(n,1;isle = false)
    x, y = pre_expo[1:d[1]], pre_expo[d[1]+1:n]
    return var_kron(x .+ reshape(x,1,d[1]), y .+ reshape(y,1,d[2]))
end

"""Returns diagonal block of the momoment matrix with partition:
    xᵃ⁺ᶜyᵇ⁺ᵈ where |a+c|,|b+d| ∈ 2N """
function get_ℝ_block_diag(d::Tuple{Int,Int},t::Tuple{Int,Int})
    mom_mat = Moments.make_mon_expo(sum(d),t)
    # return Dict("Default" => mom_mat)
    mom_vec = mom_mat[1,:]
    halfsum(arr) = [sum(arr[1:d[1]]),sum(arr[(d[1]+1):end])]
    tf(arr)      = map(x -> (isodd(x[1]),isodd(x[2])), halfsum.(arr))
    tf2(arr)     = map(p -> p[1], findall(arr))

    isevev(pair) = pair == (false,false)
    isodod(pair) = pair == (true,true)
    isevod(pair) = pair == (false,true)
    isodev(pair) = pair == (true,false)

    ee = tf2(isevev.(tf(mom_vec)))
    oo = tf2(isodod.(tf(mom_vec)))
    eo = tf2(isevod.(tf(mom_vec)))
    oe = tf2(isodev.(tf(mom_vec)))

    ℝ_block_diag = Dict("ee" => mom_mat[ee,ee],
                        "oo" => mom_mat[oo,oo],
                        "eo" => mom_mat[eo,eo],
                        "oe" => mom_mat[oe,oe])

    for key in ["ee","oo","eo","oe"]
         isempty(ℝ_block_diag[key]) ? delete!(ℝ_block_diag,key) : continue
    end
    return ℝ_block_diag
end

## Complex

""" ℜe(xx*⊗yy*), ℑm(xx*⊗yy*)
xx*⊗yy* = (xᵣₑ + i xᵢₘ)(xᵣₑ - i xᵢₘ)ᵀ ⊗ (yᵣₑ + i yᵢₘ)(yᵣₑ - i yᵢₘ)ᵀ
= (xᵣₑxᵣₑᵀ + xᵢₘxᵢₘᵀ + i (xᵢₘxᵣₑᵀ - xᵣₑxᵢₘᵀ)) ⊗ (yᵣₑyᵣₑᵀ + yᵢₘyᵢₘᵀ + i (yᵢₘyᵣₑᵀ - yᵣₑyᵢₘᵀ))
⟹ ℜe(xx*⊗yy*) = (xᵣₑxᵣₑᵀ + xᵢₘxᵢₘᵀ ) ⊗ (yᵣₑyᵣₑᵀ + yᵢₘyᵢₘᵀ ) - (xᵢₘxᵣₑᵀ - xᵣₑxᵢₘᵀ) ⊗ (yᵢₘyᵣₑᵀ - yᵣₑyᵢₘᵀ)
⟹ ℑm(xx*⊗yy*) = (xᵣₑxᵣₑᵀ + xᵢₘxᵢₘᵀ) ⊗ (yᵢₘyᵣₑᵀ - yᵣₑyᵢₘᵀ) + (xᵢₘxᵣₑᵀ - xᵣₑxᵢₘᵀ)⊗(yᵣₑyᵣₑᵀ + yᵢₘyᵢₘᵀ) """
U(vec,dₐ) = vec .+ reshape(vec,1,dₐ)
W(vec1,vec2,dₐ) = vec1 .+ reshape(vec2,1,dₐ)
function make_xx̄ᵀ_tens_yȳᵀ(d)
    d₁,d₂ = d
    n = sum(2 .* d)
    pre_expo = make_mon_expo(n,1;isle =false)
    xᵣₑ, xᵢₘ = pre_expo[1:d₁], pre_expo[d₁+1:2*d₁]
    yᵣₑ, yᵢₘ = pre_expo[2*d₁+1:2*d₁+d₂], pre_expo[2*d₁+d₂+1:2*d₁+2*d₂]

    xᵣₑxᵣₑᵀ, xᵢₘxᵢₘᵀ, yᵣₑyᵣₑᵀ, yᵢₘyᵢₘᵀ = U(xᵣₑ,d₁), U(xᵢₘ,d₁), U(yᵣₑ,d₂), U(yᵢₘ,d₂)
    xᵢₘxᵣₑᵀ, xᵣₑxᵢₘᵀ, yᵢₘyᵣₑᵀ, yᵣₑyᵢₘᵀ = W(xᵢₘ,xᵣₑ,d₁), W(xᵣₑ,xᵢₘ,d₁), W(yᵢₘ,yᵣₑ,d₂), W(yᵣₑ,yᵢₘ,d₂)

    Real_dict =  Dict(("+ 1")  => var_kron(xᵣₑxᵣₑᵀ,yᵣₑyᵣₑᵀ),
                      ("+ 2")  => var_kron(xᵣₑxᵣₑᵀ,yᵢₘyᵢₘᵀ),
                      ("+ 3")  => var_kron(xᵢₘxᵢₘᵀ,yᵣₑyᵣₑᵀ),
                      ("+ 4")  => var_kron(xᵢₘxᵢₘᵀ,yᵢₘyᵢₘᵀ),
                      ("- 5")  => var_kron(xᵢₘxᵣₑᵀ,yᵢₘyᵣₑᵀ),
                      ("+ 6")  => var_kron(xᵢₘxᵣₑᵀ,yᵣₑyᵢₘᵀ),
                      ("+ 7")  => var_kron(xᵣₑxᵢₘᵀ,yᵢₘyᵣₑᵀ),
                      ("- 8")  => var_kron(xᵣₑxᵢₘᵀ,yᵣₑyᵢₘᵀ))

    Imag_dict =  Dict(("+ 1")  => var_kron(xᵣₑxᵣₑᵀ,yᵢₘyᵣₑᵀ),
                      ("- 2") => var_kron(xᵣₑxᵣₑᵀ,yᵣₑyᵢₘᵀ),
                      ("+ 3")  => var_kron(xᵢₘxᵢₘᵀ,yᵢₘyᵣₑᵀ),
                      ("- 4") => var_kron(xᵢₘxᵢₘᵀ,yᵣₑyᵢₘᵀ),
                      ("+ 5")  => var_kron(xᵢₘxᵣₑᵀ,yᵣₑyᵣₑᵀ),
                      ("+ 6")  => var_kron(xᵢₘxᵣₑᵀ,yᵢₘyᵢₘᵀ),
                      ("- 7") => var_kron(xᵣₑxᵢₘᵀ,yᵣₑyᵣₑᵀ),
                      ("- 8") => var_kron(xᵣₑxᵢₘᵀ,yᵢₘyᵢₘᵀ))

    return Dict("real" => Real_dict,
                "imag" => Imag_dict)
end




## Block diagonalization
# The complex blocks
"""Splits a vector of exponents into the register components ααᶥββᶥ →  (α,αᶥ,β,βᶥ)"""
split_expo(ααᶥββᶥ,d) =     (ααᶥββᶥ[1:d[1]],
                            ααᶥββᶥ[d[1]+1:2*d[1]],
                            ααᶥββᶥ[1+2*d[1]:2*d[1]+d[2]],
                            ααᶥββᶥ[1+2*d[1]+d[2]:end])

"""Gets the principal block of the moment matrix after block diagonalization"""
function get_xx̄yȳMM_blocks(d,t)
    Iᵗ = Moments.make_mon_expo(d,t[1])
    Iᵗ_temp = map(x ->  sum.(split_expo(x,d)),Iᵗ)
    r_list  = map(v -> v[1]-v[2],Iᵗ_temp) ; s_list  = map(v -> v[3]-v[4],Iᵗ_temp)

    Iᵗ_d = Dict()
    for r ∈ -t[1]:t[1], s ∈ -t[1]:t[1]
        T = Iᵗ[(r_list .== r) .& (s_list .== s) ]
        isempty(T) ? continue : Iᵗ_d[r,s] = T
    end
    return Dict(zip(keys(Iᵗ_d),[Iᵗ_d[k] .+ reshape(Iᵗ_d[k],1,:) for k in keys(Iᵗ_d)]))
end

## The real analogue
"""Returns all  (γ,γᶥ,ζ,ζᶥ),(δ,δᶥ,η,ηᶥ) ∈ (ℕᵈ)⁴ s.t.
(γ,γᶥ,ζ,ζᶥ)+(δ,δᶥ,η,ηᶥ)=(α,αᶥ,β,βᶥ)  ∀  ααᶥββᶥ ∈ Moment matrixₜ
"""
function get_γγᶥζζᶥ_δδᶥηηᶥ(d,t)
    M = Moments.make_mon_expo(d,t)
    M_vec = M[:,1]
    γγᶥζζᶥ_δδᶥηηᶥ_dict = Dict()
    for k in Moments.make_mon_expo(d,2*t[1])# Define variables in the moment matrix.
        fd = hcat(map(x -> [x[1],x[2]],findall(M .== [k]))...)
        γγᶥζζᶥ_δδᶥηηᶥ_dict[k] = [M_vec[fd[1,:]], M_vec[fd[2,:]]]
    end
    return γγᶥζζᶥ_δδᶥηηᶥ_dict
end

"""The coefficients function"""
pf(x) = prod(factorial.(x))
function get_coef(part1,part2,d)
    γ,γᶥ,ζ,ζᶥ = split_expo(part1,d) ; δ,δᶥ,η,ηᶥ = split_expo(part2,d)
    a = ((-1.0)^sum(δᶥ+ηᶥ))*((1.0im)^sum(δ+δᶥ+η+ηᶥ))*prod(pf.([[γ,γᶥ,ζ,ζᶥ]+[δ,δᶥ,η,ηᶥ]...]))
    return a/prod(pf.([γ,γᶥ,ζ,ζᶥ,δ,δᶥ,η,ηᶥ]))
end
get_coef(γ_δ_pair,d) = map((x,y)->get_coef(x,y,d),γ_δ_pair...)


""" get_expo(γγᶥδδᶥ,ζζᶥηηᶥ,d) = [γ+γᶥ,δ+δᶥ,η+ηᶥ,ζ+ζᶥ] """
function get_expo(γ_,δ_,d)
    γ,γᶥ,ζ,ζᶥ = split_expo(γ_,d) ; δ,δᶥ,η,ηᶥ = split_expo(δ_,d)
    return vcat(γ+γᶥ,δ+δᶥ,ζ+ζᶥ,η+ηᶥ)
end
get_expo(γ_δ_,d) = map((x,y)->get_expo(x,y,d),γ_δ_...)


"""f"""
function get_ℜℑααᶥββᶥᴿ(d,ααᶥββᶥ::Array{Int64,1},γdict)
    γ_δ_ = γdict[ααᶥββᶥ]
    return   get_expo(γ_δ_,d),get_coef(γ_δ_,d)
end

function get_ℜℑααᶥββᶥᴿ(d,B,γdict)
    MMex = map(x->get_ℜℑααᶥββᶥᴿ(d,x,γdict)[1],B)
    MMCoef = map(x->get_ℜℑααᶥββᶥᴿ(d,x,γdict)[2],B)

    real_mask = isreal.(MMCoef)

    RMMex = copy(MMex) ; IMMex = copy(MMex)
    RMMex[.!real_mask] .= [[]]
    IMMex[real_mask] .= [[]]

    RMMCoef = copy(MMCoef) ; IMMCoef = copy(MMCoef)
    RMMCoef[.!real_mask] .= [[0.0]] ;RMMCoef = real.(RMMCoef)
    IMMCoef[real_mask] .= [[0.0]] ; IMMCoef = imag.(IMMCoef)

    MMexᴿ = hcat(vcat(RMMex,IMMex),vcat(IMMex,RMMex))
    MMCoefᴿ = hcat(vcat(RMMCoef,IMMCoef),vcat((-1.0)*IMMCoef,RMMCoef))
    return MMexᴿ,MMCoefᴿ
end

""""""
function get_ℂ_block_diag(d,t;noBlock = false)
    if noBlock
        MMexᴿ = Dict("Default" => map(x->[x],make_mon_expo(d,t)))
        MMCoefᴿ = Dict("Default" => fill([1.0],size(MMexᴿ["Default"])...))
        return MMexᴿ,MMCoefᴿ
    end
    xx̄yȳMM_blocks = get_xx̄yȳMM_blocks(d,t)
    γγᶥζζᶥ_δδᶥηηᶥ = get_γγᶥζζᶥ_δδᶥηηᶥ(d,t)

    MMexᴿ= Dict() ; MMCoefᴿ = Dict()
    for B in keys(xx̄yȳMM_blocks)
        MMexᴿ[B], MMCoefᴿ[B] = get_ℜℑααᶥββᶥᴿ(d,xx̄yȳMM_blocks[B],γγᶥζζᶥ_δδᶥηηᶥ)
    end
    return MMexᴿ,MMCoefᴿ
end





end
