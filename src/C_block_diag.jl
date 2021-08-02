module C_block_diag

include(dirname(@__FILE__)*"\\Moments\\Moments.jl")
using .Moments

export  split_expo,
        get_all_partitions,
        get_all_blocks,
        get_coef,
        get_γγᶥδδᶥζζᶥηηᶥ,
        get_ℜℑααᶥββᶥᴿ,
        get_real_blocks

## The complex blocks
"""Splits a vector of exponents into the register components ααᶥββᶥ →  (α,αᶥ,β,βᶥ)"""
split_expo(ααᶥββᶥ,d) =     (ααᶥββᶥ[1:d[1]],
                            ααᶥββᶥ[d[1]+1:2*d[1]],
                            ααᶥββᶥ[1+2*d[1]:2*d[1]+d[2]],
                            ααᶥββᶥ[1+2*d[1]+d[2]:end])

#@test split_expo([1,1,1,2,2,2,3,3,3,4,4,4],(3,3)) == ([1, 1, 1], [2, 2, 2], [3, 3, 3], [4, 4, 4])

"""Partitions Iᵗ such that M( Iᵗ , Iᵗ ) would be block diagonal"""
get_r(v) = v[1] - v[2]
get_s(v) = v[3] - v[4]
function get_all_partitions(d,t)
    Iᵗ = Moments.make_mon_expo(d,t[1])
    Iᵗ_temp = map(x ->  sum.(split_expo(x,d)),Iᵗ)

    r_list = map(v -> get_r(v),Iᵗ_temp)
    s_list = map(v -> get_s(v),Iᵗ_temp)

    Iᵗ_dict = Dict()
    Iᵗ_perm = []
    for r ∈ -t[1]:t[1], s ∈ -t[1]:t[1]
        T = Iᵗ[(r_list .== r) .& (s_list .== s) ]
        isempty(T) ? continue : Iᵗ_dict[r,s] = T
    end
    return Iᵗ_dict
end


"""Gets the principal block of the moment matrix after block diagonalization"""
get_all_blocks(Iᵗd) = Dict(zip(keys(Iᵗd),[Iᵗd[k] .+ reshape(Iᵗd[k],1,:) for k in keys(Iᵗd)]))



## The real analogue

"""Returns all  γγᶥδδᶥ,ζζᶥηηᶥ ∈ (ℕᵈ)⁴ s.t.
ααᶥββᶥ =  γγᶥδδᶥ,ζζᶥηηᶥ
∀  ααᶥββᶥ ∈ Moment matrixₜ
"""
function get_γγᶥδδᶥζζᶥηηᶥ(d,t)
    M = Moments.make_mon_expo(d,t)
    M_vec = M[:,1]
    lok = Moments.make_mon_expo(d,2*t[1]) # Define variables in the moment matrix.
    nara = Dict()
    for k in lok
        fd = hcat(map(x -> [x[1],x[2]],findall(M .== [k]))...)
        nara[k] = [M_vec[fd[1,:]], M_vec[fd[2,:]]]
    end
    return nara
end

"""The coefficients function"""
pf(x) = prod(factorial.(x))
function get_coef(γγᶥδδᶥ,ηηᶥζζᶥ,d)
    γ,γᶥ,δ,δᶥ = split_expo(γγᶥδδᶥ,d)
    η,ηᶥ,ζ,ζᶥ = split_expo(ηηᶥζζᶥ,d)
    a = ((-1)^sum(δᶥ + ζᶥ) * (im)^sum(δ+δᶥ+ζ+ζᶥ))* prod(pf.([γ+δ γᶥ+δᶥ η+ζ ηᶥ+ζᶥ]))
    return a/prod(pf.([γ,γᶥ,δ,δᶥ,η,ηᶥ,ζ,ζᶥ]))
end

"""f"""
function get_ℜℑααᶥββᶥᴿ(d,ααᶥββᶥ,γdict,K)
    Tar = γdict[ααᶥββᶥ]
    Res_vec = []
    Res_coef = []
    for k ∈ 1:K
        if k < length(Tar[1])
            γ_,η_ = Tar[1][k],Tar[2][k]
            push!(Res_coef,get_coef(γ_,η_,d))

            γ,γᶥ,δ,δᶥ = split_expo(γ_,d)
            η,ηᶥ,ζ,ζᶥ = split_expo(η_,d)
            push!(Res_vec,vcat(γ+γᶥ,δ+δᶥ,η+ηᶥ,ζ+ζᶥ))
        else
            push!(Res_coef,0.0)
            push!(Res_vec,Int.(zeros(length(ααᶥββᶥ))))
        end
    end
    return  Res_vec,Res_coef
end

"""Get the exponents and coefficients in matrix form"""
function get_real_block(d,B,γdic,K)
    Res_vec = Dict()
    Res_coef = Dict()
    for i in 1:size(B)[1], j in i:size(B)[1]
            Res_vec[i,j], Res_coef[i,j] = get_ℜℑααᶥββᶥᴿ(d,B[i,j],γdic,K)
            Res_vec[j,i]  = Res_vec[i,j]
            Res_coef[j,i] = Res_coef[i,j]
    end

    Expo_mat_dict = Dict()
    Coef_mat_dict = Dict()
    for k in 1:K
        Expo_mat_dict[k] = [Res_vec[i,j][k] for i = 1:size(B)[1], j = 1:size(B)[1]]
        Coef_mat_dict[k] = [Res_coef[i,j][k] for i = 1:size(B)[1], j = 1:size(B)[1]]
    end
    return clean_real_blocks(Coef_mat_dict, Expo_mat_dict)
end

function clean_real_blocks(Coef_mat_dict, Expo_mat_dict)
    for k in keys(Coef_mat_dict)
        if iszero(Coef_mat_dict[k]) ## Get rid of zeros
            delete!(Coef_mat_dict,k)
            delete!(Expo_mat_dict,k)
            continue
        end

        if isreal(Coef_mat_dict[k]) ## Reduce to real
            Coef_mat_dict[k] = real(Coef_mat_dict[k])
        # elseif imag(Coef_mat_dict[k]) == Coef_mat_dict[k]
        #     Coef_mat_dict[k] = imag(Coef_mat_dict[k])
        end
    end
    return Coef_mat_dict, Expo_mat_dict
end

function get_real_blocks(d,t)
    Iᵗ_dict = C_block_diag.get_all_partitions(d,t)
    block_dict = C_block_diag.get_all_blocks(Iᵗ_dict)
    γγᶥδδᶥζζᶥηηᶥ_dict = C_block_diag.get_γγᶥδδᶥζζᶥηηᶥ(d,t)
    K = maximum([length(x[1]) for x in values(γγᶥδδᶥζζᶥηηᶥ_dict)])
    Coef_mat_dict = Dict()
    Expo_mat_dict = Dict()
    for k in keys(block_dict)
        Coef_mat_dict[k], Expo_mat_dict[k] = C_block_diag.get_real_block(d,block_dict[k],γγᶥδδᶥζζᶥηηᶥ_dict,K)
    end
    return Coef_mat_dict, Expo_mat_dict
end



end
