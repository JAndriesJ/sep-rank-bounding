## Step 1: Tell Julia in which folder you are:
    sep_rank_proj_path = dirname(dirname(@__FILE__))
    srcDir  = sep_rank_proj_path*"\\src\\";
    testDir = sep_rank_proj_path*"\\test\\";

## Examples
include(srcDir*"Examples\\Examples.jl")
using .Examples

examples     = Examples.get_examples()
examples_all = merge(examples["ent"],examples["sep"])
# Get information on the examples in a tabular form, used the DataFrames package.
df           = Examples.get_example_overview(false)
show(df)
# Loading a specific example
ρ = examples_all["Eq25HHH96"]
# Routine test to see if the examples are still doing what I think they should.
# include(testDir*"testExamples.jl")
## Moments
include(srcDir*"Moments\\Moments.jl")
using .Moments
# Routine test to see if the moments code is still doing what I think they should.
# look at the file to see what test are conducted.
# include(testDir*"testMoments.jl")
##  Constraints
include(srcDir*"Constraints\\Utils_cons.jl")
# include(srcDir*"Constraints\\R_constraints.jl")
include(srcDir*"Constraints\\C_constraints.jl")
using .Utils_cons
# using .R_constraints
using .C_constraints
const uc = Utils_cons

### CHECK IF THIS IS ALL WHAT IT SHOULD BE!!!
using JuMP
#
t = (2,2)
d = (2,2)

model = Model()
mom_list = C_constraints.make_mon_expo_keys(d,t[1]) # Define variables in the moment matrix.
@variable(model, Lx[mom_list] ) #

# PSD_MM_con  = R_constraints.make_PSD_MM_con(d,t,Lx)
# ord4_con    = R_constraints.make_ord4_con(d,Lx)
L_xx̄ᵀ_tens_yȳᵀ  = Utils_cons.idx2varxx̄ᵀtyȳᵀ(Lx,Moments.make_xx̄ᵀ_tens_yȳᵀ(d))



# @assert ord4_con == [Lx[[2, 0, 2, 0]]  Lx[[2, 0, 1, 1]]  Lx[[1, 1, 2, 0]]  Lx[[1, 1, 1, 1]]
#                      Lx[[2, 0, 1, 1]]  Lx[[2, 0, 0, 2]]  Lx[[1, 1, 1, 1]]  Lx[[1, 1, 0, 2]]
#                      Lx[[1, 1, 2, 0]]  Lx[[1, 1, 1, 1]]  Lx[[0, 2, 2, 0]]  Lx[[0, 2, 1, 1]]
#                      Lx[[1, 1, 1, 1]]  Lx[[1, 1, 0, 2]]  Lx[[0, 2, 1, 1]]  Lx[[0, 2, 0, 2]]]
#
# loc_cons_S₁     = R_constraints.make_loc_cons_S₁(ρ,d,t,Lx)
# loc_cons_S₂     = R_constraints.make_loc_cons_S₂(ρ,d,t,Lx)
# loc_cons_S₃, g₂ = R_constraints.make_loc_cons_S₃(ρ,d,t,Lx)
# G_con           = R_constraints.make_G_con(ρ,d,t,Lx)

## Model
# include(srcDir*"Model\\Utils_Model.jl")
include(srcDir*"Model\\R_sep_Model.jl")
# include(srcDir*"Model\\C_sep_Model.jl")
using .R_sep_Model
# using .C_sep_Model
# using .Utils_Model
ρ = examples_all["T1CD12i"]
d = (2,2)
t = (2,2)
Lx,sep_mod = R_sep_Model.Modelξₜˢᵉᵖ(ρ,d,t)

## Compute
include(srcDir*"sep_Compute.jl")
using .sep_Compute
# using JuMP
# using MosekTools
sep_Compute.Computeξₜˢᵉᵖ(Lx,sep_mod)
# using JuMP
# V = value.(Lx).data
# Lx.data
# Lx_dict = Dict(zip(Lx.data,V))
# value.(ord4_con[1,1])





## Batch model, sep_Compute and post process.
if true
    t = (2,2)
    save_dir = sep_rank_proj_path*"\\PreRunModels\\C-models\\t=$(t[1])\\"
    Utils_Model.batch_model(t,examples_all,df,save_dir,["S₁sG","S₂sG","S₃sG"],false)
    sep_Compute.batch_Computeξₜˢᵉᵖ(save_dir)

    ## unstacking the constraints
    include(srcDir *"PostProcTools.jl")
    using .PostProcTools
    PostProcTools.unstack_constraints(save_dir,df)
end

## individual case.
# ex = "p12HK14b"
ex = "T1CD12i"
t = (2,2)
d = (2,2)
LxR,opt_modR = sep_Compute.quick_run_spec_example(ex,t,"S₁sG",true)
LxC,opt_modC = sep_Compute.quick_run_spec_example(ex,t,"S₁sG",false)
## sadf
d = (3,3)
t = (2,2)
n = sum(2 .*d)
m1 = Moments.make_mon_expo_vect(n,2, true)
m[90]
C_block_diag.split_expo(m[90],d)

m1 = Moments.make_mon_expo_mat(n,t,true)
m2 = Moments.make_mon_expo_mat(d,t,true)
m1 == m2

## function get_all_blocks()
# n = sum(2 .* d)
# model = JuMP.Model()
# ## Create variables
# list_of_keys = C_constraints.make_mom_expo_keys(n,t) # Define variables in the moment matrix.
# @variable(model, Lx[list_of_keys] ) #
#      PSD Moment matrix blocks
# Coef_mat_dict, Expo_mat_dict = C_block_diag.get_real_blocks(d,t)
# PSD_MM_con1 = Dict()
# for block in keys(Coef_mat_dict)
#     PSD_MM_con1[block] = Utils_cons.index_to_var(Lx, Expo_mat_dict[block],Coef_mat_dict[block])
# end
# sum([size(PSD_MM_con1[block])[1] for b in keys(PSD_MM_con1)])
#

## Test C_block_diag.jl
include(srcDir *"C_block_diag.jl")
using .C_block_diag

@testset "C_block_diag.get_all_partitions" begin
    t = (2,2)
    d = (2,2)
    Iᵗ = Moments.make_mon_expo_mat(sum(2 .*d),(0,t[1]),true)
    Iᵗ_part = C_block_diag.get_all_partitions(d,t)
    @test sum([length(Iᵗ_part[b]) for b in keys(Iᵗ_part)]) == length(Iᵗ)
    get_rs(v) = v[1] - v[2],v[3] - v[4]
    kit = [ unique(map(v->get_rs(v),(map(x->sum.(C_block_diag.split_expo(x,d)), Iᵗ_part[b]))))[1] for b in keys(Iᵗ_part)]
    @test kit == [keys(Iᵗ_part)...]


    t = (3,3)
    d = (2,2)
    Iᵗ = Moments.make_mon_expo_mat(sum(2 .*d),(0,t[1]),true)
    Iᵗ_part = C_block_diag.get_all_partitions(d,t)
    @test sum([length(Iᵗ_part[b]) for b in keys(Iᵗ_part)]) == length(Iᵗ)
    kit = [ unique(map(v->get_rs(v),(map(x->sum.(C_block_diag.split_expo(x,d)), Iᵗ_part[b]))))[1] for b in keys(Iᵗ_part)]
    @test kit == [keys(Iᵗ_part)...]

    t = (2,2)
    d = (3,3)
    Iᵗ = Moments.make_mon_expo_mat(sum(2 .*d),(0,t[1]),true)
    Iᵗ_part = C_block_diag.get_all_partitions(d,t)
    @test sum([length(Iᵗ_part[b]) for b in keys(Iᵗ_part)]) == length(Iᵗ)
    kit = [ unique(map(v->get_rs(v),(map(x->sum.(C_block_diag.split_expo(x,d)), Iᵗ_part[b]))))[1] for b in keys(Iᵗ_part)]
    @test kit == [keys(Iᵗ_part)...]
end

@testset "C_block_diag.get_all_blocks" begin
    t_i,d_i = rand(2:4,2)
    t = (t_i,t_i)
    d = (d_i ,d_i)
    Iᵗ_part = C_block_diag.get_all_partitions(d,t)
    Blocks = C_block_diag.get_all_blocks(Iᵗ_part)

    Iᵗ = Moments.make_mon_expo_mat(sum(2 .*d),(0,t[1]),true)
    @test sum([size(Blocks[b])[1] for b in keys(Blocks)]) == length(Iᵗ)
end

@testset "C_block_diag.get_γγᶥδδᶥζζᶥηηᶥ" begin
    t_i,d_i = rand(2:4,2)
    t = (t_i,t_i)
    d = (d_i ,d_i)
    γγᶥδδᶥζζᶥηηᶥ_dict = C_block_diag.get_γγᶥδδᶥζζᶥηηᶥ(d,t)
    @test length(unique(keys(γγᶥδδᶥζζᶥηηᶥ_dict))) == length(keys(γγᶥδδᶥζζᶥηηᶥ_dict))
    @test [ keys(γγᶥδδᶥζζᶥηηᶥ_dict)...] == [ unique(γγᶥδδᶥζζᶥηηᶥ_dict[k][1] + γγᶥδδᶥζζᶥηηᶥ_dict[k][2])[1]   for  k in keys(γγᶥδδᶥζζᶥηηᶥ_dict)]
end

@testset "C_block_diag.get_γγᶥδδᶥζζᶥηηᶥ" begin
    t_i,d_i = rand(2:3,2)
    t = (t_i,t_i)
    d = (d_i ,d_i)
    γγᶥδδᶥζζᶥηηᶥ_dict = C_block_diag.get_γγᶥδδᶥζζᶥηηᶥ(d,t)
    pf(x) = prod(factorial.(x))
    for k in keys(γγᶥδδᶥζζᶥηηᶥ_dict)
        γγᶥδδᶥ = γγᶥδδᶥζζᶥηηᶥ_dict[k][1][1]
        ζζᶥηηᶥ = γγᶥδδᶥζζᶥηηᶥ_dict[k][2][1]

        γ,γᶥ,δ,δᶥ = C_block_diag.split_expo(γγᶥδδᶥ,d)
        η,ηᶥ,ζ,ζᶥ = C_block_diag.split_expo(ζζᶥηηᶥ,d)

        a = (-1)^sum(δ + ζ) * (im)^sum(δ+δᶥ+ζ+ζᶥ) * pf(γ + δ) * pf(γᶥ + δᶥ) * pf(η + ζ) * pf(ηᶥ + ζᶥ)
        b = pf(γ)*pf(γᶥ)*pf(δ)*pf(δᶥ)*pf(η)*pf(ηᶥ)*pf(ζ)*pf(ζᶥ)
        @test C_block_diag.get_coef(γγᶥδδᶥ,ζζᶥηηᶥ,d) == a/b
    end
end

@testset "C_block_diag.get_ℜℑααᶥββᶥᴿ" begin
    t = (2,2)
    d = (2,2)
    Iᵗ_dict = C_block_diag.get_all_partitions(d,t)
    block_dict = C_block_diag.get_all_blocks(Iᵗ_dict)
    γγᶥδδᶥζζᶥηηᶥ_dict = C_block_diag.get_γγᶥδδᶥζζᶥηηᶥ(d,t)


    B = block_dict[(-2, 0)]
    ααᶥββᶥ = B[2,3]

    k = 1
    γγᶥδδᶥ = Tar[1][k]
    ηηᶥζζᶥ = Tar[2][k]

    γ,γᶥ,δ,δᶥ = C_block_diag.split_expo(γγᶥδδᶥ,d)
    η,ηᶥ,ζ,ζᶥ = C_block_diag.split_expo(ηηᶥζζᶥ,d)

    K = maximum([length(x[1]) for x in values(γγᶥδδᶥζζᶥηηᶥ_dict)])
    Res_vec,Res_coef = C_block_diag.get_ℜℑααᶥββᶥᴿ(d,ααᶥββᶥ,γγᶥδδᶥζζᶥηηᶥ_dict,K)

    @test Res_coef[k] == get_coef(γγᶥδδᶥ,ηηᶥζζᶥ,d)
    @test Res_vec[k]  == vcat(γ+γᶥ,δ+δᶥ,η+ηᶥ,ζ+ζᶥ)
end

@testset "C_block_diag.get_ℜℑααᶥββᶥᴿ" begin
    t = (4,4)
    d = (2,2)
    Iᵗ_dict = C_block_diag.get_all_partitions(d,t)
    block_dict = C_block_diag.get_all_blocks(Iᵗ_dict)
    γγᶥδδᶥζζᶥηηᶥ_dict = C_block_diag.get_γγᶥδδᶥζζᶥηηᶥ(d,t)
    K = maximum([length(x[1]) for x in values(γγᶥδδᶥζζᶥηηᶥ_dict)])

    for b ∈ keys(block_dict)
        B = block_dict[b]
        i_ind,j_ind =  size(B)
        for i in 1:i_ind
            for j in 1:j_ind
                ααᶥββᶥ = B[i,j]
                Res_vec,Res_coef = C_block_diag.get_ℜℑααᶥββᶥᴿ(d,ααᶥββᶥ,γγᶥδδᶥζζᶥηηᶥ_dict,K)
                @test isreal(Res_coef)
            end
        end
    end
end

## Test
# d = (3,3)
# n = sum(2 .* d)
# t = (2,2)
# Iᵗ_dict = C_block_diag.get_all_partitions(d,t)
# block_dict = C_block_diag.get_all_blocks(Iᵗ_dict)
# B = block_dict[(1,0)]
#
# γγᶥδδᶥζζᶥηηᶥ_dict = C_block_diag.get_γγᶥδδᶥζζᶥηηᶥ(d,t)
# K = maximum([length(x[1]) for x in values(γγᶥδδᶥζζᶥηηᶥ_dict)])
# # see if the coefficients are ever complex
## Test
GAH1 = Utils_cons.index_to_var(Lx, Expo_mat_dict[1,0][1],Coef_mat_dict[1,0][1])
GAH2 = Utils_cons.index_to_var(Lx, Expo_mat_dict[1,0][2],Coef_mat_dict[1,0][2])
GAH3 = Utils_cons.index_to_var(Lx, Expo_mat_dict[1,0][3],Coef_mat_dict[1,0][3])

GAHlist = Utils_cons.index_to_var(Lx, Expo_mat_dict[1,0],Coef_mat_dict[1,0])
GAHlist == GAH1 + GAH2 + GAH3
## d


Coef_mat_dict, Expo_mat_dict = C_block_diag.get_real_blocks(d,t)


model = JuMP.Model()
list_of_keys = C_constraints.make_mom_expo_keys(n,t) # Define variables in the moment matrix.
@variable(model, Lx[list_of_keys] ) #









## Gen constraints and see if they hold.
# isReal = false
# d = (3,3)
# if isReal
#     n = sum(d)
# else
#     n = sum(2 .* d)
# end
# model = JuMP.Model()
# list_of_keys = C_constraints.make_mom_expo_keys(n,t) # Define variables in the moment matrix.
# @variable(model, Lx[list_of_keys] ) #
#
# PSD_MM_con     = R_constraints.make_PSD_MM_con(d,t,Lx)
# L_xxᵀ_tens_yyᵀ = R_constraints.make_ord4_con(d,Lx)
# loc_con_S₁     = R_constraints.make_loc_cons_S₁(ρ,d,t,Lx)
# loc_con_S₂     = R_constraints.make_loc_cons_S₂(ρ,d,t,Lx)
# loc_con_S₃     = R_constraints.make_loc_cons_S₃(ρ,d,t,Lx)
# G_con          = R_constraints.make_G_con(ρ,d,t,Lx)
#
# PSD_MM_con     = C_constraints.make_PSD_MM_con(d,t,Lx)
# L_xx̄ᵀ_tens_yȳᵀ = C_constraints.make_ord4_con(d,Lx)
# loc_con_S₁     = C_constraints.make_loc_cons_S₁(ρ,d,t,Lx)
# loc_con_S₂     = C_constraints.make_loc_cons_S₂(ρ,d,t,Lx)
# loc_con_S₃     = C_constraints.make_loc_cons_S₃(ρ,d,t,Lx)
# G_con          = C_constraints.make_Gℝ_con(ρ,d,t,Lx)
#
#
# values.(PSD_MM_con["Default"])
#
# for key in keys(PSD_MM_con)
#
# end





function rec_LMM_rank(d,t,Lx,isReal=true)
    if isReal
        PSD_MM_con = R_constraints.make_PSD_MM_con(d,t,Lx)
    else
        PSD_MM_con = C_constraints.make_PSD_MM_con(d,t,Lx)
    end
    T = 0
    for key in keys(PSD_MM_con)
        T = T + rank(value.(PSD_MM_con[key]))
    end
    return T
end

rec_LMM_rank(d,t,Lx,false)

function rec_Lxx̄ᵀ_tens_yȳᵀ(d,isReal=true)
    if isReal
        xx̄ᵀ_tens_yȳᵀ = Moments.make_xxᵀ_tens_yyᵀ(d)
        return value.(xx̄ᵀ_tens_yȳᵀ)
    else
        xx̄ᵀ_tens_yȳᵀ = Moments.make_xx̄ᵀ_tens_yȳᵀ(d)
        T =  Dict()
        for k in ["real", "imag"]
            sleut = [key for key in keys(xx̄ᵀ_tens_yȳᵀ[k])]
            T[k] =   value.(xx̄ᵀ_tens_yȳᵀ[k][popfirst!(sleut)])
            for key in sleut
                if key[1] == '+'
                    T[k] = T[k] + value.(xx̄ᵀ_tens_yȳᵀ[k][key])
                elseif key[1] == '-'
                    T[k] = T[k] - value.(xx̄ᵀ_tens_yȳᵀ[k][key])
                end
            end
        end
        return T["real"] + im T["imag"]
    end

    return value.(xx̄ᵀ_tens_yȳᵀ)
end

Lxx̄ᵀ_tens_yȳᵀ = rec_Lxx̄ᵀ_tens_yȳᵀ(d,true)

loc_cons_S₁ = C_constraints.make_loc_cons_S₁(ρ,d,t,Lx)




eigvals(value.(loc_cons_S₁[("Default", 1, 1)]))











# end
# end
## Moments
# include(testDir*"\\testMoments.jl")
## Constraints
 ## R_Constraints
 ## C_Constraints


# ## Selecting all complex examples
# df_imag = filter(:isReal => ==(false),df)
# df_imag_RAND = df_imag[contains.(df_imag.ex,"RAND"),:]
#
# examples_imag = Dict()
# examples_imag_rand = Dict()
# for exa in df_imag.ex
#     if contains(exa,"RAND")
#         examples_imag_rand[exa] = examples_all[exa]
#     end
#     examples_imag[exa] = examples_all[exa]
# end

## Batch Model
t        = (2,2)
ex_set = examples_real_NOTrand
meta_df = df_real_NOTRAND
save_dir = sep_rank_proj_path*"\\PreRunModels\\R-models\\t=$(t[1])\\"
Utils_Model.batch_model(t,ex_set,meta_df,save_dir,["S₁sG","S₂sG","S₃sG"],true)

batch_Computeξₜˢᵉᵖ(save_dir)




batch_Computeξₜˢᵉᵖ(save_dir)
### Batch Run
