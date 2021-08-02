

include(dirname(@__FILE__)*"\\Moments\\Moments.jl")
include(dirname(@__FILE__)*"\\Constraints\\Utils_cons.jl")
using .Moments
using .Utils_cons

d = (2,2)
t = (2,2)
n = sum(d.*2)

xx̄yȳMM_blocks = get_xx̄yȳMM_blocks(d,t)
γγᶥζζᶥ_δδᶥηηᶥ = get_γγᶥζζᶥ_δδᶥηηᶥ(d,t)


I²₀₋₁ =  xx̄yȳMM_blocks[(0,-1)]
display(I²₀₋₁)
split_expo(ααᶥββᶥ,d) =     (ααᶥββᶥ[1:d[1]],
                            ααᶥββᶥ[d[1]+1:2*d[1]],
                            ααᶥββᶥ[1+2*d[1]:2*d[1]+d[2]],
                            ααᶥββᶥ[1+2*d[1]+d[2]:end])
display(map(x->split_expo(x,d),I²₀₋₁))
display(γγᶥζζᶥ_δδᶥηηᶥ[I²₀₋₁[1,1]])


pf(x) = prod(factorial.(x))
function get_coef(part1,part2,d)
    γ,γᶥ,ζ,ζᶥ = split_expo(part1,d) ; δ,δᶥ,η,ηᶥ = split_expo(part2,d)
    a = ((-1.0)^sum(δᶥ+ηᶥ))*((1.0im)^sum(δ+δᶥ+η+ηᶥ))*prod(pf.([[γ,γᶥ,ζ,ζᶥ]+[δ,δᶥ,η,ηᶥ]...]))
    return a/prod(pf.([γ,γᶥ,ζ,ζᶥ,δ,δᶥ,η,ηᶥ]))
end
get_coef(γ_δ_pair,d) = map((x,y)->get_coef(x,y,d),γ_δ_pair...)

get_coef(γγᶥζζᶥ_δδᶥηηᶥ[I²₀₋₁[1,1]],d)


MMexᴿ,MMCoefᴿ = Moments.get_ℜℑααᶥββᶥᴿ(d,I²₀₋₁,γγᶥζζᶥ_δδᶥηηᶥ)

MMex = map(x->get_ℜℑααᶥββᶥᴿ(d,x,γγᶥζζᶥ_δδᶥηηᶥ)[1],I²₀₋₁)
MMCoef = map(x->get_ℜℑααᶥββᶥᴿ(d,x,γγᶥζζᶥ_δδᶥηηᶥ)[2],I²₀₋₁)















using JuMP
model = JuMP.Model()
@variable(model, Lx[Moments.make_mon_expo(d,t[1]*2)] ) ## Create variables


MMexᴿ,MMCoefᴿ = Moments.get_ℂ_block_diag(d,t .- 1)
PSD_con = Utils_cons.idx2var(Lx,MMexᴿ,MMCoefᴿ)










B = MMexᴿ[(-1,0)]
C = MMCoefᴿ[(-1,0)]

b = B[1,1]
c = C[1,1]

Utils_cons.idx2var(Lx,b)
Utils_cons.idx2var(Lx,b,c)
Utils_cons.idx2var(Lx,b,c,[Utils_cons.eᵢ(n,0)])

Gar  = Utils_cons.idx2var_arr(Lx,B,C)
Gary = Utils_cons.idx2var_arr(Lx,B,C,[Utils_cons.eᵢ(n,1)])

DGar  = Utils_cons.idx2var(Lx,MMexᴿ,MMCoefᴿ)
DGary = Utils_cons.idx2var(Lx,MMexᴿ,MMCoefᴿ,[Utils_cons.eᵢ(n,1)])

DGar[(-1,0)] == Gar
DGary[(-1,0)] == Gary


# @test split_expo([1,1,1,2,2,2,3,3,3,4,4,4],(3,3)) == ([1, 1, 1], [2, 2, 2], [3, 3, 3], [4, 4, 4])
# γ,γᶥ,δ,δᶥ = split_expo([2,0,0,0,0,0,1,0,0,0,0,0],(3,3))
# η,ηᶥ,ζ,ζᶥ = split_expo([0,0,0,0,0,0,2,0,0,0,0,0],(3,3))
# a = ((-1)^sum(δᶥ + ζᶥ) * (im)^sum(δ+δᶥ+ζ+ζᶥ))* pf(γ + δ) * pf(γᶥ + δᶥ) * pf(η + ζ) * pf(ηᶥ + ζᶥ)
# b = pf(γ)*pf(γᶥ)*pf(δ)*pf(δᶥ)*pf(η)*pf(ηᶥ)*pf(ζ)*pf(ζᶥ)
# @test get_coef([2,0,0,0,0,0,1,0,0,0,0,0],[0,0,0,0,0,0,2,0,0,0,0,0],(3,3)) ==  a/b
# asdf = [keys(γγᶥδδᶥζζᶥηηᶥ)...]
# @test length(asdf) == binomial(sum(2 .* d)+sum(t),sum(t))

# MMB = get_xx̄yȳMM_blocks(d,t)
# γγᶥδδᶥζζᶥηηᶥ = get_γγᶥδδᶥζζᶥηηᶥ(d,t)
#
# MMexᴿ,MMCoefᴿ = get_ℜℑααᶥββᶥᴿ(d,MMB[(-1,-1)],γγᶥδδᶥζζᶥηηᶥ)
#
# γγᶥδδᶥζζᶥηηᶥ[MMB[(-1,-1)][1,1]]
#



include(srcDir*"Examples\\Utils_states.jl")
using .Utils_states
using LinearAlgebra

ρ =    [2 3 4 3 4 5 4 5 6
        3 5 7 4 6 8 5 7 9
        4 7 10 5 8 11 6 9 12
        3 4 5 5 6 7 7 8 9
        4 6 8 6 8 10 8 10 12
        5 8 11 7 10 13 9 12 15
        4 5 6 7 8 9 10 11 12
        5 7 9 8 10 12 11 13 15
        6 9 12 9 12 15 12 15 18]
ρᵀᵇ =  Utils_states.take_Pᵀ(ρ,1,(3,3))

m = diagm(sqrt.(1 ./ diag(ρ)))
mρm = m*ρ*m
mρmᵀᵇ =  Utils_states.take_Pᵀ(mρm,1,(3,3))

eigvals(ρ)[1]
eigvals(ρᵀᵇ)[1]

eigvals(m*ρ*m)[1]
eigvals(m*ρᵀᵇ*m)[1]

eigvals(mρmᵀᵇ)[1]
