module testExamples
include(pwd()*"\\src\\data.jl")
using .Examples
import Test

@testset "random sep states" begin
    d = rand(3:5)
    r = rand(2:6)
    ρ_rand = gen_rand_state(d,r,343)
    @test size(ρ_rand) == (d^2,d^2)
    @test tr(ρ_rand) == 1
    @test issymmetric(ρ_rand)
    @test ishermitian(ρ_rand)

    ρ_rand_batch = generate_random_states()
    key_samp = rand(keys(ρ_rand_batch),5)
    for key in key_samp
        ρ_rand = ρ_rand_batch[key]
        d = key[1]
        @test size(ρ_rand) == (d^2,d^2)
        @test tr(ρ_rand) == 1
        @test issymmetric(ρ_rand)
        @test ishermitian(ρ_rand)
    end
end

@testset "lit example utils" begin
    i₀ = rand(1:d); i₁ = rand(1:d)
    e₀ = get_std_base_vec(n,i₀)
    e₁ = get_std_base_vec(n,i₁)
    @test ψ(i₀,i₁,d) == kron(e₀,e₁)
    e = e₀ + e₁
    @test ψ(d) == kron(e,e)

    sq(ϕ)  ==
    sqs(ϕ₁,ϕ₂) ==
    psep(i₀::Int,i₁::Int,n::Int)  ==
end

@testset "lit example sep states" begin
    get_sep_example(d)
end

@testset "lit example ent states" begin
    get_ent_example(d)
end

end
