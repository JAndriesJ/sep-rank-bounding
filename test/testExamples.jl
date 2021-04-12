module testExamples
include(pwd()*"\\src\\Examples.jl")
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
    sep_examples = Examples.get_sep_example("tr")


    for key in keys(sep_examples)
        ρ = sep_examples[key]
        @test ishermitian(ρ)
        d = Int(sqrt(size(ρ)[1]))
        ρ_1  =  Examples.partial_transpose_per(ρ, 1, [d,d])
        @test ishermitian(ρ_1)
        ρ_2  =  Examples.partial_transpose_per(ρ, 2, [d,d])
        @test ishermitian(ρ_2)
    end

    H_1 = sep_examples[2, 4, "sep6"]
    λ₁ = λ₂ = 42 ;
    u¹₁ = u²₂ = [sqrt(14)/14, sqrt(14)/7, 3/sqrt(14)] ;
    u²₁ = u¹₂ = [sqrt(3)/3, sqrt(3)/3, sqrt(3)/3] ;
    H_2 = λ₁ * kron(u¹₁ *transpose(u¹₁) , u²₁*transpose(u²₁))  + λ₂ * kron(u¹₂ *transpose(u¹₂), u²₂*transpose(u²₂))
    @test H_1 == H_2
end

@testset "lit example ent states" begin
    get_ent_example(d)
end

end
