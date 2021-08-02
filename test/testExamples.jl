module testExamples
using LinearAlgebra
using Test

srcDir = dirname(dirname(@__FILE__))*"\\src\\Examples\\"
include(srcDir *"Examples_sep.jl")
include(srcDir *"Examples_ent.jl")
include(srcDir *"Examples.jl")

using .Examples_sep
using .Examples_ent
using .Examples

@testset "Test the separable examples" begin
    ρ_sep  = Examples_sep.get_sep_example()
end

@testset "Test the entangeled examples" begin
    ρ_ent  = Examples_ent.get_ent_example()
end

@testset "Test the  example tools" begin
    ρ = Examples.get_examples()
    dfρ = Examples.get_example_overview()
    println(dfρ)
end


end
