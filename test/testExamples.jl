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
    ρ_sep = Dict()
    ρ_sep["Randd₁4d₂4"]   = Examples_sep.gen_ρ_rand(4, 20)
    ρ_sep["CDsep1d₁2d₂2"] = Examples_sep.get_ρ_CDsep(1)
    ρ_sep["CDsep2d₁2d₂2"] = Examples_sep.get_ρ_CDsep(2)
    ρ_sep["CDsep3d₁3d₂3"] = Examples_sep.get_ρ_CDsep(3)
    ρ_sep["CDsep4d₁3d₂3"] = Examples_sep.get_ρ_CDsep(4)
    ρ_sep["CDsep5d₁3d₂3"] = Examples_sep.get_ρ_CDsep(5)
    ρ_sep["DNY1d₁3d₂3"]   = Examples_sep.get_ρ_DNY1()
    ρ_sep["HHH1d₁2d₂2"]   = Examples_sep.get_ρ_HHH1(0.5)
end

@testset "Test the entangeled examples" begin
    ρ_ent = Dict()
    ρ_ent["wiki1d₁2d₂2"] = Examples_ent.get_ρ_wiki(1)
    ρ_ent["wiki2d₁3d₂3"] = Examples_ent.get_ρ_wiki(2)
    ρ_ent["HHH1d₁2d₂2"]  = Examples_ent.get_ρ_HHH1(0.75)
    ρ_ent["HHH2d₁2d₂2"]  = Examples_ent.get_ρ_HHH2(0.75)
    ρ_ent["DNY2d₁2d₂2"]  = Examples_ent.get_ρ_DNY2()
    ρ_ent["BPd₁3d₂3"]    = Examples_ent.get_ρ_BP()
    ρ_ent["HKd₁3d₂3"]    = Examples_ent.get_ρ_HK()
    ρ_ent["CD1d₁2d₂3"]   = Examples_ent.get_ρ_CD1(3)
    ρ_ent["CD1d₁2d₂4"]   = Examples_ent.get_ρ_CD1(4)
    ρ_ent["CD1d₁2d₂5"]   = Examples_ent.get_ρ_CD1(5)
    ρ_ent["CD2d₁3d₂4"]   = Examples_ent.get_ρ_CD2()
end

@testset "Test the  example tools" begin
    ρ = Examples.get_examples()
    dfρ = Examples.get_example_overview()
    println(dfρ)
end


end
