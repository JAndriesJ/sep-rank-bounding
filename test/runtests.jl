module runTests
sep_rank_proj_path = dirname(dirname(@__FILE__))
# Pkg.activate(sep_rank_proj_path)

testDir = dirname(@__FILE__)
include(testDir*"\\testUtils_states.jl")
include(testDir*"\\testExamples.jl")
include(testDir*"\\testUtils.jl")
include(testDir*"\\testMoments.jl")
# include(testDir*"\\testSep_Constriants.jl")
include(testDir*"\\testSep_Model.jl")
include(testDir*"\\testSep_Compute.jl")

end
