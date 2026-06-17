using Test
using Chemfiles
using EntMix

@testset "EntMix.jl" begin
    include("test_smearing.jl")
    include("test_entropy.jl")
    include("test_molecules.jl")
end
