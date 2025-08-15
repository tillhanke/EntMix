module EntMix
include("Entropy.jl")

module Molecule
include("moldetect.jl")
end

export entropy, dens, entropy_distribution, Molecule

end # module EntMix
