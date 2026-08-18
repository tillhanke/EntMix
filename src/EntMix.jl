module EntMix
include("Entropy.jl")


module Molecules
include("moldetect.jl")
end

export entropy, dens, entropy_distribution, Molecules
end # module EntMix
