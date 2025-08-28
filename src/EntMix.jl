module EntMix
include("Entropy.jl")

module Molecule
include("moldetect.jl")
export generate_smiles, molecule_smiles_dict, get_molecules, mol_dictionary, mol_types, type_from_name!
end

export entropy, dens, entropy_distribution, Molecule

end # module EntMix
