using Chemfiles: Frame, UnitCell, positions, covalent_radius, lengths, guess_bonds!, set_type!, name
import Chemfiles

"""
Assign types based on Atom names or list of atomtypes

type_from_name!(frame::Frame)
Assign types based on Atom names

type_from_name!(frame::Frame, types::Vector{String})
Assign types based on list of atomtypes
"""
function type_from_name!(frame::Frame; force=false)
    if Chemfiles.type(frame[1]) != "" && !force
        @warn "Atom types already set, doing nothing"
        return
    end
    for i in 0:length(frame)-1
        set_type!(view(frame, i), name(view(frame, i)))  # copy name (which is 'C', 'H', etc.) into type
    end
end

function type_from_name!(frame::Frame, types::Vector{String})
    for i in 0:length(frame)-1
        set_type!(view(frame, i), types[i+1])  
    end
end

"""
Generate list of molecules from Frame

Args:
- frame: Frame - frame to generate molecules from

Returns:
- molecules: Vector{Vector{Int}} - list of molecules with atom indices in Chemfiles format (starting from 0)
"""
function get_molecules(frame::Frame)
    try
        guess_bonds!(frame)
    catch e
        if contains(string(e), "missing Van der Waals radius for")
            @warn "try setting Atom types first"
        end
        rethrow(e)
    end
    topo = Chemfiles.Topology(frame)
    bonds = Array{Int}(Chemfiles.bonds(topo))
    adjacent = [[] for _ in 1:size(frame)]
    for bond in eachcol(bonds)
        push!(adjacent[bond[1]+1], bond[2])
        push!(adjacent[bond[2]+1], bond[1])
    end
    visited = zeros(Bool, size(frame))
    molecules = Vector{Vector{Int}}()
    for i in 0:size(frame)-1
        if !visited[i+1]
            mol = Int[]
            queue = [i]
            visited[i+1] = true
            while !isempty(queue)
                atom = pop!(queue)
                push!(mol, atom)
                for neighbor in adjacent[atom+1]
                    if !visited[neighbor+1]
                        push!(queue, neighbor)
                        visited[neighbor+1] = true
                    end
                end
            end
            push!(molecules, mol)
        end
    end
    return molecules
end

"""
Generate dictionary of atom indices to molecule indices

Args:
- frame: Frame - frame to generate dictionary from
- molecules: Vector{Vector{Int}} - list of molecules with atom indices in Chemfiles format (starting from 0)

Returns:
- at_to_mol: Dict{Int, Int} - dictionary of atom indices to molecule indices
"""
function mol_dictionary(frame::Frame, molecules::Vector{Vector{Int}})
    at_to_mol = Dict{Int, Int}()
    for (id, mol) in enumerate(molecules)
        for atom in mol
            at_to_mol[atom] = id
        end
    end
    return at_to_mol
end

"""
Generate dictionary of molecule types

Args:
- frame: Frame - frame to generate dictionary from
- molecules: Vector{Vector{Int}} - list of molecules with atom indices in Chemfiles format (starting from 0)

Returns:
- moltypes: Dict{String, Array{Int}} - dictionary of molecule types
"""
function mol_types(frame::Frame, molecules::Vector{Vector{Int}})
    moltypes = Dict{String, Array{Int}}()
    for (molid, mol) in enumerate(molecules)
        molname = ""
        for atid in mol
            molname *= String(Chemfiles.type(frame[atid]))
        end
        molname = join(sort(collect(molname)))     
        if haskey(moltypes, molname)
            push!(moltypes[molname], molid)
        else
            moltypes[molname] = [molid]
        end
    end 
    return moltypes
end

"""
Generate SMILES code for a molecule containing the specified atom

Args:
- frame: Frame - frame containing the molecular structure
- atom_id: Int - atom ID (in Chemfiles format, starting from 0) to find molecule for
- adjacent: Vector{Vector{Int}} - optional bond graph containing all bonds between atoms

Returns:
- smiles: String - SMILES code for the molecule containing the specified atom
"""
function generate_smiles(frame::Frame, atom_id::Int, adjacent::Union{Vector{Vector{Int}}, Nothing}=nothing)
    # Get molecules and find which one contains the atom
    if adjacent === nothing
        molecules = get_molecules(frame)
        # Build adjacency from detected bonds
        try
            guess_bonds!(frame)
        catch e
            if contains(string(e), "missing Van der Waals radius for")
                @warn "try setting Atom types first"
            end
            rethrow(e)
        end
        topo = Chemfiles.Topology(frame)
        bonds = Array{Int}(Chemfiles.bonds(topo))
        adj = [Int[] for _ in 1:size(frame)]
        for bond in eachcol(bonds)
            push!(adj[bond[1]+1], bond[2])
            push!(adj[bond[2]+1], bond[1])
        end
        adjacent = adj
    else
        # Use provided bond graph to find molecules
        visited = zeros(Bool, size(frame))
        molecules = Vector{Vector{Int}}()
        for i in 0:size(frame)-1
            if !visited[i+1]
                mol = Int[]
                queue = [i]
                visited[i+1] = true
                while !isempty(queue)
                    atom = pop!(queue)
                    push!(mol, atom)
                    for neighbor in adjacent[atom+1]
                        if !visited[neighbor+1]
                            push!(queue, neighbor)
                            visited[neighbor+1] = true
                        end
                    end
                end
                push!(molecules, mol)
            end
        end
    end
    
    # Find which molecule contains the atom_id
    target_molecule = nothing
    for mol in molecules
        if atom_id in mol
            target_molecule = mol
            break
        end
    end
    
    if target_molecule === nothing
        throw(ArgumentError("Atom ID $atom_id not found in any molecule"))
    end
    
    # Generate SMILES for this molecule
    return _generate_molecule_smiles(frame, target_molecule, adjacent)
end

"""
Generate dictionary of SMILES codes mapping to molecule atom collections

Args:
- frame: Frame - frame to generate SMILES dictionary from

Returns:
- smiles_dict: Dict{String, Vector{Vector{Int}}} - dictionary mapping SMILES codes to vectors of atom ID vectors for each molecule instance
"""
function molecule_smiles_dict(frame::Frame)
    molecules = get_molecules(frame)
    
    # Build adjacency list from bonds
    try
        guess_bonds!(frame)
    catch e
        if contains(string(e), "missing Van der Waals radius for")
            @warn "try setting Atom types first"
        end
        rethrow(e)
    end
    topo = Chemfiles.Topology(frame)
    bonds = Array{Int}(Chemfiles.bonds(topo))
    adjacent = [Int[] for _ in 1:size(frame)]
    for bond in eachcol(bonds)
        push!(adjacent[bond[1]+1], bond[2])
        push!(adjacent[bond[2]+1], bond[1])
    end
    
    smiles_dict = Dict{String, Vector{Vector{Int}}}()
    
    for mol in molecules
        smiles = _generate_molecule_smiles(frame, mol, adjacent)
        if haskey(smiles_dict, smiles)
            push!(smiles_dict[smiles], mol)
        else
            smiles_dict[smiles] = [mol]
        end
    end
    
    return smiles_dict
end

"""
Internal function to generate SMILES string for a molecule

Args:
- frame: Frame - frame containing the molecular structure
- molecule: Vector{Int} - list of atom indices in the molecule
- adjacent: Vector{Vector{Int}} - adjacency list for bonds

Returns:
- smiles: String - SMILES string for the molecule
"""
function _generate_molecule_smiles(frame::Frame, molecule::Vector{Int}, adjacent::Vector{Vector{Int}})
    if isempty(molecule)
        return ""
    end
    
    # For a simple implementation, we'll create a canonical molecular formula
    # This is not a true SMILES but serves as a unique identifier for the molecule type
    
    # Count atoms by type
    atom_counts = Dict{String, Int}()
    for atom_idx in molecule
        atom_type = String(Chemfiles.type(frame[atom_idx]))
        if atom_type == ""
            atom_type = String(Chemfiles.name(frame[atom_idx]))
        end
        atom_counts[atom_type] = get(atom_counts, atom_type, 0) + 1
    end
    
    # Sort atom types for canonical representation
    sorted_types = sort(collect(keys(atom_counts)))
    
    # Build molecular formula string (simplified SMILES-like representation)
    formula_parts = String[]
    for atom_type in sorted_types
        count = atom_counts[atom_type]
        if count == 1
            push!(formula_parts, atom_type)
        else
            push!(formula_parts, "$atom_type$count")
        end
    end
    
    # For a more sophisticated approach, we could implement proper SMILES generation
    # with graph traversal, but for now we use molecular formula as a unique identifier
    return join(formula_parts, "")
end

