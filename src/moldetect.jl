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
    if size(frame) == 0
        return  # Empty frame, nothing to do
    end
    if Chemfiles.type(frame[0]) != "" && !force
        @warn "Atom types already set, doing nothing"
        return
    end
    for i in 0:size(frame)-1
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
    # Get molecules first to determine what atoms exist
    if adjacent === nothing
        try
            molecules = get_molecules(frame)
            
            # Check if frame is empty after molecule detection
            if isempty(molecules)
                throw(ArgumentError("No molecules found in frame"))
            end
            
            # Check if atom_id is valid based on molecules found
            max_atom_id = maximum(maximum(mol) for mol in molecules if !isempty(mol))
            min_atom_id = minimum(minimum(mol) for mol in molecules if !isempty(mol))
            
            if atom_id < min_atom_id || atom_id > max_atom_id
                throw(ArgumentError("Atom ID $atom_id not found. Valid atom IDs: $min_atom_id to $max_atom_id"))
            end
        catch e
            # If get_molecules fails, try to work with manual frame construction
            if adjacent === nothing
                throw(ArgumentError("Cannot detect molecules in frame and no adjacency provided: $e"))
            end
        end
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
    
    # Check if frame is empty
    if isempty(molecules)
        return Dict{String, Vector{Vector{Int}}}()
    end
    
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
    
    # For small molecules, try to generate a more SMILES-like representation
    if length(molecule) <= 10
        return _generate_simple_smiles(frame, molecule, adjacent)
    else
        # For larger molecules, fall back to molecular formula
        return _generate_molecular_formula(frame, molecule)
    end
end

"""
Generate a simple SMILES-like representation for small molecules
"""
function _generate_simple_smiles(frame::Frame, molecule::Vector{Int}, adjacent::Vector{Vector{Int}})
    # Create a subgraph for this molecule
    atom_to_local = Dict{Int, Int}()
    for (local_idx, atom_idx) in enumerate(molecule)
        atom_to_local[atom_idx] = local_idx
    end
    
    # Build local adjacency list
    local_adjacent = [Int[] for _ in 1:length(molecule)]
    for (local_idx, atom_idx) in enumerate(molecule)
        for neighbor in adjacent[atom_idx + 1]  # Convert to 1-indexed
            if neighbor in keys(atom_to_local)
                push!(local_adjacent[local_idx], atom_to_local[neighbor])
            end
        end
    end
    
    # Find a good starting atom (prefer carbons, then heteroatoms, then hydrogens)
    start_atom = _find_canonical_start_atom(frame, molecule)
    start_local = atom_to_local[start_atom]
    
    # Generate SMILES using DFS traversal
    visited = falses(length(molecule))
    smiles_parts = String[]
    
    function dfs(local_idx::Int, parent::Int = -1)
        if visited[local_idx]
            return
        end
        
        visited[local_idx] = true
        atom_idx = molecule[local_idx]
        atom_symbol = String(Chemfiles.type(frame[atom_idx]))
        if atom_symbol == ""
            atom_symbol = String(Chemfiles.name(frame[atom_idx]))
        end
        
        # Count unvisited neighbors
        unvisited_neighbors = Int[]
        for neighbor_local in local_adjacent[local_idx]
            if neighbor_local != parent && !visited[neighbor_local]
                push!(unvisited_neighbors, neighbor_local)
            end
        end
        
        # Add atom symbol
        push!(smiles_parts, atom_symbol)
        
        # Add branches
        if length(unvisited_neighbors) > 1
            for (i, neighbor) in enumerate(unvisited_neighbors)
                if i > 1
                    push!(smiles_parts, "(")
                end
                dfs(neighbor, local_idx)
                if i > 1
                    push!(smiles_parts, ")")
                end
            end
        elseif length(unvisited_neighbors) == 1
            dfs(unvisited_neighbors[1], local_idx)
        end
    end
    
    dfs(start_local)
    
    # If not all atoms visited, add remaining components
    for local_idx in 1:length(molecule)
        if !visited[local_idx]
            push!(smiles_parts, ".")
            dfs(local_idx)
        end
    end
    
    smiles = join(smiles_parts, "")
    
    # If the SMILES is getting too complex, fall back to molecular formula
    if length(smiles) > 50 || count(==('.'), smiles) > 0
        return _generate_molecular_formula(frame, molecule)
    end
    
    return smiles
end

"""
Find the most suitable starting atom for SMILES generation
"""
function _find_canonical_start_atom(frame::Frame, molecule::Vector{Int})
    # Priority: C > N > O > others > H
    priority_order = ["C", "N", "O", "S", "P", "F", "Cl", "Br", "I"]
    
    for symbol in priority_order
        for atom_idx in molecule
            atom_type = String(Chemfiles.type(frame[atom_idx]))
            if atom_type == ""
                atom_type = String(Chemfiles.name(frame[atom_idx]))
            end
            if atom_type == symbol
                return atom_idx
            end
        end
    end
    
    # If no priority atoms found, return first non-hydrogen
    for atom_idx in molecule
        atom_type = String(Chemfiles.type(frame[atom_idx]))
        if atom_type == ""
            atom_type = String(Chemfiles.name(frame[atom_idx]))
        end
        if atom_type != "H"
            return atom_idx
        end
    end
    
    # If only hydrogens, return the first one
    return molecule[1]
end

"""
Generate molecular formula as fallback
"""
function _generate_molecular_formula(frame::Frame, molecule::Vector{Int})
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
    
    # Build molecular formula string
    formula_parts = String[]
    for atom_type in sorted_types
        count = atom_counts[atom_type]
        if count == 1
            push!(formula_parts, atom_type)
        else
            push!(formula_parts, "$atom_type$count")
        end
    end
    
    return join(formula_parts, "")
end

