using Chemfiles: Frame, UnitCell, positions, covalent_radius, lengths, guess_bonds!, set_type!, name, type
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

function adjacent(frame::Frame)
    adj = Vector{Vector{Int}}(undef, size(frame))
    bonds = Chemfiles.bonds(Chemfiles.Topology(frame))
    if isempty(bonds)
        try
            guess_bonds!(frame)
            bonds = Chemfiles.bonds(Chemfiles.Topology(frame))
        catch e
            if contains(string(e), "missing Van der Waals radius for")
                @warn "try setting Atom types first"
            end
            rethrow(e)
        end
    end

    for bond in eachcol(bonds)
        if !isassigned(adj, bond[1]+1)
            adj[bond[1]+1] = Int[]
        end
        push!(adj[bond[1]+1], bond[2])
        if !isassigned(adj, bond[2]+1)
            adj[bond[2]+1] = Int[]
        end
        push!(adj[bond[2]+1], bond[1])
    end
    return adj
end

"""
Generate a simplified SMILES code for Molecule
### Args:
- frame Chemfiles Frame containing information of atoms
- mol:: Vector{Int}  Vector of atom indices (in Chemfiles format, starting at 0) within the molecule
### Returns:
- smiles:: String containing the Smiles code
"""
function smiles(frame, mol::Vector{Int})
    smi = ""
    visited = zeros(Bool, size(frame))
    adj = adjacent(frame)
    function _trav_(atid, prev) 
        # filter neighbors by Hydrogen
        neighbors = filter(n -> Chemfiles.type(frame[n]) != "H", adj[atid+1])
        # remove prev atom from neighbords
        @info neighbors
        deleteat!(neighbors, neighbors .== prev)
        visited[atid+1] = true
        if length(neighbors) == 1
            if visited[neighbors[1]+1]
                # something with loops
                @info "found loop"
            else
                return Chemfiles.type(frame[atid]) * _trav_(neighbors[1], atid)
            end
        elseif length(neighbors) == 0
            # found end of molecule
            return Chemfiles.type(frame[atid])
        end
        smi = Chemfiles.type(frame[atid])
        # manage splits
        for nb in neighbors
            if visited[nb+1]
                # something with loops
                @info "found loop"
            else
                smi *= "(" * _trav_(nb, atid) * ")"
            end
        end
        return smi
    end
    atid = mol[1] 
    while Chemfiles.type(frame[atid]) == "H"
        atid += 1
    end
    smi = _trav_(atid, -1)
    return smi
end

"""
Generate list of molecules from Frame

Args:
- frame: Frame - frame to generate molecules from

Returns:
- molecules: Vector{Vector{Int}} - list of molecules with atom indices in Chemfiles format (starting from 0)
"""
function get_molecules(frame::Frame)
    if isempty(Chemfiles.bonds(Chemfiles.Topology(frame)))
        try
            guess_bonds!(frame)
        catch e
            if contains(string(e), "missing Van der Waals radius for")
                @warn "try setting Atom types first"
            end
            rethrow(e)
        end
    end
    adj = adjacent(frame)
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
                for neighbor in adj[atom+1]
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
function mol_dictionary(molecules::Vector{Vector{Int}})
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
        molname = smiles(frame, mol)    
        if haskey(moltypes, molname)
            push!(moltypes[molname], molid)
        else
            moltypes[molname] = [molid]
        end
    end 
    return moltypes
end

