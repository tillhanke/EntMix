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


"""
Function returns a Dict of SMILES that contain lists of molecules with the corresponding type.
### Args:
- frame: Frame - frame of trajectory containing the atoms and their positions and types
- molecule: Vector{Int} - atom indices in Chemfiles format (starting from 0)
### Returns
- smiles: String - canonical SMILES representation of the molecule
"""
function smiles(frame::Frame, molecule::Vector{Int})
    path, io = mktemp()
    amats = length(molecule)
    write(io, "$amats\n\n")
    for atid in molecule
        write(io, "$(type(frame[atid]))  $(positions(frame)[1, atid+1]) $(positions(frame)[2, atid+1]) $(positions(frame)[3, atid+1])\n")
    end
    close(io)
    
    smiles = split(read(`sh -c "obabel -i xyz -o can $path 2> /dev/null"`, String), "\t")[1]
    return smiles
end
