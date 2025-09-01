using Chemfiles: Frame, UnitCell, positions, covalent_radius, lengths, guess_bonds!, set_type!, name, type
using StaticArrays
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
#### Optional:
- attypes:: Vector{String}  Vector with length of frame with all atomtypes
### Returns:
- smiles:: String containing the Smiles code
"""
function smiles(frame, mol::Vector{Int}; startid=1, adj=nothing, attypes=nothing)
    if attypes === nothing
        smis = Chemfiles.type.(frame)
    else
        smis = attypes
    end
    visited = zeros(Bool, size(frame))
    if adj === nothing
        adj = adjacent(frame)
    end
    loop_label = Dict{Tuple{Int, Int}, String}()
    next_loopid = 1
    """
    A function that searches through the tree and breaks all loops.
    """
    function _breakloop(atid, prev)
        nbs = filter(n-> (Chemfiles.type(frame[n]) != "H") && (n !=prev), adj[atid+1])
        visited[atid+1] = true
        @debug "breakloop: $atid from $prev"
        for nb in nbs
            if visited[nb+1]
                if ! haskey(loop_label, (atid, nb))
                    if next_loopid > 9
                        loop_label[(atid, nb)], loop_label[(nb, atid)] = "%$next_loopid", "%$next_loopid"
                    else
                        loop_label[(atid, nb)], loop_label[(nb, atid)] = "$next_loopid", "$next_loopid"
                    end
                    next_loopid += 1
                    smis[atid+1] *= loop_label[(atid, nb)]
                    smis[nb+1] *= loop_label[(atid, nb)]
                    deleteat!(adj[atid+1], adj[atid+1] .== nb) # remove loop from adj_graph
                    deleteat!(adj[nb+1], adj[nb+1] .== atid) # remove loop from adj_graph
                end
            else
                _breakloop(nb, atid)
            end
        end
    end
    function _smi(atid, prev)
        @debug "smi: $atid from $prev"
        nbs = filter(n-> Chemfiles.type(frame[n]) != "H", adj[atid+1])
        deleteat!(nbs, nbs .== prev)
        nbsmis = ""
        while length(nbs) > 1
            nb = pop!(nbs)
            nbsmis *= "(" * _smi(nb, atid) * ")"
        end
        while length(nbs) > 0 
            nb = pop!(nbs)
            nbsmis *= _smi(nb, atid)
        end
        return smis[atid+1] * nbsmis
    end

    # find starting atom
    atid = mol[startid] 
    ind = startid 
    while Chemfiles.type(frame[atid]) == "H"
        ind = ((ind+1)%length(mol)) + 1
        atid = mol[ind]
    end
    _breakloop(atid, -1)
    return _smi(atid, -1)
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

function _score_smiles(smi)
    score = 0
    score -= length(collect(eachmatch(r"\)", smi)))*1.1  # i dont like brakets
    anyloop = match(r"[0-9]", smi)
    firstloop = match(r"^.[0-9][0-9]*", smi) 
    if (anyloop !== nothing) 
        if (firstloop !== nothing)
            score += 1  # prefer loops at the beginning
        else 
            score -= 1
        end
    end
    if firstloop !== nothing
        score -= parse(Int, firstloop.match[2:end])
    end
    return score
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
    attypes = Chemfiles.type.(frame)
    adj = adjacent(frame)
    all_possible_smis = Dict{Set{String}, String}()  # Set of all possible SMILES and Identification SMILES
    for (molid, mol) in enumerate(molecules)
        molname = smiles(frame, mol; adj=adj, attypes=attypes)    
        found = false
        for set in keys(all_possible_smis)
            if molname in set
                molname = all_possible_smis[set]
                found = true
                break
            end
        end
        if !found
            bestscore = -Inf
            smis = Set{String}([molname])
            for i in Vector(1:length(mol))[Chemfiles.type.(frame)[mol.+1] .!= "H"] 
                # generate all possible smiles by iterating over all possible start atoms
                smi = smiles(frame, mol; startid=i, adj=adj, attypes=attypes)
                score = _score_smiles(smi)
                if bestscore < score
                    bestscore = score
                    molname = smi
                end
                push!(smis, smi)
            end
            all_possible_smis[smis] = molname
        end
        if haskey(moltypes, molname)
            push!(moltypes[molname], molid)
        else
            moltypes[molname] = [molid]
        end
    end 
    return moltypes
end

