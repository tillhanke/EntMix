using Chemfiles: Frame, UnitCell, positions, covalent_radius, lengths, vdw_radius
using HCubature: hcubature
using LinearAlgebra
using StaticArrays
include("Smearing.jl")

"""
Calculate configurational entropy of mixing for molecules from Trajectory using Chemfiles.jl

Methods:
entropy_distribution(
    r::SVector{3, Float64}: space point to calculate entropy at
    frame::Frame: frame to calculate entropy from
    atomcollections::Vector{Vector{Int}}: list of atom collections to be identified as different molecule species (atom indices in Chemfiles format)
    sigma::Float64: atom broadening factor multiplied to the atomic covalent radius
    dfunc=slater: density distribution function
    )
entropy_distribution(
    r::SVector{3, Float64}: space point to calculate entropy at
    atoms_positions_collections::Vector{Vector{SVector{3, Float64}}}: list of atom positions collections (atom positions in space)
    scaled_sigma_collections::Vector{Vector{Float64}}}: list of scaled sigma collections 
    box::SVector{3, Float64}: box dimensions 
    dfunc=slater: density distribution function
    )

returns:
- s: Float64 - configurational entropy of mixing at point r
"""
function entropy_distribution(
    r::SVector{3, Float64}, 
    frame::Frame, 
    atomcollections::Vector{Vector{Int}},
    sigma::Float64;
    dfunc=slater,
    natoms=ones(length(atoms_positions_collections))  # Number of atoms in each molecule 
    )
    atoms_positions_collections = [
        [SVector{3, Float64}(positions(frame)[:,i]) for i in atoms.+1] 
        for atoms in atomcollections
    ]
    scaled_sigma_collections = [
        sigma.*covalent_radius.([@view frame[i] for i in atoms]) for atoms in atomcollections
    ]
    box = SVector{3}(lengths(UnitCell(frame)))
    return entropy_distribution(r, atoms_positions_collections, scaled_sigma_collections, box; dfunc=dfunc, natoms=natoms)
end
function entropy_distribution(
    r::SVector{3, Float64},
    atoms_positions_collections::Vector{Vector{SVector{3, Float64}}},
    scaled_sigma_collections::Vector{Vector{Float64}},
    box::SVector{3, Float64};
    dfunc=slater,
    natoms=ones(length(atoms_positions_collections))  # Number of atoms in each molecule 
)
    density_collections = zeros(length(atoms_positions_collections))
    for i in 1:length(atoms_positions_collections)
        density_collections[i] = dens(
            r, 
            atoms_positions_collections[i], 
            scaled_sigma_collections[i], 
            box; 
            dfunc=dfunc
           )/natoms[i]
    end
    total_density = sum(density_collections)
    if total_density == 0
        return 0
    end
    s = 0
    for d in density_collections
        if d != 0
            s -= d/total_density * log(d/total_density)
        end
    end
    return s
end

"""
Calculate periodic boundary distance to r for atomid in frame

Args:
- r: SVector{3, Float64} - point in space
- frame: Frame - single frame from Trajectory
- atomid: Int - atom index

Returns:
- dr_pbc: SVector{3, Float64} - periodic boundary distance to r
"""
function pbc_distance(r::SVector{3, Float64}, frame::Frame, atomid::Int)
    box = SVector{3}(lengths(UnitCell(frame)))
    r_atom = positions(frame)[atomid+1]
    return pbc_distance(r, r_atom, box)
end
function pbc_distance(atid1::Int, atid2::Int, frame::Frame)
    box = SVector{3}(lengths(UnitCell(frame)))
    r1 = SVector{3, Float64}(positions(frame)[:, atid1 + 1])
    r2 = SVector{3, Float64}(positions(frame)[:, atid2 + 1])
    return pbc_distance(r1, r2, box)
end
function pbc_distance!(
    r1::SVector{3, Float64},
    r2::SVector{3, Float64},
    box::SVector{3, Float64},
    dr_pbc::Vector{Float64}
)
    for ind in 1:3
        dr_pbc[ind] = mod((r1[ind] - r2[ind] + box[ind]/2) , box[ind]) - box[ind]/2
    end
    return dr_pbc
end
function pbc_distance(
    r1::SVector{3, Float64},
    r2::SVector{3, Float64},
    box::SVector{3, Float64}
)
    dr_pbc = zeros(3)
    pbc_distance!(r1, r2, box, dr_pbc)
    return dr_pbc
end


"""
Calculate density distribution at point r for collection of atoms in frame

Args:
- r: SVector{3, Float64} - point in space
- frame: Frame - single frame from Trajectory
- atoms: Vector{Int} - atom indices
- sigma: Float64 - atom broadening factor multiplied to the atomic covalent radius

Kwargs:
- dfunc: Function - density distribution function

returns:
- d: Float64 - density distribution at point r
"""
function dens(r:: SVector{3, Float64}, frame::Frame, atoms::Vector{Int}, sigma::Float64; dfunc=slater)
    d = 0
    atom_positions = [SVector{3, Float64}(positions(frame)[:,i]) for i in atoms.+1] 
    scaled_sigma = sigma.*covalent_radius.([@view frame[i] for i in atoms])
    return dens(r, atom_positions, scaled_sigma, SVector{3}(lengths(UnitCell(frame))); dfunc=dfunc)
end
function dens(
    r::SVector{3, Float64},
    atom_positions::Vector{SVector{3, Float64}},
    scaled_sigma::Vector{Float64},
    box::SVector{3, Float64};
    dfunc=slater
)
    d = 0
    dist = zeros(3)
    for index in eachindex(scaled_sigma)
        pbc_distance!(r, atom_positions[index], box, dist)
        d += dfunc(norm(dist), scaled_sigma[index])
    end
    return d
end

"""
Calculate total entropy for the box

Args:
- frame: Frame - single frame from Trajectory
- atomcollections: Vector{Vector{Int}} - atom collections to be identified as different molecule species
- sigma: Float64 - atom broadening factor multiplied to the atomic covalent radius

Kwargs:
- baselength: String - "Covalent", "VDW", "homo" - base length for the atom broadening factor with homo being 1 for all atoms
"""
function entropy(frame::Frame, atomcollections::Vector{Vector{Int}}, sigma::Float64; baselength="Covalent", dfunc=slater, natoms=ones(length(atomcollections)))
    baselength = lowercase(baselength)
    if startswith(baselength, "cov") 
        scaled_sigma_collections = [
            sigma.*covalent_radius.([@view frame[i] for i in atoms]) for atoms in atomcollections
        ]
    elseif startswith(baselength, "vdw")
        scaled_sigma_collections = [
            sigma.*vdw_radius.([@view frame[i] for i in atoms]) for atoms in atomcollections
        ]
    elseif startswith(baselength, "homo")
        scaled_sigma_collections = [
            [sigma for i in atoms] for atoms in atomcollections
        ]
    else
        @error "Invalid baselength: $baselength"
    end
    args = (
        atoms_positions_collections = [
            [SVector{3, Float64}(positions(frame)[:,i]) for i in atoms.+1] 
            for atoms in atomcollections
        ],
        scaled_sigma_collections = scaled_sigma_collections,
        box = SVector{3}(lengths(UnitCell(frame)))
    )
    integral = hcubature(
        (r) -> entropy_distribution(r, args...; dfunc=dfunc, natoms),
        SVector{3, Float64}(zeros(3)),
        lengths(UnitCell(frame)),
        rtol=1e-2, maxevals=20000, initdiv=7
    )
    return integral[1]/prod(lengths(UnitCell(frame)))
end
function entropy(atoms::Vector{Vector{SVector{3, Float64}}}, sigmas, dfunc, box::SVector{3, Float64})
    args = (
        atoms_positions_collections = atoms,
        scaled_sigma_collections = sigmas,
        box = box
    )
    integral = hcubature(
        (r) -> entropy_distribution(r, args...; dfunc=dfunc),
        SVector{3, Float64}(zeros(3)),
        box,
        rtol=1e-2, maxevals=20000, initdiv=7
    )
    return integral[1]/prod(box)
end
