using Integrals

BOHR = 1.8897259886
# Van der Waals radii
# Values from: https://en.wikipedia.org/wiki/Atomic_radii_of_the_elements_(data_page)
VDWradii = Dict(
    "H"=>1.2 , 
    "C"=>1.70, 
    "F"=>1.47,
    "O"=>1.52,
    "N"=>1.55,
    "S"=>1.80,
    "Cl"=>1.75,
    "Br"=>1.85,
    "Si"=>2.10,
)

"""
Calculate entropy distribution at point r

Args:
- r: Array{Float64,1} - point in space
- atoms1: Array{Float64,2} - array of atom coordinates for first molecule
- sigma1: Array{Float64,1} - array of atom broadening values for first molecule
- atoms2: Array{Float64,2} - array of atom coordinates for second molecule
- sigma2: Array{Float64,1} - array of atom broadening values for second molecule

Optional:
- dfunc: Function - density distribution function
   - must take 3 arguments: r, r0, sigma
   - r: Array{Float64,2} - point in space
   - r0: Array{Float64,2} - atom coordinates
   - sigma: Float64 - atom broadening value
- na1: Int - number of atoms in the first molecule
- na2: Int - number of atoms in the second molecule

returns:
- s: Float64 - entropy at point r
"""
function entropy_distribution(r, atoms1, sigma1, atoms2, sigma2; dfunc=slater, na1=1, na2=1)
    d_a = dens(r, dfunc=dfunc, sigma=sigma1, atoms=atoms1)/na1
    d_b = dens(r, dfunc=dfunc, sigma=sigma2, atoms=atoms2)/na2
    s = 0
    if d_a != 0
        xa = d_a/(d_a+d_b)
        s -= xa*log(xa) 
    end
    if d_b != 0
        xb = d_b/(d_a+d_b)
        s -= xb*log(xb) 
    end
    return s
end

"""
Integral Problem compatible
"""
function entropy_distribution(r; dfunc=slater, atoms1=[0 0 0], sigma1=[.1], atoms2=[1 1 1], sigma2=[.1], na1=1, na2=1)
    return entropy_distribution(r, atoms1, sigma1, atoms2, sigma2; dfunc=slater, na1=1, na2=1)
end

"""
Calculate molecule density at point r

Args:
- r: Array{Float64,1} - point in space
- dfunc: Function - density distribution function
  - must take 3 arguments: r, r0, sigma
  - r: Array{Float64,2} - point in space
  - r0: Array{Float64,2} - atom coordinates
  - sigma: Float64 - atom broadening value
- atoms: Array{Float64,2} - array of atom coordinates
- sigma: Array{Float64,1} - array of atom broadening values
returns:
- d: Float64 - density at point r
"""
function dens(r; dfunc=slater, atoms=[0 0 0], sigma=[.1])
    d = 0
    for (r0, sig) in zip(eachrow(atoms), sigma)
        d += dfunc([r[1] r[2] r[3]], r0, sig)
    end
    return d
end

"""
Slater type density distribution function

Args:
- r: Array{Float64,2} - point in space
- r0: Array{Float64,2} - atom coordinates
- sigma: Float64 - atom broadening value
- n: Int - power of the function
returns:
- d: Float64 - density at point r for atom at r0
"""
function slater(r, r0, sigma; n=1)
    norm = 1/(factorial(n+1)*sigma^(n+2)*4*pi)
    dist = broadcast(-,r,transpose(r0)) 
    dist = sum(dist.*dist, dims=2).^(1/2)
    if size(dist) == (1,1)
        dist = dist[1]
    end
    return @. norm* dist^(n-1)*exp(-dist/sigma)
end

"""
Gaussian density distribution function

Args:
- r: Array{Float64,2} - point in space
- r0: Array{Float64,2} - atom coordinates
- sigma: Float64 - atom broadening value
returns:
- g: Float64 - density at point r for atom at r0
"""
function gaus(r, r0, sigma)
    dist = broadcast(-,r,transpose(r0))   
    g = 1/(sigma^3*(pi*2)^(3/2))*exp.(-sum(dist .* dist, dims=2) / (2*sigma^2))
    if size(g) == (1,1)
        g = g[1]
    end
    return g
end

"""
Calculate configurational Entropy of mixing for two molecules

Args:
* (coords_m1, elems_m1): Tuple{Array{Float64,2}, Array{String,1}} - coordinates and elements of the first molecule
* (coords_m2, elems_m2): Tuple{Array{Float64,2}, Array{String,1}} - coordinates and elements of the second molecule
* radial_factor: Float64 - scaling factor for Van der Waals radii
* dfunc: Function - density distribution function
  + must take 3 arguments: r, r0, sigma
  + r: Array{Float64,2} - point in space
  + r0: Array{Float64,2} - atom coordinates
  + sigma: Float64 - atom broadening value
* periodic: Bool - whether to use periodic boundary conditions
* box: Array{Float64,2} - box dimensions
* na1: Int - number of atoms in the first molecule
* na2: Int - number of atoms in the second molecule

Returns:
* Float64 - entropy of mixing
"""
function entropy(m1, m2; radial_factor=0.6, dfunc=slater, periodic=false, box=Nothing, na1=1, na2=1)

    scaled_vdw = copy(VDWradii)
    map!(x->x*BOHR*radial_factor, values(scaled_vdw))
    @debug "Scaled VdW radii:" scaled_vdw 

    if periodic
        if box == Nothing
            error("Periodic boundary conditions require box dimensions")
        end
        m1_p = add_periodic(m1[1], m1[2], box; delta=radial_factor)
        m2_p = add_periodic(m2[1], m2[2], box; delta=radial_factor)
        if size(m1_p[1])[1] != 0
            m1_coords = [ 
                         m1[1]
                         m1_p[1] 
                        ]
            m1_elems = [
                        m1[2]
                        m1_p[2]
                       ]
        end
        if size(m2_p[1])[1] != 0
            m2_coords = [
                         m2[1]
                         m2_p[1]
                        ]
            m2_elems = [
                        m2[2]
                        m2_p[2]
                       ]
        end
        @debug "Added periodic atoms" size(m1_p[1]), size(m2_p[1])
        
    else
        @debug "Molecule 1" m1

        box = [minimum([m1[1]
                        m2[1]
                       ], dims=1)
               maximum([m1[1]
                        m2[1]
                       ], dims=1)]
        m1_coords = m1[1]
        m1_elems = m1[2]
        m2_coords = m2[1]
        m2_elems = m2[2]
    end
    @debug  box

    # minmax = (minimum(frame[2], dims=1), maximum(frame[2], dims=1))
    @debug "Integral box:" [box[1, i] for i in 1:3], [box[2, i] for i in 1:3] 

    prob = IntegralProblem(
        (r, kwargs)-> entropy_distribution(
            r; 
            kwargs...
        ), 
        [box[1, i] for i in 1:3], 
        [box[2, i] for i in 1:3], 
        (
            dfunc=dfunc, 
            sigma1=[scaled_vdw[el] for el in m1_elems],
            atoms1=m1_coords, 
            sigma2=[scaled_vdw[el] for el in m2_elems],
            atoms2=m2_coords,
            na1=na1,
            na2=na2
        )
    )
    @debug "Solving"
    @debug "Integral volume" prod(box[2,:]-box[1,:])
    solution=solve(prob, HCubatureJL();reltol=1e-1, abstol=1e-1)    
    @debug "Solution" solution
    return solution.u/prod(box[2,:]-box[1,:])
end

"""
Add atoms near the periodic boundary over at the other side

Args:
- atoms: Array{Float64,2} - array of atom coordinates
- elems: Array{String,1} - array of atom elements
- box: Array{Float64,1} - box dimensions
- delta: Float64 - distance from the boundary to wrap atoms

Returns:
additional atoms and their elements
- atoms: Array{Float64,2} - array of atom coordinates
- elems: Array{String,1} - array of atom elements
"""
function add_periodic(atoms, elems, box; delta=0.1)
    new_ats = []
    new_elems = []
    basevec = box[2,:]-box[1,:]
    for (element, coord) in zip(elems, eachrow(atoms))
        negs = coord.-box[1,:] .< delta
        pos = abs.(coord.-box[2,:]) .< delta 
        drs = []
        for (n, dr) in zip(negs, eachrow([1 0 0; 0 1 0; 0 0 1]))
            if n
                push!(drs, dr.*basevec) 
            end
        end
        for i in 1:length(drs)
            push!(new_ats, coord+drs[i])
            push!(new_elems, element)
            for j in i+1:length(drs)
                push!(new_ats, coord+drs[i]+drs[j])
                push!(new_elems, element)
                for k in j+1:length(drs)
                    push!(new_ats, coord+drs[i]+drs[j]+drs[k])
                    push!(new_elems, element)
                end
            end
        end
        drs = []
        for (n, dr) in zip(pos, eachrow([-1 0 0; 0 -1 0; 0 0 -1]))
            if n
                push!(drs, dr.*basevec) 
            end
        end
        for i in 1:length(drs)
            push!(new_ats, coord+drs[i])
            push!(new_elems, element)
            for j in i+1:length(drs)
                push!(new_ats, coord+drs[i]+drs[j])
                push!(new_elems, element)
                for k in j+1:length(drs)
                    push!(new_ats, coord+drs[i]+drs[j]+drs[k])
                    push!(new_elems, element)
                end
            end
        end
    end
    if length(new_ats) == 0
        return [], []
    end
    new_ats = stack(new_ats, dims=1)
    return new_ats, new_elems
end

"""
Wrap atoms into the box

Args:
- atoms: Array{Float64,2} - array of atom coordinates
- box: Array{Float64,1} - box dimensions in the form [x_min y_min z_min; x_max y_max z_max]

Returns:
- atoms: Array{Float64,2} - updated array of atom coordinates inside the box
"""
function wrap!(atoms, box)
    for i in 1:size(atoms,1)
        for j in 1:3
            while atoms[i,j] < box[1,j]
                atoms[i,j] += box[2,j]-box[1,j]
                # if atoms[i,j] < box[1,j]
                #     error("Atom $i is still outside the box")
                # end
            end
            while atoms[i,j] > box[2,j]
                atoms[i,j] -= box[2,j]-box[1,j]
                # if atoms[i,j] > box[2,j]
                #     error("Atom $i is still outside the box")
                # end
            end
        end
    end
    return atoms
end

"""
Wrap atoms into the box

Args:
- atoms: Array{Float64,2} - array of atom coordinates
- box: Array{Float64,1} - box dimensions in the form [x_min y_min z_min; x_max y_max z_max]

Returns:
- atoms: Array{Float64,2} - updated array of atom coordinates inside the box
"""
function wrap(atoms, box)
    new_atoms = copy(atoms)
    return wrap!(new_atoms, box)
end

