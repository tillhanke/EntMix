using Plots
using LinearAlgebra
using Integrals

include("MyJlFuncs.jl")

# --------------------------
# INPUT
# --------------------------
# 
# Files of Molecular xyz
# FILE1 = "220k-sep-mol1.xyz"
# FILE2 = "220k-sep-mol2.xyz"

# Radial factor for the smearing. It will be multiplied by covalent radii of each atom
# Tested values are: slater .3
radial_factor = .6

if length(ARGS) == 2
    println("Using input files from command line")
    FILE1 = ARGS[1]
    FILE2 = ARGS[2]
elseif length(ARGS) == 3
    println("Using input files from command line")
    FILE1 = ARGS[1]
    FILE2 = ARGS[2]
    radial_factor = parse(Float64, ARGS[3])
    println("Using radial factor: ", radial_factor)
else
    println("Usage: julia CalcEntropy.jl FILE1 FILE2 [radial_factor]")
    println("FILE1: xyz file of the first molecule")
    println("FILE2: xyz file of the second molecule")
    println("radial_factor: factor to multiply the VdW radii of the atoms (default .6)")
    println("Exiting")
    exit()
end
# Van der Waals radii
# Values from: https://en.wikipedia.org/wiki/Atomic_radii_of_the_elements_(data_page)
bohr = 1.8897259886

radii = Dict(
    "H"=>1.2*bohr*radial_factor, 
    "C"=>1.70*bohr*radial_factor, 
    "F"=>1.47*bohr*radial_factor,
    "O"=>1.52*bohr*radial_factor,
    "N"=>1.55*bohr*radial_factor,
    "S"=>1.80*bohr*radial_factor,
    "Cl"=>1.75*bohr*radial_factor,
    "Br"=>1.85*bohr*radial_factor,
)
println("Scaled VdW radii: ", radii)

# Function to use for smearing. Available options are: slater, gauss, rect
smearing = slater
# --------------------------
# END INPUT
# --------------------------

# Define the function
# Density based on smearing

function dens(r; dfunc=slater, atoms=[0 0 0], sigma=[.1])
    d = 0
    for (r0, sig) in zip(eachrow(atoms), sigma)
        d += dfunc([r[1] r[2] r[3]], r0, sig)
    end
    return d
end

# Entropy
function s(r; dfunc=slater, atoms1=[0 0 0], sigma1=[.1], atoms2=[1 1 1], sigma2=[.1])
    d_a = dens(r, dfunc=dfunc, sigma=sigma1, atoms=atoms1)
    d_b = dens(r, dfunc=dfunc, sigma=sigma2, atoms=atoms2)
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

# Load files

mol1 = parse_xyz(FILE1)
mol2 = parse_xyz(FILE2)

mol1_sig = [radii[el] for el in mol1[1]]
mol2_sig = [radii[el] for el in mol2[1]]

minmax = (minimum(mol1[2][:,1]), maximum(mol1[2][:,1]))

prob = IntegralProblem(
    (r,args)-> s(
        r; 
        dfunc=args[1], 
        sigma1=args[2], 
        atoms1=args[3], 
        sigma2=args[4],
        atoms2=args[5] 
        ) , 
    ones(3)*minmax[1], 
    ones(3)*minmax[2], 
    (smearing, mol1_sig, mol1[2], mol2_sig, mol2[2])
)


println("")
println("Calculating Entropy")
println("Entropy:\t", solve(prob, HCubatureJL();reltol=1e-1, abstol=1e-1).u/((minmax[2]-minmax[1])^3))
