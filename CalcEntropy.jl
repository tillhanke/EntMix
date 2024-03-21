using Plots
using LinearAlgebra
using Integrals
using ArgParse

include("MyJlFuncs.jl")

argparser = ArgParseSettings()
@add_arg_table! argparser begin
    "--radial_factor", "-f"
        help="Factor to multiply the VdW radii of the atoms (default .6)"
        default=1
        arg_type=Float64
    "--maxstep", "-m"
        help="Maximum number of steps to use"
        default=typemax(Int)
        arg_type=Int
    "xyzfile"
        help="xyz file to use"
        required=true
        arg_type=String
    "amatoms"
        help="Number of atoms in the first molecule"
        required=true
        arg_type=Int
end

args = parse_args(argparser)
radial_factor = args["radial_factor"]
file = args["xyzfile"]
n_atoms = args["amatoms"]
maxstep = args["maxstep"]

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
smearing = gauss
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

traj_data = parse_traj(file, frames=maxstep) 
if length(traj_data) != 1
    println("Found $(length(traj_data)) frames in the file.")
end
println("Step\tEntropy:")
for (i, frame) in enumerate(traj_data)
    mol1_sig = [radii[el] for el in frame[1][1:n_atoms]]
    mol2_sig = [radii[el] for el in frame[1][n_atoms+1:end]]

    minmax = (minimum(frame[2], dims=1), maximum(frame[2], dims=1))
    prob = IntegralProblem(
        (r,args)-> s(
            r; 
            dfunc=args[1], 
            sigma1=args[2], 
            atoms1=args[3], 
            sigma2=args[4],
            atoms2=args[5] 
            ) , 
        [minmax[1][i] for i in 1:3], 
        [minmax[2][i] for i in 1:3], 
        (smearing, mol1_sig, frame[2][1:n_atoms,:], mol2_sig, frame[2][n_atoms+1:end,:])
    )
    println(i,"\t", solve(prob, HCubatureJL();reltol=1e-1, abstol=1e-1).u/prod(minmax[2]-minmax[1]))
end

