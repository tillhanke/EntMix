using LinearAlgebra
using Integrals
using ArgParse

include("MyJlFuncs.jl")

argparser = ArgParseSettings()
@add_arg_table! argparser begin
    "--radial_factor", "-f"
        help="Factor to multiply the VdW radii of the atoms (default .6)"
        default=1.
        arg_type=Float64
    "--maxstep", "-m"
        help="Maximum number of steps to use (default 0, means all steps)"
        default=0
        arg_type=Int
    "--startstep", "-s"
        help="Initial step to use"
        default=1
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
startstep = args["startstep"]

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
@debug println("Scaled VdW radii: ", radii)

# Function to use for smearing. Available options are: slater, gaus, rect
smearing = gaus
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

# traj_data = parse_traj(file, frames=maxstep)  # maxstep == 0 means all frames
# if length(traj_data) != 1
#     println("Found $(length(traj_data)) frames in the file.")
# end
println("Using $file with $n_atoms atoms in the first molecule")
println("Step\tEntropy:")
if startstep != 1
    open(file, "r") do fio
        global am_at = parse(Int, readline(fio))
        seekstart(fio)
        for i in 1:(am_at+2)*(startstep-1)
            readline(fio)
        end
        global offset = position(fio)
    end
else
    open(file, "r") do fio
        global am_at = parse(Int, readline(fio))
    end
    offset = 0
end
for step in startstep:maxstep  
    @debug "Reading step $step"
    open(file, "r") do fio
        global offset
        seek(fio, offset)
        tslines = [readline(fio) for i in 1:am_at+2]
        global frame = parse_xyz_list(tslines)
        global offset = position(fio)
        global eofio = eof(fio) 
    end


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
    @debug "Solving"
    println(step,"\t", solve(prob, HCubatureJL();reltol=1e-1, abstol=1e-1).u/prod(minmax[2]-minmax[1]))
    global eofio
    if eofio
        break
    end
end

