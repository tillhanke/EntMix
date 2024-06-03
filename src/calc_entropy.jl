#!/usr/bin/env julia
include("Functions.jl")
include("lammpstrj_parser.jl")
include("xyz_parser.jl")
include("main.jl")

using ArgParse
using Base.Threads
using DelimitedFiles

argparser = ArgParseSettings()
@add_arg_table! argparser begin
    "--radial_factor", "-f"
        help="Factor to multiply the VdW radii of the atoms"
        default=.6
        arg_type=Float64
    "--maxstep", "-m"
        help="Maximum step to use (default all steps)"
        default=nothing
        arg_type=Int
    "--startstep", "-s"
        help="Initial step to use default to first step in file"
        default=nothing
        arg_type=Int
    "--densfunc", "-d" 
        help="Function to use for smearing. Available options are: slater, gaus, rect"
        default="slater"
        arg_type=String
    "--nonperiodic", "-n"
        help="If set, the box will be treated as non periodic. Only applicable for lammpstrj files since xyz files are always non periodic."
        action="store_true"
    "--outfile", "-o"
        help="Output file to save the entropy values, if not set, the values will be printed to the console.
Only applicable for lammpstrj files. If set threads can be used with julia -t <nthreads>."
        arg_type=String
    "--debug", "-v"
        help="Print debug information"
        action="store_true"
    "trajfile"
        help="xyz or lammpstrj file to use. 
Using a xyz file will result in a non periodic box calculation. 
Lammpstrj files will be treated as periodic if not defined otherwise with --nonperiodic."
        required=true
        arg_type=String
    "amatoms"
        help="Number of atoms in the first molecule"
        required=true
        arg_type=Int
end

args = parse_args(argparser)
if args["debug"]
    ENV["JULIA_DEBUG"] = Main
end
radial_factor = args["radial_factor"]
file = args["trajfile"]
if split(file, ".")[end] == "lammpstrj"
    trajtype = "lammpstrj"
elseif split(file, ".")[end] == "xyz"
    trajtype = "xyz"
else
    println("Invalid file type")
    println("Only .xyz and .lammpstrj files are supported")
    exit(1)
end
n_atoms = args["amatoms"]
maxstep = args["maxstep"]
startstep = args["startstep"]
outfile = args["outfile"]
if args["densfunc"] == "slater"
    dfunc = slater
elseif args["densfunc"] == "gaus"
    dfunc = gaus
elseif args["densfunc"] == "rect"
    dfunc = rect
else
    println("Invalid density function")
    exit(1)
end

# run Main function
if trajtype == "lammpstrj"
    lammpstrj_entropy(file, n_atoms, maxstep, startstep; radial_factor=radial_factor, dfunc=dfunc, outfile=outfile)
elseif trajtype == "xyz"
    xyz_entropy(file, n_atoms, maxstep, startstep; radial_factor=radial_factor, dfunc=dfunc)
end

