# Mixing Entropy Calculation based on MD Trajectory

This module evaluates trajectory data from md simulations and calculates the configurational entropy distribution.

Currently only binary mixtures are supported.

## Installation
Clone the repository and install the julia module:
```bash
git clone [repo]
julia
```
then install the module from the directory
```julia
using Pkg
Pkg.add("./EntMix")
```

## Usage
Using the provided script allows to calculate entropy from lammpstrj or xyz files directly:
```bash
julia EntMix/src/calc_entropy.jl -m 100 traj.lammpstrj 300
```
This will calculate the mixing entropy of the first 100 timesteps within `traj.lammpstrj`. 
The two molecules must be sorted identically within each timestep and the first 300 atoms are all part of the first molecule type.

A complete usage is provided via
```bash
julia calc_entropy.jl --help
  > usage: calc_entropy.jl [-f RADIAL_FACTOR] [-m MAXSTEP] [-s STARTSTEP]
  >                        [-d DENSFUNC] [-n] [-o OUTFILE] [-v] [-h]
  >                        trajfile amatoms
  > 
  > positional arguments:
  >   trajfile              xyz or lammpstrj file to use. Using a xyz file
  >                         will result in a non periodic box calculation.
  >                         Lammpstrj files will be treated as periodic if
  >                         not defined otherwise with --nonperiodic.
  >   amatoms               Number of atoms in the first molecule (type:
  >                         Int64)
  > 
  > optional arguments:
  >   -f, --radial_factor RADIAL_FACTOR
  >                         Factor to multiply the VdW radii of the atoms
  >                         (type: Float64, default: 0.6)
  >   -m, --maxstep MAXSTEP
  >                         Maximum step to use (default all steps) (type:
  >                         Int64)
  >   -s, --startstep STARTSTEP
  >                         Initial step to use default to first step in
  >                         file (type: Int64)
  >   -d, --densfunc DENSFUNC
  >                         Function to use for smearing. Available
  >                         options are: slater, gaus, rect (default:
  >                         "slater")
  >   -n, --nonperiodic     If set, the box will be treated as non
  >                         periodic. Only applicable for lammpstrj files
  >                         since xyz files are always non periodic.
  >   -o, --outfile OUTFILE
  >                         Output file to save the entropy values, if not
  >                         set, the values will be printed to the
  >                         console. Only applicable for lammpstrj files.
  >                         If set threads can be used with julia -t
  >                         <nthreads>.
  >   -v, --debug           Print debug information
  >   -h, --help            show this help message and exit 
```
