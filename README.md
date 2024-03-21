# Mixing Entropy Calculation based on MD Trajectory

This script calculates mixing entropy from MD trajectories in xyz format.

Currently only binary mixtures are supported.
Usage example:
```bash
julia CalcEntropy.jl -m 100 traj.xyz 300
```
This will calculate the mixing entropy of the first 100 timesteps within `traj.xyz`. 
The two molecules are sorted within each timestep and the first 300 atoms are all part of the first molecule type.

A complete usage is provided via
```bash
julia CalcEntropy.jl --help
  > usage: CalcEntropy.jl [-f RADIAL_FACTOR] [-m MAXSTEP] [-h] xyzfile
  >                      amatoms
  >
  > positional arguments:
  >   xyzfile               xyz file to use
  >   amatoms               Number of atoms in the first molecule (type:
  >                         Int64)
  > 
  > optional arguments:
  >   -f, --radial_factor RADIAL_FACTOR
  >                         Factor to multiply the VdW radii of the atoms
  >                         (default .6) (type: Float64, default: 1.0)
  >   -m, --maxstep MAXSTEP
  >                         Maximum number of steps to use (type: Int64,
  >                         default: 9223372036854775807)
  >   -h, --help            show this help message and exit
  > 
```
