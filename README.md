# Entropy of Mixing Calculation based on MD Trajectory

This module evaluates trajectory data from md simulations and calculates the configurational entropy distribution.

Currently only binary mixtures are supported.

## Installation
### Install Julia
It is recommendet to use Juliaup to install julia and fix the default version to 1.12.5 (which is the version, the module was developed with).

Juliaup can be found [here](https://github.com/JuliaLang/juliaup)

### Install Module
Clone the repository and install the julia module:
```bash
git clone [repo]
julia
```
then install the module from the directory by entering the Pkg manager with `]`
```julia
dev "./EntMix"
```
this adds the directory as module directory and makes it editable. So a update of the module can be done by simply pulling the repository.

## Usage
The simplest way to use the functionality of this package is to make use of the `calc_entropy.jl` script in the [EntMixScripts](https://github.com/tillhanke/EntMixScripts/tree/master/src) repository.

### Run Tests:
```bash 
julia --project=./ -e "using Pkg;Pkg.test()"
```

# Related publications:
- https://doi.org/10.1021/acs.jpclett.4c02819
- https://doi.org/10.1063/5.0304658

