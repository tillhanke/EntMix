#!/usr/bin/env julia
include("Functions.jl")
include("lammpstrj_parser.jl")
include("xyz_parser.jl")

using Base.Threads
using DelimitedFiles

"""
Calculate configurational Entropy of mixing for a given xyz file
Printing out all entropies for each step from startstep to maxstep

Args:
* file: String - path to the xyz file containing both molecules in each step
* n_atoms: Int - number of atoms in the first molecule
* maxstep: Int - maximum step to calculate entropy for
* startstep: Int - step to start calculation from
* kwargs: Dict - additional arguments for the entropy function
"""
function xyz_entropy(file, n_atoms, maxstep, startstep; kwargs...)
    if maxstep == nothing
        maxstep = -1  # default for xyz
    end
    if startstep == nothing
        startstep = 1  # default for xyz
    end
    elements, positions = parse_xyz(file, start=startstep, stop=maxstep)
    maxstep = length(elements)
    @info "Not using periodic boundary conditions"
    println("Step\tEntropy")
    for (el, coo, step) in zip(elements, positions, startstep:maxstep)
        println(
                step, 
                "\t", 
                entropy(
                    (coo[1:n_atoms,:], el[1:n_atoms]),
                    (coo[n_atoms+1:end,:], el[n_atoms+1:end]);
                    kwargs...
                    )
               )
    end
end

"""
Calculate configurational Entropy of mixing for a given lammpstrj file
Printing out all entropies for each step from startstep to maxstep or saving to a file

Args:
* file: String - path to the xyz file containing both molecules in each step
* n_atoms: Int - number of atoms in the first molecule
* maxstep: Int - maximum step to calculate entropy for
* startstep: Int - step to start calculation from
* outfile: String - path to the file to save the entropies to
* kwargs: Dict - additional arguments for the entropy function
"""
function lammpstrj_entropy(file, n_atoms, maxstep=nothing, startstep=nothing; outfile=nothing, kwargs...)
    if maxstep == nothing 
        maxstep = Inf  # default for lammpstrj
    end
    if startstep == nothing
        startstep = 0  # default for lammpstrj
    end
    elements, positions, boxes, steps = parse_lammpstrj(file; start=startstep, stop=maxstep, ret_steps=true)
    
    starttime = time()
    if outfile != nothing
        outarray = Array{Float64}(undef, size(elements)[1], 2)
        @info "Using $(Threads.nthreads()) threads"
        # @threads for (el, coo, box, step) in zip(
        #                                           eachslice(elements, dims=1),
        #                                           eachslice(positions, dims=1),
        #                                           eachslice(boxes, dims=1),
        #                                           [i for i in startstep:skippedsteps:maxstep]
        #                                          )
        @threads for i in 1:size(elements)[1]
            outarray[i, 1] = steps[i]
            outarray[i, 2] = entropy(
                (positions[i,1:n_atoms,:], elements[i,1:n_atoms]),
                (positions[i,n_atoms+1:end,:], elements[i,n_atoms+1:end]);
                periodic=true,
                box=boxes[i,:,:],
                kwargs...
                )
        end
        writedlm(outfile, outarray)
        @info "Time for $(size(steps)[1]) steps taken: ", time()-starttime
        return outarray
    else
        println("# Step\tEntropy")

        for (el, coo, box, step) in zip(
                                      eachslice(elements, dims=1),
                                      eachslice(positions, dims=1),
                                      eachslice(boxes, dims=1),
                                      steps
                                     )
            println(
                    step,
                    "\t",
                    entropy(
                        (coo[1:n_atoms,:], el[1:n_atoms]),
                        (coo[n_atoms+1:end,:], el[n_atoms+1:end]);
                        periodic=true,
                        box=box,
                        kwargs...
                        )
                   )
        end
    end
    @info "Time for $(size(steps)[1]) steps taken: ", time()-starttime
end


