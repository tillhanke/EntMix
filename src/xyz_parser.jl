"""
Find the offset of the step in the xyz file 

returns:
- offset: Int - offset of the start of step in the file
"""
function find_xyz_offset(file, step)
    open(file, "r") do fio
        am_at = parse(Int, readline(fio))
        seekstart(fio)
        for i in 1:(am_at+2)*(step-1)
            readline(fio)
            if eof(fio)
                error("Step $step not found")
            end
        end
        global offset = position(fio)
    end
    return offset
end

"""
Parse single xyz timestep in form of a list of strings

returns:
- elements: Array{String} - array of elements
- atom_coords: Array{Float64,2} - array of atom coordinates
"""
function parse_xyz_list(xyzlist::Vector{String})
    amount_atoms = parse(Int, xyzlist[1])
    if length(xyzlist) != amount_atoms + 2
        error("Amount of atoms does not match the length of list")
    end

    atom_coords = zeros(amount_atoms,3)
    elements = Array{String}(undef, amount_atoms)
    for (at_ind, line) in enumerate(xyzlist[3:end])
        el, x, y, z = parse_xyz_line(line)
        elements[at_ind] = string(el)
        atom_coords[at_ind, 1:3] = [parse(Float64, x), parse(Float64, y), parse(Float64, z)]
    end
    return elements, atom_coords
end

"""
Parse single line of xyz file

returns:
- el: String - element
- x: String - x coordinate
- y: String - y coordinate
- z: String - z coordinate
"""
function parse_xyz_line(xyzline::String)
    pattern = r" *([^ ]*) *([^ ]*) *([^ ]*) *([^ ]*)" 
    return match(pattern, xyzline)
end

"""
Parse whole trajectory file, with option to define start and end step

Args:
- file: String - path to the file
- start: Int - start step
- end: Int - end step

returns:
- elements: Array{String,2} - array of elements per step
- atom_coords: Array{Float64,3} - array of atom coordinates
"""
function parse_xyz(file; start=1, stop=0)
    if start == nothing
        start = 1
    end
    offset = find_xyz_offset(file, start)
    elements = []
    atom_coords = []
    open(file, "r") do fio
        amats = parse(Int, readline(fio))
        seek(fio, offset)
        tstep = start
        if stop != 0
            @debug "Reading steps from $start to $stop"
            if stop <= start
                error("End step is before or the same as start step")
            end
        else
            @debug "Reading steps from $start to end"
        end
        while tstep != stop
            @debug "Reading step" tstep
            onestep = [readline(fio) for i in 1:amats+2]
            el, ac = parse_xyz_list(onestep)
            push!(elements, el)
            push!(atom_coords, ac)
            tstep += 1
            if eof(fio)
                break
            end
        end
    end
    return elements, atom_coords
end

