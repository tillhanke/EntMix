"""
Parse a lammpstrj file with elements and 3 coordinates per atom
#### Args:
- filename: String  path to lammpstrj file
- start: Int  timestep to start reading from (default 0)
- stop: Int  timestep to stop reading at
#### Returns:
- elements: Array{String, 2}  elements[i, j] is the element of atom j in timestep i
- coords: Array{Float64, 3}  coords[i, j, k] is the kth coordinate of atom j in timestep i
- boxs: Array{Float64, 3}  boxs[i, j, k] is boxbounds in timestep i 
"""
function parse_lammpstrj(filename::String; start=0, stop=Inf, ret_steps=false)
    # read head of file (9 lines)
    head = Dict()
    open(filename) do file
        head = read_lammpstrj_head(file)
        firststep = head["timestep"]
        @debug "Head contains: $(head)"
        if start == 0
            start = firststep
        end
        if head["timestep"] > start
            error("Start timestep is before first timestep in file")
        end
        @debug "Skipping to head of 2nd timestep"
        for l in 1:(head["n_atoms"]+9)
            readline(file)
        end
        if !eof(file)
            @debug "getting deltats"
            head2 = read_lammpstrj_head(file)
            # delta ts per timestep written to lammpstrj file
            delts = head2["timestep"] - head["timestep"]
        # approx amount of bytes per timestep
            ts_byte = position(file) 
            @debug delts, ts_byte
        else
            @info "Only found one timestep in file"
            delts = 1
            ts_byte = 0
            stop = firststep
        end
        if stop != Inf
            @debug "Stopping at timestep $(stop)"
        else
            stop, last_pos = seekend_lammpstrj(file)
            ts_byte = last_pos÷((stop-firststep)÷delts-1)
            @debug "Set tsbyte by approx from final step to $(ts_byte)"
        end
        elements = Array{String}(undef, (stop-start)÷delts+1,  head["n_atoms"])
        coords = Array{Float64}(undef, (stop-start)÷delts+1,  head["n_atoms"], 3)
        boxs = Array{Float64}(undef, (stop-start)÷ delts+1,  2, 3)
        @debug "initialized arrays:", size(elements), size(coords), size(boxs)

        # skip over timesteps until start
        seekstart(file)
        if start > head["timestep"]
            start_byte = find_lammpstrj_timestep(file, start, approx_byte=ts_byte, delts=delts)
            @debug "found starting timestep at byte $(start_byte)"
        end
        @debug begin
            head = read_lammpstrj_head(file) 
            "Starting at timestep $(head["timestep"])"
        end
        @debug "file position: $(position(file))"

        # read data
        @debug "Reading total amount of $((stop-start÷delts)+1) timesteps"
        starttime = time()
        if !("x" in keys(head["columns"]) )
            if "xu" in keys(head["columns"])
                @info "loading unwrapped coordinates"
            else
                error("No x column found in timestep")
            end
        end


        for tstep in 1:(stop-start)÷delts+1
            if tstep%1000 == 0
                @debug "Reading $(tstep)th at time: $(time() - starttime)"
            end
            head = read_lammpstrj_head!(file)
            if !("x" in keys(head["columns"]) )
                if "xu" in keys(head["columns"])
                    head["columns"]["x"] = head["columns"]["xu"]
                    head["columns"]["y"] = head["columns"]["yu"]
                    head["columns"]["z"] = head["columns"]["zu"]
                end
            end
            xc, yc, zc = head["columns"]["x"], head["columns"]["y"], head["columns"]["z"]
            elc = head["columns"]["element"]
            for l in 1:head["n_atoms"]
                line = readline(file)
                line = split(line)  # split into array of Strings
                elements[tstep, l] = line[elc]
                coords[tstep, l, :] = [parse(Float64, x) for x in line[xc:zc]]
            end
            boxs[tstep, :, :] = head["box_bounds"]
        end
        if ret_steps
            return elements, coords, boxs, [i for i in start:delts:stop] 
        end
        return elements, coords, boxs
    end
end

"""
Search for the last timestep and return the timestep number
#### Args:
- file: IOStream of a lammpstrj file at any position 
#### Returns:
- Int  timestep number of last timestep
- Int  position of last timestep in file
"""
function seekend_lammpstrj(file::IOStream)
    startpos = position(file)
    seekend(file)
    @debug "Searching for final timestep"
    laststep = 0
    last_pos = 0
    s = ""
    while true
        skip(file, -2)
        c = read(file, Char)
        s = string(c, s)
        if c == '\n'
            s = ""
        end
        if endswith(s, "STEP")
            readline(file)
            last_pos = position(file) - 15  # position before "ITEM: TIMESTEP"
            laststep = parse(Int, readline(file))
            @debug "Final Timestep:", laststep
            break
        end
    end
    seek(file, startpos)
    return laststep, last_pos   
end

"""
Read the head of a lammpstrj file timestep
#### Args:
- file: IOStream at position of head start
- offset: Int  offset from current position to start of head
Returns:
- Dict with keys:
  + timestep: Int
  + n_atoms: Int
  + box_bounds: Array{Float64, 2}  [xmin ymin zmin; xmax ymax zmax]
"""
function read_lammpstrj_head(file::IOStream; offset=0)
    pos = position(file)
    head = read_lammpstrj_head!(file; offset=offset)
    seek(file, pos)
    return head
end

"""
Read the head of a lammpstrj file timestep and update the file position
#### Args:
- file: IOStream at position of head start
- offset: Int  offset from current position to start of head
Returns:
- Dict with keys:
  + timestep: Int
  + n_atoms: Int
  + box_bounds: Array{Float64, 2}  [xmin ymin zmin; xmax ymax zmax]
"""
function read_lammpstrj_head!(file::IOStream; offset=0)
    head = Dict()
    skip(file, offset)
    readline(file) # ITEM: TIMESTEP
    head["timestep"] = parse(Int, readline(file))
    readline(file) # ITEM: NUMBER OF ATOMS
    head["n_atoms"] = parse(Int, readline(file))
    readline(file) # ITEM: BOX BOUNDS
    head["box_bounds"] = [parse(Float64, x) for x in split(readline(file))]
    head["box_bounds"] = [head["box_bounds"] [parse(Float64, x) for x in split(readline(file))]]
    head["box_bounds"] = [head["box_bounds"] [parse(Float64, x) for x in split(readline(file))]]
    cols = split(readline(file))[3:end] # ITEM: ATOMS id type x y z
    head["columns"] = Dict()
    for (i, c) in enumerate(cols)
        head["columns"][c] = i
    end
    return head
end

"""
Searches for start of defined Timestep in lammpstrj file
#### Args:
- file: IOStream at position of head start
- timestep: Int  timestep to search for
- approx_byte: Int  approximate bytes per timestep 
- delts: Int  delta timestep per timestep written to lammpstrj file
Returns:
- position of timestep in file
"""
function find_lammpstrj_timestep(file::IOStream, timestep::Int; approx_byte=nothing, delts=nothing)
    pos = position(file)
    if approx_byte == nothing || delts == nothing
        @debug "Reading head of first two timesteps to get approx_byte and delts"
        head = read_lammpstrj_head!(file)
        for l in 1:head["n_atoms"]
            readline(file)
        end
        if delts == nothing
            @debug "Reading head of second timestep to get delts"
            delts = read_lammpstrj_head(file)["timestep"] - head["timestep"]
            @debug delts
        end
        if approx_byte == nothing
            approx_byte = position(file) - pos  # approx bytes per timestep
            @debug approx_byte 
        end
    end
    curr_head = read_lammpstrj_head(file)
    curr_ts = curr_head["timestep"]
    skip_bytes = approx_byte*((timestep - curr_ts)÷delts-5)
    seek_pos = position(file) + skip_bytes
    if seek_pos < 0
        seek_pos = 0
    end
    seek(file, seek_pos)  # start approx 5 steps before position
    @debug "Searching for timestep $(timestep)"
    while !eof(file)
        l = readline(file)
        if endswith(l, "STEP")
            pos = position(file)
            step = parse(Int, readline(file))
            if step > timestep
                @debug "Approximation too far at $(step), seeking back"
                tofar = approx_byte*((step-timestep)÷delts+1)
                seek(file, pos - tofar)
                continue
            end
            if step == timestep
                seek(file, pos-length(l)-1)
                return pos - length(l) - 1
            end
        end
    end
    error("Did not find timestep in file")
end






