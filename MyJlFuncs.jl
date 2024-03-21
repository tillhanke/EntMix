using LinearAlgebra

function parse_traj(path; frames=maxtype(Int))
    println("loading ", path)
    io = open(path, "r")
    raw = read(io, String)
    close(io)
    println("I am loaded")
    offset = 0
    fsize = length(raw)
    frameid = 0
    list_of_frames = []
    while offset < fsize
        # println("Parsing frame ", frameid)
        # println("current Offset: ", offset)

        additionaloffset, framedata = parse_xyz_string(raw[offset+1:end]; offset=true)
        frameid += 1
        offset += additionaloffset
        push!(list_of_frames, framedata)
        if frameid == frames
            break
        end
    end
    println("Read through $frameid frames")
    return list_of_frames
end

function parse_xyz(path)
    println("Parsing ", path)
    io = open(path, "r")
    raw = read(io, String)
    close(io)
    println("File opened")
    println("Filesize: ", length(raw))
    return parse_xyz_string(raw)
end

function parse_xyz_string(inp; offset=false)
    amount_atoms = parse(Int, match(r"( *)?([0-9]*)", inp)[2])
    # println("Amount of atoms: ", amount_atoms)
    off = match(r"\n", inp).offset +1
    off = match(r"\n", inp, off).offset +1
    atom_coords = zeros(amount_atoms,3)
    elements = Array{String}(undef, amount_atoms)
    for at_ind in 1:amount_atoms
        el, x, y, z = (match(r"([^ ]*) *(-?\d*\.\d*) *(-?\d*\.\d*) *(-?\d*\.\d*)", inp, off)[i] for i in 1:4)
        elements[at_ind] = string(el)
        atom_coords[at_ind, 1:3] = [parse(Float64, x), parse(Float64, y), parse(Float64, z)]
        off = match(r"\n", inp, off).offset+1
    end
    if offset
        return off, (elements, atom_coords)
    end
    return elements, atom_coords
end

function constructGrid(origin, maxval; dist=0, dsize=0)
    maxval = maxval .- origin
    @assert xor(dsize != 0 , dist != 0) "Either size or dist must be given"
    if dist == 0
        if size(dsize) == ()
            dsize = [dsize, dsize, dsize]
        end
        dist = (maxval./dsize)
    else
        if size(dist) == ()
            dist = [dist, dist, dist]
        end
        dsize = @. Int(ceil(maxval/dist))
    end

    r = Matrix{Float64}(undef, prod(dsize), 3)
    for x = 0:dsize[1]-1
        for y = 0:dsize[2]-1
            for z = 0:dsize[3]-1
                r[dsize[1]*dsize[2]*z+dsize[1]*y+x+1,1:3] = origin .+ ([x,y,z] .* dist)
            end
        end
    end
    return r
end

function triangle(r, r0, sigma)
    height = 1/(sigma^3*(4*pi/3-pi))
    dist = broadcast(-,r,transpose(r0)) 
    p = height.- sum(dist.*dist, dims=2).^(1/2).*height/sigma
    for i in 1:length(p)
        if p[i] < 0
            p[i] = 0
        end
    end
    return p
end

function slater(r, r0, sigma; n=1)
    norm = 1/(factorial(n+1)*sigma^(n+2)*4*pi)
    dist = broadcast(-,r,transpose(r0)) 
    dist = sum(dist.*dist, dims=2).^(1/2)
    if size(dist) == (1,1)
        dist = dist[1]
    end
    return @. norm* dist^(n-1)*exp(-dist/sigma)
end

function rect(r, r0, sigma)
    d = broadcast(-,r,transpose(r0)) 
    dist = sum(d.*d, dims=2)
    inds = findall(x -> x < sigma, dist)
    dens = zeros(size(r)[1])
    dens[inds] .= 1
    if size(dens) == (1,)
        dens = dens[1]
    end
    return dens
end

function gaus(r, r0, sigma)
    dist = broadcast(-,r,transpose(r0))   
    g = 1/(sigma^3*(pi*2)^(3/2))*exp.(-sum(dist .* dist, dims=2) / (2*sigma^2))
    if size(g) == (1,1)
        g = g[1]
    end
    return g
end

# --------------------------
# Old Functions
# --------------------------

function calc_dens(r, atoms, elements, radii, dens_func)
    dens = zeros(size(r)[1])
    for at_ind in 1:size(atoms)[1]
        dens += dens_func(r, atoms[at_ind, 1:3], radii[elements[at_ind]])
    end
    return dens
end

function entropy_from_dens(dens1, dens2, r)
    both = dens1 .+ dens2
    inds1 = findall(!iszero, dens1)
    inds2 = findall(!iszero, dens2)
    frac1 = dens1[inds1]./both[inds1]
    frac2 = dens2[inds2]./both[inds2]
    entrop = sum(frac1.*log.(frac1))/size(r)[1] + sum(frac2.*log.(frac2))/size(r)[1]
    return entrop
end

function entropy_from_dens(dens, dv; norm=false)
    if norm
        norm = 1/(sum(dens)*dv)
    else 
        norm = 1
    end
    ids = findall(!iszero, dens)
    return sum(@. norm* dens[ids]*log(2, norm* dens[ids]))*dv
end

function entropy(r, atom1, atom2, elem1, elem2, radii; debug=false, dens_func=gaus)
    println("Calculating density of Typ 1:")
    dens1 = calc_dens(r, atom1, elem1, radii, dens_func)
    println("Calculating density of Typ 2:")
    dens2 = calc_dens(r, atom1, elem1, radii, dens_func)
    entrop = entropy_from_dens(dens1, dens2, r)
    if debug
        return entrop, dens1, dens2
    else
        return entrop
    end
end


