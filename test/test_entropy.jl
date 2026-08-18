using Test
using EntMix
using Chemfiles
using LinearAlgebra: norm
using StaticArrays

const pbc_distance = EntMix.pbc_distance
const pbc_distance! = EntMix.pbc_distance!
const wrap = EntMix.wrap
const dens = EntMix.dens
const entropy_distribution = EntMix.entropy_distribution

@testset "wrap" begin
    box = [10.0, 10.0, 10.0]
    # points already inside the box are unchanged
    @test wrap([1.0, 2.0, 3.0], box) ≈ [1.0, 2.0, 3.0]
    # points outside are folded back in
    @test wrap([11.0, -1.0, 25.0], box) ≈ [1.0, 9.0, 5.0]
    # result is always within [0, box)
    w = wrap([-0.001, 123.4, -55.5], box)
    @test all(0 .<= w .< box)
end

@testset "pbc_distance" begin
    box = [10.0, 10.0, 10.0]

    @testset "minimum image convention" begin
        # atoms near opposite faces are actually close across the boundary
        d = pbc_distance([0.0, 0.0, 0.0], [9.0, 0.0, 0.0], box)
        @test norm(d) ≈ 1.0
        # within the box the distance is the plain difference
        d2 = pbc_distance([1.0, 1.0, 1.0], [2.0, 3.0, 1.0], box)
        @test d2 ≈ [-1.0, -2.0, 0.0]
        # each component stays within [-box/2, box/2)
        d3 = pbc_distance([0.0, 0.0, 0.0], [4.9, 5.1, 7.5], box)
        @test all(-box ./ 2 .<= d3 .< box ./ 2)
    end

    @testset "in-place variant matches allocating variant" begin
        r1 = [0.3, 8.2, 4.0]
        r2 = [9.5, 1.1, 0.2]
        out = zeros(3)
        pbc_distance!(r1, r2, box, out)
        @test out ≈ pbc_distance(r1, r2, box)
    end

    @testset "frame-based methods" begin
        frame = Frame()
        cell = UnitCell([10.0, 10.0, 10.0])
        set_cell!(frame, cell)
        add_atom!(frame, Atom("H"), [0.0, 0.0, 0.0])
        add_atom!(frame, Atom("H"), [9.0, 0.0, 0.0])
        # distance from a point to an atom (chemfiles 0-based index)
        @test norm(pbc_distance([0.0, 0.0, 0.0], frame, 1)) ≈ 1.0
        # distance between two atoms
        @test norm(pbc_distance(0, 1, frame)) ≈ 1.0
    end
end

@testset "dens" begin
    box = [100.0, 100.0, 100.0]  # large box so PBC images don't interfere

    @testset "single atom matches the kernel" begin
        pos = [SVector{3,Float64}(50.0, 50.0, 50.0)]
        sigma = [1.0]
        r = [50.0, 50.0, 50.0]
        @test dens(r, pos, sigma, box; dfunc=EntMix.gaus) ≈ EntMix.gaus(0.0, 1.0)
        # a point one unit away
        r2 = [51.0, 50.0, 50.0]
        @test dens(r2, pos, sigma, box; dfunc=EntMix.gaus) ≈ EntMix.gaus(1.0, 1.0)
    end

    @testset "density is additive over atoms" begin
        positions = [SVector{3,Float64}(40.0, 50.0, 50.0),
                     SVector{3,Float64}(60.0, 50.0, 50.0)]
        sigma = [1.0, 1.0]
        r = [50.0, 50.0, 50.0]
        single = dens(r, [positions[1]], [1.0], box; dfunc=EntMix.gaus)
        both = dens(r, positions, sigma, box; dfunc=EntMix.gaus)
        @test both ≈ 2 * single
    end

    @testset "frame method agrees with vector method" begin
        frame = Frame()
        set_cell!(frame, UnitCell(box))
        add_atom!(frame, Atom("C"), [50.0, 50.0, 50.0])
        atoms = [0]
        r = [50.5, 50.0, 50.0]
        d_frame = dens(r, frame, atoms, 1.0; dfunc=EntMix.gaus)
        pos = [SVector{3,Float64}(50.0, 50.0, 50.0)]
        scaled = [1.0 * covalent_radius(frame[0])]
        d_vec = dens(r, pos, scaled, SVector{3}(box); dfunc=EntMix.gaus)
        @test d_frame ≈ d_vec
    end
end

@testset "entropy_distribution" begin
    box = [100.0, 100.0, 100.0]
    r = [50.0, 50.0, 50.0]

    @testset "perfectly mixed point gives maximal entropy" begin
        # two species with identical single atoms at the query point
        atom = SVector{3,Float64}(50.0, 50.0, 50.0)
        collections = [[atom], [atom]]
        sigmas = [[1.0], [1.0]]
        s = entropy_distribution(r, collections, sigmas, box; dfunc=EntMix.gaus)
        @test s ≈ log(2)

        # three identical species -> log(3)
        collections3 = [[atom], [atom], [atom]]
        sigmas3 = [[1.0], [1.0], [1.0]]
        s3 = entropy_distribution(r, collections3, sigmas3, box; dfunc=EntMix.gaus)
        @test s3 ≈ log(3)
    end

    @testset "single species gives zero mixing entropy" begin
        atom = SVector{3,Float64}(50.0, 50.0, 50.0)
        s = entropy_distribution(r, [[atom]], [[1.0]], box; dfunc=EntMix.gaus)
        @test s == 0
    end

    @testset "empty region (compact kernels) gives zero" begin
        far = SVector{3,Float64}(10.0, 10.0, 10.0)
        collections = [[far], [far]]
        sigmas = [[1.0], [1.0]]
        # query point is well outside every compact-support kernel
        s = entropy_distribution(r, collections, sigmas, box; dfunc=EntMix.constant)
        @test s == 0
    end

    @testset "entropy is bounded by log(number of species)" begin
        a = SVector{3,Float64}(50.0, 50.0, 50.0)
        b = SVector{3,Float64}(51.5, 50.0, 50.0)  # slightly displaced second species
        collections = [[a], [b]]
        sigmas = [[1.0], [1.0]]
        s = entropy_distribution(r, collections, sigmas, box; dfunc=EntMix.gaus)
        @test 0 <= s <= log(2) + 1e-12
    end

    @testset "frame method agrees with vector method" begin
        frame = Frame()
        set_cell!(frame, UnitCell(box))
        add_atom!(frame, Atom("C"), [50.0, 50.0, 50.0])
        add_atom!(frame, Atom("C"), [50.0, 50.0, 50.0])
        collections = [[0], [1]]
        s_frame = entropy_distribution(r, frame, collections, 1.0; dfunc=EntMix.gaus)
        @test s_frame ≈ log(2)
    end
end
