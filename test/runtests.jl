using Test
using Chemfiles
using EntMix

const Molecules = EntMix.Molecules

@testset "EntMix.jl" begin
    @testset "Molecules.geometric_center" begin
        frame = Frame()
        add_atom!(frame, Atom("H"), [0.0, 0.0, 0.0])
        add_atom!(frame, Atom("H"), [2.0, 0.0, 0.0])
        add_atom!(frame, Atom("H"), [1.0, 3.0, 0.0])
        mol = [0, 1, 2]

        gc = Molecules.geometric_center(frame, mol)
        @test length(gc) == 3
        @test gc ≈ [1.0, 1.0, 0.0]

        # a single atom: its geometric center is its own position
        single = Frame()
        add_atom!(single, Atom("C"), [4.0, -1.0, 7.0])
        @test Molecules.geometric_center(single, [0]) ≈ [4.0, -1.0, 7.0]
    end

    @testset "Molecules.center_of_mass" begin
        # equal masses -> center of mass coincides with geometric center
        frame = Frame()
        add_atom!(frame, Atom("H"), [0.0, 0.0, 0.0])
        add_atom!(frame, Atom("H"), [2.0, 0.0, 0.0])
        add_atom!(frame, Atom("H"), [1.0, 3.0, 0.0])
        mol = [0, 1, 2]

        com_equal = Molecules.center_of_mass(frame, mol)
        @test length(com_equal) == 3
        @test com_equal ≈ Molecules.geometric_center(frame, mol)

        # unequal masses -> center of mass is mass weighted, not geometric
        diatomic = Frame()
        add_atom!(diatomic, Atom("O"), [0.0, 0.0, 0.0])
        add_atom!(diatomic, Atom("H"), [10.0, 0.0, 0.0])
        mO = mass(Atom("O"))
        mH = mass(Atom("H"))

        com = Molecules.center_of_mass(diatomic, [0, 1])
        @test com[1] ≈ mH * 10.0 / (mO + mH)
        @test com[2] ≈ 0.0 atol = 1e-12
        @test com[3] ≈ 0.0 atol = 1e-12
    end
end
