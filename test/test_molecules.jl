using Test
using Chemfiles
using EntMix

const Molecules = EntMix.Molecules

"""
Build a frame with two water molecules (O H H, O H H) and explicit bonds so that
molecule detection does not rely on `guess_bonds!`.
"""
function two_waters()
    frame = Frame()
    set_cell!(frame, UnitCell([20.0, 20.0, 20.0]))
    add_atom!(frame, Atom("O"), [0.0, 0.0, 0.0])
    add_atom!(frame, Atom("H"), [0.96, 0.0, 0.0])
    add_atom!(frame, Atom("H"), [-0.24, 0.93, 0.0])
    add_atom!(frame, Atom("O"), [10.0, 10.0, 10.0])
    add_atom!(frame, Atom("H"), [10.96, 10.0, 10.0])
    add_atom!(frame, Atom("H"), [9.76, 10.93, 10.0])
    Chemfiles.add_bond!(frame, 0, 1)
    Chemfiles.add_bond!(frame, 0, 2)
    Chemfiles.add_bond!(frame, 3, 4)
    Chemfiles.add_bond!(frame, 3, 5)
    return frame
end

@testset "Molecules detection" begin
    @testset "adjacent builds the bond graph" begin
        frame = two_waters()
        adj = Molecules.adjacent(frame)
        @test length(adj) == 6
        @test Set(adj[1]) == Set([1, 2])   # O of water 1 bonded to its two H
        @test adj[2] == [0]                 # H bonded back to its O
        @test adj[3] == [0]
        @test Set(adj[4]) == Set([4, 5])    # O of water 2
        @test adj[5] == [3]
        @test adj[6] == [3]
    end

    @testset "get_molecules separates connected components" begin
        frame = two_waters()
        mols = Molecules.get_molecules(frame)
        @test length(mols) == 2
        # each molecule has three atoms; together they cover all six atoms
        @test all(length(m) == 3 for m in mols)
        @test Set(vcat(mols...)) == Set(0:5)
        # the two oxygens end up in different molecules
        component = Dict(a => i for (i, m) in enumerate(mols) for a in m)
        @test component[0] != component[3]
        @test component[0] == component[1] == component[2]
    end

    @testset "mol_dictionary maps atoms to molecule indices" begin
        molecules = [[0, 1, 2], [3, 4, 5]]
        d = Molecules.mol_dictionary(molecules)
        @test d[0] == 1
        @test d[1] == 1
        @test d[2] == 1
        @test d[3] == 2
        @test d[5] == 2
        @test length(d) == 6
    end
end

@testset "Molecules.type_from_name!" begin
    @testset "explicit type list" begin
        frame = Frame()
        add_atom!(frame, Atom("A"), [0.0, 0.0, 0.0])
        add_atom!(frame, Atom("B"), [1.0, 0.0, 0.0])
        Molecules.type_from_name!(frame, ["Na", "Cl"])
        @test Chemfiles.type(frame[0]) == "Na"
        @test Chemfiles.type(frame[1]) == "Cl"
    end

    @testset "copy names into types with force" begin
        frame = Frame()
        add_atom!(frame, Atom("C"), [0.0, 0.0, 0.0])
        # deliberately overwrite the type, then restore it from the name
        Chemfiles.set_type!(view(frame, 0), "wrong")
        @test Chemfiles.type(frame[0]) == "wrong"
        Molecules.type_from_name!(frame; force=true)
        @test Chemfiles.type(frame[0]) == name(frame[0])
    end
end

@testset "Molecules.smiles" begin
    function molecule(symbols, bonds)
        frame = Frame()
        for s in symbols
            add_atom!(frame, Atom(s), [0.0, 0.0, 0.0])
        end
        for (i, j) in bonds
            Chemfiles.add_bond!(frame, i, j)
        end
        return frame
    end

    @testset "water reduces to its heavy atom" begin
        frame = molecule(["O", "H", "H"], [(0, 1), (0, 2)])
        @test Molecules.smiles(frame, [0, 1, 2]) == "O"
    end

    @testset "methane reduces to its single carbon" begin
        frame = molecule(["C", "H", "H", "H", "H"], [(0, 1), (0, 2), (0, 3), (0, 4)])
        @test Molecules.smiles(frame, [0, 1, 2, 3, 4]) == "C"
    end

    @testset "ethane keeps both carbons" begin
        frame = molecule(["C", "C", "H", "H", "H", "H", "H", "H"],
                         [(0, 1), (0, 2), (0, 3), (0, 4), (1, 5), (1, 6), (1, 7)])
        @test Molecules.smiles(frame, collect(0:7)) == "CC"
    end
end

@testset "Molecules.mol_types groups identical molecules" begin
    frame = two_waters()
    molecules = [[0, 1, 2], [3, 4, 5]]
    types = Molecules.mol_types(frame, molecules)
    # both molecules are water -> a single species containing both molecule ids
    @test length(types) == 1
    key = first(keys(types))
    @test Set(types[key]) == Set([1, 2])
end

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
