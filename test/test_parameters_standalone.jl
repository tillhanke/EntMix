using Test

# Load the source file at module level to avoid world age issues
include("../src/parameters.jl")

@testset "Parameters Standalone" begin
    @testset "Constants Tests" begin
        @test BOHR ≈ 1.8897259886
        @test typeof(BOHR) == Float64
    end
    
    @testset "VDW Radii Dictionary Tests" begin
        # Test that dictionary exists and has expected type
        @test isa(VDWradii, Dict{String, Float64})
        
        # Test that common elements are present
        @test haskey(VDWradii, "H")
        @test haskey(VDWradii, "C")
        @test haskey(VDWradii, "N")
        @test haskey(VDWradii, "O")
        @test haskey(VDWradii, "F")
        @test haskey(VDWradii, "S")
        @test haskey(VDWradii, "Cl")
        @test haskey(VDWradii, "Br")
        @test haskey(VDWradii, "Si")
        
        # Test specific values from Wikipedia data
        @test VDWradii["H"] ≈ 1.2
        @test VDWradii["C"] ≈ 1.70
        @test VDWradii["N"] ≈ 1.55
        @test VDWradii["O"] ≈ 1.52
        @test VDWradii["F"] ≈ 1.47
        @test VDWradii["S"] ≈ 1.80
        @test VDWradii["Cl"] ≈ 1.75
        @test VDWradii["Br"] ≈ 1.85
        @test VDWradii["Si"] ≈ 2.10
        
        # Test that all values are positive
        for (element, radius) in VDWradii
            @test radius > 0.0
            @test isfinite(radius)
        end
        
        # Test that dictionary has expected number of entries
        @test length(VDWradii) == 9
    end
    
    @testset "Atomic Numbers Dictionary Tests" begin
        # Test that dictionary exists and has expected type
        @test isa(atomic_numbers, Dict{String, Int})
        
        # Test first period elements
        @test atomic_numbers["H"] == 1
        @test atomic_numbers["He"] == 2
        
        # Test second period elements
        @test atomic_numbers["Li"] == 3
        @test atomic_numbers["Be"] == 4
        @test atomic_numbers["B"] == 5
        @test atomic_numbers["C"] == 6
        @test atomic_numbers["N"] == 7
        @test atomic_numbers["O"] == 8
        @test atomic_numbers["F"] == 9
        @test atomic_numbers["Ne"] == 10
        
        # Test third period elements
        @test atomic_numbers["Na"] == 11
        @test atomic_numbers["Mg"] == 12
        @test atomic_numbers["Al"] == 13
        @test atomic_numbers["Si"] == 14
        @test atomic_numbers["P"] == 15
        @test atomic_numbers["S"] == 16
        @test atomic_numbers["Cl"] == 17
        @test atomic_numbers["Ar"] == 18
        
        # Test some transition metals
        @test atomic_numbers["Fe"] == 26
        @test atomic_numbers["Co"] == 27
        @test atomic_numbers["Ni"] == 28
        @test atomic_numbers["Cu"] == 29
        @test atomic_numbers["Zn"] == 30
        
        # Test some heavy elements  
        @test atomic_numbers["Au"] == 79
        @test atomic_numbers["Hg"] == 80
        @test atomic_numbers["Pb"] == 82
        @test atomic_numbers["U"] == 92
        
        # Test superheavy elements
        @test atomic_numbers["Uuo"] == 118
        
        # Test that all atomic numbers are positive
        for (element, number) in atomic_numbers
            @test number > 0
            @test number <= 118  # Currently known elements
        end
        
        # Test that dictionary has expected number of entries (118 elements)
        @test length(atomic_numbers) == 118
        
        # Test that all keys are strings and values are integers
        for (element, number) in atomic_numbers
            @test isa(element, String)
            @test isa(number, Int)
            @test length(element) >= 1
            @test length(element) <= 3  # Max 3 characters for element symbols
        end
        
        # Test some specific edge cases
        @test haskey(atomic_numbers, "I")  # Iodine
        @test haskey(atomic_numbers, "Te") # Tellurium
        @test atomic_numbers["I"] == 53
        @test atomic_numbers["Te"] == 52
    end
    
    @testset "Dictionary Consistency Tests" begin
        # Test that VDW radii elements are also in atomic numbers
        for element in keys(VDWradii)
            if element != "Si"  # Si might not be in all lists
                @test haskey(atomic_numbers, element) || element ∈ ["Si"]
            end
        end
        
        # Test common elements are in both dictionaries
        common_elements = ["H", "C", "N", "O", "F", "S", "Cl", "Br"]
        for element in common_elements
            @test haskey(VDWradii, element)
            @test haskey(atomic_numbers, element)
        end
    end
end