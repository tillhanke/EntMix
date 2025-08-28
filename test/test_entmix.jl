using Test
using EntMix
using LinearAlgebra  
using StaticArrays

function test_entmix_module()
    @testset "EntMix Module Tests" begin
        @testset "Module Loading" begin
            # Test that EntMix module can be loaded
            @test isdefined(Main, :EntMix)
            
            # Test that main module exists
            @test isa(EntMix, Module)
        end
        
        @testset "Exported Functions" begin
            # Test that expected functions are exported
            exported_names = names(EntMix)
            
            @test :entropy in exported_names
            @test :dens in exported_names  
            @test :entropy_distribution in exported_names
            @test :Molecule in exported_names
            
            # Test that exported functions are callable
            @test isa(EntMix.entropy, Function)
            @test isa(EntMix.dens, Function)
            @test isa(EntMix.entropy_distribution, Function)
            @test isa(EntMix.Molecule, Module)
        end
        
        @testset "Molecule Submodule" begin
            # Test that Molecule submodule exists and is accessible
            @test isa(EntMix.Molecule, Module)
            
            # Test that Molecule submodule has expected functions
            molecule_names = names(EntMix.Molecule, all=true)
            
            # Check for key functions (they should be defined even if not exported)
            expected_functions = [:get_molecules, :mol_types, :mol_dictionary, :type_from_name!]
            
            for func in expected_functions
                @test func in molecule_names || isdefined(EntMix.Molecule, func)
            end
        end
        
        @testset "Function Signatures" begin
            # Test that exported functions have expected signatures
            
            # entropy function should have multiple methods
            @test length(methods(EntMix.entropy)) >= 2
            
            # dens function should have multiple methods  
            @test length(methods(EntMix.dens)) >= 2
            
            # entropy_distribution should have multiple methods
            @test length(methods(EntMix.entropy_distribution)) >= 2
        end
        
        @testset "Module Structure" begin
            # Test that module includes expected source files
            # We can't directly test file inclusion, but we can test that 
            # functions from different files are available
            
            # Functions from Entropy.jl should be available
            @test isdefined(EntMix, :entropy)
            @test isdefined(EntMix, :dens)
            @test isdefined(EntMix, :entropy_distribution)
            
            # Smearing functions should be available (included via Entropy.jl)
            # These are not exported but should be accessible
            @test isdefined(EntMix, :slater) || @test_nowarn EntMix.eval(:(slater(1.0, 1.0)))
        end
        
        @testset "Basic Functionality Test" begin
            # Test that we can call basic functions with minimal parameters
            # This verifies the module loads correctly and functions are callable
            
            using StaticArrays
            using LinearAlgebra
            
            # Test basic entropy_distribution call
            r = SVector{3, Float64}([1.0, 1.0, 1.0])
            atoms_collections = [
                [SVector{3, Float64}([0.0, 0.0, 0.0])],
                [SVector{3, Float64}([2.0, 2.0, 2.0])]
            ]
            sigma_collections = [[1.0], [1.0]]
            box = SVector{3, Float64}([10.0, 10.0, 10.0])
            
            # This should not error (though result might be zero)
            @test_nowarn result = EntMix.entropy_distribution(r, atoms_collections, sigma_collections, box)
            result = EntMix.entropy_distribution(r, atoms_collections, sigma_collections, box)
            @test isfinite(result)
            @test result >= 0.0
            
            # Test basic dens call
            @test_nowarn density = EntMix.dens(r, atoms_collections[1], sigma_collections[1], box)
            density = EntMix.dens(r, atoms_collections[1], sigma_collections[1], box)
            @test isfinite(density)
            @test density >= 0.0
        end
    end
end