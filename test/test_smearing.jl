using Test
using LinearAlgebra
using StaticArrays

function test_smearing_functions()
    # Include the source files here within the function
    include("../src/Smearing.jl")
    @testset "Slater Function Tests" begin
        # Test 1D slater function
        @testset "slater(r, sigma)" begin
            # Test basic functionality
            result = slater(1.0, 2.0)
            @test result > 0.0  # Should be positive
            @test isfinite(result)
            
            # Test n parameter
            result_n1 = slater(1.0, 2.0, n=1)
            result_n2 = slater(1.0, 2.0, n=2) 
            @test result_n1 != result_n2
            
            # Test symmetry (distance should be abs(r))
            @test slater(1.0, 2.0) == slater(-1.0, 2.0)
            
            # Test at zero distance
            result_zero = slater(0.0, 1.0)
            # For n=1, r^(n-1) = r^0 = 1, so function is non-zero at r=0
            @test result_zero > 0.0  
            
            # Test sigma scaling
            result_sigma1 = slater(1.0, 1.0)
            result_sigma2 = slater(1.0, 2.0)
            @test result_sigma1 != result_sigma2
        end
        
        # Test 3D slater function  
        @testset "slater(r, r0, sigma)" begin
            r = [1.0, 2.0, 3.0]
            r0 = [0.0, 0.0, 0.0]
            sigma = 1.0
            
            result = slater(r, r0, sigma)
            @test result > 0.0
            @test isfinite(result)
            
            # Test with same position
            result_same = slater(r0, r0, sigma)
            # Distance is 0, but for n=1, r^(n-1) = r^0 = 1, so non-zero
            @test result_same > 0.0
            
            # Test symmetry
            result1 = slater(r, r0, sigma)  
            result2 = slater(r0, r, sigma)
            @test result1 == result2
            
            # Test with different positions
            r1 = [1.0, 0.0, 0.0]
            r2 = [2.0, 0.0, 0.0]
            result_close = slater(r1, r0, sigma)
            result_far = slater(r2, r0, sigma)
            @test result_close > result_far  # Closer should have higher density
        end
    end
    
    @testset "Gaussian Function Tests" begin
        # Test 1D gaussian function
        @testset "gaus(r, sigma)" begin
            result = gaus(1.0, 2.0)
            @test result > 0.0
            @test isfinite(result)
            
            # Test symmetry
            @test gaus(1.0, 2.0) == gaus(-1.0, 2.0)
            
            # Test at zero
            result_zero = gaus(0.0, 1.0)
            @test result_zero > 0.0  # Gaussian is positive everywhere
            
            # Test maximum at r=0
            result_zero = gaus(0.0, 1.0)
            result_nonzero = gaus(1.0, 1.0) 
            @test result_zero > result_nonzero
        end
        
        # Test 3D gaussian function
        @testset "gaus(r, r0, sigma)" begin
            r = [1.0, 2.0, 3.0]
            r0 = [0.0, 0.0, 0.0]
            sigma = 1.0
            
            result = gaus(r, r0, sigma)
            @test result > 0.0
            @test isfinite(result)
            
            # Test maximum at same position
            result_same = gaus(r0, r0, sigma)
            result_diff = gaus(r, r0, sigma)
            @test result_same > result_diff
            
            # Test symmetry
            result1 = gaus(r, r0, sigma)
            result2 = gaus(r0, r, sigma)
            @test result1 == result2
        end
    end
    
    @testset "Cauchy Function Tests" begin
        # Test 1D cauchy function
        @testset "cauchy(r, sigma)" begin
            result = cauchy(1.0, 2.0)
            @test result > 0.0
            @test isfinite(result)
            
            # Test symmetry
            @test cauchy(1.0, 2.0) == cauchy(-1.0, 2.0)
            
            # Test at zero
            result_zero = cauchy(0.0, 1.0)
            @test result_zero > 0.0
            
            # Test maximum at r=0
            result_zero = cauchy(0.0, 1.0)
            result_nonzero = cauchy(1.0, 1.0)
            @test result_zero > result_nonzero
        end
        
        # Test 3D cauchy function
        @testset "cauchy(r, r0, sigma)" begin
            r = [1.0, 2.0, 3.0] 
            r0 = [0.0, 0.0, 0.0]
            sigma = 1.0
            
            result = cauchy(r, r0, sigma)
            @test result > 0.0
            @test isfinite(result)
            
            # Test symmetry
            result1 = cauchy(r, r0, sigma)
            result2 = cauchy(r0, r, sigma)
            @test result1 == result2
            
            # Test maximum at same position
            result_same = cauchy(r0, r0, sigma)
            result_diff = cauchy(r, r0, sigma)
            @test result_same > result_diff
        end
    end
    
    @testset "Algebraic Function Tests" begin
        # Test 1D algebraic function
        @testset "algebraic(r, sigma)" begin
            result = algebraic(1.0, 2.0)
            @test result > 0.0
            @test isfinite(result)
            
            # Test symmetry
            @test algebraic(1.0, 2.0) == algebraic(-1.0, 2.0)
            
            # Test alpha parameter
            result_alpha3 = algebraic(1.0, 2.0, alpha=3)
            result_alpha4 = algebraic(1.0, 2.0, alpha=4)
            @test result_alpha3 != result_alpha4
            
            # Test at zero
            result_zero = algebraic(0.0, 1.0)
            @test result_zero > 0.0
        end
        
        # Test 3D algebraic function  
        @testset "algebraic(r, r0, sigma)" begin
            r = [1.0, 2.0, 3.0]
            r0 = [0.0, 0.0, 0.0]
            sigma = 1.0
            
            result = algebraic(r, r0, sigma)
            @test result > 0.0
            @test isfinite(result)
            
            # Test symmetry
            result1 = algebraic(r, r0, sigma)
            result2 = algebraic(r0, r, sigma)
            @test result1 == result2
            
            # Test alpha parameter
            result_alpha3 = algebraic(r, r0, sigma, alpha=3)
            result_alpha4 = algebraic(r, r0, sigma, alpha=4)
            @test result_alpha3 != result_alpha4
            
            # Test maximum at same position
            result_same = algebraic(r0, r0, sigma)
            result_diff = algebraic(r, r0, sigma)
            @test result_same > result_diff
        end
    end
end