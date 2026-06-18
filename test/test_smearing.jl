using Test
using EntMix
using HCubature: hcubature
using LinearAlgebra: norm

# Density kernels live in the EntMix namespace (included via Entropy.jl)
const slater = EntMix.slater
const gaus = EntMix.gaus
const cauchy = EntMix.cauchy
const linear = EntMix.linear
const epanechnikov = EntMix.epanechnikov
const cubicspline = EntMix.cubicspline
const wendland = EntMix.wendland
const gengauss = EntMix.gengauss
const lorentz = EntMix.lorentz
const constant = EntMix.constant
const cosine = EntMix.cosine

"""
Integrate a spherically symmetric kernel f(r, sigma) over all of R^3 via the
radial integral ∫₀^∞ 4π r² f(r) dr. For compact-support kernels pass the cutoff,
otherwise the half-open domain is mapped to [0, 1) with r = t/(1-t).
"""
function radial_norm(f, sigma; cutoff=nothing, rtol=1e-6)
    if cutoff !== nothing
        val, _ = hcubature(x -> 4π * x[1]^2 * f(x[1], sigma), (0.0,), (cutoff,); rtol=rtol)
        return val
    end
    val, _ = hcubature((0.0,), (1.0,); rtol=rtol) do x
        t = x[1]
        t >= 1 && return 0.0
        r = t / (1 - t)
        jac = 1 / (1 - t)^2
        return 4π * r^2 * f(r, sigma) * jac
    end
    return val
end

@testset "Smearing kernels" begin
    sigmas = (0.5, 1.0, 2.3)

    @testset "normalization (integrates to 1 over R^3)" begin
        for sigma in sigmas
            # compact support kernels: integrate up to their cutoff
            @test radial_norm(linear, sigma; cutoff=sigma) ≈ 1.0 atol = 1e-4
            @test radial_norm(epanechnikov, sigma; cutoff=sigma) ≈ 1.0 atol = 1e-4
            @test radial_norm(cubicspline, sigma; cutoff=2 * sigma) ≈ 1.0 atol = 1e-4
            @test radial_norm(wendland, sigma; cutoff=sigma) ≈ 1.0 atol = 1e-4
            @test radial_norm(constant, sigma; cutoff=sigma) ≈ 1.0 atol = 1e-4
            @test radial_norm(cosine, sigma; cutoff=sigma) ≈ 1.0 atol = 1e-4

            # infinite-support kernels: integrate over the whole half-line
            @test radial_norm(slater, sigma) ≈ 1.0 atol = 1e-3
            @test radial_norm(gaus, sigma) ≈ 1.0 atol = 1e-3
            @test radial_norm(gengauss, sigma) ≈ 1.0 atol = 1e-3
            @test radial_norm(lorentz, sigma) ≈ 1.0 atol = 1e-3
            # cauchy has a slow 1/r^2 radial tail; looser tolerance
            @test radial_norm(cauchy, sigma; rtol=1e-7) ≈ 1.0 atol = 5e-2
        end
    end

    @testset "compact support kernels vanish beyond cutoff" begin
        sigma = 1.5
        @test linear(sigma + 1e-6, sigma) == 0
        @test epanechnikov(sigma + 1e-6, sigma) == 0.0
        @test cubicspline(2 * sigma + 1e-6, sigma) == 0.0
        @test wendland(sigma + 1e-6, sigma) == 0.0
        @test constant(sigma + 1e-6, sigma) == 0.0
        @test cosine(sigma + 1e-6, sigma) == 0.0
        # constant kernel really is constant inside its support
        @test constant(0.0, sigma) == constant(0.9 * sigma, sigma)
    end

    @testset "non-negativity and monotone decay" begin
        sigma = 1.0
        rs = range(0.0, 4.0; length=40)
        for f in (slater, gaus, cauchy, linear, epanechnikov,
                  cubicspline, wendland, gengauss, lorentz, constant, cosine)
            @test all(f(r, sigma) >= 0 for r in rs)
        end
        # smooth bell-shaped kernels should be maximal at the centre
        for f in (gaus, cauchy, epanechnikov, cubicspline, wendland, gengauss, lorentz, constant)
            @test f(0.0, sigma) >= f(0.5, sigma) >= f(1.0, sigma)
        end
        # cosine decreases monotonically to zero at its cutoff
        @test cosine(0.0, sigma) >= cosine(0.5, sigma) >= cosine(1.0, sigma)
    end

    @testset "scalar and vector forms agree" begin
        sigma = 1.3
        r0 = [1.0, -2.0, 0.5]
        for (fvec, fscal) in ((slater, slater), (gaus, gaus), (cauchy, cauchy),
                              (epanechnikov, epanechnikov), (cubicspline, cubicspline),
                              (wendland, wendland), (gengauss, gengauss),
                              (lorentz, lorentz), (constant, constant),
                              (cosine, cosine))
            for d in (0.0, 0.7, 1.5, 3.0)
                r = r0 .+ [d, 0.0, 0.0]
                @test fvec(r, r0, sigma) ≈ fscal(d, sigma)
            end
        end
        # linear takes AbstractFloat vectors
        @test linear([1.0, 1.0, 1.0], [1.0, 1.0, 1.0], 1.0) ≈ linear(0.0, 1.0)
    end

    @testset "argument validation" begin
        # 3D vector forms assert dimension 3
        @test_throws AssertionError slater([1.0, 0.0], [0.0, 0.0], 1.0)
        @test_throws AssertionError gengauss([1.0, 0.0], [0.0, 0.0], 1.0)
        @test_throws AssertionError cauchy([1.0], [0.0], 1.0)
        @test_throws AssertionError cosine([1.0, 0.0], [0.0, 0.0], 1.0)
        # linear asserts a non-negative distance
        @test_throws AssertionError linear(-1.0, 1.0)
    end

    @testset "gengauss reduces toward a Gaussian shape for p=2" begin
        sigma = 1.0
        # p=2 gives exp(-(r/sigma)^2); the ratio at r=1 and r=2 is exp(-1)/exp(-4)
        a, b = gengauss(1.0, sigma; p=2), gengauss(2.0, sigma; p=2)
        @test a / b ≈ exp(3)
        @test radial_norm((r, s) -> gengauss(r, s; p=2), sigma) ≈ 1.0 atol = 1e-3
    end
end
