using LinearAlgebra
using SpecialFunctions: gamma

"""
Slater type density distribution function

Args:
- r - points in space
- r0 - atom coordinates
- sigma: Float64 - atom broadening value
- n: Int - power of the function
returns:
- d: Float64 - density at point r for atom at r0
"""
function slater(r:: Float64, sigma:: Float64; n=1)
    fnorm = 1/(factorial(n+1)*sigma^(n+2)*4*pi)
    dist = abs(r)
    return fnorm * dist^(n-1)*exp(-dist/sigma)
end
function slater(r::AbstractVector{<:Real}, r0::AbstractVector{<:Real}, sigma:: Float64; n=1)
    @assert length(r) == length(r0) == 3 "r and r0 must be dimension 3"
    return slater(norm(r .- r0), sigma; n=n)
end

"""
Gaussian density distribution function
checked and found to be normalized

Args:
- r - points in space
- r0 - atom coordinates
- sigma: Float64 - atom broadening value
returns:
- g: Float64 - density at point r for atom at r0
"""
function gaus(r::AbstractVector{<:Real}, r0::AbstractVector{<:Real}, sigma:: Float64)
    dist = r .- r0   
    g = 1/(sigma^3*(pi*2)^(3/2))*exp(-sum(dist .* dist) / (2*sigma^2))
    return g
end
function gaus(r:: Float64, sigma:: Float64)
    dist = abs(r)    
    g = 1/(sigma^3*(pi*2)^(3/2))*exp(-dist^2 / (2*sigma^2))
    return g
end

"""
Cauchy distribution function
Attention: this is a very broad function, that converges very slow to zero
Args:
- r - points in space
- r0 - atom coordinates
- sigma: Float64 - atom broadening value
returns:
- c: Float64 - density at point r for atom at r0
"""
function cauchy(r::AbstractVector{<:Real}, r0::AbstractVector{<:Real}, sigma:: Float64)
    @assert length(r) == length(r0) == 3 "dimension of vectors must be 3"
    # @error "This function does not converge"
    dr = 0.

    for i in eachindex(r)
        dr += (r[i] - r0[i])^2
    end
    return cauchy(sqrt(dr), sigma) 
end
function cauchy(r:: Float64, sigma:: Float64)
    # @error "This function does not converge"
    dist = abs(r)
    n = 4*pi*pi*sigma^3/2^(3/2)
    c = 1/(1 + (dist/sigma)^4)
    return c/n
end

"""
Linear triangle function in 1 and 3d
checked and found to be normalized

Args:
- r  position in space
- r0 center of smearing
- sigma  broadening value
"""
function linear(r::AbstractFloat, sigma::AbstractFloat)
    @assert r >= 0 "Distance r=$r should be greater or equal than 0";
    if r > sigma 
        return 0
    else 
        n = pi*sigma^2/3
        return (-1/sigma^2 * r + 1/sigma)/n
    end
end
function linear(r::AbstractVector{<:AbstractFloat}, r0::AbstractVector{<:AbstractFloat}, sigma::AbstractFloat)
    dr = 0.
    @assert length(r) == length(r0) == 3 "r and r0 must be dimension 3"
    for i in eachindex(r)
        dr += (r[i] - r0[i])^2
    end
    dr = sqrt(dr)
    return linear(dr, sigma)
end

"""
Epanechnikov density distribution function
compact support: zero for distances larger than sigma
normalized in 3D (checked analytically)

Args:
- r - points in space
- r0 - atom coordinates
- sigma: Float64 - atom broadening value (cutoff radius)
returns:
- e: Float64 - density at point r for atom at r0
"""
function epanechnikov(r:: Float64, sigma:: Float64)
    dist = abs(r)
    if dist > sigma
        return 0.0
    end
    norm_const = 15/(8*pi*sigma^3)
    return norm_const * (1 - (dist/sigma)^2)
end
function epanechnikov(r::AbstractVector{<:Real}, r0::AbstractVector{<:Real}, sigma:: Float64)
    @assert length(r) == length(r0) == 3 "r and r0 must be dimension 3"
    return epanechnikov(norm(r .- r0), sigma)
end

"""
Cubic spline (Monaghan SPH) density distribution function
compact support: zero for distances larger than 2*sigma
C2 smooth, normalized in 3D (checked analytically)

Args:
- r - points in space
- r0 - atom coordinates
- sigma: Float64 - atom broadening value (smoothing length)
returns:
- c: Float64 - density at point r for atom at r0
"""
function cubicspline(r:: Float64, sigma:: Float64)
    q = abs(r)/sigma
    norm_const = 1/(pi*sigma^3)
    if q <= 1
        return norm_const * (1 - 1.5*q^2 + 0.75*q^3)
    elseif q <= 2
        return norm_const * 0.25*(2 - q)^3
    else
        return 0.0
    end
end
function cubicspline(r::AbstractVector{<:Real}, r0::AbstractVector{<:Real}, sigma:: Float64)
    @assert length(r) == length(r0) == 3 "r and r0 must be dimension 3"
    return cubicspline(norm(r .- r0), sigma)
end

"""
Wendland C2 density distribution function
compact support: zero for distances larger than sigma
smooth and positive definite, normalized in 3D (checked analytically)

Args:
- r - points in space
- r0 - atom coordinates
- sigma: Float64 - atom broadening value (cutoff radius)
returns:
- w: Float64 - density at point r for atom at r0
"""
function wendland(r:: Float64, sigma:: Float64)
    q = abs(r)/sigma
    if q > 1
        return 0.0
    end
    norm_const = 21/(2*pi*sigma^3)
    return norm_const * (1 - q)^4 * (1 + 4*q)
end
function wendland(r::AbstractVector{<:Real}, r0::AbstractVector{<:Real}, sigma:: Float64)
    @assert length(r) == length(r0) == 3 "r and r0 must be dimension 3"
    return wendland(norm(r .- r0), sigma)
end

"""
Generalized (super-)Gaussian density distribution function exp(-(r/sigma)^p)
p=2 recovers a Gaussian shape, larger p gives a flatter top and sharper edge
normalized in 3D (checked analytically)

Args:
- r - points in space
- r0 - atom coordinates
- sigma: Float64 - atom broadening value
- p: power of the exponent (default 4)
returns:
- g: Float64 - density at point r for atom at r0
"""
function gengauss(r:: Float64, sigma:: Float64; p=4)
    dist = abs(r)
    norm_const = p/(4*pi*sigma^3*gamma(3/p))
    return norm_const * exp(-(dist/sigma)^p)
end
function gengauss(r::AbstractVector{<:Real}, r0::AbstractVector{<:Real}, sigma:: Float64; p=4)
    @assert length(r) == length(r0) == 3 "r and r0 must be dimension 3"
    return gengauss(norm(r .- r0), sigma; p=p)
end

"""
Lorentzian density distribution function 1/(1 + (r/sigma)^2)^3
like cauchy but with a faster 1/r^6 tail so it converges quickly
normalized in 3D (checked analytically)

Args:
- r - points in space
- r0 - atom coordinates
- sigma: Float64 - atom broadening value
returns:
- l: Float64 - density at point r for atom at r0
"""
function lorentz(r:: Float64, sigma:: Float64)
    dist = abs(r)
    norm_const = 4/(pi^2*sigma^3)
    return norm_const / (1 + (dist/sigma)^2)^3
end
function lorentz(r::AbstractVector{<:Real}, r0::AbstractVector{<:Real}, sigma:: Float64)
    @assert length(r) == length(r0) == 3 "r and r0 must be dimension 3"
    return lorentz(norm(r .- r0), sigma)
end

"""
Constant (uniform sphere / top-hat) density distribution function
compact support: constant value for distances up to sigma, zero beyond
normalized in 3D (checked analytically)

Args:
- r - points in space
- r0 - atom coordinates
- sigma: Float64 - atom broadening value (cutoff radius)
returns:
- c: Float64 - density at point r for atom at r0
"""
function constant(r:: Float64, sigma:: Float64)
    if abs(r) > sigma
        return 0.0
    end
    return 3/(4*pi*sigma^3)
end
function constant(r::AbstractVector{<:Real}, r0::AbstractVector{<:Real}, sigma:: Float64)
    @assert length(r) == length(r0) == 3 "r and r0 must be dimension 3"
    return constant(norm(r .- r0), sigma)
end



