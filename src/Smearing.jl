using LinearAlgebra

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
function slater(r:: AbstractVector, r0:: AbstractVector, sigma:: Float64; n=1)
    fnorm = 1/(factorial(n+1)*sigma^(n+2)*4*pi)
    dist = norm(r-r0)
    val = fnorm * dist^(n-1)*exp(-dist/sigma)
    return val
end

"""
Gaussian density distribution function

Args:
- r - points in space
- r0 - atom coordinates
- sigma: Float64 - atom broadening value
returns:
- g: Float64 - density at point r for atom at r0
"""
function gaus(r:: AbstractVector , r0:: AbstractVector, sigma:: Float64)
    dist = r .- r0   
    g = 1/(sigma^3*(pi*2)^(3/2))*exp(-sum(dist .* dist) / (2*sigma^2))
    return g
end
function gaus(r:: Float64, sigma:: Float64)
    dist = abs(r)    
    g = 1/(sigma^3*(pi*2)^(3/2))*exp(-(dist * dist) / (2*sigma^2))
    return g
end

"""
Cauchy distribution function

Args:
- r - points in space
- r0 - atom coordinates
- sigma: Float64 - atom broadening value
returns:
- c: Float64 - density at point r for atom at r0
"""
function cauchy(r:: AbstractVector, r0:: AbstractVector, sigma:: Float64)
    dist = r .- r0
    c = 1/(pi*sigma*(1 + sum(dist .* dist)/sigma^2))
    return c
end
function cauchy(r:: Float64, sigma:: Float64)
    dist = abs(r)
    c = 1/(pi*sigma*(1 + dist * dist/sigma^2))
    return c
end

"""
Algebraic distribution function

Args:
- r - points in space
- r0 - atom coordinates
- sigma: Float64 - atom broadening value
returns:
- c: Float64 - density at point r for atom at r0
"""
function algebraic(r:: AbstractVector, r0:: AbstractVector, sigma:: Float64; alpha=3)
    dist = r .- r0
    dist = sum(dist .* dist)^(1/2)
    fnorm = alpha*(sigma^(3-2/alpha))/(4*pi^2)
    return fnorm/(sum(dist)^(alpha) + sigma^alpha)
end
function algebraic(r:: Float64, sigma:: Float64; alpha=3)
    dist = abs(r)
    fnorm = alpha*sin(pi/alpha)*sigma^(1-1/alpha)/(2*pi)
    return fnorm/(dist^alpha + sigma^alpha)
end

