abstract type AbstractParams end

"""
# FieldNames
`sigma` # Surface tension (physical parameter)\\
`grav` # gravity (physical parameter)\\
`rneedle` # Needle radius (physical parameter)\\
`volume0` # Prescribed volume (physical parameter)\\
`deltaRho` # Density difference (physical parameter)\\
`worthingtonNumber` # Worthington number (needed for initial shape guess, physical parameter)\\
`area0` # Area calcu\\
"""
@kwdef mutable struct ParamsPhys <: AbstractParams
    sigma::Float64 = 4
    grav::Float64 = 1.2
    rneedle::Float64 = 1.4
    volume0::Float64 = 16
    deltarho::Float64 = 1.1
    Wo::Float64 = deltarho*grav*volume0/(2*pi*sigma*rneedle)
    area0::Float64 = 1.0
end

"""
# FieldNames
`N` # Number of gridpoints (numeric parameter)\\
`nMaxIter` # Maximum amount of iterations (numeric parameter)\\
`epsilon` # Forward Convergance Criterion (numeric parameter)\\
"""
@kwdef mutable struct ParamsNum <: AbstractParams
    N::Int64 = 40
    nMaxIter::Int64 = 100
    epsilon::Float64 = 1e-12
end

@kwdef mutable struct VarsShape <: AbstractParams
    r = undef
    z = undef
    s = undef
end

"""
# FieldNames
`D` # The first-order differentiation matrix. \\
`w` # The integration weights. \\
`s` # The Chebyshev points within the specified domain. \\
`N` # The number of points in the grid (copied from params_num.N). \\
`C` # Scaling factor for s = s0/C
"""
@kwdef mutable struct VarsNum <: AbstractParams
    D = undef
    #DD = undef
    #wmat = undef
    w = undef
    s0 = undef
    D0 = undef
    w0 = undef
    #wmat0 = undef
    N = undef
    ws = undef
    Ds = undef
    s = undef
    #wsmat = undef
    C = undef
end

"""
# FieldNames
`r` # Droplet radius points \\
`z` # Droplet height point. \\
`psi` # Droplet angle points \\
`C` # Sclaing factor for s = s0/C\\
`p0` # Pressure
"""
@kwdef mutable struct VarsSol <: AbstractParams
    r = undef
    z = undef
    psi = undef
    C = undef
    p0 = undef
    #sigmas = undef
    #sigmap = undef
end