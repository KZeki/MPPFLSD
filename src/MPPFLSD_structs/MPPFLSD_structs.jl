abstract type AbstractParams end

"""
# FieldNames
`sigma` # Surface tension (physical parameter)\\
`grav` # gravity (physical parameter)\\
`rneedle` # Needle radius (physical parameter)\\
`volume0` # Prescribed volume (physical parameter)\\
`deltaRho` # Density difference (physical parameter)\\
`worthingtonNumber` # Worthington number (needed for initial shape guess, physical parameter)\\
`area0` # Area calculation\\
"""
@kwdef mutable struct ParamsPhys <: AbstractParams
    sigma::Float64 = 4
    grav::Float64 = 1.2
    rneedle::Float64 = 10
    volume0::Float64 = 0
    deltarho::Float64 = -1.1
    Wo::Float64 = deltarho*grav*volume0/(2*pi*sigma*rneedle)
    area0::Float64 = 1.0
    puffParams::Dict = Dict(:impose_contact_angle=>true, :contact_angle=>0,
        :puff_exp_a=>-40, :puff_exp_b=>-0.2,
        :puff_sig_pmin=>-0.002, :puff_sig_pmax=>10, :puff_sig_a=>-4.49, :puff_sig_b=>-10)
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
    r::Vector{Float64} = Float64[]
    z::Vector{Float64} = Float64[]
    s::Vector{Float64} = Float64[]
end

"""
# FieldNames
`N` # The number of points in the grid (copied from params_num.N) \\
`C` # Linear scaling factor between numerical (s0) and actual (s) grid given by C = L0/L, where L0 is guess of numerical grid length and L is actual length \\
`s` # The Chebyshev points on the actual domain linearly mapped by s=s0/C \\
`s0` # The Chebyshev points on numerical grid \\
`w0` # The integration weights optained from function `introw` \\
`ws` # (ws = w/C, calcualted on last numeric grid update) \\
`D0` # The first-order differentiation matrix. \\
`Ds` # (Ds = D*s, calcualted on last numeric grid update) \\
"""
@kwdef mutable struct VarsNum <: AbstractParams
    N::Integer = 40;
    C::Float64 = 0.75
    s::Vector{Float64} = Float64[]
    s0::Vector{Float64} = Float64[]
    w0::Vector{Float64} = Float64[]
    ws::Vector{Float64} = Float64[]
    #wmat = undef
    #wmat0 = undef
    #wsmat = undef
    D0::Matrix{Float64} = Float64[;;]
    Ds::Matrix{Float64} = Float64[;;]
    #DD = undef
end

"""
# FieldNames
`r` # Droplet radius points \\
`z` # Droplet height point. \\
`psi` # Droplet angle points \\
`C` # Linear scaling factor between numerical (s0) and actual (s) grid given by C = L0/L, where L0 is guess of numerical grid length and L is actual length \\
`p0` # Pressure \\
"""
@kwdef mutable struct VarsSol <: AbstractParams
    C::Float64 = 0.75
    p0::Float64 = 1.0
    r::Vector{Float64} = Float64[]
    z::Vector{Float64} = Float64[]
    psi::Vector{Float64} = Float64[]
    #sigmas::Vector{Float64} = Float64[]
    #sigmap::Vector{Float64} = Float64[]
end