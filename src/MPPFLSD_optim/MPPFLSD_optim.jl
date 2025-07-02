using ..MPPFLSD_structs: ParamsPhys, ParamsNum, SimOutput
using ..MPPFLSD_util: interpolate_model, optim_params_to_ParamsPhys
using ..MPPFLSD_main: gen_single_drop
using ..Distributions: Uniform
using ..Optim: Adam, optimize, Options, minimizer
using ..StableRNGs: AbstractRNG, StableRNG



include("compute_profile_error.jl")

include("optim_params_to_error_calc.jl")

include("run_profile_optimization.jl")