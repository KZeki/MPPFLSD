using LinearAlgebra: norm, I, Diagonal
using DataFrames
using Interpolations: LinearInterpolation, Line 
using Trapz: trapz
using ..MPPFLSD_structs: ParamsPhys

include("MPPFLSD_util_shape_gen.jl")

include("compute_vertical_shift_for_zero_on_volume.jl")
include("fit_experiment_data.jl")

include("interpolate_model.jl")
include("optim_params_to_ParamsPhys.jl")