
function interpolate_model(r_sim::Vector, z_sim::Vector, r_exp::Vector)
    # Linear interpolation of the model's r-z surface
    itp = LinearInterpolation(r_sim, z_sim, extrapolation_bc=Line())
    z_model_interp = itp(r_exp)
    return z_model_interp
end
