"""
# compute\\_profile\\_error

# Inputs:
`r_sim` - r coordinates of simulation data\\
`z_sim` - z coordinates of simulation data\\
`r_exp` - r coordinates of experimental data\\
`z_exp` - z coordinates of experimental data\\
`def_error_func` - Defaults `:H1`, Choose which error function to use, `:Hmean` for the mean absolute error, `:H1` for mean of squared error, `:H2` for mean of squared error divided by (1+r)\\
`scale_sim` - Scale `r_sim` and `z_sim` by a scalar amount, for example to turn m into mm, so sim and exp data use the same unit\\
`scale_exp` - Scale `r_exp` and `z_exp` by a scalar amount, for example to turn m into mm, so sim and exp data use the same unit\\
`interpolate_on` - Default `:sim`, whether to linearly interpolate on `:sim` or `:exp` data before computing the error,
    interpolating on exp may be preferred in cases with a large steady state to cavity ratio, as this can reduce the effect of overfitting since the mean is taken \\

# Outputs
`error` - computed error between sim and exp profiles, based on the error function chosen in `def_error_func`\\
"""
function compute_profile_error( r_sim::Vector{Float64},
                                z_sim::Vector{Float64},
                                r_exp::Vector{Float64},
                                z_exp::Vector{Float64};
                                def_error_func::Symbol=:H1,
                                scale_sim::Float64=1e0,
                                scale_exp::Float64=1e0, 
                                interpolate_on::Symbol=:sim)
    if interpolate_on == :sim
        r_compare = r_exp*scale_exp
        z_model_interp = interpolate_model(r_sim*scale_sim, z_sim*scale_sim, r_compare)
        z_compare = z_exp*scale_exp
    elseif interpolate_on == :exp
        r_compare = r_sim*scale_sim
        z_model_interp = interpolate_model(r_exp*scale_exp, z_exp*scale_exp, r_compare)
        z_compare = z_sim*scale_sim
    else
        error("kwarg for interpolate_on = :$interpolate_on not found, try `:sim` or `:exp`")
    end 

    if def_error_func == :Hmean
        # Mean squared error
        error = sum(abs.(z_model_interp .- z_compare))/length(z_compare)
        return error
    elseif def_error_func == :H1
        # Mean squared error
        error = sum((z_model_interp .- z_compare).^2)/length(z_compare)
        return error
    elseif def_error_func == :H2
        # Mean squared error with range scale
        error = sum((z_model_interp .- z_compare).^2 ./(1 .+ r_compare))/length(z_compare)
        return error
    elseif def_error_func == :H3
        # Mean squared error with quadratic range scale
        error = sum((z_model_interp .- z_compare).^2 ./(1 .+ r_compare.^2))/length(z_compare)

        return error
    else
        @error ("No def_error_func of `$(def_error_func)` found try `:Hmean`, `:H1`, `:H2`, `:H3`") 
    end
end
