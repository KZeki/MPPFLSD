"""
# optim\\_params\\_to\\_error\\_calc
# Inputs:
`optim_params_vec::Vector{Float64}` - Simple parameter vector passed to optim, example [43.47;-4.49;-729.56]\\
`r_exp::Vector{Float64}` - Vector with r coordinates of experimental data\\
`z_exp::Vector{Float64};` - Vector with rz coordinates of experimental data\\
# Kwargs:
`default_params_phys::ParamsPhys` - Pass the default ParamsPhys that should be used with the parameters\\
`params_num::ParamsNum=ParamsNum(;N=40)` - Pass the ParamsNum, i.e., numeric parameters that should be used\\
`def_error_func::Symbol` - Check `compute_profile_error`, should either be `:Hmean`, `:H1` ,`:H2`\\
`scale_sim::Float64=1e0` - Scales the simulation data by this amount, usefull if sim data is in m, 
    exp data is in mm and you want to calculate error in mm,
    has effect on how the `1/(1+r)` of `:H2` scales\\
`scale_exp::Float64=1e0` - Scales the experiment data by this amount, usefull if exp data is in mm and should be in m\\
`interpolate_on::Symbol` - Either `:exp` or `:sim` in order to indicate which data should be linearly interpolated when calculating error\\
`return_profiles::Bool=false` - Whether to return simulated profiles\\

# Outputs:
`error` - computed error between experimetnal data and simulated deformed surface\\
or \\
`(error,vars_sol.r,vars_sol.z)` - if `return_profiles` is set to true also returns the r and z profiles of the simulation\\

"""
function optim_params_to_error_calc(
        optim_params_vec::Vector{Float64},
        r_exp::Vector{Float64},
        z_exp::Vector{Float64};
        default_params_phys::ParamsPhys,
        params_num::ParamsNum=ParamsNum(;N=40),
        def_error_func::Symbol=:H1,
        scale_sim::Float64=1e0,
        scale_exp::Float64=1e0,
        interpolate_on::Symbol=:exp,
        return_profiles::Bool=false)::Union{Float64,Tuple{Float64,Vector{Float64},Vector{Float64}}}
    puff_params = optim_params_to_ParamsPhys(optim_params_vec; default_ParamsPhys=default_params_phys)
    vars_sol,_,_ = gen_single_drop(puff_params,params_num;verbose=false)
    error = compute_profile_error(vars_sol.r, vars_sol.z, r_exp, z_exp; def_error_func=def_error_func, scale_sim=scale_sim, scale_exp=scale_exp,interpolate_on=interpolate_on)

    if return_profiles == false
        return error
    else
        return (error,vars_sol.r,vars_sol.z)
    end
end
