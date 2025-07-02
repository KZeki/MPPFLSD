
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
