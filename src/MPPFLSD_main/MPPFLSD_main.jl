using ..MPPFLSD_structs
using ..MPPFLSD_util
using ..LinearAlgebra
import ..Interpolations
import ..Plots

include("guess_shape.jl")
include("numerical_grid.jl")
include("jacobian_rhs_simple.jl")
include("solve_forward_young_laplace.jl")
include("update_numerical_grid.jl")


"""
(vars_sol,vars_num,params_phys) = gen_single_drop(params_phys::ParamsPhys,params_num::ParamsNum; verbose=true)::Tuple{VarsSol,VarsNum,ParamsPhys}
"""
function gen_single_drop(
    params_phys::ParamsPhys,
    params_num::ParamsNum;
    verbose=true)::Tuple{VarsSol,VarsNum,ParamsPhys}

    shape_guess = guess_shape(params_phys, 1000);
    shape_guess.z = zeros(1000);
    shape_guess.r = collect(LinRange(0,params_phys.rneedle,1000));
    shape_guess.s = collect(LinRange(0,params_phys.rneedle,1000));
    

    vars_num = numerical_grid(params_num, [0, shape_guess.s[end]]);

    vars_sol = solve_forward_young_laplace(params_phys, params_num, shape_guess, vars_num; verbose=verbose);
    
    vars_num = update_numerical_grid(vars_sol, vars_num; elastic = false)

    # store the converged value of area0 in the parameters
    _, params_phys.area0 = calculate_volume_area(vars_sol, vars_num; verbose=false);

    return (vars_sol,vars_num,params_phys)
end

"""
(volume, area) = calculate_volume_area(vars_sol::VarsSol, vars_num::VarsNum; verbose::Bool=true)::Tuple{Float64,Float64}
"""
function calculate_volume_area(
    vars_sol::VarsSol,
    vars_num::VarsNum;
    verbose::Bool=true)::Tuple{Float64,Float64}

    # calculate the volume and the area
    volume = pi*dot(vars_num.ws, (vars_sol.r.^2 .* sin.(vars_sol.psi) ) );
    area = pi*2*dot(vars_num.ws,vars_sol.r);

    if verbose
        println("volume = ", volume);
        println("area = ", area);
        println("pressure = ", vars_sol.p0);
    end
    return (volume, area)
end

"""
(kappas,kappap) = find_curvature(vars_sol::VarsSol, vars_num::VarsNum)::Tuple{Any,Any}
"""
function find_curvature(
    vars_sol::VarsSol,
    vars_num::VarsNum)::Tuple{Any,Any}
    # determine the curvatures
    # NOTE: kappap = sin(psi)/r, which is problematic for r=0. This is
    # solved here by taking kappap(0) = kappas(0)
    kappas = vars_num.Ds*vars_sol.psi;
    kappap = deepcopy(kappas);
    kappap[2:end] = sin.(vars_sol.psi[2:end])./vars_sol.r[2:end];

    return (kappas,kappap)
end