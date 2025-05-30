"""
SOLVE_FORWARD_YOUNG_LAPLACE Solves the Young-Laplace equation for the
shape of a drop.

This function solves the Young-Laplace equation for a constant surface 
tension, using a Newton-Raphson iteration.

# Usage:
  [vars_sol, vars_num] = `solve_forward_young_laplace`(params_phys, params_num)

# Inputs:
  `params_phys` - A structure containing physical parameters:
                deltarho (density difference), grav (gravity), volume0 (initial volume),
                sigma (surface tension), rneedle (radius of the needle), 
                Wo (Worhtington number).
  `params_num` - A structure containing numerical parameters:
               N (number of points), eps_fw (tolerance for convergence),
               maxiter (maximum iterations), alpha (relaxation parameter).

# Outputs:
  `vars_sol` - A structure containing solution variables:
             r, z (coordinates of the interface), psi (angle of the tangent to the interface),
             C (stretch parameter), p0 (pressure).
  `vars_num` - A structure containing numerical variables used in the solution:
             Differentiation and integration matrices, Chebyshev points, and modified
             matrices post-solution (ws, Ds).
"""
function solve_forward_young_laplace(
    params_phys::ParamsPhys,
    params_num::ParamsNum,
    shape_guess::VarsShape,
    vars_num::VarsNum;
    verbose=true)::VarsSol
      
    # interpolate the shape in the Chebyshev points
    #r = interp1(shape_guess.s,shape_guess.r,vars_num.s0);
    #z = interp1(shape_guess.s,shape_guess.z,vars_num.s0);
    r = Interpolations.linear_interpolation(shape_guess.s, shape_guess.r, extrapolation_bc=Interpolations.Line())(vars_num.s0)
    z = Interpolations.linear_interpolation(shape_guess.s, shape_guess.z, extrapolation_bc=Interpolations.Line())(vars_num.s0)

    psi = atan.(vars_num.D0*z,vars_num.D0*r);   # intial psi value 
    C = 1;                                # initial stretch parameter
    p0 = 2*params_phys.sigma/params_phys.rneedle;   # initial pressure
    u = ones(3*params_num.N+2,1);             # initial solution vector
    
    # store some variables for the iteration
    iter = 0;
    vars_sol = VarsSol();
    vars_sol.r = r;
    vars_sol.z = z;
    vars_sol.psi = psi;
    vars_sol.C = C;
    vars_sol.p0 = p0; 

    rhs_version = :puff_exp
    try
        rhs_version = params_phys.puffParams[:puff_rhs_method]
    catch
        @warn("Dictionary `puffParams` (a field of the struct params_phys::ParamsPhys (argument of function `jacobian_rhs_simple`)) does not have a key `:puff_rhs_method` defined, default value `$rhs_version` is used`")
    end

    # start the Newton-Raphson iteration
    while rms(u) > params_num.epsilon
    
        iter = iter + 1;
        
        if iter > params_num.nMaxIter
            #error("Newton-Raphson iteration of forward young-laplace did not converge in nMaxIter=$(params_num.nMaxIter) steps!")
            @warn("Newton-Raphson iteration of forward young-laplace did not converge in nMaxIter=$(params_num.nMaxIter) steps!")
            return vars_sol
        end    
    
        # build the Jacobian and RHS

        A, b = jacobian_rhs_simple(params_phys,vars_sol,vars_num;rhs_version=rhs_version);

        # solve the system of equations
        u = A\b;
        
        # update variables
        vars_sol.r   = vars_sol.r   + u[1:params_num.N];
        vars_sol.z   = vars_sol.z   + u[params_num.N+1:2*params_num.N];
        vars_sol.psi = vars_sol.psi + u[2*params_num.N+1:3*params_num.N];
        vars_sol.C   = vars_sol.C   + u[3*params_num.N+1];
        vars_sol.p0  = vars_sol.p0  + u[3*params_num.N+2];
        
        if verbose
            println("iter $iter: rms(u) = $(rms(u)) \n");
        end

    end

    if any( vars_sol.r .< -1e-6 )
        error("Negative r-coordinates encountered")
    end

    # Following not used in simple case
    # vars_sol.sigmas = params_phys.sigma*ones(vars_num.N,1);
    # vars_sol.sigmap = params_phys.sigma*ones(vars_num.N,1);

    return vars_sol
end