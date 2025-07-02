"""
NUMERICAL_GRID Generates numerical grids and differentiation matrices 
using Chebyshev points.

This function creates differentiation matrices and integration weights 
for a given domain using Chebyshev points. It is designed to work with 
the Chebfun package.

# Usage:
  `vars_num` = `numerical_grid`(params_num, domain)\\

# Inputs:
  `params_num` - A structure containing the number of points in the grid.
               It should have a field 'N' representing this number.\\
  `domain` - The left and right boundaries of the 1D domain, specified as
           a vector [left\\_boundary, right\\_boundary].

# Output:
  `vars_num` - A structure containing the following fields:\\
             D: The first-order differentiation matrix.\\
             DD: The second-order differentiation matrix.\\
             w: The integration weights.\\
             s: The Chebyshev points within the specified domain.\\
             N: The number of points in the grid (copied from params_num.N).\\
"""
function numerical_grid(params_num::ParamsNum, domain::AbstractArray)::VarsNum
    
    # Initilize var_num
    vars_num = VarsNum();

    # diffmat, introw and chebpts are defined in the util module
    vars_num.D0 = diffmat(params_num.N,1,domain);
    #vars_num.DD = diffmat(params_num.N,2,domain);
    #vars_num.wmat = intmat(params_num.N,1,domain);
    vars_num.w0 = introw(params_num.N,domain);
    vars_num.s0 = chebpts(params_num.N,domain);

    vars_num.N = params_num.N; # copy for convenvience

    return vars_num
end