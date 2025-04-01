function update_numerical_grid(
    vars_sol::VarsSol,
    vars_num::VarsNum;
    elastic::Bool = false)::VarsNum

    # the integration and differentation matrices in the solution state
    if ~elastic # simple case
        vars_num.ws = vars_num.w/vars_sol.C; 
        vars_num.Ds = vars_sol.C*vars_num.D; 
        vars_num.s = vars_num.s0/vars_sol.C;
        #vars_num.wsmat = vars_num.wmat/vars_sol.C;
        vars_num.C = vars_sol.C;

    else
        # the integration and differentation matrices in the deformed state
        # NOTE: this construction of Ds is simlar to first applying D*f/C,
        # and then dividing the components by the components of (1/lams)
        vars_num.ws = vars_num.w.*vars_sol.lams'/vars_num.C; 
        vars_num.Ds = vars_num.C*vars_num.D.*repeat((1 ./vars_sol.lams); inner=(1, vars_num.N)); 
    
        # get the integration matrix in th reference state
        # NOTE: this is NOT the integration matrix for the deformed state!
        error("vars_num.wmat not implemented, needed for elastic case (requires function `intmat` first)")
        vars_num.wsstarmat = vars_num.wmat/vars_num.C;
    
        # compute the value of s in the deformed state
        vars_num.s = vars_num.wsstarmat*vars_sol.lams;
    end    

    return vars_num
end