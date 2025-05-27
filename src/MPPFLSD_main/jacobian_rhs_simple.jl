function jacobian_rhs_simple(
    params_phys::ParamsPhys,
    vars_sol::VarsSol,
    vars_num::VarsNum;
    rhs_version::Union{Symbol,Function}=:not_def)::Tuple{Matrix,Vector}
    
    if rhs_version == :puff_sigmoid
        return _puff_sigmoid(params_phys,vars_sol,vars_num)
    elseif rhs_version == :puff_exp
        return _puff_exp(params_phys,vars_sol,vars_num)
    elseif rhs_version == :puff_test
        return _puff_test(params_phys,vars_sol,vars_num)
    elseif rhs_version == :pendant_drop
        return _pendant_drop(params_phys,vars_sol,vars_num)
    elseif typeof(rhs_version) == Function
        return rhs_version(params_phys,vars_sol,vars_num)
    else
        error("Kwarg `rhs_version` of function `jacobian_rhs_simple` with value `$rhs_version` is not defined or found.
            Try passing a Function that solves for A and b or try using a valid Symbol `:puff_sigmoid`, `:puff_exp`, `:puff_test` or `:pendant_drop`")
    end
end

function _puff_test(
    params_phys::ParamsPhys,
    vars_sol::VarsSol,
    vars_num::VarsNum)::Tuple{Matrix,Vector}

    error("no puff function is currently tested, this is just a placeholder")

    return (A, b)
end


#
# Puff with sigmoid form
#
function _puff_sigmoid(
    params_phys::ParamsPhys,
    vars_sol::VarsSol,
    vars_num::VarsNum)::Tuple{Matrix,Vector}

    D = vars_num.D0;
    w = vars_num.w0;
    N = vars_num.N;

    r = vars_sol.r;
    z = vars_sol.z;
    psi = vars_sol.psi;
    C = vars_sol.C;
    p0 = vars_sol.p0;

    g = params_phys.grav;
    sigma = params_phys.sigma;
    rho = params_phys.deltarho;
    V = params_phys.volume0;
    
    # initialize some variables
    Z = zeros(N,N);            # matrix filled with zeros
    IDL = [1 zeros(1,N-1)]; # line with single one and rest zeros
    ZL = zeros(1,N);         # line completely filled with zeros
    b = ones(3*N+2,1); # solution vector and right hand side
    
    # determine r from psi
    A11 = C*D;
    A13 = Diagonal(sin.(psi));
    A14 = D*r;
    b1 = -(C*D*r-cos.(psi));
    
    # determine z from psi 
    A22 = C*D;
    A23 = Diagonal(-cos.(psi));
    A24 = D*z;
    b2 = -(C*D*z-sin.(psi));
    
    # determine psi from Laplace law
    a = params_phys.puffParams[:puff_sigmoid_a];
    b = params_phys.puffParams[:puff_sigmoid_b];
    #a = -4.49
    #b = -10
    pmin = params_phys.puffParams[:puff_sigmoid_pmin];
    pmax = params_phys.puffParams[:puff_sigmoid_pmax];
    #pmax = 10
    #pmin = -0.002
    #A31 = -sigma*Diagonal(sin.(psi)./r.^2) .- Diagonal((a*r.*exp.(-r.^2/(2*b^2)))/b.^2);
    A31 = -Diagonal( (b*exp.(a.-b*r).*(pmax-pmin)) ./ ( exp.(a.-b*r) .+ 1 ).^2 ) - Diagonal(sigma*sin.(psi)./r.^2 )
    A32 = g*rho*Diagonal(ones(N));
    A33 = C*sigma*D + sigma*Diagonal(cos.(psi)./r);
    A34 = sigma*(D*psi);
    A35 = -ones(N, 1);

    b3 = p0 .- g*rho*z .- sigma*(C*D*psi .+ sin.(psi)./r) .+ (pmin .+ (pmax-pmin)./(exp.(a .- b*r) .+ 1))


    # impose the needle radius as a BC (imposes the domain length)
    A41 = reverse(IDL);
    b4 = (params_phys.rneedle-r[end]);
    
    # determine pressure - use volume
    A51 = 2*pi * w' .* (z .* cos.(psi))'
    A52 = 2*pi * w' .* (r .* cos.(psi))'
    A53 = -2*pi * w' .* (r .* z .* sin.(psi))'
    A54 = -V

    b5 = +C*V - 2*pi * w' * (r .* z .* cos.(psi))

    
    # boundary condition r(0) = 0
    A11[1,:] = IDL;
    A13[1,:] = ZL;
    A14[1] = 0;
    b1[1] = -r[1];

    #impose_contact_angle = params_phys.impose_contact_angle;
    impose_contact_angle = params_phys.puffParams[:impose_contact_angle];
    contact_angle = params_phys.puffParams[:contact_angle];
    if impose_contact_angle    
        # boundary condition psi(s0) = contact_angle
        A31[2,:] = ZL; 
        A32[2,:] = ZL; 
        A33[2,:] = reverse(IDL); 
        A34[2,:] .= 0; 
        A35[2,:] .= 0;
        b3[2] = contact_angle .- psi[end];
    else
        # boundary condition z(s0) = 0
        A22[1,:] = reverse(IDL); 
        A23[1,:] = ZL; 
        A24[1] .= 0;
        b2[1] = -z[end];
    end
    
    # boundary condition phi(0) = 0
    A31[1,:] = ZL;
    A32[1,:] = ZL;

    A33[1,:] = IDL; 
 
    A34[1,:] .= 0; 
    A35[1,:] .= 0;
    b3[1] = -psi[1];

    
    # assemble matrices
    Z1 = zeros(N,1);

    A =[[A11   Z A13 A14  Z1];
        [  Z A22 A23 A24  Z1];
        [A31 A32 A33 A34 A35];
        [A41  ZL  ZL   0   0];
        [A51 A52 A53 A54   0]];
    b = [b1;b2;b3;b4;b5];
    
    return (A,b)
end

#
# Puff with exponential form
#
function _puff_exp(
    params_phys::ParamsPhys,
    vars_sol::VarsSol,
    vars_num::VarsNum)::Tuple{Matrix,Vector}
    
    D = vars_num.D0;
    w = vars_num.w0;
    N = vars_num.N;

    r = vars_sol.r;
    z = vars_sol.z;
    psi = vars_sol.psi;
    C = vars_sol.C;
    p0 = vars_sol.p0;

    g = params_phys.grav;
    sigma = params_phys.sigma;
    rho = params_phys.deltarho;
    V = params_phys.volume0;
    
    # initialize some variables
    Z = zeros(N,N);            # matrix filled with zeros
    IDL = [1 zeros(1,N-1)]; # line with single one and rest zeros
    ZL = zeros(1,N);         # line completely filled with zeros
    b = ones(3*N+2,1); # solution vector and right hand side
    
    # determine r from psi
    A11 = C*D;
    A13 = Diagonal(sin.(psi));
    A14 = D*r;
    b1 = -(C*D*r-cos.(psi));
    
    # determine z from psi 
    A22 = C*D;
    A23 = Diagonal(-cos.(psi));
    A24 = D*z;
    b2 = -(C*D*z-sin.(psi));
    
    # determine psi from Laplace law
    a = params_phys.puffParams[:puff_exp_a];
    b = params_phys.puffParams[:puff_exp_b];
    #a = -40
    #b = 0.2
    A31 = -sigma*Diagonal(sin.(psi)./r.^2) .- Diagonal((a*r.*exp.(-r.^2/(2*b^2)))/b^2);
    A32 = g*rho*Diagonal(ones(N));
    A33 = C*sigma*D + sigma*Diagonal(cos.(psi)./r);
    A34 = sigma*(D*psi);
    A35 = -ones(N);
    b3 = p0 .- g*rho*z .- sigma*(C*D*psi .+sin.(psi)./r) .- a*exp.(-r.^2/(2*b^2));

    # impose the needle radius as a BC (imposes the domain length)
    A41 = reverse(IDL);
    b4 = (params_phys.rneedle-r[end]);
    
    A51 = 2*pi * w' .* (z .* cos.(psi))'
    A52 = 2*pi * w' .* (r .* cos.(psi))'
    A53 = -2*pi * w' .* (r .* z .* sin.(psi))'
    A54 = -V
    b5 = -(-C*V + 2*pi * w' * (r .* z .* cos.(psi)))
    
    # boundary condition r(0) = 0
    A11[1,:] = IDL;
    A13[1,:] = ZL;
    A14[1] = 0.0;
    b1[1] = -r[1];

    impose_contact_angle = params_phys.puffParams[:impose_contact_angle];
    contact_angle = params_phys.puffParams[:contact_angle];
    if impose_contact_angle    
        # boundary condition psi(s0) = contact_angle
        A31[2,:] = ZL; 
        A32[2,:] = ZL; 
        A33[2,:] = reverse(IDL); 
        A34[2,:] .= 0.0; 
        A35[2,:] .= 0.0;
        b3[2] = contact_angle .- psi[end];
    else
        # boundary condition z(s0) = 0
        A22[1,:] = reverse(IDL); 
        A23[1,:] = ZL; 
        A24[1] = 0.0;
        b2[1] = -z[end];
    end
    
    # boundary condition phi(0) = 0
    A31[1,:] = ZL;
    A32[1,:] = ZL;

    A33[1,:] = IDL; 
 
    A34[1,:] .= 0.0; 
    A35[1,:] .= 0.0;
    b3[1] = -psi[1];

    
    # assemble matrices
    Z1 = zeros(N,1);
    
    A =[[A11   Z A13 A14  Z1];
        [  Z A22 A23 A24  Z1];
        [A31 A32 A33 A34 A35];
        [A41  ZL  ZL   0   0];
        [A51 A52 A53 A54   0]];
    b = [b1;b2;b3;b4;b5];
 
    return (A,b)
end


#
# Pendant drop, no puff
#
function _pendant_drop(
    params_phys::ParamsPhys,
    vars_sol::VarsSol,
    vars_num::VarsNum)::Tuple{Matrix,Vector}
    
    D = vars_num.D0;
    w = vars_num.w0;
    N = vars_num.N;

    r = vars_sol.r;
    z = vars_sol.z;
    psi = vars_sol.psi;
    C = vars_sol.C;
    p0 = vars_sol.p0;

    g = params_phys.grav;
    sigma = params_phys.sigma;
    rho = params_phys.deltarho;
    V = params_phys.volume0;
    
    # initialize some variables 
    Z = zeros(N,N);            # matrix filled with zeros
    IDL = [1 zeros(1,N-1)]; # line with single one and rest zeros
    ZL = zeros(1,N);         # line completely filled with zeros 
    b = ones(3*N+2,1); # solution vector and right hand side
    A11 = C*D;
    A13 = Diagonal(sin.(psi));
    A14 = D*r;
    b1 = -(C*D*r-cos.(psi));
    # determine r from psi
    A11 = C*D;
    A13 = Diagonal(sin.(psi));
    A14 = D*r;
    b1 = -(C*D*r-cos.(psi));
    
    # determine z from psi 
    A22 = C*D;
    A23 = Diagonal(-cos.(psi));
    A24 = D*z;
    b2 = -(C*D*z-sin.(psi));
    
    # determine psi from Laplace law
    A31 = -sigma*Diagonal(sin.(psi)./r.^2);
    A32 = g*rho*Diagonal(ones(N));
    A33 = C*sigma*D + sigma*Diagonal(cos.(psi)./r);
    A34 = sigma*(D*psi);
    A35 = -ones(N,1);
    b3 = p0.-g*rho*z-sigma*(C*D*psi+sin.(psi)./r);
    
    # impose the needle radius as a BC (imposes the domain length)
    A41 = reverse(IDL);
    b4 = (params_phys.rneedle-r[end]);
    
    # determine pressure - use volume
    A51 = pi*2*(w.*r.*sin.(psi))';
    A53 = pi*(w.*r.^2 .*cos.(psi))';
    A54 = -V;
    b5 = -(pi*w'*(r.^2 .*sin.(psi)).-C*V);
    
    # boundary condition r(0) = 0
    A11[1,:] = IDL;
    A13[1,:] = ZL; 
    A14[1] = 0;
    b1[1] = -r[1];
    
    # boundary condition z(s0) = 0
    A22[1,:] = reverse(IDL); 
    A23[1,:] = ZL; 
    A24[1] = 0;
    b2[1] = -z[end];
    
    # boundary condition phi(0) = 0
    A31[1,:] = ZL;
    A32[1,:] = ZL;

    A33[1,:] = IDL; 
 
    A34[1,:] .= 0; 
    A35[1,:] .= 0;
    b3[1] = -psi[1];

    
    # assemble matrices
    Z1 = zeros(N,1);
     
    A = [[A11   Z A13 A14  Z1];
       [  Z A22 A23 A24  Z1];
       [A31 A32 A33 A34 A35];
       [A41  ZL  ZL   0   0];
       [A51 Z1' A53 A54   0]];
     
    b = [b1;b2;b3;b4;b5];

    return (A,b)
end