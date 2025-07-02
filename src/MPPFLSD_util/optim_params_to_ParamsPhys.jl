function optim_params_to_ParamsPhys(
    vec_params::Vector{Float64};
    default_ParamsPhys::ParamsPhys=ParamsPhys(;
        rneedle = 52e-3,
        volume0=0,
        deltarho=-(1050.0 - 1.225),
        grav=9.8,
        sigma = 0.02212,
        puffParams = Dict(
            :puff_rhs_method=>:puff_sigmoid,
            :impose_contact_angle=>true,
            :contact_angle=>0) ) )::ParamsPhys
    
    puff_sigmoid_pmax = vec_params[1]#*0 + 43.47
    puff_sigmoid_a = vec_params[2]#*0 + -4.49
    puff_sigmoid_b = vec_params[3]#*0 -729.56

    default_ParamsPhys.puffParams[:puff_sigmoid_pmax] = puff_sigmoid_pmax
    default_ParamsPhys.puffParams[:puff_sigmoid_a] = puff_sigmoid_a
    default_ParamsPhys.puffParams[:puff_sigmoid_b] = puff_sigmoid_b

    return default_ParamsPhys
end