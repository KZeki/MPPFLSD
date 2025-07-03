

function run_profile_optimization(
    experiment_r::Vector{Float64},
    experiment_z::Vector{Float64};
    use_error_func::Symbol = :H1,
    scale_sim::Float64 = 1e0,
    scale_exp::Float64 = 1e0,
    interpolate_on::Symbol = :exp,
    initial_vec::Vector{Float64} = [43.47;-4.49;-729.56],
    lower_bounds::Vector{Float64} = initial_vec .- abs.(initial_vec*0.10),
    upper_bounds::Vector{Float64} = initial_vec .+ abs.(initial_vec*0.10),
    stable_rng::AbstractRNG = StableRNG(1),
    N_random_initial::Int64 = 3,
    N_samples_per_initial::Int64 = 1,
    N_max_iter::Int64 = 1_000,
    default_ParamsPhys::ParamsPhys = ParamsPhys(;
            rneedle = 52e-3,
            volume0=0,
            deltarho=-(1050.0 - 1.225),
            grav=9.8,
            sigma = 0.02212,
            puffParams = Dict(
                :puff_rhs_method=>:puff_sigmoid,
                :puff_sigmoid_pmin=>-0.002,
                :impose_contact_angle=>true,
                :contact_angle=>0) ),
    params_num::ParamsNum = ParamsNum(;N=40),
    tune_Adam_alpha::Float64 = 0.001,
    tune_Adam_beta_mean::Float64 = 0.1,
    tune_Adam_beta_var::Float64 = 0.9,
    tune_Adam_epsilon::Float64 =  1e-8)
    #
    # Begin with defining storage size for inputs
    #

    # Matrix containing the intial parameters 
    # There is one row for every intial vec that should be tested, columns contain its parameters size of "initial_vec"
    initial_params_matrix = fill(NaN,N_random_initial,length(initial_vec)); 
    # Storing the initial puffParams for performing "gen_single_drop"
    initial_puffParams_vec = fill(deepcopy(default_ParamsPhys),N_random_initial);
    
    #
    # Begin with defining storage size for outputs
    # Note that there are going to be N_random_initial*N_samples_per_initial simulations,
    # Each returning an optimized profile with beloning parameters
    #

    # Store the optimal parameters from every simulation
    output_best_fit_params_simulation_matrix = fill([NaN],N_random_initial,N_samples_per_initial);
    # Store the "optimized" error value of every simulation
    output_best_fit_error_simulation_matrix = fill(NaN,N_random_initial,N_samples_per_initial);
    # Store the optimal profile from every simulation
    output_best_fit_r_profiles_simulation_matrix = fill([NaN],N_random_initial,N_samples_per_initial);
    output_best_fit_z_profiles_simulation_matrix = fill([NaN],N_random_initial,N_samples_per_initial);
    

    # Fill in the the columns of the intial parameters matrix with random starting parameters

    random_vec = rand(stable_rng,Uniform(0,1),N_random_initial,length(initial_vec))
    for i in 1:length(initial_vec)
        initial_params_matrix[:,i] = abs(upper_bounds[i]-lower_bounds[i])*random_vec[:,i] .+ lower_bounds[i]
    end
    
    # Run simulation for every initial parameter set
    for ith_random_initial in 1:N_random_initial
        initial_puffParams_vec[ith_random_initial] = optim_params_to_ParamsPhys(
            initial_params_matrix[ith_random_initial,:];
            default_ParamsPhys = default_ParamsPhys)

        for jth_sample in 1:N_samples_per_initial

            result_j = optimize(
                x-> optim_params_to_error_calc(
                        x,
                        experiment_r,
                        experiment_z;
                        default_params_phys = deepcopy(default_ParamsPhys),
                        def_error_func = use_error_func,
                        scale_sim = scale_sim,
                        scale_exp = scale_exp,
                        interpolate_on = interpolate_on),
                initial_params_matrix[ith_random_initial,:],
                Adam(; 
                        alpha = tune_Adam_alpha, 
                        beta_mean = tune_Adam_beta_mean,
                        beta_var = tune_Adam_beta_var,
                        epsilon = tune_Adam_epsilon),
                Options(iterations=N_max_iter) )

            # Save jth_sample optimal parameters
            output_best_fit_params_simulation_matrix[ith_random_initial,jth_sample] = minimizer(result_j)

            # Save jth_sample error and, profiles r and z
            output_best_fit_error_simulation_matrix[ith_random_initial,jth_sample], 
            output_best_fit_r_profiles_simulation_matrix[ith_random_initial,jth_sample],
            output_best_fit_z_profiles_simulation_matrix[ith_random_initial,jth_sample] = optim_params_to_error_calc(
                output_best_fit_params_simulation_matrix[ith_random_initial,jth_sample],
                experiment_r,
                experiment_z;
                default_params_phys = deepcopy(default_ParamsPhys),
                params_num = params_num,
                def_error_func = use_error_func,
                scale_sim = scale_sim,
                scale_exp = scale_exp,
                interpolate_on = interpolate_on,
                return_profiles = true)
        end
    end


    return SimOutput(
        initial_params = initial_params_matrix,
        output_best_fit_params = output_best_fit_params_simulation_matrix,
        output_best_fit_error = output_best_fit_error_simulation_matrix,
        output_best_fit_r_profiles = output_best_fit_r_profiles_simulation_matrix,
        output_best_fit_z_profiles = output_best_fit_z_profiles_simulation_matrix)
end