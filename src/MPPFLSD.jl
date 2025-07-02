module MPPFLSD
    #######################
    ### Begin - Imports ###

    using LinearAlgebra
    import Distributions
    import DataFrames
    import Interpolations
    import Optim
    import Plots
    import Statistics
    import StableRNGs
    import Trapz
    #import ProtoStructs
    #using ApproxFun: points, Chebyshev, ..
    #import Statistics: mean
    
    ###   End - Imports ###
    #######################    

    #######################
    ### Begin - Exports ###
    
    export ParamsPhys, ParamsNum
    export diffmat, chebpts, introw, clencurt
    export fit_experiment_data
    export gen_single_drop, calculate_volume_area, find_curvature
    export compute_profile_error, optim_params_to_error_calc, run_profile_optimization
    export plot_shape, plot_curvature
    export example_simple
    
    ###   End - Exports ###
    #######################

    ####################
    ### Begin - MAIN ###

    module MPPFLSD_structs
        export AbstractParams, ParamsPhys, ParamsNum
        export VarsNum, VarsShape, VarsSol
        export SimOutput
        include(raw"MPPFLSD_structs/MPPFLSD_structs.jl")
    end
    import .MPPFLSD_structs: ParamsPhys, ParamsNum

    module MPPFLSD_util
        export rms, diffmat, chebpts, introw, clencurt
        export fit_experiment_data
        include(raw"MPPFLSD_util/MPPFLSD_util.jl")
    end
    import .MPPFLSD_util: diffmat, chebpts, introw, clencurt
    import .MPPFLSD_util: fit_experiment_data

    module MPPFLSD_main
        export gen_single_drop, calculate_volume_area, find_curvature
        include(raw"MPPFLSD_main/MPPFLSD_main.jl")
    end
    import .MPPFLSD_main: gen_single_drop, calculate_volume_area, find_curvature

    module MPPFLSD_optim
        export compute_profile_error, optim_params_to_error_calc, run_profile_optimization
        include(raw"MPPFLSD_optim/MPPFLSD_optim.jl")
    end
    import .MPPFLSD_optim: compute_profile_error, optim_params_to_error_calc, run_profile_optimization

    module MPPFLSD_plot
        export plot_shape, plot_curvature
        include(raw"MPPFLSD_plot/MPPFLSD_plot.jl")
    end
    import .MPPFLSD_plot: plot_shape, plot_curvature

    module MPPFLSD_example
        export example_simple
        include(raw"MPPFLSD_example/MPPFLSD_example.jl")
    end
    import .MPPFLSD_example: example_simple

    ###  End - MAIN ###
    ###################
end