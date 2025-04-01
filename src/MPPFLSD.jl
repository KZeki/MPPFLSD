module MPPFLSD
    #######################
    ### Begin - Imports ###

    using LinearAlgebra
    import Interpolations
    import Plots
    #import ProtoStructs
    #using ApproxFun: points, Chebyshev, ..
    #import Statistics: mean
    
    ###   End - Imports ###
    #######################    

    #######################
    ### Begin - Exports ###
    
    export ParamsPhys, ParamsNum
    export diffmat, chebpts, introw, clencurt
    export gen_single_drop, calculate_volume_area, find_curvature
    export plot_shape, plot_curvature
    export example_simple
    
    ###   End - Exports ###
    #######################

    ####################
    ### Begin - MAIN ###

    module MPPFLSD_structs
        export AbstractParams, ParamsPhys, ParamsNum
        export VarsNum, VarsShape, VarsSol
        include(raw"MPPFLSD_structs/MPPFLSD_structs.jl")
    end
    import .MPPFLSD_structs: ParamsPhys, ParamsNum

    module MPPFLSD_util
        export rms, diffmat, chebpts, introw, clencurt
        include(raw"MPPFLSD_util/MPPFLSD_util.jl")
    end
    import .MPPFLSD_util: diffmat, chebpts, introw, clencurt

    module MPPFLSD_main
        export gen_single_drop, calculate_volume_area, find_curvature
        include(raw"MPPFLSD_main/MPPFLSD_main.jl")
    end
    import .MPPFLSD_main: gen_single_drop, calculate_volume_area, find_curvature

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