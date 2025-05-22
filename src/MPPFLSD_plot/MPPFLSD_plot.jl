import ..Plots: plot, plot!, xlabel!, ylabel!, xlims!, ylims!

function plot_shape(r, z; show::Bool=true, label::String="")
        
    # PLOT_SHAPE Plots the shape profile in the (r, z) plane.
    #
    # INPUTS:
    #   r       - Radial coordinates.
    #   z       - Axial coordinates.
   
    #plt = plot(r, z; markershape=:circle, aspect_ratio=1, label=label, grid=false, xlim=[-0.01, r[end]], ylim=[minimum(z),maximum(z)], show=show)
    plt = plot(r, z; markershape=:circle, label=label, grid=false, xlim=[-0.01, r[end]], ylim=[minimum(z),maximum(z)], show=show)
    
    xlabel!("r")
    ylabel!("z")
    
    return plt
end

function plot_curvature(z, kappas, kappap; show::Bool=true)
    # PLOT_CURVATURE Plots curvatures versus the z-coordinate.
    #
    # INPUTS:
    #   z       - z-coordinates.
    #   kappas  - Meridional curvature.
    #   kappap  - Azimuthal curvature.

    ## plot the curvatures versus the z-coordinate

    plt = plot(z, kappas; label="κₛ", linewidth=2, aspect_ratio=1, grid=false, show=show)
    xlabel!("z")
    ylabel!("κ")
    plot!(z, kappap; label="κᵩ",linewidth=2)
    plot!(z, kappas+kappap; label="κₛ + κᵩ", linewidth=2)
    xlims!(z[1],0)

    return plt
end