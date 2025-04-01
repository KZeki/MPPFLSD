import ..MPPFLSD: ParamsPhys,
                  ParamsNum,
                  gen_single_drop,
                  calculate_volume_area,
                  find_curvature,
                  plot_shape,
                  plot_curvature
import ..Plots: plot

"""
    # Function for running complete set of simple example + plotting

    params_phys = ParamsPhys();
    params_num = ParamsNum();

    # Solve for the droplet shape (Young-Laplace)
    vars_sol, vars_num = gen_single_drop(params_phys, params_num; verbose=false)

    # Post processing and plotting
    volume, area = calculate_volume_area(vars_sol, vars_num; verbose=true);
    kappas, kappap = find_curvature(vars_sol, vars_num);

    # Plot results
    plot_shape(vars_sol.r, vars_sol.z, 1);
    plot_curvature(vars_sol.z, kappas, kappap, 2);


"""
function example_simple()

    params_phys = ParamsPhys();
    params_num = ParamsNum();

    # Solve for the droplet shape (Young-Laplace)
    vars_sol, vars_num = gen_single_drop(params_phys, params_num; verbose=false)

    # Post processing and plotting
    volume, area = calculate_volume_area(vars_sol, vars_num; verbose=true);
    kappas, kappap = find_curvature(vars_sol, vars_num);

    # Plot results
    shape_plt = plot_shape(vars_sol.r, vars_sol.z; show=false);
    curv_plt = plot_curvature(vars_sol.z, kappas, kappap; show=false);
    combined_plt = plot(shape_plt, curv_plt, layout=(1,2); show=true)

    return combined_plt
end