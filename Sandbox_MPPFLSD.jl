using Pkg
Pkg.activate("BEP",shared=true)
# same as doing ] activate @BEP in REPL

# # Execute once then comment out
# path_dev_package_abslolute = raw"C:\Users\path_this_repo_was_coppied_to\MPPFLSD"
# Pkg.develop(path_dev_package_abslolute);
# Pkg.add(["Revise"; "Plots"]); # needs added packages

using Revise # allows for changes in packges to reflect without recompiling
using Plots
using MPPFLSD
using Trapz

params_phys = ParamsPhys(;rneedle = 10, volume0=0, deltarho=-1.1,grav=9.8,
    puffParams = Dict(:puff_rhs_method=>:puff_sigmoid, :impose_contact_angle=>true, :contact_angle=>0,
        :puff_exp_a=>-40, :puff_exp_b=>0.2,
        :puff_sig_pmin=>-0.002, :puff_sig_pmax=>10, :puff_sig_a=>-4.49, :puff_sig_b=>-10 ));
params_num = ParamsNum();

# solve for the droplet shape (Young-Laplace)

# solve for the droplet shape (Young-Laplace)
vars_sol, vars_num, _ = gen_single_drop(params_phys, params_num; verbose=true);

# post processing and plotting

volume, area = calculate_volume_area(vars_sol, vars_num; verbose=true);

kappas, kappap = find_curvature(vars_sol, vars_num);

#plot(vars_sol.r,vars_sol.z; xaxis=:log,yaxis=:log)
shape_plt = plot_shape(vars_sol.r, vars_sol.z, label="y1")
plot!([0; params_phys.rneedle], repeat([-params_phys.volume0/(pi*params_phys.rneedle)], 2), label="-H", color=:green)
plot!([params_phys.rneedle; params_phys.rneedle], [0; -params_phys.volume0/(pi*params_phys.rneedle)], label="-Rd", color=:red )
xlims!(-0.01, params_phys.rneedle*1.5)
volume_between_0_and_z_of_r = trapz(vars_sol.r.^2*pi,vars_sol.z)
volume_between_bottom_dish_and_z_of_r = trapz(vars_sol.r.^2*pi,(vars_sol.z .- -params_phys.volume0/(pi*params_phys.rneedle)))
titletext = "volume between z(r) and 0 = $(round(volume_between_0_and_z_of_r;digits=4)) \n volume below z(r) = $(round(volume_between_bottom_dish_and_z_of_r;digits=4))"
plot!(shape_plt,legend=:bottomright, plot_titlefontsize=8, plot_title=titletext)


curv_plt = plot_curvature(vars_sol.z, kappas, kappap);
combined_plt = plot(shape_plt, curv_plt, layout=(2); show=true)
