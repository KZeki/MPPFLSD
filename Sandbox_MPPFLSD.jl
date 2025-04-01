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

params_phys = ParamsPhys();
params_num = ParamsNum();

# solve for the droplet shape (Young-Laplace)

vars_sol, vars_num, _ = gen_single_drop(params_phys, params_num; verbose=false);

# post processing and plotting

volume, area = calculate_volume_area(vars_sol, vars_num; verbose=true);

kappas, kappap = find_curvature(vars_sol, vars_num);

shape_plt = plot_shape(vars_sol.r, vars_sol.z);
curv_plt = plot_curvature(vars_sol.z, kappas, kappap);
combined_plt = plot(shape_plt, curv_plt, layout=(1,2); show=true)