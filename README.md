
# Measuring Puff Pressure From Liquid Surface Deformations
This package ports the functionality to run the simple case of TensioMet2 (Matlab library: https://github.com/NickJaensson/tensiomet2) into Julia, and is used to generate pendant droplets by solving the Young-Laplace equation using spectral methods.

This is then used to generate the (steady state) Lidquid Surface Deformation caused by that applying a pressure to a fluid. 

### Example running simple case of pendant drop

```julia
using Pkg
Pkg.activate("BEP",shared=true)
# same as doing ] activate @BEP in REPL

# # Execute once then comment out
# path_dev_package_abslolute = raw"C:\Users\path_this_repo_was_coppied_to\MPPFLSD"
# Pkg.develop(path_dev_package_abslolute);
# Pkg.add(["Revise"; "Plots]); # might want to add Revise for developing, and additional 

using Revise # allows for changes in packges to reflect without recompiling
using Plots
using MPPFLSD

# Assign structs containing the parameters you would want to solve for
params_phys = ParamsPhys();
params_num = ParamsNum();

# Solve for the droplet shape (Young-Laplace)
vars_sol, vars_num = gen_single_drop(params_phys, params_num; verbose=false)

# Post processing and plotting
volume, area = calculate_volume_area(vars_sol, vars_num; verbose=true);
kappas, kappap = find_curvature(vars_sol, vars_num);

# Plot results
shape_plt = plot_shape(vars_sol.r, vars_sol.z);
curv_plt = plot_curvature(vars_sol.z, kappas, kappap);
combined_plt = plot(shape_plt, curv_plt, layout=(1,2); show=true)

# Could also have run the function: example_simple()

```