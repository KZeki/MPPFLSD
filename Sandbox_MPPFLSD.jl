using Pkg
Pkg.activate("BEP",shared=true)
# same as doing ] activate @BEP in REPL

# # Execute once then comment out
# path_dev_package_abslolute = raw"C:\Users\path_this_repo_was_coppied_to\MPPFLSD"
# Pkg.develop(path_dev_package_abslolute);
# Pkg.add(["Revise"; "Plots"]); # needs added packages

using Revise # allows for changes in packges to reflect without recompiling
#using Plots
using GLMakie
pt = 4/3
inch = 96
cm = inch/2.54 
fontsize_theme = Theme(fontsize = 10pt)
set_theme!(fontsize_theme)
using MPPFLSD
using Trapz
#using Debugger
#break_on(:error)

puff_params_Pinv = Dict(:puff_sigmoid_pmin=>-0.002, :puff_sigmoid_pmax=>43.47, :puff_sigmoid_a=>-4.49, :puff_sigmoid_b=>-729.56);

puff_params_Pc = Dict(:puff_sigmoid_pmin=>-0.61, :puff_sigmoid_pmax=>54.07, :puff_sigmoid_a=>-1.45, :puff_sigmoid_b=>-418.06);

puff_params_Pmeas = Dict(:puff_sigmoid_pmin=>-1.36, :puff_sigmoid_pmax=>46.31, :puff_sigmoid_a=>-2.55, :puff_sigmoid_b=>-512.02);
   

params_phys = ParamsPhys(;rneedle = 52e-3, volume0=0, deltarho=-(1050.0 - 1.225), grav=9.8, sigma = 0.02212,
    puffParams = Dict(:puff_rhs_method=>:puff_sigmoid, :impose_contact_angle=>true, :contact_angle=>0,
        :puff_exp_a=>-40, :puff_exp_b=>0.2,
        puff_params_Pinv... ));

params_num = ParamsNum(N=40);


# path = raw"C:\Users\zekig\OneDrive\Documenten\SchoolWork\BEP"
# using CSV, DataFrames
# tempA_names = ["tempA$i.csv" for i in 1:6 ]
# tempb_names = ["tempb$i.csv" for i in 1:6 ]
# A_df = CSV.read.(joinpath.(path,tempA_names),DataFrame,header=false);
# b_df = CSV.read.(joinpath.(path,tempb_names),DataFrame,header=false);
# global A = Matrix.(A_df);
# global b = vec.(Matrix.(b_df));
# global iteration = 1;

# solve for the droplet shape (Young-Laplace)
vars_sol, vars_num, _ = gen_single_drop(params_phys, params_num; verbose=true);

# post processing and plotting

volume, area = calculate_volume_area(vars_sol, vars_num; verbose=true);

kappas, kappap = find_curvature(vars_sol, vars_num);

begin
    #using GLMakie
    fig = Figure()
    ax = Axis(fig[1, -1:1], xlabel="r", ylabel="z")
    # Plot the main shape
    lines!(ax, vars_sol.r, vars_sol.z, label="z(r)", color=:blue)
    # Plot the green line from (0, y) to (rneedle, y)
    hline_y = -params_phys.volume0 / (π * params_phys.rneedle)
    lines!(ax, [0, params_phys.rneedle], [hline_y, hline_y], label="-H", color=:green)
    # Plot the red vertical line at rneedle
    lines!(ax, [params_phys.rneedle, params_phys.rneedle], [0, hline_y], label="-Rd", color=:red)
    # Set x limits
    xlims!(ax, -0.01, params_phys.rneedle * 1.5)
    # Compute volumes
    volume_between_0_and_z_of_r = trapz(vars_sol.r .^ 2 .* π, vars_sol.z)
    volume_between_bottom_dish_and_z_of_r = trapz(vars_sol.r .^ 2 .* π, vars_sol.z .- hline_y)
    # Title with computed volumes
    titletext = "volume between z(r) and 0 = $(round(volume_between_0_and_z_of_r; digits=4))\n" *
                "volume below z(r) = $(round(volume_between_bottom_dish_and_z_of_r; digits=4))"
    # Add title and legend
    axislegend(ax, position=:rb)
    Label(fig[0, :], titletext, fontsize=10)
    fig
end


#plot(vars_sol.r,vars_sol.z; xaxis=:log,yaxis=:log)
shape_plt = plot_shape(vars_sol.r, vars_sol.z, label="y1")
plot(vars_sol.r,vars_sol.z)
plot!([0; params_phys.rneedle], repeat([-params_phys.volume0/(pi*params_phys.rneedle)], 2), label="-H", color=:green)
plot!([params_phys.rneedle; params_phys.rneedle], [0; -params_phys.volume0/(pi*params_phys.rneedle)], label="-Rd", color=:red )
xlims!(-0.01, params_phys.rneedle*1.5)
volume_between_0_and_z_of_r = trapz(vars_sol.r.^2*pi,vars_sol.z)
volume_between_bottom_dish_and_z_of_r = trapz(vars_sol.r.^2*pi,(vars_sol.z .- -params_phys.volume0/(pi*params_phys.rneedle)))
titletext = "volume between z(r) and 0 = $(round(volume_between_0_and_z_of_r;digits=4)) \n volume below z(r) = $(round(volume_between_bottom_dish_and_z_of_r;digits=4))"
plot!(shape_plt,legend=:bottomright, plot_titlefontsize=8, plot_title=titletext)


curv_plt = plot_curvature(vars_sol.z, kappas, kappap);
#layout = Layout(width=800, height=600, margin=attr(l=65, r=50, b=65, t=90))
combined_plt = plot(shape_plt, curv_plt, layout=(2); show=true)
gr()
plotlyjs()


#p_r(r,Rd) = (Rd-r);

sigmoid(x) = 1 / (1 + exp(-x))
pmax = [ 43.47;54.07;46.31 ]
pmin = [-0.002;-0.61;1.36]
a = [-4.49;-1.45;-2.55]
b = [-729.56;-418.06;-512.02]
pr2(r) = pmin' .+ (pmax'.-pmin')./(1 .+exp.(a'.-b'.*r))
#pr(r) = 10 - (10-1)/params_phys.rneedle*r

fig_sigmoid = Figure()
ax_sigmoid = Axis(fig_sigmoid[1,1];xlabel="radius [mm]",ylabel="Pressure [Pa]",title=rich("Pressure profile p", subscript("puff"), "(r) for different parameters"))
sigmoid_r = LinRange(0,52e-3,1000)
lines_sigmoid_pinv = lines!(ax_sigmoid,sigmoid_r*1e3, pr2(sigmoid_r)[:,1]) 
lines_sigmoid_pc = lines!(ax_sigmoid,sigmoid_r*1e3, pr2(sigmoid_r)[:,2]) 
lines_sigmoid_pmeas = lines!(ax_sigmoid,sigmoid_r*1e3, pr2(sigmoid_r)[:,3])
axislegend(ax_sigmoid,[lines_sigmoid_pinv;lines_sigmoid_pc;lines_sigmoid_pmeas],[L"P_{inv}(r)";L"P_{c}(r)";L"P_{meas}(r)"];)