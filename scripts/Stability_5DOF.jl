using JLD2, GLMakie, Julianim, Statistics, MathTeXEngine
include("..\\src\\DynamicModels.jl")

path = abspath(@__DIR__, "..")
## Get data from static equilibrium
static_data = jldopen(joinpath(path, "data\\static\\Equilibrium_5DOF_final_reduced.jld2"), "r")
sol_static = static_data["solutions"]
pm_static = static_data["pm"]
close(static_data)
# Figure directory and general settings
f_name = "Stability_5DOF_final"
mkpath(joinpath(path, "figures\\" * f_name))
set_publication_theme!()
update_theme!(
fonts = Attributes(
        :bold => texfont(:bold),
        :bolditalic => texfont(:bolditalic),
        :italic => texfont(:italic),
        :regular => texfont(:regular)
    )
)

# Extract the data from static solution
T1_static = [sol.T1 for sol in sol_static]
T2_static = [sol.T2 for sol in sol_static]
N1_static = [sol.N1 for sol in sol_static]
N2_static = [sol.N2 for sol in sol_static]
xt1_static = [sol.xt1 for sol in sol_static]
xn1_static = [sol.xn1 for sol in sol_static]
wt1_static = [sol.wt1 for sol in sol_static]
xt2_static = [sol.xt2 for sol in sol_static]
xn2_static = [sol.xn2 for sol in sol_static]
wt2_static = [sol.wt2 for sol in sol_static]
n_static = length(T1_static)
u_static = zeros(5, n_static)
u_static[1, :] = [sol.psi1 for sol in sol_static]
u_static[2, :] = [sol.psi2 for sol in sol_static]
u_static[3, :] = [sol.x_D for sol in sol_static]
u_static[4, :] = [sol.y_D for sol in sol_static]
u_static[5, :] = [sol.psi_D for sol in sol_static]

Δ_T1 = (maximum(T1_static) - minimum(T1_static))/20
Δ_T2 = (maximum(T2_static) - minimum(T2_static))/20

# Parameters for the dynamic model
# ζ_\psi = 5.0e-2 / (2 * sqrt(4.0e3 * 6.348e-4))
n_steps = 200_001 # Number of time steps

pm = pm_static # Copy the static parameters to the dynamic model
pm.friction_law_1.wₜ = 0.0 # Reset wₜ so that previous history is lost
pm.friction_law_2.wₜ = 0.0
pm.T₁, pm.T₂, pm.N₁, pm.N₂, pm.xₜ₁, pm.xₜ₂, pm.xₙ₁, pm.xₙ₂ = ntuple(_ -> zeros(n_steps), 8)
pm.index = Ref(1) # Initialize index for the friction law

# Time vector
tf = 200.0
t = range(0.0, tf, length = n_steps)
Δt = t[2] - t[1]

# Mode
#    1 - stability -> Time integration for a series of static equilibrium solutions
#    2 - error -> Error (final value of time integration) analysis for all static 
#                 equilibrium solutions, reduced solutions recommended
#    3 - evolution -> Time integration with a perturbation or a forcing term for
#                     a series of static equilibrium solutions
mode = 1

# For stability and evolution modes
static_indices = [7 91791 533652 734134 1024500] # Indices for the static equilibrium solutions
# static_indices = 1:length(T1_static) # WARNING!!! This selects all static equilibrium solutions
# It is recommended to use the "plotting/Point_selection.jl" script to get
# the indices of the desired static equilibrium solutions

# For evolution mode
ϵ₀ = 0.1 # Perturbation for the initial condition
Mₑ = nothing # Forcing momentum

if mode == 1 # Stability analysis
    n_iter = length(static_indices)

    u = zeros(10, n_steps, n_iter)
    N₁ = zeros(n_steps-1, n_iter)
    T₁ = zeros(n_steps-1, n_iter)
    N₂ = zeros(n_steps-1, n_iter)
    T₂ = zeros(n_steps-1, n_iter)

    # Initialize figures for plotting
    fig_psi1 = Figure(size = (1280, 720))
    ax_psi1 = Axis(fig_psi1[1, 1], 
        # title=L"\psi_1 \text{ vs time}", 
        xlabel=L"\tau", 
        ylabel=L"\psi_1"
    )
    lines!(ax_psi1, 
        [t[1], t[end]], 
        [NaN, NaN]; 
        color = :black,
        linestyle = :dash,
        linewidth = 5,
        label = "Static solution"
    )
    lines!(ax_psi1, 
        [t[1], t[end]], 
        [NaN, NaN]; 
        color = :black,
        label = "Time integration"
    )

    fig_psi2 = Figure(size = (1280, 720))
    ax_psi2 = Axis(fig_psi2[1, 1], 
        # title=L"\psi_2 \text{ vs time}", 
        xlabel=L"\tau", 
        ylabel=L"\psi_2"; 
    )
    lines!(ax_psi2, 
        [t[1], t[end]], 
        [NaN, NaN]; 
        color = :black,
        linestyle = :dash,
        linewidth = 5,
        label = "Static solution"
    )
    lines!(ax_psi2, 
        [t[1], t[end]], 
        [NaN, NaN]; 
        color = :black,
        label = "Time integration"
    )

    fig_x_D = Figure(size = (1280, 720))
    ax_x_D = Axis(fig_x_D[1, 1], 
        # title=L"x_D \text{ vs time}", 
        xlabel=L"\tau", 
        ylabel=L"x_D"; 
    )
    lines!(ax_x_D, 
        [t[1], t[end]], 
        [NaN, NaN]; 
        color = :black,
        linestyle = :dash,
        linewidth = 5,
        label = "Static solution"
    )
    lines!(ax_x_D, 
        [t[1], t[end]], 
        [NaN, NaN]; 
        color = :black,
        label = "Time integration"
    )

    fig_y_D = Figure(size = (1280, 720))
    ax_y_D = Axis(fig_y_D[1, 1], 
        # title=L"y_D \text{ vs time}", 
        xlabel=L"\tau", 
        ylabel=L"y_D"; 
    )
    lines!(ax_y_D, 
        [t[1], t[end]], 
        [NaN, NaN]; 
        color = :black,
        linestyle = :dash,
        linewidth = 5,
        label = "Static solution"
    )
    lines!(ax_y_D, 
        [t[1], t[end]], 
        [NaN, NaN]; 
        color = :black,
        label = "Time integration"
    )

    fig_psi_D = Figure(size = (1280, 720))
    ax_psi_D = Axis(fig_psi_D[1, 1], 
        # title=L"\psi_D \text{ vs time}", 
        xlabel=L"\tau", 
        ylabel=L"\psi_D"; 
    )
    lines!(ax_psi_D, 
        [t[1], t[end]], 
        [NaN, NaN]; 
        color = :black,
        linestyle = :dash,
        linewidth = 5,
        label = "Static solution"
    )
    lines!(ax_psi_D, 
        [t[1], t[end]], 
        [NaN, NaN]; 
        color = :black,
        label = "Time integration"
    )

    fig_N1 = Figure(size = (1280, 720))
    ax_N1 = Axis(fig_N1[1, 1], 
        # title=L"N_1 \text{ vs time}", 
        xlabel=L"\tau", 
        ylabel=L"N_1"; 
    )
    lines!(ax_N1, 
        [t[1], t[end]], 
        [NaN, NaN]; 
        color = :black,
        linestyle = :dash,
        linewidth = 5,
        label = "Static solution"
    )
    lines!(ax_N1, 
        [t[1], t[end]], 
        [NaN, NaN]; 
        color = :black,
        label = "Time integration"
    )

    fig_T1 = Figure(size = (1280, 720))
    ax_T1 = Axis(fig_T1[1, 1], 
        # title=L"T_1 \text{ vs time}", 
        xlabel=L"\tau", 
        ylabel=L"T_1"; 
    )
    lines!(ax_T1, 
        [t[1], t[end]], 
        [NaN, NaN]; 
        color = :black,
        linestyle = :dash,
        linewidth = 5,
        label = "Static solution"
    )
    lines!(ax_T1, 
        [t[1], t[end]], 
        [NaN, NaN]; 
        color = :black,
        label = "Time integration"
    )

    fig_N2 = Figure(size = (1280, 720))
    ax_N2 = Axis(fig_N2[1, 1], 
        # title=L"N_2 \text{ vs time}", 
        xlabel=L"\tau", 
        ylabel=L"N_2"; 
    )
    lines!(ax_N2, 
        [t[1], t[end]], 
        [NaN, NaN]; 
        color = :black,
        linestyle = :dash,
        linewidth = 5,
        label = "Static solution"
    )
    lines!(ax_N2, 
        [t[1], t[end]], 
        [NaN, NaN]; 
        color = :black,
        label = "Time integration"
    )

    fig_T2 = Figure(size = (1280, 720))
    ax_T2 = Axis(fig_T2[1, 1], 
        # title=L"T_1 \text{ vs time}", 
        xlabel=L"\tau", 
        ylabel=L"T_2"; 
    )
    lines!(ax_T2, 
        [t[1], t[end]], 
        [NaN, NaN]; 
        color = :black,
        linestyle = :dash,
        linewidth = 5,
        label = "Static solution"
    )
    lines!(ax_T2, 
        [t[1], t[end]], 
        [NaN, NaN]; 
        color = :black,
        label = "Time integration"
    )

    # Initial condition plot
    fig_initial_cond = Figure()
    ax_initial_cond = Axis(fig_initial_cond[1, 1], 
        # title=L"Static solution positions", 
        xlabel=L"T_1", 
        ylabel=L"T_2", 
        aspect = 1
    )
    co_initial_cond = contourf!(ax_initial_cond, 
        T1_static, 
        T2_static, 
        ones(n_static);
        colormap = :grays,
    )
    limits!(ax_initial_cond, 
        minimum(T1_static)-Δ_T1, maximum(T1_static)+Δ_T1, 
        minimum(T2_static)-Δ_T2, maximum(T2_static)+Δ_T2
    )

    for (i, idx_static) in enumerate(static_indices)
        # Set initial condition from static equilibrium solution
        u[:, 1, i] = vcat(u_static[:, idx_static], zeros(5)) 
        
        # Set wₜ so that the initial T is consistent with the static equilibrium solution
        pm.friction_law_1.wₜ = wt1_static[idx_static]
        pm.friction_law_2.wₜ = wt2_static[idx_static]
        for j in 1:(n_steps-1)
            u[:, j+1, i] = rk4(five_dof_model, u[:, j, i], pm, t[j], Δt)
        end
        N₁[:, i] = pm.N₁[1:end-1]
        T₁[:, i] = pm.T₁[1:end-1]
        N₂[:, i] = pm.N₂[1:end-1]
        T₂[:, i] = pm.T₂[1:end-1]
        pm.index = Ref(1) # Reset the index for the next iteration
        pm.friction_law_1.wₜ = 0.0 # Reset wₜ so that previous history is lost
        pm.friction_law_2.wₜ = 0.0 

        # Plot the results
        lines!(ax_psi1, 
            t, 
            u[1, :, i]; 
            color = COLORS[i],
            label = "Time integration"
        )
        lines!(ax_psi1, 
            [t[1], t[end]], 
            u_static[1, idx_static].*[1, 1]; 
            color = COLORS[i], 
            linestyle = :dash,
            linewidth = 5,
            label = "Static solution"
        )

        lines!(ax_psi2, 
            t, 
            u[2, :, i]; 
            color = COLORS[i],
            label = "Time integration"
        )
        lines!(ax_psi2, 
            [t[1], t[end]],
            u_static[2, idx_static].*[1, 1]; 
            color = COLORS[i], 
            linestyle = :dash,
            linewidth = 5,
            label = "Static solution"
        )

        lines!(ax_x_D, 
            t, 
            u[3, :, i]; 
            color = COLORS[i],
            label = "Time integration"
        )
        lines!(ax_x_D, 
            [t[1], t[end]], 
            u_static[3, idx_static].*[1, 1]; 
            color = COLORS[i], 
            linestyle = :dash,
            linewidth = 5,
            label = "Static solution"
        )

        lines!(ax_y_D, 
            t, 
            u[4, :, i]; 
            color = COLORS[i],
            label = "Time integration"
        )
        lines!(ax_y_D, 
            [t[1], t[end]], 
            u_static[4, idx_static].*[1, 1]; 
            color = COLORS[i], 
            linestyle = :dash,
            linewidth = 5,
            label = "Static solution"
        )

        lines!(ax_psi_D, 
            t, 
            u[5, :, i]; 
            color = COLORS[i],
            label = "Time integration"
        )
        lines!(ax_psi_D, 
            [t[1], t[end]], 
            u_static[5, idx_static].*[1, 1]; 
            color = COLORS[i], 
            linestyle = :dash,
            linewidth = 5,
            label = "Static solution"
        )

        lines!(ax_N1, 
            t[1:end-1], 
            N₁[:, i]; 
            color = COLORS[i],
            label = "Time integration"
        )
        lines!(ax_N1, 
            [t[1], t[end-1]], 
            N1_static[idx_static].*[1, 1]; 
            color = COLORS[i], 
            linestyle = :dash,
            linewidth = 5,
            label = "Static solution"
        )

        lines!(ax_T1, 
            t[1:end-1], 
            T₁[:, i]; 
            color = COLORS[i],
            label = "Time integration"
        )
        lines!(ax_T1, 
            [t[1], t[end-1]], 
            T1_static[idx_static].*[1, 1]; 
            color = COLORS[i], 
            linestyle = :dash,
            linewidth = 5,
            label = "Static solution"
        )

        lines!(ax_N2, 
            t[1:end-1], 
            N₂[:, i]; 
            color = COLORS[i],
            label = "Time integration"
        )
        lines!(ax_N2, 
            [t[1], t[end-1]], 
            N2_static[idx_static].*[1, 1]; 
            color = COLORS[i], 
            linestyle = :dash,
            linewidth = 5,
            label = "Static solution"
        )

        lines!(ax_T2, 
            t[1:end-1], 
            T₂[:, i]; 
            color = COLORS[i],
            label = "Time integration"
        )
        lines!(ax_T2, 
            [t[1], t[end-1]], 
            T2_static[idx_static].*[1, 1]; 
            color = COLORS[i], 
            linestyle = :dash,
            linewidth = 5,
            label = "Static solution"
        )

        scatter!(ax_initial_cond, 
            T1_static[idx_static], 
            T2_static[idx_static]; 
            color = COLORS[i], 
            marker = :circle,
            markersize = 25,
        )
    end
    legend_psi1 = axislegend(ax_psi1, unique = true)
    translate!(legend_psi1.blockscene, Vec3f(0, -40, 0))
    legend_psi2 = axislegend(ax_psi2, unique = true)
    translate!(legend_psi2.blockscene, Vec3f(0, -40, 0))
    legend_x_D = axislegend(ax_x_D, unique = true)
    translate!(legend_x_D.blockscene, Vec3f(0, -40, 0))
    legend_y_D = axislegend(ax_y_D, unique = true)
    translate!(legend_y_D.blockscene, Vec3f(0, -40, 0))
    legend_psi_D = axislegend(ax_psi_D, unique = true)
    translate!(legend_psi_D.blockscene, Vec3f(0, -40, 0))
    legend_N1 = axislegend(ax_N1, unique = true)
    translate!(legend_N1.blockscene, Vec3f(0, -40, 0))
    legend_T1 = axislegend(ax_T1, unique = true)
    translate!(legend_T1.blockscene, Vec3f(0, -40, 0))
    legend_N2 = axislegend(ax_N2, unique = true)
    translate!(legend_N2.blockscene, Vec3f(0, -40, 0))
    legend_T2 = axislegend(ax_T2, unique = true)
    translate!(legend_T2.blockscene, Vec3f(0, -40, 0))

    # Save figures
    save(joinpath(path, "figures\\" * f_name * "\\psi1.png"), fig_psi1, px_per_unit = 4)
    save(joinpath(path, "figures\\" * f_name * "\\psi2.png"), fig_psi2, px_per_unit = 4)
    save(joinpath(path, "figures\\" * f_name * "\\x_D.png"), fig_x_D, px_per_unit = 4)
    save(joinpath(path, "figures\\" * f_name * "\\y_D.png"), fig_y_D, px_per_unit = 4)
    save(joinpath(path, "figures\\" * f_name * "\\psi_D.png"), fig_psi_D, px_per_unit = 4)
    save(joinpath(path, "figures\\" * f_name * "\\N1.png"), fig_N1, px_per_unit = 4)
    save(joinpath(path, "figures\\" * f_name * "\\T1.png"), fig_T1, px_per_unit = 4)
    save(joinpath(path, "figures\\" * f_name * "\\N2.png"), fig_N2, px_per_unit = 4)
    save(joinpath(path, "figures\\" * f_name * "\\T2.png"), fig_T2, px_per_unit = 4)
    save(joinpath(path, "figures\\" * f_name * "\\initial_cond.png"), fig_initial_cond, px_per_unit = 4)

elseif mode == 2 # Error analysis
    ϵ = zeros(5, n_static)
    RMS = zeros(n_static)
    u = zeros(10, n_steps)
    for i in 1:n_static
        # Set initial condition from static equilibrium solution
        u[:, 1] = vcat(u_static[:, i], zeros(5))

        # Set wₜ so that the initial T is consistent with the static equilibrium solution
        pm.friction_law_1.wₜ = wt1_static[i]
        pm.friction_law_2.wₜ = wt2_static[i]
        for j in 1:(n_steps-1)
            u[:, j+1] = rk4(five_dof_model, u[:, j], pm, t[j], Δt)
        end
        pm.index = Ref(1) # Reset the index for the next iteration
        pm.friction_law_1.wₜ = 0.0 # Reset wₜ so that previous history is lost
        pm.friction_law_2.wₜ = 0.0 

        # Get the errors and RMS when converged
        ϵ2 = 0.0
        for k in 1:5
            avg_u = mean(u[k, round(Int, 0.9 * n_steps):end])
            ϵ[k, i] = avg_u - u_static[k, i]
            ϵ2 += ϵ[k, i]^2
        end
        RMS[i] = sqrt(ϵ2 / 5)
    end
    # Initialize figures for plotting
    fig_error_psi1 = Figure()
    ax_error_psi1 = Axis(fig_error_psi1[1, 1], 
        # title=L"|ϵ_\psi_1| vs T₁ & T₂", 
        xlabel=L"T₁", 
        ylabel=L"T₂", 
        aspect = 1
    )
    co_error_psi1 = contourf!(ax_error_psi1, 
        T1_static, 
        T2_static, 
        abs.(ϵ[1, 1:end-1])
    )
    Colorbar(fig_error_psi1[1, 2], co_error_psi1)
    save(joinpath(path, "figures\\" * f_name * "\\error_psi1.png"), fig_error_psi1)

    fig_error_psi2 = Figure()
    ax_error_psi2 = Axis(fig_error_psi2[1, 1], 
        # title=L"|ϵ_\psi₂| vs T₁ & T₂", 
        xlabel=L"T₁", 
        ylabel=L"T₂", 
        aspect = 1
    )
    co_error_psi2 = contourf!(ax_error_psi2,
        T1_static, 
        T2_static, 
        abs.(ϵ[2, 1:end-1])
    )
    Colorbar(fig_error_psi2[1, 2], co_error_psi2)
    save(joinpath(path, "figures\\" * f_name * "\\error_psi2.png"), fig_error_psi2)

    fig_error_x_D = Figure()
    ax_error_x_D = Axis(fig_error_x_D[1, 1], 
        # title=L"|ϵ_x_D| vs T₁ & T₂", 
        xlabel=L"T₁", 
        ylabel=L"T₂", 
        aspect = 1
    )
    co_error_x_D = contourf!(ax_error_x_D, 
        T1_static, 
        T2_static, 
        abs.(ϵ[3, 1:end-1])
    )
    Colorbar(fig_error_x_D[1, 2], co_error_x_D)
    save(joinpath(path, "figures\\" * f_name * "\\error_x_D.png"), fig_error_x_D)

    fig_error_y_D = Figure()
    ax_error_y_D = Axis(fig_error_y_D[1, 1], 
        # title=L"|ϵ_y_D| vs T₁ & T₂", 
        xlabel=L"T₁", 
        ylabel=L"T₂", 
        aspect = 1
    )
    co_error_y_D = contourf!(ax_error_y_D, 
        T1_static, 
        T2_static, 
        abs.(ϵ[4, 1:end-1])
    )
    Colorbar(fig_error_y_D[1, 2], co_error_y_D)
    save(joinpath(path, "figures\\" * f_name * "\\error_y_D.png"), fig_error_y_D)

    fig_error_psi_D = Figure()
    ax_error_psi_D = Axis(fig_error_psi_D[1, 1], 
        # title=L"|\varepsilon_\psi_D| vs T₁ & T₂", 
        xlabel=L"T₁", 
        ylabel=L"T₂", 
        aspect = 1
    )
    co_error_psi_D = contourf!(ax_error_psi_D, 
        T1_static, 
        T2_static, 
        abs.(ϵ[5, 1:end-1])
    )
    Colorbar(fig_error_psi_D[1, 2], co_error_psi_D)
    save(joinpath(path, "figures\\" * f_name * "\\error_psi_D.png"), fig_error_psi_D)

    fig_rms = Figure()
    ax_rms = Axis(fig_rms[1, 1], 
        # title=L"RMS vs T₁ & T₂", 
        xlabel=L"T₁", 
        ylabel=L"T₂", 
        aspect = 1
    )
    co_rms = contourf!(ax_rms, 
        T1_static, 
        T2_static, 
        RMS[1:end-1]
    )
    Colorbar(fig_rms[1, 2], co_rms)
    save(joinpath(path, "figures\\" * f_name * "\\RMS.png"), fig_rms)

elseif mode == 3 # Evolution analysis
    # Initialize figures for plotting
    fig_evolution = Figure()
    ax_evolution = Axis(fig_evolution[1, 1], 
        # title=L"\text{Evolution of } T_1 \text{ & } T_2", 
        xlabel=L"T_1", 
        ylabel=L"T_2", 
        aspect = 1,
    )
    co_evolution = contourf!(ax_evolution ,
        T1_static, 
        T2_static, 
        0.7.*ones(n_static); 
        colormap = :grays, 
        levels = range(0, 1; length = 10)   
    )
    limits!(ax_evolution, 
        minimum(T1_static)-Δ_T1, maximum(T1_static)+Δ_T1, 
        minimum(T2_static)-Δ_T2, maximum(T2_static)+Δ_T2
    )

    fig_arrows = Figure()
    ax_arrows = Axis(fig_arrows[1, 1], 
        # title=L"\text{Evolution of } T_1 \text{ & } T_2", 
        xlabel=L"T_1", 
        ylabel=L"T_2", 
        aspect = 1,
    )
    co_arrows = contourf!(ax_arrows ,
        T1_static, 
        T2_static, 
        0.7.*ones(n_static); 
        colormap = :grays, 
        levels = range(0, 1; length = 10)   
    )
    limits!(ax_arrows, 
        minimum(T1_static)-Δ_T1, maximum(T1_static)+Δ_T1, 
        minimum(T2_static)-Δ_T2, maximum(T2_static)+Δ_T2
    )

    u = zeros(10, n_steps)
    n_final = round(Int, 0.9 * n_steps)

    pm.index = Ref(1)
    for idx_static in static_indices
        # Set initial condition from static equilibrium solution´
        u₀ = u_static[:, idx_static]
        if ϵ₀ !== nothing # Perturbation
            u₀ .+= ϵ₀ .* u₀ .* (2 .* rand(length(u₀)) .- 1)
        end

        # Forcing term Mₑ
        if Mₑ !== nothing
            pm.Mₑ = Mₑ
            pm.feq = 1.0
        end

        u[:, 1] = vcat(u₀, zeros(5))

        # Set wₜ so that the initial T is consistent with the static equilibrium solution
        pm.friction_law_1.wₜ = wt1_static[idx_static]
        pm.friction_law_2.wₜ = wt2_static[idx_static]
        for j in 1:(n_steps-1)
            u[:, j+1] = rk4(five_dof_model, u[:, j], pm, t[j], Δt)
        end
        pm.index = Ref(1) # Reset the index for the next iteration
        pm.friction_law_1.wₜ = 0.0 # Reset wₜ so that previous history is lost
        pm.friction_law_2.wₜ = 0.0 

        T₁_final = mean(pm.T₁[n_final:end-1])
        T₂_final = mean(pm.T₂[n_final:end-1])

        # Plot the results
        lines!(ax_evolution,
            pm.T₁[1:end-1], 
            pm.T₂[1:end-1]; 
            color = :black
        )
        scatter!(ax_evolution, 
            T₁_final, 
            T₂_final; 
            color = :red, 
            markersize = 10,
            marker = :circle
        )

        ini_pt = Point2f(T1_static[idx_static], T2_static[idx_static])
        fin_pt = Point2f(T₁_final,T₂_final)
        dir = fin_pt - ini_pt
        arrows!(ax_arrows, 
            [ini_pt], 
            [dir];
            color = :black
        )
    end
    display(GLMakie.Screen(), fig_evolution)
    display(GLMakie.Screen(), fig_arrows)
else 
    error("Invalid mode $mode. Please choose 1 (stability) / 2 (error) / 3 (evolution).")
end