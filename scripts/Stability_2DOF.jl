using JLD2, Parameters, GLMakie
include("..\\src\\DynamicModels.jl")

path = abspath(@__DIR__, "..")
## Get data from static equilibrium
freetext = "_test_storage"
static_data = jldopen(joinpath(path, "data\\static\\Equilibrium_2DOF" * freetext *".jld2"), "r")
sol_static = static_data["solutions"]
pm_static = static_data["pm"]
close(static_data)

# Extract the data for plotting
psi_static = [sol.psi for sol in sol_static]
y_D_static = [sol.y_D for sol in sol_static]
N_static = [sol.N for sol in sol_static]
T_static = [sol.T for sol in sol_static]
xt_static = [sol.xt for sol in sol_static]
xn_static = [sol.xn for sol in sol_static]
wt_static = [sol.wt for sol in sol_static]

# Parameters for the dynamic model
n_steps = 2_000_001 # Number of time steps

pm = pm_static # Copy the static parameters to the dynamic model
pm.friction_law.wₜ = 0.0 # Reset wₜ so that previous history is lost
pm.T, pm.N, pm.xₜ, pm.xₙ = ntuple(_ -> zeros(n_steps), 4) # Initialize T, N, xₜ, xₙ
pm.index = Ref(1) # Initialize index for the friction law

n_static::Int  = length(T_static)
distance::Int = 100 # Distance between the static equilibrium solutions
n_iter::Int = floor(Int, ((n_static - 1)  / distance) + 1)

# Time vector
tf = 200.0
t = range(0.0, tf, length = n_steps)
dt = t[2] - t[1]
u = zeros(4, n_steps, n_iter)
T = zeros(n_steps-1, n_iter)
N = zeros(n_steps-1, n_iter)
# Initialize figures for plotting
fig_psi = Figure(size = (1280, 720))
ax_psi = Axis(fig_psi[1, 1], title="ψ vs time", xlabel="t", ylabel="ψ"; xgridvisible = true, ygridvisible = true, xgridcolor = :gray, ygridcolor = :gray)
fig_y_D = Figure(size = (1280, 720))
ax_y_D = Axis(fig_y_D[1, 1], title="y_D vs time", xlabel="t", ylabel="y_D"; xgridvisible = true, ygridvisible = true, xgridcolor = :gray, ygridcolor = :gray)
fig_N = Figure(size = (1280, 720))
ax_N = Axis(fig_N[1, 1], title="N vs time", xlabel="t", ylabel="N"; xgridvisible = true, ygridvisible = true, xgridcolor = :gray, ygridcolor = :gray)
fig_T = Figure(size = (1280, 720))
ax_T = Axis(fig_T[1, 1], title="T vs time", xlabel="t", ylabel="T"; xgridvisible = true, ygridvisible = true, xgridcolor = :gray, ygridcolor = :gray)
# Initial condition plot
fig_initial_cond = Figure(size = (1280, 720))
ax_initial_cond = Axis(fig_initial_cond[1, 1], title="ψ vs T", xlabel="T", ylabel="ψ")
lines!(ax_initial_cond, T_static, psi_static, linewidth = 2)
for i in 1:n_iter
    index_static::Int = round(Int, (i - 1) * distance + 1)
    if index_static > n_static
        break
    end

    # Set initial condition from static equilibrium solution
    u[:, 1, i] = [psi_static[index_static], y_D_static[index_static], 0.0, 0.0]

    # Set wₜ so that the initial T is consistent with the static equilibrium solution
    pm.friction_law.wₜ = wt_static[index_static]
    for j in 1:(n_steps-1)
        u[:, j+1, i] = rk4(two_dof_model, u[:, j, i], pm, t[j], dt)
    end
    N[:, i] = pm.N[1:end-1]
    T[:, i] = pm.T[1:end-1]
    pm.index = Ref(1) # Reset the index for the next iteration
    pm.friction_law.wₜ = 0.0 # Reset wₜ so that previous history is lost


    # Plot the results
    lines!(ax_psi, t, u[1, :, i]; linewidth = 2, color = :blue)
    lines!(ax_psi, [t[1], t[end]], psi_static[index_static].*[1,1]; color = :red, linestyle = :dash)
    lines!(ax_y_D, t, u[2, :, i]; linewidth = 2, color = :blue)
    lines!(ax_y_D, [t[1], t[end]], y_D_static[index_static].*[1,1]; color = :red, linestyle = :dash)
    lines!(ax_N, t[1:end-1], N[:, i]; linewidth = 2, color = :blue)
    lines!(ax_N, [t[1], t[end-1]], N_static[index_static].*[1,1]; color = :red, linestyle = :dash)
    lines!(ax_T, t[1:end-1], T[:, i]; linewidth = 2, color = :blue)
    lines!(ax_T, [t[1], t[end-1]], T_static[index_static].*[1,1]; color = :red, linestyle = :dash)
    scatter!(ax_initial_cond, T_static[index_static], psi_static[index_static]; markersize = 10, color = :red)
end

# Save figures
mkpath(joinpath(path, "figures/Stability_2DOF" * freetext))
save(joinpath(path, "figures/Stability_2DOF" * freetext * "/psi.png"), fig_psi)
save(joinpath(path, "figures/Stability_2DOF" * freetext * "/y_D.png"), fig_y_D)
save(joinpath(path, "figures/Stability_2DOF" * freetext * "/N.png"), fig_N)
save(joinpath(path, "figures/Stability_2DOF" * freetext * "/T.png"), fig_T)
save(joinpath(path, "figures/Stability_2DOF" * freetext * "/initial_cond.png"), fig_initial_cond)