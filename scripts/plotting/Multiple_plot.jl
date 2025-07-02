using JLD2, GLMakie, Statistics
include("..\\..\\src\\Common.jl")

path = abspath(@__DIR__, "..\\..")
file_name = "TBMS_coulomb_progressive_loading"
mkpath("julia/figures/"*file_name)
# Load files
file_list = ["TBMS_tf200_dt1e-4_rk4_coulomb_new_parameters2", "TBMS_tf200_dt1e-4_rk4_coulomb_progressive_FC2", "TBMS_tf200_dt1e-4_rk4_coulomb_progressive_FC1"]
file_descriptions = ["ref", "tanh(t/50)", "tanh(t/100)"]
@assert length(file_list) == length(file_descriptions)

## DOF Plots
# Create a figure for each variable
fig_psi1 = Figure(size = (1920, 720))
ax_psi1 = Axis(fig_psi1[1, 1], title="ψ₁ vs time", xlabel="t", ylabel="ψ₁"; xgridvisible = true, ygridvisible = true, xgridcolor = :gray, ygridcolor = :gray)
fig_psi2 = Figure(size = (1920, 720))
ax_psi2 = Axis(fig_psi2[1, 1], title="ψ₂ vs time", xlabel="t", ylabel="ψ₂"; xgridvisible = true, ygridvisible = true, xgridcolor = :gray, ygridcolor = :gray)
fig_xD = Figure(size = (1920, 720))
ax_xD = Axis(fig_xD[1, 1], title="x_D vs time", xlabel="t", ylabel="x_D"; xgridvisible = true, ygridvisible = true, xgridcolor = :gray, ygridcolor = :gray)
fig_yD = Figure(size = (1920, 720))
ax_yD = Axis(fig_yD[1, 1], title="y_D vs time", xlabel="t", ylabel="y_D"; xgridvisible = true, ygridvisible = true, xgridcolor = :gray, ygridcolor = :gray)

for file in file_list
    data = load("julia/data/"*file*".jld2")
    u = data["u"]
    t = data["t"]

    # Plot ψ₁
    lines!(ax_psi1, t, u[1, :], label = file_descriptions[file_list .== file][1])
    # Plot ψ₂
    lines!(ax_psi2, t, u[2, :], label = file_descriptions[file_list .== file][1])
    # Plot x_D
    lines!(ax_xD, t, u[3, :], label = file_descriptions[file_list .== file][1])
    # Plot y_D
    lines!(ax_yD, t, u[4, :], label = file_descriptions[file_list .== file][1])
end

# Save the figures
axislegend(ax_psi1)
save("julia/figures/"*file_name*"/psi_1.png", fig_psi1)
axislegend(ax_psi2)
save("julia/figures/"*file_name*"/psi_2.png", fig_psi2)
axislegend(ax_xD)
save("julia/figures/"*file_name*"/x_D.png", fig_xD)
axislegend(ax_yD)
save("julia/figures/"*file_name*"/y_D.png", fig_yD)

## SCALED DOF Plots
# Create a figure for each variable
fig_psi1_scaled = Figure(size = (1920, 720))
ax_psi1_scaled = Axis(fig_psi1_scaled[1, 1], title="ψ₁/ψ₁_conv vs time", xlabel="t", ylabel="ψ₁"; xgridvisible = true, ygridvisible = true, xgridcolor = :gray, ygridcolor = :gray)
fig_psi2_scaled = Figure(size = (1920, 720))
ax_psi2_scaled = Axis(fig_psi2_scaled[1, 1], title="ψ₂/ψ₂_conv vs time", xlabel="t", ylabel="ψ₂"; xgridvisible = true, ygridvisible = true, xgridcolor = :gray, ygridcolor = :gray)
fig_xD_scaled = Figure(size = (1920, 720))
ax_xD_scaled = Axis(fig_xD_scaled[1, 1], title="x_D/x_D_conv vs time", xlabel="t", ylabel="x_D"; xgridvisible = true, ygridvisible = true, xgridcolor = :gray, ygridcolor = :gray)
fig_yD_scaled = Figure(size = (1920, 720))
ax_yD_scaled = Axis(fig_yD_scaled[1, 1], title="y_D/y_D_conv vs time", xlabel="t", ylabel="y_D"; xgridvisible = true, ygridvisible = true, xgridcolor = :gray, ygridcolor = :gray)

for file in file_list
    data = load("julia/data/"*file*".jld2")
    u = data["u"]
    t = data["t"]
    
    # Get the final value of each variable
    conv_psi1 = mean(u[1, round(Int, 0.9 * size(u, 2)):end])
    conv_psi2 = mean(u[2, round(Int, 0.9 * size(u, 2)):end])
    conv_xD = mean(u[3, round(Int, 0.9 * size(u, 2)):end])
    conv_yD = mean(u[4, round(Int, 0.9 * size(u, 2)):end])

    # Plot ψ₁
    lines!(ax_psi1_scaled, t, u[1, :]./conv_psi1, label = file_descriptions[file_list .== file][1])
    # Plot ψ₂
    lines!(ax_psi2_scaled, t, u[2, :]./conv_psi2, label = file_descriptions[file_list .== file][1])
    # Plot x_D
    lines!(ax_xD_scaled, t, u[3, :]./conv_xD, label = file_descriptions[file_list .== file][1])
    # Plot y_D
    lines!(ax_yD_scaled, t, u[4, :]./conv_yD, label = file_descriptions[file_list .== file][1])
end

# Save the figures
axislegend(ax_psi1_scaled)
save("julia/figures/"*file_name*"/psi_1_scaled.png", fig_psi1_scaled)
axislegend(ax_psi2_scaled)
save("julia/figures/"*file_name*"/psi_2_scaled.png", fig_psi2_scaled)
axislegend(ax_xD_scaled)
save("julia/figures/"*file_name*"/x_D_scaled.png", fig_xD_scaled)
axislegend(ax_yD_scaled)
save("julia/figures/"*file_name*"/y_D_scaled.png", fig_yD_scaled)