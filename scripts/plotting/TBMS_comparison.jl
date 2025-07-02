using JLD2, GLMakie
include("..\\..\\src\\Common.jl")

path = abspath(@__DIR__, "..\\..")
f_name_1 = "TBMS_tf200_dt1e-4_rk4_linear_ref1"
label_1 = "ref"
f_name_2 = "TBMS_tf200_dt1e-4_rk4_linear_progressive_FC1"
label_2 = "progressive loading"

data_1 = load(joinpath(path, "data\\" * f_name_1 * ".jld2"))
data_2 = load(joinpath(path, "data\\" * f_name_2 * ".jld2"))

u_1 = data_1["u"]
u_2 = data_2["u"]
t_1 = data_1["t"]
t_2 = data_2["t"]
pm_1 = data_1["pm"]
pm_2 = data_2["pm"]

f_name = "hola"

mkpath(joinpath(path, "figures\\" * f_name))

print("Plot DOF? [Y/n]: ")
response = readline() |> strip
if response == "Y" || response == "y"
    fig1 = Figure(size = (1920, 720))
    ax1 = Axis(fig1[1, 1], title="ψ₁ vs time", xlabel="t", ylabel="ψ₁"; xgridvisible = true, ygridvisible = true, xgridcolor = :gray, ygridcolor = :gray)
    lines!(t_1, u_1[1, :], label = label_1)
    lines!(t_2, u_2[1, :], label = label_2)
    axislegend(ax1)
    # limits!(ax1, minimum(t), maximum(t), minimum(u[1, :]) - 1e-7*abs(minimum(u[1, :])), maximum(u[1, :]) + 1e-7*abs(maximum(u[1, :])))
    save(joinpath(path, "figures\\" * f_name * "\\psi_1.png"), fig1)

    fig2 = Figure(size = (1920, 720))
    ax2 = Axis(fig2[1, 1], title="ψ₂ vs time", xlabel="t", ylabel="ψ₂"; xgridvisible = true, ygridvisible = true, xgridcolor = :gray, ygridcolor = :gray)
    lines!(t_1, u_1[2, :], label = label_1)
    lines!(t_2, u_2[2, :], label = label_2)
    axislegend(ax2)
    # limits!(ax2, minimum(t), maximum(t), minimum(u[2, :]) - 1e-7*abs(minimum(u[2, :])), maximum(u[2, :]) + 1e-7*abs(maximum(u[2, :])))
    save(joinpath(path, "figures\\" * f_name * "\\psi_2.png"), fig2)

    fig3 = Figure(size = (1920, 720))
    ax3 = Axis(fig3[1, 1], title="x_D vs time", xlabel="t", ylabel="x_D"; xgridvisible = true, ygridvisible = true, xgridcolor = :gray, ygridcolor = :gray)
    lines!(t_1, u_1[3, :], label = label_1)
    lines!(t_2, u_2[3, :], label = label_2)
    axislegend(ax3)
    # limits!(ax3, minimum(t), maximum(t), minimum(u[3, :]) - 1e-7*abs(minimum(u[3, :])), maximum(u[3, :]) + 1e-7*abs(maximum(u[3, :])))
    save(joinpath(path, "figures\\" * f_name * "\\x_D.png"), fig3)

    fig4 = Figure(size = (1920, 720))
    ax4 = Axis(fig4[1, 1], title="y_D vs time", xlabel="t", ylabel="y_D"; xgridvisible = true, ygridvisible = true, xgridcolor = :gray, ygridcolor = :gray)
    lines!(t_1, u_1[4, :], label = label_1)
    lines!(t_2, u_2[4, :], label = label_2)
    axislegend(ax4)
    # limits!(ax4, minimum(t), maximum(t), minimum(u[4, :]) - 1e-7*abs(minimum(u[4, :])), maximum(u[4, :]) + 1e-7*abs(maximum(u[4, :])))
    save(joinpath(path, "figures\\" * f_name * "\\y_D.png"), fig4)
else 
    println("[info] DOF not plotted.")
end

print("Plot DOF difference? [Y/n]: ")
response = readline() |> strip
if response == "Y" || response == "y"
    fig1_diff = Figure(size = (1920, 720))
    ax1_diff = Axis(fig1_diff[1, 1], title="ψ₁ ($label_1 - $label_2) vs time", xlabel="t", ylabel="ψ₁"; xgridvisible = true, ygridvisible = true, xgridcolor = :gray, ygridcolor = :gray)
    lines!(t_1, u_1[1, :] .- u_2[1, :])
    # limits!(ax1, minimum(t), maximum(t), minimum(u[1, :]) - 1e-7*abs(minimum(u[1, :])), maximum(u[1, :]) + 1e-7*abs(maximum(u[1, :])))
    save(joinpath(path, "figures\\" * f_name * "\\psi_1_diff.png"), fig1_diff)

    fig2_diff = Figure(size = (1920, 720))
    ax2_diff = Axis(fig2_diff[1, 1], title="ψ₂ ($label_1 - $label_2) vs time", xlabel="t", ylabel="ψ₂"; xgridvisible = true, ygridvisible = true, xgridcolor = :gray, ygridcolor = :gray)
    lines!(t_1, u_1[2, :] .- u_2[2, :])
    # limits!(ax2, minimum(t), maximum(t), minimum(u[2, :]) - 1e-7*abs(minimum(u[2, :])), maximum(u[2, :]) + 1e-7*abs(maximum(u[2, :])))
    save(joinpath(path, "figures\\" * f_name * "\\psi_2_diff.png"), fig2_diff)

    fig3_diff = Figure(size = (1920, 720))
    ax3_diff = Axis(fig3_diff[1, 1], title="x_D ($label_1 - $label_2) vs time", xlabel="t", ylabel="x_D"; xgridvisible = true, ygridvisible = true, xgridcolor = :gray, ygridcolor = :gray)
    lines!(t_1, u_1[3, :] .- u_2[3, :])
    # limits!(ax3, minimum(t), maximum(t), minimum(u[3, :]) - 1e-7*abs(minimum(u[3, :])), maximum(u[3, :]) + 1e-7*abs(maximum(u[3, :])))
    save(joinpath(path, "figures\\" * f_name * "\\x_D_diff.png"), fig3_diff)

    fig4_diff = Figure(size = (1920, 720))
    ax4_diff = Axis(fig4_diff[1, 1], title="y_D ($label_1 - $label_2) vs time", xlabel="t", ylabel="y_D"; xgridvisible = true, ygridvisible = true, xgridcolor = :gray, ygridcolor = :gray)
    lines!(t_1, u_1[4, :] .- u_2[4, :])
    # limits!(ax4, minimum(t), maximum(t), minimum(u[4, :]) - 1e-7*abs(minimum(u[4, :])), maximum(u[4, :]) + 1e-7*abs(maximum(u[4, :])))
    save(joinpath(path, "figures\\" * f_name * "\\y_D_diff.png"), fig4_diff)
else 
    println("[info] DOF difference not plotted.")
end

print("Plot zoomed DOF end? [Y/n]: ")
response = readline() |> strip
if response == "Y" || response == "y"
    begin_index = round(Int, 0.95 * length(t))
    fig1_end = Figure(size = (1920, 720))
    ax1_end = Axis(fig1_end[1, 1], title="ψ₁ vs time", xlabel="t", ylabel="ψ₁"; xgridvisible = true, ygridvisible = true, xgridcolor = :gray, ygridcolor = :gray)
    lines!(t_1[begin_index:end], u_1[1, begin_index:end], label = label_1)
    lines!(t_2[begin_index:end], u_2[1, begin_index:end], label = label_2)
    axislegend(ax1_end)
    save(joinpath(path, "figures\\" * f_name * "\\psi_1_end.png"), fig1_end)

    fig4_end = Figure(size = (1920, 720))
    ax4_end = Axis(fig4_end[1, 1], title="y_D vs time", xlabel="t", ylabel="y_D"; xgridvisible = true, ygridvisible = true, xgridcolor = :gray, ygridcolor = :gray)
    lines!(t_1[begin_index:end], u_1[4, begin_index:end], label = label_1)
    lines!(t_2[begin_index:end], u_2[4, begin_index:end], label = label_2)
    axislegend(ax4_end)
    save(joinpath(path, "figures\\" * f_name * "\\y_D_end.png"), fig4_end)
else 
    println("[info] zoomed DOF end not plotted.")
end

print("Plot zoomed DOF difference end? [Y/n]: ")
response = readline() |> strip
if response == "Y" || response == "y"
    begin_index = round(Int, 0.95 * length(t))
    fig1_diff_end = Figure(size = (1920, 720))
    ax1_diff_end = Axis(fig1_diff_end[1, 1], title="ψ₁ ($label_1 - $label_2) vs time", xlabel="t", ylabel="ψ₁"; xgridvisible = true, ygridvisible = true, xgridcolor = :gray, ygridcolor = :gray)
    lines!(t_1[begin_index:end], u_1[1, begin_index:end] .- u_2[1, begin_index:end])
    save(joinpath(path, "figures\\" * f_name * "\\psi_1_diff_end.png"), fig1_end)

    fig4_diff_end = Figure(size = (1920, 720))
    ax4_diff_end = Axis(fig4_diff_end[1, 1], title="y_D ($label_1 - $label_2) vs time", xlabel="t", ylabel="y_D"; xgridvisible = true, ygridvisible = true, xgridcolor = :gray, ygridcolor = :gray)
    lines!(t_1[begin_index:end], u_1[4, begin_index:end] .- u_2[4, begin_index:end])
    save(joinpath(path, "figures\\" * f_name * "\\y_D_diff_end.png"), fig4_diff_end)
else 
    println("[info] zoomed DOF difference end not plotted.")
end