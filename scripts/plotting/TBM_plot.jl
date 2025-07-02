using JLD2, GLMakie
include("..\\..\\src\\Common.jl")

path = abspath(@__DIR__, "..\\..")
f_name = "TBM_tf200_dt1e-4_rk4_coulomb_test_matlab2"
@load joinpath(path, "data\\" * f_name * ".jld2") u t f pm

mkpath(joinpath(path, "figures/"*f_name))

## DOF PLOTS
print("Plot DOF? [Y/n]: ")
response = readline() |> strip
if response == "Y" || response == "y"
    fig1 = Figure(size = (1920, 720))
    ax1 = Axis(fig1[1, 1], title="ψ₁ vs time", xlabel="t", ylabel="ψ₁"; xgridvisible = true, ygridvisible = true, xgridcolor = :gray, ygridcolor = :gray)
    lines!(t, u[1, :])
    save(joinpath(path, "figures\\" * f_name * "\\psi_1.png"), fig1)

    fig2 = Figure(size = (1920, 720))
    ax2 = Axis(fig2[1, 1], title="ψ₂ vs time", xlabel="t", ylabel="ψ₂"; xgridvisible = true, ygridvisible = true, xgridcolor = :gray, ygridcolor = :gray)
    lines!(t, u[2, :])
    save(joinpath(path, "figures\\" * f_name * "\\psi_2.png"), fig2)

    fig3 = Figure(size = (1920, 720))
    ax3 = Axis(fig3[1, 1], title="x_D vs time", xlabel="t", ylabel="x_D"; xgridvisible = true, ygridvisible = true, xgridcolor = :gray, ygridcolor = :gray)
    lines!(t, u[3, :])
    save(joinpath(path, "figures\\" * f_name * "\\x_D.png"), fig3)

    fig4 = Figure(size = (1920, 720))
    ax4 = Axis(fig4[1, 1], title="y_D vs time", xlabel="t", ylabel="y_D"; xgridvisible = true, ygridvisible = true, xgridcolor = :gray, ygridcolor = :gray)
    lines!(t, u[4, :])
    save(joinpath(path, "figures\\" * f_name * "\\y_D.png"), fig4)

    fig5 = Figure(size = (1920, 720))
    ax5 = Axis(fig5[1, 1], title="ψ_D vs time", xlabel="t", ylabel="ψ_D"; xgridvisible = true, ygridvisible = true, xgridcolor = :gray, ygridcolor = :gray)
    lines!(t, u[5, :])
    save(joinpath(path, "figures\\" * f_name * "\\psi_D.png"), fig5)
else 
    println("[info] DOF not plotted.")
end

print("Plot zoomed DOF end? [Y/n]: ")
response = readline() |> strip
if response == "Y" || response == "y"
    begin_index = round(Int, 0.95 * length(t))
    fig1_end = Figure(size = (1920, 720))
    ax1_end = Axis(fig1_end[1, 1], title="ψ₁ vs time", xlabel="t", ylabel="ψ₁"; xgridvisible = true, ygridvisible = true, xgridcolor = :gray, ygridcolor = :gray)
    lines!(t[begin_index:end], u[1, begin_index:end])
    save(joinpath(path, "figures\\" * f_name * "\\psi_1_end.png"), fig1_end)

    fig4_end = Figure(size = (1920, 720))
    ax4_end = Axis(fig4_end[1, 1], title="y_D vs time", xlabel="t", ylabel="y_D"; xgridvisible = true, ygridvisible = true, xgridcolor = :gray, ygridcolor = :gray)
    lines!(t[begin_index:end], u[4, begin_index:end])
    save(joinpath(path, "figures\\" * f_name * "\\y_D_end.png"), fig4_end)
else 
    println("[info] zoomed DOF end not plotted.")
end

## FORCE PLOT
μ = 0.5
print("Plot force? [Y/n]: ")
response = readline() |> strip
if response == "Y" || response == "y"
    fig_T1 = Figure(size = (1920, 720))
    ax_T1 = Axis(fig_T1[1, 1], title="T₁ vs time", xlabel="t", ylabel="T₁"; xgridvisible = true, ygridvisible = true, xgridcolor = :gray, ygridcolor = :gray)
    lines!(t[1:end-1], pm.T₁[1:end-1])
    lines!(t[1:end-1], μ*pm.N₁[1:end-1], color = :red, linestyle = :dash) 
    lines!(t[1:end-1], -μ*pm.N₁[1:end-1], color = :red, linestyle = :dash) 
    save(joinpath(path, "figures\\" * f_name * "\\T1.png"), fig_T1)

    fig_N1 = Figure(size = (1920, 720))
    ax_N1 = Axis(fig_N1[1, 1], title="N₁ vs time", xlabel="t", ylabel="N₁"; xgridvisible = true, ygridvisible = true, xgridcolor = :gray, ygridcolor = :gray)
    lines!(t[1:end-1], pm.N₁[1:end-1])
    save(joinpath(path, "figures\\" * f_name * "\\N1.png"), fig_N1)

    fig_T2 = Figure(size = (1920, 720))
    ax_T2 = Axis(fig_T2[1, 1], title="T₂ vs time", xlabel="t", ylabel="T₂"; xgridvisible = true, ygridvisible = true, xgridcolor = :gray, ygridcolor = :gray)
    lines!(t[1:end-1], pm.T₂[1:end-1])
    lines!(t[1:end-1], μ*pm.N₂[1:end-1], color = :red, linestyle = :dash) 
    lines!(t[1:end-1], -μ*pm.N₂[1:end-1], color = :red, linestyle = :dash) 
    save(joinpath(path, "figures\\" * f_name * "\\T2.png"), fig_T2)

    fig_N2 = Figure(size = (1920, 720))
    ax_N2 = Axis(fig_N2[1, 1], title="N₂ vs time", xlabel="t", ylabel="N₂"; xgridvisible = true, ygridvisible = true, xgridcolor = :gray, ygridcolor = :gray)
    lines!(t[1:end-1], pm.N₂[1:end-1])
    save(joinpath(path, "figures\\" * f_name * "\\N2.png"), fig_N2)

else
    println("[info] Force not plotted.")
end

## CONTACT POSITION PLOT
print("Plot contact position? [Y/n]: ")
response = readline() |> strip
if response == "Y" || response == "y"
    print("Generate gifs for contact position? [Y/n]: ")
    response = readline() |> strip

    fig_contact_pos1 = Figure(size = (1080, 1080))
    ax_contact1 = Axis(fig_contact_pos1[1, 1], title="Position of the left contact", xlabel="x₁", ylabel="y₁"; xgridvisible = true, ygridvisible = true, xgridcolor = :gray, ygridcolor = :gray)
    lines!(pm.xₜ₁[1:end-1].*cos(pm.θ) .+ pm.xₙ₁[1:end-1].*sin(pm.θ), pm.xₜ₁[1:end-1].*sin(pm.θ) .- pm.xₙ₁[1:end-1].*cos(pm.θ))
    save(joinpath(path, "figures\\" * f_name * "\\left_contact.png"), fig_contact_pos1)
    if response == "Y" || response == "y"
        # Initialize scatter plot with a single placeholder point
        scatter_pos1 = scatter!(ax_contact1, [Point2f(0, 0)], color=:red, markersize=15)
        # Generate animation
        record(fig_contact_pos1, joinpath(path, "figures\\" * f_name * "\\left_contact.gif"), 1:500:length(t)-1) do i
            # Compute new coordinates
            new_x = pm.xₜ₁[i] * cos(pm.θ) + pm.xₙ₁[i] * sin(pm.θ)
            new_y = pm.xₜ₁[i] * sin(pm.θ) - pm.xₙ₁[i] * cos(pm.θ)
            # Directly modify scatter plot data
            scatter_pos1[1][] = [Point2f(new_x, new_y)]
        end
    end

    fig_contact_pos2 = Figure(size = (1080, 1080))
    ax_contact2 = Axis(fig_contact_pos2[1, 1], title="Position of the right contact", xlabel="x₂", ylabel="y₂"; xgridvisible = true, ygridvisible = true, xgridcolor = :gray, ygridcolor = :gray)
    lines!(-pm.xₜ₂[1:end-1].*cos(pm.θ) .- pm.xₙ₂[1:end-1].*sin(pm.θ), pm.xₜ₂[1:end-1].*sin(pm.θ) .- pm.xₙ₂[1:end-1].*cos(pm.θ))
    save(joinpath(path, "figures\\" * f_name * "\\right_contact.png"), fig_contact_pos2)
    if response == "Y" || response == "y"
        # Initialize scatter plot with a single placeholder point
        scatter_pos2 = scatter!(ax_contact2, [Point2f(0, 0)], color=:red, markersize=15)
        # Generate animation
        record(fig_contact_pos2, joinpath(path, "figures\\" * f_name * "\\right_contact.gif"), 1:500:length(t)-1) do i
            # Compute new coordinates
            new_x = -pm.xₜ₂[i] * cos(pm.θ) - pm.xₙ₂[i] * sin(pm.θ)
            new_y = pm.xₜ₂[i] * sin(pm.θ) - pm.xₙ₂[i] * cos(pm.θ)
            # Directly modify scatter plot data
            scatter_pos2[1][] = [Point2f(new_x, new_y)]
        end
    end

    # Check contact points position
    fig_contact_points1 = Figure(size = (1080, 1080))
    sin_θ = sin(pm.θ)
    cos_θ = cos(pm.θ)
    sin_ψ₁ = sin.(u[1, :])
    cos_ψ₁ = cos.(u[1, :])
    sin_ψ₂ = sin.(u[2, :])
    cos_ψ₂ = cos.(u[2, :])
    sin_θ_p_ψ = sin.(pm.θ .+ u[5, :])
    cos_θ_p_ψ = cos.(pm.θ .+ u[5, :])
    sin_θ_m_ψ = sin.(pm.θ .- u[5, :])
    cos_θ_m_ψ = cos.(pm.θ .- u[5, :])

    # Position of the contact points
    x_B = sin_ψ₁ .- pm.b*(1 .- cos_ψ₁)
    y_B = .-(1 .- cos_ψ₁) .- pm.b*sin_ψ₁
    x_E = u[3, :] .- pm.dₜ*(cos_θ_m_ψ .- cos_θ) .- pm.dₙ*(sin_θ_m_ψ .- sin_θ)
    y_E = u[4, :] .+ pm.dₙ*(cos_θ_m_ψ .- cos_θ) .- pm.dₜ*(sin_θ_m_ψ .- sin_θ)
    x_A = sin_ψ₂ .+ pm.b*(1 .- cos_ψ₂)
    y_A = .-(1 .-cos_ψ₂) .+ pm.b*sin_ψ₂
    x_F = u[3, :] .+ pm.dₜ*(cos_θ_p_ψ .- cos_θ) .+ pm.dₙ*(sin_θ_p_ψ .- sin_θ)
    y_F = u[4, :] .+ pm.dₙ*(cos_θ_p_ψ .- cos_θ) .- pm.dₜ*(sin_θ_p_ψ .- sin_θ)

    # Relative position of the contact
    xₜ₁ = (x_B .- x_E).*cos_θ .+ (y_B .- y_E).*sin_θ
    xₙ₁ = (x_B .- x_E).*sin_θ .- (y_B .- y_E).*cos_θ
    xₜ₂ = .-(x_A .- x_F).*cos_θ .+ (y_A .- y_F).*sin_θ
    xₙ₂ = .-(x_A .- x_F).*sin_θ .- (y_A .- y_F).*cos_θ

    fig_contact_points1 = Figure(size = (1080, 1080))
    ax_contact_points1 = Axis(fig_contact_points1[1, 1], title="Points of the left contact", xlabel="x", ylabel="y"; xgridvisible = true, ygridvisible = true, xgridcolor = :gray, ygridcolor = :gray)
    lines!(ax_contact_points1, x_B, y_B, label = "B")
    lines!(ax_contact_points1, x_E, y_E, label = "E")
    axislegend(ax_contact_points1)
    save(joinpath(path, "figures\\" * f_name * "\\left_contact_points.png"), fig_contact_points1)

    fig_contact_points2 = Figure(size = (1080, 1080))
    ax_contact_points2 = Axis(fig_contact_points2[1, 1], title="Points of the right contact", xlabel="x", ylabel="y"; xgridvisible = true, ygridvisible = true, xgridcolor = :gray, ygridcolor = :gray)
    lines!(ax_contact_points2, x_F, y_F, label = "F")
    lines!(ax_contact_points2, x_A, y_A, label = "A")
    axislegend(ax_contact_points2)
    save(joinpath(path, "figures\\" * f_name * "\\right_contact_points.png"), fig_contact_points2)

    fig_check_contact = Figure(size = (1080, 1080))
    ax_check_contact1 = Axis(fig_check_contact[1, 1], title="Tangential left")
    lines!(ax_check_contact1, 1:length(xₜ₁)-1, xₜ₁[1:end-1] .- pm.xₜ₁[1:end-1])
    ax_check_contact2 = Axis(fig_check_contact[1, 2], title="Normal left")
    lines!(ax_check_contact2, 1:length(xₙ₁)-1, xₙ₁[1:end-1] .- pm.xₙ₁[1:end-1] )
    ax_check_contact3 = Axis(fig_check_contact[2, 1], title="Tangential right")
    lines!(ax_check_contact3, 1:length(xₜ₂)-1, xₜ₂[1:end-1] .- pm.xₜ₂[1:end-1] )
    ax_check_contact4 = Axis(fig_check_contact[2, 2], title="Normal right")
    lines!(ax_check_contact4, 1:length(xₙ₂)-1, xₙ₂[1:end-1] .- pm.xₙ₂[1:end-1] )
    display(GLMakie.Screen(), fig_check_contact)
    
else
    println("[info] Contact position not plotted.")
end
 
## FRICTION CYCLES PLOT
print("Plot friction cycles? [Y/n]: ")
response = readline() |> strip
if response == "Y" || response == "y"
    print("Generate gifs for friction cycles? [Y/n]: ")
    response = readline() |> strip
    fig_cycle_left = Figure(size = (1080, 1080))
    ax_cycle1 = Axis(fig_cycle_left[1, 1], title="Left contact", xlabel="xₜ₁", ylabel="T₁"; xgridvisible = true, ygridvisible = true, xgridcolor = :gray, ygridcolor = :gray)
    lines!(pm.xₜ₁[1:end-1], pm.T₁[1:end-1])
    save(joinpath(path, "figures\\" * f_name * "\\cycle_left.png"), fig_cycle_left)
    if response == "Y" || response == "y"
        # Initialize scatter plot with a single placeholder point
        scatter_pos1 = scatter!(ax_cycle1, [Point2f(0, 0)], color=:red, markersize=15)
        # Generate animation
        record(fig_cycle_left, joinpath(path, "figures\\" * f_name * "\\cycle_left.gif"), 1:500:length(t)-1) do i
            # Compute new coordinates
            new_x = pm.xₜ₁[i]
            new_y = pm.T₁[i]
            # Directly modify scatter plot data
            scatter_pos1[1][] = [Point2f(new_x, new_y)]
        end
    end
    
    fig_cycle_right = Figure(size = (1080, 1080))
    ax_cycle2 = Axis(fig_cycle_right[1, 1], title="Right contact", xlabel="xₜ₂", ylabel="T₂"; xgridvisible = true, ygridvisible = true, xgridcolor = :gray, ygridcolor = :gray)
    lines!(pm.xₜ₂[1:end-1], pm.T₂[1:end-1])
    save(joinpath(path, "figures\\" * f_name * "\\cycle_right.png"), fig_cycle_right)
    if response == "Y" || response == "y"
        # Initialize scatter plot with a single placeholder point
        scatter_pos1 = scatter!(ax_cycle2, [Point2f(0, 0)], color=:red, markersize=15)
        # Generate animation
        record(fig_cycle_right, joinpath(path, "figures\\" * f_name * "\\cycle_right.gif"), 1:500:length(t)-1) do i
            # Compute new coordinates
            new_x = pm.xₜ₂[i]
            new_y = pm.T₂[i]
            # Directly modify scatter plot data
            scatter_pos1[1][] = [Point2f(new_x, new_y)]
        end
    end
else
    println("[info] Friction cycles not plotted.")
end