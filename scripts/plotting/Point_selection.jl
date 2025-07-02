using LinearAlgebra, JLD2, GLMakie, Statistics, Colors

path = abspath(@__DIR__, "..\\..")
## INITIALIZATION
begin
    # Get data from static equilibrium
    freetext = "_final"
    static_data = jldopen(joinpath(path, "\\data\\static\\Equilibrium_5DOF" * freetext *".jld2"), "r")
    sol_static = static_data["solutions"]
    pm_static = static_data["pm"]
    close(static_data)

    # Extract the data from static solution
    global T1_static = [sol.T1 for sol in sol_static]
    global T2_static = [sol.T2 for sol in sol_static]
    global wt1_static = [sol.wt1 for sol in sol_static]
    global wt2_static = [sol.wt2 for sol in sol_static]
    global n_static = length(T1_static)
    global u_static = zeros(5, n_static)
    u_static[1, :] = [sol.psi1 for sol in sol_static]
    u_static[2, :] = [sol.psi2 for sol in sol_static]
    u_static[3, :] = [sol.x_D for sol in sol_static]
    u_static[4, :] = [sol.y_D for sol in sol_static]
    u_static[5, :] = [sol.psi_D for sol in sol_static]
end

begin 
    # Plot the contourf
    fig_static_solution = Figure(size = (600, 600))
    ax_static_solution = Axis(fig_static_solution[1, 1], 
        title="Static equilibrium positions", 
        xlabel="T₁", 
        ylabel="T₂", 
        aspect = 1,
        xgridvisible = true, ygridvisible = true, 
        xgridcolor = :gray, ygridcolor = :gray
    )
    co_static_solution = contourf!(ax_static_solution,
        T1_static, 
        T2_static, 
        ones(n_static); # u_static[1, :]; 
        colormap = :grays
    )
    Colorbar(fig_static_solution[1, 2], 
        co_static_solution, 
        # colorrange = (minimum(u_static[1, :]) + 0.1 * (maximum(u_static[1, :]) - minimum(u_static[1, :])), maximum(u_static[1, :]))
    )
    Δ_T1 = (maximum(T1_static) - minimum(T1_static))/20
    Δ_T2 = (maximum(T2_static) - minimum(T2_static))/20
    limits!(ax_static_solution, 
        minimum(T1_static)-Δ_T1, maximum(T1_static)+Δ_T1, 
        minimum(T2_static)-Δ_T2, maximum(T2_static)+Δ_T2
    )

    # Add a scatter marker to show selected point
    pt = Observable(Point2f(0, 0))
    scatter!(ax_static_solution, pt, color=:red, markersize=10)

    function get_idx_at(p; tol=0.51)
        # Find nearest grid indices
        # Tolerance is relative to the distance between T values, value of 0.51 should always ensure that the point is found
        ΔT₁ =  abs(T1_static[findfirst(x -> abs(x - T1_static[1]) > 1e-8, T1_static[2:end]) + 1] - T1_static[1])
        println("ΔT₁ = $ΔT₁")
        for idx in eachindex(T1_static)
             if abs(T1_static[idx] - p[1]) < tol*ΔT₁ && abs(T2_static[idx] - p[2]) < tol*ΔT₁
                return idx
            end
        end
        println("Point ($(p[1]), $(p[2])) not found in the surface.")
        return nothing
    end

    on(events(ax_static_solution).mousebutton) do event
        if (event.button == Mouse.middle || event.button == Mouse.button_5) && event.action == Mouse.press
            p = mouseposition(ax_static_solution.scene) 
            i = get_idx_at(p, tol=0.51)
            if i !== nothing
                pt[] = Point2f(T1_static[i], T2_static[i])  # Update the point position
                println("Clicked at i=$i, T₁=$(round(T1_static[i], sigdigits=3)), T₂=$(round(T2_static[i], sigdigits=3)), u=$(round(u_static[1, i], sigdigits=3))")
            end
        end
    end

    display(GLMakie.Screen(), fig_static_solution)
end