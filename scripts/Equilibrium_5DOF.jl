using GLMakie, JLD2, MAT
include("..\\src\\StaticModels.jl")

path = abspath(@__DIR__, "..")
## INITIALIZATION
begin
    pm = TBM_5DOF_Parameters(
        # Parameters given by Zara et al.
        I_B = 6.348e-4, # kg m² / rad
        k_ψ = 4.0e3, # N m / rad
        c_ψ = 5.0e-2, # N m s / rad (5.0e-3)
        I_D = 4.558e-8, # kg m²
        m_D = 1.0e-5, # kg (1e-2)
        k_Dx = 2e-1, # N / m (20)
        k_Dψ = 2e-3, # N m / rad (20) 
        θ = pi/4, # rad
        FC = 160.0, # N
        μ = 0.25, # 0.5
        kₙ = 2.0e5, # N / m (3.0e5)
        kₜ = 2.0e5, # N / m (3.0e5)
        # Parameters not given by Zara et al.
        f = 0.1, # m  
        b = 0.05, # m
        dₙ = 0.007, # m
        dₜ = 0.001 # m
    )

    reduced = true # Determine only half of the possible equilibrium positions
    verbose = false # Print information about the solver convergence
    interactive = false
    freetext = "_final_even_lower_def"

    if reduced
        freetext = freetext * "_reduced"
    end
end

# Static range of solutions
begin
    T_range = range(-0.0008, 0.001, length = 19) # T values
    # T_range = range(0.0, -0.0025, length = 251) # T values
    
    x_D_limit = [-5e-2, 5e-2]
    y_D_limit = [-5e-2, 5e-2]
    ψ_D_limit = [-pi/10, pi/10]

    # Variable where to store solutions
    solutions = Vector{NamedTuple{
        (:psi1, :psi2, :x_D, :y_D, :psi_D, :T1, :T2, :N1, :N2, :xt1, :xn1, :wt1, :xt2, :xn2, :wt2), 
        Tuple{Float64, Float64, Float64, Float64, Float64, Float64, Float64, Float64, Float64, Float64, Float64, Float64, Float64, Float64, Float64}
    }}()
    contact_limit = Vector{NamedTuple{
        (:T1, :T2), 
        Tuple{Float64, Float64}
    }}()

    start_time = time()
    # Loop over the range of T values
    local u₀ = zeros(5)
    local u₀₋₁ = zeros(5)
    for (idx_T1, T₁) in enumerate(T_range), (idx_T2, T₂) in enumerate(T_range)
        if reduced && T₁ > T₂
            continue # Skip the second half of the solutions if reduced is true
        end

        # Update the initial guess for the equilibrium position from different T₁ value
        if (idx_T1 > 1 && abs(T_range[idx_T1-1] - T₁) > 1e-12)
            u₀ = u₀₋₁
        end
        # Solve nonlinear system for equilibrium position
        eq_pos = equilibrium_pos_5DOF(pm, T₁, T₂, u₀)
        if eq_pos === nothing
            verbose && println("[warning] T₁ = $T₁, T₂ = $T₂ : No solution found")
            continue
        end
        # Save result for the equilibrium position from different T₁ value
        if (idx_T1 > 1 && abs(T_range[idx_T1-1] - T₁) > 1e-12)
            u₀₋₁ = eq_pos
        end
        u₀ = eq_pos

        verbose && println("[info] T₁ = $T₁, T₂ = $T₂ : Solution found")

        # Relative displacement of the contact points
        xₜ₁, xₙ₁, xₜ₂, xₙ₂ = contact_points_displacements(eq_pos, pm)

        # Normal forces
        kₙ = pm.friction_law_1.kₙ
        kₜ = pm.friction_law_1.kₜ
        μ = pm.friction_law_1.μ
        N₁ = kₙ*xₙ₁
        N₂ = kₙ*xₙ₂

        if (abs(T₁) < abs(μ*N₁) && abs(T₂) < abs(μ*N₂) && 
            N₁ > 0.0 && N₂ > 0.0) # Contact forces conditions are satisfied
            wₜ₁ = xₜ₁ - T₁/kₜ
            wₜ₂ = xₜ₂ - T₂/kₜ
            push!(solutions, (psi1 = eq_pos[1], psi2 = eq_pos[2], x_D = eq_pos[3], y_D = eq_pos[4], psi_D = eq_pos[5], T1 = T₁, T2 = T₂, N1 = N₁, N2 = N₂, xt1 = xₜ₁, xn1 = xₙ₁, wt1 = wₜ₁, xt2 = xₜ₂, xn2 = xₙ₂, wt2 = wₜ₂))
        else
            push!(contact_limit, (T1 = T₁, T2 = T₂))
        end 

    end
    elapsed_time = time() - start_time
    println("[info] Elapsed time: $elapsed_time seconds")

    # Extract the data for plotting
    psi1_values = [sol.psi1 for sol in solutions]
    psi2_values = [sol.psi2 for sol in solutions]
    x_D_values = [sol.x_D for sol in solutions]
    y_D_values = [sol.y_D for sol in solutions]
    psi_D_values = [sol.psi_D for sol in solutions]
    T1_values = [sol.T1 for sol in solutions]
    T2_values = [sol.T2 for sol in solutions]
    N1_values = [sol.N1 for sol in solutions]
    N2_values = [sol.N2 for sol in solutions]
    xt1_values = [sol.xt1 for sol in solutions]
    xn1_values = [sol.xn1 for sol in solutions]
    wt1_values = [sol.wt1 for sol in solutions]
    xt2_values = [sol.xt2 for sol in solutions]
    xn2_values = [sol.xn2 for sol in solutions]
    wt2_values = [sol.wt2 for sol in solutions]

end

# Data plotting and storage
begin
    # Plotting
    if interactive
        print("Generate figures? [Y/n]: ")
        response = readline() |> strip
    else
        response = "Y"
    end
    if response == "Y" || response == "y"

        mkpath(joinpath(path, "figures\\Equilibrium_5DOF" * freetext))

        Δ_T1 = (maximum(T1_values) - minimum(T1_values))/20
        Δ_T2 = (maximum(T2_values) - minimum(T2_values))/20

        fig_psi1 = Figure(size = (720, 720))
        ax_psi1 = Axis(fig_psi1[1, 1], 
            title="ψ₁ vs T₁ & T₂", 
            xlabel="T₁", 
            ylabel="T₂", 
            aspect = 1
        )
        limits!(ax_psi1, 
            minimum(T1_values)-Δ_T1, maximum(T1_values)+Δ_T1, 
            minimum(T2_values)-Δ_T2, maximum(T2_values)+Δ_T2
        )
        co = contourf!(ax_psi1, 
            T1_values, 
            T2_values, 
            psi1_values
        )
        Colorbar(fig_psi1[1, 2], co)
        save(joinpath(path, "figures\\Equilibrium_5DOF" * freetext * "\\psi_1.png"), fig_psi1)

        fig_psi2 = Figure(size = (720, 720))
        ax_psi2 = Axis(fig_psi2[1, 1], 
            title="ψ₂ vs T₁ & T₂", 
            xlabel="T₁", 
            ylabel="T₂", 
            aspect = 1
        )
        limits!(ax_psi2,
             minimum(T1_values)-Δ_T1, maximum(T1_values)+Δ_T1, 
             minimum(T2_values)-Δ_T2, maximum(T2_values)+Δ_T2
        )
        co = contourf!(ax_psi2, 
            T1_values,
            T2_values, 
            psi2_values
        )
        Colorbar(fig_psi2[1, 2], co)
        save(joinpath(path, "figures\\Equilibrium_5DOF" * freetext * "\\psi_2.png"), fig_psi2)

        fig_x_D = Figure(size = (720, 720))
        ax_x_D = Axis(fig_x_D[1, 1], 
            title="x_D vs T₁ & T₂", 
            xlabel="T₁", 
            ylabel="T₂", 
            aspect = 1
        )
        limits!(ax_x_D, 
            minimum(T1_values)-Δ_T1, maximum(T1_values)+Δ_T1, 
            minimum(T2_values)-Δ_T2, maximum(T2_values)+Δ_T2
        )
        co = contourf!(ax_x_D, 
            T1_values, 
            T2_values, 
            x_D_values
        )
        Colorbar(fig_x_D[1, 2], co)
        save(joinpath(path, "figures\\Equilibrium_5DOF" * freetext * "\\x_D.png"), fig_x_D)

        fig_y_D = Figure(size = (720, 720))
        ax_y_D = Axis(fig_y_D[1, 1], 
            title="y_D vs T₁ & T₂", 
            xlabel="T₁", 
            ylabel="T₂", 
            aspect = 1
        )
        limits!(ax_y_D, 
            minimum(T1_values)-Δ_T1, maximum(T1_values)+Δ_T1, 
            minimum(T2_values)-Δ_T2, maximum(T2_values)+Δ_T2
        )
        co = contourf!(ax_y_D, 
            T1_values, 
            T2_values, 
            y_D_values
        )
        Colorbar(fig_y_D[1, 2], co)
        save(joinpath(path, "figures\\Equilibrium_5DOF" * freetext * "\\y_D.png"), fig_y_D)

        fig_psi_D = Figure(size = (720, 720))
        ax_psi_D = Axis(fig_psi_D[1, 1], 
            title="ψ_D vs T₁ & T₂", 
            xlabel="T₁", 
            ylabel="T₂", 
            aspect = 1
        )
        limits!(ax_psi_D, 
            minimum(T1_values)-Δ_T1, maximum(T1_values)+Δ_T1, 
            minimum(T2_values)-Δ_T2, maximum(T2_values)+Δ_T2
        )
        co = contourf!(ax_psi_D, 
            T1_values, 
            T2_values, 
            psi_D_values
        )
        Colorbar(fig_psi_D[1, 2], co)
        save(joinpath(path, "figures\\Equilibrium_5DOF" * freetext * "\\psi_D.png"), fig_psi_D)

        fig_N1 = Figure(size = (720, 720))
        ax_N1 = Axis(fig_N1[1, 1], 
            title="N₁ vs T₁ & T₂", 
            xlabel="T₁", 
            ylabel="T₂", 
            aspect = 1
        )
        limits!(ax_N1, 
            minimum(T1_values)-Δ_T1, maximum(T1_values)+Δ_T1, 
            minimum(T2_values)-Δ_T2, maximum(T2_values)+Δ_T2
        )
        co = contourf!(ax_N1, 
            T1_values, T2_values, N1_values
        )
        Colorbar(fig_N1[1, 2], co)
        save(joinpath(path, "figures\\Equilibrium_5DOF" * freetext * "\\N_1.png"), fig_N1)

        fig_N2 = Figure(size = (720, 720))
        ax_N2 = Axis(fig_N2[1, 1], 
            title="N₂ vs T₁ & T₂", 
            xlabel="T₁", 
            ylabel="T₂", 
            aspect = 1
        )
        limits!(ax_N2, 
            minimum(T1_values)-Δ_T1, maximum(T1_values)+Δ_T1, 
            minimum(T2_values)-Δ_T2, maximum(T2_values)+Δ_T2
        )
        co = contourf!(ax_N2, 
            T1_values, T2_values, N2_values
        )
        Colorbar(fig_N2[1, 2], co)
        save(joinpath(path, "figures\\Equilibrium_5DOF" * freetext * "\\N_2.png"), fig_N2)

        fig_xt1 = Figure(size = (720, 720))
        ax_xt1 = Axis(fig_xt1[1, 1], 
            title="xₜ₁ vs T₁ & T₂", 
            xlabel="T₁", 
            ylabel="T₂", 
            aspect = 1
        )
        limits!(ax_xt1, 
            minimum(T1_values)-Δ_T1, maximum(T1_values)+Δ_T1, 
            minimum(T2_values)-Δ_T2, maximum(T2_values)+Δ_T2
        )
        co = contourf!(ax_xt1, 
            T1_values, 
            T2_values, 
            xt1_values
        )
        Colorbar(fig_xt1[1, 2], co)
        save(joinpath(path, "figures\\Equilibrium_5DOF" * freetext * "\\x_t_1.png"), fig_xt1)

        fig_xn1 = Figure(size = (720, 720))
        ax_xn1 = Axis(fig_xn1[1, 1], 
            title="xₙ₁ vs T₁ & T₂", 
            xlabel="T₁", 
            ylabel="T₂", 
            aspect = 1
        )
        limits!(ax_xn1, 
            minimum(T1_values)-Δ_T1, maximum(T1_values)+Δ_T1, 
            minimum(T2_values)-Δ_T2, maximum(T2_values)+Δ_T2
        )
        co = contourf!(ax_xn1, 
            T1_values, 
            T2_values, 
            xn1_values
        )
        Colorbar(fig_xn1[1, 2], co)
        save(joinpath(path, "figures\\Equilibrium_5DOF" * freetext * "\\x_n_1.png"), fig_xn1)

        fig_wt1 = Figure(size = (720, 720))
        ax_wt1 = Axis(fig_wt1[1, 1], 
            title="wₜ₁ vs T₁ & T₂", 
            xlabel="T₁", 
            ylabel="T₂", 
            aspect = 1
        )
        limits!(ax_wt1, 
            minimum(T1_values)-Δ_T1, maximum(T1_values)+Δ_T1, 
            minimum(T2_values)-Δ_T2, maximum(T2_values)+Δ_T2
        )
        co = contourf!(ax_wt1, 
            T1_values, 
            T2_values, 
            wt1_values
        )
        Colorbar(fig_wt1[1, 2], co)
        save(joinpath(path, "figures\\Equilibrium_5DOF" * freetext * "\\w_t_1.png"), fig_wt1)

        fig_xt2 = Figure(size = (720, 720))
        ax_xt2 = Axis(fig_xt2[1, 1], 
            title="xₜ₂ vs T₂ & T₂", 
            xlabel="T₂", 
            ylabel="T₂", 
            aspect = 1
        )
        limits!(ax_xt2, 
            minimum(T1_values)-Δ_T1, maximum(T1_values)+Δ_T1, 
            minimum(T2_values)-Δ_T2, maximum(T2_values)+Δ_T2
        )
        co = contourf!(ax_xt2, 
            T1_values, 
            T2_values, 
            xt2_values
        )
        Colorbar(fig_xt2[1, 2], co)
        save(joinpath(path, "figures\\Equilibrium_5DOF" * freetext * "\\x_t_2.png"), fig_xt2)

        fig_xn2 = Figure(size = (720, 720))
        ax_xn2 = Axis(fig_xn2[1, 1], 
            title="xₙ₂ vs T₂ & T₂", 
            xlabel="T₂", 
            ylabel="T₂", 
            aspect = 1
        )
        limits!(ax_xn2, 
            minimum(T1_values)-Δ_T1, maximum(T1_values)+Δ_T1, 
            minimum(T2_values)-Δ_T2, maximum(T2_values)+Δ_T2
        )
        co = contourf!(ax_xn2, 
            T1_values, 
            T2_values, 
            xn2_values
        )
        Colorbar(fig_xn2[1, 2], co)
        save(joinpath(path, "figures\\Equilibrium_5DOF" * freetext * "\\x_n_2.png"), fig_xn2)

        fig_wt2 = Figure(size = (720, 720))
        ax_wt2 = Axis(fig_wt2[1, 1], 
            title="wₜ₂ vs T₂ & T₂", 
            xlabel="T₂", 
            ylabel="T₂", 
            aspect = 1
        )
        limits!(ax_wt2, 
            minimum(T1_values)-Δ_T1, maximum(T1_values)+Δ_T1, 
            minimum(T2_values)-Δ_T2, maximum(T2_values)+Δ_T2
        )
        co = contourf!(ax_wt2, 
            T1_values, 
            T2_values, 
            wt2_values
        )
        Colorbar(fig_wt2[1, 2], co)
        save(joinpath(path, "figures\\Equilibrium_5DOF" * freetext * "\\w_t_2.png"), fig_wt2)

        #=
        fig_lim_contact = Figure(size = (720, 720))
        ax_lim_contact = Axis(fig_lim_contact[1, 1], 
            title="Contact limit", 
            xlabel="T₂", 
            ylabel="T₂", 
            aspect = 1
        )
        limits!(ax_lim_contact, 
            minimum(T1_values)-Δ_T1, maximum(T1_values)+Δ_T1, 
            minimum(T2_values)-Δ_T2, maximum(T2_values)+Δ_T2
        )
        co = contourf!(ax_lim_contact, 
            T1_cont, 
            T2_cont, 
            ones(length(T1_cont))
        )
        save(joinpath(path, "figures\\Equilibrium_5DOF" * freetext * "\\contact_limit.png"), fig_lim_contact)

        fig_lim_position = Figure(size = (720, 720))
        ax_lim_position = Axis(fig_lim_position[1, 1], 
            title="Position limit", 
            xlabel="T₂", 
            ylabel="T₂", 
            aspect = 1
        )
        limits!(ax_lim_position, 
            minimum(T1_values)-Δ_T1, maximum(T1_values)+Δ_T1, 
            minimum(T2_values)-Δ_T2, maximum(T2_values)+Δ_T2
        )
        co = contourf!(ax_lim_position, 
            T1_post, 
            T2_post, 
            ones(length(T1_post))
        )
        save(joinpath(path, "figures\\Equilibrium_5DOF" * freetext * "\\position_limit.png"), fig_lim_position)
        =#

    else
        println("[info] Figures not generated")
    end

    f_path = joinpath(path, "data\\static\\Equilibrium_5DOF" * freetext * ".jld2")
    if interactive
        print("Save results to file? [Y/n]: ")
        response = readline() |> strip
    else
        response = "Y"
    end
    if response == "Y" || response == "y"
        if isfile(f_path)
            println("[warning] File $f_path already exists. Choose a different name or delete the file")
        else
            mkpath(joinpath(path, "data\\static"))
            @save f_path solutions pm
            println("[info] Results saved to $f_path")
        end
    else
        println("[info] Results not saved")
    end

    if interactive
        print("Plot results? [Y/n]: ")
        response = readline() |> strip
    else
        response = "n"
    end
    if response == "Y" || response == "y"
        mat_file = replace(f_path, ".jld2" => ".mat")
        matwrite(mat_file, Dict(
                                "T1" => collect(T1_values),
                                "T2" => collect(T2_values),
                                "psi1" => collect(psi1_values),
                                "psi2" => collect(psi2_values),
                                "x_D" => collect(x_D_values),
                                "y_D" => collect(y_D_values),
                                "psi_D" => collect(psi_D_values),
                                "N_1" => collect(N1_values),
                                "N_2" => collect(N2_values),
                                "x_t1" => collect(xt1_values),
                                "x_n1" => collect(xn1_values),
                                "w_t1" => collect(wt1_values),
                                "x_t2" => collect(xt2_values),
                                "x_n2" => collect(xn2_values),
                                "w_t2" => collect(wt2_values),
                                "zeta" => pm.ζ_ψ,
                                "alpha" => pm.α,
                                "beta" => pm.β,
                                "f_x" => pm.fₓ,
                                "f_psi" => pm.f_ψ,
                                "theta" => pm.θ,
                                "b" => pm.b,
                                "d_n" => pm.dₙ,
                                "d_t" => pm.dₜ
        ))
        println("[info] Results also saved to MATLAB file in $mat_file")
    else
        println("[info] Results not saved to MATLAB")
    end
end
   

#=
## AUXILIARY ANALYSISd
# Linear approximation
function linear_aproximation(pm::TBM_5DOF_Parameters, T₁::Float64, T₂::Float64)
    @unpack θ, b, dₙ, dₜ, fₓ, f_ψ, α, β, FC, friction_law_1 = pm
    kₙ = friction_law_1.kₙ

    A1 = -(kₙ^2*b*α*FC*fₓ^2*(b^2*f_ψ^2 + 2*β*dₜ^2)*cos_θ^4 + 3*(b^2*f_ψ^2 + (2*dₜ^2*β)/3)*fₓ^2*α*FC*kₙ^2*sin_θ*cos_θ^3 - (-2*α*FC*kₙ*b*((dₜ^2*β + (3*f_ψ^2)/2)*fₓ^2 + α*f_ψ^2)*sin_θ^2 - b*f_ψ^2*α*FC*fₓ^2)*kₙ*cos_θ^2 + (2*α*FC*kₙ^2*((dₜ^2*β + f_ψ^2/2)*fₓ^2 + α*f_ψ^2)*sin_θ^3 + f_ψ^2*α*sin_θ*FC*fₓ^2*kₙ)*cos_θ)/((fₓ^2*kₙ*(b^2*f_ψ^2 + 2*β*dₜ^2)*cos_θ^2 + 2*sin_θ*cos_θ*b*f_ψ^2*fₓ^2*kₙ + 2*kₙ*((dₜ^2*β + f_ψ^2/2)*fₓ^2 + α*f_ψ^2)*sin_θ^2 + f_ψ^2*fₓ^2)*(kₙ*(b^2*fₓ^2 + 2*α)*cos_θ^2 + 2*sin_θ*cos_θ*b*fₓ^2*kₙ + (kₙ*sin_θ^2 + 1)*fₓ^2))
    B1 = -(fₓ^2*(β*(b*dₙ + dₜ)*b^2*dₜ*fₓ^2 + α*(b^2*f_ψ^2 + 2*b*β*dₙ*dₜ + 4*β*dₜ^2))*kₙ^2*cos_θ^5 - kₙ^2*((b^2*dₜ - 3*b*dₙ - 2*dₜ)*β*b*dₜ*fₓ^4 + (f_ψ^2*b^3 + 2*(β*dₜ^2 - f_ψ^2)*b - 2*β*dₙ*dₜ)*α*fₓ^2 + 2*b*f_ψ^2*α^2)*sin_θ*cos_θ^4 + (((b^3*dₙ - b^2*dₜ + 3*b*dₙ + dₜ)*β*dₜ*fₓ^4 - (b^2*f_ψ^2 - 2*b*β*dₙ*dₜ - 6*β*dₜ^2 - f_ψ^2)*α*fₓ^2 + 2*f_ψ^2*α^2)*kₙ*sin_θ^2 + fₓ^2*((b^2*f_ψ^2 + b*β*dₙ*dₜ + 2*β*dₜ^2)*fₓ^2 + 2*α*f_ψ^2))*kₙ*cos_θ^3 - (kₙ*(β*(b^3*dₜ - 3*b^2*dₙ - b*dₜ - dₙ)*dₜ*fₓ^4 + α*(f_ψ^2*b^3 + (2*β*dₜ^2 - f_ψ^2)*b - 2*β*dₙ*dₜ)*fₓ^2 + 2*b*f_ψ^2*α^2)*sin_θ^3 + fₓ^2*((f_ψ^2*b^3 - 2*(-β*dₜ^2 + f_ψ^2)*b - β*dₙ*dₜ)*fₓ^2 + 2*b*f_ψ^2*α)*sin_θ)*kₙ*cos_θ^2 + (-2*(β*(b^2*dₜ - 3/2*dₙ*b - 1/2*dₜ)*dₜ*fₓ^4 + (b^2*f_ψ^2 - dₜ^2*β - 1/2*f_ψ^2)*α*fₓ^2 - f_ψ^2*α^2)*kₙ^2*sin_θ^4 - 2*fₓ^2*((b^2*f_ψ^2 - 1/2*dₙ*dₜ*β*b - dₜ^2*β - 1/2*f_ψ^2)*fₓ^2 - α*f_ψ^2)*kₙ*sin_θ^2 + f_ψ^2*fₓ^4)*cos_θ - fₓ^2*sin_θ*(((b*dₜ - dₙ)*β*dₜ*fₓ^2 + b*f_ψ^2*α)*kₙ^2*sin_θ^4 + 2*((((2*β*dₜ^2 + f_ψ^2)*b - β*dₙ*dₜ)*fₓ^2)/2 + b*f_ψ^2*α)*kₙ*sin_θ^2 + b*f_ψ^2*fₓ^2))/((fₓ^2*kₙ*(b^2*f_ψ^2 + 2*β*dₜ^2)*cos_θ^2 + 2*sin_θ*cos_θ*b*f_ψ^2*fₓ^2*kₙ + 2*kₙ*((dₜ^2*β + f_ψ^2/2)*fₓ^2 + α*f_ψ^2)*sin_θ^2 + f_ψ^2*fₓ^2)*(kₙ*(b^2*fₓ^2 + 2*α)*cos_θ^2 + 2*sin_θ*cos_θ*b*fₓ^2*kₙ + (kₙ*sin_θ^2 + 1)*fₓ^2))
    C1 = -(fₓ^2*(-β*(b*dₙ + dₜ)*b^2*dₜ*fₓ^2 + α*(b^2*f_ψ^2 - 2*b*β*dₙ*dₜ))*kₙ^2*cos_θ^5 - (-(b^2*dₜ - 3*b*dₙ - 2*dₜ)*β*b*dₜ*fₓ^4 + (-f_ψ^2*b^3 + 2*(-β*dₜ^2 - f_ψ^2)*b + 2*β*dₙ*dₜ)*α*fₓ^2 - 2*b*f_ψ^2*α^2)*sin_θ*kₙ^2*cos_θ^4 + ((-(b^3*dₙ - b^2*dₜ + 3*b*dₙ + dₜ)*β*dₜ*fₓ^4 - α*(-3*b^2*f_ψ^2 + 2*b*β*dₙ*dₜ - 2*β*dₜ^2 - f_ψ^2)*fₓ^2 + 2*f_ψ^2*α^2)*kₙ*sin_θ^2 - β*dₜ*fₓ^4*dₙ*b)*kₙ*cos_θ^3 - kₙ*((-β*(b^3*dₜ - 3*b^2*dₙ - b*dₜ - dₙ)*dₜ*fₓ^4 + α*(-f_ψ^2*b^3 + (-2*β*dₜ^2 - 3*f_ψ^2)*b + 2*β*dₙ*dₜ)*fₓ^2 - 2*b*f_ψ^2*α^2)*kₙ*sin_θ^3 + fₓ^2*(β*dₙ*dₜ*fₓ^2 - 2*α*b*f_ψ^2)*sin_θ)*cos_θ^2 + (-2*(-β*(b^2*dₜ - 3/2*dₙ*b - 1/2*dₜ)*dₜ*fₓ^4 + α*(-b^2*f_ψ^2 - dₜ^2*β - 1/2*f_ψ^2)*fₓ^2 - f_ψ^2*α^2)*kₙ^2*sin_θ^4 - 2*(dₙ*dₜ*β*b*fₓ^2/2 - α*f_ψ^2)*fₓ^2*kₙ*sin_θ^2)*cos_θ - sin_θ*fₓ^2*(-((b*dₜ - dₙ)*β*dₜ*fₓ^2 + b*f_ψ^2*α)*kₙ^2*sin_θ^4 + kₙ*sin_θ^2*β*dₙ*dₜ*fₓ^2))/((fₓ^2*kₙ*(b^2*f_ψ^2 + 2*β*dₜ^2)*cos_θ^2 + 2*sin_θ*cos_θ*b*f_ψ^2*fₓ^2*kₙ + 2*kₙ*((dₜ^2*β + f_ψ^2/2)*fₓ^2 + α*f_ψ^2)*sin_θ^2 + f_ψ^2*fₓ^2)*(kₙ*(b^2*fₓ^2 + 2*α)*cos_θ^2 + 2*sin_θ*cos_θ*b*fₓ^2*kₙ + (kₙ*sin_θ^2 + 1)*fₓ^2))
    psi1 = A1 .+ B1*T₁ .+ C1*T₂

    A2 = -(-kₙ^2*b*α*FC*fₓ^2*(b^2*f_ψ^2 + 2*β*dₜ^2)*cos_θ^4 - 3*(b^2*f_ψ^2 + (2*dₜ^2*β)/3)*fₓ^2*α*FC*kₙ^2*sin_θ*cos_θ^3 - (2*α*FC*kₙ*b*((dₜ^2*β + (3*f_ψ^2)/2)*fₓ^2 + α*f_ψ^2)*sin_θ^2 + b*f_ψ^2*α*FC*fₓ^2)*kₙ*cos_θ^2 + (-2*α*FC*kₙ^2*((dₜ^2*β + f_ψ^2/2)*fₓ^2 + α*f_ψ^2)*sin_θ^3 - f_ψ^2*α*sin_θ*FC*fₓ^2*kₙ)*cos_θ)/((fₓ^2*kₙ*(b^2*f_ψ^2 + 2*β*dₜ^2)*cos_θ^2 + 2*sin_θ*cos_θ*b*f_ψ^2*fₓ^2*kₙ + 2*kₙ*((dₜ^2*β + f_ψ^2/2)*fₓ^2 + α*f_ψ^2)*sin_θ^2 + f_ψ^2*fₓ^2)*(kₙ*(b^2*fₓ^2 + 2*α)*cos_θ^2 + 2*sin_θ*cos_θ*b*fₓ^2*kₙ + (kₙ*sin_θ^2 + 1)*fₓ^2))
    B2 = -(fₓ^2*(β*(b*dₙ + dₜ)*b^2*dₜ*fₓ^2 - α*(b^2*f_ψ^2 - 2*b*β*dₙ*dₜ))*kₙ^2*cos_θ^5 - ((b^2*dₜ - 3*b*dₙ - 2*dₜ)*β*b*dₜ*fₓ^4 + (f_ψ^2*b^3 + 2*(β*dₜ^2 + f_ψ^2)*b - 2*β*dₙ*dₜ)*α*fₓ^2 + 2*b*f_ψ^2*α^2)*sin_θ*kₙ^2*cos_θ^4 + (kₙ*((b^3*dₙ - b^2*dₜ + 3*b*dₙ + dₜ)*β*dₜ*fₓ^4 - 3*(b^2*f_ψ^2 - 2/3*dₙ*dₜ*β*b + 2/3*dₜ^2*β + 1/3*f_ψ^2)*α*fₓ^2 - 2*f_ψ^2*α^2)*sin_θ^2 + β*dₜ*fₓ^4*dₙ*b)*kₙ*cos_θ^3 - ((β*(b^3*dₜ - 3*b^2*dₙ - b*dₜ - dₙ)*dₜ*fₓ^4 + (f_ψ^2*b^3 + (2*β*dₜ^2 + 3*f_ψ^2)*b - 2*β*dₙ*dₜ)*α*fₓ^2 + 2*b*f_ψ^2*α^2)*kₙ*sin_θ^3 - fₓ^2*(β*dₙ*dₜ*fₓ^2 - 2*α*b*f_ψ^2)*sin_θ)*kₙ*cos_θ^2 + (-2*(β*(b^2*dₜ - 3/2*dₙ*b - 1/2*dₜ)*dₜ*fₓ^4 + (b^2*f_ψ^2 + dₜ^2*β + 1/2*f_ψ^2)*α*fₓ^2 + f_ψ^2*α^2)*kₙ^2*sin_θ^4 + 2*(dₙ*dₜ*β*b*fₓ^2/2 - α*f_ψ^2)*fₓ^2*kₙ*sin_θ^2)*cos_θ - fₓ^2*sin_θ*(((b*dₜ - dₙ)*β*dₜ*fₓ^2 + b*f_ψ^2*α)*kₙ^2*sin_θ^4 - kₙ*sin_θ^2*β*dₙ*dₜ*fₓ^2))/((fₓ^2*kₙ*(b^2*f_ψ^2 + 2*β*dₜ^2)*cos_θ^2 + 2*sin_θ*cos_θ*b*f_ψ^2*fₓ^2*kₙ + 2*kₙ*((dₜ^2*β + f_ψ^2/2)*fₓ^2 + α*f_ψ^2)*sin_θ^2 + f_ψ^2*fₓ^2)*(kₙ*(b^2*fₓ^2 + 2*α)*cos_θ^2 + 2*sin_θ*cos_θ*b*fₓ^2*kₙ + (kₙ*sin_θ^2 + 1)*fₓ^2))
    C2 = -(fₓ^2*(-β*(b*dₙ + dₜ)*b^2*dₜ*fₓ^2 - α*(b^2*f_ψ^2 + 2*b*β*dₙ*dₜ + 4*β*dₜ^2))*kₙ^2*cos_θ^5 - (-(b^2*dₜ - 3*b*dₙ - 2*dₜ)*β*b*dₜ*fₓ^4 + (-f_ψ^2*b^3 + 2*(-β*dₜ^2 + f_ψ^2)*b + 2*β*dₙ*dₜ)*α*fₓ^2 - 2*b*f_ψ^2*α^2)*sin_θ*kₙ^2*cos_θ^4 + ((-(b^3*dₙ - b^2*dₜ + 3*b*dₙ + dₜ)*β*dₜ*fₓ^4 - 3*α*(-1/3*b^2*f_ψ^2 + 2/3*dₙ*dₜ*β*b + 2*dₜ^2*β + 1/3*f_ψ^2)*fₓ^2 - 2*f_ψ^2*α^2)*kₙ*sin_θ^2 - fₓ^2*((b^2*f_ψ^2 + b*β*dₙ*dₜ + 2*β*dₜ^2)*fₓ^2 + 2*α*f_ψ^2))*kₙ*cos_θ^3 - kₙ*((-β*(b^3*dₜ - 3*b^2*dₙ - b*dₜ - dₙ)*dₜ*fₓ^4 + α*(-f_ψ^2*b^3 + (-2*β*dₜ^2 + f_ψ^2)*b + 2*β*dₙ*dₜ)*fₓ^2 - 2*b*f_ψ^2*α^2)*kₙ*sin_θ^3 - fₓ^2*((f_ψ^2*b^3 - 2*(-β*dₜ^2 + f_ψ^2)*b - β*dₙ*dₜ)*fₓ^2 + 2*b*f_ψ^2*α)*sin_θ)*cos_θ^2 + (-2*(-β*(b^2*dₜ - 3/2*dₙ*b - 1/2*dₜ)*dₜ*fₓ^4 + α*(-b^2*f_ψ^2 + dₜ^2*β + 1/2*f_ψ^2)*fₓ^2 + f_ψ^2*α^2)*kₙ^2*sin_θ^4 + 2*fₓ^2*((b^2*f_ψ^2 - 1/2*dₙ*dₜ*β*b - dₜ^2*β - 1/2*f_ψ^2)*fₓ^2 - α*f_ψ^2)*kₙ*sin_θ^2 - f_ψ^2*fₓ^4)*cos_θ - sin_θ*fₓ^2*(-((b*dₜ - dₙ)*β*dₜ*fₓ^2 + b*f_ψ^2*α)*kₙ^2*sin_θ^4 - 2*((((2*β*dₜ^2 + f_ψ^2)*b - β*dₙ*dₜ)*fₓ^2)/2 + b*f_ψ^2*α)*kₙ*sin_θ^2 - b*f_ψ^2*fₓ^2))/((fₓ^2*kₙ*(b^2*f_ψ^2 + 2*β*dₜ^2)*cos_θ^2 + 2*sin_θ*cos_θ*b*f_ψ^2*fₓ^2*kₙ + 2*kₙ*((dₜ^2*β + f_ψ^2/2)*fₓ^2 + α*f_ψ^2)*sin_θ^2 + f_ψ^2*fₓ^2)*(kₙ*(b^2*fₓ^2 + 2*α)*cos_θ^2 + 2*sin_θ*cos_θ*b*fₓ^2*kₙ + (kₙ*sin_θ^2 + 1)*fₓ^2))
    psi2 = A2 .+ B2*T₁ .+ C2*T₂

    A3 = 0
    B3 = (((b^2*f_ψ^2 + 2*β*dₜ^2)*kₙ + f_ψ^2)*cos_θ + sin_θ*kₙ*(b*f_ψ^2 + 2*β*dₙ*dₜ))*α/(cos_θ^2*b^2*f_ψ^2*fₓ^2*kₙ + 2*sin_θ*cos_θ*b*f_ψ^2*fₓ^2*kₙ + 2*(fₓ^2/2 + α)*f_ψ^2*kₙ*sin_θ^2 + 2*β*dₜ^2*fₓ^2*kₙ + f_ψ^2*fₓ^2)
    C3 = -(((b^2*f_ψ^2 + 2*β*dₜ^2)*kₙ + f_ψ^2)*cos_θ + sin_θ*kₙ*(b*f_ψ^2 + 2*β*dₙ*dₜ))*α/(cos_θ^2*b^2*f_ψ^2*fₓ^2*kₙ + 2*sin_θ*cos_θ*b*f_ψ^2*fₓ^2*kₙ + 2*(fₓ^2/2 + α)*f_ψ^2*kₙ*sin_θ^2 + 2*β*dₜ^2*fₓ^2*kₙ + f_ψ^2*fₓ^2)
    x_D = A3 .+ B3*T₁ .+ C3*T₂

    A4 = α*(FC*cos_θ^2*b^2*kₙ + 2*FC*sin_θ*cos_θ*b*kₙ + FC*kₙ*sin_θ^2 + FC)/(sin_θ^2*fₓ^2*kₙ + 2*sin_θ*cos_θ*b*fₓ^2*kₙ + kₙ*(b^2*fₓ^2 + 2*α)*cos_θ^2 + fₓ^2)
    B4 = α*((kₙ + 1)*sin_θ + b*kₙ*cos_θ)/(sin_θ^2*fₓ^2*kₙ + 2*sin_θ*cos_θ*b*fₓ^2*kₙ + kₙ*(b^2*fₓ^2 + 2*α)*cos_θ^2 + fₓ^2)
    C4 = α*((kₙ + 1)*sin_θ + b*kₙ*cos_θ)/(sin_θ^2*fₓ^2*kₙ + 2*sin_θ*cos_θ*b*fₓ^2*kₙ + kₙ*(b^2*fₓ^2 + 2*α)*cos_θ^2 + fₓ^2)
    y_D = A4 .+ B4*T₁ .+ C4*T₂

    A5 = 0
    B5 = (2*dₙ*kₙ*(fₓ^2/2 + α)*sin_θ^2 - kₙ*((b^2*dₜ - 2*b*dₙ - dₜ)*fₓ^2 - 2*α*dₜ)*cos_θ*sin_θ + (kₙ*(b*dₙ + 2*dₜ)*b*cos_θ^2 - b*dₜ*kₙ + dₙ)*fₓ^2)*β/(2*(fₓ^2/2 + α)*f_ψ^2*kₙ*sin_θ^2 + 2*sin_θ*cos_θ*b*f_ψ^2*fₓ^2*kₙ + fₓ^2*(cos_θ^2*b^2*f_ψ^2*kₙ + 2*kₙ*dₜ^2*β + f_ψ^2))
    C5 = -(2*dₙ*kₙ*(fₓ^2/2 + α)*sin_θ^2 - kₙ*((b^2*dₜ - 2*b*dₙ - dₜ)*fₓ^2 - 2*α*dₜ)*cos_θ*sin_θ + (kₙ*(b*dₙ + 2*dₜ)*b*cos_θ^2 - b*dₜ*kₙ + dₙ)*fₓ^2)*β/(2*(fₓ^2/2 + α)*f_ψ^2*kₙ*sin_θ^2 + 2*sin_θ*cos_θ*b*f_ψ^2*fₓ^2*kₙ + fₓ^2*(cos_θ^2*b^2*f_ψ^2*kₙ + 2*kₙ*dₜ^2*β + f_ψ^2))
    psi_D = A5 .+ B5*T₁ .+ C5*T₂

    return [psi1, psi2, x_D, y_D, psi_D]
end

begin
    fig_psi1_diff = Figure(size = (720, 720))
    ax_psi1_diff = Axis(fig_psi1_diff[1, 1], title="ψ₁-ψ₁ₗᵢₙ", xlabel="T₁", ylabel="T₂", aspect = 1)
    co = contourf!(ax_psi1_diff, T1_values, T2_values, psi1_values .- psi1_values_lin)
    Colorbar(fig_psi1_diff[1, 2], co)
    save(joinpath(path, "figures\\Equilibrium_5DOF" * freetext * "\\psi_1_diff.png"), fig_psi1_diff)

    fig_psi2_diff = Figure(size = (720, 720))
    ax_psi2_diff = Axis(fig_psi2_diff[1, 1], title="ψ₂-ψ₂ₗᵢₙ", xlabel="T₁", ylabel="T₂", aspect = 1)
    co = contourf!(ax_psi2_diff, T1_values, T2_values, psi2_values .- psi2_values_lin)
    Colorbar(fig_psi2_diff[1, 2], co)
    save(joinpath(path, "figures\\Equilibrium_5DOF" * freetext * "\\psi_2_diff.png"), fig_psi2_diff)

    fig_x_D_diff = Figure(size = (720, 720))
    ax_x_D_diff = Axis(fig_x_D_diff[1, 1], title="x_D-x_Dₗᵢₙ", xlabel="T₁", ylabel="T₂", aspect = 1)
    co = contourf!(ax_x_D_diff, T1_values, T2_values, x_D_values .- x_D_values_lin)
    Colorbar(fig_x_D_diff[1, 2], co)
    save(joinpath(path, "figures\\Equilibrium_5DOF" * freetext * "\\x_D_diff.png"), fig_x_D_diff)

    fig_y_D_diff = Figure(size = (720, 720))
    ax_y_D_diff = Axis(fig_y_D_diff[1, 1], title="y_D-y_Dₗᵢₙ", xlabel="T₁", ylabel="T₂", aspect = 1)
    co = contourf!(ax_y_D_diff, T1_values, T2_values, y_D_values .- y_D_values_lin)
    Colorbar(fig_y_D_diff[1, 2], co)
    save(joinpath(path, "figures\\Equilibrium_5DOF" * freetext * "\\y_D_diff.png"), fig_y_D_diff)

    fig_psi_D_diff = Figure(size = (720, 720))
    ax_psi_D_diff = Axis(fig_psi_D_diff[1, 1], title="ψ_D-ψ_Dₗᵢₙ", xlabel="T₁", ylabel="T₂", aspect = 1)
    co = contourf!(ax_psi_D_diff, T1_values, T2_values, psi_D_values .- psi_D_values_lin)
    Colorbar(fig_psi_D_diff[1, 2], co)
    save(joinpath(path, "figures\\Equilibrium_5DOF" * freetext * "\\psi_D_diff.png"), fig_psi_D_diff)
end

# Comparison with 2 DOF model
begin
    two_dof_data = jldopen(joinpath(path, "data\\static\\Equilibrium_2DOF.jld2"), "r")
    sol_two_dof = two_dof_data["solutions"]
    close(two_dof_data)
    psi_two_dof = [sol.psi for sol in sol_two_dof]
    T_two_dof = [sol.T for sol in sol_two_dof]     
    fig_two_dof = Figure(size = (1280, 720))
    ax_two_dof = Axis(fig_two_dof[1, 1], title="Comparison between 2 DOF & 5 DOF", xlabel="T", ylabel="ψ")
    lines!(ax_two_dof, T_two_dof, psi_two_dof, label="2 DOF", color=:blue)
    for i in eachindex(T1_values)
        if (abs(T1_values[i] - T2_values[i]) < 1e-8)
            scatter!(ax_two_dof, T1_values[i], psi1_values[i], color=:red)
        end
    end
    save(joinpath(path, "figures\\Equilibrium_5DOF" * freetext * "\\two_dof.png"), fig_two_dof)
end
=#