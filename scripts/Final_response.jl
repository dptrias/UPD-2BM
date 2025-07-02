using JLD2, GLMakie, Julianim, Statistics, MathTeXEngine
include("..\\src\\DynamicModels.jl")

path = abspath(@__DIR__, "..")
begin 
    set_publication_theme!()
    update_theme!(
    fonts = Attributes(
            :bold => texfont(:bold),
            :bolditalic => texfont(:bolditalic),
            :italic => texfont(:italic),
            :regular => texfont(:regular)
        )
    )

    # Static equilibrium file name
    eq_pos_file = "Equilibrium_5DOF_final_even_lower_def_reduced"

    # Excitation frequency
    freq_range = [0.65, 0.7, 0.75, 0.85, 0.9, 0.95, 1.05, 1.45, 1.5, 1.55]
    # Number of cycles for the frequency sweep
    n_cycles = 200
    # Not the exact timestep of the integration
    δt = 1e-3 
    # Excitation momentum
    Mₑ = 0.25e-3

    save_data = true
    save_figures = true
    freetext = "_final_less_Me"
end

function load_static_data!(eq_pos_file)
    # Get data from static equilibrium
    static_data = jldopen(joinpath(path, "data\\static\\" * eq_pos_file *".jld2"), "r")
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

    global pm = pm_static

    println("[info] Static equilibrium file loaded: $eq_pos_file")
end

begin
    println("[info] Starting final response analysis...")
    # Load the static equilibrium data
    load_static_data!(eq_pos_file)

    for freq in freq_range
        println("[info] Processing frequency: $freq")
        flush(stdout)
        
        # Create a directory for the current frequency
        working_dir = "Final_response" * freetext * "\\freq_$(replace(string(freq), "." => "-"))\\"        
        u_mean = zeros(5, n_static)
        T_1_final = zeros(n_static)
        T_2_final = zeros(n_static)

        tf = 2π*n_cycles / freq
        n_steps = round(Int, tf/δt) + 1
        Δt = tf / (n_steps - 1)

        # The last 10% of the response is considered the final (permanent) response
        n_final = round(Int, 0.9 * n_steps)
        u_final = zeros(10, n_steps - n_final + 1)
        t_final = zeros(n_steps - n_final + 1)

        # Set contacts to starting value    
        pm.feq = freq
        pm.Mₑ = Mₑ
        pm.T₁, pm.T₂, pm.N₁, pm.N₂, pm.xₜ₁, pm.xₜ₂, pm.xₙ₁, pm.xₙ₂ =  ntuple(_ -> zeros(n_steps - n_final + 1), 8)
        
        for idx_static in 1:n_static            
            pm.friction_law_1.wₜ = wt1_static[idx_static]
            pm.friction_law_2.wₜ = wt2_static[idx_static]
            
            u = zeros(10, 1)
            
            # Set initial conditions for the integration
            tᵢ = 0.0
            uᵢ₋₁ = vcat(u_static[:, idx_static], zeros(5))
            
            # Integration
            pm.index = nothing
            for i in 2:n_steps
                (i == n_final) && (pm.index = Ref(1)) # Start storing the final response
                u = rk4(five_dof_model, uᵢ₋₁, pm, tᵢ, Δt)
                tᵢ += Δt
                uᵢ₋₁ = u
                if i >= n_final
                    k = i - (n_final - 1)
                    u_final[:, k] = u
                    t_final[k] = tᵢ
                end
            end
            
            println("[info] Time integration completed for static index $idx_static.")
            flush(stdout)

            u_mean[:, idx_static] = mean(u_final[1:5, :], dims = 2)
            T_1_final[idx_static] = mean(pm.T₁)
            T_2_final[idx_static] = mean(pm.T₂)
        
        end

        if save_data
            mkpath(joinpath(path, "data\\" * working_dir))
            f_path_data = joinpath(path, "data\\" * working_dir * "statistics.jld2")
            if isfile(f_path_data)
               println("[warning] File already exists: $f_path_data. Overwriting...")
               rm(f_path_data, force=true)
            end
            # Save the results
            @save f_path_data u_mean T_1_final T_2_final
            println("[info] Data for frequency $freq saved in: $f_path_data")
        end
            
        if !save_figures
            continue # Skip plotting if neither saving nor displaying figures
        end
            
        # Plot the results
        Δ_T1 = (maximum(T1_static) - minimum(T1_static))/20
        Δ_T2 = (maximum(T2_static) - minimum(T2_static))/20
            
        fig_psi1 = Figure()
        ax_psi1 = Axis(fig_psi1[1, 1], 
            # title="ψ₁ vs T₁ & T₂", 
            xlabel=L"T_1", 
            ylabel=L"T_2", 
            aspect = 1
        )
        limits!(ax_psi1, 
            minimum(T1_static)-Δ_T1, maximum(T1_static)+Δ_T1, 
            minimum(T2_static)-Δ_T2, maximum(T2_static)+Δ_T2
        )
        co = contourf!(ax_psi1, 
            T1_static, 
            T2_static, 
            abs.(u_static[1, :] .- u_mean[1, :]),
        )
        Colorbar(fig_psi1[1, 2],
            co,
            label = L"| \psi_1 - < \psi_1 > |" 
        )

        fig_psi2 = Figure()
        ax_psi2 = Axis(fig_psi2[1, 1], 
            # title="ψ₂ vs T₁ & T₂", 
            xlabel=L"T_1", 
            ylabel=L"T_2", 
            aspect = 1
        )
        limits!(ax_psi2,
            minimum(T1_static)-Δ_T1, maximum(T1_static)+Δ_T1, 
            minimum(T2_static)-Δ_T2, maximum(T2_static)+Δ_T2
        )
        co = contourf!(ax_psi2, 
            T1_static,
            T2_static, 
            abs.(u_static[2, :] .- u_mean[2, :]),
        )
        Colorbar(fig_psi2[1, 2],
            co,
            label = L"| \psi_2 - < \psi_2 > |" 
        )
        
        fig_x_D = Figure()
        ax_x_D = Axis(fig_x_D[1, 1], 
            # title="x_D vs T₁ & T₂", 
            xlabel=L"T_1", 
            ylabel=L"T_2", 
            aspect = 1
        )
        limits!(ax_x_D, 
            minimum(T1_static)-Δ_T1, maximum(T1_static)+Δ_T1, 
            minimum(T2_static)-Δ_T2, maximum(T2_static)+Δ_T2
        )
        co = contourf!(ax_x_D, 
            T1_static, 
            T2_static, 
            abs.(u_static[3, :] .- u_mean[3, :]),
        )
        Colorbar(fig_x_D[1, 2],
            co,
            label = L"| x_D - < x_D > |" 
        )
            
        fig_y_D = Figure()
        ax_y_D = Axis(fig_y_D[1, 1], 
            # title="y_D vs T₁ & T₂", 
            xlabel=L"T_1", 
            ylabel=L"T_2", 
            aspect = 1
        )
        limits!(ax_y_D, 
            minimum(T1_static)-Δ_T1, maximum(T1_static)+Δ_T1, 
            minimum(T2_static)-Δ_T2, maximum(T2_static)+Δ_T2
        )
        co = contourf!(ax_y_D, 
            T1_static, 
            T2_static, 
            abs.(u_static[4, :] .- u_mean[4, :]),
        )
        Colorbar(fig_y_D[1, 2],
            co,
            label = L"| y_D - < y_D > |" 
        )
            
        fig_psi_D = Figure()
        ax_psi_D = Axis(fig_psi_D[1, 1], 
            # title="ψ_D vs T₁ & T₂", 
            xlabel=L"T_1", 
            ylabel=L"T_2", 
            aspect = 1
        )
        limits!(ax_psi_D, 
            minimum(T1_static)-Δ_T1, maximum(T1_static)+Δ_T1, 
            minimum(T2_static)-Δ_T2, maximum(T2_static)+Δ_T2
        )
        co = contourf!(ax_psi_D, 
            T1_static, 
            T2_static, 
            abs.(u_static[5, :] .- u_mean[5, :]),
        )
        Colorbar(fig_psi_D[1, 2],
            co,
            label = L"| \psi_D - < \psi_D > |" 
        )
            
        fig_final = Figure()
        ax_final = Axis(fig_final[1, 1], 
            # title="ψ₁ vs T₁ & T₂", 
            xlabel=L"T_1", 
            ylabel=L"T_2", 
            aspect = 1
        )
        limits!(ax_final, 
            minimum(T1_static)-Δ_T1, maximum(T1_static)+Δ_T1, 
            minimum(T2_static)-Δ_T2, maximum(T2_static)+Δ_T2
        )
        co = contourf!(ax_final, 
            T1_static, 
            T2_static, 
            colormap = :grays, 
            ones(n_static)
        )
        scatter!(ax_final, 
            T_1_final, 
            T_2_final, 
            markersize = 10,
            color = RED,
            label = L"\text{Final value}"
        )
        legend_final = axislegend(ax_final, 
            position = :rb
        )
        translate!(legend_final.blockscene, Vec3f(-20, 20, 0))
        
        mkpath(joinpath(path, "figures\\" * working_dir))
        save(joinpath(path, "figures\\" * working_dir * "psi_1.png"), fig_psi1)
        save(joinpath(path, "figures\\" * working_dir * "psi_2.png"), fig_psi2)
        save(joinpath(path, "figures\\" * working_dir * "x_D.png"), fig_x_D)
        save(joinpath(path, "figures\\" * working_dir * "y_D.png"), fig_y_D)
        save(joinpath(path, "figures\\" * working_dir * "psi_D.png"), fig_psi_D)
        save(joinpath(path, "figures\\" * working_dir * "final.png"), fig_final)
        println("[info] Figures for frequency $freq saved in: $(joinpath(path, "figures\\" * working_dir))")

        println("[info] Frequency $freq processed.")
    end
    
    println("[info] Final response analysis completed. Shutting down...")
    
end