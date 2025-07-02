using JLD2, GLMakie, Statistics, Colors, FFTW, Printf
include("..\\src\\StaticModels.jl")
include("..\\src\\DynamicModels.jl")

path = abspath(@__DIR__, "..")
## CONFIGURATION
begin 
    ## Frequency sweep parameters
    # Frequency sweep range
    freq_range = range(0.9, 1.5, length = 121) 
    # Excitation
    Mₑ = 1 / (4e3) # Dimensionless (divided by k_ψ)
    # Number of cycles for the frequency sweep
    n_cycles = 200 
    # Not the exact timestep of the integration
    δt = 1e-3 

    # Mode of the frequency sweep
    #    1 - interactive execution from equilibrium solution file
    #    2 - batch execution from equilibrium solution file
    #    3 - batch execution from user-defined initial conditions
    mode = 2
    
    # File name of the static equilibrium solution to be used for the frequency sweep
    # Required for modes 1 and 2
    eq_pos_file = "Equilibrium_5DOF_final_reduced"
    # Index of the static solution to be used for the frequency sweep
    idx_selected = 7 # Starting index for mode 1
    static_indices = [7 91791 533652 734134 1024500] # Indices mode 2

    ## Model data
    # Required for mode 3
    # Parameters
    pm_temp = TBM_5DOF_Parameters(
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
    # Tangential forces (determines initial position)
    T₁₀ = 1e-4
    T₂₀ = 1e-4 

    ## File management and plotting
    # If true, data from previous executions will be looked up
    # and loaded if available, to continue where it was terminated
    load_data = false
    # If true, the frequency sweep data is saved to a file
    save_data = true 
    # If true, the figures are saved
    save_figures = true 
    # If true, display figures
    display_figures = false
    # Freetext of the directory where the data and figures are saved
    freetext = "_final_less_Me"
    
    ## Frequency sweep settings
    # If true, the initial conditions are taken from the static solution
    # If false, the initial conditions are taken from the previous frequency sweep
    # Only available for modes 1 and 2
    from_static = false 
    # If true, the Fourier analysis is performed on the frequency sweep data
    # bear in mind that this can take a long time to compute
    fourier = true 
    # If true, the integration is carried out in a more memory-efficienty way
    # without storing the entire response in memory
    # If false, the response is also stored in a file(response.jld2) with 
    # reduced precision as long as "save_data" is true
    ram_efficient = false 
    # DOFs to plot
    dofs = [1 2 3 4 5] 
end

## FREQUENCY SWEEP FUNCTION
function frequency_sweep(u₀_input, wₜ₁₀_input, wₜ₂₀_input)
    # Copy values to ensure the originals are not modified
    # and avoid triggering the observables
    u₀ = vcat(u₀_input, zeros(5))
    wₜ₁₀ = wₜ₁₀_input
    wₜ₂₀ = wₜ₂₀_input

    amplitude = zeros(length(dofs), length(freq_range))
    avg_response = zeros(length(dofs), length(freq_range))

    save_data && mkpath(joinpath(path, "data\\" * working_dir))
    save_figures && mkpath(joinpath(path, "figures\\" * working_dir))

    f_path_data = joinpath(path, "data\\" * working_dir * "frequency_sweep.jld2")

    # Perform frequency sweep for the selected point
    
    tf₋₁ = 0.0
    if (isfile(f_path_data) && load_data)
        ref_data = load(f_path_data)
        freq_range_ref = ref_data["freq_range"]
        avg_response_ref = ref_data["avg_response"]
        amplitude_ref = ref_data["amplitude"]
        println("[info] Loaded frequency sweep data from file: $f_path_data")
        @assert length(freq_range_ref) == length(freq_range) &&
                length(avg_response_ref) == length(avg_response) && 
                length(amplitude_ref) == length(amplitude)
        avg_response = avg_response_ref
        amplitude = amplitude_ref
        starting_index = findfirst(abs.(avg_response[1, :]) .< 1e-12)
        if starting_index === nothing
            starting_index = 1
            println("[info] No frequency sweep data found, starting from the beginning...")
        else
            if !from_static
                u₀ = ref_data["u₀"]
                wₜ₁₀ = ref_data["wₜ₁₀"]
                wₜ₂₀ = ref_data["wₜ₂₀"]
            end
            tf₋₁ = ref_data["tf₋₁"]
            println("[info] Frequency sweep data found, starting from index $starting_index of $(lastindex(freq_range))")
        end
    else
        starting_index = 1
        println("[info] Performing frequency sweep from the start...")
    end
    flush(stdout)

    # Initialize figures and axes for plotting
    if !ram_efficient
        fig_frequency_sweep = Vector{Figure}(undef, length(dofs))
        ax_frequency_sweep = Vector{Axis}(undef, length(dofs))
        for (idx_dof, dof) in enumerate(dofs)
            fig_frequency_sweep[idx_dof] = Figure(size = (1280, 720))
            ax_frequency_sweep[idx_dof] = Axis(fig_frequency_sweep[idx_dof][1, 1], 
                title="$(dof_names[dof]) vs time", 
                xlabel="t", 
                ylabel="$(dof_names[dof])"; 
                xgridvisible = true, 
                ygridvisible = true, 
                xgridcolor = :gray, 
                ygridcolor = :gray
            )
        end 

        fig_contact_cycles_1 = Figure(size = (720, 720))
        ax_contact_cycles_1 = Axis(fig_contact_cycles_1[1, 1], 
            title="Contact cycles at left contact", 
            xlabel="xₜ₁", 
            ylabel="T₁"; 
            xgridvisible = true, 
            ygridvisible = true, 
            xgridcolor = :gray, 
            ygridcolor = :gray
        )

        fig_contact_cycles_2 = Figure(size = (720, 720))
        ax_contact_cycles_2 = Axis(fig_contact_cycles_2[1, 1], 
            title="Contact cycles at right contact", 
            xlabel="xₜ₂", 
            ylabel="T₂"; 
            xgridvisible = true, 
            ygridvisible = true, 
            xgridcolor = :gray, 
            ygridcolor = :gray
        )

        fig_contact_forces_1 = Figure(size = (1280, 720))
        ax_contact_forces_1 = Axis(fig_contact_forces_1[1, 1], 
            title="Contact forces at left contact", 
            xlabel="t", 
            ylabel="T₁"; 
            xgridvisible = true, 
            ygridvisible = true, 
            xgridcolor = :gray, 
            ygridcolor = :gray
        )

        fig_contact_forces_2 = Figure(size = (1280, 720))
        ax_contact_forces_2 = Axis(fig_contact_forces_2[1, 1], 
            title="Contact forces at right contact", 
            xlabel="t", 
            ylabel="T₂"; 
            xgridvisible = true, 
            ygridvisible = true, 
            xgridcolor = :gray, 
            ygridcolor = :gray
        )
    end

    if fourier
        fig_fourier = Vector{Figure}(undef, length(dofs))
        ax_fourier = Vector{Axis}(undef, length(dofs))
        for (idx_dof, dof) in enumerate(dofs)
            fig_fourier[idx_dof] = Figure(size = (1280, 720))
            ax_fourier[idx_dof] = Axis(fig_fourier[idx_dof][1, 1], 
                title="Spectrum of DOF $(dof_names[dof])", 
                xlabel="Frequency"; 
                xgridvisible = true, 
                ygridvisible = true, 
                xgridcolor = :gray, 
                ygridcolor = :gray
            )
            xlims!(ax_fourier[idx_dof], 0.0, 10.0) 
        end
    end

    # Frequency sweep loop
    for j in starting_index:lastindex(freq_range)
        tf = 2π*n_cycles / freq_range[j]
        n_steps = round(Int, tf/δt) + 1
        Δt = tf / (n_steps - 1)

        # Set contacts to starting value    
        pm.feq = freq_range[j]
        pm.friction_law_1.wₜ = wₜ₁₀
        pm.friction_law_2.wₜ = wₜ₂₀
        
        # The last 10% of the response is considered the final (permanent) response
        n_final = round(Int, 0.9 * n_steps)
        u_final = zeros(10, n_steps - n_final + 1)
        t_final = zeros(n_steps - n_final + 1)

        if !ram_efficient
            t = range(0, tf, length = n_steps)
            # Initialize internal variables of parameters structure
            pm.T₁, pm.T₂, pm.N₁, pm.N₂, pm.xₜ₁, pm.xₜ₂, pm.xₙ₁, pm.xₙ₂ =  ntuple(_ -> zeros(n_steps), 8)
            μ = pm.friction_law_1.μ
            pm.index=Ref(1)
        
            u = zeros(10, n_steps)  
            # Set initial condition for the DOFs
            u[:, 1] = u₀

            # Integration
            for i in 1:(n_steps - 1)
                u[:, i+1] = rk4(five_dof_model, u[:, i], pm, t[i], Δt)
            end

            u_final = u[:, n_final:end]
            t_final = t[n_final:end]

            # Plot the response
            idx_plot = vcat(collect(1:10:n_steps-1), n_steps) 
            for idx_dof in eachindex(dofs)
                # Plot the frequency sweep response
                lines!(ax_frequency_sweep[idx_dof], 
                    Float32.(t[idx_plot].+tf₋₁), 
                    Float32.(u[idx_dof, idx_plot]), 
                    label="f = $(round(freq_range[j], sigdigits=4))", 
                    color=pcols[j]
                )
            end

            # Plot the contact cycles
            idx_plot_contact = vcat(collect(n_final:10:n_steps-2), n_steps-1)
            lines!(ax_contact_cycles_1, 
                Float32.(pm.xₜ₁[idx_plot_contact]), 
                Float32.(pm.T₁[idx_plot_contact]), 
                label="f = $(round(freq_range[j], sigdigits=4))", 
                color=pcols[j]
            )
            lines!(ax_contact_cycles_2, 
                Float32.(pm.xₜ₂[idx_plot_contact]), 
                Float32.(pm.T₂[idx_plot_contact]), 
                label="f = $(round(freq_range[j], sigdigits=4))", 
                color=pcols[j]
            )

            # Plot contact forces
            lines!(ax_contact_forces_1,
                Float32.(t[idx_plot[1:end-1]].+tf₋₁), 
                Float32.(pm.T₁[idx_plot[1:end-1]]), 
                label="f = $(round(freq_range[j], sigdigits=4))", 
                color=pcols[j]
            )
            lines!(ax_contact_forces_1,
                Float32.(t[idx_plot[1:end-1]].+tf₋₁), 
                Float32.(μ*pm.N₁[idx_plot[1:end-1]]), 
                label="±μN₁", 
                linestyle = :dash,
                color=pcols[j]
            )
            lines!(ax_contact_forces_1,
                Float32.(t[idx_plot[1:end-1]].+tf₋₁), 
                Float32.(-μ*pm.N₁[idx_plot[1:end-1]]), 
                label="±μN₁", 
                linestyle = :dash,
                color=pcols[j]
            )

            lines!(ax_contact_forces_2,
                Float32.(t[idx_plot[1:end-1]].+tf₋₁), 
                Float32.(pm.T₂[idx_plot[1:end-1]]), 
                label="f = $(round(freq_range[j], sigdigits=4))", 
                color=pcols[j]
            )
            lines!(ax_contact_forces_2,
                Float32.(t[idx_plot[1:end-1]].+tf₋₁), 
                Float32.(μ*pm.N₂[idx_plot[1:end-1]]), 
                label="±μN₂", 
                linestyle = :dash,
                color=pcols[j]
            )
            lines!(ax_contact_forces_2,
                Float32.(t[idx_plot[1:end-1]].+tf₋₁), 
                Float32.(-μ*pm.N₂[idx_plot[1:end-1]]), 
                label="±μN₂", 
                linestyle = :dash,
                color=pcols[j]
            )
        else 
            u = zeros(10, 1)

            # Set initial conditions for the integration
            tᵢ = 0.0
            uᵢ₋₁ = u₀

            # Integration
            for i in 2:n_steps
                u = rk4(five_dof_model, uᵢ₋₁, pm, tᵢ, Δt)
                tᵢ += Δt
                uᵢ₋₁ = u
                if i >= n_final
                    k = i - (n_final - 1)
                    u_final[:, k] = u
                    t_final[k] = tᵢ
                end
            end
        end
        println("[info] Time integration for f=$(round(freq_range[j], sigdigits=4)) completed.")
        flush(stdout)

        # Final response analysis
        if fourier
            freq_fourier = 2π*rfftfreq(length(t_final), 1/Δt)
            plan = plan_rfft(u_final[1, :], flags=FFTW.ESTIMATE, timelimit=Inf)
            U = zeros(length(dofs), length(freq_fourier))
        end

        for (idx_dof, dof) in enumerate(dofs)
            avg_response[idx_dof, j] = mean(u_final[dof, :])
            amplitude[idx_dof, j] = abs(maximum(u_final[dof, :]) - minimum(u_final[dof, :]))
            
            if fourier
                U[idx_dof, :] = abs.(plan*(u_final[dof, :].-mean(u_final[dof, :])))
                lines!(ax_fourier[idx_dof], 
                    freq_fourier, 
                    U[idx_dof, :], 
                    label="f = $(round(freq_range[j], sigdigits=4))", 
                    color=pcols[j]
                )
            end
        end
        
        # Update initial conditions for the next frequency
        if !from_static
            u₀ = u[:, end]
            wₜ₁₀ = pm.friction_law_1.wₜ
            wₜ₂₀ = pm.friction_law_2.wₜ
        end
        tf₋₁ += tf
        
        # Save the frequency sweep data
        # u₀ wₜ₁₀ wₜ₂₀ are saved to allow for continuation of the frequency sweep in
        # case the script is interrupted or restarted        
        if save_data
            @save f_path_data avg_response amplitude freq_range u₀ wₜ₁₀ wₜ₂₀ tf₋₁ dofs

            if !ram_efficient
                f_path_response = joinpath(path, "data\\" * working_dir * "response.jld2")
                key_tf = "tf_f_$(round(freq_range[j], sigdigits=4))"
                key_u = "u_f_$(round(freq_range[j], sigdigits=4))"
                key_nsteps = "n_steps_f_$(round(freq_range[j], sigdigits=4))" # Needed to be capable of reconstructing the solution with reduced size
                key_nfinal = "n_final_f_$(round(freq_range[j], sigdigits=4))"
                key_xt1 = "x_t1_f_$(round(freq_range[j], sigdigits=4))"
                key_xt2 = "x_t2_f_$(round(freq_range[j], sigdigits=4))"
                key_T1 = "T1_f_$(round(freq_range[j], sigdigits=4))"
                key_T2 = "T2_f_$(round(freq_range[j], sigdigits=4))"
                key_N1 = "N1_f_$(round(freq_range[j], sigdigits=4))"
                key_N2 = "N2_f_$(round(freq_range[j], sigdigits=4))"
                jldopen(f_path_response, "a+") do file
                    file[key_tf] = tf
                    file[key_u] = Float32.(u[1:5, idx_plot]) # Save only points every 10 steps to reduce file size
                    file[key_xt1] = Float32.(pm.xₜ₁[idx_plot_contact])
                    file[key_xt2] = Float32.(pm.xₜ₂[idx_plot_contact])
                    file[key_T1] = Float32.(pm.T₁[idx_plot_contact])
                    file[key_T2] = Float32.(pm.T₂[idx_plot_contact])
                    file[key_N1] = Float32.(pm.N₁[idx_plot_contact])
                    file[key_N2] = Float32.(pm.N₂[idx_plot_contact])
                    file[key_nsteps] = n_steps
                    file[key_nfinal] = n_final
                end
            end

            if fourier
                f_path_fourier = joinpath(path, "data\\" * working_dir * "fourier.jld2")
                key_freq = "freq_fourier_f_$(round(freq_range[j], sigdigits=4))"
                key_U = "U_f_$(round(freq_range[j], sigdigits=4))"
                jldopen(f_path_fourier, "a+") do file
                    file[key_freq] = freq_fourier
                    file[key_U] = U
                end
            end

        end
        
    end
               
    # Frequency sweep figures
    fig_freq_avg_response = Vector{Figure}(undef, length(dofs))
    ax_freq_avg_response = Vector{Axis}(undef, length(dofs))

    fig_freq_ampl = Vector{Figure}(undef, length(dofs))
    ax_freq_ampl = Vector{Axis}(undef, length(dofs))

    for (idx_dof, dof) in enumerate(dofs)
        f_freetext = ""
        if dof != 1
            f_freetext = "_dof_$dof"
        end

        fig_freq_avg_response[idx_dof] = Figure(size = (1280, 720))
        ax_freq_avg_response[idx_dof] = Axis(fig_freq_avg_response[idx_dof][1, 1], 
            title= "Average response for $(dof_names[dof])",
            xlabel= "f", 
            ylabel= "<$(dof_names[dof])>"; 
            xgridvisible = true, 
            ygridvisible = true, 
            xgridcolor = :gray,
            ygridcolor = :gray
        )
        scatterlines!(ax_freq_avg_response[idx_dof], 
            freq_range, 
            avg_response[idx_dof, :], 
            color = :black, 
            markersize = 5
        )

        fig_freq_ampl[idx_dof] = Figure(size = (1280, 720))
        ax_freq_ampl[idx_dof] = Axis(fig_freq_ampl[idx_dof][1, 1], 
            title= "Response amplitude for $(dof_names[dof])", 
            xlabel= "f", 
            ylabel= "|max($(dof_names[dof]))-min($(dof_names[dof]))|"; 
            xgridvisible = true, 
            ygridvisible = true, 
            xgridcolor = :gray, 
            ygridcolor = :gray
        )
        scatterlines!(ax_freq_ampl[idx_dof], 
            freq_range, 
            amplitude[idx_dof, :], 
            color = :black, 
            markersize = 5
        )
        
        ram_efficient || axislegend(ax_frequency_sweep[idx_dof], merge = true, unique = true)
        fourier && axislegend(ax_fourier[idx_dof], merge = true, unique = true)

        if save_figures
            save(joinpath(path, "figures\\" * working_dir * "freq_avg_response" * f_freetext * ".png"), fig_freq_avg_response[idx_dof])
            save(joinpath(path, "figures\\" * working_dir * "freq_ampl" * f_freetext * ".png"), fig_freq_ampl[idx_dof])
            (mode == 1 || mode == 2) && save(joinpath(path, "figures\\" * working_dir * "static_solution" * f_freetext * ".png"), fig_static_solution)
            ram_efficient || save(joinpath(path, "figures\\" * working_dir * "frequency_sweep" * f_freetext * ".png"), fig_frequency_sweep[idx_dof])
            fourier && save(joinpath(path, "figures\\" * working_dir * "fourier" * f_freetext * ".png"), fig_fourier[idx_dof])
        end

        if display_figures
            display(GLMakie.Screen(), fig_freq_avg_response[idx_dof])
            display(GLMakie.Screen(), fig_freq_ampl[idx_dof])
            ram_efficient || display(GLMakie.Screen(), fig_frequency_sweep[idx_dof])
            fourier && display(GLMakie.Screen(), fig_fourier[idx_dof])
        end
    end
    
    # Contact cycles figures
    if !ram_efficient
        axislegend(ax_contact_cycles_1, merge = true, unique = true, position = :lt)
        axislegend(ax_contact_cycles_2, merge = true, unique = true, position = :lt)
        if save_figures
            save(joinpath(path, "figures\\" * working_dir * "contact_1.png"), fig_contact_cycles_1)
            save(joinpath(path, "figures\\" * working_dir * "contact_2.png"), fig_contact_cycles_2)
        end
        if display_figures
            display(GLMakie.Screen(), fig_contact_cycles_1)
            display(GLMakie.Screen(), fig_contact_cycles_2)
            display(GLMakie.Screen(), fig_contact_forces_1)
            display(GLMakie.Screen(), fig_contact_forces_2)
        end
    end

end

function load_data!()
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
    pm.Mₑ = Mₑ # Set the excitation momentum

    println("[info] Static equilibrium file loaded.")

    # Static equilibrium positions figure
    global fig_static_solution = Figure(size = (600, 600))
    global ax_static_solution = Axis(fig_static_solution[1, 1], 
        title="Static equilibrium positions", 
        xlabel="T₁", 
        ylabel="T₂", 
        aspect = 1,  
        xgridvisible = true, 
        ygridvisible = true, 
        xgridcolor = :gray, 
        ygridcolor = :gray
    )
    co_static_solution = contourf!(ax_static_solution, 
        T1_static, 
        T2_static, 
        ones(n_static); 
        colormap = :grays
    )
    Δ_T1 = (maximum(T1_static) - minimum(T1_static))/20
    Δ_T2 = (maximum(T2_static) - minimum(T2_static))/20
    limits!(ax_static_solution, 
        minimum(T1_static)-Δ_T1, maximum(T1_static)+Δ_T1,
        minimum(T2_static)-Δ_T2, maximum(T2_static)+Δ_T2
    )

end

function get_idx_at(point; tol=0.51)
    # Find nearest grid indices
    # Tolerance is relative to the distance between T values, value of 0.51 should always ensure that the point is found
    ΔT =  abs(T1_static[findfirst(x -> abs(x - T1_static[1]) > 1e-8, T1_static[2:end]) + 1] - T1_static[1])
    for idx in eachindex(T1_static)
            if abs(T1_static[idx] - point[1]) < tol*ΔT && abs(T2_static[idx] - point[2]) < tol*ΔT
            return idx
        end
    end
    println("[warning] Point ($(point[1]), $(point[2])) not found in the surface.")
    return nothing
end

## FREQUENCY SWEEP EXECUTION
begin
    # General utililies
    if from_static
        freetext *= "_from_static"
    end

    dof_names = ["ψ₁", "ψ₂", "x_D", "y_D", "ψ_D"]

    cols = distinguishable_colors(length(freq_range), [RGB(1,1,1), RGB(0,0,0)], dropseed=true)
    pcols = map(col -> (red(col), green(col), blue(col)), cols)

    # Execution of frequency sweep
    if mode == 1 # Interactive execution from equilibrium solution file
        println("[info] Started interactive execution from static equilibrium solution file: $eq_pos_file.jld2")
        # Load static equilibrium data
        load_data!()

        idx_static = Observable(idx_selected)  # Index of T₁ for frequency sweep
        working_dir = "Forced_response" * freetext * "\\idx_$idx_selected\\"

        # Add a scatter marker to show selected point
        pt = Observable(Point2f(T1_static[idx_selected], T2_static[idx_selected]))
        scatter!(ax_static_solution, 
            pt, 
            color=:red, 
            markersize=10
        )
        
        # Perform frequency sweep for the selected point
        @lift begin
            u₀ = u_static[:, $idx_static]
            wₜ₁₀ = wt1_static[$idx_static]
            wₜ₂₀ = wt2_static[$idx_static]
            frequency_sweep(u₀, wₜ₁₀, wₜ₂₀)
        end

        on(events(ax_static_solution).mousebutton) do event
            if (event.button == Mouse.middle || event.button == Mouse.button_5) && event.action == Mouse.press
                point = mouseposition(ax_static_solution.scene) 
                i = get_idx_at(point, tol=0.51)
                if i !== nothing
                    pt[] = Point2f(T1_static[i], T2_static[i])  # Update the point position
                    println("[info] Clicked at i=$i, T₁=$(round(T1_static[i], sigdigits=3)), T₂=$(round(T2_static[i], sigdigits=3))")
                    global working_dir = "Forced_response" * freetext * "\\idx_$i\\"
                    idx_static[] = i # Update the index for frequency sweep
                end
            end
        end
        
        display(GLMakie.Screen(), fig_static_solution) 
        println("[info] Frequency sweep completed for static index $(idx_static[])\n",
                "       On figure \"Static equilibrium positions\" click MMB or MB5 at any\n",
                "       point to perform a frequency sweep for that equilibrium position.")

    elseif mode == 2 # Batch execution from equilibrium solution file
        println("[info] Started non-interactive execution from static equilibrium solution file: $eq_pos_file.jld2")
        # Load static equilibrium data
        load_data!()

        for idx_static in static_indices
            global working_dir = "Forced_response" * freetext * "\\idx_$idx_static\\"
            local pt = (T1_static[idx_static], T2_static[idx_static])
            scatter!(ax_static_solution, pt, color=:red, markersize=10)
            
            # Perform frequency sweep for the selected point
            local u₀ = u_static[:, idx_static]
            local wₜ₁₀ = wt1_static[idx_static]
            local wₜ₂₀ = wt2_static[idx_static]
            frequency_sweep(u₀, wₜ₁₀, wₜ₂₀)

            println("[info] Frequency sweep completed for static index $idx_static")
        end
        println("[info] Non-interactive execution completed. Shutting down...")

    elseif mode == 3 # Batch execution from user-defined initial conditions
        println("[info] Started batch execution from user-defined initial conditions.")

        # Set initial conditions for the integration
        u₀ = equilibrium_pos_5DOF(pm_temp, T₁₀, T₂₀, zeros(5))
        xₜ₁₀, _, xₜ₂₀, _ = contact_points_displacements(u₀, pm_temp)
        kₜ = pm_temp.friction_law_1.kₜ
        wₜ₁₀ = xₜ₁₀ - T₁₀ / kₜ
        wₜ₂₀ = xₜ₂₀ - T₂₀ / kₜ
        pm = pm_temp
        pm.Mₑ = Mₑ

        working_dir =  "Forced_response" * freetext * "\\user_defined\\"

        frequency_sweep(u₀, wₜ₁₀, wₜ₂₀)

        println("[info] Batch execution from user-defined initial conditions completed. Shutting down...")
    else
        error("[error] Invalid mode: $mode. Use 1, 2 or 3.")
    end
    
end