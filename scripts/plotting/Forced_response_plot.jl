using JLD2, GLMakie, Statistics, Julianim, MathTeXEngine
include("..\\..\\src\\Common.jl")

path = abspath(@__DIR__, "..\\..")
## INITIALIZATION
begin
    case = "final_less_Me"
    
    # freq_plots = [1.0 1.05 1.09 1.2 1.3 1.5]
    freq_plots =[1.1 1.2 1.24 1.27 1.4]

    indices = [7 91791 533652 734134 1024500]

    set_publication_theme!()
    update_theme!(
    fonts = Attributes(
        :bold => texfont(:bold),
        :bolditalic => texfont(:bolditalic),
        :italic => texfont(:italic),
        :regular => texfont(:regular)
    ))
    dof_names = ["\\psi_1","\\psi_2", "x_D", "y_D", "\\psi_D"]
    dofs = [1 2 3 4 5] # DOFs to plot
    eq_pos_file = "Equilibrium_5DOF_final_reduced"

end

function load_static_data!(eq_pos_file)
    # Get data from static equilibrium
    static_data = jldopen(abspath(path, "..\\data\\static\\" * eq_pos_file *".jld2"), "r")
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
    global Δ_T1 = (maximum(T1_static) - minimum(T1_static))/20
    global Δ_T2 = (maximum(T2_static) - minimum(T2_static))/20

    println("[info] Static equilibrium file loaded: $eq_pos_file")
end

begin
    load_static_data!(eq_pos_file)
    fig_freq_avg_response = Vector{Figure}(undef, length(dofs))
    ax_freq_avg_response = Vector{Axis}(undef, length(dofs))
    fig_freq_ampl = Vector{Figure}(undef, length(dofs))
    ax_freq_ampl = Vector{Axis}(undef, length(dofs))
    fig_static_pos = Figure()
    ax_static_pos = Axis(fig_static_pos[1, 1], 
        # title="ψ₁ vs T₁ & T₂", 
        xlabel=L"T_1", 
        ylabel=L"T_2", 
        aspect = 1
    )
    limits!(ax_static_pos, 
        minimum(T1_static)-Δ_T1, maximum(T1_static)+Δ_T1, 
        minimum(T2_static)-Δ_T2, maximum(T2_static)+Δ_T2
    )
    co = contourf!(ax_static_pos, 
        T1_static, 
        T2_static, 
        colormap = :grays, 
        ones(n_static)
    )
    
    for (idx_dof, dof) in enumerate(dofs)
        global fig_freq_avg_response[idx_dof] = Figure(size = (1280, 720))
        label = "\$ <$(dof_names[dof])> \$"
        global ax_freq_avg_response[idx_dof] = Axis(fig_freq_avg_response[idx_dof][1, 1], 
            #title="Average response", 
            xlabel=L"f", 
            ylabel = LaTeXString(label)
        )

        global fig_freq_ampl[idx_dof] = Figure(size = (1280, 720))
        label = "\$|\\text{max}($(dof_names[dof]))-\\text{min}($(dof_names[dof]))|\$"
        global ax_freq_ampl[idx_dof] = Axis(fig_freq_ampl[idx_dof][1, 1], 
            # title="Response amplitude", 
            xlabel=L"f", 
            ylabel = LaTeXString(label)
        )
        
    end

    figure_dir = abspath(path, "..\\figures\\Forced_response_$(case)_publication\\")
    for (idx_file, idx_static) in enumerate(indices)
        local static_case = "idx_$(idx_static)"
        mkpath(figure_dir * "$static_case\\")
        local freq_sweep_file = abspath(path, "..\\", "data\\Forced_response_$case\\$static_case\\frequency_sweep.jld2")
        local response_file = abspath(path, "..\\", "data\\Forced_response_$case\\$static_case\\response.jld2")

        if !isfile(freq_sweep_file)
            error("Frequency sweep data file not found: $freq_sweep_file")
        end
        local data = load(freq_sweep_file)
        local freq_range = data["freq_range"]
        local avg_response = data["avg_response"]
        local amplitude = data["amplitude"]

        if !isfile(response_file)
            error("Response data file not found: $response_file")
        end
        response_data = load(response_file)

        
        for (idx_dof, dof) in enumerate(dofs)
            mkpath(figure_dir * "$static_case\\")
            fig_freetext = ""
            if dof != 1
                fig_freetext *= "_dof_$dof"
            end
            fig_freetext *= ".png"

            label = "\$ <$(dof_names[dof])> \$"
            fig_freq_avg_response_indv = Figure(size = (1280, 720))
            ax_freq_avg_response_indv = Axis(fig_freq_avg_response_indv[1, 1], 
                #title="Average response (individual)", 
                xlabel=L"f", 
                ylabel = LaTeXString(label)
            )
            scatterlines!(ax_freq_avg_response_indv, 
                freq_range, 
                avg_response[idx_dof, :], 
                color = :black, 
                markersize = 5
            )
            save(figure_dir * "$static_case\\freq_avg_response" * fig_freetext, fig_freq_avg_response_indv, px_per_unit = 4)

            label = "\$|\\text{max}($(dof_names[dof]))-\\text{min}($(dof_names[dof]))|\$"
            fig_freq_ampl_indv = Figure(size = (1280, 720))
            ax_freq_ampl_indv = Axis(fig_freq_ampl_indv[1, 1], 
                #title="Average response (individual)", 
                xlabel=L"f", 
                ylabel = LaTeXString(label)
            )
            scatterlines!(ax_freq_ampl_indv, 
                freq_range, 
                amplitude[idx_dof, :], 
                color = :black, 
                markersize = 5
            )
            save(figure_dir * "$static_case\\freq_ampl" * fig_freetext, fig_freq_ampl_indv, px_per_unit = 4)

            label = "\$ $(dof_names[dof]) \$"
            fig_frequency_sweep = Figure(size = (1280, 720))
            ax_frequency_sweep = Axis(fig_frequency_sweep[1, 1], 
                #title=L"$(dof_names[dof]) vs time", 
                xlabel=L"\tau", 
                ylabel = LaTeXString(label)
            )
            tf₋₁ = 0.0
            for (idx_freq, freq) in enumerate(freq_range)
                tf = response_data["tf_f_$(round(freq, sigdigits=4))"]
                n_steps = response_data["n_steps_f_$(round(freq, sigdigits=4))"]
                idx_plot = vcat(collect(1:10:n_steps-1), n_steps) 
                u::Matrix{Float32} = response_data["u_f_$(round(freq, sigdigits=4))"]
                t = range(0, tf, length = n_steps)
                lines!(ax_frequency_sweep, 
                    t[idx_plot] .+ tf₋₁, u[idx_dof, :], 
                    color=COLORS[idx_freq]
                )
                tf₋₁ += tf
            end
            # display(GLMakie.Screen(), fig_frequency_sweep)
            save(figure_dir * "$static_case\\frequency_sweep" * fig_freetext, fig_frequency_sweep, px_per_unit = 4)
            
            label = "\\text{Eq. pos. }\$$idx_file\$"
            scatterlines!(ax_freq_avg_response[idx_dof], 
                freq_range, 
                avg_response[idx_dof, :], 
                color = COLORS[idx_file], 
                markersize = 5,
                label = LaTeXString(label)
                )
            scatterlines!(ax_freq_ampl[idx_dof], 
                freq_range, 
                amplitude[idx_dof, :], 
                color = COLORS[idx_file], 
                markersize = 5,
                label = LaTeXString(label)
            )
            
        end
        scatter!(ax_static_pos,
            T1_static[idx_static], 
            T2_static[idx_static], 
            markersize = 25, 
            color = COLORS[idx_file], 
            marker = :circle
        )

        fig_contact_cycles_1 = Figure()
        ax_contact_cycles_1 = Axis(fig_contact_cycles_1[1, 1], 
            #title="Contact cycles at left contact", 
            xlabel=L"x_{t1}", 
            ylabel=L"T_1"
        )    
        fig_contact_cycles_2 = Figure()
        ax_contact_cycles_2 = Axis(fig_contact_cycles_2[1, 1], 
            #title="Contact cycles at right contact", 
            xlabel=L"x_{t2}", 
            ylabel=L"T_2"
        )

        xₜ₁_min = 0.0
        xₜ₂_min = 0.0
        xₜ₁_max = -1.0
        xₜ₂_max = -1.0
        for (idx_freq, freq) in enumerate(freq_plots)
            label = "\$f = $(round(freq, sigdigits=3))\$" 
            
            fig_T1 = Figure(size = (1280, 720))
            ax_T1 = Axis(fig_T1[1, 1],
                #title="T₁ vs time", 
                xlabel=L"\tau", 
                ylabel=L"T_1"
            )
            fig_T2 = Figure(size = (1280, 720))
            ax_T2 = Axis(fig_T2[1, 1],
                #title="T₂ vs time", 
                xlabel=L"\tau", 
                ylabel=L"T_2"
            )
            t_cycles = 6*2π/freq # Plot the last 6 cycles
            tf = response_data["tf_f_$(round(freq, sigdigits=4))"]
            n_steps = response_data["n_steps_f_$(round(freq, sigdigits=4))"]
            n_final = response_data["n_final_f_$(round(freq, sigdigits=4))"]
            idx_plot_cntc = vcat(collect(n_final:10:n_steps-2), n_steps-1)
            u::Matrix{Float32} = response_data["u_f_$(round(freq, sigdigits=4))"]
            T1::Vector{Float32} = response_data["T1_f_$(round(freq, sigdigits=4))"]
            T2::Vector{Float32} = response_data["T2_f_$(round(freq, sigdigits=4))"]
            N1::Vector{Float32} = response_data["N1_f_$(round(freq, sigdigits=4))"]
            N2::Vector{Float32} = response_data["N2_f_$(round(freq, sigdigits=4))"]
            t = range(0, tf, length = n_steps)
            t_cntc = t[idx_plot_cntc]
            idx_cycles = findmin(abs.(t_cntc .- (tf - t_cycles)))[2]
            μ = pm.friction_law_1.μ
            lines!(ax_T1, 
                t_cntc[idx_cycles:end], 
                T1[idx_cycles:end], 
                color=COLORS[idx_freq]
            )
            lines!(ax_T1, 
                t_cntc[idx_cycles:end], 
                μ*N1[idx_cycles:end], 
                color = COLORS[idx_freq],
                label = L"\pm \mu N_1",
                linestyle = :dash,
            )
            lines!(ax_T1, 
                t_cntc[idx_cycles:end], 
                -μ*N1[idx_cycles:end], 
                color = COLORS[idx_freq],
                linestyle = :dash,
            )
            axislegend(ax_T1, position = :rb)
            save(figure_dir * "$static_case\\T1_f_$(replace(string(round(freq, sigdigits=4)), "." => "-")).png", fig_T1, px_per_unit = 4)
            lines!(ax_T2, 
                t_cntc[idx_cycles:end], 
                T2[idx_cycles:end], 
                color=COLORS[idx_freq]
            )
            lines!(ax_T2, 
                t_cntc[idx_cycles:end], 
                μ*N2[idx_cycles:end], 
                color = COLORS[idx_freq],
                label = L"\pm \mu N_2",
                linestyle = :dash,
            )
            lines!(ax_T2, 
                t_cntc[idx_cycles:end], 
                -μ*N2[idx_cycles:end], 
                color = COLORS[idx_freq],
                linestyle = :dash,
            )
            axislegend(ax_T2, position = :rb)
            save(figure_dir * "$static_case\\T2_f_$(replace(string(round(freq, sigdigits=4)), "." => "-")).png", fig_T2, px_per_unit = 4)
            xₜ₁ = response_data["x_t1_f_$(round(freq, sigdigits=4))"]
            T₁ = response_data["T1_f_$(round(freq, sigdigits=4))"]
            lines!(ax_contact_cycles_1, 
                xₜ₁, T₁, 
                label = LaTeXString(label), 
                color = COLORS[idx_freq]
            )
            minimum(xₜ₁) < xₜ₁_min && (xₜ₁_min = minimum(xₜ₁))
            maximum(xₜ₁) > xₜ₁_max && (xₜ₁_max = maximum(xₜ₁))
            xₜ₂ = response_data["x_t2_f_$(round(freq, sigdigits=4))"]
            T₂ = response_data["T2_f_$(round(freq, sigdigits=4))"]
            lines!(ax_contact_cycles_2, 
                xₜ₂, T₂, 
                label = LaTeXString(label),
                color = COLORS[idx_freq]
            )
            minimum(xₜ₂) < xₜ₂_min && (xₜ₂_min = minimum(xₜ₂))
            maximum(xₜ₂) > xₜ₂_max && (xₜ₂_max = maximum(xₜ₂))
        end
        xlims!(ax_contact_cycles_1, 
            xₜ₁_min - 0.1*(xₜ₁_max - xₜ₁_min), 
            xₜ₁_max + 0.1*(xₜ₁_max - xₜ₁_min), 
        )
        xlims!(ax_contact_cycles_2, 
            xₜ₂_min - 0.1*(abs(xₜ₂_max - xₜ₂_min)), 
            xₜ₂_max + 0.1*(abs(xₜ₂_max - xₜ₂_min))
        )
        axislegend(ax_contact_cycles_1, position = :rb)
        axislegend(ax_contact_cycles_2, position = :rb)
        save(figure_dir * "$static_case\\contact_cycles_left.png", fig_contact_cycles_1, px_per_unit = 4)
        save(figure_dir * "$static_case\\contact_cycles_right.png", fig_contact_cycles_2, px_per_unit = 4)    
    end

    position_avg_response = [:rt, :rt, :rb, :rt, :rb]
    for (idx_dof, dof) in enumerate(dofs)
        fig_freetext = ""
        if dof != 1
            fig_freetext *= "_dof_$dof"
        end
        fig_freetext *= ".png"
    
        axislegend(ax_freq_avg_response[idx_dof], position = position_avg_response[idx_dof])
        axislegend(ax_freq_ampl[idx_dof])

        # display(GLMakie.Screen(), fig_freq_ampl[idx_dof])
        save(figure_dir * "freq_ampl" * fig_freetext, fig_freq_ampl[idx_dof], px_per_unit = 4)

        # display(GLMakie.Screen(), fig_freq_avg_response[idx_dof])
        save(figure_dir * "freq_avg_response" * fig_freetext, fig_freq_avg_response[idx_dof], px_per_unit = 4)
    end
    save(figure_dir * "static_pos.png", fig_static_pos, px_per_unit = 4)
    
end
