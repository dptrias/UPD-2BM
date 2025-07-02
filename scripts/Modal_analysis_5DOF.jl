using LinearAlgebra, JLD2, GLMakie, Julianim, MathTeXEngine
include("..\\src\\Common.jl")

path = abspath(@__DIR__, "..")
# Get eigenvalue data
begin
    modal_data = jldopen(joinpath(path, "data\\eigenvalues\\Eigenvalues_family_5DOF_final_even_lower_def.jld2"), "r")
    save_figs = true
    display_figs = false
    freetext = "_final"
    
    global φ = modal_data["φ"]
    global sol_static = modal_data["sol_static"]
    close(modal_data)
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
    global T1_static = [sol.T1 for sol in sol_static]
    global T2_static = [sol.T2 for sol in sol_static]
    global n_static = length(T1_static)
    global u_static = zeros(5, n_static)
    u_static[1, :] = [sol.psi1 for sol in sol_static]
    u_static[2, :] = [sol.psi2 for sol in sol_static]
    u_static[3, :] = [sol.x_D for sol in sol_static]
    u_static[4, :] = [sol.y_D for sol in sol_static]
    u_static[5, :] = [sol.psi_D for sol in sol_static]
end

begin
    η = zeros(5, n_static) # Vector of relative modal coordinates
    for i in 1:n_static
        η[:, i] = abs.(real(φ[:, [1, 3, 5, 7, 9], i]) \ u_static[:, i])
        η[:, i] ./= sum(η[:, i])
    end
    
    max_η = maximum(η)

    fig_mode_1 = Figure()
    ax_mode_1 = Axis(fig_mode_1[1, 1], 
        # title = "η₁ vs T₁ & T₂", 
        xlabel = L"T_1", 
        ylabel = L"T_2", 
        aspect = 1
    )
    co = contourf!(ax_mode_1, 
        T1_static, 
        T2_static, 
        η[1, :], 
        levels = range(0, max_η, 100), 
        colormap = :viridis
    )
    Colorbar(fig_mode_1[1, 2], 
        co,
        label = L"η_1",
    )

    fig_mode_2 = Figure()
    ax_mode_2 = Axis(fig_mode_2[1, 1], 
        # title = "η₂ vs T₁ & T₂", 
        xlabel = L"T_1", 
        ylabel = L"T_2", 
        aspect = 1
    )
    co = contourf!(ax_mode_2, 
        T1_static, 
        T2_static, 
        η[2, :], 
        levels = range(0, max_η, 100), 
        colormap = :viridis
    )
    Colorbar(fig_mode_2[1, 2], 
        co,
        label = L"η_2",
    )

    fig_mode_3 = Figure()
    ax_mode_3 = Axis(fig_mode_3[1, 1], 
        # title = "η₃ vs T₁ & T₂", 
        xlabel = L"T_1", 
        ylabel = L"T_2", 
        aspect = 1
    )
    co = contourf!(ax_mode_3, 
        T1_static, 
        T2_static, 
        η[3, :], 
        levels = range(0, max_η, 100), 
        colormap = :viridis
    )
    Colorbar(fig_mode_3[1, 2], 
        co,
        label = L"η_3",
    )

    fig_mode_4 = Figure()
    ax_mode_4 = Axis(fig_mode_4[1, 1], 
        # title = "η₄ vs T₁ & T₂", 
        xlabel = L"T_1", 
        ylabel = L"T_2", 
        aspect = 1
    )
    co = contourf!(ax_mode_4, 
        T1_static, 
        T2_static, 
        η[4, :], 
        levels = range(0, max_η, 100), 
        colormap = :viridis
    )
    Colorbar(fig_mode_4[1, 2], 
        co,
        label = L"η_4",
    )

    fig_mode_5 = Figure()
    ax_mode_5 = Axis(fig_mode_5[1, 1], 
        # title = "η₅ vs T₁ & T₂", 
        xlabel = L"T_1", 
        ylabel = L"T_2", 
        aspect = 1
    )
    co = contourf!(ax_mode_5, 
        T1_static, 
        T2_static, 
        η[5, :], 
        levels = range(0, max_η, 100), 
        colormap = :viridis
    )
    Colorbar(fig_mode_5[1, 2], 
        co,
        label = L"η_5",
    )

    if display_figs
        display(GLMakie.Screen(), fig_mode_1)
        display(GLMakie.Screen(), fig_mode_2)
        display(GLMakie.Screen(), fig_mode_3)
        display(GLMakie.Screen(), fig_mode_4)
        display(GLMakie.Screen(), fig_mode_5)
    end

    if save_figs
        mkpath(joinpath(path, "figures\\Eigenvalues_5DOF" * freetext))
        save(joinpath(path, "figures\\Eigenvalues_5DOF" * freetext * "\\mode_1.png"), fig_mode_1)
        save(joinpath(path, "figures\\Eigenvalues_5DOF" * freetext * "\\mode_2.png"), fig_mode_2)
        save(joinpath(path, "figures\\Eigenvalues_5DOF" * freetext * "\\mode_3.png"), fig_mode_3)
        save(joinpath(path, "figures\\Eigenvalues_5DOF" * freetext * "\\mode_4.png"), fig_mode_4)
        save(joinpath(path, "figures\\Eigenvalues_5DOF" * freetext * "\\mode_5.png"), fig_mode_5)
    end

end