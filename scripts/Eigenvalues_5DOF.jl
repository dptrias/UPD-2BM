using QuadEig, LinearAlgebra, JLD2, Julianim, GLMakie, MathTeXEngine, Printf
include("..\\src\\StaticModels.jl")

path = abspath(@__DIR__, "..")
## INITIALIZATION
begin
    ## Get data from static equilibrium
    freetext_static = "_final_reduced" # Freetext for static equilibrium file
    static_data = jldopen(joinpath(path, "data\\static\\Equilibrium_5DOF" * freetext_static * ".jld2"), "r")
    sol_static = static_data["solutions"]
    pm_static = static_data["pm"]
    close(static_data)

    # Extract the data from static solution
    T1_static = [sol.T1 for sol in sol_static]
    T2_static = [sol.T2 for sol in sol_static]
    n_static = length(T1_static)
    u_static = zeros(5, n_static)
    u_static[1, :] = [sol.psi1 for sol in sol_static]
    u_static[2, :] = [sol.psi2 for sol in sol_static]
    u_static[3, :] = [sol.x_D for sol in sol_static]
    u_static[4, :] = [sol.y_D for sol in sol_static]
    u_static[5, :] = [sol.psi_D for sol in sol_static]

    # Create path for figures
    freetext = "_final" # Freetext for outputs
    mkpath(joinpath(path, "figures\\Eigenvalues_5DOF" * freetext))
    mkpath(joinpath(path, "data\\eigenvalues"))
    set_publication_theme!()
    update_theme!(
    fonts = Attributes(
            :bold => texfont(:bold),
            :bolditalic => texfont(:bolditalic),
            :italic => texfont(:italic),
            :regular => texfont(:regular)
        )
    )

    # Static equilibrium positions to be analyzed
    static_indices = [7 91791 533652 734134 1024500]
    idx_for_plot = 4 # Index of the previous vector 
    #                  to be used for the modal analysis
    interactive = false 
end

# Function to normalize the eigenvectors
function normalize_φ(φ)
    for col in eachcol(φ)
        normalize!(col)
        col ./= maximum(abs.(real(col)))
        # col ./= sum(col)
    end
    return φ
end

## MAIN
# Computation of the eigenvalues and eigenvectors for the static equilibrium solutions
begin
    # Mass & Damping matrices
    M = [1 0 0             0             0;
         0 1 0             0             0;
         0 0 1/pm_static.α 0             0;
         0 0 0             1/pm_static.α 0;
         0 0 0             0             1/pm_static.β]
    C = [2pm_static.ζ_ψ 0              0 0 0;
         0              2pm_static.ζ_ψ 0 0 0;
         0              0              0 0 0;
         0              0              0 0 0;
         0              0              0 0 0]

    # Initialize eigenvalues and eigenvectors matrices
    n_analysis = length(static_indices)
    λ::Matrix{ComplexF64} = zeros(10, n_analysis)
    φ::Array{ComplexF64, 3} = zeros(5, 10, n_analysis)
    for (i, idx_static) in enumerate(static_indices)
        local K = stiffness_5DOF(u_static[:, idx_static], pm_static)

        # Quadratic pencil
        local l = QuadEig.linearize(K, C, M)
        local d = QuadEig.deflate(l)
        local sol_eigen = eigen(d.A, d.B)
        λ[:, i] = sol_eigen.values
        φ[:, :, i] = normalize_φ((l.V*sol_eigen.vectors)[1:5,:])
    end
end

function format_sci_latex(x::Float64)
    s = @sprintf("%.3g", x)  # e.g. "1.23e-04"
    if occursin('e', s)
        base, exp = split(s, 'e')
        if exp == "00"
            return base
        end
        sign = startswith(exp, "-") ? "-" : "+"
        expval = replace(exp, r"^[+-]" => "")
        return "$(base)\\times 10^{$sign$expval}"
    else
        return s
    end
end

function format_complex_latex(z::Complex)
    re_str = format_sci_latex(real(z))
    im_str = format_sci_latex(abs(imag(z)))
    sign_str = imag(z) ≥ 0 ? "\\pm" : "\\mp"
    return "$re_str $sign_str \\mathrm{i}$im_str"
end

begin
    f_eigenvalues = joinpath(path, "figures\\Eigenvalues_5DOF" * freetext * "\\eigenvalues.txt")
    f_eigenvectors = joinpath(path, "figures\\Eigenvalues_5DOF" * freetext * "\\eigenvectors.txt")

    open(f_eigenvalues, "w") do file
        for i in 1:2:10
            println(file, format_complex_latex(λ[i, idx_for_plot]) * " \\\\")
        end
    end

    open(f_eigenvectors, "w") do file
        for j in 1:2:10
            println(file, "phi_$(round(Int, (j+1)/2))")
            for i in 1:5
                println(file, format_complex_latex(φ[i, j, idx_for_plot]) * " \\\\")
            end
            println(file, "")
        end
    end
end

# Eigenvalues at chosen static equilibrium position,
# these will be considered as the reference eigenvalues
begin 
    fig_eigenvalues_zoom = Figure(size=(1250, 500))
    ax_eigenvalues_zoom = Axis(fig_eigenvalues_zoom[1, 1], 
        xlabel=L"\lambda_R", 
        ylabel=L"\lambda_I"; 
        aspect = AxisAspect(2.5)
    )
    scatter!(ax_eigenvalues_zoom, 
        real(λ[1:2:end, idx_for_plot]), 
        abs.(imag(λ[1:2:end, idx_for_plot])),
        marker = :cross,
        color = BLUE
    )
    xlims!(ax_eigenvalues_zoom, 1.1*minimum(real(λ[:, idx_for_plot])) , 0.1*abs(minimum(real(λ[:, idx_for_plot]))))
    ylims!(ax_eigenvalues_zoom, -0.1*maximum(imag(λ[:, idx_for_plot])) , 1.1*maximum(imag(λ[:, idx_for_plot])))
    textlabel!(ax_eigenvalues_zoom,
        Point2f(-1e-4, 15),
        text = L"3", padding = 0,
        shape = Circle(Point2f(0), 0.9f0),
        shape_limits = Rect2f(-sqrt(0.45), -sqrt(0.45), sqrt(1.8), sqrt(1.8)),
        keep_aspect = true
    )
    
    ax_zoom_12 = Axis(fig_eigenvalues_zoom[1, 1],
        width=Relative(0.15),
        height=Relative(0.5), 
        halign=0.2, valign=0.5,
        xlabelvisible = false,
        ylabelvisible = false,
        # xticklabelrotation = π/10,
        xaxisposition = :top,
        yaxisposition = :right,
    )
    scatter!(ax_zoom_12, 
        real(λ[[1, 3], idx_for_plot]), 
        abs.(imag(λ[[1, 3], idx_for_plot])),
        marker = :cross,
        color = BLUE
    )
    x₁ = minimum(real(λ[[1, 3], idx_for_plot]))
    x₂ = maximum(real(λ[[1, 3], idx_for_plot]))
    x₁ -= 0.01 * abs(x₁)
    x₂ += 0.01 * abs(x₂)
    ax_zoom_12.xticks = round.([x₁, x₂], sigdigits=3)
    y₁ = minimum(abs.(imag(λ[[1, 3], idx_for_plot])))
    y₂ = maximum(abs.(imag(λ[[1, 3], idx_for_plot])))
    y₁ -= 3 * abs(y₁)
    y₂ += 3 * abs(y₂)
    xlims!(ax_zoom_12, x₁, x₂)
    ylims!(ax_zoom_12, y₁, y₂)
    textlabel!(ax_zoom_12,
        Point2f(x₁ + 0.28 * (x₂ - x₁), y₁ + 0.62 * (y₂ - y₁)),
        text = L"1", padding = 0,
        shape = Circle(Point2f(0), 0.9f0),
        shape_limits = Rect2f(-sqrt(0.45), -sqrt(0.45), sqrt(1.8), sqrt(1.8)),
        keep_aspect = true
    )
    textlabel!(ax_zoom_12,
        Point2f(x₁ + 0.65 * (x₂ - x₁), y₁ + 0.25 * (y₂ - y₁)),
        text = L"2", padding = 0,
        shape = Circle(Point2f(0), 0.9f0),
        shape_limits = Rect2f(-sqrt(0.45), -sqrt(0.45), sqrt(1.8), sqrt(1.8)),
        keep_aspect = true
    )
    x₁ -= 0.005 * abs(x₁)
    x₂ += 0.005 * abs(x₂)
    y₁ -= 0.15 * abs(y₁)
    y₂ += 0.15 * abs(y₂)
    lines!(ax_eigenvalues_zoom, 
        Rect2(x₁, y₁, x₂-x₁, y₂-y₁), 
        color=:black, 
        linewidth=2
    )
    start_point = Point2f(x₁+0.0005, y₂)
    direction = Vec2f(0.0005, 9)
    arrows!(ax_eigenvalues_zoom, 
        [start_point], 
        [direction]; 
        arrowsize = 15
    )
    translate!(ax_zoom_12.blockscene, 0, 0, 150)
    
    ax_zoom_45 = Axis(fig_eigenvalues_zoom[1, 1],
        width=Relative(0.2), 
        height=Relative(0.4), 
        halign=0.7, 
        valign=0.5
    )
    ax_zoom_45.xlabelvisible = false
    ax_zoom_45.ylabelvisible = false
    scatter!(ax_zoom_45, 
        real(λ[[7, 9], idx_for_plot]), 
        abs.(imag(λ[[7, 9], idx_for_plot])),
        marker = :cross,
        color = BLUE
    )
    x₁ = minimum(real(λ[[7, 9], idx_for_plot]))
    x₂ = 0.0
    x₁ -= 0.1 * abs(x₁)
    ax_zoom_45.xticks = round.([x₁, x₂], sigdigits=4)
    y₁ = minimum(abs.(imag(λ[[7, 9], idx_for_plot])))
    y₂ = maximum(abs.(imag(λ[[7, 9], idx_for_plot])))
    y₁ -= 0.01 * abs(y₁)
    y₂ += 0.01 * abs(y₂)
    xlims!(ax_zoom_45, x₁, x₂)
    ylims!(ax_zoom_45, y₁, y₂)
    textlabel!(ax_zoom_45,
        Point2f(x₁ + 0.28 * (x₂ - x₁), y₁ + 0.62 * (y₂ - y₁)),
        text = L"4", padding = 0,
        shape = Circle(Point2f(0), 0.9f0),
        shape_limits = Rect2f(-sqrt(0.45), -sqrt(0.45), sqrt(1.8), sqrt(1.8)),
        keep_aspect = true
    )
    textlabel!(ax_zoom_45,
        Point2f(x₁ + 0.65 * (x₂ - x₁), y₁ + 0.25 * (y₂ - y₁)),
        text = L"5", padding = 0,
        shape = Circle(Point2f(0), 0.9f0),
        shape_limits = Rect2f(-sqrt(0.45), -sqrt(0.45), sqrt(1.8), sqrt(1.8)),
        keep_aspect = true
    )
    x₁ -= 200 * abs(x₁)
    x₂ += abs(x₁)
    y₁ -= 0.03 * abs(y₁)
    y₂ += 0.03 * abs(y₂)
    lines!(ax_eigenvalues_zoom, 
        Rect2(x₁, y₁, x₂-x₁, y₂-y₁), 
        color=:black, 
        linewidth=2
    )
    start_point = Point2f(x₁, y₁)
    direction = Vec2f(-0.0015, -16.5)
    arrows!(ax_eigenvalues_zoom, 
        [start_point], 
        [direction]; 
        arrowsize = 15
    )
    translate!(ax_zoom_45.blockscene, 0, 0, 150)

    display(GLMakie.Screen(), fig_eigenvalues_zoom)
    # Save figure
    save(joinpath(path, "figures\\Eigenvalues_5DOF" * freetext * "\\eigenvalues.png"), fig_eigenvalues_zoom, px_per_unit = 10)

    local f_path = joinpath(path, "data\\eigenvalues\\Modal_analysis_5DOF" * freetext * ".jld2")
    if interactive
        print("Save modal analysis results? [Y/n]: ")
        response = readline() |> strip
    else 
        response = "Y" # Automatically save results
    end
    if response == "Y" || response == "y"
        if isfile(f_path)
            println("[warning] File $f_path already exists. Choose a different name or delete the file")
        else
            @save f_path λ φ pm_static
            println("[info] Results saved to $f_path")
        end
    else
        println("[info] Results not saved")
    end
end

# Plot of eigenvalues for a series of static equilibrium solution 
begin
    fig_static_solution = Figure()
    ax_static_solution = Axis(fig_static_solution[1, 1], 
        #title="Static equilibrium positions", 
        xlabel=L"T_1", 
        ylabel=L"T_2", 
        aspect = 1
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

    # Eigenvalues plot
    fig_eig_family = Figure()
    ax_eig_family = Axis(fig_eig_family[1, 1], 
        # title="Eigenvalues", 
        xlabel=L"\lambda_R", 
        ylabel=L"\lambda_I";
        aspect = 1
    )
    xlims!(ax_eig_family, 1.1*minimum(real(λ[:, :])) , 0.1*abs(minimum(real(λ[:, :]))))
    ylims!(ax_eig_family, -0.1*maximum(imag(λ[:, :])) , 1.1*maximum(imag(λ[:, :])))
    # Eigenvalues zoom
    x₁ = -3.5e-7
    x₂ = -3e-7
    y₁ = 79.6751
    y₂ = 79.67525
    ax_inset = Axis(fig_eig_family[1, 1],
        width=Relative(0.2), 
        height=Relative(0.2), 
        halign=0.5, 
        valign=0.7, 
        title=L"\text{Mode }5"; 
    )
    ax_inset.xlabelvisible = false
    ax_inset.ylabelvisible = false
    xtick_values = [x₁, x₂]
    ax_inset.xticks = (xtick_values, string.(xtick_values))
    ytick_values = [y₁, round((y₁+y₂)/2, sigdigits=3), y₂]
    ax_inset.yticks = (ytick_values, string.(ytick_values))
    xlims!(ax_inset, x₁, x₂)
    ylims!(ax_inset, y₁, y₂)
    for (i, idx_static) in enumerate(static_indices)
        # Plot initial condition
        scatter!(ax_static_solution, 
            T1_static[idx_static], 
            T2_static[idx_static]; 
            color = COLORS[i], 
            marker = :circle,
            markersize = 25,
        )

        # Plot eigenvalues
        scatter!(ax_eig_family, 
            real(λ[1:2:10, i]), 
            abs.(imag(λ[1:2:10, i])); 
            marker = :cross, 
            color = COLORS[i]
        )
        scatter!(ax_inset, 
            real(λ[10, i]), 
            abs.(imag(λ[10, i])); 
            marker = :cross, 
            color = COLORS[i]
        )
    end
    x₁ -= 1e3 * abs(x₁)
    x₂ = abs(x₁)
    y₂ += 0.03 * abs(y₂)
    y₁ -= 0.03 * abs(y₁)
    lines!(ax_eig_family, 
        Rect2(x₁, y₁, x₂-x₁, y₂-y₁), 
        color=:black, 
        linewidth=2
    )
    # Add the arrow to the plot
    start_point = Point2f(x₁, y₁+1)
    direction = Vec2f(-0.0048, -17)
    arrows!(ax_eig_family, 
        [start_point], 
        [direction]; 
        arrowsize = 15
    )
    translate!(ax_inset.blockscene, 0, 0, 150)
    
    display(GLMakie.Screen(), fig_static_solution)
    display(GLMakie.Screen(), fig_eig_family)
    # Save figures
    save(joinpath(path, "figures\\Eigenvalues_5DOF" * freetext * "\\eigenvalues_family.png"), fig_eig_family, px_per_unit = 4)
    save(joinpath(path, "figures\\Eigenvalues_5DOF" * freetext * "\\eigenvalues_family_eq_pos.png"), fig_static_solution, px_per_unit = 4)

    f_path = joinpath(path, "data\\eigenvalues\\Eigenvalues_family_5DOF" * freetext * ".jld2")
    if interactive
        print("Save eigenvalues family results? [Y/n]: ")
        response = readline() |> strip
    else 
        response = "Y" # Automatically save results
    end
    if response == "Y" || response == "y"
        if isfile(f_path)
            println("[warning] File $f_path already exists. Choose a different name or delete the file")
        else
            @save f_path λ φ sol_static pm_static
            println("[info] Results saved to $f_path")
        end
    else
        println("[info] Results not saved")
    end
end