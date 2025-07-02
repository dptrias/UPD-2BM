using JLD2, GLMakie, Statistics
include("..\\..\\src\\Common.jl")

path = abspath(@__DIR__, "..\\..")
begin
    file_name = "Eigenvalues_comparison_stiffness"
    mkpath("julia\\figures\\Eigenvalues_5DOF_comparison")
    # Load files
    file_1 = "Modal_analysis_5DOF"
    data_1 = load("julia\\data\\eigenvalues\\" * file_1 * ".jld2")
    label_1 = "Reference"
    λ₁ = data_1["λ"]
    file_2 = "Modal_analysis_5DOF_more_stiff_cntc"
    data_2 = load("julia\\data\\eigenvalues\\" * file_2 * ".jld2")
    label_2 = "4 × Stiffness"
    λ₂ = data_2["λ"]

    fig_λ = Figure(size = (600, 600))
    ax_λ = Axis(fig_λ[1, 1], title="Eigenvalues", xlabel="Re", ylabel="Im"; 
                xgridvisible = true, ygridvisible = true, xgridcolor = :gray, ygridcolor = :gray, aspect = 1)
    scatter!(ax_λ, real(λ₁), imag(λ₁); markersize = 16, color = :blue, label = label_1)
    scatter!(ax_λ, real(λ₂), imag(λ₂); markersize = 16, color = :red, label = label_2)
    scatter!(ax_λ, real(λ₁), 2*imag(λ₁); markersize = 8, label = "2 × Reference", marker = :utriangle, strokewidth = 2, strokecolor = :blue, color = :transparent)
    axislegend(position = :lt)
    # # Eigenvalues zoom
    # x₁ = -8e-4
    # x₂ = 1e-4
    # y₁ = 32.88
    # y₂ = 32.885
    # ax_inset = Axis(fig_λ[1, 1],
    #                 width=Relative(0.2), height=Relative(0.2), halign=0.5, valign=0.9, title="Zoomed View"; 
    #                 xgridvisible = true, ygridvisible = true, xgridcolor = :gray, ygridcolor = :gray)
    #                     xlims!(ax_inset, x₁, x₂)
    # xtick_values = [x₁, 0]
    # ax_inset.xticks = (xtick_values, string.(xtick_values))
    # ylims!(ax_inset, y₁, y₂)
    # scatter!(ax_inset, real(λ₁[9:10]), imag(λ₁[9:10]); markersize = 12, color = :blue, label = label_1)
    # scatter!(ax_inset, real(λ₂[9:10]), imag(λ₂[9:10]); markersize = 12, color = :red, label = label_2)
    # translate!(ax_inset.blockscene, 0, 0, 150)
    # x₁ = -5e-2
    # x₂ = 5e-2
    # y₁ = 32
    # y₂ = 34
    # border_rect = Rect2(x₁, y₁, x₂-x₁, y₂-y₁)
    # lines!(ax_λ, border_rect, color=:black, linewidth=1)
    # start_point = Point2f(x₁, y₁+1)
    # direction = Vec2f(-0.82, -8)
    # arrows!(ax_λ, [start_point], [direction]; arrowsize = 10)
    
    # Save figure
    save(abspath(path, "..", "figures\\Eigenvalues_5DOF_comparison", file_name * ".png"), fig_λ)

    # for i=1:5
    #     idx = 2*i-1
    #     println("Im(λ_$i / λ_$(i)_ref) = ", abs(imag(λ₂[idx]))/abs(imag(λ₁[idx])))
    # end
end