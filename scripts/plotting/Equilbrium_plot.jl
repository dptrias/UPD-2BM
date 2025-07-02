using GLMakie, JLD2, Julianim, Printf, MathTeXEngine
include("..\\..\\src\\Common.jl")

path = abspath(@__DIR__, "..\\..")
function find_repetition_bounds(vec::Vector{T}) where T
    if isempty(vec)
        return Int[], Int[]
    end

    starts = [1]
    ends = Int[]

    for i in 2:length(vec)
        if vec[i] != vec[i-1]
            push!(ends, i - 1)     # previous group ended at i-1
            push!(starts, i)       # new group starts at i
        end
    end

    push!(ends, length(vec))       # close the last group
    return starts, ends
end

function split_by_pivots(vec::Vector{T}, pivots::Vector{T}) where T
    pivots_sorted = sort(pivots)
    segments = Vector{Vector{T}}()

    for i in 1:length(pivots_sorted)+1
        if i == 1
            idx_start = 1
            idx_end = findlast(==(pivots_sorted[1]), vec)
        elseif i == length(pivots_sorted)+1
            idx_start = findfirst(==(pivots_sorted[end]), vec)
            idx_end = length(vec)
        else
            idx_start = findfirst(==(pivots_sorted[i-1]), vec)
            idx_end = findlast(==(pivots_sorted[i]), vec)
        end
        push!(segments, vec[idx_start:idx_end])
    end

    return Tuple(segments)
end

begin
    file = "Equilibrium_5DOF_final"
    data = load("julia\\data\\static\\" * file * ".jld2")
    solution = data["solutions"]
    pm = data["pm"]
    μ = pm.friction_law_1.μ

    T₁ = [sol.T1 for sol in solution]
    T₂ = [sol.T2 for sol in solution]
    N₁ = [sol.N1 for sol in solution]
    N₂ = [sol.N2 for sol in solution]
    ψ₁ = [sol.psi1 for sol in solution]
    ψ₂ = [sol.psi2 for sol in solution]
    x_D = [sol.x_D for sol in solution]
    y_D = [sol.y_D for sol in solution]
    ψ_D = [sol.psi_D for sol in solution]
end

begin
    path = abspath(@__DIR__,"..\\..\\","tex\\Images\\" * file * "\\")
    mkpath(path)

    ΔT₁ = (maximum(T₁) - minimum(T₁))/20
    ΔT₂ = (maximum(T₂) - minimum(T₂))/20

    set_publication_theme!()
    update_theme!(
    fonts = Attributes(
            :bold => texfont(:bold),
            :bolditalic => texfont(:bolditalic),
            :italic => texfont(:italic),
            :regular => texfont(:regular)
        )
    )

    # ψ₁
    fig_psi1 = Figure()
    ax_psi1 = Axis(fig_psi1[1, 1], 
        xlabel = L"T_1~(\times 10^{-4})", 
        ylabel = L"T_2~(\times 10^{-4})", 
        aspect = 1,
        xtickformat = x -> string.(round.(x .* 1e4; digits = 2)),
        ytickformat = y -> string.(round.(y .* 1e4; digits = 2))
    )
    limits!(ax_psi1, minimum(T₁)-ΔT₁, maximum(T₁)+ΔT₁, minimum(T₂)-ΔT₂, maximum(T₂)+ΔT₂)

    co = contourf!(ax_psi1, T₁, T₂, ψ₁)
    Colorbar(fig_psi1[1, 2], 
        co;
        label = L"\psi_1 \text{ [rad]}"
    )
    colgap!(fig_psi1.layout, 30)

    # idx = argmin(abs.(T₁ .- (-5.655e-4)))
    idx = 91790
    lines!(ax_psi1, T₁[idx:end], T₁[idx:end], color = RED, linestyle=:dash)

    display(GLMakie.Screen(), fig_psi1)
    save(path * "\\psi_1.png", fig_psi1, px_per_unit = 4)
    
    # ψ₂
    fig_psi2 = Figure()
    ax_psi2 = Axis(fig_psi2[1, 1], 
        xlabel = L"T_1~(\times 10^{-4})", 
        ylabel = L"T_2~(\times 10^{-4})", 
        aspect = 1,
        xtickformat = x -> string.(round.(x .* 1e4; digits = 2)),
        ytickformat = y -> string.(round.(y .* 1e4; digits = 2))
    )
    limits!(ax_psi2, minimum(T₁)-ΔT₁, maximum(T₁)+ΔT₁, minimum(T₂)-ΔT₂, maximum(T₂)+ΔT₂)
    
    co = contourf!(ax_psi2, T₁, T₂, ψ₂)
    Colorbar(fig_psi2[1, 2], 
        co;
        label = L"\psi_2 \text{ [rad]}"
    )
    colgap!(fig_psi2.layout, 30)
    
    display(GLMakie.Screen(), fig_psi2)
    save(path * "\\psi_2.png", fig_psi2, px_per_unit = 4)
    
    
    # x_D
    fig_x_D = Figure()
    ax_x_D = Axis(fig_x_D[1, 1], 
        xlabel = L"T_1~(\times 10^{-4})", 
        ylabel = L"T_2~(\times 10^{-4})", 
        aspect = 1,
        xtickformat = x -> string.(round.(x .* 1e4; digits = 2)),
        ytickformat = y -> string.(round.(y .* 1e4; digits = 2))
    )
    limits!(ax_x_D, minimum(T₁)-ΔT₁, maximum(T₁)+ΔT₁, minimum(T₂)-ΔT₂, maximum(T₂)+ΔT₂)
    
    co = contourf!(ax_x_D, T₁, T₂, x_D)
    Colorbar(fig_x_D[1, 2], 
        co;
        label = L"x_D"
    )
    colgap!(fig_x_D.layout, 30)
    
    display(GLMakie.Screen(), fig_x_D)
    save(path * "\\x_D.png", fig_x_D, px_per_unit = 4)
    
    
    # y_D
    fig_y_D = Figure()
    ax_y_D = Axis(fig_y_D[1, 1], 
        xlabel = L"T_1~(\times 10^{-4})", 
        ylabel = L"T_2~(\times 10^{-4})", 
        aspect = 1,
        xtickformat = x -> string.(round.(x .* 1e4; digits = 2)),
        ytickformat = y -> string.(round.(y .* 1e4; digits = 2))
    )
    limits!(ax_y_D, minimum(T₁)-ΔT₁, maximum(T₁)+ΔT₁, minimum(T₂)-ΔT₂, maximum(T₂)+ΔT₂)
    
    co = contourf!(ax_y_D, T₁, T₂, y_D)
    Colorbar(fig_y_D[1, 2], 
        co;
        label = L"y_D"
    )
    colgap!(fig_y_D.layout, 30)
    
    display(GLMakie.Screen(), fig_y_D)
    save(path * "\\y_D.png", fig_y_D, px_per_unit = 4)
    
    
    # ψ_D
    fig_psi_D = Figure()
    ax_psi_D = Axis(fig_psi_D[1, 1], 
        xlabel = L"T_1~(\times 10^{-4})", 
        ylabel = L"T_2~(\times 10^{-4})", 
        aspect = 1,
        xtickformat = x -> string.(round.(x .* 1e4; digits = 2)),
        ytickformat = y -> string.(round.(y .* 1e4; digits = 2))
    )
    limits!(ax_psi_D, minimum(T₁)-ΔT₁, maximum(T₁)+ΔT₁, minimum(T₂)-ΔT₂, maximum(T₂)+ΔT₂)
    
    co = contourf!(ax_psi_D, T₁, T₂, ψ_D)
    Colorbar(fig_psi_D[1, 2], 
        co;
        label = L"\psi_D \text{ [rad]}"
    )
    colgap!(fig_psi_D.layout, 30)
    
    display(GLMakie.Screen(), fig_psi_D)
    save(path * "\\psi_D.png", fig_psi_D, px_per_unit = 4)

    # N₁
    fig_N_1 = Figure()
    ax_N_1 = Axis(fig_N_1[1, 1], 
        xlabel = L"T_1~(\times 10^{-4})", 
        ylabel = L"T_2~(\times 10^{-4})", 
        aspect = 1,
        xtickformat = x -> string.(round.(x .* 1e4; digits = 2)),
        ytickformat = y -> string.(round.(y .* 1e4; digits = 2))
    )
    limits!(ax_N_1, minimum(T₁)-ΔT₁, maximum(T₁)+ΔT₁, minimum(T₂)-ΔT₂, maximum(T₂)+ΔT₂)


    co = contourf!(ax_N_1, T₁, T₂, N₁)
    Colorbar(fig_N_1[1, 2], 
    co;
    label = L"N_1"
    )
    colgap!(fig_N_1.layout, 30)

    idx_start, idx_2_T₂ = find_repetition_bounds(T₁)
    idx_1_T₁, idx_1_T₂, idx_2_T₁ = split_by_pivots(idx_start, [91790, 1821102])
    lines!(ax_N_1, -μ*N₁[idx_1_T₁], T₂[idx_1_T₁], color = RED, linestyle=:solid, label = L"\pm \mu N_1", linewidth = 4)
    lines!(ax_N_1, μ*N₁[idx_2_T₁], T₂[idx_2_T₁], color = RED, linestyle=:solid, linewidth = 4)
    axislegend(ax_N_1, position = (0.05, 0.95))

    display(GLMakie.Screen(), fig_N_1)
    save(path * "\\N_1.png", fig_N_1, px_per_unit = 4)
    
    # N₂
    fig_N_2 = Figure()
    ax_N_2 = Axis(fig_N_2[1, 1], 
        xlabel = L"T_1~(\times 10^{-4})", 
        ylabel = L"T_2~(\times 10^{-4})", 
        aspect = 1,
        xtickformat = x -> string.(round.(x .* 1e4; digits = 2)),
        ytickformat = y -> string.(round.(y .* 1e4; digits = 2))
    )
    limits!(ax_N_2, minimum(T₁)-ΔT₁, maximum(T₁)+ΔT₁, minimum(T₂)-ΔT₂, maximum(T₂)+ΔT₂)

    co = contourf!(ax_N_2, T₁, T₂, N₂)
    Colorbar(fig_N_2[1, 2], 
        co;
        label = L"N_2"
    )
    colgap!(fig_N_2.layout, 30)
    
    lines!(ax_N_2, T₁[idx_1_T₂], -μ*N₂[idx_1_T₂], color = RED, linestyle=:solid, label = L"\pm \mu N_2", linewidth = 4)
    lines!(ax_N_2, T₁[idx_2_T₂], μ*N₂[idx_2_T₂], color = RED, linestyle=:solid, linewidth = 4)
    axislegend(ax_N_2, position = (0.05, 0.95))


    display(GLMakie.Screen(), fig_N_2)
    save(path * "\\N_2.png", fig_N_2, px_per_unit = 4)

end