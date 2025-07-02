using JLD2, Julianim, Printf
include("Spring.jl")

path = abspath(@__DIR__, "..\\..")
# Video parameters
begin
    FPS = 40
    time = 10
    frames = FPS * time
    periods = 2
    ts = range(0, stop=2π * periods, length=frames)
    dt = ts[2] - ts[1]
    freetext = "_final"

    t = Observable(0.0) # Observable to update automatically and everything that depends on this observable
    step!(t) = t[] += dt
end

# Get modal analysis data
begin
    modal_data = jldopen(abspath(path, "..", "data\\eigenvalues\\Modal_analysis_5DOF_final.jld2"), "r")
    λ = modal_data["λ"]
    φ = modal_data["φ"]
    close(modal_data)
end

function mode_variables(t::Observable{Float64}, λ, φ, i)
    # Calculate the mode shape for the i-th mode
    x = Observable{Float64}[]
    for j=1:5
        push!(x, @lift real(exp(λ[i]*$t) * φ[j, i]))
    end
    return x
end

modes = 1:5
for mode in modes
    t[] = 0.0
    ηᵢ = mode_variables(t, λ, φ, mode*2 - 1)

    ψ₁ = @lift $(ηᵢ[1]) / 10
    ψ₂ = @lift $(ηᵢ[2]) / 10
    x_D = @lift $(ηᵢ[3]) / 10
    y_D = @lift $(ηᵢ[4]) / 10
    ψ_D = @lift $(ηᵢ[5]) / 10

    a = 0.4 # Horizontal distance between blade and damper
    z = -0.025 # Vertical distance between blade and damper
    M_ψ_1 = @lift [cos($ψ₁) sin($ψ₁); -sin($ψ₁) cos($ψ₁)]
    M_ψ_2 = @lift [cos($ψ₂) sin($ψ₂); -sin($ψ₂) cos($ψ₂)]
    M_ψ_D = @lift [cos($ψ_D) sin($ψ_D); -sin($ψ_D) cos($ψ_D)]
    b = 0.05
    f = 0.1
    c = 0.6

    # Left blade
    #        | T_1           ^
    #        |               |
    #        |   b           |
    #        | <--->         |
    #        |               | c
    # H -----+----- B        |
    #        |         ^     |
    #        |         | f   |
    #        | C1      v     v
    x_C1, y_C1 = -a, z
    x_B = @lift x_C1 + ($M_ψ_1*[b,f])[1]
    y_B = @lift y_C1 + ($M_ψ_1*[b,f])[2]
    x_H = @lift x_C1 + ($M_ψ_1*[-b,f])[1]
    y_H = @lift y_C1 + ($M_ψ_1*[-b,f])[2]
    x_T1 = @lift x_C1 + ($M_ψ_1*[0,c])[1]
    y_T1 = @lift y_C1 + ($M_ψ_1*[0,c])[2]

    # Right blade
    #        | T_2           ^
    #        |               |
    #        |   b           |
    #        | <--->         |
    #        |               | c
    # A -----+----- C        |
    #        |         ^     |
    #        |         | f   |
    #        | C2      v     v
    x_C2, y_C2 = a, z
    x_A = @lift x_C2 + ($M_ψ_2*[-b,f])[1]
    y_A = @lift y_C2 + ($M_ψ_2*[-b,f])[2]
    x_C = @lift x_C2 + ($M_ψ_2*[b,f])[1]
    y_C = @lift y_C2 + ($M_ψ_2*[b,f])[2]
    x_T2 = @lift x_C2 + ($M_ψ_2*[0,c])[1]
    y_T2 = @lift y_C2 + ($M_ψ_2*[0,c])[2]

    # Damper
    #        D3    
    #        /\        ^
    #       /  \       | d
    #      /    \      |
    #     /   +  \     v
    #    /    G   \
    # D1 ---------- D2
    #           g
    #         <--->
    dₙ = 0.02
    dₜ = -0.001
    θ = pi/4
    d = dₙ/cos(θ)
    g = 2d/tan(θ)
    x_D1 = @lift $x_D + ($M_ψ_D*[-g, -d])[1]
    y_D1 = @lift $y_D + ($M_ψ_D*[-g, -d])[2]
    x_D2 = @lift $x_D + ($M_ψ_D*[g, -d])[1]
    y_D2 = @lift $y_D + ($M_ψ_D*[g, -d])[2]
    x_D3 = @lift $x_D + ($M_ψ_D*[0, d])[1]
    y_D3 = @lift $y_D + ($M_ψ_D*[0, d])[2]
    damper_coords = @lift Point2f[($x_D1, $y_D1), ($x_D2, $y_D2), ($x_D3, $y_D3)]

    # Damper contact points
    G_θ_D = [cos(θ) sin(θ); -sin(θ) cos(θ)]
    G_θ_I = [sin(θ) cos(θ); -cos(θ) sin(θ)]
    x_F = @lift $x_D + ($M_ψ_D*G_θ_D*[dₜ, dₙ])[1]
    y_F = @lift $y_D + ($M_ψ_D*G_θ_D*[dₜ, dₙ])[2]
    x_E = @lift $x_D - ($M_ψ_D*G_θ_I*[dₙ, dₜ])[1]
    y_E = @lift $y_D - ($M_ψ_D*G_θ_I*[dₙ, dₜ])[2]

    S = CoiledSpring(0.02, 4)
    s_coords_1 = [@lift spring_coordinates($x_B, $y_B, $x_E, $y_E, S)[i] for i in 1:2]
    s_coords_2 = [@lift spring_coordinates($x_F, $y_F, $x_A, $y_A, S)[i] for i in 1:2]

    set_publication_theme!(fullscreen=true)

    fig = Figure(size=(600, 400))
    ax = Axis(fig[1, 1])

    hidedecorations!(ax)
    hidespines!(ax)

    xlims!(ax, -0.7, 0.7)
    ylims!(ax, -0.3, 0.7)

    textlabel!(ax, Point2f(0, 0.5), text = @lift(@sprintf("t= %.2f s", $t)), fontsize = 20, text_align = (:center, :center))
    
    # Left blade
    lines!(ax, @lift([$x_H, $x_B]), @lift([$y_H, $y_B]); color=:black, linewidth=3)
    lines!(ax, @lift([x_C1, $x_T1]), @lift([y_C1, $y_T1]); color=:black, linewidth=3)

    # Right blade
    lines!(ax, @lift([$x_A, $x_C]), @lift([$y_A, $y_C]); color=:black, linewidth=3)
    lines!(ax, @lift([x_C2, $x_T2]), @lift([y_C2, $y_T2]); color=:black, linewidth=3)

    poly!(ax, damper_coords, color=:black, strokecolor=:black, strokewidth=3)

    lines!(ax, s_coords_1..., color=GREYS[2], linewidth=2)
    lines!(ax, s_coords_2..., color=GREYS[2], linewidth=2)

    format = "gif"
    local dirpath = (abspath(path,"..","figures\\Modal_analysis_5DOF" * freetext))
    if !isdir(dirpath)
        mkpath(dirpath)
    end
    local file = joinpath(dirpath, "Mode_$mode" * "." * format)
    record(fig, file, ts; framerate=FPS) do i
        step!(t)
    end
end