using QuadEig, LinearAlgebra, JLD2, Colors, GLMakie
include("..\\src\\Common.jl")

path = abspath(@__DIR__, "..")
begin
    ## Get data from static equilibrium
    freetext = ""
    static_data = jldopen(joinpath(path, "data\\static\\Equilibrium_2DOF" * freetext *".jld2"), "r")
    sol_static = static_data["solutions"]
    pm_static = static_data["pm"]
    close(static_data)

    # Extract the data 
    psi_static = [sol.psi for sol in sol_static]
    y_D_static = [sol.y_D for sol in sol_static]
    N_static = [sol.N for sol in sol_static]
    T_static = [sol.T for sol in sol_static]
    xt_static = [sol.xt for sol in sol_static]
    xn_static = [sol.xn for sol in sol_static]
    wt_static = [sol.wt for sol in sol_static]
end


function stiffness_2DOF(ψ₀, y_D₀, θ, α, fₓ, b, kₙ, kₜ)
    K = zeros(2, 2)
    sin_0 = sin(θ)
    cos_0 = cos(θ)
    sin_ψ = sin(ψ₀)
    cos_ψ = cos(ψ₀)

    K[1,1] = 1 - kₜ*((cos_ψ - b*sin_ψ)*cos_0 + (-sin_ψ - b*cos_ψ)*sin_0)*(-cos(ψ₀ + θ) + b*sin(ψ₀ + θ)) - kₜ*((sin_ψ - b*(1 - cos_ψ))*cos_0 + 
             (-1 + cos_ψ - b*sin_ψ - y_D₀)*sin_0)*(sin(ψ₀ + θ) + b*cos(ψ₀ + θ)) - 
             kₙ*((cos_ψ - b*sin_ψ)*sin_0 - (-sin_ψ - b*cos_ψ)*cos_0)*(-sin(ψ₀ + θ) - b*cos(ψ₀ + θ)) - 
             kₙ*((sin_ψ - b*(1 - cos_ψ))*sin_0 - (-1 + cos_ψ - b*sin_ψ - y_D₀)*cos_0)*(-cos(ψ₀ + θ) + b*sin(ψ₀ + θ))
    K[1,2] = kₜ*sin_0*(-cos(ψ₀ + θ) + b*sin(ψ₀ + θ)) - kₙ*cos_0*(-sin(ψ₀ + θ) - b*cos(ψ₀ + θ))
    K[1,2] = kₜ*sin_0*(-cos(ψ₀ + θ) + b*sin(ψ₀ + θ)) - kₙ*cos_0*(-sin(ψ₀ + θ) - b*cos(ψ₀ + θ))
    K[2,1] = -2*kₜ*((cos_ψ - b*sin_ψ)*cos_0 + (-sin_ψ - b*cos_ψ)*sin_0)*sin_0 + 
             2*kₙ*((cos_ψ - b*sin_ψ)*sin_0 - (-sin_ψ - b*cos_ψ)*cos_0)*cos_0
    K[2,1] = -kₜ*((cos_ψ - b*sin(ψ₀))*cos_0 + (-sin(ψ₀) - b*cos_ψ)*sin_0)*sin_0 + 
             kₙ*((cos_ψ - b*sin(ψ₀))*sin_0 - (-sin(ψ₀) - b*cos_ψ)*cos_0)*cos_0
    K[2,2] = fₓ^2/(2α) + 2*kₜ*sin_0^2 + 2*kₙ*cos_0^2

    return K
end

begin
    # Mass & Damping matrices
    M = [1 0; 
         0 1/pm_static.α]
    ζ_ψ = 5e-2 / (2*sqrt(4.0e3 * 6.348e-4)) # Damping ratio
    C = [2ζ_ψ 0; # [2*pm_static.ζ_ψ 0;
         0    0]

    n_static::Int  = length(T_static)
    distance::Int = 100 # Distance between the static equilibrium solutions
    n_iter::Int = floor(Int, ((n_static - 1)  / distance) + 1)
    eigenvalues::Matrix{ComplexF64} = zeros(4, n_iter + 1)

    # Initial condition plot
    fig_initial_cond = Figure(size = (600, 600))
    ax_initial_cond = Axis(fig_initial_cond[1, 1], title="ψ vs T", xlabel="T", ylabel="ψ")
    lines!(ax_initial_cond, T_static, psi_static, linewidth = 2)
    # Eigenvalues plot
    fig_eigenvalues = Figure(size = (600, 600))
    ax_eigenvalues = Axis(fig_eigenvalues[1, 1], title="Eigenvalues", xlabel="Re", ylabel="Im"; xgridvisible = true, ygridvisible = true, xgridcolor = :gray, ygridcolor = :gray, aspect = 1)
    # Colors for plotting
    cols = distinguishable_colors(n_iter+1, [RGB(1,1,1), RGB(0,0,0)], dropseed=true)
    pcols = map(col -> (red(col), green(col), blue(col)), cols)
    
    for i in 1:n_iter+1
        index_static::Int = round(Int, (i - 1) * distance + 1)
        if index_static > n_static
            index_static = n_static
        end
        local K = stiffness_2DOF(psi_static[index_static], y_D_static[index_static], pm_static.θ, pm_static.α, pm_static.fₓ, pm_static.b, pm_static.kₙ, pm_static.kₜ)

        # Quadratic pencil
        local l = QuadEig.linearize(K, C, M)
        local d = QuadEig.deflate(l)
        eigenvalues[:,i] = QuadEig.eigvals(d.A, d.B)

        # Plot initial condition
        scatter!(ax_initial_cond, T_static[index_static], psi_static[index_static]; markersize = 15, color = pcols[i])
        # Plot eigenvalues
        scatter!(ax_eigenvalues, real(eigenvalues[1,i]), imag(eigenvalues[1,i]); markersize = 10, color = pcols[i])
        scatter!(ax_eigenvalues, real(eigenvalues[2,i]), imag(eigenvalues[2,i]); markersize = 10, color = pcols[i])
        scatter!(ax_eigenvalues, real(eigenvalues[3,i]), imag(eigenvalues[3,i]); markersize = 10, color = pcols[i])
        scatter!(ax_eigenvalues, real(eigenvalues[4,i]), imag(eigenvalues[4,i]); markersize = 10, color = pcols[i])
    end
    save(joinpath(path, "figures\\Stability_2DOF" * freetext * "\\eigenvalues.png"), fig_eigenvalues)
    save(joinpath(path, "figures\\Stability_2DOF" * freetext * "\\eigenvalues_eq_pos.png"), fig_initial_cond)
end