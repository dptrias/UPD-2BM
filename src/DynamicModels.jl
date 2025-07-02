include("Common.jl")

## TIME INTEGRATORS
# Runge-Kutta 2
function rk2(f, u, p, t, Δt)
    k1 = f(u, p, t, true)
    k2 = f(u + 0.5Δt * k1, p, t + 0.5Δt)
    return u + Δt * k2
end

# Runge-Kutta 4
function rk4(f, u, p, t, Δt)
    k1 = f(u, p, t, true)
    k2 = f(u + 0.5Δt * k1, p, t + 0.5Δt)
    k3 = f(u + 0.5Δt * k2, p, t + 0.5Δt)
    k4 = f(u + Δt * k3, p, t + Δt)
    return u + (Δt / 6) * (k1 + 2k2 + 2k3 + k4)
end

## DYNAMIC MODELS
# Two degrees of freedom model (ψ, y_D)
function two_dof_model(u, p, t, store = false)
    @unpack feq, ζ_ψ, α, fₓ, θ, b, Mₑ, FC, T_F, friction_law = p

    x, ẋ = u[1:2], u[3:4]

    # Recurrent trigonometric functions
    sin_θ = sin(θ)
    cos_θ = cos(θ)
    sin_ψ = sin(x[1])
    cos_ψ = cos(x[1])

    # Relative position of the contact
    xₜ = (sin_ψ - b*(1 - cos_ψ))*cos_θ + (-(1 - cos_ψ) - b*sin_ψ - x[2])*sin_θ
    xₙ = (sin_ψ - b*(1 - cos_ψ))*sin_θ - (-(1 - cos_ψ) - b*sin_ψ - x[2])*cos_θ

    # Contact forces
    T, N = g(xₜ, xₙ, friction_law)

    # Storing of contact forces and relative positions
    if store && p.index !== nothing
        p.T[p.index[]] = T      # Because of the way the storing is done, for an integration with N steps,
        p.N[p.index[]] = N      # the contact variables will be stored at indices 1:N-1, so the last step will not be stored.
        p.xₜ[p.index[]] = xₜ
        p.xₙ[p.index[]] = xₙ
        p.index[] += 1
    end

    # External forces
    M, _, FC = external_forces(t, feq, Mₑ, FC, T_F)

    # Equations of motion
    ẍ = zeros(2)
    ẍ[1] = -2ζ_ψ*ẋ[1] - x[1] + M + T*( -cos(θ + x[1]) + b*sin(θ + x[1]) ) + N*( -sin(θ + x[1]) - b*cos(θ + x[1]) )
    ẍ[2] = -fₓ^2*x[2] + α*( FC + 2T*sin_θ - 2N*cos_θ )

    return [ẋ; ẍ]
end

# Five degrees of freedom model (ψ₁, ψ₂, x_D, y_D, ψ_D)
function five_dof_model(u, p, t, store = false)
    @unpack feq, ζ_ψ, α, β, fₓ, f_ψ, θ, b, dₙ, dₜ, Mₑ, FC, T_F, friction_law_1, friction_law_2 = p

    x, ẋ = u[1:5], u[6:end]

    # Relative displacement of the contact points
    xₜ₁, xₙ₁, xₜ₂, xₙ₂ = contact_points_displacements(x, p)

    # Contact forces
    T₁, N₁ = g(xₜ₁, xₙ₁, friction_law_1)
    T₂, N₂ = g(xₜ₂, xₙ₂, friction_law_2)

    # Storing of contact forces and relative positions
    if store && p.index !== nothing
        p.T₁[p.index[]] = T₁
        p.N₁[p.index[]] = N₁
        p.T₂[p.index[]] = T₂
        p.N₂[p.index[]] = N₂
        p.xₜ₁[p.index[]] = xₜ₁
        p.xₙ₁[p.index[]] = xₙ₁
        p.xₜ₂[p.index[]] = xₜ₂
        p.xₙ₂[p.index[]] = xₙ₂
        p.index[] += 1
    end

    # External forces
    M₁, M₂, FC = external_forces(t, feq, Mₑ, FC, T_F)

    # Equations of motion
    ẍ = zeros(5)
    ẍ[1] = -2ζ_ψ*ẋ[1] - x[1] + M₁ + T₁*( -cos(θ - x[5] + x[1]) + b*sin(θ - x[5] + x[1]) ) + N₁*( -sin(θ - x[5] + x[1]) - b*cos(θ - x[5] + x[1]) )
    ẍ[2] = -2ζ_ψ*ẋ[2] - x[2] + M₂ + T₂*( cos(θ + x[5] - x[2]) - b*sin(θ + x[5] - x[2]) ) + N₂*( sin(θ + x[5] - x[2]) + b*cos(θ + x[5] - x[2]) )
    ẍ[3] = -fₓ^2*x[3] + α*( T₁*cos(θ - x[5]) + N₁*sin(θ - x[5]) - T₂*cos(θ + x[5]) - N₂*sin(θ + x[5]) )
    ẍ[4] = -fₓ^2*x[4] + α*( FC + T₁*sin(θ - x[5]) - N₁*cos(θ - x[5]) + T₂*sin(θ + x[5]) - N₂*cos(θ + x[5]) )
    ẍ[5] = -f_ψ^2*x[5] + β*( T₁*(p.dₙ - xₙ₁) - N₁*(p.dₜ - xₜ₁) - T₂*(p.dₙ - xₙ₂) + N₂*(p.dₜ - xₜ₂) )
    
    return [ẋ; ẍ]
end

## AUXILIARY FUNCTIONS
# External forces
function external_forces(t, feq::Union{Float64, Nothing}, Mₑ, FC, T_F::Union{Float64, Nothing})
    if feq !== nothing
        M₁ = Mₑ * cos(feq*t)
    else 
        M₁ = 0.0
    end
    M₂ = -M₁
    if T_F !== nothing
        FC = FC * tanh(t/T_F)
    end
    return [M₁, M₂, FC]
end