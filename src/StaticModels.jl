using NonlinearSolve
include("Common.jl")

## EQUILIBRIUM POSITIONS
# 5 DOF EQUILIBRIUM
function equilibrium_pos_5DOF(p::TBM_5DOF_Parameters, T₁::Float64, T₂::Float64, u₀::Vector{Float64})::Union{Vector{Float64}, Nothing}
    f(u, p) = five_dof_equilibrium(u, p, T₁, T₂)
    problem = NonlinearProblem(f, u₀, p)
    result = solve(problem, TrustRegion(), reltol = 1e-12, abstol = 1e-12)
    
    if SciMLBase.successful_retcode(result)
        return result.u
    else
        return nothing
    end
end

function five_dof_equilibrium(x, p::TBM_5DOF_Parameters, T₁::Float64, T₂::Float64)
    @unpack θ, b, dₙ, dₜ, fₓ, f_ψ, α, β, FC, friction_law_1 = p
    kₙ = friction_law_1.kₙ

    # Relative displacement of the contact points
    xₜ₁, xₙ₁, xₜ₂, xₙ₂ = contact_points_displacements(x, p)

    # Normal forces
    N₁ = kₙ*xₙ₁
    N₂ = kₙ*xₙ₂

    F = Vector{typeof(x[1])}(undef, 5) # Initialize the vector for the static equations
    # Static equations
    F[1] = x[1] - T₁*( -cos(θ - x[5] + x[1]) + b*sin(θ - x[5] + x[1]) ) - N₁*( -sin(θ - x[5] + x[1]) - b*cos(θ - x[5] + x[1]) )
    F[2] = x[2] - T₂*( cos(θ + x[5] - x[2]) - b*sin(θ + x[5] - x[2]) ) - N₂*( sin(θ + x[5] - x[2]) + b*cos(θ + x[5] - x[2]) )
    F[3] = fₓ^2*x[3]/α - ( T₁*cos(θ - x[5]) + N₁*sin(θ - x[5]) - T₂*cos(θ + x[5]) - N₂*sin(θ + x[5]) )
    F[4] = fₓ^2*x[4]/α - (FC + T₁*sin(θ - x[5]) - N₁*cos(θ - x[5]) + T₂*sin(θ + x[5]) - N₂*cos(θ + x[5]) )
    F[5] = f_ψ^2*x[5]/β - ( T₁*(dₙ - xₙ₁) - N₁*(dₜ - xₜ₁) - T₂*(dₙ - xₙ₂) + N₂*(dₜ - xₜ₂) )
    return F
end

# 2 DOF EQUILIBRIUM
function equilibrium_pos_2DOF(p::TBM_2DOF_Parameters, T::Float64)::Union{Vector{Float64}, Nothing}

    f(u, p) = two_dof_equilibrium!(u, p, T)
    problem = NonlinearProblem(f, u₀, p)
    result = solve(problem, TrustRegion(), reltol = 1e-12, abstol = 1e-12)
    
    if SciMLBase.successful_retcode(result)
        return result.u
    else
        return nothing
    end
end

function two_dof_equilibrium!(x, p, T)
    @unpack θ, b, fₓ, α, FC, friction_law = p
    kₙ = friction_law.kₙ
    
    sin_θ = sin(θ)
    cos_θ = cos(θ)
    sin_ψ = sin(x[1])
    cos_ψ = cos(x[1])

    F = Vector{typeof(x[1])}(undef, 2) # Initialize the vector for the static equations
    xₙ = (sin_ψ - b*(1 - cos_ψ))*sin_θ - (-(1 - cos_ψ) - b*sin_ψ - x[2])*cos_θ
    F[1] = x[1] - T*( -cos(θ + x[1]) + b*sin(θ + x[1]) ) - kₙ*xₙ*( -sin(θ + x[1]) - b*cos(θ + x[1]) )
    F[2] = fₓ^2/α*x[2] - FC - 2T*sin_θ + 2kₙ*xₙ*cos_θ
    return F
end

## JACOBIANS (LINEARIZED STIFFNESS MATRICES)
function stiffness_5DOF(u₀, p)
    @unpack θ, b, dₙ, dₜ, fₓ, f_ψ, α, β, friction_law_1 = p
    kₙ = friction_law_1.kₙ
    kₜ = friction_law_1.kₜ

    ψ₁₀, ψ₂₀, x_D₀, y_D₀, ψ_D₀ = u₀[1], u₀[2], u₀[3], u₀[4], u₀[5]
    K = zeros(5, 5)
    sin_θ = sin(θ) 
    cos_θ = cos(θ)
    sin_ψ1 = sin(ψ₁₀)
    cos_ψ1 = cos(ψ₁₀)
    sin_ψ2 = sin(ψ₂₀)
    cos_ψ2 = cos(ψ₂₀)
    sin_ψD_p_θ = sin(ψ_D₀ + θ)
    cos_ψD_p_θ = cos(ψ_D₀ + θ)
    sin_ψD_m_θ = sin(ψ_D₀ - θ)
    cos_ψD_m_θ = cos(ψ_D₀ - θ)
    
    K[1,1] = 1 - kₜ*((cos_ψ1 -b*sin_ψ1)*cos_ψD_m_θ - (-sin_ψ1 -b*cos_ψ1)*sin_ψD_m_θ)*(-cos(θ - ψ_D₀ + ψ₁₀) + b*sin(θ - ψ_D₀ + ψ₁₀)) - kₜ*((sin_ψ1 -b*(1 - cos_ψ1) - x_D₀ + dₜ*(cos_ψD_m_θ - cos_θ) + dₙ*(-sin_ψD_m_θ - sin_θ))*cos_ψD_m_θ - (-1 + cos_ψ1 -b*sin_ψ1 - y_D₀ - dₙ*(cos_ψD_m_θ - cos_θ) + dₜ*(-sin_ψD_m_θ - sin_θ))*sin_ψD_m_θ)*(sin(θ - ψ_D₀ + ψ₁₀) + b*cos(θ - ψ_D₀ + ψ₁₀)) - kₙ*(-(cos_ψ1 -b*sin_ψ1)*sin_ψD_m_θ - (-sin_ψ1 -b*cos_ψ1)*cos_ψD_m_θ)*(-sin(θ - ψ_D₀ + ψ₁₀) -b*cos(θ - ψ_D₀ + ψ₁₀)) - kₙ*(-(sin_ψ1 -b*(1 - cos_ψ1) - x_D₀ + dₜ*(cos_ψD_m_θ - cos_θ) + dₙ*(-sin_ψD_m_θ - sin_θ))*sin_ψD_m_θ - (-1 + cos_ψ1 -b*sin_ψ1 - y_D₀ - dₙ*(cos_ψD_m_θ - cos_θ) + dₜ*(-sin_ψD_m_θ - sin_θ))*cos_ψD_m_θ)*(-cos(θ - ψ_D₀ + ψ₁₀) + b*sin(θ - ψ_D₀ + ψ₁₀))
    K[1,2] = 0
    K[1,3] = kₜ*cos_ψD_m_θ*(-cos(θ - ψ_D₀ + ψ₁₀) + b*sin(θ - ψ_D₀ + ψ₁₀)) - kₙ*sin_ψD_m_θ*(-sin(θ - ψ_D₀ + ψ₁₀) - b*cos(θ - ψ_D₀ + ψ₁₀))
    K[1,4] = -kₜ*sin_ψD_m_θ*(-cos(θ - ψ_D₀ + ψ₁₀) + b*sin(θ - ψ_D₀ + ψ₁₀)) - kₙ*cos_ψD_m_θ*(-sin(θ - ψ_D₀ + ψ₁₀) - b*cos(θ - ψ_D₀ + ψ₁₀))
    K[1,5] = -kₜ*((-dₜ*sin_ψD_m_θ - dₙ*cos_ψD_m_θ)*cos_ψD_m_θ - (sin_ψ1 - b*(1 - cos_ψ1) - x_D₀ + dₜ*(cos_ψD_m_θ - cos_θ) + dₙ*(-sin_ψD_m_θ - sin_θ))*sin_ψD_m_θ - (dₙ*sin_ψD_m_θ - dₜ*cos_ψD_m_θ)*sin_ψD_m_θ - (-1 + cos_ψ1 - b*sin_ψ1 - y_D₀ - dₙ*(cos_ψD_m_θ - cos_θ) + dₜ*(-sin_ψD_m_θ - sin_θ))*cos_ψD_m_θ)*(-cos(θ - ψ_D₀ + ψ₁₀) + b*sin(θ - ψ_D₀ + ψ₁₀)) - kₜ*((sin_ψ1 - b*(1 - cos_ψ1) - x_D₀ + dₜ*(cos_ψD_m_θ - cos_θ) + dₙ*(-sin_ψD_m_θ - sin_θ))*cos_ψD_m_θ - (-1 + cos_ψ1 - b*sin_ψ1 - y_D₀ - dₙ*(cos_ψD_m_θ - cos_θ) + dₜ*(-sin_ψD_m_θ - sin_θ))*sin_ψD_m_θ)*(-sin(θ - ψ_D₀ + ψ₁₀) - b*cos(θ - ψ_D₀ + ψ₁₀)) - kₙ*(-(-dₜ*sin_ψD_m_θ - dₙ*cos_ψD_m_θ)*sin_ψD_m_θ - (sin_ψ1 - b*(1 - cos_ψ1) - x_D₀ + dₜ*(cos_ψD_m_θ - cos_θ) + dₙ*(-sin_ψD_m_θ - sin_θ))*cos_ψD_m_θ - (dₙ*sin_ψD_m_θ - dₜ*cos_ψD_m_θ)*cos_ψD_m_θ + (-1 + cos_ψ1 - b*sin_ψ1 - y_D₀ - dₙ*(cos_ψD_m_θ - cos_θ) + dₜ*(-sin_ψD_m_θ - sin_θ))*sin_ψD_m_θ)*(-sin(θ - ψ_D₀ + ψ₁₀) - b*cos(θ - ψ_D₀ + ψ₁₀)) - kₙ*(-(sin_ψ1 - b*(1 - cos_ψ1) - x_D₀ + dₜ*(cos_ψD_m_θ - cos_θ) + dₙ*(-sin_ψD_m_θ - sin_θ))*sin_ψD_m_θ - (-1 + cos_ψ1 - b*sin_ψ1 - y_D₀ - dₙ*(cos_ψD_m_θ - cos_θ) + dₜ*(-sin_ψD_m_θ - sin_θ))*cos_ψD_m_θ)*(cos(θ - ψ_D₀ + ψ₁₀) - b*sin(θ - ψ_D₀ + ψ₁₀))
    K[2,1] = 0 
    K[2,2] = 1 - kₜ*(-(cos_ψ2 + b*sin_ψ2)*cos_ψD_p_θ + (-sin_ψ2 + b*cos_ψ2)*sin_ψD_p_θ)*(cos(-θ - ψ_D₀ + ψ₂₀) + b*sin(-θ - ψ_D₀ + ψ₂₀)) - kₜ*(-(sin_ψ2 + b*(1 - cos_ψ2) - x_D₀ - dₜ*(cos_ψD_p_θ - cos_θ) - dₙ*(sin_ψD_p_θ - sin_θ))*cos_ψD_p_θ + (-1 + cos_ψ2 + b*sin_ψ2 - y_D₀ - dₙ*(cos_ψD_p_θ - cos_θ) + dₜ*(sin_ψD_p_θ - sin_θ))*sin_ψD_p_θ)*(-sin(-θ - ψ_D₀ + ψ₂₀) + b*cos(-θ - ψ_D₀ + ψ₂₀)) - kₙ*(-(cos_ψ2 + b*sin_ψ2)*sin_ψD_p_θ - (-sin_ψ2 + b*cos_ψ2)*cos_ψD_p_θ)*(-sin(-θ - ψ_D₀ + ψ₂₀) + b*cos(-θ - ψ_D₀ + ψ₂₀)) - kₙ*(-(sin_ψ2 + b*(1 - cos_ψ2) - x_D₀ - dₜ*(cos_ψD_p_θ - cos_θ) - dₙ*(sin_ψD_p_θ - sin_θ))*sin_ψD_p_θ - (-1 + cos_ψ2 + b*sin_ψ2 - y_D₀ - dₙ*(cos_ψD_p_θ - cos_θ) + dₜ*(sin_ψD_p_θ - sin_θ))*cos_ψD_p_θ)*(-cos(-θ - ψ_D₀ + ψ₂₀) -b*sin(-θ - ψ_D₀ + ψ₂₀))
    K[2,3] = -kₜ*cos_ψD_p_θ*(cos(-θ - ψ_D₀ + ψ₂₀) + b*sin(-θ - ψ_D₀ + ψ₂₀)) - kₙ*sin_ψD_p_θ*(-sin(-θ - ψ_D₀ + ψ₂₀) + b*cos(-θ - ψ_D₀ + ψ₂₀))
    K[2,4] = kₜ*sin_ψD_p_θ*(cos(-θ - ψ_D₀ + ψ₂₀) + b*sin(-θ - ψ_D₀ + ψ₂₀)) - kₙ*cos_ψD_p_θ*(-sin(-θ - ψ_D₀ + ψ₂₀) + b*cos(-θ - ψ_D₀ + ψ₂₀))
    K[2,5] = -kₜ*(-(dₜ*sin_ψD_p_θ - dₙ*cos_ψD_p_θ)*cos_ψD_p_θ + (sin_ψ2 + b*(1 - cos_ψ2) - x_D₀ - dₜ*(cos_ψD_p_θ - cos_θ) - dₙ*(sin_ψD_p_θ - sin_θ))*sin_ψD_p_θ + (dₙ*sin_ψD_p_θ + dₜ*cos_ψD_p_θ)*sin_ψD_p_θ + (-1 + cos_ψ2 + b*sin_ψ2 - y_D₀ - dₙ*(cos_ψD_p_θ - cos_θ) + dₜ*(sin_ψD_p_θ - sin_θ))*cos_ψD_p_θ)*(cos(-θ - ψ_D₀ + ψ₂₀) + b*sin(-θ - ψ_D₀ + ψ₂₀)) - kₜ*(-(sin_ψ2 + b*(1 - cos_ψ2) - x_D₀ - dₜ*(cos_ψD_p_θ - cos_θ) - dₙ*(sin_ψD_p_θ - sin_θ))*cos_ψD_p_θ + (-1 + cos_ψ2 + b*sin_ψ2 - y_D₀ - dₙ*(cos_ψD_p_θ - cos_θ) + dₜ*(sin_ψD_p_θ - sin_θ))*sin_ψD_p_θ)*(sin(-θ - ψ_D₀ + ψ₂₀) - b*cos(-θ - ψ_D₀ + ψ₂₀)) - kₙ*(-(dₜ*sin_ψD_p_θ - dₙ*cos_ψD_p_θ)*sin_ψD_p_θ - (sin_ψ2 + b*(1 - cos_ψ2) - x_D₀ - dₜ*(cos_ψD_p_θ - cos_θ) - dₙ*(sin_ψD_p_θ - sin_θ))*cos_ψD_p_θ - (dₙ*sin_ψD_p_θ + dₜ*cos_ψD_p_θ)*cos_ψD_p_θ + (-1 + cos_ψ2 + b*sin_ψ2 - y_D₀ - dₙ*(cos_ψD_p_θ - cos_θ) + dₜ*(sin_ψD_p_θ - sin_θ))*sin_ψD_p_θ)*(-sin(-θ - ψ_D₀ + ψ₂₀) + b*cos(-θ - ψ_D₀ + ψ₂₀)) - kₙ*(-(sin_ψ2 + b*(1 - cos_ψ2) - x_D₀ - dₜ*(cos_ψD_p_θ - cos_θ) - dₙ*(sin_ψD_p_θ - sin_θ))*sin_ψD_p_θ - (-1 + cos_ψ2 + b*sin_ψ2 - y_D₀ - dₙ*(cos_ψD_p_θ - cos_θ) + dₜ*(sin_ψD_p_θ - sin_θ))*cos_ψD_p_θ)*(cos(-θ - ψ_D₀ + ψ₂₀) + b*sin(-θ - ψ_D₀ + ψ₂₀))
    K[3,1] = -kₜ*((cos_ψ1 - b*sin_ψ1)*cos_ψD_m_θ - (-sin_ψ1 - b*cos_ψ1)*sin_ψD_m_θ)*cos_ψD_m_θ + kₙ*(-(cos_ψ1 - b*sin_ψ1)*sin_ψD_m_θ - (-sin_ψ1 - b*cos_ψ1)*cos_ψD_m_θ)*sin_ψD_m_θ
    K[3,2] = kₜ*(-(cos_ψ2 + b*sin_ψ2)*cos_ψD_p_θ + (-sin_ψ2 + b*cos_ψ2)*sin_ψD_p_θ)*cos_ψD_p_θ + kₙ*(-(cos_ψ2 + b*sin_ψ2)*sin_ψD_p_θ - (-sin_ψ2 + b*cos_ψ2)*cos_ψD_p_θ)*sin_ψD_p_θ
    K[3,3] = fₓ^2/α + kₜ*cos_ψD_m_θ^2 + kₙ*sin_ψD_m_θ^2 + kₜ*cos_ψD_p_θ^2 + sin_ψD_p_θ^2*kₙ
    K[3,4] = -kₜ*sin_ψD_m_θ*cos_ψD_m_θ + kₙ*cos_ψD_m_θ*sin_ψD_m_θ - kₜ*sin_ψD_p_θ*cos_ψD_p_θ + kₙ*cos_ψD_p_θ*sin_ψD_p_θ
    K[3,5] = -kₜ*((-dₜ*sin_ψD_m_θ - dₙ*cos_ψD_m_θ)*cos_ψD_m_θ - (sin_ψ1 - b*(1 - cos_ψ1) - x_D₀ + dₜ*(cos_ψD_m_θ - cos_θ) + dₙ*(-sin_ψD_m_θ - sin_θ))*sin_ψD_m_θ - (dₙ*sin_ψD_m_θ - dₜ*cos_ψD_m_θ)*sin_ψD_m_θ - (-1 + cos_ψ1 - b*sin_ψ1 - y_D₀ - dₙ*(cos_ψD_m_θ - cos_θ) + dₜ*(-sin_ψD_m_θ - sin_θ))*cos_ψD_m_θ)*cos_ψD_m_θ + kₜ*((sin_ψ1 - b*(1 - cos_ψ1) - x_D₀ + dₜ*(cos_ψD_m_θ - cos_θ) + dₙ*(-sin_ψD_m_θ - sin_θ))*cos_ψD_m_θ - (-1 + cos_ψ1 - b*sin_ψ1 - y_D₀ - dₙ*(cos_ψD_m_θ - cos_θ) + dₜ*(-sin_ψD_m_θ - sin_θ))*sin_ψD_m_θ)*sin_ψD_m_θ + kₙ*(-(-dₜ*sin_ψD_m_θ - dₙ*cos_ψD_m_θ)*sin_ψD_m_θ - (sin_ψ1 - b*(1 - cos_ψ1) - x_D₀ + dₜ*(cos_ψD_m_θ - cos_θ) + dₙ*(-sin_ψD_m_θ - sin_θ))*cos_ψD_m_θ - (dₙ*sin_ψD_m_θ - dₜ*cos_ψD_m_θ)*cos_ψD_m_θ + (-1 + cos_ψ1 - b*sin_ψ1 - y_D₀ - dₙ*(cos_ψD_m_θ - cos_θ) + dₜ*(-sin_ψD_m_θ - sin_θ))*sin_ψD_m_θ)*sin_ψD_m_θ + kₙ*(-(sin_ψ1 - b*(1 - cos_ψ1) - x_D₀ + dₜ*(cos_ψD_m_θ - cos_θ) + dₙ*(-sin_ψD_m_θ - sin_θ))*sin_ψD_m_θ - (-1 + cos_ψ1 - b*sin_ψ1 - y_D₀ - dₙ*(cos_ψD_m_θ - cos_θ) + dₜ*(-sin_ψD_m_θ - sin_θ))*cos_ψD_m_θ)*cos_ψD_m_θ + kₜ*(-(dₜ*sin_ψD_p_θ - dₙ*cos_ψD_p_θ)*cos_ψD_p_θ + (sin_ψ2 + b*(1 - cos_ψ2) - x_D₀ - dₜ*(cos_ψD_p_θ - cos_θ) - dₙ*(sin_ψD_p_θ - sin_θ))*sin_ψD_p_θ + (dₙ*sin_ψD_p_θ + dₜ*cos_ψD_p_θ)*sin_ψD_p_θ + (-1 + cos_ψ2 + b*sin_ψ2 - y_D₀ - dₙ*(cos_ψD_p_θ - cos_θ) + dₜ*(sin_ψD_p_θ - sin_θ))*cos_ψD_p_θ)*cos_ψD_p_θ - kₜ*(-(sin_ψ2 + b*(1 - cos_ψ2) - x_D₀ - dₜ*(cos_ψD_p_θ - cos_θ) - dₙ*(sin_ψD_p_θ - sin_θ))*cos_ψD_p_θ + (-1 + cos_ψ2 + b*sin_ψ2 - y_D₀ - dₙ*(cos_ψD_p_θ - cos_θ) + dₜ*(sin_ψD_p_θ - sin_θ))*sin_ψD_p_θ)*sin_ψD_p_θ + kₙ*(-(dₜ*sin_ψD_p_θ - dₙ*cos_ψD_p_θ)*sin_ψD_p_θ - (sin_ψ2 + b*(1 - cos_ψ2) - x_D₀ - dₜ*(cos_ψD_p_θ - cos_θ) - dₙ*(sin_ψD_p_θ - sin_θ))*cos_ψD_p_θ - (dₙ*sin_ψD_p_θ + dₜ*cos_ψD_p_θ)*cos_ψD_p_θ + (-1 + cos_ψ2 + b*sin_ψ2 - y_D₀ - dₙ*(cos_ψD_p_θ - cos_θ) + dₜ*(sin_ψD_p_θ - sin_θ))*sin_ψD_p_θ)*sin_ψD_p_θ + kₙ*(-(sin_ψ2 + b*(1 - cos_ψ2) - x_D₀ - dₜ*(cos_ψD_p_θ - cos_θ) - dₙ*(sin_ψD_p_θ - sin_θ))*sin_ψD_p_θ - (-1 + cos_ψ2 + b*sin_ψ2 - y_D₀ - dₙ*(cos_ψD_p_θ - cos_θ) + dₜ*(sin_ψD_p_θ - sin_θ))*cos_ψD_p_θ)*cos_ψD_p_θ
    K[4,1] = kₜ*((cos_ψ1 - b*sin_ψ1)*cos_ψD_m_θ - (-sin_ψ1 - b*cos_ψ1)*sin_ψD_m_θ)*sin_ψD_m_θ + kₙ*(-(cos_ψ1 - b*sin_ψ1)*sin_ψD_m_θ - (-sin_ψ1 - b*cos_ψ1)*cos_ψD_m_θ)*cos_ψD_m_θ
    K[4,2] = -kₜ*(-(cos_ψ2 + b*sin_ψ2)*cos_ψD_p_θ + (-sin_ψ2 + b*cos_ψ2)*sin_ψD_p_θ)*sin_ψD_p_θ + kₙ*(-(cos_ψ2 + b*sin_ψ2)*sin_ψD_p_θ - (-sin_ψ2 + b*cos_ψ2)*cos_ψD_p_θ)*cos_ψD_p_θ
    K[4,3] = -kₜ*sin_ψD_m_θ*cos_ψD_m_θ + kₙ*cos_ψD_m_θ*sin_ψD_m_θ - kₜ*sin_ψD_p_θ*cos_ψD_p_θ + kₙ*cos_ψD_p_θ*sin_ψD_p_θ
    K[4,4] = fₓ^2/α + sin_ψD_m_θ^2*kₜ + cos_ψD_m_θ^2*kₙ + sin_ψD_p_θ^2*kₜ + cos_ψD_p_θ^2*kₙ
    K[4,5] = kₜ*((-dₜ*sin_ψD_m_θ - dₙ*cos_ψD_m_θ)*cos_ψD_m_θ - (sin_ψ1 - b*(1 - cos_ψ1) - x_D₀ + dₜ*(cos_ψD_m_θ - cos_θ) + dₙ*(-sin_ψD_m_θ - sin_θ))*sin_ψD_m_θ - (dₙ*sin_ψD_m_θ - dₜ*cos_ψD_m_θ)*sin_ψD_m_θ - (-1 + cos_ψ1 - b*sin_ψ1 - y_D₀ - dₙ*(cos_ψD_m_θ - cos_θ) + dₜ*(-sin_ψD_m_θ - sin_θ))*cos_ψD_m_θ)*sin_ψD_m_θ + kₜ*((sin_ψ1 - b*(1 - cos_ψ1) - x_D₀ + dₜ*(cos_ψD_m_θ - cos_θ) + dₙ*(-sin_ψD_m_θ - sin_θ))*cos_ψD_m_θ - (-1 + cos_ψ1 - b*sin_ψ1 - y_D₀ - dₙ*(cos_ψD_m_θ - cos_θ) + dₜ*(-sin_ψD_m_θ - sin_θ))*sin_ψD_m_θ)*cos_ψD_m_θ + kₙ*(-(-dₜ*sin_ψD_m_θ - dₙ*cos_ψD_m_θ)*sin_ψD_m_θ - (sin_ψ1 - b*(1 - cos_ψ1) - x_D₀ + dₜ*(cos_ψD_m_θ - cos_θ) + dₙ*(-sin_ψD_m_θ - sin_θ))*cos_ψD_m_θ - (dₙ*sin_ψD_m_θ - dₜ*cos_ψD_m_θ)*cos_ψD_m_θ + (-1 + cos_ψ1 - b*sin_ψ1 - y_D₀ - dₙ*(cos_ψD_m_θ - cos_θ) + dₜ*(-sin_ψD_m_θ - sin_θ))*sin_ψD_m_θ)*cos_ψD_m_θ - kₙ*(-(sin_ψ1 - b*(1 - cos_ψ1) - x_D₀ + dₜ*(cos_ψD_m_θ - cos_θ) + dₙ*(-sin_ψD_m_θ - sin_θ))*sin_ψD_m_θ - (-1 + cos_ψ1 - b*sin_ψ1 - y_D₀ - dₙ*(cos_ψD_m_θ - cos_θ) + dₜ*(-sin_ψD_m_θ - sin_θ))*cos_ψD_m_θ)*sin_ψD_m_θ - kₜ*(-(dₜ*sin_ψD_p_θ - dₙ*cos_ψD_p_θ)*cos_ψD_p_θ + (sin_ψ2 + b*(1 - cos_ψ2) - x_D₀ - dₜ*(cos_ψD_p_θ - cos_θ) - dₙ*(sin_ψD_p_θ - sin_θ))*sin_ψD_p_θ + (dₙ*sin_ψD_p_θ + dₜ*cos_ψD_p_θ)*sin_ψD_p_θ + (-1 + cos_ψ2 + b*sin_ψ2 - y_D₀ - dₙ*(cos_ψD_p_θ - cos_θ) + dₜ*(sin_ψD_p_θ - sin_θ))*cos_ψD_p_θ)*sin_ψD_p_θ - kₜ*(-(sin_ψ2 + b*(1 - cos_ψ2) - x_D₀ - dₜ*(cos_ψD_p_θ - cos_θ) - dₙ*(sin_ψD_p_θ - sin_θ))*cos_ψD_p_θ + (-1 + cos_ψ2 + b*sin_ψ2 - y_D₀ - dₙ*(cos_ψD_p_θ - cos_θ) + dₜ*(sin_ψD_p_θ - sin_θ))*sin_ψD_p_θ)*cos_ψD_p_θ + kₙ*(-(dₜ*sin_ψD_p_θ - dₙ*cos_ψD_p_θ)*sin_ψD_p_θ - (sin_ψ2 + b*(1 - cos_ψ2) - x_D₀ - dₜ*(cos_ψD_p_θ - cos_θ) - dₙ*(sin_ψD_p_θ - sin_θ))*cos_ψD_p_θ - (dₙ*sin_ψD_p_θ + dₜ*cos_ψD_p_θ)*cos_ψD_p_θ + (-1 + cos_ψ2 + b*sin_ψ2 - y_D₀ - dₙ*(cos_ψD_p_θ - cos_θ) + dₜ*(sin_ψD_p_θ - sin_θ))*sin_ψD_p_θ)*cos_ψD_p_θ - kₙ*(-(sin_ψ2 + b*(1 - cos_ψ2) - x_D₀ - dₜ*(cos_ψD_p_θ - cos_θ) - dₙ*(sin_ψD_p_θ - sin_θ))*sin_ψD_p_θ - (-1 + cos_ψ2 + b*sin_ψ2 - y_D₀ - dₙ*(cos_ψD_p_θ - cos_θ) + dₜ*(sin_ψD_p_θ - sin_θ))*cos_ψD_p_θ)*sin_ψD_p_θ
    K[5,1] = kₜ*((cos_ψ1 - b*sin_ψ1)*cos_ψD_m_θ - (-sin_ψ1 - b*cos_ψ1)*sin_ψD_m_θ)*(cos_ψD_m_θ + b*sin_ψD_m_θ - cos(θ - ψ_D₀ + ψ₁₀) + b*sin(θ - ψ_D₀ + ψ₁₀) + x_D₀*sin_ψD_m_θ + y_D₀*cos_ψD_m_θ + dₜ   *sin(ψ_D₀) - dₙ*cos(ψ_D₀)) + kₜ*((sin_ψ1 - b*(1 - cos_ψ1) - x_D₀ + dₜ   *(cos_ψD_m_θ - cos_θ) + dₙ*(-sin_ψD_m_θ - sin_θ))*cos_ψD_m_θ - (-1 + cos_ψ1 - b*sin_ψ1 - y_D₀ - dₙ*(cos_ψD_m_θ - cos_θ) + dₜ *(-sin_ψD_m_θ - sin_θ))*sin_ψD_m_θ)*(sin(θ - ψ_D₀ + ψ₁₀) + b*cos(θ - ψ_D₀ + ψ₁₀)) + kₙ*(-(cos_ψ1 - b*sin_ψ1)*sin_ψD_m_θ - (-sin_ψ1 - b*cos_ψ1)*cos_ψD_m_θ)*(-sin_ψD_m_θ + b*cos_ψD_m_θ - sin(θ - ψ_D₀ + ψ₁₀) - b*cos(θ - ψ_D₀ + ψ₁₀) + x_D₀*cos_ψD_m_θ - y_D₀*sin_ψD_m_θ + dₜ   *cos(ψ_D₀) + dₙ*sin(ψ_D₀)) + kₙ*(-(sin_ψ1 - b*(1 - cos_ψ1) - x_D₀ + dₜ   *(cos_ψD_m_θ - cos_θ) + dₙ*(-sin_ψD_m_θ - sin_θ))*sin_ψD_m_θ - (-1 + cos_ψ1 - b*sin_ψ1 - y_D₀ - dₙ*(cos_ψD_m_θ - cos_θ) + dₜ *(-sin_ψD_m_θ - sin_θ))*cos_ψD_m_θ)*(-cos(θ - ψ_D₀ + ψ₁₀) + b*sin(θ - ψ_D₀ + ψ₁₀))
    K[5,2] = kₜ*(-(cos_ψ2 + b*sin_ψ2)*cos(θ + ψ_D₀) + (-sin_ψ2 + b*cos_ψ2)*sin_ψD_p_θ)*(-cos(θ + ψ_D₀) + b*sin_ψD_p_θ + cos(-θ - ψ_D₀ + ψ₂₀) + b*sin(-θ - ψ_D₀ + ψ₂₀) - x_D₀*sin_ψD_p_θ - y_D₀*cos(θ + ψ_D₀) + dₜ*sin(ψ_D₀) + dₙ*cos(ψ_D₀)) + kₜ*(-(sin_ψ2 + b*(1 - cos_ψ2) - x_D₀ - dₜ   *(cos(θ + ψ_D₀) - cos_θ) - dₙ*(sin_ψD_p_θ - sin_θ))*cos(θ + ψ_D₀) + (-1 + cos_ψ2 + b*sin_ψ2 - y_D₀ - dₙ*(cos(θ + ψ_D₀) - cos_θ) + dₜ   *(sin_ψD_p_θ - sin_θ))*sin_ψD_p_θ)*(-sin(-θ - ψ_D₀ + ψ₂₀) + b*cos(-θ - ψ_D₀ + ψ₂₀)) + kₙ*(-(cos_ψ2 + b*sin_ψ2)*sin_ψD_p_θ - (-sin_ψ2 + b*cos_ψ2)*cos(θ + ψ_D₀))*(-sin_ψD_p_θ - b*cos(θ + ψ_D₀) - sin(-θ - ψ_D₀ + ψ₂₀) + b*cos(-θ - ψ_D₀ + ψ₂₀) + x_D₀*cos(θ + ψ_D₀) - y_D₀*sin_ψD_p_θ - dₜ  *cos(ψ_D₀) + dₙ*sin(ψ_D₀)) + kₙ*(-(sin_ψ2 + b*(1 - cos_ψ2) - x_D₀ - dₜ  *(cos(θ + ψ_D₀) - cos_θ) - dₙ*(sin_ψD_p_θ - sin_θ))*sin_ψD_p_θ - (-1 + cos_ψ2 + b*sin_ψ2 - y_D₀ - dₙ*(cos(θ + ψ_D₀) - cos_θ) + dₜ   *(sin_ψD_p_θ - sin_θ))*cos(θ + ψ_D₀))*(-cos(-θ - ψ_D₀ + ψ₂₀) - b*sin(-θ - ψ_D₀ + ψ₂₀))
    K[5,3] = -kₜ*cos_ψD_m_θ*(cos_ψD_m_θ + b*sin_ψD_m_θ - cos(θ - ψ_D₀ + ψ₁₀) + b*sin(θ - ψ_D₀ + ψ₁₀) + x_D₀*sin_ψD_m_θ + y_D₀*cos_ψD_m_θ + dₜ  *sin(ψ_D₀) - dₙ*cos(ψ_D₀)) + kₜ*((sin_ψ1 - b*(1 - cos_ψ1) - x_D₀ + dₜ  *(cos_ψD_m_θ - cos_θ) + dₙ*(-sin_ψD_m_θ - sin_θ))*cos_ψD_m_θ - (-1 + cos_ψ1 - b*sin_ψ1 - y_D₀ - dₙ*(cos_ψD_m_θ - cos_θ) + dₜ*(-sin_ψD_m_θ - sin_θ))*sin_ψD_m_θ)*sin_ψD_m_θ + kₜ*cos(θ + ψ_D₀)*(-cos(θ + ψ_D₀) + b*sin_ψD_p_θ + cos(-θ - ψ_D₀ + ψ₂₀) + b*sin(-θ - ψ_D₀ + ψ₂₀) - x_D₀*sin_ψD_p_θ - y_D₀*cos(θ + ψ_D₀) + dₜ  *sin(ψ_D₀) + dₙ*cos(ψ_D₀)) - kₜ*(-(sin_ψ2 + b*(1 - cos_ψ2) - x_D₀ - dₜ *(cos(θ + ψ_D₀) - cos_θ) - dₙ*(sin_ψD_p_θ - sin_θ))*cos(θ + ψ_D₀) + (-1 + cos_ψ2 + b*sin_ψ2 - y_D₀ - dₙ*(cos(θ + ψ_D₀) - cos_θ) + dₜ *(sin_ψD_p_θ - sin_θ))*sin_ψD_p_θ)*sin_ψD_p_θ + kₙ*sin_ψD_m_θ*(-sin_ψD_m_θ + b*cos_ψD_m_θ - sin(θ - ψ_D₀ + ψ₁₀) - b*cos(θ - ψ_D₀ + ψ₁₀) + x_D₀*cos_ψD_m_θ - y_D₀*sin_ψD_m_θ + dₜ  *cos(ψ_D₀) + dₙ*sin(ψ_D₀)) + kₙ*(-(sin_ψ1 - b*(1 - cos_ψ1) - x_D₀ + dₜ  *(cos_ψD_m_θ - cos_θ) + dₙ*(-sin_ψD_m_θ - sin_θ))*sin_ψD_m_θ - (-1 + cos_ψ1 - b*sin_ψ1 - y_D₀ - dₙ*(cos_ψD_m_θ - cos_θ) + dₜ*(-sin_ψD_m_θ - sin_θ))*cos_ψD_m_θ)*cos_ψD_m_θ + kₙ*sin_ψD_p_θ*(-sin_ψD_p_θ - b*cos(θ + ψ_D₀) - sin(-θ - ψ_D₀ + ψ₂₀) + b*cos(-θ - ψ_D₀ + ψ₂₀) + x_D₀*cos(θ + ψ_D₀) - y_D₀*sin_ψD_p_θ - dₜ*cos(ψ_D₀) + dₙ*sin(ψ_D₀)) + kₙ*(-(sin_ψ2 + b*(1 - cos_ψ2) - x_D₀ - dₜ*(cos(θ + ψ_D₀) - cos_θ) - dₙ*(sin_ψD_p_θ - sin_θ))*sin_ψD_p_θ - (-1 + cos_ψ2 + b*sin_ψ2 - y_D₀ - dₙ*(cos(θ + ψ_D₀) - cos_θ) + dₜ *(sin_ψD_p_θ - sin_θ))*cos(θ + ψ_D₀))*cos(θ + ψ_D₀)
    K[5,4] = kₜ*sin_ψD_m_θ*(cos_ψD_m_θ + b*sin_ψD_m_θ - cos(θ - ψ_D₀ + ψ₁₀) + b*sin(θ - ψ_D₀ + ψ₁₀) + x_D₀*sin_ψD_m_θ + y_D₀*cos_ψD_m_θ + dₜ   *sin(ψ_D₀) - dₙ*cos(ψ_D₀)) + kₜ*((sin_ψ1 - b*(1 - cos_ψ1) - x_D₀ + dₜ   *(cos_ψD_m_θ - cos_θ) + dₙ*(-sin_ψD_m_θ - sin_θ))*cos_ψD_m_θ - (-1 + cos_ψ1 - b*sin_ψ1 - y_D₀ - dₙ*(cos_ψD_m_θ - cos_θ) + dₜ *(-sin_ψD_m_θ - sin_θ))*sin_ψD_m_θ)*cos_ψD_m_θ - kₜ*sin_ψD_p_θ*(-cos(θ + ψ_D₀) + b*sin_ψD_p_θ + cos(-θ - ψ_D₀ + ψ₂₀) + b*sin(-θ - ψ_D₀ + ψ₂₀) - x_D₀*sin_ψD_p_θ - y_D₀*cos(θ + ψ_D₀) + dₜ*sin(ψ_D₀) + dₙ*cos(ψ_D₀)) - kₜ*(-(sin_ψ2 + b*(1 - cos_ψ2) - x_D₀ - dₜ   *(cos(θ + ψ_D₀) - cos_θ) - dₙ*(sin_ψD_p_θ - sin_θ))*cos(θ + ψ_D₀) + (-1 + cos_ψ2 + b*sin_ψ2 - y_D₀ - dₙ*(cos(θ + ψ_D₀) - cos_θ) + dₜ   *(sin_ψD_p_θ - sin_θ))*sin_ψD_p_θ)*cos(θ + ψ_D₀) + kₙ*cos_ψD_m_θ*(-sin_ψD_m_θ + b*cos_ψD_m_θ - sin(θ - ψ_D₀ + ψ₁₀) - b*cos(θ - ψ_D₀ + ψ₁₀) + x_D₀*cos_ψD_m_θ - y_D₀*sin_ψD_m_θ + dₜ   *cos(ψ_D₀) + dₙ*sin(ψ_D₀)) - kₙ*(-(sin_ψ1 - b*(1 - cos_ψ1) - x_D₀ + dₜ   *(cos_ψD_m_θ - cos_θ) + dₙ*(-sin_ψD_m_θ - sin_θ))*sin_ψD_m_θ - (-1 + cos_ψ1 - b*sin_ψ1 - y_D₀ - dₙ*(cos_ψD_m_θ - cos_θ) + dₜ *(-sin_ψD_m_θ - sin_θ))*cos_ψD_m_θ)*sin_ψD_m_θ + kₙ*cos(θ + ψ_D₀)*(-sin_ψD_p_θ - b*cos(θ + ψ_D₀) - sin(-θ - ψ_D₀ + ψ₂₀) + b*cos(-θ - ψ_D₀ + ψ₂₀) + x_D₀*cos(θ + ψ_D₀) - y_D₀*sin_ψD_p_θ - dₜ*cos(ψ_D₀) + dₙ*sin(ψ_D₀)) - kₙ*(-(sin_ψ2 + b*(1 - cos_ψ2) - x_D₀ - dₜ*(cos(θ + ψ_D₀) - cos_θ) - dₙ*(sin_ψD_p_θ - sin_θ))*sin_ψD_p_θ - (-1 + cos_ψ2 + b*sin_ψ2 - y_D₀ - dₙ*(cos(θ + ψ_D₀) - cos_θ) + dₜ *(sin_ψD_p_θ - sin_θ))*cos(θ + ψ_D₀))*sin_ψD_p_θ
    K[5,5] = f_ψ^2/β + kₜ*((-dₜ*sin_ψD_m_θ - dₙ*cos_ψD_m_θ)*cos_ψD_m_θ - (sin_ψ1 - b*(1 - cos_ψ1) - x_D₀ + dₜ *(cos_ψD_m_θ - cos_θ) + dₙ*(-sin_ψD_m_θ - sin_θ))*sin_ψD_m_θ - (dₙ*sin_ψD_m_θ - dₜ   *cos_ψD_m_θ)*sin_ψD_m_θ - (-1 + cos_ψ1 - b*sin_ψ1 - y_D₀ - dₙ*(cos_ψD_m_θ - cos_θ) + dₜ   *(-sin_ψD_m_θ - sin_θ))*cos_ψD_m_θ)*(cos_ψD_m_θ + b*sin_ψD_m_θ - cos(θ - ψ_D₀ + ψ₁₀) + b*sin(θ - ψ_D₀ + ψ₁₀) + x_D₀*sin_ψD_m_θ + y_D₀*cos_ψD_m_θ + dₜ*sin(ψ_D₀) - dₙ*cos(ψ_D₀)) + kₜ*((sin_ψ1 - b*(1 - cos_ψ1) - x_D₀ + dₜ*(cos_ψD_m_θ - cos_θ) + dₙ*(-sin_ψD_m_θ - sin_θ))*cos_ψD_m_θ - (-1 + cos_ψ1 - b*sin_ψ1 - y_D₀ - dₙ*(cos_ψD_m_θ - cos_θ) + dₜ  *(-sin_ψD_m_θ - sin_θ))*sin_ψD_m_θ)*(-sin_ψD_m_θ + b*cos_ψD_m_θ - sin(θ - ψ_D₀ + ψ₁₀) - b*cos(θ - ψ_D₀ + ψ₁₀) + x_D₀*cos_ψD_m_θ - y_D₀*sin_ψD_m_θ + dₜ  *cos(ψ_D₀) + dₙ*sin(ψ_D₀)) + kₜ*(-(dₜ*sin_ψD_p_θ - dₙ*cos(θ + ψ_D₀))*cos(θ + ψ_D₀) + (sin_ψ2 + b*(1 - cos_ψ2) - x_D₀ - dₜ   *(cos(θ + ψ_D₀) - cos_θ) - dₙ*(sin_ψD_p_θ - sin_θ))*sin_ψD_p_θ + (dₙ*sin_ψD_p_θ + dₜ *cos(θ + ψ_D₀))*sin_ψD_p_θ + (-1 + cos_ψ2 + b*sin_ψ2 - y_D₀ - dₙ*(cos(θ + ψ_D₀) - cos_θ) + dₜ   *(sin_ψD_p_θ - sin_θ))*cos(θ + ψ_D₀))*(-cos(θ + ψ_D₀) + b*sin_ψD_p_θ + cos(-θ - ψ_D₀ + ψ₂₀) + b*sin(-θ - ψ_D₀ + ψ₂₀) - x_D₀*sin_ψD_p_θ - y_D₀*cos(θ + ψ_D₀) + dₜ   *sin(ψ_D₀) + dₙ*cos(ψ_D₀)) + kₜ*(-(sin_ψ2 + b*(1 - cos_ψ2) - x_D₀ - dₜ  *(cos(θ + ψ_D₀) - cos_θ) - dₙ*(sin_ψD_p_θ - sin_θ))*cos(θ + ψ_D₀) + (-1 + cos_ψ2 + b*sin_ψ2 - y_D₀ - dₙ*(cos(θ + ψ_D₀) - cos_θ) + dₜ  *(sin_ψD_p_θ - sin_θ))*sin_ψD_p_θ)*(sin_ψD_p_θ + b*cos(θ + ψ_D₀) + sin(-θ - ψ_D₀ + ψ₂₀) - b*cos(-θ - ψ_D₀ + ψ₂₀) - x_D₀*cos(θ + ψ_D₀) + y_D₀*sin_ψD_p_θ + dₜ*cos(ψ_D₀) - dₙ*sin(ψ_D₀)) + kₙ*(-(-dₜ  *sin_ψD_m_θ - dₙ*cos_ψD_m_θ)*sin_ψD_m_θ - (sin_ψ1 - b*(1 - cos_ψ1) - x_D₀ + dₜ   *(cos_ψD_m_θ - cos_θ) + dₙ*(-sin_ψD_m_θ - sin_θ))*cos_ψD_m_θ - (dₙ*sin_ψD_m_θ - dₜ *cos_ψD_m_θ)*cos_ψD_m_θ + (-1 + cos_ψ1 - b*sin_ψ1 - y_D₀ - dₙ*(cos_ψD_m_θ - cos_θ) + dₜ *(-sin_ψD_m_θ - sin_θ))*sin_ψD_m_θ)*(-sin_ψD_m_θ + b*cos_ψD_m_θ - sin(θ - ψ_D₀ + ψ₁₀) - b*cos(θ - ψ_D₀ + ψ₁₀) + x_D₀*cos_ψD_m_θ - y_D₀*sin_ψD_m_θ + dₜ *cos(ψ_D₀) + dₙ*sin(ψ_D₀)) + kₙ*(-(sin_ψ1 - b*(1 - cos_ψ1) - x_D₀ + dₜ *(cos_ψD_m_θ - cos_θ) + dₙ*(-sin_ψD_m_θ - sin_θ))*sin_ψD_m_θ - (-1 + cos_ψ1 - b*sin_ψ1 - y_D₀ - dₙ*(cos_ψD_m_θ - cos_θ) + dₜ   *(-sin_ψD_m_θ - sin_θ))*cos_ψD_m_θ)*(-cos_ψD_m_θ - b*sin_ψD_m_θ + cos(θ - ψ_D₀ + ψ₁₀) - b*sin(θ - ψ_D₀ + ψ₁₀) - x_D₀*sin_ψD_m_θ - y_D₀*cos_ψD_m_θ - dₜ   *sin(ψ_D₀) + dₙ*cos(ψ_D₀)) + kₙ*(-(dₜ  *sin_ψD_p_θ - dₙ*cos(θ + ψ_D₀))*sin_ψD_p_θ - (sin_ψ2 + b*(1 - cos_ψ2) - x_D₀ - dₜ  *(cos(θ + ψ_D₀) - cos_θ) - dₙ*(sin_ψD_p_θ - sin_θ))*cos(θ + ψ_D₀) - (dₙ*sin_ψD_p_θ + dₜ   *cos(θ + ψ_D₀))*cos(θ + ψ_D₀) + (-1 + cos_ψ2 + b*sin_ψ2 - y_D₀ - dₙ*(cos(θ + ψ_D₀) - cos_θ) + dₜ*(sin_ψD_p_θ - sin_θ))*sin_ψD_p_θ)*(-sin_ψD_p_θ - b*cos(θ + ψ_D₀) - sin(-θ - ψ_D₀ + ψ₂₀) + b*cos(-θ - ψ_D₀ + ψ₂₀) + x_D₀*cos(θ + ψ_D₀) - y_D₀*sin_ψD_p_θ - dₜ *cos(ψ_D₀) + dₙ*sin(ψ_D₀)) + kₙ*(-(sin_ψ2 + b*(1 - cos_ψ2) - x_D₀ - dₜ *(cos(θ + ψ_D₀) - cos_θ) - dₙ*(sin_ψD_p_θ - sin_θ))*sin_ψD_p_θ - (-1 + cos_ψ2 + b*sin_ψ2 - y_D₀ - dₙ*(cos(θ + ψ_D₀) - cos_θ) + dₜ  *(sin_ψD_p_θ - sin_θ))*cos(θ + ψ_D₀))*(-cos(θ + ψ_D₀) + b*sin_ψD_p_θ + cos(-θ - ψ_D₀ + ψ₂₀) + b*sin(-θ - ψ_D₀ + ψ₂₀) - x_D₀*sin_ψD_p_θ - y_D₀*cos(θ + ψ_D₀) + dₜ  *sin(ψ_D₀) + dₙ*cos(ψ_D₀))
    
    return K
end

# Function to compute the static terms of the system
function static_5DOF(x, p)
    @unpack θ, b, dₙ, dₜ, fₓ, f_ψ, α, β, friction_law_1 = p
    kₙ = friction_law_1.kₙ
    kₜ = friction_law_1.kₜ

    F = zeros(5)
    
    # Recurrent trigonometric functions
    sin_θ = sin(θ)
    cos_θ = cos(θ)
    sin_ψ₁ = sin(x[1])
    cos_ψ₁ = cos(x[1])
    sin_ψ₂ = sin(x[2])
    cos_ψ₂ = cos(x[2])
    

    # Relative displacement of the contact points
    xₜ₁, xₙ₁, xₜ₂, xₙ₂ = contact_points_displacements(x, p)

    # Contact forces
    T₁ = kₜ*xₜ₁
    N₁ = kₙ*xₙ₁
    T₂ = kₜ*xₜ₂
    N₂ = kₙ*xₙ₂
    
    # Static equations
    F[1] = x[1] - T₁*( -cos(θ - x[5] + x[1]) + b*sin(θ - x[5] + x[1]) ) - N₁*( -sin(θ - x[5] + x[1]) - b*cos(θ - x[5] + x[1]) )
    F[2] = x[2] - T₂*( cos(θ + x[5] - x[2]) - b*sin(θ + x[5] - x[2]) ) - N₂*( sin(θ + x[5] - x[2]) + b*cos(θ + x[5] - x[2]) )
    F[3] = fₓ^2*x[3]/α - ( T₁*cos(θ - x[5]) + N₁*sin(θ - x[5]) - T₂*cos(θ + x[5]) - N₂*sin(θ + x[5]) )
    F[4] = fₓ^2*x[4]/α - ( T₁*sin(θ - x[5]) - N₁*cos(θ - x[5]) + T₂*sin(θ + x[5]) - N₂*cos(θ + x[5]) )
    F[5] = f_ψ^2*x[5]/β - ( T₁*(dₙ - xₙ₁) - N₁*(dₜ - xₜ₁) - T₂*(dₙ - xₙ₂) + N₂*(dₜ - xₜ₂) )
    
    return F
end

# Function to calculate the stiffness matrix numerically
function num_stiffness(x₀, F, Δx, p)
    n = length(x₀)
    K = zeros(n, n)
    
    for i in 1:n
        x₁ = copy(x₀)
        x₂ = copy(x₀)
        x₁[i] += Δx
        x₂[i] -= Δx
        
        F₁ = F(x₁, p)
        F₂ = F(x₂, p)
        
        K[:,i] = (F₁ .- F₂)./(2Δx)
    end
    
    return K
end
