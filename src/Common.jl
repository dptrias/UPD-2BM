include("Friction.jl")

## PARAMETERS STRUCTS
mutable struct TBM_5DOF_Parameters 
    ζ_ψ::Float64
    α::Float64
    β::Float64
    fₓ::Float64
    f_ψ::Float64
    θ::Float64
    b::Float64
    dₙ::Float64
    dₜ::Float64
    Mₑ::Float64
    feq::Union{Float64, Nothing}
    FC::Float64
    T_F::Union{Float64, Nothing}
    friction_law_1::Union{FrictionLaw, Nothing}
    friction_law_2::Union{FrictionLaw, Nothing}
    T₁::Union{Vector{Float64}, Nothing}
    N₁::Union{Vector{Float64}, Nothing}
    T₂::Union{Vector{Float64}, Nothing}
    N₂::Union{Vector{Float64}, Nothing}
    xₜ₁::Union{Vector{Float64}, Nothing}
    xₙ₁::Union{Vector{Float64}, Nothing}
    xₜ₂::Union{Vector{Float64}, Nothing}
    xₙ₂::Union{Vector{Float64}, Nothing}
    index::Union{Base.RefValue{Int}, Nothing}
    
    function TBM_5DOF_Parameters(;
        ζ_ψ = 0.0, α = 0.0, β = 0.0, fₓ = 0.0, f_ψ = 0.0, θ = 0.0, b = 0.0, dₙ = 0.0, dₜ = 0.0, 
        Mₑ = 0.0, feq=nothing, FC = 0.0, T_F=nothing,
        friction_law_1=Linear(kₜ = 0.0, kₙ = 0.0),
        friction_law_2=Linear(kₜ = 0.0, kₙ = 0.0),
        T₁ = nothing,
        N₁ = nothing,
        T₂ = nothing,
        N₂ = nothing,
        xₜ₁ = nothing,
        xₙ₁ = nothing,
        xₜ₂ = nothing,
        xₙ₂ = nothing,
        index = nothing
    )
        new(ζ_ψ, α, β, fₓ, f_ψ, θ, b, dₙ, dₜ, Mₑ, feq, FC, T_F, 
            friction_law_1, friction_law_2, 
            T₁, N₁, T₂, N₂, xₜ₁, xₙ₁, xₜ₂, xₙ₂, index
        )
        
    end

    function TBM_5DOF_Parameters(;
        I_B, k_ψ, c_ψ = 0.0, 
        I_D, m_D, k_Dx = 0.0, k_Dψ = 0.0, 
        θ, FC = 0.0, Mₑ = 0.0, kₙ, kₜ, μ,
        f, b, dₙ, dₜ
    )
        ω_ψ = √(k_ψ / I_B) 
        ζ_ψ = c_ψ / (2 * √(k_ψ * I_B))
        ω_Dx = √(k_Dx / m_D)
        ω_Dψ = √(k_Dψ / I_D)
        fₓ = ω_Dx / ω_ψ
        f_ψ = ω_Dψ / ω_ψ
        α = I_B / (m_D * f^2)
        β = I_B / I_D
        b = b / f
        dₙ = dₙ / f
        dₜ = dₜ / f
        FC = FC * f / k_ψ
        kₙ = kₙ * f^2 / k_ψ
        kₜ = kₜ * f^2 / k_ψ

        new(ζ_ψ, α, β, fₓ, f_ψ, θ, b, dₙ, dₜ, Mₑ, nothing, FC, nothing,
            Coulomb(kₙ = kₙ, xₙ₀ = 0.0, kₜ = kₜ, μ = μ),
            Coulomb(kₙ = kₙ, xₙ₀ = 0.0, kₜ = kₜ, μ = μ),
            nothing, nothing,
            nothing, nothing,
            nothing, nothing, 
            nothing, nothing, 
            nothing
        )
    end

end

mutable struct TBM_2DOF_Parameters 
    ζ_ψ::Float64
    α::Float64
    fₓ::Float64
    θ::Float64
    b::Float64
    Mₑ::Float64
    feq::Union{Float64, Nothing}
    FC::Float64
    T_F::Union{Float64, Nothing}
    friction_law::FrictionLaw
    T::Union{Vector{Float64}, Nothing}
    N::Union{Vector{Float64}, Nothing}
    xₜ::Union{Vector{Float64}, Nothing}
    xₙ::Union{Vector{Float64}, Nothing}
    index::Union{Base.RefValue{Int}, Nothing}
    
    function TBM_2DOF_Parameters(;
        ζ_ψ=0.0, α=0.0, fₓ=0.0, θ=0.0, b=0.0, 
        Mₑ=0.0, feq=nothing, FC=0.0, T_F = nothing, 
        friction_law=Linear(kₙ=0.0, kₜ=0.0),
        T = nothing,
        N = nothing,
        xₜ = nothing,
        xₙ = nothing,
        index = nothing
        )
        new(ζ_ψ, α, fₓ, θ, b, Mₑ, feq, FC, T_F, friction_law, T, N, xₜ, xₙ, index)
        
    end
end

## CONTACT POINT DISPLACEMENTS
function contact_points_displacements(x, p::TBM_5DOF_Parameters)
    @unpack θ, b, dₙ, dₜ = p

    # Recurrent trigonometric functions
    sin_θ = sin(θ)
    cos_θ = cos(θ)
    sin_ψ₁ = sin(x[1])
    cos_ψ₁ = cos(x[1])
    sin_ψ₂ = sin(x[2])
    cos_ψ₂ = cos(x[2])
    sin_θ_p_ψ = sin(θ + x[5])
    cos_θ_p_ψ = cos(θ + x[5])
    sin_θ_m_ψ = sin(θ - x[5])
    cos_θ_m_ψ = cos(θ - x[5])

    # Position of the contact points
    x_B = sin_ψ₁ - b*(1 - cos_ψ₁)
    y_B = -(1 - cos_ψ₁) - b*sin_ψ₁
    x_E = x[3] - dₜ*(cos_θ_m_ψ - cos_θ) - dₙ*(sin_θ_m_ψ - sin_θ)
    y_E = x[4] + dₙ*(cos_θ_m_ψ - cos_θ) - dₜ*(sin_θ_m_ψ - sin_θ)
    x_A = sin_ψ₂ + b*(1 - cos_ψ₂)
    y_A = -(1 - cos_ψ₂) + b*sin_ψ₂
    x_F = x[3] + dₜ*(cos_θ_p_ψ - cos_θ) + dₙ*(sin_θ_p_ψ - sin_θ)
    y_F = x[4] + dₙ*(cos_θ_p_ψ - cos_θ) - dₜ*(sin_θ_p_ψ - sin_θ)
    
    # Relative position of the contact
    xₜ₁ = (x_B - x_E)*cos_θ_m_ψ + (y_B - y_E)*sin_θ_m_ψ
    xₙ₁ = (x_B - x_E)*sin_θ_m_ψ - (y_B - y_E)*cos_θ_m_ψ
    xₜ₂ = -(x_A - x_F)*cos_θ_p_ψ + (y_A - y_F)*sin_θ_p_ψ
    xₙ₂ = -(x_A - x_F)*sin_θ_p_ψ - (y_A - y_F)*cos_θ_p_ψ

    return xₜ₁, xₙ₁, xₜ₂, xₙ₂
end

#=
mutable struct TBMSParameters 
    feq::Float64
    ζ_ψ::Float64
    α::Float64
    β::Float64
    fₓ::Float64
    f_ψ::Float64
    θ::Float64
    b::Float64
    dₙ::Float64
    dₜ::Float64
    Mₑ::Float64
    FC::Float64
    T₁::Vector{Float64}
    N₁::Vector{Float64}
    T₂::Vector{Float64}
    N₂::Vector{Float64}
    xₜ₁::Vector{Float64}
    xₙ₁::Vector{Float64}
    xₜ₂::Vector{Float64}
    xₙ₂::Vector{Float64}
    index::Base.RefValue{Int}
    
    function TBMSParameters(;
        feq=0.0, ζ_ψ=0.0, α=0.0, β=0.0, fₓ=0.0, f_ψ=0.0, θ=0.0, b=0.0, dₙ=0.0, dₜ=0.0, Mₑ=0.0, FC=0.0,
        T₁ = zeros(1000),
        N₁ = zeros(1000),
        T₂ = zeros(1000),
        N₂ = zeros(1000),
        xₜ₁ = zeros(1000),
        xₙ₁ = zeros(1000),
        xₜ₂ = zeros(1000),
        xₙ₂ = zeros(1000),
        index = Ref(1)
        )
        new(feq, ζ_ψ, α, β, fₓ, f_ψ, θ, b, dₙ, dₜ, Mₑ, FC, T₁, N₁, T₂, N₂, xₜ₁, xₙ₁, xₜ₂, xₙ₂, index)
        
    end
end
=#