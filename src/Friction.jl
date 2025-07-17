using Parameters: @unpack

abstract type FrictionLaw end

## 1D FRICTION LAWS
function friction_type(law::FrictionLaw)
    error("Friction law type not defined for $(law).")
end

# LINEAR FRICTION
struct Linear{T <: Real} <: FrictionLaw
    kₙ::T
    kₜ::T
end

Linear(; kₙ, kₜ) = Linear(kₙ, kₜ) 

g(xₜ, xₙ, law::Linear) = [law.kₜ * xₜ, law.kₙ * xₙ]

function friction_type(::Linear)
    return "linear"
end

# CUBIC FRICTION
struct Cubic{T <: Real} <: FrictionLaw
    kₙ::T
    kₜ::T
end

Cubic(; kₙ, kₜ) = Cubic(kₙ, kₜ) 

g(xₜ, xₙ, law::Cubic) = [law.kₜ * xₜ^3, law.kₙ * xₙ^3]

function friction_type(::Cubic)
    return "cubic"
end

# COULOMB FRICTION
mutable struct Coulomb{T <: Real} <: FrictionLaw
    const kₙ::T
    const xₙ₀::T
    const kₜ::T
    const μ::T
    wₜ::T
end

function Coulomb(kₙ, xₙ₀, kₜ, μ)
    wₜ = zero(kₙ)
    return Coulomb(kₙ, xₙ₀, kₜ, μ, wₜ)
end

Coulomb(; kₙ, xₙ₀, kₜ, μ) = Coulomb(kₙ, xₙ₀, kₜ, μ)

function g(xₜ, xₙ, law::Coulomb)
    @unpack kₙ, xₙ₀, kₜ, μ = law

    N = normal_force(xₙ, kₙ, xₙ₀)
    T, wₜ = tangential_force(xₜ, law.wₜ, kₜ, μ, N)
    law.wₜ = wₜ

    return [T, N]
end

@inline function normal_force(xₙ, kₙ, xₙ₀)
    u = xₙ - xₙ₀
    if u > 0.0
        return kₙ * u
    else
        return zero(xₙ)
    end
end

@inline function tangential_force(xₜ, wₜ, kₜ, μ, N)
    if N > 0.0
        T = kₜ * (xₜ - wₜ)
        if abs(T) < μ * N
            # Stick regime
            return T, wₜ
        else
            # Slip regime
            sg = sign(T)
            wₜ = xₜ - sg * μ * N / kₜ
            return sg * μ * N, wₜ
        end
    else
        # Lift-off
        wₜ = xₜ
        return zero(xₜ), wₜ
    end
end

function friction_type(::Coulomb)
    return "coulomb"
end