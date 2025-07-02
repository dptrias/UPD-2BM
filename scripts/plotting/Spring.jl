struct CoiledSpring{T<:Real}
    width::T
    coils::Int
end

"""
    prolate_cycloid(a, b, t)

Compute the prolate cycloid curve for a radius of `a` and a distance `b` from the center
of the circle as a function of the parameter `t`
"""
function prolate_cycloid(a, b, t)
    x = @. a * t - b * sin(t)
    y = @. -b * cos(t)
    return x, y
end


"""
    spring_coordinates(x₀, y₀, x₁, y₁, spring::CoiledSpring)

Compute the `x`, `y` coordinates of a `CoiledSpring`. The initial point of the spring is
`(x₀, y₀)` and the end point is `(x₁, y₁)`
"""
function spring_coordinates(x₀, y₀, x₁, y₁, spring::CoiledSpring)

    b, N = spring.width, spring.coils

    L = √((x₁ - x₀)^2 + (y₁ - y₀)^2)  # Compute length of the spring
    a = (L - 2b) / (2N + 1)π  # Compute the radius of the cycloid

    # Obtain the different parts of a spring between (0, L) in the x axis

    tₗ = -3π/2:0.05:-π
    xₗ, yₗ = prolate_cycloid(a, b, tₗ)

    tₘ = -π:0.05:(2N-1)π
    xₘ, yₘ = prolate_cycloid(a, b, tₘ)

    tᵣ = (2N-1)π:0.05:(2N-1/2)π
    xᵣ, yᵣ = prolate_cycloid(a, b, tᵣ)

    xₚ, yₚ = vcat(xₗ, xₘ, xᵣ) .- (a * -3π / 2 - b), vcat(yₗ, yₘ, yᵣ)

    # Rotate and add the initial point to obtain the final position of the spring

    cosα, sinα = (x₁ - x₀) / L, (y₁ - y₀) / L

    x = @. x₀ + cosα * xₚ - sinα * yₚ
    y = @. y₀ + sinα * xₚ + cosα * yₚ

    return x, y
end
