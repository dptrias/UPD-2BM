using GLMakie, JLD2
include("..\\src\\StaticModels.jl")

path = abspath(@__DIR__, "..")
f_name = "Equilibrium_2DOF_test_storage"
begin
    # Parameters given by Zara et al.
    I_B = 6.348e-4 # kg m² / rad
    k_ψ = 4.0e3 # N m / rad
    c_ψ = 5.0e-2 # N m s / rad (5.0e-3)
    m_D = 1.0e-2 # kg
    k_Dx = 20.0 # N / m
    θ = pi/4 # rad
    FC = 160.0 # N
    μ = 0.25
    kₙ = 3.0e6 # N / m (3.0e5)
    kₜ = 3.0e6 # N / m (3.0e5)
    # Parameters not given by Zara et al.
    f = 0.1 # m (0.01)     0.1   
    b = 0.05 # m (0.005)   0.05 
    
    ω_ψ = sqrt(k_ψ / I_B)
    ζ_ψ = c_ψ / (2 * sqrt(k_ψ * I_B))
    ω_Dx = sqrt(k_Dx / m_D)
    fₓ = ω_Dx / ω_ψ
    α = I_B / (m_D * f^2)
    b = b / f
    FC = FC * f / k_ψ
    kₙ = kₙ * f^2 / k_ψ
    kₜ = kₜ * f^2 / k_ψ
end

pm = TBM_2DOF_Parameters(ζ_ψ = ζ_ψ, α = α, fₓ = fₓ, θ = θ, b = b, 
                         FC = FC, T_F = nothing,
                         friction_law = Linear(kₜ=kₜ, kₙ=kₙ))

T_range = range(-2e-3, 2e-3, length = 2000) # T values  

# Store solutions
solutions = Vector{NamedTuple{
    (:psi, :y_D, :T, :N, :xt, :xn, :wt), 
    Tuple{Float64, Float64, Float64, Float64, Float64, Float64, Float64}
}}()

# Loop over the range of T values
for T in T_range

    eq_pos = equilibrium_pos_2DOF(pm, T)
    if eq_pos === nothing
        println("[warning] T = $T : No solution found")
        continue
    elseif verbose
        println("[info] T = $T : Solution found")
    end

    sin_ψ = sin(eq_pos[1])
    cos_ψ = cos(eq_pos[1])

    xₜ = (sin_ψ - pm.b*(1 - cos_ψ))*cos(pm.θ) + (-(1 - cos_ψ) - pm.b*sin_ψ - eq_pos[2])*sin(pm.θ)
    xₙ = (sin_ψ - pm.b*(1 - cos_ψ))*sin(pm.θ) - (-(1 - cos_ψ) - pm.b*sin_ψ - eq_pos[2])*cos(pm.θ)
    N = pm.kₙ*xₙ

    if (abs(T) < abs(pm.μ*N) && N > 0.0) # Below the Coulomb limit
        wₜ = xₜ - T/pm.kₜ
        push!(solutions, (psi = eq_pos[1], y_D = eq_pos[2], T = T, N = N, xt = xₜ, xn = xₙ, wt = wₜ))
    end
end

# Extract the data for plotting
T_values = [sol.T for sol in solutions]
psi_values = [sol.psi for sol in solutions]
y_D_values = [sol.y_D for sol in solutions]
N_values = [sol.N for sol in solutions]
xt_values = [sol.xt for sol in solutions]
xn_values = [sol.xn for sol in solutions]
wt_values = [sol.wt for sol in solutions]

# Linear approximation 1
# linear_sol = [two_dof_equilibrium_linear!(pm, T) for T in T_values]
# psi_values_lin = [sol[1] for sol in linear_sol]
# y_D_values_lin = [sol[2] for sol in linear_sol]

# Linear approximation 2
sin_θ = sin(θ)
cos_θ = cos(θ)
A1 = -(FC*cos_θ^2*α*b*kₙ + FC*cos_θ*sin_θ*α*kₙ)/(kₙ*(b^2*fₓ^2 + 2*α)*cos_θ^2 + 2*cos_θ*sin_θ*b*fₓ^2*kₙ + fₓ^2*(sin_θ^2*kₙ + 1))
B1 = -((2*α*kₙ + fₓ^2)*cos_θ - sin_θ*b*fₓ^2)/(kₙ*(b^2*fₓ^2 + 2*α)*cos_θ^2 + 2*cos_θ*sin_θ*b*fₓ^2*kₙ + fₓ^2*(sin_θ^2*kₙ + 1))
A2 = 2*(FC*cos_θ^2*b^2*kₙ/2 + FC*cos_θ*sin_θ*b*kₙ + FC*sin_θ^2*kₙ/2 + FC/2)*α/(kₙ*(b^2*fₓ^2 + 2*α)*cos_θ^2 + 2*cos_θ*sin_θ*b*fₓ^2*kₙ + fₓ^2*(sin_θ^2*kₙ + 1))
B2 = 2*(cos_θ*b*kₙ + (kₙ + 1)*sin_θ)*α/(kₙ*(b^2*fₓ^2 + 2*α)*cos_θ^2 + 2*cos_θ*sin_θ*b*fₓ^2*kₙ + fₓ^2*(sin_θ^2*kₙ + 1))
psi_values_lin = A1 .+ B1.*T_values
y_D_values_lin = A2 .+ B2.*T_values

mkpath(joinpath(path, "figures\\" * f_name))

fig_psi = Figure(size = (1280, 720))
ax_psi = Axis(fig_psi[1, 1], title="ψ vs T", xlabel="T", ylabel="ψ"; xgridvisible = true, ygridvisible = true, xgridcolor = :gray, ygridcolor = :gray)
lines!(T_values, psi_values, linewidth = 2, label = "Nonlinear")
lines!(T_values, psi_values_lin, label = "Linear approx.", color = :red, linestyle = :dash)
axislegend(ax_psi)
save(joinpath(path, "figures\\" * f_name * "\\psi.png"), fig_psi)

fig_psi_diff = Figure(size = (1280, 720))
ax_psi_diff = Axis(fig_psi_diff[1, 1], title="Difference between linear and non-linear models", xlabel="T", ylabel="(ψ - ψₗᵢₙ)"; xgridvisible = true, ygridvisible = true, xgridcolor = :gray, ygridcolor = :gray)
lines!(T_values, psi_values .- psi_values_lin)
save(joinpath(path, "figures\\" * f_name * "\\psi_diff.png"), fig_psi_diff)

fig_y_D = Figure(size = (1280, 720))
ax_y_D = Axis(fig_y_D[1, 1], title="y_D vs T", xlabel="T", ylabel="y_D"; xgridvisible = true, ygridvisible = true, xgridcolor = :gray, ygridcolor = :gray)
lines!(T_values, y_D_values; linewidth = 2, label = "Nonlinear")
lines!(T_values, y_D_values_lin; label = "Linear approx.", color = :red, linestyle = :dash)
axislegend(ax_y_D; position = :rb)
save(joinpath(path, "figures\\" * f_name * "\\y_D.png"), fig_y_D)
fig_y_D_diff = Figure(size = (1280, 720))
ax_y_D_diff = Axis(fig_y_D_diff[1, 1], title="Difference between linear and non-linear models", xlabel="T", ylabel="(y_D - y_Dₗᵢₙ)"; xgridvisible = true, ygridvisible = true, xgridcolor = :gray, ygridcolor = :gray)
lines!(T_values, y_D_values .- y_D_values_lin)
save(joinpath(path, "figures\\" * f_name * "\\y_D_diff.png"), fig_y_D_diff)

fig_N = Figure(size = (1280, 720))
ax_N = Axis(fig_N[1, 1], title="N vs T", xlabel="T", ylabel="N"; xgridvisible = true, ygridvisible = true, xgridcolor = :gray, ygridcolor = :gray)
lines!(T_values, N_values)
save(joinpath(path, "figures\\" * f_name * "\\N.png"), fig_N)

fig_xt = Figure(size = (1280, 720))
ax_xt = Axis(fig_xt[1, 1], title="xₜ vs T", xlabel="T", ylabel="xₜ"; xgridvisible = true, ygridvisible = true, xgridcolor = :gray, ygridcolor = :gray)
lines!(T_values, xt_values)
save(joinpath(path, "figures\\" * f_name * "\\x_t.png"), fig_xt)

fig_xn = Figure(size = (1280, 720))
ax_xn = Axis(fig_xn[1, 1], title="xₙ vs T", xlabel="T", ylabel="xₙ"; xgridvisible = true, ygridvisible = true, xgridcolor = :gray, ygridcolor = :gray)
lines!(T_values, xn_values)
save(joinpath(path, "figures\\" * f_name * "\\x_n.png"), fig_xn)

fig_wt = Figure(size = (1280, 720))
ax_wt = Axis(fig_wt[1, 1], title="wₜ vs T", xlabel="T", ylabel="wₜ"; xgridvisible = true, ygridvisible = true, xgridcolor = :gray, ygridcolor = :gray)
lines!(T_values, wt_values)
save(joinpath(path, "figures\\" * f_name * "\\w_t.png"), fig_wt)

fig_T = Figure(size = (1280, 720))
ax_T = Axis(fig_T[1, 1], title="T vs T", xlabel="T", ylabel="T"; xgridvisible = true, ygridvisible = true, xgridcolor = :gray, ygridcolor = :gray)
lines!(T_values, T_values, label = "T")
lines!(T_values, pm.μ*N_values, label = "±μN", color = :red, linestyle = :dash) 
lines!(T_values, -pm.μ*N_values, color = :red, linestyle = :dash) 
axislegend(ax_T)
save(joinpath(path, "figures\\" * f_name * "\\T.png"), fig_T)

f_path = joinpath(path, "data\\static\\" * f_name * ".jld2")
print("Save results? [Y/n]: ")
response = readline() |> strip
if response == "Y" || response == "y"
    if isfile(f_path)
        println("[warning] File $f_path already exists. Choose a different name or delete the file")
    else
        mkpath(joinpath(path, "data\\static"))
        @save f_path solutions pm
        println("[info] Results saved to $f_path")
    end
else
    println("[info] Results not saved")
end