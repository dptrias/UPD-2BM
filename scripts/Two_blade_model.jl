using Dates, GLMakie, JLD2, NLsolve, MAT
include("..\\src\\DynamicModels.jl")

path = abspath(@__DIR__, "..")
## MODEL
begin
    # Parameters given by Zara et al.
    I_B = 6.348e-4 # kg m² / rad
    k_ψ = 4.0e3 # N m / rad
    c_ψ = 5.0e-2 # N m s / rad (5.0e-3)
    I_D = 4.558e-8 # kg m²
    m_D = 1.0e-2 # kg
    k_Dx = 20.0 # N / m
    k_Dψ = 20.0 # N m / rad 
    θ = pi/4 # rad
    Mₑ = 0.0 # N m (1)
    FC = 160.0 # N
    μ = 0.5
    kₙ = 3.0e6 # N / m (3.0e5)
    kₜ = 3.0e6 # N / m (3.0e5)
    # Parameters not given by Zara et al.
    f = 0.1 # m (0.01)     0.1   
    b = 0.05 # m (0.005)   0.05 
    dₙ = 0.007 # m (0.003)  0.007
    dₜ = 0.001 # m (0.0005) 0.001
    
    ω_ψ = sqrt(k_ψ / I_B)
    ζ_ψ = c_ψ / (2 * sqrt(k_ψ * I_B))
    ω_Dx = sqrt(k_Dx / m_D)
    ω_Dψ = sqrt(k_Dψ / I_D)
    fₓ = ω_Dx / ω_ψ
    f_ψ = ω_Dψ / ω_ψ
    α = I_B / (m_D * f^2)
    β = I_B / I_D
    b = b / f
    dₙ = dₙ / f
    dₜ = dₜ / f
    FC = FC * f / k_ψ
    Mₑ = Mₑ / k_ψ
    kₙ = kₙ * f^2 / k_ψ
    kₜ = kₜ * f^2 / k_ψ
end

## FILE MANAGEMENT
function file_name(; model, Δt, tf, method, friction_type, feq = nothing,  freetext = nothing)
    exponent = floor(Int, log10(Δt)) 
    mantissa = round(Δt / 10.0^exponent)  # Round to nearest integer

    # Format time step
    formatted_dt = "$(Int(mantissa))e$(exponent)"

    f_name = model * "_"
    if feq !== nothing
        f_name *= "f$(replace(string(feq), "." => ""))_"
    end
    f_name *= "tf$(floor(Int, tf))_" *
              "dt$(formatted_dt)_" *
              method * "_" *
              friction_type
    if freetext !== nothing
        f_name *=  "_" * freetext
    end
    
    return f_name
end

## TIME INTEGRATION
start_time = time()
n_steps = 2_000_001
tf = 200.0
t = range(0.0, tf, length = n_steps)
Δt = t[2] - t[1]
u = zeros(10, n_steps)

feq = 0.0
pm = TBM_5DOF_Parameters(feq = feq, ζ_ψ = ζ_ψ, α = α, β = β, fₓ = fₓ, f_ψ = f_ψ, θ = θ, b = b, dₙ = dₙ, dₜ = dₜ, Mₑ = 0.0, FC = FC,
                         friction_law_1 = Coulomb(kₙ = kₙ, xₙ₀ = 0.0, kₜ = kₜ, μ = μ), # Linear(kₙ= kₙ, kₜ = kₜ),
                         friction_law_2 = Coulomb(kₙ = kₙ, xₙ₀ = 0.0, kₜ = kₜ, μ = μ), # Linear(kₙ= kₙ, kₜ = kₜ),
                         T₁ = zeros(n_steps), N₁ = zeros(n_steps), T₂ = zeros(n_steps), N₂ = zeros(n_steps),
                         xₜ₁ = zeros(n_steps), xₙ₁ = zeros(n_steps), xₜ₂ = zeros(n_steps), xₙ₂ = zeros(n_steps), index = Ref(1)) 


for i in 1:(length(t) - 1)
    u[:, i+1] = rk4(five_dof_model, u[:, i], pm, t[i], Δt)
end


elapsed_time = time() - start_time

# fig_time = Figure()
# ax_t = Axis(fig_time[1, 1], title="time", xlabel="i", ylabel="t")
# scatterlines!(t, linestyle = :dash, markersize=4)
# ax_Δt = Axis(fig_time[1, 2], title="Timestep vs time", xlabel="t", ylabel="dt")
# scatterlines!(t[2:end], diff(t), linestyle = :dash, markersize=4)
# display(GLMakie.Screen(), fig_time)

## PLOT RESULTS
print("Plot results? [Y/n]: ")
response = readline() |> strip
if response == "Y" || response == "y"
    fig = Figure()
    ax1 = Axis(fig[1, 1], title="ψ₁ vs time", xlabel="t", ylabel="ψ₁")
    lines!(t, u[1, :])
    ax2 = Axis(fig[1, 2], title="ψ₂ vs time", xlabel="t", ylabel="ψ₂")
    lines!(t, u[2, :])
    ax3 = Axis(fig[2, 1], title="x_D vs time", xlabel="t", ylabel="x_D")
    lines!(t, u[3, :])    
    ax4 = Axis(fig[2, 2], title="y_D vs time", xlabel="t", ylabel="y_D")
    lines!(t, u[4, :])

    ax5 = Axis(fig[1:2, 3], title="ψ_D vs time", xlabel="t", ylabel="ψ_D")
    lines!(t, u[5, :])

    display(GLMakie.Screen(), fig)
else
    println("[info] Results not plotted.")
end

## SAVE RESULTS
f_name = file_name(model = "TBM", Δt = Δt, tf = tf, method = "rk4",
                   friction_type = friction_type(pm.friction_law_1), 
                   feq = nothing, 
                   freetext = "matlab_verification") * ".jld2"
mkpath(joinpath(path, "data"))
f_path = joinpath(path, "data", f_name)
info = "Test to verify with Matlab with an intial displacement of 0.001 rad at psi_D. Elapsed time: $elapsed_time s."
print("Save results? [Y/n]: ")
response = readline() |> strip
if response == "Y" || response == "y"
    if isfile(f_path)
        println("[warning] File $f_path already exists. Choose a different name or delete the file")
    else
        @save f_path u t f pm
        println("[info] Results saved in $f_path")
        # Save info in case_info file
        open(joinpath(path, "data/case_info.txt"), "a") do io
            println(io, "File name: $f_name")
            println(io, "Date and time: $(Dates.now())")
            println(io, "Info: $info")
            println(io, "------------------------------------")
        end

        # Save matlab file
        print("Save results to matlab file? [Y/n]: ")
        response = readline() |> strip
        if response == "Y" || response == "y"    
            mat_file = joinpath(path, "data", replace(f_name, ".jld2" => ".mat"))
            matwrite(mat_file, Dict("u" => collect(u),
                                    "t" => collect(t),
                                    "feq" => pm.feq,
                                    "zeta" => pm.ζ_ψ,
                                    "alpha" => pm.α,
                                    "beta" => pm.β,
                                    "f_x" => pm.fₓ,
                                    "f_psi" => pm.f_ψ,
                                    "theta" => pm.θ,
                                    "b" => pm.b,
                                    "d_n" => pm.dₙ,
                                    "d_t" => pm.dₜ,
                                    "M_e" => pm.Mₑ,
                                    "FC" => pm.FC,
                                    "T_1" => collect(pm.T₁),
                                    "N_1" => collect(pm.N₁),
                                    "T_2" => collect(pm.T₂),
                                    "N_2" => collect(pm.N₂),
                                    "x_t1" => collect(pm.xₜ₁),
                                    "x_n1" => collect(pm.xₙ₁),
                                    "x_t2" => collect(pm.xₜ₂),
                                    "x_n2" => collect(pm.xₙ₂)
            ))
            println("[info] Results also saved to matlab file in $mat_file")
        else
            println("[info] Results not saved to matlab file.")
        end
    end
else
    println("[info] Results not saved.")
end