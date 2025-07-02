using JLD2, GLMakie, Colors
include("..\\..\\src\\Common.jl")

path = abspath(@__DIR__, "..\\..")
begin
    data_1 = load("julia\\data\\Forced_response_fourier_28_10k\\idx_1\\fourier.jld2")
    freq_10k, U_10k = data_1["freq_fourier_f_2.8"], data_1["U_f_2.8"]

    data_2 = load("julia\\data\\Forced_response_fourier_28_100k\\idx_1\\fourier.jld2")
    freq_100k, U_100k = data_2["freq_fourier_f_2.8"], data_2["U_f_2.8"]

    fig_fourier =  Figure(size = (1280, 720))
    ax_fourier = Axis(fig_fourier[1, 1], title="Spectrum", xlabel="Frequency"; xgridvisible = true, ygridvisible = true, xgridcolor = :gray, ygridcolor = :gray)

    lines!(ax_fourier, freq_10k, U_10k, label="10k", color= :red)
    lines!(ax_fourier, freq_100k, U_100k, label="100k", color= :blue)
    xlims!(ax_fourier, 0.0, 10.0)

    display(GLMakie.Screen(), fig_fourier)
end

begin
    idx = "idx_1"
    case = "Forced_response_fourier"
    data_fourier = load("julia\\data\\" * case * "\\" * idx * "\\fourier.jld2")
    data_freq = load("julia\\data\\" * case * "\\" * idx * "\\frequency_sweep.jld2")

    freq_range = data_freq["freq_range"]
    dofs = data_freq["dofs"]
    dof_names = ["ψ₁", "ψ₂", "x_D", "y_D", "ψ_D"]

    fig_fourier_sweep = Vector{Figure}(undef, length(dofs))
    ax_fourier_sweep = Vector{Axis}(undef, length(dofs))

    idx_dof = 1
    for dof in dofs
        fig_fourier_sweep[idx_dof] = Figure(size = (1280, 720))
        ax_fourier_sweep[idx_dof] = Axis(fig_fourier_sweep[idx_dof][1, 1], title="Spectrum for DOF $(dof_names[dof])", xlabel="Frequency"; xgridvisible = true, ygridvisible = true, xgridcolor = :gray, ygridcolor = :gray)
        idx_dof += 1
    end

    freq_range = 1.5:0.1:4
    cols = distinguishable_colors(length(freq_range), [RGB(1,1,1), RGB(0,0,0)], dropseed=true)
    pcols = map(col -> (red(col), green(col), blue(col)), cols)

    freq_idx = 1
    for freq in freq_range
        dof_idx = 1
        for dof in dofs
            U = data_fourier["U_f_$(round(freq, sigdigits=4))"][dof_idx, :]
            freq_fourier = data_fourier["freq_fourier_f_$(round(freq, sigdigits=4))"]
            lines!(ax_fourier_sweep[dof_idx], freq_fourier, U, label="f = $(round(freq, sigdigits=4))", color=pcols[freq_idx])
            dof_idx += 1
        end
        freq_idx += 1
    end

    for idx in eachindex(dofs)
        axislegend(ax_fourier_sweep[idx]; labelsize = 12)
        xlims!(ax_fourier_sweep[idx], 0.0, 10.0) 
        display(GLMakie.Screen(), fig_fourier_sweep[idx])
    end
end