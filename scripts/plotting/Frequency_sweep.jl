using JLD2, GLMakie, Statistics
include("..\\..\\src\\Common.jl")

path = abspath(@__DIR__, "..\\..")
begin
    file_name = "Complete_frequency_sweep"
    mkpath("julia\\figures\\Complete_frequency_sweep")
    data_dir = "julia\\data\\"
    freq_file = "\\idx_1\\frequency_sweep.jld2"
    # Load files

    data_1 = load(data_dir * "Forced_response_detailed_3" * freq_file)
    freq_1, response_1 = data_1["freq_range"], data_1["amplitude"]

    data_2 = load(data_dir * "Forced_response_detailed_25" * freq_file)
    freq_2, response_2 = data_2["freq_range"], data_2["amplitude"]

    data_3 = load(data_dir * "Forced_response_detailed_26" * freq_file)
    freq_3, response_3 = data_3["freq_range"], data_3["amplitude"]

    data_4 = load(data_dir * "Forced_response_detailed_32" * freq_file)
    freq_4, response_4 = data_4["freq_range"], data_4["amplitude"]

    data_5 = load(data_dir * "Forced_response_detailed_262" * freq_file)
    freq_5, response_5 = data_5["freq_range"], data_5["amplitude"]

    data_6 = load(data_dir * "Forced_response_detailed_266" * freq_file)
    freq_6, response_6 = data_6["freq_range"], data_6["amplitude"]

    data_7 = load(data_dir * "Forced_response_detailed_3203" * freq_file)
    freq_7, response_7 = data_7["freq_range"], data_7["amplitude"]

    data_8 = load(data_dir * "Forced_response_detailed_3204" * freq_file)
    freq_8, response_8 = data_8["freq_range"], data_8["amplitude"]

    data_9 = load(data_dir * "Forced_response_high_res" * freq_file)
    freq_9, response_9 = data_9["freq_range"], data_9["amplitude"]

    data_10 = load(data_dir * "Forced_response_detailed_3_from_static" * freq_file)
    freq_10, response_10 = data_10["freq_range"], data_10["amplitude"]
    
    data_11 = load(data_dir * "Forced_response_detailed_26_from_static" * freq_file)
    freq_11, response_11 = data_11["freq_range"], data_11["amplitude"]
    
    data_12 = load(data_dir * "Forced_response_detailed_32_from_static" * freq_file)
    freq_12, response_12 = data_12["freq_range"], data_12["amplitude"]

    data_13 = load(data_dir * "Forced_response_detailed_28" * freq_file)
    freq_13, response_13 = data_13["freq_range"], data_13["amplitude"]

    data_14 = load(data_dir * "Forced_response_detailed_28_from_static" * freq_file)
    freq_14, response_14 = data_14["freq_range"], data_14["amplitude"]

    data_15 = load(data_dir * "Forced_response_detailed_28_high_cycles" * freq_file)
    freq_15, response_15 = data_15["freq_range"], data_15["amplitude"]

    # data_16 = load(data_dir * "Forced_response_detailed_28_high_cycles_from_static" * freq_file)
    # freq_16, response_16 = data_16["freq_range"], data_16["amplitude"]

    file_eig = "Modal_analysis_5DOF"
    data_eig = load("julia\\data\\eigenvalues\\" * file_eig * ".jld2")
    λ = data_eig["λ"]
end

begin
    fig_frequency_sweep = Figure(size = (1280, 720))
    ax_frequency_sweep = Axis(fig_frequency_sweep[1, 1], title="Frequency Sweep", xlabel="Frequency", ylabel="|max(ψ₁)-min(ψ₁)|"; 
                              xgridvisible = true, ygridvisible = true, xgridcolor = :gray, ygridcolor = :gray)
    # Evolution from previous solution
    scatter!(ax_frequency_sweep, freq_1, ifelse.(response_1 .== 0.0, NaN, response_1), color = :blue, markersize = 5)
    scatter!(ax_frequency_sweep, freq_2, ifelse.(response_2 .== 0.0, NaN, response_2), color = :blue, markersize = 5)
    scatter!(ax_frequency_sweep, freq_3, ifelse.(response_3 .== 0.0, NaN, response_3), color = :blue, markersize = 5)
    scatter!(ax_frequency_sweep, freq_4, ifelse.(response_4 .== 0.0, NaN, response_4), color = :blue, markersize = 5)
    scatter!(ax_frequency_sweep, freq_5, ifelse.(response_5 .== 0.0, NaN, response_5), color = :blue, markersize = 5)
    scatter!(ax_frequency_sweep, freq_6, ifelse.(response_6 .== 0.0, NaN, response_6), color = :blue, markersize = 5)
    scatter!(ax_frequency_sweep, freq_7, ifelse.(response_7 .== 0.0, NaN, response_7), color = :blue, markersize = 5)
    scatter!(ax_frequency_sweep, freq_8, ifelse.(response_8 .== 0.0, NaN, response_8), color = :blue, markersize = 5)
    scatter!(ax_frequency_sweep, freq_9, ifelse.(response_9 .== 0.0, NaN, response_9), color = :blue, markersize = 5)
    scatter!(ax_frequency_sweep, freq_13, ifelse.(response_13 .== 0.0, NaN, response_13), color = :blue, markersize = 5)

    # Evolution from static
    scatter!(ax_frequency_sweep, freq_10, ifelse.(response_10 .== 0.0, NaN, response_10), color = :red, markersize = 5)
    scatter!(ax_frequency_sweep, freq_11, ifelse.(response_11 .== 0.0, NaN, response_11), color = :red, markersize = 5)
    scatter!(ax_frequency_sweep, freq_12, ifelse.(response_12 .== 0.0, NaN, response_12), color = :red, markersize = 5)
    scatter!(ax_frequency_sweep, freq_14, ifelse.(response_14 .== 0.0, NaN, response_14), color = :red, markersize = 5)

    # Frequency of the eigenvalues
    for i in 1:5
        lines!(ax_frequency_sweep, abs(imag(λ[round(Int, 2*i-1)])).*[1.0, 1.0], [0.0, 0.004]; linestyle = :dash, color = :black)
    end

    xlims!(ax_frequency_sweep, 1.2, 4.0)
    ylims!(ax_frequency_sweep, 0.0, 0.0005)
    display(GLMakie.Screen(), fig_frequency_sweep)
end