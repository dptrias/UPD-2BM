using JLD2, Printf, Julianim, MathTeXEngine
include("..\\..\\src\\Common.jl")

path = abspath(@__DIR__, "..\\..")
function check_existing_file(freq)
    if freq == 0.0
        return true
    end
    file_path = abspath(path, "..", "data\\Final_response" * freetext_data * "\\freq_$(replace(string(freq), "." => "-"))\\statistics.jld2")
    if isfile(file_path)
        return true
    else
        println("[info] No data found for frequency $freq Skipping...")
        return false
    end
end

begin
    frequency_range = vcat([0.0, 0.0], collect(0.6:0.05:1.6))
    freetext_data = "_final_less_Me"
    eq_pos_file = "Equilibrium_5DOF_final_even_lower_def_reduced"
    filter!(x -> check_existing_file(x), frequency_range)
end

# Video parameters
begin
    FPS = 1
    time = length(frequency_range)
    frames = FPS * time
    ts = range(1, stop=length(frequency_range), length=frames)
    freetext = "_final_less_Me"
    
    idx_freq = Observable(1)
    step!(idx_freq) = idx_freq[] += 1
end

function load_static_data!(eq_pos_file)
    # Get data from static equilibrium
    static_data = jldopen(abspath(path, "..\\data\\static\\" * eq_pos_file * ".jld2"), "r")
    sol_static = static_data["solutions"]
    pm_static = static_data["pm"]
    close(static_data)

    # Extract the data from static solution
    global T1_static = [sol.T1 for sol in sol_static]
    global T2_static = [sol.T2 for sol in sol_static]
    global wt1_static = [sol.wt1 for sol in sol_static]
    global wt2_static = [sol.wt2 for sol in sol_static]
    global n_static = length(T1_static)
    global u_static = zeros(5, n_static)
    u_static[1, :] = [sol.psi1 for sol in sol_static]
    u_static[2, :] = [sol.psi2 for sol in sol_static]
    u_static[3, :] = [sol.x_D for sol in sol_static]
    u_static[4, :] = [sol.y_D for sol in sol_static]
    u_static[5, :] = [sol.psi_D for sol in sol_static]
    
    global pm = pm_static

    println("[info] Static equilibrium file loaded: $eq_pos_file")
end

# Get modal analysis data
begin
    load_static_data!(eq_pos_file)

    set_publication_theme!()
    update_theme!(
    fonts = Attributes(
            :bold => texfont(:bold),
            :bolditalic => texfont(:bolditalic),
            :italic => texfont(:italic),
            :regular => texfont(:regular)
        )
    )

    fig = Figure()
    ax = Axis(fig[1, 1], 
        xlabel = L"T_1", 
        ylabel = L"T_2", 
        aspect = 1
    )
    Δ_T1 = (maximum(T1_static) - minimum(T1_static))/20
    Δ_T2 = (maximum(T2_static) - minimum(T2_static))/20
    limits!(ax, 
        minimum(T1_static) - Δ_T1, maximum(T1_static) + Δ_T1, 
        minimum(T2_static) - Δ_T2, maximum(T2_static) + Δ_T2
    )
    # Fondo de contorno gris con valores constantes (por ejemplo 1)
    contourf!(ax, 
        T1_static, 
        T2_static, 
        colormap = :grays, 
        ones(n_static)
    )

    idx_freq[] = 1
    freq = @lift frequency_range[$(idx_freq)]
    T_1_final = @lift begin
        if freq[] == 0.0
            T1_static
        else
            filepath = abspath(path, "..", "data", "Final_response" * freetext_data, "freq_$(replace(string($freq), "." => "-"))", "statistics.jld2")
            file = jldopen(filepath, "r")
            data = file["T_1_final"]
            close(file)
            data
        end
    end

    T_2_final = @lift begin
        if freq[] == 0.0
            T2_static
        else
            filepath = abspath(path, "..", "data", "Final_response" * freetext_data, "freq_$(replace(string($freq), "." => "-"))", "statistics.jld2")
            file = jldopen(filepath, "r")
            data = file["T_2_final"]
            close(file)
            data
        end
    end

    # Etiqueta con frecuencia actual
    textlabel!(ax,
        Point2f(0.0005, -0.0002), 
        text = @lift(LaTeXString(@sprintf("\$ f= %.2f \$", $freq))), 
        fontsize = 30, 
        text_align = (:center, :center)
    )

    # Puntos sobre la figura
    scatter!(ax, 
        @lift($T_1_final .+ 0), 
        @lift($T_2_final .+ 0), 
        color = RED,
        marker = :circle,
    )

    # Preparación para guardado del vídeo
    format = "gif"
    dirpath = abspath(path, "..", "figures", "Final_response" * freetext)
    if !isdir(dirpath)
        mkpath(dirpath)
    end
    file = joinpath(dirpath, "freq_sweep." * format)

    # Guardar el vídeo
    record(fig, file, ts[1:end-1]; framerate = FPS) do i
        println("[info] Frame: $i")
        step!(idx_freq)
        println("[info] Frequency: $(idx_freq[])")
    end
end
