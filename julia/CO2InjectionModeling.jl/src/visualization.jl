using CairoMakie
using Statistics
using SurfaceWaterIntegratedModeling: TrapStructure, numtraps, trap_states_at_timepoints, SpillEvent

export animate_multi_layer_filling, animate_multi_layer_saturation
export plot_layer_volumes_timeseries, plot_layer_fractions_timeseries


"""
    animate_multi_layer_filling(layers, seqs, domain; kwargs...)

Create an animated bird's eye view of CO2 filling in all layers simultaneously.
Displays layers in a 3x3 grid layout (for 9 layers).

Parameters:
- `layers`: Vector of Layer structs
- `seqs`: Vector{Vector{SpillEvent}} from fill_layers
- `domain`: Domain3D struct
- `output_file`: Path to save animation (default: "multi_layer_filling.gif")
- `num_frames`: Number of frames in animation (default: 30)
- `start_time`: Start time for animation (default: 0.0)
- `end_time`: End time for animation (default: auto-detect from seqs)
- `fps`: Frames per second (default: 2)
- `colormap`: Colormap for heights (default: :thermal)
- `max_CO2_height`: Maximum CO2 height for colorscale (default: 20.0)
"""
function animate_multi_layer_filling(
    layers::Vector{Layer},
    seqs::Vector{Vector{SpillEvent}},
    domain::Domain3D;
    output_file::String = "multi_layer_filling.gif",
    num_frames::Int = 30,
    start_time::Float64 = 0.0,
    end_time::Union{Float64, Nothing} = nothing,
    fps::Int = 2,
    colormap::Symbol = :thermal,
    max_CO2_height::Float64 = 20.0
)
    n_layers = length(layers)
    @assert n_layers == length(seqs) "Number of layers must match number of sequences"

    # Determine grid layout (3x3 for 9 layers, or calculate appropriate layout)
    n_cols = ceil(Int, sqrt(n_layers))
    n_rows = ceil(Int, n_layers / n_cols)

    # Determine end_time from all sequences if not provided
    if isnothing(end_time)
        max_times = Float64[]
        for seq in seqs
            if !isempty(seq)
                max_t = maximum(se.timestamp for se in seq)
                if isfinite(max_t)
                    push!(max_times, max_t)
                elseif length(seq) > 1
                    push!(max_times, seq[end-1].timestamp + 1.0)
                end
            end
        end
        end_time = isempty(max_times) ? 15.0 : maximum(max_times)
    end

    # Generate timepoints for animation
    timepoints = collect(range(start_time, stop=end_time, length=num_frames))

    # Precompute data for each layer
    println("Computing trap states for $(n_layers) layers × $(num_frames) frames...")

    layer_data = []
    for (layer_idx, layer) in enumerate(layers)
        tstruct = layer.trap_structure
        num_traps = numtraps(tstruct)
        seq = seqs[layer_idx]

        # Skip empty layers
        if isempty(seq) || num_traps == 0
            push!(layer_data, nothing)
            continue
        end

        # Get trap states at each timepoint
        tstates = trap_states_at_timepoints(tstruct, seq, timepoints; verbose=false)

        # Compute z_vol_tables for height conversion
        z_vol_tables = SurfaceWaterIntegratedModeling._compute_z_vol_tables(tstruct)

        push!(layer_data, (
            tstruct = tstruct,
            num_traps = num_traps,
            tstates = tstates,
            z_vol_tables = z_vol_tables,
            pad = layer.boundary_padding,
            name = layer.name
        ))
    end

    # Get grid size from first layer (assume all layers have same size after unpadding)
    first_valid_idx = findfirst(!isnothing, layer_data)
    if isnothing(first_valid_idx)
        error("No valid layers found")
    end

    first_layer = layer_data[first_valid_idx]
    pad = first_layer.pad
    topo_size = size(first_layer.tstruct.topography)
    nx_padded, ny_padded = topo_size
    nx = nx_padded - 2 * pad
    ny = ny_padded - 2 * pad

    # Set up figure with grid layout
    fig = Figure(size = (400 * n_cols, 350 * n_rows + 100))

    x_coords = range(0, nx * domain.dx, length=nx)
    y_coords = range(0, ny * domain.dy, length=ny)

    # Create observables and axes for each layer
    height_observables = []
    axes = []

    for layer_idx in 1:n_layers
        row = div(layer_idx - 1, n_cols) + 1
        col = mod(layer_idx - 1, n_cols) + 1

        height_data = Observable(zeros(Float64, nx, ny))
        push!(height_observables, height_data)

        layer_name = isnothing(layer_data[layer_idx]) ? "Layer $layer_idx" : layer_data[layer_idx].name

        ax = Axis(fig[row, col],
                  title = layer_name,
                  xlabel = col == 1 ? "X (m)" : "",
                  ylabel = row == n_rows ? "Y (m)" : "",
                  aspect = DataAspect(),
                  xticklabelsvisible = (row == n_rows),
                  yticklabelsvisible = (col == 1))

        hm = heatmap!(ax, x_coords, y_coords, height_data,
                      colormap = colormap,
                      colorrange = (0.0, max_CO2_height))

        push!(axes, ax)
    end

    # Add a shared colorbar
    Colorbar(fig[:, n_cols + 1], colormap = colormap, colorrange = (0.0, max_CO2_height),
             label = "CO2 Height (m)")

    # Add overall title with time
    time_label = Observable("Time: 0.0 years")
    Label(fig[0, :], time_label, fontsize = 20, font = :bold)

    println("Creating animation...")

    # Create the animation
    record(fig, output_file, eachindex(timepoints); framerate=fps) do frame_idx
        tp = timepoints[frame_idx]

        total_volume_all_layers = 0.0
        max_height_all_layers = 0.0

        for layer_idx in 1:n_layers
            ld = layer_data[layer_idx]

            if isnothing(ld)
                # Empty layer - just zeros
                height_observables[layer_idx][] = zeros(Float64, nx, ny)
                continue
            end

            filled, volumes, _ = ld.tstates[frame_idx]

            # Create height map (on padded grid, then remove padding)
            height_map_padded = zeros(Float64, nx_padded, ny_padded)

            max_height = 0.0
            for trap_id in 1:ld.num_traps
                volume = volumes[trap_id]
                if volume > 0.0 || filled[trap_id]
                    z_vol_table = ld.z_vol_tables[trap_id]
                    height = volume_to_height(volume, trap_id, z_vol_table, ld.tstruct)
                    max_height = max(max_height, height)

                    footprint = ld.tstruct.footprints[trap_id]
                    for idx in footprint
                        height_map_padded[idx] = max(height_map_padded[idx], height)
                    end
                end
            end

            # Remove padding
            height_map = height_map_padded[pad+1:end-pad, pad+1:end-pad]

            # Update observable
            height_observables[layer_idx][] = height_map

            total_volume_all_layers += sum(volumes)
            max_height_all_layers = max(max_height_all_layers, max_height)
        end

        # Update time label
        time_label[] = "Time: $(round(tp, digits=2)) years | Max height: $(round(max_height_all_layers, digits=1)) m"

        if frame_idx % 10 == 0 || frame_idx == length(timepoints)
            println("  Frame $(frame_idx)/$(length(timepoints))")
        end
    end

    println("Animation saved to: $(output_file)")
    return nothing
end


"""
    animate_multi_layer_saturation(layers, seqs, leakage_states, domain; kwargs...)

Create an animated bird's eye view of CO2 saturation in all layers simultaneously.
Displays layers in a grid layout (3x3 for 9 layers).

This visualization shows CO2 saturation (volume/capacity) instead of height,
which makes residual drainage visible: after leakage starts, the saturation
in draining traps decreases from 1.0 down to the residual saturation level.

Parameters:
- `layers`: Vector of Layer structs
- `seqs`: Vector{Vector{SpillEvent}} from fill_layers
- `leakage_states`: Vector{LeakageState} from fill_layers
- `domain`: Domain3D struct
- `output_file`: Path to save animation (default: "multi_layer_saturation.gif")
- `num_frames`: Number of frames in animation (default: 30)
- `start_time`: Start time for animation (default: 0.0)
- `end_time`: End time for animation (default: auto-detect from seqs)
- `fps`: Frames per second (default: 2)
- `colormap`: Colormap for saturation (default: :viridis)
"""
function animate_multi_layer_saturation(
    layers::Vector{Layer},
    seqs::Vector{Vector{SpillEvent}},
    leakage_states::Vector{LeakageState},
    domain::Domain3D,
    reservoir_properties::Vector{ReservoirProperties};
    output_file::String = "multi_layer_saturation.gif",
    num_frames::Int = 30,
    start_time::Float64 = 0.0,
    end_time::Union{Float64, Nothing} = nothing,
    fps::Int = 2,
    colormap::Symbol = :viridis
)
    n_layers = length(layers)
    @assert n_layers == length(seqs) "Number of layers must match number of sequences"
    @assert n_layers == length(leakage_states) "Number of layers must match number of leakage states"

    # Determine grid layout (3x3 for 9 layers, or calculate appropriate layout)
    n_cols = ceil(Int, sqrt(n_layers))
    n_rows = ceil(Int, n_layers / n_cols)

    # Determine end_time from all sequences if not provided
    if isnothing(end_time)
        max_times = Float64[]
        for seq in seqs
            if !isempty(seq)
                max_t = maximum(se.timestamp for se in seq)
                if isfinite(max_t)
                    push!(max_times, max_t)
                elseif length(seq) > 1
                    push!(max_times, seq[end-1].timestamp + 1.0)
                end
            end
        end
        end_time = isempty(max_times) ? 15.0 : maximum(max_times)
    end

    # Generate timepoints for animation
    timepoints = collect(range(start_time, stop=end_time, length=num_frames))

    # Precompute data for each layer
    println("Computing trap states for $(n_layers) layers × $(num_frames) frames...")

    layer_data = []
    for (layer_idx, layer) in enumerate(layers)
        tstruct = layer.trap_structure
        num_traps = numtraps(tstruct)
        seq = seqs[layer_idx]

        # Skip empty layers
        if isempty(seq) || num_traps == 0
            push!(layer_data, nothing)
            continue
        end

        # Get trap states at each timepoint
        tstates = trap_states_at_timepoints(tstruct, seq, timepoints; verbose=false)

        # Compute trap capacities for saturation calculation
        trap_capacities = [tstruct.trapvolumes[i] - tstruct.subvolumes[i] for i in 1:num_traps]

        push!(layer_data, (
            tstruct = tstruct,
            num_traps = num_traps,
            tstates = tstates,
            trap_capacities = trap_capacities,
            pad = layer.boundary_padding,
            name = layer.name,
            leakage_state = leakage_states[layer_idx]
        ))
    end

    # Get grid size from first layer (assume all layers have same size after unpadding)
    first_valid_idx = findfirst(!isnothing, layer_data)
    if isnothing(first_valid_idx)
        error("No valid layers found")
    end

    first_layer = layer_data[first_valid_idx]
    pad = first_layer.pad
    topo_size = size(first_layer.tstruct.topography)
    nx_padded, ny_padded = topo_size
    nx = nx_padded - 2 * pad
    ny = ny_padded - 2 * pad

    # Set up figure with grid layout
    fig = Figure(size = (400 * n_cols, 350 * n_rows + 100))

    x_coords = range(0, nx * domain.dx, length=nx)
    y_coords = range(0, ny * domain.dy, length=ny)

    # Create observables and axes for each layer
    saturation_observables = []
    axes = []

    for layer_idx in 1:n_layers
        row = div(layer_idx - 1, n_cols) + 1
        col = mod(layer_idx - 1, n_cols) + 1

        saturation_data = Observable(zeros(Float64, nx, ny))
        push!(saturation_observables, saturation_data)

        layer_name = isnothing(layer_data[layer_idx]) ? "Layer $layer_idx" : layer_data[layer_idx].name

        ax = Axis(fig[row, col],
                  title = layer_name,
                  xlabel = col == 1 ? "X (m)" : "",
                  ylabel = row == n_rows ? "Y (m)" : "",
                  aspect = DataAspect(),
                  xticklabelsvisible = (row == n_rows),
                  yticklabelsvisible = (col == 1))

        hm = heatmap!(ax, x_coords, y_coords, saturation_data,
                      colormap = colormap,
                      colorrange = (0.0, 1.0))

        push!(axes, ax)
    end

    # Add a shared colorbar
    Colorbar(fig[:, n_cols + 1], colormap = colormap, colorrange = (0.0, 1.0),
             label = "CO2 Saturation")

    # Add overall title with time
    time_label = Observable("Time: 0.0 years")
    Label(fig[0, :], time_label, fontsize = 20, font = :bold)

    println("Creating animation...")

    # Create the animation
    record(fig, output_file, eachindex(timepoints); framerate=fps) do frame_idx
        tp = timepoints[frame_idx]

        total_draining_all_layers = 0
        max_saturation_all_layers = 0.0

        for layer_idx in 1:n_layers
            ld = layer_data[layer_idx]

            if isnothing(ld)
                # Empty layer - just zeros
                saturation_observables[layer_idx][] = zeros(Float64, nx, ny)
                continue
            end

            filled, volumes, _ = ld.tstates[frame_idx]
            leakage_state = ld.leakage_state

            # Create saturation map (on padded grid, then remove padding)
            saturation_map_padded = zeros(Float64, nx_padded, ny_padded)

            max_saturation = 0.0
            num_draining = 0

            for trap_id in 1:ld.num_traps
                volume = volumes[trap_id]

                # For draining traps, compute drainage-adjusted volume
                if leakage_state.draining[trap_id]
                    drained_vol = compute_volume_with_drainage(trap_id, tp, leakage_state)
                    if !isnothing(drained_vol)
                        volume = drained_vol
                    end
                    num_draining += 1
                end

                # Compute saturation = volume / capacity
                capacity = ld.trap_capacities[trap_id]
                if capacity > 0.0
                    saturation = min(1.0, volume / capacity)
                else
                    saturation = 0.0
                end

                max_saturation = max(max_saturation, saturation)

                # Fill the trap footprint with this saturation
                if saturation > 0.0
                    footprint = ld.tstruct.footprints[trap_id]
                    for idx in footprint
                        saturation_map_padded[idx] = max(saturation_map_padded[idx], saturation)
                    end
                end
            end

            # Remove padding
            saturation_map = saturation_map_padded[pad+1:end-pad, pad+1:end-pad]

            # Rescale the saturation_map to account for the sand_irreducible_water_saturation
            saturation_map *= 1-reservoir_properties[layer_idx].sand_irreducible_water_saturation

            # Update observable
            saturation_observables[layer_idx][] = saturation_map

            total_draining_all_layers += num_draining
            max_saturation_all_layers = max(max_saturation_all_layers, max_saturation)
        end

        # Update time label
        time_label[] = "Time: $(round(tp, digits=2)) years"

        if frame_idx % 10 == 0 || frame_idx == length(timepoints)
            println("  Frame $(frame_idx)/$(length(timepoints))")
        end
    end

    println("Animation saved to: $(output_file)")
    return nothing
end


"""
    plot_layer_volumes_timeseries(snapshots; kwargs...)

Plot CO2 volumes in each layer as a function of time using subplots.

Parameters:
- `snapshots`: Vector of ReservoirSnapshot from generate_reservoir_snapshots
- `output_file`: Path to save figure (default: nothing, returns figure)
- `title`: Overall title (default: "CO2 Volume by Layer")
- `colormap`: Colormap for layers (default: :tab10)
"""
function plot_layer_volumes_timeseries(
    snapshots::Vector{ReservoirSnapshot};
    output_file::Union{String, Nothing} = nothing,
    title::String = "CO2 Volume by Layer",
    colormap::Symbol = :tab10
)
    n_snapshots = length(snapshots)
    n_layers = length(snapshots[1].layer_snapshots)

    # Extract time series data
    times = [s.timestamp for s in snapshots]
    volumes_by_layer = [
        [s.stored_by_layer_m3[i] / 1e6 for s in snapshots]  # Convert to M m³
        for i in 1:n_layers
    ]

    # Get layer names
    layer_names = [snapshots[1].layer_snapshots[i].layer_name for i in 1:n_layers]

    # Fixed 3x3 grid for 9 layers
    n_cols = 3
    n_rows = 3

    # Create figure with clean styling
    fig = Figure(
        size = (900, 800),
        backgroundcolor = :white
    )

    # Generate colors from colormap
    cmap = cgrad(colormap, n_layers, categorical = true)

    # Find global y-axis maximum for consistent scaling
    max_vol = maximum(maximum(v) for v in volumes_by_layer)
    y_max = max_vol * 1.15  # Add 15% padding

    # Create subplots for each layer
    for layer_idx in 1:n_layers
        row = div(layer_idx - 1, n_cols) + 1
        col = mod(layer_idx - 1, n_cols) + 1

        # Determine which labels to show
        show_xlabel = row == n_rows
        show_ylabel = col == 1

        ax = Axis(fig[row, col],
                  title = layer_names[layer_idx],
                  titlesize = 14,
                  titlefont = :bold,
                  xlabel = show_xlabel ? "Time (years)" : "",
                  ylabel = show_ylabel ? "Volume (M m³)" : "",
                  xlabelsize = 12,
                  ylabelsize = 12,
                  xticklabelsize = 10,
                  yticklabelsize = 10,
                  xticklabelsvisible = show_xlabel,
                  yticklabelsvisible = show_ylabel,
                  xgridvisible = true,
                  ygridvisible = true,
                  xgridcolor = (:gray, 0.3),
                  ygridcolor = (:gray, 0.3),
                  xgridstyle = :dot,
                  ygridstyle = :dot,
                  spinewidth = 1,
                  xtickwidth = 1,
                  ytickwidth = 1)

        layer_color = cmap[layer_idx]

        # Fill area under curve first (so line is on top)
        band!(ax, times, zeros(n_snapshots), volumes_by_layer[layer_idx],
              color = (layer_color, 0.4))

        # Plot volume line
        lines!(ax, times, volumes_by_layer[layer_idx],
               color = layer_color,
               linewidth = 2.5)

        # Add markers at data points
        scatter!(ax, times, volumes_by_layer[layer_idx],
                 color = layer_color,
                 markersize = 4,
                 strokewidth = 0)

        ylims!(ax, 0, y_max)
        xlims!(ax, times[1], times[end])
    end

    # Add overall title
    Label(fig[0, :], title, fontsize = 18, font = :bold, padding = (0, 0, 10, 0))

    # Adjust spacing
    colgap!(fig.layout, 15)
    rowgap!(fig.layout, 15)

    if !isnothing(output_file)
        save(output_file, fig)
        println("Figure saved to: $(output_file)")
    end

    return fig
end


"""
    plot_layer_fractions_timeseries(snapshots; kwargs...)

Plot the percentage of CO2 in each layer as a function of time on a single plot.
Uses a stacked area chart to show the distribution across layers.

Parameters:
- `snapshots`: Vector of ReservoirSnapshot from generate_reservoir_snapshots
- `output_file`: Path to save figure (default: nothing, returns figure)
- `title`: Plot title (default: "CO2 Distribution Across Layers")
- `show_leaked`: If true, shows leaked fraction as top area (default: true)
- `colormap`: Colormap to use for layers (default: :viridis)
"""
function plot_layer_fractions_timeseries(
    snapshots::Vector{ReservoirSnapshot};
    output_file::Union{String, Nothing} = nothing,
    title::String = "CO2 Distribution Across Layers",
    show_leaked::Bool = true,
    colormap::Symbol = :viridis
)
    n_snapshots = length(snapshots)
    n_layers = length(snapshots[1].layer_snapshots)

    # Extract time series data
    times = [s.timestamp for s in snapshots]

    # Compute fractions relative to total injected (to show where all CO2 went)
    fractions_by_layer = Matrix{Float64}(undef, n_snapshots, n_layers)
    leaked_fraction = zeros(n_snapshots)

    for (t_idx, s) in enumerate(snapshots)
        total_injected = s.total_injected_m3
        if total_injected > 0
            for layer_idx in 1:n_layers
                fractions_by_layer[t_idx, layer_idx] = s.stored_by_layer_m3[layer_idx] / total_injected * 100
            end
            leaked_fraction[t_idx] = s.total_leaked_m3 / total_injected * 100
        else
            fractions_by_layer[t_idx, :] .= 0.0
            leaked_fraction[t_idx] = 0.0
        end
    end

    # Get layer names
    layer_names = [snapshots[1].layer_snapshots[i].layer_name for i in 1:n_layers]

    # Create figure
    fig = Figure(size = (800, 500))

    ax = Axis(fig[1, 1],
              title = title,
              xlabel = "Time (years)",
              ylabel = "Fraction of Injected CO2 (%)",
              limits = (nothing, nothing, 0, 105))

    # Generate colors from colormap
    cmap = cgrad(colormap, n_layers + (show_leaked ? 1 : 0), categorical = true)

    # Build cumulative sums for stacked area chart
    cumulative = zeros(n_snapshots)
    bands = []
    labels = String[]

    for layer_idx in 1:n_layers
        lower = copy(cumulative)
        upper = cumulative .+ fractions_by_layer[:, layer_idx]

        b = band!(ax, times, lower, upper, color = cmap[layer_idx], label = layer_names[layer_idx])
        push!(bands, b)
        push!(labels, layer_names[layer_idx])

        cumulative = upper
    end

    # Add leaked fraction on top
    if show_leaked
        lower = copy(cumulative)
        upper = cumulative .+ leaked_fraction

        b = band!(ax, times, lower, upper, color = (:red, 0.6), label = "Leaked")
        push!(bands, b)
        push!(labels, "Leaked")
    end

    # Add legend
    Legend(fig[1, 2], bands, labels, "Location", framevisible = true, labelsize = 11)

    # Add horizontal line at 100%
    hlines!(ax, [100], color = :black, linestyle = :dash, linewidth = 1)

    if !isnothing(output_file)
        save(output_file, fig)
        println("Figure saved to: $(output_file)")
    end

    return fig
end
