using Random
using CairoMakie


"""
    animate_single_layer_filling(layer, seq, domain; kwargs...)

Create an animated bird's eye view of CO2 filling in a single layer.

Parameters:
- `layer`: Layer struct containing trap_structure
- `seq`: Vector{SpillEvent} from fill_layer
- `domain`: Domain3D struct
- `output_file`: Path to save animation (default: "layer_filling.gif")
- `num_frames`: Number of frames in animation (default: 30)
- `start_time`: Start time for animation (default: 0.0)
- `end_time`: End time for animation (default: auto-detect from seq)
- `fps`: Frames per second (default: 2)
- `colormap`: Colormap for heights (default: :thermal)
"""
function animate_single_layer_filling(
    layer::Layer,
    seq::Vector{SpillEvent},
    domain::Domain3D;
    output_file::String = "layer_filling.gif",
    num_frames::Int = 30,
    start_time::Float64 = 0.0,
    end_time::Union{Float64, Nothing} = nothing,
    fps::Int = 2,
    colormap::Symbol = :thermal,
    max_CO2_height::Float64 = 20.0
)
    tstruct = layer.trap_structure
    num_traps = numtraps(tstruct)

    # Determine end_time from sequence if not provided
    if isnothing(end_time)
        end_time = maximum(se.timestamp for se in seq)
        if !isfinite(end_time)
            end_time = seq[end-1].timestamp + 1.0
        end
    end

    # Generate timepoints for animation
    timepoints = collect(range(start_time, stop=end_time, length=num_frames))

    # Get trap states at each timepoint
    println("Computing trap states for $(num_frames) frames...")
    tstates = trap_states_at_timepoints(tstruct, seq, timepoints; verbose=false)

    # Compute z_vol_tables for height conversion
    z_vol_tables = SurfaceWaterIntegratedModeling._compute_z_vol_tables(tstruct)

    # Get grid size (remove padding if present)
    pad = layer.boundary_padding
    topo_size = size(tstruct.topography)
    nx_padded, ny_padded = topo_size
    nx = nx_padded - 2 * pad
    ny = ny_padded - 2 * pad

    # Set up figure
    fig = Figure(size = (800, 700))

    x_coords = range(0, nx * domain.dx, length=nx)
    y_coords = range(0, ny * domain.dy, length=ny)

    # Create observables for animation
    height_data = Observable(zeros(Float64, nx, ny))
    time_text = Observable("Time: 0.0 years")

    ax = Axis(fig[1, 1],
              xlabel = "X (m)",
              ylabel = "Y (m)",
              title = time_text,
              aspect = DataAspect())

    hm = heatmap!(ax, x_coords, y_coords, height_data,
                  colormap = colormap,
                  colorrange = (0.0, max_CO2_height))  # Adjust based on expected max height

    Colorbar(fig[1, 2], hm, label = "CO2 Height (m)")

    # Add layer info
    Label(fig[0, :], "Layer: $(layer.name)", fontsize = 16)

    println("Creating animation...")

    # Create the animation
    record(fig, output_file, eachindex(timepoints); framerate=fps) do frame_idx
        tp = timepoints[frame_idx]
        filled, volumes, _ = tstates[frame_idx]

        # Create height map (on padded grid, then remove padding)
        height_map_padded = zeros(Float64, nx_padded, ny_padded)

        max_height = 0.0
        for trap_id in 1:num_traps
            volume = volumes[trap_id]
            # Show height for traps with volume > 0, OR for filled parent traps
            # (parent traps with volume = 0 still have CO2 height from children)
            if volume > 0.0 || filled[trap_id]
                # Convert volume to height
                z_vol_table = z_vol_tables[trap_id]
                height = volume_to_height(volume, trap_id, z_vol_table, tstruct)
                max_height = max(max_height, height)

                # Fill the trap footprint with this height
                footprint = tstruct.footprints[trap_id]
                for idx in footprint
                    height_map_padded[idx] = max(height_map_padded[idx], height)
                end
            end
        end

        # Remove padding
        height_map = height_map_padded[pad+1:end-pad, pad+1:end-pad]

        # Update observables
        height_data[] = height_map

        total_volume = sum(volumes)
        time_text[] = "Time: $(round(tp, digits=2)) years | Max height: $(round(max_height, digits=1)) m | Total vol: $(round(total_volume, digits=0))"

        if frame_idx % 10 == 0 || frame_idx == length(timepoints)
            println("  Frame $(frame_idx)/$(length(timepoints))")
        end
    end

    println("Animation saved to: $(output_file)")
    return nothing
end


"""
    plot_single_layer_snapshot(layer, seq, domain, time; kwargs...)

Create a static plot of CO2 distribution at a specific time.
"""
function plot_single_layer_snapshot(
    layer::Layer,
    seq::Vector{SpillEvent},
    domain::Domain3D,
    time::Float64;
    output_file::Union{String, Nothing} = nothing,
    colormap::Symbol = :thermal
)
    tstruct = layer.trap_structure
    num_traps = numtraps(tstruct)

    # Get trap state at this time
    tstates = trap_states_at_timepoints(tstruct, seq, [time]; verbose=false)
    filled, volumes, _ = tstates[1]

    # Compute z_vol_tables
    z_vol_tables = SurfaceWaterIntegratedModeling._compute_z_vol_tables(tstruct)

    # Get grid size
    pad = layer.boundary_padding
    topo_size = size(tstruct.topography)
    nx_padded, ny_padded = topo_size
    nx = nx_padded - 2 * pad
    ny = ny_padded - 2 * pad

    # Create height map
    height_map_padded = zeros(Float64, nx_padded, ny_padded)

    max_height = 0.0
    for trap_id in 1:num_traps
        volume = volumes[trap_id]
        if volume > 0.0
            z_vol_table = z_vol_tables[trap_id]
            height = volume_to_height(volume, trap_id, z_vol_table, tstruct)
            max_height = max(max_height, height)

            footprint = tstruct.footprints[trap_id]
            for idx in footprint
                height_map_padded[idx] = max(height_map_padded[idx], height)
            end
        end
    end

    # Remove padding
    height_map = height_map_padded[pad+1:end-pad, pad+1:end-pad]

    # Create figure
    fig = Figure(size = (800, 700))

    x_coords = range(0, nx * domain.dx, length=nx)
    y_coords = range(0, ny * domain.dy, length=ny)

    ax = Axis(fig[1, 1],
              xlabel = "X (m)",
              ylabel = "Y (m)",
              title = "Time: $(round(time, digits=2)) years | Max height: $(round(max_height, digits=1)) m",
              aspect = DataAspect())

    hm = heatmap!(ax, x_coords, y_coords, height_map,
                  colormap = colormap,
                  colorrange = (0.0, max(max_height, 1.0)))

    Colorbar(fig[1, 2], hm, label = "CO2 Height (m)")
    Label(fig[0, :], "Layer: $(layer.name)", fontsize = 16)

    if !isnothing(output_file)
        save(output_file, fig)
        println("Figure saved to: $(output_file)")
    end

    return fig
end


function generate_injection_events(layers::Vector{Layer})
    # Define injection events
    # Make xy in the center of the domain
    trap_topo = layers[1].trap_structure.topography
    xy = CartesianIndex(div(size(trap_topo, 1), 2), div(size(trap_topo, 2), 2))

    # Some constant rate for testing
    # annual_rates_mt = fill(0.86, 15)  # Mt/year
    annual_rates_mt = fill(0.86, 15) ./ 1.0  # Mt/year

    # Add some random variation to the rates
    Random.seed!(42)  # For reproducibility
    for i in 1:length(annual_rates_mt)
        variation = randn() * 0.05  # Small random variation
        annual_rates_mt[i] += variation
        annual_rates_mt[i] = max(annual_rates_mt[i], 0.0)  # Ensure non-negative rates
    end

    println("Annual injection rates (Mt/year): ", annual_rates_mt)

    # Convert Mt/year to m³/year
    co2_density_l1 = 425.0  # The same as used in reservoir properties
    annual_rates_m3_per_year = annual_rates_mt .* 1e9 ./ co2_density_l1

    # Create injection events for bottom layer (L1)
    # Time points are cumulative: start at year 0, events mark end of each year
    n_events = length(annual_rates_mt)
    bottom_layer_events = Vector{InjectionEvent}(undef, n_events)

    # Get the grid size from the bottom layer
    grid_size = size(layers[1].trap_structure.topography)

    for (i, rate) in enumerate(annual_rates_m3_per_year)
        # Time in years (0, 1, 2, ..., 14)
        time = float(i - 1)

        # Create injection rate field (only inject at specified cell)
        injection_rate = zeros(grid_size)
        injection_rate[xy] = rate

        # if i > 5
        #     injection_rate[xy[1] + 5, xy[2] + 5] = rate
        # end

        # Optionally inject at a second location (trap 9 at (64, 4)) to test multi-site injection
        # Note: (end - 3, end - 3) was outside all leaf traps - CO2 went to runoff
        # injection_rate[64, 4] = rate
        # if i > 3
        #     injection_rate[end - 7, end - 3] = rate
        # end

        bottom_layer_events[i] = InjectionEvent(time, injection_rate)
    end

    # Create zero injection events for all other layers
    zero_injection = zeros(grid_size)
    zero_event = [InjectionEvent(0.0, zero_injection)]

    # Assemble injection events for all layers
    n_layers = length(layers)
    injection_events = Vector{Vector{InjectionEvent}}(undef, n_layers)
    injection_events[1] = bottom_layer_events  # Bottom layer (L1) has actual injection
    for i in 2:n_layers
        injection_events[i] = zero_event  # All other layers have zero injection
    end

    return injection_events
end

"""
    animate_single_layer_saturation(layer, seq, domain, leakage_state; kwargs...)

Create an animated bird's eye view of CO2 saturation in a single layer.

This visualization shows CO2 saturation (volume/capacity) instead of height,
which makes residual drainage visible: after leakage starts, the saturation
in leaking traps decreases from 1.0 down to the residual saturation level.

Parameters:
- `layer`: Layer struct containing trap_structure
- `seq`: Vector{SpillEvent} from fill_layer
- `domain`: Domain3D struct
- `leakage_state`: LeakageState from fill_layer (required for drainage calculation)
- `output_file`: Path to save animation (default: "layer_saturation.gif")
- `num_frames`: Number of frames in animation (default: 30)
- `start_time`: Start time for animation (default: 0.0)
- `end_time`: End time for animation (default: auto-detect from seq)
- `fps`: Frames per second (default: 2)
- `colormap`: Colormap for saturation (default: :viridis)
"""
function animate_single_layer_saturation(
    layer::Layer,
    seq::Vector{SpillEvent},
    domain::Domain3D,
    leakage_state::LeakageState;
    output_file::String = "layer_saturation.gif",
    num_frames::Int = 30,
    start_time::Float64 = 0.0,
    end_time::Union{Float64, Nothing} = nothing,
    fps::Int = 2,
    colormap::Symbol = :viridis
)
    tstruct = layer.trap_structure
    num_traps = numtraps(tstruct)

    # Determine end_time from sequence if not provided
    if isnothing(end_time)
        end_time = maximum(se.timestamp for se in seq)
        if !isfinite(end_time)
            end_time = seq[end-1].timestamp + 1.0
        end
    end

    # Generate timepoints for animation
    timepoints = collect(range(start_time, stop=end_time, length=num_frames))

    # Get trap states at each timepoint
    println("Computing trap states for $(num_frames) frames...")
    tstates = trap_states_at_timepoints(tstruct, seq, timepoints; verbose=false)

    # Get trap capacities for saturation calculation
    # Capacity = trapvolume - subvolume (same as fill volume)
    trap_capacities = [tstruct.trapvolumes[i] - tstruct.subvolumes[i] for i in 1:num_traps]

    # Get grid size (remove padding if present)
    pad = layer.boundary_padding
    topo_size = size(tstruct.topography)
    nx_padded, ny_padded = topo_size
    nx = nx_padded - 2 * pad
    ny = ny_padded - 2 * pad

    # Set up figure
    fig = Figure(size = (800, 700))

    x_coords = range(0, nx * domain.dx, length=nx)
    y_coords = range(0, ny * domain.dy, length=ny)

    # Create observables for animation
    saturation_data = Observable(zeros(Float64, nx, ny))
    time_text = Observable("Time: 0.0 years")

    ax = Axis(fig[1, 1],
              xlabel = "X (m)",
              ylabel = "Y (m)",
              title = time_text,
              aspect = DataAspect())

    hm = heatmap!(ax, x_coords, y_coords, saturation_data,
                  colormap = colormap,
                  colorrange = (0.0, 1.0))

    Colorbar(fig[1, 2], hm, label = "CO2 Saturation (fraction of capacity)")

    # Add layer info
    Label(fig[0, :], "Layer: $(layer.name) - CO2 Saturation", fontsize = 16)

    println("Creating animation...")

    # Create the animation
    record(fig, output_file, eachindex(timepoints); framerate=fps) do frame_idx
        tp = timepoints[frame_idx]
        filled, volumes, _ = tstates[frame_idx]

        # Create saturation map (on padded grid, then remove padding)
        saturation_map_padded = zeros(Float64, nx_padded, ny_padded)

        max_saturation = 0.0
        total_volume = 0.0
        num_leaking = 0

        for trap_id in 1:num_traps
            volume = volumes[trap_id]

            # For draining traps (leaking traps and their filled ancestors), compute drainage-adjusted volume
            # Note: we check 'draining' not 'leaking' because ancestors can drain without being at the leakage threshold
            if leakage_state.draining[trap_id]
                drained_vol = compute_volume_with_drainage(trap_id, tp, leakage_state)
                if !isnothing(drained_vol)
                    volume = drained_vol
                end
            end
            if leakage_state.leaking[trap_id]
                num_leaking += 1
            end

            total_volume += volume

            # Compute saturation = volume / capacity
            capacity = trap_capacities[trap_id]
            if capacity > 0.0
                saturation = min(1.0, volume / capacity)
            else
                saturation = 0.0
            end

            max_saturation = max(max_saturation, saturation)

            # Fill the trap footprint with this saturation
            if saturation > 0.0
                footprint = tstruct.footprints[trap_id]
                for idx in footprint
                    saturation_map_padded[idx] = max(saturation_map_padded[idx], saturation)
                end
            end
        end

        # Remove padding
        saturation_map = saturation_map_padded[pad+1:end-pad, pad+1:end-pad]

        # Update observables
        saturation_data[] = saturation_map

        time_text[] = "Time: $(round(tp, digits=2)) years | Max sat: $(round(max_saturation, digits=2)) | Leaking: $(num_leaking) traps"

        if frame_idx % 10 == 0 || frame_idx == length(timepoints)
            println("  Frame $(frame_idx)/$(length(timepoints))")
        end
    end

    println("Animation saved to: $(output_file)")
    return nothing
end