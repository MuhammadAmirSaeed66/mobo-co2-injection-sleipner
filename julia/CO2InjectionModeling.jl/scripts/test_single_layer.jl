using CO2InjectionModeling
using SurfaceWaterIntegratedModeling

# Setup similar to debugging_single_layer.jl
boundary_condition = :closed
topography = load_sleipner_topography()
domain = create_domain_from_topography(topography, 1.0)
layers = analyze_base_surfaces(topography; boundary_condition=boundary_condition)

include("debugging_utils.jl")
injection_events = generate_injection_events(layers)

# Create reservoir properties with a VERY LOW leakage height to force leakage
layer_idx = 1
original_props = generate_reservoir_properties_for_sleipner_layers()[layer_idx]

# Use a leakage height of 16 meters
leakage_height = 30.0
rprops = ReservoirProperties(
    original_props.sand_porosity,
    original_props.sand_residual_co2_saturation,
    original_props.sand_irreducible_water_saturation,
    original_props.shale_pressure_threshold,
    leakage_height,  
    original_props.residual_leakage_time
)

tstruct = layers[layer_idx].trap_structure
weather_events_layer = convert_injection_event_to_weather_event(injection_events[layer_idx], rprops, domain)

println("\n=== Forced Leakage Test ===")
println("Layer $layer_idx: $(layers[layer_idx].name)")
println("Number of traps: $(numtraps(tstruct))")

# Run fill_layer with leakage detection
println("\nRunning fill_layer with leakage detection...")
seq, leakage_state = fill_layer(
    tstruct,
    domain,
    rprops,
    weather_events_layer;
    verbose=false
)

println("\n=== Leakage State ===")
println("Number of leaking traps: $(sum(leakage_state.leaking))")
println("Number of draining traps: $(sum(leakage_state.draining))")
println("Number of leakage records: $(length(leakage_state.leakage_records))")
println("Residual saturation: $(rprops.sand_residual_co2_saturation) ($(round(100*(1-rprops.sand_residual_co2_saturation), digits=1))% will drain)")
println("Residual leakage time: $(rprops.residual_leakage_time) years")

# Compute drainage statistics
total_initial_vol_draining = sum(leakage_state.initial_volume_at_leak[i] for i in 1:numtraps(tstruct) if leakage_state.draining[i]; init=0.0)
total_residual_vol = total_initial_vol_draining * rprops.sand_residual_co2_saturation
total_drained_vol = total_initial_vol_draining * (1 - rprops.sand_residual_co2_saturation)
println("\nDrainage statistics (SWIM units):")
println("  Total initial volume in draining traps: $(round(total_initial_vol_draining, digits=1))")
println("  Total residual volume (after drainage): $(round(total_residual_vol, digits=1))")
println("  Total CO2 that will drain out: $(round(total_drained_vol, digits=1))")

if !isempty(leakage_state.leakage_records)
    println("\nLeakage events:")
    for (i, record) in enumerate(leakage_state.leakage_records)
        ancestors = get_all_parents(tstruct, record.trap_id)
        descendants = get_all_descendants(tstruct, record.trap_id)
        leaking_ancestors = filter(a -> leakage_state.leaking[a], ancestors)
        draining_descendants = filter(d -> leakage_state.draining[d], descendants)
        trap_id = record.trap_id
        initial_vol = leakage_state.initial_volume_at_leak[trap_id]
        residual_vol = initial_vol * rprops.sand_residual_co2_saturation
        drainage_end_time = record.start_time + rprops.residual_leakage_time

        # Volume in descendants
        desc_initial_vol = sum(leakage_state.initial_volume_at_leak[d] for d in draining_descendants; init=0.0)
        desc_residual_vol = desc_initial_vol * rprops.sand_residual_co2_saturation

        println("  $i. Trap $(trap_id) started leaking at time $(round(record.start_time, digits=4)) years")
        println("     Location: $(record.leakage_location)")
        println("     Initial volume at leak (SWIM units): $(round(initial_vol, digits=2))")
        println("     Residual volume after drainage: $(round(residual_vol, digits=2))")
        println("     Drainage ends at: $(round(drainage_end_time, digits=4)) years")
        println("     Leaking ancestors: $(length(leaking_ancestors)) traps")
        println("     Draining descendants: $(length(draining_descendants)) traps")
        println("       Descendants initial volume: $(round(desc_initial_vol, digits=2))")
        println("       Descendants residual volume: $(round(desc_residual_vol, digits=2))")
    end
end

# If there was leakage, generate weather events for the next layer
if !isempty(leakage_state.leakage_records)
    println("\n=== Weather Events for Next Layer ===")
    target_grid_size = size(tstruct.topography)
    upstream_weather = generate_leakage_weather_events(
        leakage_state,
        weather_events_layer,
        seq,
        tstruct,
        target_grid_size
    )
    println("Generated $(length(upstream_weather)) weather events for the overlying layer")

    println("\nFirst 5 weather events:")
    for (i, we) in enumerate(upstream_weather[1:min(5, length(upstream_weather))])
        total_rate = sum(we.rain_rate)
        println("  Event $i: timestamp=$(round(we.timestamp, digits=4)), total_leakage_rate=$(round(total_rate, digits=4))")
    end

    # Check events around year boundaries to verify rate changes
    println("\nWeather events around year boundaries:")
    for target_year in [1.0, 2.0, 3.0]
        # Find events just before and after this year
        before_events = filter(we -> we.timestamp < target_year && we.timestamp > target_year - 0.1, upstream_weather)
        after_events = filter(we -> we.timestamp >= target_year && we.timestamp < target_year + 0.1, upstream_weather)

        if !isempty(before_events)
            we = before_events[end]
            println("  Before year $(target_year): t=$(round(we.timestamp, digits=4)), rate=$(round(sum(we.rain_rate), digits=2))")
        end
        if !isempty(after_events)
            we = after_events[1]
            println("  After year $(target_year):  t=$(round(we.timestamp, digits=4)), rate=$(round(sum(we.rain_rate), digits=2))")
        end
    end

    # === Verify mass conservation at multiple timepoints ===
    println("\n=== Mass Conservation Verification Over Time ===")
    println("Checking: Injected = Stored + Leaked (from weather events)")

    # Count spill events where traps become filled
    leakage_event_time = leakage_state.leakage_records[1].start_time
    traps_filled_at_leakage = 0
    for (i, se) in enumerate(seq)
        if !isempty(se.filled) && se.filled isa Vector{IncrementalUpdate{Bool}}
            newly_filled = [u.index for u in se.filled if u.value == true]
            if se.timestamp ≈ leakage_event_time
                traps_filled_at_leakage = length(newly_filled)
            end
        end
    end
    println("\nAt leakage event (t=$(round(leakage_event_time, digits=4))): $(traps_filled_at_leakage) traps marked as filled")
    println("  (Original trap + $(traps_filled_at_leakage - 1) ancestors)")
    println()

    check_times = [1.0, 2.0, 5.0, 6.0, 6.5, 7.0, 10.0, 15.0]

    for check_time in check_times
        # 1. Compute total injected up to check_time
        injected = CO2InjectionModeling.compute_total_injected_amount(injection_events[layer_idx], check_time)

        # 2. Compute total stored in layer at check_time (SWIM units -> physical)
        # Pass leakage_state to account for residual drainage
        stored_swim = CO2InjectionModeling.compute_total_stored_volume(seq, tstruct, check_time; leakage_state=leakage_state)
        stored = CO2InjectionModeling.swim_volume_to_physical_volume(stored_swim, rprops, domain)

        # 3. Compute total leaked by integrating weather event rates
        # The weather events give rates at discrete times; integrate piecewise
        leaked = 0.0
        for i in 1:length(upstream_weather)
            we = upstream_weather[i]
            t_start = we.timestamp

            # Skip events after check_time
            if t_start >= check_time
                break
            end

            # Find end of this rate period
            t_end = (i < length(upstream_weather)) ? upstream_weather[i+1].timestamp : check_time
            t_end = min(t_end, check_time)

            # Integrate rate over [t_start, t_end]
            dt = t_end - t_start
            if dt > 0
                rate = sum(we.rain_rate)
                # Convert rate from SWIM units to physical volume
                leaked += CO2InjectionModeling.swim_volume_to_physical_volume(rate * dt, rprops, domain)
            end
        end

        # 4. Check mass balance
        balance = injected - stored - leaked
        rel_error = injected > 0 ? abs(balance) / injected * 100 : 0.0

        status = rel_error < 0.1 ? "✓" : (rel_error < 2.0 ? "~" : "✗")
        println("  Time $(check_time) years: $(status)")
        println("    Injected: $(round(injected / 1e6, digits=3)) M m³")
        println("    Stored:   $(round(stored / 1e6, digits=3)) M m³")
        println("    Leaked:   $(round(leaked / 1e6, digits=3)) M m³")
        println("    Balance:  $(round(balance / 1e6, digits=6)) M m³ ($(round(rel_error, digits=4))% error)")
    end

end

println("\nTest completed!")

# Generate visualization
println("\n=== Generating Animations ===")
println("Generating height-based animation...")
animate_single_layer_filling(
    layers[layer_idx],
    seq,
    domain;
    output_file="leakage_forced_filling.gif",
    num_frames=30,
    end_time=20.0,
    fps=3,
    # colormap=:phase,
    max_CO2_height=leakage_height
)

println("\nGenerating saturation-based animation...")
animate_single_layer_saturation(
    layers[layer_idx],
    seq,
    domain,
    leakage_state;
    output_file="leakage_forced_saturation.gif",
    num_frames=30,
    end_time=20.0,
    fps=3,
    colormap=:viridis
)
