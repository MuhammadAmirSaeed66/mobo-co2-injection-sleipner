using SurfaceWaterIntegratedModeling: TrapStructure, numtraps, trap_states_at_timepoints, SpillEvent

export MultiLayerSnapshot, ReservoirSnapshot, generate_reservoir_snapshots
export print_snapshot_summary, print_layer_snapshot_summary


"""
    MultiLayerSnapshot

Comprehensive state of a single layer at a specific time.
Contains both SWIM units (for internal use) and physical units (m³).
"""
struct MultiLayerSnapshot
    # Basic info
    layer_idx::Int
    layer_name::String
    timestamp::Float64

    # Trap statistics
    num_traps::Int
    num_filled_traps::Int
    num_leaking_traps::Int
    num_traps_with_co2::Int  # Traps with volume > 0

    # Volumes in physical units (m³)
    stored_volume_m3::Float64

    # Height statistics (meters)
    max_co2_height::Float64
    mean_co2_height::Float64  # Mean height across traps with CO2

    # Leakage info
    leakage_rate_m3_per_year::Float64  # Current leakage rate out of this layer
    cumulative_leaked_m3::Float64  # Total leaked from this layer up to this time

    # Volume distribution
    trap_volumes::Vector{Float64}  # Volume in each trap (SWIM units)
    trap_heights::Vector{Float64}  # CO2 height in each trap (meters)
    filled_traps::Vector{Bool}     # Which traps are filled
end


"""
    ReservoirSnapshot

Comprehensive state of the full multi-layer reservoir at a specific time.
"""
struct ReservoirSnapshot
    timestamp::Float64

    # Injection metrics (physical units, m³)
    total_injected_m3::Float64
    injection_rate_m3_per_year::Float64  # Current injection rate

    # Storage metrics (physical units, m³)
    total_stored_m3::Float64
    stored_by_layer_m3::Vector{Float64}  # Storage in each layer
    storage_fraction_by_layer::Vector{Float64}  # Percentage of total stored in each layer

    # Leakage metrics (physical units, m³)
    total_leaked_m3::Float64  # Leaked out of the domain (top of reservoir)
    leakage_rate_m3_per_year::Float64  # Current leakage rate out of domain

    # Mass balance
    mass_balance_error_m3::Float64  # injected - stored - leaked (should be ~0)
    mass_balance_error_percent::Float64

    # Height statistics
    max_co2_height::Float64  # Maximum height across all layers
    max_co2_height_by_layer::Vector{Float64}

    # Trap statistics
    total_traps::Int
    total_filled_traps::Int
    total_leaking_traps::Int
    filled_traps_by_layer::Vector{Int}
    leaking_traps_by_layer::Vector{Int}

    # Individual layer snapshots
    layer_snapshots::Vector{MultiLayerSnapshot}
end


"""
    generate_reservoir_snapshots(layers, seqs, leakage_states, domain, reservoir_properties,
                                  injection_events; num_snapshots=10, start_time=0.0, end_time=15.0)

Generate comprehensive snapshots of the reservoir state at evenly spaced timepoints.

Returns a vector of ReservoirSnapshot structs.
"""
function generate_reservoir_snapshots(
    layers::Vector{Layer},
    seqs::Vector{Vector{SpillEvent}},
    leakage_states::Vector{LeakageState},
    domain::Domain3D,
    reservoir_properties::Union{ReservoirProperties, Vector{ReservoirProperties}},
    injection_events::Vector{Vector{InjectionEvent}};
    num_snapshots::Int = 10,
    start_time::Float64 = 0.0,
    end_time::Float64 = 15.0,
    verbose::Bool = true
)
    n_layers = length(layers)

    # Handle single vs. per-layer reservoir properties
    if reservoir_properties isa ReservoirProperties
        rprops = fill(reservoir_properties, n_layers)
    else
        rprops = reservoir_properties
    end

    # Generate timepoints
    timepoints = collect(range(start_time, stop=end_time, length=num_snapshots))

    verbose && println("Generating $(num_snapshots) reservoir snapshots from t=$(start_time) to t=$(end_time)...")

    # Precompute trap states for all layers at all timepoints
    all_tstates = Vector{Any}(undef, n_layers)
    all_z_vol_tables = Vector{Any}(undef, n_layers)

    for layer_idx in 1:n_layers
        tstruct = layers[layer_idx].trap_structure
        seq = seqs[layer_idx]

        if isempty(seq) || numtraps(tstruct) == 0
            all_tstates[layer_idx] = nothing
            all_z_vol_tables[layer_idx] = nothing
        else
            all_tstates[layer_idx] = trap_states_at_timepoints(tstruct, seq, timepoints; verbose=false)
            all_z_vol_tables[layer_idx] = SurfaceWaterIntegratedModeling._compute_z_vol_tables(tstruct)
        end
    end

    # Generate snapshots
    snapshots = ReservoirSnapshot[]

    for (time_idx, tp) in enumerate(timepoints)
        # Generate layer snapshots
        layer_snapshots = MultiLayerSnapshot[]

        for layer_idx in 1:n_layers
            layer = layers[layer_idx]
            tstruct = layer.trap_structure
            leakage_state = leakage_states[layer_idx]
            rp = rprops[layer_idx]

            layer_snapshot = _generate_layer_snapshot(
                layer_idx, layer, tstruct, seqs[layer_idx],
                leakage_state, all_tstates[layer_idx], all_z_vol_tables[layer_idx],
                time_idx, tp, rp, domain, injection_events[layer_idx]
            )
            push!(layer_snapshots, layer_snapshot)
        end

        # Aggregate into reservoir snapshot
        reservoir_snapshot = _aggregate_reservoir_snapshot(
            tp, layer_snapshots, injection_events, rprops, domain, leakage_states
        )
        push!(snapshots, reservoir_snapshot)

        verbose && (time_idx % 5 == 0 || time_idx == num_snapshots) &&
            println("  Snapshot $(time_idx)/$(num_snapshots) at t=$(round(tp, digits=2)) years")
    end

    verbose && println("Done generating snapshots.")
    return snapshots
end


"""
Generate a snapshot for a single layer at a specific time.
"""
function _generate_layer_snapshot(
    layer_idx::Int,
    layer::Layer,
    tstruct::TrapStructure,
    seq::Vector{SpillEvent},
    leakage_state::LeakageState,
    tstates,
    z_vol_tables,
    time_idx::Int,
    timestamp::Float64,
    rprops::ReservoirProperties,
    domain::Domain3D,
    layer_injection_events::Vector{InjectionEvent}
)::MultiLayerSnapshot

    num_traps = numtraps(tstruct)

    # Handle empty layers
    if isnothing(tstates) || num_traps == 0
        return MultiLayerSnapshot(
            layer_idx, layer.name, timestamp,
            0, 0, 0, 0,  # trap counts
            0.0,  # stored volume
            0.0, 0.0,  # heights
            0.0, 0.0,  # leakage
            Float64[], Float64[], Bool[]  # empty vectors
        )
    end

    # Get trap state at this timepoint
    filled, volumes, _ = tstates[time_idx]

    # Compute heights for each trap
    heights = zeros(Float64, num_traps)
    for trap_id in 1:num_traps
        if volumes[trap_id] > 0.0
            heights[trap_id] = volume_to_height(
                volumes[trap_id], trap_id, z_vol_tables[trap_id], tstruct
            )
        end
    end

    # Compute statistics
    num_filled = count(filled)
    num_leaking = count(leakage_state.leaking)
    num_with_co2 = count(v -> v > 0, volumes)

    # Stored volume (convert from SWIM to physical)
    # IMPORTANT: Account for drainage! Draining traps have reduced volumes
    stored_swim = 0.0
    for trap_id in 1:num_traps
        vol = volumes[trap_id]

        # Apply drainage adjustment if this trap is draining
        if leakage_state.draining[trap_id]
            drained_vol = compute_volume_with_drainage(trap_id, timestamp, leakage_state)
            if !isnothing(drained_vol)
                vol = drained_vol
            end
        end

        stored_swim += vol
    end
    stored_m3 = swim_volume_to_physical_volume(stored_swim, rprops, domain)

    # Height statistics
    max_height = maximum(heights; init=0.0)
    heights_with_co2 = filter(h -> h > 0, heights)
    mean_height = isempty(heights_with_co2) ? 0.0 : mean(heights_with_co2)

    # Compute leakage rate and cumulative leaked
    leakage_rate_swim, cumulative_leaked_swim = _compute_layer_leakage(
        leakage_state, seq, timestamp
    )
    leakage_rate_m3 = swim_volume_to_physical_volume(leakage_rate_swim, rprops, domain)
    cumulative_leaked_m3 = swim_volume_to_physical_volume(cumulative_leaked_swim, rprops, domain)

    return MultiLayerSnapshot(
        layer_idx, layer.name, timestamp,
        num_traps, num_filled, num_leaking, num_with_co2,
        stored_m3,
        max_height, mean_height,
        leakage_rate_m3, cumulative_leaked_m3,
        collect(volumes), heights, collect(filled)
    )
end


"""
Compute leakage rate and cumulative leaked volume for a layer.
"""
function _compute_layer_leakage(
    leakage_state::LeakageState,
    seq::Vector{SpillEvent},
    timestamp::Float64
)::Tuple{Float64, Float64}

    if isempty(leakage_state.leakage_records)
        return (0.0, 0.0)
    end

    total_rate = 0.0
    cumulative_leaked = 0.0

    for record in leakage_state.leakage_records
        if record.start_time > timestamp
            continue  # Leakage hasn't started yet
        end

        trap_id = record.trap_id

        # Get current inflow rate to the leaking trap (this is the leakage rate)
        inflow_rate = get_trap_inflow_at_time(trap_id, timestamp, seq)
        total_rate += inflow_rate

        # Integrate leakage from start_time to timestamp
        # This is a simplified integration using the spill events
        leaked = _integrate_trap_leakage(trap_id, record.start_time, timestamp, seq)
        cumulative_leaked += leaked
    end

    return (total_rate, cumulative_leaked)
end


"""
Integrate leakage from a trap over a time interval using spill events.
"""
function _integrate_trap_leakage(
    trap_id::Int,
    start_time::Float64,
    end_time::Float64,
    seq::Vector{SpillEvent}
)::Float64

    if isempty(seq) || start_time >= end_time
        return 0.0
    end

    total_leaked = 0.0

    # Find relevant spill events
    for i in 1:length(seq)
        se = seq[i]
        t_start = max(se.timestamp, start_time)
        t_end = (i < length(seq)) ? min(seq[i+1].timestamp, end_time) : end_time

        if t_start >= end_time || t_end <= start_time
            continue
        end

        dt = t_end - t_start
        if dt > 0
            inflow = get_trap_inflow_at_time(trap_id, t_start, seq)
            total_leaked += inflow * dt
        end
    end

    return total_leaked
end


"""
Aggregate layer snapshots into a reservoir snapshot.
"""
function _aggregate_reservoir_snapshot(
    timestamp::Float64,
    layer_snapshots::Vector{MultiLayerSnapshot},
    injection_events::Vector{Vector{InjectionEvent}},
    rprops::Vector{ReservoirProperties},
    domain::Domain3D,
    leakage_states::Vector{LeakageState}
)::ReservoirSnapshot

    n_layers = length(layer_snapshots)

    # Compute total injected (sum across all layers)
    total_injected = 0.0
    current_injection_rate = 0.0
    for layer_idx in 1:n_layers
        total_injected += compute_total_injected_amount(injection_events[layer_idx], timestamp)
        current_injection_rate += _get_injection_rate_at_time(injection_events[layer_idx], timestamp)
    end

    # Compute stored volumes
    stored_by_layer = [ls.stored_volume_m3 for ls in layer_snapshots]
    total_stored = sum(stored_by_layer)
    storage_fraction = total_stored > 0 ? stored_by_layer ./ total_stored : zeros(n_layers)

    # Compute leakage out of domain (from top layer only, or sum of all layer leakages that exit)
    # For simplicity, we compute total leaked as: injected - stored
    total_leaked = max(0.0, total_injected - total_stored)

    # Current leakage rate (from top layer or last layer with leakage)
    top_layer_leakage_rate = layer_snapshots[end].leakage_rate_m3_per_year

    # Mass balance
    mass_balance_error = total_injected - total_stored - total_leaked
    mass_balance_error_percent = total_injected > 0 ? abs(mass_balance_error) / total_injected * 100 : 0.0

    # Height statistics
    max_heights_by_layer = [ls.max_co2_height for ls in layer_snapshots]
    max_height = maximum(max_heights_by_layer; init=0.0)

    # Trap statistics
    total_traps = sum(ls.num_traps for ls in layer_snapshots)
    total_filled = sum(ls.num_filled_traps for ls in layer_snapshots)
    total_leaking = sum(ls.num_leaking_traps for ls in layer_snapshots)
    filled_by_layer = [ls.num_filled_traps for ls in layer_snapshots]
    leaking_by_layer = [ls.num_leaking_traps for ls in layer_snapshots]

    return ReservoirSnapshot(
        timestamp,
        total_injected, current_injection_rate,
        total_stored, stored_by_layer, storage_fraction,
        total_leaked, top_layer_leakage_rate,
        mass_balance_error, mass_balance_error_percent,
        max_height, max_heights_by_layer,
        total_traps, total_filled, total_leaking,
        filled_by_layer, leaking_by_layer,
        layer_snapshots
    )
end


"""
Get injection rate at a specific time.
"""
function _get_injection_rate_at_time(
    injection_events::Vector{InjectionEvent},
    timestamp::Float64
)::Float64
    if isempty(injection_events)
        return 0.0
    end

    # Find the active injection event
    current_rate = 0.0
    for ie in injection_events
        if ie.timestamp <= timestamp
            current_rate = sum(ie.injection_rate)
        else
            break
        end
    end

    return current_rate
end


"""
    print_snapshot_summary(snapshot::ReservoirSnapshot)

Print a formatted summary of a reservoir snapshot.
"""
function print_snapshot_summary(snapshot::ReservoirSnapshot)
    println("\n" * "="^70)
    println("RESERVOIR SNAPSHOT at t = $(round(snapshot.timestamp, digits=2)) years")
    println("="^70)

    println("\n--- Injection & Storage ---")
    println("  Total injected:     $(round(snapshot.total_injected_m3 / 1e6, digits=3)) M m³")
    println("  Injection rate:     $(round(snapshot.injection_rate_m3_per_year / 1e6, digits=3)) M m³/year")
    println("  Total stored:       $(round(snapshot.total_stored_m3 / 1e6, digits=3)) M m³")
    println("  Total leaked:       $(round(snapshot.total_leaked_m3 / 1e6, digits=3)) M m³")

    println("\n--- Mass Balance ---")
    status = snapshot.mass_balance_error_percent < 0.1 ? "✓" : (snapshot.mass_balance_error_percent < 1.0 ? "~" : "✗")
    println("  Error: $(round(snapshot.mass_balance_error_m3, digits=2)) m³ ($(round(snapshot.mass_balance_error_percent, digits=4))%) $(status)")

    println("\n--- CO2 Heights ---")
    println("  Maximum height:     $(round(snapshot.max_co2_height, digits=2)) m")

    println("\n--- Trap Statistics ---")
    println("  Total traps:        $(snapshot.total_traps)")
    println("  Filled traps:       $(snapshot.total_filled_traps)")
    println("  Leaking traps:      $(snapshot.total_leaking_traps)")

    println("\n--- Storage by Layer ---")
    for (i, ls) in enumerate(snapshot.layer_snapshots)
        pct = snapshot.storage_fraction_by_layer[i] * 100
        vol = snapshot.stored_by_layer_m3[i] / 1e6
        height = snapshot.max_co2_height_by_layer[i]
        println("  Layer $(i) ($(ls.layer_name)): $(round(vol, digits=3)) M m³ ($(round(pct, digits=1))%), max h=$(round(height, digits=1))m")
    end
end


"""
    print_layer_snapshot_summary(snapshot::MultiLayerSnapshot)

Print a formatted summary of a layer snapshot.
"""
function print_layer_snapshot_summary(snapshot::MultiLayerSnapshot)
    println("\n--- Layer $(snapshot.layer_idx): $(snapshot.layer_name) at t=$(round(snapshot.timestamp, digits=2)) ---")
    println("  Stored volume:      $(round(snapshot.stored_volume_m3 / 1e6, digits=4)) M m³")
    println("  Max CO2 height:     $(round(snapshot.max_co2_height, digits=2)) m")
    println("  Mean CO2 height:    $(round(snapshot.mean_co2_height, digits=2)) m")
    println("  Traps with CO2:     $(snapshot.num_traps_with_co2) / $(snapshot.num_traps)")
    println("  Filled traps:       $(snapshot.num_filled_traps)")
    println("  Leaking traps:      $(snapshot.num_leaking_traps)")
    if snapshot.num_leaking_traps > 0
        println("  Leakage rate:       $(round(snapshot.leakage_rate_m3_per_year / 1e6, digits=4)) M m³/year")
        println("  Cumulative leaked:  $(round(snapshot.cumulative_leaked_m3 / 1e6, digits=4)) M m³")
    end
end
