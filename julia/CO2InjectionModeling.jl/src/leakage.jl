"""
Leakage modeling for CO2 injection simulations.

This module provides functions to detect and track leakage of CO2 from reservoir traps
when the CO2 column height exceeds a threshold (leakage_height).
"""

import Interpolations
using SurfaceWaterIntegratedModeling: TrapStructure, numtraps, subtrapsof, FilledAmount

export compute_leakage_volume, initialize_leakage_state, find_leakage_location
export volume_to_height, compute_leakage_time_estimate, generate_leakage_weather_events
export get_true_topography_bottom, get_trap_bottom_elevation
export compute_drainable_volume, compute_volume_with_drainage, compute_residual_drainage_rate


"""
    find_leakage_location(trap_id::Int, tstruct::TrapStructure) -> CartesianIndex{2}

Find the grid cell where leakage occurs for a trap.
This is the lowest point (minimum topography elevation) in the trap footprint,
where the CO2 column is thickest.
"""
function find_leakage_location(trap_id::Int, tstruct::TrapStructure)::CartesianIndex{2}
    footprint = tstruct.footprints[trap_id]
    topo_vals = tstruct.topography[footprint]
    min_idx = argmin(topo_vals)
    linear_idx = footprint[min_idx]
    return CartesianIndices(size(tstruct.topography))[linear_idx]
end


"""
    get_true_topography_bottom(trap_id::Int, tstruct::TrapStructure) -> Float64

Get the TRUE minimum topography elevation in a trap's footprint.
This is the actual lowest point in the footprint, including the area covered by child traps.
This is the correct elevation to use for computing CO2 column height for leakage detection.

For parent traps, this includes the lowest point in all child trap footprints.
"""
function get_true_topography_bottom(trap_id::Int, tstruct::TrapStructure)::Float64
    footprint = tstruct.footprints[trap_id]
    return minimum(tstruct.topography[footprint])
end


"""
    get_trap_bottom_elevation(trap_id::Int, tstruct::TrapStructure) -> Float64

Get the effective bottom elevation of a trap for volume interpolation purposes.
For parent traps, this is the maximum of the topography minimum and the child spillpoint elevation.
For leaf traps, this is simply the minimum topography elevation in the footprint.

NOTE: This is used for z_vol_table interpolation, NOT for CO2 height calculation.
For height calculation (leakage detection), use get_true_topography_bottom instead.
"""
function get_trap_bottom_elevation(trap_id::Int, tstruct::TrapStructure)::Float64
    footprint = tstruct.footprints[trap_id]
    min_base_elevation = minimum(tstruct.topography[footprint])

    # For parent traps, the effective bottom is above child spillpoints
    children = subtrapsof(tstruct, trap_id)
    if !isempty(children)
        child_spillpoint_elev = tstruct.spillpoints[children[1]].elevation
        min_base_elevation = max(min_base_elevation, child_spillpoint_elev)
    end

    return min_base_elevation
end


"""
    compute_leakage_volume(trap_id::Int, z_vol_table, tstruct::TrapStructure, leakage_height::Float64) -> Union{Float64, Nothing}

Compute the volume at which a trap reaches the leakage height threshold.

The CO2 column height is measured from the TRUE topography bottom (including child traps)
to the current water level. This correctly accounts for CO2 that fills child traps
before spilling to parent traps.

Returns:
- Float64: The volume at which CO2 height = leakage_height
- nothing: If the trap would spill before reaching leakage_height (leakage cannot occur)
"""
function compute_leakage_volume(
    trap_id::Int,
    z_vol_table::Tuple{Vector{Float64}, Vector{Float64}},
    tstruct::TrapStructure,
    leakage_height::Float64
)::Union{Float64, Nothing}

    # Get the TRUE topography bottom (not the child spillpoint!)
    # This is the actual lowest point in the trap's footprint, including child areas
    true_bottom = get_true_topography_bottom(trap_id, tstruct)

    # Leakage occurs when water level reaches true_bottom + leakage_height
    leakage_elevation = true_bottom + leakage_height

    # Get spillpoint elevation (maximum fill level)
    spillpoint_elevation = tstruct.spillpoints[trap_id].elevation

    # If leakage elevation is above spillpoint, trap spills before it can leak
    if leakage_elevation >= spillpoint_elevation
        return nothing
    end

    # Use z_vol_table to find volume at leakage_elevation
    zvals, vvals = z_vol_table

    # Handle edge cases
    if length(zvals) == 1
        # Degenerate trap - no volume
        return 0.0
    end

    if leakage_elevation <= zvals[1]
        # Leakage elevation is below the trap's z_vol_table minimum.
        # For parent traps, zvals[1] is the child spillpoint.
        # This means leakage would occur while children are still filling,
        # so this trap's own volume at leakage is 0.
        # When this trap starts filling (children full), the CO2 column is already
        # above the leakage height, so leakage starts immediately.
        return 0.0
    end

    if leakage_elevation >= zvals[end]
        # Leakage height is above trap capacity - shouldn't happen given check above
        return nothing
    end

    # Create interpolation function from z to volume
    z2v = Interpolations.linear_interpolation(zvals, vvals, extrapolation_bc=Interpolations.Line())

    return z2v(leakage_elevation)
end


"""
    volume_to_height(volume::Float64, trap_id::Int, z_vol_table, tstruct::TrapStructure) -> Float64

Convert a volume in a trap to the CO2 column height above the TRUE topography bottom.

For parent traps with filled children, the CO2 column extends from the true
topography minimum (in the child footprints) through the children up to the
current water level. This function correctly computes this total column height.
"""
function volume_to_height(
    volume::Float64,
    trap_id::Int,
    z_vol_table::Tuple{Vector{Float64}, Vector{Float64}},
    tstruct::TrapStructure
)::Float64

    # Use the TRUE topography bottom for correct CO2 column height calculation
    true_bottom = get_true_topography_bottom(trap_id, tstruct)

    zvals, vvals = z_vol_table

    if length(zvals) == 1 || volume <= 0.0
        return 0.0
    end

    # Create interpolation function from volume to z
    v2z = Interpolations.linear_interpolation(vvals, zvals, extrapolation_bc=Interpolations.Line())

    water_level = v2z(volume)
    height = max(0.0, water_level - true_bottom)

    return height
end


"""
    compute_drainable_volume(initial_volume::Float64, residual_saturation::Float64) -> Float64

Compute the volume of CO2 that can drain from a trap during residual leakage.

The drainable volume is the portion that will leak out over the residual leakage time.
The residual (non-drainable) fraction remains trapped in pore spaces.

Parameters:
- `initial_volume`: Volume at the time leakage started
- `residual_saturation`: Fraction of CO2 that remains trapped (0 to 1)

Returns:
- Drainable volume (initial_volume * (1 - residual_saturation))
"""
function compute_drainable_volume(initial_volume::Float64, residual_saturation::Float64)::Float64
    return initial_volume * (1.0 - residual_saturation)
end


"""
    compute_volume_with_drainage(trap_id::Int, current_time::Float64, leakage_state::LeakageState) -> Union{Float64, Nothing}

Compute the current stored volume in a draining trap, accounting for residual drainage.

The volume decreases linearly from initial_volume to residual_volume over the
residual_leakage_time period after drainage starts.

A trap is "draining" if it's either:
1. Directly leaking (reached leakage threshold, edge=0), or
2. A filled ancestor whose CO2 flows through a leaking trap

Parameters:
- `trap_id`: ID of the trap
- `current_time`: Current simulation time
- `leakage_state`: Current leakage state

Returns:
- For draining traps: Current volume accounting for drainage
- For non-draining traps: `nothing` (caller should use actual volume from SWIM)

If residual_leakage_time is 0 or Inf, no drainage occurs (returns initial volume).
"""
function compute_volume_with_drainage(
    trap_id::Int,
    current_time::Float64,
    leakage_state::LeakageState
)::Union{Float64, Nothing}
    # If not draining, return nothing - caller should use actual volume from SWIM
    if !leakage_state.draining[trap_id]
        return nothing
    end

    # Get leakage parameters
    t_leak = leakage_state.leakage_start_time[trap_id]

    # If current_time is before leakage started, the trap wasn't leaking yet
    # Return nothing so caller uses actual volume from SWIM
    if current_time < t_leak
        return nothing
    end

    initial_vol = leakage_state.initial_volume_at_leak[trap_id]
    residual_sat = leakage_state.residual_saturation
    residual_time = leakage_state.residual_leakage_time

    # If no drainage (residual_time is Inf or 0, or residual_sat is 1)
    if !isfinite(residual_time) || residual_time <= 0.0 || residual_sat >= 1.0
        return initial_vol
    end

    # Compute time since leakage started
    time_since_leak = current_time - t_leak

    # Compute residual volume (what remains after full drainage)
    residual_vol = initial_vol * residual_sat

    # Compute drainable volume
    drainable_vol = initial_vol - residual_vol

    # Linear drainage over residual_leakage_time
    if time_since_leak >= residual_time
        # Drainage complete
        return residual_vol
    else
        # Partial drainage
        fraction_drained = time_since_leak / residual_time
        return initial_vol - drainable_vol * fraction_drained
    end
end


"""
    compute_residual_drainage_rate(trap_id::Int, current_time::Float64, leakage_state::LeakageState) -> Float64

Compute the current residual drainage rate from a leaking trap.

This is the rate at which CO2 is draining from the trap due to residual leakage
(separate from the pass-through of new injections).

Parameters:
- `trap_id`: ID of the trap
- `current_time`: Current simulation time
- `leakage_state`: Current leakage state

Returns:
- Drainage rate (volume per time unit). Returns 0 if not draining.
"""
function compute_residual_drainage_rate(
    trap_id::Int,
    current_time::Float64,
    leakage_state::LeakageState
)::Float64
    # If not leaking, no drainage
    if !leakage_state.leaking[trap_id]
        return 0.0
    end

    # Get leakage parameters
    t_leak = leakage_state.leakage_start_time[trap_id]
    initial_vol = leakage_state.initial_volume_at_leak[trap_id]
    residual_sat = leakage_state.residual_saturation
    residual_time = leakage_state.residual_leakage_time

    # If no drainage configured
    if !isfinite(residual_time) || residual_time <= 0.0 || residual_sat >= 1.0
        return 0.0
    end

    # Check if we're still in the drainage period
    time_since_leak = current_time - t_leak
    if time_since_leak < 0.0 || time_since_leak >= residual_time
        # Not yet started or already completed
        return 0.0
    end

    # Constant drainage rate during the drainage period
    drainable_vol = initial_vol * (1.0 - residual_sat)
    return drainable_vol / residual_time
end


"""
    initialize_leakage_state(tstruct::TrapStructure, z_vol_tables, leakage_height::Float64,
                             residual_saturation::Float64, residual_leakage_time::Float64) -> LeakageState

Initialize the leakage state for a layer, precomputing leakage volumes for all traps.

Parameters:
- `tstruct`: The trap structure for the layer
- `z_vol_tables`: Volume-elevation tables for each trap
- `leakage_height`: Height threshold at which leakage occurs
- `residual_saturation`: Fraction of CO2 that remains after drainage (sand_residual_co2_saturation)
- `residual_leakage_time`: Time over which residual drainage occurs
"""
function initialize_leakage_state(
    tstruct::TrapStructure,
    z_vol_tables::Vector{Tuple{Vector{Float64}, Vector{Float64}}},
    leakage_height::Float64,
    residual_saturation::Float64,
    residual_leakage_time::Float64
)::LeakageState

    num_traps = numtraps(tstruct)

    # Precompute leakage volumes for all traps
    leakage_volumes = zeros(Float64, num_traps)
    for trap_id in 1:num_traps
        vol = compute_leakage_volume(trap_id, z_vol_tables[trap_id], tstruct, leakage_height)
        # Use Inf for traps that cannot leak (spill before reaching leakage height)
        leakage_volumes[trap_id] = isnothing(vol) ? Inf : vol
    end

    return LeakageState(
        fill(false, num_traps),           # leaking
        fill(false, num_traps),           # draining
        leakage_volumes,                   # leakage_volume
        fill(Inf, num_traps),             # leakage_start_time
        LeakageRecord[],                   # leakage_records
        leakage_height,                    # leakage_height
        fill(0.0, num_traps),             # initial_volume_at_leak (0 until leakage starts)
        residual_saturation,               # residual_saturation
        residual_leakage_time              # residual_leakage_time
    )
end


"""
    LeakageTimeEstimate

Stores the estimated time when a trap will reach its leakage volume.
Similar to SWIM's ChangeTimeEstimate.
"""
struct LeakageTimeEstimate
    trap_id::Int
    min_time::Float64  # Earliest possible leakage time
    max_time::Float64  # Latest possible leakage time
end


"""
    compute_leakage_time_estimate(trap_id, cur_amounts, cur_time, rateinfo, leakage_state, filled_traps, tstruct) -> LeakageTimeEstimate

Compute when a trap will reach its leakage threshold.
"""
function compute_leakage_time_estimate(
    trap_id::Int,
    cur_amounts::Vector{FilledAmount},
    cur_time::Float64,
    rateinfo,  # SWIM's RateInfo
    leakage_state::LeakageState,
    filled_traps::Vector{Bool},
    tstruct::TrapStructure
)::LeakageTimeEstimate

    # Already leaking - no future leakage event
    if leakage_state.leaking[trap_id]
        return LeakageTimeEstimate(trap_id, Inf, Inf)
    end

    # Get target volume for leakage
    target_vol = leakage_state.leakage_volume[trap_id]

    # If target is Inf, trap cannot leak (spills before reaching leakage height)
    if target_vol == Inf
        return LeakageTimeEstimate(trap_id, Inf, Inf)
    end

    # Get current volume
    current_vol = cur_amounts[trap_id].amount

    # If trap is not accumulating (children not filled), cannot estimate leakage
    children = subtrapsof(tstruct, trap_id)
    if !all(filled_traps[children])
        return LeakageTimeEstimate(trap_id, Inf, Inf)
    end

    # Check if already at or above target
    # Special case: if target_vol is 0, we need current_vol > 0 to trigger
    # (can't leak if there's no CO2 yet)
    if target_vol == 0.0
        if current_vol > 0.0
            # Trap has CO2 and should leak immediately (already above threshold)
            return LeakageTimeEstimate(trap_id, cur_time, cur_time)
        end
        # target_vol = 0 but no CO2 yet - will leak as soon as it gets any CO2
        # Continue to compute time estimate based on inflow
    elseif current_vol >= target_vol
        return LeakageTimeEstimate(trap_id, cur_time, cur_time)
    end

    # Compute time using inflow rate
    # Use SWIM's RateInfo accessors
    inflow = SurfaceWaterIntegratedModeling.getinflow(rateinfo, trap_id)
    smax = SurfaceWaterIntegratedModeling.getsmax(rateinfo, trap_id)
    smin = SurfaceWaterIntegratedModeling.getsmin(rateinfo, trap_id)

    # Net inflow bounds (same logic as SWIM's _compute_changetime_estimate)
    min_net_inflow = inflow - (smax - smin)
    max_net_inflow = inflow

    # If inflow is zero or negative, trap won't fill to leakage level
    if max_net_inflow <= 0
        return LeakageTimeEstimate(trap_id, Inf, Inf)
    end

    volume_needed = target_vol - current_vol

    # Time estimates (min time uses max inflow, max time uses min inflow)
    min_time = volume_needed / max_net_inflow
    max_time = (min_net_inflow > 0) ? volume_needed / min_net_inflow : Inf

    return LeakageTimeEstimate(
        trap_id,
        cur_time + min_time,
        cur_time + max_time
    )
end


"""
    set_initial_leakage_time_estimates(cur_amounts, cur_time, rateinfo, leakage_state, filled_traps, tstruct) -> Vector{LeakageTimeEstimate}

Initialize leakage time estimates for all traps.
"""
function set_initial_leakage_time_estimates(
    cur_amounts::Vector{FilledAmount},
    cur_time::Float64,
    rateinfo,
    leakage_state::LeakageState,
    filled_traps::Vector{Bool},
    tstruct::TrapStructure
)::Vector{LeakageTimeEstimate}

    return [
        compute_leakage_time_estimate(trap_id, cur_amounts, cur_time, rateinfo,
                                      leakage_state, filled_traps, tstruct)
        for trap_id in 1:numtraps(tstruct)
    ]
end


"""
    update_leakage_time_estimates!(leakage_time_est, affected_traps, cur_amounts, cur_time, rateinfo, leakage_state, filled_traps, tstruct)

Update leakage time estimates for affected traps.
"""
function update_leakage_time_estimates!(
    leakage_time_est::Vector{LeakageTimeEstimate},
    affected_traps::Vector{Int},
    cur_amounts::Vector{FilledAmount},
    cur_time::Float64,
    rateinfo,
    leakage_state::LeakageState,
    filled_traps::Vector{Bool},
    tstruct::TrapStructure
)
    for trap_id in affected_traps
        leakage_time_est[trap_id] = compute_leakage_time_estimate(
            trap_id, cur_amounts, cur_time, rateinfo,
            leakage_state, filled_traps, tstruct
        )
    end
end


"""
    find_next_leakage_event(leakage_time_est) -> Tuple{Float64, Int}

Find the next leakage event (minimum leakage time).
Returns (time, trap_id). If no leakage is pending, returns (Inf, 0).
"""
function find_next_leakage_event(leakage_time_est::Vector{LeakageTimeEstimate})::Tuple{Float64, Int}
    min_time = Inf
    min_trap = 0

    for est in leakage_time_est
        if est.min_time < min_time
            min_time = est.min_time
            min_trap = est.trap_id
        end
    end

    return (min_time, min_trap)
end


"""
    generate_leakage_weather_events(leakage_state, source_weather_events, spill_events, tstruct, target_grid_size) -> Vector{WeatherEvent}

Generate WeatherEvents for the overlying layer from leakage data.

The leakage rate at any time t equals:
    leakage_rate(t) = injection_rate(t) - d(stored)/dt + residual_drainage_rate(t)

Where:
- d(stored)/dt is the rate at which non-leaking traps are filling
- residual_drainage_rate(t) is the rate at which leaking traps are draining their existing CO2

When multiple traps are leaking, the total leakage rate is distributed among all active
leakage locations proportionally to each trap's inflow rate. This preserves mass conservation:
the sum of rates at all leakage locations equals the total leakage rate.

When multiple traps are leaking, the total leakage rate is distributed among all active
leakage locations proportionally to each trap's inflow rate. This preserves mass conservation:
the sum of rates at all leakage locations equals the total leakage rate.

This approach correctly handles:
- The initial filling period (before leakage, leakage_rate = 0)
- The transition period (some traps filling, some leaking)
- Steady state (all contributing traps filled, leakage_rate = injection_rate)
- Multiple leakage locations with proper rate distribution
- Residual drainage (existing CO2 draining from leaking traps over time)
"""
function generate_leakage_weather_events(
    leakage_state::LeakageState,
    source_weather_events::Vector{WeatherEvent},
    spill_events::Vector{SpillEvent},
    tstruct::TrapStructure,
    target_grid_size::Tuple{Int, Int}
)::Vector{WeatherEvent}

    # If no leakage occurred, return empty vector
    if isempty(leakage_state.leakage_records)
        return WeatherEvent[]
    end

    # Get the first leakage start time
    first_leakage_time = minimum(r.start_time for r in leakage_state.leakage_records)

    # Collect all timestamps where leakage rates might change:
    # 1. Leakage start times
    # 2. Weather event timestamps (injection rates change)
    # 3. Spill event timestamps (fill rates change)
    # 4. Residual drainage end times (t_leak + residual_leakage_time)
    timestamps = Set{Float64}()

    for record in leakage_state.leakage_records
        push!(timestamps, record.start_time)
        # Add the time when residual drainage ends for this trap
        if isfinite(leakage_state.residual_leakage_time)
            drainage_end_time = record.start_time + leakage_state.residual_leakage_time
            push!(timestamps, drainage_end_time)
        end
    end

    for we in source_weather_events
        if we.timestamp >= first_leakage_time
            push!(timestamps, we.timestamp)
        end
    end

    for se in spill_events
        if se.timestamp >= first_leakage_time
            push!(timestamps, se.timestamp)
        end
    end

    # Sort timestamps
    sorted_times = sort(collect(timestamps))
    if isempty(sorted_times)
        return WeatherEvent[]
    end

    # Compute stored volume at each timestamp, accounting for drainage
    stored_at_time = Dict{Float64, Float64}()
    for t in sorted_times
        stored = compute_total_stored_at_time_with_drainage(spill_events, tstruct, t, leakage_state)
        stored_at_time[t] = stored
    end

    # Get injection rate at each timestamp (sum of rain_rate)
    function get_injection_rate_at(t::Float64)
        # Find the weather event active at time t
        current_we = source_weather_events[1]
        for we in source_weather_events
            if we.timestamp <= t
                current_we = we
            else
                break
            end
        end
        if current_we.rain_rate isa Matrix
            return sum(current_we.rain_rate)
        else
            # Scalar rain_rate - uniform across domain
            return current_we.rain_rate * prod(target_grid_size)
        end
    end

    # Build WeatherEvents using mass balance: leakage = injection - d(stored)/dt
    # Note: When drainage is happening, d(stored)/dt is negative, which increases leakage_rate
    events = WeatherEvent[]

    for i in 1:length(sorted_times)
        t = sorted_times[i]

        # Only generate events after leakage has started
        if t < first_leakage_time
            continue
        end

        # Compute filling rate (d(stored)/dt)
        # Use the interval from t to the next timestamp
        # When drainage is happening, this will be negative
        if i < length(sorted_times)
            t_next = sorted_times[i + 1]
            dt = t_next - t
            if dt > 0
                stored_now = stored_at_time[t]
                stored_next = stored_at_time[t_next]
                filling_rate = (stored_next - stored_now) / dt
            else
                filling_rate = 0.0
            end
        else
            # Last timestamp - check if we're still draining
            total_drainage_rate = 0.0
            for trap_id in 1:numtraps(tstruct)
                total_drainage_rate += compute_residual_drainage_rate(trap_id, t, leakage_state)
            end
            # If draining, filling_rate is negative (stored is decreasing)
            filling_rate = -total_drainage_rate
        end

        # Leakage rate = injection rate - filling rate
        # If filling_rate is negative (drainage), this increases leakage_rate
        injection_rate = get_injection_rate_at(t)
        total_leakage_rate = max(0.0, injection_rate - filling_rate)

        # Distribute leakage among all active leakage locations
        # Active = leakage has started (start_time <= t)
        active_records = filter(r -> r.start_time <= t, leakage_state.leakage_records)

        # Create rate matrix
        rate_matrix = zeros(Float64, target_grid_size)

        if total_leakage_rate > 0 && !isempty(active_records)
            # Get inflow rate to each active leaking trap to determine distribution weights
            inflow_rates = Float64[]
            for record in active_records
                inflow = get_trap_inflow_at_time(record.trap_id, t, spill_events)
                push!(inflow_rates, max(0.0, inflow))
            end

            total_inflow = sum(inflow_rates)

            if total_inflow > 0
                # Distribute proportionally to inflow rates (preserves mass conservation)
                for (j, record) in enumerate(active_records)
                    weight = inflow_rates[j] / total_inflow
                    loc = record.leakage_location
                    rate_matrix[loc] += weight * total_leakage_rate
                end
            else
                # Fallback: equal distribution if no inflow information
                equal_rate = total_leakage_rate / length(active_records)
                for record in active_records
                    loc = record.leakage_location
                    rate_matrix[loc] += equal_rate
                end
            end
        end

        # Add event if there's any leakage
        if sum(rate_matrix) > 0
            push!(events, WeatherEvent(t, rate_matrix))
        end
    end

    # Remove duplicate timestamps and merge rates
    return merge_weather_events(events)
end


"""
    compute_total_stored_at_time(spill_events, tstruct, time) -> Float64

Compute the total stored volume at a given time using the spill event sequence.
Does not account for residual drainage from leaking traps.
"""
function compute_total_stored_at_time(
    spill_events::Vector{SpillEvent},
    tstruct::TrapStructure,
    time::Float64
)::Float64
    # Use trap_states_at_timepoints for accurate computation
    tstates = trap_states_at_timepoints(tstruct, spill_events, [time]; verbose=false)
    volumes = tstates[1][2]  # Volume in each trap
    return sum(volumes)
end


"""
    compute_total_stored_at_time_with_drainage(spill_events, tstruct, time, leakage_state) -> Float64

Compute the total stored volume at a given time, accounting for residual drainage.

For draining traps (leaking traps and their descendants), the stored volume decreases
over time as CO2 drains out. This function uses compute_volume_with_drainage to get
the correct volume for each trap.
"""
function compute_total_stored_at_time_with_drainage(
    spill_events::Vector{SpillEvent},
    tstruct::TrapStructure,
    time::Float64,
    leakage_state::LeakageState
)::Float64
    # Get base trap states
    tstates = trap_states_at_timepoints(tstruct, spill_events, [time]; verbose=false)
    volumes = tstates[1][2]  # Volume in each trap

    total = 0.0
    for trap_id in 1:numtraps(tstruct)
        vol = volumes[trap_id]
        # Check 'draining' not 'leaking' - descendants of leaking traps also drain
        if leakage_state.draining[trap_id]
            # Draining trap - use drainage-adjusted volume if available
            drained_vol = compute_volume_with_drainage(trap_id, time, leakage_state)
            if !isnothing(drained_vol)
                vol = drained_vol
            end
        end
        total += vol
    end

    return total
end


"""
    get_trap_inflow_at_time(trap_id, time, spill_events) -> Float64

Get the inflow rate to a trap at a specific time from the spill event sequence.
"""
function get_trap_inflow_at_time(
    trap_id::Int,
    time::Float64,
    spill_events::Vector{SpillEvent}
)::Float64

    # Find the spill event at or just before this time
    se_idx = 1
    for i in 1:length(spill_events)
        if spill_events[i].timestamp <= time
            se_idx = i
        else
            break
        end
    end

    se = spill_events[se_idx]

    # Get inflow from spill event
    # The inflow field can be a full vector or incremental updates
    if se.inflow isa Vector{Float64}
        return se.inflow[trap_id]
    else
        # It's a vector of IncrementalUpdate - need to reconstruct
        # For now, return 0 if we can't find it directly
        # This case requires reconstructing full state from incremental updates
        for update in se.inflow
            if update.index == trap_id
                return update.value
            end
        end
        # Fall back - search earlier events
        return _reconstruct_inflow(trap_id, se_idx, spill_events)
    end
end


"""
    _reconstruct_inflow(trap_id, end_idx, spill_events) -> Float64

Reconstruct the inflow rate for a trap from spill events.
"""
function _reconstruct_inflow(
    trap_id::Int,
    end_idx::Int,
    spill_events::Vector{SpillEvent}
)::Float64

    # Start from the beginning and apply updates
    inflow = 0.0

    for i in 1:end_idx
        se = spill_events[i]
        if se.inflow isa Vector{Float64}
            inflow = se.inflow[trap_id]
        else
            for update in se.inflow
                if update.index == trap_id
                    inflow = update.value
                end
            end
        end
    end

    return inflow
end


"""
    merge_weather_events(events::Vector{WeatherEvent}) -> Vector{WeatherEvent}

Merge WeatherEvents that have the same timestamp.
"""
function merge_weather_events(events::Vector{WeatherEvent})::Vector{WeatherEvent}
    if isempty(events)
        return events
    end

    # Sort by timestamp
    sorted_events = sort(events, by=e -> e.timestamp)

    merged = WeatherEvent[]
    current_event = sorted_events[1]

    for i in 2:length(sorted_events)
        if sorted_events[i].timestamp == current_event.timestamp
            # Merge rates
            new_rate = current_event.rain_rate .+ sorted_events[i].rain_rate
            current_event = WeatherEvent(current_event.timestamp, new_rate)
        else
            push!(merged, current_event)
            current_event = sorted_events[i]
        end
    end
    push!(merged, current_event)

    return merged
end
