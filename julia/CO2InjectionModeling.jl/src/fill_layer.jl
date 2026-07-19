using SurfaceWaterIntegratedModeling
import Interpolations
using DifferentialEquations: solve, ODEProblem, VectorContinuousCallback, terminate!
export fill_layer, InjectionEvent, get_effective_leakage_cap


"""
    fill_layer(tstruct, domain, reservoir_properties, weather_events; kwargs...)

Fill a layer with CO2 based on weather events (injection events converted to SWIM format).
Detects and tracks leakage when CO2 height exceeds the leakage threshold.

Returns:
- `seq`: Vector of SpillEvents describing the filling sequence
- `leakage_state`: LeakageState tracking which traps are leaking and when
"""
function fill_layer(tstruct::TrapStructure{<:Real},
            domain::Domain3D,
            reservoir_properties::ReservoirProperties,
            weather_events::Vector{WeatherEvent};
            time_slack::Float64=0.0,
            infiltration::Union{Matrix{<:Real}, Nothing} = nothing,
            no_leakage::Bool=false,
            verbose::Bool=false)


    @assert !isempty(weather_events)

    num_traps = numtraps(tstruct)
    if num_traps == 0
        # No traps - return empty results
        empty_leakage = LeakageState(
            Bool[],  # leaking
            Bool[],  # draining
            Float64[], Float64[], LeakageRecord[],
            reservoir_properties.leakage_height,
            Float64[],  # initial_volume_at_leak
            reservoir_properties.sand_residual_co2_saturation,
            reservoir_properties.residual_leakage_time
        )
        return Vector{SpillEvent}(), empty_leakage
    end

    # Initialize infiltration map from user input
    infiltration =
        (typeof(infiltration) == Nothing) ? zeros(size(tstruct.topography)) :
        (typeof(infiltration) <: Real)  ? ones(size(tstruct.topography)) * infiltration :
                                          infiltration

    # Compute tables to support computation of trap water volume as function of water level
    z_vol_tables = SurfaceWaterIntegratedModeling._compute_z_vol_tables(tstruct)

    # Set initial filled_traps, cur_amounts and spillgraph
    filled_traps = Vector{Bool}(tstruct.trapvolumes .== 0.0)
    cur_amounts = fill(FilledAmount(0.0, weather_events[1].timestamp), num_traps)
    sgraph = SurfaceWaterIntegratedModeling.compute_complete_spillgraph(tstruct, filled_traps)

    # Start with empty sequence
    seq = Vector{SpillEvent}()

    # The height at which leakage occurs
    leakage_height = reservoir_properties.leakage_height
    if verbose
        println("Leakage height threshold: $(leakage_height) m")
    end

    # Initialize leakage state
    leakage_state = if no_leakage
        # Create a leakage state with infinite leakage volumes (no leakage possible)
        LeakageState(
            fill(false, num_traps),  # leaking
            fill(false, num_traps),  # draining
            fill(Inf, num_traps),
            fill(Inf, num_traps),
            LeakageRecord[],
            Inf,
            fill(0.0, num_traps),  # initial_volume_at_leak
            reservoir_properties.sand_residual_co2_saturation,  # residual_saturation
            reservoir_properties.residual_leakage_time  # residual_leakage_time
        )
    else
        initialize_leakage_state(
            tstruct, z_vol_tables, leakage_height,
            reservoir_properties.sand_residual_co2_saturation,
            reservoir_properties.residual_leakage_time
        )
    end

    # Compute development within the duration of each weather event
    for (wix, we) in enumerate(weather_events)
        cur_time = we.timestamp
        end_time =
            (wix == length(weather_events)) ? Inf : weather_events[wix+1].timestamp

        @assert(all([ca.time == cur_time for ca ∈ cur_amounts]))

        # Compute inflow/runoff/infiltration rates corresponding to the fill graph and new rain rate
        rateinfo = SurfaceWaterIntegratedModeling.compute_flow(sgraph, we.rain_rate, infiltration, tstruct, verbose)

        # Compute initial time estimates for when a trap become filled, or split into subtraps
        changetimeest = SurfaceWaterIntegratedModeling._set_initial_changetime_estimates(rateinfo, cur_amounts,
                                                          cur_time, filled_traps,
                                                          tstruct)

        # Compute initial leakage time estimates
        leakage_time_est = set_initial_leakage_time_estimates(
            cur_amounts, cur_time, rateinfo, leakage_state, filled_traps, tstruct
        )

        # Register the start of this weather event as a new, fully computed, spill event
        push!(seq, SpillEvent(cur_time, copy(cur_amounts), copy(filled_traps),
                              copy(rateinfo.trap_inflow), copy(we.rain_rate),
                              copy(rateinfo.runoff)))

        # Will add new events to `seq`. `sgraph`, `rateinfo`, `changetimeest`,
        # `filled_traps`, `cur_amounts`, and `leakage_state` are modified in the process
        _fill_sequence_for_weather_event_with_leakage!(
            seq, sgraph, rateinfo, changetimeest, leakage_time_est,
            filled_traps, cur_amounts, z_vol_tables,
            tstruct, infiltration, end_time, time_slack,
            leakage_state, verbose
        )
    end

    return seq, leakage_state
end


"""
    _fill_sequence_for_weather_event_with_leakage!(...)

Modified fill sequence that handles both fill/empty events and leakage events.
"""
function _fill_sequence_for_weather_event_with_leakage!(
    seq, sgraph, rateinfo, changetimeest, leakage_time_est,
    filled_traps, cur_amounts, z_vol_tables,
    tstruct, infiltration, endtime, time_slack,
    leakage_state::LeakageState, verbose
)
    cur_time = cur_amounts[1].time
    num_traps = numtraps(tstruct)

    fill_updates = Vector{IncrementalUpdate{Bool}}()
    graph_updates = Vector{IncrementalUpdate{Int}}()

    count = 0
    while cur_time < endtime
        verbose && (mod(count+=1, 10) == 0) && println("Fill sequence iteration: ", count)

        # Find next fill/empty event (from SWIM)
        next_fill_time, fill_updates = SurfaceWaterIntegratedModeling._identify_next_status_change!(
            changetimeest, cur_amounts, rateinfo,
            filled_traps, tstruct, z_vol_tables,
            cur_time, endtime
        )

        # Find next leakage event
        next_leak_time, leak_trap = find_next_leakage_event(leakage_time_est)

        # Determine which event comes first
        if next_leak_time < next_fill_time && next_leak_time <= endtime && leak_trap > 0
            # LEAKAGE EVENT occurs first
            cur_time = next_leak_time

            verbose && println("Leakage event at time $(cur_time) in trap $(leak_trap)")

            # Mark trap as leaking
            leakage_state.leaking[leak_trap] = true

            # Only mark as draining if trap has actual volume (not pass-through)
            # Pass-through traps (volume=0) just forward CO2 to the leak point
            leakage_vol_check = leakage_state.leakage_volume[leak_trap]
            if leakage_vol_check > 0
                leakage_state.draining[leak_trap] = true
                leakage_state.leakage_start_time[leak_trap] = cur_time
            end

            # Record the leakage for upstream layer
            leakage_location = find_leakage_location(leak_trap, tstruct)
            push!(leakage_state.leakage_records, LeakageRecord(
                cur_time,
                leak_trap,
                leakage_location
            ))

            # Record the initial volume at the time leakage started (for residual drainage)
            leakage_vol = leakage_state.leakage_volume[leak_trap]
            leakage_state.initial_volume_at_leak[leak_trap] = leakage_vol

            # NOTE: We do NOT mark ancestors as draining. Ancestors are DOWNSTREAM traps
            # (traps this trap spills into). Their CO2 does not flow through this trap.
            # They simply stop receiving spillover from this trap.
            #
            # Only DESCENDANTS should drain - they are UPSTREAM traps whose CO2 flows
            # INTO this trap and then out through the leak.
            #
            # IMPORTANT: Only mark descendants as draining if this trap has actual volume.
            # Pass-through traps (volume=0) don't cause drainage - they just pass CO2 through.
            # Drainage only happens when there's actual CO2 to drain.

            # Mark all filled descendants as draining too - their CO2 flows INTO this trap
            # and will drain out through the leak over the residual_leakage_time
            # Only do this if the leaking trap has actual volume (not a pass-through)
            if leakage_vol > 0
                descendants = get_all_descendants(tstruct, leak_trap)
                for desc_id in descendants
                    if filled_traps[desc_id] && !leakage_state.draining[desc_id]
                        leakage_state.draining[desc_id] = true
                        leakage_state.leakage_start_time[desc_id] = cur_time
                        # Get descendant's volume at leak time
                        desc_vol = cur_amounts[desc_id].amount
                        leakage_state.initial_volume_at_leak[desc_id] = desc_vol
                        verbose && println("  Descendant trap $(desc_id) marked as draining (vol=$(round(desc_vol, digits=2)))")
                    end
                end
            end

            # Cap the trap amount at leakage volume
            cur_amounts[leak_trap] = FilledAmount(leakage_vol, cur_time)

            # Mark ONLY the leaking trap as filled with edge=0
            # Ancestors can continue accumulating CO2 from other sources (e.g., other injection sites)
            # When ancestors spill to this trap, the CO2 will leak out through this trap's edge=0
            #
            # NOTE: We don't mark ancestors as leaking because:
            # 1. They may receive CO2 from multiple sources (different injection sites)
            # 2. Only the CO2 path through the leaking trap should leak
            # 3. Other CO2 plumes should continue accumulating normally

            filled_traps[leak_trap] = true
            leakage_state.leaking[leak_trap] = true

            verbose && println("  Trap $(leak_trap) marked as leaking, edge set to 0 (leak out of domain)")

            # Create fill update for the leaking trap only
            leak_fill_updates = [IncrementalUpdate{Bool}(leak_trap, true)]

            # Update spillgraph for the leaking trap
            graph_updates = SurfaceWaterIntegratedModeling.update_spillgraph!(sgraph, leak_fill_updates, tstruct)

            # Set ONLY the leaking trap's edge to 0 (out of domain = leakage)
            # Ancestors keep their normal edges - they spill to this trap, which then leaks
            sgraph.edges[leak_trap] = 0

            # Update flow information with the modified spillgraph
            setsavepoint!(rateinfo)
            SurfaceWaterIntegratedModeling._update_flow!(rateinfo, graph_updates, tstruct, sgraph)

            # Create amount update for this trap
            amount_updates = [IncrementalUpdate(leak_trap, FilledAmount(leakage_vol, cur_time))]

            # Update leakage time estimate to Inf for the leaking trap (already leaking)
            leakage_time_est[leak_trap] = LeakageTimeEstimate(leak_trap, Inf, Inf)

            # Update leakage time estimates for all traps whose inflow changed
            affected_traps = unique([u.index for u in getinflowupdates(rateinfo)])
            update_leakage_time_estimates!(
                leakage_time_est, affected_traps, cur_amounts, cur_time,
                rateinfo, leakage_state, filled_traps, tstruct
            )

            # Also update SWIM's changetime estimates for affected traps and the leaking trap
            all_traps_to_update = unique(vcat(affected_traps, [leak_trap]))
            for trap in all_traps_to_update
                changetimeest[trap] = SurfaceWaterIntegratedModeling._compute_changetime_estimate(
                    trap, cur_amounts, cur_time, rateinfo, filled_traps, tstruct
                )
            end

            # Record this as a spill event
            push!(seq, SpillEvent(cur_time, amount_updates, leak_fill_updates,
                                  getinflowupdates(rateinfo), nothing,
                                  getrunoffupdates(rateinfo)))

        elseif next_fill_time <= endtime && !isempty(fill_updates)
            # FILL/EMPTY EVENT occurs first (standard SWIM logic)
            cur_time = next_fill_time

            for u in fill_updates
                filled_traps[u.index] = u.value
            end

            # Given changes in fill state, update spill graph
            graph_updates = SurfaceWaterIntegratedModeling.update_spillgraph!(sgraph, fill_updates, tstruct)

            # Given the updates to the spill graph, update flow information in `rateinfo`
            setsavepoint!(rateinfo)
            SurfaceWaterIntegratedModeling._update_flow!(rateinfo, graph_updates, tstruct, sgraph)

            # Update water amount in traps whose inflow rate is about to change, or that just filled
            amount_updates = SurfaceWaterIntegratedModeling._update_affected_amounts(
                rateinfo, cur_amounts, filled_traps,
                tstruct, z_vol_tables, cur_time
            )

            # For traps that just filled, set their amount to full capacity
            # BUT: if they're leaking, cap at leakage volume instead
            # Also mark newly filled traps as draining if they feed into a draining trap
            for tix in [u.index for u in fill_updates]
                if leakage_state.leaking[tix]
                    # Leaking trap - cap at leakage volume (or 0 if pass-through)
                    cap_vol = get_effective_leakage_cap(leakage_state, tix)
                    push!(amount_updates, IncrementalUpdate(tix, FilledAmount(cap_vol, cur_time)))
                else
                    # Normal trap - fill to capacity
                    fill_vol = tstruct.trapvolumes[tix] - tstruct.subvolumes[tix]
                    push!(amount_updates, IncrementalUpdate(tix, FilledAmount(fill_vol, cur_time)))

                    # Check if this newly filled trap feeds into a draining chain
                    # If any descendant (child, grandchild, etc.) is draining, this trap should too
                    if !leakage_state.draining[tix]
                        # Check if this trap will spill into a draining trap
                        # Use spillgraph to find where this trap spills to
                        spill_target = sgraph.edges[tix]
                        # Only check if spill_target is a valid trap index (not runoff/boundary)
                        if spill_target > 0 && spill_target <= num_traps && leakage_state.draining[spill_target]
                            # This trap spills into a draining trap - it should also drain
                            leakage_state.draining[tix] = true
                            leakage_state.leakage_start_time[tix] = cur_time
                            leakage_state.initial_volume_at_leak[tix] = fill_vol
                            verbose && println("  Newly filled trap $(tix) marked as draining (feeds into draining trap $(spill_target))")
                        end
                    end
                end
            end

            # CRITICAL FIX: Also cap amounts for leaking traps that were updated by
            # _update_affected_amounts. This handles parent traps that are in "pass-through"
            # mode because their child is leaking. Without this fix, these parent traps
            # would incorrectly report their full capacity as stored volume.
            for i in 1:length(amount_updates)
                tix = amount_updates[i].index
                if leakage_state.leaking[tix]
                    cap_vol = get_effective_leakage_cap(leakage_state, tix)
                    if amount_updates[i].value.amount > cap_vol
                        amount_updates[i] = IncrementalUpdate(tix, FilledAmount(cap_vol, cur_time))
                    end
                end
            end

            # Integrate the changes into the continuously updated `cur_amounts` vector
            SurfaceWaterIntegratedModeling._apply_updates!(cur_amounts, amount_updates)

            # Update leakage time estimates for affected traps
            affected_traps = unique([u.index for u in getinflowupdates(rateinfo)])
            update_leakage_time_estimates!(
                leakage_time_est, affected_traps, cur_amounts, cur_time,
                rateinfo, leakage_state, filled_traps, tstruct
            )

            # Add current state to result
            push!(seq, SpillEvent(cur_time, amount_updates, fill_updates,
                                  getinflowupdates(rateinfo), nothing,
                                  getrunoffupdates(rateinfo)))
        else
            # No more events
            break
        end
    end

    # Make sure all amounts are exactly computed at end
    # Set all times to endtime (or keep as-is if endtime is Inf)
    # BUT: for leaking traps, account for residual drainage
    final_time = isfinite(endtime) ? endtime : cur_time

    for (trap, cur_fill) ∈ enumerate(cur_amounts)
        if cur_fill.time < endtime
            if leakage_state.leaking[trap]
                # Leaking trap - compute volume accounting for residual drainage
                # The stored volume decreases over time as CO2 drains out
                drained_vol = compute_volume_with_drainage(trap, final_time, leakage_state)
                # For leaking traps, drained_vol should never be nothing
                final_vol = isnothing(drained_vol) ? cur_fill.amount : drained_vol
                cur_amounts[trap] = FilledAmount(final_vol, final_time)
            else
                cur_amounts[trap] = FilledAmount(
                    SurfaceWaterIntegratedModeling._compute_exact_fill(
                        rateinfo, cur_amounts, trap,
                        filled_traps, tstruct, endtime,
                        z_vol_tables, false
                    ),
                    final_time
                )
            end
        end
    end
end


"""
    get_effective_leakage_cap(leakage_state, trap_id) -> Float64

Get the effective volume cap for a leaking trap.

Returns the leakage_volume for traps that are leaking (volume at which CO2 height
reaches the leakage threshold). Returns Inf for non-leaking traps.

NOTE: Only traps that actually started leaking are marked as leaking. Ancestor traps
are NOT marked as leaking - they continue accumulating normally and spill to the
leaking trap, which then passes CO2 out of domain. This allows multiple injection
sites to work correctly when trap hierarchies share common ancestors.
"""
function get_effective_leakage_cap(leakage_state::LeakageState, trap_id::Int)::Float64
    if !leakage_state.leaking[trap_id]
        # Not leaking - no cap (return Inf)
        return Inf
    end

    # This trap is leaking - cap at its leakage volume
    return leakage_state.leakage_volume[trap_id]
end


# Keep the old function signature for backwards compatibility
function _fill_sequence_for_weather_event!(seq, sgraph, rateinfo, changetimeest,
                                           filled_traps, cur_amounts, z_vol_tables,
                                           tstruct, infiltration, endtime, time_slack, leakage_height,
                                           verbose)
    # Create a dummy leakage state with no leakage
    num_traps = numtraps(tstruct)
    leakage_state = LeakageState(
        fill(false, num_traps),  # leaking
        fill(false, num_traps),  # draining
        fill(Inf, num_traps),
        fill(Inf, num_traps),
        LeakageRecord[],
        Inf,
        fill(0.0, num_traps),  # initial_volume_at_leak
        0.0,  # residual_saturation (no residual drainage)
        Inf   # residual_leakage_time (no residual drainage)
    )
    leakage_time_est = [LeakageTimeEstimate(i, Inf, Inf) for i in 1:num_traps]

    _fill_sequence_for_weather_event_with_leakage!(
        seq, sgraph, rateinfo, changetimeest, leakage_time_est,
        filled_traps, cur_amounts, z_vol_tables,
        tstruct, infiltration, endtime, time_slack,
        leakage_state, verbose
    )
end
