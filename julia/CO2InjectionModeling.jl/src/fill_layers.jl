using SurfaceWaterIntegratedModeling
export fill_layers

function fill_layers(
        layers::Vector{Layer},
        domain::Domain3D,
        reservoir_properties::Union{ReservoirProperties, Vector{ReservoirProperties}},
        injection_events::Vector{Vector{InjectionEvent}};
        # num_snapshots::Int=10,
        # start_time::Float64=0.0,
        # end_time::Float64=15.0,
        verbose::Bool=false
)
    n_layers = length(layers)


    seqs = Vector{Vector{SpillEvent}}(undef, n_layers)
    leakage_states = Vector{LeakageState}(undef, n_layers)

    if reservoir_properties isa ReservoirProperties
        reservoir_properties = fill(reservoir_properties, n_layers)
    end

    weather_events_layer = convert_injection_event_to_weather_event(injection_events[1], reservoir_properties[1], domain)
    for layer_idx in 1:n_layers
        if verbose
            println("Filling layer $layer_idx / $n_layers")
        end

        tstruct = layers[layer_idx].trap_structure

        seq, leakage_state = fill_layer(
            tstruct,
            domain,
            reservoir_properties[layer_idx],
            weather_events_layer;
            verbose=verbose
        )

        seqs[layer_idx] = seq
        leakage_states[layer_idx] = leakage_state


        leakage_weather_events = nothing
        if !isempty(leakage_state.leakage_records)
            target_grid_size = size(tstruct.topography)
            leakage_weather_events = generate_leakage_weather_events(
                leakage_state,
                weather_events_layer,
                seq,
                tstruct,
                target_grid_size
            )
        end

        if layer_idx < n_layers
            weather_events_next_layer = convert_injection_event_to_weather_event(injection_events[layer_idx + 1], reservoir_properties[layer_idx + 1], domain)
            weather_events_layer = create_next_layer_weather_events(
                layers[layer_idx + 1].trap_structure,
                weather_events_next_layer,
                leakage_weather_events
            )

        end
    end

    return seqs, leakage_states
end

function create_next_layer_weather_events(
    next_tstruct::TrapStructure,
    next_layer_base_weather::Vector{WeatherEvent},
    leakage_weather_events::Union{Vector{WeatherEvent}, Nothing},
)::Vector{WeatherEvent}
    if isnothing(leakage_weather_events)
        return next_layer_base_weather
    end

    combined_weather_events = Vector{WeatherEvent}()

    # The next_layer_base_weather and leakage_weather_events must be sorted by time
    # Then we have to set the rain_rate to the sum of both events at each timepoint
    all_timestamps = Set{Float64}()
    push!(all_timestamps, 0.0)

    for we in next_layer_base_weather
        push!(all_timestamps, we.timestamp)
    end

    for we in leakage_weather_events
        push!(all_timestamps, we.timestamp)
    end


    sorted_timestamps = sort(collect(all_timestamps))
    nx, ny = size(next_tstruct.topography)
    for timestamp in sorted_timestamps
        # Start with base injection rate for this layer at this time
        rain_rate = zeros(Float64, nx, ny)

        # Add base injection from next_layer_base_weather
        for we in reverse(next_layer_base_weather)
            if we.timestamp <= timestamp
                if we.rain_rate isa Matrix
                    rain_rate .+= we.rain_rate
                end
                break
            end
        end

        for we in reverse(leakage_weather_events)
            if we.timestamp <= timestamp
                if we.rain_rate isa Matrix
                    rain_rate .+= we.rain_rate
                end
                break
            end
        end

        push!(combined_weather_events, WeatherEvent(timestamp, rain_rate))
    end

    return combined_weather_events
end