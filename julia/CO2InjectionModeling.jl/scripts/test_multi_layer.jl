using CO2InjectionModeling
using SurfaceWaterIntegratedModeling

# Setup similar to debugging_single_layer.jl
boundary_condition = :closed
topography = load_sleipner_topography()
domain = create_domain_from_topography(topography, 1.0)
layers = analyze_base_surfaces(topography; boundary_condition=boundary_condition)

trap_topo = layers[1].trap_structure.topography
xy = CartesianIndex(div(size(trap_topo, 1), 2), div(size(trap_topo, 2), 2))
injection_events = generate_sleipner_injection_events(layers, xy)

# Create reservoir properties with a VERY LOW leakage height to force leakage
layer_idx = 1
rprops = generate_reservoir_properties_for_sleipner_layers()

# Run fill_layer with leakage detection
seqs, leakage_states = fill_layers(
    layers,
    domain,
    rprops,
    injection_events;
    verbose=false
);

# Print drainage statistics for each layer
println("\n=== Leakage and Drainage Statistics ===")
for (layer_idx, leakage_state) in enumerate(leakage_states)
    tstruct = layers[layer_idx].trap_structure
    num_leaking = sum(leakage_state.leaking)
    num_draining = sum(leakage_state.draining)
    num_records = length(leakage_state.leakage_records)

    if num_records > 0
        total_initial_vol = sum(leakage_state.initial_volume_at_leak[i] for i in 1:numtraps(tstruct) if leakage_state.draining[i]; init=0.0)
        residual_sat = rprops[layer_idx].sand_residual_co2_saturation
        total_residual_vol = total_initial_vol * residual_sat

        println("Layer $layer_idx ($(layers[layer_idx].name)):")
        println("  Leaking traps: $(num_leaking)")
        println("  Draining traps: $(num_draining)")
        println("  Initial vol in draining: $(round(total_initial_vol, digits=1)) SWIM units")
        println("  Residual vol after drainage: $(round(total_residual_vol, digits=1)) SWIM units")
    else
        println("Layer $layer_idx ($(layers[layer_idx].name)): No leakage")
    end
end

# Generate multi-layer animations
println("\n=== Generating Multi-Layer Animations ===")
println("Generating height-based animation...")
animate_multi_layer_filling(
    layers,
    seqs,
    domain;
    output_file="multi_layer_filling.gif",
    num_frames=30,
    end_time=15.0,
    fps=3,
    max_CO2_height=rprops[1].leakage_height
)

println("\nGenerating saturation-based animation...")
animate_multi_layer_saturation(
    layers,
    seqs,
    leakage_states,
    domain,
    rprops;
    output_file="multi_layer_saturation.gif",
    num_frames=30,
    end_time=15.0,
    fps=3,
    colormap=:viridis
)

# Generate reservoir snapshots for analysis
println("\n=== Generating Reservoir Snapshots ===")
snapshots = generate_reservoir_snapshots(
    layers,
    seqs,
    leakage_states,
    domain,
    rprops,
    injection_events;
    num_snapshots=16,
    start_time=0.0,
    end_time=15.0,
    verbose=true
)

# Print summaries at key timepoints
println("\n=== Snapshot Summaries ===")
for snapshot in snapshots[1:4:end]  # Every 4th snapshot
    print_snapshot_summary(snapshot)
end

# Print detailed layer info for final snapshot
println("\n=== Final Snapshot - Layer Details ===")
final_snapshot = snapshots[end]
for layer_snapshot in final_snapshot.layer_snapshots
    print_layer_snapshot_summary(layer_snapshot)
end

# Generate timeseries plots
println("\n=== Generating Timeseries Plots ===")

# Plot CO2 volumes by layer (subplots)
println("Generating layer volumes plot...")
plot_layer_volumes_timeseries(
    snapshots;
    output_file="layer_volumes_timeseries.png",
    title="CO2 Volume by Layer Over Time"
)

# Plot CO2 fractions (stacked area chart)
println("Generating layer fractions plot...")
plot_layer_fractions_timeseries(
    snapshots;
    output_file="layer_fractions_timeseries.png",
    title="CO2 Distribution Across Layers Over Time"
)

println("\nTest completed!")
