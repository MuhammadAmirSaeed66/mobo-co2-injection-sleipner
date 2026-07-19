 # CO2 Injection Simulation Example
# Demonstrates the improved R interface for CO2InjectionModeling.jl

library(JuliaCall)

# Setup Julia and load package
julia_setup()

# REMEMBER TO SET WORKING DIR TO CO2InjectionModeling.jl
# Alternatively: Install the package globally if you want to call the lib from
# another project.
julia_command('using Pkg; Pkg.activate(".")')

# This line might take a bit of time...
julia_command('using CO2InjectionModeling')

# ============================================================================
# ADVANCED EXAMPLE: Custom leakage heights, multi-layer injection,
#                   and time-varying rates
# ============================================================================

cat("=== Advanced Example: Custom Reservoir Configuration ===\n\n")

# Step 1: Setup simulator
cat("Setting up simulator...\n")
setup_result <- julia_call("setup_simulator", boundary_condition = "closed") # "open" for open or "closed" for closed BCs
print(setup_result)


# Note: use the nx_after_bc and ny_after_bc as the closed boundary conditions pads the domain by one
nx <- setup_result$nx_after_bc
ny <- setup_result$ny_after_bc
n_layers <- setup_result$n_layers

# Step 2: Configure custom reservoir properties with layer-specific leakage heights
cat("\nConfiguring custom reservoir properties with layer-specific leakage heights...\n")

# Compute layer-specific leakage heights from density differences
# This shows how different CO2 densities at different depths affect leakage behavior
brine_density <- 1020  # kg/m³
co2_densities <- rep(425, 9)  # kg/m³, decreasing with depth
density_diffs <- brine_density - co2_densities
g <- 9.81  # m/s²
shale_threshold <- 98000.0  # Pa

# Calculate leakage heights: h = P_threshold / (Δρ * g)
# Lower density difference -> higher critical height before leakage
leakage_heights <- shale_threshold / (density_diffs * g)

# Make top layer impermeable (represents caprock)
leakage_heights[9] <- Inf

# Optionally, make layer 5 have a lower leakage height
# This demonstrates how to create heterogeneous vertical permeability
leakage_heights[5] <- leakage_heights[5] * 0.5  # 50% of normal

cat("Layer-specific leakage heights (m):\n")
print(round(leakage_heights, 2))

# Better contol over the spesific reservoir-parameters
config_result <- julia_call("configure_reservoir",
                             porosity = rep(0.4, 9),
                             residual_co2_sat = rep(0.2, 9),
                             irreducible_water_sat = rep(0.3, 9),
                             shale_pressure_threshold = rep(shale_threshold, 9),
                             leakage_height = leakage_heights,
                             residual_leakage_time = rep(1.0, 9),
                             layer_specific = TRUE)
print(config_result)

# Step 3: Setup advanced injection scenario
cat("\nSetting up multi-layer, multi-location injection scenario...\n")

# This example demonstrates:
# 1. Multiple injection wells at different locations
# 2. Injection in multiple layers simultaneously
# 3. Time-varying injection rates (ramping up, constant, ramping down)

n_times <- 15  # 20 years of injection
co2_density <- 425.0  # kg/m³ # By now this is a "must". Sorry...

# Define time-varying injection patterns
time_vector <- 0:(n_times-1)

# Well 1: Ramp up pattern (Layer 1, bottom)
well1_rates_mt <- pmin(0.5 + 0.05 * time_vector, 1.0)  # Ramps from 0.5 to 1.0 Mt/year
well1_rates_m3 <- well1_rates_mt * 1e9 / co2_density

# Well 2: Constant rate (Layer 3)
well2_rates_mt <- rep(0.8, n_times)  # Constant 0.8 Mt/year
well2_rates_m3 <- well2_rates_mt * 1e9 / co2_density

cat("\nInjection well configuration:\n")
cat("  Well 1 (Layer 1): Ramp-up pattern, 0.5 -> 1.0 Mt/year\n")
cat("  Well 2 (Layer 1): Constant rate, 0.8 Mt/year\n")

# Layer 1:
layer1_injection <- array(0, dim = c(n_times, nx, ny))
for (i in 1:n_times) {
  layer1_injection[i, nx%/%2, ny%/%2] <- well1_rates_m3[i]      # Well 1: center
}

# Layer 3:
layer3_injection <- array(0, dim = c(n_times, nx, ny))
for (i in 1:n_times) {
  layer3_injection[i, nx%/%2 + 5, ny%/%2 + 5] <- well2_rates_m3[i]  # Well 2: offset
}

# No injection array
zero_injection <- array(0, dim = c(1, nx, ny))

# Build list of injection matrices (one per layer)
injection_matrices <- list(
  layer1_injection,  # Layer 1 (bottom)
  zero_injection,    # Layer 2
  layer3_injection,  # Layer 3
  zero_injection,    # Layer 4
  zero_injection,    # Layer 5
  zero_injection,    # Layer 6
  zero_injection,    # Layer 7
  zero_injection,    # Layer 8
  zero_injection     # Layer 9 (top)
)

# Step 4: Run simulation
cat("\nRunning advanced simulation (this may take a moment)...\n")
sim_result <- julia_call("run_simulation",
                         start_time = 0.0,
                         end_time = 15.0,  # 20 years
                         time_step = 1.0,
                         injection_rate_matrices = injection_matrices,
                         verbose = FALSE)

# Step 5: Check results
if (sim_result$status == "success") {
  cat("\n=== Simulation Successful! ===\n")
  cat("Timepoints:", sim_result$timepoints, "\n")
  cat("Total CO2 volumes (m³):\n")
  print(sim_result$total_co2_volumes)
  cat("\nFinal total volume:", tail(sim_result$total_co2_volumes, 1), "m³\n")
} else {
  cat("\nSimulation failed:", sim_result$message, "\n")
  if (!is.null(sim_result$stacktrace)) {
    cat("Stacktrace:\n", sim_result$stacktrace, "\n")
  }
}

# Plot 1: Total CO2 volume over time
cat("\nPlotting total CO2 volume over time...\n")
timepoints <- sim_result$timepoints
total_volumes <- sim_result$total_co2_volumes
plot(timepoints, total_volumes / 1e6, type = "b",
     xlab = "Time (years)", ylab = "Total CO2 Volume (million m³)",
     main = "Total CO2 Volume Over Time (Multi-Well Scenario)",
     col = "blue", pch = 19)
grid()

# Plot 2: Layer-wise CO2 volumes highlighting multi-layer injection
cat("\nPlotting layer-wise CO2 volumes...\n")
layer_volumes <- sim_result$layer_co2_volumes
matplot(timepoints, layer_volumes / 1e6, type = "b",
        xlab = "Time (years)", ylab = "CO2 Volume per Layer (million m³)",
        main = "Layer-wise CO2 Volumes (Layers 1, 3, 5 have injection)",
        col = rainbow(n_layers), pch = 19, lty = 1)
legend("topleft", legend = paste("Layer", 1:n_layers),
       col = rainbow(n_layers), pch = 19, lty = 1, cex = 0.8)
grid()

# Plot 3: Injection rates over time for all wells
cat("\nPlotting injection rate schedules...\n")
plot(time_vector, well1_rates_mt, type = "l", lwd = 2, col = "red",
     xlab = "Time (years)", ylab = "Injection Rate (Mt/year)",
     main = "Time-Varying Injection Rates by Well",
     ylim = c(0, 1.2))
lines(time_vector, well2_rates_mt, lwd = 2, col = "blue")
legend("topright",
       legend = c("Well 1 (L1, ramp-up)", "Well 2 (L3, constant)"),
       col = c("red", "blue", "green", "purple"), lwd = 2, cex = 0.8)
grid()

# Step 6: Generate animations
cat("\nGenerating animation (this may take a moment)...\n")
anim_result <- julia_call("generate_birdseye_animation",
                          output_file = "co2_advanced_animation.gif",
                          num_frames = 30L,
                          fps = 2L,
                          max_CO2_height = 20.0)

if (anim_result$status == "success") {
  cat("Animation saved to:", anim_result$output_file, "\n")
} else {
  cat("Animation failed:", anim_result$message, "\n")
  if (!is.null(anim_result$stacktrace)) {
    cat("Stacktrace:\n", anim_result$stacktrace, "\n")
  }
}

