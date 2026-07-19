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
# EXAMPLE 1: Simple simulation with Sleipner defaults
# ============================================================================

cat("=== Example 1: Simple Simulation with Sleipner Defaults ===\n\n")

# Step 1: Setup simulator
cat("Setting up simulator...\n")
setup_result <- julia_call("setup_simulator", boundary_condition = "closed") # "open" for open or "closed" for closed BCs
print(setup_result)


# Note: use the nx_after_bc and ny_after_bc as the closed boundary conditions pads the domain by one
nx <- setup_result$nx_after_bc
ny <- setup_result$ny_after_bc
n_layers <- setup_result$n_layers

# Step 2: Use Sleipner default reservoir properties
cat("\nUsing Sleipner default reservoir properties...\n")
config_result <- julia_call("setup_sleipner_reservoir")
print(config_result)

# Step 3: Setup injection scenario
cat("\nSetting up injection scenario...\n")

# Historical Sleipner injection rates (1996-2010)
rates_mt <- c(0.07, 0.67, 0.85, 0.94, 0.94, 1.02,
              0.96, 0.92, 0.76, 0.87, 0.83, 0.93,
              0.82, 0.86, 0.76)  # Mt/year

# Convert to m³/year (Current implementation assumes a CO2 density of 425kg/m^3)
co2_density <- 425.0
rates_m3 <- rates_mt * 1e9 / co2_density

n_times <- length(rates_m3)

# Create injection matrix for layer 1 (bottom layer): n_times × nx × ny
layer1_injection <- array(0, dim = c(n_times, nx, ny))

# Inject at center location. Just as an example.
for (i in 1:n_times) {
  layer1_injection[i, nx%/%2, ny%/%2] <- rates_m3[i]
}

# Create zero injection for other layers (1 × nx × ny)
zero_injection <- array(0, dim = c(1, nx, ny))

# Build list of injection matrices (one per layer)
injection_matrices <- list(
  layer1_injection,  # Layer 1 (bottom)
  zero_injection,    # Layer 2
  zero_injection,    # Layer 3
  zero_injection,    # Layer 4
  zero_injection,    # Layer 5
  zero_injection,    # Layer 6
  zero_injection,    # Layer 7
  zero_injection,    # Layer 8
  zero_injection     # Layer 9 (top)
)

# Step 4: Run simulation
cat("\nRunning simulation (this may take a moment on first run)...\n")
sim_result <- julia_call("run_simulation",
                         start_time = 0.0,
                         end_time = 15.0,
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

# Plot the total CO2 volume over time using R-plotting
cat("\nPlotting total CO2 volume over time...\n")
timepoints <- sim_result$timepoints
total_volumes <- sim_result$total_co2_volumes
plot(timepoints, total_volumes / 1e6, type = "b",
     xlab = "Time (years)", ylab = "Total CO2 Volume (million m³)",
     main = "Total CO2 Volume Over Time",
     col = "blue", pch = 19)
grid()
# cat("Plot saved to 'total_co2_volume.png'\n")
# dev.copy(png, filename = "total_co2_volume.png")
# dev.off()

cat("\nPlotting layer-wise CO2 volumes over time...\n")
layer_volumes <- sim_result$layer_co2_volumes
matplot(timepoints, layer_volumes / 1e6, type = "b",
        xlab = "Time (years)", ylab = "CO2 Volume per Layer (million m³)",
        main = "Layer-wise CO2 Volumes Over Time",
        col = rainbow(n_layers), pch = 19, lty = 1)
legend("topleft", legend = paste("Layer", 1:n_layers), col = rainbow(n_layers), pch = 19, lty = 1)
grid()
# cat("Plot saved to 'layerwise_co2_volumes.png'\n")
# dev.copy(png, filename = "layerwise_co2_volumes.png")
# dev.off()

# Step 6: Generate animations
cat("\nGenerating animation (this may take a moment)...\n")
anim_result <- julia_call("generate_birdseye_animation",
                          output_file = "co2_animation.gif",
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


