# One-well R–Julia simulator smoke test.
# This example is not the complete MOBO experiment reported in the paper.

# CO2 Injection Simulation Example
# Demonstrates the improved R interface for CO2InjectionModeling.jl
library(JuliaCall)
source(file.path("R", "project_setup.R"))

# ============================================================================
# EXAMPLE 1: Simple simulation with Sleipner defaults
# ============================================================================

cat("=== Example 1: Simple Simulation with Sleipner Defaults ===\n\n")

# Step 1: Setup simulator
cat("Setting up simulator...\n")
setup_result <- julia_call("setup_simulator", boundary_condition = "closed") # "open" for open or "closed" for closed BCs
print(setup_result)

# Note: use the nx_after_bc and ny_after_bc as the closed boundary conditions pads the domain by one
nx <- setup_result$nx
ny <- setup_result$ny
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
layer1_injection <- array(0, dim = c(n_times, nx+2, ny+2))

# Inject at center location. Just as an example.
for (i in 1:n_times) {
  layer1_injection[i, nx%/%2, ny%/%2] <- rates_m3[i]
}

# Create zero injection for other layers (1 × nx × ny)
zero_injection <- array(0, dim = c(1, nx+2, ny+2))

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

cat("Plot saved to 'total_co2_volume.png'\n")
dev.copy(png, filename = "total_co2_volume.png")
dev.off()

cat("\nPlotting layer-wise CO2 volumes over time...\n")
layer_volumes <- sim_result$layer_co2_volumes
matplot(timepoints, layer_volumes / 1e6, type = "b",
        xlab = "Time (years)", ylab = "CO2 Volume per Layer (million m³)",
        main = "Layer-wise CO2 Volumes Over Time",
        col = rainbow(n_layers), pch = 19, lty = 1)
legend("topleft", legend = paste("Layer", 1:n_layers), col = rainbow(n_layers), pch = 19, lty = 1)
grid()
cat("Plot saved to 'layerwise_co2_volumes.png'\n")
dev.copy(png, filename = "layerwise_co2_volumes.png")
dev.off()

########################################
## Second option:
# ============================================================================
# EXAMPLE 1C: Single-well simulation (custom location + custom rates + MULTI-layer injection)
# ============================================================================

cat("=== Example 1C: Single-Well (Custom Location + Custom Rates + Multi-layer Injection) ===\n\n")

# Step 1: Setup simulator
cat("Setting up simulator...\n")
setup_result <- julia_call("setup_simulator", boundary_condition = "closed") # "open" or "closed"
print(setup_result)

nx <- setup_result$nx
ny <- setup_result$ny
n_layers <- setup_result$n_layers

# Step 2: Use Sleipner default reservoir properties
cat("\nUsing Sleipner default reservoir properties...\n")
config_result <- julia_call("setup_sleipner_reservoir")
print(config_result)

# Step 3: Setup single-well multi-layer injection scenario (custom)
cat("\nSetting up single-well multi-layer injection scenario (custom)...\n")

# ---- Custom injection rates (Mt/year), 15 years ----
# Example: ramp up to 1.0, short plateau, then ramp down
rates_mt <- c(
  0.20, 0.35, 0.50, 0.65, 0.80,
  0.95, 1.00, 1.00, 0.95, 0.90,
  0.80, 0.70, 0.60, 0.55, 0.50
)

# Convert to m³/year (assumes CO2 density = 425 kg/m³)
co2_density <- 425.0
rates_m3 <- rates_mt * 1e9 / co2_density
n_times <- length(rates_m3)

# ---- Custom location (single well) ----
nx_bc <- nx + 2
ny_bc <- ny + 2

# Choose a different location (offset toward SW quadrant)
ix <- (nx_bc %/% 2) - 10
iy <- (ny_bc %/% 2) - 6

# Clamp to keep inside padded domain interior
ix <- max(2, min(nx_bc - 1, ix))
iy <- max(2, min(ny_bc - 1, iy))

cat("Injecting at location (ix, iy) =", ix, iy, "within", nx_bc, "x", ny_bc, "grid.\n")

# ---- Layer split (fractions must sum to 1.0) ----
# Inject into 3 layers from one well to distribute CO2 vertically
# Layer order: 1 bottom ... 9 top
fractions <- c(L2 = 0.50, L4 = 0.30, L6 = 0.20)

cat("Layer fractions:\n")
print(fractions)

# ---- Build injection matrices for ALL layers ----
injection_matrices <- vector("list", n_layers)
for (ell in seq_len(n_layers)) {
  injection_matrices[[ell]] <- array(0, dim = c(n_times, nx_bc, ny_bc))
}

# Apply the multi-layer injection at the single location
for (i in seq_len(n_times)) {
  rate <- rates_m3[i]
  for (nm in names(fractions)) {
    ell <- as.integer(sub("L", "", nm))  # "L4" -> 4
    injection_matrices[[ell]][i, ix, iy] <- injection_matrices[[ell]][i, ix, iy] + fractions[[nm]] * rate
  }
}

# want to keep caprock layers explicitly zeroed:
injection_matrices[[8]][] <- 0
injection_matrices[[9]][] <- 0

# Step 4: Run simulation
cat("\nRunning simulation (this may take a moment on first run)...\n")
sim_result <- julia_call(
  "run_simulation",
  start_time = 0.0,
  end_time = 15.0,
  time_step = 1.0,
  injection_rate_matrices = injection_matrices,
  verbose = FALSE
)

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

# Step 6: Plot total CO2 volume over time
cat("\nPlotting total CO2 volume over time...\n")
timepoints <- sim_result$timepoints
total_volumes <- sim_result$total_co2_volumes

plot(timepoints, total_volumes / 1e6, type = "b",
     xlab = "Time (years)", ylab = "Total CO2 Volume (million m³)",
     main = "Total CO2 Volume Over Time (Single Well, Multi-layer Injection)",
     col = "blue", pch = 19)
grid()

cat("Plot saved to 'total_co2_volume_singlewell_multilayer.png'\n")
dev.copy(png, filename = "total_co2_volume_singlewell_multilayer.png")
dev.off()

# Step 7: Plot layer-wise CO2 volumes
cat("\nPlotting layer-wise CO2 volumes over time...\n")
layer_volumes <- sim_result$layer_co2_volumes

matplot(timepoints, layer_volumes / 1e6, type = "b",
        xlab = "Time (years)", ylab = "CO2 Volume per Layer (million m³)",
        main = "Layer-wise CO2 Volumes (Single Well, Multi-layer Injection)",
        col = rainbow(n_layers), pch = 19, lty = 1)
legend("topleft", legend = paste("Layer", 1:n_layers),
       col = rainbow(n_layers), pch = 19, lty = 1)
grid()

cat("Plot saved to 'layerwise_co2_volumes_singlewell_multilayer.png'\n")
dev.copy(png, filename = "layerwise_co2_volumes_singlewell_multilayer.png")
dev.off()

matplot(timepoints, layer_volumes / 1e6, type = "b",
        xlab = "Time (years)", 
        ylab = "CO2 Volume per Layer (million m³)",
        main = "Layer-wise CO2 Volumes (Single Well, Multi-layer Injection)",
        col = rainbow(n_layers), pch = 19, lty = 1,
        cex.main = 2,   # title size
        cex.lab = 1.5,    # axis label size
        cex.axis = 1.6)   # axis tick size

legend("topleft", legend = paste("Layer", 1:n_layers),
       col = rainbow(n_layers), pch = 19, lty = 1,
       cex = 1.2)

