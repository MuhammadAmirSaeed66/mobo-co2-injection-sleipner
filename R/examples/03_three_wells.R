# Three-well R–Julia simulator demonstration.
# The main paper uses the staged schedule defined in R/01_run_mobo_final.R.

library(JuliaCall)
source(file.path("R", "project_setup.R"))

cat("=== Advanced Example: Custom Reservoir Configuration (3-location injection, multi-layer) ===\n\n")

# ===============================
# Step 1: Setup simulator
# ===============================
cat("Setting up simulator...\n")
setup_result <- julia_call("setup_simulator", boundary_condition = "closed") # "open" / "closed"

print(setup_result)

nx        <- setup_result$nx
ny        <- setup_result$ny
n_layers  <- setup_result$n_layers
nx_bc     <- nx + 2  # padded for closed BC
ny_bc     <- ny + 2

# ===============================
# Step 2: Reservoir properties
# ===============================
cat("\nConfiguring custom reservoir properties with layer-specific leakage heights...\n")

brine_density   <- 1020        # kg/m³
co2_densities   <- rep(425, 9) # kg/m³
density_diffs   <- brine_density - co2_densities
g               <- 9.81        # m/s²
shale_threshold <- 98000.0     # Pa

# Leakage heights: h = P_threshold / (Δρ * g)
leakage_heights <- shale_threshold / (density_diffs * g)

# All nine numerical layers are retained as numerical reservoir intervals
# in the final paper experiment. No layer is disabled here.

# Example: make layer 5 "weaker" (lower leakage height)
leakage_heights[5] <- leakage_heights[5] * 0.5

cat("Layer-specific leakage heights (m):\n")
print(round(leakage_heights, 2))

# Here we just keep them constant for simplicity.
config_result <- julia_call(
  "configure_reservoir",
  porosity                 = rep(0.4, 9),
  residual_co2_sat         = rep(0.2, 9),
  irreducible_water_sat    = rep(0.3, 9),
  shale_pressure_threshold = rep(shale_threshold, 9),
  brine_co2_density_diff   = density_diffs,
  residual_leakage_time    = rep(1.0, 9),
  layer_specific           = TRUE
)

print(config_result)

# ===============================
# Step 3: New injection scenario
#       3 wells with different
#       injection-rate patterns
#       and multi-layer injection
# ===============================
cat("\nSetting up 3-location, multi-layer injection scenario...\n")

n_times     <- 15                # 15 time steps (≈ 15 years)
co2_density <- 425.0             # kg/m³
time_vector <- 0:(n_times - 1)

# Injection rate patterns (Mt/year, between 0.5 and 1.0) ---

# Well A: original linear ramp 0.5 -> 1.0
wellA_rates_mt <- pmin(0.5 + 0.05 * time_vector, 1.0)

# Well B: faster ramp 0.5 -> 1.0 (reaches 1 earlier)
wellB_rates_mt <- pmin(0.5 + 0.08 * time_vector, 1.0)

# Well C: delayed ramp, flat at 0.5 for first 5 years, then ramps to 1.0
wellC_rates_mt <- ifelse(
  time_vector <= 4,
  0.5,
  pmin(0.5 + 0.125 * (time_vector - 4), 1.0)
)

# Convert to m³/year
wellA_rates_m3 <- wellA_rates_mt * 1e9 / co2_density
wellB_rates_m3 <- wellB_rates_mt * 1e9 / co2_density
wellC_rates_m3 <- wellC_rates_mt * 1e9 / co2_density

cat("\nInjection well configuration:\n")
cat("  Well A (loc1): variable ramp 0.5 -> 1.0 Mt/year (baseline pattern)\n")
cat("  Well B (loc2): faster ramp 0.5 -> 1.0 Mt/year\n")
cat("  Well C (loc3): delayed ramp 0.5 -> 1.0 Mt/year\n")
cat("  Each well injects into multiple layers to increase CO2 in the whole column.\n")

# ===============================
# Well locations
# ===============================
## Elling Settings
# ix_center <- nx_bc %/% 2
# iy_center <- ny_bc %/% 2
# 
# # Three lateral positions
# loc1 <- c(ix_center,           iy_center          ) # Well A
# loc2 <- c(ix_center + 5,       iy_center + 5      ) # Well B
# loc3 <- c(ix_center - 5,       iy_center - 5      ) # Well C

# ===============================
# Alternative well locations 1:
# ===============================
## “Spread-out corners” (maximize separation)
# This tends to reduce plume interaction between wells.

# ix_center <- nx_bc %/% 2
# iy_center <- ny_bc %/% 2
# 
# margin <- 4  # keep away from padded boundary
# 
# loc1 <- c(margin,           margin)            # Well A: lower-left-ish
# loc2 <- c(nx_bc - margin,   margin)            # Well B: lower-right-ish
# loc3 <- c(margin,           ny_bc - margin)    # Well C: upper-left-ish
# 
# cat("\n[Option 1] Well locations:\n")
# cat("  Well A:", loc1, "\n")
# cat("  Well B:", loc2, "\n")
# cat("  Well C:", loc3, "\n")

# ===============================
# Alternative well locations 2:
# ===============================
## “Clustered triangle” (strong interaction)
# All three wells inject near each other but not on the same cell. This usually increases pressure/plume overlap.

ix_center <- nx_bc %/% 2
iy_center <- ny_bc %/% 2

loc1 <- c(ix_center,       iy_center)          # Well A: center
loc2 <- c(ix_center + 2,   iy_center)          # Well B: close east
loc3 <- c(ix_center,       iy_center + 2)      # Well C: close north

cat("\n[Option 2] Well locations:\n")
cat("  Well A:", loc1, "\n")
cat("  Well B:", loc2, "\n")
cat("  Well C:", loc3, "\n")

# ===============================
# Alternative well locations 3:
# ===============================

# ix_center <- nx_bc %/% 2
# iy_center <- ny_bc %/% 2
# 
# loc1 <- c(ix_center - 8,   iy_center)          # Well A: west
# loc2 <- c(ix_center,       iy_center)          # Well B: center
# loc3 <- c(ix_center + 8,   iy_center)          # Well C: east
# 
# cat("\n[Option 3] Well locations:\n")
# cat("  Well A:", loc1, "\n")
# cat("  Well B:", loc2, "\n")
# cat("  Well C:", loc3, "\n")

# ===============================
# Alternative well locations 4:
# ===============================
## “Near boundary vs interior” (boundary-condition sensitivity)
# Useful to see how “closed” vs “open” BC influences results when injection is near an edge.
# ix_center <- nx_bc %/% 2
# iy_center <- ny_bc %/% 2
# 
# margin <- 3
# 
# loc1 <- c(margin,           iy_center)          # Well A: near left boundary
# loc2 <- c(nx_bc - margin,   iy_center)          # Well B: near right boundary
# loc3 <- c(ix_center,        iy_center)          # Well C: center
# 
# cat("\n[Option 4] Well locations:\n")
# cat("  Well A:", loc1, "\n")
# cat("  Well B:", loc2, "\n")
# cat("  Well C:", loc3, "\n")

# ===============================
# Build injection arrays for ALL 9 layers
# Dimensions for each: (time, nx_bc, ny_bc)
# We distribute each well's rate into several layers.
# This way, CO2 is stored in many layers, not only the top.
# ===============================

injection_matrices <- vector("list", n_layers)
for (ell in seq_len(n_layers)) {
  injection_matrices[[ell]] <- array(0, dim = c(n_times, nx_bc, ny_bc))
}

# Well A (center): focus on upper/mid layers 1–3
fracA <- c(
  L1 = 0.6,
  L2 = 0.25,
  L3 = 0.15
)

# Well B (offset +): focus on mid/deeper layers 2,4,5
fracB <- c(
  L2 = 0.2,
  L4 = 0.4,
  L5 = 0.4
)

# Well C (offset -): focus on mid/deeper layers 3,6,7
fracC <- c(
  L3 = 0.2,
  L6 = 0.4,
  L7 = 0.4
)

for (i in seq_len(n_times)) {
  # Well A contributions
  rateA <- wellA_rates_m3[i]
  injection_matrices[[1]][i, loc1[1], loc1[2]] <- injection_matrices[[1]][i, loc1[1], loc1[2]] + fracA["L1"] * rateA
  injection_matrices[[2]][i, loc1[1], loc1[2]] <- injection_matrices[[2]][i, loc1[1], loc1[2]] + fracA["L2"] * rateA
  injection_matrices[[3]][i, loc1[1], loc1[2]] <- injection_matrices[[3]][i, loc1[1], loc1[2]] + fracA["L3"] * rateA
  
  # Well B contributions
  rateB <- wellB_rates_m3[i]
  injection_matrices[[2]][i, loc2[1], loc2[2]] <- injection_matrices[[2]][i, loc2[1], loc2[2]] + fracB["L2"] * rateB
  injection_matrices[[4]][i, loc2[1], loc2[2]] <- injection_matrices[[4]][i, loc2[1], loc2[2]] + fracB["L4"] * rateB
  injection_matrices[[5]][i, loc2[1], loc2[2]] <- injection_matrices[[5]][i, loc2[1], loc2[2]] + fracB["L5"] * rateB
  
  # Well C contributions
  rateC <- wellC_rates_m3[i]
  injection_matrices[[3]][i, loc3[1], loc3[2]] <- injection_matrices[[3]][i, loc3[1], loc3[2]] + fracC["L3"] * rateC
  injection_matrices[[6]][i, loc3[1], loc3[2]] <- injection_matrices[[6]][i, loc3[1], loc3[2]] + fracC["L6"] * rateC
  injection_matrices[[7]][i, loc3[1], loc3[2]] <- injection_matrices[[7]][i, loc3[1], loc3[2]] + fracC["L7"] * rateC
}

# No injection into layer 8 and 9 (caprock) for now
injection_matrices[[8]][] <- 0
injection_matrices[[9]][] <- 0  # keep caprock with zero injection

# ===============================
# Step 4: Run simulation
# ===============================
sim_result <- julia_call(
  "run_simulation",
  start_time              = 0.0,
  end_time                = 15.0,
  time_step               = 1.0,
  injection_rate_matrices = injection_matrices,
  verbose                 = FALSE
)

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

# ===============================
# Step 5: Plots
# ===============================

cat("\nPlotting total CO2 volume over time...\n")

timepoints    <- sim_result$timepoints
total_volumes <- sim_result$total_co2_volumes

# Plot 1: Total CO2 volume over time
plot(
  timepoints,
  total_volumes / 1e6,
  type = "b",
  xlab = "Time (years)",
  ylab = "Total CO2 Volume (million m³)",
  main = "Total CO2 Volume Over Time (3-Wells)",
  col  = "blue",
  pch  = 19
)
grid()

# Plot 2: Layer-wise CO2 volumes
layer_volumes <- sim_result$layer_co2_volumes
stopifnot(!is.null(layer_volumes), ncol(layer_volumes) == n_layers)

matplot(
  timepoints,
  layer_volumes / 1e6,
  type = "b",
  lty  = 1,
  pch  = 19,
  col  = rainbow(n_layers),
  xlab = "Time (years)",
  ylab = "CO2 Volume per Layer (million m³)",
  main = "Layer-wise CO2 Volumes (Multi-layer Injection)"
)
legend(
  "topleft",
  legend = paste("Layer", 1:n_layers),
  col    = rainbow(n_layers),
  pch    = 19,
  lty    = 1,
  cex    = 0.8
)
grid()

matplot(
  timepoints,
  layer_volumes / 1e6,
  type = "b",
  lty  = 1,
  pch  = 19,
  col  = rainbow(n_layers),
  xlab = "Time (years)",
  ylab = "CO2 Volume per Layer (million m³)",
  main = "Layer-wise CO2 Volumes (Three Wells, Multi-layer Injection)",
  cex.main = 2,   # title size
  cex.lab  = 1.5, # axis label size
  cex.axis = 1.6  # axis tick size
)

legend(
  "topleft",
  legend = paste("Layer", 1:n_layers),
  col    = rainbow(n_layers),
  pch    = 19,
  lty    = 1,
  cex    = 1.2   # match the first plot (optional)
)

# Plot 3: Injection rate patterns for each well
plot(
  time_vector,
  wellA_rates_mt,
  type = "l",
  lwd  = 2,
  col  = "red",
  xlab = "Time (years)",
  ylab = "Injection Rate (Mt/year)",
  main = "Injection Rate per Well (Mt/year)"
)
lines(time_vector, wellB_rates_mt, lwd = 2, col = "blue")
lines(time_vector, wellC_rates_mt, lwd = 2, col = "darkgreen")
grid()
legend(
  "bottomright",
  legend = c("Well A", "Well B", "Well C"),
  col    = c("red", "blue", "darkgreen"),
  lwd    = 2,
  cex    = 0.8
)

##############
## Cumulative total volume
plot(
  timepoints,
  cumsum(total_volumes) / 1e6,
  type = "l",
  lwd  = 2,
  xlab = "Time (years)",
  ylab = "Cumulative CO2 Volume (million m³·time-step)",
  main = "Cumulative CO2 Volume Over Time"
)
grid()

## Fraction of total CO₂ in each layer vs time (stacked area style)
layer_fraction <- layer_volumes / rowSums(layer_volumes)

matplot(
  timepoints,
  layer_fraction,
  type = "l",
  lty  = 1,
  lwd  = 2,
  xlab = "Time (years)",
  ylab = "Fraction of Total CO2",
  main = "Layer Fraction of Total CO2 Over Time"
)
legend(
  "right",
  legend = paste("Layer", 1:n_layers),
  col    = 1:n_layers,
  lty    = 1,
  lwd    = 2,
  cex    = 0.7
)
grid()

## Barplot of final volume per layer
final_volumes <- layer_volumes[nrow(layer_volumes), ] / 1e6  # last time step

barplot(
  final_volumes,
  names.arg = paste("L", 1:n_layers),
  xlab = "Layer",
  ylab = "CO2 Volume (million m³)",
  main = "Final CO2 Volume per Layer",
  col  = rainbow(n_layers)
)
grid(nx = NA, ny = NULL)

## pie chart
pie(
  final_volumes,
  labels = paste("L", 1:n_layers),
  main   = "Layer Contribution to Final Total CO2",
  col    = rainbow(n_layers)
)

## Total injection rate vs time (all wells combined)
total_injection_mt <- wellA_rates_mt + wellB_rates_mt + wellC_rates_mt

plot(
  time_vector,
  total_injection_mt,
  type = "l",
  lwd  = 2,
  xlab = "Time (years)",
  ylab = "Total Injection (Mt/year)",
  main = "Total Injection Rate (All Wells)"
)
grid()

## Per-layer injection fractions for each well (static barplots)
par(mfrow = c(1, 3))

barplot(
  fracA,
  main = "Well A Injection Fractions",
  xlab = "Layer",
  ylab = "Fraction",
  ylim = c(0, 1),
  col  = "lightblue"
)
grid(nx = NA, ny = NULL)

barplot(
  fracB,
  main = "Well B Injection Fractions",
  xlab = "Layer",
  ylab = "Fraction",
  ylim = c(0, 1),
  col  = "lightgreen"
)
grid(nx = NA, ny = NULL)

barplot(
  fracC,
  main = "Well C Injection Fractions",
  xlab = "Layer",
  ylab = "Fraction",
  ylim = c(0, 1),
  col  = "lightpink"
)
grid(nx = NA, ny = NULL)

par(mfrow = c(1, 1))

## Total injection per layer per well (m³)
# Helper to sum over time at a given location
sum_injection_at <- function(inj_list, ix, iy) {
  sapply(inj_list, function(arr) sum(arr[, ix, iy]))
}

injA_layer <- sum_injection_at(injection_matrices, loc1[1], loc1[2])
injB_layer <- sum_injection_at(injection_matrices, loc2[1], loc2[2])
injC_layer <- sum_injection_at(injection_matrices, loc3[1], loc3[2])

well_injections <- rbind(
  WellA = injA_layer,
  WellB = injB_layer,
  WellC = injC_layer
)

colnames(well_injections) <- paste("L", 1:n_layers)
print(well_injections)

## Compact multi-plot dashboards
par(mfrow = c(2, 2))

# 1: total volume
plot(
  timepoints, total_volumes / 1e6,
  type = "b", pch = 19,
  xlab = "Time (years)",
  ylab = "Total CO2 (million m³)",
  main = "Total CO2 Volume"
)
grid()

# 2: top 3 layers only
matplot(
  timepoints,
  layer_volumes[, 1:3] / 1e6,
  type = "l",
  lty  = 1,
  lwd  = 2,
  xlab = "Time (years)",
  ylab = "CO2 (million m³)",
  main = "Layers 1–3"
)
legend("topleft", legend = paste("L", 1:3), col = 1:3, lty = 1, lwd = 2, cex = 0.8)
grid()

# 3: injection per well
matplot(
  time_vector,
  cbind(wellA_rates_mt, wellB_rates_mt, wellC_rates_mt),
  type = "l",
  lty  = 1,
  lwd  = 2,
  xlab = "Time (years)",
  ylab = "Rate (Mt/year)",
  main = "Injection per Well"
)
legend("bottomright", legend = c("A", "B", "C"), col = 1:3, lty = 1, lwd = 2, cex = 0.8)
grid()

## Cumulative injection per well (Mt)
cumA <- cumsum(wellA_rates_mt)
cumB <- cumsum(wellB_rates_mt)
cumC <- cumsum(wellC_rates_mt)

# 4: cumulative injection
matplot(
  time_vector,
  cbind(cumA, cumB, cumC),
  type = "l",
  lty  = 1,
  lwd  = 2,
  xlab = "Time (years)",
  ylab = "Cumulative (Mt)",
  main = "Cumulative Injection"
)
legend("topleft", legend = c("A", "B", "C"), col = 1:3, lty = 1, lwd = 2, cex = 0.8)
grid()

par(mfrow = c(1, 1))

## Stacked area plot (approx) using polygon
## Stacked CO2 Volume by Layer (fixed)
cum_layer <- t(apply(layer_volumes, 1, cumsum)) / 1e6  # million m³

plot(
  NA,
  xlim = range(timepoints),
  ylim = c(0, max(cum_layer)),
  xlab = "Time (years)",
  ylab = "CO2 Volume (million m³)",
  main = "Stacked CO2 Volume by Layer"
)

cols <- rainbow(n_layers)

for (k in n_layers:1) {
  y_top    <- cum_layer[, k]
  y_bottom <- if (k == 1) rep(0, length(timepoints)) else cum_layer[, k - 1]
  
  polygon(
    c(timepoints, rev(timepoints)),
    c(y_bottom,  rev(y_top)),
    col    = cols[k],
    border = NA
  )
}

legend(
  "topleft",
  legend = paste("Layer", 1:n_layers),
  fill   = cols,
  cex    = 0.7,
  bty    = "n"
)
grid()

## Leakage heights – fixed
finite_heights <- leakage_heights
finite_heights[is.infinite(finite_heights)] <- NA

ylim_max <- max(finite_heights, na.rm = TRUE)

par(mfrow = c(1, 2))

barplot(
  finite_heights,
  names.arg = paste("L", 1:n_layers),
  xlab = "Layer",
  ylab = "Leakage Height (m)",
  main = "Leakage Heights",
  col  = "orange",
  ylim = c(0, ylim_max * 1.1)
)
grid(nx = NA, ny = NULL)

barplot(
  final_volumes,
  names.arg = paste("L", 1:n_layers),
  xlab = "Layer",
  ylab = "Final CO2 Volume (million m³)",
  main = "Final Volume per Layer",
  col  = "skyblue"
)
grid(nx = NA, ny = NULL)

par(mfrow = c(1, 1))

## safer heatmap using image() directly on layer_volumes 
## Alternative heatmap using mat directly
mat <- layer_volumes / 1e6

image(
  x = timepoints,
  y = 1:n_layers,
  z = mat,
  xlab = "Time (years)",
  ylab = "Layer",
  main = "CO2 Volume Heatmap (Layer vs Time)",
  axes = FALSE
)
axis(1)
axis(2, at = 1:n_layers, labels = 1:n_layers)
box()
grid(nx = NA, ny = NULL)

