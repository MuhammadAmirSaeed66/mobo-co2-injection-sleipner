# General R–Julia interface examples supplied with the simulator.
# These are demonstrations rather than the final paper experiment.

# CO2 Injection Simulation Example
# Demonstrates the improved R interface for CO2InjectionModeling.jl
library(JuliaCall)
source(file.path("R", "project_setup.R"))

setup_result <- julia_call("setup_simulator", boundary_condition = "closed")
print(setup_result)

# ============================================================================
# EXAMPLE 1: Simple simulation with Sleipner defaults
# ============================================================================
nx <- setup_result$nx
ny <- setup_result$ny
n_layers <- setup_result$n_layers

cat("\nUsing Sleipner default reservoir properties...\n")
config_result <- julia_call("setup_sleipner_reservoir")
print(config_result)

if (config_result$status != "success") {
  stop(paste("Config failed:", config_result$message))
}

cat("\nSetting up injection scenario...\n")

# Historical Sleipner injection rates (1996–2010), Mt/year
rates_mt <- c(
  0.07, 0.67, 0.85, 0.94, 0.94,
  1.02, 0.96, 0.92, 0.76, 0.87,
  0.83, 0.93, 0.82, 0.86, 0.76
)

co2_density <- 570.0  # kg/m^3 at bottom layer
rates_m3 <- rates_mt * 1e9 / co2_density   # m^3/year
n_times <- length(rates_m3)

# Injection array: time × nx × ny for layer 1
layer1_injection <- array(0, dim = c(n_times, nx+2, ny+2))

# Inject at (32, 59)
for (i in 1:n_times) {
  layer1_injection[i, 32, 59] <- rates_m3[i]
}

# Zero injection for other layers: 1 × nx × ny
zero_injection <- array(0, dim = c(1, nx+2, ny+2))

# 9 layers total
injection_matrices <- list(
  layer1_injection,  # layer 1
  zero_injection,    # layer 2
  zero_injection,    # layer 3
  zero_injection,    # layer 4
  zero_injection,    # layer 5
  zero_injection,    # layer 6
  zero_injection,    # layer 7
  zero_injection,    # layer 8
  zero_injection     # layer 9
)

cat("\nRunning simulation...\n")

sim_result <- julia_call(
  "run_simulation",
  start_time = 0.0,
  end_time = 15.0,
  time_step = 1.0,
  injection_rate_matrices = injection_matrices,
  verbose = FALSE
)

if (sim_result$status == "success") {
  cat("\n=== Simulation Successful! ===\n")
  cat("Timepoints:", sim_result$timepoints, "\n")
  cat("Total CO2 volumes (m³):\n")
  print(sim_result$total_co2_volumes)
  cat("\nFinal total volume:",
      tail(sim_result$total_co2_volumes, 1),
      "m³\n")
} else {
  cat("\nSimulation failed:", sim_result$message, "\n")
  if (!is.null(sim_result$stacktrace)) {
    cat("Stacktrace:\n", sim_result$stacktrace, "\n")
  }
}

# ============================================================================
# EXAMPLE 2: Custom reservoir properties with multiple injection wells
# ============================================================================

cat("\n\n=== Example 2: Custom Properties with Multiple Wells ===\n\n")

# Setup simulator (reuse from above)
cat("Configuring custom reservoir properties...\n")

# Use layer-specific density differences
brine_density <- 1020
co2_densities <- c(570, 542.5, 515, 487.5, 460, 432.5, 405, 377.5, 350)
density_diffs <- brine_density - co2_densities

custom_config <- julia_call("configure_reservoir",
                            porosity = rep(0.35, 9),
                            residual_co2_sat = rep(0.25, 9),
                            irreducible_water_sat = rep(0.3, 9),
                            shale_pressure_threshold = rep(98000.0, 9),
                            brine_co2_density_diff = density_diffs,
                            residual_leakage_time = rep(1.0, 9),
                            layer_specific = TRUE)
print(custom_config)

# Create injection scenario with two wells
cat("\nSetting up two-well injection scenario...\n")

n_times <- 10
well1_rates <- seq(0.5, 1.0, length.out = n_times) * 1e9 / 570  # Ramping up
well2_rates <- rep(0.8, n_times) * 1e9 / 570  # Constant rate

# Layer 1 injection with two wells
layer1_injection_2wells <- array(0, dim = c(n_times, nx+2, ny+2))

for (i in 1:n_times) {
  layer1_injection_2wells[i, 32, 59] <- well1_rates[i]  # Well 1
  layer1_injection_2wells[i, 35, 62] <- well2_rates[i]  # Well 2
}

# Build injection matrices
injection_matrices_2wells <- list(
  layer1_injection_2wells,
  zero_injection, zero_injection, zero_injection, zero_injection,
  zero_injection, zero_injection, zero_injection, zero_injection
)

# Run simulation
cat("\nRunning two-well simulation...\n")
sim_result_2wells <- julia_call("run_simulation",
                                start_time = 0.0,
                                end_time = 10.0,
                                time_step = 1.0,
                                injection_rate_matrices = injection_matrices_2wells,
                                verbose = FALSE)

if (sim_result_2wells$status == "success") {
  cat("\n=== Two-Well Simulation Successful! ===\n")
  cat("Final total volume:", tail(sim_result_2wells$total_co2_volumes, 1), "m³\n")
} else {
  cat("\nTwo-well simulation failed:", sim_result_2wells$message, "\n")
}

cat("\n=== All examples completed! ===\n")

# Example 1: total CO2 vs time
df_tot1 <- data.frame(
  time   = sim_result$timepoints,
  volume = sim_result$total_co2_volumes
)
library(ggplot2)
## Time series of total CO₂ volume – Example 1
ggplot(df_tot1, aes(x = time, y = volume)) +
  geom_line() +
  geom_point() +
  labs(
    title = "Example 1: Total CO2 Volume in Reservoir",
    x = "Time [years]",
    y = expression("Total CO"[2] * " volume [m"^3*"]")
  ) +
  theme_minimal()

## Time series – Example 1 vs Example 2

df_tot2 <- data.frame(
  time   = sim_result_2wells$timepoints,
  volume = sim_result_2wells$total_co2_volumes
)

df_tot1$case <- "Example 1: 1 well, Sleipner defaults"
df_tot2$case <- "Example 2: 2 wells, custom"

df_both <- rbind(df_tot1, df_tot2)

ggplot(df_both, aes(x = time, y = volume, color = case)) +
  geom_line() +
  geom_point() +
  labs(
    title = "Total CO2 Stored vs Time: Comparison of Cases",
    x = "Time [years]",
    y = expression("Total CO"[2] * " volume [m"^3*"]"),
    color = "Scenario"
  ) +
  theme_minimal()

## Plotting historical injection rates (Mt/year, Example 1)
years <- 1996:(1996 + length(rates_mt) - 1)

df_inj1 <- data.frame(
  year    = years,
  rate_Mt = rates_mt,
  rate_m3 = rates_m3
)

# Bar plot of Mt/year
ggplot(df_inj1, aes(x = year, y = rate_Mt)) +
  geom_col() +
  labs(
    title = "Sleipner Historical Injection Rates",
    x = "Year",
    y = expression("Injection rate [Mt CO"[2]*"/year]")
  ) +
  theme_minimal()

## Cumulative injected vs simulated stored (Example 1)
# Cumulative injected volume (using the same time base as rates)
cum_injected <- cumsum(rates_m3)

df_cum <- data.frame(
  time          = years - years[1],   # time in years since start, 0,1,...
  cum_injected  = cum_injected
)

# Interpolate simulation volumes at integer years to compare
stored_interp <- approx(
  x = sim_result$timepoints,
  y = sim_result$total_co2_volumes,
  xout = df_cum$time
)$y

df_cum$stored_sim <- stored_interp

# Plot both
ggplot(df_cum, aes(x = time)) +
  geom_line(aes(y = cum_injected, linetype = "Cumulative injected")) +
  geom_line(aes(y = stored_sim,   linetype = "Simulated stored"), linewidth = 0.8) +
  labs(
    title = "Cumulative Injected vs Simulated Stored CO\u2082",
    x = "Time since start [years]",
    y = expression("CO"[2] * " volume [m"^3*"]"),
    linetype = ""
  ) +
  theme_minimal()

## Two-well injection schedule (Example 2)

df_inj2 <- data.frame(
  time_years   = 0:(n_times - 1),
  well1_rate_m3 = well1_rates,
  well2_rate_m3 = well2_rates
)

# Convert to long format for ggplot
df_inj2_long <- rbind(
  data.frame(time_years = df_inj2$time_years,
             rate_m3 = df_inj2$well1_rate_m3,
             well = "Well 1"),
  data.frame(time_years = df_inj2$time_years,
             rate_m3 = df_inj2$well2_rate_m3,
             well = "Well 2")
)

ggplot(df_inj2_long, aes(x = time_years, y = rate_m3, color = well)) +
  geom_line() +
  geom_point() +
  labs(
    title = "Two-Well Injection Schedule (Example 2)",
    x = "Time [years]",
    y = expression("Injection rate [m"^3*"/year]"),
    color = "Well"
  ) +
  theme_minimal()

## Tables
## Historical injection schedule (Example 1)

library(dplyr)

tbl_injection_hist <- df_inj1 %>%
  mutate(
    cum_rate_Mt = cumsum(rate_Mt),
    cum_rate_m3 = cumsum(rate_m3)
  )

tbl_injection_hist

knitr::kable(
  tbl_injection_hist,
  digits = c(0, 2, 0, 0),
  caption = "Historical Sleipner injection schedule used as model input"
)

## Overall mass balance for Example 1
# Use df_cum from script
tbl_mass_balance <- df_cum %>%
  mutate(
    time_years = time,  # rename for clarity
    retention_fraction = stored_sim / cum_injected
  )

# Show just first and last row as summary
tbl_mass_balance_summary <- tbl_mass_balance %>%
  slice(1, n()) %>%
  select(
    time_years,
    cum_injected,
    stored_sim,
    retention_fraction
  )

tbl_mass_balance_summary

knitr::kable(
  tbl_mass_balance_summary,
  digits = c(1, 0, 0, 3),
  col.names = c(
    "Time since start [years]",
    "Cum. injected [m³]",
    "Simulated stored [m³]",
    "Retention fraction [-]"
  ),
  caption = "Cumulative injected vs simulated stored CO2 – start and end of simulation (Example 1)"
)

## Scenario comparison table (Example 1 vs Example 2)
# Total injected volumes
total_injected_ex1_m3 <- sum(rates_m3)
total_injected_ex2_m3 <- sum(df_inj2$well1_rate_m3 + df_inj2$well2_rate_m3)

# Final stored volumes
final_stored_ex1_m3 <- tail(sim_result$total_co2_volumes, 1)
final_stored_ex2_m3 <- tail(sim_result_2wells$total_co2_volumes, 1)

tbl_scenarios <- data.frame(
  scenario = c(
    "Example 1: 1 well, Sleipner defaults",
    "Example 2: 2 wells, custom"
  ),
  n_wells = c(1, 2),
  sim_duration_years = c(
    max(sim_result$timepoints) - min(sim_result$timepoints),
    max(sim_result_2wells$timepoints) - min(sim_result_2wells$timepoints)
  ),
  total_injected_m3 = c(
    total_injected_ex1_m3,
    total_injected_ex2_m3
  ),
  final_stored_m3 = c(
    final_stored_ex1_m3,
    final_stored_ex2_m3
  ),
  storage_efficiency = c(
    final_stored_ex1_m3 / total_injected_ex1_m3,
    final_stored_ex2_m3 / total_injected_ex2_m3
  )
)

tbl_scenarios

knitr::kable(
  tbl_scenarios,
  digits = c(NA, 0, 1, 0, 0, 3),
  col.names = c(
    "Scenario",
    "No. of wells",
    "Sim. duration [years]",
    "Total injected [m³]",
    "Final stored [m³]",
    "Storage efficiency [-]"
  ),
  caption = "Comparison of injection scenarios: injected volume, stored volume, and storage efficiency"
)

## Time-series table of stored CO₂ (both cases)
df_ts1 <- data.frame(
  time_years = sim_result$timepoints,
  stored_m3_ex1 = sim_result$total_co2_volumes
)

df_ts2 <- data.frame(
  time_years = sim_result_2wells$timepoints,
  stored_m3_ex2 = sim_result_2wells$total_co2_volumes
)

# Join on time (assuming same time grid or overlapping subset)
df_ts_both <- full_join(df_ts1, df_ts2, by = "time_years") %>%
  arrange(time_years)

# Optionally thin out for a concise table
df_ts_both_thinned <- df_ts_both %>%
  filter(time_years %% 3 == 0 | time_years == max(time_years))

df_ts_both_thinned

knitr::kable(
  df_ts_both_thinned,
  digits = c(1, 0, 0),
  col.names = c(
    "Time [years]",
    "Stored CO2 – Example 1 [m³]",
    "Stored CO2 – Example 2 [m³]"
  ),
  caption = "Time evolution of stored CO2 for the two injection scenarios (subset of time steps)"
)

## Two-well injection schedule (Example 2)
tbl_injection_2wells <- df_inj2 %>%
  mutate(
    total_rate_m3 = well1_rate_m3 + well2_rate_m3
  )

tbl_injection_2wells

library(tidyr)

tbl_injection_2wells_long <- df_inj2_long %>%
  rename(time_years = time_years, rate_m3_per_year = rate_m3)

knitr::kable(
  tbl_injection_2wells,
  digits = c(0, 0, 0, 0),
  col.names = c(
    "Time since start [years]",
    "Well 1 rate [m³/year]",
    "Well 2 rate [m³/year]",
    "Total rate [m³/year]"
  ),
  caption = "Two-well injection schedule for Example 2"
)

library(tidyr)

tbl_injection_2wells_long <- df_inj2_long %>%
  rename(time_years = time_years, rate_m3_per_year = rate_m3)

knitr::kable(
  tbl_injection_2wells_long,
  digits = c(0, 0),
  col.names = c(
    "Time since start [years]",
    "Injection rate [m³/year]",
    "Well"
  ),
  caption = "Layer 1 injection schedule for each well in Example 2"
)

## reservoir configuration table
tbl_reservoir_layers <- data.frame(
  layer = 1:9,
  porosity = rep(0.35, 9),
  residual_co2_sat = rep(0.25, 9),
  irreducible_water_sat = rep(0.3, 9),
  shale_pressure_threshold_Pa = rep(98000, 9),
  co2_density_kg_m3 = co2_densities,
  brine_density_kg_m3 = brine_density,
  density_difference_kg_m3 = density_diffs
)

knitr::kable(
  tbl_reservoir_layers,
  digits = c(0, 2, 2, 2, 0, 1, 0, 1),
  caption = "Layer-specific reservoir properties used in the custom two-well simulation (Example 2)"
)

