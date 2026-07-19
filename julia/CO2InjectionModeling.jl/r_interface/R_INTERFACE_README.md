# Using CO2InjectionModeling.jl from R

Use the CO2 injection simulator from R using the JuliaCall package.

## Setup

### 1. Install prerequisites

- Julia 1.6+ from [julialang.org](https://julialang.org/downloads/)
- R 4.0+ from [r-project.org](https://www.r-project.org/)
- JuliaCall R package:
  ```R
  install.packages("JuliaCall")
  ```

### 2. Install CO2InjectionModeling.jl

#### Option A: Use from the project directory

Navigate to the project directory and run:
```R
library(JuliaCall)
julia_setup()

# REMEMBER TO SET WORKING DIR TO CO2InjectionModeling.jl
julia_command('using Pkg; Pkg.activate(".")')

# This might take some time
julia_command('using CO2InjectionModeling')
```

#### Option B: Install as a package

Install the package globally so you can use it from any directory:

```bash
cd /path/to/CO2InjectionModeling.jl
julia --project=. -e 'using Pkg; Pkg.instantiate()'
julia -e 'using Pkg; Pkg.develop(path=pwd())'
```

Then in R:
```R
library(JuliaCall)
julia_setup()
julia_command('using CO2InjectionModeling')
```

## Quick Start

The simplest way to run a simulation with Sleipner defaults:

```R
# 1. Setup simulator
setup <- julia_call("setup_simulator", boundary_condition = "open")
nx <- setup$nx
ny <- setup$ny

# 2. Use Sleipner default properties
julia_call("setup_sleipner_reservoir")

# 3. Create injection scenario
# Injection matrix dimensions: n_times × nx × ny
n_times <- 15
layer1_injection <- array(0, dim = c(n_times, nx, ny))

# Inject at location (32, 59) with historical Sleipner rates
rates_mt <- c(0.07, 0.67, 0.85, 0.94, 0.94, 1.02,
              0.96, 0.92, 0.76, 0.87, 0.83, 0.93,
              0.82, 0.86, 0.76)  # Mt/year
rates_m3 <- rates_mt * 1e9 / 570  # Convert to m³/year

for (i in 1:n_times) {
  layer1_injection[i, 32, 59] <- rates_m3[i]
}

# Zero injection for other layers
zero_injection <- array(0, dim = c(1, nx, ny))

# List of injection matrices (one per layer)
injection_matrices <- list(
  layer1_injection,  # Layer 1 (bottom)
  zero_injection, zero_injection, zero_injection,
  zero_injection, zero_injection, zero_injection,
  zero_injection, zero_injection  # Layer 9 (top)
)

# 4. Run simulation
result <- julia_call("run_simulation",
                     start_time = 0.0,
                     end_time = 15.0,
                     time_step = 1.0,
                     injection_rate_matrices = injection_matrices,
                     verbose = TRUE)

# 5. Access results
print(result$timepoints)
print(result$total_co2_volumes)
```

## API Reference

### `setup_simulator()`

Setup the simulator and load topography data.

**Parameters:**
- `data_path`: Path to depth surfaces (default: `"sleipner/depth_surfaces/"`)
- `boundary_condition`: `"open"` or `"closed"` (default: `"open"`)

**Returns:** Dictionary with `status`, `n_layers`, `nx`, `ny`, `boundary_condition`

---

### `setup_sleipner_reservoir()`

Configure reservoir properties using Sleipner field default values. This is a convenience function that automatically sets:
- Porosity: 0.4
- Residual CO2 saturation: 0.2
- Irreducible water saturation: 0.3
- Shale pressure threshold: 98000.0 Pa
- Leakage heights: Computed from density differences (brine: 1020 kg/m³, CO2: 425 kg/m³)
  - Approximately 16.8 m for all layers except top
  - Top layer (L9): Inf (impermeable caprock)
- Residual leakage time: 1.0 years

**Parameters:** None (must call `setup_simulator()` first)

**Returns:** Dictionary with `status`, `n_layers`, `message`

**Example:**
```R
julia_call("setup_simulator", boundary_condition = "open")
julia_call("setup_sleipner_reservoir")
```

---

### `configure_reservoir()`

Configure custom reservoir properties.

**Parameters:**
- `porosity`: Sand porosity (0-1). Scalar or vector of length n_layers
- `residual_co2_sat`: Residual CO2 saturation (0-1). Scalar or vector
- `irreducible_water_sat`: Irreducible water saturation (0-1). Scalar or vector
- `shale_pressure_threshold`: Shale pressure threshold (Pa). Scalar or vector
- `leakage_height`: Critical CO2 height for leakage through shale (m). Scalar or vector. Use `Inf` for impermeable caprock
- `residual_leakage_time`: Residual leakage time (years). Scalar or vector
- `layer_specific`: Set to `TRUE` to provide vectors for layer-specific properties (default: `FALSE`)

**Returns:** Dictionary with `status`, `n_layers`

**Example (uniform properties):**
```R
julia_call("configure_reservoir",
           porosity = 0.35,
           residual_co2_sat = 0.2,
           irreducible_water_sat = 0.3,
           shale_pressure_threshold = 98000.0,
           leakage_height = 17.0,
           residual_leakage_time = 1.0,
           layer_specific = FALSE)
```

**Example (layer-specific properties):**
```R
# Compute layer-specific leakage heights from density differences
brine_density <- 1020  # kg/m³
co2_densities <- rep(425, 9)  # kg/m³
density_diffs <- brine_density - co2_densities
g <- 9.81  # m/s²
shale_threshold <- 98000.0  # Pa

# Calculate leakage heights: h = P_threshold / (Δρ * g)
leakage_heights <- shale_threshold / (density_diffs * g)

# Make top layer impermeable (represents caprock)
leakage_heights[9] <- Inf

julia_call("configure_reservoir",
           porosity = rep(0.4, 9),
           residual_co2_sat = rep(0.2, 9),
           irreducible_water_sat = rep(0.3, 9),
           shale_pressure_threshold = rep(98000.0, 9),
           leakage_height = leakage_heights,
           residual_leakage_time = rep(1.0, 9),
           layer_specific = TRUE)
```

---

### `run_simulation()`

Run a CO2 injection simulation.

**Parameters:**
- `start_time`: Simulation start time (years)
- `end_time`: Simulation end time (years)
- `time_step`: Time step for output snapshots (years)
- `injection_rate_matrices`: List of 3D arrays (one per layer), each with dimensions `(n_times × nx × ny)`. Each array specifies injection rates (m³/year) at each grid cell for each time point. For layers without injection, provide a `(1 × nx × ny)` array of zeros.
- `verbose`: Print progress messages (default: `FALSE`)

**Returns:** Dictionary containing:
- `status`: "success" or "error"
- `timepoints`: Vector of snapshot times (years)
- `total_co2_volumes`: Total CO2 stored in the reservoir at each timepoint (m³)
- `layer_co2_volumes`: Matrix (timepoints × layers) of volumes per layer (m³)
- `num_layers`: Number of layers
- `num_traps_per_layer`: Vector of trap counts per layer

**Example (single injection well):**
```R
# Create injection for bottom layer
n_times <- 10
layer1_injection <- array(0, dim = c(n_times, nx, ny))

# Inject at (32, 59)
rates <- rep(0.8, n_times) * 1e9 / 570  # 0.8 Mt/year constant
for (i in 1:n_times) {
  layer1_injection[i, 32, 59] <- rates[i]
}

# Zero injection for other layers
zero_injection <- array(0, dim = c(1, nx, ny))

injection_matrices <- list(
  layer1_injection,
  zero_injection, zero_injection, zero_injection,
  zero_injection, zero_injection, zero_injection,
  zero_injection, zero_injection
)

result <- julia_call("run_simulation",
                     start_time = 0.0,
                     end_time = 10.0,
                     time_step = 1.0,
                     injection_rate_matrices = injection_matrices,
                     verbose = TRUE)
```

**Example (multiple injection wells with time-varying rates):**
```R
n_times <- 10
layer1_injection <- array(0, dim = c(n_times, nx, ny))

# Well 1: Ramping up from 0.5 to 1.0 Mt/year
well1_rates <- seq(0.5, 1.0, length.out = n_times) * 1e9 / 570

# Well 2: Constant 0.8 Mt/year
well2_rates <- rep(0.8, n_times) * 1e9 / 570

for (i in 1:n_times) {
  layer1_injection[i, 32, 59] <- well1_rates[i]  # Well 1
  layer1_injection[i, 35, 62] <- well2_rates[i]  # Well 2
}
```

---

### `generate_birdseye_animation()`

Generate a bird's eye view animation showing CO2 distribution across all layers over time.

**Parameters:**
- `output_file`: Path where to save the animation (default: `"multi_layer_filling.gif"`)
- `num_frames`: Number of frames in animation (default: `30L`)
- `start_time`: Start time for animation in years (default: `0.0`)
- `end_time`: End time for animation in years, or `NULL` for auto-detect (default: `NULL`)
- `fps`: Frames per second (default: `2L`)
- `colormap`: Colormap name for CO2 heights (default: `"thermal"`)
- `max_CO2_height`: Maximum CO2 height for colorscale in meters (default: `20.0`)

**Returns:** Dictionary with `status`, `output_file`, `message`

**Example:**
```R
# Must call run_simulation() first
julia_call("generate_birdseye_animation",
           output_file = "co2_animation.gif",
           num_frames = 30L,
           fps = 2L,
           max_CO2_height = 20.0)
```

**Note:** Both `generate_cross_section_animation()` and `generate_birdseye_animation()` currently generate the same bird's eye view animation. Cross-section views are not yet available.

---

## Result Structure

All functions return a dictionary with a `status` field:
- `"success"`: Operation completed successfully
- `"error"`: Operation failed (check `message` field for details)

Simulation results include:
- `timepoints`: Times when snapshots were taken
- `total_co2_volumes`: Total CO2 stored in the reservoir at each timepoint
- `layer_co2_volumes`: CO2 volume breakdown by layer (timepoints × layers matrix)
- `num_layers`: Number of layers in the simulation
- `num_traps_per_layer`: Number of traps in each layer

## Examples

See [r_interface_example.R](r_interface_example.R) for complete examples demonstrating:
1. Simple simulation with Sleipner defaults
2. Custom reservoir properties with multiple injection wells

## Notes

- **Grid indexing**: R uses 1-based indexing, which matches the Julia interface
- **Array dimensions**: Injection matrices use dimension order `(n_times × nx × ny)`
- **Units**:
  - Injection rates: m³/year
  - Volumes: m³
  - Time: years
  - Pressure: Pa
  - Density: kg/m³
- **Performance**: First simulation run will be slower due to Julia JIT compilation. Subsequent runs are much faster.
