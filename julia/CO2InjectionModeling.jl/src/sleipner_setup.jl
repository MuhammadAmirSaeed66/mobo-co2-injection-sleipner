using NPZ
using Statistics
using CairoMakie
export SleipnerTopography, load_sleipner_topography, create_domain_from_topography
export generate_reservoir_properties_for_sleipner_layers, generate_sleipner_injection_events

struct SleipnerTopography
    surfaces::Dict{String, Any}
    top_caprock::Array{Float64,2}
    sand_layers::Vector{Dict{String, Any}}
    nx::Int
    ny::Int
    dx::Float64
    dy::Float64
    depth_min::Float64
    depth_max::Float64
end

function load_sleipner_topography(path::String = "sleipner/depth_surfaces/")
    println("\nLoading Sleipner depth surfaces...")

    # Load individual .npy files instead of .npz
    function load_surface(name::String)
        return npzread(joinpath(path, "$(name).npy"))
    end

    top_caprock = load_surface("Top_Caprock")
    top_sw = load_surface("TopSW")
    top_utsira = load_surface("TopUtsiraFm")
    base_utsira = load_surface("BaseUtsiraFm")
    thick_shale = load_surface("ThickShale")

    # Store all surfaces in a dictionary for compatibility
    surfaces = Dict{String, Any}(
        "Top_Caprock" => top_caprock,
        "TopSW" => top_sw,
        "TopUtsiraFm" => top_utsira,
        "BaseUtsiraFm" => base_utsira,
        "ThickShale" => thick_shale
    )

    sand_layers = []

    # L9: Shallowest sand (above thick shale)
    push!(sand_layers, Dict(
        "id" => 9,
        "name" => "L9",
        "top" => top_sw,
        "base" => thick_shale,
    ))

    # Load reflector surfaces
    for i in 1:7
        surfaces["Reflector$(i)"] = load_surface("Reflector$(i)")
        surfaces["Base_Reflector$(i)"] = load_surface("Base_Reflector$(i)")
    end

    # L8: First sand below thick shale
    push!(sand_layers, Dict(
        "id" => 8,
        "name" => "L8",
        "top" => top_utsira,
        "base" => surfaces["Reflector7"],
    ))

    # L7-L2: Sand layers between thin shales
    for i in 7:-1:2
        layer_num = i
        push!(sand_layers, Dict(
            "id" => layer_num,
            "name" => "L$(layer_num)",
            "top" => surfaces["Base_Reflector$(i)"],
            "base" => surfaces["Reflector$(i-1)"],
        ))
    end

    # L1: Deepest sand (below all thin shales)
    push!(sand_layers, Dict(
        "id" => 1,
        "name" => "L1",
        "top" => surfaces["Base_Reflector1"],
        "base" => base_utsira,
    ))

    # Sort by depth (shallowest to deepest by top surface mean)
    sand_layers = sort(sand_layers, by = x -> mean(x["top"]))

    nx, ny = size(top_caprock)
    depth_min = minimum(top_caprock)
    depth_max = maximum(base_utsira)

    # Bit hacky, but matches original grid spacing
    dx = 3200.0 / nx 
    dy = 5900.0 / ny

    return SleipnerTopography(
        surfaces,
        top_caprock,
        sand_layers,
        nx,
        ny,
        dx,
        dy,
        depth_min,
        depth_max,
    )
end

function create_domain_from_topography(topography::SleipnerTopography, dz::Float64)::Domain3D
    nx = topography.nx
    ny = topography.ny
    nz = Int(ceil((topography.depth_max - topography.depth_min) / dz))

    length_x = nx * topography.dx
    length_y = ny * topography.dy
    length_z = nz * dz

    Domain3D(
        nx,
        ny,
        nz,
        length_x,
        length_y,
        length_z,
        topography.depth_min,
        topography.depth_max,
    )
end


function generate_reservoir_properties_for_sleipner_layers()::Vector{ReservoirProperties}

    n_layers = 9

    # Common reservoir properties for all layers
    sand_porosity::Float64 = 0.4
    sand_residual_co2_saturation::Float64 = 0.2
    sand_irreducible_water_saturation::Float64 = 0.3
    shale_pressure_threshold::Float64 = 98000.0
    residual_leakage_time::Float64 = 2.0 # years

    # From L1 up to L9. Using values from paper.
    brine_density = 1020
    co2_densites = fill(425, n_layers) # Think it does not make sense with the current implementation of the mass-tracking to use different densities.
    brine_co2_density_differences = brine_density .- co2_densites

    # Compute leakage heights for each layer
    g = 9.81 # m/s^2
    leakage_heights = shale_pressure_threshold ./ (brine_co2_density_differences .* g)

    # Simulate impermeable caprock by setting very high leakage height at top layer
    leakage_heights[end] = Inf

    # Create ReservoirProperties for each layer
    reservoir_properties = Vector{ReservoirProperties}(undef, n_layers)
    for i in 1:n_layers
        reservoir_properties[i] = ReservoirProperties(
            sand_porosity,
            sand_residual_co2_saturation,
            sand_irreducible_water_saturation,
            shale_pressure_threshold,
            leakage_heights[i],
            residual_leakage_time
        )
    end
    return reservoir_properties
end


"""
    generate_sleipner_injection_events(layers, injection_cell::CartesianIndex)

Generate injection events for Sleipner simulation based on historical injection rates (1996-2010).
Injection occurs only in the bottom layer (L1) at the specified cell location.

# Arguments
- `layers`: Vector of layer structures
- `injection_cell`: CartesianIndex specifying the injection location in the bottom layer

# Returns
- `Vector{Vector{InjectionEvent}}`: Injection events for each layer (only L1 has non-zero injection)

# Notes
- Annual injection rates from Sleipner 2019 Benchmark (total 12.18 Mt over 14 years)
- Rates converted from Mt/year to m³/year using CO₂ density of 570 kg/m³
- Injection location is approximately at grid center (I=32, J=59) in the model
"""
function generate_sleipner_injection_events(
    layers,
    injection_cell::CartesianIndex
)::Vector{Vector{InjectionEvent}}

    # Historical annual injection rates from 1996-2010 (Mt/year)
    # Source: Sleipner 2019 Benchmark model
    annual_rates_mt = [
        0.07,  # 1996
        0.67,  # 1997
        0.85,  # 1998
        0.94,  # 1999
        0.94,  # 2000
        1.02,  # 2001
        0.96,  # 2002
        0.92,  # 2003
        0.76,  # 2004
        0.87,  # 2005
        0.83,  # 2006
        0.93,  # 2007
        0.82,  # 2008
        0.86,  # 2009
        0.76   # 2010
    ]

    # Convert Mt/year to m³/year
    co2_density_l1 = 425.0  # The same as used in reservoir properties
    annual_rates_m3_per_year = annual_rates_mt .* 1e9 ./ co2_density_l1

    # Create injection events for bottom layer (L1)
    # Time points are cumulative: start at year 0, events mark end of each year
    n_events = length(annual_rates_mt)
    bottom_layer_events = Vector{InjectionEvent}(undef, n_events)

    # Get the grid size from the bottom layer
    grid_size = size(layers[1].trap_structure.topography)

    for (i, rate) in enumerate(annual_rates_m3_per_year)
        # Time in years (0, 1, 2, ..., 14)
        time = float(i - 1)

        # Create injection rate field (only inject at specified cell)
        injection_rate = zeros(grid_size)
        injection_rate[injection_cell] = rate
        bottom_layer_events[i] = InjectionEvent(time, injection_rate)
    end

    # Create zero injection events for all other layers
    zero_injection = zeros(grid_size)
    zero_event = [InjectionEvent(0.0, zero_injection)]

    # Assemble injection events for all layers
    n_layers = length(layers)
    injection_events = Vector{Vector{InjectionEvent}}(undef, n_layers)
    injection_events[1] = bottom_layer_events  # Bottom layer (L1) has actual injection
    for i in 2:n_layers
        injection_events[i] = zero_event  # All other layers have zero injection
    end

    return injection_events
end