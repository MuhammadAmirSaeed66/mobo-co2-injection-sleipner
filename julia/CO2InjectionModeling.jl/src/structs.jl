using SurfaceWaterIntegratedModeling
export Domain3D, CellProperties, SimulationLayerSnapshot, SimulationSnapshot, ReservoirProperties
export LeakageRecord, LeakageState

struct Domain3D
    nx::Int
    ny::Int
    nz::Int
    length_x::Float64
    length_y::Float64
    length_z::Float64
    dx::Float64
    dy::Float64
    dz::Float64
    depth_min::Float64
    depth_max::Float64

    # Constructor with automatic calculation of cell sizes
    function Domain3D(nx, ny, nz, length_x, length_y, length_z, depth_min, depth_max)
        dx = length_x / nx
        dy = length_y / ny
        dz = length_z / nz
        new(nx, ny, nz, length_x, length_y, length_z, dx, dy, dz, depth_min, depth_max)
    end
end


struct InjectionEvent
    timestamp::Float64
    injection_rate::Union{Matrix{Float64}, Float64} 
end


struct CellProperties
    porosity::Array{Float64, 3}
    pressure_threshold::Array{Float64, 3}
    residual_co2_saturation::Array{Float64, 3}
    irreducible_water_saturation::Array{Float64, 3}

    function CellProperties(domain::Domain3D)
        nx, ny, nz = domain.nx, domain.ny, domain.nz
        new(zeros(nx, ny, nz), zeros(nx, ny, nz), zeros(nx, ny, nz), zeros(nx, ny, nz))
    end
end

struct ReservoirProperties
    sand_porosity::Float64
    sand_residual_co2_saturation::Float64
    sand_irreducible_water_saturation::Float64
    shale_pressure_threshold::Float64
    leakage_height::Float64
    residual_leakage_time::Float64
end

struct SimulationLayerSnapshot
    timestamp::Float64
    spill_event::SpillEvent
    filled_traps::Vector{Bool}
    injected_volume::Float64
    co2_volume::Float64
end

"""
    SimulationSnapshot

Represents the state of the entire simulation at a point in time.
Contains snapshots from all layers.
"""
struct SimulationSnapshot
    timestamp::Float64
    total_injected_volume::Float64
    total_co2_volume::Float64
    layer_snapshots::Vector{SimulationLayerSnapshot}
end

"""
    LeakageRecord

Records a leakage event for generating WeatherEvents in the overlying layer.
- `start_time`: When leakage started
- `trap_id`: The trap that is leaking
- `leakage_location`: Grid cell (CartesianIndex) where leakage occurs (trap's lowest point)
"""
struct LeakageRecord
    start_time::Float64
    trap_id::Int
    leakage_location::CartesianIndex{2}
end

"""
    LeakageState

Tracks the leakage state for all traps in a layer during simulation.
- `leaking`: Boolean vector indicating if each trap has reached leakage threshold (edge=0)
- `draining`: Boolean vector indicating if each trap is experiencing residual drainage
  (includes leaking traps AND their filled ancestors whose CO2 drains through)
- `leakage_volume`: Volume at which leakage starts for each trap (precomputed from leakage_height)
- `leakage_start_time`: When leakage/drainage started for each trap (Inf if not yet)
- `leakage_records`: Vector of LeakageRecord for generating upstream WeatherEvents
- `leakage_height`: The threshold height for leakage (from ReservoirProperties)
- `initial_volume_at_leak`: Volume in each trap when drainage started (for residual drainage)
- `residual_saturation`: Fraction of CO2 that remains after drainage (from ReservoirProperties)
- `residual_leakage_time`: Time over which residual drainage occurs (from ReservoirProperties)
"""
mutable struct LeakageState
    leaking::Vector{Bool}
    draining::Vector{Bool}
    leakage_volume::Vector{Float64}
    leakage_start_time::Vector{Float64}
    leakage_records::Vector{LeakageRecord}
    leakage_height::Float64
    # Residual leakage fields
    initial_volume_at_leak::Vector{Float64}
    residual_saturation::Float64
    residual_leakage_time::Float64
end
