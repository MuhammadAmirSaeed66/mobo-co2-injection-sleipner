using SurfaceWaterIntegratedModeling
export Layer, add_boundary_wall, analyze_base_surfaces

struct Layer
    name::String
    trap_structure::TrapStructure
    boundary_padding::Int  # Number of cells padded on each side (0 for open BC)
end

"""
Add wall padding around a topography surface for closed boundary conditions.

Parameters:
- surface: 2D array of elevations
- pad_width: Number of cells to pad on each side
- wall_height: Height to add to create wall (large value to ensure CO2 cannot escape)

Returns:
- Padded surface with walls around the edges
"""
function add_boundary_wall(surface::Matrix{<:Real}, pad_width::Int, wall_height::Float64=1000.0)
    if pad_width == 0
        return surface
    end

    nx, ny = size(surface)
    padded_nx = nx + 2 * pad_width
    padded_ny = ny + 2 * pad_width

    # Create padded surface - initialize with original data in center
    padded_surface = zeros(Float64, padded_nx, padded_ny)
    padded_surface[pad_width+1:pad_width+nx, pad_width+1:pad_width+ny] .= surface

    # Create walls by setting boundary cells to very high elevation
    # This makes them impassable for CO2
    max_elevation = maximum(surface)
    wall_elevation = max_elevation + wall_height

    # Top and bottom boundaries
    padded_surface[1:pad_width, :] .= wall_elevation
    padded_surface[end-pad_width+1:end, :] .= wall_elevation

    # Left and right boundaries
    padded_surface[:, 1:pad_width] .= wall_elevation
    padded_surface[:, end-pad_width+1:end] .= wall_elevation

    return padded_surface
end

function analyze_base_surfaces(topography::SleipnerTopography; boundary_condition::Symbol=:open)::Vector{Layer}
    @assert boundary_condition in [:open, :closed] "boundary_condition must be :open or :closed"

    # Determine padding width: 0 for open BC, 1 for closed BC
    pad_width = boundary_condition == :closed ? 1 : 0

    # Initialize empty vector (will grow as we push)
    layers = Vector{Layer}()

    # Iterate over each sand layer to create Layer structs
    for layer in reverse(topography.sand_layers)
        layer_name = layer["name"]
        base_surface = layer["top"]

        # Add boundary walls if closed BC
        padded_surface = add_boundary_wall(base_surface, pad_width)

        # Compute trap structure on padded surface
        # Adjust lengths to account for padding
        original_length_x = topography.nx * topography.dx
        original_length_y = topography.ny * topography.dy
        padded_length_x = original_length_x + 2 * pad_width * topography.dx
        padded_length_y = original_length_y + 2 * pad_width * topography.dy

        trap_structure = spillanalysis(padded_surface, lengths = (padded_length_x, padded_length_y))
        push!(layers, Layer(layer_name, trap_structure, pad_width))
    end

    return layers
end
