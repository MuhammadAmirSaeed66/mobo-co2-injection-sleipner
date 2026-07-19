export swim_volume_to_physical_volume, physical_volume_to_swim_volume
export full_volume_to_rock_volume_scaling, unit_volume_to_physical_scaling

function unit_volume_to_physical_scaling(domain::Domain3D)
    return domain.dx * domain.dy
end

function full_volume_to_rock_volume_scaling(
    reservoir_properties::ReservoirProperties,
)::Float64
    return reservoir_properties.sand_porosity * (1.0 - reservoir_properties.sand_irreducible_water_saturation)
end


function swim_volume_to_physical_volume(
    volume::Union{AbstractArray{<:Real}, Real},
    reservoir_properties::ReservoirProperties,
    domain::Domain3D
)
    scaling_factor = full_volume_to_rock_volume_scaling(reservoir_properties) * unit_volume_to_physical_scaling(domain)
    return volume .* scaling_factor
end

function physical_volume_to_swim_volume(
    volume::Union{AbstractArray{<:Real}, Real},
    reservoir_properties::ReservoirProperties,
    domain::Domain3D
)
    scaling_factor = full_volume_to_rock_volume_scaling(reservoir_properties) * unit_volume_to_physical_scaling(domain)
    return volume ./ scaling_factor
end