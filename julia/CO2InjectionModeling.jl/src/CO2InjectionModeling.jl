module CO2InjectionModeling

# Load subfiles
include("sleipner_setup.jl")
include("structs.jl")
include("layer_analysis.jl")
include("volume_conversion.jl")
include("utils.jl")  # Must come before leakage.jl (provides get_all_parents)
include("leakage.jl")


include("fill_layer.jl")
include("fill_layers.jl")
include("analysis.jl")
include("visualization.jl")



include("CO2RInterface.jl")

using .CO2RInterface
export setup_simulator, configure_reservoir, setup_sleipner_reservoir, run_simulation
export generate_cross_section_animation, generate_birdseye_animation

end # module