# Julia simulator

The simulator used by the R scripts is publicly available at:

```text
https://github.com/MuhammadAmirSaeed66/CO2InjectionModeling.jl
```

From the repository root, install it with:

```r
source("R/install_julia_simulator.R")
```

The installer downloads the simulator into:

```text
julia/CO2InjectionModeling.jl/
```

and records the source commit in:

```text
julia/SIMULATOR_COMMIT.txt
```

Expected minimum structure:

```text
CO2InjectionModeling.jl/
├── Project.toml
├── Manifest.toml        # recommended for exact dependency versions
└── src/
```

The R scripts activate this Julia project, instantiate its dependencies, and
load:

```julia
using CO2InjectionModeling
```

Retain the original simulator license and attribution files when publishing the
combined reproducibility repository.
