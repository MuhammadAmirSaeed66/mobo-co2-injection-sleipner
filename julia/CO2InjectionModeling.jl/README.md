# CO2InjectionModeling.jl

A Julia package for modeling CO2 injection in subsurface geological formations.

## Setup

### 1. Install Julia

Download and install Julia from [julialang.org](https://julialang.org/downloads/). Version 1.6 or higher is recommended.

### 2. Clone the repository

```bash
git clone git@github.com:ellingsvee/CO2InjectionModeling.jl.git
cd CO2InjectionModeling.jl
```

### 3. Install dependencies

Open Julia in the project directory and run:

```julia
using Pkg
Pkg.activate(".")
Pkg.instantiate()
```

Alternatively, run from the terminal:

```bash
julia --project=. -e 'using Pkg; Pkg.instantiate()'
```

### 4. Verify the installation

Run the main script from the terminal:

```bash
julia --project=. scripts/run.jl
```

This will execute the example simulation to verify everything is working correctly.

### TEST THAT FORK WORKS