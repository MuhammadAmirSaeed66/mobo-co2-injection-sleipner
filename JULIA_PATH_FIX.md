# Fix for “Julia is not found”

The simulator download succeeded. The remaining error means that the R package
`JuliaCall` could not locate the Julia executable on the computer.

The updated `R/project_setup.R` automatically checks:

1. `JULIA_HOME`;
2. `~/.juliaup/bin`;
3. the system `PATH`;
4. `/opt/homebrew/bin`;
5. `/usr/local/bin`;
6. standard Julia application folders on macOS.

## Run the test again

Restart R first:

```r
.rs.restartR()
```

Then run:

```r
source("R/examples/01_one_well.R")
```

For a manual one-session fix on macOS with juliaup:

```r
Sys.setenv(JULIA_HOME = path.expand("~/.juliaup/bin"))
source("R/examples/01_one_well.R")
```

Check the executable directly with:

```r
file.exists(path.expand("~/.juliaup/bin/julia"))
system2(
  path.expand("~/.juliaup/bin/julia"),
  "--version"
)
```

The warning about a missing `Manifest.toml` is not a failure.
`Project.toml` is sufficient for `Pkg.instantiate()`. A generated
`Manifest.toml` can be committed later to lock exact Julia dependencies.
