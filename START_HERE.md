# Start here — four simple steps

This package is already organized for GitHub. The Julia simulator now has a
public repository, so you do not need to copy its folder manually.

## Step 1 — download the Julia simulator automatically

1. Unzip this package.
2. Open `mobo-co2-injection-sleipner.Rproj` in RStudio.
3. Run:

```r
source("R/install_julia_simulator.R")
```

This downloads:

```text
https://github.com/MuhammadAmirSaeed66/CO2InjectionModeling.jl
```

into:

```text
julia/CO2InjectionModeling.jl/
```

The script records the exact downloaded commit in:

```text
julia/SIMULATOR_COMMIT.txt
```

It also removes the nested `.git` folder so the simulator files can be uploaded
normally with the main paper repository.


### macOS Julia path

The repository now detects a juliaup installation automatically at
`~/.juliaup/bin`. If JuliaCall still reports that Julia is not found, restart
R and run:

```r
Sys.setenv(JULIA_HOME = path.expand("~/.juliaup/bin"))
source("R/examples/01_one_well.R")
```

## Step 2 — install the R packages and run quick tests

Run:

```r
source("R/install_packages.R")
source("R/examples/01_one_well.R")
source("R/examples/03_three_wells.R")
```

You do not need to rerun the full ten-seed experiment merely to publish the
repository. The final reported results are already included.

## Step 3 — check the repository

Run:

```r
source("R/check_repository.R")
```

Confirm that it reports the Julia project, final R script, tables, and figures
as present.

## Step 4 — publish with GitHub Desktop

1. Open GitHub Desktop.
2. Choose **Add an Existing Repository from your Local Drive**.
3. Select this extracted folder.
4. Commit all files with:
   `Initial reproducibility release`.
5. Select **Publish repository**.
6. Name it:
   `mobo-co2-injection-sleipner`.
7. Make it **Public**.
8. Copy the public URL into the manuscript and response letter.

After publication, create a GitHub release called `v1.0.0`.
