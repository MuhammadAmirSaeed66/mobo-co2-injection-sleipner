# Multi-Objective Bayesian Optimization for CO2 Injection Design

This repository contains the final R analysis, R–Julia simulator examples,
reported outputs, statistical tables, sensitivity analyses, and figure files
for the manuscript:

**Multi-Objective Bayesian Optimization Framework for CO2 Injection Strategy
Design: A Sleipner-Inspired Study**

Authors: Muhammad Amir Saeed, Jo Eidsvik, and Antonio Candelieri.

## Start here

Read [`START_HERE.md`](START_HERE.md). It contains the shortest setup and
GitHub-upload instructions.

## Main files

- `R/01_run_mobo_final.R` — complete paper experiment and analyses.
- `R/examples/01_one_well.R` — one-well simulator smoke test.
- `R/examples/02_two_wells_demo.R` — two-well interface demonstration only.
- `R/examples/03_three_wells.R` — three-well simulator demonstration.
- `R/examples/r_julia_interface_demo.R` — additional interface examples.
- `R/figures/plot_sleipner_layers.R` — portable static plotting script.
- `results/reported/` — extracted tables, figures, reference fronts and archived
  analysis objects used for the revised manuscript.
- `results/archive/final_results_complete.zip` — complete submitted output
  archive, including restart/cache material.

The paper optimization compares **one-well and three-well scenarios**. The
included two-well file is an interface demonstration and is not part of the
reported 10-seed comparison.

## Final experiment settings

- Paired random seeds: 1–10
- Shared initial designs per seed: 20
- Sequential evaluations per method: 80
- Total evaluations per method and seed: 100
- Candidate designs per iteration: 30
- MNL objective bins: 5
- Predictive samples per candidate: 20
- SCAL-UCB coefficient: 2
- Acquisition methods: TS, SCAL-UCB, EHVI and EPI
- Primary IGD+ reference: 3,402-design grid-only nondominated front

## Julia simulator

The exact Julia simulator is available publicly at:

```text
https://github.com/MuhammadAmirSaeed66/CO2InjectionModeling.jl
```

Install it automatically from the repository root:

```r
source("R/install_julia_simulator.R")
```

The installer downloads the source into `julia/CO2InjectionModeling.jl/` and
records the exact source commit in `julia/SIMULATOR_COMMIT.txt`. The R workflow
then activates this project, instantiates its Julia dependencies, and loads
`CO2InjectionModeling`.

## Quick use

Open `mobo-co2-injection-sleipner.Rproj`, then run:

```r
source("R/install_packages.R")
source("R/examples/01_one_well.R")
```

After the smoke test succeeds, the complete experiment is run with:

```r
source("R/01_run_mobo_final.R")
```

The complete run is computationally intensive and uses checkpoint files to
support restart.

## Results

The final manuscript tables and figures are already provided under
`results/reported/`. A complete original results archive is provided under
`results/archive/`.

## Sleipner plotting data

The source plume-boundary data are not redistributed in this starter package.
See `data/README.md` for the expected folder and file structure.

## Citation

See `CITATION.cff`.

## Licensing

The authors should confirm the preferred open-source license before publication.
The external Julia simulator and any Sleipner source data remain governed by
their own licenses and attribution requirements.
