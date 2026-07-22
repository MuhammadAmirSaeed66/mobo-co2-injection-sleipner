<div align="center">

<a href="https://github.com/MuhammadAmirSaeed66">
  <img src="https://github.com/MuhammadAmirSaeed66.png?size=180" width="135" alt="Muhammad Amir Saeed">
</a>

# 🌍 Multi-Objective Bayesian Optimization for CO₂ Injection Design

### A reproducible R–Julia framework for mixed-variable optimization of Sleipner-inspired CO₂ storage strategies

<p>
  <a href="https://www.r-project.org/">
    <img src="https://img.shields.io/badge/R-4.4.1-276DC3?logo=r&logoColor=white" alt="R 4.4.1">
  </a>
  <a href="https://julialang.org/">
    <img src="https://img.shields.io/badge/Julia-1.12.2-9558B2?logo=julia&logoColor=white" alt="Julia 1.12.2">
  </a>
  <a href="https://github.com/MuhammadAmirSaeed66/CO2InjectionModeling.jl">
    <img src="https://img.shields.io/badge/Simulator-CO2InjectionModeling.jl-2F855A" alt="Julia simulator">
  </a>
  <img src="https://img.shields.io/badge/Reproducibility-Resources%20available-success" alt="Reproducibility resources">
  <img src="https://img.shields.io/badge/Status-Revised%20Manuscript-blue" alt="Revised manuscript">
  <img src="https://img.shields.io/badge/License-To%20be%20confirmed-orange" alt="License pending">
</p>

<p>
  <b>Companion repository for the revised manuscript submitted to<br>
  <i>Geoenergy Science and Engineering</i></b>
</p>

<p>
  <a href="#-overview">Overview</a> •
  <a href="#-key-contributions">Contributions</a> •
  <a href="#-repository-structure">Structure</a> •
  <a href="#-quick-start">Quick start</a> •
  <a href="#-reproducibility">Reproducibility</a> •
  <a href="#-citation">Citation</a>
</p>

</div>

---

## 📄 Manuscript

**Title:**  
*Multi-Objective Bayesian Optimization Framework for CO₂ Injection Strategy Design: A Sleipner-Inspired Study*

**Authors:**  
[Muhammad Amir Saeed](https://orcid.org/0000-0003-1650-8194),  
[Jo Eidsvik](https://orcid.org/0000-0002-9757-9252), and  
[Antonio Candelieri](https://orcid.org/0000-0003-1431-576X)

**Journal:**  
*Geoenergy Science and Engineering*

This repository contains the final R analysis, the Julia simulator interface, archived numerical outputs, statistical tests, sensitivity analyses, reference-front data, and figure-generation material used in the revised manuscript.

---

## 🚀 Overview

Geological CO₂ storage design involves several competing goals. Injection strategies should retain large CO₂ volumes, limit final-time unretained volume, and control full-chain cost. The design space is also mixed: injection rate is continuous, whereas well configuration and staged layer allocation are categorical.

This repository implements a three-objective **multi-objective Bayesian optimization (MOBO)** workflow using a **multinomial-logit (MNL) surrogate** and four acquisition strategies:

- **Thompson Sampling (TS)**
- **Scalarized Upper Confidence Bound (SCAL-UCB)**
- **Expected Hypervolume Improvement (EHVI)**
- **Expected Preference Improvement (EPI)**

The R optimization workflow communicates with the Julia-based `CO2InjectionModeling.jl` simulator through `JuliaCall`.

---

## ✨ Key Contributions

- Three-objective optimization of retained CO₂ volume, final-time unretained volume, and representative full-chain cost.
- Mixed categorical–continuous decision space with one- and three-well configurations.
- Joint objective-state MNL surrogate for Pareto-aware acquisition sampling.
- Paired ten-seed comparison of TS, SCAL-UCB, EHVI, and EPI.
- Exact deterministic three-dimensional hypervolume.
- Independent grid-only IGD+ reference front constructed from 3,402 designs.
- Held-out comparison of MNL with a mixed-kernel Gaussian process and random forest.
- Hyperparameter, economic, numerical-repeatability, and reference-front sensitivity analyses.
- Fully archived tables, figures, seed-level results, and restart material.

---

## 🎯 Optimization Problem

| Component | Definition |
|---|---|
| Decision vector | `X = (s, q, L_mid, L_top)` |
| Well scenario | One or three wells |
| Injection rate | 0.8–1.0 Mt yr⁻¹ per well |
| Layer decisions | Numerical reservoir intervals 1–9 |
| Objective 1 | Maximize time-integrated retained CO₂ volume |
| Objective 2 | Minimize final-time unretained CO₂ volume |
| Objective 3 | Minimize representative 2024-EUR full-chain cost |
| Simulation horizon | 15 years |
| Optimization output | Pareto-nondominated injection strategies |

> **Scope:** The simulator is used as a deterministic Sleipner-inspired numerical test environment. It is not a history-matched or field-predictive Sleipner model.

---

## 🧪 Final Experiment Settings

| Setting | Value |
|---|---:|
| Paired random seeds | 10 |
| Shared initial designs per seed | 20 |
| Sequential evaluations per method | 80 |
| Total evaluations per method and seed | 100 |
| Acquisition methods | TS, SCAL-UCB, EHVI, EPI |
| Candidate designs per iteration | 30 |
| Predictive samples per candidate | 20 |
| MNL bins per objective | 5 |
| SCAL-UCB coefficient | κ = 2 |
| Primary IGD+ reference | Grid-only nondominated front |
| Dense reference grid | 3,402 designs |
| Primary objective evaluations | 4,000 |

---

## 📊 Main Results

| Method | Final HV | Final IGD+ | Pareto size |
|---|---:|---:|---:|
| TS | 0.5917 ± 0.0277 | 0.0651 ± 0.0126 | 27.4 ± 13.1 |
| SCAL-UCB | **0.6162 ± 0.0422** | **0.0503 ± 0.0261** | **35.5 ± 10.8** |
| EHVI | 0.6075 ± 0.0374 | 0.0577 ± 0.0199 | 25.4 ± 10.7 |
| EPI | 0.6068 ± 0.0325 | 0.0623 ± 0.0141 | 14.9 ± 3.8 |

SCAL-UCB achieved the highest numerical mean hypervolume and lowest numerical mean IGD+. However, the paired Friedman tests did not detect statistically significant method effects for final hypervolume or IGD+. Statistically supported pairwise differences were limited to larger Pareto sets for TS and SCAL-UCB than for EPI.

---

## 📂 Repository Structure

```text
mobo-co2-injection-sleipner/
├── README.md
├── START_HERE.md
├── REPRODUCIBILITY.md
├── CITATION.cff
├── mobo-co2-injection-sleipner.Rproj
│
├── R/
│   ├── 01_run_mobo_final.R
│   ├── install_julia_simulator.R
│   ├── install_packages.R
│   ├── project_setup.R
│   ├── check_repository.R
│   │
│   ├── examples/
│   │   ├── 01_one_well.R
│   │   ├── 02_two_wells_demo.R
│   │   ├── 03_three_wells.R
│   │   └── r_julia_interface_demo.R
│   │
│   └── figures/
│       └── plot_sleipner_layers.R
│
├── julia/
│   └── CO2InjectionModeling.jl/
│
├── data/
│   └── README.md
│
└── results/
    ├── reported/
    │   ├── tables/
    │   ├── figures/
    │   ├── reference/
    │   ├── sensitivity/
    │   ├── surrogate_validation/
    │   ├── simulator_validation/
    │   └── cost_sensitivity/
    │
    └── archive/
        └── final_results_complete.zip
```

The paper experiment compares **one-well and three-well scenarios**. The two-well file is included only as an R–Julia interface demonstration.

---

## ⚙️ Requirements

The repository has been tested with:

- **R 4.4.1**
- **Julia 1.12.2**
- **Git**
- R packages listed in `R/install_packages.R`
- Julia dependencies defined by the simulator `Project.toml`

On macOS with `juliaup`, Julia is commonly available at:

```text
~/.juliaup/bin/julia
```

---

## ▶️ Quick Start

### 1. Clone the repository

```bash
git clone https://github.com/MuhammadAmirSaeed66/mobo-co2-injection-sleipner.git
cd mobo-co2-injection-sleipner
```

### 2. Install the Julia simulator

From R:

```r
source("R/install_julia_simulator.R")
```

This downloads the simulator from:

```text
https://github.com/MuhammadAmirSaeed66/CO2InjectionModeling.jl
```

and records the exact downloaded commit in:

```text
julia/SIMULATOR_COMMIT.txt
```

### 3. Install the R packages

```r
source("R/install_packages.R")
```

### 4. Run the one-well smoke test

```r
Sys.setenv(JULIA_HOME = path.expand("~/.juliaup/bin"))
source("R/examples/01_one_well.R")
```

### 5. Run the three-well example

```r
source("R/examples/03_three_wells.R")
```

### 6. Run the complete experiment

```r
source("R/01_run_mobo_final.R")
```

> The full ten-seed experiment is computationally intensive. Checkpoint and restart files are supported. The final reported outputs are already included under `results/reported/`.

---

## 🔁 Reproducibility

The repository provides:

- identical 20-design initialization within each paired seed;
- ten paired random seeds;
- method- and iteration-specific random-number streams;
- exact hypervolume calculations;
- an independent 3,402-design grid-only reference front;
- complete seed-level run archives;
- implementation and hyperparameter settings;
- statistical-analysis outputs;
- surrogate-validation data;
- sensitivity-analysis outputs;
- final manuscript tables and figures.

Run the repository check with:

```r
source("R/check_repository.R")
```

Detailed instructions are available in:

- [`START_HERE.md`](START_HERE.md)
- [`REPRODUCIBILITY.md`](REPRODUCIBILITY.md)
- [`data/README.md`](data/README.md)

---

## 🧩 Julia Simulator

The Julia simulator is maintained separately at:

[![Julia simulator](https://img.shields.io/badge/GitHub-CO2InjectionModeling.jl-181717?logo=github)](https://github.com/MuhammadAmirSaeed66/CO2InjectionModeling.jl)

The R workflow activates the simulator project, instantiates its dependencies, and loads:

```julia
using CO2InjectionModeling
```

The exact simulator source commit used by a local installation is recorded in `julia/SIMULATOR_COMMIT.txt`.

---

## 🗂️ Data

The source Sleipner plume-boundary data are not redistributed unless their source license permits redistribution.

See [`data/README.md`](data/README.md) for:

- source information;
- expected directory structure;
- required filenames;
- coordinate-reference information;
- plotting instructions.

---

## 📚 Citation

A machine-readable citation is provided in [`CITATION.cff`](CITATION.cff).

Until the article DOI is available, cite the repository as:

```bibtex
@software{saeed2026mobo_co2,
  author    = {Saeed, Muhammad Amir and Eidsvik, Jo and Candelieri, Antonio},
  title     = {Multi-Objective Bayesian Optimization for CO2 Injection Design},
  year      = {2026},
  version   = {1.0.0},
  url       = {https://github.com/MuhammadAmirSaeed66/mobo-co2-injection-sleipner}
}
```

---

## 👨‍💻 Lead Author

<div align="center">

<a href="https://github.com/MuhammadAmirSaeed66">
  <img src="https://github.com/MuhammadAmirSaeed66.png?size=180" width="110" alt="Muhammad Amir Saeed">
</a>

### Muhammad Amir Saeed

[![GitHub](https://img.shields.io/badge/GitHub-MuhammadAmirSaeed66-181717?logo=github)](https://github.com/MuhammadAmirSaeed66)
[![ORCID](https://img.shields.io/badge/ORCID-0000--0003--1650--8194-A6CE39?logo=orcid&logoColor=white)](https://orcid.org/0000-0003-1650-8194)

</div>

---

## 📜 License

The software license should be confirmed by all authors before the final public release. Once selected, replace the current license badge and this section with the approved license, for example MIT or BSD-3-Clause.

Third-party simulator components and external Sleipner data remain subject to their original licenses and attribution requirements.

---

## ⭐ Support

If this repository supports your research, please consider:

- citing the associated manuscript and software;
- starring the repository;
- reporting reproducibility issues through GitHub Issues.

<div align="center">


</div>
