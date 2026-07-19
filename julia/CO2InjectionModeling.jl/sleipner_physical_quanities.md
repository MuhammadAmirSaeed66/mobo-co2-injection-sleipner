## Summary of Physical Quantities from the Paper

This summary extracts key physical quantities and parameters used in the Invasion Percolation Markov Chain (IPMC) analysis of the Sleipner 2019 Benchmark model. These are based on the provided PDF content, including descriptions, tables, and appendices. Values are intended for replication in your simulation code. Where distributions are mentioned, they refer to stochastic variations for sensitivity analysis (500 realizations per scenario). Units are preserved as in the paper.

### Grid and Model Setup
- **Grid Dimensions**: 64 (I) × 118 (J) × 263 (K) cells, approximately 2 million total cells.
- **Cell Resolution**: Lateral: 50 m × 50 m. Vertical: Proportional layering with ~0.5 m in intra-shale layers, ~2 m in sandstone zones, and ~5 m in caprock.
- **Lateral Extent**: 3.2 km (E-W) × 5.9 km (N-S).
- **Number of Layers**: 9 sandstone layers (L1 deepest to L9 shallowest), separated by 8 intraformational shale layers (7 thin shales ~1.5 m thick; 1 thick shale between L8 and L9 ~7.5 m thick).
- **Injection Point**: In the deepest sandstone layer (L1), at grid center (approximately cell I=32, J=59). Located ~65 m north of the main feeder in the model.
- **Scenarios**:
  - **Continuous Shales (CS)**: Laterally continuous shale layers with no breaks.
  - **Shales with Breaks (SB)**: Includes vertical conduits (feeders) in shales:
    - Main feeder (chimney): ~100 m × 100 m lateral extent, ~65 m south of injection point, intersects all shale layers (properties set to sandstone).
    - North-Eastern feeder: One cell (50 m × 50 m) from L5 to L6, randomly selected within a low-confidence polygon (properties set to sandstone).
    - South-Western feeder: One cell (50 m × 50 m) from L7 to L8, randomly selected within a low-confidence polygon (properties set to sandstone).
- **Simulation Period**: 1996 to 2010 (14 years), with total injected CO₂ of 12.18 Mt by 2010.

### Injection Schedule
Annual CO₂ injection rates (Mt/year) from 1996 to 2010 (cumulative total: 12.18 Mt). Use these for time-series simulation:

| Year | Rate (Mt/year) | Cumulative (Mt) |
|------|----------------|-----------------|
| 1996 | 0.07           | 0.07            |
| 1997 | 0.67           | 0.74            |
| 1998 | 0.85           | 1.59            |
| 1999 | 0.94           | 2.52            |
| 2000 | 0.94           | 3.46            |
| 2001 | 1.02           | 4.48            |
| 2002 | 0.96           | 5.45            |
| 2003 | 0.92           | 6.37            |
| 2004 | 0.76           | 7.12            |
| 2005 | 0.87           | 7.99            |
| 2006 | 0.83           | 8.82            |
| 2007 | 0.93           | 9.74            |
| 2008 | 0.82           | 10.57           |
| 2009 | 0.86           | 11.42           |
| 2010 | 0.76           | 12.18           |

### Rock Properties
- **Porosity (φ)**: Modeled using Sequential Gaussian simulation with normal distributions (specific μ and σ not provided; use typical Utsira values: sandstone ~0.35–0.40; shale ~0.15–0.20).
- **Permeability (k)**: Not a direct input in IP simulations (correlated to threshold pressure via transform, e.g., Nhabanga and Ringrose, 2022). For shales: Vertical permeability (k_v) 0.08–1.5 mD (anisotropy 0.01, so horizontal k_h ~8–150 mD). Sandstone: High permeability (1–2.5 D typical for Utsira).
- **Threshold Pressure (P_th)**:
  - **Sandstone Layers**: Low values (not specified; assume <10 kPa for percolation dominance).
  - **Shale Layers**: Normal distribution.
    - Deterministic base case: μ = 98 kPa, σ = 16 kPa (same for all shales).
    - Sensitivity analysis (500 realizations): For each shale, sample μ uniformly from 49–147 kPa (50% below/above base 98 kPa), then assign normal distribution with that μ (σ likely ~16 kPa, inferred from base case and Table 3 variations ~8–14 kPa).
    - Units: kPa (convert to Pa for code: ×1000).
    - Conversion: Inputs in Hg-Air system, converted to CO₂-brine using γ = 0.0625 (likely N/m or mN/m; check Singh et al., 2010 for exact).
- **Residual CO₂ Saturation (S_gr)**: 20% (trapped CO₂ in invaded cells outside accumulations).
- **Irreducible Water Saturation (S_wir)**: 30%.
- **Maximum CO₂ Saturation (S_g max)**: 70% (in backfilled accumulations beneath baffles/seals).

### Fluid Properties
- **CO₂ Density (ρ_CO2)**: Varies vertically with depth (pressure/temperature gradient): ~350 kg/m³ in uppermost L9 to ~570 kg/m³ in lowermost L1. (Model a linear or P/T-based trend; code uses constant 700 kg/m³—needs update.)
- **Brine Density (ρ_brine)**: Not explicitly stated; infer density contrast (Δρ) ~350–400 kg/m³ from buoyancy calculations (e.g., ~120 kPa for 30 m column).
- **Interfacial Tension (γ)**: 0.0625 (units likely N/m; used for Hg-Air to CO₂-brine conversion; typical CO₂-brine values ~30–50 mN/m—verify).
- **Geothermal Gradient**: 35.6 °C/km (temperature 35 °C at 800 m depth).
- **Salinity**: Not specified (affects density and IFT; assume typical North Sea brine ~3–5% NaCl equivalent).

### Trapping Mechanisms
- **Included**: Structural, stratigraphic, and residual trapping.
- **Excluded**: Dissolution (estimated 0.85–1.8%/year at Sleipner) and mineralization.

### Reference Mass Distribution (2010)
Estimated % CO₂ mass per layer (based on seismic polygon areas; total 12.18 Mt):

| Layer | Area (km²) | Estimated % |
|-------|------------|-------------|
| L9    | 2.13       | 17          |
| L8    | 1.32       | 11          |
| L7    | 0.84       | 7           |
| L6    | 0.95       | 8           |
| L5    | 3.02       | 25          |
| L4    | 0.78       | 6           |
| L3    | 0.88       | 7           |
| L2    | 1.25       | 10          |
| L1    | 1.11       | 9           |

### Items Missing or Needing Implementation in Your Code
Your code is a good starting point but needs adjustments to replicate the paper's analysis:
- **Number of Layers**: Code loads 7 layers (Reflector1–7); update to 9 (check if `sleipner_depth_surfaces.npz` has Reflector1–9 and Base_Reflector1–9 keys).
- **Injection Schedule**: Code uses a simplified ramp (total ~20.5 Mt over 25 years); replace with annual rates from the table above for 1996–2010 accuracy.
- **Multiple Realizations**: Code runs a single simulation; implement a loop for 500 CS and 500 SB realizations, varying shale P_th per run (sample μ uniform 49000–147000 Pa, then normal dist with σ=16000 Pa).
- **Shale Breaks (SB Scenario)**: Not implemented; add logic to set P_th = sand values (or very low, e.g., 0 Pa) in specific cells/polygons for the 3 feeders (main chimney, NE, SW). Use random selection within low-confidence areas for NE/SW.
- **Variable CO₂ Density**: Code uses constant 700 kg/m³; implement depth-dependent ρ_CO2 (e.g., linear from 350 kg/m³ at top to 570 kg/m³ at bottom).
- **Threshold Pressure Units/Distributions**: Code uses fixed ranges (shale 120000–160000 Pa); update to normal dist per shale as above. Ensure independent variation per layer.
- **Residual/Irreducible Saturations**: Code varies by rock type (sand s_gr=0.12, s_wir=0.15; shale 0.08/0.30); align to paper's 20% s_gr and 30% s_wir uniformly.
- **Maximum CO₂ Saturation**: Not explicit in code; add if ips supports (70% in accumulations).
- **Simulation Time/Snapshots**: Code runs to 40 years with 2-year snapshots; limit to 14 years (1996–2010) with snapshots matching seismic (e.g., 1999, 2003, 2010).
- **Output Metrics**: Add calculation of mass per layer (Mt and %), accumulation frequency (%), and RMSE vs. reference (Table above). Generate distributions/boxplots for sensitivity.
- **Software-Specific**: Paper uses Permedia® (Hg-Air inputs converted); your ips may need equivalent conversion for IFT (γ=0.0625).
- **Depth Surfaces**: Ensure npz matches paper's topography (e.g., southern high in L1, northern migration).

If the npz lacks 9 layers or feeder polygons, download the full Sleipner 2019 Benchmark from CO2 Datashare. Let me know if you need help implementing these in code.