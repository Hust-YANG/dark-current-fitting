---
name: dark-current-fitting
description: >
  PbS CQD photodetector dark current analysis using implicit single/double-diode
  models with series resistance. Supports dual-diode (J_diff + J_rec + J_Ohm + J_TAT)
  and single-diode (J_diff + J_Ohm + J_TAT) decomposition. Global fitting via
  curve_fit with multiple initial guesses. Generates journal-ready plots
  (SVG/PDF/PNG), component data (.txt), parameters (CSV), and academic Word report (.docx).
  Auto-detects dark sweeps from mixed dark/light data.
  Trigger on: 暗电流拟合, dark current fitting, Jdark-V, diode model,
  量子点, PbS CQD, photodetector, 光电探测器.
---

# /dark-current-fitting — PbS CQD Dark Current Component Fitting

## Overview

Fits Jdark-V data using an **implicit diode model** with series resistance R_S.
V_int = V − J_D·R_S, solved via damped fixed-point iteration.

## Models

### Dual-Diode Model (8 parameters)
```
J_D = J01·[exp(A1·V_int)−1] + J02·[exp(A2·V_int)−1] + V_int/R_SH + k·V·exp(m/V_int)
      J_diff (diffusion)       J_rec (recombination)    J_Ohm (ohmic)  J_TAT (trap-assisted tunneling)

V_int = V − J_D·R_S
```

### Single-Diode Model (6 parameters, J02 ≡ 0)
```
J_D = J01·[exp(A1·V_int)−1] + V_int/R_SH + k·V·exp(m/V_int)
```

A1 = q/(n1·k_B·T), n1 ≈ 1;  A2 = q/(n2·k_B·T), n2 ≈ 2.

**Voltage convention**: V_fit = −V_raw, J_fit = −I_raw/Area. V_fit > 0 = forward bias.

## Parameters

### Dual-Diode
| Parameter | Symbol | Unit | Bounds |
|-----------|--------|------|--------|
| Diffusion saturation current | J01 | A/cm² | [1e-12, 1e-3] |
| Diffusion coefficient | A1 | V⁻¹ | [10, 50] |
| Recombination saturation current | J02 | A/cm² | [1e-12, 1e-3] |
| Recombination coefficient | A2 | V⁻¹ | [5, 30] |
| Series resistance | R_S | Ω·cm² | [0, 1e4] |
| Shunt resistance | R_SH | Ω·cm² | [1e3, 1e8] |
| TAT coefficient | k | - | [1e-10, 1e2] |
| TAT barrier | m | V | [-1.0, -0.001] |

### Single-Diode
Same but without J02, A2 (6 parameters).

## Current Components
| Component | Formula | Physics |
|-----------|---------|---------|
| J_diff | J01·[exp(A1·V_int)−1] | Diffusion (n1≈1) |
| J_rec | J02·[exp(A2·V_int)−1] | G-R recombination (n2≈2) |
| J_Ohm | V_int / R_SH | Ohmic shunt |
| J_TAT | k·V·exp(m/V_int) | Trap-assisted tunneling (m<0) |

## Output Structure
```
output/
├── model_single/
│   ├── <name>_fitting.svg / .pdf / .png / _600dpi.png
│   ├── <name>_component_fit.txt
│   ├── <name>_params.csv
│   └── <name>_report.docx
├── model_dual/
│   └── (same structure)
```

### Component .txt
**Dual**: `V(V)  J_data  J_fit  J_diff  J_rec  J_Ohm  J_TAT`
**Single**: `V(V)  J_data  J_fit  J_diff  J_Ohm  J_TAT`

## Plot Formatting
- Font: Arial bold axis labels
- Ticks: inward, no grid, no legend frame
- Colors (Wong 2011): Data #0072B2, Fit #D55E00, J_diff #009E73, J_rec #F0E442, J_Ohm #CC79A7, J_TAT #56B4E9
- Y-axis: data −1.5 to +0.6 decades
- X-axis: exact data range, no padding
- Width: 8.5 cm, export SVG+PDF+PNG(300dpi)+PNG(600dpi with --hd)

## Usage
```bash
python scripts/dark_current_fitting.py <data.txt> [options]
```

| Param | Description | Default |
|-------|-------------|---------|
| data | Input .txt data file | Required |
| -o | Output directory | ./output |
| --model | {single, dual, both} | both |
| -T | Temperature (K) | 300 |
| -a | Device area (cm²) | None |
| --vmin | Lower voltage limit (V) | -0.6 |
| --vmax | Upper voltage limit (V) | 0.2 |
| --points N | Use first N data points | all |
| --fig-width-cm | Figure width (cm) | 8.5 |
| --hd | Export 600 dpi PNG | False |
| --no-auto-dark | Disable auto dark-sweep detection | False |

## Dependencies
numpy, scipy, matplotlib, pandas, python-docx