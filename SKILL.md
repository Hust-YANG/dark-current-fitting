---
name: dark-current-fitting
description: >
  PbS CQD photodetector dark current analysis using implicit single/double-diode
  models with series resistance. Supports dual-diode (J_diff + J_rec + J_Ohm + J_tun)
  and single-diode (J_diff + J_Ohm + J_tun) decomposition. Global fitting via
  curve_fit with multiple initial guesses. Generates journal-ready plots
  (SVG/PDF/PNG), component data (.txt), parameters (CSV), and academic Word report (.docx).
  Auto-detects dark sweeps from mixed dark/light data.
  Trigger on: 暗电流拟合, dark current fitting, Jdark-V, diode model,
  量子点, PbS CQD, photodetector, 光电探测器.
---

# /dark-current-fitting — PbS CQD Dark Current Component Fitting

## Overview

Fits Jdark-V data using an **implicit diode model** with series resistance R_S.
The intrinsic junction voltage is V_int = V − J_D·R_S, making the equation implicit
(solved via damped fixed-point iteration).

## Models

### Dual-Diode Model (8 parameters)
```
J_D = J01·[exp(A1·V_int)−1] + J02·[exp(A2·V_int)−1] + V_int/R_SH + k·V_int^m
      J_diff (diffusion)       J_rec (recombination)    J_Ohm (ohmic)  J_tun (tunneling)

V_int = V − J_D·R_S
```

- A1 = q/(n1·k_B·T), n1 ≈ 1 (diffusion ideality factor)
- A2 = q/(n2·k_B·T), n2 ≈ 2 (G-R recombination ideality factor)

### Single-Diode Model (6 parameters, J02 ≡ 0)
```
J_D = J01·[exp(A1·V_int)−1] + V_int/R_SH + k·V_int^m
      J_diff (diffusion)       J_Ohm (ohmic)  J_tun (tunneling)
```

**Voltage convention**: Raw data V>0 = reverse, V<0 = forward. Script negates V & J internally
(V_fit = −V_raw, J_fit = −J_raw) so V_fit > 0 = forward bias in standard diode convention.

## Parameters

### Dual-Diode
| Parameter | Symbol | Unit | Physical Meaning | Bounds |
|-----------|--------|------|------------------|--------|
| Diffusion saturation current | J01 | A/cm² | Minority carrier diffusion | [1e-12, 1e-3] |
| Diffusion exponential coefficient | A1 | V⁻¹ | A1 = q/(n1·k·T), n1≈1 | [10, 50] |
| Recombination saturation current | J02 | A/cm² | G-R center density | [1e-12, 1e-3] |
| Recombination exponential coefficient | A2 | V⁻¹ | A2 = q/(n2·k·T), n2≈2 | [5, 30] |
| Series resistance | R_S | Ω·cm² | Contact + transport layer drop | [0, 1e4] |
| Shunt resistance | R_SH | Ω·cm² | Ohmic leakage paths | [1e3, 1e8] |
| Tunneling coefficient | k | - | ∝ trap density N_t | [1e-10, 1e2] |
| Tunneling exponent | m | - | Tunneling mechanism, typically 2–6 | [1, 6] |

### Single-Diode
Same as above but without J02 and A2 (6 parameters total).

## Current Components
| Component | Formula | Physics |
|-----------|---------|---------|
| J_diff | J01·[exp(A1·V_int)−1] | Diffusion current (n1≈1) |
| J_rec | J02·[exp(A2·V_int)−1] | G-R recombination (n2≈2), dual only |
| J_Ohm | V_int / R_SH | Ohmic shunt leakage |
| J_tun | k · V_int^m | Power-law field-assisted tunneling |

## Implicit Equation Solver

V_int = V − J_D·R_S makes the equation implicit. Solved via damped fixed-point
(Picard) iteration with convergence tolerance 1e-10 and max 200 iterations.

## Output Structure

### `--model both` (default: dual only)
```
output/
├── model_single/
│   ├── <name>_fitting.svg / .pdf / .png / _600dpi.png
│   ├── <name>_component_fit.txt
│   ├── <name>_params.csv
│   └── <name>_report.docx
├── model_dual/
│   ├── <name>_fitting.svg / .pdf / .png / _600dpi.png
│   ├── <name>_component_fit.txt
│   ├── <name>_params.csv
│   └── <name>_report.docx
```

### Component .txt format

**Dual-Diode**: `V(V)  J_data(A/cm2)  J_fit(A/cm2)  J_diff(A/cm2)  J_rec(A/cm2)  J_Ohm(A/cm2)  J_tun(A/cm2)`

**Single-Diode**: `V(V)  J_data(A/cm2)  J_fit(A/cm2)  J_diff(A/cm2)  J_Ohm(A/cm2)  J_tun(A/cm2)`

## Plot Formatting (Journal-Ready)

- **Font**: Arial bold axis labels, Arial + CJK fallback
- **Size**: Single column 8.5 cm (configurable via --fig-width-cm)
- **Ticks**: Inward, uniform size (major=4pt, width=0.8pt)
- **Legend**: No frame, transparent, non-obscuring placement
- **Grid**: None
- **Colors** (Wong 2011, color-blind friendly):
  - Data: #0072B2 (blue), Total fit: #D55E00 (orange)
  - J_diff: #009E73 (green), J_rec: #F0E442 (gold)
  - J_Ohm: #CC79A7 (purple), J_tun: #56B4E9 (sky blue)
- **Y-axis**: semilogy, data range ±0.75 decades
- **Export**: SVG + PDF (vector) + PNG (300 dpi) + PNG (600 dpi with --hd)

## Usage

```bash
python scripts/dark_current_fitting.py <data.txt> [options]
```

| Param | Description | Default |
|-------|-------------|---------|
| data | Input .txt data file | Required |
| -o, --output | Output directory | ./output |
| --model | {single, dual, both} | dual |
| -T, --temperature | Temperature (K) | 300 |
| -a, --area | Device area (cm²) | None |
| --vmin | Lower voltage limit, V_fit (V) | -0.5 |
| --vmax | Upper voltage limit, V_fit (V) | 0.2 |
| --points N | Use first N data points | all |
| --fig-width-cm | Figure width in cm | 8.5 |
| --hd | Export 600 dpi PNG | False |
| --no-auto-dark | Disable auto dark-sweep detection | False |

## Examples

```bash
# Dual-diode fitting on sample
python scripts/dark_current_fitting.py sample.txt \
  -a 0.0706858 --model dual --points 201 --hd -o ./output

# Both models
python scripts/dark_current_fitting.py sample.txt \
  -a 0.0706858 --model both --points 201 --hd -o ./output

# Single-diode only
python scripts/dark_current_fitting.py sample.txt \
  -a 0.0706858 --model single --points 201 -o ./output
```

## Dependencies

- numpy, scipy, matplotlib, pandas, python-docx