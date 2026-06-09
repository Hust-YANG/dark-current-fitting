---
name: dark-current-fitting
description: >
  PbS CQD photodetector dark current analysis using dual-model fitting:
  Model 1 (Eq.1) initial diode model and Model 2 (Eq.2) TAT model
  (J_main + J_Ohm + J_TAT/J_non). Segmented fitting with self-consistency
  validation and academic report generation. Auto-detects dark sweeps from
  mixed dark/light data. Trigger on: 暗电流拟合, dark current fitting,
  Jdark-V, TAT, trap-assisted tunneling, 光电探测器, photodiode characterization.
---

# /dark-current-fitting — PbS CQD Dark Current Component Fitting

## Overview

Fits Jdark-V data using **two models**:

### Model 1 (Eq.1) — Initial Equivalent Diode Model
```
J_dark = J0·[exp(qV/(A·kT)) - 1] + V/Rsh + k·V^m
         J_main (主二极管)       J_Ohm (欧姆漏电)  J_non (非欧姆隧穿)
```

### Model 2 (Eq.2) — Optimized General Equation (TAT Model)
```
J_dark = J0·[exp(qV/(A·kT)) - 1] + V/Rsh + B·V·exp(-c/(Vbi-V))
         J_main (主二极管)       J_Ohm (欧姆漏电)  J_TAT (陷阱辅助隧穿)
```

**Voltage convention**: Raw data V>0=reverse, V<0=forward. Script negates V & J internally (V_fit = -V_raw, J_fit = -J_raw) so standard diode convention applies.

## Model 1 Parameters (Eq.1, 5 params)

| Parameter | Symbol | Unit | Description |
|-----------|--------|------|-------------|
| Reverse saturation current density | J0 | A/cm² | Characterizes intrinsic defect recombination level |
| Ideality factor | A | - | Quality factor, ideal approaches 1.0 |
| Shunt resistance | Rsh | Ω·cm² | Film ohmic leakage pathways |
| Non-ohmic tunneling coefficient | k | - | Related to defect density |
| Non-ohmic tunneling exponent | m | - | Related to tunneling mechanism, m>1 = non-ohmic |

## Model 2 Parameters (Eq.2, 6 params)

| Parameter | Symbol | Unit | Description |
|-----------|--------|------|-------------|
| Reverse saturation current density | J0 | A/cm² | Characterizes intrinsic defect recombination level |
| Ideality factor | A | - | Quality factor, larger A = more defect recombination |
| Shunt resistance | Rsh | Ω·cm² | Film ohmic leakage pathways |
| TAT coefficient (defect density) | B | - | Proportional to trap state density Nt |
| TAT coefficient (tunneling barrier) | C | - | Proportional to tunneling barrier height |
| Built-in voltage | Vbi | V | p-n junction / Schottky junction built-in potential |

> **Naming**: Code uses lowercase `c`, reports use uppercase `C` (publication convention).

## Segmented Fitting

- **Stage 1** (|V| < V_seg, default 0.2V): J_main dominant → fit **J0, A**
- **Stage 2** (V < -V_seg): Leakage dominant → fit **remaining params** (J0, A fixed)
- **Stage 3**: Global refinement — all params co-adjusted

## Output Structure

### `--model both` (default):
```
output/
├── model_eq1/
│   ├── control_fitting.svg / .pdf / .png / _600dpi.png
│   ├── sample_fitting.svg / .pdf / .png / _600dpi.png
│   ├── j_main_comparison.svg / .pdf / .png
│   ├── leakage_comparison.svg / .pdf / .png
│   ├── control_component_fit.txt
│   ├── sample_component_fit.txt
│   └── eq1_fitting_params.csv
├── model_eq2/
│   ├── (same structure as above)
│   └── eq2_fitting_params.csv
└── fitting_report.docx (combined dual-model report)
```

## Plot Formatting (Journal-Ready)

- **Font**: Times New Roman (English) + STSong/SimSun (Chinese fallback)
- **Size**: Single column 8.5 cm, double column ~17.5 cm
- **Ticks**: Inward, uniform size
- **Colors**: Color-blind friendly (Wong 2011, Nature Methods)
- **Export**: SVG + PDF (vector) + PNG (300 dpi) + PNG (600 dpi with --hd)

## Usage

```bash
python scripts/dark_current_fitting.py <control.txt> <sample.txt> [options]
```

| Param | Description | Default |
|-------|-------------|---------|
| control | Control PD data file | Required |
| sample | Sample PD data file | Required |
| -a, --area | Device area (cm²) | None |
| --model | {both, eq1, eq48} | both |
| --points N | Use first N data points | all |
| --vmin | Lower voltage limit, V_fit (V) | -0.5 |
| --vmax | Upper voltage limit, V_fit (V) | 0.2 |
| -T | Temperature (K) | 300 |
| -o | Output directory | ./output |
| --no-auto-dark | Disable auto dark-sweep detection | False |
| --no-segmented | Use global single-fit | False |
| --vseg | Segmentation voltage threshold (V) | 0.2 |
| --no-validate | Skip self-consistency validation | False |
| --no-report | Skip report generation | False |
| --fig-width-cm | Figure width in cm | 8.5 |
| --fig-double-col | Use double-column width | False |
| --hd | Export 600 dpi PNG | False |

## Example

```bash
python scripts/dark_current_fitting.py control.txt sample.txt \
  -a 0.0706858 --model both --points 201 --hd --fig-width-cm 8.5 -o ./output
```

## Dependencies

- numpy, scipy, matplotlib, pandas, python-docx
