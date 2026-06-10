---
name: dark-current-fitting
description: >
  PbS CQD photodetector dark current analysis using dual-model fitting:
  Model 1 (Eq.1) initial diode model and Model 2 (Eq.2) TAT model
  (J_main + J_Ohm + J_TAT/J_non). Segmented fitting with journal-ready
  plots and academic report generation. Auto-detects dark sweeps.
  Trigger on: 暗电流拟合, dark current fitting, Jdark-V, TAT,
  trap-assisted tunneling, 光电探测器, photodiode characterization.
---

# /dark-current-fitting — PbS CQD Dark Current Component Fitting

## Overview

Fits Jdark-V data using **two models**:

### Model 1 (Eq.1) — Initial Equivalent Diode Model
```
J_dark = J0*[exp(qV/(A*kT)) - 1] + V/Rsh + k*V^m
```
5 params: J0, A, Rsh, k, m

### Model 2 (Eq.2) — Optimized General Equation (TAT Model)
```
J_dark = J0*[exp(qV/(A*kT)) - 1] + V/Rsh + B*V*exp(-c/(Vbi-V))
```
6 params: J0, A, Rsh, B, c, Vbi

**Voltage convention**: V_fit = -V_raw, J_fit = -J_raw (standard diode convention).

## Segmented Fitting

- **Stage 1** (|V| < V_seg): J_main dominant → fit J0, A
- **Stage 2** (V < -V_seg): Leakage dominant → fit remaining params
- **Stage 3**: Global refinement

## Plot Formatting

- **Font**: Arial (axis labels bold) + LaTeX math
- **Ticks**: inward, no top/right ticks
- **Grid**: off
- **Legend**: no frame, lower left, small font
- **Colors**: color-blind friendly (Wong 2011)
- **Export**: SVG + PDF + PNG (300 dpi) + PNG (600 dpi with --hd)

## Usage

```bash
python scripts/dark_current_fitting.py <control.txt> <sample.txt> [options]
```

| Param | Description | Default |
|-------|-------------|---------|
| -a, --area | Device area (cm2) | None |
| --model | {both, eq1, eq48} | both |
| --points N | Use first N data points | all |
| --vmin | Lower voltage limit | -0.5 |
| --vmax | Upper voltage limit | 0.2 |
| --no-segmented | Global single-fit | False |
| --no-validate | Skip validation | False |
| --no-report | Skip report | False |
| --fig-width-cm | Figure width in cm | 8.5 |
| --hd | Export 600 dpi PNG | False |

## Dependencies

- numpy, scipy, matplotlib, pandas, python-docx
