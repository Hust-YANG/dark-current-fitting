# Dark Current Fitting Skill

## Purpose

Analyze dark current J-V of PbS CQD photodetectors using implicit diode models with R_S.
V_int = V − J_D·R_S, solved via damped fixed-point iteration.

## Models

### Dual-Diode (8 params)
```
J_D = J01·[exp(A1·V_int)−1] + J02·[exp(A2·V_int)−1] + V_int/R_SH + k·V·exp(m/V_int)
      J_diff (diffusion)       J_rec (recombination)     J_Ohm      J_TAT (TAT)
```

### Single-Diode (6 params, J02=0)
```
J_D = J01·[exp(A1·V_int)−1] + V_int/R_SH + k·V·exp(m/V_int)
      J_main (diffusion+recombination)  J_Ohm      J_TAT (TAT)
```

A1=q/(n1·kT), n1≈1; A2=q/(n2·kT), n2≈2; m<0 for TAT.

## Core Functions
- `solve_implicit()` — damped fixed-point (damping=0.4, tol=1e-10, max_iter=200)
- `_compute_components()` — V_int → J_main, J_rec, J_Ohm, J_TAT
- `plot_fitting()` — dual labels J_diff, single labels J_main
- `fit_dark_current()` — 6 initial guesses, curve_fit TRF, best R² selected

## Parameter Bounds

### Dual: [J01, A1, J02, A2, R_S, R_SH, k, m]
```
lower = [1e-12, 10.0, 1e-12,  5.0,  0.0,  1e3, 1e-10, -1.0]
upper = [1e-3,  50.0, 1e-3,   30.0, 1e4,  1e8, 1e2,   -0.001]
```

### Single: [J01, A1, R_S, R_SH, k, m]
```
lower = [1e-12, 10.0, 0.0,  1e3, 1e-10, -1.0]
upper = [1e-3,  50.0, 1e4,  1e8, 1e2,   -0.001]
```

## Usage
```bash
python scripts/dark_current_fitting.py <data.txt> [options]
```
| Param | Default |
|-------|---------|
| --model | both |
| --vmin | -0.6 |
| --vmax | 0.2 |
| -T | 300 |

## Plot Formatting
- Arial bold; ticks inward; no grid; legend no frame
- Y: data −1.5/+0.6 decades; X: exact range, no padding
- Colors: Wong 2011; J_diff/J_main share #009E73
- Export: SVG+PDF+PNG(300dpi)+PNG(600dpi --hd)

## Output Files
- `<name>_fitting.svg/.pdf/.png/_600dpi.png`
- `<name>_component_fit.txt` — Dual: V, J_data, J_fit, J_diff, J_rec, J_Ohm, J_TAT; Single: V, J_data, J_fit, J_main, J_Ohm, J_TAT
- `<name>_params.csv`, `<name>_report.docx`

## Dependencies
numpy, scipy, matplotlib, pandas, python-docx