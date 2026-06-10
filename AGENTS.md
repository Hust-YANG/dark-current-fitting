# Dark Current Fitting Skill

## Purpose

Analyze dark current J-V characteristics of PbS CQD photodetectors using implicit
single and dual-diode models with series resistance R_S. V_int = V − J_D·R_S requires
solving an implicit equation via damped fixed-point iteration.

## Models

### Dual-Diode (8 parameters)
```
J_D = J01·[exp(A1·V_int)−1] + J02·[exp(A2·V_int)−1] + V_int/R_SH + k·V·exp(m/V_int)
V_int = V − J_D·R_S
```
A1 = q/(n1·k_B·T), n1 ≈ 1; A2 = q/(n2·k_B·T), n2 ≈ 2; m < 0 for TAT.

### Single-Diode (6 parameters, J02 ≡ 0)
```
J_D = J01·[exp(A1·V_int)−1] + V_int/R_SH + k·V·exp(m/V_int)
```

## Core Implementation

- `solve_implicit()` — damped fixed-point iteration (damping=0.4, tol=1e-10, max_iter=200)
- `_compute_components()` — V_int → J_diff, J_rec, J_Ohm, J_TAT
- `dark_current_dual()` / `dark_current_single()` — model functions for curve_fit
- `fit_dark_current()` — multi-guess global fitting with TRF method
- `plot_fitting()` — semilogy plot with component decomposition
- `generate_word_report()` — English academic .docx report

## Fitting Strategy

Global fitting with 6 initial guesses. curve_fit TRF method. Best R² selected.

## Parameter Bounds

### Dual-Diode: [J01, A1, J02, A2, R_S, R_SH, k, m]
```
lower = [1e-12, 10.0, 1e-12,  5.0,  0.0,  1e3, 1e-10, -1.0]
upper = [1e-3,  50.0, 1e-3,   30.0, 1e4,  1e8, 1e2,   -0.001]
```

### Single-Diode: [J01, A1, R_S, R_SH, k, m]
```
lower = [1e-12, 10.0, 0.0,  1e3, 1e-10, -1.0]
upper = [1e-3,  50.0, 1e4,  1e8, 1e2,   -0.001]
```

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
| --hd | Export additional 600 dpi PNG | False |
| --no-auto-dark | Disable auto dark-sweep detection | False |

## Plot Formatting

- Font: Arial bold axis labels
- Ticks: inward, no grid, no top/right ticks
- Legend: no frame, lower left
- Y-axis: data −1.5 decades to +0.6 decades
- X-axis: exact data range, no padding
- Colors: Wong 2011 palette
- Export: SVG + PDF + PNG(300dpi) + PNG(600dpi, --hd)

## Output Files

Per model directory:
- `<name>_fitting.svg/.pdf/.png/_600dpi.png` — Fitting plot
- `<name>_component_fit.txt` — Tab-separated component data
- `<name>_params.csv` — Parameter table
- `<name>_report.docx` — Academic Word report

Component .txt columns:
- Dual: V, J_data, J_fit, J_diff, J_rec, J_Ohm, J_TAT
- Single: V, J_data, J_fit, J_diff, J_Ohm, J_TAT

## Dependencies

numpy, scipy, matplotlib, pandas, python-docx