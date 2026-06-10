# Dark Current Fitting Skill

## Purpose

Analyze dark current (Jdark) vs voltage (V) characteristics of PbS CQD photodetectors
using implicit single and dual-diode models with series resistance R_S.
The intrinsic junction voltage V_int = V − J_D·R_S requires solving an implicit
equation via damped fixed-point iteration.
Supports both dual-diode (J_diff + J_rec + J_Ohm + J_tun) and
single-diode (J_diff + J_Ohm + J_tun) decomposition.
Global fitting via scipy curve_fit with multiple initial guesses.
Auto-detects dark sweeps from mixed dark/light data.

## Models

### Dual-Diode (8 parameters)
```
J_D = J01·[exp(A1·V_int)−1] + J02·[exp(A2·V_int)−1] + V_int/R_SH + k·V_int^m
      J_diff (diffusion)       J_rec (G-R recombination)  J_Ohm        J_tun

V_int = V − J_D·R_S
```

A1 = q/(n1·k_B·T), n1 ≈ 1 (diffusion)
A2 = q/(n2·k_B·T), n2 ≈ 2 (G-R recombination)

### Single-Diode (6 parameters, J02 ≡ 0)
```
J_D = J01·[exp(A1·V_int)−1] + V_int/R_SH + k·V_int^m
```

**Voltage convention**: V_fit = −V_raw, J_fit = −I_raw/Area.
V_fit > 0 = forward bias in all outputs.

## Core Implementation

- `solve_implicit()` — damped fixed-point iteration (damping=0.4, tol=1e-10, max_iter=200)
- `_compute_components()` — decompose V_int into J_diff, J_rec, J_Ohm, J_tun
- `dark_current_dual()` / `dark_current_single()` — model functions for curve_fit
- `fit_dark_current()` — multi-guess global fitting with TRF method
- `plot_fitting()` — semilogy plot with component decomposition
- `generate_word_report()` — English academic .docx report

## Fitting Strategy

Global fitting with 6 initial guesses covering different parameter regions.
`curve_fit` with TRF (Trust Region Reflective) method handles multi-scale parameters.
Best R² result selected.

## Parameter Bounds

### Dual-Diode: [J01, A1, J02, A2, R_S, R_SH, k, m]
```
lower = [1e-12, 10.0, 1e-12,  5.0,  0.0,  1e3, 1e-10, 1.0]
upper = [1e-3,  50.0, 1e-3,   30.0, 1e4,  1e8, 1e2,   6.0]
```

### Single-Diode: [J01, A1, R_S, R_SH, k, m]
```
lower = [1e-12, 10.0, 0.0,  1e3, 1e-10, 1.0]
upper = [1e-3,  50.0, 1e4,  1e8, 1e2,   6.0]
```

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
| --hd | Export additional 600 dpi PNG | False |
| --no-auto-dark | Disable auto dark-sweep detection | False |

## Plot Formatting (Journal-Ready)

- Font: Arial bold axis labels, Arial + CJK fallback
- Ticks: inward, major size=4pt, width=0.8pt
- No grid, no top/right ticks
- Legend: no frame, lower left
- Y-axis: data range ±0.75 decades
- Colors: Wong 2011 palette (color-blind friendly)
- Export: SVG + PDF + PNG(300dpi) + PNG(600dpi, --hd)
- Width: 8.5 cm (single column)

## Output Files

Per model directory:
- `<name>_fitting.svg` / `.pdf` / `.png` / `_600dpi.png` — Fitting plot
- `<name>_component_fit.txt` — Tab-separated component data
- `<name>_params.csv` — Parameter table
- `<name>_report.docx` — Academic Word report

Component .txt columns:
- Dual: V, J_data, J_fit, J_diff, J_rec, J_Ohm, J_tun
- Single: V, J_data, J_fit, J_diff, J_Ohm, J_tun

## Auto Dark-Sweep Detection

For multi-sweep data (Index column with "sweep/point" format):
- Computes |I@V≈0| / max|I| ratio per sweep
- Ratio < 1% → dark sweep (keep)
- Ratio ≥ 1% → light sweep (discard)

## Dependencies

numpy, scipy, matplotlib, pandas, python-docx