# Dark Current Fitting Skill

## Purpose

Analyze dark current (Jdark) vs voltage (V) characteristics of PbS CQD photodetectors
using dual-model fitting: Model 1 (Eq.1) and Model 2 (Eq.2 TAT model).
Performs segmented fitting and generates a combined academic report (.docx).
Auto-detects dark sweeps from mixed dark/light data.

## Models

### Model 1 (Eq.1) — Initial Equivalent Diode Model
```
J_dark = J0*[exp(qV/(A*kT)) - 1] + V/Rsh + k*V^m
```
5 parameters: J0, A, Rsh, k, m

### Model 2 (Eq.2) — Optimized General Equation (TAT Model)
```
J_dark = J0*[exp(qV/(A*kT)) - 1] + V/Rsh + B*V*exp(-c/(Vbi-V))
```
6 parameters: J0, A, Rsh, B, c, Vbi

**Voltage convention**: Raw data V>0=reverse, V<0=forward. Script negates V & J internally (V_fit = -V_raw, J_fit = -J_raw).

## Segmented Fitting

- **Stage 1** (|V| < V_seg, default 0.2V): J_main dominant → fit **J0, A**
- **Stage 2** (V < -V_seg): Leakage dominant → fit **remaining params** (J0, A fixed)
- **Stage 3**: Global refinement — all params co-adjusted

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

## Outputs

### Dual model (--model both):
- Per-model directory with fitting plots (SVG/PDF/PNG), component data, params CSV
- Combined .docx report (5 sections)

### Single model:
- Fitting plots (SVG/PDF/PNG) in output dir
- 2 component data .txt files
- Parameters table CSV
- Academic report (.docx)

## Plot Formatting

- Font: Times New Roman + STSong/SimSun for CJK
- Ticks: inward, uniform
- Colors: color-blind friendly (Wong 2011 palette)
- Export: SVG + PDF + PNG (300 dpi) + PNG (600 dpi with --hd)
- Width: 8.5 cm (single column) or ~17.5 cm (double column)

See SKILL.md for full documentation.
