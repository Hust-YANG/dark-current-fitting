# Dark Current Fitting — PbS CQD Photodetector

A Python tool for quantitative dark current component fitting of PbS colloidal quantum dot (CQD) infrared photodetectors.
Decomposes measured J-V characteristics into physically distinct transport components using global curve fitting
with an **implicit diode model** incorporating series resistance R<sub>S</sub>.

## Models

### Dual-Diode Model (8 parameters)

<div align="center">

**J<sub>D</sub> = J<sub>01</sub>·[exp(A<sub>1</sub>·V<sub>int</sub>) − 1] + J<sub>02</sub>·[exp(A<sub>2</sub>·V<sub>int</sub>) − 1] + V<sub>int</sub>/R<sub>SH</sub> + k·V·exp(m/V<sub>int</sub>)**

**V<sub>int</sub> = V − J<sub>D</sub>·R<sub>S</sub>** *(implicit equation, solved by damped fixed-point iteration)*

*J<sub>diff</sub> (Diffusion, n<sub>1</sub>≈1) + J<sub>rec</sub> (G-R Recombination, n<sub>2</sub>≈2) + J<sub>Ohm</sub> (Ohmic) + J<sub>TAT</sub> (Trap-Assisted Tunneling)*

</div>

**Attention:A<sub>1</sub> = q<sub>/(n<sub>1</sub>·k<sub>B</sub>·T<sub>); A<sub>2</sub> = q<sub>/(n<sub>2</sub>·k<sub>B</sub>·T<sub>)**

</div>

### Single-Diode Model (6 parameters, J<sub>02</sub> ≡ 0)

<div align="center">

**J<sub>D</sub> = J<sub>01</sub>·[exp(A<sub>1</sub>·V<sub>int</sub>) − 1] + V<sub>int</sub>/R<sub>SH</sub> + k·V·exp(m/V<sub>int</sub>)**

*J<sub>main</sub> (Diffusion + Recombination combined) + J<sub>Ohm</sub> + J<sub>TAT</sub>*

</div>

> **Voltage convention**: Raw V>0 = reverse, V<0 = forward. Script negates internally so V<sub>fit</sub>>0 = forward bias.
>
> **Implicit solver**: V<sub>int</sub> = V − J<sub>D</sub>·R<sub>S</sub> makes the equation implicit. Solved via damped fixed-point iteration (damping=0.4, tol=1×10⁻¹⁰).

## Current Components

| Component | Dual-Diode | Single-Diode | Physics |
|-----------|-----------|--------------|---------|
| Primary | **J<sub>diff</sub>** | **J<sub>main</sub>** | Diffusion (n₁≈1) / Diffusion+Recombination |
| Recombination | **J<sub>rec</sub>** | — | G-R center recombination (n₂≈2) |
| Ohmic | **J<sub>Ohm</sub>** | J<sub>Ohm</sub> | Shunt leakage through pinholes |
| Tunneling | **J<sub>TAT</sub>** | J<sub>TAT</sub> | Trap-assisted tunneling, k·V·exp(m/V<sub>int</sub>) |

## Parameters

| Parameter | Symbol | Unit | Bounds | Description |
|-----------|--------|------|--------|-------------|
| Diffusion saturation current | J<sub>01</sub> | A/cm² | [1×10⁻¹², 1×10⁻³] | Minority carrier diffusion |
| Diffusion coefficient | A<sub>1</sub> | V⁻¹ | [10, 50] | A₁ = q/(n₁·k<sub>B</sub>·T), n₁≈1 |
| Recombination saturation current | J<sub>02</sub> | A/cm² | [1×10⁻¹², 1×10⁻³] | G-R center density (dual only) |
| Recombination coefficient | A<sub>2</sub> | V⁻¹ | [5, 30] | A₂ = q/(n₂·k<sub>B</sub>·T), n₂≈2 |
| Series resistance | R<sub>S</sub> | Ω·cm² | [0, 1×10⁴] | Contact + transport layer voltage drop |
| Shunt resistance | R<sub>SH</sub> | Ω·cm² | [1×10³, 1×10⁸] | Ohmic leakage pathways |
| TAT coefficient | k | — | [1×10⁻¹⁰, 1×10²] | ∝ trap density N<sub>t</sub> |
| TAT barrier | m | V | [−1.0, −0.001] | m<0 → tunneling grows with reverse bias |

## Fitting Strategy

Global fitting via `scipy.optimize.curve_fit` (TRF method) with **6 diverse initial guesses** covering different parameter regions.
Best R² selected. No segmented fitting — fully global optimization with physically constrained parameter bounds.

## Sample Output

<p align="center">
<img width="600" alt="sample_fitting" src="https://github.com/user-attachments/assets/f7cd9e9a-ad39-47c6-8df9-a6f7d8a57dd7" />
</p>

## Installation

```bash
pip install numpy scipy matplotlib pandas python-docx
```

## Usage

```bash
# Run with both single & dual diode models (default)
python scripts/dark_current_fitting.py sample.txt -a 0.0707 --points 201 --hd

# Dual-diode only
python scripts/dark_current_fitting.py sample.txt -a 0.0707 --model dual

# Single-diode only
python scripts/dark_current_fitting.py sample.txt -a 0.0707 --model single
```

| Argument | Description | Default |
|----------|-------------|---------|
| `data` | Input data file (.txt) | **Required** |
| `-a, --area` | Device area (cm²) | None |
| `--model` | `single`, `dual`, `both` | `both` |
| `--points N` | Use first N data points | All |
| `--vmin` | Lower voltage limit (V) | −0.6 |
| `--vmax` | Upper voltage limit (V) | 0.2 |
| `-T` | Temperature (K) | 300 |
| `--fig-width-cm` | Figure width (cm) | 8.5 |
| `--hd` | Export 600 dpi PNG | False |
| `--no-auto-dark` | Disable auto dark-sweep detection | False |
| `-o` | Output directory | ./output |

### Input Format

Tab-separated text file with 3-column format: `Index  V1 (V)  I1 (A)`

Supports multi-sweep data (`1/1, 1/2, ...`) with automatic dark-sweep detection (|I@V≈0| / max|I| < 1%).

## Output Structure

```
output/
├── model_single/
│   ├── sample_fitting.svg / .pdf / .png / _600dpi.png
│   ├── sample_component_fit.txt    (V, J_data, J_fit, J_main, J_Ohm, J_TAT)
│   ├── sample_params.csv
│   └── sample_report.docx
├── model_dual/
│   ├── sample_fitting.svg / .pdf / .png / _600dpi.png
│   ├── sample_component_fit.txt    (V, J_data, J_fit, J_diff, J_rec, J_Ohm, J_TAT)
│   ├── sample_params.csv
│   └── sample_report.docx
```

## Plot Formatting (Journal-Ready)

- **Font**: Arial bold axis labels
- **Ticks**: Inward, no grid, no legend frame
- **Colors** (Wong 2011 color-blind friendly): Data `#0072B2`, Fit `#D55E00`, J<sub>diff</sub>/J<sub>main</sub> `#009E73`, J<sub>rec</sub> `#F0E442`, J<sub>Ohm</sub> `#CC79A7`, J<sub>TAT</sub> `#56B4E9`
- **Y-axis**: semilogy, data range −1.5 to +0.6 decades
- **X-axis**: Exact data range, no padding
- **Export**: SVG + PDF + PNG (300 dpi) + PNG (600 dpi with `--hd`)
- **Width**: 8.5 cm (single column)

## Academic Word Report (.docx)

Generated automatically with:
- Embedded PNG fitting figure
- Model equation
- Component-by-component physical interpretation
- Voltage-dependent dominance analysis (J<sub>diff</sub>/J<sub>main</sub> → J<sub>TAT</sub> crossover)
- Full parameter table
- Device optimization implications

## Dependencies

- **numpy** — numerical computation
- **scipy** — curve fitting (TRF method)
- **matplotlib** — publication-quality plots
- **pandas** — CSV export
- **python-docx** — Word report generation

## License

MIT
