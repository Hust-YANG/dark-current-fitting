# Dark Current Fitting — PbS CQD Photodetector

A Python tool for quantitative multi-mechanism dark current component fitting of PbS colloidal quantum dot (CQD) infrared photodetectors. Decomposes measured J-V characteristics into three physically distinct transport components using segmented curve fitting with global refinement.

## Fitting Model

<div align="center">

**J<sub>dark</sub> = J<sub>0</sub>·[exp(qV/AkT) − 1]  +  V/R<sub>sh</sub>  +  B·V·exp[−C/(V<sub>bi</sub> − V)]**

*J<sub>main</sub> (Diffusion-Recombination)  +  J<sub>Ohm</sub> (Ohmic Leakage)  +  J<sub>TAT</sub> (Trap-Assisted Tunneling)*

</div>

### Three Current Components

| Component | Equation | Physical Origin | Dominant Region |
|-----------|----------|-----------------|-----------------|
| **J<sub>main</sub>** | J<sub>0</sub>·[exp(qV/AkT) − 1] | Minority carrier diffusion + G-R recombination in depletion region | Low reverse bias (0 ~ −0.2 V) |
| **J<sub>Ohm</sub>** | V/R<sub>sh</sub> | Ohmic conduction through film pinholes and grain boundaries | Full range (typically negligible) |
| **J<sub>TAT</sub>** | B·V·exp[−C/(V<sub>bi</sub> − V)] | Thermal excitation to trap states + field-assisted tunneling | High reverse bias (< −0.2 V) |

### Fitted Parameters

| Parameter | Symbol | Unit | Description |
|-----------|--------|------|-------------|
| Reverse saturation current density | J<sub>0</sub> | A/cm² | Intrinsic recombination activity |
| Ideality factor | A | – | Quality factor (A→1 ideal, A>1 trap-mediated) |
| Shunt resistance | R<sub>sh</sub> | Ω·cm² | Film ohmic leakage pathways |
| TAT defect density coefficient | B | – | Proportional to trap state density N<sub>t</sub> |
| TAT barrier coefficient | C | – | Effective tunneling barrier height |
| Built-in voltage | V<sub>bi</sub> | V | Junction internal electric field |

## Segmented Fitting Strategy

1. **Stage 1** (|V| < 0.2 V): J<sub>main</sub> dominant → fit **J<sub>0</sub>, A**
2. **Stage 2** (V < −0.2 V): J<sub>TAT</sub> + J<sub>Ohm</sub> dominant → fix J<sub>0</sub>, A → fit **R<sub>sh</sub>, B, C, V<sub>bi</sub>**
3. **Stage 3**: Global refinement — all 6 parameters co-adjusted within physical bounds

> J<sub>BTB</sub> (band-to-band tunneling) is explicitly excluded — PbS CQD films lack the high doping and narrow depletion region required.

## Sample Output
<p align="center">
<img width="600" alt="sample_fitting_600dpi" src="https://github.com/user-attachments/assets/f7cd9e9a-ad39-47c6-8df9-a6f7d8a57dd7" />
</p>
<p align="center">
Figure 1. Dark current J-V characteristics and multi-mechanism fitting results. 
</p>


## Usage
## Installation

```bash
pip install numpy scipy matplotlib pandas python-docx
```
```bash
python scripts/dark_current_fitting.py sample.txt -a 0.0707 --points 201 --hd
```

| Argument | Description | Default |
|----------|-------------|---------|
| `data` | Input data file (.txt) | Required |
| `-a, --area` | Device area (cm²) | None |
| `--points N` | Use first N data points | All |
| `--vmin` | Lower voltage limit (V) | -0.5 |
| `--vmax` | Upper voltage limit (V) | 0.2 |
| `--no-segmented` | Use global single-fit | False |
| `--fig-width-cm` | Figure width (cm) | 8.5 |
| `--hd` | Export 600 dpi PNG | False |
| `-o` | Output directory | ./output |

### Input Format

Tab-separated text file with columns: `Index  V1 (V)  I1 (A)`

```
Index    V1 (V)    I1 (A)
1/1      -1.0000E+00      -3.63062E-004
1/2      -9.900E-001      -3.60758E-004
...
```

Supports multi-sweep data with automatic dark-sweep detection (ratio < 1% at V≈0).

## Output

```
output/
├── sample_fitting.svg          # Vector (SVG)
├── sample_fitting.pdf          # Vector (PDF)  
├── sample_fitting.png          # 300 dpi raster
├── sample_fitting_600dpi.png   # 600 dpi (--hd)
├── sample_component_fit.txt    # Component data
├── fitting_params.csv          # Parameter table
└── fitting_report.docx         # Academic report with embedded figure
```

## Dependencies

- **numpy** — numerical computation
- **scipy** — curve fitting (TRF method)
- **matplotlib** — publication-quality plots
- **pandas** — data export
- **python-docx** — Word report generation

## License

MIT
