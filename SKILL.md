---
name: dark-current-fitting
description: >
  PbS CQD photodetector dark current analysis using three-component TAT model
  (J_main + J_Ohm + J_TAT). Segmented fitting with journal-ready plots and
  academic report generation. Auto-detects dark sweeps. Trigger on:
  暗电流拟合, dark current fitting, Jdark-V, TAT, trap-assisted tunneling,
  光电探测器, photodiode characterization.
---

# /dark-current-fitting

## Overview

Fits Jdark-V data using the three-component TAT model:

```
J_dark = J0*[exp(qV/(A*kT)) - 1] + V/Rsh + B*V*exp(-c/(Vbi-V))
         J_main                   J_Ohm       J_TAT
```

6 params: J0, A, Rsh, B, C, Vbi

## Segmented Fitting

- Stage 1 (|V| < V_seg): J0, A
- Stage 2 (V < -V_seg): Rsh, B, C, Vbi
- Stage 3: Global refinement

## Usage

```bash
python scripts/dark_current_fitting.py sample.txt -a 0.0707 --points 201 --hd
```

## Output

- Fitting plots (SVG/PDF/PNG)
- Component data (.txt)
- Parameters (CSV)
- Academic report (.docx)
