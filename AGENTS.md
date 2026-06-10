# Dark Current Fitting Skill

## Purpose

Analyze Jdark-V characteristics of PbS CQD photodetectors using dual-model fitting
(Model 1 Eq.1 + Model 2 Eq.2 TAT). Segmented fitting with journal-ready plots.

## Models

### Model 1 (Eq.1)
```
J_dark = J0*[exp(qV/(A*kT)) - 1] + V/Rsh + k*V^m
```

### Model 2 (Eq.2)
```
J_dark = J0*[exp(qV/(A*kT)) - 1] + V/Rsh + B*V*exp(-c/(Vbi-V))
```

## Usage

```bash
python scripts/dark_current_fitting.py <control.txt> <sample.txt> [options]
```

See SKILL.md for full documentation.
