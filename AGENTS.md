# Dark Current Fitting Skill

Analyze Jdark-V of PbS CQD photodetectors using three-component TAT model.
Segmented fitting with journal-ready plots and academic Word report.

## Model

```
J_dark = J0*[exp(qV/(A*kT)) - 1] + V/Rsh + B*V*exp(-c/(Vbi-V))
```

## Usage

```bash
python scripts/dark_current_fitting.py sample.txt -a 0.07 --points 201 --hd
```

See SKILL.md for details.
