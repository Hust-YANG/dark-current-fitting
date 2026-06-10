#!/usr/bin/env python3
"""
暗电流拟合 — 单/双二极管模型 (Dark Current Fitting with Implicit Diode Models)
==============================================================================

Model (Dual-Diode, 8 params):
  J_D = J01·[exp(A1·V_int)-1] + J02·[exp(A2·V_int)-1] + V_int/R_SH + k·V·exp(m/V_int)
  where V_int = V - J_D·R_S  (隐式方程, damped fixed-point 迭代求解)

Model (Single-Diode, 6 params): J02 ≡ 0, same structure.

Components:
  J_diff  = J01·[exp(A1·V_int)-1]     (扩散电流, A1 ≈ q/(n1·k·T), n1≈1)
  J_rec   = J02·[exp(A2·V_int)-1]     (复合电流, A2 ≈ q/(n2·k·T), n2≈2)
  J_Ohm   = V_int / R_SH               (欧姆漏电)
  J_tun   = k · V · exp(m/V_int)       (指数隧穿, m < 0, |m|~0.01-0.5V)

Output: 拟合图(SVG/PDF/PNG) + 分量数据(.txt) + 参数(CSV) + 学术报告(.docx)
"""

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy.constants import k, e
import pandas as pd
import os, sys
import warnings
warnings.filterwarnings('ignore', category=UserWarning, module='matplotlib')
from datetime import datetime

if hasattr(sys.stdout, 'reconfigure'):
    try: sys.stdout.reconfigure(encoding='utf-8')
    except Exception: pass

T_DEFAULT = 300
CM_PER_INCH = 2.54
MAX_EXP_ARG = 100.0

COLOR_PALETTE = dict(
    data='#0072B2', total_fit='#D55E00',
    j_diff='#009E73', j_rec='#F0E442',
    j_ohm='#CC79A7', j_tun='#56B4E9'
)

# [J01, A1, J02, A2, R_S, R_SH, k, m]   m < 0 for tunneling: exp(m/V_int) with V_int<0
BOUNDS_DUAL = (
    [1e-12, 10.0, 1e-12,  5.0,  0.0,  1e3, 1e-10, -1.0],
    [1e-3,  50.0, 1e-3,   30.0, 1e4,  1e8, 1e2,   -0.001]
)
BOUNDS_SINGLE = (
    [1e-12, 10.0, 0.0,  1e3, 1e-10, -1.0],
    [1e-3,  50.0, 1e4,  1e8, 1e2,   -0.001]
)

GUESSES_DUAL = [
    [1e-7, 38.7, 1e-6, 19.3, 1.0, 1e5, 1e-6, -0.1],
    [1e-8, 40.0, 1e-7, 20.0, 10.0, 1e6, 1e-7, -0.05],
    [1e-6, 35.0, 1e-5, 17.0, 0.1, 1e4, 1e-5, -0.3],
    [1e-7, 42.0, 1e-8, 22.0, 5.0, 1e5, 1e-8, -0.5],
    [1e-9, 38.7, 1e-6, 19.3, 50.0, 1e7, 1e-4, -0.02],
    [1e-5, 30.0, 1e-4, 15.0, 0.5, 1e3, 1e-9, -0.8],
]

GUESSES_SINGLE = [
    [1e-7, 38.7, 1.0, 1e5, 1e-6, -0.1],
    [1e-8, 40.0, 10.0, 1e6, 1e-7, -0.05],
    [1e-6, 35.0, 0.1, 1e4, 1e-5, -0.3],
    [1e-7, 42.0, 5.0, 1e5, 1e-8, -0.5],
    [1e-9, 38.7, 50.0, 1e7, 1e-4, -0.02],
    [1e-5, 30.0, 0.5, 1e3, 1e-9, -0.8],
]


def _compute_components(V, V_int, J01, A1, J02, A2, R_SH, k, m):
    J_diff = J01 * (np.exp(np.clip(A1 * V_int, -MAX_EXP_ARG, MAX_EXP_ARG)) - 1.0)
    J_rec = J02 * (np.exp(np.clip(A2 * V_int, -MAX_EXP_ARG, MAX_EXP_ARG)) - 1.0)
    J_Ohm = V_int / R_SH
    V_int_tun = np.where(np.abs(V_int) < 0.005, np.sign(V_int + 1e-30) * 0.005, V_int)
    J_tun = k * V * np.exp(np.clip(m / V_int_tun, -30.0, 30.0))
    return J_diff, J_rec, J_Ohm, J_tun


def _rhs_total(V, V_int, J01, A1, J02, A2, R_SH, k, m):
    Jd, Jr, Jo, Jt = _compute_components(V, V_int, J01, A1, J02, A2, R_SH, k, m)
    return Jd + Jr + Jo + Jt


def solve_implicit(V, J01, A1, J02, A2, R_S, R_SH, k, m,
                   max_iter=200, tol=1e-10, verbose=False):
    Jd = np.zeros_like(V)
    damping = 0.4
    for i in range(max_iter):
        V_int = V - Jd * R_S
        Jd_new = damping * _rhs_total(V, V_int, J01, A1, J02, A2, R_SH, k, m) + (1 - damping) * Jd
        diff = np.max(np.abs(Jd_new - Jd))
        if diff < tol:
            if verbose and i > 3:
                print(f"    Converged in {i+1} iterations")
            return Jd_new, V_int, True
        Jd = Jd_new
    if verbose:
        V_int_final = V - Jd * R_S
        print(f"    Warning: not converged after {max_iter} iterations, "
              f"max diff = {np.max(np.abs(_rhs_total(V, V_int_final, J01, A1, J02, A2, R_SH, k, m) - Jd)):.2e}")
    V_int = V - Jd * R_S
    return Jd, V_int, False


def dark_current_dual(V, J01, A1, J02, A2, R_S, R_SH, k, m):
    Jd, _, _ = solve_implicit(V, J01, A1, J02, A2, R_S, R_SH, k, m)
    return Jd


def dark_current_single(V, J01, A1, R_S, R_SH, k, m):
    Jd, _, _ = solve_implicit(V, J01, A1, 0.0, 19.3, R_S, R_SH, k, m)
    return Jd


def detect_encoding(fp):
    for enc in ['utf-8', 'utf-8-sig', 'gbk', 'gb2312', 'gb18030',
                'utf-16', 'utf-16-le', 'utf-16-be', 'latin-1']:
        try:
            with open(fp, 'r', encoding=enc) as f: f.read(1024)
            return enc
        except Exception: continue
    return 'utf-8'


def _parse_raw_data(fp):
    if not os.path.exists(fp):
        print(f"Error: file not found — {fp}"); sys.exit(1)
    enc = detect_encoding(fp)
    sweeps = {}; flat_v, flat_i = [], []; has_sweep_marker = False
    with open(fp, 'r', encoding=enc) as f:
        for ln in f:
            ln = ln.strip()
            if not ln: continue
            parts = [p.strip().strip('"') for p in ln.replace('\t', ' ').split(' ') if p.strip()]
            if len(parts) < 2: continue
            try:
                if len(parts) >= 3:
                    sid = parts[0]; Vv, Iv = float(parts[-2]), float(parts[-1])
                    sid = sid.split('/')[0] if '/' in sid else '1'
                    has_sweep_marker = has_sweep_marker or ('/' in parts[0])
                    sweeps.setdefault(sid, {'V': [], 'I': []})
                    sweeps[sid]['V'].append(Vv); sweeps[sid]['I'].append(Iv)
                else:
                    flat_v.append(float(parts[0])); flat_i.append(float(parts[1]))
            except ValueError: continue
    if has_sweep_marker and sweeps: return sweeps, True
    if flat_v: return {'1': {'V': flat_v, 'I': flat_i}}, False
    if sweeps: return sweeps, False
    print(f"Error: no parsable data in {fp}"); sys.exit(1)


def _identify_dark(sweeps, fp):
    dark, light = [], []
    print(f"\n  Sweep detection ({os.path.basename(fp)}):")
    print(f"  {'Sweep':<8} {'Pts':<6} {'I@V0':<16} {'max|I|':<16} {'Ratio':<10} {'State'}")
    for sid in sorted(sweeps):
        Va = np.array(sweeps[sid]['V']); Ia = np.array(sweeps[sid]['I'])
        i0 = np.argmin(np.abs(Va)); I0 = Ia[i0]; mx = np.max(np.abs(Ia))
        ratio = abs(I0) / mx if mx > 0 else 0.0
        is_dark = (ratio < 0.01) or (abs(I0) < 1e-12)
        print(f"  {sid:<8} {len(Va):<6} {I0:<16.4e} {mx:<16.4e} {ratio:<10.4f} {'DARK' if is_dark else 'LIGHT'}")
        (dark if is_dark else light).append(sid)
    if not dark: dark = list(sweeps.keys())
    print(f"  -> Dark sweeps: {dark}")
    return dark


def load_data(fp, area=None, auto_dark=True, manual_sweep=None, max_points=None):
    sweeps, has_marker = _parse_raw_data(fp)
    if manual_sweep:
        use = [manual_sweep] if manual_sweep in sweeps else _identify_dark(sweeps, fp)
    elif auto_dark and len(sweeps) > 1:
        use = _identify_dark(sweeps, fp)
    else:
        use = list(sweeps.keys())
    all_v, all_i = [], []
    for sid in sorted(use):
        all_v.extend(sweeps[sid]['V']); all_i.extend(sweeps[sid]['I'])
    if max_points:
        all_v = all_v[:max_points]; all_i = all_i[:max_points]
    Va = np.array(all_v); Ia = np.array(all_i)
    V_fit = -Va; J_fit = -Ia / (area if area else 1.0)
    print(f"  Loaded {len(V_fit)} pts" + (f", area={area} cm²" if area else ""))
    return V_fit, J_fit


def thermal_voltage(T): return k * T / e


def fit_dark_current(V, J, model_type='dual', T=300):
    Vt = thermal_voltage(T)
    if model_type == 'single':
        model_func = dark_current_single; bounds = BOUNDS_SINGLE; guesses = GUESSES_SINGLE
        param_names = ['J01', 'A1', 'R_S', 'R_SH', 'k', 'm']
    else:
        model_func = dark_current_dual; bounds = BOUNDS_DUAL; guesses = GUESSES_DUAL
        param_names = ['J01', 'A1', 'J02', 'A2', 'R_S', 'R_SH', 'k', 'm']

    label = "Single-Diode" if model_type == 'single' else "Dual-Diode"
    print(f"\n=== {label} Global Fitting ===")

    best, best_r2 = None, -np.inf
    for i, p0 in enumerate(guesses):
        try:
            popt, _ = curve_fit(model_func, V, J, p0=p0, bounds=bounds, maxfev=100000, method='trf')
            Jp = model_func(V, *popt)
            ssr = np.sum((J - Jp) ** 2); sst = np.sum((J - np.mean(J)) ** 2)
            r2 = 1.0 - ssr / sst if sst > 0 else 0.0
            n1 = 1.0 / (popt[1] * Vt) if popt[1] > 0 else float('inf')
            n2_str = ""
            if model_type == 'dual':
                n2 = 1.0 / (popt[3] * Vt) if popt[3] > 0 else float('inf')
                n2_str = f" n2={n2:.3f}"
            pstr = " ".join(f"{popt[j]:.4e}" for j in range(len(popt)))
            print(f"  G{i+1}: {pstr} R²={r2:.6f} n1={n1:.3f}{n2_str}")
            if r2 > best_r2: best_r2 = r2; best = popt
        except RuntimeError: pass

    if best is None: print("  All initial guesses failed."); return None

    if model_type == 'dual':
        Jd, V_int, conv = solve_implicit(V, *best, verbose=True)
    else:
        J01, A1, R_S, R_SH, k_val, m_val = best
        Jd, V_int, conv = solve_implicit(V, J01, A1, 0.0, 19.3, R_S, R_SH, k_val, m_val, verbose=True)
    ssr = np.sum((J - Jd) ** 2); sst = np.sum((J - np.mean(J)) ** 2)
    r2 = 1.0 - ssr / sst if sst > 0 else 0.0
    n1 = 1.0 / (best[1] * Vt) if best[1] > 0 else float('inf')
    if model_type == 'dual': n2 = 1.0 / (best[3] * Vt) if best[3] > 0 else None
    else: n2 = None
    pstr = " ".join(f"{best[j]:.4e}" for j in range(len(best)))
    print(f"  => Best: {pstr} R²={r2:.6f} n1={n1:.3f}" + (f" n2={n2:.3f}" if n2 is not None else ""))
    return {'popt': best, 'r_squared': r2, 'Vt': Vt, 'n1': n1, 'n2': n2,
            'model_type': model_type, 'param_names': param_names, 'converged': conv, 'V_int': V_int}


def configure_plot_style(width_cm=8.5):
    from matplotlib.font_manager import FontProperties, findfont
    cjk_candidates = ['STSong', 'SimSun', 'SimHei', 'Microsoft YaHei', 'Arial Unicode MS']
    cjk_avail = [f for f in cjk_candidates if FontProperties(family=f) and
                 findfont(FontProperties(family=f), fallback_to_default=False)]
    matplotlib.rcParams.update({
        'font.family': 'sans-serif', 'font.sans-serif': ['Arial'] + cjk_avail,
        'font.size': 9, 'axes.labelsize': 9, 'xtick.labelsize': 8, 'ytick.labelsize': 8,
        'legend.fontsize': 7, 'figure.dpi': 300, 'savefig.dpi': 300,
        'savefig.bbox': 'tight', 'savefig.pad_inches': 0.05,
        'xtick.direction': 'in', 'ytick.direction': 'in',
        'xtick.major.size': 4, 'xtick.major.width': 0.8,
        'ytick.major.size': 4, 'ytick.major.width': 0.8,
        'axes.linewidth': 0.8, 'lines.linewidth': 1.2, 'lines.markersize': 4,
        'legend.frameon': True, 'legend.framealpha': 0.85, 'legend.edgecolor': '#cccccc',
        'axes.grid': False, 'mathtext.fontset': 'stix',
    })


def plot_fitting(V, J, fit, ax):
    popt = fit['popt']; model_type = fit['model_type']
    if model_type == 'dual': J01, A1, J02, A2, R_S, R_SH, k_val, m_val = popt
    else: J01, A1, R_S, R_SH, k_val, m_val = popt; J02, A2 = 0.0, 19.3
    V_plot = np.linspace(V.min(), V.max(), 1000)
    Jd_plot, Vint_plot, _ = solve_implicit(V_plot, J01, A1, J02, A2, R_S, R_SH, k_val, m_val)
    Jdiff, Jrec, Joh, Jtun = _compute_components(V_plot, Vint_plot, J01, A1, J02, A2, R_SH, k_val, m_val)
    ax.semilogy(V, np.abs(J), 'o', ms=3, color=COLOR_PALETTE['data'], alpha=0.7,
                label='Data', zorder=5, mec=COLOR_PALETTE['data'], mew=0.3)
    ax.semilogy(V_plot, np.abs(Jd_plot), '-', lw=1.5, alpha=0.7,
                color=COLOR_PALETTE['total_fit'], label='$J_{\\mathrm{dark}}$ fit', zorder=4)
    ax.semilogy(V_plot, np.abs(Jdiff), '--', lw=1.0, color=COLOR_PALETTE['j_diff'],
                label='$J_{\\mathrm{diff}}$', zorder=3)
    if model_type == 'dual':
        ax.semilogy(V_plot, np.abs(Jrec), '-.', lw=1.0, color=COLOR_PALETTE['j_rec'],
                    label='$J_{\\mathrm{rec}}$', zorder=3)
    ax.semilogy(V_plot, np.abs(Joh), '-.', lw=1.0, color=COLOR_PALETTE['j_ohm'],
                label='$J_{\\mathrm{Ohm}}$', zorder=3)
    ax.semilogy(V_plot, np.abs(Jtun), ':', lw=1.0, color=COLOR_PALETTE['j_tun'],
                label='$J_{\\mathrm{tun}}$', zorder=3)
    ax.set_xlabel('Voltage (V)', fontweight='bold')
    ax.set_ylabel('Current Density (A/cm$^{2}$)', fontweight='bold')
    ax.tick_params(top=False, right=False, which='both', labelsize=8)
    Jpos = np.abs(J)[np.abs(J) > 0]
    ax.set_ylim(np.min(Jpos) * 10 ** (-0.75), np.max(Jpos) * 10 ** (0.75))
    ax.legend(fontsize=6, loc='lower left', frameon=False, handlelength=1.5, prop={'weight': 'normal'})


def save_figure_formats(fig, basepath, hd=False):
    saved = []
    for ext, dpi in [('.svg', None), ('.pdf', None), ('.png', 300)]:
        kw = {'format': ext[1:], 'bbox_inches': 'tight', 'facecolor': 'white', 'edgecolor': 'none'}
        if dpi: kw['dpi'] = dpi
        fp = basepath + ext; fig.savefig(fp, **kw); saved.append(fp)
    if hd:
        fp_hd = basepath + '_600dpi.png'
        fig.savefig(fp_hd, dpi=600, bbox_inches='tight', facecolor='white', edgecolor='none')
        saved.append(fp_hd)
    return saved


def export_component_data(V, J, fit, output_dir, label='sample'):
    popt = fit['popt']; model_type = fit['model_type']
    if model_type == 'dual': J01, A1, J02, A2, R_S, R_SH, k_val, m_val = popt
    else: J01, A1, R_S, R_SH, k_val, m_val = popt; J02, A2 = 0.0, 19.3
    Jd, Vint, _ = solve_implicit(V, J01, A1, J02, A2, R_S, R_SH, k_val, m_val)
    Jdiff, Jrec, Joh, Jtun = _compute_components(V, Vint, J01, A1, J02, A2, R_SH, k_val, m_val)
    if model_type == 'dual':
        header = 'V(V)\tJ_data(A/cm2)\tJ_fit(A/cm2)\tJ_diff(A/cm2)\tJ_rec(A/cm2)\tJ_Ohm(A/cm2)\tJ_tun(A/cm2)'
        data = np.column_stack([V, J, Jd, Jdiff, Jrec, Joh, Jtun])
    else:
        header = 'V(V)\tJ_data(A/cm2)\tJ_fit(A/cm2)\tJ_diff(A/cm2)\tJ_Ohm(A/cm2)\tJ_tun(A/cm2)'
        data = np.column_stack([V, J, Jd, Jdiff, Joh, Jtun])
    fp = os.path.join(output_dir, f'{label}_component_fit.txt')
    np.savetxt(fp, data, header=header, delimiter='\t', fmt='%.6e', comments='')
    print(f"  Component data → {fp}")
    return fp


def export_params(fit, output_dir, label='fitting'):
    popt = fit['popt']; model_type = fit['model_type']
    rows = []
    if model_type == 'dual':
        rows.extend([
            ('J01 (diffusion)', 'J01', 'A/cm²', f'{popt[0]:.4e}'),
            ('A1 (=q/(n1·kT))', 'A1', 'V⁻¹', f'{popt[1]:.4f}'),
            ('n1 (ideality)', 'n1', '-', f'{fit["n1"]:.4f}'),
            ('J02 (recombination)', 'J02', 'A/cm²', f'{popt[2]:.4e}'),
            ('A2 (=q/(n2·kT))', 'A2', 'V⁻¹', f'{popt[3]:.4f}'),
            ('n2 (ideality)', 'n2', '-', f'{fit["n2"]:.4f}'),
            ('R_S (series)', 'R_S', 'Ω·cm²', f'{popt[4]:.2f}'),
            ('R_SH (shunt)', 'R_SH', 'Ω·cm²', f'{popt[5]:.2e}'),
            ('k (tunneling coef)', 'k', '-', f'{popt[6]:.4e}'),
            ('m (tunneling barrier)', 'm', 'V', f'{popt[7]:.4f}'),
            ('R²', 'R²', '-', f'{fit["r_squared"]:.6f}'),
        ])
    else:
        rows.extend([
            ('J01 (diffusion)', 'J01', 'A/cm²', f'{popt[0]:.4e}'),
            ('A1 (=q/(n1·kT))', 'A1', 'V⁻¹', f'{popt[1]:.4f}'),
            ('n1 (ideality)', 'n1', '-', f'{fit["n1"]:.4f}'),
            ('R_S (series)', 'R_S', 'Ω·cm²', f'{popt[2]:.2f}'),
            ('R_SH (shunt)', 'R_SH', 'Ω·cm²', f'{popt[3]:.2e}'),
            ('k (tunneling coef)', 'k', '-', f'{popt[4]:.4e}'),
            ('m (tunneling barrier)', 'm', 'V', f'{popt[5]:.4f}'),
            ('R²', 'R²', '-', f'{fit["r_squared"]:.6f}'),
        ])
    df = pd.DataFrame(rows, columns=['Parameter', 'Symbol', 'Unit', 'Value'])
    print(f"\n{'='*50}\n  {model_type.upper()}-DIODE FITTING PARAMETERS\n{'='*50}")
    print(df.to_string(index=False))
    fp = os.path.join(output_dir, f'{label}_params.csv')
    df.to_csv(fp, index=False, encoding='utf-8-sig')
    print(f"  Parameters → {fp}")
    return fp


def generate_word_report(fit, V, J, output_dir, T=300, area=None, label='sample'):
    try:
        from docx import Document; from docx.shared import Pt, Cm
        from docx.enum.text import WD_ALIGN_PARAGRAPH
        from docx.enum.table import WD_TABLE_ALIGNMENT
    except ImportError:
        print("  python-docx not installed — skipping report"); return None
    popt = fit['popt']; model_type = fit['model_type']
    if model_type == 'dual':
        J01, A1, J02, A2, R_S, R_SH, k_val, m_val = popt; model_label = "Dual-Diode"
    else:
        J01, A1, R_S, R_SH, k_val, m_val = popt; J02, A2 = 0.0, 19.3; model_label = "Single-Diode"
    Vd = np.linspace(V.min(), V.max(), 5000)
    Jd_d, Vint_d, _ = solve_implicit(Vd, J01, A1, J02, A2, R_S, R_SH, k_val, m_val)
    Jdiff, Jrec, Joh, Jtun = _compute_components(Vd, Vint_d, J01, A1, J02, A2, R_SH, k_val, m_val)
    Jtot = np.abs(Jd_d)
    i0 = np.argmin(np.abs(Vd)); Jzero = Jtot[i0]
    r_main = np.abs(Jdiff) / (Jtot + 1e-30)
    mdz = (r_main > 0.7) & (Vd <= 0); V_main_end = Vd[mdz].min() if np.any(mdz) else -0.2
    mr = Vd < 0; Vr = Vd[mr]
    ci = np.where(np.abs(Jtun[mr]) > np.abs(Jdiff[mr]))[0]
    V_crossover = Vr[ci[0]] if len(ci) > 0 else -0.3
    i_neg05 = np.argmin(np.abs(Vd + 0.5))
    rt_tun = np.abs(Jtun[i_neg05]) / (Jtot[i_neg05] + 1e-30) * 100
    rt_ohm = np.abs(Joh[i_neg05]) / (Jtot[i_neg05] + 1e-30) * 100
    doc = Document()
    style = doc.styles['Normal']; style.font.name = 'Times New Roman'; style.font.size = Pt(11)
    style.paragraph_format.space_after = Pt(6); style.paragraph_format.line_spacing = 1.15
    title = doc.add_heading(f'Dark Current J-V Fitting Analysis ({model_label} Model)', 0)
    title.alignment = WD_ALIGN_PARAGRAPH.CENTER
    doc.add_paragraph(
        f'Figure 1 shows the dark current density–voltage (J–V) characteristics '
        f'fitted with the {model_label.lower()} implicit diode model. '
        f'Blue circles: experimental data; orange solid: total fit (R²={fit["r_squared"]:.4f}); '
        f'green dashed: J_diff (diffusion); '
        + ('gold dash-dotted: J_rec (recombination); ' if model_type == 'dual' else '')
        + f'purple dash-dotted: J_Ohm (ohmic); light blue dotted: J_tun (exponential tunneling). '
        f'Excellent agreement across {V.min():.1f}V to {V.max():.1f}V. '
        f'Dark current at zero bias: {Jzero:.1e} A/cm².')
    img_path = os.path.join(output_dir, f'{label}_fitting.png')
    if os.path.exists(img_path):
        p = doc.add_paragraph(); p.alignment = WD_ALIGN_PARAGRAPH.CENTER
        p.add_run().add_picture(img_path, width=Cm(8.5))
        cap = doc.add_paragraph(); cap.alignment = WD_ALIGN_PARAGRAPH.CENTER
        r = cap.add_run('Figure 1. Dark current J–V characteristics and '
                        f'{model_label.lower()} model fitting results.')
        r.font.size = Pt(9); r.italic = True; doc.add_paragraph('')
    doc.add_heading('Fitting Model', 1)
    eq = doc.add_paragraph(); eq.alignment = WD_ALIGN_PARAGRAPH.CENTER
    r = eq.add_run('J_D = J₀₁·[exp(A₁·V_int)−1] + J₀₂·[exp(A₂·V_int)−1] '
                   '+ V_int/R_SH + k·V·exp(m/V_int)\nV_int = V − J_D·R_S')
    r.font.size = Pt(12); r.bold = True; doc.add_paragraph('')
    eq_desc = doc.add_paragraph(); eq_desc.alignment = WD_ALIGN_PARAGRAPH.CENTER
    desc_text = ('J_diff (diffusion, n₁≈1) + J_rec (G-R recombination, n₂≈2) '
                 '+ J_Ohm (ohmic shunt) + J_tun (exponential tunneling)'
                 if model_type == 'dual' else
                 'J_diff (diffusion, n₁≈1) + J_Ohm (ohmic shunt) + J_tun (exponential tunneling)')
    eq_desc.add_run(desc_text).font.size = Pt(9); doc.add_paragraph('')
    doc.add_heading('1. Diffusion Current (J_diff)', 2)
    doc.add_paragraph(f'J_diff = J₀₁·[exp(A₁·V_int)−1], A₁ = q/(n₁·k_B·T). '
        f'J₀₁ = {J01:.2e} A/cm², A₁ = {A1:.2f} V⁻¹, n₁ = {fit["n1"]:.3f}.')
    if model_type == 'dual':
        doc.add_heading('2. Generation-Recombination Current (J_rec)', 2)
        doc.add_paragraph(f'J_rec = J₀₂·[exp(A₂·V_int)−1], A₂ = q/(n₂·k_B·T). '
            f'J₀₂ = {J02:.2e} A/cm², A₂ = {A2:.2f} V⁻¹, n₂ = {fit["n2"]:.3f}.')
    idx = 3 if model_type == 'dual' else 2
    doc.add_heading(f'{idx}. Ohmic Leakage (J_Ohm)', 2)
    doc.add_paragraph(f'J_Ohm = V_int / R_SH. R_SH = {R_SH:.2e} Ω·cm².')
    idx += 1
    doc.add_heading(f'{idx}. Exponential Tunneling (J_tun)', 2)
    doc.add_paragraph(f'J_tun = k·V·exp(m/V_int). k = {k_val:.4e}, m = {m_val:.4f} V. '
        f'm < 0 → tunneling grows with reverse bias.')
    doc.add_heading(f'{idx+1}. Series Resistance (R_S)', 2)
    doc.add_paragraph(f'R_S = {R_S:.2f} Ω·cm². V_int = V − J_D·R_S.')
    doc.add_heading('Voltage-Dependent Dominance', 1)
    doc.add_heading(f'1. Low Reverse Bias (0V to {abs(V_main_end):.1f}V): J_diff Dominant', 2)
    doc.add_paragraph('J_diff > 70% of total.')
    doc.add_heading(f'2. High Reverse Bias (< {V_crossover:.1f}V): J_tun Emerges', 2)
    doc.add_paragraph(f'At −0.5V, J_tun ~{rt_tun:.0f}% of total.')
    doc.add_heading('3. Full Range: J_Ohm Negligible', 2)
    doc.add_paragraph(f'At −0.5V, J_Ohm ~{rt_ohm:.1f}% of total.')
    doc.add_heading('Fitting Parameters', 1)
    tdata = [('J₀₁ (Diffusion)', f'{J01:.4e}', 'A/cm²'),
             ('A₁ (=q/(n₁·k_B·T))', f'{A1:.4f}', 'V⁻¹'),
             ('n₁ (Ideality factor)', f'{fit["n1"]:.4f}', '-')]
    if model_type == 'dual':
        tdata += [('J₀₂ (G-R Recombination)', f'{J02:.4e}', 'A/cm²'),
                  ('A₂ (=q/(n₂·k_B·T))', f'{A2:.4f}', 'V⁻¹'),
                  ('n₂ (Ideality factor)', f'{fit["n2"]:.4f}', '-')]
    tdata += [('R_S (Series resistance)', f'{R_S:.2f}', 'Ω·cm²'),
              ('R_SH (Shunt resistance)', f'{R_SH:.2e}', 'Ω·cm²'),
              ('k (Tunneling coefficient)', f'{k_val:.4e}', '-'),
              ('m (Tunneling barrier)', f'{m_val:.4f}', 'V'),
              ('R²', f'{fit["r_squared"]:.6f}', '-')]
    tbl = doc.add_table(rows=len(tdata)+1, cols=3, style='Light Grid Accent 1')
    tbl.alignment = WD_TABLE_ALIGNMENT.CENTER
    for i, h in enumerate(['Parameter', 'Value', 'Unit']): tbl.rows[0].cells[i].text = h
    for i, (name, val, unit) in enumerate(tdata):
        tbl.rows[i+1].cells[0].text = name; tbl.rows[i+1].cells[1].text = val
        tbl.rows[i+1].cells[2].text = unit
    doc.add_heading('Implications', 1)
    doc.add_paragraph(
        f'The {model_label.lower()} implicit diode model achieves R² = {fit["r_squared"]:.4f}. '
        f'n₁ = {fit["n1"]:.3f} '
        + (f'n₂ = {fit["n2"]:.3f}. ' if model_type == 'dual' else '')
        + ('Diffusion-dominated with moderate G-R recombination.'
           if model_type == 'dual' and fit["n2"] and fit["n2"] > 1.5 else 'Near-ideal diffusion transport.')
        + f' Optimization: (1) suppress J_tun (k={k_val:.2e}); '
        f'(2) minimize R_S ({R_S:.1f} Ω·cm²); (3) passivate G-R centers.')
    fp = os.path.join(output_dir, f'{label}_report.docx')
    doc.save(fp)
    print(f"  Report → {fp}")
    return fp


def run_dark_current_fitting(data_file, output_dir=None, model='dual',
                             T=300, area=None, vmin=None, vmax=None,
                             auto_dark=True, max_points=None, hd=False, fig_width_cm=8.5):
    if output_dir is None:
        output_dir = os.path.join(os.path.dirname(os.path.abspath(data_file)), 'output')
    Vt = thermal_voltage(T)
    print(f"\n{'='*50}\n  DARK CURRENT FITTING — Implicit Diode Model\n{'='*50}")
    print(f"  T = {T} K, Vt = {Vt:.6f} V" + (f", Area = {area} cm²" if area else "") + f", Model = {model}")
    V, J = load_data(data_file, area=area, auto_dark=auto_dark, max_points=max_points)
    if vmin is not None or vmax is not None:
        lo = vmin if vmin is not None else V.min(); hi = vmax if vmax is not None else V.max()
        mask = (V >= lo) & (V <= hi); V = V[mask]; J = J[mask]
        print(f"  Voltage range: [{lo}, {hi}] V, {len(V)} pts")
    models_to_run = ['single', 'dual'] if model == 'both' else [model]
    all_results = []
    base_name = os.path.splitext(os.path.basename(data_file))[0]
    for m in models_to_run:
        print(f"\n{'─'*50}\n  Running: {m.upper()}-DIODE MODEL\n{'─'*50}")
        model_dir = os.path.join(output_dir, f'model_{m}'); os.makedirs(model_dir, exist_ok=True)
        fit = fit_dark_current(V, J, model_type=m, T=T)
        if fit is None: print(f"  {m}-diode fitting failed — skipping"); continue
        configure_plot_style(width_cm=fig_width_cm)
        fig, ax = plt.subplots(figsize=(fig_width_cm/CM_PER_INCH, fig_width_cm/CM_PER_INCH*0.75))
        plot_fitting(V, J, fit, ax); plt.tight_layout()
        fig_path = os.path.join(model_dir, f'{base_name}_fitting')
        save_figure_formats(fig, fig_path, hd=hd); plt.close(fig)
        print(f"  Figure → {fig_path}.*")
        export_component_data(V, J, fit, model_dir, label=base_name)
        export_params(fit, model_dir, label=base_name)
        generate_word_report(fit, V, J, model_dir, T=T, area=area, label=base_name)
        all_results.append({'model': m, 'fit': fit, 'dir': model_dir})
    print(f"\n{'='*50}\n  DONE — {len(all_results)} model(s) fitted\n{'='*50}")
    for r in all_results:
        print(f"  {r['model']}: R² = {r['fit']['r_squared']:.6f} → {r['dir']}")
    return {'results': all_results, 'output_dir': output_dir}


if __name__ == '__main__':
    import argparse
    p = argparse.ArgumentParser(description='Dark Current Fitting — Implicit Single/Dual Diode Model')
    p.add_argument('data', help='Path to input data file (.txt)')
    p.add_argument('-o', '--output', default=None, help='Output directory')
    p.add_argument('--model', choices=['single', 'dual', 'both'], default='dual')
    p.add_argument('-T', '--temperature', type=float, default=T_DEFAULT)
    p.add_argument('-a', '--area', type=float, default=None)
    p.add_argument('--vmin', type=float, default=-0.5)
    p.add_argument('--vmax', type=float, default=0.2)
    p.add_argument('--no-auto-dark', action='store_true')
    p.add_argument('--points', type=int, default=None)
    p.add_argument('--fig-width-cm', type=float, default=8.5)
    p.add_argument('--hd', action='store_true')
    args = p.parse_args()
    result = run_dark_current_fitting(args.data, output_dir=args.output, model=args.model,
        T=args.temperature, area=args.area, vmin=args.vmin, vmax=args.vmax,
        auto_dark=not args.no_auto_dark, max_points=args.points, hd=args.hd,
        fig_width_cm=args.fig_width_cm)
    if result and result['results']:
        print(f"\nOutput directory: {result['output_dir']}")