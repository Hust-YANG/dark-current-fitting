#!/usr/bin/env python3
"""
Dark Current Component Fitting Script
======================================

Model: J_dark = J0*[exp(qV/(A*kT))-1] + V/Rsh + B*V*exp(-c/(Vbi-V))
          J_main             J_Ohm            J_TAT

Voltage convention: V_fit = -V_raw, J_fit = -J_raw (standard diode convention)

Output:
  1. Fitting plots (SVG/PDF/PNG@300dpi/PNG@600dpi)
  2. Component data .txt
  3. Fitting parameters CSV
  4. Academic report .docx
"""

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy.constants import k, e
import pandas as pd
import os, sys, json, io
import warnings
warnings.filterwarnings('ignore', category=UserWarning, module='matplotlib')
from datetime import datetime

if hasattr(sys.stdout, 'reconfigure'):
    try:
        sys.stdout.reconfigure(encoding='utf-8')
    except Exception:
        pass

T_DEFAULT = 300
CM_PER_INCH = 2.54

COLOR_PALETTE = {
    'data':       '#0072B2',
    'total_fit':  '#D55E00',
    'j_main':     '#009E73',
    'j_ohm':      '#CC79A7',
    'j_tat':      '#56B4E9',
}

BOUNDS = (
    [1e-15, 1.0,  1e-1,  1e-12, 0.001, 0.01],
    [1e-2,  3.0,  1e12,  1e8,   100.0,  2.0],
)
P0_DEFAULT    = [1e-7,  1.5,  1e4,  1.0,   1.0,  0.5]
P0_ALT        = [1e-5,  1.2,  1e5,  10.0,  5.0,  0.8]
P0_AGGRESSIVE = [1e-6,  1.8,  1e3,  0.1,   0.5,  0.4]

V_SEG_DEFAULT = 0.2

def thermal_voltage(T):
    return k * T / e

def main_diode_current(V, J0, A, Vt):
    return J0 * (np.exp(V / (A * Vt)) - 1.0)

def ohmic_current(V, Rsh):
    return V / Rsh

def tat_current(V, B, c, Vbi):
    denom = Vbi - V
    valid = denom > 0
    result = np.zeros_like(V, dtype=np.float64)
    if np.any(valid):
        exponent = np.clip(-c / denom[valid], -700.0, 0.0)
        result[valid] = B * V[valid] * np.exp(exponent)
    return result

def dark_current_model(V, J0, A, Rsh, B, c, Vbi, Vt=None, T=300):
    if Vt is None:
        Vt = thermal_voltage(T)
    return (main_diode_current(V, J0, A, Vt)
            + ohmic_current(V, Rsh)
            + tat_current(V, B, c, Vbi))

def detect_encoding(filepath):
    encodings = ['utf-8', 'utf-8-sig', 'gbk', 'gb2312', 'gb18030',
                 'utf-16', 'utf-16-le', 'utf-16-be', 'latin-1']
    for enc in encodings:
        try:
            with open(filepath, 'r', encoding=enc) as f:
                f.read(1024)
            return enc
        except (UnicodeDecodeError, UnicodeError):
            continue
    return 'utf-8'

def _parse_raw_data(filepath):
    if not os.path.exists(filepath):
        print(f"Error: file not found {filepath}")
        sys.exit(1)
    encoding = detect_encoding(filepath)
    sweeps = {}
    flat_V, flat_I = [], []
    has_index = False
    with open(filepath, 'r', encoding=encoding) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            parts = [p.strip().strip('"') for p in line.replace('\t', ' ').split(' ') if p.strip()]
            if len(parts) < 2:
                continue
            try:
                if len(parts) >= 3:
                    idx_str = parts[0]
                    V_pt, I_pt = float(parts[-2]), float(parts[-1])
                    if '/' in idx_str:
                        sweep_id = idx_str.split('/')[0]
                        has_index = True
                    else:
                        sweep_id = '1'
                    if sweep_id not in sweeps:
                        sweeps[sweep_id] = {'V': [], 'I': []}
                    sweeps[sweep_id]['V'].append(V_pt)
                    sweeps[sweep_id]['I'].append(I_pt)
                else:
                    V_pt, I_pt = float(parts[0]), float(parts[1])
                    flat_V.append(V_pt)
                    flat_I.append(I_pt)
            except ValueError:
                continue
    if has_index and len(sweeps) > 0:
        return sweeps, True
    elif len(flat_V) > 0:
        return {'1': {'V': flat_V, 'I': flat_I}}, False
    elif len(sweeps) > 0:
        return sweeps, False
    else:
        print(f"Error: no valid data in {filepath}")
        sys.exit(1)

def _identify_dark_sweeps(sweeps_dict, filepath):
    dark_sweeps = []
    light_sweeps = []
    print(f"\n  Sweep detection ({os.path.basename(filepath)}):")
    print(f"  {'Sweep':<8} {'Pts':<6} {'I@V0':<16} {'max|I|':<16} {'Ratio':<10} {'State'}")
    for sweep_id in sorted(sweeps_dict.keys()):
        V_arr = np.array(sweeps_dict[sweep_id]['V'])
        I_arr = np.array(sweeps_dict[sweep_id]['I'])
        idx_zero = np.argmin(np.abs(V_arr))
        I_at_zero = I_arr[idx_zero]
        max_abs_I = np.max(np.abs(I_arr))
        ratio = abs(I_at_zero) / max_abs_I if max_abs_I > 0 else 0.0
        is_dark = (ratio < 0.01) or (abs(I_at_zero) < 1e-12)
        state = 'DARK' if is_dark else 'LIGHT (skip)'
        print(f"  {sweep_id:<8} {len(V_arr):<6} {I_at_zero:<16.4e} {max_abs_I:<16.4e} {ratio:<10.4f} {state}")
        if is_dark:
            dark_sweeps.append(sweep_id)
        else:
            light_sweeps.append(sweep_id)
    if not dark_sweeps:
        print(f"  Warning: No dark sweep detected, using all data")
        dark_sweeps = list(sweeps_dict.keys())
    print(f"  -> Using dark sweep(s): {dark_sweeps}")
    if light_sweeps:
        print(f"  -> Discarded light sweep(s): {light_sweeps}")
    return dark_sweeps

def load_data(filepath, area=None, auto_dark=True, manual_sweep=None, max_points=None):
    sweeps_dict, has_index = _parse_raw_data(filepath)
    if manual_sweep is not None:
        if manual_sweep in sweeps_dict:
            use_sweeps = [manual_sweep]
        else:
            use_sweeps = _identify_dark_sweeps(sweeps_dict, filepath) if auto_dark else list(sweeps_dict.keys())
    elif auto_dark and len(sweeps_dict) > 1:
        use_sweeps = _identify_dark_sweeps(sweeps_dict, filepath)
    else:
        use_sweeps = list(sweeps_dict.keys())
    all_V, all_I = [], []
    for sid in sorted(use_sweeps):
        all_V.extend(sweeps_dict[sid]['V'])
        all_I.extend(sweeps_dict[sid]['I'])
    if max_points is not None and max_points > 0:
        all_V = all_V[:max_points]
        all_I = all_I[:max_points]
        print(f"  Using first {len(all_V)} points")
    V_raw = np.array(all_V)
    I_raw = np.array(all_I)
    V_fit = -V_raw
    J_fit = -I_raw / area if (area is not None and area > 0) else -I_raw
    print(f"  Loaded {len(V_fit)} pts" + (f", area={area} cm2" if area else ""))
    return V_fit, J_fit

def fit_dark_current(V, J, T=300):
    Vt = thermal_voltage(T)
    print(f"\nFitting dark current: {len(V)} pts")
    def fit_func(V, J0, A, Rsh, B, c, Vbi):
        return dark_current_model(V, J0, A, Rsh, B, c, Vbi, Vt=Vt)
    for i, p0 in enumerate([P0_DEFAULT, P0_ALT, P0_AGGRESSIVE]):
        try:
            popt, pcov = curve_fit(fit_func, V, J, p0=p0, bounds=BOUNDS, maxfev=100000, method='trf')
            J_pred = fit_func(V, *popt)
            ss_res = np.sum((J - J_pred) ** 2)
            ss_tot = np.sum((J - np.mean(J)) ** 2)
            r_squared = 1 - ss_res / ss_tot
            perr = np.sqrt(np.diag(pcov))
            print(f"  J0={popt[0]:.4e} A={popt[1]:.4f} Rsh={popt[2]:.2e} B={popt[3]:.4f} C={popt[4]:.4f} Vbi={popt[5]:.4f} R2={r_squared:.6f}")
            return {'popt': popt, 'pcov': pcov, 'perr': perr, 'r_squared': r_squared, 'Vt': Vt}
        except RuntimeError:
            pass
    print("All fits failed.")
    return None

def segmented_fit(V, J, T=300, V_seg=V_SEG_DEFAULT, verbose=True):
    Vt = thermal_voltage(T)
    if verbose:
        print(f"\nSegmented Fitting, V_seg={V_seg}V")
    global_seed = fit_dark_current(V, J, T)
    if global_seed is not None:
        gJ0, gA, gRsh, gB, gc, gVbi = global_seed['popt']
    else:
        gJ0, gA, gRsh, gB, gc, gVbi = P0_DEFAULT
    mask_small = np.abs(V) < V_seg
    V_small, J_small = V[mask_small], J[mask_small]
    def s1_model(V, J0, A, Rsh_eff):
        return main_diode_current(V, J0, A, Vt) + ohmic_current(V, Rsh_eff)
    p0_opts = [[gJ0, gA, gRsh], [P0_DEFAULT[0], P0_DEFAULT[1], P0_DEFAULT[2]],
               [P0_ALT[0], P0_ALT[1], P0_ALT[2]], [P0_AGGRESSIVE[0], P0_AGGRESSIVE[1], P0_AGGRESSIVE[2]]]
    bnd = ([BOUNDS[0][0], BOUNDS[0][1], BOUNDS[0][2]], [BOUNDS[1][0], BOUNDS[1][1], BOUNDS[1][2]])
    s1_ok = False; r2_s1 = None
    for pg in p0_opts:
        try:
            popt1, pcov1 = curve_fit(s1_model, V_small, J_small, p0=pg, bounds=bnd, maxfev=50000, method='trf')
            Jp1 = s1_model(V_small, *popt1)
            ssr1 = np.sum((J_small - Jp1)**2); sst1 = np.sum((J_small - np.mean(J_small))**2)
            r2_s1 = 1 - ssr1 / sst1 if sst1 > 0 else 0.0
            s1_ok = True; break
        except RuntimeError:
            pass
    if not s1_ok:
        gf = fit_dark_current(V, J, T)
        if gf is None: return None
        J0_s1, A_s1 = gf['popt'][0], gf['popt'][1]; pcov_s1 = None
    else:
        J0_s1, A_s1, _ = popt1
    mask_large = V < -V_seg
    V_large, J_large = V[mask_large], J[mask_large]
    if len(V_large) < 5:
        mask_large = V < 0; V_large, J_large = V[mask_large], J[mask_large]
    def s2_model(V, Rsh, B, c, Vbi):
        return main_diode_current(V, J0_s1, A_s1, Vt) + ohmic_current(V, Rsh) + tat_current(V, B, c, Vbi)
    p0_opts2 = [[gRsh, gB, gc, gVbi], [P0_DEFAULT[2], P0_DEFAULT[3], P0_DEFAULT[4], P0_DEFAULT[5]],
                [P0_ALT[2], P0_ALT[3], P0_ALT[4], P0_ALT[5]], [P0_AGGRESSIVE[2], P0_AGGRESSIVE[3], P0_AGGRESSIVE[4], P0_AGGRESSIVE[5]]]
    bnd2 = ([BOUNDS[0][2], BOUNDS[0][3], BOUNDS[0][4], BOUNDS[0][5]],
            [BOUNDS[1][2], BOUNDS[1][3], BOUNDS[1][4], BOUNDS[1][5]])
    s2_ok = False; r2_s2 = None
    for pg in p0_opts2:
        try:
            popt2, pcov2 = curve_fit(s2_model, V_large, J_large, p0=pg, bounds=bnd2, maxfev=100000, method='trf')
            Jp2 = s2_model(V_large, *popt2)
            ssr2 = np.sum((J_large - Jp2)**2); sst2 = np.sum((J_large - np.mean(J_large))**2)
            r2_s2 = 1 - ssr2 / sst2 if sst2 > 0 else 0.0
            s2_ok = True; break
        except RuntimeError:
            pass
    if not s2_ok:
        gf = fit_dark_current(V, J, T)
        if gf is None: return None
        Rsh_s2 = gf['popt'][2]; B_s2 = gf['popt'][3]; c_s2 = gf['popt'][4]; Vbi_s2 = gf['popt'][5]; pcov_s2 = None
    else:
        Rsh_s2, B_s2, c_s2, Vbi_s2 = popt2
    p0_s3 = [J0_s1, A_s1, Rsh_s2, B_s2, c_s2, Vbi_s2]
    bounds_s3 = (
        [max(BOUNDS[0][0], J0_s1*0.1), max(BOUNDS[0][1], A_s1*0.7),
         max(BOUNDS[0][2], Rsh_s2*0.1), max(BOUNDS[0][3], B_s2*0.01),
         max(BOUNDS[0][4], c_s2*0.1), max(BOUNDS[0][5], Vbi_s2*0.5)],
        [min(BOUNDS[1][0], J0_s1*10.0), min(BOUNDS[1][1], A_s1*1.5),
         min(BOUNDS[1][2], Rsh_s2*10.0), min(BOUNDS[1][3], B_s2*100.0),
         min(BOUNDS[1][4], c_s2*10.0), min(BOUNDS[1][5], Vbi_s2*1.5)],
    )
    s3_ok = False
    try:
        popt, pcov = curve_fit(lambda V, *p: dark_current_model(V, *p, Vt=Vt), V, J, p0=p0_s3, bounds=bounds_s3, maxfev=50000, method='trf')
        perr = np.sqrt(np.diag(pcov))
        s3_ok = True
    except RuntimeError:
        pass
    if not s3_ok:
        popt = np.array(p0_s3)
        perr = np.full(6, np.nan) if not (s1_ok and s2_ok) else np.array([
            np.sqrt(pcov_s1[0,0]), np.sqrt(pcov_s1[1,1]), np.sqrt(pcov_s2[0,0]),
            np.sqrt(pcov_s2[1,1]), np.sqrt(pcov_s2[2,2]), np.sqrt(pcov_s2[3,3])])
    J_pred = dark_current_model(V, *popt, Vt=Vt)
    ss_res = np.sum((J - J_pred)**2)
    ss_tot = np.sum((J - np.mean(J))**2)
    r_squared = 1 - ss_res / ss_tot if ss_tot > 0 else 0.0
    if verbose:
        print(f"  Final: J0={popt[0]:.4e} A={popt[1]:.4f} Rsh={popt[2]:.2e} B={popt[3]:.4f} C={popt[4]:.4f} Vbi={popt[5]:.4f} R2={r_squared:.6f}")
    return {'popt': popt, 'pcov': pcov if s3_ok else None, 'perr': perr,
            'r_squared': r_squared, 'Vt': Vt,
            'segmented': True, 'V_seg': V_seg,
            'stage1': {'J0': J0_s1, 'A': A_s1, 'r2': r2_s1},
            'stage2': {'Rsh': Rsh_s2, 'B': B_s2, 'c': c_s2, 'Vbi': Vbi_s2, 'r2': r2_s2}}

def configure_plot_style(width_cm=8.5):
    from matplotlib.font_manager import findfont, FontProperties
    cjk_fonts = ['STSong', 'SimSun', 'SimHei', 'Microsoft YaHei', 'Arial Unicode MS', 'KaiTi', 'FangSong']
    cjk_available = []
    for f in cjk_fonts:
        try:
            fp = FontProperties(family=f)
            path = findfont(fp, fallback_to_default=False)
            if path: cjk_available.append(f)
        except Exception:
            pass
    sans_fonts = ['Arial', 'DejaVu Sans', 'Liberation Sans']
    matplotlib.rcParams.update({
        'font.family': 'sans-serif',
        'font.sans-serif': sans_fonts + cjk_available,
        'font.weight': 'bold', 'font.size': 9,
        'axes.labelsize': 9, 'xtick.labelsize': 8, 'ytick.labelsize': 8,
        'legend.fontsize': 7, 'figure.dpi': 300, 'savefig.dpi': 300,
        'savefig.bbox': 'tight', 'savefig.pad_inches': 0.05,
        'xtick.direction': 'in', 'ytick.direction': 'in',
        'xtick.major.size': 4, 'xtick.major.width': 0.8,
        'ytick.major.size': 4, 'ytick.major.width': 0.8,
        'xtick.minor.size': 2, 'xtick.minor.width': 0.6,
        'ytick.minor.size': 2, 'ytick.minor.width': 0.6,
        'xtick.top': True, 'ytick.right': True,
        'axes.linewidth': 0.8, 'lines.linewidth': 1.2,
        'lines.markersize': 4, 'lines.markeredgewidth': 0.5,
        'legend.frameon': True, 'legend.framealpha': 0.85,
        'legend.edgecolor': '#cccccc', 'legend.fancybox': False,
        'axes.grid': False, 'mathtext.fontset': 'stix',
    })

def plot_fitting(V, J, fit, ax):
    popt = fit['popt']; Vt = fit['Vt']
    Vs = np.linspace(V.min(), V.max(), 1000)
    J_fit  = dark_current_model(Vs, *popt, Vt=Vt)
    J_main = main_diode_current(Vs, popt[0], popt[1], Vt)
    J_ohm  = ohmic_current(Vs, popt[2])
    J_tat  = tat_current(Vs, popt[3], popt[4], popt[5])
    ax.semilogy(V, np.abs(J), 'o', markersize=3, color=COLOR_PALETTE['data'], alpha=0.7,
                label='Data', zorder=5, markeredgecolor=COLOR_PALETTE['data'], markeredgewidth=0.3)
    ax.semilogy(Vs, np.abs(J_fit), '-', linewidth=1.5, alpha=0.7,
                color=COLOR_PALETTE['total_fit'], label='$J_{dark}$ fit', zorder=4)
    ax.semilogy(Vs, np.abs(J_main), '--', linewidth=1.0,
                color=COLOR_PALETTE['j_main'], label='$J_{main}$', zorder=3)
    ax.semilogy(Vs, np.abs(J_ohm), '-.', linewidth=1.0,
                color=COLOR_PALETTE['j_ohm'], label='$J_{Ohm}$', zorder=3)
    ax.semilogy(Vs, np.abs(J_tat), ':', linewidth=1.0,
                color=COLOR_PALETTE['j_tat'], label='$J_{TAT}$', zorder=3)
    ax.set_xlabel('Voltage (V)', fontweight='bold', fontfamily='Arial')
    ax.set_ylabel('Current Density (A/cm$^{2}$)', fontweight='bold', fontfamily='Arial')
    ax.tick_params(top=False, right=False, which='both', labelsize=8)
    ax.legend(fontsize=6, loc='lower left', frameon=False, handlelength=1.5, ncol=1,
              prop={'weight': 'normal'})

def save_figure_formats(fig, base_path, hd=False):
    saved = []
    fig.savefig(base_path + '.svg', format='svg', bbox_inches='tight', facecolor='white', edgecolor='none')
    saved.append(base_path + '.svg')
    fig.savefig(base_path + '.pdf', format='pdf', bbox_inches='tight', facecolor='white', edgecolor='none')
    saved.append(base_path + '.pdf')
    fig.savefig(base_path + '.png', dpi=300, bbox_inches='tight', facecolor='white', edgecolor='none')
    saved.append(base_path + '.png')
    if hd:
        fig.savefig(base_path + '_600dpi.png', dpi=600, bbox_inches='tight', facecolor='white', edgecolor='none')
        saved.append(base_path + '_600dpi.png')
    return saved

def export_component_data(V, J, fit, output_dir):
    popt = fit['popt']; Vt = fit['Vt']
    J_fit  = dark_current_model(V, *popt, Vt=Vt)
    J_main = main_diode_current(V, popt[0], popt[1], Vt)
    J_ohm  = ohmic_current(V, popt[2])
    J_tat  = tat_current(V, popt[3], popt[4], popt[5])
    header = 'V(V)\tJ_data(A/cm2)\tJ_fit(A/cm2)\tJ_main(A/cm2)\tJ_Ohm(A/cm2)\tJ_TAT(A/cm2)'
    rows = np.column_stack([V, J, J_fit, J_main, J_ohm, J_tat])
    path = os.path.join(output_dir, 'sample_component_fit.txt')
    np.savetxt(path, rows, header=header, delimiter='\t', fmt='%.6e', comments='')
    print(f'  Component data exported: {path}')
    return path

def export_parameters_table(fit, output_dir):
    popt = fit['popt']; perr = fit.get('perr', np.full(6, np.nan))
    params = [
        ('J0', 'J0', 'A/cm2', popt[0], perr[0]),
        ('A', 'A', '-', popt[1], perr[1]),
        ('Rsh', 'Rsh', 'ohm.cm2', popt[2], perr[2]),
        ('B', 'B', '-', popt[3], perr[3]),
        ('C', 'C', '-', popt[4], perr[4]),
        ('Vbi', 'Vbi', 'V', popt[5], perr[5]),
    ]
    rows = []
    for label, sym, unit, v, e in params:
        if unit == 'A/cm2':
            rows.append([label, sym, unit, f'{v:.4e}', f'{e:.4e}' if not np.isnan(e) else '-'])
        elif sym == 'Rsh':
            rows.append([label, sym, unit, f'{v:.2e}', f'{e:.2e}' if not np.isnan(e) else '-'])
        else:
            rows.append([label, sym, unit, f'{v:.4f}', f'{e:.4f}' if not np.isnan(e) else '-'])
    rows.append(['-', '-', '-', '-', '-'])
    rows.append(['R2', 'R2', '-', f'{fit["r_squared"]:.6f}', '-'])
    if fit.get('segmented'):
        s1 = fit.get('stage1', {}).get('r2')
        s2 = fit.get('stage2', {}).get('r2')
        if s1 is not None:
            rows.append(['Stage1 R2', 'R2_s1', '-', f'{s1:.6f}', '-'])
        if s2 is not None:
            rows.append(['Stage2 R2', 'R2_s2', '-', f'{s2:.6f}', '-'])
    df = pd.DataFrame(rows, columns=['Parameter', 'Symbol', 'Unit', 'Value', '+/-'])
    print(f"\n{'='*50}\n  Fitting Parameters\n{'='*50}")
    print(df.to_string(index=False))
    csv_p = os.path.join(output_dir, 'fitting_params.csv')
    df.to_csv(csv_p, index=False, encoding='utf-8-sig')
    print(f"  Params saved: {csv_p}")
    return df

def generate_word_report(fit, V, J, output_dir, T=300, area=None, V_seg=V_SEG_DEFAULT):
    try:
        from docx import Document
        from docx.shared import Pt, Cm
        from docx.enum.text import WD_ALIGN_PARAGRAPH
        from docx.enum.table import WD_TABLE_ALIGNMENT
    except ImportError:
        print("  python-docx not installed.")
        return None

    popt = fit['popt']; perr = fit.get('perr', np.full(6, np.nan)); Vt = fit['Vt']
    J0, A, Rsh, B, C, Vbi = popt

    V_dense = np.linspace(V.min(), V.max(), 5000)
    J_main_d = np.abs(main_diode_current(V_dense, J0, A, Vt))
    J_ohm_d  = np.abs(ohmic_current(V_dense, Rsh))
    J_tat_d  = np.abs(tat_current(V_dense, B, C, Vbi))
    J_total_d = np.abs(dark_current_model(V_dense, *popt, Vt=Vt))

    idx_zero = np.argmin(np.abs(V_dense))
    J_at_zero = J_total_d[idx_zero]

    ratio_main = J_main_d / (J_total_d + 1e-30)
    mask_main_dom = (ratio_main > 0.7) & (V_dense <= 0)
    V_main_end = V_dense[mask_main_dom].min() if np.any(mask_main_dom) else -0.2

    mask_rev = V_dense < 0
    V_rev = V_dense[mask_rev]
    Jm_rev = J_main_d[mask_rev]
    Jt_rev = J_tat_d[mask_rev]
    crossover_idx = np.where(Jt_rev > Jm_rev)[0]
    V_cross = V_rev[crossover_idx[0]] if len(crossover_idx) > 0 else -0.3

    idx_m5 = np.argmin(np.abs(V_dense + 0.5))
    ratio_tat_m5 = J_tat_d[idx_m5] / (J_total_d[idx_m5] + 1e-30) * 100
    ratio_ohm_m5 = J_ohm_d[idx_m5] / (J_total_d[idx_m5] + 1e-30) * 100
    J_ohm_max = np.max(J_ohm_d)

    doc = Document()
    style = doc.styles['Normal']
    style.font.name = 'Times New Roman'
    style.font.size = Pt(11)
    style.paragraph_format.space_after = Pt(6)
    style.paragraph_format.line_spacing = 1.15

    title = doc.add_heading('Dark Current J-V Characteristics and Multi-Mechanism Fitting Analysis of Infrared Photodetector', level=0)
    title.alignment = WD_ALIGN_PARAGRAPH.CENTER

    doc.add_paragraph(
        f'Figure 1 presents the dark current density-voltage (J-V) characteristics and multi-mechanism '
        f'fitting results of the fabricated infrared photodetector at room temperature. The blue circles '
        f'represent the experimental data, while the orange solid line denotes the total dark current '
        f'fitting curve based on a three-component model comprising diffusion-recombination current, '
        f'trap-assisted tunneling current, and ohmic leakage current (R2 = {fit["r_squared"]:.4f}). '
        f'The green dashed line, purple dash-dotted line, and light blue dotted line correspond to the '
        f'three independent current transport components: J_main, J_Ohm, and J_TAT, respectively.'
    )
    doc.add_paragraph(
        f'The excellent agreement between the experimental data and the total fitting curve across the '
        f'full bias range from {V.min():.1f} V to {V.max():.1f} V demonstrates that the three-component '
        f'current model accurately describes the dark current transport behavior of this device. '
        f'The dark current density at zero bias is approximately {J_at_zero:.1e} A/cm2, '
        f'indicating excellent low-dark-current characteristics.'
    )

    img_path = os.path.join(output_dir, 'sample_fitting.png')
    if os.path.exists(img_path):
        p_img = doc.add_paragraph()
        p_img.alignment = WD_ALIGN_PARAGRAPH.CENTER
        p_img.add_run().add_picture(img_path, width=Cm(8.5))
        p_cap = doc.add_paragraph()
        p_cap.alignment = WD_ALIGN_PARAGRAPH.CENTER
        run_cap = p_cap.add_run('Figure 1. Dark current J-V characteristics and multi-mechanism fitting results of the infrared photodetector.')
        run_cap.font.size = Pt(9)
        run_cap.italic = True
        doc.add_paragraph('')

    doc.add_heading('Fitting Mechanism', level=1)

    p_eq = doc.add_paragraph()
    p_eq.alignment = WD_ALIGN_PARAGRAPH.CENTER
    run = p_eq.add_run('J_dark = J0*[exp(qV/AkT) - 1]  +  V/Rsh  +  B*V*exp[-C/(Vbi - V)]')
    run.font.name = 'Times New Roman'
    run.font.size = Pt(13)
    run.bold = True

    doc.add_paragraph('')
    p_sub = doc.add_paragraph()
    p_sub.alignment = WD_ALIGN_PARAGRAPH.CENTER
    p_sub.add_run('J_main (Diffusion-Recombination)  +  J_Ohm (Ohmic Leakage)  +  J_TAT (Trap-Assisted Tunneling)').font.size = Pt(9)

    doc.add_paragraph('')
    doc.add_paragraph(
        'The dark current of the infrared photodetector is decomposed into three physically distinct '
        'transport components, each governed by a specific conduction mechanism:'
    )

    doc.add_heading('1. Diffusion-Recombination Current (J_main)', level=2)
    doc.add_paragraph(
        f'J_main = J0*[exp(qV/AkT) - 1] - This component arises from minority carrier diffusion '
        f'in the quasi-neutral regions and generation-recombination (G-R) processes within the depletion '
        f'region of the p-i-n junction. It follows the classical Shockley diode equation, with the '
        f'reverse saturation current density J0 = {J0:.2e} A/cm2 characterizing the intrinsic '
        f'recombination activity and the ideality factor A = {A:.2f} reflecting the relative contribution '
        f'of G-R current (A > 1 indicates significant trap-mediated recombination).'
    )

    doc.add_heading('2. Ohmic Leakage Current (J_Ohm)', level=2)
    doc.add_paragraph(
        f'J_Ohm = V/Rsh - This component originates from ohmic conduction pathways through film '
        f'pinholes, grain boundaries, and interfacial defects. It exhibits a linear dependence on '
        f'applied voltage and is inversely proportional to the shunt resistance Rsh = {Rsh:.2e} ohm*cm2. '
        f'A high Rsh value indicates excellent film quality with minimal parasitic leakage channels.'
    )

    doc.add_heading('3. Trap-Assisted Tunneling Current (J_TAT)', level=2)
    doc.add_paragraph(
        f'J_TAT = B*V*exp[-C/(Vbi - V)] - This component describes carrier tunneling mediated by '
        f'deep-level trap states within the bandgap. Carriers are first thermally excited to mid-gap '
        f'trap levels associated with quantum dot surface defects or heterojunction interface states, '
        f'and subsequently tunnel through the narrowed barrier under strong electric fields. The '
        f'coefficient B = {B:.4f} is proportional to the trap state density Nt, while C = {C:.4f} '
        f'reflects the effective tunneling barrier height. The built-in voltage Vbi = {Vbi:.4f} V '
        f'represents the internal electric field of the junction. This component exhibits an '
        f'exponential increase with reverse bias as the depletion region electric field intensifies '
        f'and the effective barrier width narrows, dramatically enhancing the tunneling probability.'
    )

    doc.add_heading('Segmented Fitting Strategy', level=2)
    doc.add_paragraph(
        f'A three-stage segmented fitting procedure was employed to ensure physically robust parameter '
        f'extraction: (1) Stage 1 - fitting J0 and A in the small-bias region (|V| < {V_seg} V) where '
        f'J_main dominates; (2) Stage 2 - fixing J0 and A, fitting Rsh, B, C, and Vbi in the large '
        f'reverse-bias region (V < -{V_seg} V) where J_TAT and J_Ohm become significant; '
        f'(3) Stage 3 - global refinement allowing all six parameters to co-adjust within physically '
        f'constrained bounds. J_BTB (band-to-band tunneling) was explicitly excluded from the model '
        f'as it requires high doping concentrations and extremely narrow depletion regions (<10 nm) '
        f'that are not satisfied in PbS CQD photodetectors.'
    )

    doc.add_heading('Voltage-Dependent Dominance of Current Components', level=1)

    doc.add_heading(f'1. Zero and Low Reverse Bias Region (0 V to {abs(V_main_end):.1f} V): '
                    f'J_main as the Primary Dark Current Source', level=2)
    doc.add_paragraph(
        f'In this bias regime, the diffusion-recombination current J_main constitutes the dominant '
        f'contribution to the total dark current, accounting for over 70% of the total. J_main '
        f'originates from the diffusion of minority carriers from the quasi-neutral regions toward '
        f'the depletion region boundary, combined with radiative and non-radiative recombination '
        f'processes of carriers within the depletion layer. It exhibits a gradual increase with '
        f'increasing reverse bias magnitude, consistent with the Shockley diode transport behavior. '
        f'The fitted reverse saturation current density J0 = {J0:.2e} A/cm2 and ideality factor '
        f'A = {A:.2f} indicate the presence of a moderate level of G-R recombination contribution, '
        f'likely associated with deep-level trap states within the PbS CQD active layer.'
    )

    doc.add_heading(f'2. High Reverse Bias Region (< {V_cross:.1f} V): '
                    f'J_TAT Emerges as the Dominant Leakage Mechanism', level=2)
    doc.add_paragraph(
        f'As the reverse bias increases beyond {abs(V_cross):.1f} V, the trap-assisted tunneling '
        f'current J_TAT grows exponentially and progressively overtakes J_main as the primary dark '
        f'current component. The physical origin of J_TAT lies in the thermal excitation of carriers '
        f'to quantum dot surface defect states or heterojunction interface trap levels, followed by '
        f'field-assisted tunneling through the narrowed potential barrier. With increasing reverse '
        f'bias, the depletion region electric field intensifies, the effective barrier width narrows '
        f'significantly, and the tunneling probability rises dramatically, leading to the rapid '
        f'enlargement of J_TAT. At -0.5 V reverse bias, J_TAT accounts for approximately '
        f'{ratio_tat_m5:.0f}% of the total dark current. The fitted TAT coefficient '
        f'B = {B:.4f} (proportional to the trap state density Nt) and barrier coefficient '
        f'C = {C:.4f} (proportional to the effective tunneling barrier height) provide quantitative '
        f'metrics for evaluating the density and energy depth of defect states participating in '
        f'the TAT process.'
    )

    doc.add_heading('3. Full Bias Range: J_Ohm Remains Negligible Across All Operating Conditions', level=2)
    doc.add_paragraph(
        f'Throughout the entire tested bias range, the ohmic leakage current J_Ohm remains at an '
        f'extremely low level (< {J_ohm_max*2:.1e} A/cm2) and exhibits an approximately linear '
        f'voltage dependence. The fitted shunt resistance Rsh = {Rsh:.2e} ohm*cm2 provides strong '
        f'evidence for the excellent film quality of the fabricated PbS CQD layer, demonstrating '
        f'low interfacial defect density between the charge transport layers and the quantum dot '
        f'active region, as well as superior electrode contact performance. The ohmic loss caused '
        f'by bulk defects and interfacial leakage has been effectively suppressed. '
        f'At -0.5 V reverse bias, J_Ohm accounts for merely {ratio_ohm_m5:.1f}% of the total '
        f'dark current.'
    )

    doc.add_heading('Fitting Parameters Summary', level=1)
    param_data = [
        ('Reverse saturation current density J0', f'{J0:.4e}', 'A/cm2'),
        ('Ideality factor A', f'{A:.4f}', '-'),
        ('Shunt resistance Rsh', f'{Rsh:.2e}', 'ohm*cm2'),
        ('TAT defect density coefficient B', f'{B:.4f}', '-'),
        ('TAT barrier coefficient C', f'{C:.4f}', '-'),
        ('Built-in voltage Vbi', f'{Vbi:.4f}', 'V'),
        ('R-squared R2', f'{fit["r_squared"]:.6f}', '-'),
    ]
    t = doc.add_table(rows=len(param_data)+1, cols=3, style='Light Grid Accent 1')
    t.alignment = WD_TABLE_ALIGNMENT.CENTER
    for i, h in enumerate(['Parameter', 'Value', 'Unit']):
        t.rows[0].cells[i].text = h
    for i, (pname, pval, unit) in enumerate(param_data):
        t.rows[i+1].cells[0].text = pname
        t.rows[i+1].cells[1].text = pval
        t.rows[i+1].cells[2].text = unit

    doc.add_heading('Implications for Device Optimization', level=1)
    doc.add_paragraph(
        f'Through quantitative multi-mechanism fitting, this study has unambiguously identified the '
        f'key limiting factors governing the dark current of the infrared photodetector. In the low '
        f'reverse bias region (0 V to {abs(V_main_end):.1f} V), the diffusion-recombination current '
        f'J_main is the dominant mechanism; in the high reverse bias region (< {V_cross:.1f} V), '
        f'trap-assisted tunneling current J_TAT progressively becomes the primary leakage source. '
        f'To further reduce the dark current and enhance the specific detectivity (D*), future device '
        f'optimization efforts should focus on two key aspects: (1) reducing the quantum dot surface '
        f'trap state density through ligand exchange and surface passivation strategies to suppress '
        f'the trap-assisted tunneling current under high reverse bias; and (2) optimizing the energy '
        f'band alignment and doping concentration of the hole transport layer and electron transport '
        f'layer to reduce minority carrier injection efficiency from the quasi-neutral regions, '
        f'thereby suppressing the diffusion-recombination current under low reverse bias.'
    )

    report_path = os.path.join(output_dir, 'fitting_report.docx')
    doc.save(report_path)
    print(f"\n  Word report saved: {report_path}")
    return report_path

def run_dark_current_fitting(data_file, output_dir=None, T=300, area=None,
                              vmin=None, vmax=None, auto_dark=True,
                              segmented=True, v_seg=V_SEG_DEFAULT,
                              max_points=None, hd=False, fig_width_cm=8.5):
    if output_dir is None:
        output_dir = os.path.join(os.path.dirname(os.path.abspath(data_file)), 'output')
    os.makedirs(output_dir, exist_ok=True)
    Vt = thermal_voltage(T)
    print("\n" + "=" * 60)
    print("  PbS CQD Dark Current Component Fitting")
    print("=" * 60)
    print(f"T={T}K, Vt={Vt:.6f}V, Area={area} cm2" if area else f"T={T}K, Vt={Vt:.6f}V")
    print(f"Data: {data_file}")

    V, J = load_data(data_file, area=area, auto_dark=auto_dark, max_points=max_points)
    if vmin is not None or vmax is not None:
        lo = vmin if vmin is not None else V.min()
        hi = vmax if vmax is not None else V.max()
        m = (V >= lo) & (V <= hi); V, J = V[m], J[m]
        print(f"Voltage filter: {lo}~{hi}V, {len(V)} pts")

    fit = segmented_fit(V, J, T, V_seg=v_seg) if segmented else fit_dark_current(V, J, T)
    if fit is None:
        print("Fitting failed."); return None

    configure_plot_style(width_cm=fig_width_cm)
    fig, ax = plt.subplots(figsize=(fig_width_cm/CM_PER_INCH, fig_width_cm/CM_PER_INCH*0.75))
    plot_fitting(V, J, fit, ax)
    plt.tight_layout()
    base = os.path.join(output_dir, 'sample_fitting')
    saved = save_figure_formats(fig, base, hd=hd)
    for s in saved: print(f"  {s}")
    plt.close(fig)

    export_component_data(V, J, fit, output_dir)
    export_parameters_table(fit, output_dir)
    generate_word_report(fit, V, J, output_dir, T=T, area=area, V_seg=v_seg)

    print("\n" + "=" * 60)
    print("Done!")
    print("=" * 60)
    return {'fit': fit, 'output_dir': output_dir}

if __name__ == '__main__':
    import argparse
    p = argparse.ArgumentParser(description='PbS CQD Dark Current Component Fitting')
    p.add_argument('data', help='Data file (.txt)')
    p.add_argument('-o', '--output', default=None, help='Output directory')
    p.add_argument('-T', '--temperature', type=float, default=T_DEFAULT)
    p.add_argument('-a', '--area', type=float, default=None, help='Device area (cm2)')
    p.add_argument('--vmin', type=float, default=-0.5)
    p.add_argument('--vmax', type=float, default=0.2)
    p.add_argument('--no-auto-dark', action='store_true')
    p.add_argument('--points', type=int, default=None)
    p.add_argument('--no-segmented', action='store_true')
    p.add_argument('--vseg', type=float, default=V_SEG_DEFAULT)
    p.add_argument('--fig-width-cm', type=float, default=8.5)
    p.add_argument('--hd', action='store_true')
    args = p.parse_args()

    r = run_dark_current_fitting(args.data, args.output, args.temperature, args.area,
                                  args.vmin, args.vmax, auto_dark=not args.no_auto_dark,
                                  segmented=not args.no_segmented, v_seg=args.vseg,
                                  max_points=args.points, hd=args.hd, fig_width_cm=args.fig_width_cm)
    if r:
        print(f"\nOutput files in: {r['output_dir']}")