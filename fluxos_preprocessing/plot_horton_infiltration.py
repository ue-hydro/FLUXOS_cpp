#!/usr/bin/env python3
"""
FLUXOS - Horton Infiltration Module Visualization
===================================================
Generates a multi-panel figure showing the Horton decay infiltration model:
  Panel 1: Spatial map of infiltration rate at 4 timesteps
  Panel 2: Time series of mean/max infiltration rate vs theoretical Horton curve
  Panel 3: Spatial map of water depth at same timesteps

Usage:
    python plot_horton_infiltration.py --results-dir Results_regular_soil/
"""

import os
import sys
import csv
import math
import argparse
import numpy as np

# Use Agg backend for headless rendering
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.colors import Normalize, LogNorm
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.ticker as mticker


def load_regular_mesh(filepath):
    """Load regular mesh output file into structured arrays."""
    rows, cols, x, y, h, infil_rate = [], [], [], [], [], []
    with open(filepath) as f:
        reader = csv.DictReader(f)
        for row in reader:
            rows.append(int(row['irow [-]']))
            cols.append(int(row['icol [-]']))
            x.append(float(row['Xcoord']))
            y.append(float(row['Ycoord']))
            h.append(float(row['h [m]']))
            infil_rate.append(float(row['soil_infil_rate [m/s]']))

    rows = np.array(rows)
    cols = np.array(cols)
    x = np.array(x)
    y = np.array(y)
    h = np.array(h)
    infil_rate = np.array(infil_rate)

    # Build 2D grids
    rmin, rmax = rows.min(), rows.max()
    cmin, cmax = cols.min(), cols.max()
    nr = rmax - rmin + 1
    nc = cmax - cmin + 1

    h_grid = np.full((nr, nc), np.nan)
    infil_grid = np.full((nr, nc), np.nan)
    x_grid = np.full((nr, nc), np.nan)
    y_grid = np.full((nr, nc), np.nan)

    for i in range(len(rows)):
        ri = rows[i] - rmin
        ci = cols[i] - cmin
        h_grid[ri, ci] = h[i]
        infil_grid[ri, ci] = infil_rate[i]
        x_grid[ri, ci] = x[i]
        y_grid[ri, ci] = y[i]

    return {
        'h': h_grid,
        'infil_rate': infil_grid,
        'x': x_grid,
        'y': y_grid,
        'xmin': np.nanmin(x), 'xmax': np.nanmax(x),
        'ymin': np.nanmin(y), 'ymax': np.nanmax(y),
    }


def horton_rate(t_sec, f0, fc, k):
    """Compute Horton infiltration rate: f(t) = fc + (f0 - fc) * exp(-k * t)"""
    return fc + (f0 - fc) * math.exp(-k * t_sec)


def main():
    parser = argparse.ArgumentParser(description='FLUXOS Horton Infiltration Visualization')
    parser.add_argument('--results-dir', default='Results_regular_soil/',
                        help='Path to regular mesh results directory')
    parser.add_argument('--output', default='horton_infiltration_results.png',
                        help='Output image filename')
    args = parser.parse_args()

    results_dir = args.results_dir

    # Sandy Loam Horton parameters (from USDA table)
    Ks = 1.089e-5   # fc - final rate [m/s]
    f0 = 5.445e-5   # initial rate [m/s]
    k  = 2.78e-4    # decay constant [1/s]

    # ---- Discover and sort output files ----
    files = sorted(
        [f for f in os.listdir(results_dir) if f.endswith('.txt') and f[0].isdigit()],
        key=lambda f: int(f.replace('.txt', ''))
    )

    if not files:
        print(f"No output files found in {results_dir}")
        sys.exit(1)

    print(f"Found {len(files)} output files in {results_dir}")

    # ---- Load time series data ----
    times_sec = []
    mean_rates = []
    max_rates = []
    min_rates = []
    wet_counts = []

    for fname in files:
        ts = int(fname.replace('.txt', ''))
        times_sec.append(ts)
        fpath = os.path.join(results_dir, fname)
        with open(fpath) as f:
            reader = csv.DictReader(f)
            rates = [float(row['soil_infil_rate [m/s]']) for row in reader]
        rates = np.array(rates)
        wet = rates[rates > 0]
        if len(wet) > 0:
            mean_rates.append(wet.mean())
            max_rates.append(wet.max())
            min_rates.append(wet.min())
        else:
            mean_rates.append(0)
            max_rates.append(0)
            min_rates.append(0)
        wet_counts.append(len(wet))

    times_min = np.array(times_sec) / 60.0
    mean_rates = np.array(mean_rates)
    max_rates = np.array(max_rates)
    min_rates = np.array(min_rates)

    # ---- Select 4 representative timesteps for spatial maps ----
    # Pick: early, 1/3, 2/3, late
    n = len(files)
    if n >= 4:
        idx_map = [1, n // 3, 2 * n // 3, n - 1]
    else:
        idx_map = list(range(n))

    map_data = []
    for idx in idx_map:
        fpath = os.path.join(results_dir, files[idx])
        data = load_regular_mesh(fpath)
        data['time_sec'] = times_sec[idx]
        data['time_min'] = times_sec[idx] / 60.0
        map_data.append(data)
        print(f"  Loaded {files[idx]} (t={data['time_min']:.0f} min)")

    # ---- Compute theoretical Horton curve ----
    t_theory = np.linspace(0, max(times_sec), 500)
    f_theory = np.array([horton_rate(t, f0, Ks, k) for t in t_theory])
    t_theory_min = t_theory / 60.0

    # ============================================================================
    # Create multi-panel figure
    # ============================================================================
    fig = plt.figure(figsize=(20, 16))
    fig.suptitle('FLUXOS Horton Soil Infiltration Module', fontsize=20, fontweight='bold', y=0.98)
    fig.text(0.5, 0.955, r'$f(t) = f_c + (f_0 - f_c) \cdot e^{-k \cdot t_{wet}}$'
             f'   |   Sandy Loam: $f_0$={f0*1e3:.3f} mm/s, $f_c$={Ks*1e3:.4f} mm/s, $k$={k:.2e} 1/s',
             ha='center', fontsize=13, style='italic', color='#444444')

    gs = gridspec.GridSpec(3, 4, hspace=0.35, wspace=0.3,
                           left=0.06, right=0.94, top=0.92, bottom=0.06)

    # ---- Row 1: Infiltration rate spatial maps (4 timesteps) ----
    infil_max = max(d['infil_rate'][np.isfinite(d['infil_rate'])].max() for d in map_data)
    infil_min_nz = Ks * 0.5  # clip below half Ks for visualization

    for i, data in enumerate(map_data):
        ax = fig.add_subplot(gs[0, i])
        grid = data['infil_rate'].copy()
        grid[grid <= 0] = np.nan  # mask dry cells

        # Convert to mm/hr for readability
        grid_mmhr = grid * 3600 * 1000  # m/s -> mm/hr

        im = ax.imshow(grid_mmhr, origin='lower',
                       extent=[data['xmin'], data['xmax'], data['ymin'], data['ymax']],
                       cmap='YlOrRd_r',
                       vmin=Ks * 3600 * 1000 * 0.9,
                       vmax=f0 * 3600 * 1000 * 0.75,
                       aspect='equal')

        t_min = data['time_min']
        f_t = horton_rate(data['time_sec'], f0, Ks, k)
        ax.set_title(f't = {t_min:.0f} min\n$f(t)$={f_t*3600*1000:.1f} mm/hr',
                     fontsize=11, fontweight='bold')

        ax.tick_params(labelsize=7)
        ax.xaxis.set_major_formatter(mticker.FuncFormatter(lambda x, p: f'{x/1000:.1f}'))
        ax.yaxis.set_major_formatter(mticker.FuncFormatter(lambda x, p: f'{x/1e6:.4f}'))

        if i == 0:
            ax.set_ylabel('Y [×10⁶ m]', fontsize=9)

    # Colorbar for infiltration rate
    cbar_ax1 = fig.add_axes([0.95, 0.66, 0.012, 0.26])
    cbar1 = fig.colorbar(im, cax=cbar_ax1)
    cbar1.set_label('Infiltration Rate [mm/hr]', fontsize=10)
    cbar1.ax.tick_params(labelsize=8)

    # Row 1 label
    fig.text(0.02, 0.80, 'Infiltration\nRate', fontsize=12, fontweight='bold',
             ha='center', va='center', rotation=90, color='#B22222')

    # ---- Row 2: Water depth spatial maps (same 4 timesteps) ----
    for i, data in enumerate(map_data):
        ax = fig.add_subplot(gs[1, i])
        grid = data['h'].copy()
        grid[grid <= 1e-6] = np.nan  # mask dry cells

        # Convert to cm
        grid_cm = grid * 100

        im2 = ax.imshow(grid_cm, origin='lower',
                        extent=[data['xmin'], data['xmax'], data['ymin'], data['ymax']],
                        cmap='Blues',
                        vmin=0, vmax=np.nanpercentile(
                            np.concatenate([d['h'][d['h'] > 1e-6] * 100 for d in map_data]),
                            95
                        ),
                        aspect='equal')

        ax.set_title(f't = {data["time_min"]:.0f} min', fontsize=11, fontweight='bold')
        ax.tick_params(labelsize=7)
        ax.xaxis.set_major_formatter(mticker.FuncFormatter(lambda x, p: f'{x/1000:.1f}'))
        ax.yaxis.set_major_formatter(mticker.FuncFormatter(lambda x, p: f'{x/1e6:.4f}'))

        if i == 0:
            ax.set_ylabel('Y [×10⁶ m]', fontsize=9)

    # Colorbar for water depth
    cbar_ax2 = fig.add_axes([0.95, 0.37, 0.012, 0.26])
    cbar2 = fig.colorbar(im2, cax=cbar_ax2)
    cbar2.set_label('Water Depth [cm]', fontsize=10)
    cbar2.ax.tick_params(labelsize=8)

    # Row 2 label
    fig.text(0.02, 0.51, 'Water\nDepth', fontsize=12, fontweight='bold',
             ha='center', va='center', rotation=90, color='#1E3A8A')

    # ---- Row 3: Time series plots ----
    # Left half: Infiltration rate time series
    ax_ts = fig.add_subplot(gs[2, 0:2])

    # Convert to mm/hr
    ax_ts.fill_between(times_min, min_rates * 3600 * 1000, max_rates * 3600 * 1000,
                       alpha=0.15, color='steelblue', label='Min–Max range')
    ax_ts.plot(times_min, mean_rates * 3600 * 1000, 'o-', color='steelblue',
               linewidth=2, markersize=5, label='Domain mean (simulated)')
    ax_ts.plot(times_min, max_rates * 3600 * 1000, '^--', color='royalblue',
               linewidth=1, markersize=4, alpha=0.6, label='Domain max')
    ax_ts.plot(t_theory_min, f_theory * 3600 * 1000, '-', color='red',
               linewidth=2.5, label=r'Horton theory: $f_c + (f_0 - f_c)e^{-kt}$')
    ax_ts.axhline(y=Ks * 3600 * 1000, color='darkred', linestyle=':', linewidth=1.5,
                  label=f'$f_c$ = Ks = {Ks*3600*1000:.2f} mm/hr')
    ax_ts.axhline(y=f0 * 3600 * 1000, color='orange', linestyle=':', linewidth=1.5,
                  label=f'$f_0$ = {f0*3600*1000:.1f} mm/hr')

    # Mark the 4 spatial map timesteps
    for data in map_data:
        ax_ts.axvline(x=data['time_min'], color='gray', linestyle='--', alpha=0.3, linewidth=0.8)

    ax_ts.set_xlabel('Time [min]', fontsize=12)
    ax_ts.set_ylabel('Infiltration Rate [mm/hr]', fontsize=12)
    ax_ts.set_title('Infiltration Rate Decay Over Time', fontsize=13, fontweight='bold')
    ax_ts.legend(fontsize=8.5, loc='upper right', framealpha=0.9)
    ax_ts.set_xlim(0, max(times_min) * 1.02)
    ax_ts.set_ylim(0, f0 * 3600 * 1000 * 1.15)
    ax_ts.grid(True, alpha=0.3)

    # Right half: Wet cell count + water volume
    ax_wet = fig.add_subplot(gs[2, 2:4])

    color_wet = '#2E86AB'
    color_ratio = '#A23B72'

    ax_wet.bar(times_min, wet_counts, width=times_min[1] - times_min[0] if len(times_min) > 1 else 30,
               color=color_wet, alpha=0.6, label='Wet cells (infiltrating)')
    ax_wet.set_xlabel('Time [min]', fontsize=12)
    ax_wet.set_ylabel('Number of Infiltrating Cells', fontsize=12, color=color_wet)
    ax_wet.tick_params(axis='y', labelcolor=color_wet)
    ax_wet.set_title('Wetting Front Progression & Decay Ratio', fontsize=13, fontweight='bold')
    ax_wet.set_xlim(0, max(times_min) * 1.02)
    ax_wet.grid(True, alpha=0.3)

    # Secondary axis: ratio of mean rate to Ks
    ax_ratio = ax_wet.twinx()
    ratio = mean_rates / Ks
    ratio[ratio <= 0] = np.nan
    ax_ratio.plot(times_min, ratio, 's-', color=color_ratio, linewidth=2, markersize=5,
                  label='Mean rate / Ks ratio')
    ax_ratio.axhline(y=1.0, color='darkred', linestyle=':', linewidth=1.5, alpha=0.7)
    ax_ratio.set_ylabel('Mean Infil Rate / Ks', fontsize=12, color=color_ratio)
    ax_ratio.tick_params(axis='y', labelcolor=color_ratio)
    ax_ratio.set_ylim(0, max(ratio[np.isfinite(ratio)]) * 1.2 if any(np.isfinite(ratio)) else 5)

    # Combined legend
    lines1, labels1 = ax_wet.get_legend_handles_labels()
    lines2, labels2 = ax_ratio.get_legend_handles_labels()
    ax_wet.legend(lines1 + lines2, labels1 + labels2, fontsize=9, loc='upper left', framealpha=0.9)

    # ---- Save ----
    outpath = args.output
    fig.savefig(outpath, dpi=180, bbox_inches='tight', facecolor='white')
    print(f"\nVisualization saved to: {outpath}")
    print(f"  Image size: {os.path.getsize(outpath) / 1e6:.1f} MB")
    plt.close()


if __name__ == '__main__':
    main()
