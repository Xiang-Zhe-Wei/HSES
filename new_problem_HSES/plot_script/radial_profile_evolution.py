#!/usr/bin/env python
# radial_profile_evolution.py

import os
import numpy as np
import matplotlib.pyplot as plt
import yt
import argparse

# ---------- CLI ----------
ap = argparse.ArgumentParser(description='Radial density evolution vs analytic Hernquist')
ap.add_argument('-s', type=int, default=None, help='first data index (auto if None)')
ap.add_argument('-e', type=int, default=None, help='last data index (auto if None)')
ap.add_argument('-d', type=int, default=1, help='delta index  [%(default)d]')
ap.add_argument('-i', type=str, default='..', help='data path prefix [%(default)s]')
ap.add_argument('-o', type=str, default='rho_evo.png', help='output figure name')
args = ap.parse_args()

# ---------- Hernquist density ----------
def rho_hernquist(r, M=1.0, a=1.0):
    return M * a / (2 * np.pi * r * (r + a)**3)

# ---------- Data Files ----------
data_dirs = sorted([d for d in os.listdir(args.i) if d.startswith('Data_')]) # ['Data_000000', 'Data_000001', ...]
if args.s is None:
    args.s = 0
if args.e is None:
    args.e = len(data_dirs) - 1

data_dirs = data_dirs[args.s:args.e+1:args.d]

# ---------- Plotting ----------
plt.figure(figsize=(8, 6))
for idx, dname in enumerate(data_dirs):
    ds = yt.load(os.path.join(args.i, dname))
    sp = ds.sphere("c", (0.5, "unitary"))
    prof = yt.ProfilePlot(sp, "radius", ["density"], weight_field="cell_volume")
    radius = prof.profiles[0].x.value
    density = prof.profiles[0]['density'].value

    # avoid 0 happens
    mask = density > 0
    radius = radius[mask]
    density = density[mask]

    label = f"t = {ds.current_time:.3e}"
    plt.plot(radius, density, label=label, alpha=0.7)

# Plot analytic Hernquist density
radius_fine = np.linspace(radius.min(), radius.max(), 1000)
analytic_density = rho_hernquist(radius_fine)
plt.plot(radius_fine, analytic_density, 'k--', lw=2, label='Analytic (Hernquist)')

plt.xlabel('Radius (code units)')
plt.ylabel('Density (code units)')
plt.xscale('log')
plt.yscale('log')
plt.legend(loc='upper right', fontsize='small', ncol=1)
plt.title('Radial Density Evolution vs Analytic Hernquist Solution')
plt.grid(True, which='both', linestyle='--', alpha=0.5)
plt.tight_layout()
plt.savefig(args.o, dpi=300)
plt.close()
print(f"Radial density evolution plot saved to {args.o}")

# python radial_profile_evolution.py -s 0 -e 10 -d 1  -i .. -o density_plot.png