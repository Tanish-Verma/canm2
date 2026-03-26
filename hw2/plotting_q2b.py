import numpy as np
import matplotlib.pyplot as plt
import glob
import os

BASE_DIR = os.path.dirname(os.path.abspath(__file__))
OUTPUT_DIR = os.path.join(BASE_DIR, "output2")
PLOTS_DIR = os.path.join(BASE_DIR, "plots")
os.makedirs(PLOTS_DIR, exist_ok=True)

# -------------------------------------------------------
# Load all output snapshots from a folder
# Each file: two columns  x  T
# -------------------------------------------------------
def load_snapshots(folder):
    files = sorted(glob.glob(os.path.join(folder, "output_*.dat")))
    if len(files) == 0:
        raise FileNotFoundError(f"No output_*.dat files found in '{folder}'")
    snapshots = []
    for f in files:
        data = np.loadtxt(f)
        snapshots.append((data[:, 0], data[:, 1]))
    return snapshots

# -------------------------------------------------------
# Cases: (folder, r value, stable/unstable label)
# -------------------------------------------------------
cases = [
    (os.path.join(OUTPUT_DIR, "r040"), 0.40, "stable"),
    (os.path.join(OUTPUT_DIR, "r049"), 0.49, "stable"),
    (os.path.join(OUTPUT_DIR, "r051"), 0.51, "unstable"),
    (os.path.join(OUTPUT_DIR, "r070"), 0.70, "unstable"),
]

fig, axes = plt.subplots(2, 2, figsize=(12, 8))
axes = axes.flatten()

for ax, (folder, r_val, stability) in zip(axes, cases):
    snapshots = load_snapshots(folder)
    n = len(snapshots)
    colors = plt.cm.viridis(np.linspace(0.1, 0.9, n))

    for i, (x, T) in enumerate(snapshots):
        label = f't-step {i+1}' if i == 0 or i == n-1 else None
        ax.plot(x, T, color=colors[i], linewidth=1.6, label=label)

    # colorbar-style time indicator using a dummy ScalarMappable
    sm = plt.cm.ScalarMappable(cmap='viridis',
                                norm=plt.Normalize(vmin=1, vmax=n))
    sm.set_array([])
    cbar = plt.colorbar(sm, ax=ax)
    cbar.set_label('time step index', fontsize=10)

    status = "Stable" if stability == "stable" else "Unstable"
    ax.set_title(f'{status}:  r = {r_val}', fontsize=13)
    ax.set_xlabel('$x$', fontsize=12)
    ax.set_ylabel('$T$', fontsize=12)
    ax.grid(True, alpha=0.3)
    ax.tick_params(labelsize=11)

plt.suptitle('Q2(b): Forward Euler — stable vs unstable behaviour', fontsize=14)
plt.tight_layout()
outpath = os.path.join(PLOTS_DIR, 'q2b_stability.png')
plt.savefig(outpath, dpi=150, bbox_inches='tight')
plt.show()
print(f"Saved {outpath}")