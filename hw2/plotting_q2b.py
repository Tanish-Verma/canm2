import numpy as np
import matplotlib.pyplot as plt
import glob
import os

os.makedirs("plots", exist_ok=True)

LABEL_FS = 15
TITLE_FS = 15
TICK_FS = 12
C_BAR_FS = 12
LINE_LW = 2.3

fig, axes = plt.subplots(2, 2, figsize=(14, 10))

for ax, (folder, r, status) in zip(axes.flatten(), [
    ("output2/r040", 0.40, "Stable"),
    ("output2/r049", 0.49, "Stable"),
    ("output2/r051", 0.51, "Unstable"),
    ("output2/r070", 0.70, "Unstable"),
]):
    files = sorted(glob.glob(f"{folder}/output_*.dat"))
    snaps = [np.loadtxt(f)[:, :2].T for f in files]
    n = len(snaps)
    colors = plt.cm.viridis(np.linspace(0.1, 0.9, n))
    
    for i, (x, t) in enumerate(snaps):
        ax.plot(x, t, color=colors[i], lw=LINE_LW, label=f"t-step {i+1}" if i in (0, n-1) else None)
    
    sm = plt.cm.ScalarMappable(cmap="viridis", norm=plt.Normalize(vmin=1, vmax=n))
    sm.set_array([])
    cbar = plt.colorbar(sm, ax=ax)
    cbar.set_label("time step index", fontsize=C_BAR_FS)
    cbar.ax.tick_params(labelsize=TICK_FS)
    
    ax.set_title(f"{status}: r = {r}", fontsize=TITLE_FS)
    ax.set_xlabel("$x$", fontsize=LABEL_FS)
    ax.set_ylabel("$T$", fontsize=LABEL_FS)
    ax.grid(True, alpha=0.3)
    ax.tick_params(labelsize=TICK_FS)
    ax.legend(fontsize=12)

fig.suptitle("Q2(b): Forward Euler - stable vs unstable behaviour", fontsize=17)
fig.tight_layout()
outpath = "plots/q2b_stability.png"
fig.savefig(outpath, dpi=150, bbox_inches="tight")
plt.close(fig)
print(f"Saved {outpath}")