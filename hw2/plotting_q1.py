import numpy as np
import matplotlib.pyplot as plt
import glob
import os

os.makedirs("plots", exist_ok=True)

LABEL_FS = 16
TITLE_FS = 15
TICK_FS = 13
LEGEND_FS = 13
LINE_LW = 2.6
REF_LW = 2.2

for f in sorted(glob.glob("output1/error_b*.dat")):
    data = np.loadtxt(f)
    if data.ndim == 1:
        data = data.reshape(1, -1)
    nx, err = data[:, 0], data[:, 2]
    
    beta_tag = f.split("_b")[1].split(".")[0]
    beta_val = float(beta_tag) / 100
    
    nx_ref = np.array([nx[0], nx[-1]])
    err_ref = err[0] * (nx[0] / nx_ref) ** 2
    
    fig, ax = plt.subplots(figsize=(8, 6))
    ax.loglog(nx, err, "o-", lw=LINE_LW, ms=8, label=rf"$\beta = {beta_val:g}$")
    ax.loglog(nx_ref, err_ref, "k:", lw=REF_LW, label="slope $-2$")
    ax.set_xlabel("$n_x$", fontsize=LABEL_FS)
    ax.set_ylabel("$L_2$ error", fontsize=LABEL_FS)
    ax.set_title(rf"Order of accuracy - stretched grid, $\beta = {beta_val:g}$ (Q1)", fontsize=TITLE_FS)
    ax.legend(fontsize=LEGEND_FS)
    ax.grid(True, which="both", alpha=0.3)
    ax.tick_params(labelsize=TICK_FS)
    
    fig.tight_layout()
    outfile = f"plots/q1_convergence_b{beta_tag}.png"
    fig.savefig(outfile, dpi=150)
    plt.close(fig)
    print(f"Saved {outfile}")