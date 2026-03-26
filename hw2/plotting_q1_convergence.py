import numpy as np
import matplotlib.pyplot as plt
import glob
import os

os.makedirs("plots", exist_ok=True)

for f in sorted(glob.glob("output1/error_b*.dat")):
    data = np.loadtxt(f)
    if data.ndim == 1:
        data = data.reshape(1, -1)
    nx, err = data[:, 0], data[:, 2]
    
    beta_tag = f.split("_b")[1].split(".")[0]
    beta_val = float(beta_tag) / 100
    
    nx_ref = np.array([nx[0], nx[-1]])
    err_ref = err[0] * (nx[0] / nx_ref) ** 2
    
    fig, ax = plt.subplots(figsize=(7, 5))
    ax.loglog(nx, err, "o-", lw=2, ms=7, label=rf"$\beta = {beta_val:g}$")
    ax.loglog(nx_ref, err_ref, "k:", lw=1.5, label="slope $-2$")
    ax.set_xlabel("$n_x$", fontsize=14)
    ax.set_ylabel("$L_2$ error", fontsize=14)
    ax.set_title(rf"Order of accuracy - stretched grid, $\beta = {beta_val:g}$ (Q1)", fontsize=13)
    ax.legend(fontsize=12)
    ax.grid(True, which="both", alpha=0.3)
    ax.tick_params(labelsize=12)
    
    fig.tight_layout()
    outfile = f"plots/q1_convergence_b{beta_tag}.png"
    fig.savefig(outfile, dpi=150)
    plt.close(fig)
    print(f"Saved {outfile}")