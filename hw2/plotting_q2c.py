import numpy as np
import matplotlib.pyplot as plt
import os

os.makedirs("plots", exist_ok=True)

LABEL_FS = 16
TITLE_FS = 15
TICK_FS = 13
LEGEND_FS = 13
LINE_LW = 2.6
REF_LW = 2.2

def main():
    data = np.loadtxt("output2/error2.dat")
    nx, err = (data.reshape(1, -1) if data.ndim == 1 else data)[:, :2].T

    # Fit err = C * nx^-2 with fixed slope -2 
    c = np.exp(np.mean(np.log(err) + 2 * np.log(nx)))
    err_m2 = c * nx ** -2

    sat = np.where(np.diff(err) >= 0)[0]
    sat_idx = sat[0] + 1 if sat.size else None

    fig, ax = plt.subplots(figsize=(9, 6))
    ax.loglog(nx, err, "o-", lw=LINE_LW, ms=8, color="steelblue", label="$L_2$ error")
    ax.loglog(nx, err_m2, "k:", lw=REF_LW, label="slope $-2$")

    if sat_idx is not None:
        ax.loglog(nx[sat_idx], err[sat_idx], "r^", ms=11, label=f"saturation begins (nx={int(nx[sat_idx])})")

    ax.set_xlabel("$n_x$", fontsize=LABEL_FS)
    ax.set_ylabel("$L_2$ error", fontsize=LABEL_FS)
    ax.set_title("Q2(c): Convergence - forward Euler, 2nd-order central difference", fontsize=TITLE_FS)
    ax.legend(fontsize=LEGEND_FS)
    ax.grid(True, which="both", alpha=0.3)
    ax.tick_params(labelsize=TICK_FS)

    fig.tight_layout()
    out_file = "plots/q2c_convergence.png"
    fig.savefig(out_file, dpi=150)
    plt.close(fig)
    print(f"Saved {out_file}")


if __name__ == "__main__":
    main()