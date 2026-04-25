import numpy as np
import matplotlib.pyplot as plt
from glob import glob
import os

# ─── Parameters ───────────────────────────────────────────────────────────────
kappa     = 0.010   # Roll 23
nu        = 0.80    # Roll 23
c         = 1.0
Nx        = 100
dx        = 1.0 / Nx
dt        = nu * dx / c
r         = kappa * dt / (dx**2)
T         = 1.0
target_times  = [0.0, 0.25, 0.5, 0.75, 1.0]
plots_dir = "Q4/plots"
data_dir  = "Q4/data"
os.makedirs(plots_dir, exist_ok=True)

print(f"kappa={kappa}, nu={nu}, r={r:.4f}, Pe_h={nu/r:.4f}")

# ─── Actual simulation times: n*dt ────────────────────────────────────────────
def getTimes(dt, target_times):
    n_steps = np.ceil(np.array(target_times, dtype=float) / dt).astype(int)
    return n_steps * dt

actual_times = getTimes(dt, target_times)
print(f"Actual output times (n*dt): {actual_times}")

# ─── Exact solution ───────────────────────────────────────────────────────────
def exactSolution(x, t, kappa):
    return np.exp(-4.0 * np.pi**2 * kappa * t) * np.sin(2 * np.pi * (x - t))

# ─── Read all snapshots for a scheme ─────────────────────────────────────────
def readAllSnapshots(scheme, data_dir="data"):
    files = sorted(
        glob(os.path.join(data_dir, f"{scheme}_T_*.txt")),
        key=lambda p: float(os.path.basename(p).split("_T_")[1].replace(".txt", ""))
    )
    snapshots, times_found = [], []
    for f in files:
        t = float(os.path.basename(f).split("_T_")[1].replace(".txt", ""))
        data = np.loadtxt(f)
        snapshots.append((data[:, 0], data[:, 1]))
        times_found.append(t)
    return snapshots, times_found

# ─── L2 error ────────────────────────────────────────────────────────────────
def computeL2Error(u, x, t, kappa, dx):
    return np.sqrt(dx * np.sum((u - exactSolution(x, t, kappa))**2))

# ─── (b) 3-row x 5-col figure: one row per scheme, one col per output time ───
def plotAllSchemes(schemes, labels, actual_times, kappa, data_dir, plots_dir,
                   ylim=(-1.5, 1.5)):
    fig, axes = plt.subplots(3, 5, figsize=(20, 10), sharey=True)

    for row, (scheme, label) in enumerate(zip(schemes, labels)):
        snapshots, times_found = readAllSnapshots(scheme, data_dir)

        for col, t_target in enumerate(actual_times):
            ax = axes[row][col]

            idx = np.argmin([abs(tf - t_target) for tf in times_found])
            x, u = snapshots[idx]
            t_file = times_found[idx]

            # Clip blown-up values for FTCS_UW
            u_plot = np.clip(u, -10, 10)
            clipped = not np.allclose(u, u_plot)

            ax.plot(x, u_plot, 'b-', linewidth=1.8, label="Numerical")
            ax.plot(x, exactSolution(x, t_file, kappa), 'r--',
                    linewidth=1.8, label="Exact")

            ax.set_title(f"t = {t_file:.4f}s", fontsize=9)
            ax.set_xlim(0, 1)
            ax.set_ylim(ylim)
            ax.grid(True, alpha=0.3)

            if clipped:
                ax.text(0.5, 0.5, "DIVERGED", transform=ax.transAxes,
                        ha='center', va='center', fontsize=12,
                        color='red', fontweight='bold')

        # Row label on leftmost column
        axes[row][0].set_ylabel(f"{label}\nu(x,t)", fontsize=10, fontweight='bold')

    # Column titles on top row
    time_labels = ["t = 0", "t ≈ 0.25", "t ≈ 0.50", "t ≈ 0.75", "t ≈ 1.00"]
    for col, tlabel in enumerate(time_labels):
        axes[0][col].set_title(tlabel, fontsize=10, fontweight='bold')

    # x-axis label on bottom row
    for col in range(5):
        axes[2][col].set_xlabel("x", fontsize=10)

    # Legend on first subplot only
    axes[0][0].legend(fontsize=8, loc='upper right')

    fig.suptitle(f"1D Advection-Diffusion Schemes (Roll 23, κ={kappa}, ν={nu}, r={r:.2f}, Pe_h={nu/r:.2f})",
                 fontsize=13, fontweight='bold')
    plt.tight_layout()
    out_path = os.path.join(plots_dir, "Q4_AllSchemes.pdf")
    plt.savefig(out_path, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"Saved: {out_path}")

# ─── (c) L2 error table ───────────────────────────────────────────────────────
def printL2Table(schemes, labels, actual_times, kappa, dx, data_dir):
    print("\n=== L2 Error at final time ===")
    print(f"{'Scheme':<15} | {'L2 Error':>12}")
    print("-" * 32)
    for scheme, label in zip(schemes, labels):
        snapshots, times_found = readAllSnapshots(scheme, data_dir)
        if not snapshots:
            continue
        # Last snapshot
        x, u = snapshots[-1]
        t = times_found[-1]
        e = computeL2Error(u, x, t, kappa, dx)
        print(f"{label:<15} | {e:>12.6e}")

# ─── (d) CN with large nu: 3-row x 5-col figure ──────────────────────────────
def plotCNLargeNu(nu_vals, actual_times, kappa, dx, c, data_dir, plots_dir,
                  ylim=(-1.5, 1.5)):
    fig, axes = plt.subplots(3, 5, figsize=(20, 10), sharey=True)

    for row, nu_val in enumerate(nu_vals):
        dt_val = nu_val * dx / c
        r_val  = kappa * dt_val / (dx**2)
        actual_t = getTimes(dt_val, target_times)

        scheme = f"CN_nu{int(nu_val)}"
        snapshots, times_found = readAllSnapshots(scheme, data_dir)

        for col, t_target in enumerate(actual_t):
            ax = axes[row][col]

            idx = np.argmin([abs(tf - t_target) for tf in times_found])
            x, u = snapshots[idx]
            t_file = times_found[idx]

            ax.plot(x, u, 'b-', linewidth=1.8, label="CN")
            ax.plot(x, exactSolution(x, t_file, kappa), 'r--',
                    linewidth=1.8, label="Exact")
            ax.set_title(f"t = {t_file:.4f}s", fontsize=9)
            ax.set_xlim(0, 1)
            ax.set_ylim(ylim)
            ax.grid(True, alpha=0.3)

        axes[row][0].set_ylabel(f"ν={nu_val}, r={r_val:.1f}\nu(x,t)",
                                fontsize=10, fontweight='bold')

    for col in range(5):
        axes[2][col].set_xlabel("x", fontsize=10)
    axes[0][0].legend(fontsize=8, loc='upper right')

    fig.suptitle(f"Crank-Nicolson with Large ν (κ={kappa}) — Unconditional Stability ≠ Accuracy",
                 fontsize=13, fontweight='bold')
    plt.tight_layout()
    out_path = os.path.join(plots_dir, "Q4_CN_LargeNu.pdf")
    plt.savefig(out_path, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"Saved: {out_path}")

# ─── Main ─────────────────────────────────────────────────────────────────────
schemes = ["FTCS_UW", "BTCS_UW", "CN"]
labels  = ["FTCS-UW (Expl.)", "BTCS-UW (Impl.)", "Crank-Nicolson"]

plotAllSchemes(schemes, labels, actual_times, kappa, data_dir, plots_dir)
printL2Table(schemes, labels, actual_times, kappa, dx, data_dir)

nu_large = [2.0, 5.0, 10.0]
plotCNLargeNu(nu_large, actual_times, kappa, dx, c, data_dir, plots_dir)

print("\nProgram completed!")