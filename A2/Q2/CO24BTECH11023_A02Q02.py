import numpy as np
import matplotlib.pyplot as plt
from glob import glob
import os

# ─── Parameters ───────────────────────────────────────────────────────────────
A       = 0.75   # Roll 23
Nx      = 200
dx      = 1.0 / Nx
dt      = 0.8 * dx / A
t_star  = 1.0 / (2.0 * np.pi * A)
T       = 1.5 * t_star
target_times = [0.0, t_star / 2.0, t_star, T]
plots_dir = "Q2/plots"
data_dir  = "Q2/data"
os.makedirs(plots_dir, exist_ok=True)

# ─── Actual simulation times: n*dt where n = ceil(t_target / dt) ─────────────
def getTimes(dt, target_times):
    n_steps = np.ceil(np.array(target_times, dtype=float) / dt).astype(int)
    return n_steps * dt

actual_times = getTimes(dt, target_times)
print(f"t*={t_star:.6f}, T={T:.6f}")
print(f"Actual simulation times (n*dt): {actual_times}")

# ─── Read a snapshot file → (x, u) ───────────────────────────────────────────
def readSnapshot(filepath):
    data = np.loadtxt(filepath)
    return data[:, 0], data[:, 1]

# ─── Read all snapshots sorted by time ───────────────────────────────────────
def readAllSnapshots(data_dir="A2/Q2/data"):
    files = sorted(
        glob(os.path.join(data_dir, "Burgers_T_*.txt")),
        key=lambda p: float(os.path.basename(p).split("_T_")[1].replace(".txt", ""))
    )
    snapshots, times_found = [], []
    for f in files:
        t = float(os.path.basename(f).split("_T_")[1].replace(".txt", ""))
        x, u = readSnapshot(f)
        snapshots.append((x, u))
        times_found.append(t)
    return snapshots, times_found

# ─── Plot 4 subplots: t=0, t*/2, t*, T ───────────────────────────────────────
def plotBurgers(snapshots, times_found, actual_times, t_star, plots_dir):
    labels = ["t = 0", "t = t*/2", "t = t*", "t = T = 1.5t*"]

    fig, axes = plt.subplots(1, 4, figsize=(18, 4), sharey=True)

    for col, (t_actual, label) in enumerate(zip(actual_times, labels)):
        ax = axes[col]

        # Match closest snapshot
        idx = np.argmin([abs(tf - t_actual) for tf in times_found])
        x, u = snapshots[idx]
        t_file = times_found[idx]

        ax.plot(x, u, 'b-', linewidth=1.8)
        ax.set_title(f"{label}\n(t = {t_file:.4f} s)", fontsize=11, fontweight='bold')
        ax.set_xlabel("x", fontsize=12)
        ax.grid(True, alpha=0.4)
        ax.set_xlim(0, 1)

    axes[0].set_ylabel("u(x, t)", fontsize=12, fontweight='bold')
    fig.suptitle(f"Burgers Equation — Lax-Friedrichs (Roll 23, A={A})\n"
                 f"t* = {t_star:.4f} s", fontsize=13, fontweight='bold')
    plt.tight_layout()
    out_path = os.path.join(plots_dir, "Q2_Burgers.pdf")
    plt.savefig(out_path, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"Saved: {out_path}")

# ─── Main ─────────────────────────────────────────────────────────────────────
snapshots, times_found = readAllSnapshots(data_dir)
plotBurgers(snapshots, times_found, actual_times, t_star, plots_dir)

print("Program completed!")