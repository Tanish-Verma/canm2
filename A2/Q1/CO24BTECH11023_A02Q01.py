import numpy as np
import matplotlib.pyplot as plt
from glob import glob
import os

# ─── Parameters ───────────────────────────────────────────────────────────────
c             = 1.0
nu            = 0.45   # Roll 23
Nx, Ny        = 100, 100
dx            = 1.0 / Nx
dt            = nu * dx / c
target_times  = [0.0, 0.25, 0.5, 0.75, 1.0]
plots_dir     = "Q1/plots"
data_dir      = "Q1/data"
os.makedirs(plots_dir, exist_ok=True)

# ─── Actual simulation times: n*dt where n = ceil(t_target / dt) ─────────────
# This matches exactly what the C++ solver writes to disk
def getTimes(dt, target_times):
    n_steps = np.ceil(np.array(target_times, dtype=float) / dt).astype(int)
    return n_steps * dt

actual_times = getTimes(dt, target_times)  # these are the real timestamps on disk

# ─── Exact solution ───────────────────────────────────────────────────────────
def exactSolution(X, Y, t, c):
    return np.sin(2 * np.pi * (X - c * t)) * np.sin(2 * np.pi * (Y - c * t))

# ─── Read one snapshot file → (X, Y, U) meshgrids ────────────────────────────
def readSnapshot(filepath):
    data = np.loadtxt(filepath)
    x_vals = np.unique(data[:, 0])
    y_vals = np.unique(data[:, 1])
    X, Y = np.meshgrid(x_vals, y_vals, indexing='ij')
    U = data[:, 2].reshape(len(x_vals), len(y_vals))
    return X, Y, U

# ─── Read all snapshots for a scheme, sorted by time ─────────────────────────
def readAllSnapshots(scheme, data_dir="data"):
    files = sorted(
        glob(os.path.join(data_dir, f"{scheme}_T_*.txt")),
        key=lambda p: float(os.path.basename(p).split("_T_")[1].replace(".txt", ""))
    )
    snapshots, times_found = [], []
    for f in files:
        t = float(os.path.basename(f).split("_T_")[1].replace(".txt", ""))
        X, Y, U = readSnapshot(f)
        snapshots.append((X, Y, U))
        times_found.append(t)
    return snapshots, times_found

# ─── L2 error ────────────────────────────────────────────────────────────────
def computeL2Error(U, X, Y, t, c):
    dx = X[1, 0] - X[0, 0]
    dy = Y[0, 1] - Y[0, 0]
    U_ex = exactSolution(X, Y, t, c)
    return np.sqrt(dx * dy * np.sum((U - U_ex)**2))

# ─── Plot: 2 rows (Upwind, LaxF) x 5 cols (output times) ────────────────────
def plotBothSchemes(schemes, data_dir, plots_dir, actual_times, c):
    fig, axes = plt.subplots(2, len(actual_times), figsize=(20, 8))
    vmin, vmax = -1.0, 1.0

    for row, scheme in enumerate(schemes):
        snapshots, times_found = readAllSnapshots(scheme, data_dir)
        if not snapshots:
            print(f"No data found for {scheme}")
            continue

        for col, t_actual in enumerate(actual_times):
            ax = axes[row][col]

            # Match the file whose timestamp equals t_actual
            idx = np.argmin([abs(tf - t_actual) for tf in times_found])
            X, Y, U = snapshots[idx]
            t_file  = times_found[idx]

            # Exact solution evaluated at the same t_actual (= n*dt)
            U_ex = exactSolution(X, Y, t_file, c)

            cf = ax.contourf(X, Y, U, levels=20, cmap='RdBu_r', vmin=vmin, vmax=vmax)
            ax.contour(X, Y, U_ex, levels=10, colors='k', linewidths=0.8,
                       linestyles='dashed', alpha=0.7)

            ax.set_title(f"t = {t_file:.4f}s", fontsize=10)
            ax.set_xlabel("x", fontsize=9)
            ax.set_ylabel("y", fontsize=9)
            ax.set_aspect('equal')
            plt.colorbar(cf, ax=ax, shrink=0.8)

        axes[row][0].set_ylabel(f"{scheme}\ny", fontsize=10, fontweight='bold')

    scheme_labels = " (top) & ".join(schemes) + " (bottom)"
    fig.suptitle(f"2D Linear Advection — {scheme_labels}\n"
                 "Filled: Numerical | Dashed: Exact", fontsize=13, fontweight='bold')
    plt.tight_layout()
    out_path = os.path.join(plots_dir, "Q1_2D_Advection.pdf")
    plt.savefig(out_path, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"Saved: {out_path}")

# ─── Plot L2 error vs time for both schemes ───────────────────────────────────
def plotL2ErrorInTime(schemes, data_dir, plots_dir, c):
    plt.figure(figsize=(8, 5))
    for scheme in schemes:
        snapshots, times_found = readAllSnapshots(scheme, data_dir)
        errors = [computeL2Error(U, X, Y, t, c) for (X, Y, U), t in zip(snapshots, times_found)]
        plt.plot(times_found, errors, 'o-', linewidth=2, markersize=6, label=scheme)
    plt.xlabel("Time (s)", fontsize=13, fontweight='bold')
    plt.ylabel("L2 Error", fontsize=13, fontweight='bold')
    plt.title("L2 Error vs Time — 2D Advection (Roll 23, ν=0.45)", fontsize=13, fontweight='bold')
    plt.legend(fontsize=12)
    plt.grid(True, alpha=0.4)
    plt.tight_layout()
    out_path = os.path.join(plots_dir, "Q1_L2_Error.pdf")
    plt.savefig(out_path, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"Saved: {out_path}")

# ─── Main ─────────────────────────────────────────────────────────────────────
schemes = ["Upwind", "LaxF"]

print(f"Actual simulation times (n*dt): {actual_times}")

plotBothSchemes(schemes, data_dir, plots_dir, actual_times, c)
plotL2ErrorInTime(schemes, data_dir, plots_dir, c)

print("\nL2 errors at each output snapshot:")
for scheme in schemes:
    snapshots, times_found = readAllSnapshots(scheme, data_dir)
    for (X, Y, U), t in zip(snapshots, times_found):
        e = computeL2Error(U, X, Y, t, c)
        print(f"  {scheme}  t={t:.6f}  L2={e:.6e}")

print("Program completed!")