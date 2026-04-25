import numpy as np
import matplotlib.pyplot as plt
from glob import glob
import os

# ─── Parameters ───────────────────────────────────────────────────────────────
kappa0    = 0.20   # Roll 23
beta      = 0.40   # Roll 23
r0        = 0.40
Nx        = 50     # intervals => 51 grid points
dx        = 1.0 / Nx
dt        = r0 * dx * dx / kappa0
T         = 0.2
target_times = [0.0, T/3.0, 2.0*T/3.0, T]
plots_dir = "Q3/plots"
data_dir  = "Q3/data"
os.makedirs(plots_dir, exist_ok=True)

# ─── Actual simulation times: n*dt ────────────────────────────────────────────
def getTimes(dt, target_times):
    n_steps = np.ceil(np.array(target_times, dtype=float) / dt).astype(int)
    return n_steps * dt

actual_times = getTimes(dt, target_times)
print(f"Actual output times (n*dt): {actual_times}")

# ─── Read 2D snapshot → (X, Y, U) meshgrids ──────────────────────────────────
def readSnapshot(filepath):
    data = np.loadtxt(filepath)
    x_vals = np.unique(data[:, 0])
    y_vals = np.unique(data[:, 1])
    X, Y = np.meshgrid(x_vals, y_vals, indexing='ij')
    U = data[:, 2].reshape(len(x_vals), len(y_vals))
    return X, Y, U

# ─── Read all ADI snapshots sorted by time ────────────────────────────────────
def readAllSnapshots(prefix, data_dir="data"):
    files = sorted(
        glob(os.path.join(data_dir, f"{prefix}_T_*.txt")),
        key=lambda p: float(os.path.basename(p).split("_T_")[1].replace(".txt", ""))
    )
    snapshots, times_found = [], []
    for f in files:
        t = float(os.path.basename(f).split("_T_")[1].replace(".txt", ""))
        X, Y, U = readSnapshot(f)
        snapshots.append((X, Y, U))
        times_found.append(t)
    return snapshots, times_found

# ─── (b) Plot colour maps of u at 4 times ────────────────────────────────────
def plotColourMaps(snapshots, times_found, actual_times, plots_dir):
    fig, axes = plt.subplots(1, 4, figsize=(18, 4))
    titles = ["t = 0", "t = T/3", "t = 2T/3", "t = T"]

    for col, (t_actual, title) in enumerate(zip(actual_times, titles)):
        ax = axes[col]
        idx = np.argmin([abs(tf - t_actual) for tf in times_found])
        X, Y, U = snapshots[idx]
        t_file = times_found[idx]

        cf = ax.contourf(X, Y, U, levels=20, cmap='hot_r')
        ax.set_title(f"{title}\n(t = {t_file:.4f} s)", fontsize=11, fontweight='bold')
        ax.set_xlabel("x", fontsize=10)
        ax.set_ylabel("y", fontsize=10)
        ax.set_aspect('equal')
        plt.colorbar(cf, ax=ax, shrink=0.85)

    fig.suptitle(f"ADI Diffusion — u(x,y,t) (Roll 23, κ₀={kappa0}, β={beta})",
                 fontsize=13, fontweight='bold')
    plt.tight_layout()
    out_path = os.path.join(plots_dir, "Q3_ADI_ColourMaps.pdf")
    plt.savefig(out_path, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"Saved: {out_path}")

# ─── (b) Plot kappa(x,y) ─────────────────────────────────────────────────────
def plotKappa(data_dir, plots_dir, kappa0, beta):
    X, Y, K = readSnapshot(os.path.join(data_dir, "kappa.txt"))
    fig, ax = plt.subplots(figsize=(5, 4))
    cf = ax.contourf(X, Y, K, levels=20, cmap='viridis')
    ax.set_title(f"κ(x,y) = {kappa0}·(1 + {beta}·sin(πx)·sin(πy))",
                 fontsize=11, fontweight='bold')
    ax.set_xlabel("x", fontsize=10)
    ax.set_ylabel("y", fontsize=10)
    ax.set_aspect('equal')
    plt.colorbar(cf, ax=ax)
    plt.tight_layout()
    out_path = os.path.join(plots_dir, "Q3_Kappa.pdf")
    plt.savefig(out_path, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"Saved: {out_path}")

# ─── (c) Convergence plot (beta=0) ───────────────────────────────────────────
def plotConvergence(data_dir, plots_dir):
    data = np.loadtxt(os.path.join(data_dir, "convergence.txt"))
    dx_vals = data[:, 0]
    err_vals = data[:, 1]

    slope, intercept = np.polyfit(np.log(dx_vals), np.log(err_vals), 1)
    print(f"Convergence rate p = {slope:.3f}")

    dx_fit = np.linspace(dx_vals.min(), dx_vals.max(), 200)
    fit = np.exp(intercept) * dx_fit**slope

    plt.figure(figsize=(7, 5))
    plt.loglog(dx_vals, err_vals, 'bo-', linewidth=2, markersize=8, label="L2 Error (β=0)")
    plt.loglog(dx_fit, fit, 'r--', linewidth=2, label=f"Fit: O(Δx^{slope:.2f})")
    # Reference O(dx^2) line
    ref = err_vals[-1] * (dx_vals / dx_vals[-1])**2
    plt.loglog(dx_vals, ref, 'k:', linewidth=1.5, label="O(Δx²) reference")
    plt.xlabel("Δx", fontsize=13, fontweight='bold')
    plt.ylabel("L2 Error", fontsize=13, fontweight='bold')
    plt.title(f"ADI Convergence (β=0, κ₀={kappa0})", fontsize=13, fontweight='bold')
    plt.legend(fontsize=11)
    plt.grid(True, which='both', alpha=0.4)
    plt.tight_layout()
    out_path = os.path.join(plots_dir, "Q3_Convergence.pdf")
    plt.savefig(out_path, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"Saved: {out_path}")

# ─── (d) Max-u decay plot (beta = 0, 0.4, 0.8) ───────────────────────────────
def plotDecay(data_dir, plots_dir, kappa0, T):
    times = np.loadtxt(os.path.join(data_dir, "decay_time.txt"))
    betas    = [0.0, 0.4, 0.8]
    filenames = ["decay_beta0.txt", "decay_beta4.txt", "decay_beta8.txt"]
    colors   = ['b', 'g', 'r']

    # Constant-kappa0 analytical reference: exp(-2*pi^2*kappa0*t)
    t_ref = np.linspace(0, T, 300)
    u_ref = np.exp(-2 * np.pi**2 * kappa0 * t_ref)

    plt.figure(figsize=(8, 5))
    plt.plot(t_ref, u_ref, 'k--', linewidth=2, label="Exact (β=0)")

    for beta_val, fname, col in zip(betas, filenames, colors):
        max_u = np.loadtxt(os.path.join(data_dir, fname))
        plt.plot(times, max_u, color=col, linewidth=2, label=f"β = {beta_val}")

    plt.xlabel("Time (s)", fontsize=13, fontweight='bold')
    plt.ylabel("max u", fontsize=13, fontweight='bold')
    plt.title(f"Max-u Decay vs Time (κ₀={kappa0})", fontsize=13, fontweight='bold')
    plt.legend(fontsize=11)
    plt.grid(True, alpha=0.4)
    plt.tight_layout()
    out_path = os.path.join(plots_dir, "Q3_Decay.pdf")
    plt.savefig(out_path, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"Saved: {out_path}")

# ─── Main ─────────────────────────────────────────────────────────────────────
snapshots, times_found = readAllSnapshots("ADI", data_dir)

plotColourMaps(snapshots, times_found, actual_times, plots_dir)
plotKappa(data_dir, plots_dir, kappa0, beta)
plotConvergence(data_dir, plots_dir)
plotDecay(data_dir, plots_dir, kappa0, T)

print("Program completed!")