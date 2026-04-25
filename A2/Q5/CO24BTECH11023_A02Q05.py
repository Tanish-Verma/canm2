import numpy as np
import matplotlib.pyplot as plt
from glob import glob
import os

# ─── Parameters ───────────────────────────────────────────────────────────────
Pe      = 100.0   # Roll 23
kappa   = 1.0 / Pe
c       = 1.0
T_tr    = 0.5
nu_b    = 0.5
Nx_b    = 100
dx_b    = 1.0 / Nx_b
dt_b    = nu_b * dx_b / c
target_times = [0.0, 0.125, 0.25, 0.375, 0.5]
plots_dir = "Q5/plots"
data_dir  = "Q5/data"
os.makedirs(plots_dir, exist_ok=True)

# ─── Actual simulation times: n*dt ────────────────────────────────────────────
def getTimes(dt, target_times):
    n_steps = np.ceil(np.array(target_times, dtype=float) / dt).astype(int)
    return n_steps * dt

actual_times_b = getTimes(dt_b, target_times)
print(f"Pe={Pe}, kappa={kappa}")
print(f"Actual transient output times: {actual_times_b}")

# ─── Exact solutions ──────────────────────────────────────────────────────────
def steadyExact(x, Pe):
    return (np.exp(Pe * x) - 1.0) / (np.exp(Pe) - 1.0)

def transientExact(x, t, kappa):
    return np.exp(-4.0 * np.pi**2 * kappa * t) * np.sin(2.0 * np.pi * (x - t))

# ─── Read snapshot files ──────────────────────────────────────────────────────
def readAllSnapshots(scheme, data_dir="data"):
    files = sorted(
        glob(os.path.join(data_dir, f"{scheme}_T_*.txt")),
        key=lambda p: float(os.path.basename(p).split("_T_")[1].replace(".txt", ""))
    )
    snapshots, times_found = [], []
    for f in files:
        t = float(os.path.basename(f).split("_T_")[1].replace(".txt", ""))
        d = np.loadtxt(f)
        snapshots.append((d[:, 0], d[:, 1]))
        times_found.append(t)
    return snapshots, times_found

# ─── Part A(a): Steady-state plot ────────────────────────────────────────────
def plotSteadyState(data_dir, plots_dir, Pe):
    # Central differences
    dc = np.loadtxt(os.path.join(data_dir, "steady_central.txt"))
    xc, uc, uc_ex = dc[:, 0], dc[:, 1], dc[:, 2]

    # Upwind (N=50)
    du = np.loadtxt(os.path.join(data_dir, "steady_upwind_N50.txt"))
    xu, uu, uu_ex = du[:, 0], du[:, 1], du[:, 2]

    fig, axes = plt.subplots(1, 2, figsize=(13, 5))

    # Central diff
    ax = axes[0]
    ax.plot(xc, uc_ex, 'k-', linewidth=2, label="Exact")
    ax.plot(xc, uc, 'b--o', linewidth=1.8, markersize=4, label="Central diff")
    ax.set_title(f"Central Differences (Pe_h = {Pe}/50 = {Pe/50:.1f})",
                 fontsize=11, fontweight='bold')
    ax.set_xlabel("x", fontsize=12)
    ax.set_ylabel("u(x)", fontsize=12)
    ax.legend(fontsize=11)
    ax.grid(True, alpha=0.4)

    # Upwind
    ax = axes[1]
    ax.plot(xu, uu_ex, 'k-', linewidth=2, label="Exact")
    ax.plot(xu, uu, 'r--o', linewidth=1.8, markersize=4, label="Upwind")
    ax.set_title(f"Upwind Differences (Pe_h = {Pe/50:.1f})",
                 fontsize=11, fontweight='bold')
    ax.set_xlabel("x", fontsize=12)
    ax.set_ylabel("u(x)", fontsize=12)
    ax.legend(fontsize=11)
    ax.grid(True, alpha=0.4)

    fig.suptitle(f"Steady-State Advection-Diffusion (Pe={Pe:.0f}, N=50)",
                 fontsize=13, fontweight='bold')
    plt.tight_layout()
    out_path = os.path.join(plots_dir, "Q5A_SteadyState.pdf")
    plt.savefig(out_path, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"Saved: {out_path}")

# ─── Part A(c): L_inf convergence (upwind) ───────────────────────────────────
def plotSteadyConvergence(data_dir, plots_dir):
    data = np.loadtxt(os.path.join(data_dir, "steady_convergence.txt"))
    dx_vals = data[:, 0]
    err_vals = data[:, 1]

    # Fit slope using all points
    slope, intercept = np.polyfit(np.log(dx_vals), np.log(err_vals), 1)
    print(f"Steady upwind L_inf convergence rate: p = {slope:.3f}")

    dx_fit = np.linspace(dx_vals.min(), dx_vals.max(), 200)
    fit = np.exp(intercept) * dx_fit**slope

    plt.figure(figsize=(7, 5))
    plt.loglog(dx_vals, err_vals, 'ro-', linewidth=2, markersize=8, label="L∞ Error (Upwind)")
    plt.loglog(dx_fit, fit, 'b--', linewidth=2, label=f"Fit: O(Δx^{slope:.2f})")
    ref = err_vals[0] * (dx_vals / dx_vals[0])**1
    plt.loglog(dx_vals, ref, 'k:', linewidth=1.5, label="O(Δx¹) reference")
    plt.xlabel("Δx", fontsize=13, fontweight='bold')
    plt.ylabel("L∞ Error", fontsize=13, fontweight='bold')
    plt.title(f"Steady-State Upwind Convergence (Pe={Pe:.0f})",
              fontsize=13, fontweight='bold')
    plt.legend(fontsize=11)
    plt.grid(True, which='both', alpha=0.4)
    plt.tight_layout()
    out_path = os.path.join(plots_dir, "Q5A_SteadyConvergence.pdf")
    plt.savefig(out_path, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"Saved: {out_path}")

# ─── Part B(a): CN transient plot ────────────────────────────────────────────
def plotCNTransient(data_dir, plots_dir, actual_times, kappa):
    snapshots, times_found = readAllSnapshots("CN_transient", data_dir)

    fig, axes = plt.subplots(1, 5, figsize=(20, 4), sharey=True)
    time_labels = ["t = 0", "t = 0.125", "t = 0.25", "t = 0.375", "t = 0.5"]

    for col, (t_target, label) in enumerate(zip(actual_times, time_labels)):
        ax = axes[col]
        idx = np.argmin([abs(tf - t_target) for tf in times_found])
        x, u = snapshots[idx]
        t_file = times_found[idx]

        ax.plot(x, u, 'b-', linewidth=1.8, label="CN")
        ax.plot(x, transientExact(x, t_file, kappa), 'r--', linewidth=1.8, label="Exact")
        ax.set_title(f"{label}\n(t={t_file:.4f}s)", fontsize=10, fontweight='bold')
        ax.set_xlabel("x", fontsize=10)
        ax.grid(True, alpha=0.3)
        ax.set_xlim(0, 1)

    axes[0].set_ylabel("u(x,t)", fontsize=11, fontweight='bold')
    axes[0].legend(fontsize=9)
    fig.suptitle(f"CN Transient (Pe={Pe:.0f}, κ={kappa:.4f}, N=100, ν=0.5, T=0.5)",
                 fontsize=13, fontweight='bold')
    plt.tight_layout()
    out_path = os.path.join(plots_dir, "Q5B_CNTransient.pdf")
    plt.savefig(out_path, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"Saved: {out_path}")

# ─── Part B(b): Spatial convergence ──────────────────────────────────────────
def plotSpatialConvergence(data_dir, plots_dir):
    data = np.loadtxt(os.path.join(data_dir, "spatial_convergence.txt"))
    dx_vals = data[:, 0]
    err_vals = data[:, 1]

    slope, intercept = np.polyfit(np.log(dx_vals), np.log(err_vals), 1)
    print(f"Spatial convergence rate: p = {slope:.3f}")

    dx_fit = np.linspace(dx_vals.min(), dx_vals.max(), 200)
    fit = np.exp(intercept) * dx_fit**slope

    plt.figure(figsize=(7, 5))
    plt.loglog(dx_vals, err_vals, 'bo-', linewidth=2, markersize=8, label="L2 Error")
    plt.loglog(dx_fit, fit, 'r--', linewidth=2, label=f"Fit: O(Δx^{slope:.2f})")
    ref1 = err_vals[-1] * (dx_vals / dx_vals[-1])**1
    ref2 = err_vals[-1] * (dx_vals / dx_vals[-1])**2
    plt.loglog(dx_vals, ref1, 'k:', linewidth=1.5, label="O(Δx¹)")
    plt.loglog(dx_vals, ref2, 'k--', linewidth=1.5, label="O(Δx²)")
    plt.xlabel("Δx", fontsize=13, fontweight='bold')
    plt.ylabel("L2 Error", fontsize=13, fontweight='bold')
    plt.title(f"Spatial Convergence — CN (Pe={Pe:.0f}, ν=0.5, T=0.5)",
              fontsize=13, fontweight='bold')
    plt.legend(fontsize=11)
    plt.grid(True, which='both', alpha=0.4)
    plt.tight_layout()
    out_path = os.path.join(plots_dir, "Q5B_SpatialConv.pdf")
    plt.savefig(out_path, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"Saved: {out_path}")

# ─── Part B(c): Temporal convergence ─────────────────────────────────────────
def plotTemporalConvergence(data_dir, plots_dir):
    data = np.loadtxt(os.path.join(data_dir, "temporal_convergence.txt"))
    dt_vals = data[:, 0]
    err_vals = data[:, 1]

    slope, intercept = np.polyfit(np.log(dt_vals), np.log(err_vals), 1)
    print(f"Temporal convergence rate: q = {slope:.3f}")

    dt_fit = np.linspace(dt_vals.min(), dt_vals.max(), 200)
    fit = np.exp(intercept) * dt_fit**slope

    plt.figure(figsize=(7, 5))
    plt.loglog(dt_vals, err_vals, 'go-', linewidth=2, markersize=8, label="L2 Error")
    plt.loglog(dt_fit, fit, 'r--', linewidth=2, label=f"Fit: O(Δt^{slope:.2f})")
    ref2 = err_vals[-1] * (dt_vals / dt_vals[-1])**2
    plt.loglog(dt_vals, ref2, 'k--', linewidth=1.5, label="O(Δt²)")
    plt.xlabel("Δt", fontsize=13, fontweight='bold')
    plt.ylabel("L2 Error", fontsize=13, fontweight='bold')
    plt.title(f"Temporal Convergence — CN (Pe={Pe:.0f}, N=200, T=0.5)",
              fontsize=13, fontweight='bold')
    plt.legend(fontsize=11)
    plt.grid(True, which='both', alpha=0.4)
    plt.tight_layout()
    out_path = os.path.join(plots_dir, "Q5B_TemporalConv.pdf")
    plt.savefig(out_path, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"Saved: {out_path}")

# ─── Part B(d): Pe comparison table ──────────────────────────────────────────
def printPeComparison():
    # Values from C++ output
    pe_data = [(5.0, 1.260430e-03), (100.0, 5.447664e-02)]
    print("\n=== Part B(d): Pe Comparison (N=100, nu=0.5, T=0.5) ===")
    print(f"{'Pe':>6} | {'kappa':>10} | {'L2 Error':>12}")
    print("-" * 38)
    for Pe_val, err in pe_data:
        kap = 1.0 / Pe_val
        print(f"{Pe_val:>6.0f} | {kap:>10.6f} | {err:>12.6e}")
    print()

# ─── Main ─────────────────────────────────────────────────────────────────────
plotSteadyState(data_dir, plots_dir, Pe)
plotSteadyConvergence(data_dir, plots_dir)
plotCNTransient(data_dir, plots_dir, actual_times_b, kappa)
plotSpatialConvergence(data_dir, plots_dir)
plotTemporalConvergence(data_dir, plots_dir)
printPeComparison()

print("Program completed!")