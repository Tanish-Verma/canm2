import numpy as np
import matplotlib.pyplot as plt
from glob import glob
import os

def exactSolution(x,t,c,L):
    return np.sin(2*np.pi*(x - c*t)/L)

def formatApproxTime(t, decimals=2):
    return f"{round(float(t), decimals):.{decimals}f}"

def getTimes(dt, target_times):
    n_steps = np.ceil(np.array(target_times, dtype=float) / dt).astype(int)
    return n_steps * dt


def readSolutionsForMethod(method, data_dir="data"):
    files = sorted(
        glob(os.path.join(data_dir, f"{method}_T_*.txt")),
        key=lambda path: float(os.path.splitext(os.path.basename(path))[0].split("_T_")[1]),
    )

    if not files:
        return np.empty((0, 0)), np.array([])

    data = [np.loadtxt(file) for file in files]
    x = data[0][:, 0]
    u_columns = [item[:, 1] for item in data]

    return np.column_stack([x] + u_columns)

def plotSolutionsForMethod(method, data, times, c, L, xlim=None, ylim=None, plots_dir="plots"):
    x = data[:, 0]
    plt.figure(figsize=(12, 8))
    
    colors = plt.cm.Set1(np.linspace(0, 1, len(times)))

    for index, t in enumerate(times):
        u = data[:, index + 1]
        color = colors[index]
        t_label = formatApproxTime(t)

        plt.plot(x, u, color=color, linewidth=2.5, label=f"{method} T≈{t_label}")
        plt.plot(x, exactSolution(x, t, c, L), "--", color=color, linewidth=2.5, alpha=0.8, label=f"Exact T≈{t_label}")

    if xlim is not None:
        plt.xlim(xlim)
    if ylim is not None:
        plt.ylim(ylim)
        
    plt.xlabel("x", fontsize=14, fontweight='bold')
    plt.ylabel("u(x,t)", fontsize=14, fontweight='bold')
    plt.title(method, fontsize=16, fontweight='bold')
    plt.legend(fontsize=12, loc='upper left', bbox_to_anchor=(1.05, 1), borderaxespad=0)
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    # Save the separate figure
    plt.savefig(os.path.join(plots_dir, f"{method}.pdf"), dpi=150, bbox_inches='tight')
    plt.close()

def plotLaxFconvergence(data_file="error.txt", plots_dir="plots"):
    if not os.path.exists(data_file): 
        print(f"Warning: {data_file} not found.")
        return
    
    data = np.loadtxt(data_file)
    dx = data[:, 0]
    l2_error = data[:, 1]
    slope, intercept = np.polyfit(np.log(dx), np.log(l2_error), 1)
    observed_rate = slope
    print(f"Observed convergence rate: {observed_rate:.4f}")
    
    plt.figure(figsize=(8, 6))
    
    plt.loglog(dx, l2_error, 'bo-', linewidth=2, markersize=8, label="L2 Error")
    dx_fit = np.linspace(dx.min(), dx.max(), 200)
    observed_fit = np.exp(intercept) * dx_fit ** slope
    plt.loglog(dx_fit, observed_fit, 'r-', linewidth=2, label=f"Observed fit (p={observed_rate:.2f})")
    
    plt.loglog(dx, dx * (l2_error[-1]/dx[-1]), 'k--', alpha=0.6, label="O(dx)")

    plt.xlabel("Grid Spacing (dx)", fontsize=14, fontweight='bold')
    plt.ylabel("L2 Error", fontsize=14, fontweight='bold')
    plt.title("Problem 3: Lax-Friedrich Convergence", fontsize=16, fontweight='bold')
    
    plt.grid(True, which="both", ls="--", alpha=0.5)
    plt.legend(fontsize=12)
    
    plt.tight_layout()
    plt.savefig(os.path.join(plots_dir, "LaxFConvergence.pdf"), dpi=150)
    plt.close()

c = 1.0
nu = 0.5
Nx = 100
dt = nu * (1.0 / (Nx * c))
target_times = [0.0,0.25,0.5,0.75,1.0]
times = getTimes(dt, target_times)
xstart = 0.0
xend = 1.0
L = xend - xstart
methods = ["LaxF"]
data_dir = "data"
plots_dir = "plots"
os.makedirs(plots_dir, exist_ok=True)

for method in methods:
    data= readSolutionsForMethod(method, data_dir)
    plotSolutionsForMethod(method, data, times, c, L, xlim=(xstart, xend),ylim=[-1.5, 1.5] ,plots_dir=plots_dir)

plotLaxFconvergence(data_file="error.txt", plots_dir=plots_dir)
print("Program completed!")