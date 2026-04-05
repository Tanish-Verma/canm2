import numpy as np
import matplotlib.pyplot as plt
from glob import glob
import os

def exactSolution(x,t,c,L):
    return np.sin(2*np.pi*(x - c*t)/L)

def ApproxTime(t, decimals=2):
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
        t_label = ApproxTime(t)

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

times = [0.0,0.25,0.5,0.75,1.0]
nu = 0.48
Nx = 100
c = 1.0
dt = nu * (1.0 / (Nx*c))
times = getTimes(dt, times)
xstart = 0.0
xend = 1.0
L = xend - xstart
methods = ["FTBS", "FTFS"]
data_dir = "data"
plots_dir = "plots"

for method in methods:
    data= readSolutionsForMethod(method, data_dir)
    plotSolutionsForMethod(method, data, times, c, L, xlim=(xstart, xend),ylim=[-1.5, 1.5] ,plots_dir=plots_dir)

print("Program completed!")