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
    times = np.array([
        float(os.path.splitext(os.path.basename(file))[0].split("_T_")[1])
        for file in files
    ])
    u_columns = [item[:, 1] for item in data]

    return np.column_stack([x] + u_columns), times

def plotSolutionsForMethod(method, data, times, c, L, xlim=None, ylim=None, plots_dir="plots"):
    if data.size == 0 or len(times) == 0:
        return

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


def plotStableTimeSteps(method, c, L, data_dir="data/StableTimeSteps", plots_dir="plots", xlim=None, ylim=None):
    data, times = readSolutionsForMethod(method, data_dir)
    if data.size == 0 or len(times) == 0:
        return

    x = data[:, 0]
    plt.figure(figsize=(12, 8))
    colors = plt.cm.plasma(np.linspace(0, 1, len(times)))

    for i, t in enumerate(times):
        u = data[:, i + 1]
        color = colors[i]
        if i % 3 == 0:
            num_label = f"Num t={t:.4g}"
            exact_label = f"Exact t={t:.4g}"
            plt.plot(x, u, color=color, linewidth=1.8, label=num_label)
            plt.plot(x, exactSolution(x, t, c, L), "--", color=color, linewidth=2.4, alpha=0.9, label=exact_label)

    if xlim is not None:
        plt.xlim(xlim)
    if ylim is not None:
        plt.ylim(ylim)

    plt.xlabel("x", fontsize=14, fontweight='bold')
    plt.ylabel("u(x,t)", fontsize=14, fontweight='bold')
    plt.title(f"{method} StableTimeSteps", fontsize=16, fontweight='bold')
    plt.grid(True, alpha=0.3)
    plt.legend(fontsize=9, ncol=2, loc='upper left', bbox_to_anchor=(1.02, 1))
    plt.tight_layout()
    plt.savefig(os.path.join(plots_dir, f"{method}_StableTimeSteps.pdf"), dpi=150, bbox_inches='tight')
    plt.close()

c = 1.0
nu = 1.82
xstart = 0.0
xend = 1.0
Nx = 100
dt = nu * (1.0 / (Nx * c))
target_times = [0.0,0.25,0.5,0.75,1.0]
times = getTimes(dt, target_times)
L = xend - xstart
methods = ["LaxW"]
data_dir = "data"
plots_dir = "plots"

for method in methods:
    data, _ = readSolutionsForMethod(method, "data/SetTimeSteps")
    plotSolutionsForMethod(method, data, times, c, L, xlim=(xstart, xend),ylim=[-1.5, 1.5] ,plots_dir=plots_dir)
    plotStableTimeSteps(method, c, L, data_dir="data/StableTimeSteps", plots_dir=plots_dir, xlim=(xstart, xend), ylim=[-1.5, 1.5])

print("Program completed!")