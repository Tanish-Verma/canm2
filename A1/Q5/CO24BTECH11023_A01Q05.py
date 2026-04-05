import numpy as np
import matplotlib.pyplot as plt
from glob import glob
import os

def exactSolution(x,t,c,L):
    return np.sin(2*np.pi*(x - c*t)/L)

def formatApproxTime(t, decimals=2):
    return f"{round(float(t), decimals):.{decimals}f}"


def getAvailableTimes(dt, target_times):
    n_steps = np.ceil(np.array(target_times, dtype=float) / dt).astype(int)
    return n_steps * dt

def readSolutionsForMethod(method, data_dir="data"):
    files = sorted(
        glob(os.path.join(data_dir, f"{method}_T_*.txt")),
        key=lambda path: float(os.path.splitext(os.path.basename(path))[0].split("_T_")[1]),
    )

    if not files:
        raise FileNotFoundError(f"No files found for method '{method}' in '{data_dir}'.")

    data = [np.loadtxt(file) for file in files]
    x = data[0][:, 0]
    u_columns = [item[:, 1] for item in data]

    return np.column_stack([x] + u_columns)

def plotSolutionsForMethod(method, data, times, c, L, nu,xlim=None, ylim=None, plots_dir="plots"):
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
    plt.title(f"{method} at nu = {nu}", fontsize=16, fontweight='bold')
    plt.legend(fontsize=12, loc='upper left', bbox_to_anchor=(1.05, 1), borderaxespad=0)
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    # Save the separate figure
    plt.savefig(os.path.join(plots_dir, f"{method}_{nu}.pdf"), dpi=150, bbox_inches='tight')
    plt.close()

def plotSolutionsCombined(methods, all_data, times, c, L, xlim=None, ylim=None, plots_dir="plots"):
    fig, axes = plt.subplots(1, len(methods), figsize=(14, 6), sharey=True)
    if len(methods) == 1:
        axes = [axes]

    colors = plt.cm.Set1(np.linspace(0, 1, len(times)))

    for ax, method in zip(axes, methods):
        data = all_data[method]
        x = data[:, 0]

        for index, t in enumerate(times):
            u = data[:, index + 1]
            color = colors[index]
            t_label = formatApproxTime(t)
            ax.plot(x, u, color=color, linewidth=2.2, label=f"{method} T≈{t_label}")
            ax.plot(
                x,
                exactSolution(x, t, c, L),
                "--",
                color=color,
                linewidth=2.2,
                alpha=0.8,
                label=f"Exact T≈{t_label}",
            )

        if xlim is not None:
            ax.set_xlim(xlim)
        if ylim is not None:
            ax.set_ylim(ylim)

        ax.set_xlabel("x", fontsize=12, fontweight='bold')
        ax.set_title(method, fontsize=14, fontweight='bold')
        ax.grid(True, alpha=0.3)

    axes[0].set_ylabel("u(x,t)", fontsize=12, fontweight='bold')
    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, fontsize=10, loc='upper center', ncol=5, bbox_to_anchor=(0.5, 1.04))
    fig.tight_layout(rect=[0, 0, 1, 0.93])

    os.makedirs(plots_dir, exist_ok=True)
    plt.savefig(os.path.join(plots_dir, "MacC_CN.pdf"), dpi=150, bbox_inches='tight')
    plt.close()

nu = 2.14
c = 1.0
xstart = 0.0
xend = 1.0
Nx = 100.0
dx = (xend - xstart) / Nx
dt = nu * dx/c 
L = xend - xstart
methods = ["CrankNicolson", "MacCormack"]
data_dir = "data"
plots_dir = "plots"
   
all_data = {}
for method in methods:
    all_data[method] = readSolutionsForMethod(method, data_dir)

# Extract actual available times from the first method's data files
times = [0.0, 0.25, 0.5, 0.75, 1.0]
times = getAvailableTimes(dt, times)

plotSolutionsCombined(
    methods,
    all_data,
    times,
    c,
    L,
    xlim=(xstart, xend),
    ylim=[-1.5, 1.5],
    plots_dir=plots_dir,
)

nu = [2,5,10]
for nu in nu:
    data = readSolutionsForMethod("CN",f"data/nu{nu}")
    plotSolutionsForMethod("CN",data,times,c,L,xlim=(xstart,xend),ylim=[-1.5,1.5],nu=nu)

print("Program completed!")