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
        float(os.path.splitext(os.path.basename(path))[0].split("_T_")[1])
        for path in files
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


def plotFTFSBlowup(method, data_dir="data", plots_dir="plots", xlim=None, ylim=(-5, 5)):
    data, times = readSolutionsForMethod(method, data_dir)
    if data.size == 0 or len(times) == 0:
        return

    x = data[:, 0]
    
    plt.figure(figsize=(11, 7))
    
    num_steps = min(50, data.shape[1] - 1)
    
    colors = plt.cm.plasma(np.linspace(0, 1, num_steps)) 
    
    # legend_indices = {0, num_steps // 4, num_steps // 2, (3 * num_steps) // 4, num_steps - 1}
    
    for i in range(num_steps):
        u = data[:, i + 1]
        label = f"t={times[i]:.4g}" if i % 5 == 0 else None
        plt.plot(x, u, color=colors[i], alpha=0.6, linewidth=0.5, label=label)

    if xlim is not None:
        plt.xlim(xlim)
    if ylim is not None:
        plt.ylim(ylim)
        
    plt.xlabel("x", fontsize=14, fontweight='bold')
    plt.ylabel("u(x,t)", fontsize=14, fontweight='bold')
    plt.title("Problem 2: FTCS Oscillation Growth (First 50 steps)", fontsize=16, fontweight='bold')
    
    # Put legend outside to avoid covering the graph
    plt.legend(title="Sample times", fontsize=10, title_fontsize=11, loc='upper left', bbox_to_anchor=(1.02, 1))
    plt.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(os.path.join(plots_dir, "FTCS_Blowup.pdf"), dpi=150)
    plt.close()


def plotFTFSGrowthRate(data_file="data/FTCS_max_u.txt", plots_dir="plots", nsteps=200):
    file_path = data_file
    if not os.path.exists(file_path): return
    
    # Load the tracking file
    data = np.loadtxt(file_path)
    data = data[:nsteps]
    step = np.arange(1, data.shape[0] + 1)
    max_u = np.abs(data[:, 1])
    # slope, intercept = np.polyfit(np.log(max_u),step, 1)
    # observed_rate = slope
    # print(f"Observed growth rate: {observed_rate:.4f}")

    plt.figure(figsize=(8, 6))
    # Plot on a log scale 
    plt.semilogy(step, max_u, 'r-', linewidth=2.5, label="||u||_inf")
    #plotting the observed fit
    # max_u_fit = np.exp(intercept) * np.exp(observed_rate * step)
    # plt.semilogy(step, max_u_fit, 'b--', linewidth=2.5, label=f"Observed fit (p={observed_rate:.2f})")

    plt.xlabel("Time Step", fontsize=14, fontweight='bold')
    plt.ylabel("||u||_inf (log scale)", fontsize=14, fontweight='bold')
    plt.title("Problem 2: FTCS Growth of ||u||_inf (First 200 steps)", fontsize=16, fontweight='bold')
    plt.grid(True, which="both", ls="--", alpha=0.5)
    plt.legend()
    
    plt.tight_layout()
    plt.savefig(os.path.join(plots_dir, "FTCS_GrowthRate.pdf"), dpi=150)
    plt.close()


c = 1.0
xstart = 0.0
xend = 1.0
L = xend - xstart
nu = 1.35
Nx = 100
dt = nu * (1.0 / (Nx * c))
target_times = [0.0, 0.25, 0.5, 0.75, 1.0]
method = "FTCS"

set_timesteps_dir = "data/SetTimesteps"
first_50_steps_dir = "data/First50Steps"
plots_dir = "plots"

data, _ = readSolutionsForMethod(method, set_timesteps_dir)
times = getTimes(dt, target_times)
plotSolutionsForMethod(method, data, times, c, L, xlim=(xstart, xend), ylim=(-1.5, 1.5), plots_dir=plots_dir)
plotFTFSBlowup(method, first_50_steps_dir, plots_dir, xlim=(xstart, xend), ylim=(-5, 5))
plotFTFSGrowthRate("data/FTCS_max_u.txt", plots_dir, nsteps=200)

print("program completed!")