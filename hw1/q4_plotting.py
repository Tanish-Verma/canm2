import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import os

LABEL_FS = 16
TITLE_FS = 17
TICK_FS = 13
LEGEND_FS = 14
LINE_LW = 3.2
EXACT_LW = 3.6

def ensure_plots_dir():
    plots_dir = 'plots'
    if not os.path.exists(plots_dir):
        os.makedirs(plots_dir)
        print(f"Created '{plots_dir}' directory")
    return plots_dir

def save_figure(fig, filename, plots_dir='plots', dpi=300):
    """Helper function to save figures with consistent settings."""
    filepath = os.path.join(plots_dir, filename)
    fig.savefig(filepath, dpi=dpi, bbox_inches='tight')
    plt.close(fig)
    print(f"Plot saved as '{filepath}'")

def load_data(filename):
    data = np.loadtxt(filename)
    
    x = data[:, 0]
    y = data[:, 1]
    T = data[:, 2]
    Tex = data[:, 3]
    
    return x, y, T, Tex

def reshape_data(x, y, T, Tex, nx, ny):
    x2d = x.reshape(nx, ny)
    y2d = y.reshape(nx, ny)
    T2d = T.reshape(nx, ny)
    Tex2d = Tex.reshape(nx, ny)
    
    return x2d, y2d, T2d, Tex2d

def plot_contour(filename):
    if not os.path.exists(filename):
        print(f"Error: File '{filename}' not found!")
        return
    
    plots_dir = ensure_plots_dir()
    
    basename = os.path.basename(filename)
    parts = basename.replace('.dat', '').split('_')
    nx = int(parts[2])
    ny = int(parts[3])
    
    x, y, T, Tex = load_data(filename)
    X, Y, T_2d, Tex_2d = reshape_data(x, y, T, Tex, nx, ny)
    
    fig1, ax = plt.subplots(figsize=(10, 8))
    
    cs = ax.contourf(X, Y, T_2d, levels=50, cmap='viridis')
    ax.set_xlabel('x (m)', fontsize=LABEL_FS)
    ax.set_ylabel('y (m)', fontsize=LABEL_FS)
    ax.set_title('Numerical Solution (T)', fontsize=TITLE_FS, fontweight='bold')
    cbar = fig1.colorbar(cs, ax=ax)
    cbar.set_label('Temperature (K)', fontsize=LABEL_FS)
    cbar.ax.tick_params(labelsize=TICK_FS)
    ax.tick_params(labelsize=TICK_FS)
    plt.tight_layout()
    save_figure(fig1, 'contour_numerical.png', plots_dir, dpi=300)
    
def plot_error(filename,a):
    
    if not os.path.exists(filename):
        print(f"Error: File '{filename}' not found!")
        return

    plots_dir = ensure_plots_dir()

    data = np.loadtxt(filename)
    if data.ndim == 1:
        data = data.reshape(1, -1)

    nx = data[:, 0]
    l2err = data[:, 2]
    
    # Calculate grid spacing from nx (assuming domain spans [0, 1])
    grid_spacing = a / (nx - 1)
    # slope, _ = np.polyfit(np.log(grid_spacing), np.log(l2err), 1)
    slope, _ = np.polyfit(np.log(grid_spacing[2:-2]), np.log(l2err[2:-2]), 1)
    print(f"Approximate slope (order of convergence): {slope:.3f}")

    fig, ax = plt.subplots(figsize=(9, 6))
    ax.loglog(grid_spacing, l2err, marker='o', linewidth=LINE_LW, markersize=9)
    ax.set_xlabel('Grid spacing (dx)', fontsize=LABEL_FS)
    ax.set_ylabel('L2 Error Norm', fontsize=LABEL_FS)
    ax.set_title('Error Norm vs Grid Spacing (log Scale)', fontsize=TITLE_FS, fontweight='bold')
    ax.grid(True, which='both', linestyle='--', alpha=0.3, linewidth=1)
    ax.tick_params(labelsize=TICK_FS)
    
    def power_of_2_formatter(x,pos):
        log_val = np.log2(x)
        power = int(np.round(log_val))
        return f'$2^{{{power}}}$'
    
    ax.set_xticks(grid_spacing)
    ax.xaxis.set_major_formatter(ticker.FuncFormatter(power_of_2_formatter))

    plt.tight_layout()
    save_figure(fig, 'error_norm_log_log.png', plots_dir, dpi=300)

def plot_line_profiles_with_analytical_solution(grid_sizes, output_dir='output', iter_stamp=0, x_target=0.9, y_target=0.4):
    plots_dir = ensure_plots_dir()
    datasets = []

    for nx in grid_sizes:
        nx = int(nx)
        ny = (nx + 1) // 2
        fname = f"T_xy_{nx:03d}_{ny:03d}_{int(iter_stamp):04d}.dat"
        fpath = os.path.join(output_dir, fname)

        if not os.path.exists(fpath):
            print(f"Warning: File not found, skipping: {fpath}")
            continue

        x, y, t_num, t_ex = load_data(fpath)
        x2d, y2d, t2d, tex2d = reshape_data(x, y, t_num, t_ex, nx, ny)

        datasets.append({
            'nx': nx,
            'ny': ny,
            'x1d': x2d[:, 0],
            'y1d': y2d[0, :],
            'T2d': t2d,
            'Tex2d': tex2d,
        })

    if len(datasets) == 0:
        print('No valid data files found for provided grid sizes.')
        return

    datasets.sort(key=lambda d: d['nx'])
    ref = datasets[-1]

    # -------- Plot 1: T vs y at x = x_target --------
    fig_y, ax_y = plt.subplots(figsize=(9, 6))

    for d in datasets:
        if d!=ref:
            ix = int(np.argmin(np.abs(d['x1d'] - x_target)))
                ax_y.plot(d['y1d'], d['T2d'][ix, :], linewidth=LINE_LW,
                    label=f"Grid size {d['nx']}x{d['ny']}")

    ix_ref = int(np.argmin(np.abs(ref['x1d'] - x_target)))
    ax_y.plot(ref['y1d'], ref['Tex2d'][ix_ref, :], '--', linewidth=EXACT_LW, label='Exact (Tex)',color = "#000000")
    ax_y.set_xlabel('y (m)', fontsize=LABEL_FS)
    ax_y.set_ylabel('Temperature (K)', fontsize=LABEL_FS)

    x_tag = str(x_target).replace('.', 'p')
    ax_y.set_title(f'Temperature Profile Along y-axis for x = {x_target}', fontsize=TITLE_FS, fontweight='bold')
    ax_y.grid(True, alpha=0.3, linewidth=1)
    ax_y.legend(fontsize=LEGEND_FS)
    ax_y.tick_params(labelsize=TICK_FS)
    fig_y.tight_layout()
    save_figure(fig_y, f'y_profile_for_x{x_tag}.png', plots_dir,dpi=300)

    # -------- Plot 2: T vs x at y = y_target --------
    fig_x, ax_x = plt.subplots(figsize=(9, 6))

    for d in datasets:
        if d != ref:
            jy = int(np.argmin(np.abs(d['y1d'] - y_target)))
                ax_x.plot(d['x1d'], d['T2d'][:, jy], linewidth=LINE_LW,
                    label=f"Grid size {d['nx']}x{d['ny']}")

    jy_ref = int(np.argmin(np.abs(ref['y1d'] - y_target)))
    ax_x.plot(ref['x1d'], ref['Tex2d'][:, jy_ref], '--', linewidth=EXACT_LW, label='Exact (Tex)',color = "#000000")
    ax_x.set_xlabel('x (m)', fontsize=LABEL_FS)
    ax_x.set_ylabel('Temperature (K)', fontsize=LABEL_FS)
    y_tag = str(y_target).replace('.', 'p')
    ax_x.set_title(f'Temperature Profile Along x-axis for y = {y_target}', fontsize=TITLE_FS, fontweight='bold')
    ax_x.grid(True, alpha=0.3, linewidth=1)
    ax_x.legend(fontsize=LEGEND_FS)
    ax_x.tick_params(labelsize=TICK_FS)
    fig_x.tight_layout()
    save_figure(fig_x, f'x_profile_for_y{y_tag}.png', plots_dir,dpi=300)


if __name__ == "__main__":
    a = 1
    b = 0.5
    grid_sizes = [5, 9, 17, 65, 257, 513]
    x_target = 0.9*a
    y_target = 0.8 * b
    dat_file = "output/T_xy_513_257_0000.dat"
    error_file = "output/error.dat"
    plot_contour(dat_file)
    plot_error(error_file, a)

    plot_line_profiles_with_analytical_solution(grid_sizes, x_target=x_target, y_target=y_target)