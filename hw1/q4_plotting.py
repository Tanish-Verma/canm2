import numpy as np
import matplotlib.pyplot as plt
import os

def ensure_plots_dir():
    plots_dir = 'plots'
    if not os.path.exists(plots_dir):
        os.makedirs(plots_dir)
        print(f"Created '{plots_dir}' directory")
    return plots_dir

def save_figure(fig, filename, plots_dir='plots', dpi=200):
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
    
    fig1, ax1 = plt.subplots(figsize=(10, 8))
    
    cs1 = ax1.contourf(X, Y, T_2d, levels=20, cmap='viridis')
    ax1.contour(X, Y, T_2d, levels=10, colors='black', linewidths=0.5, alpha=0.3)
    ax1.set_xlabel('x (m)')
    ax1.set_ylabel('y (m)')
    ax1.set_title('Numerical Solution (T)')
    fig1.colorbar(cs1, ax=ax1)
    
    plt.tight_layout()
    save_figure(fig1, 'contour_numerical.png', plots_dir, dpi=150)
    
def plot_error(filename):
    
    if not os.path.exists(filename):
        print(f"Error: File '{filename}' not found!")
        return

    plots_dir = ensure_plots_dir()

    data = np.loadtxt(filename)
    if data.ndim == 1:
        data = data.reshape(1, -1)

    nx = data[:, 0]
    l2err = data[:, 2]

    fig, ax = plt.subplots(figsize=(8, 6))
    ax.loglog(nx, l2err, marker='o', linewidth=1.5)
    ax.set_xlabel('Grid points (nx)')
    ax.set_ylabel('L2 Error Norm')
    ax.set_title('Error Norm vs Grid Points (log Scale)')
    ax.grid(True, which='both', linestyle='--', alpha=0.4)
    
    plt.tight_layout()
    save_figure(fig, 'error_norm_log_log.png', plots_dir, dpi=150)


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
            ax_y.plot(d['y1d'], d['T2d'][ix, :], linewidth=2.5,
                    label=f"Grid size {d['nx']}x{d['ny']}")

    ix_ref = int(np.argmin(np.abs(ref['x1d'] - x_target)))
    ax_y.plot(ref['y1d'], ref['Tex2d'][ix_ref, :], '--', linewidth=3.0, label='Exact (Tex)',color = "#000000")
    ax_y.set_xlabel('y (m)')
    ax_y.set_ylabel('Temperature (k)')

    x_tag = str(x_target).replace('.', 'p')
    ax_y.set_title(f'Temperature vs y at x = {x_target:.1f}')
    ax_y.grid(True, alpha=0.3)
    ax_y.legend()
    fig_y.tight_layout()
    save_figure(fig_y, f'line_T_vs_y_x{x_tag}.png', plots_dir)

    # -------- Plot 2: T vs x at y = y_target --------
    fig_x, ax_x = plt.subplots(figsize=(9, 6))

    for d in datasets:
        if d != ref:
            jy = int(np.argmin(np.abs(d['y1d'] - y_target)))
            ax_x.plot(d['x1d'], d['T2d'][:, jy], linewidth=2.5,
                    label=f"Grid size {d['nx']}x{d['ny']}")

    jy_ref = int(np.argmin(np.abs(ref['y1d'] - y_target)))
    ax_x.plot(ref['x1d'], ref['Tex2d'][:, jy_ref], '--', linewidth=3.0, label='Exact (Tex)',color = "#000000")
    ax_x.set_xlabel('x (m)')
    ax_x.set_ylabel('Temperature (k)')
    y_tag = str(y_target).replace('.', 'p')
    ax_x.set_title(f'Temperature vs x at y = {y_target:.1f}')
    ax_x.grid(True, alpha=0.3)
    ax_x.legend()
    fig_x.tight_layout()
    save_figure(fig_x, f'line_T_vs_x_y{y_tag}.png', plots_dir)


if __name__ == "__main__":
    grid_sizes = [5, 9, 17, 65, 257, 513]
    x_target = 0.9
    y_target = 0.4
    dat_file = "output/T_xy_513_257_0000.dat"
    error_file = "output/error.dat"
    plot_contour(dat_file)
    plot_error(error_file)

    plot_line_profiles_with_analytical_solution(grid_sizes, x_target=x_target, y_target=y_target)