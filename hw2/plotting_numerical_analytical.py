import glob
import os
import re
import matplotlib.pyplot as plt
import numpy as np

FIG_W, FIG_H = 16, 6
LABEL_FS = 16
TITLE_FS = 16
TICK_FS = 13
C_BAR_FS = 13
CONTOUR_LW = 1.4

os.makedirs("plots", exist_ok=True)


def get_finest_file(pattern="output1/T_xy_*.dat"):
    files = glob.glob(pattern)
    if not files:
        raise FileNotFoundError(f"No files found with pattern: {pattern}")

    def key(path):
        m = re.search(r"T_xy_(\d+)_(\d+)_(\d+)\.dat$", os.path.basename(path))
        if m is None:
            return (-1, -1, -1)
        nx, ny, it = map(int, m.groups())
        return (nx * ny, nx, it)

    return max(files, key=key)


def parse_shape(path):
    m = re.search(r"T_xy_(\d+)_(\d+)_\d+\.dat$", os.path.basename(path))
    if m is None:
        raise ValueError(f"Could not parse nx, ny from filename: {path}")
    return int(m.group(1)), int(m.group(2))


def main():
    fpath = get_finest_file()
    nx, ny = parse_shape(fpath)

    data = np.loadtxt(fpath)
    x = data[:, 0].reshape(nx, ny)
    y = data[:, 1].reshape(nx, ny)
    t_num = data[:, 2].reshape(nx, ny)
    t_ex = data[:, 3].reshape(nx, ny)

    vmin = min(t_num.min(), t_ex.min())
    vmax = max(t_num.max(), t_ex.max())
    levels = np.linspace(vmin, vmax, 40)

    fig, axes = plt.subplots(1, 2, figsize=(FIG_W, FIG_H), constrained_layout=True)

    c0 = axes[0].contourf(x, y, t_num, levels=levels, cmap="viridis")
    axes[0].contour(x, y, t_num, levels=levels[::4], colors="white", linewidths=CONTOUR_LW, alpha=0.45)
    axes[0].set_title(f"Numerical (finest grid: {nx}x{ny})", fontsize=TITLE_FS)
    axes[0].set_xlabel("x", fontsize=LABEL_FS)
    axes[0].set_ylabel("y", fontsize=LABEL_FS)
    axes[0].tick_params(labelsize=TICK_FS)

    c1 = axes[1].contourf(x, y, t_ex, levels=levels, cmap="viridis")
    axes[1].contour(x, y, t_ex, levels=levels[::4], colors="white", linewidths=CONTOUR_LW, alpha=0.45)
    axes[1].set_title("Analytical", fontsize=TITLE_FS)
    axes[1].set_xlabel("x", fontsize=LABEL_FS)
    axes[1].set_ylabel("y", fontsize=LABEL_FS)
    axes[1].tick_params(labelsize=TICK_FS)

    cb0 = fig.colorbar(c0, ax=axes[0], shrink=0.96)
    cb1 = fig.colorbar(c1, ax=axes[1], shrink=0.96)
    cb0.set_label("Temperature", fontsize=C_BAR_FS)
    cb1.set_label("Temperature", fontsize=C_BAR_FS)
    cb0.ax.tick_params(labelsize=TICK_FS)
    cb1.ax.tick_params(labelsize=TICK_FS)

    out = "plots/q1_numerical_vs_analytical.png"
    fig.savefig(out, dpi=180, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved {out} from {fpath}")


if __name__ == "__main__":
    main()
