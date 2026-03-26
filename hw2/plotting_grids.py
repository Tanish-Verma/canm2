import numpy as np
import matplotlib.pyplot as plt
import os

NX = 33
NY = 17
LX = 1.0
LY = 0.5
BETA_VALUES = [1.5, 2.0, 5.0, 10.0, 50.0]

os.makedirs("plots", exist_ok=True)

def save_plot(fig, filename):
    fig.savefig(f"plots/{filename}", dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved plots/{filename}")


def stretched_coords(n, length, beta):
    xu = np.linspace(0.0, 1.0, n)
    lam = (1.0 + beta) / (beta - 1.0)
    ratio = lam ** (2.0 * (xu / length) - 1.0)
    xs = length * ((1.0 + beta) * ratio - beta + 1.0) / (2.0 * (1.0 + ratio))
    return xs


def uniform_coords(n, length):
    return np.linspace(0.0, length, n)


def make_2d_grid_plot():
    fig, axes = plt.subplots(1, 3, figsize=(17, 5))
    compare_betas = [None, 1.5, 10.0]

    for ax, beta in zip(axes, compare_betas):
        if beta is None:
            xs_x = uniform_coords(NX, LX)
            xs_y = uniform_coords(NY, LY)
            title = f"2D uniform grid ({NX}x{NY})"
        else:
            xs_x = stretched_coords(NX, LX, beta)
            xs_y = stretched_coords(NY, LY, beta)
            title = f"2D stretched grid beta={beta:g} ({NX}x{NY})"

        for xi in xs_x:
            ax.axvline(xi, color="steelblue", linewidth=0.7, alpha=0.8)
        for yi in xs_y:
            ax.axhline(yi, color="steelblue", linewidth=0.7, alpha=0.8)

        xx, yy = np.meshgrid(xs_x, xs_y)
        ax.plot(xx, yy, "k.", markersize=2)
        ax.set_xlim(-0.02, LX + 0.02)
        ax.set_ylim(-0.01, LY + 0.01)
        ax.set_xlabel("x", fontsize=12)
        ax.set_ylabel("y", fontsize=12)
        ax.set_title(title, fontsize=12)
        ax.set_aspect("equal")
        ax.grid(False)

    fig.tight_layout()
    save_plot(fig, "grid_2d.png")


def make_spacing_plot():
    fig, axes = plt.subplots(1, 2, figsize=(13, 5))

    for beta in BETA_VALUES:
        dx_x = np.diff(stretched_coords(NX, LX, beta))
        dx_y = np.diff(stretched_coords(NY, LY, beta))
        axes[0].plot(np.arange(dx_x.size), dx_x, linewidth=1.8, label=f"beta={beta:g}")
        axes[1].plot(np.arange(dx_y.size), dx_y, linewidth=1.8, label=f"beta={beta:g}")

    axes[0].set_title("X grid spacing for all beta")
    axes[1].set_title("Y grid spacing for all beta")
    axes[0].set_xlabel("Cell index i")
    axes[1].set_xlabel("Cell index j")
    axes[0].set_ylabel("dx")
    axes[1].set_ylabel("dy")
    axes[0].set_yscale("log")
    axes[1].set_yscale("log")
    axes[0].grid(True, linestyle="--", alpha=0.4)
    axes[1].grid(True, linestyle="--", alpha=0.4)
    axes[0].legend(fontsize=9)
    axes[1].legend(fontsize=9)

    fig.tight_layout()
    save_plot(fig, "grid_spacings_all.png")


def main():
    make_2d_grid_plot()
    make_spacing_plot()


if __name__ == "__main__":
    main()

