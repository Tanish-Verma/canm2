import numpy as np
import matplotlib.pyplot as plt
import os

os.makedirs("plots", exist_ok=True)

def save_figure(fig, filename, dpi=300):
    filepath = f"plots/{filename}"
    fig.savefig(filepath, dpi=dpi, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved {filepath}")


def main():
    kdx = np.linspace(0, 2 * np.pi, 500, endpoint=False)
    r_values = [0.4, 0.5, 0.6, 0.7]

    fig, ax = plt.subplots(figsize=(8, 5))

    for r in r_values:
        growth_factor = 1 - 4 * r * np.sin(kdx / 2) ** 2
        sigma_dt = np.full_like(growth_factor, np.nan, dtype=float)
        valid = growth_factor > 0
        sigma_dt[valid] = np.log(growth_factor[valid])
        ax.plot(kdx, sigma_dt, linewidth=2, label=f"r = {r}")

    ax.axhline(0, color="black", linestyle="--", linewidth=1.5, label="sigma = 0")
    ax.set_xlabel("k delta x", fontsize=14)
    ax.set_ylabel("sigma delta t", fontsize=14)
    ax.set_title("Von Neumann stability: growth rate vs wavenumber", fontsize=14)
    ax.set_xticks([0, np.pi / 2, np.pi, 3 * np.pi / 2, 2 * np.pi])
    ax.set_xticklabels(["0", "pi/2", "pi", "3pi/2", "2pi"], fontsize=12)
    ax.set_ylim([-5, 1])
    ax.legend(fontsize=12)
    ax.grid(True, alpha=0.3)

    fig.tight_layout()
    save_figure(fig, "stability_plot.png")


if __name__ == "__main__":
    main()