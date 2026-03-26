import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

beta_values = [1.5, 2.0, 5.0, 10.0, 50.0]
nx          = 33
ny          = 17
Lx          = 1.0
Ly          = 0.5

OUTDIR = Path(__file__).resolve().parent / 'plots'
OUTDIR.mkdir(parents=True, exist_ok=True)


def save_plot(filename):
    outpath = OUTDIR / filename
    plt.savefig(outpath, dpi=150, bbox_inches='tight')
    print(f"Saved {outpath}")

def compute_xs(xu, L, beta):
    Lambda = (1.0 + beta) / (beta - 1.0)
    R      = Lambda ** (2.0 * (xu / L) - 1.0)
    return L * ((1.0 + beta) * R - beta + 1.0) / (2.0 * (1.0 + R))

def make_grid(n, L, beta):
    xu   = np.linspace(0.0, 1.0, n)       # uniform in [0,1]
    # h = L / (n - 1)                       # uniform grid spacing in [0, L]  
    xs   = compute_xs(xu, L, beta)        # stretched physical coords on [0, L]
    dx   = np.diff(xs)
    return xu, xs, dx


def make_uniform_grid(n, L):
    xs = np.linspace(0.0, L, n)
    dx = np.diff(xs)
    return xs, dx



fig, axes = plt.subplots(1, 3, figsize=(17, 5))
compare_betas = [None,1.5 ,2.0, 10.0]  # None = uniform

for ax, beta in zip(axes, compare_betas):
    if beta is None:
        xs_x, _ = make_uniform_grid(nx, Lx)
        xs_y, _ = make_uniform_grid(ny, Ly)
        title = f'2D uniform grid ({nx}x{ny})'
    else:
        _, xs_x, _ = make_grid(nx, Lx, beta)
        _, xs_y, _ = make_grid(ny, Ly, beta)
        title = f'2D stretched grid beta={beta:g} ({nx}x{ny})'

    for xi in xs_x:
        ax.axvline(xi, color='steelblue', linewidth=0.7, alpha=0.8)
    for yi in xs_y:
        ax.axhline(yi, color='steelblue', linewidth=0.7, alpha=0.8)

    xx, yy = np.meshgrid(xs_x, xs_y)
    ax.plot(xx, yy, 'k.', markersize=2)

    ax.set_xlim(-0.02, Lx + 0.02)
    ax.set_ylim(-0.01, Ly + 0.01)
    ax.set_xlabel('x', fontsize=12)
    ax.set_ylabel('y', fontsize=12)
    ax.set_title(title, fontsize=12)
    ax.set_aspect('equal')
    ax.grid(False)

plt.tight_layout()
save_plot('grid_2d.png')
plt.close()


fig, axes = plt.subplots(1, 2, figsize=(13, 5))

for beta in beta_values:
    _, _, dx_x = make_grid(nx, Lx, beta)
    _, _, dx_y = make_grid(ny, Ly, beta)
    axes[0].plot(np.arange(dx_x.size), dx_x, linewidth=1.8, label=f'beta={beta:g}')
    axes[1].plot(np.arange(dx_y.size), dx_y, linewidth=1.8, label=f'beta={beta:g}')

axes[0].set_title('X grid spacing for all beta')
axes[1].set_title('Y grid spacing for all beta')
axes[0].set_xlabel('Cell index i')
axes[1].set_xlabel('Cell index j')
axes[0].set_ylabel('dx')
axes[1].set_ylabel('dy')
axes[0].set_yscale('log')
axes[1].set_yscale('log')
axes[0].grid(True, linestyle='--', alpha=0.4)
axes[1].grid(True, linestyle='--', alpha=0.4)
axes[0].legend(fontsize=9)
axes[1].legend(fontsize=9)

plt.tight_layout()
save_plot('grid_spacings_all.png')
plt.close()

