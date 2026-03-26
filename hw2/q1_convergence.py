import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
import re

BASE_DIR = Path(__file__).resolve().parent
ERROR_DIR = BASE_DIR / "output1"
PLOTS_DIR = BASE_DIR / "plots"
PLOTS_DIR.mkdir(parents=True, exist_ok=True)

# -------------------------------------------------------
# Read error data
# Each file has lines: nx ny l2err
# Auto-detects files matching output1/error_b*.dat
# -------------------------------------------------------

def load_error_file(filename):
    data = np.loadtxt(filename)
    if data.ndim == 1:
        data = data.reshape(1, -1)
    nx  = data[:, 0]
    err = data[:, 2]
    return nx, err

def make_slope2_ref(nx, err):
    # anchor reference line to first data point
    nx_ref  = np.array([nx[0], nx[-1]])
    err_ref = err[0] * (nx[0] / nx_ref) ** 2
    return nx_ref, err_ref

def plot_convergence(nx, err, beta_str, outfile):
    nx_ref, err_ref = make_slope2_ref(nx, err)

    fig, ax = plt.subplots(figsize=(7, 5))
    ax.loglog(nx, err, 'o-', linewidth=2, markersize=7, label=rf'$\beta = {beta_str}$')
    ax.loglog(nx_ref, err_ref, 'k:', linewidth=1.5, label='slope $-2$')

    ax.set_xlabel('$n_x$', fontsize=14)
    ax.set_ylabel('$L_2$ error', fontsize=14)
    ax.set_title(rf'Order of accuracy — stretched grid, $\beta = {beta_str}$ (Q1)', fontsize=13)
    ax.legend(fontsize=12)
    ax.grid(True, which='both', alpha=0.3)
    ax.tick_params(labelsize=12)

    plt.tight_layout()
    outpath = PLOTS_DIR / outfile
    plt.savefig(outpath, dpi=150)
    plt.show()
    print(f"Saved {outpath}")


def parse_beta(file_path):
    match = re.match(r"error_b(\d+)\.dat$", file_path.name)
    if not match:
        return None, None
    beta_tag = match.group(1)
    beta_val = int(beta_tag) / 100.0
    return beta_tag, f"{beta_val:g}"


error_files = sorted(ERROR_DIR.glob("error_b*.dat"))
if not error_files:
    raise FileNotFoundError(f"No files found matching {ERROR_DIR / 'error_b*.dat'}")

for error_file in error_files:
    beta_tag, beta_str = parse_beta(error_file)
    if beta_tag is None:
        continue

    nx, err = load_error_file(error_file)
    plot_convergence(nx, err, beta_str, f"q1_convergence_b{beta_tag}.png")