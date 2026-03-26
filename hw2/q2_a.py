import numpy as np
import matplotlib.pyplot as plt
import os

def save_figure(fig, filename, plots_dir='plots', dpi=300):
    """Helper function to save figures with consistent settings."""
    base_dir = os.path.dirname(os.path.abspath(__file__))
    target_dir = os.path.join(base_dir, plots_dir)
    os.makedirs(target_dir, exist_ok=True)
    filepath = os.path.join(target_dir, filename)
    fig.savefig(filepath, dpi=dpi, bbox_inches='tight')
    plt.close(fig)
    print(f"Plot saved as '{filepath}'")

kdx = np.linspace(0, 2 * np.pi, 500, endpoint=False)
r_values = [0.4, 0.5, 0.6, 0.7]

fig, ax = plt.subplots(figsize=(8, 5))

for r in r_values:
    sigma_dt = np.log((1 - 4 * r * np.sin(kdx / 2)**2))
    ax.plot(kdx, sigma_dt, linewidth=2, label=f'r = {r}')

ax.axhline(0, color='black', linestyle='--', linewidth=1.5, label='σ = 0')

ax.set_xlabel('k∆x', fontsize=14)
ax.set_ylabel('σ∆t', fontsize=14)
ax.set_title('Von Neumann stability: growth rate vs wavenumber', fontsize=14)
ax.set_xticks([0, np.pi/2, np.pi, 3*np.pi/2, 2*np.pi])
ax.set_xticklabels(['0', 'π/2', 'π', '3π/2', '2π'], fontsize=12)
ax.set_ylim([-5,1]) # Set lower limit for better visibility
ax.legend(fontsize=12)
ax.grid(True, alpha=0.3)

plt.tight_layout()
save_figure(plt.gcf(), 'stability_plot.png')