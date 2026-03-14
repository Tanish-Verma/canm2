import numpy as np 
import matplotlib.pyplot as plt 
import os

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

def cos_term(y,m,b):
    return np.cos((2*m+1)*np.pi*y/b)

def sinh_term(x,m,a,b):
    return np.sinh((2*m+1)*np.pi*x/b)/np.sinh((2*m+1)*np.pi*a/b)

def u(x,y,a,b,N,T0):
    u_sum = 0
    for m in range(N):
        u_sum += (1/((2*m+1)*np.pi)**2)*cos_term(y,m,b)*sinh_term(x,m,a,b)

    u_result = T0*(1 - 0.8 * u_sum)
    
    return u_result

def Plot_contour(X,Y,u_res):
    plots_dir = ensure_plots_dir()
    fig, ax = plt.subplots(figsize=(10, 8))
    cs = ax.contourf(X, Y, u_res, levels=50, cmap='viridis')
    cbar = fig.colorbar(cs, ax=ax)
    cbar.set_label('Temperature (K)', fontsize=14)
    ax.set_title('Temperature Distribution in the Plate', fontsize=16, fontweight='bold')
    ax.set_xlabel('x (m)', fontsize=14)
    ax.set_ylabel('y (m)', fontsize=14)
    ax.tick_params(labelsize=12)
    fig.tight_layout()
    save_figure(fig, 'contour_analytical.png', plots_dir, dpi=300)

def plot_yprofile(x_val,y,a,b,N,T0):
    plots_dir = ensure_plots_dir()
    fig, ax = plt.subplots(figsize=(9, 6))
    for x_i in x_val:
        u_xi = u(x_i,y,a,b,N,T0)
        ax.plot(y, u_xi, label=f'x={x_i} m', linewidth=3.0)
    ax.set_title('Temperature Profile Along y-axis for Different x Values', fontsize=16, fontweight='bold')
    ax.set_xlabel('y (m)', fontsize=14)
    ax.set_ylabel('Temperature (K)', fontsize=14)
    ax.grid(True, alpha=0.3, linewidth=1)
    ax.legend(fontsize=13, loc='best')
    ax.tick_params(labelsize=12)
    fig.tight_layout()
    save_figure(fig, 'y_profile.png', plots_dir, dpi=300)

def plot_xprofile(y_val,x,a,b,N,T0):
    plots_dir = ensure_plots_dir()
    fig, ax = plt.subplots(figsize=(9, 6))
    for y_i in y_val:
        u_yi = u(x,y_i,a,b,N,T0)
        ax.plot(x, u_yi, label=f'y={y_i} m', linewidth=3.0)
    ax.set_title('Temperature Profile Along x-axis for Different y Values', fontsize=16, fontweight='bold')
    ax.set_xlabel('x (m)', fontsize=14)
    ax.set_ylabel('Temperature (K)', fontsize=14)
    ax.grid(True, alpha=0.3, linewidth=1)
    ax.legend(fontsize=13, loc='best')
    ax.tick_params(labelsize=12)
    fig.tight_layout()
    save_figure(fig, 'x_profile.png', plots_dir, dpi=300)

N = 50
a = 1
b = 0.5
T0 = 300
nx = 513
ny = 257
x = np.linspace(0, a, nx)
y = np.linspace(0, b, ny)
X, Y = np.meshgrid(x, y)
u_res = u(X, Y, a, b, N, T0)



Plot_contour(X,Y,u_res)
x_val = a*np.array([0.1,0.3,0.5,0.7,0.9],dtype=float)
y_val = b*np.array([0.2,0.5,0.8],dtype=float)

plot_yprofile(x_val,y,a,b,N,T0)
plot_xprofile(y_val,x,a,b,N,T0)