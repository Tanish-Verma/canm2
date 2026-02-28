import numpy as np 
import matplotlib.pyplot as plt 

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
    plt.figure(figsize=(10, 8))
    plt.contourf(X, Y, u_res, levels=50, cmap='inferno')
    cbar = plt.colorbar(label='Temperature (K)')
    cbar.set_label('Temperature (K)', fontsize=14)
    plt.title('Temperature Distribution in the Plate', fontsize=16, fontweight='bold')
    plt.xlabel('x (m)', fontsize=14)
    plt.ylabel('y (m)', fontsize=14)
    plt.tick_params(labelsize=12)
    plt.savefig("plots/contour_analytical.png", dpi=300, bbox_inches='tight')
    plt.close()

def plot_yprofile(x_val,y,a,b,N,T0):
    plt.figure(figsize=(12, 8))
    for x_i in x_val:
        u_xi = u(x_i,y,a,b,N,T0)
        plt.plot(y, u_xi, label=f'x={x_i:.1f} m', linewidth=2.5)
    plt.title('Temperature Profile Along y-axis for Different x Values', fontsize=16, fontweight='bold')
    plt.xlabel('y (m)', fontsize=14)
    plt.grid(alpha=0.3, linewidth=1)
    plt.ylabel('Temperature (K)', fontsize=14)
    plt.legend(fontsize=12, loc='best')
    plt.tick_params(labelsize=12)
    plt.savefig("plots/y_profile.png", dpi=300, bbox_inches='tight')
    plt.close()

def plot_xprofile(y_val,x,a,b,N,T0):
    plt.figure(figsize=(12, 8))
    for y_i in y_val:
        u_yi = u(x,y_i,a,b,N,T0)
        plt.plot(x, u_yi, label=f'y={y_i:.1f} m', linewidth=2.5)
    plt.title('Temperature Profile Along x-axis for Different y Values', fontsize=16, fontweight='bold')
    plt.xlabel('x (m)', fontsize=14)
    plt.grid(alpha=0.3, linewidth=1)
    plt.ylabel('Temperature (K)', fontsize=14)
    plt.legend(fontsize=12, loc='best')
    plt.tick_params(labelsize=12)
    plt.savefig("plots/x_profile.png", dpi=300, bbox_inches='tight')
    plt.close()

N = 50
a = 1
b = 0.5
T0 = 300
x = np.linspace(0, a, 5)
y = np.linspace(0, b, 3)
X, Y = np.meshgrid(x, y)
u_res = u(X, Y, a, b, N, T0)

#debugging print statements
# for i in range(len(x)):
#     for j in range(len(y)):
#         print(f"x: {x[i]:.2f} m, y: {y[j]:.2f} m, u: {u_res[j,i]:.2f} K")

Plot_contour(X,Y,u_res)
x_val = a*np.array([0.1,0.3,0.5,0.7,0.9],dtype=float)
y_val = b*np.array([0.1,0.5,0.9],dtype=float)

plot_yprofile(x_val,y,a,b,N,T0)
plot_xprofile(y_val,x,a,b,N,T0)