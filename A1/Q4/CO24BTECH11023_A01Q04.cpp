#include <bits/stdc++.h>
using namespace std;

#define nu 1.82 // Courant number for Roll 23 (Problem 4)

vector<double> makeGrid(double xstart, double xend, int Nx)
{
    vector<double> grid;
    double dx = (xend - xstart) / Nx;
    for (int i = 0; i < Nx; i++)
    {
        grid.push_back(xstart + i * dx);
    }
    return grid;
}

vector<double> exactSol(const vector<double> &x, double t, double c, double L)
{
    vector<double> sol;
    for (double xi : x)
    {
        sol.push_back(sin(2 * M_PI * (xi - c * t) / L));
    }
    return sol;
}

void setinitialCondition(vector<double> &u, const vector<double> &x, double L)
{
    for (size_t i = 0; i < x.size(); i++)
    {
        u[i] = sin(2 * M_PI * x[i] / L);
    }
}

double calculatel2Norm(const vector<double> &u_num, const vector<double> &u_exact, double dx)
{
    double sum = 0.0;
    for (size_t i = 0; i < u_num.size(); ++i)
    {
        double diff = u_num[i] - u_exact[i];
        sum += diff * diff;
    }
    return sqrt(dx * sum);
}

void writesolutiontoFile(const vector<double> &u, const vector<double> &x, double t)
{
    ostringstream name;
    name << "data/LaxW_T_" << t << ".txt";
    ofstream out(name.str());
    for (size_t i = 0; i < u.size() && i < x.size(); ++i)
    {
        out << x[i] << " " << u[i] << "\n";
    }
}

void solveLaxWendroff(vector<double> &u, const vector<double> &x, double dt, int Nt, int Nx, const vector<double> &output_times)
{
    vector<double> u_new(u.size(), 0.0);

    // Write initial condition
    writesolutiontoFile(u, x, 0.0);

    int next_output_idx = 1;

    for (int n = 1; n <= Nt; n++)
    {
        // Interior points using Lax-Wendroff
        for (int i = 1; i < Nx - 1; ++i)
        {
            u_new[i] = u[i] - (nu / 2.0) * (u[i + 1] - u[i - 1]) 
                       + (nu * nu / 2.0) * (u[i + 1] - 2.0 * u[i] + u[i - 1]);
        }
        
        // Periodic Boundary Conditions
        u_new[0] = u[0] - (nu / 2.0) * (u[1] - u[Nx - 1]) 
                   + (nu * nu / 2.0) * (u[1] - 2.0 * u[0] + u[Nx - 1]);
                   
        u_new[Nx - 1] = u[Nx - 1] - (nu / 2.0) * (u[0] - u[Nx - 2]) 
                        + (nu * nu / 2.0) * (u[0] - 2.0 * u[Nx - 1] + u[Nx - 2]);
        
        // Update state
        u = u_new;

        // Output logic
        double current_time = n * dt;
        if (next_output_idx < output_times.size() && current_time >= output_times[next_output_idx])
        {
            writesolutiontoFile(u, x, output_times[next_output_idx]);
            next_output_idx++;
        }
    }
}

int main()
{
    // Extract inputs from file
    FILE *fp = fopen("input.txt", "r");
    if (!fp) {
        printf("Error: Could not open input.txt\n");
        return 1;
    }
    int Nx, c;
    double xstart, xend, tstart, tend;
    fscanf(fp, "%d %d", &Nx, &c);
    fscanf(fp, "%lf %lf %lf %lf", &xstart, &xend, &tstart, &tend);
    fclose(fp);
    
    printf("The inputs are Nx = %d, c = %d, xstart = %lf, xend = %lf, tstart = %lf, tend = %lf\n\n", Nx, c, xstart, xend, tstart, tend);

    // Build grid and calculate time step
    double dx = (xend - xstart) / Nx;
    double dt = nu * dx / c;
    int Nt = ceil((tend - tstart) / dt);
    double L = xend - xstart;
    vector<double> x = makeGrid(xstart, xend, Nx);

    printf("The values of dx and dt are %lf and %lf respectively\n\n", dx, dt);
    printf("The number of time steps Nt is %d\n\n", Nt);

    // Initialize state and output schedule
    vector<double> u(Nx, 0.0);
    setinitialCondition(u, x, L);
    vector<double> output_times = {0.0, 0.25, 0.5, 0.75, 1.0};
    
    // Solve using Lax-Wendroff
    printf("Solving using Lax-Wendroff...\n\n");
    solveLaxWendroff(u, x, dt, Nt, Nx, output_times);
    printf("Lax-Wendroff solution completed.\n\n");
    
    // Calculate final error
    double final_time = Nt * dt;
    double e = calculatel2Norm(u, exactSol(x, final_time, c, L), dx);
    printf("L2 Error for Lax-Wendroff at t=1.0: %lf\n\n", e);

    printf("Program completed.\n");
    return 0;
}