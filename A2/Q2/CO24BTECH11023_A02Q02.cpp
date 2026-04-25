#include <bits/stdc++.h>
using namespace std;

// Roll 23: A = 0.75
// Nx = 200, dx = 1/200
// dt = 0.8 * dx / A  (fixed)
// t* = 1 / (2*pi*A)
// T  = 1.5 * t*

vector<double> makeGrid(double xstart, double xend, int Nx)
{
    vector<double> g(Nx);
    double dx = (xend - xstart) / Nx;
    for (int i = 0; i < Nx; i++)
        g[i] = xstart + i * dx;
    return g;
}

void setIC(vector<double> &u, const vector<double> &x, double A, int Nx)
{
    for (int i = 0; i < Nx; i++)
        u[i] = A * sin(2 * M_PI * x[i]);
}

void writeSolution(const vector<double> &u, const vector<double> &x, double t_actual, int Nx)
{
    ostringstream fname;
    fname << "data/Burgers_T_" << fixed << setprecision(6) << t_actual << ".txt";
    ofstream out(fname.str());
    out << fixed << setprecision(8);
    for (int i = 0; i < Nx; i++)
        out << x[i] << " " << u[i] << "\n";
}

void solveBurgers(vector<double> &u, const vector<double> &x,
                  double dt, double dx, int Nx, int Nt,
                  const vector<double> &output_times)
{
    // Write IC (t=0)
    writeSolution(u, x, 0.0, Nx);

    vector<double> u_new(Nx, 0.0);
    int next_out = 1;

    for (int n = 1; n <= Nt; n++)
    {
        for (int i = 0; i < Nx; i++)
        {
            int ip1 = (i + 1) % Nx;
            int im1 = (i - 1 + Nx) % Nx;

            // Lax-Friedrichs conservative form for Burgers
            u_new[i] = 0.5 * (u[ip1] + u[im1])
                       - (dt / (2.0 * dx)) * (0.5 * u[ip1]*u[ip1] - 0.5 * u[im1]*u[im1]);
        }
        u = u_new;

        double t_now = n * dt;
        if (next_out < (int)output_times.size() && t_now >= output_times[next_out])
        {
            writeSolution(u, x, t_now, Nx);
            next_out++;
        }
    }
    printf("Burgers solve done.\n");
}

int main()
{
    const int    Nx     = 200;
    const double xstart = 0.0, xend = 1.0;
    const double A      = 0.75; // Roll 23

    double dx  = (xend - xstart) / Nx;
    double dt  = 0.8 * dx / A;
    double t_star = 1.0 / (2.0 * M_PI * A);
    double T      = 1.5 * t_star;
    int    Nt     = (int)ceil(T / dt);
    double final_t = Nt * dt;

    printf("Roll 23 — Burgers Equation (Lax-Friedrichs)\n");
    printf("A=%.4f, Nx=%d, dx=%.6f, dt=%.6f\n", A, Nx, dx, dt);
    printf("t* = %.6f s\n", t_star);
    printf("T  = 1.5 * t* = %.6f s\n", T);
    printf("Nt=%d, final_t=%.6f\n\n", Nt, final_t);

    // Verify CFL: dt * max|u0| / dx = 0.8 * dx/A * A / dx = 0.8 <= 1
    printf("CFL check: dt * max|u0| / dx = %.4f (should be <= 1)\n\n", dt * A / dx);

    vector<double> x = makeGrid(xstart, xend, Nx);
    vector<double> u(Nx, 0.0);
    setIC(u, x, A, Nx);

    // Output at t=0, t*/2, t*, T
    vector<double> output_times = {0.0, t_star / 2.0, t_star, T};

    solveBurgers(u, x, dt, dx, Nx, Nt, output_times);

    printf("Program completed.\n");
    return 0;
}