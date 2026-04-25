#include <bits/stdc++.h>
using namespace std;

// ─── Parameters (Roll 23, group 13-25) ───────────────────────────────────────
// Nu = 0.45, Nx = Ny = 100, IC = sin(2*pi*x)*sin(2*pi*y)
// T = 1s, output times = 0, 0.25, 0.5, 0.75, 1.0
// ─────────────────────────────────────────────────────────────────────────────

vector<double> makeGrid(double xstart, double xend, int N)
{
    vector<double> g(N);
    double dx = (xend - xstart) / N;
    for (int i = 0; i < N; i++)
        g[i] = xstart + i * dx;
    return g;
}

double exactSol(double x, double y, double t, double c)
{
    return sin(2 * M_PI * (x - c * t)) * sin(2 * M_PI * (y - c * t));
}

void setIC(vector<vector<double>> &u, const vector<double> &x, const vector<double> &y, int Nx, int Ny)
{
    for (int i = 0; i < Nx; i++)
        for (int j = 0; j < Ny; j++)
            u[i][j] = sin(2 * M_PI * x[i]) * sin(2 * M_PI * y[j]);
}

double computeL2Error(const vector<vector<double>> &u,
                      const vector<double> &x, const vector<double> &y,
                      double t, double c, double dx, double dy, int Nx, int Ny)
{
    double sum = 0.0;
    for (int i = 0; i < Nx; i++)
        for (int j = 0; j < Ny; j++)
        {
            double diff = u[i][j] - exactSol(x[i], y[j], t, c);
            sum += diff * diff;
        }
    return sqrt(dx * dy * sum);
}

// Write snapshot. Filename uses actual simulation time t_actual
void writeSolution(const vector<vector<double>> &u,
                   const vector<double> &x, const vector<double> &y,
                   int Nx, int Ny, double t_actual, const string &scheme)
{
    ostringstream fname;
    fname << "data/" << scheme << "_T_" << fixed << setprecision(6) << t_actual << ".txt";
    ofstream out(fname.str());
    out << fixed << setprecision(8);
    for (int i = 0; i < Nx; i++)
        for (int j = 0; j < Ny; j++)
            out << x[i] << " " << y[j] << " " << u[i][j] << "\n";
}

// ─── Scheme 1A: 2D Upwind (FTBS-FTBS) for c > 0 ─────────────────────────────
void solveUpwind(vector<vector<double>> &u,
                 const vector<double> &x, const vector<double> &y,
                 double dt, double dx, double dy,
                 double c, double nu, int Nx, int Ny, int Nt,
                 const vector<double> &output_times)
{
    writeSolution(u, x, y, Nx, Ny, 0.0, "Upwind");

    vector<vector<double>> u_new(Nx, vector<double>(Ny, 0.0));
    int next_out = 1;

    for (int n = 1; n <= Nt; n++)
    {
        for (int i = 0; i < Nx; i++)
        {
            int im1 = (i - 1 + Nx) % Nx;
            for (int j = 0; j < Ny; j++)
            {
                int jm1 = (j - 1 + Ny) % Ny;
                u_new[i][j] = u[i][j]
                              - nu * (u[i][j] - u[im1][j])
                              - nu * (u[i][j] - u[i][jm1]);
            }
        }
        u = u_new;

        // Use actual n*dt as the timestamp — this is what Python will read back
        double t_now = n * dt;
        if (next_out < (int)output_times.size() && t_now >= output_times[next_out] - 1e-12)
        {
            writeSolution(u, x, y, Nx, Ny, t_now, "Upwind");
            next_out++;
        }
    }
    printf("Upwind done.\n");
}

// ─── Scheme 1B: 2D Lax-Friedrichs ───────────────────────────────────────────
void solveLaxFriedrichs(vector<vector<double>> &u,
                        const vector<double> &x, const vector<double> &y,
                        double dt, double dx, double dy,
                        double c, double nu, int Nx, int Ny, int Nt,
                        const vector<double> &output_times)
{
    writeSolution(u, x, y, Nx, Ny, 0.0, "LaxF");

    vector<vector<double>> u_new(Nx, vector<double>(Ny, 0.0));
    int next_out = 1;

    for (int n = 1; n <= Nt; n++)
    {
        for (int i = 0; i < Nx; i++)
        {
            int ip1 = (i + 1) % Nx;
            int im1 = (i - 1 + Nx) % Nx;
            for (int j = 0; j < Ny; j++)
            {
                int jp1 = (j + 1) % Ny;
                int jm1 = (j - 1 + Ny) % Ny;
                u_new[i][j] = 0.25 * (u[ip1][j] + u[im1][j] + u[i][jp1] + u[i][jm1])
                              - (nu / 2.0) * (u[ip1][j] - u[im1][j])
                              - (nu / 2.0) * (u[i][jp1] - u[i][jm1]);
            }
        }
        u = u_new;

        // Divergence check
        double max_u = 0.0;
        for (int i = 0; i < Nx; i++)
            for (int j = 0; j < Ny; j++)
                max_u = max(max_u, fabs(u[i][j]));
        if (max_u > 1e6)
        {
            printf("LaxF diverged at step %d (t=%.6f). Stopping.\n", n, n * dt);
            return;
        }

        double t_now = n * dt;
        if (next_out < (int)output_times.size() && t_now >= output_times[next_out] - 1e-12)
        {
            writeSolution(u, x, y, Nx, Ny, t_now, "LaxF");
            next_out++;
        }
    }
    printf("LaxF done.\n");
}

int main()
{
    const int Nx = 100, Ny = 100;
    const double xstart = 0.0, xend = 1.0;
    const double ystart = 0.0, yend = 1.0;
    const double c  = 1.0;
    const double nu = 0.45; // Roll 23
    const double T  = 1.0;

    double dx = (xend - xstart) / Nx;
    double dy = (yend - ystart) / Ny;
    double dt = nu * dx / c;
    int    Nt = (int)ceil(T / dt);
    double final_t = Nt * dt;

    printf("Roll 23 — 2D Linear Advection\n");
    printf("Nx=%d, Ny=%d, dx=%.6f, dy=%.6f\n", Nx, Ny, dx, dy);
    printf("nu=%.4f, dt=%.6f, Nt=%d, final_t=%.6f\n\n", nu, dt, Nt, final_t);

    vector<double> x = makeGrid(xstart, xend, Nx);
    vector<double> y = makeGrid(ystart, yend, Ny);
    vector<double> output_times = {0.0, 0.25, 0.5, 0.75, 1.0};

    // ── Upwind ──
    vector<vector<double>> u(Nx, vector<double>(Ny, 0.0));
    setIC(u, x, y, Nx, Ny);
    printf("Solving Scheme 1A: 2D Upwind...\n");
    solveUpwind(u, x, y, dt, dx, dy, c, nu, Nx, Ny, Nt, output_times);
    double e_upwind = computeL2Error(u, x, y, final_t, c, dx, dy, Nx, Ny);
    printf("Upwind L2 error at T=%.6f: %.6e\n\n", final_t, e_upwind);

    // ── Lax-Friedrichs ──
    setIC(u, x, y, Nx, Ny);
    printf("Solving Scheme 1B: 2D Lax-Friedrichs...\n");
    solveLaxFriedrichs(u, x, y, dt, dx, dy, c, nu, Nx, Ny, Nt, output_times);
    double e_laxf = computeL2Error(u, x, y, final_t, c, dx, dy, Nx, Ny);
    printf("LaxF L2 error at T=%.6f: %.6e\n\n", final_t, e_laxf);

    printf("Program completed.\n");
    return 0;
}