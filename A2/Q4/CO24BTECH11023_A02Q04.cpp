#include <bits/stdc++.h>
using namespace std;

// Roll 23: kappa=0.010, nu=0.80, c=1, Nx=100, dx=1/100
// dt = nu*dx/c = 0.008
// r  = kappa*dt/dx^2 = 0.010*0.008/0.0001 = 0.8
// Pe_h = nu/r = 0.80/0.8 = 1.0
// IC: sin(2*pi*x), T=1s, periodic BCs
// Exact: u_ex(x,t) = exp(-4*pi^2*kappa*t) * sin(2*pi*(x-t))

// ─── Grid ─────────────────────────────────────────────────────────────────────
vector<double> makeGrid(double xstart, double xend, int Nx)
{
    vector<double> g(Nx);
    double dx = (xend - xstart) / Nx;
    for (int i = 0; i < Nx; i++)
        g[i] = xstart + i * dx;
    return g;
}

void setIC(vector<double> &u, const vector<double> &x, int Nx)
{
    for (int i = 0; i < Nx; i++)
        u[i] = sin(2 * M_PI * x[i]);
}

double exactSol(double x, double t, double kappa)
{
    return exp(-4.0 * M_PI * M_PI * kappa * t) * sin(2 * M_PI * (x - t));
}

double computeL2Error(const vector<double> &u, const vector<double> &x,
                      double t, double kappa, double dx, int Nx)
{
    double sum = 0.0;
    for (int i = 0; i < Nx; i++)
    {
        double diff = u[i] - exactSol(x[i], t, kappa);
        sum += diff * diff;
    }
    return sqrt(dx * sum);
}

void writeSolution(const vector<double> &u, const vector<double> &x,
                   double t_actual, int Nx, const string &scheme)
{
    ostringstream fname;
    fname << "data/" << scheme << "_T_" << fixed << setprecision(6) << t_actual << ".txt";
    ofstream out(fname.str());
    out << fixed << setprecision(8);
    for (int i = 0; i < Nx; i++)
        out << x[i] << " " << u[i] << "\n";
}

// ─── Periodic TDMA (Sherman-Morrison) ────────────────────────────────────────
// Solves periodic tridiagonal with uniform coefficients a, b, c
// (same a, b, c for every row — valid here since kappa and c are constant)
vector<double> solveTDMA(const vector<double> &a, const vector<double> &b,
                          const vector<double> &c, const vector<double> &d)
{
    int n = d.size();
    vector<double> c_star(n, 0.0), d_star(n, 0.0), sol(n, 0.0);
    c_star[0] = c[0] / b[0];
    d_star[0] = d[0] / b[0];
    for (int i = 1; i < n; i++)
    {
        double m = b[i] - a[i] * c_star[i - 1];
        c_star[i] = (i < n - 1) ? c[i] / m : 0.0;
        d_star[i] = (d[i] - a[i] * d_star[i - 1]) / m;
    }
    sol[n - 1] = d_star[n - 1];
    for (int i = n - 2; i >= 0; i--)
        sol[i] = d_star[i] - c_star[i] * sol[i + 1];
    return sol;
}

vector<double> solvePeriodicTDMA(double a_coeff, double b_coeff, double c_coeff,
                                  const vector<double> &d)
{
    int N = d.size();

    double gamma = -b_coeff;
    double alpha = c_coeff;
    double beta  = a_coeff;

    // Modified diagonal for T matrix
    vector<double> T_a(N, a_coeff);
    vector<double> T_b(N, b_coeff);
    vector<double> T_c(N, c_coeff);
    T_b[0]     = b_coeff - gamma;
    T_b[N - 1] = b_coeff - (alpha * beta) / gamma;

    // u vector: [gamma, 0, ..., 0, alpha]
    vector<double> u_vec(N, 0.0);
    u_vec[0]     = gamma;
    u_vec[N - 1] = alpha;

    vector<double> y = solveTDMA(T_a, T_b, T_c, d);
    vector<double> z = solveTDMA(T_a, T_b, T_c, u_vec);

    double rho = (y[0] + (beta / gamma) * y[N - 1]) /
                 (1.0 + z[0] + (beta / gamma) * z[N - 1]);

    vector<double> sol(N);
    for (int i = 0; i < N; i++)
        sol[i] = y[i] - rho * z[i];
    return sol;
}

// ─── Scheme 4A: FTCS-UW (fully explicit) ─────────────────────────────────────
void solveFTCS_UW(vector<double> &u, const vector<double> &x,
                  double dt, double dx, double nu, double r,
                  int Nx, int Nt, const vector<double> &output_times)
{
    writeSolution(u, x, 0.0, Nx, "FTCS_UW");
    vector<double> u_new(Nx, 0.0);
    int next_out = 1;

    for (int n = 1; n <= Nt; n++)
    {
        for (int i = 0; i < Nx; i++)
        {
            int im1 = (i - 1 + Nx) % Nx;
            int ip1 = (i + 1) % Nx;
            u_new[i] = u[i]
                       - nu * (u[i] - u[im1])
                       + r  * (u[ip1] - 2.0*u[i] + u[im1]);
        }
        u = u_new;

        double t_now = n * dt;
        if (next_out < (int)output_times.size() && t_now >= output_times[next_out])
        {
            writeSolution(u, x, t_now, Nx, "FTCS_UW");
            next_out++;
        }
    }
    printf("FTCS-UW done.\n");
}

// ─── Scheme 4B: BTCS-UW (fully implicit) ─────────────────────────────────────
// (1+nu+2r)*u^{n+1}_i + (-nu-r)*u^{n+1}_{i-1} + (-r)*u^{n+1}_{i+1} = u^n_i
void solveBTCS_UW(vector<double> &u, const vector<double> &x,
                  double dt, double dx, double nu, double r,
                  int Nx, int Nt, const vector<double> &output_times)
{
    writeSolution(u, x, 0.0, Nx, "BTCS_UW");

    // Uniform coefficients (same every row)
    double a_coeff = -(nu + r); // coefficient of u_{i-1}
    double b_coeff =  1.0 + nu + 2.0*r; // coefficient of u_i
    double c_coeff = -r;        // coefficient of u_{i+1}

    int next_out = 1;

    for (int n = 1; n <= Nt; n++)
    {
        // RHS = u^n
        vector<double> d(u.begin(), u.end());
        u = solvePeriodicTDMA(a_coeff, b_coeff, c_coeff, d);

        double t_now = n * dt;
        if (next_out < (int)output_times.size() && t_now >= output_times[next_out])
        {
            writeSolution(u, x, t_now, Nx, "BTCS_UW");
            next_out++;
        }
    }
    printf("BTCS-UW done.\n");
}

// ─── Scheme 4C: Crank-Nicolson ────────────────────────────────────────────────
// LHS: (1 + nu/2 + r)*u^{n+1}_i + (-nu/2 - r/2)*u^{n+1}_{i-1} + (-r/2)*u^{n+1}_{i+1}
// RHS: (1 - nu/2 - r)*u^n_i    + ( nu/2 + r/2)*u^n_{i-1}      + ( r/2)*u^n_{i+1}
void solveCN(vector<double> &u, const vector<double> &x,
             double dt, double dx, double nu, double r,
             int Nx, int Nt, const vector<double> &output_times,
             const string &scheme_name)
{
    writeSolution(u, x, 0.0, Nx, scheme_name);

    // LHS uniform coefficients
    double a_coeff = -(nu/2.0 + r/2.0); // u^{n+1}_{i-1}
    double b_coeff =  1.0 + nu/2.0 + r; // u^{n+1}_i
    double c_coeff = -r/2.0;             // u^{n+1}_{i+1}

    int next_out = 1;

    for (int n = 1; n <= Nt; n++)
    {
        // Build RHS
        vector<double> d(Nx, 0.0);
        for (int i = 0; i < Nx; i++)
        {
            int im1 = (i - 1 + Nx) % Nx;
            int ip1 = (i + 1) % Nx;
            d[i] = (1.0 - nu/2.0 - r) * u[i]
                 + (nu/2.0 + r/2.0)   * u[im1]
                 + (r/2.0)             * u[ip1];
        }
        u = solvePeriodicTDMA(a_coeff, b_coeff, c_coeff, d);

        double t_now = n * dt;
        if (next_out < (int)output_times.size() && t_now >= output_times[next_out])
        {
            writeSolution(u, x, t_now, Nx, scheme_name);
            next_out++;
        }
    }
    printf("CN done (%s).\n", scheme_name.c_str());
}

int main()
{
    const double kappa  = 0.010; // Roll 23
    const double nu_base = 0.80; // Roll 23
    const double c      = 1.0;
    const int    Nx     = 100;
    const double xstart = 0.0, xend = 1.0;
    const double T      = 1.0;

    double dx = (xend - xstart) / Nx;
    double dt = nu_base * dx / c;
    double r  = kappa * dt / (dx * dx);
    int    Nt = (int)ceil(T / dt);
    double final_t = Nt * dt;

    printf("Roll 23 — 1D Advection-Diffusion\n");
    printf("kappa=%.4f, nu=%.4f, c=%.1f\n", kappa, nu_base, c);
    printf("Nx=%d, dx=%.6f, dt=%.6f, r=%.6f\n", Nx, dx, dt, r);
    printf("Pe_h = nu/r = %.4f\n", nu_base / r);
    printf("Nt=%d, final_t=%.6f\n\n", Nt, final_t);

    vector<double> x = makeGrid(xstart, xend, Nx);
    vector<double> output_times = {0.0, 0.25, 0.5, 0.75, 1.0};
    vector<double> u(Nx, 0.0);

    // ── Scheme 4A: FTCS-UW ──
    setIC(u, x, Nx);
    printf("Solving Scheme 4A: FTCS-UW...\n");
    solveFTCS_UW(u, x, dt, dx, nu_base, r, Nx, Nt, output_times);
    double e_ftcs = computeL2Error(u, x, final_t, kappa, dx, Nx);
    printf("FTCS-UW L2 error at T=%.6f: %.6e\n\n", final_t, e_ftcs);

    // ── Scheme 4B: BTCS-UW ──
    setIC(u, x, Nx);
    printf("Solving Scheme 4B: BTCS-UW...\n");
    solveBTCS_UW(u, x, dt, dx, nu_base, r, Nx, Nt, output_times);
    double e_btcs = computeL2Error(u, x, final_t, kappa, dx, Nx);
    printf("BTCS-UW L2 error at T=%.6f: %.6e\n\n", final_t, e_btcs);

    // ── Scheme 4C: Crank-Nicolson (base nu) ──
    setIC(u, x, Nx);
    printf("Solving Scheme 4C: Crank-Nicolson (nu=%.2f)...\n", nu_base);
    solveCN(u, x, dt, dx, nu_base, r, Nx, Nt, output_times, "CN");
    double e_cn = computeL2Error(u, x, final_t, kappa, dx, Nx);
    printf("CN L2 error at T=%.6f: %.6e\n\n", final_t, e_cn);

    // ── Summary table ──
    printf("=== L2 Error Summary ===\n");
    printf("Scheme     | L2 error\n");
    printf("FTCS-UW    | %.6e\n", e_ftcs);
    printf("BTCS-UW    | %.6e\n", e_btcs);
    printf("CN         | %.6e\n\n", e_cn);

    // ── Part (d): CN with nu = 2, 5, 10 ──
    vector<double> nu_vals = {2.0, 5.0, 10.0};
    for (double nu_val : nu_vals)
    {
        double dt_val = nu_val * dx / c;
        double r_val  = kappa * dt_val / (dx * dx);
        int    Nt_val = (int)ceil(T / dt_val);
        double ft_val = Nt_val * dt_val;

        setIC(u, x, Nx);
        string sname = "CN_nu" + to_string((int)nu_val);
        printf("Solving CN (nu=%.1f, r=%.4f)...\n", nu_val, r_val);
        solveCN(u, x, dt_val, dx, nu_val, r_val, Nx, Nt_val, output_times, sname);
        double e_val = computeL2Error(u, x, ft_val, kappa, dx, Nx);
        printf("CN (nu=%.1f) L2 error at T=%.6f: %.6e\n\n", nu_val, ft_val, e_val);
    }

    printf("Program completed.\n");
    return 0;
}   