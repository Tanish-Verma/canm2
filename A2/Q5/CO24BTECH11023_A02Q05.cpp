#include <bits/stdc++.h>
using namespace std;

// Roll 23: Pe=100, kappa=1/100=0.01, c=1
// Part A: steady-state, N=50, Dirichlet BCs u(0)=0, u(1)=1
// Part B: transient CN, T=0.5, output times 0,0.125,0.25,0.375,0.5 (ALL roll groups)

// ─── Grid (includes both boundary points) ────────────────────────────────────
// N intervals => N+1 points: x[0]=0, x[1]=dx, ..., x[N]=1
vector<double> makeGrid(int N)
{
    vector<double> g(N + 1);
    double dx = 1.0 / N;
    for (int i = 0; i <= N; i++)
        g[i] = i * dx;
    return g;
}

// ─── Plain (non-periodic) TDMA ────────────────────────────────────────────────
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

// ─── Periodic TDMA (Sherman-Morrison) ────────────────────────────────────────
vector<double> solvePeriodicTDMA(double a_coeff, double b_coeff, double c_coeff,
                                  const vector<double> &d)
{
    int N = d.size();
    double gamma = -b_coeff;
    double alpha = c_coeff;
    double beta  = a_coeff;

    vector<double> T_a(N, a_coeff), T_b(N, b_coeff), T_c(N, c_coeff);
    T_b[0]     = b_coeff - gamma;
    T_b[N - 1] = b_coeff - (alpha * beta) / gamma;

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

// ─── Steady-state exact solution ─────────────────────────────────────────────
double steadyExact(double x, double Pe)
{
    return (exp(Pe * x) - 1.0) / (exp(Pe) - 1.0);
}

// ─── Part A: Solve steady-state with central differences for ux ──────────────
// c*ux = kappa*uxx  =>  on interior i=1..N-1:
// (-kappa/dx^2 - c/(2dx))*u_{i-1} + (2*kappa/dx^2)*u_i + (-kappa/dx^2 + c/(2dx))*u_{i+1} = 0
// BCs: u[0]=0, u[N]=1
void solveSteadyCentral(int N, double Pe, double kappa, double c)
{
    vector<double> x = makeGrid(N);
    double dx = 1.0 / N;
    int Ni = N - 1; // number of interior points

    double A =  -kappa / (dx * dx) - c / (2.0 * dx); // sub-diagonal
    double B =   2.0 * kappa / (dx * dx);              // main diagonal
    double C =  -kappa / (dx * dx) + c / (2.0 * dx);  // super-diagonal

    vector<double> a(Ni, A), b(Ni, B), cc(Ni, C), d(Ni, 0.0);

    // BC contribution: u[0]=0 adds nothing to first eq
    // BC contribution: u[N]=1 subtracts C from last eq RHS
    d[Ni - 1] -= C * 1.0;

    vector<double> u_int = solveTDMA(a, b, cc, d);

    // Assemble full solution including boundaries
    vector<double> u(N + 1, 0.0);
    u[0] = 0.0;
    u[N] = 1.0;
    for (int i = 1; i < N; i++)
        u[i] = u_int[i - 1];

    // Write to file
    ofstream out("data/steady_central.txt");
    out << fixed << setprecision(8);
    for (int i = 0; i <= N; i++)
        out << x[i] << " " << u[i] << " " << steadyExact(x[i], Pe) << "\n";
    printf("Steady-state (central diff) written.\n");
}

// ─── Part A: Solve steady-state with upwind differences for ux ───────────────
// c*(u_i - u_{i-1})/dx = kappa*(u_{i+1} - 2*u_i + u_{i-1})/dx^2
// => (-kappa/dx^2 - c/dx)*u_{i-1} + (2*kappa/dx^2 + c/dx)*u_i + (-kappa/dx^2)*u_{i+1} = 0
void solveSteadyUpwind(int N, double Pe, double kappa, double c)
{
    vector<double> x = makeGrid(N);
    double dx = 1.0 / N;
    int Ni = N - 1;

    double A = -kappa / (dx * dx) - c / dx; // sub-diagonal
    double B =  2.0 * kappa / (dx * dx) + c / dx; // main diagonal
    double C = -kappa / (dx * dx);           // super-diagonal

    vector<double> a(Ni, A), b(Ni, B), cc(Ni, C), d(Ni, 0.0);
    d[Ni - 1] -= C * 1.0; // BC: u[N]=1

    vector<double> u_int = solveTDMA(a, b, cc, d);

    vector<double> u(N + 1, 0.0);
    u[0] = 0.0;
    u[N] = 1.0;
    for (int i = 1; i < N; i++)
        u[i] = u_int[i - 1];

    // Write solution
    ofstream out("data/steady_upwind_N" + to_string(N) + ".txt");
    out << fixed << setprecision(8);
    for (int i = 0; i <= N; i++)
        out << x[i] << " " << u[i] << " " << steadyExact(x[i], Pe) << "\n";
}

// Compute L_inf error for upwind steady-state at a given N
double steadyUpwindLinfError(int N, double Pe, double kappa, double c)
{
    vector<double> x = makeGrid(N);
    double dx = 1.0 / N;
    int Ni = N - 1;

    double A = -kappa / (dx * dx) - c / dx;
    double B =  2.0 * kappa / (dx * dx) + c / dx;
    double C = -kappa / (dx * dx);

    vector<double> a(Ni, A), b(Ni, B), cc(Ni, C), d(Ni, 0.0);
    d[Ni - 1] -= C * 1.0;
    vector<double> u_int = solveTDMA(a, b, cc, d);

    vector<double> u(N + 1, 0.0);
    u[0] = 0.0; u[N] = 1.0;
    for (int i = 1; i < N; i++) u[i] = u_int[i - 1];

    double err = 0.0;
    for (int i = 0; i <= N; i++)
        err = max(err, fabs(u[i] - steadyExact(x[i], Pe)));
    return err;
}

// ─── Transient exact solution ─────────────────────────────────────────────────
double transientExact(double x, double t, double kappa)
{
    return exp(-4.0 * M_PI * M_PI * kappa * t) * sin(2.0 * M_PI * (x - t)); // c is 1
}

// ─── L2 error (1D periodic) ───────────────────────────────────────────────────
double computeL2Error(const vector<double> &u, const vector<double> &x,
                      double t, double kappa, double dx, int Nx)
{
    double sum = 0.0;
    for (int i = 0; i < Nx; i++)
    {
        double diff = u[i] - transientExact(x[i], t, kappa);
        sum += diff * diff;
    }
    return sqrt(dx * sum);
}

// ─── CN transient solver (periodic BCs) ──────────────────────────────────────
// Same as Q4: upwind advection + central diffusion, time-centred
// LHS: (1+nu/2+r)*u^{n+1}_i + (-nu/2-r/2)*u^{n+1}_{i-1} + (-r/2)*u^{n+1}_{i+1}
// RHS: (1-nu/2-r)*u^n_i    + ( nu/2+r/2)*u^n_{i-1}      + ( r/2)*u^n_{i+1}
void solveCN(vector<double> &u, const vector<double> &x,
             double dt, double dx, double nu, double r,
             int Nx, int Nt, const vector<double> &output_times,
             const string &scheme_name)
{
    // Write IC
    {
        ofstream out("data/" + scheme_name + "_T_0.000000.txt");
        out << fixed << setprecision(8);
        for (int i = 0; i < Nx; i++) out << x[i] << " " << u[i] << "\n";
    }

    double a_coeff = -(nu / 2.0 + r / 2.0);
    double b_coeff =  1.0 + nu / 2.0 + r;
    double c_coeff = -r / 2.0;

    int next_out = 1;

    for (int n = 1; n <= Nt; n++)
    {
        vector<double> d(Nx, 0.0);
        for (int i = 0; i < Nx; i++)
        {
            int im1 = (i - 1 + Nx) % Nx;
            int ip1 = (i + 1) % Nx;
            d[i] = (1.0 - nu / 2.0 - r) * u[i]
                 + (nu / 2.0 + r / 2.0)  * u[im1]
                 + (r / 2.0)              * u[ip1];
        }
        u = solvePeriodicTDMA(a_coeff, b_coeff, c_coeff, d);

        double t_now = n * dt;
        if (next_out < (int)output_times.size() && t_now >= output_times[next_out] - 1e-12)
        {
            ostringstream fname;
            fname << "data/" << scheme_name << "_T_" << fixed << setprecision(6) << t_now << ".txt";
            ofstream out(fname.str());
            out << fixed << setprecision(8);
            for (int i = 0; i < Nx; i++) out << x[i] << " " << u[i] << "\n";
            next_out++;
        }
    }
}

int main()
{
    const double Pe     = 100.0; // Roll 23
    const double kappa  = 1.0 / Pe;
    const double c      = 1.0;
    const double T_tr   = 0.5;   // transient T (all roll groups for Q5)

    printf("Roll 23 — Q5 Advection-Diffusion: Peclet Regimes\n");
    printf("Pe=%.0f, kappa=%.6f, c=%.1f\n\n", Pe, kappa, c);

    // ── Part A(a): Steady-state, N=50 ──
    printf("=== Part A: Steady-State ===\n");
    int N_ss = 50;
    double dx_ss = 1.0 / N_ss;
    double Peh_ss = c * dx_ss / kappa;
    printf("N=%d, dx=%.6f, Pe_h=%.4f\n", N_ss, dx_ss, Peh_ss);
    printf("Pe_h=%.2f %s 2 => central diff %s show oscillations\n",
           Peh_ss, Peh_ss > 2 ? ">" : "<=",
           Peh_ss > 2 ? "WILL" : "may not");

    solveSteadyCentral(N_ss, Pe, kappa, c);
    solveSteadyUpwind(N_ss, Pe, kappa, c);

    // ── Part A(c): Upwind L_inf convergence ──
    printf("\nUpwind L_inf convergence:\n");
    vector<int> N_vals = {10, 20, 50, 100, 200};
    ofstream conv_out("data/steady_convergence.txt");
    conv_out << fixed << setprecision(10);
    for (int N : N_vals)
    {
        double dx  = 1.0 / N;
        double err = steadyUpwindLinfError(N, Pe, kappa, c);
        printf("  N=%3d, dx=%.6f, L_inf=%.6e\n", N, dx, err);
        conv_out << dx << " " << err << "\n";
        solveSteadyUpwind(N, Pe, kappa, c); // also write solution for each N
    }
    printf("\n");

    // ── Part B(a): CN transient, N=100, nu=0.5, T=0.5 ──
    printf("=== Part B: Transient CN ===\n");
    int    Nx_b   = 100;
    double nu_b   = 0.5;
    double dx_b   = 1.0 / Nx_b;
    double dt_b   = nu_b * dx_b / c;
    double r_b    = kappa * dt_b / (dx_b * dx_b);
    int    Nt_b   = (int)ceil(T_tr / dt_b);
    double final_t_b = Nt_b * dt_b;

    printf("Part B(a): N=%d, nu=%.2f, dx=%.6f, dt=%.6f, r=%.6f\n",
           Nx_b, nu_b, dx_b, dt_b, r_b);
    printf("Nt=%d, final_t=%.6f\n", Nt_b, final_t_b);

    vector<double> x_b = makeGrid(Nx_b); // Nx_b intervals => Nx_b+1 points, but periodic uses Nx_b points
    // For periodic: use Nx_b points at x=0,dx,...,(Nx_b-1)*dx (exclude right endpoint)
    vector<double> xp(Nx_b);
    for (int i = 0; i < Nx_b; i++) xp[i] = i * dx_b;

    vector<double> output_times_b = {0.0, 0.125, 0.25, 0.375, 0.5};
    vector<double> u(Nx_b);
    for (int i = 0; i < Nx_b; i++) u[i] = sin(2.0 * M_PI * xp[i]);

    solveCN(u, xp, dt_b, dx_b, nu_b, r_b, Nx_b, Nt_b, output_times_b, "CN_transient");
    double e_ba = computeL2Error(u, xp, final_t_b, kappa, dx_b, Nx_b);
    printf("Part B(a) L2 error at T=%.4f: %.6e\n\n", final_t_b, e_ba);

    // ── Part B(b): Spatial convergence, fix nu=0.5, vary N ──
    printf("Part B(b): Spatial convergence (nu=0.5):\n");
    vector<int> N_space = {25, 50, 100, 200, 400};
    ofstream sp_out("data/spatial_convergence.txt");
    sp_out << fixed << setprecision(10);
    for (int N : N_space)
    {
        double dx  = 1.0 / N;
        double dt  = nu_b * dx / c;
        double r   = kappa * dt / (dx * dx);
        int    Nt  = (int)ceil(T_tr / dt);
        double ft  = Nt * dt;

        vector<double> xpi(N);
        for (int i = 0; i < N; i++) xpi[i] = i * dx;
        vector<double> ui(N);
        for (int i = 0; i < N; i++) ui[i] = sin(2.0 * M_PI * xpi[i]);

        vector<double> dummy_out = {0.0, T_tr};
        solveCN(ui, xpi, dt, dx, nu_b, r, N, Nt, dummy_out, "dummy_space");

        double err = computeL2Error(ui, xpi, ft, kappa, dx, N);
        printf("  N=%3d, dx=%.6f, L2=%.6e\n", N, dx, err);
        sp_out << dx << " " << err << "\n";
    }
    printf("\n");

    // ── Part B(c): Temporal convergence, fix N=200, vary nu ──
    printf("Part B(c): Temporal convergence (N=200):\n");
    int N_tc = 200;
    double dx_tc = 1.0 / N_tc;
    vector<double> nu_vals = {0.8, 0.4, 0.2, 0.1, 0.05};
    ofstream tp_out("data/temporal_convergence.txt");
    tp_out << fixed << setprecision(10);

    vector<double> xp_tc(N_tc);
    for (int i = 0; i < N_tc; i++) xp_tc[i] = i * dx_tc;

    for (double nu_val : nu_vals)
    {
        double dt  = nu_val * dx_tc / c;
        double r   = kappa * dt / (dx_tc * dx_tc);
        int    Nt  = (int)ceil(T_tr / dt);
        double ft  = Nt * dt;

        vector<double> ui(N_tc);
        for (int i = 0; i < N_tc; i++) ui[i] = sin(2.0 * M_PI * xp_tc[i]);

        vector<double> dummy_out = {0.0, T_tr};
        solveCN(ui, xp_tc, dt, dx_tc, nu_val, r, N_tc, Nt, dummy_out, "dummy_time");

        double err = computeL2Error(ui, xp_tc, ft, kappa, dx_tc, N_tc);
        printf("  nu=%.2f, dt=%.6f, L2=%.6e\n", nu_val, dt, err);
        tp_out << dt << " " << err << "\n";
    }
    printf("\n");

    // ── Part B(d): Compare Pe=100 vs Pe=5 ──
    printf("Part B(d): Pe comparison (N=100, nu=0.5, T=0.5):\n");
    vector<double> Pe_compare = {5.0, 100.0};
    for (double Pe_val : Pe_compare)
    {
        double kap = 1.0 / Pe_val;
        double dx  = 1.0 / Nx_b;
        double dt  = nu_b * dx / c;
        double r   = kap * dt / (dx * dx);
        int    Nt  = (int)ceil(T_tr / dt);
        double ft  = Nt * dt;

        vector<double> xpi(Nx_b);
        for (int i = 0; i < Nx_b; i++) xpi[i] = i * dx;
        vector<double> ui(Nx_b);
        for (int i = 0; i < Nx_b; i++) ui[i] = sin(2.0 * M_PI * xpi[i]);

        vector<double> dummy_out = {0.0, T_tr};
        solveCN(ui, xpi, dt, dx, nu_b, r, Nx_b, Nt, dummy_out, "dummy_pe");

        double err = computeL2Error(ui, xpi, ft, kap, dx, Nx_b);
        printf("  Pe=%.0f, kappa=%.6f, r=%.4f, L2=%.6e\n", Pe_val, kap, r, err);
    }

    printf("\nProgram completed.\n");
    return 0;
}