#include <bits/stdc++.h>
using namespace std;

// Roll 23: kappa0 = 0.20, beta = 0.4
// Nx = Ny = 50, dx = dy = 1/50
// r0 = kappa0 * dt / dx^2 = 0.4  =>  dt = 0.4 * dx^2 / kappa0
// T = 0.2 s
// IC: u(x,y,0) = sin(pi*x)*sin(pi*y)
// BC: u = 0 on all boundaries (Dirichlet)

// ─── Grid ─────────────────────────────────────────────────────────────────────
// Interior points: i = 1..N-2 (0-indexed), boundary at i=0 and i=N-1 are fixed = 0
// Grid includes boundary: x[i] = i*dx for i=0..N-1, with N = Nx = Ny

vector<double> makeGrid(int N)
{
    vector<double> g(N);
    double dx = 1.0 / (N - 1); // N points including both endpoints
    for (int i = 0; i < N; i++)
        g[i] = i * dx;
    return g;
}

// κ(x,y) = κ0 * (1 + β*sin(π*x)*sin(π*y))
double kappa(double x, double y, double kappa0, double beta)
{
    return kappa0 * (1.0 + beta * sin(M_PI * x) * sin(M_PI * y));
}

// Set IC: u(x,y,0) = sin(π*x)*sin(π*y)  (zero on boundary automatically)
void setIC(vector<vector<double>> &u, const vector<double> &x, const vector<double> &y, int N)
{
    for (int i = 0; i < N; i++)
        for (int j = 0; j < N; j++)
            u[i][j] = sin(M_PI * x[i]) * sin(M_PI * y[j]);
}

// Exact solution (β=0 only): u_ex = exp(-2*π²*κ0*t) * sin(π*x)*sin(π*y)
double exactSol(double x, double y, double t, double kappa0)
{
    return exp(-2.0 * M_PI * M_PI * kappa0 * t) * sin(M_PI * x) * sin(M_PI * y);
}

// 2D L2 error: sqrt(dx*dy * sum (u_ij - u_ex_ij)^2)
double computeL2Error(const vector<vector<double>> &u,
                      const vector<double> &x, const vector<double> &y,
                      double t, double kappa0, double dx, double dy, int N)
{
    double sum = 0.0;
    for (int i = 0; i < N; i++)
        for (int j = 0; j < N; j++)
        {
            double diff = u[i][j] - exactSol(x[i], y[j], t, kappa0);
            sum += diff * diff;
        }
    return sqrt(dx * dy * sum);
}

// Max of u over all interior points
double maxU(const vector<vector<double>> &u, int N)
{
    double m = 0.0;
    for (int i = 1; i < N - 1; i++)
        for (int j = 1; j < N - 1; j++)
            m = max(m, u[i][j]);
    return m;
}

// Write 2D snapshot: "x y u" for all points
void writeSolution(const vector<vector<double>> &u,
                   const vector<double> &x, const vector<double> &y,
                   int N, double t_actual, const string &prefix)
{
    ostringstream fname;
    fname << "data/" << prefix << "_T_" << fixed << setprecision(6) << t_actual << ".txt";
    ofstream out(fname.str());
    out << fixed << setprecision(8);
    for (int i = 0; i < N; i++)
        for (int j = 0; j < N; j++)
            out << x[i] << " " << y[j] << " " << u[i][j] << "\n";
}

// Write kappa field
void writeKappa(const vector<double> &x, const vector<double> &y,
                int N, double kappa0, double beta)
{
    ofstream out("data/kappa.txt");
    out << fixed << setprecision(8);
    for (int i = 0; i < N; i++)
        for (int j = 0; j < N; j++)
            out << x[i] << " " << y[j] << " " << kappa(x[i], y[j], kappa0, beta) << "\n";
}

// ─── Thomas (TDMA) algorithm ──────────────────────────────────────────────────
// Solves: a[i]*u[i-1] + b[i]*u[i] + c[i]*u[i+1] = d[i]
// Standard (non-periodic) tridiagonal — valid for Dirichlet BCs
// Only interior points passed in (size = N_int = N-2)
vector<double> solveTDMA(const vector<double> &a, const vector<double> &b,
                          const vector<double> &c, const vector<double> &d)
{
    int n = d.size();
    vector<double> c_star(n, 0.0), d_star(n, 0.0), sol(n, 0.0);

    // Forward sweep
    c_star[0] = c[0] / b[0];
    d_star[0] = d[0] / b[0];
    for (int i = 1; i < n; i++)
    {
        double m = b[i] - a[i] * c_star[i - 1];
        c_star[i] = (i < n - 1) ? c[i] / m : 0.0;
        d_star[i] = (d[i] - a[i] * d_star[i - 1]) / m;
    }

    // Backward substitution
    sol[n - 1] = d_star[n - 1];
    for (int i = n - 2; i >= 0; i--)
        sol[i] = d_star[i] - c_star[i] * sol[i + 1];

    return sol;
}

// ─── ADI solver ───────────────────────────────────────────────────────────────
// u: N×N grid (includes boundary points fixed at 0)
// Returns max_u_history: max interior u at every step
vector<double> solveADI(vector<vector<double>> &u,
                         const vector<double> &x, const vector<double> &y,
                         double dt, double dx, double dy,
                         double kappa0, double beta, int N, int Nt,
                         const vector<double> &output_times,
                         const string &prefix)
{
    int Ni = N - 2; // number of interior points per row/col

    // Precompute κ at all grid points
    vector<vector<double>> kap(N, vector<double>(N));
    for (int i = 0; i < N; i++)
        for (int j = 0; j < N; j++)
            kap[i][j] = kappa(x[i], y[j], kappa0, beta);

    // Precompute interface values (arithmetic average)
    // kap_xp[i][j] = κ_{i+1/2,j},  kap_xm[i][j] = κ_{i-1/2,j}
    // kap_yp[i][j] = κ_{i,j+1/2},  kap_ym[i][j] = κ_{i,j-1/2}
    vector<vector<double>> kap_xp(N, vector<double>(N));
    vector<vector<double>> kap_xm(N, vector<double>(N));
    vector<vector<double>> kap_yp(N, vector<double>(N));
    vector<vector<double>> kap_ym(N, vector<double>(N));
    for (int i = 0; i < N - 1; i++)
        for (int j = 0; j < N; j++)
        {
            kap_xp[i][j] = 0.5 * (kap[i][j] + kap[i + 1][j]);
            kap_xm[i + 1][j] = kap_xp[i][j]; // κ_{i+1/2} = κ_{(i+1)-1/2}
        }
    for (int i = 0; i < N; i++)
        for (int j = 0; j < N - 1; j++)
        {
            kap_yp[i][j] = 0.5 * (kap[i][j] + kap[i][j + 1]);
            kap_ym[i][j + 1] = kap_yp[i][j];
        }

    // Write IC
    writeSolution(u, x, y, N, 0.0, prefix);

    vector<double> max_u_history;
    max_u_history.push_back(maxU(u, N));

    vector<vector<double>> u_half(N, vector<double>(N, 0.0));
    int next_out = 1;

    for (int n = 1; n <= Nt; n++)
    {
        // ── X-Sweep: implicit in x, explicit in y ──────────────────────────
        // For each row j (fixed), solve tridiagonal in i=1..N-2
        // [u^{n+1/2} - u^n] / (dt/2) = Lx[u^{n+1/2}] + Ly[u^n]
        // => -rx*kxm * u^{n+1/2}_{i-1} + (1 + rx*(kxm+kxp)) * u^{n+1/2}_i - rx*kxp * u^{n+1/2}_{i+1}
        //    = u^n_i + (dt/2) * Ly[u^n]_i
        // where rx = dt/(2*dx^2)

        double rx = dt / (2.0 * dx * dx);
        double ry = dt / (2.0 * dy * dy);

        for (int j = 1; j < N - 1; j++) // each interior row
        {
            vector<double> a(Ni), b(Ni), c(Ni), d(Ni);

            for (int ii = 0; ii < Ni; ii++)
            {
                int i = ii + 1; // actual grid index

                double kxp = kap_xp[i][j];
                double kxm = kap_xm[i][j];

                // Ly[u^n]_i,j = (kyp*(u_{i,j+1}-u_{i,j}) - kym*(u_{i,j}-u_{i,j-1})) / dy^2
                double Ly_un = (kap_yp[i][j] * (u[i][j + 1] - u[i][j])
                               - kap_ym[i][j] * (u[i][j] - u[i][j - 1])) / (dy * dy);

                a[ii] = -rx * kxm;
                b[ii] =  1.0 + rx * (kxm + kxp);
                c[ii] = -rx * kxp;
                d[ii] =  u[i][j] + (dt / 2.0) * Ly_un;
            }

            // Boundary contributions (u=0 on boundary, so no correction needed)
            // a[0]*u_boundary_left + ... → since u_boundary=0, no extra term

            vector<double> sol = solveTDMA(a, b, c, d);
            for (int ii = 0; ii < Ni; ii++)
                u_half[ii + 1][j] = sol[ii];
        }
        // Boundary of u_half stays 0
        for (int i = 0; i < N; i++) { u_half[i][0] = 0.0; u_half[i][N-1] = 0.0; }
        for (int j = 0; j < N; j++) { u_half[0][j] = 0.0; u_half[N-1][j] = 0.0; }

        // ── Y-Sweep: implicit in y, explicit in x (using u_half) ──────────
        // [u^{n+1} - u^{n+1/2}] / (dt/2) = Lx[u^{n+1/2}] + Ly[u^{n+1}]
        // => -ry*kym * u^{n+1}_{i,j-1} + (1 + ry*(kym+kyp)) * u^{n+1}_{i,j} - ry*kyp * u^{n+1}_{i,j+1}
        //    = u^{n+1/2}_{i,j} + (dt/2) * Lx[u^{n+1/2}]_{i,j}

        for (int i = 1; i < N - 1; i++) // each interior column
        {
            vector<double> a(Ni), b(Ni), c(Ni), d(Ni);

            for (int jj = 0; jj < Ni; jj++)
            {
                int j = jj + 1;

                double kyp = kap_yp[i][j];
                double kym = kap_ym[i][j];

                // Lx[u^{n+1/2}]_{i,j} = (kxp*(u_half_{i+1,j}-u_half_{i,j}) - kxm*(u_half_{i,j}-u_half_{i-1,j})) / dx^2
                double Lx_uhalf = (kap_xp[i][j] * (u_half[i + 1][j] - u_half[i][j])
                                  - kap_xm[i][j] * (u_half[i][j] - u_half[i - 1][j])) / (dx * dx);

                a[jj] = -ry * kym;
                b[jj] =  1.0 + ry * (kym + kyp);
                c[jj] = -ry * kyp;
                d[jj] =  u_half[i][j] + (dt / 2.0) * Lx_uhalf;
            }

            vector<double> sol = solveTDMA(a, b, c, d);
            for (int jj = 0; jj < Ni; jj++)
                u[i][jj + 1] = sol[jj];
        }
        // Boundary of u stays 0
        for (int i = 0; i < N; i++) { u[i][0] = 0.0; u[i][N-1] = 0.0; }
        for (int j = 0; j < N; j++) { u[0][j] = 0.0; u[N-1][j] = 0.0; }

        max_u_history.push_back(maxU(u, N));

        double t_now = n * dt;
        if (next_out < (int)output_times.size() && t_now >= output_times[next_out] - 1e-12)
        {
            writeSolution(u, x, y, N, t_now, prefix);
            next_out++;
        }
    }

    printf("ADI done (beta=%.2f).\n", beta);
    return max_u_history;
}

// ─── Convergence study (beta=0) ───────────────────────────────────────────────
void convergenceStudy(double kappa0, double r0)
{
    vector<int> Nvals = {25, 50, 100}; // grid sizes (including boundary)
    ofstream out("data/convergence.txt");
    out << fixed << setprecision(10);

    for (int N : Nvals)
    {
        double dx  = 1.0 / (N - 1);
        double dy  = dx;
        double dt  = r0 * dx * dx / kappa0;
        double T   = 0.2;
        int    Nt  = (int)ceil(T / dt);
        double final_t = Nt * dt;

        vector<double> x = makeGrid(N);
        vector<double> y = makeGrid(N);
        vector<vector<double>> u(N, vector<double>(N, 0.0));
        setIC(u, x, y, N);

        // Run ADI with beta=0 (constant kappa), no file output
        vector<double> dummy_output = {0.0, T};  // only care about final state
        // We reuse solveADI with a dummy prefix that won't be read
        // but to avoid file clutter we inline a minimal loop here

        int Ni = N - 2;
        double rx = dt / (2.0 * dx * dx);
        double ry = dt / (2.0 * dy * dy);
        // With beta=0, kappa is constant = kappa0 everywhere
        // Interface values are all kappa0
        double kxp = kappa0, kxm = kappa0, kyp = kappa0, kym = kappa0;

        vector<vector<double>> u_half(N, vector<double>(N, 0.0));

        for (int n = 1; n <= Nt; n++)
        {
            // X-sweep
            for (int j = 1; j < N - 1; j++)
            {
                vector<double> a(Ni), b(Ni), c(Ni), d(Ni);
                for (int ii = 0; ii < Ni; ii++)
                {
                    int i = ii + 1;
                    double Ly_un = kappa0 * (u[i][j+1] - 2.0*u[i][j] + u[i][j-1]) / (dy*dy);
                    a[ii] = -rx * kxm;
                    b[ii] =  1.0 + rx * (kxm + kxp);
                    c[ii] = -rx * kxp;
                    d[ii] =  u[i][j] + (dt/2.0) * Ly_un;
                }
                vector<double> sol = solveTDMA(a, b, c, d);
                for (int ii = 0; ii < Ni; ii++)
                    u_half[ii+1][j] = sol[ii];
            }
            for (int i = 0; i < N; i++) { u_half[i][0]=0; u_half[i][N-1]=0; }
            for (int j = 0; j < N; j++) { u_half[0][j]=0; u_half[N-1][j]=0; }

            // Y-sweep
            for (int i = 1; i < N - 1; i++)
            {
                vector<double> a(Ni), b(Ni), c(Ni), d(Ni);
                for (int jj = 0; jj < Ni; jj++)
                {
                    int j = jj + 1;
                    double Lx_uhalf = kappa0 * (u_half[i+1][j] - 2.0*u_half[i][j] + u_half[i-1][j]) / (dx*dx);
                    a[jj] = -ry * kym;
                    b[jj] =  1.0 + ry * (kym + kyp);
                    c[jj] = -ry * kyp;
                    d[jj] =  u_half[i][j] + (dt/2.0) * Lx_uhalf;
                }
                vector<double> sol = solveTDMA(a, b, c, d);
                for (int jj = 0; jj < Ni; jj++)
                    u[i][jj+1] = sol[jj];
            }
            for (int i = 0; i < N; i++) { u[i][0]=0; u[i][N-1]=0; }
            for (int j = 0; j < N; j++) { u[0][j]=0; u[N-1][j]=0; }
        }

        double err = computeL2Error(u, x, y, final_t, kappa0, dx, dy, N);
        printf("  N=%3d, dx=%.6f, final_t=%.6f, L2 error=%.6e\n", N, dx, final_t, err);
        out << dx << " " << err << "\n";
    }
}

// ─── Max-u decay study (part d) ───────────────────────────────────────────────
void decayStudy(double kappa0, double r0, double T,
                const vector<double> &betas, const vector<double> &x, const vector<double> &y, int N)
{
    double dx = 1.0 / (N - 1);
    double dy = dx;
    double dt = r0 * dx * dx / kappa0;
    int    Nt = (int)ceil(T / dt);

    // Write time axis
    ofstream t_out("data/decay_time.txt");
    t_out << fixed << setprecision(8);
    for (int n = 0; n <= Nt; n++)
        t_out << n * dt << "\n";

    for (double beta : betas)
    {
        vector<vector<double>> u(N, vector<double>(N, 0.0));
        setIC(u, x, y, N);

        // Use full solveADI but capture max_u_history
        vector<double> dummy_out = {0.0, T};
        ostringstream prefix;
        prefix << "decay_beta" << (int)(beta * 10);

        // Inline the ADI to capture max history without extra file writes
        int Ni = N - 2;
        double rx = dt / (2.0 * dx * dx);
        double ry = dt / (2.0 * dy * dy);

        // Precompute kappa fields
        vector<vector<double>> kap(N, vector<double>(N));
        for (int i = 0; i < N; i++)
            for (int j = 0; j < N; j++)
                kap[i][j] = kappa(x[i], y[j], kappa0, beta);

        vector<vector<double>> kap_xp(N, vector<double>(N));
        vector<vector<double>> kap_xm(N, vector<double>(N));
        vector<vector<double>> kap_yp(N, vector<double>(N));
        vector<vector<double>> kap_ym(N, vector<double>(N));
        for (int i = 0; i < N-1; i++)
            for (int j = 0; j < N; j++)
            { kap_xp[i][j]=0.5*(kap[i][j]+kap[i+1][j]); kap_xm[i+1][j]=kap_xp[i][j]; }
        for (int i = 0; i < N; i++)
            for (int j = 0; j < N-1; j++)
            { kap_yp[i][j]=0.5*(kap[i][j]+kap[i][j+1]); kap_ym[i][j+1]=kap_yp[i][j]; }

        ofstream max_out("data/decay_beta" + to_string((int)(beta*10)) + ".txt");
        max_out << fixed << setprecision(8);
        max_out << maxU(u, N) << "\n";

        vector<vector<double>> u_half(N, vector<double>(N, 0.0));

        for (int n = 1; n <= Nt; n++)
        {
            // X-sweep
            for (int j = 1; j < N-1; j++)
            {
                vector<double> a(Ni), b(Ni), c(Ni), d(Ni);
                for (int ii = 0; ii < Ni; ii++)
                {
                    int i = ii+1;
                    double Ly_un = (kap_yp[i][j]*(u[i][j+1]-u[i][j])
                                  - kap_ym[i][j]*(u[i][j]-u[i][j-1])) / (dy*dy);
                    a[ii] = -rx * kap_xm[i][j];
                    b[ii] =  1.0 + rx*(kap_xm[i][j]+kap_xp[i][j]);
                    c[ii] = -rx * kap_xp[i][j];
                    d[ii] =  u[i][j] + (dt/2.0)*Ly_un;
                }
                vector<double> sol = solveTDMA(a, b, c, d);
                for (int ii = 0; ii < Ni; ii++) u_half[ii+1][j] = sol[ii];
            }
            for (int i=0;i<N;i++){u_half[i][0]=0;u_half[i][N-1]=0;}
            for (int j=0;j<N;j++){u_half[0][j]=0;u_half[N-1][j]=0;}

            // Y-sweep
            for (int i = 1; i < N-1; i++)
            {
                vector<double> a(Ni), b(Ni), c(Ni), d(Ni);
                for (int jj = 0; jj < Ni; jj++)
                {
                    int j = jj+1;
                    double Lx_uhalf = (kap_xp[i][j]*(u_half[i+1][j]-u_half[i][j])
                                     - kap_xm[i][j]*(u_half[i][j]-u_half[i-1][j])) / (dx*dx);
                    a[jj] = -ry * kap_ym[i][j];
                    b[jj] =  1.0 + ry*(kap_ym[i][j]+kap_yp[i][j]);
                    c[jj] = -ry * kap_yp[i][j];
                    d[jj] =  u_half[i][j] + (dt/2.0)*Lx_uhalf;
                }
                vector<double> sol = solveTDMA(a, b, c, d);
                for (int jj = 0; jj < Ni; jj++) u[i][jj+1] = sol[jj];
            }
            for (int i=0;i<N;i++){u[i][0]=0;u[i][N-1]=0;}
            for (int j=0;j<N;j++){u[0][j]=0;u[N-1][j]=0;}

            max_out << maxU(u, N) << "\n";
        }
        printf("Decay study done (beta=%.2f).\n", beta);
    }
}

int main()
{
    const double kappa0 = 0.20; // Roll 23
    const double beta   = 0.40; // Roll 23
    const double r0     = 0.40; // diffusion number

    // Nx=50 means 50 intervals => 51 points (0,1/50,...,50/50=1)
    // So N=51
    const int    Ngrid  = 51;
    const double dx     = 1.0 / 50.0; // = 1/Nx
    const double dy     = dx;
    const double dt     = r0 * dx * dx / kappa0;
    const double T      = 0.2;
    int          Nt     = (int)ceil(T / dt);
    double       final_t = Nt * dt;

    printf("Roll 23 — ADI Variable-Coefficient Diffusion\n");
    printf("kappa0=%.2f, beta=%.2f, N=%d (grid pts), dx=%.6f\n", kappa0, beta, Ngrid, dx);
    printf("r0=%.4f, dt=%.6f, Nt=%d, T=%.2f, final_t=%.6f\n\n", r0, dt, Nt, T, final_t);

    vector<double> x = makeGrid(Ngrid);
    vector<double> y = makeGrid(Ngrid);

    // Write kappa field
    writeKappa(x, y, Ngrid, kappa0, beta);

    // Output times: 0, T/3, 2T/3, T
    vector<double> output_times = {0.0, T/3.0, 2.0*T/3.0, T};

    // ── Main ADI solve (beta=0.4) ──
    vector<vector<double>> u(Ngrid, vector<double>(Ngrid, 0.0));
    setIC(u, x, y, Ngrid);
    printf("Solving ADI (variable kappa, beta=%.2f)...\n", beta);
    solveADI(u, x, y, dt, dx, dy, kappa0, beta, Ngrid, Nt, output_times, "ADI");
    printf("\n");

    // ── Convergence study (beta=0, N=25,50,100) ──
    printf("Running convergence study (beta=0)...\n");
    convergenceStudy(kappa0, r0);
    printf("Convergence data written to data/convergence.txt\n\n");

    // ── Decay study (beta=0, 0.4, 0.8) ──
    printf("Running decay study...\n");
    vector<double> betas = {0.0, 0.4, 0.8};
    decayStudy(kappa0, r0, T, betas, x, y, Ngrid);
    printf("\nProgram completed.\n");
    return 0;
}