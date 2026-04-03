#include <bits/stdc++.h>
using namespace std;

#define nu 2.14

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

void writesolutiontoFile(const vector<double> &u, const vector<double> &x, double t, const string &method)
{
    ostringstream name;
    name << "data/" << method << "_T_" << t << ".txt";
    ofstream out(name.str());
    for (size_t i = 0; i < u.size() && i < x.size(); ++i)
    {
        out << x[i] << " " << u[i] << "\n";
    }
}

vector<double> solveTDMA(const vector<double> &a, const vector<double> &b, const vector<double> &c, const vector<double> &d)
{
    int N = d.size();
    vector<double> c_star(N, 0.0);
    vector<double> d_star(N, 0.0);
    vector<double> sol(N, 0.0);

    // Forward elimination
    c_star[0] = c[0] / b[0];
    d_star[0] = d[0] / b[0];

    for (int i = 1; i < N; i++)
    {
        double m = 1.0 / (b[i] - a[i] * c_star[i - 1]);
        if (i < N - 1)
        {
            c_star[i] = c[i] * m;
        }
        d_star[i] = (d[i] - a[i] * d_star[i - 1]) * m;
    }

    // Backward substitution
    sol[N - 1] = d_star[N - 1];
    for (int i = N - 2; i >= 0; i--)
    {
        sol[i] = d_star[i] - c_star[i] * sol[i + 1];
    }

    return sol;
}

vector<double> solvePeriodicTDMA(double a, double b, double c, const vector<double> &d)
{
    int N = d.size();
    vector<double> x(N, 0.0);

    // 1. Choose gamma (and map alpha/beta to the periodic corner elements)
    double gamma = -b;
    double alpha = c; // Bottom-left corner
    double beta = a;  // Top-right corner

    // 2. Form T vectors (lower, main, and upper diagonals)
    vector<double> T_a(N, a); 
    vector<double> T_b(N, b); 
    vector<double> T_c(N, c); 

    // Modify the first and last main diagonal elements
    T_b[0] = b - gamma;
    T_b[N - 1] = b - (alpha * beta) / gamma;

    // 3. Form u = [gamma, 0, ..., 0, alpha]^T  <--- FIX 1: This must be alpha!
    vector<double> u(N, 0.0);
    u[0] = gamma;
    u[N - 1] = alpha; 

    // 4. Solve Ty = d using the external standard TDMA
    vector<double> y = solveTDMA(T_a, T_b, T_c, d);

    // 5. Solve Tz = u using the external standard TDMA
    vector<double> z = solveTDMA(T_a, T_b, T_c, u);

    // 6. Compute scalar: rho  <--- FIX 2: These must use beta!
    double rho = (y[0] + (beta / gamma) * y[N - 1]) /
                 (1.0 + z[0] + (beta / gamma) * z[N - 1]);

    // 7. x = y - rho*z -> final solution
    for (int i = 0; i < N; ++i)
    {
        x[i] = y[i] - rho * z[i];
    }

    return x;
}
void solveMacCormack(vector<double> &u, const vector<double> &x, double dt, int Nt, int Nx, const vector<double> &output_times)
{
    vector<double> u_new(Nx, 0.0);
    vector<double> u_pred(Nx, 0.0); // Predictor array

    writesolutiontoFile(u, x, 0.0, "MacCormack");

    int next_output_idx = 1;

    for (int n = 1; n <= Nt; n++)
    {
        // 1. Predictor Step (Forward Space)
        for (int i = 0; i < Nx - 1; ++i)
        {
            u_pred[i] = u[i] - nu * (u[i + 1] - u[i]);
        }
        u_pred[Nx - 1] = u[Nx - 1] - nu * (u[0] - u[Nx - 1]); // Periodic

        // 2. Corrector Step (Backward Space)
        for (int i = 1; i < Nx; ++i)
        {
            u_new[i] = 0.5 * (u[i] + u_pred[i] - nu * (u_pred[i] - u_pred[i - 1]));
        }
        u_new[0] = 0.5 * (u[0] + u_pred[0] - nu * (u_pred[0] - u_pred[Nx - 1])); // Periodic

        u = u_new;

        double current_time = n * dt;
        if (next_output_idx < output_times.size() && current_time >= output_times[next_output_idx])
        {
            writesolutiontoFile(u, x, output_times[next_output_idx], "MacCormack");
            next_output_idx++;
        }
    }
}

void solveCrankNicolson(vector<double> &u, const vector<double> &x, double dt, int Nt, int Nx, const vector<double> &output_times)
{
    writesolutiontoFile(u, x, 0.0, "CrankNicolson");

    int next_output_idx = 1;

    // Matrix coefficients for the LHS: u_i^{n+1} + (nu/4)u_{i+1} - (nu/4)u_{i-1}
    double a_coeff = -nu / 4.0; // Lower diagonal
    double b_coeff = 1.0;       // Main diagonal
    double c_coeff = nu / 4.0;  // Upper diagonal

    for (int n = 1; n <= Nt; n++)
    {
        vector<double> rhs(Nx, 0.0);

        // Calculate the RHS vector (b^n) from the known time step 'u'
        // RHS eq: u_i^n - (nu/4)*(u_{i+1}^n - u_{i-1}^n)
        for (int i = 1; i < Nx - 1; ++i)
        {
            rhs[i] = u[i] - (nu / 4.0) * (u[i + 1] - u[i - 1]);
        }
        // RHS Periodic bounds
        rhs[0] = u[0] - (nu / 4.0) * (u[1] - u[Nx - 1]);
        rhs[Nx - 1] = u[Nx - 1] - (nu / 4.0) * (u[0] - u[Nx - 2]);

        // Call your TDMA solver to get u^{n+1}
        u = solvePeriodicTDMA(a_coeff, b_coeff, c_coeff, rhs);

        double current_time = n * dt;
        if (next_output_idx < output_times.size() && current_time >= output_times[next_output_idx])
        {
            writesolutiontoFile(u, x, output_times[next_output_idx], "CrankNicolson");
            next_output_idx++;
        }
    }
}

int main()
{
    // Extract inputs from file
    FILE *fp = fopen("input.txt", "r");
    if (!fp)
    {
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
    double final_time = Nt * dt;
    double e;
    // Solve using MacCormack method
    printf("Solving using MacCormack method...\n");
    solveMacCormack(u, x, dt, Nt, Nx, output_times);
    e = calculatel2Norm(u, exactSol(x, final_time, c, L), dx);
    printf("L2 Error for MacCormack at t=1.0: %lf\n\n", e);
    printf("MacCormack method completed.\n");

    // Reset initial condition for Crank-Nicolson
    setinitialCondition(u, x, L);

    // Solve using Crank-Nicolson method
    printf("Solving using Crank-Nicolson method...\n");
    solveCrankNicolson(u, x, dt, Nt, Nx, output_times);
    e = calculatel2Norm(u, exactSol(x, final_time, c, L), dx);
    printf("L2 Error for Crank-Nicolson at t=1.0: %lf\n\n", e);
    printf("Crank-Nicolson method completed.\n");

    printf("Program completed.\n");
    return 0;
}