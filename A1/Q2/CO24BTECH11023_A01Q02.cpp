#include <bits/stdc++.h>
using namespace std;

vector<double> makeGrid(double xstart, double xend, int Nx)
{
    vector<double> grid;
    double dx = (xend - xstart) / (Nx);
    for (size_t i = 0; i < Nx; i++)
    {
        grid.push_back(xstart + i * dx);
    }
    // for(double xi : grid)
    // {
    //     printf("%lf ", xi);
    //     printf("\n");
    // }
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

void writesolutiontoFile(const vector<double> &u, const vector<double> &x, double t, const string &dir)
{
    ostringstream name;
    name << dir << "/FTCS_T_" << t << ".txt";
    ofstream out(name.str());
    for (size_t i = 0; i < u.size() && i < x.size(); ++i)
    {
        out << x[i] << " " << u[i] << "\n";
    }
}

void resetFTCSOutputFiles()
{
    system("rm -f data/SetTimesteps/FTCS_T_*.txt");
    system("rm -f data/First50Steps/FTCS_T_*.txt");
    system("rm -f data/FTCS_max_u.txt");
}

void solveFTCS(vector<double> &u, const vector<double> &x, double dt, double dx, double c, double nu, int Nt, int Nx, const vector<double> &output_times)
{
    std::vector<double> u_new(u.size(), 0.0);

    writesolutiontoFile(u, x, 0.0, "data/SetTimesteps");
    writesolutiontoFile(u, x, 0.0, "data/First50Steps");

    int next_output_idx = 1;

    for (int n = 1; n <= Nt; n++)
    {
        for (int i = 1; i < Nx - 1; ++i)
        {
            u_new[i] = u[i] - nu * (u[i + 1] - u[i - 1]) * 0.5;
        }
        u_new[0] = u[0] - nu * (u[1] - u[Nx - 1]) * 0.5;
        u_new[Nx - 1] = u[Nx - 1] - nu * (u[0] - u[Nx - 2]) * 0.5;
        u = u_new;

        double current_time = n * dt;
        if (calculatel2Norm(u, exactSol(x, current_time, c, x.back() - x.front()), dx) > 1.0)
        {
            printf("Warning: L2 norm exceeded 1.0 at time step %d\n", n);
        }
        if (n <= 50)
        {
            writesolutiontoFile(u, x, current_time, "data/First50Steps");
        }
        if (n <= 200)
        {
            double max_u = 0.0;
            for (double val : u)
            {
                max_u = std::max(max_u, std::abs(val));
            }
            // Open in append mode so it keeps adding lines
            ofstream max_out("data/FTCS_max_u.txt", ios_base::app);
            max_out << current_time << " " << max_u << "\n";
        }
        if (next_output_idx < output_times.size() && current_time >= output_times[next_output_idx])
        {
            writesolutiontoFile(u, x, output_times[next_output_idx], "data/SetTimesteps");
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
    int Nx;
    double c, nu;
    double xstart, xend, tstart, tend;
    int first_line_items = fscanf(fp, "%d %lf %lf", &Nx, &c, &nu);
    if (first_line_items < 2)
    {
        printf("Error: First input line must contain at least Nx and c.\n");
        fclose(fp);
        return 1;
    }
    fscanf(fp, "%lf %lf %lf %lf", &xstart, &xend, &tstart, &tend);
    fclose(fp);
    printf("The inputs are Nx = %d, c = %lf, nu = %lf, xstart = %lf, xend = %lf, tstart = %lf, tend = %lf\n\n", Nx, c, nu, xstart, xend, tstart, tend);

    // Build grid and calculate time step
    double dx = (xend - xstart) / (Nx);
    double dt = nu * dx / c;
    // int Nt = ceil((tend - tstart) / dt);
    int Nt = 200;
    double L = xend - xstart;
    vector<double> x = makeGrid(xstart, xend, Nx);

    printf("The values of dx and dt are %lf and %lf respectively\n\n", dx, dt);
    printf("The number of time steps Nt is %d\n\n", Nt);

    // Initialize state and output schedule
    resetFTCSOutputFiles();
    vector<double> u(Nx, 0.0);
    setinitialCondition(u, x, L);
    vector<double> output_times = {0.0, 0.25, 0.5, 0.75, 1.0};
    double e;
    double final_time = Nt * dt;
    // Solve using FTCS and report error
    printf("Solving using FTCS...\n\n");
    solveFTCS(u, x, dt, dx, c, nu, Nt, Nx, output_times);
    printf("FTCS solution completed.\n\n");
    e = calculatel2Norm(u, exactSol(x, final_time, c, L), dx);
    printf("L2 Error for FTCS: %lf\n\n", e);

    printf("Program completed.\n");
    return 0;
}