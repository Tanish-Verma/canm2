#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sys/stat.h>
#include <sys/types.h>

// exact solution: T(x,t) = erf((x - 0.5) / (2 * sqrt(t)))
double exact_solution(double x, double t)
{
  return erf((x - 0.5) / (2.0 * sqrt(t)));
}

void grid(int nx, double xst, double xen, double *x, double *dx)
{
  int i;

  // uniform mesh, grid points placed at cell edges
  *dx = (xen - xst) / (double)(nx - 1);
  for (i = 0; i < nx; i++)
    x[i] = xst + i * (*dx);
}

void enforce_bcs(int nx, double *x, double *T)
{
  // Dirichlet BCs: T(0) = -1, T(1) = 1
  T[0]      = -1.0;
  T[nx - 1] =  1.0;
}

void set_initial_condition(int nx, double *x, double *T, double tst)
{
  int i;

  // initialize from exact solution at t = tst
  for (i = 0; i < nx; i++)
    T[i] = exact_solution(x[i], tst);

  // enforce boundary conditions on top
  enforce_bcs(nx, x, T);
}

void get_rhs(int nx, double dx, double *x, double *T, double *rhs)
{
  int i;
  double kappa = 1.0;

  // 2nd-order central difference: d2T/dx2 ~ (T[i+1] - 2*T[i] + T[i-1]) / dx^2
  // interior points only; boundaries are fixed by BCs
  rhs[0]      = 0.0;
  rhs[nx - 1] = 0.0;
  for (i = 1; i < nx - 1; i++)
    rhs[i] = kappa * (T[i + 1] - 2.0 * T[i] + T[i - 1]) / (dx * dx);
}

void timestep_Euler(int nx, double dt, double dx, double *x, double *T, double *rhs)
{
  int i;

  // calculate rhs at current time level
  get_rhs(nx, dx, x, T, rhs);

  // forward Euler update: T^(n+1) = T^n + dt * rhs
  for (i = 1; i < nx - 1; i++)
    T[i] = T[i] + dt * rhs[i];

  // re-apply BCs after update
  enforce_bcs(nx, x, T);
}

void output_soln(int nx, int it, double tcurr, double *x, double *T)
{
  int i;
  FILE *fp;
  char fname[100];

  sprintf(fname, "output2/output_%04d.dat", it);
  fp = fopen(fname, "w");
  for (i = 0; i < nx; i++)
    fprintf(fp, "%lf %lf\n", x[i], T[i]);
  fclose(fp);

  printf(" > Written solution at it = %d, t = %lf to %s\n", it, tcurr, fname);
}

double get_l2err_norm(int nx, double *T, double *Tex)
{
  int i;
  double l2err = 0.0, val;

  for (i = 0; i < nx; i++)
  {
    val = T[i] - Tex[i];
    l2err += val * val;
  }
  l2err = l2err / (double)nx;
  l2err = sqrt(l2err);

  return l2err;
}

int main()
{
  int nx;
  double *x, *T, *rhs, *Tex;
  double tst, ten, xst, xen, dx, dt, tcurr;
  int i, it, num_time_steps, it_print;
  FILE *fp;

  // ensure output directory exists
  mkdir("output2", 0777);

  // read inputs
  fp = fopen("input2.in", "r");
  fscanf(fp, "%d\n", &nx);
  fscanf(fp, "%lf %lf\n", &xst, &xen);
  fscanf(fp, "%lf %lf\n", &tst, &ten);
  fclose(fp);

  printf("Inputs are: nx=%d xst=%lf xen=%lf tst=%lf ten=%lf\n", nx, xst, xen, tst, ten);

  x   = (double *)calloc(nx, sizeof(double));
  T   = (double *)calloc(nx, sizeof(double));
  rhs = (double *)calloc(nx, sizeof(double));
  Tex = (double *)calloc(nx, sizeof(double));

  grid(nx, xst, xen, x, &dx);              // initialize the grid

  set_initial_condition(nx, x, T, tst);    // initial condition from exact soln at t=tst

  dt = 2.9121e-4;
//   dt = 0.000001;
  num_time_steps = (int)((ten - tst) / dt) + 1;  // +1 to ensure we reach ten
  it_print = num_time_steps / 10;
  if (it_print < 1) it_print = 1;

  printf("dx = %lf, dt = %lf, r = kappa*dt/dx^2 = %lf\n", dx, dt, dt / (dx * dx));
  printf("num_time_steps = %d\n", num_time_steps);

  // time stepping loop
  for (it = 0; it < num_time_steps; it++)
  {
    tcurr = tst + (double)it * dt;

    timestep_Euler(nx, dt, dx, x, T, rhs);

    if (it % it_print == 0)
      output_soln(nx, it, tcurr, x, T);
  }

  // compute exact solution at final time and L2 error
  for (i = 0; i < nx; i++)
    Tex[i] = exact_solution(x[i], ten);

  double l2err = get_l2err_norm(nx, T, Tex);
  printf("\nnx = %d, L2 error = %9.5e\n", nx, l2err);

  // append error to file for convergence study
  fp = fopen("output2/error2.dat", "a");
  fprintf(fp, "%d %9.5e\n", nx, l2err);
  fclose(fp);

  free(rhs);
  free(Tex);
  free(T);
  free(x);

  return 0;
}