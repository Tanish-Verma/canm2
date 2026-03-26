#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#define BETA 3.0

// const double PI = 4.0*atan(1.0);
const double PI = 3.14159265358979323846;

double xu_prime(double xs, double L, double beta)
{
    double log_term = log((beta + 1.0) / (beta - 1.0));
    double L_beta = L * beta;
    double L_minus_2xs = L - 2.0 * xs;
    // Denominator base: (beta*L)^2 - (L - 2xs)^2
    double denom_base = (L_beta * L_beta) - (L_minus_2xs * L_minus_2xs);

    // dxu/dxs
    double first_deriv = (2.0 * L * L_beta) / (denom_base * log_term);
    return first_deriv;
}

double xu_double_prime(double xs, double L, double beta)
{
    double log_term = log((beta + 1.0) / (beta - 1.0));
    double L_beta = L * beta;
    double L_minus_2xs = L - 2.0 * xs;
    double denom_base = (L_beta * L_beta) - (L_minus_2xs * L_minus_2xs);

    // d2xu/dxs2
    // Note: The -8.0 coefficient correctly accounts for the inner derivative 
    // of the squared term in the denominator.
    double second_deriv = (-8.0 * L * L_beta * L_minus_2xs) / (pow(denom_base, 2) * log_term);
    return second_deriv;
}

double xs(double xu, double L, double beta)
{

  double Lambda = (1.0 + beta) / (beta - 1.0);
  double R = pow(Lambda, (2.0 * (xu / L) - 1.0));

  return L * ((1.0 + beta) * R - beta + 1.0) / (2.0 * (1.0 + R));
}
void grid(int nx, double xst, double xen, double *x, double *dx)
{
  int i;
  double dxunif;

  // uniform mesh for now;
  // can use stretching factors to place x nodes later
  dxunif = (xen - xst) / (double)(nx - 1);
  double Lx = xen - xst;

  // populate x[i] s
  for (i = 0; i < nx; i++)
    x[i] = xs(xst + i * dxunif, Lx, BETA); // stretch factor beta = BETA
  // x[i] = xs(i / (nx - 1.0), Lx, BETA); // stretch factor beta = BETA

  // dx[i] s are spacing between adjacent xs
  for (i = 0; i < nx - 1; i++)
    dx[i] = x[i + 1] - x[i];

  //// debug -- print x
  printf("--x--\n");
  for (i = 0; i < nx; i++)
  {
    printf("%d %lf xu = %lf\n", i, x[i], xst + i * dxunif);
  }

  // //// debug -- print dx
  // printf("--dx--\n");
  // for (i = 0; i < nx - 1; i++)
  //   printf("%d %lf\n", i, dx[i]);
}

void set_initial_guess(int nx, int ny, double *x, double *y, double **T)
{
  int i, j;

  for (i = 0; i < nx; i++)
    for (j = 0; j < ny; j++)
      T[i][j] = 0;
}

void calc_diffusivity(int nx, int ny, double *x, double *y, double **T, double **kdiff)
{
  int i, j;
  double A = 1.0, B = 0.0;

  for (i = 0; i < nx; i++)
    for (j = 0; j < ny; j++)
    {
      kdiff[i][j] = A;
    }
}

void calc_sources(int nx, int ny, double hx,double hy,double *x, double *y, double *dx, double *dy, double **T, double **b)
{
  int i, j;

  // calculate source (may be dependent on T)
  for (i = 0; i < nx; i++)
    for (j = 0; j < ny; j++)
    {
      b[i][j] = 0;                        // set source function here
      b[i][j] = -b[i][j] * dx[0] * dx[0]; // only works for uniform mesh
    }
}

void set_boundary_conditions(int nx, int ny, double *x, double *y, double *dx, double *dy, double **T, double **bcleft, double **bcrght, double **bctop, double **bcbot)
{

  int i, j;
  double bcE, bcW, bcN, bcS, kw, ke, kn, ks;

  // left boundary -- Homogeneous Dirichlet
  // p(y) * dTdx + q(y) * T  = r(y)
  for (j = 0; j < ny; j++)
  {
    bcleft[j][0] = 0.0;   // p(y)
    bcleft[j][1] = 1.0;   // q(y)
    bcleft[j][2] = 300.0; // r(y)
  }

  // right boundary -- Homogeneous Dirichlet
  // p(y) * dTdx + q(y) * T = r(y)
  for (j = 0; j < ny; j++)
  {
    bcrght[j][0] = 0.0;                                // p(y)
    bcrght[j][1] = 1.0;                                // q(y)
    // bcrght[j][2] = 300 * (1 + 0.2 * (2 * y[j] - 0.5)); // r(y)
    bcrght[j][2] = 300 * (1 + 0.2 * (y[j]/y[ny-1] - 0.5)); // r(y)
  }

  // bottom boundary -- Homogeneous Neumann
  // p(x) * dTdy + q(x) * T = r(x)
  for (i = 0; i < nx; i++)
  {
    bcbot[i][0] = 1.0; // p(x)
    bcbot[i][1] = 0.0; // q(x)
    bcbot[i][2] = 0.0; // r(x)
  }

  // top boundary -- Homogeneous Neumann
  // p(x) * dTdy + q(x) * T = r(x)
  for (j = 0; j < nx; j++)
  {
    bctop[j][0] = 1.0; // p(x)
    bctop[j][1] = 0.0; // q(x)
    bctop[j][2] = 0.0; // r(x)
  }
}

void get_coeffs(int nx, int ny, double hx, double hy, double *x, double *y, double *dx, double *dy, double **aP, double **aE, double **aW, double **aN, double **aS, double **b, double **T, double **kdiff, double **bcleft, double **bcrght, double **bctop, double **bcbot)
{

  int i, j;
  double pj, qj, rj;
  double hxhy_sq;
  // calculate diffusivity at [x, y]; may be dependent on T
  calc_diffusivity(nx, ny, x, y, T, kdiff);

  // calculate sources Su, Sp at [x, y]; may be dependent on T
  calc_sources(nx, ny, hx,hy,x, y, dx, dy, T, b);

  // populate values in BC arrays
  set_boundary_conditions(nx, ny, x, y, dx, dy, T, bcleft, bcrght, bctop, bcbot);

  // start populating the coefficients
  hxhy_sq = (hx * hx) / (hy * hy);
  double Lx = x[nx - 1] - x[0];
  double Ly = y[ny - 1] - y[0];

  // ------ Step 1 :: interior points ------
  for (i = 1; i < nx - 1; i++)
    for (j = 1; j < ny - 1; j++)
    {
      double Jx = xu_prime(x[i], Lx, BETA);
      double Jy = xu_prime(y[j], Ly, BETA);
      aP[i][j] = 2 * (Jx * Jx + hxhy_sq * Jy * Jy) * kdiff[i][j];
      aE[i][j] = -kdiff[i][j] * (Jx * Jx + 0.5 * hx * xu_double_prime(x[i], Lx, BETA));
      aW[i][j] = -kdiff[i][j] * (Jx * Jx - 0.5 * hx * xu_double_prime(x[i], Lx, BETA));
      aN[i][j] = -kdiff[i][j] * (hxhy_sq * Jy * Jy + hxhy_sq * 0.5 * hy * xu_double_prime(y[j], Ly, BETA));
      aS[i][j] = -kdiff[i][j] * (hxhy_sq * Jy * Jy - hxhy_sq * 0.5 * hy * xu_double_prime(y[j], Ly, BETA));
    }

  // ------ Step 2 :: left boundary ----------
  i = 0;
  double Jx = xu_prime(x[i], Lx, BETA);
  double Jxx = xu_double_prime(x[i], Lx, BETA);
  double E_int = (Jx * Jx + hx * Jxx / 2.0);
  double W_int = (Jx * Jx - hx * Jxx / 2.0);

  for (j = 0; j < ny; j++)
  {
    pj = bcleft[j][0];
    qj = bcleft[j][1];
    rj = bcleft[j][2];
    double Jy = xu_prime(y[j], Ly, BETA);
    double Jyy = xu_double_prime(y[j], Ly, BETA);
    double ptilda = pj * Jx;
    double N_int = (hxhy_sq * Jy * Jy + hxhy_sq * 0.5 * hy * Jyy);
    double S_int = (hxhy_sq * Jy * Jy - hxhy_sq * 0.5 * hy * Jyy);
    double P_int = Jx * Jx + hxhy_sq * Jy * Jy;

    aP[i][j] = 2.0 * kdiff[i][j] * (P_int * ptilda - hx * W_int * qj);
    aE[i][j] = -kdiff[i][j] * (E_int + W_int) * ptilda;
    aW[i][j] = 0.0;
    aN[i][j] = -N_int * kdiff[i][j] * ptilda;
    aS[i][j] = -S_int * kdiff[i][j] * ptilda;
    b[i][j] = b[i][j] * ptilda - 2.0 * hx * rj * kdiff[i][j] * W_int;
  }
  // ------ Step 2 :: left boundary done ---

  // ------ Step 3 :: right boundary ----------
  {
    i = nx - 1;
    Jx = xu_prime(x[i], Lx, BETA);
    Jxx = xu_double_prime(x[i], Lx, BETA);
    E_int = (Jx * Jx + hx * Jxx / 2.0);
    W_int = (Jx * Jx - hx * Jxx / 2.0);
    for (j = 0; j < ny; j++)
    {
      pj = bcrght[j][0];
      qj = bcrght[j][1];
      rj = bcrght[j][2];
      double Jy = xu_prime(y[j], Ly, BETA);
      double Jyy = xu_double_prime(y[j], Ly, BETA);
      double ptilda = pj * Jx;
      double N_int = (hxhy_sq * Jy * Jy + hxhy_sq * 0.5 * hy * Jyy);
      double S_int = (hxhy_sq * Jy * Jy - hxhy_sq * 0.5 * hy * Jyy);
      double P_int = Jx * Jx + hxhy_sq * Jy * Jy;

      aP[i][j] = 2.0 * kdiff[i][j] * (P_int * ptilda + hx * E_int * qj);
      aW[i][j] = -kdiff[i][j] * (E_int + W_int) * ptilda;
      aE[i][j] = 0.0;
      aN[i][j] = -N_int * kdiff[i][j] * ptilda;
      aS[i][j] = -S_int * kdiff[i][j] * ptilda;
      b[i][j] = b[i][j] * ptilda + 2.0 * hx * rj * kdiff[i][j] * E_int;
    }
  }
  // ------ Step 3 :: right boundary done ---

  // ------ Step 4 :: bottom boundary ----------
  j = 0;
  double Jy = xu_prime(y[j], Ly, BETA);
  double Jyy = xu_double_prime(y[j], Ly, BETA);
  double N_int = (hxhy_sq * Jy * Jy + hxhy_sq * 0.5 * hy * Jyy);
  double S_int = (hxhy_sq * Jy * Jy - hxhy_sq * 0.5 * hy * Jyy);

  for (i = 1; i < nx - 1; i++)
  {
    pj = bcbot[i][0];
    qj = bcbot[i][1];
    rj = bcbot[i][2];
    Jx = xu_prime(x[i], Lx, BETA);
    Jxx = xu_double_prime(x[i], Lx, BETA);
    double ptilda = pj * Jy;
    double E_int = (Jx * Jx + hx * Jxx / 2.0);
    double W_int = (Jx * Jx - hx * Jxx / 2.0);
    double P_int = Jx * Jx + hxhy_sq * Jy * Jy;

    aP[i][j] = 2.0 * kdiff[i][j] * (P_int * ptilda - hy * S_int * qj);
    aW[i][j] = -kdiff[i][j] * ptilda * W_int;
    aE[i][j] = -kdiff[i][j] * ptilda * E_int;
    aN[i][j] = -kdiff[i][j] * ptilda * (N_int + S_int);
    aS[i][j] = 0.0;
    b[i][j] = (b[i][j] * ptilda) - 2 * hy * rj * kdiff[i][j] * S_int;
  }
  // ------ Step 4 :: bottom boundary done ---

  // ------ Step 5 :: top boundary ----------
  j = ny - 1;
  Jy = xu_prime(y[j], Ly, BETA);
  Jyy = xu_double_prime(y[j], Ly, BETA);
  N_int = (hxhy_sq * Jy * Jy + hxhy_sq * 0.5 * hy * Jyy);
  S_int = (hxhy_sq * Jy * Jy - hxhy_sq * 0.5 * hy * Jyy);
  for (i = 1; i < nx - 1; i++)
  {
    pj = bctop[i][0];
    qj = bctop[i][1];
    rj = bctop[i][2];

    Jx = xu_prime(x[i], Lx, BETA);
    Jxx = xu_double_prime(x[i], Lx, BETA);
    double ptilda = pj * Jy;
    double E_int = (Jx * Jx + hx * Jxx / 2.0);
    double W_int = (Jx * Jx - hx * Jxx / 2.0);
    double P_int = Jx * Jx + hxhy_sq * Jy * Jy;

    aP[i][j] = 2.0 * kdiff[i][j] * (P_int * ptilda + hy * N_int * qj);
    aW[i][j] = -kdiff[i][j] * ptilda * W_int;
    aE[i][j] = -kdiff[i][j] * ptilda * E_int;
    aN[i][j] = 0.0;
    aS[i][j] = -kdiff[i][j] * ptilda * (N_int + S_int);
    b[i][j] = b[i][j] * ptilda + 2.0 * hy * rj * kdiff[i][j] * N_int;
  }
  // ------ Step 5 :: top boundary done ---

  // debug
  // for (i = 0; i < nx; i++)
  //   for (j = 0; j < ny; j++)
  //   {
  //     printf("%d %d %lf %lf %lf %lf %lf %lf\n", i, j, aP[i][j], aE[i][j], aW[i][j], aN[i][j], aS[i][j], b[i][j]);
  //   }
}

double get_max_of_array(int nx, int ny, double **arr)
{
  int i, j;
  double arrmax, val;

  arrmax = arr[0][0];
  for (i = 0; i < nx; i++)
    for (j = 0; j < ny; j++)
    {
      val = fabs(arr[i][j]);
      if (arrmax < val)
        arrmax = val;
    }
  return arrmax;
}

double get_l2err_norm(int nx, int ny, double **arr1, double **arr2)
{
  double l2err = 0.0, val;
  int i, j;

  for (i = 0; i < nx; i++)
    for (j = 0; j < ny; j++)
    {
      val = arr1[i][j] - arr2[i][j];
      l2err += val * val;
    }
  // printf("l2err = %lf\n", l2err);
  l2err = l2err / ((double)(nx * ny));
  l2err = sqrt(l2err);

  return l2err;
}

void solve_gssor(int nx, int ny, double **aP, double **aE, double **aW, double **aN, double **aS, double **b, double **T, double **Tpad, double **Tpnew, int max_iter, double tol, double relax_T)
{

  int i, j, ip, jp, iter;
  double l2err, arrmax1, arrmax2, err_ref, T_gs, rel_err;

  // copy to padded array
  for (i = 0; i < nx; i++)
    for (j = 0; j < ny; j++)
    {
      ip = i + 1;
      jp = j + 1;
      Tpad[ip][jp] = T[i][j];
    }

  // now perform iterations
  for (iter = 0; iter < max_iter; iter++)
  {
    // update Tpnew
    for (i = 0; i < nx; i++)
      for (j = 0; j < ny; j++)
      {
        ip = i + 1;
        jp = j + 1;
        T_gs = (b[i][j] - aE[i][j] * Tpad[ip + 1][jp] - aW[i][j] * Tpnew[ip - 1][jp] -
                aN[i][j] * Tpad[ip][jp + 1] - aS[i][j] * Tpnew[ip][jp - 1]) /
               aP[i][j];
        Tpnew[ip][jp] = (1.0 - relax_T) * Tpad[ip][jp] + relax_T * T_gs;
        // printf("+++%d %d %e\n", ip, jp, Tpnew[ip][jp]);
      }

    // check for convergence
    l2err = get_l2err_norm(nx + 2, ny + 2, Tpad, Tpnew);
    arrmax1 = get_max_of_array(nx + 2, ny + 2, Tpad);
    arrmax2 = get_max_of_array(nx + 2, ny + 2, Tpnew);
    err_ref = fmax(arrmax1, arrmax2);
    err_ref = fmax(err_ref, 1.0e-6);
    rel_err = l2err / err_ref;
    // printf("   > %d %9.5e  %9.5e  %9.5e\n", iter, l2err, err_ref, rel_err);
    if (rel_err < tol)
      break;

    // prepare for next iteration
    for (i = 0; i < nx; i++)
      for (j = 0; j < ny; j++)
      {
        ip = i + 1;
        jp = j + 1;
        Tpad[ip][jp] = Tpnew[ip][jp];
      }
  }

  // copy from padded array
  for (i = 0; i < nx; i++)
    for (j = 0; j < ny; j++)
    {
      ip = i + 1;
      jp = j + 1;
      T[i][j] = Tpad[ip][jp];
    }

  printf(" > Fin: %d %9.5e  %9.5e  %9.5e\n\n", iter, l2err, err_ref, rel_err);
}

void get_exact_soln(int nx, int ny, double *x, double *y, double **Tex)
{
  int i, j, k;
  double a = x[nx - 1];
  double b = y[ny - 1];
  double T0 = 300.0;

  for (i = 0; i < nx; i++)
  {
    for (j = 0; j < ny; j++)
    {
      double sum = 0.0;
      for (k = 0; k < 50; k++)
      {
        double n = 2 * k + 1;
        double lam = n * PI / b;
        sum += (1.0 / (n * n * PI * PI)) * (sinh(lam * x[i]) / sinh(lam * a)) * cos(lam * y[j]);
      }

      Tex[i][j] = T0 * (1.0 - 0.8 * sum);
    }
  }
}

void output_soln(int nx, int ny, int iter, double *x, double *y, double **T, double **Tex)
{
  int i, j;
  FILE *fp;
  char fname[100];

  sprintf(fname, "output1/T_xy_%03d_%03d_%04d.dat", nx, ny, iter);

  fp = fopen(fname, "w");
  for (i = 0; i < nx; i++)
    for (j = 0; j < ny; j++)
      fprintf(fp, "%lf %lf %lf %lf\n", x[i], y[j], T[i][j], Tex[i][j]);
  fclose(fp);

  printf(" > Done writing solution for stamp = %d to file %s\n\n", iter, fname);
}

double get_sor(int nx, int ny)
{
  double omega;
  // assuming uniform grid for calculating optimal SOR factor
  double r = (cos(PI / nx) + cos(PI / ny));
  omega = 2.0 / (1.0 + sqrt(1.0 - r * r / 4.0));
  return omega;
}

int main()
{

  int nx, ny;
  double *x, *dx, *y, *dy;
  double **aP, **aE, **aW, **aN, **aS, **b, **Sp;
  double **T, **Tex, **kdiff, **wrk1, **wrk2;
  double xst, xen, yst, yen;
  double **bcleft, **bcrght, **bctop, **bcbot;
  int i, j, max_iter;
  double relax_T, tol, l2err;
  FILE *fp;

  // read inputs
  fp = fopen("input1.in", "r");
  fscanf(fp, "%d %d\n", &nx, &ny);
  fscanf(fp, "%lf %lf\n", &xst, &xen);
  fscanf(fp, "%lf %lf\n", &yst, &yen);
  fclose(fp);

  printf("Inputs are: %d %d %lf %lf %lf %lf\n", nx, ny, xst, xen, yst, yen);

  // allocate memory
  printf("\n > Allocating Memory -- \n");
  x = (double *)calloc(nx, sizeof(double));      // grid points
  dx = (double *)calloc(nx - 1, sizeof(double)); // spacing betw grid points

  y = (double *)calloc(ny, sizeof(double));      // grid points
  dy = (double *)calloc(ny - 1, sizeof(double)); // spacing betw grid points

  printf("   >> Done allocating 1D arrays -- \n");

  // allocate 2D arrays dynamically
  // -- for T --
  T = (double **)calloc(nx, sizeof(double *));
  for (i = 0; i < nx; i++)
    T[i] = (double *)calloc(ny, sizeof(double));

  // -- for Tex --
  Tex = (double **)calloc(nx, sizeof(double *));
  for (i = 0; i < nx; i++)
    Tex[i] = (double *)calloc(ny, sizeof(double));

  // -- for aP --
  aP = (double **)calloc(nx, sizeof(double *));
  for (i = 0; i < nx; i++)
    aP[i] = (double *)calloc(ny, sizeof(double));

  // -- for aE --
  aE = (double **)calloc(nx, sizeof(double *));
  for (i = 0; i < nx; i++)
    aE[i] = (double *)calloc(ny, sizeof(double));

  // -- for aW --
  aW = (double **)calloc(nx, sizeof(double *));
  for (i = 0; i < nx; i++)
    aW[i] = (double *)calloc(ny, sizeof(double));

  // -- for aN --
  aN = (double **)calloc(nx, sizeof(double *));
  for (i = 0; i < nx; i++)
    aN[i] = (double *)calloc(ny, sizeof(double));

  // -- for aS --
  aS = (double **)calloc(nx, sizeof(double *));
  for (i = 0; i < nx; i++)
    aS[i] = (double *)calloc(ny, sizeof(double));

  // -- for b --
  b = (double **)calloc(nx, sizeof(double *));
  for (i = 0; i < nx; i++)
    b[i] = (double *)calloc(ny, sizeof(double));

  // left boundary condition
  bcleft = (double **)calloc(ny, sizeof(double *));
  for (i = 0; i < ny; i++)
    bcleft[i] = (double *)calloc(3, sizeof(double));

  // right boundary condition
  bcrght = (double **)calloc(ny, sizeof(double *));
  for (i = 0; i < ny; i++)
    bcrght[i] = (double *)calloc(3, sizeof(double));

  // bottom boundary condition
  bctop = (double **)calloc(nx, sizeof(double *));
  for (i = 0; i < nx; i++)
    bctop[i] = (double *)calloc(3, sizeof(double));

  // top boundary condition
  bcbot = (double **)calloc(nx, sizeof(double *));
  for (i = 0; i < nx; i++)
    bcbot[i] = (double *)calloc(3, sizeof(double));

  // -- for kdiff --
  kdiff = (double **)calloc(nx, sizeof(double *));
  for (i = 0; i < nx; i++)
    kdiff[i] = (double *)calloc(ny, sizeof(double));

  // -- for work arrays wrk1, wrk2 --
  wrk1 = (double **)calloc(nx + 2, sizeof(double *));
  for (i = 0; i < nx + 2; i++)
    wrk1[i] = (double *)calloc(ny + 2, sizeof(double));

  wrk2 = (double **)calloc(nx + 2, sizeof(double *));
  for (i = 0; i < nx + 2; i++)
    wrk2[i] = (double *)calloc(ny + 2, sizeof(double));

  printf("   >> Done allocating 2D arrays -- \n");
  printf(" > Done allocating memory -------- \n");

  // initialize the grid
  grid(nx, xst, xen, x, dx); // -- along x --
  grid(ny, yst, yen, y, dy); // -- along y --
  printf("\n > Done setting up grid ---------- \n");
  double hx = (xen - xst) / (nx - 1);
  double hy = (yen - yst) / (ny - 1);
  set_initial_guess(nx, ny, x, y, T); // initial condition
  printf("\n > Done setting up initial guess -- \n");

  // // ---
  get_coeffs(nx, ny, hx, hy, x, y, dx, dy,    // grid vars
             aP, aE, aW, aN, aS, b, T, kdiff, // coefficients
             bcleft, bcrght, bctop, bcbot);   // BC vars
  printf("\n > Done calculating coeffs ----- \n");

  printf("\n > Solving for T ------------- \n\n");
  max_iter = 100000;
  tol = 1.0e-10;
  relax_T = get_sor(nx, ny); // approximating optimal relaxation factor for SOR
  // relax_T = 1.5;
  clock_t start_time = clock();
  solve_gssor(nx, ny, aP, aE, aW, aN, aS, b, T, wrk1, wrk2, max_iter, tol, relax_T);
  clock_t end_time = clock();

  // // ---
  printf(" > Done solving for T ------------- \n\n");

  double runtime = ((double)(end_time - start_time)) / CLOCKS_PER_SEC;
  int minutes = (int)(runtime / 60);
  int seconds = (int)(runtime) % 60;
  int milliseconds = (int)((runtime - (int)runtime) * 1000);
  printf(" > Total runtime for solve_gssor: %d minutes, %d seconds, %d milliseconds\n\n", minutes, seconds, milliseconds);

  get_exact_soln(nx, ny, x, y, Tex);
  output_soln(nx, ny, 0, x, y, T, Tex);

  l2err = get_l2err_norm(nx, ny, T, Tex);
  printf(" > %d %d %9.5e\n", nx, ny, l2err);

  FILE *error_fp = fopen("output1/error.dat", "a");
  fprintf(error_fp, "%d %d %9.5e\n", nx, ny, l2err);
  fclose(error_fp);

  // free memory
  // ----1D arrays ---
  free(y);
  free(dy);
  free(x);
  free(dx);
  // --- Done 1D arrays ---

  // ----2D arrays ---
  for (i = 0; i < nx; i++)
    free(T[i]);
  free(T);
  for (i = 0; i < nx; i++)
    free(Tex[i]);
  free(Tex);
  for (i = 0; i < nx; i++)
    free(b[i]);
  free(b);
  for (i = 0; i < nx; i++)
    free(aP[i]);
  free(aP);
  for (i = 0; i < nx; i++)
    free(aE[i]);
  free(aE);
  for (i = 0; i < nx; i++)
    free(aW[i]);
  free(aW);
  for (i = 0; i < nx; i++)
    free(aN[i]);
  free(aN);
  for (i = 0; i < nx; i++)
    free(aS[i]);
  free(aS);
  for (i = 0; i < nx; i++)
    free(kdiff[i]);
  free(kdiff);
  for (i = 0; i < ny; i++)
    free(bcleft[i]);
  free(bcleft);
  for (i = 0; i < ny; i++)
    free(bcrght[i]);
  free(bcrght);
  for (i = 0; i < nx; i++)
    free(bctop[i]);
  free(bctop);
  for (i = 0; i < nx; i++)
    free(bcbot[i]);
  free(bcbot);
  for (i = 0; i < nx + 2; i++)
    free(wrk1[i]);
  free(wrk1);
  for (i = 0; i < nx + 2; i++)
    free(wrk2[i]);
  free(wrk2);
  // --- Done 2D arrays ---
  printf("\n > Done freeing up memory --------- \n");
  return 0;
}
