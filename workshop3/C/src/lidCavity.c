/*============================================================================*\
Solves Navier-Stokes equations for incompressible, laminar, steady flow using
artificial compressibility method on staggered grid.

The governing equations are as follows:
P_t + c^2 div[u] = 0
u_t + u . grad[u] = - grad[P] + nu div[grad[u]]

where P is p/rho and c represents artificial sounds speed.

Lid-Driven Cavity case:
Dimensions : 1x1 m
Grid size  : 128 x 128
Re number  : 100 / 1000 / 5000 / 10000
Tolerance  : 10^-8
CFL        : 0.08
c^2        : 3.8
Grid type  : Arakawa C
Boundary Conditions: u, v -> Dirichlet (as shown below)
                     p    -> Neumann (grad[p] = 0)

                                u=1, v=0
                             ---------------
                            |               |
                        u=0 |               | u=0
                        v=0 |               | v=0
                            |               |
                            |               |
                             ---------------
                                  u=0, v=0
\*============================================================================*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "functions.h"

/* Specify number of grid points in x and y directions */
#define IX 128
#define IY 128

/* Define a max function */
#ifndef max
#   define max(a,b) \
           ({ __auto_type _a = (a); \
               __auto_type _b = (b); \
             _a > _b ? _a : _b; })
#endif

/* Define a function to check if a variable is NAN */
#ifndef isnan
#   define isnan(x) x != x
#endif

/* Main program */
int main (int argc, char *argv[])
{
  double **ubufo, **ubufn, **u, **un, **uc, **utmp;
  double **vbufo, **vbufn, **v, **vn, **vc, **vtmp;
  double **pbufo, **pbufn, **p, **pn, **pc, **ptmp;

  double dx, dy, dt, Re, nu;
  double dtdx, dtdy, dtdxx, dtdyy;

  int i, j, step;
  double ut, ub, ul, ur;
  double vt, vb, vl, vr;
  double err_tot, err_u, err_v, err_p, err_d;

  /* Define first and last interior nodes number for easier iteration */
  const int xlo = 1, xhi = IX - 1, xtot = IX + 1, xm = IX/2;
  const int ylo = 1, yhi = IY - 1, ytot = IY + 1, ym = IY/2;

  /* Simulation parameters */
  /* Best for Re = 100
  const double cfl = 0.15, c2 = 5.0;*/
  /* Best for Re = 1000*/
  const double cfl = 0.20, c2 = 3.8;

  const double tol = 1.0e-8, llid = 1.0, gradp = 0.0;

  /* Getting Reynolds number */
  if(argc <= 1) {
    Re = 1000.0;
  } else {
    Re = atof(argv[1]);
  }
  printf("Re number is set to %d\n", (int) Re);

  /* Boundary condition values */
  ut = 1.0;
  ub = ul = ur = 0.0;
  vt = vb = vl = vr = 0.0;

  /* Create a log file for outputting the residuals */
  FILE *flog;
  flog = fopen("data/residual","w+t");

  /* Compute flow parameters based on inputs */
  dx = llid / (double) (IX - 1);
  dy = llid / (double) (IY - 1);
  dt = cfl * fmin(dx,dy) / ut;
  nu = ut * llid / Re;

  /* Carry out operations that their values do not change in loops */
  dtdx = dt / dx;
  dtdy = dt / dy;
  dtdxx = dt / (dx * dx);
  dtdyy = dt / (dy * dy);

  /* Generate two 2D arrays for storing old and new velocity field
     in x-direction */
  ubufo = array_2d(IX, ytot);
  ubufn = array_2d(IX, ytot);
  /* Define two pointers to the generated buffers for velocity field
     in x-direction */
  u = ubufo;
  un = ubufn;

  /* Generate two 2D arrays for storing old and new velocity field
     in y-direction */
  vbufo = array_2d(xtot, IY);
  vbufn = array_2d(xtot, IY);
  /* Define two pointers to the generated buffers for velocity field
     in y-direction */
  v = vbufo;
  vn = vbufn;

  /* Generate two 2D arrays for storing old and new pressure field*/
  pbufo = array_2d(xtot, ytot);
  pbufn = array_2d(xtot, ytot);
  /* Define two pointers to the generated buffers for pressure field*/
  p = pbufo;
  pn = pbufn;

  /* Apply initial conditions*/
  for (i = xlo; i < xhi; i++) {
    u[i][ytot - 1] = ut;
    u[i][ytot - 2] = ut;
  }

  /* Initialize error and step */
  err_tot = 1.0;
  step = 1;

  /* Start the main loop */
  while (err_tot > tol) {
    /* Solve x-momentum equation for computing u */
    #pragma omp parallel for private(i,j) schedule(auto)
    for (i = xlo; i < xhi; i++)
      for (j = ylo; j < ytot - 1; j++)
        un[i][j] = u[i][j]
                   - dtdx * (pow(avx(u, i+1, j), 2) - pow(avx(u, i, j), 2))
                   - dtdy * (avy(u, i, j+1) * avx(v, i+1, j)
                            - avy(u, i, j) * avx(v, i+1, j-1))
                   - dtdx * (p[i+1][j] - p[i][j])
                   + nu * div2(u, i, j, dtdxx, dtdyy);

    /* Dirichlet boundary conditions */
    for (j = 0; j < ytot; j++) {
      un[0][j] = ul;
      un[xhi][j] = ur;
    }

    for (i = 0; i < IX; i++) {
      un[i][0] = ub - un[i][1];
      un[i][ytot - 1] = 2.0*ut - un[i][ytot - 2];
    }

    /* Solve y-momentum for computing v */
    #pragma omp parallel for private(i,j) schedule(auto)
    for (i = xlo; i < xtot - 1; i++)
      for (j = ylo; j < yhi; j++)
        vn[i][j] = v[i][j]
                   - dtdx * (avy(u, i, j+1) * avx(v, i+1, j)
                            - avy(u, i-1, j+1) * avx(v, i, j))
                   - dtdy * (pow(avy(v, i, j+1), 2) - pow(avy(v, i, j), 2))
                   - dtdy * (p[i][j+1] - p[i][j])
                   + nu * div2(v, i, j, dtdxx, dtdyy);

    /* Dirichlet boundary conditions */
    for (j = 0; j < IY; j++) {
      vn[0][j] = vr - vn[1][j];
      vn[xtot - 1][j] = vl - vn[IX-1][j];
    }

    for (i = 0; i < xtot; i++) {
      vn[i][0] = vb;
      vn[i][yhi] = vt;
    }

    /* Solves continuity equation for computing P */
    #pragma omp parallel for private(i,j) schedule(auto)
    for (i = xlo; i < xtot - 1; i++)
      for (j = ylo; j < ytot - 1; j++)
        pn[i][j] = p[i][j] - c2 * ((un[i][j] - un[i-1][j]) * dtdx
                                   + (vn[i][j] - vn[i][j-1]) * dtdy);

    /* Neumann boundary conditions */
    for (i = 0; i < xtot; i++) {
      pn[i][0] = pn[i][1] - dy * gradp;
      pn[i][ytot - 1] = pn[i][ytot - 2]  - dy * gradp;
    }

    for (j = 0; j < ytot; j++) {
      pn[0][j] = pn[1][j] - dx * gradp;
      pn[xtot - 1][j] = pn[xtot - 2][j] - dx * gradp;
    }

    /* Compute L2-norm */
    err_tot = err_u = err_v = err_p = err_d = 0.0;
    #pragma omp parallel for private(i,j) schedule(auto) \
                             reduction(+:err_u, err_v, err_p, err_d)
    for (i = xlo; i < xhi; i++) {
      for (j = ylo; j < yhi; j++) {
        err_u += pow(un[i][j] - u[i][j], 2);
        err_v += pow(vn[i][j] - v[i][j], 2);
        err_p += pow(pn[i][j] - p[i][j], 2);
        err_d += fabs((un[i][j] - un[i-1][j]) * dy
                      + (vn[i][j] - vn[i][j-1]) * dx);
      }
    }

    err_u = sqrt(dt * dx * dy * err_u);
    err_v = sqrt(dt * dx * dy * err_v);
    err_p = sqrt(dt * dx * dy * err_p);
    err_tot = max(err_u, err_v);
    err_tot = max(err_tot, err_p);
    err_tot = max(err_tot, dt * dx * dy * err_d);

    /* Check if solution diverged */
    if (isnan(err_tot)) {
        printf("Solution Diverged after %d steps!\n", step);
        exit(EXIT_FAILURE);
    }

    /* Compute relative error */
    fprintf(flog ,"%d \t %.8e \t %.8e \t %.8e \t %.8e \t %.8e\n",
            step, err_tot, err_u, err_v, err_p, err_d);

    /* Changing pointers to point to the newly computed fields */
    utmp = u;
    u = un;
    un = utmp;

    vtmp = v;
    v = vn;
    vn = vtmp;

    ptmp = p;
    p = pn;
    pn = ptmp;

    step += 1;
  }

  printf("Converged after %d steps\n", step);
  fclose(flog);

  /* Computing new arrays for computing the fields along lines crossing
     the center of the domain in x- and y-directions */
  vc = array_2d(IX, IY);
  uc = array_2d(IX, IY);
  pc = array_2d(IX, IY);

  #pragma omp parallel for private(i,j) schedule(auto)
  for (i = 0; i < IX; i++) {
   for (j = 0; j < IY; j++) {
     uc[i][j] = 0.5 * (u[i][j] + u[i][j + 1]);
     vc[i][j] = 0.5 * (v[i][j] + v[i + 1][j]);
     pc[i][j] = 0.25 * (p[i][j] + p[i + 1][j] + p[i][j + 1] + p[i + 1][j + 1]);
   }
  }

  /* Free the memory */
  free(*ubufo);
  free(ubufo);
  free(*ubufn);
  free(ubufn);

  free(*vbufo);
  free(vbufo);
  free(*vbufn);
  free(vbufn);

  free(*pbufo);
  free(pbufo);
  free(*pbufn);
  free(pbufn);

  /* Writing fields data for post-processing */
  FILE *fuc, *fvc, *fd;

  /* Velocity field value along a line crossing the middle of x-axis */
  fuc = fopen("data/Central_U","w+t");
  fprintf(fuc, "# U, Y\n");

  for (j = 0; j < IY; j++)
    fprintf(fuc, "%5.8lf \t %5.8lf\n", 0.5 * (uc[xm][j] + uc[xm + 1][j]),
            (double) j*dy );

  fclose(fuc);

  /* Velocity field value along a line crossing the middle of y-axis */
  fvc = fopen("data/Central_V","w+t");
  fprintf(fuc, "# V, X\n");

  for (i = 0; i < IX; i++)
    fprintf(fvc, "%5.8lf \t %5.8lf\n", 0.5 * (vc[i][ym] + vc[i][ym + 1]),
            (double) i*dx );

  fclose(fvc);

  /* Writing all the field data */
  fd = fopen("data/xyuvp","w+t");
  fprintf(fd, "# X \t Y \t U \t V \t P\n");
  for (i = 0; i < IX; i++)
      for (j = 0; j < IY; j++)
         fprintf(fd, "%5.8lf \t %5.8lf \t %5.8lf \t %5.8lf \t %5.8lf\n",
                 (double) i*dx, (double) j*dy, uc[i][j], vc[i][j], pc[i][j]);
  fclose(fd);

  /* Free the memory */
  free(*uc);
  free(uc);

  free(*vc);
  free(vc);

  free(*pc);
  free(pc);

  return 0;
}
