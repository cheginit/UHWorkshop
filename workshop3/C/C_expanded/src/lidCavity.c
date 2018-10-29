/*============================================================================*\
Solves Navier-Stokes equations for incompressible, laminar, steady flow using
artificial compressibility method on staggered grid.

The governing equations are as follows:
P_t + c^2 div[u] = 0
u_t + u . grad[u] = - grad[P] + nu div[grad[u]]

where P is p/rho and c represents artificial sound's speed.

Lid-Driven Cavity case:
Dimensions : 1x1 m
Grid size  : 128 x 128
Re number  : 100 / 1000 / 5000 / 10000
Grid type  : Staggered Arakawa C
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
#include "lidCavity.h"

/* Main program */
int main(int argc, char *argv[]) {
  /* Two arrays are required for each Variable; one for old time step and one
   * for the new time step. */
  double **ubufo, **ubufn, **u, **un, **utmp;
  double **vbufo, **vbufn, **v, **vn, **vtmp;
  double **pbufo, **pbufn, **p, **pn, **ptmp;

  double dx, dy, dt, Re, nu;
  double dtdx, dtdy, dtdxx, dtdyy, dtdxdy;

  /* Boundary conditions: {top, left, bottom, right} */
  double ubc[4] = {1.0, 0.0, 0.0, 0.0};
  double vbc[4] = {0.0, 0.0, 0.0, 0.0};
  double pbc[4] = {0.0, 0.0, 0.0, 0.0};

  int i, j, count;

  double err_tot, err_u, err_v, err_p, err_d;

  /* Define first and last interior nodes number for improving
   * code readability*/
  const int xlo = 1, xhi = IX - 1, xtot = IX + 1;
  const int ylo = 1, yhi = IY - 1, ytot = IY + 1;

  /* Simulation parameters */
  double cfl, c2;
  const double tol = 1.0e-7, l_lid = 1.0;
  int itr = 1, itr_max = 1000000;

  /* Getting Reynolds number */
  Re = 100.0;
  if (argc > 1) {
    char *ptr;
    Re = strtod(argv[1], &ptr);
  }
  printf("Re number is set to %d\n", (int)Re);

  if (Re < 500.0) {
    cfl = 0.15;
    c2 = 5.0;
  } else if (Re < 2000.0) {
    cfl = 0.20;
    c2 = 5.8;
  } else {
    cfl = 0.05;
    c2 = 5.8;
  }

  /* Create a log file for outputting the residuals */
  FILE *flog;
  flog = fopen("data/residual", "w+t+e");

  /* Compute flow parameters based on inputs */
  dx = l_lid / (double)(IX - 1);
  dy = dx;
  dt = cfl * fmin(dx, dy) / ubc[0];
  nu = ubc[0] * l_lid / Re;

  /* Carry out operations that their values do not change in loops */
  dtdx = dt / dx;
  dtdy = dt / dy;
  dtdxx = dt / (dx * dx);
  dtdyy = dt / (dy * dy);
  dtdxdy = dt * dx * dy;

  /* Generate two 2D arrays for storing old and new velocity field
   * in x-direction */
  ubufo = array_2d(IX, ytot);
  ubufn = array_2d(IX, ytot);
  /* Define two pointers to the generated buffers for velocity field
   * in x-direction */
  u = ubufo;
  un = ubufn;

  /* Generate two 2D arrays for storing old and new velocity field
   * in y-direction */
  vbufo = array_2d(xtot, IY);
  vbufn = array_2d(xtot, IY);
  /* Define two pointers to the generated buffers for velocity field
   * in y-direction */
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
    u[i][ytot - 1] = ubc[0];
    u[i][ytot - 2] = ubc[0];
  }

  /* Applying boundary conditions */
  set_UBC(u, v, ubc, vbc);
  set_PBC(p, pbc, dx, dy);

  /* Start the main loop */
  do {
/* Solve x-momentum equation for computing u */
#pragma omp parallel for private(i, j) schedule(auto)
    for (i = xlo; i < xhi; i++) {
      for (j = ylo; j < ytot - 1; j++) {
        un[i][j] =
            u[i][j] -
            0.25 * dtdx *
                (pow(u[i + 1][j] + u[i][j], 2) -
                 pow(u[i][j] + u[i - 1][j], 2)) -
            0.25 * dtdy *
                ((u[i][j + 1] + u[i][j]) * (v[i + 1][j] + v[i][j]) -
                 (u[i][j] + u[i][j - 1]) * (v[i + 1][j - 1] + v[i][j - 1])) -
            dtdx * (p[i + 1][j] - p[i][j]) +
            nu * (dtdxx * (u[i + 1][j] - 2.0 * u[i][j] + u[i - 1][j]) +
                  dtdyy * (u[i][j + 1] - 2.0 * u[i][j] + u[i][j - 1]));
      }
    }

/* Solve y-momentum for computing v */
#pragma omp parallel for private(i, j) schedule(auto)
    for (i = xlo; i < xtot - 1; i++) {
      for (j = ylo; j < yhi; j++) {
        vn[i][j] =
            v[i][j] -
            0.25 * dtdx *
                ((u[i][j + 1] + u[i][j]) * (v[i + 1][j] + v[i][j]) -
                 (u[i - 1][j + 1] + u[i - 1][j]) * (v[i][j] + v[i - 1][j])) -
            0.25 * dtdy *
                (pow(v[i][j + 1] + v[i][j], 2) -
                 pow(v[i][j] + v[i][j - 1], 2)) -
            dtdy * (p[i][j + 1] - p[i][j]) +
            nu * (dtdxx * (v[i + 1][j] - 2.0 * v[i][j] + v[i - 1][j]) +
                  dtdyy * (v[i][j + 1] - 2.0 * v[i][j] + v[i][j - 1]));
      }
    }

    set_UBC(un, vn, ubc, vbc);

/* Solves continuity equation for computing P */
#pragma omp parallel for private(i, j) schedule(auto)
    for (i = xlo; i < xtot - 1; i++) {
      for (j = ylo; j < ytot - 1; j++) {
        pn[i][j] = p[i][j] - c2 * ((un[i][j] - un[i - 1][j]) * dtdx +
                                   (vn[i][j] - vn[i][j - 1]) * dtdy);
      }
    }

    set_PBC(pn, pbc, dx, dy);

    /* Compute L2-norm */
    err_u = err_v = err_p = err_d = 0.0;
#pragma omp parallel for private(i,j) schedule(auto) \
                             reduction(+:err_u, err_v, err_p, err_d)
    for (i = xlo; i < xhi; i++) {
      for (j = ylo; j < yhi; j++) {
        err_u += pow(un[i][j] - u[i][j], 2);
        err_v += pow(vn[i][j] - v[i][j], 2);
        err_p += pow(pn[i][j] - p[i][j], 2);
        err_d +=
            (un[i][j] - un[i - 1][j]) * dtdx + (vn[i][j] - vn[i][j - 1]) * dtdy;
      }
    }

    err_u = sqrt(dtdxdy * err_u);
    err_v = sqrt(dtdxdy * err_v);
    err_p = sqrt(dtdxdy * err_p);
    err_d = fabs(err_d);

    count = 4;
    err_tot = fmaxof(count, err_u, err_v, err_p, err_d);

    /* Check if solution diverged */
    if (isnan(err_tot)) {
      printf("Solution Diverged after %d iterations!\n", itr);

      /* Free the memory and terminate */
      count = 6;
      freeMem(count, ubufo, vbufo, pbufo, ubufn, vbufn, pbufn);
      exit(EXIT_FAILURE);
    }

    /* Write relative error */
    fprintf(flog, "%d \t %.8lf \t %.8lf \t %.8lf \t %.8lf \t %.8lf\n", itr,
            err_tot, err_u, err_v, err_p, err_d);

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

    itr += 1;
  } while (err_tot > tol && itr < itr_max);

  if (itr == itr_max) {
    printf("Maximum number of iterations, %d, exceeded\n", itr);

    /* Free the memory and terminate */
    count = 6;
    freeMem(count, ubufo, vbufo, pbufo, ubufn, vbufn, pbufn);

    exit(EXIT_FAILURE);
  }

  printf("Converged after %d iterations\n", itr);
  fclose(flog);

  /* Write output data */
  dump_data(u, v, p, dx, dy);
  /* Free the memory */
  count = 6;
  freeMem(count, ubufo, vbufo, pbufo, ubufn, vbufn, pbufn);
  return 0;
}

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *  *\
* ===============================Functions==================================== *
\* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *  */

/* Generate a 2D array using pointer to a pointer */
double **array_2d(const int row, const int col) {
  double **arr = (double **)0;

  arr = (double **)malloc(sizeof(double *) * row);
  arr[0] = (double *)malloc(sizeof(double) * col * row);

  for (int i = 0; i < row; i++) {
    arr[i] = (*arr + col * i);
  }

  if (!arr) {
    printf("Memory allocation error.\n");
    if (*arr) {
      free(*arr);
    }
    if (arr) {
      free(arr);
    }

    exit(EXIT_FAILURE);
  }

  for (int i = 0; i < row; i++) {
    for (int j = 0; j < col; j++) {
      arr[i][j] = 0.0;
    }
  }

  return arr;
}

/* Applying boundary conditions for velocity */
void set_UBC(double **un, double **vn, double ubc[4], double vbc[4]) {
  int i, j;

  /* Dirichlet boundary condition */
  for (i = 1; i < IX - 1; i++) {
    un[i][0] = 2.0 * ubc[2] - un[i][1];
    un[i][IY] = 2.0 * ubc[0] - un[i][IY - 1];
  }
  for (j = 0; j < IY + 1; j++) {
    un[0][j] = ubc[1];
    un[IX - 1][j] = ubc[3];
  }

  /* Dirichlet boundary condition */
  for (i = 1; i < IX; i++) {
    vn[i][0] = vbc[2];
    vn[i][IY - 1] = vbc[0];
  }
  for (j = 0; j < IY; j++) {
    vn[0][j] = 2.0 * vbc[1] - vn[1][j];
    vn[IX][j] = 2.0 * vbc[3] - vn[IX - 1][j];
  }
}

/* Applying boundary conditions for pressure */
void set_PBC(double **pn, double pbc[4], double dx, double dy) {
  /* Neumann boundary condition */
  for (int i = 1; i < IX; i++) {
    pn[i][0] = pn[i][1] - dy * pbc[2];
    pn[i][IY] = pn[i][IY - 1] - dy * pbc[0];
  }

  for (int j = 0; j < IY + 1; j++) {
    pn[0][j] = pn[1][j] - dx * pbc[1];
    pn[IX][j] = pn[IX - 1][j] - dx * pbc[3];
  }
}

/* Free the memory */
void freeMem(int count, ...) {
  double **arr;
  va_list args;

  va_start(args, count);
  for (int i = 1; i <= count; i++) {
    arr = va_arg(args, double **);
    free(*arr);
    *arr = NULL;
    free(arr);
    arr = NULL;
  }
  va_end(args);
}

/* Find mamximum of a set of float numebrs */
double fmaxof(int count, ...) {
  va_list args;
  double max;

  va_start(args, count);
  max = va_arg(args, double);

  for (int i = 2; i <= count; i++) {
    max = fmax(va_arg(args, double), max);
  }

  va_end(args);

  return max;
}

/* Save fields data to files */
void dump_data(double **u, double **v, double **p, double dx, double dy) {
  int i, j;

  /* Writing fields data for post-processing */
  FILE *fd;

  /* Computing new arrays for computing the fields along lines crossing
   the centers of the axes in x- and y-directions */
  double **u_g, **v_g, **p_g;
  v_g = array_2d(IX, IY);
  u_g = array_2d(IX, IY);
  p_g = array_2d(IX, IY);

#pragma omp parallel for private(i, j) schedule(auto)
  for (i = 0; i < IX; i++) {
    for (j = 0; j < IY; j++) {
      u_g[i][j] = 0.5 * (u[i][j + 1] + u[i][j]);
      v_g[i][j] = 0.5 * (v[i + 1][j] + v[i][j]);
      p_g[i][j] =
          0.25 * (p[i][j] + p[i + 1][j] + p[i][j + 1] + p[i + 1][j + 1]);
    }
  }

  /* Writing all the field data */
  fd = fopen("data/xyuvp", "w+t");
  fprintf(fd, "# X \t Y \t U \t V \t P\n");
  for (i = 0; i < IX; i++) {
    for (j = 0; j < IY; j++) {
      fprintf(fd, "%.8lf \t %.8lf \t %.8lf \t %.8lf \t %.8lf\n", (double)i * dx,
              (double)j * dy, u[i][j], v[i][j], p[i][j]);
    }
  }

  fclose(fd);
  /* Free the memory */
  int count = 3;
  freeMem(count, u_g, v_g, p_g);
}
