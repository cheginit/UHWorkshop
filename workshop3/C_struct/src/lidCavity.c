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
  double Re = 100.0;
  double dtdx, dtdy, dtdxx, dtdyy, dtdxdy;
  double err_tot, err_u, err_v, err_p, err_d;
  const double tol = 1.0e-7, l_lid = 1.0;

  int i, j, itr = 1;
  const int itr_max = 1000000;

  FILE *flog;

  /* Getting Reynolds number */
  if (argc > 1) {
    char *ptr;
    Re = strtod(argv[1], &ptr);
  }
  printf("Re number is set to %d\n", (int)Re);

  /* Create a log file for outputting the residuals */
  flog = fopen("data/residual", "w+t");

  g = ((struct Grid2D){.ubufo = array_2D(IX, IY + 1),
                       .ubufn = array_2D(IX, IY + 1),
                       .vbufo = array_2D(IX + 1, IY),
                       .vbufn = array_2D(IX + 1, IY),
                       .pbufo = array_2D(IX + 1, IY + 1),
                       .pbufn = array_2D(IX + 1, IY + 1),
                       .dx = l_lid / (double)(IX - 1),
                       .dy = l_lid / (double)(IY - 1)});

  f = ((struct FieldPointers){.u = g.ubufo,
                              .un = g.ubufn,
                              .v = g.vbufo,
                              .vn = g.vbufn,
                              .p = g.pbufo,
                              .pn = g.pbufn});

  /* Boundary conditions: {0:top, 1:left, 2:bottom, 3:right} */
  s = ((struct SimulationInfo){.ubc = {1.0, 0.0, 0.0, 0.0},
                               .vbc = {0.0, 0.0, 0.0, 0.0},
                               .pbc = {0.0, 0.0, 0.0, 0.0}});

  /* Set c2 and cfl according to Re based on trail and error */
  if (Re < 500.0) {
    s.cfl = 0.15;
    s.c2 = 5.0;
  } else if (Re < 2000.0) {
    s.cfl = 0.20;
    s.c2 = 5.8;
  } else {
    s.cfl = 0.05;
    s.c2 = 5.8;
  }

  s.nu = s.ubc[0] * l_lid / Re;
  s.dt = s.cfl * fmin(g.dx, g.dy) / s.ubc[0];

  /* Carry out operations that their values do not change in loops */
  dtdx = s.dt / g.dx;
  dtdy = s.dt / g.dy;
  dtdxx = s.dt / (g.dx * g.dx);
  dtdyy = s.dt / (g.dy * g.dy);
  dtdxdy = s.dt * g.dx * g.dy;

  /* Apply initial conditions*/
  for (i = 1; i < IX - 1; i++) {
    f.un[i][IY] = s.ubc[0];
    f.un[i][IY - 1] = s.ubc[0];
  }

  /* Applying boundary conditions */
  set_UBC(&f, &g, &s);
  set_PBC(&f, &g, &s);
  update(&f);

  /* Start the main loop */
  do {
/* Solve x-momentum equation for computing u */
#pragma omp parallel for private(i, j) schedule(auto)
    for (i = 1; i < IX - 1; i++) {
      for (j = 1; j < IY; j++) {
        f.un[i][j] =
            f.u[i][j] -
            0.25 * dtdx *
                (pow(f.u[i + 1][j] + f.u[i][j], 2) -
                 pow(f.u[i][j] + f.u[i - 1][j], 2)) -
            0.25 * dtdy *
                ((f.u[i][j + 1] + f.u[i][j]) * (f.v[i + 1][j] + f.v[i][j]) -
                 (f.u[i][j] + f.u[i][j - 1]) *
                     (f.v[i + 1][j - 1] + f.v[i][j - 1])) -
            dtdx * (f.p[i + 1][j] - f.p[i][j]) +
            s.nu * (dtdxx * (f.u[i + 1][j] - 2.0 * f.u[i][j] + f.u[i - 1][j]) +
                    dtdyy * (f.u[i][j + 1] - 2.0 * f.u[i][j] + f.u[i][j - 1]));
      }
    }

/* Solve y-momentum for computing v */
#pragma omp parallel for private(i, j) schedule(auto)
    for (i = 1; i < IX; i++) {
      for (j = 1; j < IY - 1; j++) {
        f.vn[i][j] =
            f.v[i][j] -
            0.25 * dtdx *
                ((f.u[i][j + 1] + f.u[i][j]) * (f.v[i + 1][j] + f.v[i][j]) -
                 (f.u[i - 1][j + 1] + f.u[i - 1][j]) *
                     (f.v[i][j] + f.v[i - 1][j])) -
            0.25 * dtdy *
                (pow(f.v[i][j + 1] + f.v[i][j], 2) -
                 pow(f.v[i][j] + f.v[i][j - 1], 2)) -
            dtdy * (f.p[i][j + 1] - f.p[i][j]) +
            s.nu * (dtdxx * (f.v[i + 1][j] - 2.0 * f.v[i][j] + f.v[i - 1][j]) +
                    dtdyy * (f.v[i][j + 1] - 2.0 * f.v[i][j] + f.v[i][j - 1]));
      }
    }

    set_UBC(&f, &g, &s);

/* Solves continuity equation for computing P */
#pragma omp parallel for private(i, j) schedule(auto)
    for (i = 1; i < IX; i++) {
      for (j = 1; j < IY; j++) {
        f.pn[i][j] = f.p[i][j] - s.c2 * ((f.un[i][j] - f.un[i - 1][j]) * dtdx +
                                         (f.vn[i][j] - f.vn[i][j - 1]) * dtdy);
      }
    }

    set_PBC(&f, &g, &s);

    /* Compute L2-norm */
    err_u = err_v = err_p = err_d = 0.0;
#pragma omp parallel for private(i,j) schedule(auto) \
                             reduction(+:err_u, err_v, err_p, err_d)
    for (i = 1; i < IX - 1; i++) {
      for (j = 1; j < IY - 1; j++) {
        err_u += pow(f.un[i][j] - f.u[i][j], 2);
        err_v += pow(f.vn[i][j] - f.v[i][j], 2);
        err_p += pow(f.pn[i][j] - f.p[i][j], 2);
        err_d += (f.un[i][j] - f.un[i - 1][j]) * dtdx +
                 (f.vn[i][j] - f.vn[i][j - 1]) * dtdy;
      }
    }

    err_u = sqrt(dtdxdy * err_u);
    err_v = sqrt(dtdxdy * err_v);
    err_p = sqrt(dtdxdy * err_p);
    err_tot = fmaxof(err_u, err_v, err_p, err_d);

    /* Check if solution diverged */
    if (isnan(err_tot)) {
      printf("Solution Diverged after %d iterations!\n", itr);

      /* Free the memory and terminate */
      freeMem(g.ubufo, g.vbufo, g.pbufo, g.ubufn, g.vbufn, g.pbufn);
      exit(EXIT_FAILURE);
    }

    /* Write relative error */
    fprintf(flog, "%d \t %.8lf \t %.8lf \t %.8lf \t %.8lf \t %.8lf\n", itr,
            err_tot, err_u, err_v, err_p, err_d);

    /* Update the fields */
    update(&f);
    itr += 1;
  } while (err_tot > tol && itr < itr_max);

  if (itr == itr_max) {
    printf("Maximum number of iterations, %d, exceeded\n", itr);

    /* Free the memory and terminate */
    freeMem(g.ubufo, g.vbufo, g.pbufo, g.ubufn, g.vbufn, g.pbufn);
    exit(EXIT_FAILURE);
  }

  printf("Converged after %d iterations\n", itr);
  fclose(flog);

  /* Write output data */
  dump_data(&g, &s);

  return 0;
}

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *  *\
* ===============================Functions==================================== *
\* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *  */

/* Generate a 2D array using pointer to a pointer */
double **array_2D(int row, int col) {
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

/* Update the fields to the new time step for the next iteration */
void update(struct FieldPointers *f) {
  double **tmp;

  tmp = f->u;
  f->u = f->un;
  f->un = tmp;

  tmp = f->v;
  f->v = f->vn;
  f->vn = tmp;

  tmp = f->p;
  f->p = f->pn;
  f->pn = tmp;
}

/* Applying boundary conditions for velocity */
void set_UBC(struct FieldPointers *f, struct Grid2D *g,
             struct SimulationInfo *s) {
  int i, j;

  /* Dirichlet boundary condition */
  for (i = 1; i < IX - 1; i++) {
    f->un[i][0] = s->ubc[2] - f->un[i][1];
    f->un[i][IY] = 2.0 * s->ubc[0] - f->un[i][IY - 1];
  }
  for (j = 0; j < IY + 1; j++) {
    f->un[0][j] = s->ubc[1];
    f->un[IX - 1][j] = s->ubc[3];
  }

  /* Dirichlet boundary condition */
  for (i = 1; i < IX; i++) {
    f->vn[i][0] = s->vbc[2];
    f->vn[i][IY - 1] = s->vbc[0];
  }
  for (j = 0; j < IY; j++) {
    f->vn[0][j] = s->vbc[1] - f->vn[1][j];
    f->vn[IX][j] = s->vbc[3] - f->vn[IX - 1][j];
  }
}

/* Applying boundary conditions for pressure */
void set_PBC(struct FieldPointers *f, struct Grid2D *g,
             struct SimulationInfo *s) {
  /* Neumann boundary condition */
  for (int i = 1; i < IX; i++) {
    f->pn[i][0] = f->pn[i][1] - g->dy * s->pbc[2];
    f->pn[i][IY] = f->pn[i][IY - 1] - g->dy * s->pbc[0];
  }

  for (int j = 0; j < IY + 1; j++) {
    f->pn[0][j] = f->pn[1][j] - g->dx * s->pbc[1];
    f->pn[IX][j] = f->pn[IX - 1][j] - g->dx * s->pbc[3];
  }
}

/* Free the memory */
void freeMem(double **phi, ...) {
  va_list args;
  va_start(args, phi);
  free(*va_arg(args, double **));
  free(va_arg(args, double **));
  va_end(args);
}

/* Find mamximum of a set of float numebrs */
double fmaxof(double errs, ...) {
  double max = errs, val;
  va_list args;

  va_start(args, errs);
  val = va_arg(args, double);
  max = fmax(val, max);
  va_end(args);
  return max;
}

/* Save fields data to files */
void dump_data(struct Grid2D *g, struct SimulationInfo *s) {
  int i, j;

  double **ug, **vg, **pg;

  FILE *fd;

  ug = array_2D(IX, IY);
  vg = array_2D(IX, IY);
  pg = array_2D(IX, IY);

#pragma omp parallel for private(i, j) schedule(auto)
  for (i = 0; i < IX; i++) {
    for (j = 0; j < IY; j++) {
      ug[i][j] = 0.5 * (f.u[i][j + 1] + f.u[i][j]);
      vg[i][j] = 0.5 * (f.v[i + 1][j] + f.v[i][j]);
      pg[i][j] = 0.25 * (f.p[i][j] + f.p[i + 1][j] + f.p[i][j + 1] +
                         f.p[i + 1][j + 1]);
    }
  }

  /* Free the memory */
  freeMem(g->ubufo, g->vbufo, g->pbufo, g->ubufn, g->vbufn, g->pbufn);

  /* Writing all the field data */
  fd = fopen("data/xyuvp", "w+t+e");
  fprintf(fd, "# X \t Y \t U \t V \t P\n");
  for (i = 0; i < IX; i++) {
    for (j = 0; j < IY; j++) {
      fprintf(fd, "%.8lf \t %.8lf \t %.8lf \t %.8lf \t %.8lf\n",
              (double)i * g->dx, (double)j * g->dy, ug[i][j], vg[i][j],
              pg[i][j]);
    }
  }
  fclose(fd);

  /* Free the memory */
  freeMem(ug, vg, pg);
}
