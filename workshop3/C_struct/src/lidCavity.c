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
#include "functions.h"

/* Main program */
int main(int argc, char *argv[]) {
  double Re, cfl;
  double dtdx, dtdy, dtdxx, dtdyy, dtdxdy;
  double err_tot, err_u, err_v, err_p, err_d;
  const double tol = 1.0e-7, l_lid = 1.0;

  int i, j, itr = 1;
  const int itr_max = 1000000;

  /* Getting Reynolds number */
  if (argc <= 1) {
    Re = 100.0;
  } else {
    char *ptr;
    Re = strtod(argv[1], &ptr);
  }
  printf("Re number is set to %d\n", (int)Re);

  /* Create a log file for outputting the residuals */
  FILE *flog;
  flog = fopen("data/residual", "w+t");

  g = ((struct Grid2D){.ubufo = array_2D(IX, IY + 1),
                       .ubufn = array_2D(IX, IY + 1),
                       .vbufo = array_2D(IX + 1, IY),
                       .vbufn = array_2D(IX + 1, IY),
                       .pbufo = array_2D(IX + 1, IY + 1),
                       .pbufn = array_2D(IX + 1, IY + 1),
                       .ubc = (double *)malloc(4 * sizeof(double)),
                       .vbc = (double *)malloc(4 * sizeof(double)),
                       .pbc = (double *)malloc(4 * sizeof(double)),
                       .dx = l_lid / (double)(IX - 1),
                       .dy = l_lid / (double)(IY - 1)});

  /* Boundary conditions: {0:top, 1:left, 2:bottom, 3:right} */
  g.ubc[0] = 1.0;
  g.ubc[1] = g.ubc[2] = g.ubc[3] = 0.0;
  g.vbc[0] = g.vbc[1] = g.vbc[2] = g.vbc[3] = 0.0;
  g.pbc[0] = g.pbc[1] = g.pbc[2] = g.pbc[3] = 0.0;

  /* Set c2 and cfl according to Re based on trail and error */
  if (Re < 500) {
    cfl = 0.15;
    s.c2 = 5.0;
  } else if (Re < 2000) {
    cfl = 0.20;
    s.c2 = 5.8;
  } else {
    cfl = 0.05;
    s.c2 = 5.8;
  }

  s = ((struct Simulation){.u = g.ubufo,
                           .un = g.ubufn,
                           .v = g.vbufo,
                           .vn = g.vbufn,
                           .p = g.pbufo,
                           .pn = g.pbufn,
                           .dt = cfl * fmin(g.dx, g.dy) / g.ubc[0],
                           .nu = g.ubc[0] * l_lid / Re,
                           .c2 = s.c2});

  /* Carry out operations that their values do not change in loops */
  dtdx = s.dt / g.dx;
  dtdy = s.dt / g.dy;
  dtdxx = s.dt / (g.dx * g.dx);
  dtdyy = s.dt / (g.dy * g.dy);
  dtdxdy = s.dt * g.dx * g.dy;

  /* Apply initial conditions*/
  for (i = 1; i < IX - 1; i++) {
    s.un[i][IY] = g.ubc[0];
    s.un[i][IY - 1] = g.ubc[0];
  }

  /* Applying boundary conditions */
  set_UBC(&s, &g);
  set_PBC(&s, &g);
  update(&s);

  /* Start the main loop */
  do {
/* Solve x-momentum equation for computing u */
#pragma omp parallel for private(i, j) schedule(auto)
    for (i = 1; i < IX - 1; i++) {
      for (j = 1; j < IY; j++) {
        s.un[i][j] =
            s.u[i][j] -
            0.25 * dtdx *
                (pow(s.u[i + 1][j] + s.u[i][j], 2) -
                 pow(s.u[i][j] + s.u[i - 1][j], 2)) -
            0.25 * dtdy *
                ((s.u[i][j + 1] + s.u[i][j]) * (s.v[i + 1][j] + s.v[i][j]) -
                 (s.u[i][j] + s.u[i][j - 1]) *
                     (s.v[i + 1][j - 1] + s.v[i][j - 1])) -
            dtdx * (s.p[i + 1][j] - s.p[i][j]) +
            s.nu * (dtdxx * (s.u[i + 1][j] - 2.0 * s.u[i][j] + s.u[i - 1][j]) +
                    dtdyy * (s.u[i][j + 1] - 2.0 * s.u[i][j] + s.u[i][j - 1]));
      }
    }

/* Solve y-momentum for computing v */
#pragma omp parallel for private(i, j) schedule(auto)
    for (i = 1; i < IX; i++) {
      for (j = 1; j < IY - 1; j++) {
        s.vn[i][j] =
            s.v[i][j] -
            0.25 * dtdx *
                ((s.u[i][j + 1] + s.u[i][j]) * (s.v[i + 1][j] + s.v[i][j]) -
                 (s.u[i - 1][j + 1] + s.u[i - 1][j]) *
                     (s.v[i][j] + s.v[i - 1][j])) -
            0.25 * dtdy *
                (pow(s.v[i][j + 1] + s.v[i][j], 2) -
                 pow(s.v[i][j] + s.v[i][j - 1], 2)) -
            dtdy * (s.p[i][j + 1] - s.p[i][j]) +
            s.nu * (dtdxx * (s.v[i + 1][j] - 2.0 * s.v[i][j] + s.v[i - 1][j]) +
                    dtdyy * (s.v[i][j + 1] - 2.0 * s.v[i][j] + s.v[i][j - 1]));
      }
    }

    set_UBC(&s, &g);

/* Solves continuity equation for computing P */
#pragma omp parallel for private(i, j) schedule(auto)
    for (i = 1; i < IX; i++) {
      for (j = 1; j < IY; j++) {
        s.pn[i][j] = s.p[i][j] - s.c2 * ((s.un[i][j] - s.un[i - 1][j]) * dtdx +
                                         (s.vn[i][j] - s.vn[i][j - 1]) * dtdy);
      }
    }

    set_PBC(&s, &g);

    /* Compute L2-norm */
    err_u = err_v = err_p = err_d = 0.0;
#pragma omp parallel for private(i,j) schedule(auto) \
                             reduction(+:err_u, err_v, err_p, err_d)
    for (i = 1; i < IX - 1; i++) {
      for (j = 1; j < IY - 1; j++) {
        err_u += pow(s.un[i][j] - s.u[i][j], 2);
        err_v += pow(s.vn[i][j] - s.v[i][j], 2);
        err_p += pow(s.pn[i][j] - s.p[i][j], 2);
        err_d += (s.un[i][j] - s.un[i - 1][j]) * dtdx +
                 (s.vn[i][j] - s.vn[i][j - 1]) * dtdy;
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

    /* Changing pointers to point to the newly computed fields */
    update(&s);
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

  g.u_g = array_2D(IX, IY);
  g.v_g = array_2D(IX, IY);
  g.p_g = array_2D(IX, IY);

#pragma omp parallel for private(i, j) schedule(auto)
  for (i = 0; i < IX; i++) {
    for (j = 0; j < IY; j++) {
      g.u_g[i][j] = 0.5 * (s.u[i][j + 1] + s.u[i][j]);
      g.v_g[i][j] = 0.5 * (s.v[i + 1][j] + s.v[i][j]);
      g.p_g[i][j] = 0.25 * (s.p[i][j] + s.p[i + 1][j] + s.p[i][j + 1] +
                            s.p[i + 1][j + 1]);
    }
  }

  /* Free the memory */
  freeMem(g.ubufo, g.vbufo, g.pbufo, g.ubufn, g.vbufn, g.pbufn);

  /* Write output data */
  dump_data(&g);

  /* Free the memory */
  freeMem(g.u_g, g.v_g, g.p_g);

  return 0;
}
