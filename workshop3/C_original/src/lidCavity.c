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
int main (int argc, char *argv[])
{
  /* 4 arrays are required for velocity fields and two for pressure.
   * Velocity: old, new, grid points, cell centers
   * Pressure: old and new */
  double **ubufo, **ubufn, **ubufg, **ubufc, **u, **un, **ug, **uc, **utmp;
  double **vbufo, **vbufn, **vbufg, **vbufc, **v, **vn, **vg, **vc,**vtmp;
  double **pbufo, **pbufn, **p, **pn, **ptmp;

  double dx, dy, dt, Re, nu;
  double dtdx, dtdy, dtdxx, dtdyy, dtdxdy;

  int i, j, itr = 1, itr_max = 1e6;
  double ut, ub, ul, ur;
  double vt, vb, vl, vr;
  double err_tot, err_u, err_v, err_p, err_d;

  /* Define first and last interior nodes number for better
   * code readability*/
  const int xlo = 1, xhi = IX - 1, xtot = IX + 1, xm = IX/2;
  const int ylo = 1, yhi = IY - 1, ytot = IY + 1, ym = IY/2;

  /* Define arrays for min and max bounds for slicing arrays.
   * First row are bounds for averaging in x direction and second
   * row for averaging in y direction. */
  const int ru[2][4] = {{xlo, IX, 0, ytot - 1},
                        {0, IX, ylo, ytot - 1}};
  const int rv[2][4] = {{xlo, xtot - 1, 0, IY},
                        {0, xtot - 1, ylo, IY}};

  const double tol = 1.0e-7, l_lid = 1.0, gradp = 0.0;

  /* Simulation parameters */
  double cfl , c2;

  /* Getting Reynolds number */
  if(argc <= 1) {
    Re = 100.0;
  } else {
    char *ptr;
    Re = strtod(argv[1], &ptr);
  }
  printf("Re number is set to %d\n", (int) Re);

  if (Re < 500) {
      cfl = 0.15;
      c2 = 5.0;
  } else if (Re < 2000) {
      cfl = 0.20;
      c2 = 5.8;
  } else {
      cfl = 0.05;
      c2 = 5.8;
  }

  /* Boundary condition values */
  ut = 1.0;
  ub = ul = ur = 0.0;
  vt = vb = vl = vr = 0.0;

  /* Create a log file for outputting the residuals */
  FILE *flog;
  flog = fopen("data/residual","w+t");

  /* Compute flow parameters based on inputs */
  dx = l_lid / (double) (IX - 1);
  dy = dx; //l_lid / (double) (IY - 1);
  dt = cfl * fmin(dx,dy) / ut;
  nu = ut * l_lid / Re;

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

  /* Generate two 2D arrays for storing grid and center velocity field
   * in x-direction */
  ubufg = array_2d(IX, ytot);
  ubufc = array_2d(IX, ytot);
  /* Define two pointers to the generated buffers for velocity field
   * in x-direction */
  ug = ubufg;
  uc = ubufc;

  /* Generate two 2D arrays for storing old and new velocity field
   * in y-direction */
  vbufo = array_2d(xtot, IY);
  vbufn = array_2d(xtot, IY);
  /* Define two pointers to the generated buffers for velocity field
   * in y-direction */
  v = vbufo;
  vn = vbufn;

  /* Generate two 2D arrays for storing grid and center velocity field
   * in y-direction */
  vbufg = array_2d(xtot, IY);
  vbufc = array_2d(xtot, IY);
  /* Define two pointers to the generated buffers for velocity field
   * in y-direction */
  vg = vbufg;
  vc = vbufc;

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

  /* Applying boundary conditions */
  /* Dirichlet boundary condition */
  for (i = xlo; i < xhi; i++) {
    u[i][0] = ub - u[i][1];
    u[i][ytot - 1] = 2.0*ut - u[i][ytot - 2];
  }
  for (j = 0; j < ytot; j++) {
    u[0][j] = ul;
    u[xhi][j] = ur;
  }

  /* Dirichlet boundary condition */
  for (i = xlo; i < xtot - 1; i++) {
    v[i][0] = vb;
    v[i][yhi] = vt;
  }
  for (j = 0; j < IY; j++) {
    v[0][j] = vl - v[1][j];
    v[xtot - 1][j] = vr - v[xtot - 2][j];
  }

  /* Neumann boundary condition */
  for (i = xlo; i < xtot - 1; i++) {
    p[i][0] = p[i][1] - dy * gradp;
    p[i][ytot - 1] = p[i][ytot - 2]  - dy * gradp;
  }
  for (j = 0; j < ytot; j++) {
    p[0][j] = p[1][j] - dx * gradp;
    p[xtot - 1][j] = p[xtot - 2][j] - dx * gradp;
  }

  /* Start the main loop */
  do {
    phi_gc(u, uc, ug, ru);
    phi_gc(v, vg, vc, rv);

    /* Solve x-momentum equation for computing u */
    #pragma omp parallel for private(i,j) schedule(auto)
    for (i = xlo; i < xhi; i++) {
      for (j = ylo; j < ytot - 1; j++) {
        un[i][j] = u[i][j]
                   - dtdx * (uc[i+1][j] * uc[i+1][j]
                             - uc[i][j] * uc[i][j])
                   - dtdy * (ug[i][j+1] * vg[i+1][j]
                             - ug[i][j] * vg[i+1][j-1])
                   - dtdx * (p[i+1][j] - p[i][j])
                   + nu * (dtdxx * (u[i+1][j] - 2.0 * u[i][j] + u[i-1][j])
                           + dtdyy * (u[i][j+1] - 2.0 * u[i][j] + u[i][j-1]));
      }
    }

    /* Solve y-momentum for computing v */
    #pragma omp parallel for private(i,j) schedule(auto)
    for (i = xlo; i < xtot - 1; i++) {
      for (j = ylo; j < yhi; j++) {
        vn[i][j] = v[i][j]
                   - dtdx * (ug[i][j+1] * vg[i+1][j]
                             - ug[i-1][j+1] * vg[i][j])
                   - dtdy * (vc[i][j+1] * vc[i][j+1]
                             - vc[i][j] * vc[i][j])
                   - dtdy * (p[i][j+1] - p[i][j])
                   + nu * (dtdxx * (v[i+1][j] - 2.0 * v[i][j] + v[i-1][j])
                           + dtdyy * (v[i][j+1] - 2.0 * v[i][j] + v[i][j-1]));
      }
    }

    /* Dirichlet boundary conditions */
    for (i = xlo; i < xhi; i++) {
      un[i][0] = ub - un[i][1];
      un[i][ytot - 1] = 2.0*ut - un[i][ytot - 2];
    }
    for (j = 0; j < ytot; j++) {
      un[0][j] = ul;
      un[xhi][j] = ur;
    }

    /* Dirichlet boundary conditions */
    for (i = xlo; i < xtot - 1; i++) {
      vn[i][0] = vb;
      vn[i][yhi] = vt;
    }
    for (j = 0; j < IY; j++) {
      vn[0][j] = vr - vn[1][j];
      vn[xtot - 1][j] = vl - vn[xtot - 2][j];
    }

    /* Solves continuity equation for computing P */
    #pragma omp parallel for private(i,j) schedule(auto)
    for (i = xlo; i < xtot - 1; i++) {
      for (j = ylo; j < ytot - 1; j++) {
        pn[i][j] = p[i][j] - c2 * ((un[i][j] - un[i-1][j]) * dtdx
                                   + (vn[i][j] - vn[i][j-1]) * dtdy);
      }
    }

    /* Neumann boundary conditions */
    for (i = xlo; i < xtot - 1; i++) {
      pn[i][0] = pn[i][1] - dy * gradp;
      pn[i][ytot - 1] = pn[i][ytot - 2]  - dy * gradp;
    }
    for (j = 0; j < ytot; j++) {
      pn[0][j] = pn[1][j] - dx * gradp;
      pn[xtot - 1][j] = pn[xtot - 2][j] - dx * gradp;
    }

    /* Compute L2-norm */
    err_u = err_v = err_p = err_d = 0.0;
    #pragma omp parallel for private(i,j) schedule(auto) \
                             reduction(+:err_u, err_v, err_p, err_d)
    for (i = xlo; i < xhi; i++) {
      for (j = ylo; j < yhi; j++) {
        err_u += pow(un[i][j] - u[i][j], 2);
        err_v += pow(vn[i][j] - v[i][j], 2);
        err_p += pow(pn[i][j] - p[i][j], 2);
        err_d += (un[i][j] - un[i-1][j]) * dtdx
                 + (vn[i][j] - vn[i][j-1]) * dtdy;
      }
    }

    err_u = sqrt(dtdxdy * err_u);
    err_v = sqrt(dtdxdy * err_v);
    err_p = sqrt(dtdxdy * err_p);
    err_tot = fmax(err_u, err_v);
    err_tot = fmax(err_tot, err_p);
    err_tot = fmax(err_tot, err_d);

    /* Check if solution diverged */
    if (isnan(err_tot)) {
      printf("Solution Diverged after %d iterations!\n", itr);

      /* Free the memory */
      free(*ubufo);
      free(ubufo);
      free(*ubufn);
      free(ubufn);
      free(*ubufg);
      free(ubufg);
      free(*ubufc);
      free(ubufc);

      free(*vbufo);
      free(vbufo);
      free(*vbufn);
      free(vbufn);
      free(*vbufg);
      free(vbufg);
      free(*vbufc);
      free(vbufc);

      free(*pbufo);
      free(pbufo);
      free(*pbufn);
      free(pbufn);

      exit(EXIT_FAILURE);
    }

    /* Write relative error */
    fprintf(flog ,"%d \t %.8lf \t %.8lf \t %.8lf \t %.8lf \t %.8lf\n",
            itr, err_tot, err_u, err_v, err_p, err_d);

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

    /* Free the memory */
    free(*ubufo);
    free(ubufo);
    free(*ubufn);
    free(ubufn);
    free(*ubufg);
    free(ubufg);
    free(*ubufc);
    free(ubufc);

    free(*vbufo);
    free(vbufo);
    free(*vbufn);
    free(vbufn);
    free(*vbufg);
    free(vbufg);
    free(*vbufc);
    free(vbufc);

    free(*pbufo);
    free(pbufo);
    free(*pbufn);
    free(pbufn);

    exit(EXIT_FAILURE);
  }

  printf("Converged after %d iterations\n", itr);
  fclose(flog);

  /* Computing new arrays for computing the fields along lines crossing
     the center of the domain in x- and y-directions */
  double **u_g, **v_g, **p_g;
  v_g = array_2d(IX, IY);
  u_g = array_2d(IX, IY);
  p_g = array_2d(IX, IY);

  #pragma omp parallel for private(i,j) schedule(auto)
  for (i = 0; i < IX; i++) {
   for (j = 0; j < IY; j++) {
     u_g[i][j] = 0.5 * (u[i][j+1] + u[i][j]);
     v_g[i][j] = 0.5 * (v[i+1][j] + v[i][j]);
     p_g[i][j] = 0.25 * (p[i][j] + p[i + 1][j] + p[i][j + 1] + p[i + 1][j + 1]);
   }
  }

  /* Free the memory */
  free(*ubufo);
  free(ubufo);
  free(*ubufn);
  free(ubufn);
  free(*ubufg);
  free(ubufg);
  free(*ubufc);
  free(ubufc);

  free(*vbufo);
  free(vbufo);
  free(*vbufn);
  free(vbufn);
  free(*vbufg);
  free(vbufg);
  free(*vbufc);
  free(vbufc);

  free(*pbufo);
  free(pbufo);
  free(*pbufn);
  free(pbufn);

  /* Writing fields data for post-processing */
  FILE *fug, *fvg, *fd;

  /* Velocity field value along a line crossing the middle of x-axis */
  fug = fopen("data/Central_U","w+t");
  fprintf(fug, "# U, Y\n");

  for (j = 0; j < IY; j++) {
    fprintf(fug, "%.8lf \t %.8lf\n", 0.5 * (u_g[xm][j] + u_g[xm + 1][j]),
            (double) j*dy );
  }

  fclose(fug);

  /* Velocity field value along a line crossing the middle of y-axis */
  fvg = fopen("data/Central_V","w+t");
  fprintf(fug, "# V, X\n");

  for (i = 0; i < IX; i++) {
    fprintf(fvg, "%.8lf \t %.8lf\n", 0.5 * (v_g[i][ym] + v_g[i][ym + 1]),
            (double) i*dx );
  }

  fclose(fvg);

  /* Writing all the field data */
  fd = fopen("data/xyuvp","w+t");
  fprintf(fd, "# X \t Y \t U \t V \t P\n");
  for (i = 0; i < IX; i++) {
    for (j = 0; j < IY; j++) {
       fprintf(fd, "%.8lf \t %.8lf \t %.8lf \t %.8lf \t %.8lf\n",
               (double) i*dx, (double) j*dy, u_g[i][j], v_g[i][j], p_g[i][j]);
    }
  }

  fclose(fd);

  /* Free the memory */
  free(*u_g);
  free(u_g);

  free(*v_g);
  free(v_g);

  free(*p_g);
  free(p_g);

  return 0;
}

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *  *\
* ===============================Functions==================================== *
\* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *  */

/* Generating a 2D array using pointer to a pointer */
double **array_2d(const int row, const int col) {
  double **arr = (double**) 0;

  arr = (double **) malloc(sizeof(double *) * row);
  arr[0] = (double *) malloc(sizeof(double) * col * row);

  for (int i = 0; i < row; i++) {
    arr[i] = (*arr + col * i);
  }

  if(!arr) {
    printf("Memory allocation error.\n");
    if(*arr) {
      free(*arr);
    }
    if(arr) {
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

/* Average in x and y-direction for computing field data on grid points */
void phi_gc(double **phi, double **phi_x, double **phi_y, const int r[2][4]) {
  int i, j;

  #pragma omp parallel for private(i,j) schedule(auto)
  for (i = r[0][0]; i < r[0][1]; i++) {
    for (j = r[0][2]; j < r[0][3]; j++) {
      phi_x[i][j] = (phi[i][j] + phi[i - 1][j]) * 0.5;
    }
  }

  #pragma omp parallel for private(i,j) schedule(auto)
  for (i = r[1][0]; i < r[1][1]; i++) {
    for (j = r[1][2]; j < r[1][3]; j++) {
      phi_y[i][j] = (phi[i][j] + phi[i][j - 1]) * 0.5;
    }
  }
}

