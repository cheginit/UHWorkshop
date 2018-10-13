#ifndef FUNCTIONS_H
#define FUNCTIONS_H

#include <math.h>
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>

/* Specify number of grid points in x and y directions */
#define IX 128
#define IY 128

/* Generate a 2D array using pointer to a pointer */
double **array_2D(const int row, const int col) {
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

struct Grid2D {
  /* Two arrays are required for each Variable; one for old time step and one
   * for the new time step. */
  double **ubufo;
  double **ubufn;
  double **vbufo;
  double **vbufn;
  double **pbufo;
  double **pbufn;

  /* Computing new arrays for computing the fields along lines crossing
   * the centers of the axes in x- and y-directions */
  double **v_g;
  double **u_g;
  double **p_g;

  /* Boundary conditions: {top, left, bottom, right} */
  double *ubc;
  double *vbc;
  double *pbc;

  double dx;
  double dy;
} g;

struct Simulation {
  /* Define two pointers to the generated buffers for each variable*/
  double **u;
  double **un;
  double **v;
  double **vn;
  double **p;
  double **pn;

  /* Flow parameters based on inputs */
  double dt;
  double nu;
  double c2;
} s;

void update(struct Simulation *s) {
  double **tmp;

  tmp = s->u;
  s->u = s->un;
  s->un = tmp;

  tmp = s->v;
  s->v = s->vn;
  s->vn = tmp;

  tmp = s->p;
  s->p = s->pn;
  s->pn = tmp;
}
void set_UBC(struct Simulation *s, struct Grid2D *g) {
  int i, j;

  /* Dirichlet boundary condition */
  for (i = 1; i < IX - 1; i++) {
    s->un[i][0] = g->ubc[2] - s->un[i][1];
    s->un[i][IY] = 2.0 * g->ubc[0] - s->un[i][IY - 1];
  }
  for (j = 0; j < IY + 1; j++) {
    s->un[0][j] = g->ubc[1];
    s->un[IX - 1][j] = g->ubc[3];
  }

  /* Dirichlet boundary condition */
  for (i = 1; i < IX; i++) {
    s->vn[i][0] = g->vbc[2];
    s->vn[i][IY - 1] = g->vbc[0];
  }
  for (j = 0; j < IY; j++) {
    s->vn[0][j] = g->vbc[1] - s->vn[1][j];
    s->vn[IX][j] = g->vbc[3] - s->vn[IX - 1][j];
  }
}

/* Applying boundary conditions for pressure */
void set_PBC(struct Simulation *s, struct Grid2D *g) {
  /* Neumann boundary condition */
  for (int i = 1; i < IX; i++) {
    s->pn[i][0] = s->pn[i][1] - g->dy * g->pbc[2];
    s->pn[i][IY] = s->pn[i][IY - 1] - g->dy * g->pbc[0];
  }

  for (int j = 0; j < IY + 1; j++) {
    s->pn[0][j] = s->pn[1][j] - g->dx * g->pbc[1];
    s->pn[IX][j] = s->pn[IX - 1][j] - g->dx * g->pbc[3];
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

void dump_data(struct Grid2D *g) {
  const int xm = IX / 2, ym = IY / 2;
  int i, j;

  /* Writing fields data for post-processing */
  FILE *fug, *fvg, *fd;

  /* Velocity field value along a line crossing the middle of x-axis */
  fug = fopen("data/Central_U", "w+t+e");
  fprintf(fug, "# U, Y\n");

  for (j = 0; j < IY; j++) {
    fprintf(fug, "%.8lf \t %.8lf\n", 0.5 * (g->u_g[xm][j] + g->u_g[xm + 1][j]),
            (double)j * g->dy);
  }

  fclose(fug);

  /* Velocity field value along a line crossing the middle of y-axis */
  fvg = fopen("data/Central_V", "w+t+e");
  fprintf(fug, "# V, X\n");

  for (i = 0; i < IX; i++) {
    fprintf(fvg, "%.8lf \t %.8lf\n", 0.5 * (g->v_g[i][ym] + g->v_g[i][ym + 1]),
            (double)i * g->dx);
  }

  fclose(fvg);

  /* Writing all the field data */
  fd = fopen("data/xyuvp", "w+t+e");
  fprintf(fd, "# X \t Y \t U \t V \t P\n");
  for (i = 0; i < IX; i++) {
    for (j = 0; j < IY; j++) {
      fprintf(fd, "%.8lf \t %.8lf \t %.8lf \t %.8lf \t %.8lf\n",
              (double)i * g->dx, (double)j * g->dy, g->u_g[i][j], g->v_g[i][j],
              g->p_g[i][j]);
    }
  }

  fclose(fd);
}

#endif /* FUNCTIONS_H */
