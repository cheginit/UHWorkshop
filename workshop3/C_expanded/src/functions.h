#ifndef FUNCTIONS_H
#define FUNCTIONS_H

#include <stdio.h>
#include <stdarg.h>
#include <stdlib.h>
#include <math.h>

/* Specify number of grid points in x and y directions */
#define IX 128
#define IY 128

/* Generate a 2D array using pointer to a pointer */
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

void set_UBC(double **un, double **vn, double ubc[4], double vbc[4]) {
    int i, j;

    /* Dirichlet boundary condition */
   for (i = 1; i < IX - 1; i++) {
      un[i][0] = ubc[2] - un[i][1];
      un[i][IY] = 2.0*ubc[0] - un[i][IY - 1];
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
      vn[0][j] = vbc[1] - vn[1][j];
      vn[IX][j] = vbc[3] - vn[IX - 1][j];
    }
}

/* Applying boundary conditions for pressure */
void set_PBC(double **pn, double pbc[4], double dx, double dy) {
    /* Neumann boundary condition */
    for (int i = 1; i < IX; i++) {
      pn[i][0] = pn[i][1] - dy * pbc[2];
      pn[i][IY] = pn[i][IY - 1]  - dy * pbc[0];
    }

    for (int j = 0; j < IY + 1; j++) {
      pn[0][j] = pn[1][j] - dx * pbc[1];
      pn[IX][j] = pn[IX - 1][j] - dx * pbc[3];
    }
}

/* Free the memory */
void freeMem(double **phi, ...) {
    va_list args;
    va_start(args, phi);
    free(*va_arg(args, double**));
    free(va_arg(args, double**));
    va_end(args);
}

/* Find mamximum of a set of float numebrs */
double fmaxof(double errs, ...){
    double max = errs, val;
    va_list args;

    va_start(args, errs);
    val = va_arg(args, double);
    max = fmax(val, max);
    va_end(args);
    return max;
}

void dump_data(double **u, double **v, double **p, double dx, double dy) {
    const int xm = IX/2, ym = IY/2;
    int i, j;

    /* Writing fields data for post-processing */
    FILE *fug, *fvg, *fd;

    /* Velocity field value along a line crossing the middle of x-axis */
    fug = fopen("data/Central_U","w+t");
    fprintf(fug, "# U, Y\n");

    for (j = 0; j < IY; j++) {
    fprintf(fug, "%.8lf \t %.8lf\n", 0.5 * (u[xm][j] + u[xm + 1][j]),
            (double) j*dy );
    }

    fclose(fug);

    /* Velocity field value along a line crossing the middle of y-axis */
    fvg = fopen("data/Central_V","w+t");
    fprintf(fug, "# V, X\n");

    for (i = 0; i < IX; i++) {
    fprintf(fvg, "%.8lf \t %.8lf\n", 0.5 * (v[i][ym] + v[i][ym + 1]),
            (double) i*dx );
    }

    fclose(fvg);

    /* Writing all the field data */
    fd = fopen("data/xyuvp","w+t");
    fprintf(fd, "# X \t Y \t U \t V \t P\n");
    for (i = 0; i < IX; i++) {
        for (j = 0; j < IY; j++) {
           fprintf(fd, "%.8lf \t %.8lf \t %.8lf \t %.8lf \t %.8lf\n",
                   (double) i*dx, (double) j*dy, u[i][j], v[i][j], p[i][j]);
        }
    }

    fclose(fd);
}

#endif /* FUNCTIONS_H */
