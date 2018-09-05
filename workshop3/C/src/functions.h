#ifndef FUNCTIONS_H
#define FUNCTIONS_H

#include <stdio.h>
#include <stdlib.h>

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

  for (i = r[0][0]; i < r[0][1]; i++) {
    for (j = r[0][2]; j < r[0][3]; j++) {
      phi_x[i][j] = (phi[i][j] + phi[i - 1][j]) * 0.5;
    }
  }

  for (i = r[1][0]; i < r[1][1]; i++) {
    for (j = r[1][2]; j < r[1][3]; j++) {
      phi_y[i][j] = (phi[i][j] + phi[i][j - 1]) * 0.5;
    }
  }
}

/* Calculate laplacian of a field at grid point i, j */
double div2(double **phi, const int r, const int c,
            const double dtdxx, const double dtdyy) {
  return dtdxx * (phi[r + 1][c] - 2.0 * phi[r][c] + phi[r - 1][c])
       + dtdyy * (phi[r][c + 1] - 2.0 * phi[r][c] + phi[r][c - 1]);
}

#endif /* FUNCTIONS_H */
