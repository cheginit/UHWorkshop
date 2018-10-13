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

#endif /* FUNCTIONS_H */
