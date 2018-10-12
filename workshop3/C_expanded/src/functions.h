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


#endif /* FUNCTIONS_H */
