#include "utilities.h"

/* Generate a 2D array using pointer to pointer */
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

/* Free the memory */
void freeMem(int count, ...) {
  va_list args;
  va_start(args, count);

  double **arr;
  for (int i = 0; i < count; i++) {
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
  va_start(args, count);

  double max = va_arg(args, double);
  for (int i = 1; i < count; i++)
    max = fmax(va_arg(args, double), max);

  va_end(args);
  return max;
}
