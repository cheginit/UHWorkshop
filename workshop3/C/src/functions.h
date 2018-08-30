/* -------------------------------Functions---------------------------------- */
/* Generating a 2D array using pointer to a pointer */
double **array_2d(const int row, const int col) {
  double **arr = (double**) 0;

  arr = (double **) malloc(sizeof(double *) * row);
  arr[0] = (double *) malloc(sizeof(double) * col * row);

  for (int i = 0; i < row; i++)
    arr[i] = (*arr + col * i);

  if(!arr) {
    printf("Memory allocation error.\n");
    if(*arr) free(*arr);
    if(arr) free(arr);
    exit(EXIT_FAILURE);
  }
  for (int i = 0; i < row; i++)
    for (int j = 0; j < col; j++)
      arr[i][j] = 0.0;

  return arr;
}

/* Average in x-direction at grid point i, j */
double avx(double **phi, const int i, const int j) {
    return (phi[i][j] + phi[i - 1][j]) * 0.5;
}

/* Average in y-direction at grid point i, j */
double avy(double **phi, const int i, const int j) {
    return (phi[i][j] + phi[i][j - 1]) * 0.5;
}

/* Calculate laplacian of a field at grid point i, j */
double div2(double **phi, const int i, const int j,
            const double dtdxx, const double dtdyy) {
  return dtdxx * (phi[i + 1][j] - 2.0 * phi[i][j] + phi[i - 1][j])
       + dtdyy * (phi[i][j + 1] - 2.0 * phi[i][j] + phi[i][j - 1]);
}
