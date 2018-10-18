#ifndef FUNCTIONS_H
#define FUNCTIONS_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

/* Specify number of grid points in x and y directions */
#define IX 128
#define IY 128

/* Generating a 2D array using pointer to a pointer */
double **array_2d(const int row, const int col);

/* Average in x and y-direction for computing field data on grid points */
void phi_gc(double **phi, double **phi_x, double **phi_y, const int r[2][4]);

#endif /* FUNCTIONS_H */
