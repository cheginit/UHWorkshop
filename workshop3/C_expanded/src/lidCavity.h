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
double **array_2d(const int row, const int col);

/* Applying boundary conditions for velocity */
void set_UBC(double **un, double **vn, double ubc[4], double vbc[4]);

/* Applying boundary conditions for pressure */
void set_PBC(double **pn, double pbc[4], double dx, double dy);

/* Free the memory */
void freeMem(int count, ...);

/* Find mamximum of a set of float numebrs */
double fmaxof(int count, ...);

/* Save fields data to files */
void dump_data(double **u, double **v, double **p, double dx, double dy);

#endif /* FUNCTIONS_H */
