#ifndef FUNCTIONS_H
#define FUNCTIONS_H

#include <math.h>
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>

/* Specify number of grid points in x and y directions */
#define IX 128
#define IY 128

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

  /* Boundary conditions: {0:top, 1:left, 2:bottom, 3:right} */
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

/* Generate a 2D array using pointer to a pointer */
double **array_2D(const int row, const int col);

/* Update the fields to the new time step for the next iteration */
void update(struct Simulation *s);

/* Applying boundary conditions for velocity */
void set_UBC(struct Simulation *s, struct Grid2D *g);

/* Applying boundary conditions for pressure */
void set_PBC(struct Simulation *s, struct Grid2D *g);

/* Free the memory */
void freeMem(double **phi, ...);

/* Find mamximum of a set of float numebrs */
double fmaxof(double errs, ...);

/* Save fields data to files */
void dump_data(struct Grid2D *g);

#endif /* FUNCTIONS_H */
