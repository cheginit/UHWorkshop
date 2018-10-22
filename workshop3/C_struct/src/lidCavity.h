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

  /* Boundary conditions */
  double *ubc;
  double *vbc;
  double *pbc;

  /* Grid Spacing */
  double dx;
  double dy;
} g;

struct FieldPointers {
  /* Two pointers to the generated buffers for each variable */
  double **u;
  double **un;
  double **v;
  double **vn;
  double **p;
  double **pn;
} f;

struct SimulationInfo {
  /* Flow parameters based on inputs */
  double dt;
  double nu;
  double c2;
  double cfl;
} s;

/* Generate a 2D array using pointer to a pointer */
double **array_2D(int row, int col);

/* Update the fields to the new time step for the next iteration */
void update(struct FieldPointers *f);

/* Applying boundary conditions for velocity */
void set_UBC(struct FieldPointers *f, struct Grid2D *g);

/* Applying boundary conditions for pressure */
void set_PBC(struct FieldPointers *f, struct Grid2D *g, struct SimulationInfo *s);

/* Free the memory */
void freeMem(double **phi, ...);

/* Find mamximum of a set of float numebrs */
double fmaxof(double errs, ...);

/* Save fields data to files */
void dump_data(struct Grid2D *g, struct SimulationInfo *s);

#endif /* FUNCTIONS_H */
