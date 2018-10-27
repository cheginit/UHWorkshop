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
  double Re;
  double l_lid;

  /* Boundary conditions */
  double ubc[4];
  double vbc[4];
  double pbc[4];

  /* Flow parameters based on inputs */
  double dt;
  double nu;
  double c2;
  double cfl;

  double dtdx;
  double dtdy;
  double dtdxx;
  double dtdyy;
  double dtdxdy;

  double err_tot;
} s;

/* Applying boundary conditions for velocity */
void initialize(struct FieldPointers *f, struct Grid2D *g,
                struct SimulationInfo *s);

/* Generate a 2D array using pointer to pointer */
double **array_2D(int row, int col);

/* Set initial condition */
void set_init(struct FieldPointers *f, struct SimulationInfo *s);

/* Update the fields to the new time step for the next iteration */
void update(struct FieldPointers *f);

/* Solve momentum for computing u and v */
void solve_U(struct FieldPointers *f, struct SimulationInfo *s);

/* Solves continuity equation for computing P */
void solve_P(struct FieldPointers *f, struct SimulationInfo *s);

/* Applying boundary conditions for velocity */
void set_BC(struct FieldPointers *f, struct Grid2D *g,
            struct SimulationInfo *s);

/* Compute L2-norm */
void l2_norm(struct FieldPointers *f, struct SimulationInfo *s, FILE *flog,
             int itr);

/* Free the memory */
void freeMem(double **phi, ...);

/* Find mamximum of a set of float numebrs */
double fmaxof(double errs, ...);

/* Save fields data to files */
void dump_data(struct Grid2D *g);

#endif /* FUNCTIONS_H */
