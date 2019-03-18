#ifndef STRUCTS_H
#define STRUCTS_H

struct Grid2D {
  /* Two arrays are required for each Variable; one for old time step and one
   * for the new time step. */
  double **ubufo;
  double **ubufn;
  double **vbufo;
  double **vbufn;
  double **pbufo;
  double **pbufn;

  /* Number of grid points */
  int nx;
  int ny;

  /* Grid Spacing */
  double dx;
  double dy;
} g;

struct FieldPointers {
  /* Pointers to the generated buffer arrays for each variable */
  double **u;
  double **un;
  double **v;
  double **vn;
  double **p;
  double **pn;
} f;

struct SimulationInfo {
  /* Reynolds number and length of the lid */
  double Re;
  double l_lid;

  /* Boundary conditions: {top, left, bottom, right} */
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

  /* Errors: {total, u err, v err, p err, div U} */
  double errs[5];
} s;

#endif /* STRUCTS_H */
