#include "simulationControls.h"

/* Initialize structs */
void initialize(struct FieldPointers *f, struct Grid2D *g,
                struct SimulationInfo *s) {

  g->ubufo = array_2D(g->nx, g->ny + 1);
  g->ubufn = array_2D(g->nx, g->ny + 1);
  g->vbufo = array_2D(g->nx + 1, g->ny);
  g->vbufn = array_2D(g->nx + 1, g->ny);
  g->pbufo = array_2D(g->nx + 1, g->ny + 1);
  g->pbufn = array_2D(g->nx + 1, g->ny + 1);

  f->u = g->ubufo;
  f->un = g->ubufn;
  f->v = g->vbufo;
  f->vn = g->vbufn;
  f->p = g->pbufo;
  f->pn = g->pbufn;

  g->dx = s->l_lid / (double)(g->nx - 1);
  g->dy = s->l_lid / (double)(g->ny - 1);

  /* Set c2 and cfl according to Re based on trail and error */
  if (s->Re < 500.0) {
    s->cfl = 0.15;
    s->c2 = 5.0;
  } else if (s->Re < 2000 - .0) {
    s->cfl = 0.20;
    s->c2 = 5.8;
  } else {
    s->cfl = 0.05;
    s->c2 = 5.8;
  }

  s->nu = s->ubc[0] * s->l_lid / s->Re;
}

/* Compute the time step based on maximum velocity in the domain  */
void set_delt(struct FieldPointers *f, struct Grid2D *g,
              struct SimulationInfo *s) {
  double umax;

  umax = fmax(fmaxarr(f->u, g->nx, g->ny + 1),
              fmaxarr(f->v, g->nx + 1, g->ny));
  s->dt = s->cfl * fmin(g->dx, g->dy) / umax;

  /* Carry out operations that their values do not change in loops */
  s->dtdx = s->dt / g->dx;
  s->dtdy = s->dt / g->dy;
  s->dtdxx = s->dt / (g->dx * g->dx);
  s->dtdyy = s->dt / (g->dy * g->dy);
  s->dtdxdy = s->dt * g->dx * g->dy;
}

/* Apply initial conditions*/
void set_init(struct FieldPointers *f, struct Grid2D *g,
              struct SimulationInfo *s) {
  for (int i = 1; i < g->nx - 1; i++) {
    f->u[i][g->ny] = s->ubc[0];
    f->u[i][g->ny - 1] = s->ubc[0];
  }
}

/* Set boundary conditions for velocity */
void set_UBC(struct FieldPointers *f, struct Grid2D *g,
            struct SimulationInfo *s) {
  int i, j;

  /* Sides */
  for (j = 0; j < g->ny + 1; j++) {
    f->un[0][j] = s->ubc[1];
    f->un[g->nx - 1][j] = s->ubc[3];
  }
  for (j = 0; j < g->ny; j++) {
    f->vn[0][j] = 2.0 * s->vbc[1] - f->vn[1][j];
    f->vn[g->nx][j] = 2.0 * s->vbc[3] - f->vn[g->nx - 1][j];
  }

  /* Bottom and top */
  for (i = 0; i < g->nx; i++) {
    f->un[i][0] = 2.0 * s->ubc[2] - f->un[i][1];
    f->un[i][g->ny] = 2.0 * s->ubc[0] - f->un[i][g->ny - 1];
  }
  for (i = 0; i < g->nx + 1; i++) {
    f->vn[i][0] = s->vbc[2];
    f->vn[i][g->ny - 1] = s->vbc[0];
  }
}

/* Set boundary conditions for pressure */
void set_PBC(struct FieldPointers *f, struct Grid2D *g,
            struct SimulationInfo *s) {
  int i, j;

  /* Sides */
  for (j = 0; j < g->ny + 1; j++) {
    f->pn[0][j] = f->pn[1][j] - g->dx * s->pbc[1];
    f->pn[g->nx][j] = f->pn[g->nx - 1][j] - g->dx * s->pbc[3];
  }

  /* Bottom and top */
  for (i = 0; i < g->nx + 1; i++) {
    f->pn[i][0] = f->pn[i][1] - g->dy * s->pbc[2];
    f->pn[i][g->ny] = f->pn[i][g->ny - 1] - g->dy * s->pbc[0];
  }
}

/* Solve momentum for computing u and v */
void solve_U(struct FieldPointers *f, struct Grid2D *g,
             struct SimulationInfo *s) {
  int i, j;

#pragma omp parallel for private(i, j) schedule(auto)
  for (i = 1; i < g->nx - 1; i++) {
    for (j = 1; j < g->ny; j++) {
      f->un[i][j] =
          f->u[i][j] -
          0.25 * s->dtdx *
              (pow(f->u[i + 1][j] + f->u[i][j], 2) -
               pow(f->u[i][j] + f->u[i - 1][j], 2)) -
          0.25 * s->dtdy *
              ((f->u[i][j + 1] + f->u[i][j]) * (f->v[i + 1][j] + f->v[i][j]) -
               (f->u[i][j] + f->u[i][j - 1]) *
                   (f->v[i + 1][j - 1] + f->v[i][j - 1])) -
          s->dtdx * (f->p[i + 1][j] - f->p[i][j]) +
          s->nu *
              (s->dtdxx * (f->u[i + 1][j] - 2.0 * f->u[i][j] + f->u[i - 1][j]) +
               s->dtdyy * (f->u[i][j + 1] - 2.0 * f->u[i][j] + f->u[i][j - 1]));
    }
  }

#pragma omp parallel for private(i, j) schedule(auto)
  for (i = 1; i < g->nx; i++) {
    for (j = 1; j < g->ny - 1; j++) {
      f->vn[i][j] =
          f->v[i][j] -
          0.25 * s->dtdx *
              ((f->u[i][j + 1] + f->u[i][j]) * (f->v[i + 1][j] + f->v[i][j]) -
               (f->u[i - 1][j + 1] + f->u[i - 1][j]) *
                   (f->v[i][j] + f->v[i - 1][j])) -
          0.25 * s->dtdy *
              (pow(f->v[i][j + 1] + f->v[i][j], 2) -
               pow(f->v[i][j] + f->v[i][j - 1], 2)) -
          s->dtdy * (f->p[i][j + 1] - f->p[i][j]) +
          s->nu *
              (s->dtdxx * (f->v[i + 1][j] - 2.0 * f->v[i][j] + f->v[i - 1][j]) +
               s->dtdyy * (f->v[i][j + 1] - 2.0 * f->v[i][j] + f->v[i][j - 1]));
    }
  }
}

/* Solves continuity equation for computing P */
void solve_P(struct FieldPointers *f, struct Grid2D *g,
             struct SimulationInfo *s) {
  int i, j;

#pragma omp parallel for private(i, j) schedule(auto)
  for (i = 1; i < g->nx; i++) {
    for (j = 1; j < g->ny; j++) {
      f->pn[i][j] =
          f->p[i][j] - s->c2 * ((f->un[i][j] - f->un[i - 1][j]) * s->dtdx +
                                (f->vn[i][j] - f->vn[i][j - 1]) * s->dtdy);
    }
  }
}

/* Compute L2-norm */
void l2_norm(struct FieldPointers *f, struct Grid2D *g,
             struct SimulationInfo *s) {
  int i, j, count;
  double err_u = 0.0, err_v = 0.0, err_p = 0.0, err_d = 0.0;

#pragma omp parallel for private(i,j) schedule(auto) \
                             reduction(+:err_u, err_v, err_p, err_d)
  for (i = 1; i < g->nx - 1; i++) {
    for (j = 1; j < g->ny - 1; j++) {
      err_u += pow(f->un[i][j] - f->u[i][j], 2);
      err_v += pow(f->vn[i][j] - f->v[i][j], 2);
      err_p += pow(f->pn[i][j] - f->p[i][j], 2);
      err_d += (f->un[i][j] - f->un[i - 1][j]) * s->dtdx +
               (f->vn[i][j] - f->vn[i][j - 1]) * s->dtdy;
    }
  }
  s->errs[1] = sqrt(s->dtdxdy * err_u);
  s->errs[2] = sqrt(s->dtdxdy * err_v);
  s->errs[3] = sqrt(s->dtdxdy * err_p);
  s->errs[4] = fabs(err_d);

  count = 4;
  s->errs[0] = fmaxof(count, s->errs[1], s->errs[2], s->errs[3], s->errs[4]);
}
