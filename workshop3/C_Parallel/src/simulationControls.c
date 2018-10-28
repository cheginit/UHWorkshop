#include "simulationControls.h"

/* Initialize structs */
void initialize(struct FieldPointers *f, struct Grid2D *g,
                struct SimulationInfo *s, int rank, int nprocs) {
  /* set neighbors */
  if (MASTER || LAST_NODE) {
    g->ghosts = 1;
  } else {
    g->ghosts = 2;
  }

  g->nx_p = g->nx_psg = (g->nx / nprocs) + g->ghosts;
  if (LAST_NODE)
    g->nx_psg = g->nx_p + 1;

  g->ubufo = array_2D(g->nx_p, g->ny + 1);
  g->ubufn = array_2D(g->nx_p, g->ny + 1);
  g->vbufo = array_2D(g->nx_psg, g->ny);
  g->vbufn = array_2D(g->nx_psg, g->ny);
  g->pbufo = array_2D(g->nx_psg, g->ny + 1);
  g->pbufn = array_2D(g->nx_psg, g->ny + 1);

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
  s->dt = s->cfl * fmin(g->dx, g->dy) / s->ubc[0];

  /* Carry out operations that their values do not change in loops */
  s->dtdx = s->dt / g->dx;
  s->dtdy = s->dt / g->dy;
  s->dtdxx = s->dt / (g->dx * g->dx);
  s->dtdyy = s->dt / (g->dy * g->dy);
  s->dtdxdy = s->dt * g->dx * g->dy;

  /* set neighbors */
  if (MASTER) {
    s->prev = MPI_PROC_NULL;
  } else {
    s->prev = rank - 1;
  }

  if (LAST_NODE) {
    s->next = MPI_PROC_NULL;
  } else {
    s->next = rank + 1;
  }
}

/* Apply initial conditions*/
void set_init(struct FieldPointers *f, struct Grid2D *g,
              struct SimulationInfo *s) {
  for (int i = 1; i < g->nx_p - 1; i++) {
    f->un[i][g->ny] = s->ubc[0];
    f->un[i][g->ny - 1] = s->ubc[0];
  }
}

/* Applying boundary conditions for velocity */
void set_BC(struct FieldPointers *f, struct Grid2D *g, struct SimulationInfo *s,
            int rank, int nprocs) {
  int i, j, tag = 0;
  MPI_Status status;

  /* Set virtual boundary conditions */
  MPI_Sendrecv(&f->un[g->nx_p - 2][1], g->ny, MPI_DOUBLE, s->next, tag,
               &f->un[g->nx_p - 1][1], g->ny, MPI_DOUBLE, s->next, tag, WORLD,
               &status);
  MPI_Sendrecv(&f->un[1][1], g->ny, MPI_DOUBLE, s->prev, tag, &f->un[0][1],
               g->ny, MPI_DOUBLE, s->prev, tag, WORLD, &status);

  MPI_Sendrecv(&f->vn[g->nx_psg - 2][1], g->ny - 1, MPI_DOUBLE, s->next, tag,
               &f->vn[g->nx_psg - 1][1], g->ny - 1, MPI_DOUBLE, s->next, tag, WORLD,
               &status);
  MPI_Sendrecv(&f->vn[1][1], g->ny - 1, MPI_DOUBLE, s->prev, tag, &f->vn[0][1], g->ny - 1,
               MPI_DOUBLE, s->prev, tag, WORLD, &status);

  MPI_Sendrecv(&f->pn[g->nx_psg - 2][1], g->ny, MPI_DOUBLE, s->next, tag,
               &f->pn[g->nx_psg - 1][1], g->ny, MPI_DOUBLE, s->next, tag,
               WORLD, &status);
  MPI_Sendrecv(&f->pn[1][1], g->ny, MPI_DOUBLE, s->prev, tag, &f->pn[0][1],
               g->ny, MPI_DOUBLE, s->prev, tag, WORLD, &status);

  /* Set physical boundary conditions */
  /* Sides */
  if (MASTER) {
    for (j = 0; j < g->ny; j++) {
      f->un[0][j] = s->ubc[1];
      f->pn[0][j] = f->pn[1][j] - g->dx * s->pbc[1];
    }
    for (j = 0; j < g->ny - 1; j++) {
      f->vn[0][j] = 2.0 * s->vbc[1] - f->vn[1][j];
    }
  }

  if (LAST_NODE) {
    for (j = 0; j < g->ny; j++) {
      f->un[g->nx_p - 1][j] = s->ubc[3];
      f->pn[g->nx_psg - 1][j] = f->pn[g->nx_psg - 2][j] - g->dx * s->pbc[3];
    }
    for (j = 0; j < g->ny - 1; j++) {
      f->vn[g->nx_psg - 1][j] = 2.0 * s->vbc[3] - f->vn[g->nx_psg - 2][j];
    }
  }

  /* Bottom */
  for (i = 0; i < g->nx_p; i++)
    f->un[i][0] = 2.0 * s->ubc[2] - f->un[i][1];
  for (i = 0; i < g->nx_psg; i++) {
    f->vn[i][0] = s->vbc[2];
    f->pn[i][0] = f->pn[i][1] - g->dy * s->pbc[2];
  }

  /* Top */
  for (i = 0; i < g->nx_p; i++)
    f->un[i][g->ny] = 2.0 * s->ubc[0] - f->un[i][g->ny - 1];
  for (i = 0; i < g->nx_psg; i++) {
    f->vn[i][g->ny - 1] = s->vbc[0];
    f->pn[i][g->ny] = f->pn[i][g->ny - 1] - g->dy * s->pbc[0];
  }
}

/* Solve momentum for computing u and v */
void solve_U(struct FieldPointers *f, struct Grid2D *g,
             struct SimulationInfo *s) {
  int i, j;

#pragma omp parallel for private(i, j) schedule(auto)
  for (i = 1; i < g->nx_p - 1; i++) {
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
  for (i = 1; i < g->nx_psg - 1; i++) {
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
  for (i = 1; i < g->nx_psg - 1; i++) {
    for (j = 1; j < g->ny; j++) {
      f->pn[i][j] =
          f->p[i][j] - s->c2 * ((f->un[i][j] - f->un[i - 1][j]) * s->dtdx +
                                (f->vn[i][j] - f->vn[i][j - 1]) * s->dtdy);
    }
  }
}

/* Compute L2-norm */
void l2_norm(struct FieldPointers *f, struct Grid2D *g,
             struct SimulationInfo *s, int itr, int rank,
             int nprocs) {
  int i, j;
  double err_u = 0.0, err_v = 0.0, err_p = 0.0, err_d = 0.0;
  
  for (i = 0; i < 5; i++)
    s->errs[i] = 0;
#pragma omp parallel for private(i,j) schedule(auto) \
                             reduction(+:err_u, err_v, err_p, err_d)
  for (i = 1; i < g->nx_p - 1; i++) {
    for (j = 1; j < g->ny - 1; j++) {
      err_u += pow(f->un[i][j] - f->u[i][j], 2);
      err_v += pow(f->vn[i][j] - f->v[i][j], 2);
      err_p += pow(f->pn[i][j] - f->p[i][j], 2);
      err_d += (f->un[i][j] - f->un[i - 1][j]) * s->dtdx +
               (f->vn[i][j] - f->vn[i][j - 1]) * s->dtdy;
    }
  }
  MPI_Reduce(&err_u, &s->errs[1], 1, MPI_DOUBLE, MPI_SUM, 0, WORLD);
  MPI_Reduce(&err_v, &s->errs[2], 1, MPI_DOUBLE, MPI_SUM, 0, WORLD);
  MPI_Reduce(&err_p, &s->errs[3], 1, MPI_DOUBLE, MPI_SUM, 0, WORLD);
  MPI_Reduce(&err_d, &s->errs[4], 1, MPI_DOUBLE, MPI_SUM, 0, WORLD);

  if (MASTER) {
    s->errs[1] = sqrt(s->dtdxdy * s->errs[1]);
    s->errs[2] = sqrt(s->dtdxdy * s->errs[2]);
    s->errs[3] = sqrt(s->dtdxdy * s->errs[3]);
    s->errs[4] = fabs(s->errs[4]);

    int count = 4;
    s->errs[0] = fmaxof(count, s->errs[1], s->errs[2], s->errs[3], s->errs[4]);
  }
}
