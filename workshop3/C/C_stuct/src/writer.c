#include "writer.h"

/* Save fields data to files */
void dump_data(struct Grid2D *g, struct FieldPointers *f) {
  int i, j, count;
  FILE *fd;
  /* Local arrays on each process for storing fields at grid points */
  double **ug, **vg, **pg;

  ug = array_2D(g->nx, g->ny);
  vg = array_2D(g->nx, g->ny);
  pg = array_2D(g->nx, g->ny);

  /* #pragma omp parallel for private(i, j) schedule(auto) */
  for (i = 0; i < g->nx; i++) {
    for (j = 0; j < g->ny; j++) {
      ug[i][j] = 0.5 * (f->u[i][j + 1] + f->u[i][j]);
      vg[i][j] = 0.5 * (f->v[i + 1][j] + f->v[i][j]);
      pg[i][j] = 0.25 * (f->p[i][j] + f->p[i + 1][j] + f->p[i][j + 1] +
                         f->p[i + 1][j + 1]);
    }
  }

  count = 6;
  freeMem(count, g->ubufo, g->vbufo, g->pbufo, g->ubufn, g->vbufn, g->pbufn);

  fd = fopen("data/xyuvp", "w+t+e");

  /* Writing the field data from MASTER */
  fprintf(fd, "# X \t Y \t U \t V \t P\n");
  for (i = 0; i < g->nx; i++) {
    for (j = 0; j < g->ny; j++) {
      fprintf(fd, "%.8lf \t %.8lf \t %.8lf \t %.8lf \t %.8lf\n",
              (double)i * g->dx, (double)j * g->dy, ug[i][j], vg[i][j],
              pg[i][j]);
    }
  }

  count = 3;
  freeMem(count, ug, vg, pg);
  fclose(fd);
}
