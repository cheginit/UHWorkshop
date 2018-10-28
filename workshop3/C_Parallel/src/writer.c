#include "writer.h"

/* Save fields data to files */
void dump_data(struct Grid2D *g, struct FieldPointers *f, int rank,
               int nprocs) {
  int i, j, offset, inner, arr_size, count;
  MPI_Status status;

  /* Local arrays on each process for storing fields at grid points */
  double **ug, **vg, **pg;

  inner = g->nx / nprocs;

  ug = array_2D(inner, g->ny);
  vg = array_2D(inner, g->ny);
  pg = array_2D(inner, g->ny);

  offset = g->ghosts / 2;

  if (LAST_NODE)
    offset = g->ghosts;

  /* #pragma omp parallel for private(i, j) schedule(auto) */
  for (i = 0; i < inner; i++) {
    for (j = 0; j < g->ny; j++) {
      ug[i][j] = 0.5 * (f->u[i + offset][j + 1] + f->u[i + offset][j]);
      vg[i][j] = 0.5 * (f->v[i + offset + 1][j] + f->v[i + offset][j]);
      pg[i][j] = 0.25 * (f->p[i + offset][j] + f->p[i + offset + 1][j] +
                         f->p[i + offset][j + 1] + f->p[i + offset + 1][j + 1]);
    }
  }

  arr_size = g->ny * inner;
  if (!MASTER) {
    MPI_Send(&(ug[0][0]), arr_size, MPI_DOUBLE, 0, 0, WORLD);
    MPI_Send(&(vg[0][0]), arr_size, MPI_DOUBLE, 0, 0, WORLD);
    MPI_Send(&(pg[0][0]), arr_size, MPI_DOUBLE, 0, 0, WORLD);

    /* Free the memory */
    count = 3;
    freeMem(count, ug, vg, pg);
  } else {
    FILE *fd;
    /* Buffer array on MASTER to get fields values from other processors */
    double **ubuff, **vbuff, **pbuff;

    ubuff = array_2D(inner, g->ny);
    vbuff = array_2D(inner, g->ny);
    pbuff = array_2D(inner, g->ny);

    fd = fopen("data/xyuvp", "w+t+e");

    /* Writing the field data from MASTER */
    fprintf(fd, "# X \t Y \t U \t V \t P\n");
    for (i = 0; i < inner; i++) {
      for (j = 0; j < g->ny; j++) {
        fprintf(fd, "%.8lf \t %.8lf \t %.8lf \t %.8lf \t %.8lf\n",
                (double)i * g->dx, (double)j * g->dy, ug[i][j], vg[i][j],
                pg[i][j]);
      }
    }

    count = 6;
    freeMem(count, g->ubufo, g->vbufo, g->pbufo, g->ubufn, g->vbufn, g->pbufn);
    /* Writing the field data from other processors */
    for (int r = 1; r < nprocs; r++) {
      MPI_Recv(&(ubuff[0][0]), arr_size, MPI_DOUBLE, r, 0, WORLD, &status);
      MPI_Recv(&(vbuff[0][0]), arr_size, MPI_DOUBLE, r, 0, WORLD, &status);
      MPI_Recv(&(pbuff[0][0]), arr_size, MPI_DOUBLE, r, 0, WORLD, &status);
      for (i = 0; i < inner; i++) {
        for (j = 0; j < g->ny; j++) {
          fprintf(fd, "%.8lf \t %.8lf \t %.8lf \t %.8lf \t %.8lf\n",
                  (double)(i + inner * r + offset) * g->dx, (double)j * g->dy,
                   ubuff[i][j], vbuff[i][j], pbuff[i][j]);
        }
      }
    }
    count = 6;
    freeMem(count, ubuff, vbuff, pbuff, ug, vg, pg);
    fclose(fd);
  }

  MPI_Finalize();
}
