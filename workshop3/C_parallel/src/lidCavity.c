/*============================================================================*\
Solves Navier-Stokes equations for incompressible, laminar, steady flow using
artificial compressibility method on staggered grid.
The governing equations are as follows:

P_t + c^2 div[u] = 0
u_t + u . grad[u] = - grad[P] + nu div[grad[u]]

where P is p/rho and c represents artificial sound's speed.

Lid-Driven Cavity case:
Dimensions : 1x1 m
Grid size  : 128 x 128
Re number  : 100 / 1000 / 5000 / 10000
Grid type  : Staggered Arakawa C

Boundary Conditions: u, v -> Dirichlet (as shown below)
                     p    -> Neumann (grad[p] = 0)
                                u=1, v=0
                             ---------------
                            |               |
                        u=0 |               | u=0
                        v=0 |               | v=0
                            |               |
                            |               |
                             ---------------
                                  u=0, v=0
\*============================================================================*/
#include "simulationControls.h"
#include "writer.h"

int main(int argc, char *argv[]) {
  int itr = 1, count = 6;
  const double tol = 1.0e-6;
  const int itr_max = 1000000;

  FILE *flog;

  int rank, nprocs, provided;

  MPI_Init_thread(&argc, &argv, MPI_THREAD_FUNNELED, &provided);
  MPI_Comm_size(WORLD, &nprocs);
  MPI_Comm_rank(WORLD, &rank);

  g.nx = IX;
  g.ny = IY;

  /* Boundary conditions: {top, left, bottom, right} */
  s = ((struct SimulationInfo){.ubc = {1.0, 0.0, 0.0, 0.0},
                               .vbc = {0.0, 0.0, 0.0, 0.0},
                               .pbc = {0.0, 0.0, 0.0, 0.0}});

  s.Re = 100.0;
  s.l_lid = 1.0;

  /* Getting Reynolds number */
  if (MASTER) {
    if (argc > 1) {
      char *ptr;
      s.Re = strtod(argv[1], &ptr);
    }
    printf("Re number is set to %d\n", (int)s.Re);

    /* Create a log file for outputting the residuals */
    flog = fopen("data/residual", "w+t+e");
  }
  MPI_Bcast(&s.Re, 1, MPI_DOUBLE, 0, WORLD);

  initialize(&f, &g, &s, rank, nprocs);
  set_init(&f, &g, &s);
  set_BC(&f, &g, &s, rank, nprocs);
  update(&f);

  /* Start the main loop */
  do {
    solve_U(&f, &g, &s);
    solve_P(&f, &g, &s);
    set_BC(&f, &g, &s, rank, nprocs);
    l2_norm(&f, &g, &s, rank, nprocs);

    if (MASTER) {
      /* Check if solution diverged */
      if (isnan(s.errs[0])) {
        printf("Solution Diverged after %d iterations!\n", itr);
        /* Free the memory and terminate */
        count = 6;
        freeMem(count, g.ubufo, g.vbufo, g.pbufo, g.ubufn, g.vbufn, g.pbufn);
        exit(EXIT_FAILURE);
      }
      fprintf(flog, "%d\t%.8lf\t%.8lf\t%.8lf\t%.8lf\t%.8lf\n", itr, s.errs[0],
              s.errs[1], s.errs[2], s.errs[3], s.errs[4]);
      itr += 1;
    }
    MPI_Bcast(&s.errs[0], 1, MPI_DOUBLE, 0, WORLD);
    MPI_Bcast(&itr, 1, MPI_INT, 0, WORLD);

    /* Update the fields */
    update(&f);

  } while (s.errs[0] > tol && itr < itr_max);

  if (MASTER) {
    if (itr == itr_max) {
      printf("Maximum number of iterations, %d, exceeded\n", itr);

      /* Free the memory and terminate */
      count = 6;
      freeMem(count, g.ubufo, g.vbufo, g.pbufo, g.ubufn, g.vbufn, g.pbufn);
      exit(EXIT_FAILURE);
    }

    printf("Converged after %d iterations\n", itr);
    fclose(flog);
  }

  /* Write output data */
  dump_data(&g, &f, rank, nprocs);
  return 0;
}
