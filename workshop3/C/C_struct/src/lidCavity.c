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
  int itr = 1, count;
  const double tol = 1.0e-7;
  const int itr_max = 1000000;

  FILE *flog;

  g.nx = IX;
  g.ny = IY;

  /* Boundary conditions: {top, left, bottom, right} */
  s = ((struct SimulationInfo){.ubc = {1.0, 0.0, 0.0, 0.0},
                               .vbc = {0.0, 0.0, 0.0, 0.0},
                               .pbc = {0.0, 0.0, 0.0, 0.0}});

  s.Re = 100.0;
  s.l_lid = 1.0;

  /* Getting Reynolds number */
  if (argc > 1) {
    char *ptr;
    s.Re = strtod(argv[1], &ptr);
  }
  printf("Re number is set to %d\n", (int)s.Re);

  /* Create a log file for outputting the residuals */
  flog = fopen("data/residual", "w+t+e");
  fprintf(flog, "# iteration\ttotal\tu\tv\tp\tdivergence\n");

  initialize(&f, &g, &s);
  set_init(&f, &g, &s);

  /* Start the main loop */
  do {
    solve_U(&f, &g, &s);
    set_UBC(&f, &g, &s);
    solve_P(&f, &g, &s);
    set_PBC(&f, &g, &s);
    l2_norm(&f, &g, &s);

    /* Check if solution diverged */
    if (isnan(s.errs[0])) {
      printf("Solution Diverged after %d iterations!\n", itr);

      /* Free the memory and terminate */
      count = 6;
      freeMem(count, g.ubufo, g.vbufo, g.pbufo, g.ubufn, g.vbufn, g.pbufn);
      fclose(flog);
      exit(EXIT_FAILURE);
    }
    fprintf(flog, "%d\t%.8lf\t%.8lf\t%.8lf\t%.8lf\t%.8lf\n", itr, s.errs[0],
            s.errs[1], s.errs[2], s.errs[3], s.errs[4]);

    /* Update the fields */
    update(&f);
    itr += 1;
  } while (s.errs[0] > tol && itr < itr_max);

  if (itr == itr_max) {
    printf("Maximum number of iterations (%d) exceeded\n", itr);
  } else {
    printf("Converged after %d iterations\n", itr);
  }

  fclose(flog);

  /* Write output data */
  dump_data(&g, &f);
  return 0;
}
