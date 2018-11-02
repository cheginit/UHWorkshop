#ifndef SIMULATIONCONTROLS_H
#define SIMULATIONCONTROLS_H

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "globals.h"
#include "structs.h"
#include "utilities.h"

/* Applying boundary conditions for velocity */
void initialize(struct FieldPointers *f, struct Grid2D *g,
                struct SimulationInfo *s);

/* Set initial condition */
void set_init(struct FieldPointers *f, struct Grid2D *g,
              struct SimulationInfo *s);

/* Applying boundary conditions for velocity */
void set_BC(struct FieldPointers *f, struct Grid2D *g,
            struct SimulationInfo *s);

/* Update the fields to the new time step for the next iteration */
void update(struct FieldPointers *f);

/* Solve momentum for computing u and v */
void solve_U(struct FieldPointers *f, struct Grid2D *g,
             struct SimulationInfo *s);

/* Solves continuity equation for computing P */
void solve_P(struct FieldPointers *f, struct Grid2D *g,
             struct SimulationInfo *s);

/* Compute L2-norm */
void l2_norm(struct FieldPointers *f, struct Grid2D *g,
             struct SimulationInfo *s);

#endif /* SIMULATIONCONTROLS_H */
