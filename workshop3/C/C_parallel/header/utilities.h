#ifndef UTILITITES_H
#define UTILITITES_H

#include <math.h>
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>

#include "structs.h"

/* Generate a 2D array using pointer to pointer */
double **array_2D(int row, int col);

/* Update the fields to the new time step for the next iteration */
void update(struct FieldPointers *f);

/* Free the memory */
void freeMem(int count, ...);

/* Find mamximum of a set of float numebrs */
double fmaxof(int count, ...);

#endif /* UTILITITES_H */
