#ifndef GLOBALS_H
#define GLOBALS_H

/* Specify number of grid points in x and y directions */
#define IX 128
#define IY 128

/* MPI variables */
#define MASTER (rank == 0)
#define NODE (rank != 0)
#define LAST_NODE (rank == nprocs - 1)
#define WORLD MPI_COMM_WORLD

#endif /* GLOBALS_H */
