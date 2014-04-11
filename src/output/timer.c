#if defined(TIMER)

#include "standard_include.h"


#if defined(PARALLEL)

#include <mpi.h>

/*==========================================================================*/
void timer_start(char *label) {
/*==========================================================================*/

  int rank;

  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  if (rank == 0) {
    fprintf(stderr, "timer start | %.9f | %s\n", timer_time(), label);
  }

}


/*==========================================================================*/
void timer_stop(char *label) {
/*==========================================================================*/

  int rank;

  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  if (rank == 0) {
    fprintf(stderr, "timer stop  | %.9f | %s\n", timer_time(), label);
  }

}

#endif


/*==========================================================================*/
double timer_resolution(void) {
/*==========================================================================*/

  struct timespec ts;

  clock_getres(CLOCK_MONOTONIC, &ts);

  return (double)ts.tv_sec + 1e-9 * (double)ts.tv_nsec;
}


/*==========================================================================*/
double timer_time(void) {
/*==========================================================================*/

  struct timespec ts;

  clock_gettime(CLOCK_MONOTONIC, &ts);

  return (double)ts.tv_sec + 1e-9 * (double)ts.tv_nsec;
}

#endif
