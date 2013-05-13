#define HP_VECLIB
#define NO_PRAGMA
#define PARALLEL
#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <cps.h>
#include <sys/cnx_ail.h>
#ifdef PARALLEL
#include "mpi.h"
#else
#include "../typ_defs/mpi_f.h"
#endif
#include "../typ_defs/defines.h"




