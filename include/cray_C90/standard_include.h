#ifndef _CRAY
#define _CRAY
#endif
#define C90
#define SIMP_NINT
#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#ifdef PARALLEL
#include "mpi.h"
#else
#include "../typ_defs/mpi_f.h"
#endif
#include "../typ_defs/defines.h"




