/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
/*                                                                          */
/*                         PI_MD:                                           */
/*             The future of simulation technology                          */
/*             ------------------------------------                         */
/*                   Module: communicate_utils_pimd.c                       */
/*                                                                          */
/* This routine communicates energy terms                                   */
/*                                                                          */
/*                                                                          */
/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

#include "standard_include.h"
#include "../typ_defs/typedefs_class.h"
#include "../typ_defs/typedefs_gen.h"
#include "../typ_defs/typedefs_bnd.h"
#include "../proto_defs/proto_integrate_pimd_local.h"
#include "../proto_defs/proto_communicate_wrappers.h"

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void communicate_utils_pimd(double *kinet, COMMUNICATE *communicate)

/*==========================================================================*/
/*         Begin Routine                                                    */ 
     {/*begin routine*/
/*==========================================================================*/

#include "../typ_defs/typ_mask.h"


     double kinet_temp;
     MPI_Comm comm_beads;

     comm_beads = communicate->comm_beads;
     kinet_temp = *kinet;

     Allreduce(&kinet_temp,kinet, 1, MPI_DOUBLE, MPI_SUM,0, comm_beads);

}/*end routine*/





