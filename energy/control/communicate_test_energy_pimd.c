/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
/*                                                                          */
/*                         PI_MD:                                           */
/*             The future of simulation technology                          */
/*             ------------------------------------                         */
/*                   Module: communicate_test_energy_pimd.c                 */
/*                                                                          */
/* This routine communicates energy terms                                   */
/*                                                                          */
/*                                                                          */
/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

#include "standard_include.h"
#include "../typ_defs/typedefs_gen.h"
#include "../typ_defs/typedefs_class.h"
#include "../typ_defs/typedefs_bnd.h"
#include "../typ_defs/typedefs_par.h"
#include "../typ_defs/typedefs_cp.h"
#include "../proto_defs/proto_energy_ctrl_local.h"
#include "../proto_defs/proto_communicate_wrappers.h"

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void communicate_test_energy_pimd(double *vgen,double *vgen2,double *vgen3,
                                      double *cpu, MPI_Comm world)

/*==========================================================================*/
/*         Begin Routine                                                    */ 
     {/*begin routine*/
/*==========================================================================*/

#include "../typ_defs/typ_mask.h"

     int iii;
     double vgen_temp,vgen2_temp,vgen3_temp,cpu_temp;

     cpu_temp = *cpu;
     vgen_temp = *vgen;
     vgen2_temp = *vgen2;
     vgen3_temp = *vgen3;
     Reduce(&cpu_temp,cpu, 1, MPI_DOUBLE, MPI_MAX, 0,world);
     Allreduce(&vgen_temp,vgen, 1, MPI_DOUBLE, MPI_SUM,0, world);
     Allreduce(&vgen2_temp,vgen2, 1, MPI_DOUBLE, MPI_SUM,0, world);
     Allreduce(&vgen3_temp,vgen3, 1, MPI_DOUBLE, MPI_SUM,0, world);
     Barrier(world);



}/*end routine*/





