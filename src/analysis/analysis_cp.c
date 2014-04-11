/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
/*                                                                          */
/*                         PI_MD:                                           */
/*             The future of simulation technology                          */
/*             ------------------------------------                         */
/*                     Module: analysis_cp                                  */
/*                                                                          */
/* This subprogram performs on the fly analysis of CP data                  */
/*                                                                          */
/*                                                                          */
/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

#include "standard_include.h"
#include "../typ_defs/typedefs_gen.h"
#include "../typ_defs/typedefs_class.h"
#include "../typ_defs/typedefs_bnd.h"
#include "../typ_defs/typedefs_cp.h"
#include "../typ_defs/typedefs_stat.h"
#include "../proto_defs/proto_analysis_cp_entry.h"
#include "../proto_defs/proto_analysis_md_entry.h"
#include "../proto_defs/proto_analysis_local_entry.h"
#include "../proto_defs/proto_communicate_wrappers.h"

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void analysis_cp(CLASS *class,BONDED *bonded,GENERAL_DATA *general_data,
                 CP *cp,ANALYSIS *analysis)

/*=======================================================================*/
{ /*begin routine*/

/*=======================================================================*/
/*            Local variable declarations                                */
/*=======================================================================*/

  int myid           = class->communicate.myid;
  int num_proc       = class->communicate.np;
  MPI_Comm comm_stat = class->communicate.comm_states;
  MPI_Comm world     = class->communicate.world;

 /*=========================================================================*/
 /* I- Gather the positions and the velocities on processor 0               */
 /*=========================================================================*/


 /*=========================================================================*/
 /* II- Do some preliminary stuff on processor 0                            */
 /*=========================================================================*/

  if(myid==0)
  {
   if (general_data->timeinfo.itime==1)
   {
    prelim_analysis(class,general_data,analysis);
   }; /* endif */
  }; /* endif */

 /*=========================================================================*/
 /* III- Calculate the radial distribution functions on processor 0         */
 /*=========================================================================*/

  if(myid==0)
  {
   if (analysis->rdf.calcul_gr_on==1)
   {
    calcul_gr(class,general_data,analysis);
   }; /* endif */
  }; /* endif */

 /*=========================================================================*/
 /* IV- Diagonalize the Hessian and compute frequencies on proc. 0         */
 /*=========================================================================*/

  if(myid==0)
  {
   if (analysis->harmonic_analysis.calcul_freq_on==1)
   {
     calcul_freqs(class,general_data,analysis);
   }; /* endif */
  }; /* endif */

 /*=========================================================================*/
 /* V- Special preliminary stuff....                                       */
 /*=========================================================================*/

 /*=========================================================================*/
 /* V- Waiting each others.....                                             */
 /*=========================================================================*/

 if(num_proc>1){  Barrier(world); }

/*=======================================================================*/
}/*end routine*/
/*==========================================================================*/


