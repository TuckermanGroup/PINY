/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
/*                                                                          */
/*                         PI_MD:                                           */
/*             The future of simulation technology                          */
/*             ------------------------------------                         */
/*                     Module: analysis_md                                  */
/*                                                                          */
/* This subprogram performs on the fly analysis of MD data                  */
/*                                                                          */
/*                                                                          */
/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

#include "standard_include.h"
#include "../typ_defs/typedefs_gen.h"
#include "../typ_defs/typedefs_class.h"
#include "../typ_defs/typedefs_bnd.h"
#include "../typ_defs/typedefs_stat.h"
#include "../proto_defs/proto_analysis_md_entry.h"
#include "../proto_defs/proto_communicate_wrappers.h"

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void analysis_pimd(CLASS *class,GENERAL_DATA *general_data,
		   BONDED *bonded, ANALYSIS *analysis)
/*==========================================================================*/
{ /*begin routine*/
 /*=========================================================================*/
 /*              Local variable declarations                                */
 /*=========================================================================*/

   MPI_Comm world = class->communicate.world;
  int num_proc       = class->communicate.np;

 /*=========================================================================*/

 /*=========================================================================*/
 /* I- Do some preliminary stuff on each processor                          */
 /*=========================================================================*/

  if (general_data->timeinfo.itime==1)
  {
   prelim_analysis(class,general_data,analysis);
  }; /* endif */

 /*=========================================================================*/
 /* II- Caculate the radial distribution functions on each processor        */
 /*=========================================================================*/
  
  if (analysis->rdf.calcul_gr_on==1)
  {
   calcul_gr(class,general_data,analysis);
  }; /* endif */
  
 /*=========================================================================*/
 /* III- Calculate the velocity correlation functions on processor 0        */
 /*      using classical operator (centroid = bead No 1) limit              */
 /*=========================================================================*/

  if (class->communicate.myid==0)
  {
   if (analysis->velocorel.calcul_atm_on==1)
   {
    calcul_vovt_atm(class,general_data,analysis);
   }; /* endif */
  }; /* endif */

 /*=========================================================================*/
 /* IV- Calculate the R, V, R**2, and V**2 semiclassical operators          */
 /*     and their respa-pimd time average                                   */
 /*=========================================================================*/

#ifdef TOTO
   calcul_simp_sc_op(class,general_data,analysis);
#endif

 /*=========================================================================*/
 /* IV- Calculate the R, V, R**2, and V**2 semiclassical operators          */
 /*     and their respa-pimd time average                                   */
 /*=========================================================================*/

#ifdef OK
   calcul_vovt_bead(class,general_data,analysis);
#endif

 /*=========================================================================*/
 /* V- Calculate the incoherent scattering functions on processor 0         */
 /*=========================================================================*/

  if(class->communicate.myid==0)
  {
   if (analysis->iikt_iso_corel.calcul_on==1)
   {
    calcul_iikt_iso_cmd(class,general_data,analysis);
   }; /* endif */
  }; /* endif */

 /*=========================================================================*/
 /* IV- Waiting each others.....                                            */
 /*=========================================================================*/

  if(num_proc>1){ Barrier(world); }

 /*=========================================================================*/
} /*end routine*/
/*==========================================================================*/

