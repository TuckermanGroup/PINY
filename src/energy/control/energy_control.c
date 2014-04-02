/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
/*                                                                          */
/*                         PI_MD:                                           */
/*             The future of simulation technology                          */
/*             ------------------------------------                         */
/*                   Module: energy_control.c                               */
/*                                                                          */
/* This routine calls the required force and PE routines                    */
/*                                                                          */
/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

#include "standard_include.h"
#include "../typ_defs/typedefs_gen.h"
#include "../typ_defs/typedefs_class.h"
#include "../typ_defs/typedefs_bnd.h"
#include "../proto_defs/proto_energy_ctrl_entry.h"
#include "../proto_defs/proto_intra_entry.h"
#include "../proto_defs/proto_real_space_entry.h"
#include "../proto_defs/proto_recip3d_entry.h"
#include "../proto_defs/proto_energy_ctrl_local.h"
#include "../proto_defs/proto_intra_con_entry.h"
#include "../proto_defs/proto_math.h"


/*========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*========================================================================*/

void energy_control(CLASS *class, BONDED *bonded, 
                    GENERAL_DATA *general_data)

/*========================================================================*/
  { /* Begin Routine */
/*========================================================================*/

  int    iii,i;
  double pvten_temp[10];
  int    np_forc     = class->communicate.np_forc;
  int    myid        = class->communicate.myid;
  MPI_Comm comm_forc = class->communicate.comm_forc;

/*======================================================================*/
/* I) Initialize energy variables                                       */
  
  energy_control_initial(class,bonded,general_data);

/*======================================================================*/
/* II) Get intermolecular real space force and  PE                      */

  energy_control_inter_real(class,bonded,general_data);

/*======================================================================*/
/* III) Get intermolecular recip space force and  PE                    */

  energy_control_inter_recip(class,bonded,general_data); 

/*======================================================================*/
/* IV) Get surface PE                                                   */

  energy_control_surf(class,bonded,general_data);

/*======================================================================*/
/* V) Get intramolecular bond force and  PE                             */

  energy_control_intra(class,bonded,general_data);

/*======================================================================*/
/* VI) Finish energy routine                                            */

  energy_control_final(class,bonded,general_data);

/*------------------------------------------------------------------------*/
    }/*end routine */
/*========================================================================*/



#ifdef SAVE
  for(i=1;i<=9;i++){pvten_temp[i]=0.0;}
  Allreduce(&(general_data->ptens.pvten_tot[1]), &(pvten_temp[1]),9,MPI_DOUBLE,
            MPI_SUM,0,class->communicate.world);
  for(iii=0;iii<=np_forc;iii++){
   Barrier(class->communicate.world);
   if(iii==myid){
     for(i=1;i<=9;i++){printf("after intra %d %.12g\n",i,pvten_temp[i]);}}
  }
  if(myid==0){scanf("%d",&iii);}
#endif






