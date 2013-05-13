/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
/*                                                                          */
/*                         PI_MD:                                           */
/*             The future of simulation technology                          */
/*             ------------------------------------                         */
/*                   Module: cp_energy_control.c                            */
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
#include "../typ_defs/typedefs_cp.h"
#include "../proto_defs/proto_energy_ctrl_entry.h"
#include "../proto_defs/proto_energy_ctrl_local.h"
#include "../proto_defs/proto_energy_ctrl_cp_entry.h"
#include "../proto_defs/proto_energy_cp_entry.h"
#include "../proto_defs/proto_energy_ctrl_cp_local.h"
#include "../proto_defs/proto_intra_entry.h"
#include "../proto_defs/proto_real_space_entry.h"
#include "../proto_defs/proto_recip3d_entry.h"
#include "../proto_defs/proto_communicate_wrappers.h"

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void cp_energy_control(CLASS *class, BONDED *bonded, 
                       GENERAL_DATA *general_data, CP *cp)

/*=======================================================================*/

{/*Begin Routine*/

/*=======================================================================*/
/*         Local Variable declarations                                   */
#include "../typ_defs/typ_mask.h"

  int iii;
  int myid_state       = class->communicate.myid_state;
  MPI_Comm world       = class->communicate.world;
  int simopt_cp        = general_data->simopts.cp;
  int simopt_cp_min    = general_data->simopts.cp_min;
  int np_forc          = class->communicate.np_forc; /*DY*/
  int nproc            = class->communicate.np;

/*======================================================================*/
/* I)   Initialize energies and forces                                  */

  energy_control_initial(class,bonded,general_data);

/*======================================================================*/
/* II)  Get electron contribution to the atm force and total E          */
/*      Get forces on basis set parameters                              */

  energy_control_elec(class,bonded,general_data,cp);


/*======================================================================*/
/*======================================================================*/
/*  ONLY GET CLASSICAL FORCES AND ENERGY IF NECESSARY                   */

  if(((simopt_cp+simopt_cp_min)==1) &&((myid_state==0)||(np_forc > 1))){

/*======================================================================*/
/* III) calculate inter_real energy                                     */

    energy_control_inter_real(class,bonded,general_data);

/*======================================================================*/
/* IV) Get surface PE                                                   */

    energy_control_surf(class,bonded,general_data);

/*======================================================================*/
/* V) calculate intra energy                                          */

    energy_control_intra(class,bonded,general_data);

  }/*endif: Classical stuff */
  if(nproc>1){Barrier(world);}
/*======================================================================*/
/*======================================================================*/
/* VI) Finalize energies                                               */

  energy_control_final(class,bonded,general_data);

/*------------------------------------------------------------------------*/
   }/*endroutine*/
/*=======================================================================*/








