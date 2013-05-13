/*===================================================================*/
/*ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*===================================================================*/
/*                                                                   */
/*                         PI_MD:                                    */
/*             The future of simulation technology                   */
/*             ------------------------------------                  */
/*                Module: control_vnhc_smpl.c                        */
/*                                                                   */
/* Control routine for velocity resampling                           */
/*                                                                   */
/*===================================================================*/
/*ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*===================================================================*/

#include "standard_include.h"
#include "../typ_defs/typedefs_gen.h"
#include "../typ_defs/typedefs_class.h"
#include "../typ_defs/typedefs_bnd.h"
#include "../typ_defs/typedefs_par.h"
#include "../proto_defs/proto_vel_sampl_class_entry.h"
#include "../proto_defs/proto_vel_sampl_class_local.h"
#include "../proto_defs/proto_communicate_wrappers.h"

/*===================================================================*/
/*ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*===================================================================*/

void control_vnhc_smpl(CLASS *class,GENERAL_DATA *general_data)

/*===================================================================*/
{/*begin routine*/

  int error_check_on  = general_data->error_check_on;
  MPI_Comm comm_beads = class->communicate.comm_beads;

/*===================================================================*/
/* 0) Write to screen                                                */

if(error_check_on==1){
  PRINT_LINE_STAR;
  printf("Sampling extended classical velocities\n");
  PRINT_LINE_DASH;printf("\n");
}/*endif*/

/*===================================================================*/
/* I) Sample extended class velocities                               */

   sampl_vnhc(&(class->therm_info_bead),(class->therm_bead), 
              &(class->therm_info_class),&(class->therm_class), 
              &(general_data->baro),&(general_data->par_rahman),
	      &(general_data->ensopts),&(general_data->statepoint), 
              &(class->int_scr),(class->clatoms_info.pi_beads), 
	      &(class->vel_samp_class.iseed), &(class->vel_samp_class.iseed2),
	      &(class->vel_samp_class.qseed),(general_data->cell.iperd),
              (class->communicate.myid),(class->clatoms_info.pi_beads_proc),
              (general_data->cell.hmat_int_typ),
              (general_data->cell.hmat_cons_typ));

/*====================================================================*/
/* II) Write to screen                                                */

if(error_check_on==1){
  PRINT_LINE_DASH;
  printf("Extended classical velocity sampling complete\n");
  PRINT_LINE_STAR;printf("\n");
}/*endif*/

/*====================================================================*/
  }/*end routine*/
/*====================================================================*/







