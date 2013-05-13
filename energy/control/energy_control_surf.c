/*========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*========================================================================*/
/*                                                                        */
/*                       PI_MD:                                           */
/*           The future of simulation technology                          */
/*           ------------------------------------                         */
/*                 Module: energy_control_surf.c                          */
/*                                                                        */
/* This routine calls the required force and PE routines                  */
/*                                                                        */
/*========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*========================================================================*/

#include "standard_include.h"
#include "../typ_defs/typedefs_gen.h"
#include "../typ_defs/typedefs_class.h"
#include "../typ_defs/typedefs_bnd.h"
#include "../proto_defs/proto_energy_ctrl_local.h"
#include "../proto_defs/proto_surf_entry.h"

/*========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*========================================================================*/

void energy_control_surf(CLASS *class, BONDED *bonded, 
                         GENERAL_DATA *general_data)

/*========================================================================*/
{/*Begin Routine*/
/*=======================================================================*/
/*         Local Variable declarations                                   */

  int ip,iver_get;
  double vsurf;
  
  int iget_full_inter   = class->energy_ctrl.iget_full_inter;
  int iget_full_intra   = class->energy_ctrl.iget_full_intra;

  int pi_beads          = class->clatoms_info.pi_beads;  
  int pi_beads_proc     = class->clatoms_info.pi_beads_proc;  

  int isurf_on          = class->surface.isurf_on; 

  int iget_pe_real_inter= class->energy_ctrl.iget_pe_real_inter;
  int iget_pv_real_inter= class->energy_ctrl.iget_pv_real_inter;

  int pimd_on;
  pimd_on               =  general_data->simopts.pimd
                          +general_data->simopts.debug_pimd
                          +general_data->simopts.cp_pimd
                          +general_data->simopts.cp_wave_pimd
                          +general_data->simopts.cp_wave_min_pimd
                          +general_data->simopts.debug_cp_pimd;

/*======================================================================*/
/* I) Zero stuff and set some constants */

  vsurf    = 0.0;

  iver_get = 0;
  if((pimd_on==1)&&(iget_full_inter==1)){iver_get = 1;}

/*======================================================================*/
/* VII) Get intramolecular tors force and  PE                           */

  if(iget_full_intra==1){

    if(isurf_on != 0){
     for(ip=1;ip<=pi_beads_proc;ip++){
       surf_pot(&(class->clatoms_info),&(class->clatoms_pos[ip]),
                &(class->atommaps),&(class->surface),&(general_data->cell),
                &(bonded->intra_scr),&(general_data->ptens),&vsurf,iver_get,
                &(class->class_comm_forc_pkg),
                iget_pv_real_inter,iget_pe_real_inter);
     }/*endfor*/
     vsurf/=pi_beads;
    }/*endif*/

  }/*endif*/

/*======================================================================*/
/* IX) Collect and store PE                                           */
  
  if( iget_full_intra==1 ){

    if(iget_pe_real_inter==1){
      general_data->stat_avg.vvdw       += vsurf;
      general_data->stat_avg.vintert    += (vsurf);
    }/*endif*/

    general_data->stat_avg.vsurft        = vsurf;

  }/*endif*/

/*-----------------------------------------------------------------------*/
   }/*end routine */
/*==========================================================================*/
