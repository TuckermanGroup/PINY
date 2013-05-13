/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
/*                                                                          */
/*                         PI_MD:                                           */
/*             The future of simulation technology                          */
/*             ------------------------------------                         */
/*                     Module: mall_coord                                   */
/*                                                                          */
/* This subprogram reads atm-atm_NHC vol-vol_NHC input for a MD on a        */ 
/* LD-classical potential energy surface (LD-PES)                           */
/*                                                                          */
/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

#include "standard_include.h"
#include "../typ_defs/typedefs_class.h"
#include "../typ_defs/typedefs_par.h"
#include "../typ_defs/typedefs_gen.h"
#include "../typ_defs/typedefs_bnd.h"
#include "../proto_defs/proto_coords_local.h"
#include "../proto_defs/proto_friend_lib_entry.h"

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void mall_coord(CLASS *class,GENERAL_DATA *general_data)

/*======================================================================*/
/*                Begin Routine */
{   /*begin routine */
/*======================================================================*/
/*               Local variable declarations                            */

  int i,pimd_on,iii;
  int pi_beads      = class->clatoms_info.pi_beads;
  int pi_beads_proc = class->clatoms_info.pi_beads_proc;
  int natm_mall     = class->clatoms_info.natm_mall;
  int natm_tot      = class->clatoms_info.natm_tot;
  int nhess_mall; 

/*========================================================================*/
/* I) Malloc the vector structures                                        */

  pimd_on = general_data->simopts.pimd 
          + general_data->simopts.cp_pimd 
          + general_data->simopts.cp_wave_pimd 
          + general_data->simopts.cp_wave_min_pimd 
          + general_data->simopts.debug_pimd 
          + general_data->simopts.debug_cp_pimd;

  if(pi_beads>1||pimd_on==1){   
   class->therm_bead      = (THERM_POS *)cmalloc(pi_beads *
                              sizeof(THERM_POS))-1;
  }/*endif*/
  class->clatoms_pos      = (CLATOMS_POS *)cmalloc(pi_beads *
                             sizeof(CLATOMS_POS))-1;

/*========================================================================*/
/* I) Malloc atom positions, velocities, and forces                       */


  for(i=1;i<=pi_beads;i++){
  (class->clatoms_pos)[i].x   = (double *)cmalloc(natm_mall*sizeof(double))-1;
  (class->clatoms_pos)[i].y   = (double *)cmalloc(natm_mall*sizeof(double))-1;
  (class->clatoms_pos)[i].z   = (double *)cmalloc(natm_mall*sizeof(double))-1;
  (class->clatoms_pos)[i].vx  = (double *)cmalloc(natm_mall*sizeof(double))-1;
  (class->clatoms_pos)[i].vy  = (double *)cmalloc(natm_mall*sizeof(double))-1;
  (class->clatoms_pos)[i].vz  = (double *)cmalloc(natm_mall*sizeof(double))-1;

   if(pi_beads>1||pimd_on==1){
    (class->clatoms_pos)[i].fxt  = 
                                 (double *)cmalloc(natm_mall*sizeof(double))-1;
    (class->clatoms_pos)[i].fyt  = 
                                 (double *)cmalloc(natm_mall*sizeof(double))-1;
    (class->clatoms_pos)[i].fzt  = 
                                 (double *)cmalloc(natm_mall*sizeof(double))-1;
    (class->clatoms_pos)[i].fxm  = 
                                 (double *)cmalloc(natm_mall*sizeof(double))-1;
    (class->clatoms_pos)[i].fym  = 
                                 (double *)cmalloc(natm_mall*sizeof(double))-1;
    (class->clatoms_pos)[i].fzm  = 
                                 (double *)cmalloc(natm_mall*sizeof(double))-1;
    (class->clatoms_pos)[i].mass = 
                                 (double *)cmalloc(natm_mall*sizeof(double))-1;
   }
  (class->clatoms_pos)[i].fx  = (double *)cmalloc(natm_mall*sizeof(double))-1;
  (class->clatoms_pos)[i].fy  = (double *)cmalloc(natm_mall*sizeof(double))-1;
  (class->clatoms_pos)[i].fz  = (double *)cmalloc(natm_mall*sizeof(double))-1;
  }/*endfor*/

/*========================================================================*/
/* II) Malloc atomic hessian if necessary                                 */

  nhess_mall = (class->clatoms_info.hess_calc > 1 ? natm_tot*natm_tot : natm_tot);
  
  if(class->clatoms_info.hess_calc > 0){

    class->clatoms_pos[1].hess_xx = 
               (double *)cmalloc(nhess_mall*sizeof(double))-1;
    class->clatoms_pos[1].hess_xy = 
               (double *)cmalloc(nhess_mall*sizeof(double))-1;
    class->clatoms_pos[1].hess_xz = 
               (double *)cmalloc(nhess_mall*sizeof(double))-1;
    class->clatoms_pos[1].hess_yy = 
               (double *)cmalloc(nhess_mall*sizeof(double))-1;
    class->clatoms_pos[1].hess_yz = 
               (double *)cmalloc(nhess_mall*sizeof(double))-1;
    class->clatoms_pos[1].hess_zz = 
               (double *)cmalloc(nhess_mall*sizeof(double))-1;
  }/* end if */

/*========================================================================*/
    } /* end routine */
/*==========================================================================*/









