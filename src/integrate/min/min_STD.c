/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
/*                                                                          */
/*                         PI_MD:                                           */
/*             The future of simulation technology                          */
/*             ------------------------------------                         */
/*                     Module: min_CG                                       */
/*                                                                          */
/* This subprogram minimizes a classical class using steepest descent      */
/*                                                                          */
/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/


#include "standard_include.h"
#include "../typ_defs/typedefs_class.h"
#include "../typ_defs/typedefs_bnd.h"
#include "../typ_defs/typedefs_gen.h"
#include "../proto_defs/proto_integrate_min_entry.h"
#include "../proto_defs/proto_integrate_min_local.h"
#include "../proto_defs/proto_intra_con_entry.h"
#include "../proto_defs/proto_energy_ctrl_entry.h"

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
void min_STD(CLASS *class,BONDED *bonded,GENERAL_DATA *general_data)
/*========================================================================*/
{/*begin routine*/
/*========================================================================*/
/*             Local variable declarations                                */
    int ipart,ifirst=1;
    double dt,tol_glob,dt22,dt2;

/* Local Pointers         */

   int min_atm_com_fix_opt = general_data->minopts.min_atm_com_fix_opt;

   double *clatoms_x    = class->clatoms_pos[1].x;
   double *clatoms_y    = class->clatoms_pos[1].y;
   double *clatoms_z    = class->clatoms_pos[1].z;
   double *clatoms_fx    = class->clatoms_pos[1].fx;
   double *clatoms_fy    = class->clatoms_pos[1].fy;
   double *clatoms_fz    = class->clatoms_pos[1].fz;
   double *clatoms_mass  = class->clatoms_info.mass;
   double *clatoms_xold = class->clatoms_info.xold;
   double *clatoms_yold = class->clatoms_info.yold;
   double *clatoms_zold = class->clatoms_info.zold;
/*==========================================================================*/
/* 0) Useful constants                                                      */

    (class->energy_ctrl.int_res_tra)    = 0;
    (class->energy_ctrl.int_res_ter)    = 0;
    dt   = (general_data->timeinfo.dt);
    dt2  = dt/2.0;
    dt22 = dt*dt/2.0;
    general_data->stat_avg.iter_shake = 0;
    general_data->stat_avg.iter_ratl = 0;

/*==========================================================================*/
/* I)Save positions                                                         */

    for(ipart=1;ipart<=(class->clatoms_info.natm_tot);ipart++){
      clatoms_xold[ipart] = clatoms_x[ipart];
      clatoms_yold[ipart] = clatoms_y[ipart];
      clatoms_zold[ipart] = clatoms_z[ipart];
    /*endfor*/}

/*==========================================================================*/
/* II) Get new energy and forces                                            */

    (class->energy_ctrl.iget_full_inter)= 1;
    (class->energy_ctrl.iget_res_inter) = 0;
    (class->energy_ctrl.iget_full_intra)= 1;
    (class->energy_ctrl.iget_res_intra) = 0;
    energy_control(class,bonded,general_data);

    if(min_atm_com_fix_opt==1){
       proj_com_out(class->clatoms_info.natm_tot,
                    clatoms_fx,clatoms_fy,clatoms_fz);
    }/*endif*/

/*==========================================================================*/
/* III) Evolve positions                                                     */


    for(ipart=1;ipart<=(class->clatoms_info.natm_tot);ipart++){
      clatoms_x[ipart] = clatoms_xold[ipart] + clatoms_fx[ipart]*dt
                                 /clatoms_mass[ipart];
      clatoms_y[ipart] = clatoms_yold[ipart] + clatoms_fy[ipart]*dt
                                 /clatoms_mass[ipart];
      clatoms_z[ipart] = clatoms_zold[ipart] + clatoms_fz[ipart]*dt
                                 /clatoms_mass[ipart];
    /*endfor*/}

/*==========================================================================*/
/* IV) Shake if necessary: Put particles on surface and add in constraint  */
/*                          force so the force test in control_min          */
/*                          will converge.                                  */

    if((bonded->constrnt.iconstrnt)==1){
         (bonded->constrnt.iroll) = 0;
         shake_control(bonded,&(class->clatoms_info),&(class->clatoms_pos[1]), 
                   &(general_data->cell),&(general_data->ptens),
                   &(general_data->statepoint),
                   &(general_data->baro),&(general_data->par_rahman),
                   &(general_data->stat_avg),dt,&tol_glob,ifirst,
                   &(class->class_comm_forc_pkg),&(class->ewd_scr));
         for(ipart=1;ipart<=class->clatoms_info.natm_tot;ipart++){
           clatoms_fx[ipart] = (clatoms_x[ipart]-clatoms_xold[ipart])
                             *clatoms_mass[ipart]/dt;
           clatoms_fy[ipart] = (clatoms_y[ipart]-clatoms_yold[ipart])
                             *clatoms_mass[ipart]/dt;
           clatoms_fz[ipart] = (clatoms_z[ipart]-clatoms_zold[ipart])
                             *clatoms_mass[ipart]/dt;
         }
    /*endif*/}  

/*==========================================================================*/
/* V) Get positions of ghost atoms                                          */

   get_ghost_pos(&(class->clatoms_info),&(class->clatoms_pos[1]),
                 &(class->ghost_atoms));

/*--------------------------------------------------------------------------*/
/*end routine*/}
/*==========================================================================*/














