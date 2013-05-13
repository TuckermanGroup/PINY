/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
/*                                                                          */
/*                         PI_MD:                                           */
/*             The future of simulation technology                          */
/*             ------------------------------------                         */
/*               Module: energy_control_pimd.c                              */
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
#include "../proto_defs/proto_pimd_entry.h"
#include "../proto_defs/proto_pimd_local.h"
#include "../proto_defs/proto_math.h"
#include "../proto_defs/proto_communicate_wrappers.h"

#define DEBUG_OFF

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void energy_control_pimd(CLASS *class, BONDED *bonded, 
                         GENERAL_DATA *general_data)

/*==========================================================================*/
{/*Begin Routine*/

#include "../typ_defs/typ_mask.h"

  /*=======================================================================*/
  /*         Local Variable declarations                                   */
  
  int    i,iii,ip,ipart;
  double akvirx;
  double vbond_free,vbend_free,vtors_free;
  double *fx,*fy,*fz,*fxm,*fym,*fzm,*fxt,*fyt,*fzt;

  int    pi_beads        = class->clatoms_info.pi_beads; 
  int    pi_beads_proc   = class->clatoms_info.pi_beads_proc;
  int    myid            = class->communicate.myid;
  int    *recv_count_atm = class->class_comm_forc_pkg.recv_count_atm;
  double *fxtemp         = class->ewd_scr.x;
  double *fytemp         = class->ewd_scr.y;
  double *fztemp         = class->ewd_scr.z;
  double *pvten          = general_data->ptens.pvten;
  double *pvten_tot      = general_data->ptens.pvten_tot;
  int    natm_tot        = class->clatoms_info.natm_tot;
  int    myatm_start     = class->clatoms_info.myatm_start;
  int    myatm_end       = class->clatoms_info.myatm_end;
  int    np_forc         = class->communicate.np_forc;
  MPI_Comm comm_beads    = class->communicate.comm_beads;
  MPI_Comm comm_forc     = class->communicate.comm_forc;
  int iget_full_inter    = class->energy_ctrl.iget_full_inter;
  int iget_res_inter     = class->energy_ctrl.iget_res_inter;
  int iget_full_intra    = class->energy_ctrl.iget_full_intra;
  int iget_res_intra     = class->energy_ctrl.iget_res_intra;
  double volume          = general_data->cell.vol;
  double t_ext           = general_data->statepoint.t_ext;
  int bond_free_num      = bonded->bond_free.num;
  int bend_free_num      = bonded->bend_free.num;
  int tors_free_num      = bonded->tors_free.num;

/*========================================================================*/
/* 0) Zero stuff      */

  akvirx     = 0.0;
  vbond_free = 0.0; 
  vbend_free = 0.0; 
  vtors_free = 0.0; 

/*========================================================================*/
/* I) Initialize energy variables                                       */
  
  energy_control_initial(class,bonded,general_data);

/*======================================================================*/
/* II) Get intermolecular real space force and  PE   */
  
  energy_control_inter_real(class,bonded,general_data);

/*======================================================================*/
/* III) Get intermolecular recip space force and  PE */
  
  energy_control_inter_recip(class,bonded,general_data);

/*======================================================================*/
/* III) Get surface force and  PE */
  
  energy_control_surf(class,bonded,general_data);

/*======================================================================*/
/* IV) Get intramolecular bond force and  PE         */

  energy_control_intra(class,bonded,general_data);

/*======================================================================*/
/* V) Construct staging/centroid forces                                 */

  for(ip=1; ip<= pi_beads_proc; ip++){
   fx = class->clatoms_pos[ip].fx;
   fy = class->clatoms_pos[ip].fy;
   fz = class->clatoms_pos[ip].fz;
   for(i=1;i<= natm_tot;i++){
     fx[i] /= pi_beads;
     fy[i] /= pi_beads;
     fz[i] /= pi_beads;
   }/*endfor*/
  }/*endfor*/

  if( ((iget_full_inter)==1)){
   if(pi_beads>1){
    for(ip=1;ip<= pi_beads_proc;ip++){
     fxt = class->clatoms_pos[ip].fxt;
     fyt = class->clatoms_pos[ip].fyt;
     fzt = class->clatoms_pos[ip].fzt;
     for(i=1;i<= natm_tot;i++){
       fxt[i] /= pi_beads;
       fyt[i] /= pi_beads;
       fzt[i] /= pi_beads;
     }/*endfor*/
    }/*endfor*/
   }/*endif*/
  }/*endif*/
  for(i=1;i<=9;i++){
   pvten[i]     /= pi_beads;
   pvten_tot[i] /= pi_beads;
  }/*endfor*/

  /*----------------------------------------------------------------------*/
  /*  Add in mode forces and correct pressure tensor                      */

  if(pi_beads>1){pimd_pvten_forc_corr(class,general_data,&akvirx);}

  if(np_forc>1){
    for(ip=1;ip<= pi_beads_proc;ip++){
      fxt = class->clatoms_pos[ip].fxt;
      fyt = class->clatoms_pos[ip].fyt;
      fzt = class->clatoms_pos[ip].fzt;

      Allreduce(&(fxt[1]),&(fxtemp[1]),natm_tot,MPI_DOUBLE,
                           MPI_SUM,0,comm_forc);
      Allreduce(&(fyt[1]),&(fytemp[1]),natm_tot,MPI_DOUBLE,
                           MPI_SUM,0,comm_forc);
      Allreduce(&(fzt[1]),&(fztemp[1]),natm_tot,MPI_DOUBLE,
                           MPI_SUM,0,comm_forc);
      for(ipart=1;ipart<=natm_tot;ipart++){
        fxt[ipart] = fxtemp[ipart];
        fyt[ipart] = fytemp[ipart];
        fzt[ipart] = fztemp[ipart];
      }/*endfor : ipart*/
    }/*endfor : ip*/
  }/*endif*/


  /*----------------------------------------------------------------------*/
  /* Allreduce forces if force level par is on                            */


  if( (iget_full_inter+iget_res_inter)>=1){

    if(np_forc>1){
      for(ip=1;ip<=pi_beads_proc;ip++){
        fx = class->clatoms_pos[ip].fx;
        fy = class->clatoms_pos[ip].fy;
        fz = class->clatoms_pos[ip].fz;
        Reduce_scatter(&fx[1],&fxtemp[myatm_start],&recv_count_atm[1],
                       MPI_DOUBLE,MPI_SUM,comm_forc);
        Reduce_scatter(&fy[1],&fytemp[myatm_start],&recv_count_atm[1],
                       MPI_DOUBLE,MPI_SUM,comm_forc);
        Reduce_scatter(&fz[1],&fztemp[myatm_start],&recv_count_atm[1],
                       MPI_DOUBLE,MPI_SUM,comm_forc);
        for(ipart=myatm_start;ipart<=myatm_end;ipart++){
          fx[ipart] = fxtemp[ipart];
          fy[ipart] = fytemp[ipart];
          fz[ipart] = fztemp[ipart];
        }/*endfor : ipart*/
      }/*endfor : ip*/
    }/*endif : np_forc>1*/
   
  }/*endif : inter RESPA step*/

  /*----------------------------------------------------------------------*/
  /* Transform forces                                                     */

  if(((iget_full_intra)==1)|| ((iget_res_intra)==1)){
    control_pimd_trans_force(class,general_data);
  }/*endif*/

  for(ip=1;ip<=pi_beads_proc;ip++){
    fx = class->clatoms_pos[ip].fx;
    fy = class->clatoms_pos[ip].fy;
    fz = class->clatoms_pos[ip].fz;
    fxm = class->clatoms_pos[ip].fxm;
    fym = class->clatoms_pos[ip].fym;
    fzm = class->clatoms_pos[ip].fzm;
    for(i=myatm_start;i<=myatm_end;i++){
      fx[i] += fxm[i];
      fy[i] += fym[i];
      fz[i] += fzm[i];
    }/*endfor*/
  }/*endfor*/
  

  /*----------------------------------------------------------------------*/
  /*   Get kinetic contribution to the pressure                    */

   general_data->stat_avg.press_kin =  (natm_tot*t_ext)/(volume*BOLTZ);

/*=========================================================================*/
/* Mode free energy calculations */

  if(myid==0){
   if( (iget_full_intra==1) ||  (iget_res_intra==1) ){

     if((bond_free_num)!=0){
      bond_free_mode(&(class->clatoms_info),&(class->clatoms_pos[1]),
                &(bonded->bond_free),&(general_data->cell),
                &(bonded->intra_scr),&(general_data->ptens),&vbond_free,
                &(class->energy_ctrl),class->communicate.np_forc);
     }/*endif*/
     if((bend_free_num)!=0){
      bend_free_mode(&(class->clatoms_info),&(class->clatoms_pos[1]),
                &(bonded->bend_free),&(general_data->cell),
                &(bonded->intra_scr),&(general_data->ptens),&vbend_free,
                &(class->energy_ctrl),class->communicate.np_forc);
     }/*endif*/
     if((tors_free_num)!=0){
      tors_free_mode(&(class->clatoms_info),&(class->clatoms_pos[1]),
                &(bonded->tors_free),&(general_data->cell),
                &(bonded->intra_scr),&(general_data->ptens),&vtors_free,
                &(class->energy_ctrl),class->communicate.np_forc);
     }/*endif*/
     if((bond_free_num+bend_free_num+tors_free_num)!=0){
      distrib_ghost_force_mode(&(class->clatoms_info),&(class->clatoms_pos[1]),
                               &(class->ghost_atoms));
     }/*endif*/
   }/*endif*/
  }/*endif:myid=0*/

/*======================================================================*/
/* XI) Save the relevant energy terms                                   */

  if(iget_full_inter==1){

    general_data->stat_avg.vbond_free = vbond_free;
    general_data->stat_avg.vbend_free = vbend_free;
    general_data->stat_avg.vtors_free = vtors_free;
    general_data->stat_avg.vintrat   += (vbond_free+vbend_free+vtors_free);
    general_data->stat_avg.pi_ke_vir  = akvirx; 

 }/*endif*/

/*======================================================================*/
/* X) Finish energy routine                                             */

  energy_control_final(class,bonded,general_data);

/*-----------------------------------------------------------------------*/
   }/*end routine */
/*==========================================================================*/









