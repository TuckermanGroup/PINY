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
#include "../proto_defs/proto_pimd_entry.h"
#include "../proto_defs/proto_math.h"
/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void energy_control_inter_recip(CLASS *class, BONDED *bonded, 
                                GENERAL_DATA *general_data)

/*==========================================================================*/
   {/*Begin Routine*/
/*=======================================================================*/
/*         Local Variable declarations                                   */
  
  double vrecip,vself,vbgr,vecor,wght_now,vdummy;
  double press_inter;
  int i,iii,iver_get,irespa,ip;
  int pimd_on,ewald_on;

  /* Local pointers */
  int iget_full_intra    = class->energy_ctrl.iget_full_intra;
  int iget_res_intra     = class->energy_ctrl.iget_res_intra;
  int iget_full_inter    = class->energy_ctrl.iget_full_inter;
  int iget_res_inter     = class->energy_ctrl.iget_res_inter;
  int iget_pe_real_inter = class->energy_ctrl.iget_pe_real_inter;
  int iget_pv_real_inter = class->energy_ctrl.iget_pv_real_inter;
  int pi_beads           = class->clatoms_info.pi_beads;  
  int pi_beads_proc      = class->clatoms_info.pi_beads_proc;  
  int nchrg              = class->clatoms_info.nchrg;
  int np_tot             = class->communicate.np;
  int np_beads           = class->communicate.np_beads;
  int np_forc            = class->communicate.np_forc;  
  MPI_Comm comm_forc     = class->communicate.comm_forc;
  int myid               = class->communicate.myid;
  int myid_bead          = class->communicate.myid_bead;
  int pme_on             = class->part_mesh.pme_on;
  int pme_res_on         = class->part_mesh.pme_res_on;
  double alp_ewd         = general_data->ewald.alp_ewd;
  double alp_clus        = general_data->ewald.alp_clus;
  int nktot              = general_data->ewald.nktot;
  int nktot_res          = general_data->ewald.nktot_res;
  int *kastr             = general_data->ewald.kastr;
  int *kbstr             = general_data->ewald.kbstr;
  int *kcstr             = general_data->ewald.kcstr;
  int *kastr_res         = general_data->ewald.kastr_res;
  int *kbstr_res         = general_data->ewald.kbstr_res;
  int *kcstr_res         = general_data->ewald.kcstr_res;
  int *ibrk1             = general_data->ewald.ibrk1;
  int *ibrk2             = general_data->ewald.ibrk2;
  int *ibrk3             = general_data->ewald.ibrk3;
  int *ibrk1_res         = general_data->ewald.ibrk1_res;
  int *ibrk2_res         = general_data->ewald.ibrk2_res;
  double *clus_corr_r    = general_data->ewald.clus_corr_r;
  double *dclus_corr_r   = general_data->ewald.dclus_corr_r;
  int error_check_on     = general_data->error_check_on;
  double pext            = general_data->statepoint.pext;
  double vol             = general_data->cell.vol;
  int iperd              = general_data->cell.iperd;
  double *pvten_tot      = general_data->ptens.pvten_tot;
  double wght_ter_res    = class->for_scr.wght_ter_res;
  double wght_ter        = class->for_scr.wght_ter;
  double wght_tra_res    = bonded->intra_scr.wght_tra_res;
  int necor              = bonded->ecor.num;
  int np_send_self_bgr;
  double wght_diff;
  wght_diff              = wght_ter_res - wght_ter;

/*======================================================================*/
/* 0) Zero stuff and assign flags*/

  pimd_on  =  general_data->simopts.pimd
             +general_data->simopts.debug_pimd
             +general_data->simopts.cp_pimd
             +general_data->simopts.cp_wave_pimd
             +general_data->simopts.cp_wave_min_pimd
             +general_data->simopts.debug_cp_pimd;

  np_send_self_bgr = ((pimd_on == 1 && (np_beads + np_forc) == np_tot) 
                      ? np_tot : np_forc);

  iver_get = 0; if( (pimd_on==1) && (iget_full_inter==1) ){iver_get = 1;}
  ewald_on = 0; if( (nchrg>0) && (iperd>0) ){ewald_on=1;}

  vself       = 0.0;
  vbgr        = 0.0;   
  vecor       = 0.0;
  vrecip      = 0.0;

/*=================================================================*/
/* I) Self term calculation                                        */

  if(ewald_on == 1){
    if( (myid_bead==0) || (np_forc > 1) ){
      if( (iget_full_intra==1) || (iget_res_intra==0) ){ /* tra is correct*/
        ewald3d_selfbgr(&(class->clatoms_info),&general_data->ewald,
                        &(general_data->ptens),vol,wght_tra_res,
                        &vself,&vbgr,np_send_self_bgr,iget_pv_real_inter,iperd);
      }/* endif : intra respa */ 
    }/* endif : myid_bead = 0 */
  }/*endif:ewald_on*/

/*=================================================================*/
/* II) Reciprocal space calculation                                */

  if(ewald_on == 1){
   for(ip=1;ip<=pi_beads_proc;ip++){
  /*---------------------------------------------------------------*/
  /*  i) PME off                                                   */
   if(pme_on==0 && pme_res_on==0){
     if( (iget_res_inter == 0) && (iget_full_inter == 1) ){
#ifdef DEBUG
printf("ewald 1\n");
#endif
        ewald3d_recip(&(class->clatoms_info),&(class->clatoms_pos[ip]),
                   &(general_data->cell),&(general_data->ptens),alp_ewd,alp_clus,
                   nktot,kastr,kbstr,kcstr,ibrk1,ibrk2,
                   &(class->ewd_scr),&vrecip,wght_ter,iver_get,
                   &(class->class_comm_forc_pkg),iget_pv_real_inter, 
                   clus_corr_r,dclus_corr_r);
     }/*endif*/
     if( (iget_res_inter == 1) && (iget_full_inter == 0) && (nktot_res > 0)){
#ifdef DEBUG
printf("ewald 2\n");
#endif
        ewald3d_recip(&(class->clatoms_info),&(class->clatoms_pos[ip]),
                   &(general_data->cell),&(general_data->ptens),alp_ewd,alp_clus,
                   nktot_res,kastr_res,kbstr_res,kcstr_res,ibrk1_res,ibrk2_res,
                   &(class->ewd_scr),&vrecip,wght_ter_res,iver_get,
                   &(class->class_comm_forc_pkg),iget_pv_real_inter, 
                   clus_corr_r,dclus_corr_r);
     }/*endif*/
     if( (iget_res_inter == 1) && (iget_full_inter == 1) && (nktot_res > 0)){
#ifdef DEBUG
printf("ewald 3\n");
#endif
        ewald3d_recip_both(&(class->clatoms_info),&(class->clatoms_pos[ip]),
                   &(general_data->cell),&(general_data->ptens),alp_ewd,
                   nktot,kastr,kbstr,kcstr,ibrk1,ibrk2,ibrk3,
                   &(class->ewd_scr),&vrecip,wght_ter,wght_ter_res,
                   iver_get,&(class->class_comm_forc_pkg),iget_pv_real_inter, 
                   clus_corr_r,dclus_corr_r);
     }/*endif*/
     if( (iget_res_inter == 1) && (iget_full_inter == 1) && (nktot_res == 0)){
#ifdef DEBUG
printf("ewald 4\n");
#endif
        ewald3d_recip(&(class->clatoms_info),&(class->clatoms_pos[ip]),
                   &(general_data->cell),&(general_data->ptens),alp_ewd,alp_clus,
                   nktot,kastr,kbstr,kcstr,ibrk1,ibrk2,
                   &(class->ewd_scr),&vrecip,wght_ter,iver_get,
                   &(class->class_comm_forc_pkg),iget_pv_real_inter, 
                   clus_corr_r,dclus_corr_r);
     }/*endif*/
   }/*endif:standard ewald*/
  /*---------------------------------------------------------------*/
  /*  ii) PME on : PME_res off                                      */
   if(pme_on==1 && pme_res_on==0){
     if( (iget_res_inter == 1) && nktot_res > 0){
#ifdef DEBUG
printf("ewald 5\n");
#endif
        irespa = 0;
        ewald3d_recip_pme(&(class->clatoms_info),&(class->clatoms_pos[ip]),
                   &(general_data->cell),&(general_data->ptens),alp_ewd,
                   nktot,kastr,kbstr,kcstr,
                   &(class->ewd_scr),&vrecip,wght_ter_res,iver_get,
                   &(class->part_mesh),irespa,
                   &(class->class_comm_forc_pkg),iget_pv_real_inter,
                   &(general_data->pme_fft_pkg),&(class->for_scr), 
                   clus_corr_r,dclus_corr_r);
     }else{
      if( iget_full_inter == 1 ){
#ifdef DEBUG
printf("ewald 6\n");
#endif
        irespa = 0;
        ewald3d_recip_pme(&(class->clatoms_info),&(class->clatoms_pos[ip]),
                   &(general_data->cell),&(general_data->ptens),alp_ewd,
                   nktot,kastr,kbstr,kcstr,
                   &(class->ewd_scr),&vrecip,wght_ter,iver_get,
                   &(class->part_mesh),irespa,
                   &(class->class_comm_forc_pkg),iget_pv_real_inter,
                   &(general_data->pme_fft_pkg),&(class->for_scr), 
                   clus_corr_r,dclus_corr_r);
      }/*endif*/
     }/*endif*/
   }/*endif:pme:no_pme_res*/
  /*---------------------------------------------------------------*/
  /*  iii) PME on : PME_res on :                      */
   if(pme_on==1 && pme_res_on==1 ){
     if( (iget_res_inter == 1) && (iget_full_inter == 0) && (nktot_res > 0)){
#ifdef DEBUG
printf("ewald 7\n");
#endif
        irespa = 1;
        ewald3d_recip_pme(&(class->clatoms_info),&(class->clatoms_pos[ip]),
                   &(general_data->cell),&(general_data->ptens),alp_ewd,
                   nktot_res,kastr_res,kbstr_res,kcstr_res,
                   &(class->ewd_scr),&vrecip,wght_ter_res,iver_get,
                   &(class->part_mesh),irespa,
                   &(class->class_comm_forc_pkg),iget_pv_real_inter,
                   &(general_data->pme_res_fft_pkg),&(class->for_scr), 
                   clus_corr_r,dclus_corr_r);
     }/*endif*/
     if( (iget_res_inter == 1) && (iget_full_inter == 1) && (nktot_res > 0)){
#ifdef DEBUG
printf("ewald 8\n");
#endif
         irespa = 1;vdummy = 0.0;
         ewald3d_recip_pme(&(class->clatoms_info),&(class->clatoms_pos[ip]),
                   &(general_data->cell),&(general_data->ptens),alp_ewd,
                   nktot_res,kastr_res,kbstr_res,kcstr_res,
                   &(class->ewd_scr),&vdummy,wght_diff,iver_get,
                   &(class->part_mesh),irespa,
                   &(class->class_comm_forc_pkg),iget_pv_real_inter,
                   &(general_data->pme_res_fft_pkg),&(class->for_scr), 
                   clus_corr_r,dclus_corr_r);
         irespa = 0;
         ewald3d_recip_pme(&(class->clatoms_info),&(class->clatoms_pos[ip]),
                  &(general_data->cell),&(general_data->ptens),alp_ewd,
                  nktot,kastr,kbstr,kcstr,
                  &(class->ewd_scr),&vrecip,wght_ter,iver_get,
                  &(class->part_mesh),irespa,
                  &(class->class_comm_forc_pkg),iget_pv_real_inter,
                  &(general_data->pme_fft_pkg),&(class->for_scr), 
                  clus_corr_r,dclus_corr_r);
     }/*endif*/
     if( (iget_res_inter == 1) && (iget_full_inter == 1) && (nktot_res == 0)){
#ifdef DEBUG
printf("ewald 9\n");
#endif
         irespa = 0;
         ewald3d_recip_pme(&(class->clatoms_info),&(class->clatoms_pos[ip]),
                 &(general_data->cell),&(general_data->ptens),alp_ewd,
                 nktot,kastr,kbstr,kcstr,
                 &(class->ewd_scr),&vrecip,wght_ter,iver_get,
                 &(class->part_mesh),irespa,
                 &(class->class_comm_forc_pkg),iget_pv_real_inter,
                 &(general_data->pme_fft_pkg),&(class->for_scr), 
                 clus_corr_r,dclus_corr_r);
     }/*endif*/
   }/*endif:pme:pme_res*/
  /*---------------------------------------------------------------*/
   }/*endfor:ip*/
   vrecip /= (double) pi_beads;
  }/*endif:ewald_on*/

/*=================================================================*/
/* III) Ecorr calculation                                          */

  if( (ewald_on==1) && (necor > 0) ){
   for(ip=1;ip<=pi_beads_proc;ip++){
  /*---------------------------------------------------------------*/
  /*  i) nktot_res == 0 */
    if( (nktot_res == 0) && (iget_full_inter == 1) ){
#ifdef DEBUG
printf("ecor 1\n");
#endif
       ecor(&(class->clatoms_info),&(class->clatoms_pos[ip]),
            &(bonded->ecor),&(general_data->cell),&(bonded->intra_scr),
            &(general_data->ptens),&vecor,iver_get,
            &(class->class_comm_forc_pkg),
            iget_pe_real_inter,iget_pv_real_inter);
    }/*endif*/
  /*---------------------------------------------------------------*/
  /*  ii) nktot_res > 0 */
    if(nktot_res > 0 ){
     if( (iget_res_inter == 1) && (iget_full_inter == 0) ){
#ifdef DEBUG
printf("ecor 2\n");
#endif
       ecor(&(class->clatoms_info),&(class->clatoms_pos[ip]),
            &(bonded->ecor),&(general_data->cell),&(bonded->intra_scr),
            &(general_data->ptens),&vecor,iver_get,
            &(class->class_comm_forc_pkg),
            iget_pe_real_inter,iget_pv_real_inter);
     }/*endif*/
     if( (iget_res_inter == 1) && (iget_full_inter == 1) ){
#ifdef DEBUG
printf("ecor 3\n");
#endif
       ecor_both(&(class->clatoms_info),&(class->clatoms_pos[ip]),
                 &(bonded->ecor),&(general_data->cell),&(bonded->intra_scr),
                 &(general_data->ptens),&vecor,iver_get,
                 &(class->class_comm_forc_pkg),
                 iget_pe_real_inter,iget_pv_real_inter);
     }/*endif*/
    }/*endif*/
  /*---------------------------------------------------------------*/
   }/*endfor:pi_beads*/
   vecor /= (double) pi_beads;
  }/*endif:ewald_on*/

/*======================================================================*/
/* IV)  Get intermolecular contribution to the pressure                  */
  
  if(iget_full_inter==1){
    if(pimd_on==0){
     if(iget_pv_real_inter==1){
       general_data->stat_avg.press_inter = (pvten_tot[1]
                                           + pvten_tot[5]
                                           + pvten_tot[9])/(3.0*vol);
     }/*endif*/
    }else{
       get_vir_press(class,general_data,&press_inter);
       general_data->stat_avg.press_inter =  (press_inter)/(vol);
    }/*endif*/
  }/*endif*/

/*======================================================================*/
/* V) Collect and store PE                                           */
  
  if( (iget_full_intra==1) || (iget_full_inter==1) ){

    general_data->stat_avg.vrecip       = vrecip+vself+vbgr+vecor;

    if(iget_pe_real_inter==1){
      general_data->stat_avg.vintert    += (vrecip+vself+vbgr+vecor);
      general_data->stat_avg.vcoul      += (vrecip+vself+vbgr+vecor);
    }/*endif*/

  }/*endif*/

/*-----------------------------------------------------------------------*/
    }/*end routine */
/*==========================================================================*/





