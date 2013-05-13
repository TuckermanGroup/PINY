/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
/*                                                                          */
/*                         PI_MD:                                           */
/*             The future of simulation technology                          */
/*             ------------------------------------                         */
/*                     Module: force_npol                                   */
/*                                                                          */
/* This routine uses splines to calculate the real space intermolecular PE  */
/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

#include "standard_include.h"
#include "../typ_defs/typedefs_class.h"
#include "../typ_defs/typedefs_bnd.h"
#include "../typ_defs/typedefs_gen.h"
#include "../proto_defs/proto_real_space_local.h"
#include "../proto_defs/proto_energy_ctrl_entry.h"
#include "../proto_defs/proto_math.h"


/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void force_npol(CLATOMS_INFO *clatoms_info,CLATOMS_POS *clatoms_pos,
                FOR_SCR *for_scr,ATOMMAPS *atommaps,
                CELL *cell,PTENS *ptens,INTERACT *interact,
                ENERGY_CTRL *energy_ctrl,INTRA_SCR *intra_scr,double *vtert,
                double *vvdw,double *vcoul, int num_call,int lst_typ)

/*=====================================================================*/
/*           Begin Routine                                             */
      {/*Begin Routine*/
/*=====================================================================*/
/*           Local variable declarations                                */

  int i, iii;
  int imap = 0;
  int isum,icoul,ibrk;
  double wght_full, wght_res, wght_lnk;
  int intact_now,num_brk,intact_use;
  int energy_ctrl_iget_pe_real_inter;
  int energy_ctrl_iget_pv_real_inter;

/*========================================================================*/
/*             Store the Nasty Pointers for the dumb compiler             */

  int nchrg                       = clatoms_info->nchrg;
  int natm_typ                    = atommaps->natm_typ;
  double wght_ter                 = (for_scr->wght_ter);
  double wght_ter_res             = (for_scr->wght_ter_res); 
  int iperd                       = cell->iperd;
  int energy_ctrl_iget_full_inter = energy_ctrl->iget_full_inter;
  int energy_ctrl_iget_res_inter  = energy_ctrl->iget_res_inter;
  int energy_ctrl_isep_vvdw       = energy_ctrl->isep_vvdw;
  int energy_ctrl_int_res_ter     = energy_ctrl->int_res_ter;
  int energy_ctrl_itime           = energy_ctrl->itime;
  int energy_ctrl_iget_pe_real_inter_freq=energy_ctrl->iget_pe_real_inter_freq;
  int ishave_opt                  = interact->ishave_opt;
  int pi_beads                    = clatoms_info->pi_beads;
  int diele_opt                   = interact->dielectric_opt;
  double *xmod                    = clatoms_info->xmod;
  double *ymod                    = clatoms_info->ymod;
  double *zmod                    = clatoms_info->zmod;
  double *intra_scr_fx1           = intra_scr->fx1;
  double *intra_scr_fy1           = intra_scr->fy1;
  double *intra_scr_fz1           = intra_scr->fz1;
  double *intra_scr_fx2           = intra_scr->fx2;
  double *intra_scr_fy2           = intra_scr->fy2;
  double *intra_scr_fz2           = intra_scr->fz2;
  double *intra_scr_p11           = intra_scr->p11;
  double *intra_scr_p22           = intra_scr->p22;
  double *intra_scr_p33           = intra_scr->p33;
  double *intra_scr_p12           = intra_scr->p12;
  double *intra_scr_p21           = intra_scr->p21;
  double *intra_scr_p13           = intra_scr->p13;
  double *intra_scr_p31           = intra_scr->p31;
  double *intra_scr_p23           = intra_scr->p23;
  double *intra_scr_p32           = intra_scr->p32;
  double *clatoms_x               = clatoms_pos->x;
  double *clatoms_y               = clatoms_pos->y;
  double *clatoms_z               = clatoms_pos->z;
  double *clatoms_fx              = clatoms_pos->fx;
  double *clatoms_fy              = clatoms_pos->fy;
  double *clatoms_fz              = clatoms_pos->fz;
  double *clatoms_fxt             = clatoms_pos->fxt;
  double *clatoms_fyt             = clatoms_pos->fyt;
  double *clatoms_fzt             = clatoms_pos->fzt;
  double *clatoms_q               = clatoms_info->q;
  double *intra_scr_dx12          = intra_scr->dx12;
  double *intra_scr_dy12          = intra_scr->dy12;
  double *intra_scr_dz12          = intra_scr->dz12;
  double *intra_scr_dx56          = intra_scr->dx56;
  double *intra_scr_dy56          = intra_scr->dy56;
  double *intra_scr_dz56          = intra_scr->dz56;
  double *intra_scr_q2            = intra_scr->q2;
  double *interact_vcut_coul      = interact->vcut_coul;
  double *intra_scr_vcut_coul     = intra_scr->vcut_coul;
  double *interact_cutoff         = interact->cutoff;
  double *interact_cutoff_res     = interact->cutoff_res;
  double *intra_scr_cutoff        = intra_scr->cutoff;
  double *intra_scr_cutoff_res    = intra_scr->cut;
  double *interact_rmin_spl       = interact->rmin_spl;
  double *intra_scr_dr            = intra_scr->dr;
  double *intra_scr_dri           = intra_scr->dri;
  double *intra_scr_r             = intra_scr->r;
  double *intra_scr_del_r         = intra_scr->del_r;
  double *intra_scr_del_rc        = intra_scr->del_rc;
  double *intra_scr_swit          = intra_scr->swit;
  double *intra_scr_swit_hard     = intra_scr->dz43;
  double *intra_scr_vpot          = intra_scr->vpot;
  double *intra_scr_spl_tmp       = intra_scr->spl_tmp;
  double *intra_scr_c_0           = intra_scr->c_0;
  double *interact_cv0            = interact->cv0;
  double *interact_cv0_c          = interact->cv0_c;
  double *interact_cdv0           = interact->cdv0;
  double *interact_cdv0_c         = interact->cdv0_c;
  double *interact_dr_spl         = interact->dr_spl;
  double *interact_dri_spl        = interact->dri_spl;
  double *ptens_pvten             = ptens->pvten;
  double *ptens_pvten_tot         = ptens->pvten_tot;
  double *ptens_pvten_tmp         = ptens->pvten_tmp;
  int *for_scr_i_index            = for_scr->i_index;
  int *for_scr_j_index            = for_scr->j_index;
  int *for_scr_num_brk_i          = for_scr->num_brk_i;
  int *intra_scr_iatm_typ         = intra_scr->iatm_typ;
  int *intra_scr_jatm_typ         = intra_scr->jatm_typ;
  int *atommaps_iatm_atm_typ      = atommaps->iatm_atm_typ;
  int *for_scr_i_indext           = for_scr->i_indext;

  int *for_scr_i_indext2          = for_scr->i_indext2;
  int *inter_map_index            = interact->inter_map_index;
  int interact_nsplin             = interact->nsplin;
  int interact_nsplin_m2          = interact->nsplin-2;
  int hess_calc                   = clatoms_info->hess_calc;
  double interact_rheal_res       = interact->rheal_res;
  double interact_rheal_resi;
  interact_rheal_resi             = 1.0/interact->rheal_res;

/*========================================================================*/
/* Initialize variables */

  energy_ctrl_iget_pe_real_inter = 0;
  energy_ctrl_iget_pv_real_inter = 0;
  if(energy_ctrl_iget_full_inter==1){
   energy_ctrl_iget_pe_real_inter = 1;
   energy_ctrl_iget_pv_real_inter = 1;
  }/*endif*/
  if(energy_ctrl_itime > 1){
   if((energy_ctrl_itime % energy_ctrl_iget_pe_real_inter_freq)!=0){
    energy_ctrl_iget_pe_real_inter = 0;
    if(energy_ctrl_int_res_ter ==1){energy_ctrl_iget_pv_real_inter = 0;}
   }/*endif*/
  }/*endif*/
  energy_ctrl->iget_pe_real_inter = energy_ctrl_iget_pe_real_inter;
  energy_ctrl->iget_pv_real_inter = energy_ctrl_iget_pv_real_inter;

/*========================================================================*/
/* I) Local varibales and Error check                                     */

  intact_now = for_scr->intact;
  num_brk    = for_scr->num_brk;
  wght_lnk   = (for_scr->wght_lnk); 
  wght_full  = wght_ter*wght_lnk;
  wght_res   = wght_ter_res*wght_lnk;

  isum = 0;
  for (ibrk = 1; ibrk <= num_brk; ++ibrk) {
    isum += for_scr_num_brk_i[ibrk];
  }/*endfor*/
  if(isum != intact_now) {
    printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
    printf("Internal Error:\n ");
    printf("number of interactions inconsistent with\n");
    printf("the breakpoints: in force_npol\n");
    printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
    fflush(stdout);
    exit(1);
  }/*endif*/

/*========================================================================*/
/* II) Gather the displacements/charges/interaction types in this set     */

  npol_posdata_fetch(intra_scr_dx12,intra_scr_dy12,intra_scr_dz12,
                     for_scr_i_index,for_scr_j_index,for_scr_i_indext,
                     intra_scr_dx56,intra_scr_dy56,intra_scr_dz56,
                     clatoms_x,clatoms_y,clatoms_z,atommaps_iatm_atm_typ,
                     xmod,ymod,zmod,intra_scr_q2,clatoms_q,
                     num_brk,for_scr_num_brk_i,iperd,pi_beads,lst_typ,&icoul,
                     intact_now,nchrg,natm_typ);
     
/*========================================================================*/
/*  III) Distance Calculatation                                           */

  npol_dist_calc(intra_scr_dx12,intra_scr_dy12,intra_scr_dz12,
                 intra_scr_dx56,intra_scr_dy56,intra_scr_dz56,
                 intra_scr_r,iperd,pi_beads,intact_now,cell);

/*========================================================================*/
/* IV) Shave the interactions down using the cutoff                       */
  
   intact_use = intact_now;
   if(ishave_opt==1){
    npol_shave_skin(&intact_use,intra_scr_dx12,intra_scr_dy12,intra_scr_dz12,
                    intra_scr_r,intra_scr_q2,
                    for_scr_i_index,for_scr_j_index,
                    num_brk,for_scr_num_brk_i,
                    intra_scr_cutoff,interact_cutoff,
                    intra_scr_cutoff_res,interact_cutoff_res,
                    for_scr_i_indext,icoul,
                    energy_ctrl_iget_res_inter,
                    energy_ctrl_iget_full_inter,
                    energy_ctrl_int_res_ter);
   }/*endif*/

/*========================================================================*/
/* V) Gather ye rosebud spline parameters of i,j interaction while ye may */


  npol_spldata_fetch(intra_scr_r,
                     intra_scr_del_r,interact_rmin_spl,
                     interact_dri_spl,interact_dr_spl,
                     intra_scr_swit_hard,intra_scr_swit,
                     intra_scr_cutoff,intra_scr_cutoff_res,
                     interact_cutoff,interact_cutoff_res,
                     intra_scr_vcut_coul,
                     interact_vcut_coul,
                     for_scr_i_indext,
                     inter_map_index,
                     interact_rheal_res,interact_rheal_resi,
                     interact_nsplin,interact_nsplin_m2,
                     energy_ctrl_iget_res_inter,
                     iperd,icoul,diele_opt,intact_now,
                     ishave_opt);


/*========================================================================*/
/* VII) Get the PE using the spline parameters if necessary               */

  if ((energy_ctrl_iget_full_inter    == 1) && 
      (energy_ctrl_iget_pe_real_inter == 1)) {

   npol_pot_calc(intact_use,intra_scr_del_r,intra_scr_vpot,
                 for_scr_i_indext,
                 interact_cv0,
                 intra_scr_swit_hard,interact_cv0_c,
                 intra_scr_spl_tmp,
                 energy_ctrl_isep_vvdw,wght_lnk,
                 icoul,iperd,diele_opt, intra_scr_q2,
                 intra_scr_r,intra_scr_vcut_coul,
                 vvdw,vtert,vcoul);

  }/*endif:get the PE when needed */

/*=====================================================================*/
/* VIII)  Get the forces using the spline                              */

   npol_force_calc(intact_use,intra_scr_del_r,
                   intra_scr_spl_tmp,
                   for_scr_i_indext,
                   interact_cdv0,intra_scr_swit_hard,
                   intra_scr_fx1,intra_scr_dx12,
                   intra_scr_fy1,intra_scr_dy12,
                   intra_scr_fz1,intra_scr_dz12,
                   intra_scr_fx2,intra_scr_fy2,
                   intra_scr_fz2,
                   intra_scr_q2,intra_scr_r,
                   interact_cdv0_c,intra_scr_c_0,intra_scr_swit,
                   iperd,icoul,
                   diele_opt,wght_res,wght_full,
                   energy_ctrl_iget_full_inter,
                   energy_ctrl_iget_res_inter);

/*=====================================================================*/
/* IX)  If neeced, get the Hessian using the spline                               */

  if(hess_calc == 3){
     npol_hess_calc(clatoms_info,clatoms_pos,cell,intra_scr_cutoff);
  }/* endif */

/*========================================================================*/
/* X) Sum the forces                                                     */

   npol_force_reduc(num_brk,for_scr_num_brk_i,lst_typ,
                    clatoms_fx,intra_scr_fx2,
                    clatoms_fy,intra_scr_fy2,
                    clatoms_fz,intra_scr_fz2,
                    for_scr_i_index,for_scr_j_index);
       
/*========================================================================*/
/* XI) Construct and reduce the pressure tensors                           */

  if( iperd > 0 && iperd <= 3){
   npol_press_reduc(intact_use,iperd,ptens_pvten,
           ptens_pvten_tmp,
           intra_scr_p11,intra_scr_p22,intra_scr_p33,
           intra_scr_p12,intra_scr_p21,intra_scr_p13,
           intra_scr_p31,intra_scr_p23,intra_scr_p32,
           intra_scr_fx1,intra_scr_fy1,intra_scr_fz1,
           intra_scr_fx2,intra_scr_fy2,intra_scr_fz2,
           intra_scr_dx12,intra_scr_dy12,
           intra_scr_dz12,ptens_pvten_tot,
           energy_ctrl_iget_full_inter,
           energy_ctrl_int_res_ter,energy_ctrl_iget_pv_real_inter,
           wght_lnk,wght_ter);
  }/*endif*/

/*========================================================================*/
/* XII) Construct and reduce the virial estimator                          */

  if( (energy_ctrl_iget_full_inter   ==1)&&(pi_beads>1)&&
      (energy_ctrl_iget_pe_real_inter==1) ){

   npol_vir_reduc(energy_ctrl_int_res_ter,intact_use,
                  intra_scr_fx1,intra_scr_fy1,
                  intra_scr_fz1,wght_lnk,
                  num_brk,for_scr_num_brk_i,
                  clatoms_fxt,clatoms_fyt,
                  clatoms_fzt,for_scr_i_index,
                  for_scr_j_index);
       
  }/*endif*/
  
/*--------------------------------------------------------------------------*/
   }/* end routine force_npol */
/*==========================================================================*/





/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void npol_vspl_fetch(int n,double del[], double spl_out[],
                    int i_index[],double c0_data[],double swit[])

/*=======================================================================*/
/*              Begin Routine                                            */
    {/*Begin Routine*/
/*=======================================================================*/
/*            Local variable declarations                                */

  int i,ktemp0,ktemp1,ktemp2,ktemp3,iii;
  double del_tmp,c0,c1,c2,c3;
  double f0,fp1,fm1,fm2;
  static double oneth = (1.0/3.0);
  static double onesi = (1.0/6.0);
  /*========================================================================*/
  /* I) Fetch the coefs */

  for(i=1;i<=n;i++){
    ktemp0 = i_index[i];
    ktemp1 = ktemp0 + 1;
    ktemp2 = ktemp0 - 1;
    ktemp3 = ktemp0 - 2;
    c0  = c0_data[ktemp0];
    fp1 = c0_data[ktemp1];
    fm1 = c0_data[ktemp2];
    fm2 = c0_data[ktemp3];
    c1 = oneth*fp1+0.5*c0-fm1+onesi*fm2;
    c2 = 0.5*fp1-c0+0.5*fm1;
    c3 = onesi*fp1-0.5*c0+0.5*fm1-onesi*fm2;
    del_tmp = del[i];
    spl_out[i]  =  swit[i]*(c0 + del_tmp*(c1 + del_tmp* 
                                     (c2 + del_tmp*c3)));

  }/*endfor*/

/*--------------------------------------------------------------------------*/
}/* end routine vspl_fetch */
/*==========================================================================*/





/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void pten_sum(int n,int iperd, double pvten[],double pvten_tmp[],
              double p11[],  double p22[],  double p33[],
              double p12[],  double p21[],  double p13[],
              double p31[],  double p23[],  double p32[],
              double fx[],  double fy[],  double fz[],
              double dx12[], double dy12[], double dz12[])

/*=======================================================================*/
/*              Begin Routine                                            */
    {/*Begin Routine*/
/*=======================================================================*/
/*            Local variable declarations                                */

  int i;

/*=======================================================================*/
/* I) construct and sum the pressure tensor */
  
 for(i=1;i <= 9; ++i) {
  pvten_tmp[i] = 0.0;
 }
  if(iperd == 2 || iperd == 3) {
    for(i=1;i <= n; ++i) {
      p11[i] = dx12[i]*fx[i];
      p22[i] = dy12[i]*fy[i];
      p33[i] = dz12[i]*fz[i];
      p12[i] = dx12[i]*fy[i];
    }/*endfor*/
    for(i=1;i <= n; ++i) {
      pvten_tmp[1] += p11[i];
      pvten_tmp[5] += p22[i];
      pvten_tmp[9] += p33[i];
      pvten_tmp[2] += p12[i];
    }/*endfor*/
  }/*endif*/
  if(iperd == 3) {
    for(i=1;i <= n; ++i) {
      p13[i] = dx12[i]*fz[i];
      p23[i] = dy12[i]*fz[i];
    }/*endfor*/
    for(i=1;i <= n; ++i) {
      pvten_tmp[3] += p13[i];
      pvten_tmp[6] += p23[i];
    }/*endfor*/
  }/*endif*/
 pvten_tmp[4] = pvten_tmp[2];     
 pvten_tmp[7] = pvten_tmp[3];     
 pvten_tmp[8] = pvten_tmp[6];     
 for(i=1;i <= 9; ++i) {
  pvten[i] += pvten_tmp[i];
 }

/*--------------------------------------------------------------------------*/
}/* end routine */
/*==========================================================================*/




/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void pten_sum_roll(int n,int iperd, double pvten[],double pvten_tmp[],
                   double p11[],  double p22[],  double p33[],
                   double p12[],  double p21[],  double p13[],
                   double p23[], 
                   double fx[],  double fy[],  double fz[],
                   double dx12[], double dy12[], double dz12[])

/*=======================================================================*/
/*              Begin Routine                                            */
    {/*Begin Routine*/
/*=======================================================================*/
/*            Local variable declarations                                */

  int i,n2;

/*=======================================================================*/
/* I) construct and sum the pressure tensor */
  
 for(i=1;i <= 9; ++i) {
  pvten_tmp[i] = 0.0;
 }
  if(iperd == 2 || iperd == 3) {
    for(i=1;i <= n; i++) {
      p11[i] = dx12[i]*fx[i];
      p22[i] = dy12[i]*fy[i];
      p33[i] = dz12[i]*fz[i];
      p12[i] = dx12[i]*fy[i];
    }/*endfor*/
    n2 = n-4;
    for(i=1;i <= n2; i+=5) {
      pvten_tmp[1] += (p11[i]+p11[(i+1)]+p11[(i+2)]+p11[(i+3)]+p11[(i+4)]);
      pvten_tmp[5] += (p22[i]+p22[(i+1)]+p22[(i+2)]+p22[(i+3)]+p22[(i+4)]);
      pvten_tmp[9] += (p33[i]+p33[(i+1)]+p33[(i+2)]+p33[(i+3)]+p33[(i+4)]);
      pvten_tmp[2] += (p12[i]+p12[(i+1)]+p12[(i+2)]+p12[(i+3)]+p12[(i+4)]);
    }/*endfor*/
    n2 = i;
    for(i=n2;i <= n; i++) {
      pvten_tmp[1] += (p11[i]);
      pvten_tmp[5] += (p22[i]);
      pvten_tmp[9] += (p33[i]);
      pvten_tmp[2] += (p12[i]);
    }
  }/*endif*/
  if(iperd == 3) {
    for(i=1;i <= n; ++i) {
      p13[i] = dx12[i]*fz[i];
      p23[i] = dy12[i]*fz[i];
    }/*endfor*/
    n2 = n-4;
    for(i=1;i <= n2; i+=5) {
      pvten_tmp[3] += (p13[i]+p13[(i+1)]+p13[(i+2)]+p13[(i+3)]+p13[(i+4)]);
      pvten_tmp[6] += (p23[i]+p23[(i+1)]+p23[(i+2)]+p23[(i+3)]+p23[(i+4)]);
    }/*endfor*/
    n2 = i;
    for(i=n2;i <= n; i++) {
      pvten_tmp[3] += p13[i];
      pvten_tmp[6] += p23[i];
    }
  }/*endif*/
 pvten_tmp[4] = pvten_tmp[2];     
 pvten_tmp[7] = pvten_tmp[3];     
 pvten_tmp[8] = pvten_tmp[6];     
 for(i=1;i <= 9; ++i) {
  pvten[i] += pvten_tmp[i];
 }

/*--------------------------------------------------------------------------*/
}/* end routine */
/*==========================================================================*/




/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void npol_posdata_fetch(double *dx12,double *dy12,double *dz12,
                        int *i_index,int *j_index,int *jatm_typ,
                        double *dx56,double *dy56,double *dz56,
                        double *clatoms_x,double *clatoms_y,double *clatoms_z,
                        int *clatoms_atm_typ,
                        double *xmod,double *ymod,double *zmod,
                        double *q2,double *clatoms_q,
                        int num_brk,int *num_brk_i,int iperd,int pi_beads,
                        int lst_typ,int *icoul_ret,int intact_now,int nchrg,
                        int natm_typ)

/*=======================================================================*/
/*              Begin Routine                                            */
    {/*Begin Routine*/
/*=======================================================================*/
/*            Local variable declarations                                */

  int i,ibrk,lower;
  int jpart,ktemp,icoul;
  int upper,n2;
  double qsum,qtemp;
  double x1t,y1t,z1t;
  double x2t,y2t,z2t;
  double x5t,y5t,z5t;
  double x6t,y6t,z6t;
  double q1t,q2t;
  int natm1_typ2,i2,i3,jjtemp,iitemp,kktemp;
  int iatm_typt,jatm_typt;


/*=======================================================================*/
/* I) Fetch the positions, get displacements and interaction number      */

   natm1_typ2 = (natm_typ+1) * 2;
  /*---------------------------------------------------------------*/
  /* For a link list */
if(lst_typ == 3){
    for (jpart = 1; jpart <= intact_now; jpart++) {
      x1t         = clatoms_x[i_index[jpart]];
      y1t         = clatoms_y[i_index[jpart]];
      z1t         = clatoms_z[i_index[jpart]];
      x2t         = clatoms_x[j_index[jpart]];
      y2t         = clatoms_y[j_index[jpart]];
      z2t         = clatoms_z[j_index[jpart]];
      dx12[jpart] = x1t-x2t;
      dy12[jpart] = y1t-y2t;
      dz12[jpart] = z1t-z2t;
      jatm_typt   = clatoms_atm_typ[j_index[jpart]];
      iatm_typt   = clatoms_atm_typ[i_index[jpart]];
      i2          = iatm_typt;
      i3          = jatm_typt;
      jjtemp      = MIN(i2,i3); 
      iitemp      = i2+i3-jjtemp;  /*MAX(i2,i3)*/
      kktemp      = ((jjtemp - 1) * (natm1_typ2 - jjtemp))/2; 
      jatm_typ[jpart] = (kktemp + iitemp - jjtemp + 1); 
    }/*endfor*/
    if(iperd > 0 && pi_beads > 1){     
      for (jpart = 1; jpart <= intact_now; jpart++) {
        x5t         = xmod[i_index[jpart]];
        y5t         = ymod[i_index[jpart]];
        z5t         = zmod[i_index[jpart]];
        x6t         = xmod[j_index[jpart]];
        y6t         = ymod[j_index[jpart]];
        z6t         = zmod[j_index[jpart]];
        dx56[jpart] = x5t-x6t;
        dy56[jpart] = y5t-y6t;
        dz56[jpart] = z5t-z6t;
      }/*endfor*/
    }/*endif*/
}/*endif*/

  /*---------------------------------------------------------------*/
  /* For a nolist or verlet list */

if(lst_typ != 3){
   lower = 1;
   for (ibrk = 1; ibrk <= num_brk; ibrk++) {
     upper = num_brk_i[ibrk] + lower - 1;
     ktemp = i_index[lower];
#ifndef NO_PRAGMA
#pragma ivdep
#endif
     for(jpart=lower; jpart<=upper;jpart++){
      x1t         = clatoms_x[ktemp];
      y1t         = clatoms_y[ktemp];
      z1t         = clatoms_z[ktemp];
      x2t         = clatoms_x[j_index[jpart]];
      y2t         = clatoms_y[j_index[jpart]];
      z2t         = clatoms_z[j_index[jpart]];
      dx12[jpart] = x1t-x2t;
      dy12[jpart] = y1t-y2t;
      dz12[jpart] = z1t-z2t;
      jatm_typt   = clatoms_atm_typ[j_index[jpart]];
      iatm_typt   = clatoms_atm_typ[ktemp];
      i2          = iatm_typt;
      i3          = jatm_typt;
      jjtemp      = MIN(i2,i3); 
      iitemp      = i2+i3-jjtemp;  /*MAX(i2,i3)*/
      kktemp      = ((jjtemp - 1) * (natm1_typ2 - jjtemp))/2; 
      jatm_typ[jpart] = (kktemp + iitemp - jjtemp + 1); 
     }/*endfor*/
     if(iperd > 0 && pi_beads>1){     
       for(jpart=lower; jpart<=upper;jpart++){
        x5t         = xmod[ktemp];
        y5t         = ymod[ktemp];
        z5t         = zmod[ktemp];
        x6t         = xmod[j_index[jpart]];
        y6t         = ymod[j_index[jpart]];
        z6t         = zmod[j_index[jpart]];
        dx56[jpart] = x5t-x6t;
        dy56[jpart] = y5t-y6t;
        dz56[jpart] = z5t-z6t;
       }/*endfor*/
     }/*endif*/
     lower += num_brk_i[ibrk];
    }/*endfor*/
}/*endif lst_typ */ 


/*=======================================================================*/
/*  Find out if charged atoms are around and store q1*q2                 */

 icoul = 0;
 if(nchrg>0){
     icoul = 1;
     for (jpart = 1; jpart <= intact_now; jpart++) {
       q1t         = clatoms_q[i_index[jpart]];
       q2t         = clatoms_q[j_index[jpart]];
       q2[jpart]  = q1t*q2t;
     }/*endfor*/
     n2 = intact_now-4;
     qsum = 0.0;
     for (jpart = 1; jpart <= n2; jpart+=5) {
       qtemp = (q2[jpart] * q2[jpart] + 
                q2[(jpart+1)] * q2[(jpart+1)] +
                q2[(jpart+2)] * q2[(jpart+2)] +
                q2[(jpart+3)] * q2[(jpart+3)] +
                q2[(jpart+4)] * q2[(jpart+4)]);
       qsum += qtemp;
     }/*endfor*/
     n2 = jpart;
     for (jpart = n2; jpart <= intact_now; jpart++) {
      qsum += (q2[jpart]*q2[jpart]);
     }/*endfor*/

     if(qsum == 0.0) {icoul = 0;}
 }/*endif*/
 *icoul_ret = icoul;

/*--------------------------------------------------------------------------*/
    }/* end routine */
/*==========================================================================*/



/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void npol_spldata_fetch(double *intra_scr_r,
                        double *intra_scr_del_r,double *interact_rmin_spl,
                        double *interact_dri_spl,double *interact_dr_spl,
                        double *intra_scr_swit_hard,double *intra_scr_swit,
                        double *intra_scr_cutoff,double *intra_scr_cutoff_res,
                        double *interact_cutoff,double *interact_cutoff_res,
                        double *intra_scr_vcut_coul,
                        double *interact_vcut_coul,
                        int *for_scr_i_indext,
                        int *inter_map_index,
                        double interact_rheal_res,double interact_rheal_resi,
                        int interact_nsplin,int interact_nsplin_m2,
                        int energy_ctrl_iget_res_inter,
                        int iperd,int icoul,int diele_opt,int intact_now,
                        int ishave_opt)


/*=======================================================================*/
/*              Begin Routine                                            */
    {/*Begin Routine*/
/*=======================================================================*/
/*            Local variable declarations                                */

  int jpart,ktemp,iii;
  int kktemp,jjtemp;
  int for_scr_i_indext_jpart;
  double rmin_now, dr,dri,r;
  double dr_cut,rsp;
  double br;

/*=======================================================================*/
/* I) Get the spline stuff                                              */


 if((iperd == 0) && (icoul != 0) && (diele_opt==0) ) {
    for (jpart = 1; jpart <= intact_now; jpart++) {
      for_scr_i_indext_jpart = for_scr_i_indext[jpart];
      ktemp =  inter_map_index[(for_scr_i_indext_jpart)];
      intra_scr_vcut_coul[jpart] = interact_vcut_coul[ktemp];
    }/*endfor*/
 }/*endif*/


switch (ishave_opt){

 case 0:
  if(energy_ctrl_iget_res_inter == 0) {
    for (jpart = 1; jpart <= intact_now; jpart++) {
      for_scr_i_indext_jpart = for_scr_i_indext[jpart];
      ktemp =  inter_map_index[(for_scr_i_indext_jpart)];
      rmin_now         = interact_rmin_spl[ktemp];
      dr               = interact_dr_spl[ktemp];
      dri              = interact_dri_spl[ktemp];

      intra_scr_cutoff[jpart] = interact_cutoff[for_scr_i_indext_jpart];

      r                = intra_scr_r[jpart];
      kktemp           = (int) (((r - rmin_now)*dri)+0.5) + 3;
      kktemp           = MIN(kktemp,interact_nsplin_m2);
      kktemp           = MAX(kktemp,3);
      rsp              =((double)(kktemp-3))*dr+rmin_now;
      intra_scr_del_r[jpart] = (r - rsp)*dri;
      jjtemp  =  kktemp + (inter_map_index[for_scr_i_indext_jpart]-1)*interact_nsplin ;

      for_scr_i_indext[jpart] = jjtemp;
      dr_cut           = (intra_scr_cutoff[jpart]-r);
      intra_scr_swit_hard[jpart] = (dr_cut >= 0.0 ? 1.0:0.0);
    }/*endfor*/
  }else{

    for (jpart = 1; jpart <= intact_now; jpart++) {
      for_scr_i_indext_jpart = for_scr_i_indext[jpart];
      ktemp =  inter_map_index[(for_scr_i_indext_jpart)];
      rmin_now         = interact_rmin_spl[ktemp];
      dr               = interact_dr_spl[ktemp];
      dri              = interact_dri_spl[ktemp];

      intra_scr_cutoff[jpart]    = interact_cutoff[for_scr_i_indext_jpart];
      intra_scr_cutoff_res[jpart]= interact_cutoff_res[for_scr_i_indext_jpart];

      r                = intra_scr_r[jpart];
      kktemp           = (int) (((r - rmin_now)*dri)+0.5) + 3;
      kktemp           = MIN(kktemp,interact_nsplin_m2);
      kktemp           = MAX(kktemp,3);
      rsp              =((double)(kktemp-3))*dr+rmin_now;
      intra_scr_del_r[jpart] = (r - rsp)*dri;
      jjtemp  = kktemp + (inter_map_index[for_scr_i_indext_jpart]-1)*interact_nsplin;

      for_scr_i_indext[jpart] = jjtemp;
      dr_cut           = (intra_scr_cutoff[jpart]-r);
      intra_scr_swit_hard[jpart] = (dr_cut >= 0.0 ? 1.0:0.0);
    }/*endfor*/
  }/*endif*/
 break;

 case 1:
    for (jpart = 1; jpart <= intact_now; jpart++) {
      for_scr_i_indext_jpart = for_scr_i_indext[jpart];
      ktemp = inter_map_index[(for_scr_i_indext_jpart)];
      rmin_now         = interact_rmin_spl[ktemp];
      dr               = interact_dr_spl[ktemp];
      dri              = interact_dri_spl[ktemp];
      r                = intra_scr_r[jpart];
      kktemp           = (int) (((r - rmin_now)*dri)+0.5) + 3;
      kktemp           = MIN(kktemp,interact_nsplin_m2);
      kktemp           = MAX(kktemp,3);
      rsp              =((double)(kktemp-3))*dr+rmin_now;
      intra_scr_del_r[jpart] = (r - rsp)*dri;
      jjtemp  = kktemp + (inter_map_index[for_scr_i_indext_jpart]-1)*interact_nsplin;
      for_scr_i_indext[jpart] = jjtemp;
      dr_cut           = (intra_scr_cutoff[jpart]-r);
      intra_scr_swit_hard[jpart] = (dr_cut >= 0.0 ? 1.0:0.0);
    }/*endfor*/
     
 break;

}/*end switch*/


/*========================================================================*/
/* II) Calculate the RESPA switching function                              */

  if(energy_ctrl_iget_res_inter == 1) {
    for (jpart = 1; jpart <= intact_now; jpart++) {
      r = intra_scr_r[jpart];
      br = (r - intra_scr_cutoff_res[jpart] 
            + interact_rheal_res)*interact_rheal_resi;
      br = MAX(br,0.0);
      br = MIN(br,1.0);
      intra_scr_swit[jpart] = br * br * (br * 2.0 - 3.0) + 1.0;
    }/*endfor*/
  }/*endif*/

/*--------------------------------------------------------------------------*/
    }/* end routine */
/*==========================================================================*/




/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void npol_dist_calc(double *dx12,double *dy12,double *dz12,
                    double *dx56,double *dy56,double *dz56,
                    double *r,
                    int iperd,int pi_beads,int intact_now,CELL *cell)

/*=======================================================================*/
/*              Begin Routine                                            */
    {/*Begin Routine*/
/*=======================================================================*/
/*            Local variable declarations                                */

int jpart;
double r2;

  /*========================================================================*/
  /* VIII) Impose boundary conditions                                       */

  if(iperd>0){
#define DEBUG_OFF
#ifdef DEBUG
    period(intact_now,dx12,dy12,dz12,cell);
#endif
#ifdef DEBUG_OFF
   if(pi_beads==1){     
    period(intact_now,dx12,dy12,dz12,cell);
   }else{
    period_pimd(intact_now,dx12,dy12,dz12,dx56,dy56,dz56,cell);
   }/*endif*/
#endif
  }/*endif*/

  /*========================================================================*/
  /* VIII) Distance                                                         */

  for (jpart = 1; jpart <= intact_now; jpart++) {
    r2 = (dx12[jpart]*dx12[jpart]+ dy12[jpart]*dy12[jpart]
        + dz12[jpart]*dz12[jpart]);
    r[jpart] = sqrt(r2);
  }/*endfor*/

/*--------------------------------------------------------------------------*/
    }/* end routine */
/*==========================================================================*/





/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void npol_pot_calc(int intact_now,double *intra_scr_del_r,
                   double *intra_scr_vpot,
                   int *for_scr_i_indext,
                   double *interact_cv0,double *intra_scr_swit_hard,
                   double *interact_cv0_c,double *intra_scr_spl_tmp,
                   int energy_ctrl_isep_vvdw,double wght_lnk,
                   int icoul,int iperd,int diele_opt, double *intra_scr_q2,
                   double *intra_scr_r,double *intra_scr_vcut_coul,
                   double *vvdw,double *vtert,double *vcoul)

/*=======================================================================*/
/*              Begin Routine                                            */
    {/*Begin Routine*/
/*=======================================================================*/
/*            Local variable declarations                                */
    double vvdw_temp,vtemp;
    int jpart;
    int iii;
/*=======================================================================*/
/*  i) Vanderwaals                                                       */

    npol_vspl_fetch(intact_now,intra_scr_del_r,intra_scr_vpot,
                   for_scr_i_indext,interact_cv0,intra_scr_swit_hard);
    if(energy_ctrl_isep_vvdw == 1){
      vvdw_temp = dsum1_npol(intact_now,intra_scr_vpot)*wght_lnk;
      *(vvdw) += vvdw_temp;
    }/*endif:sep_vdw*/
    
/*=======================================================================*/
/*  ii) Coulomb                                                          */

    if((icoul != 0)) {
      if(iperd == 0 && (diele_opt==0)) {
        for (jpart = 1; jpart <= intact_now; jpart++) {
            intra_scr_vpot[jpart] += 
              intra_scr_q2[jpart]*intra_scr_swit_hard[jpart]*
              (1.0/intra_scr_r[jpart]-intra_scr_vcut_coul[jpart]);
        }/*endfor*/
      }/*endif*/

      if(iperd > 0 || (diele_opt==1)) {
        npol_vspl_fetch(intact_now,intra_scr_del_r,intra_scr_spl_tmp,
                       for_scr_i_indext,interact_cv0_c,intra_scr_swit_hard);
        for (jpart = 1; jpart <= intact_now; jpart++) {
          vtemp = intra_scr_q2[jpart]*intra_scr_spl_tmp[jpart];
          intra_scr_vpot[jpart] += vtemp;
        }/*endfor*/
      }/*endif*/
    }/*endif:coulomb interactions*/

/*=======================================================================*/
/*  iii) Sum the PE                                                      */

    vtemp = dsum1_npol(intact_now,intra_scr_vpot)*wght_lnk;
    *vtert += vtemp;
    if(energy_ctrl_isep_vvdw == 1){ *vcoul += vtemp - vvdw_temp;}

/*--------------------------------------------------------------------------*/
    }/* end routine */
/*==========================================================================*/





/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void npol_force_calc(int intact_now,double *intra_scr_del_r,
                        double *intra_scr_spl_tmp,
                        int *for_scr_i_indext,
                        double *interact_cdv0,
                        double *intra_scr_swit_hard,
                        double *intra_scr_fx1,double *intra_scr_dx12,
                        double *intra_scr_fy1,double *intra_scr_dy12,
                        double *intra_scr_fz1,double *intra_scr_dz12,
                        double *intra_scr_fx2,double *intra_scr_fy2,
                        double *intra_scr_fz2,
                        double *intra_scr_q2,double *intra_scr_r,
                        double *interact_cdv0_c,double *intra_scr_spl_tmp2,
                        double *intra_scr_swit,
                        int iperd,int icoul,
                        int diele_opt,double wght_res,double wght_full,
                        int energy_ctrl_iget_full_inter,
                        int energy_ctrl_iget_res_inter)

/*=======================================================================*/
/*              Begin Routine                                            */
    {/*Begin Routine*/
/*=======================================================================*/
/*            Local variable declarations                                */
    
    int jpart,iii;
    double fmult,dvtemp,intra_scr_r3;

/*=======================================================================*/
/* I) Vanderwaals spline                                                 */

  npol_vspl_fetch(intact_now,intra_scr_del_r,intra_scr_spl_tmp,
             for_scr_i_indext,interact_cdv0,intra_scr_swit_hard);

/*=======================================================================*/
/* I) Vanderwaals only */

if(icoul==0){

  for (jpart = 1; jpart <= intact_now; jpart++) {
    dvtemp = intra_scr_spl_tmp[jpart];
    intra_scr_fx1[jpart] = intra_scr_dx12[jpart]*dvtemp;
    intra_scr_fy1[jpart] = intra_scr_dy12[jpart]*dvtemp;
    intra_scr_fz1[jpart] = intra_scr_dz12[jpart]*dvtemp;
  }/*endfor*/

}/*endfor*/

/*=======================================================================*/
/*  II) Vanderwalls plus  Coulomb                                       */

if((icoul != 0)) {
    if(iperd == 0 && diele_opt==0) {
      for (jpart = 1; jpart <= intact_now; jpart++) {

        intra_scr_r3 = intra_scr_r[jpart]*
                       intra_scr_r[jpart]*intra_scr_r[jpart];
        dvtemp=(intra_scr_q2[jpart]/intra_scr_r3*
               intra_scr_swit_hard[jpart] + intra_scr_spl_tmp[jpart]);
        intra_scr_fx1[jpart] = intra_scr_dx12[jpart]*dvtemp;
        intra_scr_fy1[jpart] = intra_scr_dy12[jpart]*dvtemp;
        intra_scr_fz1[jpart] = intra_scr_dz12[jpart]*dvtemp;

      }/*endfor*/
    }/*endif*/
    if(iperd > 0 || diele_opt==1) {
      npol_vspl_fetch(intact_now,intra_scr_del_r,intra_scr_spl_tmp2,
                 for_scr_i_indext,interact_cdv0_c,intra_scr_swit_hard);
      for (jpart = 1; jpart <= intact_now; jpart++) {
        dvtemp = (intra_scr_q2[jpart]*intra_scr_spl_tmp2[jpart]
                 +intra_scr_spl_tmp[jpart]);
        intra_scr_fx1[jpart] = intra_scr_dx12[jpart]*dvtemp;
        intra_scr_fy1[jpart] = intra_scr_dy12[jpart]*dvtemp;
        intra_scr_fz1[jpart] = intra_scr_dz12[jpart]*dvtemp;
      }/*endfor*/
    }/*endif*/
}/*endif:coulomb interactions*/

/*========================================================================*/
/* III) RESPA forces                                                      */

  /*------------------------------------------------------------------------*/
  /* i) Inter Respa only */

  if((energy_ctrl_iget_full_inter==0)&&(energy_ctrl_iget_res_inter==1)){
    for (jpart = 1; jpart <= intact_now; jpart++) {
      fmult = intra_scr_swit[jpart]*wght_res;
      intra_scr_fx2[jpart] = intra_scr_fx1[jpart]*fmult;
      intra_scr_fy2[jpart] = intra_scr_fy1[jpart]*fmult;
      intra_scr_fz2[jpart] = intra_scr_fz1[jpart]*fmult;
    }/*endfor*/
  }/*endif*/
  /*------------------------------------------------------------------------*/
  /* ii) Inter Respa + correction */
  if((energy_ctrl_iget_full_inter==1)&&(energy_ctrl_iget_res_inter==1)){
    for (jpart = 1; jpart <= intact_now; jpart++) {
      fmult = (intra_scr_swit[jpart]*wght_res  
               +  (1.0-intra_scr_swit[jpart])*wght_full);
      intra_scr_fx2[jpart] = intra_scr_fx1[jpart]*fmult;
      intra_scr_fy2[jpart] = intra_scr_fy1[jpart]*fmult;
      intra_scr_fz2[jpart] = intra_scr_fz1[jpart]*fmult;
    }/*endfor*/
  }/*endif*/
  /*----------------------------------------------------------------------*/
  /* iii) Full only */
  if((energy_ctrl_iget_full_inter==1)&&(energy_ctrl_iget_res_inter==0)){
    for (jpart = 1; jpart <= intact_now; jpart++) {
      intra_scr_fx2[jpart] = intra_scr_fx1[jpart]*wght_full;
      intra_scr_fy2[jpart] = intra_scr_fy1[jpart]*wght_full;
      intra_scr_fz2[jpart] = intra_scr_fz1[jpart]*wght_full;
    }/*endfor*/
  }/*endif*/
  
/*--------------------------------------------------------------------------*/
    }/* end routine */
/*==========================================================================*/




/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void npol_force_reduc(int num_brk,int *for_scr_num_brk_i,
                      int lst_typ,
                      double *clatoms_fx,double *intra_scr_fx2,
                      double *clatoms_fy,double *intra_scr_fy2,
                      double *clatoms_fz,double *intra_scr_fz2,
                      int *for_scr_i_index,int *for_scr_j_index)
       

/*=======================================================================*/
/*              Begin Routine                                            */
    {/*Begin Routine*/
/*=======================================================================*/
/*            Local variable declarations                                */

   int lower,ibrk,upper,jpart,ktemp,n2;

/*=======================================================================*/


   lower = 1;
   for (ibrk = 1; ibrk <= num_brk; ibrk++) {
     upper = for_scr_num_brk_i[ibrk] + lower - 1;
  /*------------------------------------------------------------------------*/
  /*  i) ith particle */
    if(lst_typ == 3){ 
#ifndef NO_PRAGMA
#pragma ivdep
#endif
    for(jpart=lower; jpart<=upper;jpart++){
      ktemp = for_scr_i_index[jpart];
      clatoms_fx[ktemp] += intra_scr_fx2[jpart];
      clatoms_fy[ktemp] += intra_scr_fy2[jpart];
      clatoms_fz[ktemp] += intra_scr_fz2[jpart];
     }/*endfor*/
    }else{
     ktemp = for_scr_i_index[lower];
     n2 = upper-4;
     for(jpart=lower; jpart<=n2;jpart+=5){
       clatoms_fx[ktemp] += (intra_scr_fx2[jpart]
                            +intra_scr_fx2[(jpart+1)]
                            +intra_scr_fx2[(jpart+2)]
                            +intra_scr_fx2[(jpart+3)]
                            +intra_scr_fx2[(jpart+4)]);
       clatoms_fy[ktemp] += (intra_scr_fy2[jpart]
                            +intra_scr_fy2[(jpart+1)]
                            +intra_scr_fy2[(jpart+2)]
                            +intra_scr_fy2[(jpart+3)]
                            +intra_scr_fy2[(jpart+4)]);
       clatoms_fz[ktemp] += (intra_scr_fz2[jpart]
                            +intra_scr_fz2[(jpart+1)]
                            +intra_scr_fz2[(jpart+2)]
                            +intra_scr_fz2[(jpart+3)]
                            +intra_scr_fz2[(jpart+4)]);
     }/*endfor*/
     n2 = jpart;
     for(jpart=n2; jpart<=upper;jpart++){
       clatoms_fx[ktemp] += intra_scr_fx2[jpart];
       clatoms_fy[ktemp] += intra_scr_fy2[jpart];
       clatoms_fz[ktemp] += intra_scr_fz2[jpart];
     }/*endfor*/

    }/*endif lst_typ */
  /*------------------------------------------------------------------------*/
  /*  ii) jth particle */
#ifndef NO_PRAGMA
#pragma ivdep
#endif
    for(jpart=lower; jpart<=upper;jpart++){
      ktemp = for_scr_j_index[jpart];
      clatoms_fx[ktemp] -= intra_scr_fx2[jpart];
      clatoms_fy[ktemp] -= intra_scr_fy2[jpart];
      clatoms_fz[ktemp] -= intra_scr_fz2[jpart];
    }/*endfor*/
    lower += for_scr_num_brk_i[ibrk];
  }/*endfor*/
  

/*--------------------------------------------------------------------------*/
    }/* end routine */
/*==========================================================================*/




/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void npol_press_reduc(int intact_now,int iperd,double *ptens_pvten,
           double *ptens_pvten_tmp,
           double *intra_scr_p11,double *intra_scr_p22,double *intra_scr_p33,
           double *intra_scr_p12,double *intra_scr_p21,double *intra_scr_p13,
           double *intra_scr_p31,double *intra_scr_p23,double *intra_scr_p32,
           double *intra_scr_fx1,double *intra_scr_fy1,double *intra_scr_fz1,
           double *intra_scr_fx2,double *intra_scr_fy2,double *intra_scr_fz2,
           double *intra_scr_dx12,double *intra_scr_dy12,
           double *intra_scr_dz12,double *ptens_pvten_tot,
           int energy_ctrl_iget_full_inter,int energy_ctrl_int_res_ter,
           int energy_ctrl_iget_pv_real_inter,
           double wght_lnk,double wght_ter)
       
/*=======================================================================*/
/*              Begin Routine                                            */
    {/*Begin Routine*/
/*=======================================================================*/
/*            Local variable declarations                                */

    int jpart,i;

/*=======================================================================*/

  pten_sum_roll(intact_now,iperd,ptens_pvten,ptens_pvten_tmp,
           intra_scr_p11,intra_scr_p22,intra_scr_p33,
           intra_scr_p12,intra_scr_p21,intra_scr_p13,
           intra_scr_p23,
           intra_scr_fx2,intra_scr_fy2,intra_scr_fz2,
           intra_scr_dx12,intra_scr_dy12,intra_scr_dz12);

  if((energy_ctrl_iget_full_inter   ==1)&&(energy_ctrl_int_res_ter==1)&&
     (energy_ctrl_iget_pv_real_inter==1) ){
    for (jpart = 1; jpart <= intact_now; jpart++) {
      intra_scr_fx1[jpart] *= wght_lnk;
      intra_scr_fy1[jpart] *= wght_lnk;
      intra_scr_fz1[jpart] *= wght_lnk;
    }/*endfor*/
    pten_sum(intact_now,iperd,ptens_pvten_tot,ptens_pvten_tmp,
             intra_scr_p11,intra_scr_p22,intra_scr_p33,
             intra_scr_p12,intra_scr_p21,intra_scr_p13,
             intra_scr_p31,intra_scr_p23,intra_scr_p32,
             intra_scr_fx1,intra_scr_fy1,intra_scr_fz1,
             intra_scr_dx12,intra_scr_dy12,intra_scr_dz12);
  }/*endif*/

  if((energy_ctrl_iget_full_inter==1)&&(energy_ctrl_int_res_ter == 0)){
   for(i=1;i<=9;++i){
     ptens_pvten_tot[i] += (ptens_pvten_tmp[i]*wght_lnk/wght_ter);
   }/*endfor*/
  }/*endif*/

/*--------------------------------------------------------------------------*/
    }/* end routine */
/*==========================================================================*/



/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void npol_vir_reduc(int energy_ctrl_int_res_ter,int intact_now,
                    double *intra_scr_fx1,double *intra_scr_fy1,
                    double *intra_scr_fz1,double wght_lnk,
                    int num_brk,int *for_scr_num_brk_i,
                    double *clatoms_fxt,double *clatoms_fyt,
                    double *clatoms_fzt,int *for_scr_i_index,
                    int *for_scr_j_index)

/*=======================================================================*/
/*              Begin Routine                                            */
    {/*Begin Routine*/
/*=======================================================================*/
/*            Local variable declarations                                */

    int ibrk,upper,lower,jpart,ktemp;

/*=======================================================================*/

    if((energy_ctrl_int_res_ter == 0)){
     for (jpart = 1; jpart <= intact_now; jpart++) {
      intra_scr_fx1[jpart] *= wght_lnk;
      intra_scr_fy1[jpart] *= wght_lnk;
      intra_scr_fz1[jpart] *= wght_lnk;
    }/*endfor*/
   }/*endif*/

/*=======================================================================*/

   lower = 1;
   for (ibrk = 1; ibrk <= num_brk; ibrk++) {
    upper = for_scr_num_brk_i[ibrk] + lower - 1;
  /*------------------------------------------------------------------------*/
  /*  i) ith particle */

#ifndef NO_PRAGMA
#pragma ivdep
#endif
    for(jpart=lower;jpart<=upper;jpart++){
     ktemp = for_scr_i_index[jpart];
     clatoms_fxt[ktemp] += intra_scr_fx1[jpart];
     clatoms_fyt[ktemp] += intra_scr_fy1[jpart];
     clatoms_fzt[ktemp] += intra_scr_fz1[jpart];
    }/*endfor*/

  /*------------------------------------------------------------------------*/
  /*  ii) jth particle */
#ifndef NO_PRAGMA
#pragma ivdep
#endif
    for(jpart=lower; jpart<=upper;jpart++){
      ktemp = for_scr_j_index[jpart];
      clatoms_fxt[ktemp] -= intra_scr_fx1[jpart];
      clatoms_fyt[ktemp] -= intra_scr_fy1[jpart];
      clatoms_fzt[ktemp] -= intra_scr_fz1[jpart];
     }/*endfor*/
    lower += for_scr_num_brk_i[ibrk];
   }/*endfor*/


/*--------------------------------------------------------------------------*/
    }/* end routine */
/*==========================================================================*/



/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void npol_shave_skin(int *intact_use,double *intra_scr_dx12,
                    double *intra_scr_dy12,double *intra_scr_dz12,
                    double *intra_scr_r,double *intra_scr_q2,
                    int *for_scr_i_index,int *for_scr_j_index,
                    int num_brk,int *for_scr_num_brk_i,
                    double *intra_scr_cutoff,double *interact_cutoff,
                    double *intra_scr_cutoff_res,double *interact_cutoff_res,
                    int *for_scr_i_indext,int icoul,
                    int energy_ctrl_iget_res_inter,
                    int energy_ctrl_iget_full_inter,
                    int energy_ctrl_int_res_ter)


/*=======================================================================*/
/*              Begin Routine                                            */
    {/*Begin Routine*/
/*=======================================================================*/
/*            Local variable declarations                                */

    int ibrk,upper,lower,jpart,ktemp,ic_old,ic_tot;
    int intact_now = *intact_use;

/*=======================================================================*/
/* I) Fetch the appropriate cutoffs */
 
  if(energy_ctrl_iget_full_inter==1){
    for (jpart = 1; jpart <= intact_now; jpart++) {
      ktemp = for_scr_i_indext[jpart];
      intra_scr_cutoff[jpart]  = interact_cutoff[ktemp]; 
    }/*endfor*/
  }else{
    for (jpart = 1; jpart <= intact_now; jpart++) {
      ktemp = for_scr_i_indext[jpart];
      intra_scr_cutoff[jpart]  = interact_cutoff_res[ktemp]; 
    }/*endfor*/
  }/*endif*/

/*=======================================================================*/
/* II) Shave                                                             */

   lower = 1;
   ic_tot = 0;
   for (ibrk = 1; ibrk <= num_brk; ibrk++) {
    upper  = for_scr_num_brk_i[ibrk] + lower - 1;
    ic_old = ic_tot;    
    if(icoul==0){
     for(jpart=lower;jpart<=upper;jpart++){
      if(intra_scr_r[jpart]<intra_scr_cutoff[jpart]){
        ic_tot++;
        intra_scr_dx12[ic_tot]   = intra_scr_dx12[jpart];
        intra_scr_dy12[ic_tot]   = intra_scr_dy12[jpart];
        intra_scr_dz12[ic_tot]   = intra_scr_dz12[jpart];
        intra_scr_r[ic_tot]      = intra_scr_r[jpart];
        intra_scr_cutoff[ic_tot] = intra_scr_cutoff[jpart];
        for_scr_i_index[ic_tot]  = for_scr_i_index[jpart];
        for_scr_j_index[ic_tot]  = for_scr_j_index[jpart];
        for_scr_i_indext[ic_tot] = for_scr_i_indext[jpart];
      }/*endif*/
     }/*endfor*/
    }else{
     for(jpart=lower;jpart<=upper;jpart++){
      if(intra_scr_r[jpart]<intra_scr_cutoff[jpart]){
        ic_tot++;
        intra_scr_dx12[ic_tot]   = intra_scr_dx12[jpart];
        intra_scr_dy12[ic_tot]   = intra_scr_dy12[jpart];
        intra_scr_dz12[ic_tot]   = intra_scr_dz12[jpart];
        intra_scr_r[ic_tot]      = intra_scr_r[jpart];
        intra_scr_q2[ic_tot]     = intra_scr_q2[jpart];
        intra_scr_cutoff[ic_tot] = intra_scr_cutoff[jpart];
        for_scr_i_index[ic_tot]  = for_scr_i_index[jpart];
        for_scr_j_index[ic_tot]  = for_scr_j_index[jpart];
        for_scr_i_indext[ic_tot] = for_scr_i_indext[jpart];
      }/*endif*/
     }/*endfor*/
    }/*endif*/
    lower += for_scr_num_brk_i[ibrk];
    for_scr_num_brk_i[ibrk] = ic_tot-ic_old;
   }/*endfor*/
   *intact_use = ic_tot;

/*=======================================================================*/
/* III) Assign the respa cutoff if necessary                             */

  if(energy_ctrl_int_res_ter==1&&energy_ctrl_iget_full_inter==0){
    for (jpart = 1; jpart <= ic_tot; jpart++) {
     intra_scr_cutoff_res[jpart]  = intra_scr_cutoff[jpart]; 
    }/*endfor*/
  }/*endif*/

  if(energy_ctrl_int_res_ter==1&&energy_ctrl_iget_full_inter==1){
    for (jpart = 1; jpart <= ic_tot; jpart++) {
      ktemp = for_scr_i_indext[jpart];
      intra_scr_cutoff_res[jpart]  = interact_cutoff_res[ktemp]; 
    }/*endfor*/
  }/*endif*/


/*--------------------------------------------------------------------------*/
    }/* end routine */
/*==========================================================================*/





/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
  double dsum1_npol(int n,double *a)
/*==========================================================================*/
{
  int i,j,n2;
  double dsum1;
/*==========================================================================*/

  dsum1 = 0.;
  n2 = n-4;
  for(i=1; i<=n2 ;i+=5){
    dsum1 += (a[i]+a[(i+1)]+a[(i+2)]+a[(i+3)]+a[(i+4)]);
  }
  n2 = i;
  for(i=n2; i<=n ;i++){
    dsum1 += a[i];
  }
  return dsum1;
/*==========================================================================*/
 }
/*==========================================================================*/



/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
void npol_hess_calc(CLATOMS_INFO *clatoms_info,CLATOMS_POS *clatoms_pos,CELL *cell,
                    double *rcut)
/*==========================================================================*/
{/* begin function */

/*--------------------------------------------------------------------------*/
/*  Local variables */

  int i,j,k,l;
  int hess_ind;
  double dx,dy,dz,r,r2,r3,r4,r5;
  double arg,rrpi;
  double alp3;
  double f1,f2,eee;
  int iperd = cell->iperd;
  int npart = clatoms_info->natm_tot;
  double alp_ewd = clatoms_info->alp_ewd;
  double *hxx = clatoms_pos->hess_xx;
  double *hxy = clatoms_pos->hess_xy;
  double *hxz = clatoms_pos->hess_xz;
  double *hyy = clatoms_pos->hess_yy;
  double *hyz = clatoms_pos->hess_yz;
  double *hzz = clatoms_pos->hess_zz;
  double *x = clatoms_pos->x;
  double *y = clatoms_pos->y;
  double *z = clatoms_pos->z;
  double *q = clatoms_info->q;

/*--------------------------------------------------------------------------*/
/*  Dumb calculation first */

  rrpi = 1.0/sqrt(M_PI);
  alp3 = alp_ewd*alp_ewd*alp_ewd;

/*--------------------------------------------------------------------------*/
/*  Periodicity 0 */

  if(iperd == 0){

      for(i=1;i<=npart;i++){
       for(j=1;j<=npart;j++){
         hess_ind = (i-1)*npart + j;
         if(i!=j){
           dx = x[i]-x[j];
           dy = y[i]-y[j];
           dz = z[i]-z[j];
           r = sqrt(dx*dx + dy*dy + dz*dz);
           r2 = r*r;
           r3 = r*r*r;
           hxx[hess_ind] += (q[i]*q[j]/r3)*
                                  (1.0-3.0*dx*dx/r2);
           hyy[hess_ind] += (q[i]*q[j]/r3)*
                                  (1.0-3.0*dy*dy/r2);
           hzz[hess_ind] += (q[i]*q[j]/r3)*
                                 (1.0-3.0*dz*dz/r2);

           hxy[hess_ind] -= 3.0*(q[i]*q[j]/r3)*dx*dy/r2;
           hxz[hess_ind] -= 3.0*(q[i]*q[j]/r3)*dx*dz/r2;
           hyz[hess_ind] -= 3.0*(q[i]*q[j]/r3)*dy*dz/r2;
         } else {
	  for(k=1;k<=npart;k++){
	    if(k != j){
             dx = x[k]-x[j];
             dy = y[k]-y[j];
             dz = z[k]-z[j];
             r = sqrt(dx*dx + dy*dy + dz*dz);
             r2 = r*r;
             r3 = r*r*r;
             hxx[hess_ind] -= (q[k]*q[j]/r3)*
                                    (1.0-3.0*dx*dx/r2);
             hyy[hess_ind] -= (q[k]*q[j]/r3)*
                                    (1.0-3.0*dy*dy/r2);
             hzz[hess_ind] -= (q[k]*q[j]/r3)*
                                    (1.0-3.0*dz*dz/r2);
             hxy[hess_ind] += 3.0*(q[k]*q[j]/r3)*dx*dy/r2;
             hxz[hess_ind] += 3.0*(q[k]*q[j]/r3)*dx*dz/r2;
             hyz[hess_ind] += 3.0*(q[k]*q[j]/r3)*dy*dz/r2;
	    }/* endif */
          }/* endfor */
        }/* endif */
      }/* endfor */
    }/* endfor */

  }/*endif iperd */

/*--------------------------------------------------------------------------*/
/*  Periodicity > 0 */

  if(iperd != 0){

      for(i=1;i<=npart;i++){
       for(j=1;j<=npart;j++){
         hess_ind = (i-1)*npart + j;
         if(i!=j){
           dx = x[i]-x[j];
           dy = y[i]-y[j];
           dz = z[i]-z[j];
           period(1,&dx,&dy,&dz,cell);
           r = sqrt(dx*dx + dy*dy + dz*dz);
           if(r<rcut[j]){
             arg = r*alp_ewd;
             eee = exp(-arg*arg);
             r2 = r*r;
             r3 = r2*r;
             r4 = r3*r;
             r5 = r3*r2;
             f1 = gerfc(arg)/r3 + 2.0*alp_ewd*rrpi*eee/r2;
             f2 = 3.0*gerfc(arg)/r5 + 6.0*alp_ewd*rrpi*eee/r4 + 4.0*alp3*rrpi*eee/r2;
             hxx[hess_ind] += q[i]*q[j]*
                                    (f1-f2*dx*dx);
             hyy[hess_ind] += q[i]*q[j]*
                                    (f1-f2*dy*dy);
             hzz[hess_ind] += q[i]*q[j]*
                                    (f1-f2*dz*dz);

             hxy[hess_ind] -= f2*q[i]*q[j]*dx*dy;
             hxz[hess_ind] -= f2*q[i]*q[j]*dx*dz;
             hyz[hess_ind] -= f2*q[i]*q[j]*dy*dz;
           }/* endif */
         } else {
	  for(k=1;k<=npart;k++){
	    if(k != j){
             dx = x[k]-x[j];
             dy = y[k]-y[j];
             dz = z[k]-z[j];
             period(1,&dx,&dy,&dz,cell);
             r = sqrt(dx*dx + dy*dy + dz*dz);
             if(r<rcut[j]){
               arg = r*alp_ewd;
               eee = exp(-arg*arg);
               r2 = r*r;
               r3 = r2*r;
               r4 = r3*r;
               r5 = r3*r2;
               f1 = gerfc(arg)/r3 + 2.0*alp_ewd*rrpi*eee/r2;
               f2 = 3.0*gerfc(arg)/r5 + 6.0*alp_ewd*rrpi*eee/r4 + 4.0*alp3*rrpi*eee/r2;
               hxx[hess_ind] -= q[k]*q[j]*
                                      (f1-f2*dx*dx);
               hyy[hess_ind] -= q[k]*q[j]*
                                      (f1-f2*dy*dy);
               hzz[hess_ind] -= q[k]*q[j]*
                                      (f1-f2*dz*dz);
               hxy[hess_ind] += f2*q[k]*q[j]*dx*dy;
               hxz[hess_ind] += f2*q[k]*q[j]*dx*dz;
               hyz[hess_ind] += f2*q[k]*q[j]*dy*dz;
             }/* endif */
	    }/* endif */
          }/* endfor */
        }/* endif */
      }/* endfor */
    }/* endfor */

  }/* endif iperd */

/*==========================================================================*/
}/* end function */
/*==========================================================================*/










