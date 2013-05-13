/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
/*                                                                          */
/*                         PI_MD:                                           */
/*             The future of simulation technology                          */
/*             ------------------------------------                         */
/*                     Module: ecorr.c                                      */
/*                                                                          */
/* This routine computes the energy and forces from                         */ 
/* the intramolecular ewald corrections.                                    */
/*                                                                          */
/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

#include "standard_include.h"
#include "../typ_defs/typedefs_class.h"
#include "../typ_defs/typedefs_bnd.h"
#include "../typ_defs/typedefs_gen.h"
#include "../proto_defs/proto_intra_entry.h"
#include "../proto_defs/proto_energy_ctrl_entry.h"

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void ecor(CLATOMS_INFO *clatoms_info,CLATOMS_POS *clatoms_pos,
          ECOR *ecor,CELL *cell, 
          INTRA_SCR *intra_scr,PTENS *ptens, double *vecort, int iver_get,
          CLASS_COMM_FORC_PKG *class_comm_forc_pkg,
          int iget_pe_real_inter,int iget_pv_real_inter)

/*==========================================================================*/
/*               Begin subprogram:                                          */
 {/*begin routine*/
/*==========================================================================*/
/*               Local variable declarations                                */
  int iecor,ibig,ibig1,ioff,iend,nnow;     /* Indices and counters   */
  int kktemp;

  double rsp,dvtemp;
  double r122,r12;
  int i,ktemp,ieoff,mtemp,iii;
  int ktemp1,ktemp2,ktemp5,ktemp6,iboff;

/* Define local pointers                                                 */
  int *intra_scr_i_index  = intra_scr->iatm_typ;
  double *intra_scr_q12   = intra_scr->q1;
  double *intra_scr_fx1   = intra_scr->fx1;
  double *intra_scr_fy1   = intra_scr->fy1;
  double *intra_scr_fz1   = intra_scr->fz1;
  double *intra_scr_dx12  = intra_scr->dx12;
  double *intra_scr_dy12  = intra_scr->dy12;
  double *intra_scr_dz12  = intra_scr->dz12;
  double *intra_scr_dx56  = intra_scr->dx56;
  double *intra_scr_dy56  = intra_scr->dy56;
  double *intra_scr_dz56  = intra_scr->dz56;
  double *intra_scr_spl_tmp  = intra_scr->vpot;
  double *intra_scr_p11   = intra_scr->p11;
  double *intra_scr_p22   = intra_scr->p22;
  double *intra_scr_p33   = intra_scr->p33;
  double *intra_scr_p21   = intra_scr->p21;
  double *intra_scr_p12   = intra_scr->p12;
  double *intra_scr_p13   = intra_scr->p13;
  double *intra_scr_p23   = intra_scr->p23;
  double *intra_scr_del_r = intra_scr->del_r;
  double *clatoms_x       = clatoms_pos->x;
  double *clatoms_y       = clatoms_pos->y;
  double *clatoms_z       = clatoms_pos->z;
  double *clatoms_q       = clatoms_info->q;
  double *clatoms_fx      = clatoms_pos->fx;
  double *clatoms_fy      = clatoms_pos->fy;
  double *clatoms_fz      = clatoms_pos->fz;
  double *clatoms_fxt     = clatoms_pos->fxt;
  double *clatoms_fyt     = clatoms_pos->fyt;
  double *clatoms_fzt     = clatoms_pos->fzt;
  double *ptens_pvten_tmp = ptens->pvten_tmp;
  double *xmod            = clatoms_info->xmod;
  double *ymod            = clatoms_info->ymod;
  double *zmod            = clatoms_info->zmod;
  int *ecor_j1            = ecor->j1;
  int *ecor_j2            = ecor->j2;
  int nsplin_m2           = ecor->nsplin_m2;
  int np_forc             = class_comm_forc_pkg->num_proc;
  int myid_forc           = class_comm_forc_pkg->myid;
  double dr_spl           = ecor->dr_spl;
  double dri_spl          = ecor->dri_spl;
  double rmin_spl         = ecor->rmin_spl;
  int iperd               = cell->iperd;
  int intra_perds         = cell->intra_perds;
  int pi_beads            = clatoms_info->pi_beads;
  int nlen_use            = intra_scr->nlen;
  int ntot                = ecor->num;
  double *pvten           = ptens->pvten;
  double *pvten_tot       = ptens->pvten_tot;
  double *ecor_cv0        = ecor->cv0;
  double *ecor_cdv0       = ecor->cdv0;
  int nktot_res           = ecor->nktot_res;
  double wfor             = intra_scr->wght_ter;


  if(nktot_res>0){
    ecor_cdv0       = ecor->cdv0_res;
    wfor            = intra_scr->wght_ter_res;
  }/*endif*/

/*==========================================================================*/
/* I) loop over all the ecors in steps of nlen to save memory               */

  for(i=1;i<=9;i++){(ptens_pvten_tmp)[i]=0;}


  for(ibig=1;ibig <= ntot;ibig += nlen_use) {

/*--------------------------------------------------------------------------*/
/*  A) Offsets to save some memory by doing only nlen ecors at a time       */

    ibig1 = ibig-1;
    ioff = -ibig1;
    iend = MIN(ntot,ibig1+nlen_use);
    nnow = iend-ibig1;

/*--------------------------------------------------------------------------*/
/*  B) Gather displacements                                                 */

    for(ieoff=1,iecor=ibig;iecor <= iend; ++iecor,++ieoff) {       
       ktemp1 = ecor_j1[iecor];
       ktemp2 = ecor_j2[iecor];
       intra_scr_dx12[ieoff] = clatoms_x[ktemp1]-clatoms_x[ktemp2];
       intra_scr_dy12[ieoff] = clatoms_y[ktemp1]-clatoms_y[ktemp2];
       intra_scr_dz12[ieoff] = clatoms_z[ktemp1]-clatoms_z[ktemp2];
    }/*endfor*/

    if(intra_perds==1&&pi_beads>1){     
       for(ieoff=1,iecor=ibig;iecor <= iend; ++iecor,++ieoff) {       
         ktemp5 = ecor_j1[iecor];
         ktemp6 = ecor_j2[iecor];
         intra_scr_dx56[ieoff] = xmod[ktemp5]-xmod[ktemp6];
         intra_scr_dy56[ieoff] = ymod[ktemp5]-ymod[ktemp6];
         intra_scr_dz56[ieoff] = zmod[ktemp5]-zmod[ktemp6];
       }/*endfor*/
    }/*endif*/

/*--------------------------------------------------------------------------*/
/*  C) Periodic boundary conditions                                         */

    if(intra_perds == 1) {
      if(pi_beads==1){     
        period(nnow,intra_scr_dx12,intra_scr_dy12,intra_scr_dz12,cell);
      }else{
        period_pimd(nnow,intra_scr_dx12,intra_scr_dy12,intra_scr_dz12,
                    intra_scr_dx56,intra_scr_dy56,intra_scr_dz56,cell);
      }/*endif*/
    }/*endif*/

/*--------------------------------------------------------------------------*/
/*  D) Spline info                                                          */

    for(iecor=1,iboff=ibig;iecor <= nnow; ++iecor,++iboff) {
       r122         = intra_scr_dx12[iecor]*intra_scr_dx12[iecor]
                    + intra_scr_dy12[iecor]*intra_scr_dy12[iecor]
                    + intra_scr_dz12[iecor]*intra_scr_dz12[iecor];
       r12          = sqrt(r122);
       kktemp       = (int) (((r12 - rmin_spl)*dri_spl)+0.5) + 3;
       kktemp       = MIN(kktemp,nsplin_m2);
       kktemp       = MAX(kktemp,3);
       rsp          = ((double)(kktemp-3))*dr_spl+rmin_spl;
       intra_scr_del_r[iecor]   = (r12 - rsp)*dri_spl;
       intra_scr_i_index[iecor] = kktemp;
       intra_scr_q12[iecor]     = 
                       clatoms_q[ecor_j1[iboff]]*clatoms_q[ecor_j2[iboff]];
    }/*endfor iecor*/

/*--------------------------------------------------------------------------*/
/*  E) Sum the potential energy                                             */

    if(iget_pe_real_inter==1){
     ecor_vspl_fetch(nnow,intra_scr_del_r,intra_scr_spl_tmp,
                     intra_scr_i_index,ecor_cv0);
     for(iecor=1;iecor<=nnow;++iecor)
           {*vecort += (intra_scr_spl_tmp)[iecor]*intra_scr_q12[iecor];}
    }/*endif*/

/*--------------------------------------------------------------------------*/
/*  F) Get the forces                                                       */

    ecor_vspl_fetch(nnow,intra_scr_del_r,intra_scr_spl_tmp,
                    intra_scr_i_index,ecor_cdv0);
    for(iecor = 1; iecor <= nnow; iecor++) {
       dvtemp = intra_scr_spl_tmp[iecor]*intra_scr_q12[iecor];
       intra_scr_fx1[iecor] = intra_scr_dx12[iecor]*dvtemp;
       intra_scr_fy1[iecor] = intra_scr_dy12[iecor]*dvtemp;
       intra_scr_fz1[iecor] = intra_scr_dz12[iecor]*dvtemp;
    }/*endfor*/

/*--------------------------------------------------------------------------*/
/*  G) Pressure tensor, etc.                                                */

    ecor_pvten_sum_roll(nnow,intra_scr_dx12,intra_scr_dy12,intra_scr_dz12,
                       intra_scr_fx1,intra_scr_fy1,intra_scr_fz1,
                       intra_scr_p11,intra_scr_p22,intra_scr_p33,intra_scr_p12,
                       intra_scr_p13,intra_scr_p23,ptens_pvten_tmp,iperd);

/*--------------------------------------------------------------------------*/
/*  H) Scatter the forces                                                   */

    if(iver_get==1){

      for(iecor=ibig;iecor <= iend; ++iecor) {
         ktemp = ecor_j1[iecor];
         mtemp = iecor+ioff;
         clatoms_fxt[ktemp] += intra_scr_fx1[mtemp];
         clatoms_fyt[ktemp] += intra_scr_fy1[mtemp];
         clatoms_fzt[ktemp] += intra_scr_fz1[mtemp];
       }/*endfor*/

       for(iecor=ibig;iecor <= iend; ++iecor) {
         ktemp = ecor_j2[iecor];
         mtemp = iecor+ioff;
         clatoms_fxt[ktemp] -= intra_scr_fx1[mtemp];
         clatoms_fyt[ktemp] -= intra_scr_fy1[mtemp];
         clatoms_fzt[ktemp] -= intra_scr_fz1[mtemp];
       }/*endfor*/

    }/*endif*/

    for(ieoff=1;ieoff <= nnow; ++ieoff) {
       (intra_scr_fx1)[ieoff] *= wfor;
       (intra_scr_fy1)[ieoff] *= wfor;
       (intra_scr_fz1)[ieoff] *= wfor;
    }/*endfor*/

    for(iecor=ibig;iecor <= iend; ++iecor) {
       ktemp = ecor_j1[iecor];
       mtemp = iecor+ioff;
       clatoms_fx[ktemp] += intra_scr_fx1[mtemp];
       clatoms_fy[ktemp] += intra_scr_fy1[mtemp];
       clatoms_fz[ktemp] += intra_scr_fz1[mtemp];
    }/*endfor*/

    for(iecor=ibig;iecor <= iend; ++iecor) {
       ktemp = ecor_j2[iecor];
       mtemp = iecor+ioff;
       clatoms_fx[ktemp] -= intra_scr_fx1[mtemp];
       clatoms_fy[ktemp] -= intra_scr_fy1[mtemp];
       clatoms_fz[ktemp] -= intra_scr_fz1[mtemp];
    }/*endfor*/

  }/* endfor ibig */

/*==========================================================================*/
/* Increment the Pressure tensor */

  ptens_pvten_tmp[4] = ptens_pvten_tmp[2];
  ptens_pvten_tmp[7] = ptens_pvten_tmp[3];
  ptens_pvten_tmp[8] = ptens_pvten_tmp[6];
  for(i=1;i<=9;i++){pvten[i]     += ptens_pvten_tmp[i]*wfor;}

  if(iget_pv_real_inter==1){    
    for(i=1;i<=9;i++){pvten_tot[i] += ptens_pvten_tmp[i];}
  }/*endif*/

/*--------------------------------------------------------------------------*/
  }/* end routine */
/*==========================================================================*/





/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void ecor_both(CLATOMS_INFO *clatoms_info,CLATOMS_POS *clatoms_pos,
               ECOR *ecor,CELL *cell, 
               INTRA_SCR *intra_scr,PTENS *ptens, double *vecort, int iver_get,
               CLASS_COMM_FORC_PKG *class_comm_forc_pkg,
               int iget_pe_real_inter,int iget_pv_real_inter)

/*==========================================================================*/
/*               Begin subprogram:                                          */
    {/*begin routine*/
/*==========================================================================*/
/*               Local variable declarations                                */
  int iecor,ibig,ibig1,ioff,iend,nnow;     /* Indices and counters   */
  int kktemp;

  double rsp,dvtemp;
  double r122,r12;
  int i,ktemp,ieoff,mtemp,iii;
  int ktemp1,ktemp2,ktemp5,ktemp6,iboff;

/* Define local pointers                                                 */
  int *intra_scr_i_index  = intra_scr->iatm_typ;
  double *intra_scr_q12   = intra_scr->q1;
  double *intra_scr_fx1   = intra_scr->fx1;
  double *intra_scr_fy1   = intra_scr->fy1;
  double *intra_scr_fz1   = intra_scr->fz1;
  double *intra_scr_fx2   = intra_scr->fx2;
  double *intra_scr_fy2   = intra_scr->fy2;
  double *intra_scr_fz2   = intra_scr->fz2;
  double *intra_scr_dx12  = intra_scr->dx12;
  double *intra_scr_dy12  = intra_scr->dy12;
  double *intra_scr_dz12  = intra_scr->dz12;
  double *intra_scr_dx56  = intra_scr->dx56;
  double *intra_scr_dy56  = intra_scr->dy56;
  double *intra_scr_dz56  = intra_scr->dz56;
  double *intra_scr_spl_tmp  = intra_scr->vpot;
  double *intra_scr_p11   = intra_scr->p11;
  double *intra_scr_p22   = intra_scr->p22;
  double *intra_scr_p33   = intra_scr->p33;
  double *intra_scr_p21   = intra_scr->p21;
  double *intra_scr_p12   = intra_scr->p12;
  double *intra_scr_p13   = intra_scr->p13;
  double *intra_scr_p23   = intra_scr->p23;
  double *intra_scr_del_r = intra_scr->del_r;
  double *clatoms_x       = clatoms_pos->x;
  double *clatoms_y       = clatoms_pos->y;
  double *clatoms_z       = clatoms_pos->z;
  double *clatoms_q       = clatoms_info->q;
  double *clatoms_fx      = clatoms_pos->fx;
  double *clatoms_fy      = clatoms_pos->fy;
  double *clatoms_fz      = clatoms_pos->fz;
  double *clatoms_fxt     = clatoms_pos->fxt;
  double *clatoms_fyt     = clatoms_pos->fyt;
  double *clatoms_fzt     = clatoms_pos->fzt;
  double *ptens_pvten_tmp = ptens->pvten_tmp;
  double *ptens_pvten_tmp_res= ptens->pvten_tmp_res;
  double *xmod            = clatoms_info->xmod;
  double *ymod            = clatoms_info->ymod;
  double *zmod            = clatoms_info->zmod;
  int *ecor_j1            = ecor->j1;
  int *ecor_j2            = ecor->j2;
  int nsplin_m2           = ecor->nsplin_m2;
  int np_forc             = class_comm_forc_pkg->num_proc;
  int myid_forc           = class_comm_forc_pkg->myid;
  double dr_spl           = ecor->dr_spl;
  double dri_spl          = ecor->dri_spl;
  double rmin_spl         = ecor->rmin_spl;
  int iperd               = cell->iperd;
  int intra_perds         = cell->intra_perds;
  int pi_beads            = clatoms_info->pi_beads;
  int nlen_use            = (intra_scr->nlen);
  int ntot                = ecor->num;
  double *pvten           = ptens->pvten;
  double *pvten_tot       = ptens->pvten_tot;
  double wfor             = (intra_scr->wght_ter);
  double wfor_dif         = ((intra_scr->wght_ter_res)
                            -(intra_scr->wght_ter));
  double *ecor_cv0        = ecor->cv0;
  double *ecor_cdv0       = ecor->cdv0;
  double *ecor_cdv0_res   = ecor->cdv0_res;



/*==========================================================================*/
/* I) loop over all the ecors in steps of nlen to save memory               */

  for(i=1;i<=9;i++){ptens_pvten_tmp[i]=0;}
  for(i=1;i<=9;i++){ptens_pvten_tmp_res[i]=0;}


  for(ibig=1;ibig <= ntot;ibig += nlen_use) {

/*--------------------------------------------------------------------------*/
/*  A) Offsets to save some memory by doing only nlen ecors at a time       */

    ibig1 = ibig-1;
    ioff = -ibig1;
    iend = MIN(ntot,ibig1+nlen_use);
    nnow = iend-ibig1;

/*--------------------------------------------------------------------------*/
/*  B) Gather displacements                                                 */

    for(ieoff=1,iecor=ibig;iecor <= iend; ++iecor,++ieoff) {       
       ktemp1 = ecor_j1[iecor];
       ktemp2 = ecor_j2[iecor];
       intra_scr_dx12[ieoff] = clatoms_x[ktemp1]-clatoms_x[ktemp2];
       intra_scr_dy12[ieoff] = clatoms_y[ktemp1]-clatoms_y[ktemp2];
       intra_scr_dz12[ieoff] = clatoms_z[ktemp1]-clatoms_z[ktemp2];
    }/*endfor*/

    if(intra_perds==1&&pi_beads>1){     
       for(ieoff=1,iecor=ibig;iecor <= iend; ++iecor,++ieoff) {       
         ktemp5 = ecor_j1[iecor];
         ktemp6 = ecor_j2[iecor];
         intra_scr_dx56[ieoff] = xmod[ktemp5]-xmod[ktemp6];
         intra_scr_dy56[ieoff] = ymod[ktemp5]-ymod[ktemp6];
         intra_scr_dz56[ieoff] = zmod[ktemp5]-zmod[ktemp6];
       }/*endfor*/
    }/*endif*/

/*--------------------------------------------------------------------------*/
/*  C) Periodic boundary conditions                                         */

    if(intra_perds == 1) {
      if(pi_beads==1){     
        period(nnow,intra_scr_dx12,intra_scr_dy12,intra_scr_dz12,cell);
      }else{
        period_pimd(nnow,intra_scr_dx12,intra_scr_dy12,intra_scr_dz12,
                    intra_scr_dx56,intra_scr_dy56,intra_scr_dz56,cell);
      }/*endif*/
    }/*endif*/

/*--------------------------------------------------------------------------*/
/*  D) Spline info                                                          */

    for(iecor=1,iboff=ibig;iecor <= nnow; ++iecor,++iboff) {
       r122         = intra_scr_dx12[iecor]*intra_scr_dx12[iecor]
                    + intra_scr_dy12[iecor]*intra_scr_dy12[iecor]
                    + intra_scr_dz12[iecor]*intra_scr_dz12[iecor];
       r12          = sqrt(r122);
       kktemp       = (int) (((r12 - rmin_spl)*dri_spl)+0.5) + 3;
       kktemp       = MIN(kktemp,nsplin_m2);
       kktemp       = MAX(kktemp,3);
       rsp          = ((double)(kktemp-3))*dr_spl+rmin_spl;
       intra_scr_del_r[iecor]   = (r12 - rsp)*dri_spl;
       intra_scr_i_index[iecor] = kktemp;
       intra_scr_q12[iecor]     = 
                       clatoms_q[ecor_j1[iboff]]*clatoms_q[ecor_j2[iboff]];
    }/*endfor iecor*/

/*--------------------------------------------------------------------------*/
/*  E) Sum the potential energy                                             */

    if(iget_pe_real_inter==1){
     ecor_vspl_fetch(nnow,intra_scr_del_r,intra_scr_spl_tmp,
                     intra_scr_i_index,ecor_cv0);
     for(iecor=1;iecor<=nnow;++iecor)
           {*vecort += (intra_scr_spl_tmp)[iecor]*intra_scr_q12[iecor];}
    }/*endif*/

/*--------------------------------------------------------------------------*/
/*  F) Get the forces                                                       */

    ecor_vspl_fetch(nnow,intra_scr_del_r,intra_scr_spl_tmp,
                    intra_scr_i_index,ecor_cdv0);
    for(iecor = 1; iecor <= nnow; iecor++) {
       dvtemp = intra_scr_spl_tmp[iecor]*intra_scr_q12[iecor];
       intra_scr_fx1[iecor] = intra_scr_dx12[iecor]*dvtemp;
       intra_scr_fy1[iecor] = intra_scr_dy12[iecor]*dvtemp;
       intra_scr_fz1[iecor] = intra_scr_dz12[iecor]*dvtemp;
    }/*endfor*/

    ecor_vspl_fetch(nnow,intra_scr_del_r,intra_scr_spl_tmp,
                    intra_scr_i_index,ecor_cdv0_res);
    for(iecor = 1; iecor <= nnow; iecor++) {
       dvtemp = intra_scr_spl_tmp[iecor]*intra_scr_q12[iecor];
       intra_scr_fx2[iecor] = intra_scr_dx12[iecor]*dvtemp;
       intra_scr_fy2[iecor] = intra_scr_dy12[iecor]*dvtemp;
       intra_scr_fz2[iecor] = intra_scr_dz12[iecor]*dvtemp;
    }/*endfor*/

/*--------------------------------------------------------------------------*/
/*  G) Pressure tensor, etc.                                                */

    ecor_pvten_sum_roll(nnow,intra_scr_dx12,intra_scr_dy12,intra_scr_dz12,
                       intra_scr_fx1,intra_scr_fy1,intra_scr_fz1,
                       intra_scr_p11,intra_scr_p22,intra_scr_p33,intra_scr_p12,
                       intra_scr_p13,intra_scr_p23,ptens_pvten_tmp,iperd);

    ecor_pvten_sum_roll(nnow,intra_scr_dx12,intra_scr_dy12,intra_scr_dz12,
                       intra_scr_fx2,intra_scr_fy2,intra_scr_fz2,
                       intra_scr_p11,intra_scr_p22,intra_scr_p33,intra_scr_p12,
                       intra_scr_p13,intra_scr_p23,ptens_pvten_tmp_res,iperd);

/*--------------------------------------------------------------------------*/
/*  H) Scatter the forces                                                   */

    if(iver_get==1){

      for(iecor=ibig;iecor <= iend; ++iecor) {
         ktemp = ecor_j1[iecor];
         mtemp = iecor+ioff;
         clatoms_fxt[ktemp] += intra_scr_fx1[mtemp];
         clatoms_fyt[ktemp] += intra_scr_fy1[mtemp];
         clatoms_fzt[ktemp] += intra_scr_fz1[mtemp];
       }/*endfor*/

       for(iecor=ibig;iecor <= iend; ++iecor) {
         ktemp = ecor_j2[iecor];
         mtemp = iecor+ioff;
         clatoms_fxt[ktemp] -= intra_scr_fx1[mtemp];
         clatoms_fyt[ktemp] -= intra_scr_fy1[mtemp];
         clatoms_fzt[ktemp] -= intra_scr_fz1[mtemp];
       }/*endfor*/

    }/*endif: iverget*/

    for(iecor=1;iecor <= nnow; ++iecor) {
       intra_scr_fx1[iecor] = wfor*intra_scr_fx1[iecor]  
                            + wfor_dif*intra_scr_fx2[iecor];
       intra_scr_fy1[iecor] = wfor*intra_scr_fy1[iecor]  
                            + wfor_dif*intra_scr_fy2[iecor];
       intra_scr_fz1[iecor] = wfor*intra_scr_fz1[iecor]  
                            + wfor_dif*intra_scr_fz2[iecor];
    }/*endfor*/

    for(iecor=ibig;iecor <= iend; ++iecor) {
       ktemp = ecor_j1[iecor];
       mtemp = iecor+ioff;
       clatoms_fx[ktemp] += intra_scr_fx1[mtemp];
       clatoms_fy[ktemp] += intra_scr_fy1[mtemp];
       clatoms_fz[ktemp] += intra_scr_fz1[mtemp];
    }/*endfor*/

    for(iecor=ibig;iecor <= iend; ++iecor) {
       ktemp = ecor_j2[iecor];
       mtemp = iecor+ioff;
       clatoms_fx[ktemp] -= intra_scr_fx1[mtemp];
       clatoms_fy[ktemp] -= intra_scr_fy1[mtemp];
       clatoms_fz[ktemp] -= intra_scr_fz1[mtemp];
    }/*endfor*/

  }/* endfor ibig */

/*==========================================================================*/
/* Increment the Pressure tensor */

  ptens_pvten_tmp[4]     = ptens_pvten_tmp[2];
  ptens_pvten_tmp[7]     = ptens_pvten_tmp[3];
  ptens_pvten_tmp[8]     = ptens_pvten_tmp[6];
  ptens_pvten_tmp_res[4] = ptens_pvten_tmp_res[2];
  ptens_pvten_tmp_res[7] = ptens_pvten_tmp_res[3];
  ptens_pvten_tmp_res[8] = ptens_pvten_tmp_res[6];
  for(i=1;i<=9;i++){
    pvten[i] += (ptens_pvten_tmp[i]*wfor + ptens_pvten_tmp_res[i]*wfor_dif);
  }/*endfor*/

  if(iget_pv_real_inter==1){    
    for(i=1;i<=9;i++){pvten_tot[i] += ptens_pvten_tmp[i];}
  }/*endif*/

/*--------------------------------------------------------------------------*/
  }/* end routine */
/*==========================================================================*/






/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void ecor_vspl_fetch(int n,double del[], double spl_out[],
                    int i_index[],double c0_data[])

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
    spl_out[i]  =  (c0 + del_tmp*(c1 + del_tmp* (c2 + del_tmp*c3)));
  }/*endfor*/


/*--------------------------------------------------------------------------*/
}/* end routine vspl_fetch */
/*==========================================================================*/







/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void ecor_pvten_sum_roll(int n,double *dx12, double *dy12, double *dz12,
                         double *fx,double *fy,double *fz,
                         double *p11,double *p22,double *p33,double *p12,
                         double *p13,double *p23,double *pvten_tmp,int iperd)

/*=======================================================================*/
/*              Begin Routine                                            */
    {/*Begin Routine*/
/*=======================================================================*/
/*            Local variable declarations                                */

int i,n2;
  
/*=======================================================================*/

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

/*--------------------------------------------------------------------------*/
}/* end routine */
/*==========================================================================*/




