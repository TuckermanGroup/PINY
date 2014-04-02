/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
/*                                                                          */
/*                         PI_MD:                                           */
/*             The future of simulation technology                          */
/*             ------------------------------------                         */
/*                     Module: bond_both.c                                  */
/*                                                                          */
/* This routine computes the energy and forces from                         */
/* the intramolecular bond potential.                                       */
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

void bond_both(CLATOMS_INFO *clatoms_info,CLATOMS_POS *clatoms_pos,
          BOND *bond,CELL *cell,
          INTRA_SCR *intra_scr,PTENS *ptens, double *vbondt, int iver_get,
          CLASS_COMM_FORC_PKG *class_comm_forc_pkg,int iget_pv_real_inter)

/*==========================================================================*/
{/*begin routine*/
  /*=======================================================================*/
  /*            Local variable declarations                                */
  int ibond,ibig,ibig1,ioff,iend,nnow;     /* Indices and counters   */
  double r122,r12;                         /* (r-r_0)^2 and (r-r_0)  */
  double vbond;                            /* Bond potential         */
  double dvbond;                           /* Derivative of bond pot */
  double pre_res;                          /* RESPA Force prefactor        */
  double dvbond_res;                       /* Derivative of RESPA bond pot */
  double pre;                              /* Force prefactor        */
  double rr0,rr0_res;                      /* (r-r_0)                */
  double wfor,wfor_res;
  double temp;
  double req,req_res;
  int i,ktemp,iboff,iii;
  int nlen_use,nlen_now;
  int ntot;
  int ngo,irem,istart_big;

  /* Define local pointers                                                 */
  double *intra_scr_x1     = intra_scr->x1;
  double *intra_scr_y1     = intra_scr->y1;
  double *intra_scr_z1     = intra_scr->z1;
  double *intra_scr_x2     = intra_scr->x2;
  double *intra_scr_y2     = intra_scr->y2;
  double *intra_scr_z2     = intra_scr->z2;
  double *intra_scr_x5     = intra_scr->x5;
  double *intra_scr_y5     = intra_scr->y5;
  double *intra_scr_z5     = intra_scr->z5;
  double *intra_scr_x6     = intra_scr->x6;
  double *intra_scr_y6     = intra_scr->y6;
  double *intra_scr_z6     = intra_scr->z6;
  double *intra_scr_fx1    = intra_scr->fx1;
  double *intra_scr_fy1    = intra_scr->fy1;
  double *intra_scr_fz1    = intra_scr->fz1;
  double *intra_scr_fx2    = intra_scr->fx2;
  double *intra_scr_fy2    = intra_scr->fy2;
  double *intra_scr_fz2    = intra_scr->fz2;
  double *intra_scr_dx12   = intra_scr->dx12;
  double *intra_scr_dy12   = intra_scr->dy12;
  double *intra_scr_dz12   = intra_scr->dz12;
  double *intra_scr_dx56   = intra_scr->dx56;
  double *intra_scr_dy56   = intra_scr->dy56;
  double *intra_scr_dz56   = intra_scr->dz56;
  double *intra_scr_vpot   = intra_scr->vpot;
  double *intra_scr_eq     = intra_scr->eq;
  double *intra_scr_eq_res = intra_scr->fx4;
  double *intra_scr_c_0    = intra_scr->c_0;
  double *intra_scr_c_1    = intra_scr->c_1;
  double *intra_scr_c_2    = intra_scr->c_2;
  double *intra_scr_c_3    = intra_scr->c_3;
  double *intra_scr_c_4    = intra_scr->c_4;
  double *intra_scr_c_5    = intra_scr->c_5;
  double *intra_scr_c_6    = intra_scr->c_6;
  double *intra_scr_dc_0   = intra_scr->dc_0;
  double *intra_scr_dc_1   = intra_scr->dc_1;
  double *intra_scr_dc_2   = intra_scr->dc_2;
  double *intra_scr_dc_3   = intra_scr->dc_3;
  double *intra_scr_dc_4   = intra_scr->dc_4;
  double *intra_scr_dc_5   = intra_scr->dc_5;
  double *intra_scr_dc_6   = intra_scr->dc_6;
  double *intra_scr_p11    = intra_scr->p11;
  double *intra_scr_p22    = intra_scr->p22;
  double *intra_scr_p33    = intra_scr->p33;
  double *intra_scr_p21    = intra_scr->p21;
  double *intra_scr_p12    = intra_scr->p12;
  double *intra_scr_p31    = intra_scr->p31;
  double *intra_scr_p13    = intra_scr->p13;
  double *intra_scr_p32    = intra_scr->p32;
  double *intra_scr_p23    = intra_scr->p23;
  double *intra_scr_p11_res    = intra_scr->dx13;
  double *intra_scr_p22_res    = intra_scr->dy13;
  double *intra_scr_p33_res    = intra_scr->dz13;
  double *intra_scr_p21_res    = intra_scr->dx23;
  double *intra_scr_p12_res    = intra_scr->dy23;
  double *intra_scr_p31_res    = intra_scr->dz23;
  double *intra_scr_p13_res    = intra_scr->dx42;
  double *intra_scr_p32_res    = intra_scr->dy42;
  double *intra_scr_p23_res    = intra_scr->dz42;
  double *clatoms_x        = clatoms_pos->x;
  double *clatoms_y        = clatoms_pos->y;
  double *clatoms_z        = clatoms_pos->z;
  double *clatoms_fx       = clatoms_pos->fx;
  double *clatoms_fy       = clatoms_pos->fy;
  double *clatoms_fz       = clatoms_pos->fz;
  double *clatoms_fxt      = clatoms_pos->fxt;
  double *clatoms_fyt      = clatoms_pos->fyt;
  double *clatoms_fzt      = clatoms_pos->fzt;
  double *bond_eq_pow      = bond->eq_pow;
  double *bond_eq_pow_res  = bond->eq_pow_res;
  double *bond_c_0         = bond->c_0;
  double *bond_c_1         = bond->c_1;
  double *bond_c_2         = bond->c_2;
  double *bond_c_3         = bond->c_3;
  double *bond_c_4         = bond->c_4;
  double *bond_c_5         = bond->c_5;
  double *bond_c_6         = bond->c_6;
  double *bond_dc_0        = bond->dc_0;
  double *bond_dc_1        = bond->dc_1;
  double *bond_dc_2        = bond->dc_2;
  double *bond_dc_3        = bond->dc_3;
  double *bond_dc_4        = bond->dc_4;
  double *bond_dc_5        = bond->dc_5;
  double *bond_dc_6        = bond->dc_6;
  double *ptens_pvten_tmp      = ptens->pvten_tmp;
  double *ptens_pvten_tmp_res  = ptens->pvten_tmp_res;
  double *xmod            = clatoms_info->xmod;
  double *ymod            = clatoms_info->ymod;
  double *zmod            = clatoms_info->zmod;
  int *bond_j1_pow         = bond->j1_pow;
  int *bond_j2_pow         = bond->j2_pow;
  int *bond_jtyp_pow       = bond->jtyp_pow;
  int *iblock_size         = bond->iblock_pow_size;
  int *iblock_conflict_1   = bond->iblock_pow_conflict_1;
  int *iblock_conflict_2   = bond->iblock_pow_conflict_2;
  int np_forc             = class_comm_forc_pkg->num_proc;
  int myid_forc           = class_comm_forc_pkg->myid;

  /*=======================================================================*/
  /* I) loop over all the bonds in steps of nlen to save memory            */

  wfor_res    = (intra_scr->wght_tra_res);
  if(intra_scr->int_res_inter == 1){
    wfor        = (intra_scr->wght_ter_res);
  }else{
    wfor        = (intra_scr->wght_ter);
  }/*endif*/
  for(i=1;i<=9;i++){(ptens_pvten_tmp)[i]=0.0;}
  for(i=1;i<=9;i++){(ptens_pvten_tmp_res)[i]=0.0;}

  nlen_use = (intra_scr->nlen);
  ntot = bond->npow;

  for(ibig=1;ibig <= ntot;ibig += nlen_use) {
    
    /*----------------------------------------------------------------------*/
    /*  A) Offsets to save some memory by doing only nlen bonds at a time   */
    
    ibig1 = ibig-1;
    ioff = -ibig1;
    nlen_now = nlen_use;
    iend = MIN(ntot,ibig1+nlen_now);
    nnow = iend-ibig1;
    
    /*----------------------------------------------------------------------*/
    /*  B) Gather positions                                                 */

    for(iboff=1,ibond=ibig;ibond <= iend; ++ibond,++iboff) { 
       ktemp = bond_j1_pow[ibond];
       intra_scr_x1[iboff] = clatoms_x[ktemp];
       intra_scr_y1[iboff] = clatoms_y[ktemp];
       intra_scr_z1[iboff] = clatoms_z[ktemp];
    }/*endfor*/
    for(iboff=1,ibond=ibig;ibond <= iend; ++ibond,++iboff) { 
      ktemp = bond_j2_pow[ibond];
      intra_scr_x2[iboff] = clatoms_x[ktemp];
      intra_scr_y2[iboff] = clatoms_y[ktemp];
      intra_scr_z2[iboff] = clatoms_z[ktemp];
    }/*endfor*/

  if(cell->intra_perds==1&&clatoms_info->pi_beads>1){     
    for(iboff=1,ibond=ibig;ibond <= iend; ++ibond,++iboff) { 
       ktemp = bond_j1_pow[ibond];
       intra_scr_x5[iboff] = xmod[ktemp];
       intra_scr_y5[iboff] = ymod[ktemp];
       intra_scr_z5[iboff] = zmod[ktemp];
    }/*endfor*/
    for(iboff=1,ibond=ibig;ibond <= iend; ++ibond,++iboff) { 
      ktemp = bond_j2_pow[ibond];
      intra_scr_x6[iboff] = xmod[ktemp];
      intra_scr_y6[iboff] = ymod[ktemp];
      intra_scr_z6[iboff] = zmod[ktemp];
    }/*endfor*/
  }/*endif*/


    /*----------------------------------------------------------------------*/
    /*  C) Gather the power series coefficients and equilibrium positions   */
    
    for(iboff=1,ibond=ibig;ibond <= iend; ++ibond,++iboff) {       
      ktemp = bond_jtyp_pow[ibond];

      intra_scr_eq[iboff]     = bond_eq_pow[ktemp];
      intra_scr_eq_res[iboff] = bond_eq_pow_res[ktemp];
      intra_scr_c_0[iboff] = bond_c_0[ktemp];
      intra_scr_c_1[iboff] = bond_c_1[ktemp];
      intra_scr_c_2[iboff] = bond_c_2[ktemp];
      intra_scr_c_3[iboff] = bond_c_3[ktemp];
      intra_scr_c_4[iboff] = bond_c_4[ktemp];
      intra_scr_c_5[iboff] = bond_c_5[ktemp];
      intra_scr_c_6[iboff] = bond_c_6[ktemp];
      
      intra_scr_dc_0[iboff] = bond_dc_0[ktemp];
      intra_scr_dc_1[iboff] = bond_dc_1[ktemp];
      intra_scr_dc_2[iboff] = bond_dc_2[ktemp];
      intra_scr_dc_3[iboff] = bond_dc_3[ktemp];
      intra_scr_dc_4[iboff] = bond_dc_4[ktemp];
      intra_scr_dc_5[iboff] = bond_dc_5[ktemp];
      intra_scr_dc_6[iboff] = bond_dc_6[ktemp];
    }/*endfor*/
    
    /*----------------------------------------------------------------------*/
    /*  D) Calculate the basis vectors r12                                  */
    
    for(ibond=1;ibond <= nnow; ++ibond) {
      intra_scr_dx12[ibond] = intra_scr_x1[ibond] - intra_scr_x2[ibond];
      intra_scr_dy12[ibond] = intra_scr_y1[ibond] - intra_scr_y2[ibond];
      intra_scr_dz12[ibond] = intra_scr_z1[ibond] - intra_scr_z2[ibond];
    }/*endfor*/
    
  if(cell->intra_perds==1&&clatoms_info->pi_beads>1){     
    for(ibond=1;ibond <= nnow; ++ibond) {
      intra_scr_dx56[ibond] = intra_scr_x5[ibond] - intra_scr_x6[ibond];
      intra_scr_dy56[ibond] = intra_scr_y5[ibond] - intra_scr_y6[ibond];
      intra_scr_dz56[ibond] = intra_scr_z5[ibond] - intra_scr_z6[ibond];
    }/*endfor*/
  }/*endif*/
    
    /*---------------------------------------------------------------------*/
    /*  E) Periodic boundary conditions                                    */
    
    if(cell->intra_perds == 1) {
      if(clatoms_info->pi_beads==1){     
        period(nnow,intra_scr_dx12,intra_scr_dy12,intra_scr_dz12,cell);
      }else{
        period_pimd(nnow,intra_scr_dx12,intra_scr_dy12,intra_scr_dz12,
                         intra_scr_dx56,intra_scr_dy56,intra_scr_dz56,cell);
      }
    }
    
    for(ibond=1;ibond <= nnow; ++ibond) {
      r122 = intra_scr_dx12[ibond]*intra_scr_dx12[ibond]
        + intra_scr_dy12[ibond]*intra_scr_dy12[ibond]
          + intra_scr_dz12[ibond]*intra_scr_dz12[ibond];
      r12  = sqrt(r122);
/*      if(bond_jtyp_pow[ibond]==1){
        rbar_1 += r12;
        ntyp1++;
      }
      if(bond_jtyp_pow[ibond]==2){
        rbar_2 += r12;
        ntyp2++;
      }*/

      req     = intra_scr_eq[ibond];
      req_res = intra_scr_eq_res[ibond];      

      rr0     = r12 - req;
      rr0_res = r12 - req_res;

      /*--------------------------------------------------------------------*/
      /*  F) Get the bonding potential energy using Horner's method         */
      
      vbond =  ((((( intra_scr_c_6[ibond]
                    *rr0 + intra_scr_c_5[ibond])
                   *rr0 + intra_scr_c_4[ibond])
                  *rr0 + intra_scr_c_3[ibond])
                 *rr0 + intra_scr_c_2[ibond])
                *rr0 + intra_scr_c_1[ibond])
        *rr0 + intra_scr_c_0[ibond];
      
      intra_scr_vpot[ibond] = vbond;
      
      /*---------------------------------------------------------------------*/
      /*  G) Get the force on the atoms using the chain rule */
      
      dvbond =   (((( intra_scr_dc_6[ibond]
                     *rr0 + intra_scr_dc_5[ibond])
                    *rr0 + intra_scr_dc_4[ibond])
                   *rr0 + intra_scr_dc_3[ibond])
                  *rr0 + intra_scr_dc_2[ibond])
        *rr0 + intra_scr_dc_1[ibond];
      
      dvbond_res =       (((( intra_scr_dc_6[ibond]
                         *rr0_res + intra_scr_dc_5[ibond])
                        *rr0_res + intra_scr_dc_4[ibond])
                       *rr0_res + intra_scr_dc_3[ibond])
                      *rr0_res + intra_scr_dc_2[ibond])
            *rr0_res + intra_scr_dc_1[ibond];
      
      pre        = -dvbond/r12;
      pre_res    = -dvbond_res/r12;
      intra_scr_fx1[ibond] = intra_scr_dx12[ibond]*pre;
      intra_scr_fy1[ibond] = intra_scr_dy12[ibond]*pre;
      intra_scr_fz1[ibond] = intra_scr_dz12[ibond]*pre;
      intra_scr_fx2[ibond] = intra_scr_dx12[ibond]*pre_res;
      intra_scr_fy2[ibond] = intra_scr_dy12[ibond]*pre_res;
      intra_scr_fz2[ibond] = intra_scr_dz12[ibond]*pre_res;
      /*endfor ibond*/}
    
    /*----------------------------------------------------------------------*/
    /*  H) Sum the potential energy                                         */
    
    for(ibond=1;ibond <= nnow; ++ibond){
      *vbondt += intra_scr_vpot[ibond];
    }
    
    /*----------------------------------------------------------------------*/
    /*  J) Pressure tensor, etc.                                            */
    
    if(cell->iperd == 2 || cell->iperd == 3) {
      for(ibond=1;ibond <= nnow; ++ibond) {
        intra_scr_p11[ibond] = intra_scr_dx12[ibond]*intra_scr_fx1[ibond];
        intra_scr_p22[ibond] = intra_scr_dy12[ibond]*intra_scr_fy1[ibond];
        intra_scr_p33[ibond] = intra_scr_dz12[ibond]*intra_scr_fz1[ibond];
        intra_scr_p12[ibond] = intra_scr_dx12[ibond]*intra_scr_fy1[ibond];
        intra_scr_p11_res[ibond] = intra_scr_dx12[ibond]*intra_scr_fx2[ibond];
        intra_scr_p22_res[ibond] = intra_scr_dy12[ibond]*intra_scr_fy2[ibond];
        intra_scr_p33_res[ibond] = intra_scr_dz12[ibond]*intra_scr_fz2[ibond];
        intra_scr_p12_res[ibond] = intra_scr_dx12[ibond]*intra_scr_fy2[ibond];
      }/*endfor*/
      for(ibond=1;ibond <= nnow; ++ibond) {
/*        if(bond_jtyp_pow[ibond]==1){
         pbond_1+=intra_scr_p11[ibond]+
                  intra_scr_p22[ibond]+
                  intra_scr_p33[ibond];
        }else{
         pbond_2+=intra_scr_p11[ibond]+
                  intra_scr_p22[ibond]+
                  intra_scr_p33[ibond];
        }*/
        ptens_pvten_tmp[1] += intra_scr_p11[ibond];
        ptens_pvten_tmp[5] += intra_scr_p22[ibond];
        ptens_pvten_tmp[9] += intra_scr_p33[ibond];
        ptens_pvten_tmp[2] += intra_scr_p12[ibond];
        ptens_pvten_tmp_res[1] += intra_scr_p11_res[ibond];
        ptens_pvten_tmp_res[5] += intra_scr_p22_res[ibond];
        ptens_pvten_tmp_res[9] += intra_scr_p33_res[ibond];
        ptens_pvten_tmp_res[2] += intra_scr_p12_res[ibond];
      }/*endfor*/
    /*endif*/}
    if((cell->iperd) == 3) {
      for(ibond=1;ibond <= nnow; ++ibond) {
        intra_scr_p13[ibond] = intra_scr_dx12[ibond]*intra_scr_fz1[ibond];
        intra_scr_p23[ibond] = intra_scr_dy12[ibond]*intra_scr_fz1[ibond];
        intra_scr_p13_res[ibond] = intra_scr_dx12[ibond]*intra_scr_fz2[ibond];
        intra_scr_p23_res[ibond] = intra_scr_dy12[ibond]*intra_scr_fz2[ibond];
      }/*endfor*/
      for(ibond=1;ibond <= nnow; ++ibond) {
        ptens_pvten_tmp[3] += intra_scr_p13[ibond];
        ptens_pvten_tmp[6] += intra_scr_p23[ibond];
        ptens_pvten_tmp_res[3] += intra_scr_p13_res[ibond];
        ptens_pvten_tmp_res[6] += intra_scr_p23_res[ibond];
      }/*endfor*/
    /*endif*/}

    /*----------------------------------------------------------------------*/
    /*  I) Virial estimator                                               */

    if(iver_get==1){

      for(iboff=1,ibond=ibig;ibond <= iend; ++ibond,++iboff) { 
       ktemp = bond_j1_pow[ibond];
       clatoms_fxt[ktemp] += intra_scr_fx1[iboff];
       clatoms_fyt[ktemp] += intra_scr_fy1[iboff];
       clatoms_fzt[ktemp] += intra_scr_fz1[iboff];
      }/*endfor*/

      for(iboff=1,ibond=ibig;ibond <= iend; ++ibond,++iboff) { 
       ktemp = bond_j2_pow[ibond];
       clatoms_fxt[ktemp] += -intra_scr_fx1[iboff];
       clatoms_fyt[ktemp] += -intra_scr_fy1[iboff];
       clatoms_fzt[ktemp] += -intra_scr_fz1[iboff];
      }/*endfor*/
    }/*endif : iver_get*/

    /*----------------------------------------------------------------------*/
    /*  I) Scatter the forces                                               */

    for(iboff=1;iboff <= iend+ioff; ++iboff) { 
     intra_scr_fx1[iboff] = wfor_res*intra_scr_fx2[iboff]
                         + wfor*(intra_scr_fx1[iboff]-intra_scr_fx2[iboff]);
     intra_scr_fy1[iboff] = wfor_res*intra_scr_fy2[iboff]
                         + wfor*(intra_scr_fy1[iboff]-intra_scr_fy2[iboff]);
     intra_scr_fz1[iboff] = wfor_res*intra_scr_fz2[iboff]
                         + wfor*(intra_scr_fz1[iboff]-intra_scr_fz2[iboff]);
    }/*endfor*/

     for(iboff=1,ibond=ibig;ibond <= iend; ++ibond,++iboff) { 
      ktemp = bond_j1_pow[ibond];
      clatoms_fx[ktemp] += intra_scr_fx1[iboff];
      clatoms_fy[ktemp] += intra_scr_fy1[iboff];
      clatoms_fz[ktemp] += intra_scr_fz1[iboff];
     }/*endfor*/

     for(iboff=1,ibond=ibig;ibond <= iend; ++ibond,++iboff) { 
      ktemp = bond_j2_pow[ibond];
      clatoms_fx[ktemp] += -intra_scr_fx1[iboff];
      clatoms_fy[ktemp] += -intra_scr_fy1[iboff];
      clatoms_fz[ktemp] += -intra_scr_fz1[iboff];
     }/*endfor*/
  }/* endfor ibig */
  /*=======================================================================*/
  /* II) Increment the Pressure tensor */

  ptens_pvten_tmp[4]     = ptens_pvten_tmp[2];
  ptens_pvten_tmp[7]     = ptens_pvten_tmp[3];
  ptens_pvten_tmp[8]     = ptens_pvten_tmp[6];
  ptens_pvten_tmp_res[4] = ptens_pvten_tmp_res[2];
  ptens_pvten_tmp_res[7] = ptens_pvten_tmp_res[3];
  ptens_pvten_tmp_res[8] = ptens_pvten_tmp_res[6];

  for(i=1;i<=9;i++){
    ptens->pvten[i]     += (wfor_res*ptens_pvten_tmp_res[i]
                  + wfor*(ptens_pvten_tmp[i]-ptens_pvten_tmp_res[i]));
  }/*endfor*/

  if(iget_pv_real_inter==1){    
   for(i=1;i<=9;i++){
    ptens->pvten_tot[i] += ptens_pvten_tmp[i];
   }/*endfor*/ 
  }/*endif*/


/*--------------------------------------------------------------------------*/
}/* end routine */
/*==========================================================================*/






