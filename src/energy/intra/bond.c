/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
/*                                                                          */
/*                         PI_MD:                                           */
/*             The future of simulation technology                          */
/*             ------------------------------------                         */
/*                     Module: bond.c                                       */
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

void bond(CLATOMS_INFO *clatoms_info,CLATOMS_POS *clatoms_pos,
          BOND *bond,CELL *cell,
          INTRA_SCR *intra_scr,PTENS *ptens, double *vbondt, int iver_get,
          int irespa,CLASS_COMM_FORC_PKG *class_comm_forc_pkg,
          int iget_pv_real_inter)

/*==========================================================================*/
{/*begin routine*/
  /*=======================================================================*/
  /*            Local variable declarations                                */
  int ibond,ibig,ibig1,ioff,iend,nnow;     /* Indices and counters   */
  double r122,r12;                         /* (r-r_0)^2 and (r-r_0)  */
  double vbond;                            /* Bond potential         */
  double dvbond;                           /* Derivative of bond pot */
  double pre;                              /* Force prefactor        */
  double rr0;                              /* (r-r_0)                */
  double wfor;
  double temp;
  int i,ktemp,iboff,iii;
  int nlen_use,nlen_now;
  int ktemp1,ktemp2;
  int ntot;
  int ngo,irem,istart_big;

  /* Define local pointers                                                 */
  double *intra_scr_fx1   = intra_scr->fx1;
  double *intra_scr_fy1   = intra_scr->fy1;
  double *intra_scr_fz1   = intra_scr->fz1;
  double *intra_scr_dx12  = intra_scr->dx12;
  double *intra_scr_dy12  = intra_scr->dy12;
  double *intra_scr_dz12  = intra_scr->dz12;
  double *intra_scr_dx56  = intra_scr->dx56;
  double *intra_scr_dy56  = intra_scr->dy56;
  double *intra_scr_dz56  = intra_scr->dz56;
  double *intra_scr_vpot  = intra_scr->vpot;
  double *intra_scr_p11   = intra_scr->p11;
  double *intra_scr_p22   = intra_scr->p22;
  double *intra_scr_p33   = intra_scr->p33;
  double *intra_scr_p21   = intra_scr->p21;
  double *intra_scr_p12   = intra_scr->p12;
  double *intra_scr_p31   = intra_scr->p31;
  double *intra_scr_p13   = intra_scr->p13;
  double *intra_scr_p32   = intra_scr->p32;
  double *intra_scr_p23   = intra_scr->p23;
  double *clatoms_x       = clatoms_pos->x;
  double *clatoms_y       = clatoms_pos->y;
  double *clatoms_z       = clatoms_pos->z;
  double *clatoms_fx      = clatoms_pos->fx;
  double *clatoms_fy      = clatoms_pos->fy;
  double *clatoms_fz      = clatoms_pos->fz;
  double *clatoms_fxt      = clatoms_pos->fxt;
  double *clatoms_fyt      = clatoms_pos->fyt;
  double *clatoms_fzt      = clatoms_pos->fzt;
  double *bond_eq_pow;
  double *bond_c_0        = bond->c_0;
  double *bond_c_1        = bond->c_1;
  double *bond_c_2        = bond->c_2;
  double *bond_c_3        = bond->c_3;
  double *bond_c_4        = bond->c_4;
  double *bond_c_5        = bond->c_5;
  double *bond_c_6        = bond->c_6;
  double *bond_dc_0       = bond->dc_0;
  double *bond_dc_1       = bond->dc_1;
  double *bond_dc_2       = bond->dc_2;
  double *bond_dc_3       = bond->dc_3;
  double *bond_dc_4       = bond->dc_4;
  double *bond_dc_5       = bond->dc_5;
  double *bond_dc_6       = bond->dc_6;
  double *ptens_pvten_tmp = ptens->pvten_tmp;
  double *xmod            = clatoms_info->xmod;
  double *ymod            = clatoms_info->ymod;
  double *zmod            = clatoms_info->zmod;
  int *bond_j1_pow        = bond->j1_pow;
  int *bond_j2_pow        = bond->j2_pow;
  int *bond_jtyp_pow      = bond->jtyp_pow;
  int *iblock_size        = bond->iblock_pow_size;
  int *iblock_conflict_1  = bond->iblock_pow_conflict_1;
  int *iblock_conflict_2  = bond->iblock_pow_conflict_2;
  int np_forc             = class_comm_forc_pkg->num_proc;
  int myid_forc           = class_comm_forc_pkg->myid;

  /*=======================================================================*/
  /* I) loop over all the bonds in steps of nlen to save memory            */

  if(irespa==1){
     bond_eq_pow  = bond->eq_pow_res;
   }else{
     bond_eq_pow  = bond->eq_pow;
   }/*endif*/
  wfor    = (intra_scr->wght_tra_res);
  for(i=1;i<=9;i++){(ptens_pvten_tmp)[i]=0.0;}

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
       ktemp1 = bond_j1_pow[ibond];
       ktemp2 = bond_j2_pow[ibond];
       intra_scr_dx12[iboff] = clatoms_x[ktemp1] - clatoms_x[ktemp2];
       intra_scr_dy12[iboff] = clatoms_y[ktemp1] - clatoms_y[ktemp2];
       intra_scr_dz12[iboff] = clatoms_z[ktemp1] - clatoms_z[ktemp2];
    }/*endfor*/

  if(cell->intra_perds==1&&clatoms_info->pi_beads>1){     
    for(iboff=1,ibond=ibig;ibond <= iend; ++ibond,++iboff) { 
       ktemp1 = bond_j1_pow[ibond];
       ktemp2 = bond_j2_pow[ibond];
       intra_scr_dx56[iboff] = xmod[ktemp1] - xmod[ktemp2];
       intra_scr_dy56[iboff] = ymod[ktemp1] - ymod[ktemp2];
       intra_scr_dz56[iboff] = zmod[ktemp1] - zmod[ktemp2];
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

    for(ibond=1,iboff=ibig;ibond <= nnow; ++ibond,++iboff) {

      ktemp = bond_jtyp_pow[iboff];

      r122 = intra_scr_dx12[ibond]*intra_scr_dx12[ibond]
        + intra_scr_dy12[ibond]*intra_scr_dy12[ibond]
          + intra_scr_dz12[ibond]*intra_scr_dz12[ibond];
      r12  = sqrt(r122);

      rr0 = r12 - bond_eq_pow[ktemp];
      /*--------------------------------------------------------------------*/
      /*  F) Get the bonding potential energy using Horner's method         */
      
      vbond =  ((((( bond_c_6[ktemp]
                    *rr0 + bond_c_5[ktemp])
                   *rr0 + bond_c_4[ktemp])
                  *rr0 + bond_c_3[ktemp])
                 *rr0 + bond_c_2[ktemp])
                *rr0 + bond_c_1[ktemp])
        *rr0 + bond_c_0[ktemp];
      
      intra_scr_vpot[ibond] = vbond;
      
      /*---------------------------------------------------------------------*/
      /*  G) Get the force on the atoms using the chain rule */
      
      dvbond =   (((( bond_dc_6[ktemp]
                     *rr0 + bond_dc_5[ktemp])
                    *rr0 + bond_dc_4[ktemp])
                   *rr0 + bond_dc_3[ktemp])
                  *rr0 + bond_dc_2[ktemp])
        *rr0 + bond_dc_1[ktemp];
      
      pre    = -dvbond/r12;
      intra_scr_fx1[ibond] = intra_scr_dx12[ibond]*pre;
      intra_scr_fy1[ibond] = intra_scr_dy12[ibond]*pre;
      intra_scr_fz1[ibond] = intra_scr_dz12[ibond]*pre;
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
      }/*endfor*/
    /*endif*/}
    if((cell->iperd) == 3) {
      for(ibond=1;ibond <= nnow; ++ibond) {
        intra_scr_p13[ibond] = intra_scr_dx12[ibond]*intra_scr_fz1[ibond];
        intra_scr_p23[ibond] = intra_scr_dy12[ibond]*intra_scr_fz1[ibond];
      }/*endfor*/
      for(ibond=1;ibond <= nnow; ++ibond) {
        ptens_pvten_tmp[3] += intra_scr_p13[ibond];
        ptens_pvten_tmp[6] += intra_scr_p23[ibond];
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
    }/*endif*/

    /*----------------------------------------------------------------------*/
    /*  I) Scatter the forces                                               */

    for(iboff=1;iboff <= iend+ioff; ++iboff) { 
     intra_scr_fx1[iboff] *= wfor;
     intra_scr_fy1[iboff] *= wfor;
     intra_scr_fz1[iboff] *= wfor;
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

  ptens_pvten_tmp[4] = ptens_pvten_tmp[2];
  ptens_pvten_tmp[7] = ptens_pvten_tmp[3];
  ptens_pvten_tmp[8] = ptens_pvten_tmp[6];

  for(i=1;i<=9;i++){
    ptens->pvten[i]     += ptens_pvten_tmp[i]*wfor;
  }/*endfor*/

  if(iget_pv_real_inter==1){
   for(i=1;i<=9;i++){
    ptens->pvten_tot[i] += ptens_pvten_tmp[i];
   }/*endfor*/
  }/*endif*/

/*--------------------------------------------------------------------------*/
}/* end routine */
/*==========================================================================*/






/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void bond_free(CLATOMS_INFO *clatoms_info,CLATOMS_POS *clatoms_pos,
               BOND_FREE *bond_free,CELL *cell,
               INTRA_SCR *intra_scr,PTENS *ptens, double *vbondt,
               ENERGY_CTRL *energy_ctrl,int np_forc)

/*==========================================================================*/
{/*begin routine*/
  /*======================================================================*/
  /*           Local variable declarations                                */
  int ibin,nnow;
  int pi_beads = clatoms_info->pi_beads; 
  double r122,r12;                         /* (r-r_0)^2 and (r-r_0)  */
  double dvbond;                           /* Derivative of bond pot */
  double pre;                              /* Force prefactor        */
  double rr0;                              /* (r-r_0)                */
  double apow,apow1,temp;
  double wfor;
  double *xmod            = clatoms_info->xmod;
  double *ymod            = clatoms_info->ymod;
  double *zmod            = clatoms_info->zmod;
  int iget_pv_real_inter = energy_ctrl->iget_pv_real_inter;

  /*=======================================================================*/
  /* I) Get the two particles */
  
  wfor             = intra_scr->wght_tra_res;
  intra_scr->x1[1] = clatoms_pos->x[bond_free->j1];
  intra_scr->y1[1] = clatoms_pos->y[bond_free->j1];
  intra_scr->z1[1] = clatoms_pos->z[bond_free->j1];
  intra_scr->x2[1] = clatoms_pos->x[bond_free->j2];
  intra_scr->y2[1] = clatoms_pos->y[bond_free->j2];
  intra_scr->z2[1] = clatoms_pos->z[bond_free->j2];

  if(cell->intra_perds==1&&clatoms_info->pi_beads>1){     
    intra_scr->x5[1] = xmod[bond_free->j1];
    intra_scr->y5[1] = ymod[bond_free->j1];
    intra_scr->z5[1] = zmod[bond_free->j1];
    intra_scr->x6[1] = xmod[bond_free->j2];
    intra_scr->y6[1] = ymod[bond_free->j2];
    intra_scr->z6[1] = zmod[bond_free->j2];
  }/*endif*/  

  /*------------------------------------------------------------------------*/
  /* III) Calculate the basis vectors r12                                   */
  
  intra_scr->dx12[1] = intra_scr->x1[1] - intra_scr->x2[1];
  intra_scr->dy12[1] = intra_scr->y1[1] - intra_scr->y2[1];
  intra_scr->dz12[1] = intra_scr->z1[1] - intra_scr->z2[1];

  if(cell->intra_perds==1&&clatoms_info->pi_beads>1){     
    intra_scr->dx56[1] = intra_scr->x5[1] - intra_scr->x6[1];
    intra_scr->dy56[1] = intra_scr->y5[1] - intra_scr->y6[1];
    intra_scr->dz56[1] = intra_scr->z5[1] - intra_scr->z6[1];
  }/*endif*/

  /*------------------------------------------------------------------------*/
  /*  E) Periodic boundary conditions                                       */
  
  nnow = 1;
    if(cell->intra_perds == 1) {
      if(clatoms_info->pi_beads==1){     
        period(nnow,intra_scr->dx12,intra_scr->dy12,intra_scr->dz12,cell);
      }else{
        period_pimd(nnow,intra_scr->dx12,intra_scr->dy12,intra_scr->dz12,
                         intra_scr->dx56,intra_scr->dy56,intra_scr->dz56,cell);
      }
    }

  /*------------------------------------------------------------------------*/
  /*  F) Get the bonding potential energy and fill histogram                */
  
  r122      = (intra_scr->dx12[1]*intra_scr->dx12[1]
               + intra_scr->dy12[1]*intra_scr->dy12[1]
               + intra_scr->dz12[1]*intra_scr->dz12[1]);
  r12       = sqrt(r122);
  rr0       = r12 - (bond_free->eq);
  apow      = (double)bond_free->npow;
  apow1     = apow-1.0;
  temp      = pow(rr0,apow);
  *vbondt   = (bond_free->fk)*temp/apow;

  if(energy_ctrl->iget_full_inter==1){
   ibin      = (int) ((r12-bond_free->rmin)/(bond_free->del) + 1);
   ibin      = MIN(bond_free->nhist,ibin);
   ibin      = MAX(1,ibin);
   bond_free->hist[ibin] += 1.0;
  }/*endif*/  
  /*------------------------------------------------------------------------*/
  /*  G) Get the force on the atoms using the chain rule                    */
  
  dvbond     = bond_free->fk*pow(rr0,apow1);
  
  pre        = -dvbond/r12;
  intra_scr->fx1[1] = intra_scr->dx12[1]*pre;
  intra_scr->fy1[1] = intra_scr->dy12[1]*pre;
  intra_scr->fz1[1] = intra_scr->dz12[1]*pre;
  /*------------------------------------------------------------------------*/
  /*  J) Pressure tensor, etc.                                              */
  
  if(cell->iperd == 2 || cell->iperd == 3) {
    intra_scr->p11[1] = intra_scr->dx12[1]*intra_scr->fx1[1];
    intra_scr->p22[1] = intra_scr->dy12[1]*intra_scr->fy1[1];
    intra_scr->p33[1] = intra_scr->dz12[1]*intra_scr->fz1[1];
    intra_scr->p12[1] = intra_scr->dx12[1]*intra_scr->fy1[1];
    intra_scr->p21[1] = intra_scr->dy12[1]*intra_scr->fx1[1];
    
    ptens->pvten[1] += intra_scr->p11[1]*wfor*pi_beads;
    ptens->pvten[5] += intra_scr->p22[1]*wfor*pi_beads;
    ptens->pvten[9] += intra_scr->p33[1]*wfor*pi_beads;
    ptens->pvten[2] += intra_scr->p12[1]*wfor*pi_beads;
    ptens->pvten[4] += intra_scr->p21[1]*wfor*pi_beads;
   if(iget_pv_real_inter==1){    
    ptens->pvten_tot[1] += intra_scr->p11[1]*pi_beads;
    ptens->pvten_tot[5] += intra_scr->p22[1]*pi_beads;
    ptens->pvten_tot[9] += intra_scr->p33[1]*pi_beads;
    ptens->pvten_tot[2] += intra_scr->p12[1]*pi_beads;
    ptens->pvten_tot[4] += intra_scr->p21[1]*pi_beads;
   }/*endif*/
  }/*endif*/
  if((cell->iperd) == 3) {
    intra_scr->p13[1] = intra_scr->dx12[1]*intra_scr->fz1[1];
    intra_scr->p31[1] = intra_scr->dz12[1]*intra_scr->fx1[1];
    intra_scr->p23[1] = intra_scr->dy12[1]*intra_scr->fz1[1];
    intra_scr->p32[1] = intra_scr->dz12[1]*intra_scr->fy1[1];
    
    ptens->pvten[3] += intra_scr->p13[1]*wfor*pi_beads;
    ptens->pvten[7] += intra_scr->p31[1]*wfor*pi_beads;
    ptens->pvten[6] += intra_scr->p23[1]*wfor*pi_beads;
    ptens->pvten[8] += intra_scr->p32[1]*wfor*pi_beads;
    
   if(iget_pv_real_inter==1){    
    ptens->pvten_tot[3] += intra_scr->p13[1]*pi_beads;
    ptens->pvten_tot[7] += intra_scr->p31[1]*pi_beads;
    ptens->pvten_tot[6] += intra_scr->p23[1]*pi_beads;
    ptens->pvten_tot[8] += intra_scr->p32[1]*pi_beads;
   }/*endif*/
  }/*endif*/
  /*========================================================================*/
  /*  I) Scatter the forces                                                 */
  
  clatoms_pos->fx[bond_free->j1] +=  intra_scr->fx1[1]*wfor*pi_beads;
  clatoms_pos->fy[bond_free->j1] +=  intra_scr->fy1[1]*wfor*pi_beads;
  clatoms_pos->fz[bond_free->j1] +=  intra_scr->fz1[1]*wfor*pi_beads;
  clatoms_pos->fx[bond_free->j2] += -intra_scr->fx1[1]*wfor*pi_beads;
  clatoms_pos->fy[bond_free->j2] += -intra_scr->fy1[1]*wfor*pi_beads;
  clatoms_pos->fz[bond_free->j2] += -intra_scr->fz1[1]*wfor*pi_beads;
  
  
  /*-----------------------------------------------------------------------*/
}/* end routine */
/*==========================================================================*/



/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void bond_free_mode(CLATOMS_INFO *clatoms_info,CLATOMS_POS *clatoms_pos,
               BOND_FREE *bond_free,CELL *cell,
               INTRA_SCR *intra_scr,PTENS *ptens, double *vbondt,
               ENERGY_CTRL *energy_ctrl,int np_forc)

/*==========================================================================*/
{/*begin routine*/
  /*======================================================================*/
  /*           Local variable declarations                                */
  int ibin,nnow;
  int pi_beads = clatoms_info->pi_beads; 
  double r122,r12;                         /* (r-r_0)^2 and (r-r_0)  */
  double dvbond;                           /* Derivative of bond pot */
  double pre;                              /* Force prefactor        */
  double rr0;                              /* (r-r_0)                */
  double apow,apow1;
  double wfor;
  int iget_pv_real_inter = energy_ctrl->iget_pv_real_inter;

  /*=======================================================================*/
  /* I) Get the two particles */
  
  wfor    = intra_scr->wght_tra_res;
  intra_scr->x1[1] = clatoms_info->xmod[bond_free->j1];
  intra_scr->y1[1] = clatoms_info->ymod[bond_free->j1];
  intra_scr->z1[1] = clatoms_info->zmod[bond_free->j1];
  intra_scr->x2[1] = clatoms_info->xmod[bond_free->j2];
  intra_scr->y2[1] = clatoms_info->ymod[bond_free->j2];
  intra_scr->z2[1] = clatoms_info->zmod[bond_free->j2];
  
  /*------------------------------------------------------------------------*/
  /* III) Calculate the basis vectors r12                                   */
  
  intra_scr->dx12[1] = intra_scr->x1[1] - intra_scr->x2[1];
  intra_scr->dy12[1] = intra_scr->y1[1] - intra_scr->y2[1];
  intra_scr->dz12[1] = intra_scr->z1[1] - intra_scr->z2[1];
  
  /*------------------------------------------------------------------------*/
  /*  E) Periodic boundary conditions                                       */
  
  nnow = 1;
     period(nnow,intra_scr->dx12,intra_scr->dy12,intra_scr->dz12,cell);
  
  /*------------------------------------------------------------------------*/
  /*  F) Get the bonding potential energy and fill histogram                */
  
  r122      = (intra_scr->dx12[1]*intra_scr->dx12[1]
               + intra_scr->dy12[1]*intra_scr->dy12[1]
               + intra_scr->dz12[1]*intra_scr->dz12[1]);
  r12       = sqrt(r122);
  rr0       = r12 - bond_free->eq;
  apow      = (double)bond_free->npow;
  apow1     = apow-1.0;
  *vbondt   = bond_free->fk*pow(rr0,apow)/apow;

  if(energy_ctrl->iget_full_inter==1){
   ibin      = (int) ((r12-bond_free->rmin)/(bond_free->del) + 1);
   ibin      = MIN(bond_free->nhist,ibin);
   ibin      = MAX(1,ibin);
   bond_free->hist[ibin] += 1.0;
  }/*endif*/  
  /*------------------------------------------------------------------------*/
  /*  G) Get the force on the atoms using the chain rule                    */
  
  dvbond     = bond_free->fk*pow(rr0,apow1);
  
  pre        = -dvbond/r12;
  intra_scr->fx1[1] = intra_scr->dx12[1]*pre;
  intra_scr->fy1[1] = intra_scr->dy12[1]*pre;
  intra_scr->fz1[1] = intra_scr->dz12[1]*pre;
  /*------------------------------------------------------------------------*/
  /*  J) Pressure tensor, etc.                                              */
  
  if(cell->iperd == 2 || cell->iperd == 3) {
    intra_scr->p11[1] = intra_scr->dx12[1]*intra_scr->fx1[1];
    intra_scr->p22[1] = intra_scr->dy12[1]*intra_scr->fy1[1];
    intra_scr->p33[1] = intra_scr->dz12[1]*intra_scr->fz1[1];
    intra_scr->p12[1] = intra_scr->dx12[1]*intra_scr->fy1[1];
    intra_scr->p21[1] = intra_scr->dy12[1]*intra_scr->fx1[1];
    
    ptens->pvten[1] += intra_scr->p11[1]*wfor*pi_beads;
    ptens->pvten[5] += intra_scr->p22[1]*wfor*pi_beads;
    ptens->pvten[9] += intra_scr->p33[1]*wfor*pi_beads;
    ptens->pvten[2] += intra_scr->p12[1]*wfor*pi_beads;
    ptens->pvten[4] += intra_scr->p21[1]*wfor*pi_beads;
    
   if(iget_pv_real_inter==1){    
    ptens->pvten_tot[1] += intra_scr->p11[1]*pi_beads;
    ptens->pvten_tot[5] += intra_scr->p22[1]*pi_beads;
    ptens->pvten_tot[9] += intra_scr->p33[1]*pi_beads;
    ptens->pvten_tot[2] += intra_scr->p12[1]*pi_beads;
    ptens->pvten_tot[4] += intra_scr->p21[1]*pi_beads;
   }/*endif*/
  }/*endif*/
  if((cell->iperd) == 3) {
    intra_scr->p13[1] = intra_scr->dx12[1]*intra_scr->fz1[1];
    intra_scr->p31[1] = intra_scr->dz12[1]*intra_scr->fx1[1];
    intra_scr->p23[1] = intra_scr->dy12[1]*intra_scr->fz1[1];
    intra_scr->p32[1] = intra_scr->dz12[1]*intra_scr->fy1[1];
    
    ptens->pvten[3] += intra_scr->p13[1]*wfor*pi_beads;
    ptens->pvten[7] += intra_scr->p31[1]*wfor*pi_beads;
    ptens->pvten[6] += intra_scr->p23[1]*wfor*pi_beads;
    ptens->pvten[8] += intra_scr->p32[1]*wfor*pi_beads;
    
   if(iget_pv_real_inter==1){    
    ptens->pvten_tot[3] += intra_scr->p13[1]*pi_beads;
    ptens->pvten_tot[7] += intra_scr->p31[1]*pi_beads;
    ptens->pvten_tot[6] += intra_scr->p23[1]*pi_beads;
    ptens->pvten_tot[8] += intra_scr->p32[1]*pi_beads;
   }/*endif*/
  }/*endif*/
  /*========================================================================*/
  /*  I) Scatter the forces                                                 */
  
  clatoms_pos->fx[bond_free->j1] +=  intra_scr->fx1[1]*wfor;
  clatoms_pos->fy[bond_free->j1] +=  intra_scr->fy1[1]*wfor;
  clatoms_pos->fz[bond_free->j1] +=  intra_scr->fz1[1]*wfor;
  clatoms_pos->fx[bond_free->j2] += -intra_scr->fx1[1]*wfor;
  clatoms_pos->fy[bond_free->j2] += -intra_scr->fy1[1]*wfor;
  clatoms_pos->fz[bond_free->j2] += -intra_scr->fz1[1]*wfor;
  
  
  /*-----------------------------------------------------------------------*/
}/* end routine */
/*==========================================================================*/
