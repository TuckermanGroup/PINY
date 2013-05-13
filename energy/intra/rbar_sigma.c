/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
/*                                                                          */
/*                         PI_MD:                                           */
/*             The future of simulation technology                          */
/*             ------------------------------------                         */
/*                     Module: rbar_sigma.c                                 */
/*                                                                          */
/* This routine computes the energy and forces from                         */ 
/* the intramolecular biasing potential.                                    */
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

void rbar_sig_free(CLATOMS_INFO *clatoms_info,CLATOMS_POS *clatoms_pos,
                   RBAR_SIG_FREE *rbar_sig_free,CELL *cell,
                   INTRA_SCR *intra_scr,PTENS *ptens, double *vbar_free, 
                   ENERGY_CTRL *energy_ctrl,int np_forc)

/*==========================================================================*/
   {/*begin routine*/
/*=======================================================================*/
/*            Local variable declarations                                */

  int ibond,ibig,ibig1,ioff,iend,nnow;     /* Indices and counters   */
  double r122;                             /* (r-r_0)^2 and (r-r_0)  */
  double pre;                              /* Force prefactor        */
  int iii;
  int ibin_bar,ibin_sig,ibin;
  double const1,const3,const4,r_bar,sigma,rr0_bar;
  double rr0_sig;
  double r2_bar;

 /* Define local pointers    */

  double *intra_scr_x1    = intra_scr->x1;
  double *intra_scr_y1    = intra_scr->y1;
  double *intra_scr_z1    = intra_scr->z1;
  double *intra_scr_x2    = intra_scr->x2;
  double *intra_scr_y2    = intra_scr->y2;
  double *intra_scr_z2    = intra_scr->z2;
  double *intra_scr_fx1   = intra_scr->fx1;
  double *intra_scr_fy1   = intra_scr->fy1;
  double *intra_scr_fz1   = intra_scr->fz1;
  double *intra_scr_dx12  = intra_scr->dx12;
  double *intra_scr_dy12  = intra_scr->dy12;
  double *intra_scr_dz12  = intra_scr->dz12;
  double *intra_scr_vpot  = intra_scr->vpot;
  double *intra_scr_eq    = intra_scr->eq;
  double *intra_scr_p11   = intra_scr->p11;
  double *intra_scr_p22   = intra_scr->p22;
  double *intra_scr_p33   = intra_scr->p33;
  double *intra_scr_p21   = intra_scr->p21;
  double *intra_scr_p12   = intra_scr->p12;
  double *intra_scr_p31   = intra_scr->p31;
  double *intra_scr_p13   = intra_scr->p13;
  double *intra_scr_p32   = intra_scr->p32;
  double *intra_scr_p23   = intra_scr->p23;
  double *intra_scr_r     = intra_scr->r;
  double *intra_scr_r2    = intra_scr->r3;

  int pi_beads            = clatoms_info->pi_beads;
  double *clatoms_x       = clatoms_pos->x;
  double *clatoms_y       = clatoms_pos->y;
  double *clatoms_z       = clatoms_pos->z;
  double *clatoms_pos_fx  = clatoms_pos->fx;
  double *clatoms_pos_fy  = clatoms_pos->fy;
  double *clatoms_pos_fz  = clatoms_pos->fz;
  double *clatoms_fxt     = clatoms_pos->fxt;
  double *clatoms_fyt     = clatoms_pos->fyt;
  double *clatoms_fzt     = clatoms_pos->fzt;
  double *ptens_pvten_tmp = ptens->pvten_tmp;
  double *ptens_pvten     = ptens->pvten;
  double *ptens_pvten_tot = ptens->pvten_tot;

  int nfree               = rbar_sig_free->nfree;
  double rnfree           = rbar_sig_free->rnfree;
  double eq_bar           = rbar_sig_free->eq_bar;
  double eq_sigma         = rbar_sig_free->eq_sigma;
  double fk_bar           = rbar_sig_free->fk_bar;
  double fk_sigma         = rbar_sig_free->fk_sigma;
  int nhist_bar           = rbar_sig_free->nhist_bar;
  int nhist_sig           = rbar_sig_free->nhist_sig;
  double smin             = rbar_sig_free->smin;
  double smax             = rbar_sig_free->smax;
  double rmin             = rbar_sig_free->rmin;
  double rmax             = rbar_sig_free->rmax;
  double del_bar          = rbar_sig_free->del_bar;
  double del_sig          = rbar_sig_free->del_sig;
  double **hist           = rbar_sig_free->hist;
  double **hist_rn      = rbar_sig_free->hist_rn;
  int *j1                 = rbar_sig_free->j1;
  int *j2                 = rbar_sig_free->j2;

  int energy_ctrl_iget_full_inter = energy_ctrl->iget_full_inter;
  int intra_perds                 = cell->intra_perds;
  int iperd                       = cell->iperd;
  double wfor                     = intra_scr->wght_tra_res;

/*=======================================================================*/
/* I) Get the pairs of particles */
 
  
   for(ibond=1;ibond<=nfree;ibond++){
      intra_scr_x1[ibond] = clatoms_x[j1[ibond]];
      intra_scr_y1[ibond] = clatoms_y[j1[ibond]];
      intra_scr_z1[ibond] = clatoms_z[j1[ibond]];
      intra_scr_x2[ibond] = clatoms_x[j2[ibond]];
      intra_scr_y2[ibond] = clatoms_y[j2[ibond]];
      intra_scr_z2[ibond] = clatoms_z[j2[ibond]];
   }/*endfor*/

/*=======================================================================*/
/* II) Calculate the basis vectors                                       */
  
   for(ibond=1;ibond<=nfree;ibond++){
      intra_scr_dx12[ibond] = intra_scr_x1[ibond] - intra_scr_x2[ibond];
      intra_scr_dy12[ibond] = intra_scr_y1[ibond] - intra_scr_y2[ibond];
      intra_scr_dz12[ibond] = intra_scr_z1[ibond] - intra_scr_z2[ibond];
   }/*endfor*/
  
/*=======================================================================*/
/* III) Periodic boundary conditions                                     */
  
   if(intra_perds == 1) {
     period(nfree,intra_scr_dx12,intra_scr_dy12,intra_scr_dz12,cell);
   }/*endif*/

/*=======================================================================*/
/* IV) Get the biasing potential energy                                  */
  
   r_bar = 0.0;
   r2_bar = 0.0;
   for(ibond=1;ibond<=nfree;ibond++){
      r122                = (intra_scr_dx12[ibond]*intra_scr_dx12[ibond]
                            +intra_scr_dy12[ibond]*intra_scr_dy12[ibond]
                            +intra_scr_dz12[ibond]*intra_scr_dz12[ibond]);
      intra_scr_r[ibond]  = sqrt(r122);
      r_bar              += intra_scr_r[ibond];
      r2_bar             += r122;
   }/*endfor*/
   r_bar     /= rnfree;
   r2_bar    /= rnfree;
   sigma      = sqrt(fabs(r2_bar-r_bar*r_bar));

   rr0_bar     = (r_bar - eq_bar);
   *vbar_free  = 0.5*fk_bar*rr0_bar*rr0_bar;

   rr0_sig     = (sigma - eq_sigma);
   *vbar_free += 0.5*fk_sigma*rr0_sig*rr0_sig;

/*=======================================================================*/
/* V) Fill the histogram                                                 */

   if(energy_ctrl_iget_full_inter==1){

     ibin_sig    = (sigma-smin)/del_sig + 1.0;
     ibin_sig    = MIN(nhist_sig,ibin_sig);
     ibin_sig    = MAX(1,ibin_sig);

     ibin_bar    = (r_bar-rmin)/del_bar + 1.0;
     ibin_bar    = MIN(nhist_bar,ibin_bar);
     ibin_bar    = MAX(1,ibin_bar);

     hist[ibin_bar][ibin_sig] += 1.0;

     for(ibond=1;ibond<=nfree;ibond++){
        ibin   = (intra_scr_r[ibond]-rmin)/del_bar + 1.0;
        ibin   = MIN(nhist_bar,ibin);
        ibin   = MAX(1,ibin_bar);
        hist_rn[ibond][ibin] += 1.0;
     }/*endfor*/

   }/*endif*/  

/*=======================================================================*/  
/*  VI) Get the force on the atoms using the chain rule                  */
  
   const3 = -fk_bar*(r_bar-eq_bar)/rnfree;
   const4 = -fk_sigma*(sigma-eq_sigma)/(sigma*rnfree);

   for(ibond=1;ibond<=nfree;ibond++){

      const1 = (intra_scr_r[ibond]-r_bar)*const4;
      pre    = (const1+const3)/intra_scr_r[ibond];

      intra_scr_fx1[ibond] = intra_scr_dx12[ibond]*pre;
      intra_scr_fy1[ibond] = intra_scr_dy12[ibond]*pre;
      intra_scr_fz1[ibond] = intra_scr_dz12[ibond]*pre;

     }/*endfor*/

/*=======================================================================*/  
/*  VII) Pressure tensor                                                 */
  
  if(iperd == 2 || iperd == 3) {

    for(ibond=1;ibond<=nfree;ibond++){
      intra_scr_p11[ibond] = intra_scr_dx12[ibond]*intra_scr_fx1[ibond];
      intra_scr_p22[ibond] = intra_scr_dy12[ibond]*intra_scr_fy1[ibond];
      intra_scr_p33[ibond] = intra_scr_dz12[ibond]*intra_scr_fz1[ibond];
      intra_scr_p12[ibond] = intra_scr_dx12[ibond]*intra_scr_fy1[ibond];
      intra_scr_p21[ibond] = intra_scr_dy12[ibond]*intra_scr_fx1[ibond];
     }/*endfor*/
    
     for(ibond=1;ibond<=nfree;ibond++){
       ptens_pvten[1] += intra_scr_p11[ibond]*wfor;
       ptens_pvten[5] += intra_scr_p22[ibond]*wfor;
       ptens_pvten[9] += intra_scr_p33[ibond]*wfor;
       ptens_pvten[2] += intra_scr_p12[ibond]*wfor;
       ptens_pvten[4] += intra_scr_p21[ibond]*wfor;
    
       ptens_pvten_tot[1] += intra_scr_p11[ibond];
       ptens_pvten_tot[5] += intra_scr_p22[ibond];
       ptens_pvten_tot[9] += intra_scr_p33[ibond];
       ptens_pvten_tot[2] += intra_scr_p12[ibond];
       ptens_pvten_tot[4] += intra_scr_p21[ibond];
     }/*endfor*/

   }/*endif*/

   if(iperd == 3) {

     for(ibond=1;ibond<=nfree;ibond++){
       intra_scr_p13[ibond] = intra_scr_dx12[ibond]*intra_scr_fz1[ibond];
       intra_scr_p31[ibond] = intra_scr_dz12[ibond]*intra_scr_fx1[ibond];
       intra_scr_p23[ibond] = intra_scr_dy12[ibond]*intra_scr_fz1[ibond];
       intra_scr_p32[ibond] = intra_scr_dz12[ibond]*intra_scr_fy1[ibond];
     }/*endfor*/
    
     for(ibond=1;ibond<=nfree;ibond++){
       ptens_pvten[3] += intra_scr_p13[ibond]*wfor;
       ptens_pvten[7] += intra_scr_p31[ibond]*wfor;
       ptens_pvten[6] += intra_scr_p23[ibond]*wfor;
       ptens_pvten[8] += intra_scr_p32[ibond]*wfor;
    
       ptens_pvten_tot[3] += intra_scr_p13[ibond];
       ptens_pvten_tot[7] += intra_scr_p31[ibond];
       ptens_pvten_tot[6] += intra_scr_p23[ibond];
       ptens_pvten_tot[8] += intra_scr_p32[ibond];
     }/*endfor*/

   }/*endif*/

/*========================================================================*/
/*  VIII) Scatter the forces                                              */
  
   for(ibond=1;ibond<=nfree;ibond++){
     clatoms_pos_fx[j1[ibond]] +=  intra_scr_fx1[ibond]*wfor;
     clatoms_pos_fy[j1[ibond]] +=  intra_scr_fy1[ibond]*wfor;
     clatoms_pos_fz[j1[ibond]] +=  intra_scr_fz1[ibond]*wfor;
     clatoms_pos_fx[j2[ibond]] += -intra_scr_fx1[ibond]*wfor;
     clatoms_pos_fy[j2[ibond]] += -intra_scr_fy1[ibond]*wfor;
     clatoms_pos_fz[j2[ibond]] += -intra_scr_fz1[ibond]*wfor;
   }/*endfor*/
  
/*--------------------------------------------------------------------------*/
     }/* end routine */
/*==========================================================================*/

