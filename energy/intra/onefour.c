/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
/*                                                                          */
/*                         PI_MD:                                           */
/*             The future of simulation technology                          */
/*             ------------------------------------                         */
/*                     Module: onefour.c                                    */
/*                                                                          */
/* This routine computes the energy and forces from                         */ 
/* the intramolecular onefour potential.                                    */
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

void onfo(CLATOMS_INFO *clatoms_info,CLATOMS_POS *clatoms_pos,
          ONFO *onfo,CELL *cell,
          INTRA_SCR *intra_scr,PTENS *ptens, double *vonfot,
	  double *vonfo_vdw,double *vonfo_coul, int iver_get,
          CLASS_COMM_FORC_PKG *class_comm_forc_pkg,
          int iget_pv_real_inter)

/*==========================================================================*/
/*               Begin subprogram:                                          */
      {/*begin routine*/
/*==========================================================================*/
/*               Local variable declarations                                */
      int ionfo,ibig,ibig1,ioff,iend,nnow;     /* Indices and counters   */
      double r122i,r12i;                         /* (r-r_0)^2 and (r-r_0)  */
      double vonfo,vonfoe;                       /* Onfo potential         */
      double dvonfo,xxx;                       /* Derivative of onfo pot */
      double pre;                              /* Force prefactor        */
      double q12s;
      double wfor;
      int i,ktemp,iooff;
  int ktemp1,ktemp2,ktemp5,ktemp6,iboff,iii;
  int nlen_use,nlen_now;
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
  double *intra_scr_p13   = intra_scr->p13;
  double *intra_scr_p23   = intra_scr->p23;
  double *clatoms_x       = clatoms_pos->x;
  double *clatoms_y       = clatoms_pos->y;
  double *clatoms_z       = clatoms_pos->z;
  double *clatoms_q       = clatoms_info->q;
  double *clatoms_fx      = clatoms_pos->fx;
  double *clatoms_fy      = clatoms_pos->fy;
  double *clatoms_fz      = clatoms_pos->fz;
  double *clatoms_fxt      = clatoms_pos->fxt;
  double *clatoms_fyt      = clatoms_pos->fyt;
  double *clatoms_fzt      = clatoms_pos->fzt;
  double *onfo_sc         = onfo->sc;
  double *onfo_s6         = onfo->s6;
  double *onfo_feps       = onfo->feps;
  double *ptens_pvten_tmp = ptens->pvten_tmp;
  double *xmod            = clatoms_info->xmod;
  double *ymod            = clatoms_info->ymod;
  double *zmod            = clatoms_info->zmod;
  int *onfo_j1            = onfo->j1;
  int *onfo_j2            = onfo->j2;
  int *onfo_jtyp         = onfo->jtyp;
  double *intra_scr_vpot_e= intra_scr->dx13;
  int *iblock_size        = onfo->iblock_size;
  int *iblock_conflict_1  = onfo->iblock_conflict_1;
  int *iblock_conflict_2  = onfo->iblock_conflict_2;
  int np_forc             = class_comm_forc_pkg->num_proc;
  int myid_forc           = class_comm_forc_pkg->myid;

/*==========================================================================*/
/* I) loop over all the onfos in steps of nlen to save memory               */


      wfor    = (intra_scr->wght_tra_res);
      for(i=1;i<=9;i++){(ptens_pvten_tmp)[i]=0;}


  nlen_use = (intra_scr->nlen);
  ntot = (onfo->num);

  ntot = onfo->num;

for(ibig=1;ibig <= ntot;ibig += nlen_use) {

/*--------------------------------------------------------------------------*/
/*  A) Offsets to save some memory by doing only nlen onfos at a time       */

    ibig1 = ibig-1;
    ioff = -ibig1;
    nlen_now = nlen_use;
    iend = MIN(ntot,ibig1+nlen_now);
    nnow = iend-ibig1;

/*--------------------------------------------------------------------------*/
/*  B) Gather positions and charges                                         */

       for(iooff=1,ionfo=ibig;ionfo <= iend; ++ionfo,++iooff) {       
         ktemp1 = onfo_j1[ionfo];
         ktemp2 = onfo_j2[ionfo];
         intra_scr_dx12[iooff] = clatoms_x[ktemp1]-clatoms_x[ktemp2];
         intra_scr_dy12[iooff] = clatoms_y[ktemp1]-clatoms_y[ktemp2];
         intra_scr_dz12[iooff] = clatoms_z[ktemp1]-clatoms_z[ktemp2];
       }/*endfor*/

       if(cell->intra_perds==1&&clatoms_info->pi_beads>1){     
         for(iooff=1,ionfo=ibig;ionfo <= iend; ++ionfo,++iooff) {       
           ktemp5 = onfo_j1[ionfo];
           ktemp6 = onfo_j2[ionfo];
           intra_scr_dx56[iooff] = xmod[ktemp5]-xmod[ktemp6];
           intra_scr_dy56[iooff] = ymod[ktemp5]-ymod[ktemp6];
           intra_scr_dz56[iooff] = zmod[ktemp5]-zmod[ktemp6];
         }/*endfor*/
       }/*endif*/

/*--------------------------------------------------------------------------*/
/*  E) Periodic boundary conditions                                         */

     if(cell->intra_perds == 1) {
       if(clatoms_info->pi_beads==1){     
         period(nnow,intra_scr_dx12,intra_scr_dy12,intra_scr_dz12,cell);
       }else{
         period_pimd(nnow,intra_scr_dx12,intra_scr_dy12,intra_scr_dz12,
                           intra_scr_dx56,intra_scr_dy56,intra_scr_dz56,cell);
       }/*endif*/
     }/*endif*/

       for(ionfo=1,iboff=ibig;ionfo <= nnow; ++ionfo,++iboff) {
         r122i = 1.0/(intra_scr_dx12[ionfo]*intra_scr_dx12[ionfo]
                    + intra_scr_dy12[ionfo]*intra_scr_dy12[ionfo]
                    + intra_scr_dz12[ionfo]*intra_scr_dz12[ionfo]);
         r12i  = sqrt(r122i);

         ktemp = onfo_jtyp[iboff];
/*--------------------------------------------------------------------------*/
/*  F) Get the onfoing potential energy                                     */

         q12s   = clatoms_q[onfo_j1[iboff]]*clatoms_q[onfo_j2[iboff]]
                 *onfo_sc[ktemp];
         xxx    = onfo_s6[ktemp]*r122i*r122i*r122i;
         vonfo  = onfo_feps[ktemp]*xxx*(xxx-1.0);
         vonfoe = q12s*r12i;

         (intra_scr_vpot)[ionfo] = vonfo+vonfoe;
	 intra_scr_vpot_e[ionfo] = vonfoe;

/*--------------------------------------------------------------------------*/
/*  G) Get the force on the atoms using the chain rule                      */

         dvonfo = 6.0*onfo_feps[ktemp]*xxx*(2.0*xxx-1.0)*r122i
                + q12s*r12i*r122i;
         pre    = dvonfo;

         intra_scr_fx1[ionfo] = intra_scr_dx12[ionfo]*pre;
         intra_scr_fy1[ionfo] = intra_scr_dy12[ionfo]*pre;
         intra_scr_fz1[ionfo] = intra_scr_dz12[ionfo]*pre;

       /*endfor ionfo*/}

/*--------------------------------------------------------------------------*/
/*  H) Sum the potential energy                                             */

       for(ionfo=1;ionfo <= nnow; ++ionfo) {
	 *vonfot     += intra_scr_vpot[ionfo];
	 *vonfo_coul += intra_scr_vpot_e[ionfo];
       }

/*--------------------------------------------------------------------------*/
/*  J) Pressure tensor, etc.                                                */

       if(cell->iperd == 2 || cell->iperd == 3) {
          for(ionfo=1;ionfo <= nnow; ++ionfo) {
            intra_scr_p11[ionfo] = intra_scr_dx12[ionfo]
                                     *intra_scr_fx1[ionfo];
            intra_scr_p22[ionfo] = intra_scr_dy12[ionfo]
                                     *intra_scr_fy1[ionfo];
            intra_scr_p33[ionfo] = intra_scr_dz12[ionfo]
                                     *intra_scr_fz1[ionfo];
            intra_scr_p12[ionfo] = intra_scr_dx12[ionfo]
                                     *intra_scr_fy1[ionfo];
            intra_scr_p21[ionfo] = intra_scr_dy12[ionfo]
                                     *intra_scr_fx1[ionfo];
          /*endfor*/}
          for(ionfo=1;ionfo <= nnow; ++ionfo) {
            ptens_pvten_tmp[1] += intra_scr_p11[ionfo];
            ptens_pvten_tmp[5] += intra_scr_p22[ionfo];
            ptens_pvten_tmp[9] += intra_scr_p33[ionfo];
            ptens_pvten_tmp[2] += intra_scr_p12[ionfo];
          /*endfor*/}
       /*endif*/}
       if((cell->iperd) == 3) {
          for(ionfo=1;ionfo <= nnow; ++ionfo) {
            intra_scr_p13[ionfo] = intra_scr_dx12[ionfo]
                                  *intra_scr_fz1[ionfo];
            intra_scr_p23[ionfo] = intra_scr_dy12[ionfo]
                                  *intra_scr_fz1[ionfo];
          /*endfor*/}
          for(ionfo=1;ionfo <= nnow; ++ionfo) {
             ptens_pvten_tmp[3] += intra_scr_p13[ionfo];
             ptens_pvten_tmp[6] += intra_scr_p23[ionfo];
          /*endfor*/}
       /*endif*/}
/*--------------------------------------------------------------------------*/
/*  I) Scatter the forces                                                   */
      if(iver_get==1){

       for(iooff=1,ionfo=ibig;ionfo <= iend; ++ionfo,++iooff) {       
         ktemp = onfo_j1[ionfo];
         clatoms_fxt[ktemp] += intra_scr_fx1[iooff];
         clatoms_fyt[ktemp] += intra_scr_fy1[iooff];
         clatoms_fzt[ktemp] += intra_scr_fz1[iooff];
       }

       for(iooff=1,ionfo=ibig;ionfo <= iend; ++ionfo,++iooff) {       
         ktemp = onfo_j2[ionfo];
         clatoms_fxt[ktemp] -= intra_scr_fx1[iooff];
         clatoms_fyt[ktemp] -= intra_scr_fy1[iooff];
         clatoms_fzt[ktemp] -= intra_scr_fz1[iooff];
       /*endfor*/}
    }/*endif*/

       for(iooff=1;iooff <= iend+ioff; ++iooff) {
          intra_scr_fx1[iooff] *= wfor;
          intra_scr_fy1[iooff] *= wfor;
          intra_scr_fz1[iooff] *= wfor;
       /*endfor*/}

       for(iooff=1,ionfo=ibig;ionfo <= iend; ++ionfo,++iooff) {       
         ktemp = onfo_j1[ionfo];
         clatoms_fx[ktemp] += intra_scr_fx1[iooff];
         clatoms_fy[ktemp] += intra_scr_fy1[iooff];
         clatoms_fz[ktemp] += intra_scr_fz1[iooff];
       }

       for(iooff=1,ionfo=ibig;ionfo <= iend; ++ionfo,++iooff) {       
         ktemp = onfo_j2[ionfo];
         clatoms_fx[ktemp] -= intra_scr_fx1[iooff];
         clatoms_fy[ktemp] -= intra_scr_fy1[iooff];
         clatoms_fz[ktemp] -= intra_scr_fz1[iooff];
       }/*endfor*/

  } /* endfor ibig */
/*==========================================================================*/
/* Increment the Pressure tensor */

   ptens_pvten_tmp[4] = ptens_pvten_tmp[2];
   ptens_pvten_tmp[7] = ptens_pvten_tmp[3];
   ptens_pvten_tmp[8] = ptens_pvten_tmp[6];
   for(i=1;i<=9;i++){
       (ptens->pvten)[i]     += (ptens_pvten_tmp)[i]*wfor;
   }/*endfor*/

   if(iget_pv_real_inter==1){    
     for(i=1;i<=9;i++){
       (ptens->pvten_tot)[i] += (ptens_pvten_tmp)[i];
     }/*endfor*/
   }/*endif*/

   *vonfo_vdw = *vonfot - *vonfo_coul;

/*--------------------------------------------------------------------------*/
/* end routine */}
/*==========================================================================*/
