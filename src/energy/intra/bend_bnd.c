/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
/*                                                                          */
/*                         PI_MD:                                           */
/*             The future of simulation technology                          */
/*             ------------------------------------                         */
/*                     Module: bend_bnd.c                                   */
/*                                                                          */
/* This routine computes the energy and forces from                         */
/* the intramolecular bend_bnd potential.                                   */
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

void bend_bnd(CLATOMS_INFO *clatoms_info,CLATOMS_POS *clatoms_pos,
              BEND_BND *bend_bnd,CELL *cell,
              INTRA_SCR *intra_scr,PTENS *ptens,double *vbend_bndt,
              double *vbend_bnd_bend,double *vbend_bnd_bond, int iver_get,
              CLASS_COMM_FORC_PKG *class_comm_forc_pkg,
              int iget_pv_real_inter)

/*==========================================================================*/
/*   Begin Routine                                                          */
      {/*begin routine*/
/*=======================================================================*/
/*            Local variable declarations                                */

  int iii;
  int ibend_bnd,ibig,ibig1,ioff,iend,nnow; /* Indices and counters       */
  int nlen_use,nlen_now;
  int ntot;
  int ngo,irem,istart_big;

  double r122,r12;                         /* (r-r_0)^2 and (r-r_0)      */
  double r322,r32;
  double cost,sint,seps;
  double vbendc,vbends;            /* Bend_bnd potential         */
  double dvbendc,dvbends;          /* Derivative of bend_bnd pot */
  double vbond, dvbond;
  double sisum;                             /* sine sum                   */
  double pre;                              /* Force prefactor            */
  double rpmag;                            /* 1/(r12*r32)                */
  double cos122,cos322;                    /* cos/r122, cos/r322         */
  double rr0;
  double wfor;
  int i,ktemp,iboff;
  double fxtmp,fytmp,fztmp;
  double x1,y1,z1,x2,y2,z2,x3,y3,z3,x5,y5,z5,x6,y6,z6;
  double x7,y7,z7;
  int ktemp1,ktemp2,ktemp3,ktemp5,ktemp6,ktemp7;

  /* Define local pointers                                                 */
  double *intra_scr_fx1    = intra_scr->fx1;
  double *intra_scr_fy1    = intra_scr->fy1;
  double *intra_scr_fz1    = intra_scr->fz1;
  double *intra_scr_fx2    = intra_scr->fx2;
  double *intra_scr_fy2    = intra_scr->fy2;
  double *intra_scr_fz2    = intra_scr->fz2;
  double *intra_scr_fx3    = intra_scr->fx3;
  double *intra_scr_fy3    = intra_scr->fy3;
  double *intra_scr_fz3    = intra_scr->fz3;
  double *intra_scr_dx12   = intra_scr->dx12;
  double *intra_scr_dy12   = intra_scr->dy12;
  double *intra_scr_dz12   = intra_scr->dz12;
  double *intra_scr_dx23   = intra_scr->dx23;
  double *intra_scr_dy23   = intra_scr->dy23;
  double *intra_scr_dz23   = intra_scr->dz23;
  double *intra_scr_dx43   = intra_scr->dx43;
  double *intra_scr_dy43   = intra_scr->dy43;
  double *intra_scr_dz43   = intra_scr->dz43;
  double *intra_scr_dx56   = intra_scr->dx56;
  double *intra_scr_dy56   = intra_scr->dy56;
  double *intra_scr_dz56   = intra_scr->dz56;
  double *intra_scr_dx67   = intra_scr->dx67;
  double *intra_scr_dy67   = intra_scr->dy67;
  double *intra_scr_dz67   = intra_scr->dz67;
  double *intra_scr_dx87   = intra_scr->dx87;
  double *intra_scr_dy87   = intra_scr->dy87;
  double *intra_scr_dz87   = intra_scr->dz87;
  double *intra_scr_vpot   = intra_scr->vpot;
  double *intra_scr_p11    = intra_scr->p11;
  double *intra_scr_p22    = intra_scr->p22;
  double *intra_scr_p33    = intra_scr->p33;
  double *intra_scr_p12    = intra_scr->p12;
  double *intra_scr_p13    = intra_scr->p13;
  double *intra_scr_p23    = intra_scr->p23;
  double *bend_bnd_eq_bond = bend_bnd->eq_bond;
  double *bend_bnd_cbend_0 = bend_bnd->cbend_0;
  double *bend_bnd_cbend_1 = bend_bnd->cbend_1;
  double *bend_bnd_cbend_2 = bend_bnd->cbend_2;
  double *bend_bnd_cbend_3 = bend_bnd->cbend_3;
  double *bend_bnd_cbend_4 = bend_bnd->cbend_4;
  double *bend_bnd_cbend_5 = bend_bnd->cbend_5;
  double *bend_bnd_cbend_6 = bend_bnd->cbend_6;
  double *bend_bnd_sbend_0 = bend_bnd->sbend_0;
  double *bend_bnd_sbend_1 = bend_bnd->sbend_1;
  double *bend_bnd_sbend_2 = bend_bnd->sbend_2;
  double *bend_bnd_sbend_3 = bend_bnd->sbend_3;
  double *bend_bnd_sbend_4 = bend_bnd->sbend_4;
  double *bend_bnd_sbend_5 = bend_bnd->sbend_5;
  double *bend_bnd_sbend_6 = bend_bnd->sbend_6;
  double *bend_bnd_dcbend_0= bend_bnd->dcbend_0;
  double *bend_bnd_dcbend_1= bend_bnd->dcbend_1;
  double *bend_bnd_dcbend_2= bend_bnd->dcbend_2;
  double *bend_bnd_dcbend_3= bend_bnd->dcbend_3;
  double *bend_bnd_dcbend_4= bend_bnd->dcbend_4;
  double *bend_bnd_dcbend_5= bend_bnd->dcbend_5;
  double *bend_bnd_dcbend_6= bend_bnd->dcbend_6;
  double *bend_bnd_dsbend_0= bend_bnd->dsbend_0;
  double *bend_bnd_dsbend_1= bend_bnd->dsbend_1;
  double *bend_bnd_dsbend_2= bend_bnd->dsbend_2;
  double *bend_bnd_dsbend_3= bend_bnd->dsbend_3;
  double *bend_bnd_dsbend_4= bend_bnd->dsbend_4;
  double *bend_bnd_dsbend_5= bend_bnd->dsbend_5;
  double *bend_bnd_dsbend_6= bend_bnd->dsbend_6;
  double *bend_bnd_cbond_0 = bend_bnd->cbond_0;
  double *bend_bnd_cbond_1 = bend_bnd->cbond_1;
  double *bend_bnd_cbond_2 = bend_bnd->cbond_2;
  double *bend_bnd_cbond_3 = bend_bnd->cbond_3;
  double *bend_bnd_cbond_4 = bend_bnd->cbond_4;
  double *bend_bnd_cbond_5 = bend_bnd->cbond_5;
  double *bend_bnd_cbond_6 = bend_bnd->cbond_6;
  double *bend_bnd_dcbond_0 = bend_bnd->dcbond_0;
  double *bend_bnd_dcbond_1 = bend_bnd->dcbond_1;
  double *bend_bnd_dcbond_2 = bend_bnd->dcbond_2;
  double *bend_bnd_dcbond_3 = bend_bnd->dcbond_3;
  double *bend_bnd_dcbond_4 = bend_bnd->dcbond_4;
  double *bend_bnd_dcbond_5 = bend_bnd->dcbond_5;
  double *bend_bnd_dcbond_6 = bend_bnd->dcbond_6;
  double *clatoms_x         = clatoms_pos->x;
  double *clatoms_y         = clatoms_pos->y;
  double *clatoms_z         = clatoms_pos->z;
  double *clatoms_fx        = clatoms_pos->fx;
  double *clatoms_fy        = clatoms_pos->fy;
  double *clatoms_fz        = clatoms_pos->fz;
  double *clatoms_fxt        = clatoms_pos->fxt;
  double *clatoms_fyt        = clatoms_pos->fyt;
  double *clatoms_fzt        = clatoms_pos->fzt;
  double *ptens_pvten_tmp   = ptens->pvten_tmp;
  double *xmod            = clatoms_info->xmod;
  double *ymod            = clatoms_info->ymod;
  double *zmod            = clatoms_info->zmod;
  int *bend_bnd_j1          = bend_bnd->j1;
  int *bend_bnd_j2          = bend_bnd->j2;
  int *bend_bnd_j3          = bend_bnd->j3;
  int *bend_bnd_jtyp        = bend_bnd->jtyp;
  int *iblock_size        = bend_bnd->iblock_size;
  int *iblock_conflict_1  = bend_bnd->iblock_conflict_1;
  int *iblock_conflict_2  = bend_bnd->iblock_conflict_2;
  int *iblock_conflict_3  = bend_bnd->iblock_conflict_3;
  int np_forc             = class_comm_forc_pkg->num_proc;
  int myid_forc           = class_comm_forc_pkg->myid;

/*=======================================================================*/
/* 0) Initialize                                                         */

  seps = 1.0e-8;
  wfor    = (intra_scr->wght_tra_res);
  for(i=1;i<=9;i++){(ptens_pvten_tmp)[i]=0;}

/*=======================================================================*/
/* LOOP OVER ALL THE BEND_BNDS IN STEPS OF NLEN TO SAVE MEMORY            */


  nlen_use = (intra_scr->nlen);
  ntot = (bend_bnd->num);

  for(ibig=1;ibig <= ntot;ibig += nlen_use) {
/*=======================================================================*/
/*  I) Offsets to save some memory by doing only nlen bend_bnds at a time   */

    ibig1 = ibig-1;
    ioff = -ibig1;
    nlen_now = nlen_use;
    iend = MIN(ntot,(ibig1+nlen_now));
    nnow = iend-ibig1;

/*=======================================================================*/
/*  II) Gather positions                                                */
    
    for(iboff=ibig+ioff,ibend_bnd=ibig;ibend_bnd <= iend; 
                                                ++ibend_bnd,++iboff){
      ktemp1 = bend_bnd_j1[ibend_bnd];
      ktemp2 = bend_bnd_j2[ibend_bnd];
      ktemp3 = bend_bnd_j3[ibend_bnd];
      x1 = clatoms_x[ktemp1];
      y1 = clatoms_y[ktemp1];
      z1 = clatoms_z[ktemp1];
      x2 = clatoms_x[ktemp2];
      y2 = clatoms_y[ktemp2];
      z2 = clatoms_z[ktemp2];
      x3 = clatoms_x[ktemp3];
      y3 = clatoms_y[ktemp3];
      z3 = clatoms_z[ktemp3];
      intra_scr_dx12[iboff] = x1 - x2;
      intra_scr_dy12[iboff] = y1 - y2;
      intra_scr_dz12[iboff] = z1 - z2;
      intra_scr_dx23[iboff] = x3 - x2;
      intra_scr_dy23[iboff] = y3 - y2;
      intra_scr_dz23[iboff] = z3 - z2;
      intra_scr_dx43[iboff] = x1 - x3;
      intra_scr_dy43[iboff] = y1 - y3;
      intra_scr_dz43[iboff] = z1 - z3;
    }/*endfor*/

   if(cell->intra_perds==1&&clatoms_info->pi_beads>1){     
    for(iboff=ibig+ioff,ibend_bnd=ibig;ibend_bnd <= iend; ++ibend_bnd,++iboff){
        ktemp5 = bend_bnd_j1[ibend_bnd];
        ktemp6 = bend_bnd_j2[ibend_bnd];
        ktemp7 = bend_bnd_j3[ibend_bnd];
        x5 = xmod[ktemp5];
        y5 = ymod[ktemp5];
        z5 = zmod[ktemp5];
        x6 = xmod[ktemp6];
        y6 = ymod[ktemp6];
        z6 = zmod[ktemp6];
        x7 = xmod[ktemp7];
        y7 = ymod[ktemp7];
        z7 = zmod[ktemp7];
        intra_scr_dx56[iboff] = x5 - x6;
        intra_scr_dy56[iboff] = y5 - y6;
        intra_scr_dz56[iboff] = z5 - z6;
        intra_scr_dx67[iboff] = x7 - x6;
        intra_scr_dy67[iboff] = y7 - y6;
        intra_scr_dz67[iboff] = z7 - z6;
        intra_scr_dx87[iboff] = x5 - x7;
        intra_scr_dy87[iboff] = y5 - y7;
        intra_scr_dz87[iboff] = z5 - z7;
      }/*endfor*/
    }/*endif*/     

/*=======================================================================*/
/* IV) Periodic boundary conditions                                      */
    
   if(cell->intra_perds == 1) {
     if(clatoms_info->pi_beads==1){     
       period(nnow,intra_scr_dx12,intra_scr_dy12,intra_scr_dz12,cell);
       period(nnow,intra_scr_dx23,intra_scr_dy23,intra_scr_dz23,cell);
       period(nnow,intra_scr_dx43,intra_scr_dy43,intra_scr_dz43,cell);
     }else{
       period_pimd(nnow,intra_scr_dx12,intra_scr_dy12,intra_scr_dz12,
                         intra_scr_dx56,intra_scr_dy56,intra_scr_dz56,cell);
       period_pimd(nnow,intra_scr_dx23,intra_scr_dy23,intra_scr_dz23,
                         intra_scr_dx67,intra_scr_dy67,intra_scr_dz67,cell);
       period_pimd(nnow,intra_scr_dx43,intra_scr_dy43,intra_scr_dz43,
                         intra_scr_dx87,intra_scr_dy87,intra_scr_dz87,cell);
     }/*endif*/
   }/*endif*/

/*=======================================================================*/
/* VI) Get the bend energies                                             */

    for(ibend_bnd=1,iboff=ibig;ibend_bnd <= nnow; ++ibend_bnd,++iboff) {

      r122  = (intra_scr_dx12[ibend_bnd]*intra_scr_dx12[ibend_bnd]
               + intra_scr_dy12[ibend_bnd]*intra_scr_dy12[ibend_bnd]
               + intra_scr_dz12[ibend_bnd]*intra_scr_dz12[ibend_bnd]);
      r12    = sqrt(r122);
      
      r322   = (intra_scr_dx23[ibend_bnd]*intra_scr_dx23[ibend_bnd]
                + intra_scr_dy23[ibend_bnd]*intra_scr_dy23[ibend_bnd]
                + intra_scr_dz23[ibend_bnd]*intra_scr_dz23[ibend_bnd]);
      r32    = sqrt(r322);
      rpmag  = 1.0/(r12*r32);
      
      /*-------------------------------------------------------------------*/
      /*  A) Calculate the cosine and sine of the angle between them       */
      /*     (cos(theta_{123}), sine(theta_{123}))                         */
      
      cost = (intra_scr_dx12[ibend_bnd]*intra_scr_dx23[ibend_bnd]
              +  intra_scr_dy12[ibend_bnd]*intra_scr_dy23[ibend_bnd]
              +  intra_scr_dz12[ibend_bnd]*intra_scr_dz23[ibend_bnd])
              /(r12*r32);
      
      cost   = (cost < 1.0 ? cost:1.0);
      cost   = (cost > -1.0 ? cost:-1.0);
      
      sint   = sqrt(1.0 - cost*cost);
      sint   = (sint > seps ? sint:seps);
      
      cos122 = cost/r122;
      cos322  = cost/r322;
      

      ktemp = bend_bnd_jtyp[iboff];

      /*--------------------------------------------------------------------*/
      /*  B) Get the bend_bnding potential energy via cosine and sine power */
      /*     series using Horner's method                                   */
      
      vbendc = ((((((bend_bnd_cbend_6[ktemp]
                     *cost + bend_bnd_cbend_5[ktemp])
                    *cost + bend_bnd_cbend_4[ktemp])
                   *cost + bend_bnd_cbend_3[ktemp])
                  *cost + bend_bnd_cbend_2[ktemp])
                 *cost + bend_bnd_cbend_1[ktemp])
                *cost + bend_bnd_cbend_0[ktemp]);
      
      sisum = (((((bend_bnd_sbend_6[ktemp]
                  *sint + bend_bnd_sbend_5[ktemp])
                 *sint + bend_bnd_sbend_4[ktemp])
                *sint + bend_bnd_sbend_3[ktemp])
               *sint + bend_bnd_sbend_2[ktemp])
              *sint);
      
      vbends = sisum*cost + bend_bnd_sbend_1[ktemp]*sint;
      
      intra_scr_vpot[ibend_bnd] = vbendc + vbends;

      /*-------------------------------------------------------------------*/
      /*  C) Get the force on the atoms using the chain rule and Horner's  */
      
      dvbendc =  (((((bend_bnd_dcbend_6[ktemp]
                      *cost + bend_bnd_dcbend_5[ktemp])
                     *cost + bend_bnd_dcbend_4[ktemp])
                    *cost + bend_bnd_dcbend_3[ktemp])
                   *cost + bend_bnd_dcbend_2[ktemp])
                  *cost + bend_bnd_dcbend_1[ktemp]);
      
      dvbends =  ((((bend_bnd_dsbend_6[ktemp]
                     *sint + bend_bnd_dsbend_5[ktemp])
                    *sint + bend_bnd_dsbend_4[ktemp])
                   *sint + bend_bnd_dsbend_3[ktemp])
                  *sint + bend_bnd_dsbend_2[ktemp]);
      
      dvbends = dvbends*cost + bend_bnd_dsbend_1[ktemp];
      
      pre     = -(dvbendc + sisum - dvbends*cost/sint);

      intra_scr_fx1[ibend_bnd] = (intra_scr_dx23[ibend_bnd]*rpmag
                                 -  intra_scr_dx12[ibend_bnd]*cos122)*pre;
      intra_scr_fy1[ibend_bnd] = (intra_scr_dy23[ibend_bnd]*rpmag
                                 -  intra_scr_dy12[ibend_bnd]*cos122)*pre;
      intra_scr_fz1[ibend_bnd] = (intra_scr_dz23[ibend_bnd]*rpmag
                                 -  intra_scr_dz12[ibend_bnd]*cos122)*pre;
      intra_scr_fx3[ibend_bnd] = (intra_scr_dx12[ibend_bnd]*rpmag
                                 -  intra_scr_dx23[ibend_bnd]*cos322)*pre;
      intra_scr_fy3[ibend_bnd] = (intra_scr_dy12[ibend_bnd]*rpmag
                                 -  intra_scr_dy23[ibend_bnd]*cos322)*pre;
      intra_scr_fz3[ibend_bnd] = (intra_scr_dz12[ibend_bnd]*rpmag
                                 -  intra_scr_dz23[ibend_bnd]*cos322)*pre;
      intra_scr_fx2[ibend_bnd] =-(intra_scr_fx1[ibend_bnd]
                                 +  intra_scr_fx3[ibend_bnd]); 
      intra_scr_fy2[ibend_bnd] =-(intra_scr_fy1[ibend_bnd]
                                 +  intra_scr_fy3[ibend_bnd]);
      intra_scr_fz2[ibend_bnd] =-(intra_scr_fz1[ibend_bnd]        
                                 +  intra_scr_fz3[ibend_bnd]);
      intra_scr_p11[ibend_bnd]=intra_scr_dx12[ibend_bnd]
                               *intra_scr_fx1[ibend_bnd]
                               +intra_scr_dx23[ibend_bnd]
                               *intra_scr_fx3[ibend_bnd];
      intra_scr_p22[ibend_bnd]=intra_scr_dy12[ibend_bnd]
                               *intra_scr_fy1[ibend_bnd]
                               +intra_scr_dy23[ibend_bnd]
                               *intra_scr_fy3[ibend_bnd];
      intra_scr_p33[ibend_bnd]=intra_scr_dz12[ibend_bnd]
                               *intra_scr_fz1[ibend_bnd]
                               +intra_scr_dz23[ibend_bnd]
                               *intra_scr_fz3[ibend_bnd];
      intra_scr_p12[ibend_bnd]=intra_scr_dx12[ibend_bnd]
                               *intra_scr_fy1[ibend_bnd]
                               +intra_scr_dx23[ibend_bnd]
                               *intra_scr_fy3[ibend_bnd];
      intra_scr_p13[ibend_bnd]=intra_scr_dx12[ibend_bnd]
                               *intra_scr_fz1[ibend_bnd]
                               +intra_scr_dx23[ibend_bnd]
                               *intra_scr_fz3[ibend_bnd];
      intra_scr_p23[ibend_bnd]=intra_scr_dy12[ibend_bnd]
                               *intra_scr_fz1[ibend_bnd]
                               +intra_scr_dy23[ibend_bnd]
                               *intra_scr_fz3[ibend_bnd];
    }/*endfor ibend_bnd */

    for(ibend_bnd=1;ibend_bnd <= nnow; ++ibend_bnd){
      *vbend_bnd_bend += intra_scr_vpot[ibend_bnd];
    }/*endfor*/
    
/*==========================================================================*/
/* VIII) Get the bond energy                                                */

    for(ibend_bnd=1,iboff=ibig;ibend_bnd <= nnow; ++ibend_bnd,++iboff) {
      ktemp = bend_bnd_jtyp[iboff];

      r122 = intra_scr_dx43[ibend_bnd]*intra_scr_dx43[ibend_bnd]
           + intra_scr_dy43[ibend_bnd]*intra_scr_dy43[ibend_bnd]
           + intra_scr_dz43[ibend_bnd]*intra_scr_dz43[ibend_bnd];
      r12  = sqrt(r122);
      rr0 = r12 - bend_bnd_eq_bond[ktemp];

      /*--------------------------------------------------------------------*/
      /*  A) Get the bend_bnding potential energy using Horner's method     */
      
      vbond =  ((((( bend_bnd_cbond_6[ktemp]
                    *rr0 + bend_bnd_cbond_5[ktemp])
                   *rr0 + bend_bnd_cbond_4[ktemp])
                  *rr0 + bend_bnd_cbond_3[ktemp])
                 *rr0 + bend_bnd_cbond_2[ktemp])
                *rr0 + bend_bnd_cbond_1[ktemp])
        *rr0 + bend_bnd_cbond_0[ktemp];
      
      intra_scr_vpot[ibend_bnd] += vbond;
      
      /*--------------------------------------------------------------------*/
      /*  B) Get the force on the atoms using the chain rule */
      
      dvbond =   (((( bend_bnd_dcbond_6[ktemp]
                     *rr0 + bend_bnd_dcbond_5[ktemp])
                    *rr0 + bend_bnd_dcbond_4[ktemp])
                   *rr0 + bend_bnd_dcbond_3[ktemp])
                  *rr0 + bend_bnd_dcbond_2[ktemp])
        *rr0 + bend_bnd_dcbond_1[ktemp];
      
      pre    = -dvbond/r12;
      fxtmp = intra_scr_dx43[ibend_bnd]*pre;
      fytmp = intra_scr_dy43[ibend_bnd]*pre;
      fztmp = intra_scr_dz43[ibend_bnd]*pre;
      intra_scr_fx1[ibend_bnd] += fxtmp;
      intra_scr_fy1[ibend_bnd] += fytmp;
      intra_scr_fz1[ibend_bnd] += fztmp;
      intra_scr_fx3[ibend_bnd] -= fxtmp;
      intra_scr_fy3[ibend_bnd] -= fytmp;
      intra_scr_fz3[ibend_bnd] -= fztmp;
      intra_scr_p11[ibend_bnd] += intra_scr_dx43[ibend_bnd]*fxtmp;
      intra_scr_p22[ibend_bnd] += intra_scr_dy43[ibend_bnd]*fytmp;
      intra_scr_p33[ibend_bnd] += intra_scr_dz43[ibend_bnd]*fztmp;
      intra_scr_p12[ibend_bnd] += intra_scr_dx43[ibend_bnd]*fytmp;
      intra_scr_p13[ibend_bnd] += intra_scr_dx43[ibend_bnd]*fztmp;
      intra_scr_p23[ibend_bnd] += intra_scr_dy43[ibend_bnd]*fztmp;
    /*endfor ibend_bnd*/}
    
/*==========================================================================*/
/*  IX) Sum the potential energy                                            */
    
    for(ibend_bnd=1;ibend_bnd <= nnow; ++ibend_bnd){
      *vbend_bndt += intra_scr_vpot[ibend_bnd];
    }/*endfor*/

/*==========================================================================*/
/* X) Sum the pressure tensor                                               */

    if(cell->iperd == 2 || cell->iperd == 3) {
      for(ibend_bnd=1;ibend_bnd <= nnow; ++ibend_bnd){
        ptens_pvten_tmp[1] += intra_scr_p11[ibend_bnd];
        ptens_pvten_tmp[5] += intra_scr_p22[ibend_bnd];
        ptens_pvten_tmp[9] += intra_scr_p33[ibend_bnd];
        ptens_pvten_tmp[2] += intra_scr_p12[ibend_bnd];
      }/*endfor*/
    }/*endif*/
    if(cell->iperd == 3) {
      for(ibend_bnd=1;ibend_bnd <= nnow; ++ibend_bnd){
        ptens_pvten_tmp[3] += intra_scr_p13[ibend_bnd];
        ptens_pvten_tmp[6] += intra_scr_p23[ibend_bnd];
      }/*endfor*/
    }/*endif*/

/*==========================================================================*/
/* XI) Scatter the forces                                                   */

if(iver_get==1){

    for(iboff=ibig+ioff,ibend_bnd=ibig;ibend_bnd <= iend; ++ibend_bnd,++iboff){
      ktemp = bend_bnd_j1[ibend_bnd];
      clatoms_fxt[ktemp] += intra_scr_fx1[iboff];
      clatoms_fyt[ktemp] += intra_scr_fy1[iboff];
      clatoms_fzt[ktemp] += intra_scr_fz1[iboff];
    }/*endfor*/

    for(iboff=ibig+ioff,ibend_bnd=ibig;ibend_bnd <= iend; ++ibend_bnd,++iboff){
      ktemp = bend_bnd_j2[ibend_bnd];
      clatoms_fxt[ktemp] += intra_scr_fx2[iboff];
      clatoms_fyt[ktemp] += intra_scr_fy2[iboff];
      clatoms_fzt[ktemp] += intra_scr_fz2[iboff];
    }/*endfor*/

    for(iboff=ibig+ioff,ibend_bnd=ibig;ibend_bnd <= iend; ++ibend_bnd,++iboff){
      ktemp = bend_bnd_j3[ibend_bnd];
      clatoms_fxt[ktemp] += intra_scr_fx3[iboff];
      clatoms_fyt[ktemp] += intra_scr_fy3[iboff];
      clatoms_fzt[ktemp] += intra_scr_fz3[iboff];
    }/*endfor*/

  }/*endif*/

    for(iboff=ibig+ioff;iboff <= iend+ioff; ++iboff) {
      intra_scr_fx1[iboff]  *= wfor;
      intra_scr_fy1[iboff]  *= wfor;
      intra_scr_fz1[iboff]  *= wfor;
      intra_scr_fx2[iboff]  *= wfor;
      intra_scr_fy2[iboff]  *= wfor;
      intra_scr_fz2[iboff]  *= wfor;
      intra_scr_fx3[iboff]  *= wfor;
      intra_scr_fy3[iboff]  *= wfor;
      intra_scr_fz3[iboff]  *= wfor;
    }/*endfor*/

    for(iboff=ibig+ioff,ibend_bnd=ibig;ibend_bnd <= iend; ++ibend_bnd,++iboff){
      ktemp = bend_bnd_j1[ibend_bnd];
      clatoms_fx[ktemp] += intra_scr_fx1[iboff];
      clatoms_fy[ktemp] += intra_scr_fy1[iboff];
      clatoms_fz[ktemp] += intra_scr_fz1[iboff];
    }/*endfor*/

    for(iboff=ibig+ioff,ibend_bnd=ibig;ibend_bnd <= iend; ++ibend_bnd,++iboff){
      ktemp = bend_bnd_j2[ibend_bnd];
      clatoms_fx[ktemp] += intra_scr_fx2[iboff];
      clatoms_fy[ktemp] += intra_scr_fy2[iboff];
      clatoms_fz[ktemp] += intra_scr_fz2[iboff];
    }/*endfor*/

    for(iboff=ibig+ioff,ibend_bnd=ibig;ibend_bnd <= iend; ++ibend_bnd,++iboff){
      ktemp = bend_bnd_j3[ibend_bnd];
      clatoms_fx[ktemp] += intra_scr_fx3[iboff];
      clatoms_fy[ktemp] += intra_scr_fy3[iboff];
      clatoms_fz[ktemp] += intra_scr_fz3[iboff];
    }/*endfor*/

  }/* endfor ibig */

/* END BIG LOOP */
/*======================================================================*/
/*======================================================================*/

/*======================================================================*/
/* XII) Increment the Pressure tensor */

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

  *vbend_bnd_bond = (*vbend_bndt)-(*vbend_bnd_bend);

  /*----------------------------------------------------------------------*/
}/*end routine*/
/*==========================================================================*/


