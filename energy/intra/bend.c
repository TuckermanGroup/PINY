/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
/*                                                                          */
/*                         PI_MD:                                           */
/*             The future of simulation technology                          */
/*             ------------------------------------                         */
/*                     Module: bend.c                                       */
/*                                                                          */
/* This routine computes the energy and forces from                         */ 
/* the intramolecular bend potential.                                       */
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

void bend(CLATOMS_INFO *clatoms_info,CLATOMS_POS *clatoms_pos,
          BEND *bend,CELL *cell,
          INTRA_SCR *intra_scr,PTENS *ptens,double *vbendt, int iver_get,
          CLASS_COMM_FORC_PKG *class_comm_forc_pkg,int iget_pv_real_inter)

/*==========================================================================*/
{/*begin routine*/
  /*=======================================================================*/
  /*            Local variable declarations                                */
  int iii;
  int ibend,ibig,ibig1,ioff,iend,nnow;     /* Indices and counters   */
  int iloop;

  double r122i,r12i;                       /* (r-r_0)^2 and (r-r_0)  */
  double r322i,r32i;
  double cost,sint,seps;
  double vbendc,vbends;                    /* Bond potential         */
  double dvbendc,dvbends;                  /* Derivative of bend pot */
  double sisum;                            /* sine sum               */
  double pre;                              /* Force prefactor        */
  double rpmag;                            /* 1/(r12*r32)            */
  double cos122,cos322;                    /* cos/r122, cos/r322     */
  double wfor;
  double x2,y2,z2;
  int i,ktemp,iboff,ktemp1,ktemp2,ktemp3,ktemp5,ktemp6,ktemp7;
  int nlen_use,nlen_now;
  int ntot;
  int ngo,irem,istart_big;

  /* Define local pointers                                                 */
  double *intra_scr_fx1   = intra_scr->fx1;
  double *intra_scr_fy1   = intra_scr->fy1;
  double *intra_scr_fz1   = intra_scr->fz1;
  double *intra_scr_fx2   = intra_scr->fx2;
  double *intra_scr_fy2   = intra_scr->fy2;
  double *intra_scr_fz2   = intra_scr->fz2;
  double *intra_scr_fx3   = intra_scr->fx3;
  double *intra_scr_fy3   = intra_scr->fy3;
  double *intra_scr_fz3   = intra_scr->fz3;
  double *intra_scr_dx12  = intra_scr->dx12;
  double *intra_scr_dy12  = intra_scr->dy12;
  double *intra_scr_dz12  = intra_scr->dz12;
  double *intra_scr_dx23  = intra_scr->dx23;
  double *intra_scr_dy23  = intra_scr->dy23;
  double *intra_scr_dz23  = intra_scr->dz23;
  double *intra_scr_dx56  = intra_scr->dx56;
  double *intra_scr_dy56  = intra_scr->dy56;
  double *intra_scr_dz56  = intra_scr->dz56;
  double *intra_scr_dx67  = intra_scr->dx67;
  double *intra_scr_dy67  = intra_scr->dy67;
  double *intra_scr_dz67  = intra_scr->dz67;
  double *intra_scr_vpot  = intra_scr->vpot;
  double *intra_scr_p11   = intra_scr->p11;
  double *intra_scr_p22   = intra_scr->p22;
  double *intra_scr_p33   = intra_scr->p33;
  double *intra_scr_p12   = intra_scr->p12;
  double *intra_scr_p13   = intra_scr->p13;
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
  double *bend_c_0        = bend->c_0;
  double *bend_c_1        = bend->c_1;
  double *bend_c_2        = bend->c_2;
  double *bend_c_3        = bend->c_3;
  double *bend_c_4        = bend->c_4;
  double *bend_c_5        = bend->c_5;
  double *bend_c_6        = bend->c_6;
  double *bend_s_0        = bend->s_0;
  double *bend_s_1        = bend->s_1;
  double *bend_s_2        = bend->s_2;
  double *bend_s_3        = bend->s_3;
  double *bend_s_4        = bend->s_4;
  double *bend_s_5        = bend->s_5;
  double *bend_s_6        = bend->s_6;
  double *bend_dc_0       = bend->dc_0;
  double *bend_dc_1       = bend->dc_1;
  double *bend_dc_2       = bend->dc_2;
  double *bend_dc_3       = bend->dc_3;
  double *bend_dc_4       = bend->dc_4;
  double *bend_dc_5       = bend->dc_5;
  double *bend_dc_6       = bend->dc_6;
  double *bend_ds_0       = bend->ds_0;
  double *bend_ds_1       = bend->ds_1;
  double *bend_ds_2       = bend->ds_2;
  double *bend_ds_3       = bend->ds_3;
  double *bend_ds_4       = bend->ds_4;
  double *bend_ds_5       = bend->ds_5;
  double *bend_ds_6       = bend->ds_6;
  double *ptens_pvten_tmp = ptens->pvten_tmp;
  double *xmod            = clatoms_info->xmod;
  double *ymod            = clatoms_info->ymod;
  double *zmod            = clatoms_info->zmod;
  int *bend_j1_pow        = bend->j1_pow;
  int *bend_j2_pow        = bend->j2_pow;
  int *bend_j3_pow        = bend->j3_pow;
  int *bend_jtyp_pow      = bend->jtyp_pow;
  int *iblock_size        = bend->iblock_pow_size;
  int *iblock_conflict_1  = bend->iblock_pow_conflict_1;
  int *iblock_conflict_2  = bend->iblock_pow_conflict_2;
  int *iblock_conflict_3  = bend->iblock_pow_conflict_3;
  int np_forc             = class_comm_forc_pkg->num_proc;
  int myid_forc           = class_comm_forc_pkg->myid;

  /*=======================================================================*/
  /* 0) Lower cutoff for sin(theta)                                        */

  seps = 1.0e-8;

  /* I) loop over all the bends in steps of nlen to save memory            */
  wfor    = (intra_scr->wght_tra_res);
  for(i=1;i<=9;i++){(ptens_pvten_tmp)[i]=0;}

/*=======================================================================*/
/*=======================================================================*/
/* LOOP OVER ALL THE BENDS IN STEPS OF NLEN TO SAVE MEMORY               */

  nlen_use = (intra_scr->nlen);
  ntot     = bend->npow;

  for(ibig=1;ibig <= ntot;ibig += nlen_use) {
    
    /*----------------------------------------------------------------------*/
    /*  A) Offsets to save some memory by doing only nlen bends at a time   */

    ibig1 = ibig-1;
    ioff = -ibig1;
    nlen_now = nlen_use;
    iend = MIN(ntot,ibig1+nlen_now);
    nnow = iend-ibig1;

    /*---------------------------------------------------------------------*/
    /*  B) Gather positions                                                */
    
    for(iboff=1,ibend=ibig;ibend <= iend; ++ibend,++iboff) {

      ktemp1 = bend_j1_pow[ibend];
      ktemp2 = bend_j2_pow[ibend];
      ktemp3 = bend_j3_pow[ibend];
      x2 = clatoms_x[ktemp2];
      y2 = clatoms_y[ktemp2];
      z2 = clatoms_z[ktemp2];
      intra_scr_dx12[iboff] = clatoms_x[ktemp1]-x2;
      intra_scr_dy12[iboff] = clatoms_y[ktemp1]-y2;
      intra_scr_dz12[iboff] = clatoms_z[ktemp1]-z2;
      intra_scr_dx23[iboff] = clatoms_x[ktemp3]-x2;
      intra_scr_dy23[iboff] = clatoms_y[ktemp3]-y2;
      intra_scr_dz23[iboff] = clatoms_z[ktemp3]-z2;

    }/*endfor*/

  if(cell->intra_perds==1&&clatoms_info->pi_beads>1){     

    for(iboff=1,ibend=ibig;ibend <= iend; ++ibend,++iboff) {
      ktemp5 = bend_j1_pow[ibend];
      ktemp6 = bend_j2_pow[ibend];
      ktemp7 = bend_j3_pow[ibend];
      intra_scr_dx56[iboff] = xmod[ktemp5]-xmod[ktemp6];
      intra_scr_dy56[iboff] = ymod[ktemp5]-ymod[ktemp6];
      intra_scr_dz56[iboff] = zmod[ktemp5]-zmod[ktemp6];
      intra_scr_dx67[iboff] = xmod[ktemp7]-xmod[ktemp6];
      intra_scr_dy67[iboff] = ymod[ktemp7]-ymod[ktemp6];
      intra_scr_dz67[iboff] = zmod[ktemp7]-zmod[ktemp6];
    }/*endfor*/

  }/*endif*/

    /*---------------------------------------------------------------------*/
    /*  E) Periodic boundary conditions                                    */
    
    if(cell->intra_perds == 1) {

      if(clatoms_info->pi_beads==1){     
        period(nnow,intra_scr_dx12,intra_scr_dy12,intra_scr_dz12,cell);
        period(nnow,intra_scr_dx23,intra_scr_dy23,intra_scr_dz23,cell);
      }else{
        period_pimd(nnow,intra_scr_dx12,intra_scr_dy12,intra_scr_dz12,
                         intra_scr_dx56,intra_scr_dy56,intra_scr_dz56,cell);
        period_pimd(nnow,intra_scr_dx23,intra_scr_dy23,intra_scr_dz23,
                         intra_scr_dx67,intra_scr_dy67,intra_scr_dz67,cell);
      }/*endif*/

    }/*endif*/
    
    for(ibend=1,iboff=ibig;ibend <= nnow; ++ibend,++iboff) {
      r122i  = 1.0/(intra_scr_dx12[ibend]*intra_scr_dx12[ibend]
                  + intra_scr_dy12[ibend]*intra_scr_dy12[ibend]
                  + intra_scr_dz12[ibend]*intra_scr_dz12[ibend]);
      r12i    = sqrt(r122i);
      
      r322i   = 1.0/(intra_scr_dx23[ibend]*intra_scr_dx23[ibend]
                   + intra_scr_dy23[ibend]*intra_scr_dy23[ibend]
                   + intra_scr_dz23[ibend]*intra_scr_dz23[ibend]);
      r32i    = sqrt(r322i);
      rpmag  = r12i*r32i;
      
      /*-------------------------------------------------------------------*/
      /*  F) Calculate the cosine and sine of the angle between them       */
      /*     (cos(theta_{123}), sine(theta_{123}))                         */
      
      cost = (intra_scr_dx12[ibend]*intra_scr_dx23[ibend]
              +  intra_scr_dy12[ibend]*intra_scr_dy23[ibend]
              +  intra_scr_dz12[ibend]*intra_scr_dz23[ibend])*rpmag;
      
      cost   = (cost < 1.0 ? cost:1.0);
      cost   = (cost > -1.0 ? cost:-1.0);
      
      sint   = sqrt(1.0 - cost*cost);
      sint   = (sint > seps ? sint:seps);
      
      cos122 = cost*r122i;
      cos322  = cost*r322i;
      
      ktemp1 = bend_jtyp_pow[iboff];

      /*------------------------------------------------------------------*/
      /*  H) Get the bending potential energy from the cosine and sine power*/
      /*     series using Horner's method                                   */
      
      vbendc = ((((((bend_c_6[ktemp1]
                     *cost + bend_c_5[ktemp1])
                    *cost + bend_c_4[ktemp1])
                   *cost + bend_c_3[ktemp1])
                  *cost + bend_c_2[ktemp1])
                 *cost + bend_c_1[ktemp1])
                *cost + bend_c_0[ktemp1]);
      
      sisum = (((((bend_s_6[ktemp1]
                  *sint + bend_s_5[ktemp1])
                 *sint + bend_s_4[ktemp1])
                *sint + bend_s_3[ktemp1])
               *sint + bend_s_2[ktemp1])
              *sint);
      
      vbends = sisum*cost + bend_s_1[ktemp1]*sint;
      
      intra_scr_vpot[ibend] = vbendc + vbends;

      /*-------------------------------------------------------------------*/
      /*  I) Get the force on the atoms using the chain rule and Horner's  */
      
      dvbendc =  (((((bend_dc_6[ktemp1]
                      *cost + bend_dc_5[ktemp1])
                     *cost + bend_dc_4[ktemp1])
                    *cost + bend_dc_3[ktemp1])
                   *cost + bend_dc_2[ktemp1])
                  *cost + bend_dc_1[ktemp1]);
      
      dvbends =  ((((bend_ds_6[ktemp1]
                     *sint + bend_ds_5[ktemp1])
                    *sint + bend_ds_4[ktemp1])
                   *sint + bend_ds_3[ktemp1])
                  *sint + bend_ds_2[ktemp1]);
      
      dvbends = dvbends*cost + bend_ds_1[ktemp1];
      
      pre     = -(dvbendc + sisum - dvbends*cost/sint);

      intra_scr_fx1[ibend] = (intra_scr_dx23[ibend]*rpmag
                           -  intra_scr_dx12[ibend]*cos122)*pre;
      intra_scr_fy1[ibend] = (intra_scr_dy23[ibend]*rpmag
                           -  intra_scr_dy12[ibend]*cos122)*pre;
      intra_scr_fz1[ibend] = (intra_scr_dz23[ibend]*rpmag
                           -  intra_scr_dz12[ibend]*cos122)*pre;
      intra_scr_fx3[ibend] = (intra_scr_dx12[ibend]*rpmag
                           -  intra_scr_dx23[ibend]*cos322)*pre;
      intra_scr_fy3[ibend] = (intra_scr_dy12[ibend]*rpmag
                           -  intra_scr_dy23[ibend]*cos322)*pre;
      intra_scr_fz3[ibend] = (intra_scr_dz12[ibend]*rpmag
                           -  intra_scr_dz23[ibend]*cos322)*pre;
      intra_scr_fx2[ibend] =-(intra_scr_fx1[ibend]
                           +  intra_scr_fx3[ibend]); 
      intra_scr_fy2[ibend] =-(intra_scr_fy1[ibend]
                           +  intra_scr_fy3[ibend]);
      intra_scr_fz2[ibend] =-(intra_scr_fz1[ibend]        
                           +  intra_scr_fz3[ibend]);
    }/*endfor ibend */
    
    /*----------------------------------------------------------------------*/
    /*  J) Get the pressure tensor                                          */
    
    if(cell->iperd == 2 || cell->iperd == 3) {
      for(ibend=1;ibend <= nnow; ++ibend) {
        
        intra_scr_p11[ibend]  = intra_scr_dx12[ibend]*intra_scr_fx1[ibend]
                              + intra_scr_dx23[ibend]*intra_scr_fx3[ibend];
        intra_scr_p22[ibend]  = intra_scr_dy12[ibend]*intra_scr_fy1[ibend]
                              + intra_scr_dy23[ibend]*intra_scr_fy3[ibend];
        intra_scr_p33[ibend]  = intra_scr_dz12[ibend]*intra_scr_fz1[ibend]
                              + intra_scr_dz23[ibend]*intra_scr_fz3[ibend];
        intra_scr_p12[ibend]  = intra_scr_dx12[ibend]*intra_scr_fy1[ibend]
                              + intra_scr_dx23[ibend]*intra_scr_fy3[ibend];
      }/*endfor*/
    }/*endif*/
    if(cell->iperd == 3) {
      for(ibend=1;ibend <= nnow; ++ibend) {
        intra_scr_p13[ibend] = intra_scr_dx12[ibend]*intra_scr_fz1[ibend]
                             + intra_scr_dx23[ibend]*intra_scr_fz3[ibend];
        intra_scr_p23[ibend] = intra_scr_dy12[ibend]*intra_scr_fz1[ibend]
                             + intra_scr_dy23[ibend]*intra_scr_fz3[ibend];
      }/*endfor*/
    }/*endif*/
    /*---------------------------------------------------------------------*/
    /*  K) Sum the potential energy                                        */
    
    for(ibend=1;ibend <= nnow; ++ibend){
      *vbendt += intra_scr_vpot[ibend];
    }/*endfor*/
    /*---------------------------------------------------------------------*/
    /*  L) Sum the pressure tensor                                         */

    if(cell->iperd == 2 || cell->iperd == 3) {
      for(ibend=1;ibend <= nnow; ++ibend){
        ptens_pvten_tmp[1] += intra_scr_p11[ibend];
        ptens_pvten_tmp[5] += intra_scr_p22[ibend];
        ptens_pvten_tmp[9] += intra_scr_p33[ibend];
        ptens_pvten_tmp[2] += intra_scr_p12[ibend];
      }/*endfor*/
    }/*endif*/
    if(cell->iperd == 3) {
      for(ibend=1;ibend <= nnow; ++ibend){
        ptens_pvten_tmp[3] += intra_scr_p13[ibend];
        ptens_pvten_tmp[6] += intra_scr_p23[ibend];
      }/*endfor*/
    }/*endif*/
    /*---------------------------------------------------------------------*/
    /*  M) Scatter the forces                                              */
    if(iver_get==1){

      for(iboff=1,ibend=ibig;ibend <= iend; ++ibend,++iboff) {
       ktemp = bend_j1_pow[ibend];
       clatoms_fxt[ktemp] += intra_scr_fx1[iboff];
       clatoms_fyt[ktemp] += intra_scr_fy1[iboff];
       clatoms_fzt[ktemp] += intra_scr_fz1[iboff];
      }

      for(iboff=1,ibend=ibig;ibend <= iend; ++ibend,++iboff) {
       ktemp = bend_j2_pow[ibend];
       clatoms_fxt[ktemp] += intra_scr_fx2[iboff];
       clatoms_fyt[ktemp] += intra_scr_fy2[iboff];
       clatoms_fzt[ktemp] += intra_scr_fz2[iboff];
      }/*endfor*/

      for(iboff=1,ibend=ibig;ibend <= iend; ++ibend,++iboff) {
        ktemp = bend_j3_pow[ibend];
        clatoms_fxt[ktemp] += intra_scr_fx3[iboff];
        clatoms_fyt[ktemp] += intra_scr_fy3[iboff];
        clatoms_fzt[ktemp] += intra_scr_fz3[iboff];
      }/*endfor*/
    }/*endif*/

    for(iboff=1;iboff <= iend+ioff; ++iboff) {
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

     for(iboff=1,ibend=ibig;ibend <= iend; ++ibend,++iboff) {
       ktemp = bend_j1_pow[ibend];
       clatoms_fx[ktemp] += intra_scr_fx1[iboff];
       clatoms_fy[ktemp] += intra_scr_fy1[iboff];
       clatoms_fz[ktemp] += intra_scr_fz1[iboff];
     }/*endfor*/

     for(iboff=1,ibend=ibig;ibend <= iend; ++ibend,++iboff) {
      ktemp = bend_j2_pow[ibend];
      clatoms_fx[ktemp] += intra_scr_fx2[iboff];
      clatoms_fy[ktemp] += intra_scr_fy2[iboff];
      clatoms_fz[ktemp] += intra_scr_fz2[iboff];
     }/*endfor*/

     for(iboff=1,ibend=ibig;ibend <= iend; ++ibend,++iboff) {
      ktemp = bend_j3_pow[ibend];
      clatoms_fx[ktemp] += intra_scr_fx3[iboff];
      clatoms_fy[ktemp] += intra_scr_fy3[iboff];
      clatoms_fz[ktemp] += intra_scr_fz3[iboff];
     }/*endfor*/
  }/* endfor ibig */
  /*======================================================================*/
  /* Increment the Pressure tensor */
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

  /*----------------------------------------------------------------------*/
}/*end routine*/
/*==========================================================================*/


/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void bend_free(CLATOMS_INFO *clatoms_info,CLATOMS_POS *clatoms_pos,
               BEND_FREE *bend_free,CELL *cell,
               INTRA_SCR *intra_scr,PTENS *ptens,double *vbendt,
               ENERGY_CTRL *energy_ctrl,int np_forc)

/*==========================================================================*/
{/*begin routine*/
/*==========================================================================*/
  /*         Local variable declarations                                */
  int ibin,nnow;                           /* Indices and counters   */
  int pi_beads = clatoms_info->pi_beads; 
  double r122,r12;                         /* (r-r_0)^2 and (r-r_0)  */
  double r322,r32;
  double cost,sint,theta,seps;
  double dvbend;                          /* Derivative of bend pot */
  double pre;                              /* Force prefactor        */
  double rpmag;                            /* 1/(r12*r32)            */
  double cos122,cos322;                    /* cos/r122, cos/r322     */
  double apow,apow1;
  double theta0;
  double wfor;
  double *xmod            = clatoms_info->xmod;
  double *ymod            = clatoms_info->ymod;
  double *zmod            = clatoms_info->zmod;
  int iget_pv_real_inter = energy_ctrl->iget_pv_real_inter;

  /*======================================================================*/
  /* 0) Lower cutoff for sin(theta)                                       */
  
  seps = 1.0e-8;
  /*----------------------------------------------------------------------*/
  /*  I) Gather positions                                                 */
  
  wfor    = intra_scr->wght_tra_res;
  intra_scr->x1[1] = clatoms_pos->x[bend_free->j1];
  intra_scr->y1[1] = clatoms_pos->y[bend_free->j1];
  intra_scr->z1[1] = clatoms_pos->z[bend_free->j1];
  intra_scr->x2[1] = clatoms_pos->x[bend_free->j2];
  intra_scr->y2[1] = clatoms_pos->y[bend_free->j2];
  intra_scr->z2[1] = clatoms_pos->z[bend_free->j2];
  intra_scr->x3[1] = clatoms_pos->x[bend_free->j3];
  intra_scr->y3[1] = clatoms_pos->y[bend_free->j3];
  intra_scr->z3[1] = clatoms_pos->z[bend_free->j3];
  
  if(cell->intra_perds==1&&clatoms_info->pi_beads>1){     
    intra_scr->x5[1] = xmod[bend_free->j1];
    intra_scr->y5[1] = ymod[bend_free->j1];
    intra_scr->z5[1] = zmod[bend_free->j1];
    intra_scr->x6[1] = xmod[bend_free->j2];
    intra_scr->y6[1] = ymod[bend_free->j2];
    intra_scr->z6[1] = zmod[bend_free->j2];
    intra_scr->x7[1] = xmod[bend_free->j3];
    intra_scr->y7[1] = ymod[bend_free->j3];
    intra_scr->z7[1] = zmod[bend_free->j3];
  }/*endif*/

  /*----------------------------------------------------------------------*/
  /* II) Calculate the basis vectors r12,r32                              */

  intra_scr->dx12[1] = intra_scr->x1[1] - intra_scr->x2[1];
  intra_scr->dy12[1] = intra_scr->y1[1] - intra_scr->y2[1];
  intra_scr->dz12[1] = intra_scr->z1[1] - intra_scr->z2[1];
  intra_scr->dx23[1] = intra_scr->x3[1] - intra_scr->x2[1];
  intra_scr->dy23[1] = intra_scr->y3[1] - intra_scr->y2[1];
  intra_scr->dz23[1] = intra_scr->z3[1] - intra_scr->z2[1];

  if(cell->intra_perds==1&&clatoms_info->pi_beads>1){     
    intra_scr->dx56[1] = intra_scr->x5[1] - intra_scr->x6[1];
    intra_scr->dy56[1] = intra_scr->y5[1] - intra_scr->y6[1];
    intra_scr->dz56[1] = intra_scr->z5[1] - intra_scr->z6[1];
    intra_scr->dx67[1] = intra_scr->x7[1] - intra_scr->x6[1];
    intra_scr->dy67[1] = intra_scr->y7[1] - intra_scr->y6[1];
    intra_scr->dz67[1] = intra_scr->z7[1] - intra_scr->z6[1];
  }/*endif*/

  /*---------------------------------------------------------------------*/
  /*  III) Periodic boundary conditions                                  */

  nnow = 1;
    if(cell->intra_perds == 1) {
      if(clatoms_info->pi_beads==1){     
        period(nnow,intra_scr->dx12,intra_scr->dy12,intra_scr->dz12,cell);
        period(nnow,intra_scr->dx23,intra_scr->dy23,intra_scr->dz23,cell);
      }else{
        period_pimd(nnow,intra_scr->dx12,intra_scr->dy12,intra_scr->dz12,
                      intra_scr->dx56,intra_scr->dy56,intra_scr->dz56,cell);
        period_pimd(nnow,intra_scr->dx23,intra_scr->dy23,intra_scr->dz23,
                   intra_scr->dx67,intra_scr->dy67,intra_scr->dz67,cell);
      }
   }
  r122   = intra_scr->dx12[1]*intra_scr->dx12[1]
    + intra_scr->dy12[1]*intra_scr->dy12[1]
      + intra_scr->dz12[1]*intra_scr->dz12[1];
  r12    = sqrt(r122);
  
  r322   = intra_scr->dx23[1]*intra_scr->dx23[1]
    + intra_scr->dy23[1]*intra_scr->dy23[1]
      + intra_scr->dz23[1]*intra_scr->dz23[1];
  r32    = sqrt(r322);
  rpmag  = 1.0/(r12*r32);
  
  /*----------------------------------------------------------------------*/
  /*  IV) Calculate the cosine and sine of the angle between them         */
  /*     (cos(theta_{123}), sine(theta_{123}))                            */
  
  cost     = (intra_scr->dx12[1]*intra_scr->dx23[1]
              +  intra_scr->dy12[1]*intra_scr->dy23[1]
              +  intra_scr->dz12[1]*intra_scr->dz23[1])/(r12*r32);
  
  cost   = (cost < 1.0 ? cost:1.0);
  cost   = (cost > -1.0 ? cost:-1.0);
  
  theta  = acos(cost);
  sint   = sin(theta);
  theta0 = theta - bend_free->eq;
  
  cos122 = cost/r122;
  cos322  = cost/r322;
  rpmag  = 1.0/(r12*r32);

  /*-----------------------------------------------------------------------*/
  /*  IV) Calculate PE and fill bins                                       */
  
  apow      = (double)(bend_free->npow);
  apow1     = apow-1.0;
  
  *vbendt = bend_free->fk*pow(theta0,apow)/apow;
  
  if(energy_ctrl->iget_full_inter==1){
   ibin      = (int) (theta/(bend_free->del) + 1);
   ibin      = MIN((bend_free->nhist),ibin);
   ibin      = MAX(1,ibin);
   bend_free->hist[ibin] += 1.0;
  }/*endif*/  
  /*-----------------------------------------------------------------------*/
  /*  I) Get the force on the atoms using the chain rule                   */
  dvbend = bend_free->fk*pow(theta0,apow1);
  pre    = dvbend/sint;
  intra_scr->fx1[1] = (intra_scr->dx23[1]*rpmag-intra_scr->dx12[1]*cos122)*pre;
  intra_scr->fy1[1] = (intra_scr->dy23[1]*rpmag-intra_scr->dy12[1]*cos122)*pre;
  intra_scr->fz1[1] = (intra_scr->dz23[1]*rpmag-intra_scr->dz12[1]*cos122)*pre;
  intra_scr->fx3[1] = (intra_scr->dx12[1]*rpmag-intra_scr->dx23[1]*cos322)*pre;
  intra_scr->fy3[1] = (intra_scr->dy12[1]*rpmag-intra_scr->dy23[1]*cos322)*pre;
  intra_scr->fz3[1] = (intra_scr->dz12[1]*rpmag-intra_scr->dz23[1]*cos322)*pre;
  intra_scr->fx2[1] =-(intra_scr->fx1[1] + intra_scr->fx3[1]); 
  intra_scr->fy2[1] =-(intra_scr->fy1[1] + intra_scr->fy3[1]);
  intra_scr->fz2[1] =-(intra_scr->fz1[1] + intra_scr->fz3[1]);
  
  /*----------------------------------------------------------------------*/
  /*  II) Get the pressure tensor                                         */
  
  if(cell->iperd == 2 || cell->iperd == 3) {
    intra_scr->p11[1] = intra_scr->dx12[1]*intra_scr->fx1[1]
                      + intra_scr->dx23[1]*intra_scr->fx3[1];
    intra_scr->p22[1] = intra_scr->dy12[1]*intra_scr->fy1[1]
                      + intra_scr->dy23[1]*intra_scr->fy3[1];
    intra_scr->p33[1] = intra_scr->dz12[1]*intra_scr->fz1[1]
                      + intra_scr->dz23[1]*intra_scr->fz3[1];
    intra_scr->p12[1] = intra_scr->dx12[1]*intra_scr->fy1[1]
                      + intra_scr->dx23[1]*intra_scr->fy3[1];
    intra_scr->p21[1] =  intra_scr->dy12[1]*intra_scr->fx1[1]
                       + intra_scr->dy23[1]*intra_scr->fx3[1];
  }/*endif*/
  if(cell->iperd == 3) {
    intra_scr->p13[1] = intra_scr->dx12[1]*intra_scr->fz1[1]
                      + intra_scr->dx23[1]*intra_scr->fz3[1];
    intra_scr->p31[1] = intra_scr->dz12[1]*intra_scr->fx1[1]
                      + intra_scr->dz23[1]*intra_scr->fx3[1];
    intra_scr->p23[1] = intra_scr->dy12[1]*intra_scr->fz1[1]
                      + intra_scr->dy23[1]*intra_scr->fz3[1];
    intra_scr->p32[1] = intra_scr->dz12[1]*intra_scr->fy1[1]
                      + intra_scr->dz23[1]*intra_scr->fy3[1];
  }/*endif*/
  /*-----------------------------------------------------------------------*/
  /*  L) Sum the pressure tensor                                           */
  
  if(cell->iperd == 2 || cell->iperd == 3) {
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
  if(cell->iperd == 3) {
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
  /*-----------------------------------------------------------------------*/
  /*  M) Scatter the forces                                                */
  clatoms_pos->fx[bend_free->j1]  += intra_scr->fx1[1]*wfor*pi_beads;
  clatoms_pos->fy[bend_free->j1]  += intra_scr->fy1[1]*wfor*pi_beads;
  clatoms_pos->fz[bend_free->j1]  += intra_scr->fz1[1]*wfor*pi_beads;
  clatoms_pos->fx[bend_free->j2]  += intra_scr->fx2[1]*wfor*pi_beads;
  clatoms_pos->fy[bend_free->j2]  += intra_scr->fy2[1]*wfor*pi_beads;
  clatoms_pos->fz[bend_free->j2]  += intra_scr->fz2[1]*wfor*pi_beads;
  clatoms_pos->fx[bend_free->j3]  += intra_scr->fx3[1]*wfor*pi_beads;
  clatoms_pos->fy[bend_free->j3]  += intra_scr->fy3[1]*wfor*pi_beads;
  clatoms_pos->fz[bend_free->j3]  += intra_scr->fz3[1]*wfor*pi_beads;
  /*-----------------------------------------------------------------------*/
}/*end routine*/
/*==========================================================================*/



/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void bend_free_mode(CLATOMS_INFO *clatoms_info,CLATOMS_POS *clatoms_pos,
               BEND_FREE *bend_free,CELL *cell,
               INTRA_SCR *intra_scr,PTENS *ptens,double *vbendt,
               ENERGY_CTRL *energy_ctrl,int np_forc)

/*==========================================================================*/
{/*begin routine*/
/*==========================================================================*/
  /*         Local variable declarations                                */
  int ibin,nnow;                           /* Indices and counters   */
  int pi_beads = clatoms_info->pi_beads; 
  double r122,r12;                         /* (r-r_0)^2 and (r-r_0)  */
  double r322,r32;
  double cost,sint,theta,seps;
  double dvbend;                          /* Derivative of bend pot */
  double pre;                              /* Force prefactor        */
  double rpmag;                            /* 1/(r12*r32)            */
  double cos122,cos322;                    /* cos/r122, cos/r322     */
  double apow,apow1;
  double theta0;
  double wfor;
  int iget_pv_real_inter = energy_ctrl->iget_pv_real_inter;

  /*======================================================================*/
  /* 0) Lower cutoff for sin(theta)                                       */
  
  seps = 1.0e-8;
  /*----------------------------------------------------------------------*/
  /*  I) Gather positions                                                 */
  
  wfor    = intra_scr->wght_tra_res;
  intra_scr->x1[1] = clatoms_info->xmod[bend_free->j1];
  intra_scr->y1[1] = clatoms_info->ymod[bend_free->j1];
  intra_scr->z1[1] = clatoms_info->zmod[bend_free->j1];
  intra_scr->x2[1] = clatoms_info->xmod[bend_free->j2];
  intra_scr->y2[1] = clatoms_info->ymod[bend_free->j2];
  intra_scr->z2[1] = clatoms_info->zmod[bend_free->j2];
  intra_scr->x3[1] = clatoms_info->xmod[bend_free->j3];
  intra_scr->y3[1] = clatoms_info->ymod[bend_free->j3];
  intra_scr->z3[1] = clatoms_info->zmod[bend_free->j3];
  
  /*----------------------------------------------------------------------*/
  /* II) Calculate the basis vectors r12,r32                              */

  intra_scr->dx12[1] = intra_scr->x1[1] - intra_scr->x2[1];
  intra_scr->dy12[1] = intra_scr->y1[1] - intra_scr->y2[1];
  intra_scr->dz12[1] = intra_scr->z1[1] - intra_scr->z2[1];
  intra_scr->dx23[1] = intra_scr->x3[1] - intra_scr->x2[1];
  intra_scr->dy23[1] = intra_scr->y3[1] - intra_scr->y2[1];
  intra_scr->dz23[1] = intra_scr->z3[1] - intra_scr->z2[1];

  /*---------------------------------------------------------------------*/
  /*  III) Periodic boundary conditions                                  */

  nnow = 1;
     period(nnow,intra_scr->dx12,intra_scr->dy12,intra_scr->dz12,cell);
     period(nnow,intra_scr->dx23,intra_scr->dy23,intra_scr->dz23,cell);

  r122   = intra_scr->dx12[1]*intra_scr->dx12[1]
    + intra_scr->dy12[1]*intra_scr->dy12[1]
      + intra_scr->dz12[1]*intra_scr->dz12[1];
  r12    = sqrt(r122);
  
  r322   = intra_scr->dx23[1]*intra_scr->dx23[1]
    + intra_scr->dy23[1]*intra_scr->dy23[1]
      + intra_scr->dz23[1]*intra_scr->dz23[1];
  r32    = sqrt(r322);
  rpmag  = 1.0/(r12*r32);
  
  /*----------------------------------------------------------------------*/
  /*  IV) Calculate the cosine and sine of the angle between them         */
  /*     (cos(theta_{123}), sine(theta_{123}))                            */
  
  cost     = (intra_scr->dx12[1]*intra_scr->dx23[1]
              +  intra_scr->dy12[1]*intra_scr->dy23[1]
              +  intra_scr->dz12[1]*intra_scr->dz23[1])/(r12*r32);
  
  cost   = (cost < 1.0 ? cost:1.0);
  cost   = (cost > -1.0 ? cost:-1.0);
  
  theta  = acos(cost);
  sint   = sin(theta);
  theta0 = theta - bend_free->eq;
  
  cos122 = cost/r122;
  cos322  = cost/r322;
  rpmag  = 1.0/(r12*r32);

  /*-----------------------------------------------------------------------*/
  /*  IV) Calculate PE and fill bins                                       */
  
  apow      = (double)(bend_free->npow);
  apow1     = apow-1.0;
  
  *vbendt = bend_free->fk*pow(theta0,apow)/apow;
  
  if(energy_ctrl->iget_full_inter==1){
   ibin      = (int) (theta/(bend_free->del) + 1);
   ibin      = MIN((bend_free->nhist),ibin);
   ibin      = MAX(1,ibin);
   bend_free->hist[ibin] += 1.0;
  }/*endif*/  
  /*-----------------------------------------------------------------------*/
  /*  I) Get the force on the atoms using the chain rule                   */
  dvbend = bend_free->fk*pow(theta0,apow1);
  pre    = dvbend/sint;
  intra_scr->fx1[1] = (intra_scr->dx23[1]*rpmag-intra_scr->dx12[1]*cos122)*pre;
  intra_scr->fy1[1] = (intra_scr->dy23[1]*rpmag-intra_scr->dy12[1]*cos122)*pre;
  intra_scr->fz1[1] = (intra_scr->dz23[1]*rpmag-intra_scr->dz12[1]*cos122)*pre;
  intra_scr->fx3[1] = (intra_scr->dx12[1]*rpmag-intra_scr->dx23[1]*cos322)*pre;
  intra_scr->fy3[1] = (intra_scr->dy12[1]*rpmag-intra_scr->dy23[1]*cos322)*pre;
  intra_scr->fz3[1] = (intra_scr->dz12[1]*rpmag-intra_scr->dz23[1]*cos322)*pre;
  intra_scr->fx2[1] =-(intra_scr->fx1[1] + intra_scr->fx3[1]); 
  intra_scr->fy2[1] =-(intra_scr->fy1[1] + intra_scr->fy3[1]);
  intra_scr->fz2[1] =-(intra_scr->fz1[1] + intra_scr->fz3[1]);
  
  /*----------------------------------------------------------------------*/
  /*  II) Get the pressure tensor                                         */
  
  if(cell->iperd == 2 || cell->iperd == 3) {
    intra_scr->p11[1] = intra_scr->dx12[1]*intra_scr->fx1[1]
                      + intra_scr->dx23[1]*intra_scr->fx3[1];
    intra_scr->p22[1] = intra_scr->dy12[1]*intra_scr->fy1[1]
                      + intra_scr->dy23[1]*intra_scr->fy3[1];
    intra_scr->p33[1] = intra_scr->dz12[1]*intra_scr->fz1[1]
                      + intra_scr->dz23[1]*intra_scr->fz3[1];
    intra_scr->p12[1] = intra_scr->dx12[1]*intra_scr->fy1[1]
                      + intra_scr->dx23[1]*intra_scr->fy3[1];
    intra_scr->p21[1] =  intra_scr->dy12[1]*intra_scr->fx1[1]
                       + intra_scr->dy23[1]*intra_scr->fx3[1];
  }/*endif*/
  if(cell->iperd == 3) {
    intra_scr->p13[1] = intra_scr->dx12[1]*intra_scr->fz1[1]
                      + intra_scr->dx23[1]*intra_scr->fz3[1];
    intra_scr->p31[1] = intra_scr->dz12[1]*intra_scr->fx1[1]
                      + intra_scr->dz23[1]*intra_scr->fx3[1];
    intra_scr->p23[1] = intra_scr->dy12[1]*intra_scr->fz1[1]
                      + intra_scr->dy23[1]*intra_scr->fz3[1];
    intra_scr->p32[1] = intra_scr->dz12[1]*intra_scr->fy1[1]
                      + intra_scr->dz23[1]*intra_scr->fy3[1];
  }/*endif*/
  /*-----------------------------------------------------------------------*/
  /*  L) Sum the pressure tensor                                           */
  
  if(cell->iperd == 2 || cell->iperd == 3) {
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
  if(cell->iperd == 3) {
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
  /*-----------------------------------------------------------------------*/
  /*  M) Scatter the forces                                                */
  clatoms_pos->fx[bend_free->j1]  += intra_scr->fx1[1]*wfor;
  clatoms_pos->fy[bend_free->j1]  += intra_scr->fy1[1]*wfor;
  clatoms_pos->fz[bend_free->j1]  += intra_scr->fz1[1]*wfor;
  clatoms_pos->fx[bend_free->j2]  += intra_scr->fx2[1]*wfor;
  clatoms_pos->fy[bend_free->j2]  += intra_scr->fy2[1]*wfor;
  clatoms_pos->fz[bend_free->j2]  += intra_scr->fz2[1]*wfor;
  clatoms_pos->fx[bend_free->j3]  += intra_scr->fx3[1]*wfor;
  clatoms_pos->fy[bend_free->j3]  += intra_scr->fy3[1]*wfor;
  clatoms_pos->fz[bend_free->j3]  += intra_scr->fz3[1]*wfor;
  /*-----------------------------------------------------------------------*/
}/*end routine*/
/*==========================================================================*/








