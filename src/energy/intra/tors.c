/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
/*                                                                          */
/*                         PI_MD:                                           */
/*             The future of simulation technology                          */
/*             ------------------------------------                         */
/*                     Module: tors.c                                       */
/*                                                                          */
/* This routine computes the energy and forces from                         */ 
/* the intramolecular tors potential.                                       */
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

void tors(CLATOMS_INFO *clatoms_info,CLATOMS_POS *clatoms_pos,
          TORS *tors,CELL *cell,
          INTRA_SCR *intra_scr,PTENS *ptens, double *vtorst, int iver_get,
          CLASS_COMM_FORC_PKG *class_comm_forc_pkg,
          int iget_pv_real_inter)

/*==========================================================================*/
{/*begin routine*/
  /*=======================================================================*/
  /*            Local variable declarations                                */
  int itors,ibig,ibig1,ioff,iend,nnow;       /* Counters,offsets, etc  */
  int nlen;
  int ktemp,itoff;
  int nlen_now;
  int ngo,irem,istart_big;

  double r122,r12;
  double r132,r13;
  double r232,r23;
  double r422,r42;
  double r432,r43;
  double temp1p,temp1q,temp2q,temp2p;
  double temp3p,temp3q,temp4q,temp4p;
                                               /* Distances              */
  double dxp,dyp,dzp,rp2i,rpi,rpqi;
  double dxq,dyq,dzq,rq2i,rqi;     
  double dp13,dp42,dq13,dq42;
  double cost,sint,seps,vtorsc,vtorss,dvtorsc,dvtorss,pre;          
  double sisum;
  double costrp2,costrq2,dp13crp2,dq42crq2;
  double dp13rpq,dq13rpq,dp42rpq,dq42rpq;
  double d1323i,palp,palp1;
  double d4223i,qalp,qalp1; 
  double fpx1,fpy1,fpz1;
  double fpx2,fpy2,fpz2;  
  double fqx4,fqy4,fqz4;
  double fqx2,fqy2,fqz2;  
  int i;
  int ktemp1,ktemp2,ktemp3,ktemp4,ktemp5,ktemp6,ktemp7,ktemp8;
  int iboff;
  double x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,x5,y5,z5,x6,y6,z6;
  double x7,y7,z7,x8,y8,z8;


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
  double *intra_scr_fx4   = intra_scr->fx4;
  double *intra_scr_fy4   = intra_scr->fy4;
  double *intra_scr_fz4   = intra_scr->fz4;
  double *intra_scr_dx12  = intra_scr->dx12;
  double *intra_scr_dy12  = intra_scr->dy12;
  double *intra_scr_dz12  = intra_scr->dz12;
  double *intra_scr_dx13  = intra_scr->dx13;
  double *intra_scr_dy13  = intra_scr->dy13;
  double *intra_scr_dz13  = intra_scr->dz13;
  double *intra_scr_dx42  = intra_scr->dx42;
  double *intra_scr_dy42  = intra_scr->dy42;
  double *intra_scr_dz42  = intra_scr->dz42;
  double *intra_scr_dx43  = intra_scr->dx43;
  double *intra_scr_dy43  = intra_scr->dy43;
  double *intra_scr_dz43  = intra_scr->dz43;
  double *intra_scr_dx23  = intra_scr->dx23;
  double *intra_scr_dy23  = intra_scr->dy23;
  double *intra_scr_dz23  = intra_scr->dz23;
  double *intra_scr_dx56  = intra_scr->dx56;
  double *intra_scr_dy56  = intra_scr->dy56;
  double *intra_scr_dz56  = intra_scr->dz56;
  double *intra_scr_dx57  = intra_scr->dx57;
  double *intra_scr_dy57  = intra_scr->dy57;
  double *intra_scr_dz57  = intra_scr->dz57;
  double *intra_scr_dx86  = intra_scr->dx86;
  double *intra_scr_dy86  = intra_scr->dy86;
  double *intra_scr_dz86  = intra_scr->dz86;
  double *intra_scr_dx87  = intra_scr->dx87;
  double *intra_scr_dy87  = intra_scr->dy87;
  double *intra_scr_dz87  = intra_scr->dz87;
  double *intra_scr_dx67  = intra_scr->dx67;
  double *intra_scr_dy67  = intra_scr->dy67;
  double *intra_scr_dz67  = intra_scr->dz67;
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
  double *tors_eq_pow     = tors->eq_pow;
  double *tors_c_0        = tors->c_0;
  double *tors_c_1        = tors->c_1;
  double *tors_c_2        = tors->c_2;
  double *tors_c_3        = tors->c_3;
  double *tors_c_4        = tors->c_4;
  double *tors_c_5        = tors->c_5;
  double *tors_c_6        = tors->c_6;
  double *tors_s_0        = tors->s_0;
  double *tors_s_1        = tors->s_1;
  double *tors_s_2        = tors->s_2;
  double *tors_s_3        = tors->s_3;
  double *tors_s_4        = tors->s_4;
  double *tors_s_5        = tors->s_5;
  double *tors_s_6        = tors->s_6;
  double *tors_dc_0       = tors->dc_0;
  double *tors_dc_1       = tors->dc_1;
  double *tors_dc_2       = tors->dc_2;
  double *tors_dc_3       = tors->dc_3;
  double *tors_dc_4       = tors->dc_4;
  double *tors_dc_5       = tors->dc_5;
  double *tors_dc_6       = tors->dc_6;
  double *tors_ds_0       = tors->ds_0;
  double *tors_ds_1       = tors->ds_1;
  double *tors_ds_2       = tors->ds_2;
  double *tors_ds_3       = tors->ds_3;
  double *tors_ds_4       = tors->ds_4;
  double *tors_ds_5       = tors->ds_5;
  double *tors_ds_6       = tors->ds_6;
  double *xmod            = clatoms_info->xmod;
  double *ymod            = clatoms_info->ymod;
  double *zmod            = clatoms_info->zmod;
  double wfor             = (intra_scr->wght_tra);

  double *ptens_pvten_tmp = ptens->pvten_tmp;
  double *pvten           = ptens->pvten;
  double *pvten_tot       = ptens->pvten_tot;

  int *tors_j1_pow        = tors->j1_pow;
  int *tors_j2_pow        = tors->j2_pow;
  int *tors_j3_pow        = tors->j3_pow;
  int *tors_j4_pow        = tors->j4_pow;
  int *tors_jtyp_pow      = tors->jtyp_pow;
  int *iblock_size        = tors->iblock_pow_size;
  int *iblock_conflict_1  = tors->iblock_pow_conflict_1;
  int *iblock_conflict_2  = tors->iblock_pow_conflict_2;
  int *iblock_conflict_3  = tors->iblock_pow_conflict_3;
  int *iblock_conflict_4  = tors->iblock_pow_conflict_4;
  int np_forc             = class_comm_forc_pkg->num_proc;
  int myid_forc           = class_comm_forc_pkg->myid;
  int intra_perds         = cell->intra_perds;
  int iperd               = cell->iperd;
  int pi_beads            = clatoms_info->pi_beads;

  int nlen_use = (intra_scr->nlen);
  int ntot     = tors->npow;


  /*=======================================================================*/
  /* 0) Lower cutoff for sin(theta)                                        */

  seps = 1.0e-8;
  for(i=1;i<=9;i++){(ptens_pvten_tmp)[i]=0;}

  /*=======================================================================*/
  /* I) loop over all the tors in steps of nlen to save memory             */

  for(ibig=1;ibig <= ntot;ibig += nlen_use) {
    /*---------------------------------------------------------------------*/
    /*  A) Offsets to save some memory by doing only nlen bonds at a time  */
    
    ibig1 = ibig-1;
    ioff = -ibig1;
    nlen_now = nlen_use;
    iend = MIN(ntot,ibig1+nlen_now);
    nnow = iend-ibig1;

    /*---------------------------------------------------------------------*/
    /*  B) Gather positions                                                */
    for(itoff=1,itors=ibig;itors <= iend; ++itors,++itoff) {       
      ktemp1 = tors_j1_pow[itors];
      ktemp2 = tors_j2_pow[itors];
      ktemp3 = tors_j3_pow[itors];
      ktemp4 = tors_j4_pow[itors];
      x1 = clatoms_x[ktemp1];
      y1 = clatoms_y[ktemp1];
      z1 = clatoms_z[ktemp1];
      x2 = clatoms_x[ktemp2];
      y2 = clatoms_y[ktemp2];
      z2 = clatoms_z[ktemp2];
      x3 = clatoms_x[ktemp3];
      y3 = clatoms_y[ktemp3];
      z3 = clatoms_z[ktemp3];
      x4 = clatoms_x[ktemp4];
      y4 = clatoms_y[ktemp4];
      z4 = clatoms_z[ktemp4];
      intra_scr_dx12[itoff] = x1-x2;
      intra_scr_dy12[itoff] = y1-y2;
      intra_scr_dz12[itoff] = z1-z2;
      intra_scr_dx13[itoff] = x1-x3;
      intra_scr_dy13[itoff] = y1-y3;
      intra_scr_dz13[itoff] = z1-z3;
      intra_scr_dx42[itoff] = x4-x2;
      intra_scr_dy42[itoff] = y4-y2;
      intra_scr_dz42[itoff] = z4-z2;
      intra_scr_dx43[itoff] = x4-x3;
      intra_scr_dy43[itoff] = y4-y3;
      intra_scr_dz43[itoff] = z4-z3;
      intra_scr_dx23[itoff] = x2-x3;
      intra_scr_dy23[itoff] = y2-y3;
      intra_scr_dz23[itoff] = z2-z3;
    }/*endfor*/

  if(intra_perds==1&&pi_beads>1){     
    for(itoff=1,itors=ibig;itors <= iend; ++itors,++itoff) {       
      ktemp5 = tors_j1_pow[itors];
      ktemp6 = tors_j2_pow[itors];
      ktemp7 = tors_j3_pow[itors];
      ktemp8 = tors_j4_pow[itors];
      x5 = xmod[ktemp5];
      y5 = ymod[ktemp5];
      z5 = zmod[ktemp5];
      x6 = xmod[ktemp6];
      y6 = ymod[ktemp6];
      z6 = zmod[ktemp6];
      x7 = xmod[ktemp7];
      y7 = ymod[ktemp7];
      z7 = zmod[ktemp7];
      x8 = xmod[ktemp8];
      y8 = ymod[ktemp8];
      z8 = zmod[ktemp8];
      intra_scr_dx56[itoff] = x5-x6;
      intra_scr_dy56[itoff] = y5-y6;
      intra_scr_dz56[itoff] = z5-z6;
      intra_scr_dx57[itoff] = x5-x7;
      intra_scr_dy57[itoff] = y5-y7;
      intra_scr_dz57[itoff] = z5-z7;
      intra_scr_dx86[itoff] = x8-x6;
      intra_scr_dy86[itoff] = y8-y6;
      intra_scr_dz86[itoff] = z8-z6;
      intra_scr_dx87[itoff] = x8-x7;
      intra_scr_dy87[itoff] = y8-y7;
      intra_scr_dz87[itoff] = z8-z7;
      intra_scr_dx67[itoff] = x6-x7;
      intra_scr_dy67[itoff] = y6-y7;
      intra_scr_dz67[itoff] = z6-z7;
    }/*endfor*/
  }/*endif*/

    /*----------------------------------------------------------------------*/
    /*  E) Periodic boundary conditions                                     */

    if(intra_perds == 1) {
      if(pi_beads==1){     
        period(nnow,intra_scr_dx12,intra_scr_dy12,intra_scr_dz12,cell);
        period(nnow,intra_scr_dx13,intra_scr_dy13,intra_scr_dz13,cell);
        period(nnow,intra_scr_dx42,intra_scr_dy42,intra_scr_dz42,cell);
        period(nnow,intra_scr_dx43,intra_scr_dy43,intra_scr_dz43,cell);
        period(nnow,intra_scr_dx23,intra_scr_dy23,intra_scr_dz23,cell);
      }else{
        period_pimd(nnow,intra_scr_dx12,intra_scr_dy12,intra_scr_dz12,
                         intra_scr_dx56,intra_scr_dy56,intra_scr_dz56,cell);
        period_pimd(nnow,intra_scr_dx13,intra_scr_dy13,intra_scr_dz13,
                         intra_scr_dx57,intra_scr_dy57,intra_scr_dz57,cell);
        period_pimd(nnow,intra_scr_dx42,intra_scr_dy42,intra_scr_dz42,
                         intra_scr_dx86,intra_scr_dy86,intra_scr_dz86,cell);
        period_pimd(nnow,intra_scr_dx43,intra_scr_dy43,intra_scr_dz43,
                         intra_scr_dx87,intra_scr_dy87,intra_scr_dz87,cell);
        period_pimd(nnow,intra_scr_dx23,intra_scr_dy23,intra_scr_dz23,
                         intra_scr_dx67,intra_scr_dy67,intra_scr_dz67,cell);

      }
    }

    for(itors=1,iboff=ibig;itors <= nnow; ++itors,++iboff) {
      r122   = ((intra_scr_dx12)[itors]*(intra_scr_dx12)[itors]
                + (intra_scr_dy12)[itors]*(intra_scr_dy12)[itors]
                + (intra_scr_dz12)[itors]*(intra_scr_dz12)[itors]);
      r12    = sqrt(r122);
      r132   = ((intra_scr_dx13)[itors]*(intra_scr_dx13)[itors]
                + (intra_scr_dy13)[itors]*(intra_scr_dy13)[itors]
                + (intra_scr_dz13)[itors]*(intra_scr_dz13)[itors]);
      r13    = sqrt(r132);
      r422   = ((intra_scr_dx42)[itors]*(intra_scr_dx42)[itors]
                + (intra_scr_dy42)[itors]*(intra_scr_dy42)[itors]
                + (intra_scr_dz42)[itors]*(intra_scr_dz42)[itors]);
      r42    = sqrt(r422);
      r432   = ((intra_scr_dx43)[itors]*(intra_scr_dx43)[itors]
                + (intra_scr_dy43)[itors]*(intra_scr_dy43)[itors]
                + (intra_scr_dz43)[itors]*(intra_scr_dz43)[itors]);
      r43    = sqrt(r432);
      r232   = ((intra_scr_dx23)[itors]*(intra_scr_dx23)[itors]
                + (intra_scr_dy23)[itors]*(intra_scr_dy23)[itors]
                + (intra_scr_dz23)[itors]*(intra_scr_dz23)[itors]);
      r23    = sqrt(r232);


      /*------------------------------------------------------------------*/
      /*F) Construct vector in plane of 123 perp to r23 call it(dxp,dyp,dzp) */

      d1323i  = 1.0/((intra_scr_dx13)[itors]*(intra_scr_dx23)[itors]
                + (intra_scr_dy13)[itors]*(intra_scr_dy23)[itors]
                + (intra_scr_dz13)[itors]*(intra_scr_dz23)[itors]);
      palp   =-((intra_scr_dx12)[itors]*(intra_scr_dx23)[itors]
                +  (intra_scr_dy12)[itors]*(intra_scr_dy23)[itors]
                +  (intra_scr_dz12)[itors]*(intra_scr_dz23)[itors])*d1323i;
      palp1  = 1.0+palp;
      dxp    = palp*(intra_scr_dx13)[itors]+(intra_scr_dx12)[itors];
      dyp    = palp*(intra_scr_dy13)[itors]+(intra_scr_dy12)[itors];
      dzp    = palp*(intra_scr_dz13)[itors]+(intra_scr_dz12)[itors];
      rp2i   = 1.0/(dxp*dxp+dyp*dyp+dzp*dzp);
      rpi    = sqrt(rp2i);
      /*-----------------------------------------------------------------*/
      /*  G) Get dot product of this vector (dxp,dyp,dzp) with r13 and r42 */
      
      dp13  = (dxp*(intra_scr_dx13)[itors]+dyp*(intra_scr_dy13)[itors]
               +dzp*(intra_scr_dz13)[itors]);
      dp42  = (dxp*(intra_scr_dx42)[itors]+dyp*(intra_scr_dy42)[itors]
               +dzp*(intra_scr_dz42)[itors]);
      
      /*-----------------------------------------------------------------*/
      /*  H) Get derivative of parameter palp used to construct (dxp,dyp,dzp)*/
      /*     with respect to atoms (1-4)                                     */
      
      fpx1   = -palp1*d1323i*(intra_scr_dx23)[itors];
      fpy1   = -palp1*d1323i*(intra_scr_dy23)[itors];
      fpz1   = -palp1*d1323i*(intra_scr_dz23)[itors];
      fpx2   = ((intra_scr_dx23)[itors]-dxp)*d1323i;
      fpy2   = ((intra_scr_dy23)[itors]-dyp)*d1323i;
      fpz2   = ((intra_scr_dz23)[itors]-dzp)*d1323i;

      /*------------------------------------------------------------------*/
      /*  I) Construct vector in plane of 234 perp to r23 (dxq,dyq,dzq) */
      
      d4223i  = 1.0/((intra_scr_dx42)[itors]*(intra_scr_dx23)[itors]
                + (intra_scr_dy42)[itors]*(intra_scr_dy23)[itors]
                + (intra_scr_dz42)[itors]*(intra_scr_dz23)[itors]);
      qalp   =-((intra_scr_dx43)[itors]*(intra_scr_dx23)[itors]
                +  (intra_scr_dy43)[itors]*(intra_scr_dy23)[itors]
                +  (intra_scr_dz43)[itors]*(intra_scr_dz23)[itors])*d4223i;
      qalp1  = 1.0+qalp;
      dxq    = qalp*(intra_scr_dx42)[itors]+(intra_scr_dx43)[itors];
      dyq    = qalp*(intra_scr_dy42)[itors]+(intra_scr_dy43)[itors];
      dzq    = qalp*(intra_scr_dz42)[itors]+(intra_scr_dz43)[itors];
      rq2i   = 1.0/(dxq*dxq+dyq*dyq+dzq*dzq);
      rqi     = sqrt(rq2i);
      /*------------------------------------------------------------------*/
      /*  J) Get dot product of this vector (dxq,dyq,dzq) with r13 and r42*/
      
      dq13  = (dxq*(intra_scr_dx13)[itors]+dyq*(intra_scr_dy13)[itors]
               +dzq*(intra_scr_dz13)[itors]);
      dq42  = (dxq*(intra_scr_dx42)[itors]+dyq*(intra_scr_dy42)[itors]
               +dzq*(intra_scr_dz42)[itors]);
      /*-------------------------------------------------------------------*/
      /*  K) Get derivative of parameter qalp used to construct (dxq,dyq,dzq)*/
      /*     with respect to atoms (1-4)                                   */
      
      fqx4   = -qalp1*d4223i*(intra_scr_dx23)[itors];
      fqy4   = -qalp1*d4223i*(intra_scr_dy23)[itors];
      fqz4   = -qalp1*d4223i*(intra_scr_dz23)[itors];
      fqx2   = (qalp*(intra_scr_dx23)[itors]-dxq)*d4223i;
      fqy2   = (qalp*(intra_scr_dy23)[itors]-dyq)*d4223i;
      fqz2   = (qalp*(intra_scr_dz23)[itors]-dzq)*d4223i;
      /*--------------------------------------------------------------------*/
      /*  L) Get cosine  and sine of the angle                              */

      rpqi = rpi*rqi;
      cost   = (dxp*dxq+dyp*dyq+dzp*dzq)*rpqi;
      cost   = (cost < 1.0 ? cost:1.0);
      cost   = (cost > -1.0 ? cost:-1.0);
      
      sint = sqrt(1.0 - cost*cost);
      sint = (sint > seps ? sint:seps);


      ktemp = tors_jtyp_pow[iboff];

      /*---------------------------------------------------------------------*/
      /*  M) Calculate the potential energy from cosine and sine power       */
      /*     series using Horner's method                                    */
      
      vtorsc = ((((((tors_c_6[ktemp]
		     *cost + tors_c_5[ktemp])
		    *cost + tors_c_4[ktemp])
		   *cost + tors_c_3[ktemp])
		  *cost + tors_c_2[ktemp])
		 *cost + tors_c_1[ktemp])
		*cost + tors_c_0[ktemp]);

      
      sisum = (((((tors_s_6[ktemp]
		  *sint + tors_s_5[ktemp])
		 *sint + tors_s_4[ktemp])
		*sint + tors_s_3[ktemp])
	       *sint + tors_s_2[ktemp])
	      *sint);
      
      vtorss = sisum*cost + (tors_s_1)[ktemp]*sint;
      
      
      intra_scr_vpot[itors] = vtorsc + vtorss;

      /*-------------------------------------------------------------------*/
      /*  N) Get the force on the atoms using the chain rule */
      
      
      dvtorsc = (((((tors_dc_6[ktemp]
		     *cost + tors_dc_5[ktemp])
		    *cost + tors_dc_4[ktemp])
		   *cost + tors_dc_3[ktemp])
		  *cost + tors_dc_2[ktemp])
		 *cost + tors_dc_1[ktemp]);

      dvtorss =  ((((tors_ds_6[ktemp]
		     *sint + tors_ds_5[ktemp])
		    *sint + tors_ds_4[ktemp])
		   *sint + tors_ds_3[ktemp])
		  *sint + tors_ds_2[ktemp]);


      dvtorss = dvtorss*cost + (tors_ds_1)[ktemp];
      
      pre     = -(dvtorsc + sisum - dvtorss*cost/sint);


      costrp2 = cost*rp2i;
      costrq2 = cost*rq2i;
      dp13crp2 = dp13*costrp2;
      dq42crq2 = dq42*costrq2;
      dp13rpq = dp13*rpqi;
      dq13rpq = dq13*rpqi;
      dp42rpq = dp42*rpqi;
      dq42rpq = dq42*rpqi;
      temp1p  = (dq13rpq-dp13crp2)*pre;
      temp1q  = (dp42rpq-dq42crq2)*pre;
      temp2p  = (costrp2*palp1)*pre;
      temp2q  = (costrq2*qalp1)*pre;
      temp3q  = (rpqi*palp1)*pre;  
      temp3p  = (rpqi*qalp1)*pre;  
      temp4p  = (costrp2-qalp*rpqi)*pre;
      temp4q  = (qalp*costrq2-rpqi)*pre;

      (intra_scr_fx1)[itors] = -dxp*temp2p+dxq*temp3q+temp1p*fpx1;
      (intra_scr_fy1)[itors] = -dyp*temp2p+dyq*temp3q+temp1p*fpy1;
      (intra_scr_fz1)[itors] = -dzp*temp2p+dzq*temp3q+temp1p*fpz1;
      (intra_scr_fx4)[itors] = -dxq*temp2q+dxp*temp3p+temp1q*fqx4;
      (intra_scr_fy4)[itors] = -dyq*temp2q+dyp*temp3p+temp1q*fqy4;
      (intra_scr_fz4)[itors] = -dzq*temp2q+dzp*temp3p+temp1q*fqz4;
      (intra_scr_fx2)[itors] =  dxp*temp4p+dxq*temp4q+temp1p*fpx2+temp1q*fqx2;
      (intra_scr_fy2)[itors] =  dyp*temp4p+dyq*temp4q+temp1p*fpy2+temp1q*fqy2;
      (intra_scr_fz2)[itors] =  dzp*temp4p+dzq*temp4q+temp1p*fpz2+temp1q*fqz2;
      (intra_scr_fx3)[itors] = -(intra_scr_fx1[itors]+intra_scr_fx2[itors]
                                +intra_scr_fx4[itors]);
      (intra_scr_fy3)[itors] = -(intra_scr_fy1[itors]+intra_scr_fy2[itors]
                                +intra_scr_fy4[itors]);
      (intra_scr_fz3)[itors] = -(intra_scr_fz1[itors]+intra_scr_fz2[itors]
                                +intra_scr_fz4[itors]);

    /*--------------------------------------------------------------------*/
    /*  O) Get the pressure tensor                                        */
      
      (intra_scr_p11)[itors] = intra_scr_dx13[itors]*intra_scr_fx1[itors]
                             + intra_scr_dx23[itors]*intra_scr_fx2[itors] 
                             + intra_scr_dx43[itors]*intra_scr_fx4[itors];
      (intra_scr_p22)[itors] = intra_scr_dy13[itors]*intra_scr_fy1[itors]
                             + intra_scr_dy23[itors]*intra_scr_fy2[itors] 
                             + intra_scr_dy43[itors]*intra_scr_fy4[itors] ;
      (intra_scr_p33)[itors] = intra_scr_dz13[itors]*intra_scr_fz1[itors]
                             + intra_scr_dz23[itors]*intra_scr_fz2[itors] 
                             + intra_scr_dz43[itors]*intra_scr_fz4[itors] ;
      (intra_scr_p12)[itors] = intra_scr_dx13[itors]*intra_scr_fy1[itors]
                             + intra_scr_dx23[itors]*intra_scr_fy2[itors] 
                             + intra_scr_dx43[itors]*intra_scr_fy4[itors] ;
      (intra_scr_p21)[itors] = intra_scr_dy13[itors]*intra_scr_fx1[itors]
                             + intra_scr_dy23[itors]*intra_scr_fx2[itors] 
                             + intra_scr_dy43[itors]*intra_scr_fx4[itors] ;
      (intra_scr_p13)[itors] = intra_scr_dx13[itors]*intra_scr_fz1[itors]
                             + intra_scr_dx23[itors]*intra_scr_fz2[itors] 
                             + intra_scr_dx43[itors]*intra_scr_fz4[itors] ;
      (intra_scr_p31)[itors] = intra_scr_dz13[itors]*intra_scr_fx1[itors]
                             + intra_scr_dz23[itors]*intra_scr_fx2[itors] 
                             + intra_scr_dz43[itors]*intra_scr_fx4[itors] ;
      (intra_scr_p23)[itors] = intra_scr_dy13[itors]*intra_scr_fz1[itors]
                             + intra_scr_dy23[itors]*intra_scr_fz2[itors] 
                             + intra_scr_dy43[itors]*intra_scr_fz4[itors] ;
      (intra_scr_p32)[itors] = intra_scr_dz13[itors]*intra_scr_fy1[itors]
                             + intra_scr_dz23[itors]*intra_scr_fy2[itors] 
                             + intra_scr_dz43[itors]*intra_scr_fy4[itors] ; 
    }/*endfor itors*/

/*--------------------------------------------------------------------------*/
/*  P) Sum the potential energy                                             */

          for(itors=1;itors <= nnow; ++itors){
              (*vtorst) += (intra_scr_vpot)[itors];
          /*endfor*/}
/*--------------------------------------------------------------------------*/
/*  Q) Sum the pressure tensor                                              */

          if(iperd == 2 || iperd == 3) {
              for(itors=1;itors <= nnow; ++itors){
                 (ptens_pvten_tmp)[1] += (intra_scr_p11)[itors];
                 (ptens_pvten_tmp)[5] += (intra_scr_p22)[itors];
                 (ptens_pvten_tmp)[9] += (intra_scr_p33)[itors];
                 (ptens_pvten_tmp)[2] += (intra_scr_p12)[itors];
                 (ptens_pvten_tmp)[4] += (intra_scr_p21)[itors];
              /*endfor*/}
          /*endif*/}
          if(iperd == 3) {
              for(itors=1;itors <= nnow; ++itors){
                 (ptens_pvten_tmp)[3] += (intra_scr_p13)[itors];
                 (ptens_pvten_tmp)[7] += (intra_scr_p31)[itors];
                 (ptens_pvten_tmp)[6] += (intra_scr_p23)[itors];
                 (ptens_pvten_tmp)[8] += (intra_scr_p32)[itors];
             /*endfor*/}
          /*endif*/}

/*--------------------------------------------------------------------------*/
/*  R) Scatter the forces                                                   */
      if(iver_get==1){

        for(itoff=1,itors=ibig;itors <= iend; ++itors,++itoff) {       
         ktemp = tors_j1_pow[itors];
         clatoms_fxt[ktemp]  +=  intra_scr_fx1[itoff];
         clatoms_fyt[ktemp]  +=  intra_scr_fy1[itoff];
         clatoms_fzt[ktemp]  +=  intra_scr_fz1[itoff];
        }/*endfor*/

        for(itoff=1,itors=ibig;itors <= iend; ++itors,++itoff) {       
         ktemp = tors_j2_pow[itors];
         clatoms_fxt[ktemp]  +=  intra_scr_fx2[itoff];
         clatoms_fyt[ktemp]  +=  intra_scr_fy2[itoff];
         clatoms_fzt[ktemp]  +=  intra_scr_fz2[itoff];
        }/*endfor*/

        for(itoff=1,itors=ibig;itors <= iend; ++itors,++itoff) {       
         ktemp = tors_j3_pow[itors];
         clatoms_fxt[ktemp]  +=  intra_scr_fx3[itoff];
         clatoms_fyt[ktemp]  +=  intra_scr_fy3[itoff];
         clatoms_fzt[ktemp]  +=  intra_scr_fz3[itoff];
        }/*endfor*/

        for(itoff=1,itors=ibig;itors <= iend; ++itors,++itoff) {       
         ktemp = tors_j4_pow[itors];
         clatoms_fxt[ktemp]  +=  intra_scr_fx4[itoff];
         clatoms_fyt[ktemp]  +=  intra_scr_fy4[itoff];
         clatoms_fzt[ktemp]  +=  intra_scr_fz4[itoff];
        }/*endfor*/

     }/*endif*/

     for(itoff=1;itoff <= iend+ioff; ++itoff) {
          intra_scr_fx1[itoff] *= wfor;
          intra_scr_fy1[itoff] *= wfor;
          intra_scr_fz1[itoff] *= wfor;
          intra_scr_fx2[itoff] *= wfor;
          intra_scr_fy2[itoff] *= wfor;
          intra_scr_fz2[itoff] *= wfor;
          intra_scr_fx3[itoff] *= wfor;
          intra_scr_fy3[itoff] *= wfor;
          intra_scr_fz3[itoff] *= wfor;
          intra_scr_fx4[itoff] *= wfor;
          intra_scr_fy4[itoff] *= wfor;
          intra_scr_fz4[itoff] *= wfor;
      }/*endfor*/

      for(itoff=1,itors=ibig;itors <= iend; ++itors,++itoff) {       
        ktemp = tors_j1_pow[itors];
        clatoms_fx[ktemp]  +=  intra_scr_fx1[itoff];
        clatoms_fy[ktemp]  +=  intra_scr_fy1[itoff];
        clatoms_fz[ktemp]  +=  intra_scr_fz1[itoff];
      }/*endfor*/

      for(itoff=1,itors=ibig;itors <= iend; ++itors,++itoff) {       
        ktemp = tors_j2_pow[itors];
        clatoms_fx[ktemp]  +=  intra_scr_fx2[itoff];
        clatoms_fy[ktemp]  +=  intra_scr_fy2[itoff];
        clatoms_fz[ktemp]  +=  intra_scr_fz2[itoff];
      }/*endfor*/

      for(itoff=1,itors=ibig;itors <= iend; ++itors,++itoff) {       
        ktemp = tors_j3_pow[itors];
        clatoms_fx[ktemp]  +=  intra_scr_fx3[itoff];
        clatoms_fy[ktemp]  +=  intra_scr_fy3[itoff];
        clatoms_fz[ktemp]  +=  intra_scr_fz3[itoff];
      }/*endfor*/

      for(itoff=1,itors=ibig;itors <= iend; ++itors,++itoff) {       
        ktemp = tors_j4_pow[itors];
        clatoms_fx[ktemp]  +=  intra_scr_fx4[itoff];
        clatoms_fy[ktemp]  +=  intra_scr_fy4[itoff];
        clatoms_fz[ktemp]  +=  intra_scr_fz4[itoff];
      }/*endfor*/

 }/*endfor ibig*/

/*==========================================================================*/
/* Increment the Pressure tensor */

  for(i=1;i<=9;i++){
       (pvten)[i]     += (ptens_pvten_tmp)[i]*wfor;
  }/*endfor*/

  if(iget_pv_real_inter==1){    
   for(i=1;i<=9;i++){
       (pvten_tot)[i] += (ptens_pvten_tmp)[i];
   }/*endfor*/
  }/*endif*/

/*--------------------------------------------------------------------------*/
   /* end routine */}
/*==========================================================================*/







/*==========================================================================*/
/*               Header:                                                    */
/*==========================================================================*/

void tors_free(CLATOMS_INFO *clatoms_info,CLATOMS_POS *clatoms_pos,
               TORS_FREE *tors_free,CELL *cell,
               INTRA_SCR *intra_scr,PTENS *ptens, double *vtorst,
               ENERGY_CTRL *energy_ctrl,int np_forc)

/*==========================================================================*/
/*               Begin subprogram:                                          */
    {/*begin routine*/
/*==========================================================================*/
/*               Local variable declarations                                */
   int ibin,jbin,nnow,i;  
   double r122,r12,sign_theta;
   double r132,r13,arg,tpi;
   double r232,r23;
   double r422,r42;
   double r432,r43;
                                                 /* Distances              */
   double dxp,dyp,dzp,rp2i,rpi,rpqi;
   double dxq,dyq,dzq,rq2i,rqi;     
   double dp13,dp42,dq13,dq42;
   double cost,sint,seps,pre;          
   double costrp2,costrq2,dp13crp2,dq42crq2;
   double dp13rpq,dq13rpq,dp42rpq,dq42rpq;
   double d1323i,palp,palp1;
   double d4223i,qalp,qalp1; 
   double fpx1,fpy1,fpz1;
   double fpx2,fpy2,fpz2;  
   double fqx4,fqy4,fqz4;
   double fqx2,fqy2,fqz2;  
   double apow,apow1,theta0,theta,dvtors;
                                                /* Force temporaries       */
   double temp1p,temp1q,temp2q,temp2p;
   double temp3p,temp3q,temp4q,temp4p;

   double xcrs_crs,ycrs_crs,zcrs_crs,ssdot;
   double theta_now[4];

/*-------------------------------------------------------------------------*/
/*               Local pointer declarations                                */

   double *xmod            = clatoms_info->xmod;
   double *ymod            = clatoms_info->ymod;
   double *zmod            = clatoms_info->zmod;
   int pi_beads            = clatoms_info->pi_beads; 

   int num                = tors_free->num;
   int *tors_free_j1      = tors_free->j1;
   int *tors_free_j2      = tors_free->j2;
   int *tors_free_j3      = tors_free->j3;
   int *tors_free_j4      = tors_free->j4;
   int tors_free_nhist    = tors_free->nhist;

   double *tors_free_eq  = tors_free->eq;
   double tors_free_npow = tors_free->npow;
   double tors_free_fk   = tors_free->fk;
   double tors_free_del   = tors_free->del;
   double *tors_free_hist = tors_free->hist;
   double **tors_free_hist_2d = tors_free->hist_2d;

   double *clatoms_pos_x  = clatoms_pos->x;
   double *clatoms_pos_y  = clatoms_pos->y;
   double *clatoms_pos_z  = clatoms_pos->z;
   double *clatoms_pos_fx = clatoms_pos->fx;
   double *clatoms_pos_fy = clatoms_pos->fy;
   double *clatoms_pos_fz = clatoms_pos->fz;

   double *intra_scr_x1 =  intra_scr->x1;
   double *intra_scr_y1 =  intra_scr->y1; 
   double *intra_scr_z1 =  intra_scr->z1; 
   double *intra_scr_x2 =  intra_scr->x2; 
   double *intra_scr_y2 =  intra_scr->y2; 
   double *intra_scr_z2 =  intra_scr->z2; 
   double *intra_scr_x3 =  intra_scr->x3; 
   double *intra_scr_y3 =  intra_scr->y3; 
   double *intra_scr_z3 =  intra_scr->z3;  
   double *intra_scr_x4 =  intra_scr->x4;  
   double *intra_scr_y4 =  intra_scr->y4;  
   double *intra_scr_z4 =  intra_scr->z4;  

   double *intra_scr_fx1 =  intra_scr->fx1;
   double *intra_scr_fy1 =  intra_scr->fy1; 
   double *intra_scr_fz1 =  intra_scr->fz1; 
   double *intra_scr_fx2 =  intra_scr->fx2; 
   double *intra_scr_fy2 =  intra_scr->fy2; 
   double *intra_scr_fz2 =  intra_scr->fz2; 
   double *intra_scr_fx3 =  intra_scr->fx3; 
   double *intra_scr_fy3 =  intra_scr->fy3; 
   double *intra_scr_fz3 =  intra_scr->fz3;  
   double *intra_scr_fx4 =  intra_scr->fx4;  
   double *intra_scr_fy4 =  intra_scr->fy4;  
   double *intra_scr_fz4 =  intra_scr->fz4;  

   double *intra_scr_x5=  intra_scr->x5;
   double *intra_scr_y5=  intra_scr->y5;
   double *intra_scr_z5=  intra_scr->z5;
   double *intra_scr_x6=  intra_scr->x6;
   double *intra_scr_y6=  intra_scr->y6;
   double *intra_scr_z6=  intra_scr->z6;
   double *intra_scr_x7=  intra_scr->x7;
   double *intra_scr_y7=  intra_scr->y7;
   double *intra_scr_z7=  intra_scr->z7;
   double *intra_scr_x8=  intra_scr->x8;
   double *intra_scr_y8=  intra_scr->y8; 
   double *intra_scr_z8=  intra_scr->z8;

   double *intra_scr_dx12 = intra_scr->dx12;
   double *intra_scr_dy12 = intra_scr->dy12;
   double *intra_scr_dz12 = intra_scr->dz12;
   double *intra_scr_dx13 = intra_scr->dx13;
   double *intra_scr_dy13 = intra_scr->dy13;
   double *intra_scr_dz13 = intra_scr->dz13;
   double *intra_scr_dx42 = intra_scr->dx42;
   double *intra_scr_dy42 = intra_scr->dy42;
   double *intra_scr_dz42 = intra_scr->dz42;
   double *intra_scr_dx43 = intra_scr->dx43;
   double *intra_scr_dy43 = intra_scr->dy43;
   double *intra_scr_dz43 = intra_scr->dz43;
   double *intra_scr_dx23 = intra_scr->dx23;
   double *intra_scr_dy23 = intra_scr->dy23;
   double *intra_scr_dz23 = intra_scr->dz23;

   double *intra_scr_dx56 = intra_scr->dx56;
   double *intra_scr_dy56 = intra_scr->dy56;
   double *intra_scr_dz56 = intra_scr->dz56;
   double *intra_scr_dx57 = intra_scr->dx57;
   double *intra_scr_dy57 = intra_scr->dy57;
   double *intra_scr_dz57 = intra_scr->dz57;
   double *intra_scr_dx86 = intra_scr->dx86;
   double *intra_scr_dy86 = intra_scr->dy86;
   double *intra_scr_dz86 = intra_scr->dz86;
   double *intra_scr_dx87 = intra_scr->dx87;
   double *intra_scr_dy87 = intra_scr->dy87;
   double *intra_scr_dz87 = intra_scr->dz87;
   double *intra_scr_dx67 = intra_scr->dx67;
   double *intra_scr_dy67 = intra_scr->dy67;
   double *intra_scr_dz67 = intra_scr->dz67;

   double *intra_scr_p11 = intra_scr->p11;
   double *intra_scr_p22 = intra_scr->p22;
   double *intra_scr_p33 = intra_scr->p33;
   double *intra_scr_p12 = intra_scr->p12;
   double *intra_scr_p21 = intra_scr->p21;
   double *intra_scr_p31 = intra_scr->p31;
   double *intra_scr_p13 = intra_scr->p13;
   double *intra_scr_p32 = intra_scr->p32;
   double *intra_scr_p23 = intra_scr->p23;

   double wfor           = intra_scr->wght_tra;

   double *ptens_pvten     = ptens->pvten;
   double *ptens_pvten_tot = ptens->pvten_tot;

   int iperd              = cell->iperd;
   int intra_perds        = cell->intra_perds;

   int iget_pv_real_inter = energy_ctrl->iget_pv_real_inter;
   int iget_full_inter    = energy_ctrl->iget_full_inter;
 
/*==========================================================================*/
/* 0) Lower cutoff for sin(theta)                                           */
 
  /* printf("In the tors free routine\n"); */

   seps = 1.0e-8;
   tpi  = 2.0*M_PI;

/* I) loop over all the tors in steps of nlen to save memory                */

   (*vtorst) = 0.0;
/*--------------------------------------------------------------------------*/
/*  B) Gather positions                                                     */

  for(i=1;i<=num;i++){
   intra_scr_x1[i] = clatoms_pos_x[tors_free_j1[i]];
   intra_scr_y1[i] = clatoms_pos_y[tors_free_j1[i]];
   intra_scr_z1[i] = clatoms_pos_z[tors_free_j1[i]];
   intra_scr_x2[i] = clatoms_pos_x[tors_free_j2[i]];
   intra_scr_y2[i] = clatoms_pos_y[tors_free_j2[i]];
   intra_scr_z2[i] = clatoms_pos_z[tors_free_j2[i]];
   intra_scr_x3[i] = clatoms_pos_x[tors_free_j3[i]];
   intra_scr_y3[i] = clatoms_pos_y[tors_free_j3[i]];
   intra_scr_z3[i] = clatoms_pos_z[tors_free_j3[i]];
   intra_scr_x4[i] = clatoms_pos_x[tors_free_j4[i]];
   intra_scr_y4[i] = clatoms_pos_y[tors_free_j4[i]];
   intra_scr_z4[i] = clatoms_pos_z[tors_free_j4[i]];
  }/*endfor*/

  if(intra_perds==1&&pi_beads>1){     
    for(i=1;i<=num;i++){
     intra_scr_x5[i] = xmod[tors_free_j1[i]];
     intra_scr_y5[i] = ymod[tors_free_j1[i]];
     intra_scr_z5[i] = zmod[tors_free_j1[i]];
     intra_scr_x6[i] = xmod[tors_free_j2[i]];
     intra_scr_y6[i] = ymod[tors_free_j2[i]];
     intra_scr_z6[i] = zmod[tors_free_j2[i]];
     intra_scr_x7[i] = xmod[tors_free_j3[i]];
     intra_scr_y7[i] = ymod[tors_free_j3[i]];
     intra_scr_z7[i] = zmod[tors_free_j3[i]];
     intra_scr_x8[i] = xmod[tors_free_j4[i]];
     intra_scr_y8[i] = ymod[tors_free_j4[i]];
     intra_scr_z8[i] = zmod[tors_free_j4[i]];
    }/*endfor*/
   }/*endif*/
 
/*--------------------------------------------------------------------------*/
/*  D) Calculate the basis vectors (r12,r13,r43,r42)                       */

   for(i=1;i<=num;i++){
    intra_scr_dx12[i] = intra_scr_x1[i] - intra_scr_x2[i];
    intra_scr_dy12[i] = intra_scr_y1[i] - intra_scr_y2[i];
    intra_scr_dz12[i] = intra_scr_z1[i] - intra_scr_z2[i];
    intra_scr_dx13[i] = intra_scr_x1[i] - intra_scr_x3[i];
    intra_scr_dy13[i] = intra_scr_y1[i] - intra_scr_y3[i];
    intra_scr_dz13[i] = intra_scr_z1[i] - intra_scr_z3[i];
    intra_scr_dx42[i] = intra_scr_x4[i] - intra_scr_x2[i];
    intra_scr_dy42[i] = intra_scr_y4[i] - intra_scr_y2[i];
    intra_scr_dz42[i] = intra_scr_z4[i] - intra_scr_z2[i];
    intra_scr_dx43[i] = intra_scr_x4[i] - intra_scr_x3[i];
    intra_scr_dy43[i] = intra_scr_y4[i] - intra_scr_y3[i];
    intra_scr_dz43[i] = intra_scr_z4[i] - intra_scr_z3[i];
    intra_scr_dx23[i] = intra_scr_x2[i] - intra_scr_x3[i];
    intra_scr_dy23[i] = intra_scr_y2[i] - intra_scr_y3[i];
    intra_scr_dz23[i] = intra_scr_z2[i] - intra_scr_z3[i];
   }/*endfor*/
 
   if(intra_perds==1&&pi_beads>1){     
    for(i=1;i<=num;i++){
     intra_scr_dx56[i] = intra_scr_x5[i] - intra_scr_x6[i];
     intra_scr_dy56[i] = intra_scr_y5[i] - intra_scr_y6[i];
     intra_scr_dz56[i] = intra_scr_z5[i] - intra_scr_z6[i];
     intra_scr_dx57[i] = intra_scr_x5[i] - intra_scr_x7[i];
     intra_scr_dy57[i] = intra_scr_y5[i] - intra_scr_y7[i];
     intra_scr_dz57[i] = intra_scr_z5[i] - intra_scr_z7[i];
     intra_scr_dx86[i] = intra_scr_x8[i] - intra_scr_x6[i];
     intra_scr_dy86[i] = intra_scr_y8[i] - intra_scr_y6[i];
     intra_scr_dz86[i] = intra_scr_z8[i] - intra_scr_z6[i];
     intra_scr_dx87[i] = intra_scr_x8[i] - intra_scr_x7[i];
     intra_scr_dy87[i] = intra_scr_y8[i] - intra_scr_y7[i];
     intra_scr_dz87[i] = intra_scr_z8[i] - intra_scr_z7[i];
     intra_scr_dx67[i] = intra_scr_x6[i] - intra_scr_x7[i];
     intra_scr_dy67[i] = intra_scr_y6[i] - intra_scr_y7[i];
     intra_scr_dz67[i] = intra_scr_z6[i] - intra_scr_z7[i];
    }/*endfor*/
   }/*endif*/
/*--------------------------------------------------------------------------*/
/*  E) Periodic boundary conditions                                         */

   nnow = num;
   if(pi_beads==1){     
        period(nnow,intra_scr_dx12,intra_scr_dy12,intra_scr_dz12,cell);
        period(nnow,intra_scr_dx13,intra_scr_dy13,intra_scr_dz13,cell);
        period(nnow,intra_scr_dx42,intra_scr_dy42,intra_scr_dz42,cell);
        period(nnow,intra_scr_dx43,intra_scr_dy43,intra_scr_dz43,cell);
        period(nnow,intra_scr_dx23,intra_scr_dy23,intra_scr_dz23,cell);
   }else{
        period_pimd(nnow,intra_scr_dx12,intra_scr_dy12,intra_scr_dz12,
                         intra_scr_dx56,intra_scr_dy56,intra_scr_dz56,cell);
        period_pimd(nnow,intra_scr_dx13,intra_scr_dy13,intra_scr_dz13,
                         intra_scr_dx57,intra_scr_dy57,intra_scr_dz57,cell);
        period_pimd(nnow,intra_scr_dx42,intra_scr_dy42,intra_scr_dz42,
                         intra_scr_dx86,intra_scr_dy86,intra_scr_dz86,cell);
        period_pimd(nnow,intra_scr_dx43,intra_scr_dy43,intra_scr_dz43,
                         intra_scr_dx87,intra_scr_dy87,intra_scr_dz87,cell);
        period_pimd(nnow,intra_scr_dx23,intra_scr_dy23,intra_scr_dz23,
                         intra_scr_dx67,intra_scr_dy67,intra_scr_dz67,cell);

   }/*endif*/

  for(i=1;i<=num;i++){
   r122   = (intra_scr_dx12[i]*intra_scr_dx12[i]
           +  intra_scr_dy12[i]*intra_scr_dy12[i]
           +  intra_scr_dz12[i]*intra_scr_dz12[i]);
   r12    = sqrt(r122);
   r132   = (intra_scr_dx13[i]*intra_scr_dx13[i]
           +  intra_scr_dy13[i]*intra_scr_dy13[i]
           +  intra_scr_dz13[i]*intra_scr_dz13[i]);
   r13    = sqrt(r132); 
   r422   = (intra_scr_dx42[i]*intra_scr_dx42[i]
           +  intra_scr_dy42[i]*intra_scr_dy42[i]
           +  intra_scr_dz42[i]*intra_scr_dz42[i]);
   r42    = sqrt(r422);
   r432   = (intra_scr_dx43[i]*intra_scr_dx43[i]
           +  intra_scr_dy43[i]*intra_scr_dy43[i]
           +  intra_scr_dz43[i]*intra_scr_dz43[i]);
   r43    = sqrt(r432);
   r232   = (intra_scr_dx23[i]*intra_scr_dx23[i]
           +  intra_scr_dy23[i]*(intra_scr_dy23)[i]
           +  intra_scr_dz23[i]*intra_scr_dz23[i]);
   r23    = sqrt(r232);
 
/*------------------------------------------------------------------*/
/*F) Construct vector in plane of 123 perp to r23 call it(dxp,dyp,dzp) */

   d1323i  = 1.0/(intra_scr_dx13[i]*intra_scr_dx23[i]
                + intra_scr_dy13[i]*intra_scr_dy23[i]
                + intra_scr_dz13[i]*intra_scr_dz23[i]);
   palp   =-(intra_scr_dx12[i]*intra_scr_dx23[i]
          +  intra_scr_dy12[i]*intra_scr_dy23[i]
          +  intra_scr_dz12[i]*intra_scr_dz23[i])*d1323i;
   palp1  = 1.0+palp;
   dxp    = palp*intra_scr_dx13[i]+intra_scr_dx12[i];
   dyp    = palp*intra_scr_dy13[i]+intra_scr_dy12[i];
   dzp    = palp*intra_scr_dz13[i]+intra_scr_dz12[i];
   rp2i   = 1.0/(dxp*dxp+dyp*dyp+dzp*dzp);
   rpi    = sqrt(rp2i);

/*-----------------------------------------------------------------*/
/*  G) Get dot product of this vector (dxp,dyp,dzp) with r13 and r42 */
      
   dp13  = (dxp*intra_scr_dx13[i]+dyp*intra_scr_dy13[i]
           +dzp*intra_scr_dz13[i]);
   dp42  = (dxp*intra_scr_dx42[i]+dyp*intra_scr_dy42[i]
           +dzp*intra_scr_dz42[i]);
      
/*-----------------------------------------------------------------*/
/*  H) Get derivative of parameter palp used to construct (dxp,dyp,dzp)*/
/*     with respect to atoms (1-4)                                     */
      
   fpx1   = -palp1*d1323i*intra_scr_dx23[i];
   fpy1   = -palp1*d1323i*intra_scr_dy23[i];
   fpz1   = -palp1*d1323i*intra_scr_dz23[i];
   fpx2   = (intra_scr_dx23[i]-dxp)*d1323i;
   fpy2   = (intra_scr_dy23[i]-dyp)*d1323i;
   fpz2   = (intra_scr_dz23[i]-dzp)*d1323i;

/*------------------------------------------------------------------*/
/*  I) Construct vector in plane of 234 perp to r23 (dxq,dyq,dzq) */
      
   d4223i  = 1.0/(intra_scr_dx42[i]*intra_scr_dx23[i]
                + intra_scr_dy42[i]*intra_scr_dy23[i]
                + intra_scr_dz42[i]*intra_scr_dz23[i]);
   qalp   =-(intra_scr_dx43[i]*intra_scr_dx23[i]
          +  intra_scr_dy43[i]*intra_scr_dy23[i]
          +  intra_scr_dz43[i]*intra_scr_dz23[i])*d4223i;
   qalp1  = 1.0+qalp;
   dxq    = qalp*intra_scr_dx42[i]+intra_scr_dx43[i];
   dyq    = qalp*intra_scr_dy42[i]+intra_scr_dy43[i];
   dzq    = qalp*intra_scr_dz42[i]+intra_scr_dz43[i];
   rq2i   = 1.0/(dxq*dxq+dyq*dyq+dzq*dzq);
   rqi    = sqrt(rq2i);
/*------------------------------------------------------------------*/
/*  J) Get dot product of this vector (dxq,dyq,dzq) with r13 and r42*/
      
   dq13  = (dxq*intra_scr_dx13[i]+dyq*intra_scr_dy13[i]
           +dzq*intra_scr_dz13[i]);
   dq42  = (dxq*intra_scr_dx42[i]+dyq*intra_scr_dy42[i]
           +dzq*intra_scr_dz42[i]);
/*-------------------------------------------------------------------*/
/*  K) Get derivative of parameter qalp used to construct (dxq,dyq,dzq)*/
/*     with respect to atoms (1-4)                                   */
      
   fqx4   = -qalp1*d4223i*intra_scr_dx23[i];
   fqy4   = -qalp1*d4223i*intra_scr_dy23[i];
   fqz4   = -qalp1*d4223i*intra_scr_dz23[i];
   fqx2   = (qalp*intra_scr_dx23[i]-dxq)*d4223i;
   fqy2   = (qalp*intra_scr_dy23[i]-dyq)*d4223i;
   fqz2   = (qalp*intra_scr_dz23[i]-dzq)*d4223i;

/*--------------------------------------------------------------------*/
/*  L) Get cosine  and sine of the angle                              */

   rpqi   = rpi*rqi;

   xcrs_crs = (dyq*dzp - dyp*dzq);
   ycrs_crs = (dzq*dxp - dzp*dxq);
   zcrs_crs = (dxq*dyp - dxp*dyq);
   sint     = (xcrs_crs*intra_scr_dx23[i]
              +ycrs_crs*intra_scr_dy23[i]
              +zcrs_crs*intra_scr_dz23[i])*rpqi/r23;
   sint   = (sint < 1.0 ? sint:1.0);
   sint   = (sint > -1.0 ? sint:-1.0);

   cost   = (dxp*dxq+dyp*dyq+dzp*dzq)*rpqi;
   cost   = (cost < 1.0 ? cost:1.0);
   cost   = (cost > -1.0 ? cost:-1.0);

   theta      = acos(cost);
   sign_theta = (sint >= 0 ? 1.0: -1.0);
   theta     *= sign_theta;
   sint       = (sign_theta*sint > seps ? sint:sign_theta*seps);

   theta_now[i] = theta;
   theta0  = theta - (tors_free_eq[i]);
   arg     = theta0/tpi;
   theta0 -= (tpi*NINT(arg));

/*--------------------------------------------------------------------------*/
/*  M) Calculate the potential energy and bin the stuff                     */

   apow      = (double)(tors_free_npow);
   apow1     = apow-1.0;

   (*vtorst) += (tors_free_fk)*pow(theta0,apow)/apow;

/*--------------------------------------------------------------------------*/
/*  N) Get the force on the atoms using the chain rule and Horner's method  */

   dvtors = (tors_free_fk)*pow(theta0,apow1);
   pre    = dvtors/sint;

   costrp2 = cost*rp2i;
   costrq2 = cost*rq2i;
   dp13crp2 = dp13*costrp2;
   dq42crq2 = dq42*costrq2;
   dp13rpq = dp13*rpqi;
   dq13rpq = dq13*rpqi;
   dp42rpq = dp42*rpqi;
   dq42rpq = dq42*rpqi;
   temp1p  = (dq13rpq-dp13crp2)*pre;
   temp1q  = (dp42rpq-dq42crq2)*pre;
   temp2p  = (costrp2*palp1)*pre;
   temp2q  = (costrq2*qalp1)*pre;
   temp3q  = (rpqi*palp1)*pre;  
   temp3p  = (rpqi*qalp1)*pre;  
   temp4p  = (costrp2-qalp*rpqi)*pre;
   temp4q  = (qalp*costrq2-rpqi)*pre;

   intra_scr_fx1[i] = -dxp*temp2p+dxq*temp3q+temp1p*fpx1;
   intra_scr_fy1[i] = -dyp*temp2p+dyq*temp3q+temp1p*fpy1;
   intra_scr_fz1[i] = -dzp*temp2p+dzq*temp3q+temp1p*fpz1;

   intra_scr_fx4[i] = -dxq*temp2q+dxp*temp3p+temp1q*fqx4;
   intra_scr_fy4[i] = -dyq*temp2q+dyp*temp3p+temp1q*fqy4;
   intra_scr_fz4[i] = -dzq*temp2q+dzp*temp3p+temp1q*fqz4;

   intra_scr_fx2[i] =  dxp*temp4p+dxq*temp4q+temp1p*fpx2+temp1q*fqx2;
   intra_scr_fy2[i] =  dyp*temp4p+dyq*temp4q+temp1p*fpy2+temp1q*fqy2;
   intra_scr_fz2[i] =  dzp*temp4p+dzq*temp4q+temp1p*fpz2+temp1q*fqz2;

   intra_scr_fx3[i] = -(intra_scr_fx1[i]+intra_scr_fx2[i]
                      +intra_scr_fx4[i]);
   intra_scr_fy3[i] = -(intra_scr_fy1[i]+intra_scr_fy2[i]
                      +intra_scr_fy4[i]);
   intra_scr_fz3[i] = -(intra_scr_fz1[i]+intra_scr_fz2[i]
                      +intra_scr_fz4[i]);

/*--------------------------------------------------------------------------*/
/*  O) Get the pressure tensor (this form is ok for periodic systems):      */

   intra_scr_p11[i] = intra_scr_dx13[i]*intra_scr_fx1[i]
                       + intra_scr_dx23[i]*intra_scr_fx2[i] 
                       + intra_scr_dx43[i]*intra_scr_fx4[i];
   intra_scr_p22[i] = intra_scr_dy13[i]*intra_scr_fy1[i]
                       + intra_scr_dy23[i]*intra_scr_fy2[i] 
                       + intra_scr_dy43[i]*intra_scr_fy4[i];
   intra_scr_p33[i] = intra_scr_dz13[i]*intra_scr_fz1[i]
                       + intra_scr_dz23[i]*intra_scr_fz2[i] 
                       + intra_scr_dz43[i]*intra_scr_fz4[i];
   intra_scr_p12[i] = intra_scr_dx13[i]*intra_scr_fy1[i]
                       + intra_scr_dx23[i]*intra_scr_fy2[i] 
                       + intra_scr_dx43[i]*intra_scr_fy4[i];
   intra_scr_p21[i] = intra_scr_dy13[i]*intra_scr_fx1[i]
                       + intra_scr_dy23[i]*intra_scr_fx2[i] 
                       + intra_scr_dy43[i]*intra_scr_fx4[i];
   intra_scr_p13[i] = intra_scr_dx13[i]*intra_scr_fz1[i]
                       + intra_scr_dx23[i]*intra_scr_fz2[i] 
                       + intra_scr_dx43[i]*intra_scr_fz4[i];
   intra_scr_p31[i] = intra_scr_dz13[i]*intra_scr_fx1[i]
                       + intra_scr_dz23[i]*intra_scr_fx2[i]
                       + intra_scr_dz43[i]*intra_scr_fx4[i];
   intra_scr_p23[i] = intra_scr_dy13[i]*intra_scr_fz1[i]
                       + intra_scr_dy23[i]*intra_scr_fz2[i] 
                       + intra_scr_dy43[i]*intra_scr_fz4[i];
   intra_scr_p32[i] = intra_scr_dz13[i]*intra_scr_fy1[i]
                       + intra_scr_dz23[i]*intra_scr_fy2[i] 
                       + intra_scr_dz43[i]*intra_scr_fy4[i]; 
 }/*endfor*/
 
/*--------------------------------------------------------------------------*/
/*  Q) Sum the pressure tensor                                              */

   if(iperd == 2 || iperd == 3) {
    for(i=1;i<=num;i++){
     ptens_pvten[1] += intra_scr_p11[i]*wfor*pi_beads;
     ptens_pvten[5] += intra_scr_p22[i]*wfor*pi_beads;
     ptens_pvten[9] += intra_scr_p33[i]*wfor*pi_beads;
     ptens_pvten[2] += intra_scr_p12[i]*wfor*pi_beads;
     ptens_pvten[4] += intra_scr_p21[i]*wfor*pi_beads;
    }/*endfor*/
    if(iget_pv_real_inter==1){    
     for(i=1;i<=num;i++){
       ptens_pvten_tot[1] += intra_scr_p11[i]*pi_beads;
       ptens_pvten_tot[5] += intra_scr_p22[i]*pi_beads;
       ptens_pvten_tot[9] += intra_scr_p33[i]*pi_beads;
       ptens_pvten_tot[2] += intra_scr_p12[i]*pi_beads;
       ptens_pvten_tot[4] += intra_scr_p21[i]*pi_beads;
     }/*endfor*/
    }/*endif*/

   }/*endif*/
 
   if(iperd == 3) {
    for(i=1;i<=num;i++){
     ptens_pvten[3] += intra_scr_p13[i]*wfor*pi_beads;
     ptens_pvten[7] += intra_scr_p31[i]*wfor*pi_beads;
     ptens_pvten[6] += intra_scr_p23[i]*wfor*pi_beads;
     ptens_pvten[8] += intra_scr_p32[i]*wfor*pi_beads;
    }/*endfor*/

    if(iget_pv_real_inter==1){    
     for(i=1;i<=num;i++){
      ptens_pvten_tot[3] += intra_scr_p13[i]*pi_beads;
      ptens_pvten_tot[7] += intra_scr_p31[i]*pi_beads;
      ptens_pvten_tot[6] += intra_scr_p23[i]*pi_beads;
      ptens_pvten_tot[8] += intra_scr_p32[i]*pi_beads;
     }/*endfor*/
    }/*endif*/

   }/*endif*/
 
/*--------------------------------------------------------------------------*/
/*  R) Scatter the forces                                                   */

  for(i=1;i<=num;i++){
   clatoms_pos_fx[tors_free_j1[i]]  += intra_scr_fx1[i]*wfor*pi_beads;
   clatoms_pos_fy[tors_free_j1[i]]  += intra_scr_fy1[i]*wfor*pi_beads;
   clatoms_pos_fz[tors_free_j1[i]]  += intra_scr_fz1[i]*wfor*pi_beads;

   clatoms_pos_fx[tors_free_j2[i]]  += intra_scr_fx2[i]*wfor*pi_beads;
   clatoms_pos_fy[tors_free_j2[i]]  += intra_scr_fy2[i]*wfor*pi_beads;
   clatoms_pos_fz[tors_free_j2[i]]  += intra_scr_fz2[i]*wfor*pi_beads;

   clatoms_pos_fx[tors_free_j3[i]]  += intra_scr_fx3[i]*wfor*pi_beads;
   clatoms_pos_fy[tors_free_j3[i]]  += intra_scr_fy3[i]*wfor*pi_beads;
   clatoms_pos_fz[tors_free_j3[i]]  += intra_scr_fz3[i]*wfor*pi_beads;

   clatoms_pos_fx[tors_free_j4[i]]  += intra_scr_fx4[i]*wfor*pi_beads;
   clatoms_pos_fy[tors_free_j4[i]]  += intra_scr_fy4[i]*wfor*pi_beads;
   clatoms_pos_fz[tors_free_j4[i]]  += intra_scr_fz4[i]*wfor*pi_beads;
 }/*endfor*/

/*--------------------------------------------------------------------------*/
/* Histogram */

   if(iget_full_inter==1){
    if(num==1){
     ibin      = (int) ( ((theta_now[1]+M_PI)/tors_free_del) + 1.0 );
     ibin      = MIN((tors_free_nhist),ibin);
     ibin      = MAX(1,ibin);
     tors_free_hist[ibin] += 1.0;
    }/*endif*/
    if(num==2){
     jbin      = (int) ( ((theta_now[1]+M_PI)/tors_free_del) + 1.0 );
     jbin      = MIN((tors_free_nhist),jbin);
     jbin      = MAX(1,jbin);
     ibin      = (int) ( ((theta_now[2]+M_PI)/tors_free_del) + 1.0 );
     ibin      = MIN((tors_free_nhist),ibin);
     ibin      = MAX(1,ibin);
     tors_free_hist_2d[jbin][ibin] += 1.0;
    }/*endif*/
  }/*endif*/

/*--------------------------------------------------------------------------*/
   }/* end routine */
/*==========================================================================*/





/*==========================================================================*/
/*               Header:                                                    */
/*==========================================================================*/

void tors_free_mode(CLATOMS_INFO *clatoms_info,CLATOMS_POS *clatoms_pos,
                    TORS_FREE *tors_free,CELL *cell,
                    INTRA_SCR *intra_scr,PTENS *ptens, double *vtorst,
                    ENERGY_CTRL *energy_ctrl,int np_forc)

/*==========================================================================*/
/*               Begin subprogram:                                          */
{/*begin routine*/
/*==========================================================================*/
/*               Local variable declarations                                */
   int ibin,jbin,nnow,i;  
   double r122,r12,sign_theta;
   double r132,r13,arg,tpi;
   double r232,r23;
   double r422,r42;
   double r432,r43;
                                                 /* Distances              */
   double dxp,dyp,dzp,rp2i,rpi,rpqi;
   double dxq,dyq,dzq,rq2i,rqi;     
   double dp13,dp42,dq13,dq42;
   double cost,sint,seps,pre;          
   double costrp2,costrq2,dp13crp2,dq42crq2;
   double dp13rpq,dq13rpq,dp42rpq,dq42rpq;
   double d1323i,palp,palp1;
   double d4223i,qalp,qalp1; 
   double fpx1,fpy1,fpz1;
   double fpx2,fpy2,fpz2;  
   double fqx4,fqy4,fqz4;
   double fqx2,fqy2,fqz2;  
   double apow,apow1,theta0,theta,dvtors;
                                                /* Force temporaries       */
   double temp1p,temp1q,temp2q,temp2p;
   double temp3p,temp3q,temp4q,temp4p;

   double xcrs_crs,ycrs_crs,zcrs_crs,ssdot;
   double theta_now[4];

/*-------------------------------------------------------------------------*/
/*               Local pointer declarations                                */

   double *xmod            = clatoms_info->xmod;
   double *ymod            = clatoms_info->ymod;
   double *zmod            = clatoms_info->zmod;
   int pi_beads            = clatoms_info->pi_beads; 

   int num                = tors_free->num;
   int *tors_free_j1      = tors_free->j1;
   int *tors_free_j2      = tors_free->j2;
   int *tors_free_j3      = tors_free->j3;
   int *tors_free_j4      = tors_free->j4;
   int tors_free_nhist    = tors_free->nhist;

   double *tors_free_eq  = tors_free->eq;
   double tors_free_npow = tors_free->npow;
   double tors_free_fk   = tors_free->fk;
   double tors_free_del   = tors_free->del;
   double *tors_free_hist = tors_free->hist;
   double **tors_free_hist_2d = tors_free->hist_2d;

   double *clatoms_pos_fx = clatoms_pos->fx;
   double *clatoms_pos_fy = clatoms_pos->fy;
   double *clatoms_pos_fz = clatoms_pos->fz;

   double *intra_scr_x1 =  intra_scr->x1;
   double *intra_scr_y1 =  intra_scr->y1; 
   double *intra_scr_z1 =  intra_scr->z1; 
   double *intra_scr_x2 =  intra_scr->x2; 
   double *intra_scr_y2 =  intra_scr->y2; 
   double *intra_scr_z2 =  intra_scr->z2; 
   double *intra_scr_x3 =  intra_scr->x3; 
   double *intra_scr_y3 =  intra_scr->y3; 
   double *intra_scr_z3 =  intra_scr->z3;  
   double *intra_scr_x4 =  intra_scr->x4;  
   double *intra_scr_y4 =  intra_scr->y4;  
   double *intra_scr_z4 =  intra_scr->z4;  

   double *intra_scr_fx1 =  intra_scr->fx1;
   double *intra_scr_fy1 =  intra_scr->fy1; 
   double *intra_scr_fz1 =  intra_scr->fz1; 
   double *intra_scr_fx2 =  intra_scr->fx2; 
   double *intra_scr_fy2 =  intra_scr->fy2; 
   double *intra_scr_fz2 =  intra_scr->fz2; 
   double *intra_scr_fx3 =  intra_scr->fx3; 
   double *intra_scr_fy3 =  intra_scr->fy3; 
   double *intra_scr_fz3 =  intra_scr->fz3;  
   double *intra_scr_fx4 =  intra_scr->fx4;  
   double *intra_scr_fy4 =  intra_scr->fy4;  
   double *intra_scr_fz4 =  intra_scr->fz4;  

   double *intra_scr_x5=  intra_scr->x5;
   double *intra_scr_y5=  intra_scr->y5;
   double *intra_scr_z5=  intra_scr->z5;
   double *intra_scr_x6=  intra_scr->x6;
   double *intra_scr_y6=  intra_scr->y6;
   double *intra_scr_z6=  intra_scr->z6;
   double *intra_scr_x7=  intra_scr->x7;
   double *intra_scr_y7=  intra_scr->y7;
   double *intra_scr_z7=  intra_scr->z7;
   double *intra_scr_x8=  intra_scr->x8;
   double *intra_scr_y8=  intra_scr->y8; 
   double *intra_scr_z8=  intra_scr->z8;

   double *intra_scr_dx12 = intra_scr->dx12;
   double *intra_scr_dy12 = intra_scr->dy12;
   double *intra_scr_dz12 = intra_scr->dz12;
   double *intra_scr_dx13 = intra_scr->dx13;
   double *intra_scr_dy13 = intra_scr->dy13;
   double *intra_scr_dz13 = intra_scr->dz13;
   double *intra_scr_dx42 = intra_scr->dx42;
   double *intra_scr_dy42 = intra_scr->dy42;
   double *intra_scr_dz42 = intra_scr->dz42;
   double *intra_scr_dx43 = intra_scr->dx43;
   double *intra_scr_dy43 = intra_scr->dy43;
   double *intra_scr_dz43 = intra_scr->dz43;
   double *intra_scr_dx23 = intra_scr->dx23;
   double *intra_scr_dy23 = intra_scr->dy23;
   double *intra_scr_dz23 = intra_scr->dz23;

   double *intra_scr_dx56 = intra_scr->dx56;
   double *intra_scr_dy56 = intra_scr->dy56;
   double *intra_scr_dz56 = intra_scr->dz56;
   double *intra_scr_dx57 = intra_scr->dx57;
   double *intra_scr_dy57 = intra_scr->dy57;
   double *intra_scr_dz57 = intra_scr->dz57;
   double *intra_scr_dx86 = intra_scr->dx86;
   double *intra_scr_dy86 = intra_scr->dy86;
   double *intra_scr_dz86 = intra_scr->dz86;
   double *intra_scr_dx87 = intra_scr->dx87;
   double *intra_scr_dy87 = intra_scr->dy87;
   double *intra_scr_dz87 = intra_scr->dz87;
   double *intra_scr_dx67 = intra_scr->dx67;
   double *intra_scr_dy67 = intra_scr->dy67;
   double *intra_scr_dz67 = intra_scr->dz67;

   double *intra_scr_p11 = intra_scr->p11;
   double *intra_scr_p22 = intra_scr->p22;
   double *intra_scr_p33 = intra_scr->p33;
   double *intra_scr_p12 = intra_scr->p12;
   double *intra_scr_p21 = intra_scr->p21;
   double *intra_scr_p31 = intra_scr->p31;
   double *intra_scr_p13 = intra_scr->p13;
   double *intra_scr_p32 = intra_scr->p32;
   double *intra_scr_p23 = intra_scr->p23;

   double wfor           = intra_scr->wght_tra;

   double *ptens_pvten     = ptens->pvten;
   double *ptens_pvten_tot = ptens->pvten_tot;

   int iperd              = cell->iperd;
   int intra_perds        = cell->intra_perds;

   int iget_pv_real_inter = energy_ctrl->iget_pv_real_inter;
   int iget_full_inter    = energy_ctrl->iget_full_inter;
 
/*==========================================================================*/
/* 0) Lower cutoff for sin(theta)                                           */

   seps = 1.0e-8;
   tpi  = 2.0*M_PI;

/* I) loop over all the tors in steps of nlen to save memory                */

   (*vtorst) = 0.0;
/*--------------------------------------------------------------------------*/
/*  B) Gather positions                                                     */

  for(i=1;i<=num;i++){
   intra_scr_x1[i] = xmod[tors_free_j1[i]];
   intra_scr_y1[i] = ymod[tors_free_j1[i]];
   intra_scr_z1[i] = zmod[tors_free_j1[i]];
   intra_scr_x2[i] = xmod[tors_free_j2[i]];
   intra_scr_y2[i] = ymod[tors_free_j2[i]];
   intra_scr_z2[i] = zmod[tors_free_j2[i]];
   intra_scr_x3[i] = xmod[tors_free_j3[i]];
   intra_scr_y3[i] = ymod[tors_free_j3[i]];
   intra_scr_z3[i] = zmod[tors_free_j3[i]];
   intra_scr_x4[i] = xmod[tors_free_j4[i]];
   intra_scr_y4[i] = ymod[tors_free_j4[i]];
   intra_scr_z4[i] = zmod[tors_free_j4[i]];
  }/*endfor*/

/*--------------------------------------------------------------------------*/
/*  D) Calculate the basis vectors (r12,r13,r43,r42)                       */

   for(i=1;i<=num;i++){
    intra_scr_dx12[i] = intra_scr_x1[i] - intra_scr_x2[i];
    intra_scr_dy12[i] = intra_scr_y1[i] - intra_scr_y2[i];
    intra_scr_dz12[i] = intra_scr_z1[i] - intra_scr_z2[i];
    intra_scr_dx13[i] = intra_scr_x1[i] - intra_scr_x3[i];
    intra_scr_dy13[i] = intra_scr_y1[i] - intra_scr_y3[i];
    intra_scr_dz13[i] = intra_scr_z1[i] - intra_scr_z3[i];
    intra_scr_dx42[i] = intra_scr_x4[i] - intra_scr_x2[i];
    intra_scr_dy42[i] = intra_scr_y4[i] - intra_scr_y2[i];
    intra_scr_dz42[i] = intra_scr_z4[i] - intra_scr_z2[i];
    intra_scr_dx43[i] = intra_scr_x4[i] - intra_scr_x3[i];
    intra_scr_dy43[i] = intra_scr_y4[i] - intra_scr_y3[i];
    intra_scr_dz43[i] = intra_scr_z4[i] - intra_scr_z3[i];
    intra_scr_dx23[i] = intra_scr_x2[i] - intra_scr_x3[i];
    intra_scr_dy23[i] = intra_scr_y2[i] - intra_scr_y3[i];
    intra_scr_dz23[i] = intra_scr_z2[i] - intra_scr_z3[i];
   }/*endfor*/
 
/*--------------------------------------------------------------------------*/
/*  E) Periodic boundary conditions                                         */

   nnow = num;
   if(intra_perds==1){
     period(nnow,intra_scr_dx12,intra_scr_dy12,intra_scr_dz12,cell);
     period(nnow,intra_scr_dx13,intra_scr_dy13,intra_scr_dz13,cell);
     period(nnow,intra_scr_dx42,intra_scr_dy42,intra_scr_dz42,cell);
     period(nnow,intra_scr_dx43,intra_scr_dy43,intra_scr_dz43,cell);
     period(nnow,intra_scr_dx23,intra_scr_dy23,intra_scr_dz23,cell);
   }/*endif*/

  for(i=1;i<=num;i++){
   r122   = (intra_scr_dx12[i]*intra_scr_dx12[i]
           +  intra_scr_dy12[i]*intra_scr_dy12[i]
           +  intra_scr_dz12[i]*intra_scr_dz12[i]);
   r12    = sqrt(r122);
   r132   = (intra_scr_dx13[i]*intra_scr_dx13[i]
           +  intra_scr_dy13[i]*intra_scr_dy13[i]
           +  intra_scr_dz13[i]*intra_scr_dz13[i]);
   r13    = sqrt(r132); 
   r422   = (intra_scr_dx42[i]*intra_scr_dx42[i]
           +  intra_scr_dy42[i]*intra_scr_dy42[i]
           +  intra_scr_dz42[i]*intra_scr_dz42[i]);
   r42    = sqrt(r422);
   r432   = (intra_scr_dx43[i]*intra_scr_dx43[i]
           +  intra_scr_dy43[i]*intra_scr_dy43[i]
           +  intra_scr_dz43[i]*intra_scr_dz43[i]);
   r43    = sqrt(r432);
   r232   = (intra_scr_dx23[i]*intra_scr_dx23[i]
           +  intra_scr_dy23[i]*(intra_scr_dy23)[i]
           +  intra_scr_dz23[i]*intra_scr_dz23[i]);
   r23    = sqrt(r232);
 
/*------------------------------------------------------------------*/
/*F) Construct vector in plane of 123 perp to r23 call it(dxp,dyp,dzp) */

   d1323i  = 1.0/(intra_scr_dx13[i]*intra_scr_dx23[i]
                + intra_scr_dy13[i]*intra_scr_dy23[i]
                + intra_scr_dz13[i]*intra_scr_dz23[i]);
   palp   =-(intra_scr_dx12[i]*intra_scr_dx23[i]
          +  intra_scr_dy12[i]*intra_scr_dy23[i]
          +  intra_scr_dz12[i]*intra_scr_dz23[i])*d1323i;
   palp1  = 1.0+palp;
   dxp    = palp*intra_scr_dx13[i]+intra_scr_dx12[i];
   dyp    = palp*intra_scr_dy13[i]+intra_scr_dy12[i];
   dzp    = palp*intra_scr_dz13[i]+intra_scr_dz12[i];
   rp2i   = 1.0/(dxp*dxp+dyp*dyp+dzp*dzp);
   rpi    = sqrt(rp2i);

/*-----------------------------------------------------------------*/
/*  G) Get dot product of this vector (dxp,dyp,dzp) with r13 and r42 */
      
   dp13  = (dxp*intra_scr_dx13[i]+dyp*intra_scr_dy13[i]
           +dzp*intra_scr_dz13[i]);
   dp42  = (dxp*intra_scr_dx42[i]+dyp*intra_scr_dy42[i]
           +dzp*intra_scr_dz42[i]);
      
/*-----------------------------------------------------------------*/
/*  H) Get derivative of parameter palp used to construct (dxp,dyp,dzp)*/
/*     with respect to atoms (1-4)                                     */
      
   fpx1   = -palp1*d1323i*intra_scr_dx23[i];
   fpy1   = -palp1*d1323i*intra_scr_dy23[i];
   fpz1   = -palp1*d1323i*intra_scr_dz23[i];
   fpx2   = (intra_scr_dx23[i]-dxp)*d1323i;
   fpy2   = (intra_scr_dy23[i]-dyp)*d1323i;
   fpz2   = (intra_scr_dz23[i]-dzp)*d1323i;

/*------------------------------------------------------------------*/
/*  I) Construct vector in plane of 234 perp to r23 (dxq,dyq,dzq) */
      
   d4223i  = 1.0/(intra_scr_dx42[i]*intra_scr_dx23[i]
                + intra_scr_dy42[i]*intra_scr_dy23[i]
                + intra_scr_dz42[i]*intra_scr_dz23[i]);
   qalp   =-(intra_scr_dx43[i]*intra_scr_dx23[i]
          +  intra_scr_dy43[i]*intra_scr_dy23[i]
          +  intra_scr_dz43[i]*intra_scr_dz23[i])*d4223i;
   qalp1  = 1.0+qalp;
   dxq    = qalp*intra_scr_dx42[i]+intra_scr_dx43[i];
   dyq    = qalp*intra_scr_dy42[i]+intra_scr_dy43[i];
   dzq    = qalp*intra_scr_dz42[i]+intra_scr_dz43[i];
   rq2i   = 1.0/(dxq*dxq+dyq*dyq+dzq*dzq);
   rqi    = sqrt(rq2i);
/*------------------------------------------------------------------*/
/*  J) Get dot product of this vector (dxq,dyq,dzq) with r13 and r42*/
      
   dq13  = (dxq*intra_scr_dx13[i]+dyq*intra_scr_dy13[i]
           +dzq*intra_scr_dz13[i]);
   dq42  = (dxq*intra_scr_dx42[i]+dyq*intra_scr_dy42[i]
           +dzq*intra_scr_dz42[i]);
/*-------------------------------------------------------------------*/
/*  K) Get derivative of parameter qalp used to construct (dxq,dyq,dzq)*/
/*     with respect to atoms (1-4)                                   */
      
   fqx4   = -qalp1*d4223i*intra_scr_dx23[i];
   fqy4   = -qalp1*d4223i*intra_scr_dy23[i];
   fqz4   = -qalp1*d4223i*intra_scr_dz23[i];
   fqx2   = (qalp*intra_scr_dx23[i]-dxq)*d4223i;
   fqy2   = (qalp*intra_scr_dy23[i]-dyq)*d4223i;
   fqz2   = (qalp*intra_scr_dz23[i]-dzq)*d4223i;

/*--------------------------------------------------------------------*/
/*  L) Get cosine  and sine of the angle                              */

   rpqi   = rpi*rqi;

   xcrs_crs = (dyq*dzp - dyp*dzq);
   ycrs_crs = (dzq*dxp - dzp*dxq);
   zcrs_crs = (dxq*dyp - dxp*dyq);
   sint     = (xcrs_crs*intra_scr_dx23[i]
              +ycrs_crs*intra_scr_dy23[i]
              +zcrs_crs*intra_scr_dz23[i])*rpqi/r23;
   sint   = (sint < 1.0 ? sint:1.0);
   sint   = (sint > -1.0 ? sint:-1.0);

   cost   = (dxp*dxq+dyp*dyq+dzp*dzq)*rpqi;
   cost   = (cost < 1.0 ? cost:1.0);
   cost   = (cost > -1.0 ? cost:-1.0);

   theta      = acos(cost);
   sign_theta = (sint >= 0 ? 1.0: -1.0);
   theta     *= sign_theta;
   sint       = (sign_theta*sint > seps ? sint:sign_theta*seps);

   theta_now[i] = theta;
   theta0  = theta - (tors_free_eq[i]);
   arg     = theta0/tpi;
   theta0 -= (tpi*NINT(arg));

/*--------------------------------------------------------------------------*/
/*  M) Calculate the potential energy and bin the stuff                     */

   apow      = (double)(tors_free_npow);
   apow1     = apow-1.0;

   (*vtorst) += (tors_free_fk)*pow(theta0,apow)/apow;

/*--------------------------------------------------------------------------*/
/*  N) Get the force on the atoms using the chain rule and Horner's method  */

   dvtors = (tors_free_fk)*pow(theta0,apow1);
   pre    = dvtors/sint;

   costrp2 = cost*rp2i;
   costrq2 = cost*rq2i;
   dp13crp2 = dp13*costrp2;
   dq42crq2 = dq42*costrq2;
   dp13rpq = dp13*rpqi;
   dq13rpq = dq13*rpqi;
   dp42rpq = dp42*rpqi;
   dq42rpq = dq42*rpqi;
   temp1p  = (dq13rpq-dp13crp2)*pre;
   temp1q  = (dp42rpq-dq42crq2)*pre;
   temp2p  = (costrp2*palp1)*pre;
   temp2q  = (costrq2*qalp1)*pre;
   temp3q  = (rpqi*palp1)*pre;  
   temp3p  = (rpqi*qalp1)*pre;  
   temp4p  = (costrp2-qalp*rpqi)*pre;
   temp4q  = (qalp*costrq2-rpqi)*pre;

   intra_scr_fx1[i] = -dxp*temp2p+dxq*temp3q+temp1p*fpx1;
   intra_scr_fy1[i] = -dyp*temp2p+dyq*temp3q+temp1p*fpy1;
   intra_scr_fz1[i] = -dzp*temp2p+dzq*temp3q+temp1p*fpz1;

   intra_scr_fx4[i] = -dxq*temp2q+dxp*temp3p+temp1q*fqx4;
   intra_scr_fy4[i] = -dyq*temp2q+dyp*temp3p+temp1q*fqy4;
   intra_scr_fz4[i] = -dzq*temp2q+dzp*temp3p+temp1q*fqz4;

   intra_scr_fx2[i] =  dxp*temp4p+dxq*temp4q+temp1p*fpx2+temp1q*fqx2;
   intra_scr_fy2[i] =  dyp*temp4p+dyq*temp4q+temp1p*fpy2+temp1q*fqy2;
   intra_scr_fz2[i] =  dzp*temp4p+dzq*temp4q+temp1p*fpz2+temp1q*fqz2;

   intra_scr_fx3[i] = -(intra_scr_fx1[i]+intra_scr_fx2[i]
                      +intra_scr_fx4[i]);
   intra_scr_fy3[i] = -(intra_scr_fy1[i]+intra_scr_fy2[i]
                      +intra_scr_fy4[i]);
   intra_scr_fz3[i] = -(intra_scr_fz1[i]+intra_scr_fz2[i]
                      +intra_scr_fz4[i]);

/*--------------------------------------------------------------------------*/
/*  O) Get the pressure tensor (this form is ok for periodic systems):      */

   intra_scr_p11[i] = intra_scr_dx13[i]*intra_scr_fx1[i]
                       + intra_scr_dx23[i]*intra_scr_fx2[i] 
                       + intra_scr_dx43[i]*intra_scr_fx4[i];
   intra_scr_p22[i] = intra_scr_dy13[i]*intra_scr_fy1[i]
                       + intra_scr_dy23[i]*intra_scr_fy2[i] 
                       + intra_scr_dy43[i]*intra_scr_fy4[i];
   intra_scr_p33[i] = intra_scr_dz13[i]*intra_scr_fz1[i]
                       + intra_scr_dz23[i]*intra_scr_fz2[i] 
                       + intra_scr_dz43[i]*intra_scr_fz4[i];
   intra_scr_p12[i] = intra_scr_dx13[i]*intra_scr_fy1[i]
                       + intra_scr_dx23[i]*intra_scr_fy2[i] 
                       + intra_scr_dx43[i]*intra_scr_fy4[i];
   intra_scr_p21[i] = intra_scr_dy13[i]*intra_scr_fx1[i]
                       + intra_scr_dy23[i]*intra_scr_fx2[i] 
                       + intra_scr_dy43[i]*intra_scr_fx4[i];
   intra_scr_p13[i] = intra_scr_dx13[i]*intra_scr_fz1[i]
                       + intra_scr_dx23[i]*intra_scr_fz2[i] 
                       + intra_scr_dx43[i]*intra_scr_fz4[i];
   intra_scr_p31[i] = intra_scr_dz13[i]*intra_scr_fx1[i]
                       + intra_scr_dz23[i]*intra_scr_fx2[i]
                       + intra_scr_dz43[i]*intra_scr_fx4[i];
   intra_scr_p23[i] = intra_scr_dy13[i]*intra_scr_fz1[i]
                       + intra_scr_dy23[i]*intra_scr_fz2[i] 
                       + intra_scr_dy43[i]*intra_scr_fz4[i];
   intra_scr_p32[i] = intra_scr_dz13[i]*intra_scr_fy1[i]
                       + intra_scr_dz23[i]*intra_scr_fy2[i] 
                       + intra_scr_dz43[i]*intra_scr_fy4[i]; 
 }/*endfor*/
 
/*--------------------------------------------------------------------------*/
/*  Q) Sum the pressure tensor                                              */

   if(iperd == 2 || iperd == 3) {
    for(i=1;i<=num;i++){
     ptens_pvten[1] += intra_scr_p11[i]*wfor*pi_beads;
     ptens_pvten[5] += intra_scr_p22[i]*wfor*pi_beads;
     ptens_pvten[9] += intra_scr_p33[i]*wfor*pi_beads;
     ptens_pvten[2] += intra_scr_p12[i]*wfor*pi_beads;
     ptens_pvten[4] += intra_scr_p21[i]*wfor*pi_beads;
    }/*endfor*/
    if(iget_pv_real_inter==1){    
     for(i=1;i<=num;i++){
       ptens_pvten_tot[1] += intra_scr_p11[i]*pi_beads;
       ptens_pvten_tot[5] += intra_scr_p22[i]*pi_beads;
       ptens_pvten_tot[9] += intra_scr_p33[i]*pi_beads;
       ptens_pvten_tot[2] += intra_scr_p12[i]*pi_beads;
       ptens_pvten_tot[4] += intra_scr_p21[i]*pi_beads;
     }/*endfor*/
    }/*endif*/

   }/*endif*/
 
   if(iperd == 3) {
    for(i=1;i<=num;i++){
     ptens_pvten[3] += intra_scr_p13[i]*wfor*pi_beads;
     ptens_pvten[7] += intra_scr_p31[i]*wfor*pi_beads;
     ptens_pvten[6] += intra_scr_p23[i]*wfor*pi_beads;
     ptens_pvten[8] += intra_scr_p32[i]*wfor*pi_beads;
    }/*endfor*/

    if(iget_pv_real_inter==1){    
     for(i=1;i<=num;i++){
      ptens_pvten_tot[3] += intra_scr_p13[i]*pi_beads;
      ptens_pvten_tot[7] += intra_scr_p31[i]*pi_beads;
      ptens_pvten_tot[6] += intra_scr_p23[i]*pi_beads;
      ptens_pvten_tot[8] += intra_scr_p32[i]*pi_beads;
     }/*endfor*/
    }/*endif*/

   }/*endif*/
 
/*--------------------------------------------------------------------------*/
/*  R) Scatter the forces                                                   */

  for(i=1;i<=num;i++){
   clatoms_pos_fx[tors_free_j1[i]]  += intra_scr_fx1[i]*wfor;
   clatoms_pos_fy[tors_free_j1[i]]  += intra_scr_fy1[i]*wfor;
   clatoms_pos_fz[tors_free_j1[i]]  += intra_scr_fz1[i]*wfor;

   clatoms_pos_fx[tors_free_j2[i]]  += intra_scr_fx2[i]*wfor;
   clatoms_pos_fy[tors_free_j2[i]]  += intra_scr_fy2[i]*wfor;
   clatoms_pos_fz[tors_free_j2[i]]  += intra_scr_fz2[i]*wfor;

   clatoms_pos_fx[tors_free_j3[i]]  += intra_scr_fx3[i]*wfor;
   clatoms_pos_fy[tors_free_j3[i]]  += intra_scr_fy3[i]*wfor;
   clatoms_pos_fz[tors_free_j3[i]]  += intra_scr_fz3[i]*wfor;

   clatoms_pos_fx[tors_free_j4[i]]  += intra_scr_fx4[i]*wfor;
   clatoms_pos_fy[tors_free_j4[i]]  += intra_scr_fy4[i]*wfor;
   clatoms_pos_fz[tors_free_j4[i]]  += intra_scr_fz4[i]*wfor;
 }/*endfor*/

/*--------------------------------------------------------------------------*/
/* Histogram */

   if(iget_full_inter==1){
    if(num==1){
     ibin      = (int) ( ((theta_now[1]+M_PI)/tors_free_del) + 1.0 );
     ibin      = MIN((tors_free_nhist),ibin);
     ibin      = MAX(1,ibin);
     tors_free_hist[ibin] += 1.0;
    }/*endif*/
    if(num==2){
     jbin      = (int) ( ((theta_now[1]+M_PI)/tors_free_del) + 1.0 );
     jbin      = MIN((tors_free_nhist),jbin);
     jbin      = MAX(1,jbin);
     ibin      = (int) ( ((theta_now[2]+M_PI)/tors_free_del) + 1.0 );
     ibin      = MIN((tors_free_nhist),ibin);
     ibin      = MAX(1,ibin);
     tors_free_hist_2d[jbin][ibin] += 1.0;
    }/*endif*/
  }/*endif*/

/*--------------------------------------------------------------------------*/
  }/* end routine */
/*==========================================================================*/







