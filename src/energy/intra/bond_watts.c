/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
/*                                                                          */
/*                         PI_MD:                                           */
/*             The future of simulation technology                          */
/*             ------------------------------------                         */
/*                     Module: bond_watts.c                                 */
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

void bond_watts_33(CLATOMS_INFO *clatoms_info,CLATOMS_POS *clatoms_pos,
          GRP_BOND_WATTS *grp_bond_watts,CELL *cell,
          INTRA_SCR *intra_scr,PTENS *ptens, double *vbendt, 
          double *vbondt, double *vtot, int iver_get,
          CLASS_COMM_FORC_PKG *class_comm_forc_pkg,
          int iget_pv_real_inter)

/*==========================================================================*/
{/*begin routine*/
  /*=======================================================================*/
  /*            Local variable declarations                                */

  int iii;
  int ibend_bnd,ibig,ibig1,ioff,iend,nnow; /* Indices and counters       */
  int nlen_use,nlen_now;
  int ntot;
  int ngo,irem,istart_big;

  double r122,r12;                         /* (r-r_0)^2 and (r-r_0)      */
  double r132,r13;
  double theta;
  double cost,sint,seps;
  double vbendc,vbends;            /* Bend_bnd potential         */
  double dvbendc,dvbends;          /* Derivative of bend_bnd pot */
  double dvbond;
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
  double dr12_dx12,dr13_dx13;
  double dr12_dy12,dr13_dy13;
  double dr12_dz12,dr13_dz13;
  double cost_2,sint_2;
  double sin_use,cos_use;
  double dcost_dx12,dcost_dx13;
  double dcost_dy12,dcost_dy13;
  double dcost_dz12,dcost_dz13,r_red,ds1,ds2,ds3;
  double vbond1,vbond2,vbend;
  double dvbond1,dvbond2;
  double dvbond1_dcost,dvbond2_dcost;
  double dvbond1_dr12,dvbond1_dr13,dvbond2_dr13;
  double dvbend,dvbend_dcost;
  double dvbend_dr12,dvbend_dr13;
  int ktemp1,ktemp2,ktemp3,ktemp5,ktemp6,ktemp7,iwatt;


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
  double *intra_scr_dx13   = intra_scr->dx23;
  double *intra_scr_dy13   = intra_scr->dy23;
  double *intra_scr_dz13   = intra_scr->dz23;
  double *intra_scr_dx56   = intra_scr->dx56;
  double *intra_scr_dy56   = intra_scr->dy56;
  double *intra_scr_dz56   = intra_scr->dz56;
  double *intra_scr_dx57   = intra_scr->dx67;
  double *intra_scr_dy57   = intra_scr->dy67;
  double *intra_scr_dz57   = intra_scr->dz67;
  double *intra_scr_vpot   = intra_scr->vpot;
  double *intra_scr_p11    = intra_scr->p11;
  double *intra_scr_p22    = intra_scr->p22;
  double *intra_scr_p33    = intra_scr->p33;
  double *intra_scr_p12    = intra_scr->p12;
  double *intra_scr_p13    = intra_scr->p13;
  double *intra_scr_p23    = intra_scr->p23;
  double *intra_scr_vbond  = intra_scr->vpot;
  double *intra_scr_vbend  = intra_scr->dx87;
  double **eq  = grp_bond_watts->eq_33;
  double **c_0 = grp_bond_watts->c_0_33;
  double **c_1 = grp_bond_watts->c_1_33;
  double **c_2 = grp_bond_watts->c_2_33;
  double **c_3 = grp_bond_watts->c_3_33;
  double **c_4 = grp_bond_watts->c_4_33;
  double **c_5 = grp_bond_watts->c_5_33;
  double **c_6 = grp_bond_watts->c_6_33;
  double **dc_0= grp_bond_watts->dc_0_33;
  double **dc_1= grp_bond_watts->dc_1_33;
  double **dc_2= grp_bond_watts->dc_2_33;
  double **dc_3= grp_bond_watts->dc_3_33;
  double **dc_4= grp_bond_watts->dc_4_33;
  double **dc_5= grp_bond_watts->dc_5_33;
  double **dc_6= grp_bond_watts->dc_6_33;
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
  double *cos_thet0_2     = grp_bond_watts->cos_thet0_2;
  double *sin_thet0_2     = grp_bond_watts->sin_thet0_2;
  int *j1          = grp_bond_watts->j1_33;
  int *j2          = grp_bond_watts->j2_33;
  int *j3          = grp_bond_watts->j3_33;
  int *jtyp        = grp_bond_watts->jtyp_33;
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
  ntot = (grp_bond_watts->num_33);

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
    
    for(iboff=ibig+ioff,iwatt=ibig;iwatt <= iend; 
                                                ++iwatt,++iboff){
      ktemp1 = j1[iwatt];
      ktemp2 = j2[iwatt];
      ktemp3 = j3[iwatt];
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
      intra_scr_dx13[iboff] = x1 - x3;
      intra_scr_dy13[iboff] = y1 - y3;
      intra_scr_dz13[iboff] = z1 - z3;
    }/*endfor*/

   if(cell->intra_perds==1&&clatoms_info->pi_beads>1){     
    for(iboff=ibig+ioff,iwatt=ibig;iwatt <= iend; ++iwatt,++iboff){
        ktemp5 = j1[iwatt];
        ktemp6 = j2[iwatt];
        ktemp7 = j3[iwatt];
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
        intra_scr_dx57[iboff] = x5 - x7;
        intra_scr_dy57[iboff] = y5 - y7;
        intra_scr_dz57[iboff] = z5 - z7;
      }/*endfor*/
    }/*endif*/     

/*=======================================================================*/
/* IV) Periodic boundary conditions                                      */
    
   if(cell->intra_perds == 1) {
     if(clatoms_info->pi_beads==1){     
       period(nnow,intra_scr_dx12,intra_scr_dy12,intra_scr_dz12,cell);
       period(nnow,intra_scr_dx13,intra_scr_dy13,intra_scr_dz13,cell);
     }else{
       period_pimd(nnow,intra_scr_dx12,intra_scr_dy12,intra_scr_dz12,
                         intra_scr_dx56,intra_scr_dy56,intra_scr_dz56,cell);
       period_pimd(nnow,intra_scr_dx13,intra_scr_dy13,intra_scr_dz13,
                         intra_scr_dx57,intra_scr_dy57,intra_scr_dz57,cell);
     }/*endif*/
   }/*endif*/

/*=======================================================================*/
/* VI) Get the Watts energies                                            */

    for(iwatt=1,iboff=ibig;iwatt <= nnow; ++iwatt,++iboff) {

      ktemp = jtyp[iboff];

      /*-------------------------------------------------------------------*/
      /*  A) Get the bond lengths. */
      

      r122  = (intra_scr_dx12[iwatt]*intra_scr_dx12[iwatt]
               + intra_scr_dy12[iwatt]*intra_scr_dy12[iwatt]
               + intra_scr_dz12[iwatt]*intra_scr_dz12[iwatt]);
      r12    = sqrt(r122);
      
      r132  = (intra_scr_dx13[iwatt]*intra_scr_dx13[iwatt]
               + intra_scr_dy13[iwatt]*intra_scr_dy13[iwatt]
               + intra_scr_dz13[iwatt]*intra_scr_dz13[iwatt]);
      r13    = sqrt(r132);
      
      dr12_dx12 = intra_scr_dx12[iwatt]/r12;
      dr13_dx13 = intra_scr_dx13[iwatt]/r13;
      dr12_dy12 = intra_scr_dy12[iwatt]/r12;
      dr13_dy13 = intra_scr_dy13[iwatt]/r13;
      dr12_dz12 = intra_scr_dz12[iwatt]/r12;
      dr13_dz13 = intra_scr_dz13[iwatt]/r13;
      
      /*-------------------------------------------------------------------*/
      /*  B) Calculate the cosine and sine of the angle*/
      
      rpmag  = 1.0/(r12*r13);
      cost = (intra_scr_dx12[iwatt]*intra_scr_dx13[iwatt]
              +  intra_scr_dy12[iwatt]*intra_scr_dy13[iwatt]
              +  intra_scr_dz12[iwatt]*intra_scr_dz13[iwatt])
              *rpmag;
      
      cost   = (cost < 1.0 ? cost:1.0);
      cost   = (cost > -1.0 ? cost:-1.0);
      cost_2 = sqrt(0.5*(1.0+cost));
      
      sint   = sqrt(1.0 - cost*cost);
      sint   = (sint > seps ? sint:seps);
      sint_2 = sqrt(0.5*(1.0-cost));

      cos_use = cos_thet0_2[ktemp]*cost_2 + sint_2*sin_thet0_2[ktemp];
      sin_use = cos_thet0_2[ktemp]*sint_2 - cost_2*sin_thet0_2[ktemp];

      dcost_dx12 = (rpmag*intra_scr_dx13[iwatt]
                          -cost*intra_scr_dx12[iwatt]/(r122));
      dcost_dx13 = (rpmag*intra_scr_dx12[iwatt]
                          -cost*intra_scr_dx13[iwatt]/(r132));
      dcost_dy12 = (rpmag*intra_scr_dy13[iwatt]
                          -cost*intra_scr_dy12[iwatt]/(r122));
      dcost_dy13 = (rpmag*intra_scr_dy12[iwatt]
                          -cost*intra_scr_dy13[iwatt]/(r132));
      dcost_dz12 = (rpmag*intra_scr_dz13[iwatt]
                          -cost*intra_scr_dz12[iwatt]/(r122));
      dcost_dz13 = (rpmag*intra_scr_dz12[iwatt]
                          -cost*intra_scr_dz13[iwatt]/(r132));

       theta = acos(cost)*180.0/M_PI;
      /*-------------------------------------------------------------------*/
      /*  C) Construct the local modes */
      
      ds1 = r12*cos_use - eq[1][ktemp];
      ds2 = r13*cos_use - eq[2][ktemp];
      r_red = r12/eq[1][ktemp] + r13/eq[2][ktemp];
      ds3 = r_red*sin_use;

      /*--------------------------------------------------------------------*/
      /*  D) Get the Watts potential energy */
      
      vbond1 =  ((((( c_6[1][ktemp]
                    *ds1 + c_5[1][ktemp])
                   *ds1 + c_4[1][ktemp])
                  *ds1 + c_3[1][ktemp])
                 *ds1 + c_2[1][ktemp])
                *ds1 + c_1[1][ktemp])
        *ds1 + c_0[1][ktemp];
      
      vbond2 =  ((((( c_6[2][ktemp]
                    *ds2 + c_5[2][ktemp])
                   *ds2 + c_4[2][ktemp])
                  *ds2 + c_3[2][ktemp])
                 *ds2 + c_2[2][ktemp])
                *ds2 + c_1[2][ktemp])
        *ds2 + c_0[2][ktemp];
      
      intra_scr_vbond[iwatt] = vbond1+vbond2;
      
      vbend =  ((((( c_6[3][ktemp]
                    *ds3 + c_5[3][ktemp])
                   *ds3 + c_4[3][ktemp])
                  *ds3 + c_3[3][ktemp])
                 *ds3 + c_2[3][ktemp])
                *ds3 + c_1[3][ktemp])
        *ds3 + c_0[3][ktemp];

      
      intra_scr_vbend[iwatt] = vbend;



      /*-------------------------------------------------------------------*/
      /*  C) Get the force on the atoms using the chain rule and Horner's  */
      
      dvbond1 =  ((((( dc_6[1][ktemp]
                    *ds1 + dc_5[1][ktemp])
                   *ds1 + dc_4[1][ktemp])
                  *ds1 + dc_3[1][ktemp])
                 *ds1 + dc_2[1][ktemp])
                *ds1 + dc_1[1][ktemp]);
      
      dvbond1_dr12 = cos_use*dvbond1;
      dvbond1_dcost = r12*(0.5*sin_use/sint)*dvbond1;

      dvbond2 =  ((((( dc_6[2][ktemp]
                    *ds2 + dc_5[2][ktemp])
                   *ds2 + dc_4[2][ktemp])
                  *ds2 + dc_3[2][ktemp])
                 *ds2 + dc_2[2][ktemp])
                *ds2 + dc_1[2][ktemp]);
  
      dvbond2_dr13 = cos_use*dvbond2;      
      dvbond2_dcost = r13*(0.5*sin_use/sint)*dvbond2;

      
      dvbend =  ((((( dc_6[3][ktemp]
                    *ds3 + dc_5[3][ktemp])
                   *ds3 + dc_4[3][ktemp])
                  *ds3 + dc_3[3][ktemp])
                 *ds3 + dc_2[3][ktemp])
                *ds3 + dc_1[3][ktemp]);

      dvbend_dr12 = sin_use*dvbend/eq[1][ktemp];
      dvbend_dr13 = sin_use*dvbend/eq[2][ktemp];
      dvbend_dcost = -r_red*(0.5*cos_use/sint)*dvbend;

      
      intra_scr_fx2[iwatt]  = (dr12_dx12*(dvbond1_dr12+dvbend_dr12)+
                dcost_dx12*(dvbond1_dcost+dvbond2_dcost+dvbend_dcost));
      intra_scr_fy2[iwatt]  = (dr12_dy12*(dvbond1_dr12+dvbend_dr12)+
                dcost_dy12*(dvbond1_dcost+dvbond2_dcost+dvbend_dcost));
      intra_scr_fz2[iwatt]  = (dr12_dz12*(dvbond1_dr12+dvbend_dr12)+
                dcost_dz12*(dvbond1_dcost+dvbond2_dcost+dvbend_dcost));
      intra_scr_fx3[iwatt]  = (dr13_dx13*(dvbond2_dr13+dvbend_dr13)+
                 dcost_dx13*(dvbond1_dcost+dvbond2_dcost+dvbend_dcost));
      intra_scr_fy3[iwatt]  = (dr13_dy13*(dvbond2_dr13+dvbend_dr13)+
                 dcost_dy13*(dvbond1_dcost+dvbond2_dcost+dvbend_dcost));
      intra_scr_fz3[iwatt]  = (dr13_dz13*(dvbond2_dr13+dvbend_dr13)+
                 dcost_dz13*(dvbond1_dcost+dvbond2_dcost+dvbend_dcost));

      intra_scr_fx1[iwatt]  = -(intra_scr_fx2[iwatt]
                              + intra_scr_fx3[iwatt]);
      intra_scr_fy1[iwatt]  = -(intra_scr_fy2[iwatt]
                              + intra_scr_fy3[iwatt]);
      intra_scr_fz1[iwatt]  = -(intra_scr_fz2[iwatt]
                              + intra_scr_fz3[iwatt]);


      intra_scr_p11[iwatt]=-(intra_scr_dx12[iwatt]
                               *intra_scr_fx2[iwatt]
                               +intra_scr_dx13[iwatt]
                               *intra_scr_fx3[iwatt]);
      intra_scr_p22[iwatt]=-(intra_scr_dy12[iwatt]
                               *intra_scr_fy2[iwatt]
                               +intra_scr_dy13[iwatt]
                               *intra_scr_fy3[iwatt]);
      intra_scr_p33[iwatt]=-(intra_scr_dz12[iwatt]
                               *intra_scr_fz2[iwatt]
                               +intra_scr_dz13[iwatt]
                               *intra_scr_fz3[iwatt]);
      intra_scr_p12[iwatt]=-(intra_scr_dx12[iwatt]
                               *intra_scr_fy2[iwatt]
                               +intra_scr_dx13[iwatt]
                               *intra_scr_fy3[iwatt]);
      intra_scr_p13[iwatt]=-(intra_scr_dx12[iwatt]
                               *intra_scr_fz2[iwatt]
                               +intra_scr_dx13[iwatt]
                               *intra_scr_fz3[iwatt]);
      intra_scr_p23[iwatt]=-(intra_scr_dy12[iwatt]
                               *intra_scr_fz2[iwatt]
                               +intra_scr_dy13[iwatt]
                               *intra_scr_fz3[iwatt]);
    }/*endfor iwatt */

    for(iwatt=1;iwatt <= nnow; ++iwatt){
      *vbendt += intra_scr_vbend[iwatt];
      *vbondt += intra_scr_vbond[iwatt];
    }/*endfor*/
    
/*==========================================================================*/
/* X) Sum the pressure tensor                                               */

    if(cell->iperd == 2 || cell->iperd == 3) {
      for(iwatt=1;iwatt <= nnow; ++iwatt){
        ptens_pvten_tmp[1] += intra_scr_p11[iwatt];
        ptens_pvten_tmp[5] += intra_scr_p22[iwatt];
        ptens_pvten_tmp[9] += intra_scr_p33[iwatt];
        ptens_pvten_tmp[2] += intra_scr_p12[iwatt];
      }/*endfor*/
    }/*endif*/
    if(cell->iperd == 3) {
      for(iwatt=1;iwatt <= nnow; ++iwatt){
        ptens_pvten_tmp[3] += intra_scr_p13[iwatt];
        ptens_pvten_tmp[6] += intra_scr_p23[iwatt];
      }/*endfor*/
    }/*endif*/

/*==========================================================================*/
/* XI) Scatter the forces                                                   */

if(iver_get==1){

    for(iboff=ibig+ioff,iwatt=ibig;iwatt <= iend; ++iwatt,++iboff){
      ktemp = j1[iwatt];
      clatoms_fxt[ktemp] += intra_scr_fx1[iboff];
      clatoms_fyt[ktemp] += intra_scr_fy1[iboff];
      clatoms_fzt[ktemp] += intra_scr_fz1[iboff];
    }/*endfor*/

    for(iboff=ibig+ioff,iwatt=ibig;iwatt <= iend; ++iwatt,++iboff){
      ktemp = j2[iwatt];
      clatoms_fxt[ktemp] += intra_scr_fx2[iboff];
      clatoms_fyt[ktemp] += intra_scr_fy2[iboff];
      clatoms_fzt[ktemp] += intra_scr_fz2[iboff];
    }/*endfor*/

    for(iboff=ibig+ioff,iwatt=ibig;iwatt <= iend; ++iwatt,++iboff){
      ktemp = j3[iwatt];
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

    for(iboff=ibig+ioff,iwatt=ibig;iwatt <= iend; ++iwatt,++iboff){
      ktemp = j1[iwatt];
      clatoms_fx[ktemp] += intra_scr_fx1[iboff];
      clatoms_fy[ktemp] += intra_scr_fy1[iboff];
      clatoms_fz[ktemp] += intra_scr_fz1[iboff];
    }/*endfor*/

    for(iboff=ibig+ioff,iwatt=ibig;iwatt <= iend; ++iwatt,++iboff){
      ktemp = j2[iwatt];
      clatoms_fx[ktemp] += intra_scr_fx2[iboff];
      clatoms_fy[ktemp] += intra_scr_fy2[iboff];
      clatoms_fz[ktemp] += intra_scr_fz2[iboff];
    }/*endfor*/

    for(iboff=ibig+ioff,iwatt=ibig;iwatt <= iend; ++iwatt,++iboff){
      ktemp = j3[iwatt];
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

  *vtot = *vbondt + *vbendt;

  /*----------------------------------------------------------------------*/
}/*end routine*/
/*==========================================================================*/








