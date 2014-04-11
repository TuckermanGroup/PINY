/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
/*                                                                          */
/*                         PI_MD:                                           */
/*             The future of simulation technology                          */
/*             ------------------------------------                         */
/*                   Module: grp_bond_con.c                                 */
/*                                                                          */
/* This routine controls the atom based group constraint routines           */
/*                                                                          */
/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

#include "standard_include.h"
#include "../typ_defs/typedefs_class.h"
#include "../typ_defs/typedefs_gen.h"
#include "../typ_defs/typedefs_bnd.h"
#include "../proto_defs/proto_intra_con_local.h"
#include "../proto_defs/proto_friend_lib_entry.h"
#include "../proto_defs/proto_math.h"

#define NCON_23 2
#define NAT_23 3


/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
void shake_23(GRP_BOND_CON *grp_bond_con,
              CLATOMS_INFO *clatoms_info,CLATOMS_POS *clatoms_pos,
              PTENS *ptens,double dt,double *aiter,
              CLASS_COMM_FORC_PKG *class_comm_forc_pkg)


/*==========================================================================*/
/*        Begin Routine                                                     */
{/*Begin Routine*/
/*=======================================================================*/
/*         Local Variable declarations                                   */
  
 double xl0[NCON_23+1];
 double rmu[NCON_23+1];
 double amat[NCON_23+1][NCON_23+1],ainv[NCON_23+1][NCON_23+1];
 double r12s,r13s;
 double dxn12,dyn12,dzn12,dxn13,dyn13,dzn13;
 double dlmax,dlmax1,dlmax2,rdet_a;
 double rmu12,rmu22,rm12;

 double dts;
 
 int iter,igrp,*ind1,*ind2,*ind3,jtyp;
 int iii,ktemp;

 double *xlam1,*xlam2;    
 double *avec1,*avec2;
 double *rm1,*rm2,*rm3;
 double *rmm11,*rmm12,*rmm21,*rmm22;
 double *dxl1,*dxl2;
 double *dx12,*dy12,*dz12,*dx13,*dy13,*dz13;
 double *dxo12,*dyo12,*dzo12,*dxo13,*dyo13,*dzo13;
 double **x,**y,**z,**xo,**yo,**zo;
 double *dij1,*dij2;
 double *p11,*p12,*p13,*p23,*p33,*p22;
/* Local pointers */
  double *clatoms_mass         = clatoms_info->mass;
  double *clatoms_x            = clatoms_pos->x;
  double *clatoms_y            = clatoms_pos->y;
  double *clatoms_z            = clatoms_pos->z;
  double *clatoms_vx           = clatoms_pos->vx;
  double *clatoms_vy           = clatoms_pos->vy;
  double *clatoms_vz           = clatoms_pos->vz;
  double *clatoms_xold         = clatoms_info->xold;
  double *clatoms_yold         = clatoms_info->yold;
  double *clatoms_zold         = clatoms_info->zold;
  int *grp_bond_con_j1_23      = grp_bond_con->j1_23;
  int *grp_bond_con_j2_23      = grp_bond_con->j2_23;
  int *grp_bond_con_j3_23      = grp_bond_con->j3_23;
  int *grp_bond_con_jtyp_23    = grp_bond_con->jtyp_23;
  double **grp_bond_con_eq_23  = grp_bond_con->eq_23;
  double **grp_bond_con_al_23  = grp_bond_con->al_23;
  double *ptens_pvten_inc      = ptens->pvten_inc;
  double pnorm;


  int ngrp,irem,igrp_off;
  int ngrp_tot                 = grp_bond_con->num_23;
  int np_forc                  = class_comm_forc_pkg->num_proc;
  int myid_forc                = class_comm_forc_pkg->myid;

/*=======================================================================*/

  ngrp = (ngrp_tot);
  igrp_off = 0;

/*=======================================================================*/

     xlam1= dvector(1,ngrp);
     xlam2= dvector(1,ngrp);

     avec1= dvector(1,ngrp);
     avec2= dvector(1,ngrp);
 
     rm1 = dvector(1,ngrp);
     rm2 = dvector(1,ngrp);
     rm3 = dvector(1,ngrp);

     rmm11 = dvector(1,ngrp);
     rmm12 = dvector(1,ngrp);
     rmm21 = dvector(1,ngrp);
     rmm22 = dvector(1,ngrp);

     dxl1= dvector(1,ngrp);
     dxl2= dvector(1,ngrp);
 
     dx12= dvector(1,ngrp);
     dy12= dvector(1,ngrp);
     dz12= dvector(1,ngrp);

     dx13= dvector(1,ngrp);
     dy13= dvector(1,ngrp);
     dz13= dvector(1,ngrp);

     dxo12= dvector(1,ngrp);
     dyo12= dvector(1,ngrp);
     dzo12= dvector(1,ngrp);

     dxo13= dvector(1,ngrp);
     dyo13= dvector(1,ngrp);
     dzo13= dvector(1,ngrp);

      dij1= dvector(1,ngrp);
      dij2= dvector(1,ngrp);

         x= dmatrix(1,3,1,ngrp);
         y= dmatrix(1,3,1,ngrp);
         z= dmatrix(1,3,1,ngrp);

        xo= dmatrix(1,3,1,ngrp);
        yo= dmatrix(1,3,1,ngrp);
        zo= dmatrix(1,3,1,ngrp);

       p11= dvector(1,ngrp);
       p12= dvector(1,ngrp);
       p13= dvector(1,ngrp);
       p22= dvector(1,ngrp);
       p23= dvector(1,ngrp);
       p33= dvector(1,ngrp);

      ind1 = (int *) cmalloc(ngrp*sizeof(int))-1;
      ind2 = (int *) cmalloc(ngrp*sizeof(int))-1;
      ind3 = (int *) cmalloc(ngrp*sizeof(int))-1;
/*=======================================================================*/
/*=======================================================================*/

 dts = dt*dt;
 pnorm = 2.0/dts;
 *aiter = 0.0;
/* assign positions and masses */
 for(igrp=1;igrp <= ngrp; igrp++) {
   ind1[igrp] = grp_bond_con_j1_23[(igrp+igrp_off)];
   ind2[igrp] = grp_bond_con_j2_23[(igrp+igrp_off)];
   ind3[igrp] = grp_bond_con_j3_23[(igrp+igrp_off)];
 }

 for(igrp=1;igrp <= ngrp; igrp++) {
   ktemp = ind1[igrp];
   x[1][igrp] = clatoms_x[ktemp];
   y[1][igrp] = clatoms_y[ktemp];
   z[1][igrp] = clatoms_z[ktemp];
 }

 for(igrp=1;igrp <= ngrp; igrp++) {
   ktemp = ind2[igrp];
   x[2][igrp] = clatoms_x[ktemp];
   y[2][igrp] = clatoms_y[ktemp];
   z[2][igrp] = clatoms_z[ktemp];
 }

 for(igrp=1;igrp <= ngrp; igrp++) {
   ktemp = ind3[igrp];
   x[3][igrp] = clatoms_x[ktemp];
   y[3][igrp] = clatoms_y[ktemp];
   z[3][igrp] = clatoms_z[ktemp];
 }

 for(igrp=1;igrp <= ngrp; igrp++) {
   ktemp = ind1[igrp];
   xo[1][igrp] = clatoms_xold[ktemp];
   yo[1][igrp] = clatoms_yold[ktemp];
   zo[1][igrp] = clatoms_zold[ktemp];
   rm1[igrp] = 1.0/clatoms_mass[ktemp];
 }

 for(igrp=1;igrp <= ngrp; igrp++) {
   ktemp = ind2[igrp];
   xo[2][igrp] = clatoms_xold[ktemp];
   yo[2][igrp] = clatoms_yold[ktemp];
   zo[2][igrp] = clatoms_zold[ktemp];
   rm2[igrp] = 1.0/clatoms_mass[ktemp];
 }

 for(igrp=1;igrp <= ngrp; igrp++) {
   ktemp = ind3[igrp];
   xo[3][igrp] = clatoms_xold[ktemp];
   yo[3][igrp] = clatoms_yold[ktemp];
   zo[3][igrp] = clatoms_zold[ktemp];
   rm3[igrp] = 1.0/clatoms_mass[ktemp];
 }


 for(igrp=1;igrp <= ngrp; igrp++) {
   jtyp = grp_bond_con_jtyp_23[(igrp+igrp_off)];
   dij1[igrp] = grp_bond_con_eq_23[1][jtyp];
   dij2[igrp] = grp_bond_con_eq_23[2][jtyp];
 }/* end for*/

/* Initial Guess for multipliers */
 for(igrp=1;igrp <= ngrp; igrp++) {
    rmm11[igrp] = -(rm1[igrp]+rm2[igrp]);
    rmm12[igrp] = -rm1[igrp];
    rmm21[igrp] = -rm1[igrp];
    rmm22[igrp] = -(rm1[igrp]+rm3[igrp]);
 }
 for(igrp=1;igrp <= ngrp; igrp++) {
    dx12[igrp] = x[1][igrp]-x[2][igrp];
    dy12[igrp] = y[1][igrp]-y[2][igrp];
    dz12[igrp] = z[1][igrp]-z[2][igrp];
  }

 for(igrp=1;igrp <= ngrp; igrp++) {
    dx13[igrp] = x[1][igrp]-x[3][igrp];
    dy13[igrp] = y[1][igrp]-y[3][igrp];
    dz13[igrp] = z[1][igrp]-z[3][igrp];
  }
 for(igrp=1;igrp <= ngrp; igrp++) {
    dxo12[igrp] = xo[1][igrp]-xo[2][igrp];
    dyo12[igrp] = yo[1][igrp]-yo[2][igrp];
    dzo12[igrp] = zo[1][igrp]-zo[2][igrp];
  }

 for(igrp=1;igrp <= ngrp; igrp++) {
    dxo13[igrp] = xo[1][igrp]-xo[3][igrp];
    dyo13[igrp] = yo[1][igrp]-yo[3][igrp];
    dzo13[igrp] = zo[1][igrp]-zo[3][igrp];
  }

 for(igrp=1;igrp <= ngrp; igrp++) {

   amat[1][1] = 2.0*rmm11[igrp]*(dx12[igrp]*dxo12[igrp]
              + dy12[igrp]*dyo12[igrp]
              + dz12[igrp]*dzo12[igrp]);
   amat[1][2] = 2.0*rmm12[igrp]*(dx12[igrp]*dxo13[igrp]
              + dy12[igrp]*dyo13[igrp]
              + dz12[igrp]*dzo13[igrp]);
   amat[2][1] = 2.0*rmm21[igrp]*(dx13[igrp]*dxo12[igrp]
              + dy13[igrp]*dyo12[igrp]
              + dz13[igrp]*dzo12[igrp]);
   amat[2][2] = 2.0*rmm22[igrp]*(dx13[igrp]*dxo13[igrp]
              + dy13[igrp]*dyo13[igrp]
              + dz13[igrp]*dzo13[igrp]);

   r12s = dx12[igrp]*dx12[igrp] + dy12[igrp]*dy12[igrp]
        + dz12[igrp]*dz12[igrp];
   r13s = dx13[igrp]*dx13[igrp] + dy13[igrp]*dy13[igrp]
        + dz13[igrp]*dz13[igrp];
   avec1[igrp] = dij1[igrp]*dij1[igrp] - r12s;
   avec2[igrp] = dij2[igrp]*dij2[igrp] - r13s;

   rdet_a = 1.0/(amat[1][1]*amat[2][2] - amat[1][2]*amat[2][1]);
   ainv[1][1] =  amat[2][2]*rdet_a;
   ainv[1][2] = -amat[1][2]*rdet_a;
   ainv[2][1] = -amat[2][1]*rdet_a;
   ainv[2][2] =  amat[1][1]*rdet_a;

   xlam1[igrp] = ainv[1][1]*avec1[igrp] + ainv[1][2]*avec2[igrp];
   xlam2[igrp] = ainv[2][1]*avec1[igrp] + ainv[2][2]*avec2[igrp];


 } /* first loop over groups */

/* Iterative do loop for multiplier */
 dlmax = 1.0;
 iter = 0;
 do {
  ++iter;
  if(iter > grp_bond_con->max_iter) {
    printf("$$$$$$$$$$$$$$$$$$$$_WARNING_$$$$$$$$$$$$$$$$$$$$\n");
    printf("Shake not converged after %d iterations.\n",
            grp_bond_con->max_iter);
    printf("The present tolerance is %g \n",dlmax);
    printf("The desired tolerance is %g \n",grp_bond_con->tol);
    printf("$$$$$$$$$$$$$$$$$$$$_WARNING_$$$$$$$$$$$$$$$$$$$$\n");
    fflush(stdout);
    break;
  }/*endif*/
/* Get elements of matrix */
 for(igrp=1;igrp <= ngrp; igrp++) {

  dxn12 = 2.0*dx12[igrp] 
        + rmm11[igrp]*xlam1[igrp]*dxo12[igrp]
        + rmm12[igrp]*xlam2[igrp]*dxo13[igrp];
  dyn12 = 2.0*dy12[igrp] 
        + rmm11[igrp]*xlam1[igrp]*dyo12[igrp]
        + rmm12[igrp]*xlam2[igrp]*dyo13[igrp];
  dzn12 = 2.0*dz12[igrp] 
        + rmm11[igrp]*xlam1[igrp]*dzo12[igrp]
        + rmm12[igrp]*xlam2[igrp]*dzo13[igrp];
  dxn13 = 2.0*dx13[igrp] 
        + rmm21[igrp]*xlam1[igrp]*dxo12[igrp]
        + rmm22[igrp]*xlam2[igrp]*dxo13[igrp];
  dyn13 = 2.0*dy13[igrp] 
        + rmm21[igrp]*xlam1[igrp]*dyo12[igrp]
        + rmm22[igrp]*xlam2[igrp]*dyo13[igrp];
  dzn13 = 2.0*dz13[igrp] 
        + rmm21[igrp]*xlam1[igrp]*dzo12[igrp]
        + rmm22[igrp]*xlam2[igrp]*dzo13[igrp];

   amat[1][1] = rmm11[igrp]*(dxn12*dxo12[igrp]
              + dyn12*dyo12[igrp]
              + dzn12*dzo12[igrp]);
   amat[1][2] = rmm12[igrp]*(dxn12*dxo13[igrp]
              + dyn12*dyo13[igrp]
              + dzn12*dzo13[igrp]);
   amat[2][1] = rmm21[igrp]*(dxn13*dxo12[igrp]
              + dyn13*dyo12[igrp]
              + dzn13*dzo12[igrp]);
   amat[2][2] = rmm22[igrp]*(dxn13*dxo13[igrp]
              + dyn13*dyo13[igrp]
              + dzn13*dzo13[igrp]);
 
  rdet_a = 1.0/(amat[1][1]*amat[2][2] - amat[1][2]*amat[2][1]);
  ainv[1][1] =  amat[2][2]*rdet_a;
  ainv[1][2] = -amat[1][2]*rdet_a;
  ainv[2][1] = -amat[2][1]*rdet_a;
  ainv[2][2] =  amat[1][1]*rdet_a;

  xl0[1] = xlam1[igrp];
  xl0[2] = xlam2[igrp];

  xlam1[igrp] = ainv[1][1]*avec1[igrp] + ainv[1][2]*avec2[igrp];
  xlam2[igrp] = ainv[2][1]*avec1[igrp] + ainv[2][2]*avec2[igrp];

  dxl1[igrp] = fabs(xlam1[igrp]-xl0[1]);
  dxl2[igrp] = fabs(xlam2[igrp]-xl0[2]);
 } /* end loop over groups */

/* test for largest */

   dlmax1= dxl1[1];
   dlmax2= dxl2[1];
 for(igrp=2;igrp <= ngrp; igrp++) {
     dlmax1= (dlmax1 > dxl1[igrp] ? dlmax1:dxl1[igrp]);
     dlmax2= (dlmax2 > dxl2[igrp] ? dlmax2:dxl2[igrp]);
 }
   dlmax = (dlmax1 > dlmax2 ? dlmax1:dlmax2);

 } while(dlmax > grp_bond_con->tol);
 *aiter += (double) iter;

/* position update */
#ifndef NO_PRAGMA
#pragma IVDEP
#endif
 for(igrp=1;igrp <= ngrp; igrp++) {
  ktemp = ind1[igrp];
  clatoms_x[ktemp] -= (xlam1[igrp]*dxo12[igrp] + xlam2[igrp]*dxo13[igrp])*rm1[igrp];
  clatoms_y[ktemp] -= (xlam1[igrp]*dyo12[igrp] + xlam2[igrp]*dyo13[igrp])*rm1[igrp];
  clatoms_z[ktemp] -= (xlam1[igrp]*dzo12[igrp] + xlam2[igrp]*dzo13[igrp])*rm1[igrp];
 }

#ifndef NO_PRAGMA
#pragma IVDEP
#endif
 for(igrp=1;igrp <= ngrp; igrp++) {
  ktemp = ind2[igrp];
  clatoms_x[ktemp] += xlam1[igrp]*dxo12[igrp]*rm2[igrp];
  clatoms_y[ktemp] += xlam1[igrp]*dyo12[igrp]*rm2[igrp];
  clatoms_z[ktemp] += xlam1[igrp]*dzo12[igrp]*rm2[igrp];
 }

#ifndef NO_PRAGMA
#pragma IVDEP
#endif
 for(igrp=1;igrp <= ngrp; igrp++) {
  ktemp = ind3[igrp];
  clatoms_x[ktemp] += xlam2[igrp]*dxo13[igrp]*rm3[igrp];
  clatoms_y[ktemp] += xlam2[igrp]*dyo13[igrp]*rm3[igrp];
  clatoms_z[ktemp] += xlam2[igrp]*dzo13[igrp]*rm3[igrp];
 }

/* Velocity update */

#ifndef NO_PRAGMA
#pragma IVDEP
#endif
 for(igrp=1;igrp <= ngrp; igrp++) {
  ktemp = ind1[igrp];
  clatoms_vx[ktemp]-=(xlam1[igrp]*dxo12[igrp]+xlam2[igrp]*dxo13[igrp])*rm1[igrp]/dt;
  clatoms_vy[ktemp]-=(xlam1[igrp]*dyo12[igrp]+xlam2[igrp]*dyo13[igrp])*rm1[igrp]/dt;
  clatoms_vz[ktemp]-=(xlam1[igrp]*dzo12[igrp]+xlam2[igrp]*dzo13[igrp])*rm1[igrp]/dt;
 }

#ifndef NO_PRAGMA
#pragma IVDEP
#endif
 for(igrp=1;igrp <= ngrp; igrp++) {
  ktemp = ind2[igrp];
  clatoms_vx[ktemp] += xlam1[igrp]*dxo12[igrp]*rm2[igrp]/dt;
  clatoms_vy[ktemp] += xlam1[igrp]*dyo12[igrp]*rm2[igrp]/dt;
  clatoms_vz[ktemp] += xlam1[igrp]*dzo12[igrp]*rm2[igrp]/dt;
 }

#ifndef NO_PRAGMA
#pragma IVDEP
#endif
 for(igrp=1;igrp <= ngrp; igrp++) {
  ktemp = ind3[igrp];
  clatoms_vx[ktemp] += xlam2[igrp]*dxo13[igrp]*rm3[igrp]/dt;
  clatoms_vy[ktemp] += xlam2[igrp]*dyo13[igrp]*rm3[igrp]/dt;
  clatoms_vz[ktemp] += xlam2[igrp]*dzo13[igrp]*rm3[igrp]/dt;
 }

/* Pressure tensor update */
/* Compute difference vectors */
 for(igrp=1;igrp <= ngrp; igrp++) {
   dxn12 = dxo12[igrp];
   dyn12 = dyo12[igrp];
   dzn12 = dzo12[igrp];
   dxn13 = dxo13[igrp];
   dyn13 = dyo13[igrp];
   dzn13 = dzo13[igrp];

   p11[igrp] = xlam1[igrp]*dxn12*dxn12 + xlam2[igrp]*dxn13*dxn13;
   p22[igrp] = xlam1[igrp]*dyn12*dyn12 + xlam2[igrp]*dyn13*dyn13;
   p33[igrp] = xlam1[igrp]*dzn12*dzn12 + xlam2[igrp]*dzn13*dzn13;
   p12[igrp] = xlam1[igrp]*dxn12*dyn12 + xlam2[igrp]*dxn13*dyn13;
   p13[igrp] = xlam1[igrp]*dxn12*dzn12 + xlam2[igrp]*dxn13*dzn13;
   p23[igrp] = xlam1[igrp]*dyn12*dzn12 + xlam2[igrp]*dyn13*dzn13;
 }/*end for*/

#ifndef NO_PRAGMA
#pragma IVDEP
#endif
 for(igrp=1;igrp <= ngrp; igrp++) {
    ptens_pvten_inc[1] -= (p11[igrp]*pnorm); 
    ptens_pvten_inc[2] -= (p12[igrp]*pnorm);
    ptens_pvten_inc[3] -= (p13[igrp]*pnorm);
    ptens_pvten_inc[4] -= (p12[igrp]*pnorm);
    ptens_pvten_inc[5] -= (p22[igrp]*pnorm);
    ptens_pvten_inc[6] -= (p23[igrp]*pnorm);
    ptens_pvten_inc[7] -= (p13[igrp]*pnorm);
    ptens_pvten_inc[8] -= (p23[igrp]*pnorm);
    ptens_pvten_inc[9] -= (p33[igrp]*pnorm);
 }/*end for*/

/* Save multiplier */

 for(igrp=1;igrp <= ngrp; igrp++) {
   grp_bond_con_al_23[1][(igrp+igrp_off)] = xlam1[igrp];
   grp_bond_con_al_23[2][(igrp+igrp_off)] = xlam2[igrp];
 } /* end for igrp */

/* free locally assigned memory */
     free_dvector(xlam1,1,ngrp);
     free_dvector(xlam2,1,ngrp);

     free_dvector(avec1,1,ngrp);
     free_dvector(avec2,1,ngrp);
 
     free_dvector(rm1,1,ngrp);
     free_dvector(rm2,1,ngrp);
     free_dvector(rm3,1,ngrp);

     free_dvector(rmm11,1,ngrp);
     free_dvector(rmm12,1,ngrp);
     free_dvector(rmm21,1,ngrp);
     free_dvector(rmm22,1,ngrp);

     free_dvector(dxl1,1,ngrp);
     free_dvector(dxl2,1,ngrp);

     free_dvector(dx12,1,ngrp);
     free_dvector(dy12,1,ngrp);
     free_dvector(dz12,1,ngrp);

     free_dvector(dx13,1,ngrp);
     free_dvector(dy13,1,ngrp);
     free_dvector(dz13,1,ngrp);

     free_dvector(dxo12,1,ngrp);
     free_dvector(dyo12,1,ngrp);
     free_dvector(dzo12,1,ngrp);

     free_dvector(dxo13,1,ngrp);
     free_dvector(dyo13,1,ngrp);
     free_dvector(dzo13,1,ngrp);

    free_dvector(dij1,1,ngrp);
    free_dvector(dij2,1,ngrp);

    free_dmatrix(x,1,3,1,ngrp); 
    free_dmatrix(y,1,3,1,ngrp); 
    free_dmatrix(z,1,3,1,ngrp); 

    free_dmatrix(xo,1,3,1,ngrp); 
    free_dmatrix(yo,1,3,1,ngrp); 
    free_dmatrix(zo,1,3,1,ngrp); 

    free_dvector(p11,1,ngrp);
    free_dvector(p12,1,ngrp);
    free_dvector(p13,1,ngrp);
    free_dvector(p22,1,ngrp);
    free_dvector(p23,1,ngrp);
    free_dvector(p33,1,ngrp);


    cfree(&(ind1[1]));
    cfree(&(ind2[1]));
    cfree(&(ind3[1]));


/*=======================================================================*/
} /* end routine */
/*=======================================================================*/


/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void rattle_23(GRP_BOND_CON *grp_bond_con,
              CLATOMS_INFO *clatoms_info,CLATOMS_POS *clatoms_pos,
              PTENS *ptens,double dt,
              CLASS_COMM_FORC_PKG *class_comm_forc_pkg)

/*==========================================================================*/
/*        Begin Routine                                                     */
{/* Begin routine */
/*=======================================================================*/
/*         Local Variable declarations                                   */

  double avec[NCON_23+1];
  double amat[NCON_23+1][NCON_23+1],ainv[NCON_23+1][NCON_23+1];
  double rmu1,rmu2;


  double dvx12,dvy12,dvz12,dvx13,dvy13,dvz13;
  double r12s,r13s,dot23,dot2t2,dot3t3;
  double rdet_a;

  int igrp,*ind1,*ind2,*ind3,jtyp;
  int ktemp;

  double **x,**y,**z;
  double **vx,**vy,**vz;
  double *p11,*p12,*p13,*p22,*p23,*p33;
  double *rm1,*rm2,*rm3;
  double *dx12,*dy12,*dz12,*dx13,*dy13,*dz13;
  double **xlam;

/* Local pointers */
  double *clatoms_mass         = clatoms_info->mass;
  double *clatoms_x            = clatoms_pos->x;
  double *clatoms_y            = clatoms_pos->y;
  double *clatoms_z            = clatoms_pos->z;
  double *clatoms_vx           = clatoms_pos->vx;
  double *clatoms_vy           = clatoms_pos->vy;
  double *clatoms_vz           = clatoms_pos->vz;
  int *grp_bond_con_j1_23      = grp_bond_con->j1_23;
  int *grp_bond_con_j2_23      = grp_bond_con->j2_23;
  int *grp_bond_con_j3_23      = grp_bond_con->j3_23;
  int *grp_bond_con_jtyp_23    = grp_bond_con->jtyp_23;
  double *ptens_pvten_inc      = ptens->pvten_inc;
  double pnorm;

  int ngrp,irem,igrp_off;
  int ngrp_tot                 = grp_bond_con->num_23;
  int np_forc                  = class_comm_forc_pkg->num_proc;
  int myid_forc                = class_comm_forc_pkg->myid;

/*=======================================================================*/

  ngrp = (ngrp_tot);
  igrp_off = 0;

/*=======================================================================*/
 
/* assign local memory */
    x= dmatrix(1,3,1,ngrp);
    y= dmatrix(1,3,1,ngrp);
    z= dmatrix(1,3,1,ngrp);
   vx= dmatrix(1,3,1,ngrp);
   vy= dmatrix(1,3,1,ngrp);
   vz= dmatrix(1,3,1,ngrp);
  p11= dvector(1,ngrp);
  p12= dvector(1,ngrp);
  p13= dvector(1,ngrp);
  p22= dvector(1,ngrp);
  p23= dvector(1,ngrp);
  p33= dvector(1,ngrp);
  rm1= dvector(1,ngrp);
  rm2= dvector(1,ngrp);
  rm3= dvector(1,ngrp);
  dx12= dvector(1,ngrp);
  dy12= dvector(1,ngrp);
  dz12= dvector(1,ngrp);
  dx13= dvector(1,ngrp);
  dy13= dvector(1,ngrp);
  dz13= dvector(1,ngrp);
 xlam= dmatrix(1,2,1,ngrp);
 ind1 = (int *) cmalloc(ngrp*sizeof(int))-1;
 ind2 = (int *) cmalloc(ngrp*sizeof(int))-1;
 ind3 = (int *) cmalloc(ngrp*sizeof(int))-1;
/*=======================================================================*/
/*=======================================================================*/

 pnorm = 2.0/dt;

 for(igrp=1;igrp <= ngrp; igrp++) {
 
  ind1[igrp] = grp_bond_con_j1_23[(igrp+igrp_off)];
  ind2[igrp] = grp_bond_con_j2_23[(igrp+igrp_off)];
  ind3[igrp] = grp_bond_con_j3_23[(igrp+igrp_off)];
}

 for(igrp=1;igrp <= ngrp; igrp++) {
  rm1[igrp] = 1.0/clatoms_mass[ind1[igrp]];
  rm2[igrp] = 1.0/clatoms_mass[ind2[igrp]];
  rm3[igrp] = 1.0/clatoms_mass[ind3[igrp]];
}

 for(igrp=1;igrp <= ngrp; igrp++) {
  ktemp = ind1[igrp];
  x[1][igrp] = clatoms_x[ktemp];
  y[1][igrp] = clatoms_y[ktemp];
  z[1][igrp] = clatoms_z[ktemp];
} /*end for*/

 for(igrp=1;igrp <= ngrp; igrp++) {
  ktemp = ind2[igrp];
  x[2][igrp] = clatoms_x[ktemp];
  y[2][igrp] = clatoms_y[ktemp];
  z[2][igrp] = clatoms_z[ktemp];
} /*end for*/

 for(igrp=1;igrp <= ngrp; igrp++) {
  ktemp = ind3[igrp];
  x[3][igrp] = clatoms_x[ktemp];
  y[3][igrp] = clatoms_y[ktemp];
  z[3][igrp] = clatoms_z[ktemp];
} /*end for*/

 for(igrp=1;igrp <= ngrp; igrp++) {
  ktemp = ind1[igrp];
  vx[1][igrp] = clatoms_vx[ktemp];
  vy[1][igrp] = clatoms_vy[ktemp];
  vz[1][igrp] = clatoms_vz[ktemp];
} /*end for*/

 for(igrp=1;igrp <= ngrp; igrp++) {
  ktemp = ind2[igrp];
  vx[2][igrp] = clatoms_vx[ktemp];
  vy[2][igrp] = clatoms_vy[ktemp];
  vz[2][igrp] = clatoms_vz[ktemp];
} /*end for*/

 for(igrp=1;igrp <= ngrp; igrp++) {
  ktemp = ind3[igrp];
  vx[3][igrp] = clatoms_vx[ktemp];
  vy[3][igrp] = clatoms_vy[ktemp];
  vz[3][igrp] = clatoms_vz[ktemp];
} /*end for*/

 for(igrp=1;igrp <= ngrp; igrp++) {
   rmu1 = rm1[igrp] + rm2[igrp];
   rmu2 = rm1[igrp] + rm3[igrp];

/* Define useful constants */
  
 dx12[igrp] = x[1][igrp]-x[2][igrp];
 dy12[igrp] = y[1][igrp]-y[2][igrp];
 dz12[igrp] = z[1][igrp]-z[2][igrp];
 r12s = dx12[igrp]*dx12[igrp] + dy12[igrp]*dy12[igrp] + dz12[igrp]*dz12[igrp];

 dx13[igrp] = x[1][igrp]-x[3][igrp];
 dy13[igrp] = y[1][igrp]-y[3][igrp];
 dz13[igrp] = z[1][igrp]-z[3][igrp];
 r13s = dx13[igrp]*dx13[igrp] + dy13[igrp]*dy13[igrp]+ dz13[igrp]*dz13[igrp];

 dvx12 = vx[1][igrp]-vx[2][igrp];
 dvy12 = vy[1][igrp]-vy[2][igrp];
 dvz12 = vz[1][igrp]-vz[2][igrp];

 dvx13 = vx[1][igrp]-vx[3][igrp];
 dvy13 = vy[1][igrp]-vy[3][igrp];
 dvz13 = vz[1][igrp]-vz[3][igrp];

/* Get elements of vector */

 dot2t2 = dx12[igrp]*dvx12 + dy12[igrp]*dvy12 + dz12[igrp]*dvz12;
 dot3t3 = dx13[igrp]*dvx13 + dy13[igrp]*dvy13 + dz13[igrp]*dvz13;
 dot23 = dx12[igrp]*dx13[igrp] + dy12[igrp]*dy13[igrp] + dz12[igrp]*dz13[igrp];

 avec[1] = dot2t2;
 avec[2] = dot3t3;

/* Get elements of matrix */

 amat[1][1] = rmu1*r12s;
 amat[1][2] = rm1[igrp]*dot23;
 amat[2][1] = amat[1][2];
 amat[2][2] = rmu2*r13s;

 rdet_a = 1.0/(amat[1][1]*amat[2][2] - amat[1][2]*amat[2][1]);
 ainv[1][1] =  amat[2][2]*rdet_a;
 ainv[1][2] = -amat[1][2]*rdet_a;
 ainv[2][1] = -amat[2][1]*rdet_a;
 ainv[2][2] =  amat[1][1]*rdet_a;

 xlam[1][igrp] = ainv[1][1]*avec[1] + ainv[1][2]*avec[2];
 xlam[2][igrp] = ainv[2][1]*avec[1] + ainv[2][2]*avec[2];
}/* endfor */

#ifndef NO_PRAGMA
#pragma IVDEP
#endif
 for(igrp=1;igrp <= ngrp; igrp++) {
   double xlam1,xlam2;
   double dx_12,dx_13;
   double dy_12,dy_13;
   double dz_12,dz_13;
   int ktemp1,ktemp2,ktemp3;

    xlam1= xlam[1][igrp];
    xlam2= xlam[2][igrp];

    dx_12= dx12[igrp]; dx_13= dx13[igrp];
    dy_12= dy12[igrp]; dy_13= dy13[igrp];
    dz_12= dz12[igrp]; dz_13= dz13[igrp];
    ktemp1=ind1[igrp]; ktemp2=ind2[igrp]; ktemp3=ind3[igrp];

   clatoms_vx[ktemp1] -= (xlam1*dx_12 + xlam2*dx_13)*rm1[igrp];
   clatoms_vy[ktemp1] -= (xlam1*dy_12 + xlam2*dy_13)*rm1[igrp];
   clatoms_vz[ktemp1] -= (xlam1*dz_12 + xlam2*dz_13)*rm1[igrp];

   clatoms_vx[ktemp2] += xlam1*dx_12*rm2[igrp];
   clatoms_vy[ktemp2] += xlam1*dy_12*rm2[igrp];
   clatoms_vz[ktemp2] += xlam1*dz_12*rm2[igrp];

   clatoms_vx[ktemp3] += xlam2*dx_13*rm3[igrp];
   clatoms_vy[ktemp3] += xlam2*dy_13*rm3[igrp];
   clatoms_vz[ktemp3] += xlam2*dz_13*rm3[igrp];

/* Pressure Tensor update */

    p11[igrp] = xlam1*dx_12*dx_12 + xlam2*dx_13*dx_13;
    p22[igrp] = xlam1*dy_12*dy_12 + xlam2*dy_13*dy_13;
    p33[igrp] = xlam1*dz_12*dz_12 + xlam2*dz_13*dz_13;
    p12[igrp] = xlam1*dx_12*dy_12 + xlam2*dx_13*dy_13;
    p13[igrp] = xlam1*dx_12*dz_12 + xlam2*dx_13*dz_13;
    p23[igrp] = xlam1*dy_12*dz_12 + xlam2*dy_13*dz_13;
}/*end for*/

#ifndef NO_PRAGMA
#pragma IVDEP
#endif
 for(igrp=1;igrp <= ngrp; igrp++) {
    ptens_pvten_inc[1] -= (p11[igrp]*pnorm); 
    ptens_pvten_inc[2] -= (p12[igrp]*pnorm);
    ptens_pvten_inc[3] -= (p13[igrp]*pnorm);
    ptens_pvten_inc[4] -= (p12[igrp]*pnorm);
    ptens_pvten_inc[5] -= (p22[igrp]*pnorm);
    ptens_pvten_inc[6] -= (p23[igrp]*pnorm);
    ptens_pvten_inc[7] -= (p13[igrp]*pnorm);
    ptens_pvten_inc[8] -= (p23[igrp]*pnorm);
    ptens_pvten_inc[9] -= (p33[igrp]*pnorm);

 }/* end for igrp */
/* free locally assigned memory */  
   free_dmatrix(x,1,3,1,ngrp);
   free_dmatrix(y,1,3,1,ngrp);
   free_dmatrix(z,1,3,1,ngrp);

   free_dmatrix(vx,1,3,1,ngrp);
   free_dmatrix(vy,1,3,1,ngrp);
   free_dmatrix(vz,1,3,1,ngrp);

   free_dvector(p11,1,ngrp);
   free_dvector(p12,1,ngrp);
   free_dvector(p13,1,ngrp);
   free_dvector(p22,1,ngrp);
   free_dvector(p23,1,ngrp);
   free_dvector(p33,1,ngrp);

   free_dvector(rm1,1,ngrp);
   free_dvector(rm2,1,ngrp);
   free_dvector(rm3,1,ngrp);

   free_dvector(dx12,1,ngrp); 
   free_dvector(dy12,1,ngrp); 
   free_dvector(dz12,1,ngrp); 

   free_dvector(dx13,1,ngrp); 
   free_dvector(dy13,1,ngrp); 
   free_dvector(dz13,1,ngrp); 

   free_dmatrix(xlam,1,2,1,ngrp); 


   cfree(&(ind1[1]));
   cfree(&(ind2[1]));
   cfree(&(ind3[1]));



/*=======================================================================*/
} /* end routine */
/*=======================================================================*/




