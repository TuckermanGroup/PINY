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

#define NCON_33 3
#define NAT_33 3


/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
void shake_33(GRP_BOND_CON *grp_bond_con,
              CLATOMS_INFO *clatoms_info,CLATOMS_POS *clatoms_pos,
              PTENS *ptens,double dt,double *aiter,
              CLASS_COMM_FORC_PKG *class_comm_forc_pkg)

/*==========================================================================*/
/*        Begin Routine                                                     */
{/*Begin Routine*/
/*=======================================================================*/
/*         Local Variable declarations                                   */
 double xl0[NCON_33+1];


 double dxn[NCON_33+1],dyn[NCON_33+1],dzn[NCON_33+1];
 double amat[NCON_33+1][NCON_33+1],ainv[NCON_33+1][NCON_33+1];
 double dlmax,dlmax1,dlmax2,dlmax3,det,rdet_a;
 double rms1,rms2,rms3;

 double dts;
 int i,ktemp,ktemp1,ktemp2,ktemp3;
 int iter,igrp,*ind1,*ind2,*ind3,jtyp;
/* Local arrays   */
   double **dx,**dy,**dz;
   double **dxt,**dyt,**dzt;
   double *rm1, *rm2, *rm3;
   double **xlam,**dxl;
   double **avec,**dij;
   double *rmm11,*rmm12, *rmm13, *rmm21, *rmm22, *rmm23;
   double *rmm31, *rmm32, *rmm33;
   double *p11,*p22,*p33,*p12,*p13,*p23;

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
  int *grp_bond_con_j1_33      = grp_bond_con->j1_33;
  int *grp_bond_con_j2_33      = grp_bond_con->j2_33;
  int *grp_bond_con_j3_33      = grp_bond_con->j3_33;
  int *grp_bond_con_jtyp_33    = grp_bond_con->jtyp_33;
  double **grp_bond_con_eq_33  = grp_bond_con->eq_33;
  double **grp_bond_con_al_33  = grp_bond_con->al_33;
  double *ptens_pvten_inc      = ptens->pvten_inc;
  double pnorm;
  double x_1,x_2,x_3,y_1,y_2,y_3,z_1,z_2,z_3,xo_1,xo_2,xo_3;
  double yo_1,yo_2,yo_3,zo_1,zo_2,zo_3;

  int ngrp,irem,igrp_off;
  int ngrp_tot = grp_bond_con->num_33;
  int np_forc                  = class_comm_forc_pkg->num_proc;
  int myid_forc                = class_comm_forc_pkg->myid;

/*=======================================================================*/
  ngrp = (ngrp_tot);
  igrp_off = 0;
/*=======================================================================*/
 
/* assign local arrays */ 
    dx= dmatrix(1,3,1,ngrp);
    dy= dmatrix(1,3,1,ngrp);
    dz= dmatrix(1,3,1,ngrp);
   dxt= dmatrix(1,3,1,ngrp);
   dyt= dmatrix(1,3,1,ngrp);
   dzt= dmatrix(1,3,1,ngrp);
  xlam= dmatrix(1,3,1,ngrp);
  avec= dmatrix(1,3,1,ngrp);
   dxl= dmatrix(1,3,1,ngrp);
   dij= dmatrix(1,3,1,ngrp); 
    rm1= dvector(1,ngrp); rm2= dvector(1,ngrp); rm3= dvector(1,ngrp);
   rmm11= dvector(1,ngrp); rmm12= dvector(1,ngrp); rmm13= dvector(1,ngrp);
   rmm21= dvector(1,ngrp); rmm22= dvector(1,ngrp); rmm23= dvector(1,ngrp);
   rmm31= dvector(1,ngrp); rmm32= dvector(1,ngrp); rmm33= dvector(1,ngrp);
    p11= dvector(1,ngrp);
    p12= dvector(1,ngrp);
    p13= dvector(1,ngrp);
    p22= dvector(1,ngrp);
    p23= dvector(1,ngrp);
    p33= dvector(1,ngrp);

/*=======================================================================*/

/* Set zeta matrix and reciprocal mass matrix */

 dts = dt*dt;
 pnorm = 2.0/dts;
 *aiter = 0.0;
/* I assign positions and masses */

 for(igrp=1;igrp <= ngrp; igrp++) {

/* assign jtyps  */

    jtyp = grp_bond_con_jtyp_33[igrp];
    dij[1][igrp] = grp_bond_con_eq_33[1][jtyp];
    dij[2][igrp] = grp_bond_con_eq_33[2][jtyp];
    dij[3][igrp] = grp_bond_con_eq_33[3][jtyp];

 }/*endfor*/

 for(igrp=1;igrp <= ngrp; igrp++) {

    rm1[igrp]  = 1.0/clatoms_mass[grp_bond_con_j1_33[(igrp+igrp_off)]];
    rm2[igrp]  = 1.0/clatoms_mass[grp_bond_con_j2_33[(igrp+igrp_off)]];
    rm3[igrp]  = 1.0/clatoms_mass[grp_bond_con_j3_33[(igrp+igrp_off)]];

 }/*endfor*/

 for(igrp=1;igrp <= ngrp; igrp++) {

    ktemp1 = grp_bond_con_j1_33[(igrp+igrp_off)];
    x_1  = clatoms_x[ktemp1];
    y_1  = clatoms_y[ktemp1];
    z_1  = clatoms_z[ktemp1];

    ktemp2 = grp_bond_con_j2_33[(igrp+igrp_off)];
    x_2  = clatoms_x[ktemp2];
    y_2  = clatoms_y[ktemp2];
    z_2  = clatoms_z[ktemp2];

    ktemp3 = grp_bond_con_j3_33[(igrp+igrp_off)];
    x_3  = clatoms_x[ktemp3];
    y_3  = clatoms_y[ktemp3];
    z_3  = clatoms_z[ktemp3];

/* Compute difference vectors */
    dxt[1][igrp]= x_1 - x_2;
    dxt[2][igrp]= x_1 - x_3;
    dxt[3][igrp]= x_2 - x_3;
    dyt[1][igrp]= y_1 - y_2;
    dyt[2][igrp]= y_1 - y_3;
    dyt[3][igrp]= y_2 - y_3;
    dzt[1][igrp]= z_1 - z_2;
    dzt[2][igrp]= z_1 - z_3;
    dzt[3][igrp]= z_2 - z_3;

 }/*endfor*/

 for(igrp=1;igrp <= ngrp; igrp++) {

    ktemp1 = grp_bond_con_j1_33[(igrp+igrp_off)];
    xo_1 = clatoms_xold[ktemp1];
    yo_1 = clatoms_yold[ktemp1];
    zo_1 = clatoms_zold[ktemp1];

    ktemp2 = grp_bond_con_j2_33[(igrp+igrp_off)];
    xo_2 = clatoms_xold[ktemp2];
    yo_2 = clatoms_yold[ktemp2];
    zo_2 = clatoms_zold[ktemp2];

    ktemp3 = grp_bond_con_j3_33[(igrp+igrp_off)];
    xo_3 = clatoms_xold[ktemp3];
    yo_3 = clatoms_yold[ktemp3];
    zo_3 = clatoms_zold[ktemp3];

    dx[1][igrp]= xo_1 - xo_2;
    dx[2][igrp]= xo_1 - xo_3;
    dx[3][igrp]= xo_2 - xo_3;
    dy[1][igrp]= yo_1 - yo_2;
    dy[2][igrp]= yo_1 - yo_3;
    dy[3][igrp]= yo_2 - yo_3;
    dz[1][igrp]= zo_1 - zo_2;
    dz[2][igrp]= zo_1 - zo_3;
    dz[3][igrp]= zo_2 - zo_3;

 }/*endfor*/

/* compute the reduced mass matrix elements */

 for(igrp=1;igrp <= ngrp; igrp++) {
    rms1 = rm1[igrp];
    rms2 = rm2[igrp];
    rms3 = rm3[igrp];
    rmm11[igrp] = -(rms1+rms2); rmm12[igrp] = -rms1; rmm13[igrp] = rms2;
    rmm21[igrp] = -rms1; rmm22[igrp] = -(rms1+rms3); rmm23[igrp] = -rms3;
    rmm31[igrp] = rms2; rmm32[igrp] = -rms3; rmm33[igrp] = -(rms2+rms3);
 }/*endfor*/


/* Get initial guess for lambda */
 for(igrp=1;igrp <= ngrp; igrp++) {
    amat[1][1] = 2.0*rmm11[igrp]* (dxt[1][igrp]*dx[1][igrp] 
               + dyt[1][igrp]*dy[1][igrp] + dzt[1][igrp]*dz[1][igrp]);
    amat[1][2] = 2.0*rmm12[igrp]* (dxt[1][igrp]*dx[2][igrp] 
               + dyt[1][igrp]*dy[2][igrp] + dzt[1][igrp]*dz[2][igrp]);
    amat[1][3] = 2.0*rmm13[igrp]* (dxt[1][igrp]*dx[3][igrp] 
               + dyt[1][igrp]*dy[3][igrp] + dzt[1][igrp]*dz[3][igrp]);

    amat[2][1] = 2.0*rmm21[igrp]* (dxt[2][igrp]*dx[1][igrp] 
               + dyt[2][igrp]*dy[1][igrp] + dzt[2][igrp]*dz[1][igrp]);
    amat[2][2] = 2.0*rmm22[igrp]* (dxt[2][igrp]*dx[2][igrp] 
               + dyt[2][igrp]*dy[2][igrp] + dzt[2][igrp]*dz[2][igrp]);
    amat[2][3] = 2.0*rmm23[igrp]* (dxt[2][igrp]*dx[3][igrp] 
               + dyt[2][igrp]*dy[3][igrp] + dzt[2][igrp]*dz[3][igrp]);

    amat[3][1] = 2.0*rmm31[igrp]* (dxt[3][igrp]*dx[1][igrp] 
               + dyt[3][igrp]*dy[1][igrp] + dzt[3][igrp]*dz[1][igrp]);
    amat[3][2] = 2.0*rmm32[igrp]* (dxt[3][igrp]*dx[2][igrp] 
               + dyt[3][igrp]*dy[2][igrp] + dzt[3][igrp]*dz[2][igrp]);
    amat[3][3] = 2.0*rmm33[igrp]* (dxt[3][igrp]*dx[3][igrp] 
               + dyt[3][igrp]*dy[3][igrp] + dzt[3][igrp]*dz[3][igrp]);


    avec[1][igrp]= dij[1][igrp]*dij[1][igrp] - (dxt[1][igrp]*dxt[1][igrp]
                  +dyt[1][igrp]*dyt[1][igrp] + dzt[1][igrp]*dzt[1][igrp]);

    avec[2][igrp]= dij[2][igrp]*dij[2][igrp] - (dxt[2][igrp]*dxt[2][igrp]
                  +dyt[2][igrp]*dyt[2][igrp] + dzt[2][igrp]*dzt[2][igrp]);

    avec[3][igrp]= dij[3][igrp]*dij[3][igrp] - (dxt[3][igrp]*dxt[3][igrp]
                  +dyt[3][igrp]*dyt[3][igrp] + dzt[3][igrp]*dzt[3][igrp]);
        

/* Get inverse of matrix */

  det = (amat[1][1] * (amat[2][2] * amat[3][3] - amat[3][2] * amat[2][3]) + 
	 amat[2][1] * (amat[3][2] * amat[1][3] - amat[1][2] * amat[3][3]) + 
	 amat[3][1] * (amat[1][2] * amat[2][3] - amat[2][2] * amat[1][3]));
  rdet_a = 1.0/det;
  ainv[1][1] = (amat[2][2] * amat[3][3] - amat[3][2] * amat[2][3])*rdet_a;
  ainv[2][2] = (amat[1][1] * amat[3][3] - amat[3][1] * amat[1][3])*rdet_a;
  ainv[3][3] = (amat[1][1] * amat[2][2] - amat[2][1] * amat[1][2])*rdet_a;
  ainv[2][1] = (amat[3][1] * amat[2][3] - amat[2][1] * amat[3][3])*rdet_a;
  ainv[1][2] = (amat[1][3] * amat[3][2] - amat[1][2] * amat[3][3])*rdet_a;
  ainv[3][1] = (amat[2][1] * amat[3][2] - amat[3][1] * amat[2][2])*rdet_a;
  ainv[1][3] = (amat[1][2] * amat[2][3] - amat[1][3] * amat[2][2])*rdet_a;
  ainv[3][2] = (amat[3][1] * amat[1][2] - amat[3][2] * amat[1][1])*rdet_a;
  ainv[2][3] = (amat[1][3] * amat[2][1] - amat[2][3] * amat[1][1])*rdet_a;

/* initial guess for multipliers */

  xlam[1][igrp] = ainv[1][1]*avec[1][igrp] + ainv[1][2]*avec[2][igrp]
                 +ainv[1][3]*avec[3][igrp];
        
  xlam[2][igrp] = ainv[2][1]*avec[1][igrp] + ainv[2][2]*avec[2][igrp]
                 +ainv[2][3]*avec[3][igrp];
        
  xlam[3][igrp] = ainv[3][1]*avec[1][igrp] + ainv[3][2]*avec[2][igrp]
                 +ainv[3][3]*avec[3][igrp];
        
} /*end for*/

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

#ifndef NO_PRAGMA
#pragma IVDEP
#endif
   for(igrp=1;igrp <= ngrp; igrp++) {
/* Set up guess of difference vectors */

    dxn[1] = 2.0*dxt[1][igrp]; 
    dyn[1] = 2.0*dyt[1][igrp];
    dzn[1] = 2.0*dzt[1][igrp];

    dxn[2] = 2.0*dxt[2][igrp]; 
    dyn[2] = 2.0*dyt[2][igrp];
    dzn[2] = 2.0*dzt[2][igrp];

    dxn[3] = 2.0*dxt[3][igrp]; 
    dyn[3] = 2.0*dyt[3][igrp];
    dzn[3] = 2.0*dzt[3][igrp];

    dxn[1] += (rmm11[igrp]*xlam[1][igrp]*dx[1][igrp]+ 
               rmm12[igrp]*xlam[2][igrp]*dx[2][igrp]+
               rmm13[igrp]*xlam[3][igrp]*dx[3][igrp]);

    dxn[2] += (rmm21[igrp]*xlam[1][igrp]*dx[1][igrp]+ 
               rmm22[igrp]*xlam[2][igrp]*dx[2][igrp]+
               rmm23[igrp]*xlam[3][igrp]*dx[3][igrp]);

    dxn[3] += (rmm31[igrp]*xlam[1][igrp]*dx[1][igrp]+ 
               rmm32[igrp]*xlam[2][igrp]*dx[2][igrp]+
               rmm33[igrp]*xlam[3][igrp]*dx[3][igrp]);

    dyn[1] += (rmm11[igrp]*xlam[1][igrp]*dy[1][igrp]+ 
               rmm12[igrp]*xlam[2][igrp]*dy[2][igrp]+
               rmm13[igrp]*xlam[3][igrp]*dy[3][igrp]);

    dyn[2] += (rmm21[igrp]*xlam[1][igrp]*dy[1][igrp]+ 
               rmm22[igrp]*xlam[2][igrp]*dy[2][igrp]+
               rmm23[igrp]*xlam[3][igrp]*dy[3][igrp]);

    dyn[3] += (rmm31[igrp]*xlam[1][igrp]*dy[1][igrp]+ 
               rmm32[igrp]*xlam[2][igrp]*dy[2][igrp]+
               rmm33[igrp]*xlam[3][igrp]*dy[3][igrp]);

    dzn[1] += (rmm11[igrp]*xlam[1][igrp]*dz[1][igrp]+ 
               rmm12[igrp]*xlam[2][igrp]*dz[2][igrp]+
               rmm13[igrp]*xlam[3][igrp]*dz[3][igrp]);

    dzn[2] += (rmm21[igrp]*xlam[1][igrp]*dz[1][igrp]+ 
               rmm22[igrp]*xlam[2][igrp]*dz[2][igrp]+
               rmm23[igrp]*xlam[3][igrp]*dz[3][igrp]);

    dzn[3] += (rmm31[igrp]*xlam[1][igrp]*dz[1][igrp]+ 
               rmm32[igrp]*xlam[2][igrp]*dz[2][igrp]+
               rmm33[igrp]*xlam[3][igrp]*dz[3][igrp]);

/* Construct A-matrix */

      amat[1][1] = rmm11[igrp]* (dxn[1]*dx[1][igrp] 
               + dyn[1]*dy[1][igrp] + dzn[1]*dz[1][igrp]);
      amat[1][2] = rmm12[igrp]* (dxn[1]*dx[2][igrp] 
               + dyn[1]*dy[2][igrp] + dzn[1]*dz[2][igrp]);
      amat[1][3] = rmm13[igrp]* (dxn[1]*dx[3][igrp] 
               + dyn[1]*dy[3][igrp] + dzn[1]*dz[3][igrp]);

      amat[2][1] = rmm21[igrp]* (dxn[2]*dx[1][igrp] 
               + dyn[2]*dy[1][igrp] + dzn[2]*dz[1][igrp]);
      amat[2][2] = rmm22[igrp]* (dxn[2]*dx[2][igrp] 
               + dyn[2]*dy[2][igrp] + dzn[2]*dz[2][igrp]);
      amat[2][3] = rmm23[igrp]* (dxn[2]*dx[3][igrp] 
               + dyn[2]*dy[3][igrp] + dzn[2]*dz[3][igrp]);

      amat[3][1] = rmm31[igrp]* (dxn[3]*dx[1][igrp] 
               + dyn[3]*dy[1][igrp] + dzn[3]*dz[1][igrp]);
      amat[3][2] = rmm32[igrp]* (dxn[3]*dx[2][igrp] 
               + dyn[3]*dy[2][igrp] + dzn[3]*dz[2][igrp]);
      amat[3][3] = rmm33[igrp]* (dxn[3]*dx[3][igrp] 
               + dyn[3]*dy[3][igrp] + dzn[3]*dz[3][igrp]);

/* Get inverse of matrix */

    det = (amat[1][1] * (amat[2][2] * amat[3][3] - amat[3][2] * amat[2][3]) + 
	   amat[2][1] * (amat[3][2] * amat[1][3] - amat[1][2] * amat[3][3]) + 
	   amat[3][1] * (amat[1][2] * amat[2][3] - amat[2][2] * amat[1][3]));
    rdet_a = 1.0/det;
    ainv[1][1] = (amat[2][2] * amat[3][3] - amat[3][2] * amat[2][3])*rdet_a;
    ainv[2][2] = (amat[1][1] * amat[3][3] - amat[3][1] * amat[1][3])*rdet_a;
    ainv[3][3] = (amat[1][1] * amat[2][2] - amat[2][1] * amat[1][2])*rdet_a;
    ainv[2][1] = (amat[3][1] * amat[2][3] - amat[2][1] * amat[3][3])*rdet_a;
    ainv[1][2] = (amat[1][3] * amat[3][2] - amat[1][2] * amat[3][3])*rdet_a;
    ainv[3][1] = (amat[2][1] * amat[3][2] - amat[3][1] * amat[2][2])*rdet_a;
    ainv[1][3] = (amat[1][2] * amat[2][3] - amat[1][3] * amat[2][2])*rdet_a;
    ainv[3][2] = (amat[3][1] * amat[1][2] - amat[3][2] * amat[1][1])*rdet_a;
    ainv[2][3] = (amat[1][3] * amat[2][1] - amat[2][3] * amat[1][1])*rdet_a;

  xl0[1] = xlam[1][igrp];
  xl0[2] = xlam[2][igrp];
  xl0[3] = xlam[3][igrp];

 xlam[1][igrp] = ainv[1][1]*avec[1][igrp]  
               + ainv[1][2]*avec[2][igrp] + ainv[1][3]*avec[3][igrp];
 xlam[2][igrp] = ainv[2][1]*avec[1][igrp]  
               + ainv[2][2]*avec[2][igrp] + ainv[2][3]*avec[3][igrp];
 xlam[3][igrp] = ainv[3][1]*avec[1][igrp]  
               + ainv[3][2]*avec[2][igrp] + ainv[3][3]*avec[3][igrp];


 dxl[1][igrp] = fabs(xlam[1][igrp]-xl0[1]);
 dxl[2][igrp] = fabs(xlam[2][igrp]-xl0[2]);
 dxl[3][igrp] = fabs(xlam[3][igrp]-xl0[3]);
}  /* end for loop over groups */
     dlmax1= dxl[1][1];
     dlmax2= dxl[2][1];
     dlmax3= dxl[3][1];

  for(igrp=2;igrp<=ngrp;igrp++){
     dlmax1= dlmax1 > dxl[1][igrp] ? dlmax1: dxl[1][igrp];
     dlmax2= dlmax2 > dxl[2][igrp] ? dlmax2: dxl[2][igrp];
     dlmax3= dlmax3 > dxl[3][igrp] ? dlmax3: dxl[3][igrp];
  }

   dlmax = (dlmax1 > dlmax2 ? dlmax1:dlmax2);
   dlmax = (dlmax3 > dlmax ? dlmax3:dlmax);
 } while(dlmax > grp_bond_con->tol);
 *aiter += (double) iter;

/* position update */

#ifndef NO_PRAGMA
#pragma IVDEP
#endif
 for(igrp=1;igrp <= ngrp; igrp++) {
    double xlam1,xlam2,xlam3;
    double rms1,rms2,rms3;
    double dx1,dx2,dx3;
    double dy1,dy2,dy3;
    double dz1,dz2,dz3;
    int ktemp1,ktemp2,ktemp3;

     rms1=rm1[igrp];
     rms2=rm2[igrp];
     rms3=rm3[igrp];

     xlam1=xlam[1][igrp];
     xlam2=xlam[2][igrp];
     xlam3=xlam[3][igrp];
    
     dx1=dx[1][igrp]; dx2=dx[2][igrp]; dx3=dx[3][igrp];
     dy1=dy[1][igrp]; dy2=dy[2][igrp]; dy3=dy[3][igrp];
     dz1=dz[1][igrp]; dz2=dz[2][igrp]; dz3=dz[3][igrp];

     ktemp1 = grp_bond_con_j1_33[(igrp+igrp_off)];
     ktemp2 = grp_bond_con_j2_33[(igrp+igrp_off)]; 
     ktemp3 = grp_bond_con_j3_33[(igrp+igrp_off)];

    clatoms_x[ktemp1] -= (xlam1*dx1 + xlam2*dx2)*rms1;
    clatoms_y[ktemp1] -= (xlam1*dy1 + xlam2*dy2)*rms1;
    clatoms_z[ktemp1] -= (xlam1*dz1 + xlam2*dz2)*rms1;

    clatoms_vx[ktemp1] -= (xlam1*dx1 + xlam2*dx2)*rms1/dt;
    clatoms_vy[ktemp1] -= (xlam1*dy1 + xlam2*dy2)*rms1/dt;
    clatoms_vz[ktemp1] -= (xlam1*dz1 + xlam2*dz2)*rms1/dt;

    clatoms_x[ktemp2] -= (-xlam1*dx1 + xlam3*dx3)*rms2;
    clatoms_y[ktemp2] -= (-xlam1*dy1 + xlam3*dy3)*rms2;
    clatoms_z[ktemp2] -= (-xlam1*dz1 + xlam3*dz3)*rms2;

    clatoms_vx[ktemp2] -= (-xlam1*dx1 + xlam3*dx3)*rms2/dt;
    clatoms_vy[ktemp2] -= (-xlam1*dy1 + xlam3*dy3)*rms2/dt;
    clatoms_vz[ktemp2] -= (-xlam1*dz1 + xlam3*dz3)*rms2/dt;

    clatoms_x[ktemp3] -= (-xlam2*dx2 - xlam3*dx3)*rms3;
    clatoms_y[ktemp3] -= (-xlam2*dy2 - xlam3*dy3)*rms3;
    clatoms_z[ktemp3] -= (-xlam2*dz2 - xlam3*dz3)*rms3;

    clatoms_vx[ktemp3] -= (-xlam2*dx2 - xlam3*dx3)*rms3/dt;
    clatoms_vy[ktemp3] -= (-xlam2*dy2 - xlam3*dy3)*rms3/dt;
    clatoms_vz[ktemp3] -= (-xlam2*dz2 - xlam3*dz3)*rms3/dt;
  }

/* Pressure tensor update */

 for(igrp=1;igrp <= ngrp; igrp++) {
    double xlam1,xlam2,xlam3;
    double rms1,rms2,rms3;
    double dx1,dx2,dx3;
    double dy1,dy2,dy3;
    double dz1,dz2,dz3;

     rms1=rm1[igrp];
     rms2=rm2[igrp];
     rms3=rm3[igrp];

     xlam1=xlam[1][igrp];
     xlam2=xlam[2][igrp];
     xlam3=xlam[3][igrp];
    
     dx1=dx[1][igrp]; dx2=dx[2][igrp]; dx3=dx[3][igrp];
     dy1=dy[1][igrp]; dy2=dy[2][igrp]; dy3=dy[3][igrp];
     dz1=dz[1][igrp]; dz2=dz[2][igrp]; dz3=dz[3][igrp];

    p11[igrp] = xlam1*dx1*dx1 + xlam2*dx2*dx2 + xlam3*dx3*dx3;
    p22[igrp] = xlam1*dy1*dy1 + xlam2*dy2*dy2 + xlam3*dy3*dy3;
    p33[igrp] = xlam1*dz1*dz1 + xlam2*dz2*dz2 + xlam3*dz3*dz3;
    p12[igrp] = xlam1*dx1*dy1 + xlam2*dx2*dy2 + xlam3*dx3*dy3;
    p13[igrp] = xlam1*dx1*dz1 + xlam2*dx2*dz2 + xlam3*dx3*dz3;
    p23[igrp] = xlam1*dy1*dz1 + xlam2*dy2*dz2 + xlam3*dy3*dz3;
 }/*endfor */

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
    grp_bond_con_al_33[1][(igrp+igrp_off)] = xlam[1][igrp];
    grp_bond_con_al_33[2][(igrp+igrp_off)] = xlam[2][igrp];
    grp_bond_con_al_33[3][(igrp+igrp_off)] = xlam[3][igrp];
 } /* end for */

/* free all the memory allocated locally */
     free_dmatrix(dx,1,3,1,ngrp);
     free_dmatrix(dy,1,3,1,ngrp);
     free_dmatrix(dz,1,3,1,ngrp);

     free_dmatrix(dxt,1,3,1,ngrp);
     free_dmatrix(dyt,1,3,1,ngrp);
     free_dmatrix(dzt,1,3,1,ngrp);

     free_dmatrix(xlam,1,3,1,ngrp);
     free_dmatrix(avec,1,3,1,ngrp);
     free_dmatrix(dxl,1,3,1,ngrp);
     free_dmatrix(dij,1,3,1,ngrp);

     free_dvector(rm1,1,ngrp);
     free_dvector(rm2,1,ngrp);
     free_dvector(rm3,1,ngrp);

     free_dvector(rmm11,1,ngrp);
     free_dvector(rmm12,1,ngrp);
     free_dvector(rmm13,1,ngrp);
     free_dvector(rmm21,1,ngrp);
     free_dvector(rmm22,1,ngrp);
     free_dvector(rmm23,1,ngrp);
     free_dvector(rmm31,1,ngrp);
     free_dvector(rmm32,1,ngrp);
     free_dvector(rmm33,1,ngrp);

     free_dvector(p11,1,ngrp);
     free_dvector(p12,1,ngrp);
     free_dvector(p13,1,ngrp);
     free_dvector(p22,1,ngrp);
     free_dvector(p23,1,ngrp);
     free_dvector(p33,1,ngrp);

/*=======================================================================*/
} /* end routine */
/*=======================================================================*/

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void rattle_33(GRP_BOND_CON *grp_bond_con,
              CLATOMS_INFO *clatoms_info,CLATOMS_POS *clatoms_pos,
              PTENS *ptens,double dt,
              CLASS_COMM_FORC_PKG *class_comm_forc_pkg)

/*==========================================================================*/
/*        Begin Routine                                                     */
{/*Begin Routine*/
/*=======================================================================*/
/*         Local Variable declarations                                   */
  
  double avec[NCON_33+1];
  double amat[NCON_33+1][NCON_33+1],ainv[NCON_33+1][NCON_33+1];
/*  double x[NAT_33+1],y[NAT_33+1],z[NAT_33+1];
  double vx[NAT_33+1],vy[NAT_33+1],vz[NAT_33+1]; */

/*  double dx[NCON_33+1],dy[NCON_33+1],dz[NCON_33+1];
  double dvx[NCON_33+1],dvy[NCON_33+1],dvz[NCON_33+1]; */
  double rms1,rms2,rms3;


  double det,rdet_a;

  int i,j,k;
  int ktemp,ktemp1,ktemp2,ktemp3;
  int igrp,jtyp;

   double *p11,*p12,*p13,*p23,*p22,*p33;
   double *rm1,*rm2,*rm3;
   double **dx,**dy,**dz;
   double **dvx,**dvy,**dvz;
   double ***rmm;
   double **xlam;

/* Local pointers */
  double *clatoms_mass         = clatoms_info->mass;
  double *clatoms_x            = clatoms_pos->x;
  double *clatoms_y            = clatoms_pos->y;
  double *clatoms_z            = clatoms_pos->z;
  double *clatoms_vx           = clatoms_pos->vx;
  double *clatoms_vy           = clatoms_pos->vy;
  double *clatoms_vz           = clatoms_pos->vz;
  int *grp_bond_con_j1_33      = grp_bond_con->j1_33;
  int *grp_bond_con_j2_33      = grp_bond_con->j2_33;
  int *grp_bond_con_j3_33      = grp_bond_con->j3_33;
  int *grp_bond_con_jtyp_33    = grp_bond_con->jtyp_33;
  double *ptens_pvten_inc      = ptens->pvten_inc;
  double pnorm;
  double x_1,x_2,x_3,y_1,y_2,y_3,z_1,z_2,z_3;
  double vx_1,vx_2,vx_3,vy_1,vy_2,vy_3,vz_1,vz_2,vz_3;

  int ngrp,irem,igrp_off;
  int ngrp_tot = grp_bond_con->num_33;
  int np_forc                  = class_comm_forc_pkg->num_proc;
  int myid_forc                = class_comm_forc_pkg->myid;

/*=======================================================================*/

  ngrp = ngrp_tot;
  igrp_off = 0;

/*=======================================================================*/

  p11= dvector(1,ngrp);
  p12= dvector(1,ngrp);
  p13= dvector(1,ngrp);
  p22= dvector(1,ngrp);
  p23= dvector(1,ngrp);
  p33= dvector(1,ngrp);
  rm1= dvector(1,ngrp);
  rm2= dvector(1,ngrp);
  rm3= dvector(1,ngrp);

  dx= dmatrix(1,3,1,ngrp);
  dy= dmatrix(1,3,1,ngrp);
  dz= dmatrix(1,3,1,ngrp);
 dvx= dmatrix(1,3,1,ngrp);
 dvy= dmatrix(1,3,1,ngrp);
 dvz= dmatrix(1,3,1,ngrp);
 xlam= dmatrix(1,3,1,ngrp);
 rmm= d3tensor(1,3,1,3,1,ngrp);
  
/*=======================================================================*/

 pnorm = 2.0/dt;

/* Set zeta matrix and reciprocal mass matrix */

 for(igrp=1;igrp <= ngrp; igrp++) {

  rm1[igrp] = 1.0/clatoms_mass[grp_bond_con_j1_33[(igrp+igrp_off)]];
  rm2[igrp] = 1.0/clatoms_mass[grp_bond_con_j2_33[(igrp+igrp_off)]];
  rm3[igrp] = 1.0/clatoms_mass[grp_bond_con_j3_33[(igrp+igrp_off)]];

 }/*endfor*/

 for(igrp=1;igrp <= ngrp; igrp++) {

  ktemp1 = grp_bond_con_j1_33[(igrp+igrp_off)];
  x_1 = clatoms_x[ktemp1];
  y_1 = clatoms_y[ktemp1];
  z_1 = clatoms_z[ktemp1];
  ktemp2 = grp_bond_con_j2_33[(igrp+igrp_off)];
  x_2 = clatoms_x[ktemp2];
  y_2 = clatoms_y[ktemp2];
  z_2 = clatoms_z[ktemp2];
  ktemp3 = grp_bond_con_j3_33[(igrp+igrp_off)];
  x_3 = clatoms_x[ktemp3];
  y_3 = clatoms_y[ktemp3];
  z_3 = clatoms_z[ktemp3];

/* Compute difference vectors */

  dx[1][igrp] = x_1-x_2;
  dx[2][igrp] = x_1-x_3;
  dx[3][igrp] = x_2-x_3;

  dy[1][igrp] = y_1-y_2;
  dy[2][igrp] = y_1-y_3;
  dy[3][igrp] = y_2-y_3;

  dz[1][igrp] = z_1-z_2;
  dz[2][igrp] = z_1-z_3;
  dz[3][igrp] = z_2-z_3;

}/*endfor*/

 for(igrp=1;igrp <= ngrp; igrp++) {

  ktemp1 = grp_bond_con_j1_33[(igrp+igrp_off)];
  vx_1 = clatoms_vx[ktemp1];
  vy_1 = clatoms_vy[ktemp1];
  vz_1 = clatoms_vz[ktemp1];
  ktemp2 = grp_bond_con_j2_33[(igrp+igrp_off)];
  vx_2 = clatoms_vx[ktemp2];
  vy_2 = clatoms_vy[ktemp2];
  vz_2 = clatoms_vz[ktemp2];
  ktemp3 = grp_bond_con_j3_33[(igrp+igrp_off)];
  vx_3 = clatoms_vx[ktemp3];
  vy_3 = clatoms_vy[ktemp3];
  vz_3 = clatoms_vz[ktemp3];

/* Compute difference vectors */

   dvx[1][igrp] = vx_1-vx_2;
   dvx[2][igrp] = vx_1-vx_3;
   dvx[3][igrp] = vx_2-vx_3;

   dvy[1][igrp] = vy_1-vy_2;
   dvy[2][igrp] = vy_1-vy_3;
   dvy[3][igrp] = vy_2-vy_3;

   dvz[1][igrp] = vz_1-vz_2;
   dvz[2][igrp] = vz_1-vz_3;
   dvz[3][igrp] = vz_2-vz_3;

 } /*end for*/

 for(igrp=1;igrp <= ngrp; igrp++) {
  rms1 = rm1[igrp];
  rms2 = rm2[igrp];
  rms3 = rm3[igrp];

  rmm[1][1][igrp] = -(rms1+rms2); rmm[1][2][igrp] = -rms1; 
                      rmm[1][3][igrp] = rms2;
  rmm[2][1][igrp] = -rms1; rmm[2][2][igrp] = -(rms1+rms3); 
                     rmm[2][3][igrp] = -rms3;
  rmm[3][1][igrp] = rms2; rmm[3][2][igrp] = -rms3; rmm[3][3][igrp] = 
                    -(rms2+rms3);
 }/* end for */


/* Get lambda */
#ifndef NO_PRAGMA
#pragma IVDEP
#endif
 for(igrp=1;igrp <= ngrp; igrp++) {
   amat[1][1] =-rmm[1][1][igrp]*
               (dx[1][igrp]*dx[1][igrp] + dy[1][igrp]*dy[1][igrp] + dz[1][igrp]*dz[1][igrp]);
   amat[1][2] =-rmm[1][2][igrp]*
               (dx[1][igrp]*dx[2][igrp] + dy[1][igrp]*dy[2][igrp] + dz[1][igrp]*dz[2][igrp]);
   amat[1][3] =-rmm[1][3][igrp]*
               (dx[1][igrp]*dx[3][igrp] + dy[1][igrp]*dy[3][igrp] + dz[1][igrp]*dz[3][igrp]);

   amat[2][1] =-rmm[2][1][igrp]*
               (dx[2][igrp]*dx[1][igrp] + dy[2][igrp]*dy[1][igrp] + dz[2][igrp]*dz[1][igrp]);
   amat[2][2] =-rmm[2][2][igrp]*
               (dx[2][igrp]*dx[2][igrp] + dy[2][igrp]*dy[2][igrp] + dz[2][igrp]*dz[2][igrp]);
   amat[2][3] =-rmm[2][3][igrp]*
               (dx[2][igrp]*dx[3][igrp] + dy[2][igrp]*dy[3][igrp] + dz[2][igrp]*dz[3][igrp]);

   amat[3][1] =-rmm[3][1][igrp]*
        (dx[3][igrp]*dx[1][igrp] + dy[3][igrp]*dy[1][igrp] + dz[3][igrp]*dz[1][igrp]);
   amat[3][2] =-rmm[3][2][igrp]*
       (dx[3][igrp]*dx[2][igrp] + dy[3][igrp]*dy[2][igrp] + dz[3][igrp]*dz[2][igrp]);
   amat[3][3] =-rmm[3][3][igrp]*
       (dx[3][igrp]*dx[3][igrp] + dy[3][igrp]*dy[3][igrp] + dz[3][igrp]*dz[3][igrp]);

  avec[1] = dvx[1][igrp]*dx[1][igrp] + dvy[1][igrp]*dy[1][igrp] + dvz[1][igrp]*dz[1][igrp];
  avec[2] = dvx[2][igrp]*dx[2][igrp] + dvy[2][igrp]*dy[2][igrp] + dvz[2][igrp]*dz[2][igrp];
  avec[3] = dvx[3][igrp]*dx[3][igrp] + dvy[3][igrp]*dy[3][igrp] + dvz[3][igrp]*dz[3][igrp];


/* Get inverse of matrix */

    det = (amat[1][1] * (amat[2][2] * amat[3][3] - amat[3][2] * amat[2][3]) + 
	   amat[2][1] * (amat[3][2] * amat[1][3] - amat[1][2] * amat[3][3]) + 
	   amat[3][1] * (amat[1][2] * amat[2][3] - amat[2][2] * amat[1][3]));
    rdet_a = 1.0/det;
    ainv[1][1] = (amat[2][2] * amat[3][3] - amat[3][2] * amat[2][3])*rdet_a;
    ainv[2][2] = (amat[1][1] * amat[3][3] - amat[3][1] * amat[1][3])*rdet_a;
    ainv[3][3] = (amat[1][1] * amat[2][2] - amat[2][1] * amat[1][2])*rdet_a;
    ainv[2][1] = (amat[3][1] * amat[2][3] - amat[2][1] * amat[3][3])*rdet_a;
    ainv[1][2] = (amat[1][3] * amat[3][2] - amat[1][2] * amat[3][3])*rdet_a;
    ainv[3][1] = (amat[2][1] * amat[3][2] - amat[3][1] * amat[2][2])*rdet_a;
    ainv[1][3] = (amat[1][2] * amat[2][3] - amat[1][3] * amat[2][2])*rdet_a;
    ainv[3][2] = (amat[3][1] * amat[1][2] - amat[3][2] * amat[1][1])*rdet_a;
    ainv[2][3] = (amat[1][3] * amat[2][1] - amat[2][3] * amat[1][1])*rdet_a;

/* Get initial guess for multipliers */

 xlam[1][igrp] = ainv[1][1]*avec[1] + ainv[1][2]*avec[2] + ainv[1][3]*avec[3];
 xlam[2][igrp] = ainv[2][1]*avec[1] + ainv[2][2]*avec[2] + ainv[2][3]*avec[3];
 xlam[3][igrp] = ainv[3][1]*avec[1] + ainv[3][2]*avec[2] + ainv[3][3]*avec[3];
}/* end for */

/* update velocities */
#ifndef NO_PRAGMA
#pragma IVDEP
#endif
 for(igrp=1;igrp <= ngrp; igrp++) {
  double xlam1,xlam2,xlam3;
  double dx1,dx2,dx3;
  double dy1,dy2,dy3;
  double dz1,dz2,dz3;
  int ktemp1,ktemp2,ktemp3;

    dx1= dx[1][igrp]; dx2= dx[2][igrp]; dx3= dx[3][igrp];
    dy1= dy[1][igrp]; dy2= dy[2][igrp]; dy3= dy[3][igrp];
    dz1= dz[1][igrp]; dz2= dz[2][igrp]; dz3= dz[3][igrp];

    xlam1= xlam[1][igrp];
    xlam2= xlam[2][igrp];
    xlam3= xlam[3][igrp];

    ktemp1=grp_bond_con_j1_33[(igrp+igrp_off)];
    ktemp2=grp_bond_con_j2_33[(igrp+igrp_off)]; 
    ktemp3=grp_bond_con_j3_33[(igrp+igrp_off)];

    clatoms_vx[ktemp1] -= (xlam1*dx1 + xlam2*dx2)*rm1[igrp];
    clatoms_vy[ktemp1] -= (xlam1*dy1 + xlam2*dy2)*rm1[igrp];
    clatoms_vz[ktemp1] -= (xlam1*dz1 + xlam2*dz2)*rm1[igrp];

    clatoms_vx[ktemp2] += (xlam1*dx1 - xlam3*dx3)*rm2[igrp];
    clatoms_vy[ktemp2] += (xlam1*dy1 - xlam3*dy3)*rm2[igrp];
    clatoms_vz[ktemp2] += (xlam1*dz1 - xlam3*dz3)*rm2[igrp];

    clatoms_vx[ktemp3] += (xlam2*dx2 + xlam3*dx3)*rm3[igrp];
    clatoms_vy[ktemp3] += (xlam2*dy2 + xlam3*dy3)*rm3[igrp];
    clatoms_vz[ktemp3] += (xlam2*dz2 + xlam3*dz3)*rm3[igrp];

/* Pressure tensor update */

    p11[igrp] = xlam1*dx1*dx1 + xlam2*dx2*dx2 + xlam3*dx3*dx3;
    p22[igrp] = xlam1*dy1*dy1 + xlam2*dy2*dy2 + xlam3*dy3*dy3;
    p33[igrp] = xlam1*dz1*dz1 + xlam2*dz2*dz2 + xlam3*dz3*dz3;
    p12[igrp] = xlam1*dx1*dy1 + xlam2*dx2*dy2 + xlam3*dx3*dy3;
    p13[igrp] = xlam1*dx1*dz1 + xlam2*dx2*dz2 + xlam3*dx3*dz3;
    p23[igrp] = xlam1*dy1*dz1 + xlam2*dy2*dz2 + xlam3*dy3*dz3;
}/* end for */

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
 } /* end for igrp */

 /* free locally assigned memory */ 

   free_dvector(p11,1,ngrp);
   free_dvector(p12,1,ngrp);
   free_dvector(p13,1,ngrp);
   free_dvector(p22,1,ngrp);
   free_dvector(p23,1,ngrp);
   free_dvector(p33,1,ngrp);

   free_dvector(rm1,1,ngrp);
   free_dvector(rm2,1,ngrp);
   free_dvector(rm3,1,ngrp);

   free_dmatrix(dx,1,3,1,ngrp);
   free_dmatrix(dy,1,3,1,ngrp);
   free_dmatrix(dz,1,3,1,ngrp);

   free_dmatrix(dvx,1,3,1,ngrp);
   free_dmatrix(dvy,1,3,1,ngrp);
   free_dmatrix(dvz,1,3,1,ngrp);

   free_dmatrix(xlam,1,3,1,ngrp);

   free_d3tensor(rmm,1,3,1,3,1,ngrp);

/*=======================================================================*/
} /* end routine */
/*=======================================================================*/

