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
#include "../proto_defs/proto_communicate_wrappers.h"

#define NCON_43 3
#define NAT_43 4


/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
void shake_43_rolli(GRP_BOND_CON *grp_bond_con,
              CLATOMS_INFO *clatoms_info,CLATOMS_POS *clatoms_pos,
              PTENS *ptens,double dt,double *aiter,
              BARO *baro,int ifirst,
              CLASS_COMM_FORC_PKG *class_comm_forc_pkg)


/*==========================================================================*/
/*        Begin Routine                                                     */
{/*Begin Routine*/
/*=======================================================================*/
/*         Local Variable declarations                                   */

#include "../typ_defs/typ_mask.h"

 double xl0[NCON_43+1];
 double dxn[NCON_43+1],dyn[NCON_43+1],dzn[NCON_43+1];
 double amat[NCON_43+1][NCON_43+1],ainv[NCON_43+1][NCON_43+1];
 double dlmax,dlmax1,dlmax2,dlmax3,det,rdet_a;
 double rms1,rms2,rms3,rms4;

 double dts;
 double ftemp; 
 int i,ktemp;
 int iter,igrp,*ind1,*ind2,*ind3,*ind4,jtyp;
  
/* Local arrays   */
   double **dx,**dy,**dz;
   double **dxo,**dyo,**dzo;
   double **dxt,**dyt,**dzt;
   double **xlam,**dxl;
   double **avec,**dij;
   double *rmm11,*rmm12, *rmm13, *rmm21, *rmm22, *rmm23;
   double *rmm31, *rmm32, *rmm33;
   double *rm1, *rm2, *rm3,*rm4;
   double **x,**y,**z,**xo,**yo,**zo;
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
  int *grp_bond_con_j1_43      = grp_bond_con->j1_43;
  int *grp_bond_con_j2_43      = grp_bond_con->j2_43;
  int *grp_bond_con_j3_43      = grp_bond_con->j3_43;
  int *grp_bond_con_j4_43      = grp_bond_con->j4_43;
  int *grp_bond_con_jtyp_43    = grp_bond_con->jtyp_43;
  double **grp_bond_con_eq_43  = grp_bond_con->eq_43;
  double **grp_bond_con_al_43  = grp_bond_con->al_43;
  double *ptens_pvten_inc      = ptens->pvten_inc;
  double *ptens_pvten_tmp      = ptens->pvten_tmp;
  double *ptens_pvten_tmp2     = ptens->pvten_tmp_res;
  double pnorm;
  double baro_roll_scv         = baro->roll_scv;

  int ngrp,irem,igrp_off;
  int ngrp_tot                 = grp_bond_con->num_43;
  int np_forc                  = class_comm_forc_pkg->num_proc;
  int myid_forc                = class_comm_forc_pkg->myid;
  MPI_Comm comm_forc           = class_comm_forc_pkg->comm;

/*=======================================================================*/
  ngrp = (ngrp_tot);
  igrp_off = 0;
/*=======================================================================*/
 
/* assign local arrays */ 
  if(ngrp > 0){
    dx= dmatrix(1,3,1,ngrp);
    dy= dmatrix(1,3,1,ngrp);
    dz= dmatrix(1,3,1,ngrp);

    dxo= dmatrix(1,3,1,ngrp);
    dyo= dmatrix(1,3,1,ngrp);
    dzo= dmatrix(1,3,1,ngrp);

     x= dmatrix(1,4,1,ngrp);
     y= dmatrix(1,4,1,ngrp);
     z= dmatrix(1,4,1,ngrp);
    xo= dmatrix(1,4,1,ngrp);
    yo= dmatrix(1,4,1,ngrp);
    zo= dmatrix(1,4,1,ngrp);

    rm1= dvector(1,ngrp); 
    rm2= dvector(1,ngrp); 
    rm3= dvector(1,ngrp); 
    rm4= dvector(1,ngrp); 

   dxt= dmatrix(1,3,1,ngrp);
   dyt= dmatrix(1,3,1,ngrp);
   dzt= dmatrix(1,3,1,ngrp);

  xlam= dmatrix(1,3,1,ngrp);
  avec= dmatrix(1,3,1,ngrp);
   dxl= dmatrix(1,3,1,ngrp);
   dij= dmatrix(1,3,1,ngrp); 
    
   rmm11= dvector(1,ngrp); rmm12= dvector(1,ngrp); rmm13= dvector(1,ngrp);
   rmm21= dvector(1,ngrp); rmm22= dvector(1,ngrp); rmm23= dvector(1,ngrp);
   rmm31= dvector(1,ngrp); rmm32= dvector(1,ngrp); rmm33= dvector(1,ngrp);

    p11= dvector(1,ngrp);
    p12= dvector(1,ngrp);
    p13= dvector(1,ngrp);
    p22= dvector(1,ngrp);
    p23= dvector(1,ngrp);
    p33= dvector(1,ngrp);
    ind1 = (int *) calloc((ngrp+1),sizeof(int));
    ind2 = (int *) calloc((ngrp+1),sizeof(int));
    ind3 = (int *) calloc((ngrp+1),sizeof(int));
    ind4 = (int *) calloc((ngrp+1),sizeof(int));

  }/*endif*/

/*=======================================================================*/

/* Set zeta matrix and reciprocal mass matrix */

 dts = dt*dt;
 pnorm = 2.0/dts;
 *aiter = 0.0;
  ptens_pvten_tmp[1] = 0.0;
  ptens_pvten_tmp[2] = 0.0;
  ptens_pvten_tmp[3] = 0.0;
  ptens_pvten_tmp[4] = 0.0;
  ptens_pvten_tmp[5] = 0.0;
  ptens_pvten_tmp[6] = 0.0;
  ptens_pvten_tmp[7] = 0.0;
  ptens_pvten_tmp[8] = 0.0;
  ptens_pvten_tmp[9] = 0.0;
 
 if(ifirst == 2){
   for(igrp=1;igrp <= ngrp ; igrp++) {
     grp_bond_con_al_43[1][(igrp+igrp_off)] = 0.0;
     grp_bond_con_al_43[2][(igrp+igrp_off)] = 0.0;
     grp_bond_con_al_43[3][(igrp+igrp_off)] = 0.0;
   }
 }

/* I assign positions and masses */

 for(igrp=1;igrp <= ngrp ; igrp++) {
    ind1[igrp] = grp_bond_con_j1_43[(igrp+igrp_off)];
    ind2[igrp] = grp_bond_con_j2_43[(igrp+igrp_off)];
    ind3[igrp] = grp_bond_con_j3_43[(igrp+igrp_off)];
    ind4[igrp] = grp_bond_con_j4_43[(igrp+igrp_off)];
  }
 for(igrp=1;igrp <= ngrp ; igrp++) {
    ktemp = ind1[igrp];
    x[1][igrp] = clatoms_x[ktemp];
    y[1][igrp] = clatoms_y[ktemp];
    z[1][igrp] = clatoms_z[ktemp];
  }

 for(igrp=1;igrp <= ngrp ; igrp++) {
    ktemp = ind2[igrp];
    x[2][igrp] = clatoms_x[ktemp];
    y[2][igrp] = clatoms_y[ktemp];
    z[2][igrp] = clatoms_z[ktemp];
  }

 for(igrp=1;igrp <= ngrp ; igrp++) {
    ktemp = ind3[igrp];
    x[3][igrp] = clatoms_x[ktemp];
    y[3][igrp] = clatoms_y[ktemp];
    z[3][igrp] = clatoms_z[ktemp];
  }

 for(igrp=1;igrp <= ngrp ; igrp++) {
    ktemp = ind4[igrp];
    x[4][igrp] = clatoms_x[ktemp];
    y[4][igrp] = clatoms_y[ktemp];
    z[4][igrp] = clatoms_z[ktemp];
  }


 for(igrp=1;igrp <= ngrp ; igrp++) {
    ktemp = ind1[igrp];
    xo[1][igrp] = clatoms_xold[ktemp];
    yo[1][igrp] = clatoms_yold[ktemp];
    zo[1][igrp] = clatoms_zold[ktemp];
    rm1[igrp] = 1.0/clatoms_mass[ktemp];
  }

 for(igrp=1;igrp <= ngrp ; igrp++) {
    ktemp = ind2[igrp];
    xo[2][igrp] = clatoms_xold[ktemp];
    yo[2][igrp] = clatoms_yold[ktemp];
    zo[2][igrp] = clatoms_zold[ktemp];
    rm2[igrp] = 1.0/clatoms_mass[ktemp];
  }

 for(igrp=1;igrp <= ngrp ; igrp++) {
    ktemp = ind3[igrp];
    xo[3][igrp] = clatoms_xold[ktemp];
    yo[3][igrp] = clatoms_yold[ktemp];
    zo[3][igrp] = clatoms_zold[ktemp];
    rm3[igrp] = 1.0/clatoms_mass[ktemp];
  }

 for(igrp=1;igrp <= ngrp ; igrp++) {
    ktemp = ind4[igrp];
    xo[4][igrp] = clatoms_xold[ktemp];
    yo[4][igrp] = clatoms_yold[ktemp];
    zo[4][igrp] = clatoms_zold[ktemp];
    rm4[igrp] = 1.0/clatoms_mass[ktemp];
  }


/* assign jtyps  */
 for(igrp=1;igrp <= ngrp ; igrp++) {
    jtyp = grp_bond_con_jtyp_43[(igrp+igrp_off)];
    dij[1][igrp] = grp_bond_con_eq_43[1][jtyp];
    dij[2][igrp] = grp_bond_con_eq_43[2][jtyp];
    dij[3][igrp] = grp_bond_con_eq_43[3][jtyp];
 }/*endfor*/

/* compute the reduced mass matrix elements */
 for(igrp=1;igrp <= ngrp ; igrp++) {

    rms1 = rm1[igrp];  rms2 = rm2[igrp]; 
    rms3 = rm3[igrp];  rms4 = rm4[igrp];

    rmm11[igrp] = (rms1+rms2); rmm12[igrp] = rmm13[igrp] = rms1;
    rmm21[igrp] = rms1; rmm22[igrp] = (rms1+rms3); rmm23[igrp] = rms1;
    rmm31[igrp] = rms1; rmm32[igrp] = rms1; rmm33[igrp] = (rms1+rms4);
 }/*endfor*/


/* Compute difference vectors */
 for(igrp=1;igrp <= ngrp ; igrp++) {
    dxt[1][igrp]= x[1][igrp] - x[2][igrp];
    dxt[2][igrp]= x[1][igrp] - x[3][igrp];
    dxt[3][igrp]= x[1][igrp] - x[4][igrp];
  }
     
 for(igrp=1;igrp <= ngrp ; igrp++) {
    dyt[1][igrp]= y[1][igrp] - y[2][igrp];
    dyt[2][igrp]= y[1][igrp] - y[3][igrp];
    dyt[3][igrp]= y[1][igrp] - y[4][igrp];
  }
     
 for(igrp=1;igrp <= ngrp ; igrp++) {
    dzt[1][igrp]= z[1][igrp] - z[2][igrp];
    dzt[2][igrp]= z[1][igrp] - z[3][igrp];
    dzt[3][igrp]= z[1][igrp] - z[4][igrp];
  }
     
 for(igrp=1;igrp <= ngrp ; igrp++) {
    dxo[1][igrp]= (xo[1][igrp] - xo[2][igrp]);
    dxo[2][igrp]= (xo[1][igrp] - xo[3][igrp]);
    dxo[3][igrp]= (xo[1][igrp] - xo[4][igrp]);
  }
     
 for(igrp=1;igrp <= ngrp ; igrp++) {
    dx[1][igrp]= dxo[1][igrp]*baro_roll_scv;
    dx[2][igrp]= dxo[2][igrp]*baro_roll_scv;
    dx[3][igrp]= dxo[3][igrp]*baro_roll_scv;
 }

 for(igrp=1;igrp <= ngrp ; igrp++) {
    dyo[1][igrp]= (yo[1][igrp] - yo[2][igrp]);
    dyo[2][igrp]= (yo[1][igrp] - yo[3][igrp]);
    dyo[3][igrp]= (yo[1][igrp] - yo[4][igrp]);
  }
     
 for(igrp=1;igrp <= ngrp ; igrp++) {
    dy[1][igrp]= dyo[1][igrp]*baro_roll_scv;
    dy[2][igrp]= dyo[2][igrp]*baro_roll_scv;
    dy[3][igrp]= dyo[3][igrp]*baro_roll_scv;
 }

 for(igrp=1;igrp <= ngrp ; igrp++) {
    dzo[1][igrp]= (zo[1][igrp] - zo[2][igrp]);
    dzo[2][igrp]= (zo[1][igrp] - zo[3][igrp]);
    dzo[3][igrp]= (zo[1][igrp] - zo[4][igrp]);
 }

 for(igrp=1;igrp <= ngrp ; igrp++) {
    dz[1][igrp]= dzo[1][igrp]*baro_roll_scv;
    dz[2][igrp]= dzo[2][igrp]*baro_roll_scv;
    dz[3][igrp]= dzo[3][igrp]*baro_roll_scv;
 }


 for(igrp=1;igrp <= ngrp ; igrp++) {
    avec[1][igrp]= dij[1][igrp]*dij[1][igrp] - (dxt[1][igrp]*dxt[1][igrp]
                  +dyt[1][igrp]*dyt[1][igrp] + dzt[1][igrp]*dzt[1][igrp]);

    avec[2][igrp]= dij[2][igrp]*dij[2][igrp] - (dxt[2][igrp]*dxt[2][igrp]
                  +dyt[2][igrp]*dyt[2][igrp] + dzt[2][igrp]*dzt[2][igrp]);

    avec[3][igrp]= dij[3][igrp]*dij[3][igrp] - (dxt[3][igrp]*dxt[3][igrp]
                  +dyt[3][igrp]*dyt[3][igrp] + dzt[3][igrp]*dzt[3][igrp]);
 }


/* Get initial guess for lambda */

 if(ifirst == 2 || ifirst == 0){
#ifndef NO_PRAGMA
#pragma IVDEP
#endif
 for(igrp=1;igrp <= ngrp ; igrp++) {
    amat[1][1] = -2.0*rmm11[igrp]* (dxt[1][igrp]*dx[1][igrp] 
               + dyt[1][igrp]*dy[1][igrp] + dzt[1][igrp]*dz[1][igrp]);
    amat[1][2] = -2.0*rmm12[igrp]* (dxt[1][igrp]*dx[2][igrp] 
               + dyt[1][igrp]*dy[2][igrp] + dzt[1][igrp]*dz[2][igrp]);
    amat[1][3] = -2.0*rmm13[igrp]* (dxt[1][igrp]*dx[3][igrp] 
               + dyt[1][igrp]*dy[3][igrp] + dzt[1][igrp]*dz[3][igrp]);

    amat[2][1] = -2.0*rmm21[igrp]* (dxt[2][igrp]*dx[1][igrp] 
               + dyt[2][igrp]*dy[1][igrp] + dzt[2][igrp]*dz[1][igrp]);
    amat[2][2] = -2.0*rmm22[igrp]* (dxt[2][igrp]*dx[2][igrp] 
               + dyt[2][igrp]*dy[2][igrp] + dzt[2][igrp]*dz[2][igrp]);
    amat[2][3] = -2.0*rmm23[igrp]* (dxt[2][igrp]*dx[3][igrp] 
               + dyt[2][igrp]*dy[3][igrp] + dzt[2][igrp]*dz[3][igrp]);

    amat[3][1] = -2.0*rmm31[igrp]* (dxt[3][igrp]*dx[1][igrp] 
               + dyt[3][igrp]*dy[1][igrp] + dzt[3][igrp]*dz[1][igrp]);
    amat[3][2] = -2.0*rmm32[igrp]* (dxt[3][igrp]*dx[2][igrp] 
               + dyt[3][igrp]*dy[2][igrp] + dzt[3][igrp]*dz[2][igrp]);
    amat[3][3] = -2.0*rmm33[igrp]* (dxt[3][igrp]*dx[3][igrp] 
               + dyt[3][igrp]*dy[3][igrp] + dzt[3][igrp]*dz[3][igrp]);


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
 } else {
 for(igrp=1;igrp <= ngrp ; igrp++) {
   xlam[1][igrp] = grp_bond_con_al_43[1][(igrp+igrp_off)];
   xlam[2][igrp] = grp_bond_con_al_43[2][(igrp+igrp_off)];
   xlam[3][igrp] = grp_bond_con_al_43[3][(igrp+igrp_off)];
   grp_bond_con_al_43[1][(igrp+igrp_off)] = 0.0;
   grp_bond_con_al_43[2][(igrp+igrp_off)] = 0.0;
   grp_bond_con_al_43[3][(igrp+igrp_off)] = 0.0;
 }/*endfor*/
}


/* Iterative do loop for multiplier */

 if(ngrp > 0){

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
   for(igrp=1;igrp <= ngrp ; igrp++) {
/* Set up guess of difference vectors */

    dxn[1] = -2.0*dxt[1][igrp]; 
    dyn[1] = -2.0*dyt[1][igrp];
    dzn[1] = -2.0*dzt[1][igrp];

    dxn[2] = -2.0*dxt[2][igrp]; 
    dyn[2] = -2.0*dyt[2][igrp];
    dzn[2] = -2.0*dzt[2][igrp];

    dxn[3] = -2.0*dxt[3][igrp]; 
    dyn[3] = -2.0*dyt[3][igrp];
    dzn[3] = -2.0*dzt[3][igrp];

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

 }/*endif for ngrp > 0*/

/* position update */

#ifndef NO_PRAGMA
#pragma IVDEP
#endif
 for(igrp=1;igrp <= ngrp ; igrp++) {
    double xlam1,xlam2,xlam3;
    double rms1,rms2,rms3,rms4;
    double dx1,dx2,dx3;
    double dy1,dy2,dy3;
    double dz1,dz2,dz3;
    int ktemp1,ktemp2,ktemp3,ktemp4;

     rms1=rm1[igrp];
     rms2=rm2[igrp];
     rms3=rm3[igrp];
     rms4=rm4[igrp];

     xlam1=xlam[1][igrp];
     xlam2=xlam[2][igrp];
     xlam3=xlam[3][igrp];
    
     dx1=dxo[1][igrp]; dx2=dxo[2][igrp]; dx3=dxo[3][igrp]; 
     dy1=dyo[1][igrp]; dy2=dyo[2][igrp]; dy3=dyo[3][igrp];
     dz1=dzo[1][igrp]; dz2=dzo[2][igrp]; dz3=dzo[3][igrp];

     ktemp1 = ind1[igrp]; ktemp2 = ind2[igrp]; 
     ktemp3 = ind3[igrp]; ktemp4 = ind4[igrp];

    clatoms_x[ktemp1] -= (xlam1*dx1 + xlam2*dx2 + xlam3*dx3)*rms1*baro_roll_scv;
    clatoms_y[ktemp1] -= (xlam1*dy1 + xlam2*dy2 + xlam3*dy3)*rms1*baro_roll_scv;
    clatoms_z[ktemp1] -= (xlam1*dz1 + xlam2*dz2 + xlam3*dz3)*rms1*baro_roll_scv;

    clatoms_vx[ktemp1] -= (xlam1*dx1 + xlam2*dx2 + xlam3*dx3)*rms1/dt;
    clatoms_vy[ktemp1] -= (xlam1*dy1 + xlam2*dy2 + xlam3*dy3)*rms1/dt;
    clatoms_vz[ktemp1] -= (xlam1*dz1 + xlam2*dz2 + xlam3*dz3)*rms1/dt;

    clatoms_x[ktemp2] += (xlam1*dx1)*rms2*baro_roll_scv;
    clatoms_y[ktemp2] += (xlam1*dy1)*rms2*baro_roll_scv;
    clatoms_z[ktemp2] += (xlam1*dz1)*rms2*baro_roll_scv;

    clatoms_vx[ktemp2] += (xlam1*dx1)*rms2/dt;
    clatoms_vy[ktemp2] += (xlam1*dy1)*rms2/dt;
    clatoms_vz[ktemp2] += (xlam1*dz1)*rms2/dt;

    clatoms_x[ktemp3] += (xlam2*dx2)*rms3*baro_roll_scv;
    clatoms_y[ktemp3] += (xlam2*dy2)*rms3*baro_roll_scv;
    clatoms_z[ktemp3] += (xlam2*dz2)*rms3*baro_roll_scv;

    clatoms_vx[ktemp3] += (xlam2*dx2)*rms3/dt;
    clatoms_vy[ktemp3] += (xlam2*dy2)*rms3/dt;
    clatoms_vz[ktemp3] += (xlam2*dz2)*rms3/dt;

    clatoms_x[ktemp4] += (xlam3*dx3)*rms4*baro_roll_scv;
    clatoms_y[ktemp4] += (xlam3*dy3)*rms4*baro_roll_scv;
    clatoms_z[ktemp4] += (xlam3*dz3)*rms4*baro_roll_scv;

    clatoms_vx[ktemp4] += (xlam3*dx3)*rms4/dt;
    clatoms_vy[ktemp4] += (xlam3*dy3)*rms4/dt;
    clatoms_vz[ktemp4] += (xlam3*dz3)*rms4/dt;

  }


/* Pressure tensor update */
/* unscaled old distances */

 for(igrp=1;igrp <= ngrp ; igrp++) {
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
    
     dx1=dxo[1][igrp]; dx2=dxo[2][igrp]; dx3=dxo[3][igrp];
     dy1=dyo[1][igrp]; dy2=dyo[2][igrp]; dy3=dyo[3][igrp];
     dz1=dzo[1][igrp]; dz2=dzo[2][igrp]; dz3=dzo[3][igrp];

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
 for(igrp=1;igrp <= ngrp ; igrp++) {
     ptens_pvten_tmp[1] -= (p11[igrp]*pnorm); 
     ptens_pvten_tmp[2] -= (p12[igrp]*pnorm);
     ptens_pvten_tmp[3] -= (p13[igrp]*pnorm);
     ptens_pvten_tmp[4] -= (p12[igrp]*pnorm);
     ptens_pvten_tmp[5] -= (p22[igrp]*pnorm);
     ptens_pvten_tmp[6] -= (p23[igrp]*pnorm);
     ptens_pvten_tmp[7] -= (p13[igrp]*pnorm);
     ptens_pvten_tmp[8] -= (p23[igrp]*pnorm);
     ptens_pvten_tmp[9] -= (p33[igrp]*pnorm);
 }/*end for*/

/* Save multiplier */

 for(igrp=1;igrp <= ngrp ; igrp++) {
    grp_bond_con_al_43[1][(igrp+igrp_off)] += xlam[1][igrp];
    grp_bond_con_al_43[2][(igrp+igrp_off)] += xlam[2][igrp];
    grp_bond_con_al_43[3][(igrp+igrp_off)] += xlam[3][igrp];
 } /* end for */

/*=======================================================================*/
/*  IV)Allreduce pvten_tmp     */

  if(np_forc > 1 ){
   for(i=1;i<=9;i++){
    ptens_pvten_tmp2[i] = ptens_pvten_tmp[i];
   }/*endfor*/
   Allreduce(&(ptens_pvten_tmp2[1]), &(ptens_pvten_tmp[1]),9,MPI_DOUBLE,
                   MPI_SUM,0,comm_forc);
  }/*endif*/

 for(i=1; i<=9; i++){
  ptens_pvten_inc[i] += ptens_pvten_tmp[i];
 }
 if(ifirst == 0){
  ftemp   = (ptens_pvten_tmp[1]+ptens_pvten_tmp[5]+ptens_pvten_tmp[9]);
  baro->f_lnv_p += ftemp;
  baro->v_lnv   += 0.5*ftemp*(baro->roll_scg)*dt/(baro->mass_lnv);
 }

/* free all the memory allocated locally */
 if(ngrp > 0){
     free_dmatrix(dx,1,3,1,ngrp);
     free_dmatrix(dy,1,3,1,ngrp);
     free_dmatrix(dz,1,3,1,ngrp);

     free_dmatrix(dxo,1,3,1,ngrp);
     free_dmatrix(dyo,1,3,1,ngrp);
     free_dmatrix(dzo,1,3,1,ngrp);

     free_dmatrix(x,1,4,1,ngrp); 
     free_dmatrix(y,1,4,1,ngrp); 
     free_dmatrix(z,1,4,1,ngrp); 

     free_dmatrix(xo,1,4,1,ngrp);
     free_dmatrix(yo,1,4,1,ngrp);
     free_dmatrix(zo,1,4,1,ngrp);

     free_dvector(rm1,1,ngrp);
     free_dvector(rm2,1,ngrp);
     free_dvector(rm3,1,ngrp);
     free_dvector(rm4,1,ngrp);

     free_dmatrix(dxt,1,3,1,ngrp);
     free_dmatrix(dyt,1,3,1,ngrp);
     free_dmatrix(dzt,1,3,1,ngrp);

     free_dmatrix(xlam,1,3,1,ngrp);
     free_dmatrix(avec,1,3,1,ngrp);
     free_dmatrix(dxl,1,3,1,ngrp);
     free_dmatrix(dij,1,3,1,ngrp);

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

    free(ind1);
    free(ind2);
    free(ind3);
    free(ind4);
 }/*endif*/

/*=======================================================================*/
} /* end routine */
/*=======================================================================*/


/*==========================================================================*/

void rattle_43_rolli(GRP_BOND_CON *grp_bond_con,
              CLATOMS_INFO *clatoms_info,CLATOMS_POS *clatoms_pos,
              PTENS *ptens,double dt,BARO *baro,int ifirst,
              CLASS_COMM_FORC_PKG *class_comm_forc_pkg)

/*==========================================================================*/
/*        Begin Routine                                                     */
{/*Begin Routine*/
/*=======================================================================*/
/*         Local Variable declarations                                   */

#include "../typ_defs/typ_mask.h"

  double avec[NCON_43+1];
  double amat[NCON_43+1][NCON_43+1],ainv[NCON_43+1][NCON_43+1];
  double rms1,rms2,rms3,rms4;
  double det,rdet_a;
  double roll_sci,f_lnv_inc;
  double ftemp;

  int i,j,k,iii;
  int ktemp,ktemp4;
  int igrp,*ind1,*ind2,*ind3,*ind4,jtyp;

   double **x,**y,**z;
   double **vx,**vy,**vz;
   double *p11,*p12,*p13,*p23,*p22,*p33;
   double *rm1,*rm2,*rm3,*rm4;
   double **dx,**dy,**dz;
   double **dvx,**dvy,**dvz;
   double *rmm11,*rmm12,*rmm13,*rmm21,*rmm22,*rmm23,*rmm31,*rmm32,*rmm33;
   double **xlam;

/* Local pointers */
  double *clatoms_mass         = clatoms_info->mass;
  double *clatoms_x            = clatoms_pos->x;
  double *clatoms_y            = clatoms_pos->y;
  double *clatoms_z            = clatoms_pos->z;
  double *clatoms_vx           = clatoms_pos->vx;
  double *clatoms_vy           = clatoms_pos->vy;
  double *clatoms_vz           = clatoms_pos->vz;
  int *grp_bond_con_j1_43      = grp_bond_con->j1_43;
  int *grp_bond_con_j2_43      = grp_bond_con->j2_43;
  int *grp_bond_con_j3_43      = grp_bond_con->j3_43;
  int *grp_bond_con_j4_43      = grp_bond_con->j4_43;
  int *grp_bond_con_jtyp_43    = grp_bond_con->jtyp_43;
  double *ptens_pvten_inc      = ptens->pvten_inc;
  double *ptens_pvten_tmp      = ptens->pvten_tmp;
  double *ptens_pvten_tmp2      = ptens->pvten_tmp_res;
  double pnorm;
  double *clatoms_roll_sc = clatoms_info->roll_sc;
  double baro_v_lnv_g = baro->v_lnv_g;


  int ngrp,irem,igrp_off;
  int ngrp_tot                 = grp_bond_con->num_43;
  int np_forc                  = class_comm_forc_pkg->num_proc;
  int myid_forc                = class_comm_forc_pkg->myid;
  MPI_Comm comm_forc           = class_comm_forc_pkg->comm;

/*=======================================================================*/

  ngrp = (ngrp_tot);
  igrp_off = 0;

/*=======================================================================*/

  if(ngrp > 0){
   x= dmatrix(1,4,1,ngrp);
   y= dmatrix(1,4,1,ngrp);
   z= dmatrix(1,4,1,ngrp);
  vx= dmatrix(1,4,1,ngrp);
  vy= dmatrix(1,4,1,ngrp);
  vz= dmatrix(1,4,1,ngrp);

  p11= dvector(1,ngrp);
  p12= dvector(1,ngrp);
  p13= dvector(1,ngrp);
  p22= dvector(1,ngrp);
  p23= dvector(1,ngrp);
  p33= dvector(1,ngrp);

  rm1= dvector(1,ngrp);
  rm2= dvector(1,ngrp);
  rm3= dvector(1,ngrp);
  rm4= dvector(1,ngrp);

  dx= dmatrix(1,3,1,ngrp);
  dy= dmatrix(1,3,1,ngrp);
  dz= dmatrix(1,3,1,ngrp);
 dvx= dmatrix(1,3,1,ngrp);
 dvy= dmatrix(1,3,1,ngrp);
 dvz= dmatrix(1,3,1,ngrp);
 xlam= dmatrix(1,3,1,ngrp);

 ind1 = (int *) calloc((ngrp+1),sizeof(int));
 ind2 = (int *) calloc((ngrp+1),sizeof(int));
 ind3 = (int *) calloc((ngrp+1),sizeof(int));
 ind4 = (int *) calloc((ngrp+1),sizeof(int));

  rmm11 = dvector(1,ngrp);
  rmm12 = dvector(1,ngrp);
  rmm13 = dvector(1,ngrp);
  rmm21 = dvector(1,ngrp);
  rmm22 = dvector(1,ngrp);
  rmm23 = dvector(1,ngrp);
  rmm31 = dvector(1,ngrp);
  rmm32 = dvector(1,ngrp);
  rmm33 = dvector(1,ngrp);
 }/*endif*/
/*=======================================================================*/
 pnorm = 2.0/dt;

/* Set zeta matrix and reciprocal mass matrix */

 ptens_pvten_tmp[1] = 0.0;
 ptens_pvten_tmp[2] = 0.0;
 ptens_pvten_tmp[3] = 0.0;
 ptens_pvten_tmp[4] = 0.0;
 ptens_pvten_tmp[5] = 0.0;
 ptens_pvten_tmp[6] = 0.0;
 ptens_pvten_tmp[7] = 0.0;
 ptens_pvten_tmp[8] = 0.0;
 ptens_pvten_tmp[9] = 0.0;

 for(igrp=1;igrp <= ngrp ; igrp++) {

  ind1[igrp] = grp_bond_con_j1_43[(igrp+igrp_off)];
  ind2[igrp] = grp_bond_con_j2_43[(igrp+igrp_off)];
  ind3[igrp] = grp_bond_con_j3_43[(igrp+igrp_off)];
  ind4[igrp] = grp_bond_con_j4_43[(igrp+igrp_off)];
 }

 for(igrp=1;igrp <= ngrp ; igrp++) {
  ktemp = ind1[igrp];
  x[1][igrp] = clatoms_x[ktemp];
  y[1][igrp] = clatoms_y[ktemp];
  z[1][igrp] = clatoms_z[ktemp];
  rm1[igrp] = 1.0/clatoms_mass[ktemp];
 }

 for(igrp=1;igrp <= ngrp ; igrp++) {
  ktemp = ind2[igrp];
  x[2][igrp] = clatoms_x[ktemp];
  y[2][igrp] = clatoms_y[ktemp];
  z[2][igrp] = clatoms_z[ktemp];
  rm2[igrp] = 1.0/clatoms_mass[ktemp];
 }

 for(igrp=1;igrp <= ngrp ; igrp++) {
  ktemp = ind3[igrp];
  x[3][igrp] = clatoms_x[ktemp];
  y[3][igrp] = clatoms_y[ktemp];
  z[3][igrp] = clatoms_z[ktemp];
  rm3[igrp] = 1.0/clatoms_mass[ktemp];
}

 for(igrp=1;igrp <= ngrp ; igrp++) {
  ktemp = ind4[igrp];
  x[4][igrp] = clatoms_x[ktemp];
  y[4][igrp] = clatoms_y[ktemp];
  z[4][igrp] = clatoms_z[ktemp];
  rm4[igrp] = 1.0/clatoms_mass[ktemp];
}


 for(igrp=1;igrp <= ngrp ; igrp++) {
  ktemp = ind1[igrp];
  ktemp4 = ind4[igrp];
  roll_sci=1.0/clatoms_roll_sc[ktemp4];/*all roll scales the same in same cons*/
  vx[1][igrp] = clatoms_vx[ktemp] + x[1][igrp]*baro_v_lnv_g*roll_sci;
  vy[1][igrp] = clatoms_vy[ktemp] + y[1][igrp]*baro_v_lnv_g*roll_sci;
  vz[1][igrp] = clatoms_vz[ktemp] + z[1][igrp]*baro_v_lnv_g*roll_sci;
 }

 for(igrp=1;igrp <= ngrp ; igrp++) {
  ktemp = ind2[igrp];
  ktemp4 = ind4[igrp];
  roll_sci=1.0/clatoms_roll_sc[ktemp4];/*all roll scales the same in same cons*/
  vx[2][igrp] = clatoms_vx[ktemp] + x[2][igrp]*baro_v_lnv_g*roll_sci;
  vy[2][igrp] = clatoms_vy[ktemp] + y[2][igrp]*baro_v_lnv_g*roll_sci;
  vz[2][igrp] = clatoms_vz[ktemp] + z[2][igrp]*baro_v_lnv_g*roll_sci;
 }

 for(igrp=1;igrp <= ngrp ; igrp++) {
  ktemp = ind3[igrp];
  ktemp4 = ind4[igrp];
  roll_sci=1.0/clatoms_roll_sc[ktemp4];/*all roll scales the same in same cons*/
  vx[3][igrp] = clatoms_vx[ktemp] + x[3][igrp]*baro_v_lnv_g*roll_sci;
  vy[3][igrp] = clatoms_vy[ktemp] + y[3][igrp]*baro_v_lnv_g*roll_sci;
  vz[3][igrp] = clatoms_vz[ktemp] + z[3][igrp]*baro_v_lnv_g*roll_sci;
 }

 for(igrp=1;igrp <= ngrp ; igrp++) {
  ktemp = ind4[igrp];
  roll_sci=1.0/clatoms_roll_sc[ktemp];/*all roll scales the same in same cons*/
  vx[4][igrp] = clatoms_vx[ktemp] + x[4][igrp]*baro_v_lnv_g*roll_sci;
  vy[4][igrp] = clatoms_vy[ktemp] + y[4][igrp]*baro_v_lnv_g*roll_sci;
  vz[4][igrp] = clatoms_vz[ktemp] + z[4][igrp]*baro_v_lnv_g*roll_sci;
 }



/* Compute difference vectors */
 for(igrp=1;igrp <= ngrp ; igrp++) {

   dvx[1][igrp] = vx[1][igrp]-vx[2][igrp];
   dvx[2][igrp] = vx[1][igrp]-vx[3][igrp];
   dvx[3][igrp] = vx[1][igrp]-vx[4][igrp];

   dvy[1][igrp] = vy[1][igrp]-vy[2][igrp];
   dvy[2][igrp] = vy[1][igrp]-vy[3][igrp];
   dvy[3][igrp] = vy[1][igrp]-vy[4][igrp];

   dvz[1][igrp] = vz[1][igrp]-vz[2][igrp];
   dvz[2][igrp] = vz[1][igrp]-vz[3][igrp];
   dvz[3][igrp] = vz[1][igrp]-vz[4][igrp];

   dx[1][igrp] = x[1][igrp]-x[2][igrp];
   dx[2][igrp] = x[1][igrp]-x[3][igrp];
   dx[3][igrp] = x[1][igrp]-x[4][igrp];

   dy[1][igrp] = y[1][igrp]-y[2][igrp];
   dy[2][igrp] = y[1][igrp]-y[3][igrp];
   dy[3][igrp] = y[1][igrp]-y[4][igrp];

   dz[1][igrp] = z[1][igrp]-z[2][igrp];
   dz[2][igrp] = z[1][igrp]-z[3][igrp];
   dz[3][igrp] = z[1][igrp]-z[4][igrp];

} 

 for(igrp=1;igrp <= ngrp ; igrp++) {
  rms1 = rm1[igrp];
  rms2 = rm2[igrp];
  rms3 = rm3[igrp];
  rms4 = rm4[igrp];

  rmm11[igrp] = (rms1+rms2); rmm12[igrp] = rms1; rmm13[igrp] = rms1;
  rmm21[igrp] = rms1; rmm22[igrp] = (rms1+rms3); rmm23[igrp] = rms1;
  rmm31[igrp] = rms1; rmm32[igrp] = rms1; rmm33[igrp] = (rms1+rms4);
 }


/* Get lambda */
#ifndef NO_PRAGMA
#pragma IVDEP
#endif
 for(igrp=1;igrp <= ngrp ; igrp++) {
  amat[1][1] = rmm11[igrp]*
               (dx[1][igrp]*dx[1][igrp] + dy[1][igrp]*dy[1][igrp] + dz[1][igrp]*dz[1][igrp]);
  amat[1][2] = rmm12[igrp]*
               (dx[1][igrp]*dx[2][igrp] + dy[1][igrp]*dy[2][igrp] + dz[1][igrp]*dz[2][igrp]);
  amat[1][3] = rmm13[igrp]*
               (dx[1][igrp]*dx[3][igrp] + dy[1][igrp]*dy[3][igrp] + dz[1][igrp]*dz[3][igrp]);

   amat[2][1] = rmm21[igrp]*
               (dx[2][igrp]*dx[1][igrp] + dy[2][igrp]*dy[1][igrp] + dz[2][igrp]*dz[1][igrp]);
  amat[2][2] = rmm22[igrp]*
               (dx[2][igrp]*dx[2][igrp] + dy[2][igrp]*dy[2][igrp] + dz[2][igrp]*dz[2][igrp]);
  amat[2][3] = rmm23[igrp]*
               (dx[2][igrp]*dx[3][igrp] + dy[2][igrp]*dy[3][igrp] + dz[2][igrp]*dz[3][igrp]);

   amat[3][1] = rmm31[igrp]*
        (dx[3][igrp]*dx[1][igrp] + dy[3][igrp]*dy[1][igrp] + dz[3][igrp]*dz[1][igrp]);
   amat[3][2] = rmm32[igrp]*
       (dx[3][igrp]*dx[2][igrp] + dy[3][igrp]*dy[2][igrp] + dz[3][igrp]*dz[2][igrp]);
   amat[3][3] = rmm33[igrp]*
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
 for(igrp=1;igrp <= ngrp ; igrp++) {
  double xlam1,xlam2,xlam3;
  double dx1,dx2,dx3;
  double dy1,dy2,dy3;
  double dz1,dz2,dz3;
  int ktemp1,ktemp2,ktemp3,ktemp4;


    dx1= dx[1][igrp]; dx2= dx[2][igrp]; dx3= dx[3][igrp];
    dy1= dy[1][igrp]; dy2= dy[2][igrp]; dy3= dy[3][igrp];
    dz1= dz[1][igrp]; dz2= dz[2][igrp]; dz3= dz[3][igrp];

    xlam1= xlam[1][igrp];
    xlam2= xlam[2][igrp];
    xlam3= xlam[3][igrp];

    ktemp1=ind1[igrp]; ktemp2=ind2[igrp];
    ktemp3=ind3[igrp]; ktemp4=ind4[igrp];

    rms1 = rm1[igrp]; rms2 = rm2[igrp];
    rms3 = rm3[igrp]; rms4 = rm4[igrp];


    clatoms_vx[ktemp1] -= (xlam1*dx1 + xlam2*dx2 + xlam3*dx3)*rms1;
    clatoms_vy[ktemp1] -= (xlam1*dy1 + xlam2*dy2 + xlam3*dy3)*rms1;
    clatoms_vz[ktemp1] -= (xlam1*dz1 + xlam2*dz2 + xlam3*dz3)*rms1;

    clatoms_vx[ktemp2] += (xlam1*dx1)*rms2;
    clatoms_vy[ktemp2] += (xlam1*dy1)*rms2;
    clatoms_vz[ktemp2] += (xlam1*dz1)*rms2;

    clatoms_vx[ktemp3] += (xlam2*dx2)*rms3;
    clatoms_vy[ktemp3] += (xlam2*dy2)*rms3;
    clatoms_vz[ktemp3] += (xlam2*dz2)*rms3;

    clatoms_vx[ktemp4] += (xlam3*dx3)*rms4;
    clatoms_vy[ktemp4] += (xlam3*dy3)*rms4;
    clatoms_vz[ktemp4] += (xlam3*dz3)*rms4;

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
 for(igrp=1;igrp <= ngrp ; igrp++) {
    ptens_pvten_tmp[1] -= (p11[igrp]*pnorm);
    ptens_pvten_tmp[2] -= (p12[igrp]*pnorm);
    ptens_pvten_tmp[3] -= (p13[igrp]*pnorm);
    ptens_pvten_tmp[4] -= (p12[igrp]*pnorm);
    ptens_pvten_tmp[5] -= (p22[igrp]*pnorm);
    ptens_pvten_tmp[6] -= (p23[igrp]*pnorm);
    ptens_pvten_tmp[7] -= (p13[igrp]*pnorm);
    ptens_pvten_tmp[8] -= (p23[igrp]*pnorm);
    ptens_pvten_tmp[9] -= (p33[igrp]*pnorm);
 } /* end for igrp */

/*=======================================================================*/
/*  IV)Allreduce pvten_tmp     */

  if(np_forc > 1 ){
   for(i=1;i<=9;i++){
    ptens_pvten_tmp2[i] = ptens_pvten_tmp[i];
   }/*endfor*/
   Allreduce(&(ptens_pvten_tmp2[1]), &(ptens_pvten_tmp[1]),9,MPI_DOUBLE,
                   MPI_SUM,0,comm_forc);
  }/*endif*/

  ptens_pvten_inc[1] += ptens_pvten_tmp[1];
  ptens_pvten_inc[2] += ptens_pvten_tmp[2];
  ptens_pvten_inc[3] += ptens_pvten_tmp[3];
  ptens_pvten_inc[4] += ptens_pvten_tmp[4];
  ptens_pvten_inc[5] += ptens_pvten_tmp[5];
  ptens_pvten_inc[6] += ptens_pvten_tmp[6];
  ptens_pvten_inc[7] += ptens_pvten_tmp[7];
  ptens_pvten_inc[8] += ptens_pvten_tmp[8];
  ptens_pvten_inc[9] += ptens_pvten_tmp[9];


 if(ifirst == 0){
  f_lnv_inc   = (ptens_pvten_tmp[1]+ptens_pvten_tmp[5]
               +ptens_pvten_tmp[9]);
  baro->f_lnv_p += f_lnv_inc;
  baro->v_lnv_g += f_lnv_inc*(baro->roll_scg)*0.5*dt/(baro->mass_lnv);
 }

 /* free locally assigned memory */
 if(ngrp > 0){
   free_dmatrix(x,1,4,1,ngrp);
   free_dmatrix(y,1,4,1,ngrp);
   free_dmatrix(z,1,4,1,ngrp);

   free_dmatrix(vx,1,4,1,ngrp);
   free_dmatrix(vy,1,4,1,ngrp);
   free_dmatrix(vz,1,4,1,ngrp);

   free_dvector(p11,1,ngrp);
   free_dvector(p12,1,ngrp);
   free_dvector(p13,1,ngrp);
   free_dvector(p22,1,ngrp);
   free_dvector(p23,1,ngrp);
   free_dvector(p33,1,ngrp);

   free_dvector(rm1,1,ngrp);
   free_dvector(rm2,1,ngrp);
   free_dvector(rm3,1,ngrp);
   free_dvector(rm4,1,ngrp);

   free_dmatrix(dx,1,3,1,ngrp);
   free_dmatrix(dy,1,3,1,ngrp);
   free_dmatrix(dz,1,3,1,ngrp);

   free_dmatrix(dvx,1,3,1,ngrp);
   free_dmatrix(dvy,1,3,1,ngrp);
   free_dmatrix(dvz,1,3,1,ngrp);

   free_dmatrix(xlam,1,3,1,ngrp);

   free(ind1);
   free(ind2);
   free(ind3);
   free(ind4);

   free_dvector(rmm11,1,ngrp); 
   free_dvector(rmm12,1,ngrp); 
   free_dvector(rmm13,1,ngrp); 
   free_dvector(rmm21,1,ngrp); 
   free_dvector(rmm22,1,ngrp); 
   free_dvector(rmm23,1,ngrp); 
   free_dvector(rmm31,1,ngrp); 
   free_dvector(rmm32,1,ngrp); 
   free_dvector(rmm33,1,ngrp); 
 }/*endif*/

/*=======================================================================*/
} /* end routine */
/*=======================================================================*/
