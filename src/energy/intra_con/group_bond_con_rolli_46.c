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
/* ===================================================================== */
#include "standard_include.h"
#include "../typ_defs/typedefs_class.h"
#include "../typ_defs/typedefs_gen.h"
#include "../typ_defs/typedefs_bnd.h"
#include "../proto_defs/proto_intra_con_local.h"
#include "../proto_defs/proto_math.h"
#include "../proto_defs/proto_friend_lib_entry.h"
#include "../proto_defs/proto_communicate_wrappers.h"

#define NCON_46 6
#define NAT_46 4


/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
void shake_46_rolli(GRP_BOND_CON *grp_bond_con,
              CLATOMS_INFO *clatoms_info,CLATOMS_POS *clatoms_pos,
              PTENS *ptens,double dt,double *aiter,
              BARO *baro, int ifirst,
              CLASS_COMM_FORC_PKG *class_comm_forc_pkg)


/*==========================================================================*/
/*        Begin Routine                                                     */
{/*Begin Routine*/
/*=======================================================================*/
/*         Local Variable declarations                                   */

#include "../typ_defs/typ_mask.h"

 double xl0[NCON_46+1],dmax;


 double rms1,rms2,rms3,rms4;
 double ftemp;
 double dts;
 int i,j,iii;
 int iter,igrp,*ind1,*ind2,*ind3,*ind4,jtyp;
 int ktemp,ktemp1,ktemp2,ktemp3,ktemp4;
 int na,job,info,ipvt[NCON_46+1];       /* For dgefa and dgesl */
/* AAA */

 double *rmass1,*rmass2,*rmass3,*rmass4,*dlmax,*txlam,*tamat;
 double **dx,**dy,**dz;
 double **dxt,**dyt,**dzt;
 double **dxn,**dyn,**dzn;
 double **xlam,**avec,**dxl,**dij,**amat;
 double ***rmassm;
 double **x,**y,**z;
 double **xo,**yo,**zo;
 double *p11,*p12,*p13,*p22,*p23,*p33; 

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
  int *grp_bond_con_j1_46      = grp_bond_con->j1_46;
  int *grp_bond_con_j2_46      = grp_bond_con->j2_46;
  int *grp_bond_con_j3_46      = grp_bond_con->j3_46;
  int *grp_bond_con_j4_46      = grp_bond_con->j4_46;
  int *grp_bond_con_jtyp_46    = grp_bond_con->jtyp_46;
  double **grp_bond_con_eq_46  = grp_bond_con->eq_46;
  double **grp_bond_con_al_46  = grp_bond_con->al_46;
  double *ptens_pvten_inc      = ptens->pvten_inc;
  double *ptens_pvten_tmp      = ptens->pvten_tmp;
  double *ptens_pvten_tmp2     = ptens->pvten_tmp_res;
  double pnorm;
  double baro_roll_scv         = baro->roll_scv;

  int ngrp,irem,igrp_off;
  int ngrp_tot                 = grp_bond_con->num_46;
  int np_forc                  = class_comm_forc_pkg->num_proc;
  int myid_forc                = class_comm_forc_pkg->myid;
  MPI_Comm comm_forc           = class_comm_forc_pkg->comm;

/*=======================================================================*/

  ngrp = (ngrp_tot);
  igrp_off = 0;

/*=======================================================================*/

  if(ngrp > 0){
     dlmax= dvector(1,6); 
     txlam= dvector(1,6); 
    rmassm= d3tensor(1,6,1,6,1,ngrp);
     amat= dmatrix(1,36,1,ngrp);
     xlam= dmatrix(1,6,1,ngrp);  
     avec= dmatrix(1,6,1,ngrp);  
     dxl= dmatrix(1,6,1,ngrp);  
     dij= dmatrix(1,6,1,ngrp);  
      dx= dmatrix(1,6,1,ngrp);
      dy= dmatrix(1,6,1,ngrp);
      dz= dmatrix(1,6,1,ngrp);

      dxt= dmatrix(1,6,1,ngrp);
      dyt= dmatrix(1,6,1,ngrp);
      dzt= dmatrix(1,6,1,ngrp);

      dxn= dmatrix(1,6,1,ngrp);
      dyn= dmatrix(1,6,1,ngrp);
      dzn= dmatrix(1,6,1,ngrp);
   rmass1= dvector(1,ngrp);
   rmass2= dvector(1,ngrp);
   rmass3= dvector(1,ngrp);
   rmass4= dvector(1,ngrp);
        x= dmatrix(1,4,1,ngrp);
        y= dmatrix(1,4,1,ngrp);
        z= dmatrix(1,4,1,ngrp);
       xo= dmatrix(1,4,1,ngrp);
       yo= dmatrix(1,4,1,ngrp);
       zo= dmatrix(1,4,1,ngrp);
      p11= dvector(1,ngrp);
      p12= dvector(1,ngrp);
      p13= dvector(1,ngrp);
      p22= dvector(1,ngrp);
      p23= dvector(1,ngrp);
      p33= dvector(1,ngrp);
   tamat= dvector(1,36);
 
    ind1= (int *)calloc((ngrp+1),sizeof(int));
    ind2= (int *)calloc((ngrp+1),sizeof(int));
    ind3= (int *)calloc((ngrp+1),sizeof(int));
    ind4= (int *)calloc((ngrp+1),sizeof(int));
  }/*endif*/

/*=======================================================================*/
/* Malloc up some vectors and matrices */
 na = NCON_46;
/* AA */

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
  for(igrp=1;igrp <= ngrp; igrp++) {
   grp_bond_con_al_46[1][(igrp+igrp_off)] = 0.0;
   grp_bond_con_al_46[2][(igrp+igrp_off)] = 0.0;
   grp_bond_con_al_46[3][(igrp+igrp_off)] = 0.0;
   grp_bond_con_al_46[4][(igrp+igrp_off)] = 0.0;
   grp_bond_con_al_46[5][(igrp+igrp_off)] = 0.0;
   grp_bond_con_al_46[6][(igrp+igrp_off)] = 0.0;
  }/*endif*/
 }/*endif*/

 for(igrp=1;igrp <= ngrp; igrp++) {

  ind1[igrp] = grp_bond_con_j1_46[(igrp+igrp_off)];
  ind2[igrp] = grp_bond_con_j2_46[(igrp+igrp_off)];
  ind3[igrp] = grp_bond_con_j3_46[(igrp+igrp_off)];
  ind4[igrp] = grp_bond_con_j4_46[(igrp+igrp_off)];
}

 for(igrp=1;igrp <= ngrp; igrp++) {
     ktemp= ind1[igrp];
    x[1][igrp] = clatoms_x[ktemp];
    y[1][igrp] = clatoms_y[ktemp];
    z[1][igrp] = clatoms_z[ktemp];
   rmass1[igrp]  = 1.0/clatoms_mass[ktemp];
 }

 for(igrp=1;igrp <= ngrp; igrp++) {
     ktemp= ind2[igrp];
    x[2][igrp] = clatoms_x[ktemp];
    y[2][igrp] = clatoms_y[ktemp];
    z[2][igrp] = clatoms_z[ktemp];
   rmass2[igrp]  = 1.0/clatoms_mass[ktemp];
 }

 for(igrp=1;igrp <= ngrp; igrp++) {
     ktemp= ind3[igrp];
    x[3][igrp] = clatoms_x[ktemp];
    y[3][igrp] = clatoms_y[ktemp];
    z[3][igrp] = clatoms_z[ktemp];
   rmass3[igrp]  = 1.0/clatoms_mass[ktemp];
 }

 for(igrp=1;igrp <= ngrp; igrp++) {
     ktemp= ind4[igrp];
    x[4][igrp] = clatoms_x[ktemp];
    y[4][igrp] = clatoms_y[ktemp];
    z[4][igrp] = clatoms_z[ktemp];
   rmass4[igrp]  = 1.0/clatoms_mass[ktemp];
 }

 for(igrp=1;igrp <= ngrp; igrp++) {
     ktemp= ind1[igrp];
    xo[1][igrp] = clatoms_xold[ktemp];
    yo[1][igrp] = clatoms_yold[ktemp];
    zo[1][igrp] = clatoms_zold[ktemp];
 }

 for(igrp=1;igrp <= ngrp; igrp++) {
     ktemp= ind2[igrp];
    xo[2][igrp] = clatoms_xold[ktemp];
    yo[2][igrp] = clatoms_yold[ktemp];
    zo[2][igrp] = clatoms_zold[ktemp];
 }

 for(igrp=1;igrp <= ngrp; igrp++) {
     ktemp= ind3[igrp];
    xo[3][igrp] = clatoms_xold[ktemp];
    yo[3][igrp] = clatoms_yold[ktemp];
    zo[3][igrp] = clatoms_zold[ktemp];
 }

 for(igrp=1;igrp <= ngrp; igrp++) {
     ktemp= ind4[igrp];
    xo[4][igrp] = clatoms_xold[ktemp];
    yo[4][igrp] = clatoms_yold[ktemp];
    zo[4][igrp] = clatoms_zold[ktemp];
 }



/* ============================================================================= */
/* Gather the equilibrium bond lengths */
 for(igrp=1;igrp <= ngrp; igrp++) {
  jtyp = grp_bond_con_jtyp_46[(igrp+igrp_off)];

   dij[1][igrp] = grp_bond_con_eq_46[1][jtyp];
   dij[2][igrp] = grp_bond_con_eq_46[2][jtyp];
   dij[3][igrp] = grp_bond_con_eq_46[3][jtyp];
   dij[4][igrp] = grp_bond_con_eq_46[4][jtyp];
   dij[5][igrp] = grp_bond_con_eq_46[5][jtyp];
   dij[6][igrp] = grp_bond_con_eq_46[6][jtyp];
 }/*end for*/

/* ============================================================================= */
/* Calculate the recip mass tensor and bond distances                            */

 for(igrp=1;igrp <= ngrp; igrp++) {
    rms1 = rmass1[igrp];
    rms2 = rmass2[igrp];
    rms3 = rmass3[igrp];
    rms4 = rmass4[igrp];

   rmassm[1][1][igrp] = -(rms1+rms2);
   rmassm[1][2][igrp] = rmassm[1][3][igrp] = -rms1;
   rmassm[1][4][igrp] = rmassm[1][5][igrp] = rms2; rmassm[1][6][igrp] = 0.0;
 
   rmassm[2][1][igrp] = -rms1; rmassm[2][2][igrp] = -(rms1+rms3); 
                          rmassm[2][3][igrp] = -rms1;
   rmassm[2][4][igrp] = -rms3; rmassm[2][5][igrp] = 0.0; 
   rmassm[2][6][igrp] = rms3;

   rmassm[3][1][igrp] = rmassm[3][2][igrp] = -rms1;
   rmassm[3][3][igrp] = -(rms1+rms4);
   rmassm[3][4][igrp] = 0.0; rmassm[3][5][igrp] = rmassm[3][6][igrp] = -rms4;

   rmassm[4][1][igrp] = rms2;    rmassm[4][2][igrp] = -rms3;
   rmassm[4][3][igrp] = 0.0;
   rmassm[4][4][igrp] = -(rms2+rms3); rmassm[4][5][igrp] = -rms2; 
                                   rmassm[4][6][igrp] = rms3;

   rmassm[5][1][igrp] = rms2;  rmassm[5][2][igrp] = 0.0; 
   rmassm[5][3][igrp] = -rms4;
   rmassm[5][4][igrp] = -rms2; rmassm[5][5][igrp] = -(rms2+rms4); 
                          rmassm[5][6][igrp] = -rms4;

   rmassm[6][1][igrp] = 0.0;    rmassm[6][2][igrp] = rms3;
   rmassm[6][3][igrp] = -rms4;
   rmassm[6][4][igrp] = rms3; rmassm[6][5][igrp] = -rms4; 
                         rmassm[6][6][igrp] = -(rms3+rms4);
 }

/* Compute difference vectors : Old distances scaled*/
 for(igrp=1;igrp <= ngrp; igrp++) {
    dxt[1][igrp] = x[1][igrp]-x[2][igrp];
    dxt[2][igrp] = x[1][igrp]-x[3][igrp];
    dxt[3][igrp] = x[1][igrp]-x[4][igrp];
    dxt[4][igrp] = x[2][igrp]-x[3][igrp];
    dxt[5][igrp] = x[2][igrp]-x[4][igrp];
    dxt[6][igrp] = x[3][igrp]-x[4][igrp];
 
    dyt[1][igrp] = y[1][igrp]-y[2][igrp];
    dyt[2][igrp] = y[1][igrp]-y[3][igrp];
    dyt[3][igrp] = y[1][igrp]-y[4][igrp];
    dyt[4][igrp] = y[2][igrp]-y[3][igrp];
    dyt[5][igrp] = y[2][igrp]-y[4][igrp];
    dyt[6][igrp] = y[3][igrp]-y[4][igrp];
 
    dzt[1][igrp] = z[1][igrp]-z[2][igrp];
    dzt[2][igrp] = z[1][igrp]-z[3][igrp];
    dzt[3][igrp] = z[1][igrp]-z[4][igrp];
    dzt[4][igrp] = z[2][igrp]-z[3][igrp];
    dzt[5][igrp] = z[2][igrp]-z[4][igrp];
    dzt[6][igrp] = z[3][igrp]-z[4][igrp];

    dx[1][igrp] = (xo[1][igrp]-xo[2][igrp])*baro_roll_scv;
    dx[2][igrp] = (xo[1][igrp]-xo[3][igrp])*baro_roll_scv;
    dx[3][igrp] = (xo[1][igrp]-xo[4][igrp])*baro_roll_scv;
    dx[4][igrp] = (xo[2][igrp]-xo[3][igrp])*baro_roll_scv;
    dx[5][igrp] = (xo[2][igrp]-xo[4][igrp])*baro_roll_scv;
    dx[6][igrp] = (xo[3][igrp]-xo[4][igrp])*baro_roll_scv;

    dy[1][igrp] = (yo[1][igrp]-yo[2][igrp])*baro_roll_scv;
    dy[2][igrp] = (yo[1][igrp]-yo[3][igrp])*baro_roll_scv;
    dy[3][igrp] = (yo[1][igrp]-yo[4][igrp])*baro_roll_scv;
    dy[4][igrp] = (yo[2][igrp]-yo[3][igrp])*baro_roll_scv;
    dy[5][igrp] = (yo[2][igrp]-yo[4][igrp])*baro_roll_scv;
    dy[6][igrp] = (yo[3][igrp]-yo[4][igrp])*baro_roll_scv;

    dz[1][igrp] = (zo[1][igrp]-zo[2][igrp])*baro_roll_scv;
    dz[2][igrp] = (zo[1][igrp]-zo[3][igrp])*baro_roll_scv;
    dz[3][igrp] = (zo[1][igrp]-zo[4][igrp])*baro_roll_scv;
    dz[4][igrp] = (zo[2][igrp]-zo[3][igrp])*baro_roll_scv;
    dz[5][igrp] = (zo[2][igrp]-zo[4][igrp])*baro_roll_scv;
    dz[6][igrp] = (zo[3][igrp]-zo[4][igrp])*baro_roll_scv;

} /* end loop over groups */
/* =========================================================================== */
/* Get initial guess for lambda                                                */

  for(igrp=1;igrp <= ngrp; igrp++) {
   avec[1][igrp] = dij[1][igrp]*dij[1][igrp] - (dxt[1][igrp]*dxt[1][igrp]
                 + dyt[1][igrp]*dyt[1][igrp] +  dzt[1][igrp]*dzt[1][igrp]);
   avec[2][igrp] = dij[2][igrp]*dij[2][igrp] - (dxt[2][igrp]*dxt[2][igrp]
                 + dyt[2][igrp]*dyt[2][igrp] +  dzt[2][igrp]*dzt[2][igrp]);
   avec[3][igrp] = dij[3][igrp]*dij[3][igrp] - (dxt[3][igrp]*dxt[3][igrp]
                 + dyt[3][igrp]*dyt[3][igrp] +  dzt[3][igrp]*dzt[3][igrp]);
   avec[4][igrp] = dij[4][igrp]*dij[4][igrp] - (dxt[4][igrp]*dxt[4][igrp]
                 + dyt[4][igrp]*dyt[4][igrp] +  dzt[4][igrp]*dzt[4][igrp]);
   avec[5][igrp] = dij[5][igrp]*dij[5][igrp] - (dxt[5][igrp]*dxt[5][igrp]
                 + dyt[5][igrp]*dyt[5][igrp] +  dzt[5][igrp]*dzt[5][igrp]);
   avec[6][igrp] = dij[6][igrp]*dij[6][igrp] - (dxt[6][igrp]*dxt[6][igrp]
                 + dyt[6][igrp]*dyt[6][igrp] +  dzt[6][igrp]*dzt[6][igrp]);
  }/* endfor */
 if(ifirst == 2 || ifirst == 0){
    iii = 0;
  for(i=1; i <= NCON_46; i++){
   for(j=1; j <= NCON_46; j++){
      iii++;
     for(igrp=1;igrp <= ngrp; igrp++) {
         amat[iii][igrp] = 2.0*rmassm[i][j][igrp]*
               (dxt[i][igrp]*dx[j][igrp] + dyt[i][igrp]*dy[j][igrp]
              + dzt[i][igrp]*dz[j][igrp]);
    }/*endfor igrp*/
   }/*endfor j*/
  }/*endfor i*/

  for(igrp=1;igrp <= ngrp; igrp++) {

   xlam[1][igrp] = avec[1][igrp];
   xlam[2][igrp] = avec[2][igrp];
   xlam[3][igrp] = avec[3][igrp];
   xlam[4][igrp] = avec[4][igrp];
   xlam[5][igrp] = avec[5][igrp];
   xlam[6][igrp] = avec[6][igrp];

      txlam[1]=xlam[1][igrp];
      txlam[2]=xlam[2][igrp];
      txlam[3]=xlam[3][igrp];
      txlam[4]=xlam[4][igrp];
      txlam[5]=xlam[5][igrp];
      txlam[6]=xlam[6][igrp];

   for(i=1; i<=36 ; i++) {
      tamat[i]= amat[i][igrp];
   }

/* Solve linear system A xlam = avec */
   ipvt[1]=0; ipvt[2]=0; ipvt[3]=0;
   ipvt[4]=0; ipvt[5]=0; ipvt[6]=0;

#ifdef IBM_ESSL
  dgef(&(tamat[1]),&na,&na,&(ipvt[1]));
#else
  DGEFA(&(tamat[1]),&na,&na,&(ipvt[1]),&info);
#endif
  job = 1;
#ifdef IBM_ESSL
  job = 1;  /*changed from 0 */
  dges(&(tamat[1]),&na,&na,&(ipvt[1]),&(txlam[1]),&job);
#else
  DGESL(&(tamat[1]),&na,&na,&(ipvt[1]),&(txlam[1]),&job);
#endif
     xlam[1][igrp] = txlam[1];
     xlam[2][igrp] = txlam[2];
     xlam[3][igrp] = txlam[3];
     xlam[4][igrp] = txlam[4];
     xlam[5][igrp] = txlam[5];
     xlam[6][igrp] = txlam[6];

 } /* end loop over groups */
 } else {
  for(igrp=1;igrp <= ngrp; igrp++) {
   xlam[1][igrp] = grp_bond_con_al_46[1][(igrp+igrp_off)];
   xlam[2][igrp] = grp_bond_con_al_46[2][(igrp+igrp_off)];
   xlam[3][igrp] = grp_bond_con_al_46[3][(igrp+igrp_off)];
   xlam[4][igrp] = grp_bond_con_al_46[4][(igrp+igrp_off)];
   xlam[5][igrp] = grp_bond_con_al_46[5][(igrp+igrp_off)];
   xlam[6][igrp] = grp_bond_con_al_46[6][(igrp+igrp_off)];

   grp_bond_con_al_46[1][(igrp+igrp_off)] = 0.0;
   grp_bond_con_al_46[2][(igrp+igrp_off)] = 0.0;
   grp_bond_con_al_46[3][(igrp+igrp_off)] = 0.0;
   grp_bond_con_al_46[4][(igrp+igrp_off)] = 0.0;
   grp_bond_con_al_46[5][(igrp+igrp_off)] = 0.0;
   grp_bond_con_al_46[6][(igrp+igrp_off)] = 0.0;
  } /* end for */
 }

/* Iterative loop to convergence */

 if(ngrp > 0){
  dmax = 1.0;
  iter = 0;
  do {
   ++iter;
   if(iter > grp_bond_con->max_iter) {
    printf("$$$$$$$$$$$$$$$$$$$$_WARNING_$$$$$$$$$$$$$$$$$$$$\n");
    printf("Group constraint Shake not converged after %d iterations.\n",
            grp_bond_con->max_iter);
    printf("The present tolerance is %g \n",dmax);
    printf("The desired tolerance is %g \n",grp_bond_con->tol);
    printf("$$$$$$$$$$$$$$$$$$$$_WARNING_$$$$$$$$$$$$$$$$$$$$\n");
    fflush(stdout);
    break;
   }/*endif*/
  for(igrp=1;igrp <= ngrp; igrp++) {
/* Set up guess of difference vectors */
    dxn[1][igrp] = 2.0*dxt[1][igrp]; dxn[2][igrp] = 2.0*dxt[2][igrp];
    dxn[3][igrp] = 2.0*dxt[3][igrp]; dxn[4][igrp] = 2.0*dxt[4][igrp];
    dxn[5][igrp] = 2.0*dxt[5][igrp]; dxn[6][igrp] = 2.0*dxt[6][igrp];
 
    dyn[1][igrp] = 2.0*dyt[1][igrp]; dyn[2][igrp] = 2.0*dyt[2][igrp];
    dyn[3][igrp] = 2.0*dyt[3][igrp]; dyn[4][igrp] = 2.0*dyt[4][igrp];
    dyn[5][igrp] = 2.0*dyt[5][igrp]; dyn[6][igrp] = 2.0*dyt[6][igrp];

    dzn[1][igrp] = 2.0*dzt[1][igrp]; dzn[2][igrp] = 2.0*dzt[2][igrp];
    dzn[3][igrp] = 2.0*dzt[3][igrp]; dzn[4][igrp] = 2.0*dzt[4][igrp];
    dzn[5][igrp] = 2.0*dzt[5][igrp]; dzn[6][igrp] = 2.0*dzt[6][igrp];

   }/*endfor*/
   for(i=1; i <= NCON_46; i++) {
    for(j=1; j <= NCON_46; j++) {
      for(igrp=1;igrp <= ngrp; igrp++) {
        dxn[i][igrp] += rmassm[i][j][igrp]*xlam[j][igrp]*dx[j][igrp];
        dyn[i][igrp] += rmassm[i][j][igrp]*xlam[j][igrp]*dy[j][igrp];
        dzn[i][igrp] += rmassm[i][j][igrp]*xlam[j][igrp]*dz[j][igrp];
      }/*endfor*/
    }/*endfor*/
   }/*endfor*/

/* Construct A-matrix */

   iii = 0;
   for(i=1; i <= NCON_46; i++) {
    for(j=1; j <= NCON_46; j++) {
      iii++;
      for(igrp=1;igrp <= ngrp; igrp++) {
          amat[iii][igrp] = rmassm[i][j][igrp]*
                (dxn[i][igrp]*dx[j][igrp] + dyn[i][igrp]*dy[j][igrp]
               + dzn[i][igrp]*dz[j][igrp]);
      }/*endfor*/
    }/*endfor*/
   }/*endfor*/

  for(igrp=1;igrp <= ngrp; igrp++) {
    xl0[1] = xlam[1][igrp];
    xl0[2] = xlam[2][igrp];
    xl0[3] = xlam[3][igrp];
    xl0[4] = xlam[4][igrp];
    xl0[5] = xlam[5][igrp];
    xl0[6] = xlam[6][igrp];

    xlam[1][igrp] = avec[1][igrp];
    xlam[2][igrp] = avec[2][igrp];
    xlam[3][igrp] = avec[3][igrp];
    xlam[4][igrp] = avec[4][igrp];
    xlam[5][igrp] = avec[5][igrp];
    xlam[6][igrp] = avec[6][igrp];
 
     txlam[1]= xlam[1][igrp];
     txlam[2]= xlam[2][igrp];
     txlam[3]= xlam[3][igrp];
     txlam[4]= xlam[4][igrp];
     txlam[5]= xlam[5][igrp];
     txlam[6]= xlam[6][igrp];

    for(i=1;i<=36; i++ ) {
      tamat[i]= amat[i][igrp];
    }

/* Solve linear system A xlam = avec */
   ipvt[1] = 0; ipvt[2] = 0; ipvt[3] = 0;
   ipvt[4] = 0; ipvt[5] = 0; ipvt[6] = 0;

#ifdef IBM_ESSL
  dgef(&(tamat[1]),&na,&na,&(ipvt[1]));
#else
   DGEFA(&(tamat[1]),&na,&na,&(ipvt[1]),&info);
#endif
   job = 1;
#ifdef IBM_ESSL
  job = 1;  /*changed from 0 */
  dges(&(tamat[1]),&na,&na,&(ipvt[1]),&(txlam[1]),&job);
#else
   DGESL(&(tamat[1]),&na,&na,&(ipvt[1]),&(txlam[1]),&job);
#endif
    xlam[1][igrp] = txlam[1];
    xlam[2][igrp] = txlam[2];
    xlam[3][igrp] = txlam[3];
    xlam[4][igrp] = txlam[4];
    xlam[5][igrp] = txlam[5];
    xlam[6][igrp] = txlam[6];

    dxl[1][igrp] = fabs(xlam[1][igrp]-xl0[1]);
    dxl[2][igrp] = fabs(xlam[2][igrp]-xl0[2]);
    dxl[3][igrp] = fabs(xlam[3][igrp]-xl0[3]);
    dxl[4][igrp] = fabs(xlam[4][igrp]-xl0[4]);
    dxl[5][igrp] = fabs(xlam[5][igrp]-xl0[5]);
    dxl[6][igrp] = fabs(xlam[6][igrp]-xl0[6]);
 }/*endfor*/
/* Convergence criteria */
      dlmax[1]= dxl[1][1];
      dlmax[2]= dxl[2][1];
      dlmax[3]= dxl[3][1];
      dlmax[4]= dxl[4][1];
      dlmax[5]= dxl[5][1];
      dlmax[6]= dxl[6][1];
  for(igrp=2;igrp <= ngrp; igrp++) {
    dlmax[1]= (dlmax[1] > dxl[1][igrp] ? dlmax[1]: dxl[1][igrp]);
    dlmax[2]= (dlmax[2] > dxl[2][igrp] ? dlmax[2]: dxl[2][igrp]);
    dlmax[3]= (dlmax[3] > dxl[3][igrp] ? dlmax[3]: dxl[3][igrp]);
    dlmax[4]= (dlmax[4] > dxl[4][igrp] ? dlmax[4]: dxl[4][igrp]);
    dlmax[5]= (dlmax[5] > dxl[5][igrp] ? dlmax[5]: dxl[5][igrp]);
    dlmax[6]= (dlmax[6] > dxl[6][igrp] ? dlmax[6]: dxl[6][igrp]);
   }/*end loop over groups */
    dmax=dlmax[1];
   for(i=2;i <= NCON_46; i++) {
    dmax = (dmax > dlmax[i] ? dmax : dlmax[i]);
   }/*endfor*/
  } while(dmax > grp_bond_con->tol);
  *aiter += (double) iter;

 }/*endif for ngrp > 0*/

/* Position update */

#ifndef NO_PRAGMA
#pragma IVDEP
#endif
 for(igrp=1;igrp <= ngrp; igrp++) {
   double xlam1,xlam2,xlam3,xlam4,xlam5,xlam6;
   double dx1,dx2,dx3,dx4,dx5,dx6;
   double dy1,dy2,dy3,dy4,dy5,dy6;
   double dz1,dz2,dz3,dz4,dz5,dz6;


    ktemp1 =ind1[igrp];
    ktemp2 =ind2[igrp];
    ktemp3 =ind3[igrp];
    ktemp4 =ind4[igrp];

    xlam1= xlam[1][igrp];
    xlam2= xlam[2][igrp];
    xlam3= xlam[3][igrp];
    xlam4= xlam[4][igrp];
    xlam5= xlam[5][igrp];
    xlam6= xlam[6][igrp];

    rms1 = rmass1[igrp];
    rms2 = rmass2[igrp];
    rms3 = rmass3[igrp];
    rms4 = rmass4[igrp];


  x[1][igrp] = clatoms_x[ktemp1]; y[1][igrp] = clatoms_y[ktemp1]; 
  z[1][igrp] = clatoms_z[ktemp1];
  x[2][igrp] = clatoms_x[ktemp2]; y[2][igrp] = clatoms_y[ktemp2]; 
  z[2][igrp] = clatoms_z[ktemp2];
  x[3][igrp] = clatoms_x[ktemp3]; y[3][igrp] = clatoms_y[ktemp3]; 
  z[3][igrp] = clatoms_z[ktemp3];
  x[4][igrp] = clatoms_x[ktemp4]; y[4][igrp] = clatoms_y[ktemp4]; 
  z[4][igrp] = clatoms_z[ktemp4];

  xo[1][igrp] = clatoms_xold[ktemp1]; yo[1][igrp] = clatoms_yold[ktemp1]; 
  zo[1][igrp] = clatoms_zold[ktemp1];
  xo[2][igrp] = clatoms_xold[ktemp2]; yo[2][igrp] = clatoms_yold[ktemp2]; 
  zo[2][igrp] = clatoms_zold[ktemp2];
  xo[3][igrp] = clatoms_xold[ktemp3]; yo[3][igrp] = clatoms_yold[ktemp3]; 
  zo[3][igrp] = clatoms_zold[ktemp3];
  xo[4][igrp] = clatoms_xold[ktemp4]; yo[4][igrp] = clatoms_yold[ktemp4]; 
  zo[4][igrp] = clatoms_zold[ktemp4];

    dx[1][igrp] = dx1 = (xo[1][igrp]-xo[2][igrp]);
    dx[2][igrp] = dx2 = (xo[1][igrp]-xo[3][igrp]);
    dx[3][igrp] = dx3 = (xo[1][igrp]-xo[4][igrp]);
    dx[4][igrp] = dx4 = (xo[2][igrp]-xo[3][igrp]);
    dx[5][igrp] = dx5 = (xo[2][igrp]-xo[4][igrp]);
    dx[6][igrp] = dx6 = (xo[3][igrp]-xo[4][igrp]);

    dy[1][igrp] = dy1 = (yo[1][igrp]-yo[2][igrp]);
    dy[2][igrp] = dy2 = (yo[1][igrp]-yo[3][igrp]);
    dy[3][igrp] = dy3 = (yo[1][igrp]-yo[4][igrp]);
    dy[4][igrp] = dy4 = (yo[2][igrp]-yo[3][igrp]);
    dy[5][igrp] = dy5 = (yo[2][igrp]-yo[4][igrp]);
    dy[6][igrp] = dy6 = (yo[3][igrp]-yo[4][igrp]);

    dz[1][igrp] = dz1 = (zo[1][igrp]-zo[2][igrp]);
    dz[2][igrp] = dz2 = (zo[1][igrp]-zo[3][igrp]);
    dz[3][igrp] = dz3 = (zo[1][igrp]-zo[4][igrp]);
    dz[4][igrp] = dz4 = (zo[2][igrp]-zo[3][igrp]);
    dz[5][igrp] = dz5 = (zo[2][igrp]-zo[4][igrp]);
    dz[6][igrp] = dz6 = (zo[3][igrp]-zo[4][igrp]);


  clatoms_x[ktemp1] -=  ( xlam1*dx1 + xlam2*dx2 + xlam3*dx3)*rms1 *baro_roll_scv;
  clatoms_y[ktemp1] -=  ( xlam1*dy1 + xlam2*dy2 + xlam3*dy3)*rms1 *baro_roll_scv;
  clatoms_z[ktemp1] -=  ( xlam1*dz1 + xlam2*dz2 + xlam3*dz3)*rms1 *baro_roll_scv;

  clatoms_x[ktemp2] -=  (-xlam1*dx1 + xlam4*dx4 + xlam5*dx5)*rms2 *baro_roll_scv;
  clatoms_y[ktemp2] -=  (-xlam1*dy1 + xlam4*dy4 + xlam5*dy5)*rms2 *baro_roll_scv;
  clatoms_z[ktemp2] -=  (-xlam1*dz1 + xlam4*dz4 + xlam5*dz5)*rms2 *baro_roll_scv;


  clatoms_x[ktemp3] -=  (-xlam2*dx2 - xlam4*dx4 + xlam6*dx6)*rms3 *baro_roll_scv;
  clatoms_y[ktemp3] -=  (-xlam2*dy2 - xlam4*dy4 + xlam6*dy6)*rms3 *baro_roll_scv;
  clatoms_z[ktemp3] -=  (-xlam2*dz2 - xlam4*dz4 + xlam6*dz6)*rms3 *baro_roll_scv;

  clatoms_x[ktemp4] -=  (-xlam3*dx3 - xlam5*dx5 - xlam6*dx6)*rms4 *baro_roll_scv;
  clatoms_y[ktemp4] -=  (-xlam3*dy3 - xlam5*dy5 - xlam6*dy6)*rms4 *baro_roll_scv;
  clatoms_z[ktemp4] -=  (-xlam3*dz3 - xlam5*dz5 - xlam6*dz6)*rms4 *baro_roll_scv;

/* Velocity update */

  clatoms_vx[ktemp1] -=  ( xlam1*dx1 + xlam2*dx2 + xlam3*dx3) *rms1/dt;
  clatoms_vy[ktemp1] -=  ( xlam1*dy1 + xlam2*dy2 + xlam3*dy3) *rms1/dt;
  clatoms_vz[ktemp1] -=  ( xlam1*dz1 + xlam2*dz2 + xlam3*dz3) *rms1/dt;

  clatoms_vx[ktemp2] -=  (-xlam1*dx1 + xlam4*dx4 + xlam5*dx5) *rms2/dt;
  clatoms_vy[ktemp2] -=  (-xlam1*dy1 + xlam4*dy4 + xlam5*dy5) *rms2/dt;
  clatoms_vz[ktemp2] -=  (-xlam1*dz1 + xlam4*dz4 + xlam5*dz5) *rms2/dt;

  clatoms_vx[ktemp3] -=  (-xlam2*dx2 - xlam4*dx4 + xlam6*dx6) *rms3/dt;
  clatoms_vy[ktemp3] -=  (-xlam2*dy2 - xlam4*dy4 + xlam6*dy6) *rms3/dt;
  clatoms_vz[ktemp3] -=  (-xlam2*dz2 - xlam4*dz4 + xlam6*dz6) *rms3/dt;

  clatoms_vx[ktemp4] -=  (-xlam3*dx3 - xlam5*dx5 - xlam6*dx6) *rms4/dt;
  clatoms_vy[ktemp4] -=  (-xlam3*dy3 - xlam5*dy5 - xlam6*dy6) *rms4/dt;
  clatoms_vz[ktemp4] -=  (-xlam3*dz3 - xlam5*dz5 - xlam6*dz6) *rms4/dt;


/* Pressure tensor update */
/* Compute difference vectors: use unscaled old distances */

     p11[igrp]= xlam1*dx1*dx1 + xlam2*dx2*dx2 + xlam3*dx3*dx3
              + xlam4*dx4*dx4 + xlam5*dx5*dx5 + xlam6*dx6*dx6;

     p22[igrp]= xlam1*dy1*dy1 + xlam2*dy2*dy2 + xlam3*dy3*dy3
              + xlam4*dy4*dy4 + xlam5*dy5*dy5 + xlam6*dy6*dy6;

     p33[igrp]= xlam1*dz1*dz1 + xlam2*dz2*dz2 + xlam3*dz3*dz3
              + xlam4*dz4*dz4 + xlam5*dz5*dz5 + xlam6*dz6*dz6;

     p12[igrp]= xlam1*dx1*dy1 + xlam2*dx2*dy2 + xlam3*dx3*dy3
              + xlam4*dx4*dy4 + xlam5*dx5*dy5 + xlam6*dx6*dy6;

     p13[igrp]= xlam1*dx1*dz1 + xlam2*dx2*dz2 + xlam3*dx3*dz3
              + xlam4*dx4*dz4 + xlam5*dx5*dz5 + xlam6*dx6*dz6;

     p23[igrp]= xlam1*dy1*dz1 + xlam2*dy2*dz2 + xlam3*dy3*dz3
              + xlam4*dy4*dz4 + xlam5*dy5*dz5 + xlam6*dy6*dz6;

}/*end for*/

#ifndef NO_PRAGMA
#pragma IVDEP
#endif
 for(igrp=1;igrp <= ngrp; igrp++) {
    ptens_pvten_tmp[1] -= (p11[igrp]*pnorm); 
    ptens_pvten_tmp[2] -= (p12[igrp]*pnorm);
    ptens_pvten_tmp[3] -= (p13[igrp]*pnorm);
    ptens_pvten_tmp[4] -= (p12[igrp]*pnorm);
    ptens_pvten_tmp[5] -= (p22[igrp]*pnorm);
    ptens_pvten_tmp[6] -= (p23[igrp]*pnorm);
    ptens_pvten_tmp[7] -= (p13[igrp]*pnorm);
    ptens_pvten_tmp[8] -= (p23[igrp]*pnorm);
    ptens_pvten_tmp[9] -= (p33[igrp]*pnorm);
 } /* end for  */


/* Save multiplier */

 for(igrp=1;igrp <= ngrp; igrp++) {
    grp_bond_con_al_46[1][(igrp+igrp_off)] += xlam[1][igrp];
    grp_bond_con_al_46[2][(igrp+igrp_off)] += xlam[2][igrp];
    grp_bond_con_al_46[3][(igrp+igrp_off)] += xlam[3][igrp];
    grp_bond_con_al_46[4][(igrp+igrp_off)] += xlam[4][igrp];
    grp_bond_con_al_46[5][(igrp+igrp_off)] += xlam[5][igrp];
    grp_bond_con_al_46[6][(igrp+igrp_off)] += xlam[6][igrp];

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
  ftemp   = (ptens_pvten_tmp[1]+ptens_pvten_tmp[5]+ptens_pvten_tmp[9]);
  baro->f_lnv_p += ftemp;
  baro->v_lnv   += 0.5*ftemp*(baro->roll_scg)*dt/(baro->mass_lnv);
 }
/* free locally assigned memory */
 if(ngrp > 0){
    free_dvector(dlmax,1,6);
    free_dvector(txlam,1,6);
    free_d3tensor(rmassm,1,6,1,6,1,ngrp);
    free_dmatrix(amat,1,36,1,ngrp); 
    free_dmatrix(xlam,1,6,1,ngrp);
    free_dmatrix(avec,1,6,1,ngrp);
    free_dmatrix(dxl,1,6,1,ngrp);
    free_dmatrix(dij,1,6,1,ngrp);

    free_dmatrix(dx,1,6,1,ngrp);
    free_dmatrix(dy,1,6,1,ngrp);
    free_dmatrix(dz,1,6,1,ngrp);

    free_dmatrix(dxt,1,6,1,ngrp);
    free_dmatrix(dyt,1,6,1,ngrp);
    free_dmatrix(dzt,1,6,1,ngrp);

    free_dmatrix(dxn,1,6,1,ngrp);
    free_dmatrix(dyn,1,6,1,ngrp);
    free_dmatrix(dzn,1,6,1,ngrp);

    free_dvector(rmass1,1,ngrp);
    free_dvector(rmass2,1,ngrp);
    free_dvector(rmass3,1,ngrp);
    free_dvector(rmass4,1,ngrp);

    free_dmatrix(x,1,ngrp,1,36); 
    free_dmatrix(y,1,ngrp,1,36); 
    free_dmatrix(z,1,ngrp,1,36); 

    free_dmatrix(xo,1,ngrp,1,36); 
    free_dmatrix(yo,1,ngrp,1,36); 
    free_dmatrix(zo,1,ngrp,1,36); 

    free_dvector(p11,1,ngrp);
    free_dvector(p12,1,ngrp);
    free_dvector(p13,1,ngrp);
    free_dvector(p22,1,ngrp);
    free_dvector(p23,1,ngrp);
    free_dvector(p33,1,ngrp);

    free_dvector(tamat,1,36); 

    free(ind1); free(ind2); free(ind3); free(ind4);
 }/*endif*/

/*==========================================================================*/
} /* end routine */
/*==========================================================================*/

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void rattle_46_rolli(GRP_BOND_CON *grp_bond_con,
              CLATOMS_INFO *clatoms_info,CLATOMS_POS *clatoms_pos,
              PTENS *ptens,double dt,BARO *baro,int ifirst,
              CLASS_COMM_FORC_PKG *class_comm_forc_pkg)

/*==========================================================================*/
/*        Begin Routine                                                     */
{/*Begin Routine*/
/*=======================================================================*/
/*         Local Variable declarations                                   */
  
#include "../typ_defs/typ_mask.h"

 double rmass_1,rmass_2,rmass_3,rmass_4;
 double avec[NCON_46+1];
 double pnorm;

 double roll_sci,dlam;
 double f_lnv_inc;
 int i,j,k,iii;
 int igrp,*ind1,*ind2,*ind3,*ind4,jtyp;
 int ktemp,ktemp1,ktemp2,ktemp3,ktemp4;
 int na,job,info,ipvt[NCON_46+1];       /* For dgefa and dgesl */

  double *rmass1,*rmass2,*rmass3,*rmass4;
  double **x,**y,**z;
  double **vx,**vy,**vz;
  double *p11,*p22,*p33,*p12,*p13,*p23;
  double ***rmassm;
  double **dvx,**dvy,**dvz;
  double **dx,**dy,**dz;
  double **amat,**xlam;
  double *txlam,*tamat;

/* Local pointers */

  double *clatoms_mass         = clatoms_info->mass;
  double *clatoms_x            = clatoms_pos->x;
  double *clatoms_y            = clatoms_pos->y;
  double *clatoms_z            = clatoms_pos->z;
  double *clatoms_vx           = clatoms_pos->vx;
  double *clatoms_vy           = clatoms_pos->vy;
  double *clatoms_vz           = clatoms_pos->vz;
  int *grp_bond_con_j1_46      = grp_bond_con->j1_46;
  int *grp_bond_con_j2_46      = grp_bond_con->j2_46;
  int *grp_bond_con_j3_46      = grp_bond_con->j3_46;
  int *grp_bond_con_j4_46      = grp_bond_con->j4_46;
  int *grp_bond_con_jtyp_46    = grp_bond_con->jtyp_46;
  double *ptens_pvten_inc      = ptens->pvten_inc;
  double *ptens_pvten_tmp      = ptens->pvten_tmp;
  double *ptens_pvten_tmp2      = ptens->pvten_tmp_res;
  double *clatoms_roll_sc = clatoms_info->roll_sc;
  double baro_v_lnv_g = baro->v_lnv_g;

  int ngrp,irem,igrp_off;
  int ngrp_tot                 = grp_bond_con->num_46;
  int np_forc                  = class_comm_forc_pkg->num_proc;
  int myid_forc                = class_comm_forc_pkg->myid;
  MPI_Comm comm_forc           = class_comm_forc_pkg->comm;

/*=======================================================================*/
 
  ngrp = (ngrp_tot);
  igrp_off = 0;

/*=======================================================================*/

 if(ngrp > 0){
  rmass1 = dvector(1,ngrp);
  rmass2 = dvector(1,ngrp);
  rmass3 = dvector(1,ngrp);
  rmass4 = dvector(1,ngrp);
       x = dmatrix(1,4,1,ngrp);
       y = dmatrix(1,4,1,ngrp);
       z = dmatrix(1,4,1,ngrp);

       vx = dmatrix(1,4,1,ngrp);
       vy = dmatrix(1,4,1,ngrp);
       vz = dmatrix(1,4,1,ngrp);
      p11 = dvector(1,ngrp);
      p12 = dvector(1,ngrp);
      p13 = dvector(1,ngrp);
      p22 = dvector(1,ngrp);
      p23 = dvector(1,ngrp);
      p33 = dvector(1,ngrp);
   rmassm = d3tensor(1,6,1,6,1,ngrp);
      dvx = dmatrix(1,6,1,ngrp);
      dvy = dmatrix(1,6,1,ngrp);
      dvz = dmatrix(1,6,1,ngrp);
       dx = dmatrix(1,6,1,ngrp);
       dy = dmatrix(1,6,1,ngrp);
       dz = dmatrix(1,6,1,ngrp);
    txlam = dvector(1,6);
     xlam = dmatrix(1,6,1,ngrp);
     amat = dmatrix(1,36,1,ngrp);
    tamat = dvector(1,36);
    ind1 = (int *)calloc((ngrp+1),sizeof(int));
    ind2 = (int *)calloc((ngrp+1),sizeof(int));
    ind3 = (int *)calloc((ngrp+1),sizeof(int));
    ind4 = (int *)calloc((ngrp+1),sizeof(int));
 }/*endif*/

/*=======================================================================*/

/* Malloc up some vectors and matrices */
 na = NCON_46;
 pnorm = 2.0/dt;

 ptens_pvten_tmp[1] = 0.0;
 ptens_pvten_tmp[2] = 0.0;
 ptens_pvten_tmp[3] = 0.0;
 ptens_pvten_tmp[4] = 0.0;
 ptens_pvten_tmp[5] = 0.0;
 ptens_pvten_tmp[6] = 0.0;
 ptens_pvten_tmp[7] = 0.0;
 ptens_pvten_tmp[8] = 0.0;
 ptens_pvten_tmp[9] = 0.0;

/* Gather masses, positions and velocities of atoms */
 for(igrp=1;igrp <= ngrp; igrp++) {
  ind1[igrp] = grp_bond_con_j1_46[(igrp+igrp_off)];
  ind2[igrp] = grp_bond_con_j2_46[(igrp+igrp_off)];
  ind3[igrp] = grp_bond_con_j3_46[(igrp+igrp_off)];
  ind4[igrp] = grp_bond_con_j4_46[(igrp+igrp_off)];
 }


 for(igrp=1;igrp <= ngrp; igrp++) {
     ktemp= ind1[igrp];
  x[1][igrp] = clatoms_x[ktemp];
  y[1][igrp] = clatoms_y[ktemp];
  z[1][igrp] = clatoms_z[ktemp];
  rmass1[igrp] = 1.0/clatoms_mass[ktemp];
 }

 for(igrp=1;igrp <= ngrp; igrp++) {
     ktemp= ind2[igrp];
  x[2][igrp] = clatoms_x[ktemp];
  y[2][igrp] = clatoms_y[ktemp];
  z[2][igrp] = clatoms_z[ktemp];
  rmass2[igrp] = 1.0/clatoms_mass[ktemp];
 }

 for(igrp=1;igrp <= ngrp; igrp++) {
     ktemp= ind3[igrp];
  x[3][igrp] = clatoms_x[ktemp];
  y[3][igrp] = clatoms_y[ktemp];
  z[3][igrp] = clatoms_z[ktemp];
  rmass3[igrp] = 1.0/clatoms_mass[ktemp];
 }

 for(igrp=1;igrp <= ngrp; igrp++) {
     ktemp= ind4[igrp];
  x[4][igrp] = clatoms_x[ktemp];
  y[4][igrp] = clatoms_y[ktemp];
  z[4][igrp] = clatoms_z[ktemp];
  rmass4[igrp] = 1.0/clatoms_mass[ktemp];
 }


 for(igrp=1;igrp <= ngrp; igrp++) {
   ktemp= ind1[igrp];
   ktemp3= ind3[igrp];
  roll_sci=1.0/clatoms_roll_sc[ktemp3];/*all roll scales the same in same cons*/
  vx[1][igrp] = clatoms_vx[ktemp]+x[1][igrp]*baro_v_lnv_g*roll_sci;
  vy[1][igrp] = clatoms_vy[ktemp]+y[1][igrp]*baro_v_lnv_g*roll_sci;
  vz[1][igrp] = clatoms_vz[ktemp]+z[1][igrp]*baro_v_lnv_g*roll_sci;
 }

 for(igrp=1;igrp <= ngrp; igrp++) {
   ktemp= ind2[igrp];
   ktemp3= ind3[igrp];
  roll_sci=1.0/clatoms_roll_sc[ktemp3];/*all roll scales the same in same cons*/
  vx[2][igrp] = clatoms_vx[ktemp]+x[2][igrp]*baro_v_lnv_g*roll_sci;
  vy[2][igrp] = clatoms_vy[ktemp]+y[2][igrp]*baro_v_lnv_g*roll_sci;
  vz[2][igrp] = clatoms_vz[ktemp]+z[2][igrp]*baro_v_lnv_g*roll_sci;
 }

 for(igrp=1;igrp <= ngrp; igrp++) {
   ktemp= ind3[igrp];
  roll_sci=1.0/clatoms_roll_sc[ktemp];/*all roll scales the same in same cons*/
  vx[3][igrp] = clatoms_vx[ktemp]+x[3][igrp]*baro_v_lnv_g*roll_sci;
  vy[3][igrp] = clatoms_vy[ktemp]+y[3][igrp]*baro_v_lnv_g*roll_sci;
  vz[3][igrp] = clatoms_vz[ktemp]+z[3][igrp]*baro_v_lnv_g*roll_sci;
 }

 for(igrp=1;igrp <= ngrp; igrp++) {
   ktemp= ind4[igrp];
   ktemp3= ind3[igrp];
  roll_sci=1.0/clatoms_roll_sc[ktemp3];/*all roll scales the same in same cons*/
  vx[4][igrp] = clatoms_vx[ktemp]+x[4][igrp]*baro_v_lnv_g*roll_sci;
  vy[4][igrp] = clatoms_vy[ktemp]+y[4][igrp]*baro_v_lnv_g*roll_sci;
  vz[4][igrp] = clatoms_vz[ktemp]+z[4][igrp]*baro_v_lnv_g*roll_sci;
 }
 

/* Set reciprocal mass matrix */
 for(igrp=1;igrp <= ngrp; igrp++) {
    rmass_1 = rmass1[igrp];
    rmass_2 = rmass2[igrp];
    rmass_3 = rmass3[igrp];
    rmass_4 = rmass4[igrp];


  rmassm[1][1][igrp] = -(rmass_1+rmass_2);
  rmassm[1][2][igrp] = rmassm[1][3][igrp] = -rmass_1;
  rmassm[1][4][igrp] = rmassm[1][5][igrp] = rmass_2;  
  rmassm[1][6][igrp] = 0.0;

  rmassm[2][1][igrp] = -rmass_1; rmassm[2][2][igrp] = -(rmass_1+rmass_3); 
                          rmassm[2][3][igrp] = -rmass_1;
  rmassm[2][4][igrp] = -rmass_3; rmassm[2][5][igrp] = 0.0;
  rmassm[2][6][igrp] = rmass_3;

  rmassm[3][1][igrp] = rmassm[3][2][igrp] = -rmass_1;
  rmassm[3][3][igrp] = -(rmass_1+rmass_4);
  rmassm[3][4][igrp] = 0.0; rmassm[3][5][igrp] = rmassm[3][6][igrp] = -rmass_4;

  rmassm[4][1][igrp] = rmass_2; rmassm[4][2][igrp] = -rmass_3;
  rmassm[4][3][igrp] = 0.0;
  rmassm[4][4][igrp] = -(rmass_2+rmass_3); rmassm[4][5][igrp] = -rmass_2; 
                                   rmassm[4][6][igrp] = rmass_3;

  rmassm[5][1][igrp] = rmass_2; rmassm[5][2][igrp] = 0.0;
  rmassm[5][3][igrp] = -rmass_4;
  rmassm[5][4][igrp] = -rmass_2; rmassm[5][5][igrp] = -(rmass_2+rmass_4); 
                          rmassm[5][6][igrp] = -rmass_4;

  rmassm[6][1][igrp] = 0.0; rmassm[6][2][igrp] = rmass_3;
  rmassm[6][3][igrp] = -rmass_4;
  rmassm[6][4][igrp] = rmass_3; rmassm[6][5][igrp] = -rmass_4; 
                         rmassm[6][6][igrp] = -(rmass_3+rmass_4);
 }


 for(igrp=1;igrp <= ngrp; igrp++) {
    dvx[1][igrp] = vx[1][igrp]-vx[2][igrp];
    dvx[2][igrp] = vx[1][igrp]-vx[3][igrp];
    dvx[3][igrp] = vx[1][igrp]-vx[4][igrp];
    dvx[4][igrp] = vx[2][igrp]-vx[3][igrp];
    dvx[5][igrp] = vx[2][igrp]-vx[4][igrp];
    dvx[6][igrp] = vx[3][igrp]-vx[4][igrp];
 }

 for(igrp=1;igrp <= ngrp; igrp++) {
    dvy[1][igrp] = vy[1][igrp]-vy[2][igrp];
    dvy[2][igrp] = vy[1][igrp]-vy[3][igrp];
    dvy[3][igrp] = vy[1][igrp]-vy[4][igrp];
    dvy[4][igrp] = vy[2][igrp]-vy[3][igrp];
    dvy[5][igrp] = vy[2][igrp]-vy[4][igrp];
    dvy[6][igrp] = vy[3][igrp]-vy[4][igrp];
 }

 for(igrp=1;igrp <= ngrp; igrp++) {
    dvz[1][igrp] = vz[1][igrp]-vz[2][igrp];
    dvz[2][igrp] = vz[1][igrp]-vz[3][igrp];
    dvz[3][igrp] = vz[1][igrp]-vz[4][igrp];
    dvz[4][igrp] = vz[2][igrp]-vz[3][igrp];
    dvz[5][igrp] = vz[2][igrp]-vz[4][igrp];
    dvz[6][igrp] = vz[3][igrp]-vz[4][igrp];
 }

 for(igrp=1;igrp <= ngrp; igrp++) {
    dx[1][igrp] = x[1][igrp]-x[2][igrp];
    dx[2][igrp] = x[1][igrp]-x[3][igrp];
    dx[3][igrp] = x[1][igrp]-x[4][igrp];
    dx[4][igrp] = x[2][igrp]-x[3][igrp];
    dx[5][igrp] = x[2][igrp]-x[4][igrp];
    dx[6][igrp] = x[3][igrp]-x[4][igrp];
 }

 for(igrp=1;igrp <= ngrp; igrp++) {
    dy[1][igrp] = y[1][igrp]-y[2][igrp];
    dy[2][igrp] = y[1][igrp]-y[3][igrp];
    dy[3][igrp] = y[1][igrp]-y[4][igrp];
    dy[4][igrp] = y[2][igrp]-y[3][igrp];
    dy[5][igrp] = y[2][igrp]-y[4][igrp];
    dy[6][igrp] = y[3][igrp]-y[4][igrp];
 }

 for(igrp=1;igrp <= ngrp; igrp++) {
    dz[1][igrp] = z[1][igrp]-z[2][igrp];
    dz[2][igrp] = z[1][igrp]-z[3][igrp];
    dz[3][igrp] = z[1][igrp]-z[4][igrp];
    dz[4][igrp] = z[2][igrp]-z[3][igrp];
    dz[5][igrp] = z[2][igrp]-z[4][igrp];
    dz[6][igrp] = z[3][igrp]-z[4][igrp];

 }/*endfor*/

/* ========================================================================== */
/* Get initial guess for lambda */
  iii = 0;
  for(i=1; i <= NCON_46; i++){
   for(j=1; j <= NCON_46; j++){
      iii++;
    for(igrp=1;igrp <= ngrp; igrp++) {
       amat[iii][igrp] =-rmassm[i][j][igrp]*
                 (dx[i][igrp]*dx[j][igrp] + dy[i][igrp]*dy[j][igrp] 
                + dz[i][igrp]*dz[j][igrp]);
    }/*endfor*/
   }/*endfor*/
  }/*endfor*/


  for(i=1; i <= NCON_46; i++){
    for(igrp=1;igrp <= ngrp; igrp++) {
      avec[i] = dvx[i][igrp]*dx[i][igrp] + dvy[i][igrp]*dy[i][igrp]  
              + dvz[i][igrp]*dz[i][igrp];
      xlam[i][igrp] = avec[i];
    }/*endfor*/
  }/*endfor*/

/* Solve linear system A xlam = avec */
  for(igrp=1;igrp <= ngrp; igrp++) {

   ipvt[1] = 0; ipvt[2] = 0; ipvt[3] = 0;
   ipvt[4] = 0; ipvt[5] = 0; ipvt[6] = 0;
     
    txlam[1] = xlam[1][igrp];
    txlam[2] = xlam[2][igrp];
    txlam[3] = xlam[3][igrp];
    txlam[4] = xlam[4][igrp];
    txlam[5] = xlam[5][igrp];
    txlam[6] = xlam[6][igrp];

    for(i=1;i<=36;i++) {
      tamat[i] = amat[i][igrp];
    }
 
#ifdef IBM_ESSL
  dgef(&(tamat[1]),&na,&na,&(ipvt[1]));
#else
  DGEFA(&(tamat[1]),&na,&na,&(ipvt[1]),&info);
#endif
  job = 1;
#ifdef IBM_ESSL
  job = 1;  /*changed from 0 */
  dges(&(tamat[1]),&na,&na,&(ipvt[1]),&(txlam[1]),&job);
#else
  DGESL(&(tamat[1]),&na,&na,&(ipvt[1]),&(txlam[1]),&job);
#endif

    xlam[1][igrp] =  txlam[1];
    xlam[2][igrp] =  txlam[2];
    xlam[3][igrp] =  txlam[3];
    xlam[4][igrp] =  txlam[4];
    xlam[5][igrp] =  txlam[5];
    xlam[6][igrp] =  txlam[6];
}/*end for*/

/* ======================================================================  */
/* Velocity update */

#ifndef NO_PRAGMA
#pragma IVDEP
#endif
 for(igrp=1;igrp <= ngrp; igrp++) {
   double xlam1,xlam2,xlam3,xlam4,xlam5,xlam6; 
   double dx1,dx2,dx3,dx4,dx5,dx6;
   double dy1,dy2,dy3,dy4,dy5,dy6;
   double dz1,dz2,dz3,dz4,dz5,dz6;
  
    ktemp1 = ind1[igrp];
    ktemp2 = ind2[igrp];
    ktemp3 = ind3[igrp];
    ktemp4 = ind4[igrp];

     dx1 = dx[1][igrp]; dx2 = dx[2][igrp];
     dx3 = dx[3][igrp]; dx4 = dx[4][igrp];
     dx5 = dx[5][igrp]; dx6 = dx[6][igrp];

     dy1 = dy[1][igrp]; dy2 = dy[2][igrp]; 
     dy3 = dy[3][igrp]; dy4 = dy[4][igrp]; 
     dy5 = dy[5][igrp]; dy6 = dy[6][igrp]; 

     dz1 = dz[1][igrp]; dz2 = dz[2][igrp]; 
     dz3 = dz[3][igrp]; dz4 = dz[4][igrp]; 
     dz5 = dz[5][igrp]; dz6 = dz[6][igrp]; 
  
     xlam1 = xlam[1][igrp];
     xlam2 = xlam[2][igrp];
     xlam3 = xlam[3][igrp];
     xlam4 = xlam[4][igrp];
     xlam5 = xlam[5][igrp];
     xlam6 = xlam[6][igrp];

     rmass_1 = rmass1[igrp];
     rmass_2 = rmass2[igrp];
     rmass_3 = rmass3[igrp];
     rmass_4 = rmass4[igrp];

  clatoms_vx[ktemp1] -=  ( xlam1*dx1 + xlam2*dx2 + xlam3*dx3)*rmass_1;
  clatoms_vy[ktemp1] -=  ( xlam1*dy1 + xlam2*dy2 + xlam3*dy3)*rmass_1;
  clatoms_vz[ktemp1] -=  ( xlam1*dz1 + xlam2*dz2 + xlam3*dz3)*rmass_1;

  clatoms_vx[ktemp2] -=  (-xlam1*dx1 + xlam4*dx4 + xlam5*dx5)*rmass_2;
  clatoms_vy[ktemp2] -=  (-xlam1*dy1 + xlam4*dy4 + xlam5*dy5)*rmass_2;
  clatoms_vz[ktemp2] -=  (-xlam1*dz1 + xlam4*dz4 + xlam5*dz5)*rmass_2;

  clatoms_vx[ktemp3] -=  (-xlam2*dx2 - xlam4*dx4 + xlam6*dx6)*rmass_3;
  clatoms_vy[ktemp3] -=  (-xlam2*dy2 - xlam4*dy4 + xlam6*dy6)*rmass_3;
  clatoms_vz[ktemp3] -=  (-xlam2*dz2 - xlam4*dz4 + xlam6*dz6)*rmass_3;

  clatoms_vx[ktemp4] -=  (-xlam3*dx3 - xlam5*dx5 - xlam6*dx6)*rmass_4;
  clatoms_vy[ktemp4] -=  (-xlam3*dy3 - xlam5*dy5 - xlam6*dy6)*rmass_4;
  clatoms_vz[ktemp4] -=  (-xlam3*dz3 - xlam5*dz5 - xlam6*dz6)*rmass_4;

/* Pressure tensor update */

     p11[igrp] = xlam1*dx1*dx1 + xlam2*dx2*dx2 + xlam3*dx3*dx3
               + xlam4*dx4*dx4 + xlam5*dx5*dx5 + xlam6*dx6*dx6;
     p22[igrp] = xlam1*dy1*dy1 + xlam2*dy2*dy2 + xlam3*dy3*dy3
               + xlam4*dy4*dy4 + xlam5*dy5*dy5 + xlam6*dy6*dy6;
     p33[igrp] = xlam1*dz1*dz1 + xlam2*dz2*dz2 + xlam3*dz3*dz3
               + xlam4*dz4*dz4 + xlam5*dz5*dz5 + xlam6*dz6*dz6;
     p12[igrp] = xlam1*dx1*dy1 + xlam2*dx2*dy2 + xlam3*dx3*dy3
               + xlam4*dx4*dy4 + xlam5*dx5*dy5 + xlam6*dx6*dy6;
     p13[igrp] = xlam1*dx1*dz1 + xlam2*dx2*dz2 + xlam3*dx3*dz3
               + xlam4*dx4*dz4 + xlam5*dx5*dz5 + xlam6*dx6*dz6;
     p23[igrp] = xlam1*dy1*dz1 + xlam2*dy2*dz2 + xlam3*dy3*dz3
               + xlam4*dy4*dz4 + xlam5*dy5*dz5 + xlam6*dy6*dz6;
  }/*endfor*/

#ifndef NO_PRAGMA
#pragma IVDEP
#endif
 for(igrp=1;igrp <= ngrp; igrp++) {
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
   free_dvector(rmass1,1,ngrp);
   free_dvector(rmass2,1,ngrp);
   free_dvector(rmass3,1,ngrp);
   free_dvector(rmass4,1,ngrp);

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

   free_d3tensor(rmassm,1,6,1,6,1,ngrp);

   free_dmatrix(dvx,1,6,1,ngrp);
   free_dmatrix(dvy,1,6,1,ngrp);
   free_dmatrix(dvz,1,6,1,ngrp);

   free_dmatrix(dx,1,6,1,ngrp);
   free_dmatrix(dy,1,6,1,ngrp);
   free_dmatrix(dz,1,6,1,ngrp);

   free_dvector(txlam,1,6);
   free_dmatrix(xlam,1,6,1,ngrp);
   free_dmatrix(amat,1,36,1,ngrp);
   free_dvector(tamat,1,36);

   free(ind1); free(ind2); free(ind3); free(ind4);
 }/*endif*/
   
/*=======================================================================*/
/*=======================================================================*/
} /* end routine */
/*=======================================================================*/


