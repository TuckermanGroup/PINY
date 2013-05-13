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
#include "../proto_defs/proto_math.h"
#include "../proto_defs/proto_friend_lib_entry.h"
#include "../proto_defs/proto_communicate_wrappers.h"
#include "../proto_defs/proto_energy_ctrl_entry.h"

#define NCON_46 6
#define NAT_46 4

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void shake_46_rollf(GRP_BOND_CON *grp_bond_con,
                    CLATOMS_INFO *clatoms_info,CLATOMS_POS *clatoms_pos,
                    PTENS *ptens,double dt,double *aiter,
                    PAR_RAHMAN *par_rahman, int ifirst,CELL *cell,
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
 double prrm1, prrm2, prrm3, prrm4, prrm5, prrm6, prrm7, prrm8, prrm9;
 int i,j,k,iii;
 int iter,igrp,*ind1,*ind2,*ind3,*ind4,jtyp;
 int na,job,info,ipvt[NCON_46+1];       /* For dgefa and dgesl */
 int ktemp,ktemp1,ktemp2,ktemp3,ktemp4;

   double dx1,dx2,dx3,dx4,dx5,dx6;
   double dy1,dy2,dy3,dy4,dy5,dy6;
   double dz1,dz2,dz3,dz4,dz5,dz6;
   double dx01,dx02,dx03,dx04,dx05,dx06;
   double dy01,dy02,dy03,dy04,dy05,dy06;
   double dz01,dz02,dz03,dz04,dz05,dz06;
   double xlam1,xlam2,xlam3,xlam4,xlam5,xlam6;


/* TOP */
 double *rm1,*rm2,*rm3,*rm4; 
 double ***rmm;
 double **dx,**dy,**dz;
 double **dx0,**dy0,**dz0;
 double **dxt,**dyt,**dzt;
 double **dxn,**dyn,**dzn;
 double **avec,**amat,**dij,**dxl,**xlam;
 double *dlmax,*txlam;
 double **x,**y,**z;
 double **xo,**yo,**zo;
 double *p11,*p22,*p33,*p12,*p13,*p23;
 double *tamat;
 double pnorm;

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

  int ngrp                     = grp_bond_con->num_46;
  int max_iter                 = grp_bond_con->max_iter;
  double tol                   = grp_bond_con->tol;
  int *grp_bond_con_j1_46      = grp_bond_con->j1_46;
  int *grp_bond_con_j2_46      = grp_bond_con->j2_46;
  int *grp_bond_con_j3_46      = grp_bond_con->j3_46;
  int *grp_bond_con_j4_46      = grp_bond_con->j4_46;
  int *grp_bond_con_jtyp_46    = grp_bond_con->jtyp_46;
  double **grp_bond_con_eq_46  = grp_bond_con->eq_46;
  double **grp_bond_con_al_46  = grp_bond_con->al_46;

  double *pvten_inc            = ptens->pvten_inc;
  double *pvten_tmp            = ptens->pvten_tmp;
  double *pvten_tmp2           = ptens->pvten_tmp_res;

  double *roll_mtv             = par_rahman->roll_mtv;
  double *fgmat_p              = par_rahman->fgmat_p;
  double *vgmat                = par_rahman->vgmat;
  double roll_scg              = par_rahman->roll_scg;
  double mass_hm               = par_rahman->mass_hm;

  int iperd                    = cell->iperd;
  int hmat_cons_typ            = cell->hmat_cons_typ;
  int hmat_int_typ             = cell->hmat_int_typ;

  int np_forc                  = class_comm_forc_pkg->num_proc;
  int myid_forc                = class_comm_forc_pkg->myid;
  MPI_Comm comm_forc           = class_comm_forc_pkg->comm;

/*=======================================================================*/
  if(ngrp > 0){
     rm1= dvector(1,ngrp);
     rm2= dvector(1,ngrp);
     rm3= dvector(1,ngrp);
     rm4= dvector(1,ngrp);
      dlmax= dvector(1,6);
      txlam= dvector(1,6);
     rmm= d3tensor(1,6,1,6,1,ngrp); 
         dx= dmatrix(1,6,1,ngrp);
         dy= dmatrix(1,6,1,ngrp);
         dz= dmatrix(1,6,1,ngrp);
        dxt= dmatrix(1,6,1,ngrp);
        dyt= dmatrix(1,6,1,ngrp);
        dzt= dmatrix(1,6,1,ngrp);
        dxn= dmatrix(1,6,1,ngrp);
        dyn= dmatrix(1,6,1,ngrp);
        dzn= dmatrix(1,6,1,ngrp);
        dx0= dmatrix(1,6,1,ngrp);
        dy0= dmatrix(1,6,1,ngrp);
        dz0= dmatrix(1,6,1,ngrp);
        dij= dmatrix(1,6,1,ngrp);
        dxl= dmatrix(1,6,1,ngrp);
       amat= dmatrix(1,36,1,ngrp);
       tamat= dvector(1,36);
       avec= dmatrix(1,6,1,ngrp);
       xlam= dmatrix(1,6,1,ngrp);
          x= dmatrix(1,4,1,ngrp);
          y= dmatrix(1,4,1,ngrp);
          z= dmatrix(1,4,1,ngrp);
          xo= dmatrix(1,4,1,ngrp);
          yo= dmatrix(1,4,1,ngrp);
          zo= dmatrix(1,4,1,ngrp);
     ind1= (int *)calloc((ngrp+1),sizeof(int));
     ind2= (int *)calloc((ngrp+1),sizeof(int));
     ind3= (int *)calloc((ngrp+1),sizeof(int));
     ind4= (int *)calloc((ngrp+1),sizeof(int));
         p11= dvector(1,ngrp);
         p12= dvector(1,ngrp);
         p13= dvector(1,ngrp);
         p22= dvector(1,ngrp);
         p23= dvector(1,ngrp);
         p33= dvector(1,ngrp);
  }/*endif*/

/*=======================================================================*/
/* Malloc up some vectors and matrices */
 na = NCON_46;
/* AA */

 dts = dt*dt;
 pnorm = 2.0/dts;
 *aiter = 0.0;
 pvten_tmp[1] = 0.0;
 pvten_tmp[2] = 0.0;
 pvten_tmp[3] = 0.0;
 pvten_tmp[4] = 0.0;
 pvten_tmp[5] = 0.0;
 pvten_tmp[6] = 0.0;
 pvten_tmp[7] = 0.0;
 pvten_tmp[8] = 0.0;
 pvten_tmp[9] = 0.0;


 if(ifirst == 2){
  for(igrp=1;igrp <= ngrp; igrp++) {
   grp_bond_con_al_46[1][igrp] = 0.0;
   grp_bond_con_al_46[2][igrp] = 0.0;
   grp_bond_con_al_46[3][igrp] = 0.0;
   grp_bond_con_al_46[4][igrp] = 0.0;
   grp_bond_con_al_46[5][igrp] = 0.0;
   grp_bond_con_al_46[6][igrp] = 0.0;
  }/*endif*/
 }/*endif*/

 for(igrp=1;igrp <= ngrp; igrp++) {
    ind1[igrp] = grp_bond_con_j1_46[igrp];
    ind2[igrp] = grp_bond_con_j2_46[igrp];
    ind3[igrp] = grp_bond_con_j3_46[igrp];
    ind4[igrp] = grp_bond_con_j4_46[igrp];
 }/*end for*/

 for(igrp=1;igrp <= ngrp; igrp++) {
    rm1[igrp] = 1.0/clatoms_mass[ind1[igrp]];
    rm2[igrp] = 1.0/clatoms_mass[ind2[igrp]];
    rm3[igrp] = 1.0/clatoms_mass[ind3[igrp]];
    rm4[igrp] = 1.0/clatoms_mass[ind4[igrp]];
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
    ktemp = ind4[igrp];
    x[4][igrp] = clatoms_x[ktemp];
    y[4][igrp] = clatoms_y[ktemp];
    z[4][igrp] = clatoms_z[ktemp];
  }
 for(igrp=1;igrp <= ngrp; igrp++) {
    ktemp = ind1[igrp];
    xo[1][igrp] = clatoms_xold[ktemp];
    yo[1][igrp] = clatoms_yold[ktemp];
    zo[1][igrp] = clatoms_zold[ktemp];
 }
for(igrp=1;igrp <= ngrp; igrp++) {
    ktemp = ind2[igrp];
    xo[2][igrp] = clatoms_xold[ktemp];
    yo[2][igrp] = clatoms_yold[ktemp];
    zo[2][igrp] = clatoms_zold[ktemp];
 }
 for(igrp=1;igrp <= ngrp; igrp++) {
    ktemp = ind3[igrp];
    xo[3][igrp] = clatoms_xold[ktemp];
    yo[3][igrp] = clatoms_yold[ktemp];
    zo[3][igrp] = clatoms_zold[ktemp];
 }
 for(igrp=1;igrp <= ngrp; igrp++) {
    ktemp = ind4[igrp];
    xo[4][igrp] = clatoms_xold[ktemp];
    yo[4][igrp] = clatoms_yold[ktemp];
    zo[4][igrp] = clatoms_zold[ktemp];
 }
/* II collect jtyps  */
 for(igrp=1;igrp <= ngrp; igrp++) {
    jtyp = grp_bond_con_jtyp_46[igrp];

    dij[1][igrp] = grp_bond_con_eq_46[1][jtyp];
    dij[2][igrp] = grp_bond_con_eq_46[2][jtyp];
    dij[3][igrp] = grp_bond_con_eq_46[3][jtyp];
    dij[4][igrp] = grp_bond_con_eq_46[4][jtyp];
    dij[5][igrp] = grp_bond_con_eq_46[5][jtyp];
    dij[6][igrp] = grp_bond_con_eq_46[6][jtyp];

 } /*end for*/

/* III calculate mass tensor */
 for(igrp=1;igrp <= ngrp; igrp++) {
    rms1=rm1[igrp];
    rms2=rm2[igrp];
    rms3=rm3[igrp];
    rms4=rm4[igrp];

  rmm[1][1][igrp] = -(rms1+rms2); rmm[1][2][igrp] = rmm[1][3][igrp] = -rms1;
  rmm[1][4][igrp] = rmm[1][5][igrp] = rms2; rmm[1][6][igrp] = 0.0;

  rmm[2][1][igrp] = -rms1; rmm[2][2][igrp] = -(rms1+rms3); rmm[2][3][igrp] = -rms1;
  rmm[2][4][igrp] = -rms3; rmm[2][5][igrp] = 0.0;   rmm[2][6][igrp] = rms3;

  rmm[3][1][igrp] = rmm[3][2][igrp] = -rms1; rmm[3][3][igrp] = -(rms1+rms4);
  rmm[3][4][igrp] = 0.0; rmm[3][5][igrp] = rmm[3][6][igrp] = -rms4;

  rmm[4][1][igrp] = rms2; rmm[4][2][igrp] = -rms3; rmm[4][3][igrp] = 0.0;
  rmm[4][4][igrp] = -(rms2+rms3); rmm[4][5][igrp] = -rms2; rmm[4][6][igrp] = rms3;

  rmm[5][1][igrp] = rms2; rmm[5][2][igrp] = 0.0; rmm[5][3][igrp] = -rms4;
  rmm[5][4][igrp] = -rms2; rmm[5][5][igrp] = -(rms2+rms4); rmm[5][6][igrp] = -rms4;
  rmm[6][1][igrp] = 0.0; rmm[6][2][igrp] = rms3; rmm[6][3][igrp] = -rms4;
  rmm[6][4][igrp] = rms3; rmm[6][5][igrp] = -rms4; rmm[6][6][igrp] = -(rms3+rms4);
}/*endfor*/


/* Compute difference vectors */
 for(igrp=1;igrp <= ngrp; igrp++) {
     dxt[1][igrp]= x[1][igrp]-x[2][igrp];
     dxt[2][igrp]= x[1][igrp]-x[3][igrp];
     dxt[3][igrp]= x[1][igrp]-x[4][igrp];
     dxt[4][igrp]= x[2][igrp]-x[3][igrp];
     dxt[5][igrp]= x[2][igrp]-x[4][igrp];
     dxt[6][igrp]= x[3][igrp]-x[4][igrp];
   }
 for(igrp=1;igrp <= ngrp; igrp++) {
     dyt[1][igrp]= y[1][igrp]-y[2][igrp];
     dyt[2][igrp]= y[1][igrp]-y[3][igrp];
     dyt[3][igrp]= y[1][igrp]-y[4][igrp];
     dyt[4][igrp]= y[2][igrp]-y[3][igrp];
     dyt[5][igrp]= y[2][igrp]-y[4][igrp];
     dyt[6][igrp]= y[3][igrp]-y[4][igrp];
   }
 for(igrp=1;igrp <= ngrp; igrp++) {
     dzt[1][igrp]= z[1][igrp]-z[2][igrp];
     dzt[2][igrp]= z[1][igrp]-z[3][igrp];
     dzt[3][igrp]= z[1][igrp]-z[4][igrp];
     dzt[4][igrp]= z[2][igrp]-z[3][igrp];
     dzt[5][igrp]= z[2][igrp]-z[4][igrp];
     dzt[6][igrp]= z[3][igrp]-z[4][igrp];
   }
 for(igrp=1;igrp <= ngrp; igrp++) {
     dx0[1][igrp]= xo[1][igrp]-xo[2][igrp];
     dx0[2][igrp]= xo[1][igrp]-xo[3][igrp];
     dx0[3][igrp]= xo[1][igrp]-xo[4][igrp];
     dx0[4][igrp]= xo[2][igrp]-xo[3][igrp];
     dx0[5][igrp]= xo[2][igrp]-xo[4][igrp];
     dx0[6][igrp]= xo[3][igrp]-xo[4][igrp];
   }
 for(igrp=1;igrp <= ngrp; igrp++) {

     dy0[1][igrp]= yo[1][igrp]-yo[2][igrp];
     dy0[2][igrp]= yo[1][igrp]-yo[3][igrp];
     dy0[3][igrp]= yo[1][igrp]-yo[4][igrp];
     dy0[4][igrp]= yo[2][igrp]-yo[3][igrp];
     dy0[5][igrp]= yo[2][igrp]-yo[4][igrp];
     dy0[6][igrp]= yo[3][igrp]-yo[4][igrp];
   }
 for(igrp=1;igrp <= ngrp; igrp++) {
     dz0[1][igrp]= zo[1][igrp]-zo[2][igrp];
     dz0[2][igrp]= zo[1][igrp]-zo[3][igrp];
     dz0[3][igrp]= zo[1][igrp]-zo[4][igrp];
     dz0[4][igrp]= zo[2][igrp]-zo[3][igrp];
     dz0[5][igrp]= zo[2][igrp]-zo[4][igrp];
     dz0[6][igrp]= zo[3][igrp]-zo[4][igrp];
   }

    prrm1= roll_mtv[1]; prrm2= roll_mtv[2];
    prrm3= roll_mtv[3]; prrm4= roll_mtv[4];
    prrm5= roll_mtv[5]; prrm6= roll_mtv[6];
    prrm7= roll_mtv[7]; prrm8= roll_mtv[8];
    prrm9= roll_mtv[9];

 for(igrp=1;igrp <= ngrp; igrp++) {
    dx[1][igrp]= dx0[1][igrp]*prrm1 + dy0[1][igrp]*prrm2 
                +dz0[1][igrp]*prrm3;
    dx[2][igrp]= dx0[2][igrp]*prrm1 + dy0[2][igrp]*prrm2 
               + dz0[2][igrp]*prrm3;
    dx[3][igrp]= dx0[3][igrp]*prrm1 + dy0[3][igrp]*prrm2 
               + dz0[3][igrp]*prrm3;
    dx[4][igrp]= dx0[4][igrp]*prrm1 + dy0[4][igrp]*prrm2 
               + dz0[4][igrp]*prrm3;
    dx[5][igrp]= dx0[5][igrp]*prrm1 + dy0[5][igrp]*prrm2 
               + dz0[5][igrp]*prrm3;
    dx[6][igrp]= dx0[6][igrp]*prrm1 + dy0[6][igrp]*prrm2 
               + dz0[6][igrp]*prrm3;
}/*end for*/

 for(igrp=1;igrp <= ngrp; igrp++) {
    dy[1][igrp]= dx0[1][igrp]*prrm4 + dy0[1][igrp]*prrm5 
               + dz0[1][igrp]*prrm6;
    dy[2][igrp]= dx0[2][igrp]*prrm4 + dy0[2][igrp]*prrm5 
               + dz0[2][igrp]*prrm6;
    dy[3][igrp]= dx0[3][igrp]*prrm4 + dy0[3][igrp]*prrm5 
               + dz0[3][igrp]*prrm6;
    dy[4][igrp]= dx0[4][igrp]*prrm4 + dy0[4][igrp]*prrm5 
               + dz0[4][igrp]*prrm6;
    dy[5][igrp]= dx0[5][igrp]*prrm4 + dy0[5][igrp]*prrm5 
               + dz0[5][igrp]*prrm6;
    dy[6][igrp]= dx0[6][igrp]*prrm4 + dy0[6][igrp]*prrm5 
               + dz0[6][igrp]*prrm6;
 }/*end for*/

 for(igrp=1;igrp <= ngrp; igrp++) {
    dz[1][igrp]= dx0[1][igrp]*prrm7 + dy0[1][igrp]*prrm8 
               + dz0[1][igrp]*prrm9;
    dz[2][igrp]= dx0[2][igrp]*prrm7 + dy0[2][igrp]*prrm8 
               + dz0[2][igrp]*prrm9;
    dz[3][igrp]= dx0[3][igrp]*prrm7 + dy0[3][igrp]*prrm8 
               + dz0[3][igrp]*prrm9;
    dz[4][igrp]= dx0[4][igrp]*prrm7 + dy0[4][igrp]*prrm8 
               + dz0[4][igrp]*prrm9;
    dz[5][igrp]= dx0[5][igrp]*prrm7 + dy0[5][igrp]*prrm8 
               + dz0[5][igrp]*prrm9;
    dz[6][igrp]= dx0[6][igrp]*prrm7 + dy0[6][igrp]*prrm8 
               + dz0[6][igrp]*prrm9;
}/* end for*/


 for(igrp=1;igrp <= ngrp; igrp++) {
   avec[1][igrp] = dij[1][igrp]*dij[1][igrp] -(dxt[1][igrp]*dxt[1][igrp]
                 + dyt[1][igrp]*dyt[1][igrp] + dzt[1][igrp]*dzt[1][igrp]);

   avec[2][igrp] = dij[2][igrp]*dij[2][igrp] -(dxt[2][igrp]*dxt[2][igrp]
                 + dyt[2][igrp]*dyt[2][igrp] + dzt[2][igrp]*dzt[2][igrp]);

   avec[3][igrp] = dij[3][igrp]*dij[3][igrp] -(dxt[3][igrp]*dxt[3][igrp]
                 + dyt[3][igrp]*dyt[3][igrp] + dzt[3][igrp]*dzt[3][igrp]);

   avec[4][igrp] = dij[4][igrp]*dij[4][igrp] -(dxt[4][igrp]*dxt[4][igrp]
                 + dyt[4][igrp]*dyt[4][igrp] + dzt[4][igrp]*dzt[4][igrp]);

   avec[5][igrp] = dij[5][igrp]*dij[5][igrp] -(dxt[5][igrp]*dxt[5][igrp]
                 + dyt[5][igrp]*dyt[5][igrp] + dzt[5][igrp]*dzt[5][igrp]);

   avec[6][igrp] = dij[6][igrp]*dij[6][igrp] -(dxt[6][igrp]*dxt[6][igrp]
                 + dyt[6][igrp]*dyt[6][igrp] + dzt[6][igrp]*dzt[6][igrp]);
 }/*end for*/

/* Get initial guess for lambda */
 if(ifirst == 2 || ifirst == 0){
 
  iii = 0;
  for(i=1; i <= NCON_46; i++){
   for(j=1; j <= NCON_46; j++){
    iii++;
 for(igrp=1;igrp <= ngrp; igrp++) {
    amat[iii][igrp] = 2.0*rmm[i][j][igrp]*
               (dxt[i][igrp]*dx[j][igrp] + dyt[i][igrp]*dy[j][igrp]
              + dzt[i][igrp]*dz[j][igrp]);
    }/*endfor*/
   }/*endfor*/
  }/*endfor*/
 for(igrp=1;igrp <= ngrp; igrp++) {

   txlam[1] = xlam[1][igrp] = avec[1][igrp];
   txlam[2] = xlam[2][igrp] = avec[2][igrp];
   txlam[3] = xlam[3][igrp] = avec[3][igrp];
   txlam[4] = xlam[4][igrp] = avec[4][igrp];
   txlam[5] = xlam[5][igrp] = avec[5][igrp];
   txlam[6] = xlam[6][igrp] = avec[6][igrp];

     for(i=1;i<=36;i++) {
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
   xlam[1][igrp] = txlam[1] ;
   xlam[2][igrp] = txlam[2] ;
   xlam[3][igrp] = txlam[3] ;
   xlam[4][igrp] = txlam[4] ;
   xlam[5][igrp] = txlam[5] ;
   xlam[6][igrp] = txlam[6] ;
         
}/* end loop over groups */
 } else {
 for(igrp=1;igrp <= ngrp; igrp++) {
   xlam[1][igrp] = grp_bond_con_al_46[1][igrp];
   xlam[2][igrp] = grp_bond_con_al_46[2][igrp];
   xlam[3][igrp] = grp_bond_con_al_46[3][igrp];
   xlam[4][igrp] = grp_bond_con_al_46[4][igrp];
   xlam[5][igrp] = grp_bond_con_al_46[5][igrp];
   xlam[6][igrp] = grp_bond_con_al_46[6][igrp];

   grp_bond_con_al_46[1][igrp] = 0.0;
   grp_bond_con_al_46[2][igrp] = 0.0;
   grp_bond_con_al_46[3][igrp] = 0.0;
   grp_bond_con_al_46[4][igrp] = 0.0;
   grp_bond_con_al_46[5][igrp] = 0.0;
   grp_bond_con_al_46[6][igrp] = 0.0;

 }/* end loop over groups */
 }/* end if */

/* Iterative loop to convergence */

 if(ngrp > 0){

  dmax = 1.0;
  iter = 0;
  do {
   ++iter;
   if(iter > max_iter) {
    printf("$$$$$$$$$$$$$$$$$$$$_WARNING_$$$$$$$$$$$$$$$$$$$$\n");
    printf("Group constraint Shake not converged after %d iterations.\n",
            max_iter);
    printf("The present tolerance is %g \n",dmax);
    printf("The desired tolerance is %g \n",tol);
    printf("$$$$$$$$$$$$$$$$$$$$_WARNING_$$$$$$$$$$$$$$$$$$$$\n");
    fflush(stdout);
    break;
   }/*endif*/
/* Set up guess of difference vectors */
   for(i=1; i <= NCON_46; i++) {
    for(igrp=1;igrp <= ngrp; igrp++) {
    dxn[i][igrp] = 2.0*dxt[i][igrp];
    dyn[i][igrp] = 2.0*dyt[i][igrp];
    dzn[i][igrp] = 2.0*dzt[i][igrp];
     }/*end loop over groups */
   }/*endfor*/
   for(i=1; i <= NCON_46; i++) {
    for(j=1; j <= NCON_46; j++) {
     for(igrp=1;igrp <= ngrp; igrp++) {
     dxn[i][igrp] += rmm[i][j][igrp]*xlam[j][igrp]*dx[j][igrp];
     dyn[i][igrp] += rmm[i][j][igrp]*xlam[j][igrp]*dy[j][igrp];
     dzn[i][igrp] += rmm[i][j][igrp]*xlam[j][igrp]*dz[j][igrp];
      }/*end loop over group */
    }/*endfor*/
   }/*endfor*/

/* Construct A-matrix */

   iii = 0;
   for(i=1; i <= NCON_46; i++) {
    for(j=1; j <= NCON_46; j++) {
       iii++;
 for(igrp=1;igrp <= ngrp; igrp++) {
      amat[iii][igrp] = rmm[i][j][igrp]*
                (dxn[i][igrp]*dx[j][igrp] + dyn[i][igrp]*dy[j][igrp]
               + dzn[i][igrp]*dz[j][igrp]);
      }/*end loop over groups */
    }/*endfor*/
   }/*endfor*/

 for(igrp=1;igrp <= ngrp; igrp++) {
    xl0[1] = xlam[1][igrp];
    xl0[2] = xlam[2][igrp];
    xl0[3] = xlam[3][igrp];
    xl0[4] = xlam[4][igrp];
    xl0[5] = xlam[5][igrp];
    xl0[6] = xlam[6][igrp];

    txlam[1] = xlam[1][igrp] = avec[1][igrp];
    txlam[2] = xlam[2][igrp] = avec[2][igrp];
    txlam[3] = xlam[3][igrp] = avec[3][igrp];
    txlam[4] = xlam[4][igrp] = avec[4][igrp];
    txlam[5] = xlam[5][igrp] = avec[5][igrp];
    txlam[6] = xlam[6][igrp] = avec[6][igrp];

     for(i=1;i<=36;i++) {
       tamat[i]= amat[i][igrp];
     }

/* Solve linear system A xlam = avec */

   for(i=1; i <=na;i++){ ipvt[i] = 0;}
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
    xlam[1][igrp] = txlam[1]; xlam[2][igrp] = txlam[2];
    xlam[3][igrp] = txlam[3]; xlam[4][igrp] = txlam[4];
    xlam[5][igrp] = txlam[5]; xlam[6][igrp] = txlam[6];


    dxl[1][igrp] = fabs(xlam[1][igrp]-xl0[1]);
    dxl[2][igrp] = fabs(xlam[2][igrp]-xl0[2]);
    dxl[3][igrp] = fabs(xlam[3][igrp]-xl0[3]);
    dxl[4][igrp] = fabs(xlam[4][igrp]-xl0[4]);
    dxl[5][igrp] = fabs(xlam[5][igrp]-xl0[5]);
    dxl[6][igrp] = fabs(xlam[6][igrp]-xl0[6]);
 }/*end loop over groups */

    dlmax[1]= dxl[1][1];
    dlmax[2]= dxl[2][1];
    dlmax[3]= dxl[3][1];
    dlmax[4]= dxl[4][1];
    dlmax[5]= dxl[5][1];
    dlmax[6]= dxl[6][1];

 for(igrp=2;igrp <= ngrp; igrp++) {
    dlmax[1]= (dlmax[1] > dxl[1][igrp] ? dlmax[1] : dxl[1][igrp]);
    dlmax[2]= (dlmax[2] > dxl[2][igrp] ? dlmax[2] : dxl[2][igrp]);
    dlmax[3]= (dlmax[3] > dxl[3][igrp] ? dlmax[3] : dxl[3][igrp]);
    dlmax[4]= (dlmax[4] > dxl[4][igrp] ? dlmax[4] : dxl[4][igrp]);
    dlmax[5]= (dlmax[5] > dxl[5][igrp] ? dlmax[5] : dxl[5][igrp]);
    dlmax[6]= (dlmax[6] > dxl[6][igrp] ? dlmax[6] : dxl[6][igrp]);
   }/*end loop over groups */
     dmax= dlmax[1];
    for(i=2;i <= 6; i++) {
    dmax = (dmax > dlmax[i] ? dmax:dlmax[i]);
   }/*endfor*/
  } while(dmax > tol);
  *aiter += (double) iter;
}/*endif for ngrp > 0*/

/* Position update */
#ifndef NO_PRAGMA
#pragma IVDEP
#endif
 for(igrp=1;igrp <= ngrp; igrp++) {
    ktemp1 = ind1[igrp];
    ktemp2 = ind2[igrp];
    ktemp3 = ind3[igrp];
    ktemp4 = ind4[igrp];

  dx1= dx[1][igrp]; dx2= dx[2][igrp]; dx3= dx[3][igrp];
  dx4= dx[4][igrp]; dx5= dx[5][igrp]; dx6= dx[6][igrp];

  dy1= dy[1][igrp]; dy2= dy[2][igrp]; dy3= dy[3][igrp];
  dy4= dy[4][igrp]; dy5= dy[5][igrp]; dy6= dy[6][igrp];

  dz1= dz[1][igrp]; dz2= dz[2][igrp]; dz3= dz[3][igrp];
  dz4= dz[4][igrp]; dz5= dz[5][igrp]; dz6= dz[6][igrp];


  dx01= dx0[1][igrp]; dx02= dx0[2][igrp]; dx03= dx0[3][igrp];
  dx04= dx0[4][igrp]; dx05= dx0[5][igrp]; dx06= dx0[6][igrp];

  dy01= dy0[1][igrp]; dy02= dy0[2][igrp]; dy03= dy0[3][igrp];
  dy04= dy0[4][igrp]; dy05= dy0[5][igrp]; dy06= dy0[6][igrp];

  dz01= dz0[1][igrp]; dz02= dz0[2][igrp]; dz03= dz0[3][igrp];
  dz04= dz0[4][igrp]; dz05= dz0[5][igrp]; dz06= dz0[6][igrp];


  xlam1= xlam[1][igrp]; xlam2= xlam[2][igrp];
  xlam3= xlam[3][igrp]; xlam4= xlam[4][igrp];
  xlam5= xlam[5][igrp]; xlam6= xlam[6][igrp];

  clatoms_x[ktemp1] -= ( xlam1*dx1 + xlam2*dx2 + xlam3*dx3)*rm1[igrp];
  clatoms_y[ktemp1] -= ( xlam1*dy1 + xlam2*dy2 + xlam3*dy3)*rm1[igrp];
  clatoms_z[ktemp1] -= ( xlam1*dz1 + xlam2*dz2 + xlam3*dz3)*rm1[igrp];

  clatoms_x[ktemp2] -=  (-xlam1*dx1 + xlam4*dx4 + xlam5*dx5)*rm2[igrp];
  clatoms_y[ktemp2] -=  (-xlam1*dy1 + xlam4*dy4 + xlam5*dy5)*rm2[igrp];
  clatoms_z[ktemp2] -=  (-xlam1*dz1 + xlam4*dz4 + xlam5*dz5)*rm2[igrp];

  clatoms_x[ktemp3] -=  (-xlam2*dx2 - xlam4*dx4 + xlam6*dx6)*rm3[igrp];
  clatoms_y[ktemp3] -=  (-xlam2*dy2 - xlam4*dy4 + xlam6*dy6)*rm3[igrp];
  clatoms_z[ktemp3] -=  (-xlam2*dz2 - xlam4*dz4 + xlam6*dz6)*rm3[igrp];

  clatoms_x[ktemp4] -=  (-xlam3*dx3 - xlam5*dx5 - xlam6*dx6)*rm4[igrp];
  clatoms_y[ktemp4] -=  (-xlam3*dy3 - xlam5*dy5 - xlam6*dy6)*rm4[igrp];
  clatoms_z[ktemp4] -=  (-xlam3*dz3 - xlam5*dz5 - xlam6*dz6)*rm4[igrp];


/* Velocity update */

  clatoms_vx[ktemp1] -=  ( xlam1*dx01 + xlam2*dx02 + xlam3*dx03)
                       *rm1[igrp]/dt;
  clatoms_vy[ktemp1] -=  ( xlam1*dy01 + xlam2*dy02 + xlam3*dy03)
                       *rm1[igrp]/dt;
  clatoms_vz[ktemp1] -=  ( xlam1*dz01 + xlam2*dz02 + xlam3*dz03)
                       *rm1[igrp]/dt;

  clatoms_vx[ktemp2] -=  (-xlam1*dx01 + xlam4*dx04 + xlam5*dx05)
                       *rm2[igrp]/dt;
  clatoms_vy[ktemp2] -=  (-xlam1*dy01 + xlam4*dy04 + xlam5*dy05)
                       *rm2[igrp]/dt;
  clatoms_vz[ktemp2] -=  (-xlam1*dz01 + xlam4*dz04 + xlam5*dz05)
                       *rm2[igrp]/dt;

  clatoms_vx[ktemp3] -=  (-xlam2*dx02 - xlam4*dx04 + xlam6*dx06)
                       *rm3[igrp]/dt;
  clatoms_vy[ktemp3] -=  (-xlam2*dy02 - xlam4*dy04 + xlam6*dy06)
                       *rm3[igrp]/dt;
  clatoms_vz[ktemp3] -=  (-xlam2*dz02 - xlam4*dz04 + xlam6*dz06)
                       *rm3[igrp]/dt;

  clatoms_vx[ktemp4] -=  (-xlam3*dx03 - xlam5*dx05 - xlam6*dx06)
                       *rm4[igrp]/dt;
  clatoms_vy[ktemp4] -=  (-xlam3*dy03 - xlam5*dy05 - xlam6*dy06)
                       *rm4[igrp]/dt;
  clatoms_vz[ktemp4] -=  (-xlam3*dz03 - xlam5*dz05 - xlam6*dz06)
                       *rm4[igrp]/dt;



/* Pressure tensor update */
/* Compute difference vectors: use unscaled old distances */
     p11[igrp] = (xlam1*dx01*dx01 + xlam2*dx02*dx02 + xlam3*dx03*dx03
            +xlam4*dx04*dx04 + xlam5*dx05*dx05 + xlam6*dx06*dx06);
     p22[igrp] = (xlam1*dy01*dy01 + xlam2*dy02*dy02 + xlam3*dy03*dy03
           + xlam4*dy04*dy04 + xlam5*dy05*dy05 + xlam6*dy06*dy06);
     p33[igrp] = (xlam1*dz01*dz01 + xlam2*dz02*dz02 + xlam3*dz03*dz03
           + xlam4*dz04*dz04 + xlam5*dz05*dz05 + xlam6*dz06*dz06 );
     p12[igrp] = (xlam1*dx01*dy01 + xlam2*dx02*dy02 + xlam3*dx03*dy03
           + xlam4*dx04*dy04 + xlam5*dx05*dy05 + xlam6*dx06*dy06);
     p13[igrp] = (xlam1*dx01*dz01 + xlam2*dx02*dz02 + xlam3*dx03*dz03
           + xlam4*dx04*dz04 + xlam5*dx05*dz05 + xlam6*dx06*dz06);
     p23[igrp] = (xlam1*dy01*dz01 + xlam2*dy02*dz02 + xlam3*dy03*dz03
           + xlam4*dy04*dz04 + xlam5*dy05*dz05 + xlam6*dy06*dz06);
       
 }/*end for*/
  

 
#ifndef NO_PRAGMA
#pragma IVDEP
#endif
 for(igrp=1;igrp <= ngrp; igrp++) {
    pvten_tmp[1] -= (p11[igrp]*pnorm);
    pvten_tmp[2] -= (p12[igrp]*pnorm);
    pvten_tmp[3] -= (p13[igrp]*pnorm);
    pvten_tmp[4] -= (p12[igrp]*pnorm);
    pvten_tmp[5] -= (p22[igrp]*pnorm);
    pvten_tmp[6] -= (p23[igrp]*pnorm);
    pvten_tmp[7] -= (p13[igrp]*pnorm);
    pvten_tmp[8] -= (p23[igrp]*pnorm);
    pvten_tmp[9] -= (p33[igrp]*pnorm);
  }/*end for*/

/* Save multiplier */
 for(igrp=1;igrp <= ngrp; igrp++) {
    grp_bond_con_al_46[1][igrp] = xlam[1][igrp];
    grp_bond_con_al_46[2][igrp] = xlam[2][igrp];
    grp_bond_con_al_46[3][igrp] = xlam[3][igrp];
    grp_bond_con_al_46[4][igrp] = xlam[4][igrp];
    grp_bond_con_al_46[5][igrp] = xlam[5][igrp];
    grp_bond_con_al_46[6][igrp] = xlam[6][igrp];
 } /* end for */
  

/*=======================================================================*/
/*  IV)Allreduce pvten_tmp     */

  if(np_forc > 1 ){
   for(i=1;i<=9;i++){
    pvten_tmp2[i] = pvten_tmp[i];
   }/*endfor*/
   Allreduce(&(pvten_tmp2[1]), &(pvten_tmp[1]),9,MPI_DOUBLE,
                   MPI_SUM,0,comm_forc);
  }/*endif*/

  pvten_inc[1] += pvten_tmp[1];
  pvten_inc[2] += pvten_tmp[2];
  pvten_inc[3] += pvten_tmp[3];
  pvten_inc[4] += pvten_tmp[4];
  pvten_inc[5] += pvten_tmp[5];
  pvten_inc[6] += pvten_tmp[6];
  pvten_inc[7] += pvten_tmp[7];
  pvten_inc[8] += pvten_tmp[8];
  pvten_inc[9] += pvten_tmp[9];
 
 if(ifirst == 0){
   constr_cell_mat(iperd,hmat_cons_typ,hmat_int_typ,pvten_tmp);
   for(i=1;i<=9;i++){
     fgmat_p[i]+= pvten_tmp[i];    
     vgmat[i]  += (pvten_tmp[i]*roll_scg*0.5*dt/mass_hm);
   }/*endfor*/
 }/*endif*/
/* free memory assigned */
 if(ngrp > 0){
   free_dvector(rm1,1,ngrp); free_dvector(rm2,1,ngrp);
   free_dvector(rm3,1,ngrp); free_dvector(rm4,1,ngrp);

   free_dvector(dlmax,1,6);
   free_dvector(txlam,1,6);

   free_d3tensor(rmm,1,6,1,6,1,ngrp);

   free_dmatrix(dx,1,6,1,ngrp);
   free_dmatrix(dy,1,6,1,ngrp);
   free_dmatrix(dz,1,6,1,ngrp);

   free_dmatrix(dxt,1,6,1,ngrp);
   free_dmatrix(dyt,1,6,1,ngrp);
   free_dmatrix(dzt,1,6,1,ngrp);

   free_dmatrix(dxn,1,6,1,ngrp);
   free_dmatrix(dyn,1,6,1,ngrp);
   free_dmatrix(dzn,1,6,1,ngrp);

   free_dmatrix(dx0,1,6,1,ngrp);
   free_dmatrix(dy0,1,6,1,ngrp);
   free_dmatrix(dz0,1,6,1,ngrp);

   free_dmatrix(dij,1,6,1,ngrp);
   free_dmatrix(dxl,1,6,1,ngrp);
   free_dmatrix(amat,1,36,1,ngrp);
   free_dvector(tamat,1,36);
   free_dmatrix(avec,1,6,1,ngrp);
   free_dmatrix(xlam,1,6,1,ngrp);

   free_dmatrix(x,1,4,1,ngrp);
   free_dmatrix(y,1,4,1,ngrp);
   free_dmatrix(z,1,4,1,ngrp);

   free_dmatrix(xo,1,4,1,ngrp);
   free_dmatrix(yo,1,4,1,ngrp);
   free_dmatrix(zo,1,4,1,ngrp);

   free(ind1); free(ind2); free(ind3); free(ind4);

   free_dvector(p11,1,ngrp);
   free_dvector(p12,1,ngrp);
   free_dvector(p13,1,ngrp);
   free_dvector(p22,1,ngrp);
   free_dvector(p23,1,ngrp);
   free_dvector(p33,1,ngrp);

 }/*endif*/

/*==========================================================================*/
} /* end routine */
/*==========================================================================*/



/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void rattle_46_rollf(GRP_BOND_CON *grp_bond_con,
                     CLATOMS_INFO *clatoms_info,CLATOMS_POS *clatoms_pos,
                     PTENS *ptens,double dt,
                     PAR_RAHMAN *par_rahman,int ifirst, CELL *cell,
                     CLASS_COMM_FORC_PKG *class_comm_forc_pkg)

/*==========================================================================*/
/*        Begin Routine                                                     */
   {/*Begin Routine*/
/*=======================================================================*/
/*         Local Variable declarations                                   */
#include "../typ_defs/typ_mask.h"

 double pnorm;
 double rms1,rms2,rms3,rms4;
 int i,j,k,iii,ktemp;
 int igrp,*ind1,*ind2,*ind3,*ind4,jtyp;
 int na,job,info,ipvt[NCON_46+1];       /* For dgefa and dgesl */

 double roll_sci,dlam;
 int n;
 int ktemp1,ktemp2,ktemp3,ktemp4;

 double *rm1,*rm2,*rm3,*rm4;
 double *p11,*p22,*p33,*p12,*p13,*p23;
 double **x,**y,**z,**vx,**vy,**vz;
 double **avec;
 double ***rmm;
 double **dx,**dy,**dz;
 double **dvx,**dvy,**dvz;
 double **amat,*tamat;
 double **xlam;
 double *txlam;
 double roll_mtvvi[10],rolli_by_vg[10];
 double det_roll_mtvv;

/* Local pointers */

  double *clatoms_roll_sc      = clatoms_info->roll_sc;
  double *clatoms_mass         = clatoms_info->mass;
  double *clatoms_x            = clatoms_pos->x;
  double *clatoms_y            = clatoms_pos->y;
  double *clatoms_z            = clatoms_pos->z;
  double *clatoms_vx           = clatoms_pos->vx;
  double *clatoms_vy           = clatoms_pos->vy;
  double *clatoms_vz           = clatoms_pos->vz;

  int ngrp                     = grp_bond_con->num_46;
  int *grp_bond_con_j1_46      = grp_bond_con->j1_46;
  int *grp_bond_con_j2_46      = grp_bond_con->j2_46;
  int *grp_bond_con_j3_46      = grp_bond_con->j3_46;
  int *grp_bond_con_j4_46      = grp_bond_con->j4_46;
  int *grp_bond_con_jtyp_46    = grp_bond_con->jtyp_46;

  double *pvten_inc            = ptens->pvten_inc;
  double *pvten_tmp            = ptens->pvten_tmp;
  double *pvten_tmp2           = ptens->pvten_tmp_res;

  double *roll_mtvv            = par_rahman->roll_mtvv;
  double *vgmat_g              = par_rahman->vgmat_g;
  double *fgmat_p              = par_rahman->fgmat_p;
  double roll_scg              = par_rahman->roll_scg;
  double mass_hm               = par_rahman->mass_hm;

  int iperd                    = cell->iperd;
  int hmat_cons_typ            = cell->hmat_cons_typ;
  int hmat_int_typ             = cell->hmat_int_typ;

  int np_forc                  = class_comm_forc_pkg->num_proc;
  int myid_forc                = class_comm_forc_pkg->myid;
  MPI_Comm comm_forc           = class_comm_forc_pkg->comm;

/*=======================================================================*/

  if(ngrp > 0){
    rm1= dvector(1,ngrp);
    rm2= dvector(1,ngrp);
    rm3= dvector(1,ngrp);
    rm4= dvector(1,ngrp);
    p11= dvector(1,ngrp);
    p12= dvector(1,ngrp);
    p13= dvector(1,ngrp);
    p22= dvector(1,ngrp);
    p23= dvector(1,ngrp);
    p33= dvector(1,ngrp);
      x= dmatrix(1,4,1,ngrp);
      y= dmatrix(1,4,1,ngrp);
      z= dmatrix(1,4,1,ngrp);
     vx= dmatrix(1,4,1,ngrp);
     vy= dmatrix(1,4,1,ngrp);
     vz= dmatrix(1,4,1,ngrp);
    rmm= d3tensor(1,6,1,6,1,ngrp);
     dx= dmatrix(1,6,1,ngrp);
     dy= dmatrix(1,6,1,ngrp);
     dz= dmatrix(1,6,1,ngrp);
    dvx= dmatrix(1,6,1,ngrp);
    dvy= dmatrix(1,6,1,ngrp);
    dvz= dmatrix(1,6,1,ngrp);
   amat= dmatrix(1,36,1,ngrp);
  tamat= (double *)cmalloc(36*sizeof(double))-1;
   avec= dmatrix(1,6,1,ngrp);
   xlam= dmatrix(1,6,1,ngrp);
  txlam= dvector(1,6);
  ind1= (int *)calloc((ngrp+1),sizeof(int));
  ind2= (int *)calloc((ngrp+1),sizeof(int));
  ind3= (int *)calloc((ngrp+1),sizeof(int));
  ind4= (int *)calloc((ngrp+1),sizeof(int));
 }/*endif*/

/*=======================================================================*/

/* Malloc up some vectors and matrices */
 na = NCON_46;
 pnorm = 2.0/dt;

 n=3;
 gethinv(roll_mtvv,roll_mtvvi,&det_roll_mtvv,n);
 matmul_tt(roll_mtvvi,vgmat_g,rolli_by_vg,n);

 pvten_tmp[1] = 0.0;
 pvten_tmp[2] = 0.0;
 pvten_tmp[3] = 0.0;
 pvten_tmp[4] = 0.0;
 pvten_tmp[5] = 0.0;
 pvten_tmp[6] = 0.0;
 pvten_tmp[7] = 0.0;
 pvten_tmp[8] = 0.0;
 pvten_tmp[9] = 0.0;


 for(igrp=1;igrp <= ngrp; igrp++) {
    ind1[igrp] = grp_bond_con_j1_46[igrp];
    ind2[igrp] = grp_bond_con_j2_46[igrp];
    ind3[igrp] = grp_bond_con_j3_46[igrp];
    ind4[igrp] = grp_bond_con_j4_46[igrp];
}/*end for*/

 for(igrp=1;igrp <= ngrp; igrp++) {
  rm1[igrp] = 1.0/clatoms_mass[ind1[igrp]];
  rm2[igrp] = 1.0/clatoms_mass[ind2[igrp]];
  rm3[igrp] = 1.0/clatoms_mass[ind3[igrp]];
  rm4[igrp] = 1.0/clatoms_mass[ind4[igrp]];
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
  ktemp = ind4[igrp];
  x[4][igrp] = clatoms_x[ktemp];
  y[4][igrp] = clatoms_y[ktemp];
  z[4][igrp] = clatoms_z[ktemp];
 }



 for(igrp=1;igrp <= ngrp; igrp++) {
    ktemp= ind1[igrp];
    ktemp3= ind3[igrp];
    roll_sci=1.0/clatoms_roll_sc[ktemp3];/*all roll scales the same in same cons*/
    vx[1][igrp] = clatoms_vx[ktemp]
        + (x[1][igrp]*rolli_by_vg[1] +y[1][igrp]*rolli_by_vg[2] 
         + z[1][igrp]*rolli_by_vg[3]) *roll_sci;
    vy[1][igrp] = clatoms_vy[ktemp]
        + (x[1][igrp]*rolli_by_vg[4] +y[1][igrp]*rolli_by_vg[5] 
         + z[1][igrp]*rolli_by_vg[6]) *roll_sci;
    vz[1][igrp] = clatoms_vz[ktemp]
        + (x[1][igrp]*rolli_by_vg[7] +y[1][igrp]*rolli_by_vg[8]
         + z[1][igrp]*rolli_by_vg[9]) *roll_sci;
  }/*end for*/

 for(igrp=1;igrp <= ngrp; igrp++) {
    ktemp= ind2[igrp];
    ktemp3= ind3[igrp];
    roll_sci=1.0/clatoms_roll_sc[ktemp3];/*all roll scales the same in same cons*/
    vx[2][igrp] = clatoms_vx[ktemp]
        + (x[2][igrp]*rolli_by_vg[1] +y[2][igrp]*rolli_by_vg[2] 
         + z[2][igrp]*rolli_by_vg[3]) *roll_sci;
    vy[2][igrp] = clatoms_vy[ktemp]
        + (x[2][igrp]*rolli_by_vg[4] +y[2][igrp]*rolli_by_vg[5] 
         + z[2][igrp]*rolli_by_vg[6]) *roll_sci;
    vz[2][igrp] = clatoms_vz[ktemp]
        + (x[2][igrp]*rolli_by_vg[7] +y[2][igrp]*rolli_by_vg[8] 
         + z[2][igrp]*rolli_by_vg[9]) *roll_sci;
  }/*end for*/

 for(igrp=1;igrp <= ngrp; igrp++) {
    ktemp= ind3[igrp];
    roll_sci=1.0/clatoms_roll_sc[ktemp];/*all roll scales the same in same cons*/
    vx[3][igrp] = clatoms_vx[ktemp]
        + (x[3][igrp]*rolli_by_vg[1] +y[3][igrp]*rolli_by_vg[2] 
         + z[3][igrp]*rolli_by_vg[3]) *roll_sci;
    vy[3][igrp] = clatoms_vy[ktemp]
        + (x[3][igrp]*rolli_by_vg[4] +y[3][igrp]*rolli_by_vg[5] 
        + z[3][igrp]*rolli_by_vg[6]) *roll_sci;
    vz[3][igrp] = clatoms_vz[ktemp]
        + (x[3][igrp]*rolli_by_vg[7] +y[3][igrp]*rolli_by_vg[8] 
        + z[3][igrp]*rolli_by_vg[9]) *roll_sci;
  }/*end for*/

 for(igrp=1;igrp <= ngrp; igrp++) {
    ktemp= ind4[igrp];
    ktemp3= ind3[igrp];
    roll_sci=1.0/clatoms_roll_sc[ktemp3];/*all roll scales the same in same cons*/
    vx[4][igrp] = clatoms_vx[ktemp]
        + (x[4][igrp]*rolli_by_vg[1] +y[4][igrp]*rolli_by_vg[2] 
        + z[4][igrp]*rolli_by_vg[3]) *roll_sci;
    vy[4][igrp] = clatoms_vy[ktemp]
        + (x[4][igrp]*rolli_by_vg[4] +y[4][igrp]*rolli_by_vg[5] 
        + z[4][igrp]*rolli_by_vg[6]) *roll_sci;
    vz[4][igrp] = clatoms_vz[ktemp]
        + (x[4][igrp]*rolli_by_vg[7] +y[4][igrp]*rolli_by_vg[8] 
         + z[4][igrp]*rolli_by_vg[9]) *roll_sci;
  }/*end for*/

/* Set reciprocal mass matrix */
 for(igrp=1;igrp <= ngrp; igrp++) {
   rms1= rm1[igrp];
   rms2= rm2[igrp];
   rms3= rm3[igrp];
   rms4= rm4[igrp];

  rmm[1][1][igrp] = -(rms1+rms2); rmm[1][2][igrp] = rmm[1][3][igrp] = -rms1;
  rmm[1][4][igrp] = rmm[1][5][igrp] = rms2; rmm[1][6][igrp] = 0.0;

  rmm[2][1][igrp] = -rms1; rmm[2][2][igrp] = -(rms1+rms3); rmm[2][3][igrp] = -rms1;
  rmm[2][4][igrp] = -rms3; rmm[2][5][igrp] = 0.0; rmm[2][6][igrp] = rms3;

  rmm[3][1][igrp] = rmm[3][2][igrp] = -rms1; rmm[3][3][igrp] = -(rms1+rms4);
  rmm[3][4][igrp] = 0.0; rmm[3][5][igrp] = rmm[3][6][igrp] = -rms4;

  rmm[4][1][igrp] = rms2; rmm[4][2][igrp] = -rms3; rmm[4][3][igrp] = 0.0;
  rmm[4][4][igrp] = -(rms2+rms3); rmm[4][5][igrp] = -rms2; rmm[4][6][igrp] = rms3;

  rmm[5][1][igrp] = rms2; rmm[5][2][igrp] = 0.0; rmm[5][3][igrp] = -rms4;
  rmm[5][4][igrp] = -rms2; rmm[5][5][igrp] = -(rms2+rms4); rmm[5][6][igrp] = -rms4;

  rmm[6][1][igrp] = 0.0; rmm[6][2][igrp] = rms3; rmm[6][3][igrp] = -rms4;
  rmm[6][4][igrp] = rms3; rmm[6][5][igrp] = -rms4; rmm[6][6][igrp] = -(rms3+rms4);

 } /*end for*/

/* Compute difference vectors */

 for(igrp=1;igrp <= ngrp; igrp++) {
    dvx[1][igrp]= vx[1][igrp] - vx[2][igrp];
    dvx[2][igrp]= vx[1][igrp] - vx[3][igrp];
    dvx[3][igrp]= vx[1][igrp] - vx[4][igrp];
    dvx[4][igrp]= vx[2][igrp] - vx[3][igrp];
    dvx[5][igrp]= vx[2][igrp] - vx[4][igrp];
    dvx[6][igrp]= vx[3][igrp] - vx[4][igrp];
 }/*end for*/

 for(igrp=1;igrp <= ngrp; igrp++) {
    dvy[1][igrp]= vy[1][igrp] - vy[2][igrp];
    dvy[2][igrp]= vy[1][igrp] - vy[3][igrp];
    dvy[3][igrp]= vy[1][igrp] - vy[4][igrp];
    dvy[4][igrp]= vy[2][igrp] - vy[3][igrp];
    dvy[5][igrp]= vy[2][igrp] - vy[4][igrp];
    dvy[6][igrp]= vy[3][igrp] - vy[4][igrp];
 }/*end for*/

 for(igrp=1;igrp <= ngrp; igrp++) {
    dvz[1][igrp]= vz[1][igrp] - vz[2][igrp];
    dvz[2][igrp]= vz[1][igrp] - vz[3][igrp];
    dvz[3][igrp]= vz[1][igrp] - vz[4][igrp];
    dvz[4][igrp]= vz[2][igrp] - vz[3][igrp];
    dvz[5][igrp]= vz[2][igrp] - vz[4][igrp];
    dvz[6][igrp]= vz[3][igrp] - vz[4][igrp];
 }/*end for*/

 for(igrp=1;igrp <= ngrp; igrp++) {
    dx[1][igrp]= x[1][igrp] - x[2][igrp];
    dx[2][igrp]= x[1][igrp] - x[3][igrp];
    dx[3][igrp]= x[1][igrp] - x[4][igrp];
    dx[4][igrp]= x[2][igrp] - x[3][igrp];
    dx[5][igrp]= x[2][igrp] - x[4][igrp];
    dx[6][igrp]= x[3][igrp] - x[4][igrp];
 }/*end for*/

 for(igrp=1;igrp <= ngrp; igrp++) {
    dy[1][igrp]= y[1][igrp] - y[2][igrp];
    dy[2][igrp]= y[1][igrp] - y[3][igrp];
    dy[3][igrp]= y[1][igrp] - y[4][igrp];
    dy[4][igrp]= y[2][igrp] - y[3][igrp];
    dy[5][igrp]= y[2][igrp] - y[4][igrp];
    dy[6][igrp]= y[3][igrp] - y[4][igrp];
 }/*end for*/

 for(igrp=1;igrp <= ngrp; igrp++) {
    dz[1][igrp]= z[1][igrp] - z[2][igrp];
    dz[2][igrp]= z[1][igrp] - z[3][igrp];
    dz[3][igrp]= z[1][igrp] - z[4][igrp];
    dz[4][igrp]= z[2][igrp] - z[3][igrp];
    dz[5][igrp]= z[2][igrp] - z[4][igrp];
    dz[6][igrp]= z[3][igrp] - z[4][igrp];
 }

/* Get initial guess for lambda */

  iii = 0;
  for(i=1; i <= NCON_46; i++){
   for(j=1; j <= NCON_46; j++){
    iii++;
    for(igrp=1;igrp <= ngrp; igrp++) {
      amat[iii][igrp] = -rmm[i][j][igrp]*
               (dx[i][igrp]*dx[j][igrp] + dy[i][igrp]*dy[j][igrp] + dz[i][igrp]*dz[j][igrp]);
    }/*endfor*/
   }/*endfor*/
  }/*endfor*/

for(igrp=1;igrp <= ngrp; igrp++) {

   avec[1][igrp] = dvx[1][igrp]*dx[1][igrp] + dvy[1][igrp]*dy[1][igrp]
                   + dvz[1][igrp]*dz[1][igrp];
   avec[2][igrp] = dvx[2][igrp]*dx[2][igrp] + dvy[2][igrp]*dy[2][igrp]
                   + dvz[2][igrp]*dz[2][igrp];
   avec[3][igrp] = dvx[3][igrp]*dx[3][igrp] + dvy[3][igrp]*dy[3][igrp]
                   + dvz[3][igrp]*dz[3][igrp];
   avec[4][igrp] = dvx[4][igrp]*dx[4][igrp] + dvy[4][igrp]*dy[4][igrp]
                   + dvz[4][igrp]*dz[4][igrp];
   avec[5][igrp] = dvx[5][igrp]*dx[5][igrp] + dvy[5][igrp]*dy[5][igrp]
                   + dvz[5][igrp]*dz[5][igrp];
   avec[6][igrp] = dvx[6][igrp]*dx[6][igrp] + dvy[6][igrp]*dy[6][igrp]
                   + dvz[6][igrp]*dz[6][igrp];
 }
for(igrp=1;igrp <= ngrp; igrp++) {
   txlam[1]= xlam[1][igrp] = avec[1][igrp];
   txlam[2]= xlam[2][igrp] = avec[2][igrp];
   txlam[3]= xlam[3][igrp] = avec[3][igrp];
   txlam[4]= xlam[4][igrp] = avec[4][igrp];
   txlam[5]= xlam[5][igrp] = avec[5][igrp];
   txlam[6]= xlam[6][igrp] = avec[6][igrp];

   for(i=1; i<=36; i++){
    tamat[i] = amat[i][igrp];
   }


/* Solve linear system A xlam = avec */

  for(i=1; i <=na;i++) ipvt[i] = 0;
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

 }/*endfor*/


/* Velocity update */
#ifndef NO_PRAGMA
#pragma IVDEP
#endif
 for(igrp=1;igrp <= ngrp; igrp++) {
   double xlam1, xlam2, xlam3, xlam4, xlam5, xlam6;
   double dx1,dx2,dx3,dx4,dx5,dx6;
   double dy1,dy2,dy3,dy4,dy5,dy6;
   double dz1,dz2,dz3,dz4,dz5,dz6;
 
     ktemp1= ind1[igrp];
     ktemp2= ind2[igrp];
     ktemp3= ind3[igrp];
     ktemp4= ind4[igrp];

    xlam1= xlam[1][igrp]; xlam2= xlam[2][igrp];
    xlam3= xlam[3][igrp]; xlam4= xlam[4][igrp];
    xlam5= xlam[5][igrp]; xlam6= xlam[6][igrp];

    dx1= dx[1][igrp]; dx2= dx[2][igrp]; dx3= dx[3][igrp];
    dx4= dx[4][igrp]; dx5= dx[5][igrp]; dx6= dx[6][igrp];

    dy1= dy[1][igrp]; dy2= dy[2][igrp]; dy3= dy[3][igrp];
    dy4= dy[4][igrp]; dy5= dy[5][igrp]; dy6= dy[6][igrp];

    dz1= dz[1][igrp]; dz2= dz[2][igrp]; dz3= dz[3][igrp];
    dz4= dz[4][igrp]; dz5= dz[5][igrp]; dz6= dz[6][igrp];


  clatoms_vx[ktemp1] -=  ( xlam1*dx1 + xlam2*dx2 + xlam3*dx3)*rm1[igrp];
  clatoms_vy[ktemp1] -=  ( xlam1*dy1 + xlam2*dy2 + xlam3*dy3)*rm1[igrp];
  clatoms_vz[ktemp1] -=  ( xlam1*dz1 + xlam2*dz2 + xlam3*dz3)*rm1[igrp];

  clatoms_vx[ktemp2] -=  (-xlam1*dx1 + xlam4*dx4 + xlam5*dx5)*rm2[igrp];
  clatoms_vy[ktemp2] -=  (-xlam1*dy1 + xlam4*dy4 + xlam5*dy5)*rm2[igrp];
  clatoms_vz[ktemp2] -=  (-xlam1*dz1 + xlam4*dz4 + xlam5*dz5)*rm2[igrp];

  clatoms_vx[ktemp3] -=  (-xlam2*dx2 - xlam4*dx4 + xlam6*dx6)*rm3[igrp];
  clatoms_vy[ktemp3] -=  (-xlam2*dy2 - xlam4*dy4 + xlam6*dy6)*rm3[igrp];
  clatoms_vz[ktemp3] -=  (-xlam2*dz2 - xlam4*dz4 + xlam6*dz6)*rm3[igrp];

  clatoms_vx[ktemp4] -=  (-xlam3*dx3 - xlam5*dx5 - xlam6*dx6)*rm4[igrp];
  clatoms_vy[ktemp4] -=  (-xlam3*dy3 - xlam5*dy5 - xlam6*dy6)*rm4[igrp];
  clatoms_vz[ktemp4] -=  (-xlam3*dz3 - xlam5*dz5 - xlam6*dz6)*rm4[igrp];

/* Pressure tensor update */

     p11[igrp]= xlam1*dx1*dx1 +xlam2*dx2*dx2 +xlam3*dx3*dx3
               +xlam4*dx4*dx4 +xlam5*dx5*dx5 +xlam6*dx6*dx6;

     p22[igrp]= xlam1*dy1*dy1 +xlam2*dy2*dy2 +xlam3*dy3*dy3
               +xlam4*dy4*dy4 +xlam5*dy5*dy5 +xlam6*dy6*dy6;

     p33[igrp]= xlam1*dz1*dz1 +xlam2*dz2*dz2 +xlam3*dz3*dz3
               +xlam4*dz4*dz4 +xlam5*dz5*dz5 +xlam6*dz6*dz6;

     p12[igrp]= xlam1*dx1*dy1 +xlam2*dx2*dy2 +xlam3*dx3*dy3
               +xlam4*dx4*dy4 +xlam5*dx5*dy5 +xlam6*dx6*dy6;

     p13[igrp]= xlam1*dx1*dz1 +xlam2*dx2*dz2 +xlam3*dx3*dz3
               +xlam4*dx4*dz4 +xlam5*dx5*dz5 +xlam6*dx6*dz6;

     p23[igrp]= xlam1*dy1*dz1 +xlam2*dy2*dz2 +xlam3*dy3*dz3
               +xlam4*dy4*dz4 +xlam5*dy5*dz5 +xlam6*dy6*dz6;
  }/*endfor*/

#ifndef NO_PRAGMA
#pragma IVDEP
#endif
 for(igrp=1;igrp <= ngrp; igrp++) {
    pvten_tmp[1] -= (p11[igrp]*pnorm);
    pvten_tmp[2] -= (p12[igrp]*pnorm);
    pvten_tmp[3] -= (p13[igrp]*pnorm);
    pvten_tmp[4] -= (p12[igrp]*pnorm);
    pvten_tmp[5] -= (p22[igrp]*pnorm);
    pvten_tmp[6] -= (p23[igrp]*pnorm);
    pvten_tmp[7] -= (p13[igrp]*pnorm);
    pvten_tmp[8] -= (p23[igrp]*pnorm);
    pvten_tmp[9] -= (p33[igrp]*pnorm);

 } /* end for igrp */

/*=======================================================================*/
/*  IV)Allreduce pvten_tmp     */

  if(np_forc > 1 ){
   for(i=1;i<=9;i++){
    pvten_tmp2[i] = pvten_tmp[i];
   }/*endfor*/
   Allreduce(&(pvten_tmp2[1]), &(pvten_tmp[1]),9,MPI_DOUBLE,
                   MPI_SUM,0,comm_forc);
  }/*endif*/

 for(i=1; i<=9; i++){ 
  pvten_inc[i] += pvten_tmp[i];
 }

 if(ifirst == 0){
   constr_cell_mat(iperd,hmat_cons_typ,hmat_int_typ,pvten_tmp);
   for(i=1;i<=9;i++){      
     fgmat_p[i] += pvten_tmp[i];
     vgmat_g[i] += (pvten_tmp[i]*roll_scg*0.5*dt/mass_hm);
  }/*endfor*/
 }/*endif*/
/* free locally assigned memmory */
 if(ngrp > 0){
    free_dvector(rm1,1,ngrp);
    free_dvector(rm2,1,ngrp);
    free_dvector(rm3,1,ngrp);
    free_dvector(rm4,1,ngrp);

    free_dvector(p11,1,ngrp);
    free_dvector(p12,1,ngrp);
    free_dvector(p13,1,ngrp);
    free_dvector(p22,1,ngrp);
    free_dvector(p23,1,ngrp);
    free_dvector(p33,1,ngrp);

    free_dmatrix(x,1,4,1,ngrp);
    free_dmatrix(y,1,4,1,ngrp);
    free_dmatrix(z,1,4,1,ngrp);

    free_dmatrix(vx,1,4,1,ngrp);
    free_dmatrix(vy,1,4,1,ngrp);
    free_dmatrix(vz,1,4,1,ngrp);

    free_d3tensor(rmm,1,6,1,6,1,ngrp);

    free_dmatrix(dx,1,4,1,ngrp);
    free_dmatrix(dy,1,4,1,ngrp);
    free_dmatrix(dz,1,4,1,ngrp);

    free_dmatrix(dvx,1,4,1,ngrp);
    free_dmatrix(dvy,1,4,1,ngrp);
    free_dmatrix(dvz,1,4,1,ngrp);

    free_dmatrix(amat,1,36,1,ngrp);
    cfree(&(tamat[1]));
    free_dmatrix(avec,1,6,1,ngrp);
    free_dmatrix(xlam,1,6,1,ngrp);
    free_dvector(txlam,1,6);


    free(ind1); free(ind2);
    free(ind3); free(ind4);
 }/*endif*/

/*=======================================================================*/
} /* end routine */
/*=======================================================================*/

