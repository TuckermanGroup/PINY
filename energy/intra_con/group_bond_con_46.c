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

#define NCON_46 6
#define NAT_46 4


/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
void shake_46(GRP_BOND_CON *grp_bond_con,
              CLATOMS_INFO *clatoms_info,CLATOMS_POS *clatoms_pos,
              PTENS *ptens,double dt,double *aiter,
              CLASS_COMM_FORC_PKG *class_comm_forc_pkg)

/*==========================================================================*/
/*        Begin Routine                                                     */
{/*Begin Routine*/
/*=======================================================================*/
/*         Local Variable declarations                                   */
  
 double txlam[NCON_46+1],xl0[NCON_46+1],dlmax[NCON_46+1],dmax;
 double rms1,rms2,rms3,rms4;
 double dts;
 int i,j,iii,ktemp;
 int iter,igrp,*ind1,*ind2,*ind3,*ind4,jtyp;
 int na,job,info,ipvt[NCON_46+1];       /* For dgefa and dgesl */


  double *rm1,*rm2,*rm3,*rm4;
  double **xlam, **avec, **amat, **dxl,**dij;
  double *tamat;
  double **dx,**dy,**dz,**dxt,**dyt,**dzt,**dxn,**dyn,**dzn; 
  double **x,**y,**z,**xo,**yo,**zo;
  double *p11,*p22,*p33,*p12,*p13,*p23;
  double ***rmm;
  
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
  double pnorm;

  int ngrp,irem,igrp_off;
  int ngrp_tot = grp_bond_con->num_46;
  int np_forc                  = class_comm_forc_pkg->num_proc;
  int myid_forc                = class_comm_forc_pkg->myid;

/*=======================================================================*/

  ngrp = ngrp_tot;
  igrp_off = 0;

/*=======================================================================*/

/* Assign local memory   */
       rmm= d3tensor(1,6,1,6,1,ngrp);
      xlam= dmatrix(1,6,1,ngrp);
      avec= dmatrix(1,6,1,ngrp);
      amat= dmatrix(1,36,1,ngrp);
     tamat= (double *)cmalloc(36*sizeof(double))-1;
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
       rm1= dvector(1,ngrp); 
       rm2= dvector(1,ngrp); 
       rm3= dvector(1,ngrp); 
       rm4= dvector(1,ngrp); 
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
      ind1 = (int *) cmalloc(ngrp*sizeof(int))-1;
      ind2 = (int *) cmalloc(ngrp*sizeof(int))-1;
      ind3 = (int *) cmalloc(ngrp*sizeof(int))-1;
      ind4 = (int *) cmalloc(ngrp*sizeof(int))-1;
  
  

/*=======================================================================*/ 
/* Malloc up some vectors and matrices */
 na = NCON_46;  /* = 6 */


 dts = dt*dt;
 pnorm = 2.0/dts;
 *aiter = 0.0;
/* I collect masses and positions */
 for(igrp=1;igrp <= ngrp; igrp++) {
    ind1[igrp] = grp_bond_con_j1_46[(igrp+igrp_off)];
    ind2[igrp] = grp_bond_con_j2_46[(igrp+igrp_off)];
    ind3[igrp] = grp_bond_con_j3_46[(igrp+igrp_off)];
    ind4[igrp] = grp_bond_con_j4_46[(igrp+igrp_off)];
  }

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
    jtyp = grp_bond_con_jtyp_46[(igrp+igrp_off)];

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

     dx[1][igrp]= xo[1][igrp]-xo[2][igrp];
     dx[2][igrp]= xo[1][igrp]-xo[3][igrp];
     dx[3][igrp]= xo[1][igrp]-xo[4][igrp];
     dx[4][igrp]= xo[2][igrp]-xo[3][igrp];
     dx[5][igrp]= xo[2][igrp]-xo[4][igrp];
     dx[6][igrp]= xo[3][igrp]-xo[4][igrp];
   }
 for(igrp=1;igrp <= ngrp; igrp++) {

     dy[1][igrp]= yo[1][igrp]-yo[2][igrp];
     dy[2][igrp]= yo[1][igrp]-yo[3][igrp];
     dy[3][igrp]= yo[1][igrp]-yo[4][igrp];
     dy[4][igrp]= yo[2][igrp]-yo[3][igrp];
     dy[5][igrp]= yo[2][igrp]-yo[4][igrp];
     dy[6][igrp]= yo[3][igrp]-yo[4][igrp];
   }

 for(igrp=1;igrp <= ngrp; igrp++) {
     dz[1][igrp]= zo[1][igrp]-zo[2][igrp];
     dz[2][igrp]= zo[1][igrp]-zo[3][igrp];
     dz[3][igrp]= zo[1][igrp]-zo[4][igrp];
     dz[4][igrp]= zo[2][igrp]-zo[3][igrp];
     dz[5][igrp]= zo[2][igrp]-zo[4][igrp];
     dz[6][igrp]= zo[3][igrp]-zo[4][igrp];
   }


/* Solve for initial lambda */

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
  for(i=1; i <= NCON_46; i++){
      for(igrp=1;igrp <= ngrp; igrp++) {
          avec[i][igrp] = dij[i][igrp]*dij[i][igrp] - 
                         (dxt[i][igrp]*dxt[i][igrp]+ dyt[i][igrp]*dyt[i][igrp]                               + dzt[i][igrp]*dzt[i][igrp]);
            xlam[i][igrp] = avec[i][igrp];
       }/*endfor*/
  }/*endfor*/

/* Solve linear system A xlam = avec */

 for(igrp=1;igrp <= ngrp; igrp++) {
      txlam[1]=xlam[1][igrp]; txlam[2]=xlam[2][igrp];
      txlam[3]=xlam[3][igrp]; txlam[4]=xlam[4][igrp];
      txlam[5]=xlam[5][igrp]; txlam[6]=xlam[6][igrp];

   for(i=1; i<=36; i++){
    tamat[i] = amat[i][igrp];
   }

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
      xlam[1][igrp]=txlam[1]; xlam[2][igrp]=txlam[2];
      xlam[3][igrp]=txlam[3]; xlam[4][igrp]=txlam[4];
      xlam[5][igrp]=txlam[5]; xlam[6][igrp]=txlam[6];
} /*end loop over group of atoms */
       
/* Iterative loop to convergence */

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
/* Set up guess of difference vectors */
 for(igrp=1;igrp <= ngrp; igrp++) {

   dxn[1][igrp] = 2.0*dxt[1][igrp]; dyn[1][igrp] = 2.0*dyt[1][igrp]; 
   dxn[2][igrp] = 2.0*dxt[2][igrp]; dyn[2][igrp] = 2.0*dyt[2][igrp];
   dxn[3][igrp] = 2.0*dxt[3][igrp]; dyn[3][igrp] = 2.0*dyt[3][igrp];

   dxn[4][igrp] = 2.0*dxt[4][igrp]; dyn[4][igrp] = 2.0*dyt[4][igrp]; 
   dxn[5][igrp] = 2.0*dxt[5][igrp]; dyn[5][igrp] = 2.0*dyt[5][igrp]; 
   dxn[6][igrp] = 2.0*dxt[6][igrp]; dyn[6][igrp] = 2.0*dyt[6][igrp]; 

   dzn[1][igrp] = 2.0*dzt[1][igrp];
   dzn[2][igrp] = 2.0*dzt[2][igrp];
   dzn[3][igrp] = 2.0*dzt[3][igrp];
   dzn[4][igrp] = 2.0*dzt[4][igrp];
   dzn[5][igrp] = 2.0*dzt[5][igrp];
   dzn[6][igrp] = 2.0*dzt[6][igrp];
 }

   for(i=1; i <= NCON_46; i++) {
     for(j=1; j <= NCON_46; j++) {
       for(igrp=1;igrp <= ngrp; igrp++) {
         dxn[i][igrp] += rmm[i][j][igrp]*xlam[j][igrp]*dx[j][igrp];
         dyn[i][igrp] += rmm[i][j][igrp]*xlam[j][igrp]*dy[j][igrp];
         dzn[i][igrp] += rmm[i][j][igrp]*xlam[j][igrp]*dz[j][igrp];
       }/*endfor*/
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
                +dzn[i][igrp]*dz[j][igrp]);
           }/*endfor*/
        }/*endfor*/
   }/*endfor*/

 for(igrp=1;igrp <= ngrp; igrp++) {

    xl0[1] = xlam[1][igrp]; xl0[2] = xlam[2][igrp];
    xl0[3] = xlam[3][igrp]; xl0[4] = xlam[4][igrp];
    xl0[5] = xlam[5][igrp]; xl0[6] = xlam[6][igrp];

    xlam[1][igrp] = avec[1][igrp]; xlam[2][igrp] = avec[2][igrp];
    xlam[3][igrp] = avec[3][igrp]; xlam[4][igrp] = avec[4][igrp];
    xlam[5][igrp] = avec[5][igrp]; xlam[6][igrp] = avec[6][igrp];

    txlam[1]=xlam[1][igrp]; txlam[2]=xlam[2][igrp];
    txlam[3]=xlam[3][igrp]; txlam[4]=xlam[4][igrp];
    txlam[5]=xlam[5][igrp]; txlam[6]=xlam[6][igrp];

    for(i=1; i<=36; i++){
     tamat[i] = amat[i][igrp];
    }

  ipvt[1]=0; ipvt[2]=0; ipvt[3]=0;
  ipvt[4]=0; ipvt[5]=0; ipvt[6]=0;

/* Solve linear system A xlam = avec */
   
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

    xlam[1][igrp]=txlam[1]; xlam[2][igrp]=txlam[2];
    xlam[3][igrp]=txlam[3]; xlam[4][igrp]=txlam[4];
    xlam[5][igrp]=txlam[5]; xlam[6][igrp]=txlam[6];

    dxl[1][igrp] = fabs(xlam[1][igrp]-xl0[1]);
    dxl[2][igrp] = fabs(xlam[2][igrp]-xl0[2]);
    dxl[3][igrp] = fabs(xlam[3][igrp]-xl0[3]);
    dxl[4][igrp] = fabs(xlam[4][igrp]-xl0[4]);
    dxl[5][igrp] = fabs(xlam[5][igrp]-xl0[5]);
    dxl[6][igrp] = fabs(xlam[6][igrp]-xl0[6]);
   
 } /*endfor */


/* test  for convergence over all the groups */
/* assign initial values for dlmaxs */
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
  }/*endfor*/
        dmax= dlmax[1];
     for(i=2;i <= NCON_46; i++) {
        dmax = (dmax > dlmax[i] ? dmax:dlmax[i]);
     }/*endfor*/
} while(dmax > grp_bond_con->tol);
  *aiter += (double) iter;

/* Position update */

#ifndef NO_PRAGMA
#pragma IVDEP
#endif
 for(igrp=1;igrp <= ngrp; igrp++) {
   double dx1,dx2,dx3,dx4,dx5,dx6;
   double dy1,dy2,dy3,dy4,dy5,dy6;
   double dz1,dz2,dz3,dz4,dz5,dz6;
   double xlam1,xlam2,xlam3,xlam4,xlam5,xlam6;
   double t1,t2,t3,t4,t5,t6;
   int ktemp1,ktemp2,ktemp3,ktemp4;

     xlam1=xlam[1][igrp];
     xlam2=xlam[2][igrp];
     xlam3=xlam[3][igrp];
     xlam4=xlam[4][igrp];
     xlam5=xlam[5][igrp];
     xlam6=xlam[6][igrp];
    ktemp1=ind1[igrp]; ktemp2=ind2[igrp]; 
    ktemp3=ind3[igrp]; ktemp4=ind4[igrp];

    dx1= dx[1][igrp]; dx2= dx[2][igrp]; dx3= dx[3][igrp];
    dx4= dx[4][igrp]; dx5= dx[5][igrp]; dx6= dx[6][igrp];

    dy1= dy[1][igrp]; dy2= dy[2][igrp]; dy3= dy[3][igrp];
    dy4= dy[4][igrp]; dy5= dy[5][igrp]; dy6= dy[6][igrp];

    dz1= dz[1][igrp]; dz2= dz[2][igrp]; dz3= dz[3][igrp];
    dz4= dz[4][igrp]; dz5= dz[5][igrp]; dz6= dz[6][igrp];

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

  clatoms_vx[ktemp1] -=  ( xlam1*dx1 + xlam2*dx2 + xlam3*dx3)*rm1[igrp]/dt;
  clatoms_vy[ktemp1] -=  ( xlam1*dy1 + xlam2*dy2 + xlam3*dy3)*rm1[igrp]/dt;
  clatoms_vz[ktemp1] -=  ( xlam1*dz1 + xlam2*dz2 + xlam3*dz3)*rm1[igrp]/dt;

  clatoms_vx[ktemp2] -=  (-xlam1*dx1 + xlam4*dx4 + xlam5*dx5)*rm2[igrp]/dt;
  clatoms_vy[ktemp2] -=  (-xlam1*dy1 + xlam4*dy4 + xlam5*dy5)*rm2[igrp]/dt;
  clatoms_vz[ktemp2] -=  (-xlam1*dz1 + xlam4*dz4 + xlam5*dz5)*rm2[igrp]/dt;

  clatoms_vx[ktemp3] -=  (-xlam2*dx2 - xlam4*dx4 + xlam6*dx6)*rm3[igrp]/dt;
  clatoms_vy[ktemp3] -=  (-xlam2*dy2 - xlam4*dy4 + xlam6*dy6)*rm3[igrp]/dt;
  clatoms_vz[ktemp3] -=  (-xlam2*dz2 - xlam4*dz4 + xlam6*dz6)*rm3[igrp]/dt;

  clatoms_vx[ktemp4] -=  (-xlam3*dx3 - xlam5*dx5 - xlam6*dx6)*rm4[igrp]/dt;
  clatoms_vy[ktemp4] -=  (-xlam3*dy3 - xlam5*dy5 - xlam6*dy6)*rm4[igrp]/dt;
  clatoms_vz[ktemp4] -=  (-xlam3*dz3 - xlam5*dz5 - xlam6*dz6)*rm4[igrp]/dt;

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
    grp_bond_con_al_46[1][(igrp+igrp_off)] = xlam[1][igrp];
    grp_bond_con_al_46[2][(igrp+igrp_off)] = xlam[2][igrp];
    grp_bond_con_al_46[3][(igrp+igrp_off)] = xlam[3][igrp];
    grp_bond_con_al_46[4][(igrp+igrp_off)] = xlam[4][igrp];
    grp_bond_con_al_46[5][(igrp+igrp_off)] = xlam[5][igrp];
    grp_bond_con_al_46[6][(igrp+igrp_off)] = xlam[6][igrp];
 } /* end for */


/* free all the memory assigned in this subroutine */
   free_d3tensor(rmm,1,6,1,6,1,ngrp);
   free_dmatrix(xlam,1,6,1,ngrp);
   free_dmatrix(avec,1,6,1,ngrp);
   free_dmatrix(amat,1,36,1,ngrp);
   cfree(&(tamat[1]));
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
   free_dvector(rm1,1,ngrp);
   free_dvector(rm2,1,ngrp);
   free_dvector(rm3,1,ngrp);
   free_dvector(rm4,1,ngrp);

   free_dmatrix(x,1,4,1,ngrp);
   free_dmatrix(y,1,4,1,ngrp);
   free_dmatrix(z,1,4,1,ngrp);

   free_dmatrix(xo,1,4,1,ngrp); 
   free_dmatrix(yo,1,4,1,ngrp); 
   free_dmatrix(zo,1,4,1,ngrp); 

   free_dvector(p11,1,ngrp);
   free_dvector(p12,1,ngrp);
   free_dvector(p13,1,ngrp);
   free_dvector(p22,1,ngrp);
   free_dvector(p23,1,ngrp);
   free_dvector(p33,1,ngrp);

   cfree(&(ind1[1]));
   cfree(&(ind2[1]));
   cfree(&(ind3[1]));
   cfree(&(ind4[1]));

/*==========================================================================*/
} /* end routine */
/*==========================================================================*/



/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
void rattle_46(GRP_BOND_CON *grp_bond_con,
               CLATOMS_INFO *clatoms_info,CLATOMS_POS *clatoms_pos,
               PTENS *ptens,double dt,
              CLASS_COMM_FORC_PKG *class_comm_forc_pkg)

/*==========================================================================*/
/*        Begin Routine                                                     */
{/*Begin Routine*/
/*=======================================================================*/
/*         Local Variable declarations                                   */
  
 double pnorm;
 double rms1,rms2,rms3,rms4;
 int i,j,k,iii,ktemp;
 int igrp,*ind1,*ind2,*ind3,*ind4,jtyp;
 int na,job,info,ipvt[NCON_46+1];       /* For dgefa and dgesl */

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

  int ngrp,irem,igrp_off;
  int ngrp_tot = grp_bond_con->num_46;
  int np_forc                  = class_comm_forc_pkg->num_proc;
  int myid_forc                = class_comm_forc_pkg->myid;

/*=======================================================================*/

  ngrp = ngrp_tot;
  igrp_off = 0;

/*=======================================================================*/

/* assign local arrays */
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
  ind1 = (int *) cmalloc(ngrp*sizeof(int))-1;
  ind2 = (int *) cmalloc(ngrp*sizeof(int))-1;
  ind3 = (int *) cmalloc(ngrp*sizeof(int))-1;
  ind4 = (int *) cmalloc(ngrp*sizeof(int))-1;

/*=======================================================================*/

/* Malloc up some vectors and matrices */
 na = NCON_46;
 pnorm = 2.0/dt;

 for(igrp=1;igrp <= ngrp; igrp++) {
  ind1[igrp] = grp_bond_con_j1_46[(igrp+igrp_off)];
  ind2[igrp] = grp_bond_con_j2_46[(igrp+igrp_off)];
  ind3[igrp] = grp_bond_con_j3_46[(igrp+igrp_off)];
  ind4[igrp] = grp_bond_con_j4_46[(igrp+igrp_off)];
 }

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
  vx[1][igrp] = clatoms_vx[ktemp]; 
  vy[1][igrp] = clatoms_vy[ktemp]; 
  vz[1][igrp] = clatoms_vz[ktemp]; 
 }

 for(igrp=1;igrp <= ngrp; igrp++) {
  ktemp = ind2[igrp];
  vx[2][igrp] = clatoms_vx[ktemp];
  vy[2][igrp] = clatoms_vy[ktemp];
  vz[2][igrp] = clatoms_vz[ktemp];
 }

 for(igrp=1;igrp <= ngrp; igrp++) {
  ktemp = ind3[igrp];
  vx[3][igrp] = clatoms_vx[ktemp]; 
  vy[3][igrp] = clatoms_vy[ktemp]; 
  vz[3][igrp] = clatoms_vz[ktemp]; 
 }

 for(igrp=1;igrp <= ngrp; igrp++) {
  ktemp = ind4[igrp];
  vx[4][igrp] = clatoms_vx[ktemp];
  vy[4][igrp] = clatoms_vy[ktemp];
  vz[4][igrp] = clatoms_vz[ktemp];
 }

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
    dvx[1][igrp] = vx[1][igrp]-vx[2][igrp]; dvx[2][igrp] = vx[1][igrp]-vx[3][igrp];
    dvx[3][igrp] = vx[1][igrp]-vx[4][igrp]; dvx[4][igrp] = vx[2][igrp]-vx[3][igrp];
    dvx[5][igrp] = vx[2][igrp]-vx[4][igrp]; dvx[6][igrp] = vx[3][igrp]-vx[4][igrp];

    dvy[1][igrp] = vy[1][igrp]-vy[2][igrp]; dvy[2][igrp] = vy[1][igrp]-vy[3][igrp];
    dvy[3][igrp] = vy[1][igrp]-vy[4][igrp]; dvy[4][igrp] = vy[2][igrp]-vy[3][igrp];
    dvy[5][igrp] = vy[2][igrp]-vy[4][igrp]; dvy[6][igrp] = vy[3][igrp]-vy[4][igrp];

    dvz[1][igrp] = vz[1][igrp]-vz[2][igrp]; dvz[2][igrp] = vz[1][igrp]-vz[3][igrp];
    dvz[3][igrp] = vz[1][igrp]-vz[4][igrp]; dvz[4][igrp] = vz[2][igrp]-vz[3][igrp];
    dvz[5][igrp] = vz[2][igrp]-vz[4][igrp]; dvz[6][igrp] = vz[3][igrp]-vz[4][igrp];

    dx[1][igrp] = x[1][igrp]-x[2][igrp]; dx[2][igrp] = x[1][igrp]-x[3][igrp];
    dx[3][igrp] = x[1][igrp]-x[4][igrp]; dx[4][igrp] = x[2][igrp]-x[3][igrp];
    dx[5][igrp] = x[2][igrp]-x[4][igrp]; dx[6][igrp] = x[3][igrp]-x[4][igrp];

    dy[1][igrp] = y[1][igrp]-y[2][igrp]; dy[2][igrp] = y[1][igrp]-y[3][igrp];
    dy[3][igrp] = y[1][igrp]-y[4][igrp]; dy[4][igrp] = y[2][igrp]-y[3][igrp];
    dy[5][igrp] = y[2][igrp]-y[4][igrp]; dy[6][igrp] = y[3][igrp]-y[4][igrp];

    dz[1][igrp] = z[1][igrp]-z[2][igrp]; dz[2][igrp] = z[1][igrp]-z[3][igrp];
    dz[3][igrp] = z[1][igrp]-z[4][igrp]; dz[4][igrp] = z[2][igrp]-z[3][igrp];
    dz[5][igrp] = z[2][igrp]-z[4][igrp]; dz[6][igrp] = z[3][igrp]-z[4][igrp];

 }/*endfor*/

/* Get initial guess for lambda */
  iii = 0;
  for(i=1; i <= NCON_46; i++){
   for(j=1; j <= NCON_46; j++){
       iii++;
     for(igrp=1;igrp <= ngrp; igrp++) {
        amat[iii][igrp] =-rmm[i][j][igrp]*
               (dx[i][igrp]*dx[j][igrp] + dy[i][igrp]*dy[j][igrp]                                               + dz[i][igrp]*dz[j][igrp]);
     }/*endfor*/
   }/*endfor*/
  }/*endfor*/

for(igrp=1;igrp <= ngrp; igrp++) {

   avec[1][igrp] = dvx[1][igrp]*dx[1][igrp] + dvy[1][igrp]*dy[1][igrp]                                            + dvz[1][igrp]*dz[1][igrp];
   avec[2][igrp] = dvx[2][igrp]*dx[2][igrp] + dvy[2][igrp]*dy[2][igrp]                                            + dvz[2][igrp]*dz[2][igrp];
   avec[3][igrp] = dvx[3][igrp]*dx[3][igrp] + dvy[3][igrp]*dy[3][igrp]                                            + dvz[3][igrp]*dz[3][igrp];
   avec[4][igrp] = dvx[4][igrp]*dx[4][igrp] + dvy[4][igrp]*dy[4][igrp]                                            + dvz[4][igrp]*dz[4][igrp];
   avec[5][igrp] = dvx[5][igrp]*dx[5][igrp] + dvy[5][igrp]*dy[5][igrp]                                            + dvz[5][igrp]*dz[5][igrp];
   avec[6][igrp] = dvx[6][igrp]*dx[6][igrp] + dvy[6][igrp]*dy[6][igrp]                                            + dvz[6][igrp]*dz[6][igrp];
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
   ipvt[1]= 0; ipvt[2]= 0; ipvt[3]= 0;
   ipvt[4]= 0; ipvt[5]= 0; ipvt[6]= 0;

#ifdef IBM_ESSL
  dgef(&(tamat[1]),&na,&na,&(ipvt[1]));
#else
  DGEFA(&(tamat[1]),&na,&na,&(ipvt[1]),&info);
#endif
  job = 1;
#ifdef IBM_ESSL
  job = 1;  /* changed from 0 */
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

   ind1[igrp] = grp_bond_con_j1_46[(igrp+igrp_off)];
   ind2[igrp] = grp_bond_con_j2_46[(igrp+igrp_off)];
   ind3[igrp] = grp_bond_con_j3_46[(igrp+igrp_off)];
   ind4[igrp] = grp_bond_con_j4_46[(igrp+igrp_off)];

     xlam1=xlam[1][igrp];
     xlam2=xlam[2][igrp];
     xlam3=xlam[3][igrp];
     xlam4=xlam[4][igrp];
     xlam5=xlam[5][igrp];
     xlam6=xlam[6][igrp];

    dx1= dx[1][igrp]; dx2= dx[2][igrp]; dx3= dx[3][igrp];
    dx4= dx[4][igrp]; dx5= dx[5][igrp]; dx6= dx[6][igrp];

    dy1= dy[1][igrp]; dy2= dy[2][igrp]; dy3= dy[3][igrp];
    dy4= dy[4][igrp]; dy5= dy[5][igrp]; dy6= dy[6][igrp];

    dz1= dz[1][igrp]; dz2= dz[2][igrp]; dz3= dz[3][igrp];
    dz4= dz[4][igrp]; dz5= dz[5][igrp]; dz6= dz[6][igrp];

  clatoms_vx[ind1[igrp]] -=  ( xlam1*dx1 + xlam2*dx2 + xlam3*dx3)*rm1[igrp];
  clatoms_vy[ind1[igrp]] -=  ( xlam1*dy1 + xlam2*dy2 + xlam3*dy3)*rm1[igrp];
  clatoms_vz[ind1[igrp]] -=  ( xlam1*dz1 + xlam2*dz2 + xlam3*dz3)*rm1[igrp];

  clatoms_vx[ind2[igrp]] -=  (-xlam1*dx1 + xlam4*dx4 + xlam5*dx5)*rm2[igrp];
  clatoms_vy[ind2[igrp]] -=  (-xlam1*dy1 + xlam4*dy4 + xlam5*dy5)*rm2[igrp];
  clatoms_vz[ind2[igrp]] -=  (-xlam1*dz1 + xlam4*dz4 + xlam5*dz5)*rm2[igrp];

  clatoms_vx[ind3[igrp]] -=  (-xlam2*dx2 - xlam4*dx4 + xlam6*dx6)*rm3[igrp];
  clatoms_vy[ind3[igrp]] -=  (-xlam2*dy2 - xlam4*dy4 + xlam6*dy6)*rm3[igrp];
  clatoms_vz[ind3[igrp]] -=  (-xlam2*dz2 - xlam4*dz4 + xlam6*dz6)*rm3[igrp];

  clatoms_vx[ind4[igrp]] -=  (-xlam3*dx3 - xlam5*dx5 - xlam6*dx6)*rm4[igrp];
  clatoms_vy[ind4[igrp]] -=  (-xlam3*dy3 - xlam5*dy5 - xlam6*dy6)*rm4[igrp];
  clatoms_vz[ind4[igrp]] -=  (-xlam3*dz3 - xlam5*dz5 - xlam6*dz6)*rm4[igrp];

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

/* free memory locally assigned */
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


    cfree(&(ind1[1]));
    cfree(&(ind2[1]));
    cfree(&(ind3[1]));
    cfree(&(ind4[1]));

/*=======================================================================*/
} /* end routine */
/*=======================================================================*/

