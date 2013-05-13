/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
/*                                                                          */
/*                         PI_MD:                                           */
/*             The future of simulation technology                          */
/*             ------------------------------------                         */
/*                     Module: recip_pme.c                                  */
/*                                                                          */
/* This subprogram set compute the reciprocal space part of the ewald sum   */
/* using the Particle Mesh Ewald approximation                              */
/*                                                                          */
/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
/*               Includes:                                                  */

#include "standard_include.h"
#include "../typ_defs/typedefs_gen.h"
#include "../typ_defs/typedefs_class.h"
#include "../proto_defs/proto_math.h"
#include "../proto_defs/proto_recip3d_entry.h"
#include "../proto_defs/proto_communicate_wrappers.h"
#include "../proto_defs/proto_recip3d_local.h"

#define DEBUG_PME_OFF
#define FFT
#define COMM_GRID

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
/* PME controller */
/*==========================================================================*/

void ewald3d_recip_pme(CLATOMS_INFO *clatoms_info,
                       CLATOMS_POS *clatoms_pos,
                       CELL *cell,PTENS *ptens,
                       double alp_ewd,int nktot,
                       int kastr[],int kbstr[],int kcstr[],
                       EWD_SCR *ewd_scr,double *vrecip,
                       double wght_now, int iver_get,PART_MESH *part_mesh,
                       int irespa,CLASS_COMM_FORC_PKG *class_comm_forc_pkg,
                       int iget_pv_real_inter,PARA_FFT_PKG3D *para_fft_pkg3d,
                       FOR_SCR *for_scr, 
                       double *clus_corr_r, double *dclus_corr_r)

/*==========================================================================*/
     {/*begin routine*/ 
/*==========================================================================*/
/*    Local Variables                                                       */

   int pme_para_opt = part_mesh->pme_para_opt;

/*==========================================================================*/

  switch(pme_para_opt){

    case 0: ewald3d_recip_pme_none(clatoms_info,clatoms_pos,cell,ptens,
                          alp_ewd,nktot,kastr,kbstr,kcstr,ewd_scr,vrecip,
                          wght_now,iver_get,part_mesh,irespa,
                          class_comm_forc_pkg,iget_pv_real_inter,
                          para_fft_pkg3d,clus_corr_r,dclus_corr_r);
                          break;

    case 1: ewald3d_recip_pme_hybr(clatoms_info,clatoms_pos,cell,ptens,
                          alp_ewd,nktot,kastr,kbstr,kcstr,ewd_scr,vrecip,
                          wght_now,iver_get,part_mesh,irespa,
                          class_comm_forc_pkg,iget_pv_real_inter,
                          para_fft_pkg3d,clus_corr_r,dclus_corr_r);
                          break;

    case 2: ewald3d_recip_pme_full_g(clatoms_info,clatoms_pos,cell,ptens,
                          alp_ewd,nktot,kastr,kbstr,kcstr,ewd_scr,vrecip,
                          wght_now,iver_get,part_mesh,irespa,
                          class_comm_forc_pkg,iget_pv_real_inter,
                          para_fft_pkg3d,for_scr,clus_corr_r,dclus_corr_r);
                          break;

  }/*end switch*/

/*--------------------------------------------------------------------------*/
     }/*end routine*/
/*==========================================================================*/




/*==========================================================================*/
/* No-parallel */
/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void ewald3d_recip_pme_none(CLATOMS_INFO *clatoms_info,
                          CLATOMS_POS *clatoms_pos,
                          CELL *cell,PTENS *ptens,
                          double alp_ewd,int nktot,
                          int kastr[],int kbstr[],int kcstr[],
                          EWD_SCR *ewd_scr,double *vrecip,
                          double wght_now, int iver_get,PART_MESH *part_mesh,
                          int irespa,CLASS_COMM_FORC_PKG *class_comm_forc_pkg,
                          int iget_pv_real_inter,
                          PARA_FFT_PKG3D *para_fft_pkg3d,
                          double *clus_corr_r, double *dclus_corr_r)
               
/*==========================================================================*/
     {/*begin routine*/ 
/*==========================================================================*/
/*    Local Variables                                                       */

#include "../typ_defs/typ_mask.h"
  int i,j,k,iatm,ktemp,n,j2,k1,j1,inow;
  int ia,ib,ic;
  int ja,jb,jc;
  int iii;
  int nk1,nk2,nk3,init,ireal;   
  double mn_a_tmp,mn_b_tmp,mn_c_tmp,prep,vnow;
  double atemp,btemp,ctemp;
  double qgrid_dx,qgrid_dy,qgrid_dz;
  double falp2,vol,tpi,pivol,g2,preg,rvol;
  double aka,akb,akc,q_sum1;
  double xk,yk,zk,smag;
  double qgrid_real,qgrid_imag,vnow_tmp;
  int ngrid_a,ngrid_b,ngrid_c,nrecip_grid;
  int ngrid_ab;
  int ngrid_bc,nedge,n_interp;
  double grid_a,grid_b,grid_c;
  int iatm1,iend,nnow,itemp;
  int nlen_pme;
  double pvten_temp[10];
/* Define local pointers                                                  */

  int nchrg               = clatoms_info->nchrg;
  int natm_tot            = clatoms_info->natm_tot;
  double *clatoms_x       = clatoms_pos->x;
  double *clatoms_y       = clatoms_pos->y;
  double *clatoms_z       = clatoms_pos->z;
  double *clatoms_q       = clatoms_info->q;
  double *clatoms_fx      = clatoms_pos->fx;
  double *clatoms_fy      = clatoms_pos->fy;
  double *clatoms_fz      = clatoms_pos->fz;
  double *clatoms_fxt     = clatoms_pos->fxt;
  double *clatoms_fyt     = clatoms_pos->fyt;
  double *clatoms_fzt     = clatoms_pos->fzt;
  double *cell_hmati      = cell->hmati;
  double *cell_hmat      = cell->hmat;
  double *pvten     = ptens->pvten;
  double *pvten_tot = ptens->pvten_tot;
  double *ewd_scr_x       = ewd_scr->x;
  double *ewd_scr_y       = ewd_scr->y;
  double *ewd_scr_z       = ewd_scr->z;
  double *ewd_scr_q       = ewd_scr->q;
  double *ewd_scr_fx      = ewd_scr->fx;
  double *ewd_scr_fy      = ewd_scr->fy;
  double *ewd_scr_fz      = ewd_scr->fz; 
  double *pvten_tmp = ptens->pvten_tmp;
  int *clatoms_ichrg      = clatoms_info->ichrg;

  int *isphere_real;
  int *isphere_imag;
  int *isphere_real_c;
  int *isphere_imag_c;
  double *bweight_tot;
  double *qgrid          = part_mesh->qgrid;
  double *qgrid_tmp      = part_mesh->qgrid_scr;
  double *aj             = part_mesh->aj;
  double *rn             = part_mesh->rn;
  double *rn1            = part_mesh->rn1;
  int *iatemp            = part_mesh->iatemp;
  int *ibtemp            = part_mesh->ibtemp;
  int *ictemp            = part_mesh->ictemp;
  double *frac_a         = part_mesh->frac_a;
  double *frac_b         = part_mesh->frac_b;
  double *frac_c         = part_mesh->frac_c;
  double *qgrid_tmp_real = part_mesh->qgrid_tmp_real;
  double *qgrid_tmp_imag = part_mesh->qgrid_tmp_imag;
  double **qgrid_now     = part_mesh->qgrid_now;
  double **ua            = part_mesh->ua;
  double **ub            = part_mesh->ub;
  double **uc            = part_mesh->uc;
  double **mn_a          = part_mesh->mn_a;
  double **mn_b          = part_mesh->mn_b;
  double **mn_c          = part_mesh->mn_c;
  double **dmn_a         = part_mesh->dmn_a;
  double **dmn_b         = part_mesh->dmn_b;
  double **dmn_c         = part_mesh->dmn_c;
  int **igrid_a          = part_mesh->igrid_a;
  int **igrid_b          = part_mesh->igrid_b;
  int **igrid_c          = part_mesh->igrid_c;
  int **igrid_now        = part_mesh->igrid_now;
  int np_forc            = class_comm_forc_pkg->num_proc;
  int myid_forc          = class_comm_forc_pkg->myid;
  MPI_Comm comm_forc = class_comm_forc_pkg->comm;
  double dnp_forc_i;
  int iperd               = cell->iperd; 

  dnp_forc_i = 1.0/(double) (np_forc);

/*==========================================================================*/
/* I) Construct some useful constants                                       */

   nlen_pme       = part_mesh->nlen_pme;
   if(irespa==0){
    ngrid_a        = part_mesh->ngrid_a;
    ngrid_b        = part_mesh->ngrid_b;
    ngrid_c        = part_mesh->ngrid_c;
    n_interp       = part_mesh->n_interp;
    bweight_tot    = part_mesh->bweight_tot;
   }else{
    ngrid_a        = part_mesh->ngrid_a_res;
    ngrid_b        = part_mesh->ngrid_b_res;
    ngrid_c        = part_mesh->ngrid_c_res;
    n_interp       = part_mesh->n_interp_res;
    bweight_tot    = part_mesh->bweight_tot_res;
   }/*endif*/
   nrecip_grid    = 2*ngrid_a*ngrid_b*ngrid_c;
   grid_a         = (double) ngrid_a;
   grid_b         = (double) ngrid_b;
   grid_c         = (double) ngrid_c;
   ngrid_bc       = ngrid_b*ngrid_c;
   ngrid_ab       = ngrid_a*ngrid_b;
   for(j=1;j<=n_interp;j++){
     aj[j] = (double) (j-1);
     rn[j] = (double) (j);
     if(j > 1){rn1[j] = 1.0/((double)(j-1));}
   }/*endfor*/
     rn1[1] = 0.0;

/*==========================================================================*/
/* II) Find charged atoms, get their scaled imaged particle coordinates     */


   for(iatm = 1;iatm <= nchrg;++iatm){
      ktemp = clatoms_ichrg[iatm];
      ewd_scr_x[iatm] = clatoms_x[ktemp];
      ewd_scr_y[iatm] = clatoms_y[ktemp];
      ewd_scr_z[iatm] = clatoms_z[ktemp];
      ewd_scr_q[iatm] = clatoms_q[ktemp];
   }/*endfor*/
   for(iatm = 1;iatm <= nchrg;++iatm){
     atemp = ewd_scr_x[iatm]*cell_hmati[1]
           + ewd_scr_y[iatm]*cell_hmati[4]
           + ewd_scr_z[iatm]*cell_hmati[7];
     btemp = ewd_scr_x[iatm]*cell_hmati[2]
           + ewd_scr_y[iatm]*cell_hmati[5]
           + ewd_scr_z[iatm]*cell_hmati[8];
     ctemp = ewd_scr_x[iatm]*cell_hmati[3]
           + ewd_scr_y[iatm]*cell_hmati[6]
           + ewd_scr_z[iatm]*cell_hmati[9];
     atemp = atemp - NINT((atemp-0.5));
     btemp = btemp - NINT((btemp-0.5));
     ctemp = ctemp - NINT((ctemp-0.5));
     ewd_scr_x[iatm] = atemp*grid_a;
     ewd_scr_y[iatm] = btemp*grid_b;
     ewd_scr_z[iatm] = ctemp*grid_c;
   }/*endfor*/
/*==========================================================================*/
/* III) Calculate the Cardinal B spline interpolation functions of the      */
/*     charge weighted density on the real space grid                       */

   for(i=1;i<=nrecip_grid;i++){
     qgrid[i]=0.0;
   }/*endfor*/

   for(iatm=1;iatm<=nchrg;iatm+=nlen_pme){
     iatm1 = iatm-1;
     iend = MIN(nchrg,iatm1+nlen_pme);
     nnow = iend-iatm1;
/*--------------------------------------------------------------------------*/ 
/* A) Using current fraction, find the grid points on which M_n is non-zero */
     for(i=1;i<=nnow;i++){
      itemp     = i+iatm1;
      iatemp[i] = (int) (ewd_scr_x[itemp]);
      ibtemp[i] = (int) (ewd_scr_y[itemp]);
      ictemp[i] = (int) (ewd_scr_z[itemp]);
      frac_a[i] = ewd_scr_x[itemp] - (double) (iatemp[i]);
      frac_b[i] = ewd_scr_y[itemp] - (double) (ibtemp[i]);
      frac_c[i] = ewd_scr_z[itemp] - (double) (ictemp[i]);
     }/*endfor*/
     for(j=1;j<=n_interp;j++){
      for(i=1;i<=nnow;i++){
       ua[j][i]    = frac_a[i] + aj[j];
       ub[j][i]    = frac_b[i] + aj[j];
       uc[j][i]    = frac_c[i] + aj[j];
       j2       = j-2;
       ia       = iatemp[i] - j2;
       ib       = ibtemp[i] - j2;
       ic       = ictemp[i] - j2;
       ia       = (ia>0 ? ia:ngrid_a+ia);
       ib       = (ib>0 ? ib:ngrid_b+ib);
       ic       = (ic>0 ? ic:ngrid_c+ic);

       igrid_a[j][i] = 2*ia - 1;
       igrid_b[j][i] = 2*(ib - 1)*ngrid_a;
       igrid_c[j][i] = 2*(ic - 1)*ngrid_ab;

      }/*endfor*/
     }/*endfor*/
/*--------------------------------------------------------------------------*/ 
/* B) Initialize M2 and get the Mn's using the recursion relation           */ 
/*    Note: M_n is defined on 0 to n. Since frac between 0 and 1, the       */
/*          calculation is performed in an order that takes advantage of it */
     for(i=1;i<=nnow;i++){
      mn_a[1][i] = 1.0 - fabs(ua[1][i]-1.0);
      mn_b[1][i] = 1.0 - fabs(ub[1][i]-1.0);
      mn_c[1][i] = 1.0 - fabs(uc[1][i]-1.0);
      mn_a[2][i] = 1.0 - fabs(ua[2][i]-1.0);
      mn_b[2][i] = 1.0 - fabs(ub[2][i]-1.0);
      mn_c[2][i] = 1.0 - fabs(uc[2][i]-1.0);
     }/*endfor*/
     for(j=3;j<=n_interp;j++){
      for(i=1;i<=nnow;i++){
       mn_a[j][i]   = 0.0;
       mn_b[j][i]   = 0.0;
       mn_c[j][i]   = 0.0;
      }/*endfor*/
     }/*endfor*/

     for(n=3;n<=n_interp;n++){
       for(j=n;j>=2;j--){
        j1 = j-1;
        for(i=1;i<=nnow;i++){
         mn_a_tmp  = (ua[j][i]*mn_a[j][i]+(rn[n]-ua[j][i])*mn_a[j1][i])*rn1[n];
         mn_b_tmp  = (ub[j][i]*mn_b[j][i]+(rn[n]-ub[j][i])*mn_b[j1][i])*rn1[n];
         mn_c_tmp  = (uc[j][i]*mn_c[j][i]+(rn[n]-uc[j][i])*mn_c[j1][i])*rn1[n];
         mn_a[j][i] = mn_a_tmp;
         mn_b[j][i] = mn_b_tmp;
         mn_c[j][i] = mn_c_tmp;
        }/*end for: i*/
       }/*end for: j*/
       for(i=1;i<=nnow;i++){
        mn_a[1][i] = ua[1][i]*mn_a[1][i]*rn1[n];
        mn_b[1][i] = ub[1][i]*mn_b[1][i]*rn1[n];
        mn_c[1][i] = uc[1][i]*mn_c[1][i]*rn1[n];
      }/*endfor */
     }/*end for: k*/

/*--------------------------------------------------------------------------*/ 
/* C) Stuff Mn's on the grid:Sum is ordered to remove recursive dependencies*/ 
/*                           in qgrid. igrid_now is unique for fixed i      */ 
     for(jc=1;jc<=n_interp;jc++){
     for(jb=1;jb<=n_interp;jb++){
     for(ja=1;ja<=n_interp;ja++){
     for(i=1;i<=nnow;i++){
       itemp     = i+iatm1;
       igrid_now[ja][i] = igrid_a[ja][i]+igrid_b[jb][i]+igrid_c[jc][i];
       qgrid_now[ja][i] = mn_a[ja][i]*mn_b[jb][i]*mn_c[jc][i]*ewd_scr_q[itemp];
     }}/*endfor*/
     for(i=1;i<=nnow;i++){
#ifndef NO_PRAGMA
#pragma IVDEP
#endif
     for(ja=1;ja<=n_interp;ja++){
       qgrid[igrid_now[ja][i]] += qgrid_now[ja][i];
     }}}}/*end for: j*/
   }/*end for: iatm*/

/*==========================================================================*/
/* IV) Fourier Transform qgrid                                              */

   para_fft_gen3d_bck_to_g(qgrid,qgrid_tmp,para_fft_pkg3d); 
   sngl_upack_coef(qgrid_tmp_real,qgrid_tmp_imag,qgrid,para_fft_pkg3d);

/*==========================================================================*/
/* VI) Compute the potential energy on the spherically cutoff grid.         */

   falp2 = 4.0*alp_ewd*alp_ewd;
   vol = getdeth(cell->hmat);
   tpi = 2.0*M_PI;
   pivol = vol/(4.0*M_PI);
   rvol  = 1.0/vol;  

   for(i=1;i<=9;i++){pvten_tmp[i]=0.0;}

   vnow = 0.0;
   for(i=1;i <= nktot; ++i) {
     aka = tpi*( (double) kastr[i] );
     akb = tpi*( (double) kbstr[i] );
     akc = tpi*( (double) kcstr[i] );
     xk = (aka*cell_hmati[1] + akb*cell_hmati[2] + akc*cell_hmati[3]);
     yk = (aka*cell_hmati[4] + akb*cell_hmati[5] + akc*cell_hmati[6]);
     zk = (aka*cell_hmati[7] + akb*cell_hmati[8] + akc*cell_hmati[9]);
     g2 = xk*xk + yk*yk + zk*zk;
     preg = exp((-g2/falp2))/(g2*pivol);
     prep = -2.0*preg*((g2/falp2)+1.0)/g2;
     if(iperd != 3){
       preg += clus_corr_r[i]*rvol;
       prep += dclus_corr_r[i]*rvol;
     }/*endif*/
     smag       = (qgrid_tmp_real[i]*qgrid_tmp_real[i]
                  +qgrid_tmp_imag[i]*qgrid_tmp_imag[i])*bweight_tot[i];
     vnow += (smag*preg);
     prep    *= smag;
     pvten_tmp[1] += prep*xk*xk;
     pvten_tmp[5] += prep*yk*yk;
     pvten_tmp[9] += prep*zk*zk;
     pvten_tmp[2] += prep*xk*yk;
     pvten_tmp[3] += prep*xk*zk;
     pvten_tmp[6] += prep*yk*zk;
     qgrid_tmp_real[i] *= (preg*bweight_tot[i]);
     qgrid_tmp_imag[i] *= (preg*bweight_tot[i]);
   }/*endfor*/
   qgrid_tmp_real[(nktot+1)] = 0.0;
   qgrid_tmp_imag[(nktot+1)] = 0.0;

   *vrecip += vnow*dnp_forc_i;
   if(iperd!=3) {
     q_sum1 = dsum1(natm_tot,clatoms_q,1);
     *vrecip += (0.5*q_sum1*q_sum1*clus_corr_r[(nktot+1)]*rvol*dnp_forc_i);
   }/*endif*/

   pvten_tmp[4]  = pvten_tmp[2];
   pvten_tmp[7]  = pvten_tmp[3];
   pvten_tmp[8]  = pvten_tmp[6];
   pvten_tmp[1] += vnow;
   pvten_tmp[5] += vnow;
   pvten_tmp[9] += vnow;


   for(i=1;i<=9;i++){
    pvten[i]     += (pvten_tmp[i]*wght_now*dnp_forc_i);
   }/*endfor*/

   if(iget_pv_real_inter==1&&irespa==0){
    for(i=1;i<=9;i++){   
     pvten_tot[i] += (pvten_tmp[i]*dnp_forc_i);
    }/*endfor*/
   }/*endif*/

/*==========================================================================*/
/* VII) Fourier Transform qgrid                                             */

   pme_sngl_pack_coef(qgrid_tmp_real,qgrid_tmp_imag,qgrid,para_fft_pkg3d);
   para_fft_gen3d_fwd_to_r(qgrid,qgrid_tmp,para_fft_pkg3d); 

/*==========================================================================*/
/* VIII) Calculate the force                                                */

   for(iatm = 1;iatm <= nchrg;++iatm){
      ewd_scr_fx[iatm] = 0.0;
      ewd_scr_fy[iatm] = 0.0;
      ewd_scr_fz[iatm] = 0.0;
   }/*endfor*/
   for(iatm=1;iatm<=nchrg;iatm+=nlen_pme){
     iatm1 = iatm-1;
     iend = MIN(nchrg,iatm1+nlen_pme);
     nnow = iend-iatm1;
/*--------------------------------------------------------------------------*/ 
/* A) Using current fraction, find the grid points on which M_n is non-zero */
     for(i=1;i<=nnow;i++){
      itemp     = i+iatm1;
      iatemp[i] = (int) (ewd_scr_x[itemp]);
      ibtemp[i] = (int) (ewd_scr_y[itemp]);
      ictemp[i] = (int) (ewd_scr_z[itemp]);
      frac_a[i] = ewd_scr_x[itemp] - (double) (iatemp[i]);
      frac_b[i] = ewd_scr_y[itemp] - (double) (ibtemp[i]);
      frac_c[i] = ewd_scr_z[itemp] - (double) (ictemp[i]);
     }/*endfor*/
     for(j=1;j<=n_interp;j++){
      for(i=1;i<=nnow;i++){
       ua[j][i]    = frac_a[i] + aj[j];
       ub[j][i]    = frac_b[i] + aj[j];
       uc[j][i]    = frac_c[i] + aj[j];
       j2       = j-2;
       ia       = iatemp[i] - j2;
       ib       = ibtemp[i] - j2;
       ic       = ictemp[i] - j2;
       ia       = (ia>0 ? ia:ngrid_a+ia);
       ib       = (ib>0 ? ib:ngrid_b+ib);
       ic       = (ic>0 ? ic:ngrid_c+ic);

       igrid_a[j][i] = 2*ia - 1;
       igrid_b[j][i] = 2*(ib - 1)*ngrid_a;
       igrid_c[j][i] = 2*(ic - 1)*ngrid_ab;

      }/*endfor*/
     }/*endfor*/
/*--------------------------------------------------------------------------*/ 
/* B) Initialize M2 and get the Mn's using the recursion relation           */ 
/*    Note: M_n is defined on 0 to n. Since frac between 0 and 1, the       */
/*          calculation is performed in an order that takes advantage of it */
     for(i=1;i<=nnow;i++){
      mn_a[1][i] = 1.0 - fabs(ua[1][i]-1.0);
      mn_b[1][i] = 1.0 - fabs(ub[1][i]-1.0);
      mn_c[1][i] = 1.0 - fabs(uc[1][i]-1.0);
      mn_a[2][i] = 1.0 - fabs(ua[2][i]-1.0);
      mn_b[2][i] = 1.0 - fabs(ub[2][i]-1.0);
      mn_c[2][i] = 1.0 - fabs(uc[2][i]-1.0);
     }/*endfor*/
     for(j=3;j<=n_interp;j++){
      for(i=1;i<=nnow;i++){
       mn_a[j][i]   = 0.0;
       mn_b[j][i]   = 0.0;
       mn_c[j][i]   = 0.0;
      }/*endfor*/
     }/*endfor*/

     for(n=3;n<=n_interp;n++){
       for(j=n;j>=2;j--){
        j1 = j-1;
        for(i=1;i<=nnow;i++){
         mn_a_tmp  = (ua[j][i]*mn_a[j][i]+(rn[n]-ua[j][i])*mn_a[j1][i])*rn1[n];
         mn_b_tmp  = (ub[j][i]*mn_b[j][i]+(rn[n]-ub[j][i])*mn_b[j1][i])*rn1[n];
         mn_c_tmp  = (uc[j][i]*mn_c[j][i]+(rn[n]-uc[j][i])*mn_c[j1][i])*rn1[n];
         mn_a[j][i] = mn_a_tmp;
         mn_b[j][i] = mn_b_tmp;
         mn_c[j][i] = mn_c_tmp;
        }/*end for: i*/
       }/*end for: j*/
       for(i=1;i<=nnow;i++){
        mn_a[1][i] = ua[1][i]*mn_a[1][i]*rn1[n];
        mn_b[1][i] = ub[1][i]*mn_b[1][i]*rn1[n];
        mn_c[1][i] = uc[1][i]*mn_c[1][i]*rn1[n];
       }/*endfor */
       if(n==(n_interp-1)){
        for(i=1;i<=nnow;i++){
         dmn_a[1][i] = mn_a[1][i];
         dmn_b[1][i] = mn_b[1][i];
         dmn_c[1][i] = mn_c[1][i];
        }/*endfor*/
        for(j=2;j<=n_interp;j++){    
         j1 = j-1;
         for(i=1;i<=nnow;i++){
          dmn_a[j][i] = mn_a[j][i] - mn_a[j1][i];
          dmn_b[j][i] = mn_b[j][i] - mn_b[j1][i];
          dmn_c[j][i] = mn_c[j][i] - mn_c[j1][i];
         }/*endfor*/
        }/*endfor*/
       }/*endif*/
     }/*end for: n*/
/*--------------------------------------------------------------------------*/ 
/* C) Calculate the force                                                   */ 
     for(jc=1;jc<=n_interp;jc++){
     for(jb=1;jb<=n_interp;jb++){
     for(ja=1;ja<=n_interp;ja++){
     for(i=1;i<=nnow;i++){
       igrid_now[ja][i] = igrid_a[ja][i]+igrid_b[jb][i]+igrid_c[jc][i];
     }}/*endfor*/
     for(ja=1;ja<=n_interp;ja++){
     for(i=1;i<=nnow;i++){
       qgrid_now[ja][i] = qgrid[igrid_now[ja][i]];
     }}/*endfor*/
     for(ja=1;ja<=n_interp;ja++){
     for(i=1;i<=nnow;i++){
        itemp = i+iatm1;
        atemp = dmn_a[ja][i]* mn_b[jb][i]* mn_c[jc][i]*grid_a;
        btemp =  mn_a[ja][i]*dmn_b[jb][i]* mn_c[jc][i]*grid_b;
        ctemp =  mn_a[ja][i]* mn_b[jb][i]*dmn_c[jc][i]*grid_c;
        qgrid_dx = atemp*cell_hmati[1]+btemp*cell_hmati[2]+ctemp*cell_hmati[3];
        qgrid_dy = atemp*cell_hmati[4]+btemp*cell_hmati[5]+ctemp*cell_hmati[6];
        qgrid_dz = atemp*cell_hmati[7]+btemp*cell_hmati[8]+ctemp*cell_hmati[9];
        qgrid_now[ja][i] *= ewd_scr_q[itemp];
        ewd_scr_fx[itemp] -= (qgrid_dx*qgrid_now[ja][i]);
        ewd_scr_fy[itemp] -= (qgrid_dy*qgrid_now[ja][i]);
        ewd_scr_fz[itemp] -= (qgrid_dz*qgrid_now[ja][i]); 
     }}}}/*end for: j*/
   }/*end for: iatm*/

/*==========================================================================*/
/* IX) Scatter the force back and get virial estimator if PIMD on           */

    if(iver_get==1){
#ifndef NO_PRAGMA
#pragma IVDEP
#endif
     for(iatm = 1;iatm <= nchrg;++iatm){
      ktemp = clatoms_ichrg[iatm];
      clatoms_fxt[ktemp] += ewd_scr_fx[iatm];
      clatoms_fyt[ktemp] += ewd_scr_fy[iatm];
      clatoms_fzt[ktemp] += ewd_scr_fz[iatm];
     }/*endfor*/
    }/*endif*/
    for(iatm = 1;iatm <= nchrg;++iatm){
      ewd_scr_fx[iatm] *= wght_now*dnp_forc_i;
      ewd_scr_fy[iatm] *= wght_now*dnp_forc_i;
      ewd_scr_fz[iatm] *= wght_now*dnp_forc_i;
    }/*endfor*/
#ifndef NO_PRAGMA
#pragma IVDEP
#endif
    for(iatm = 1;iatm <= nchrg;++iatm){
      ktemp = clatoms_ichrg[iatm];
      clatoms_fx[ktemp] += ewd_scr_fx[iatm];
      clatoms_fy[ktemp] += ewd_scr_fy[iatm];
      clatoms_fz[ktemp] += ewd_scr_fz[iatm];
    }/*endfor*/

/*-----------------------------------------------------------------------*/ 
}/*end routine*/
/*=======================================================================*/





/*==========================================================================*/
/* No-parallel */
/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void ewald3d_recip_pme_hybr(CLATOMS_INFO *clatoms_info,
                          CLATOMS_POS *clatoms_pos,
                          CELL *cell,PTENS *ptens,
                          double alp_ewd,int nktot,
                          int kastr[],int kbstr[],int kcstr[],
                          EWD_SCR *ewd_scr,double *vrecip,
                          double wght_now, int iver_get,PART_MESH *part_mesh,
                          int irespa,CLASS_COMM_FORC_PKG *class_comm_forc_pkg,
                          int iget_pv_real_inter,
                          PARA_FFT_PKG3D *para_fft_pkg3d,
                          double *clus_corr_r, double *dclus_corr_r)

/*==========================================================================*/
     {/*begin routine*/ 
/*==========================================================================*/
/*    Local Variables                                                       */

#include "../typ_defs/typ_mask.h"
  int i,j,k,iatm,ktemp,n,j2,k1,j1,inow,iproc;
  int ia,ib,ic;
  int ja,jb,jc;
  int iii;
  int nk1,nk2,nk3,init,ireal;   
  double mn_a_tmp,mn_b_tmp,mn_c_tmp,prep,vnow;
  double atemp,btemp,ctemp;
  double qgrid_dx,qgrid_dy,qgrid_dz;
  double falp2,vol,tpi,pivol,g2,preg,rvol;
  double aka,akb,akc,q_sum1;
  double xk,yk,zk,smag;
  double qgrid_real,qgrid_imag,vnow_tmp;
  int ngrid_a,ngrid_b,ngrid_c,nrecip_grid;
  int ngrid_ab;
  int ngrid_bc,nedge,n_interp;
  double grid_a,grid_b,grid_c;
  int iatm1,iend,nnow,itemp;
  int nlen_pme;
  double pvten_temp[10];
  int idiv,irem,natm_proc,iatm_str,iatm_end;

/* Define local pointers                                                  */

  int nchrg               = clatoms_info->nchrg;
  int natm_tot            = clatoms_info->natm_tot;
  double *clatoms_x       = clatoms_pos->x;
  double *clatoms_y       = clatoms_pos->y;
  double *clatoms_z       = clatoms_pos->z;
  double *clatoms_q       = clatoms_info->q;
  double *clatoms_fx      = clatoms_pos->fx;
  double *clatoms_fy      = clatoms_pos->fy;
  double *clatoms_fz      = clatoms_pos->fz;
  double *clatoms_fxt     = clatoms_pos->fxt;
  double *clatoms_fyt     = clatoms_pos->fyt;
  double *clatoms_fzt     = clatoms_pos->fzt;
  double *cell_hmati      = cell->hmati;
  double *cell_hmat      = cell->hmat;
  double *pvten     = ptens->pvten;
  double *pvten_tot = ptens->pvten_tot;
  double *ewd_scr_x       = ewd_scr->x;
  double *ewd_scr_y       = ewd_scr->y;
  double *ewd_scr_z       = ewd_scr->z;
  double *ewd_scr_q       = ewd_scr->q;
  double *ewd_scr_fx      = ewd_scr->fx;
  double *ewd_scr_fy      = ewd_scr->fy;
  double *ewd_scr_fz      = ewd_scr->fz; 
  double *pvten_tmp = ptens->pvten_tmp;
  int *clatoms_ichrg      = clatoms_info->ichrg;

  int *isphere_real;
  int *isphere_imag;
  int *isphere_real_c;
  int *isphere_imag_c;
  double *bweight_tot;
  double *qgrid          = part_mesh->qgrid;
  double *qgrid_tmp      = part_mesh->qgrid_scr;
  double *aj             = part_mesh->aj;
  double *rn             = part_mesh->rn;
  double *rn1            = part_mesh->rn1;
  int *iatemp            = part_mesh->iatemp;
  int *ibtemp            = part_mesh->ibtemp;
  int *ictemp            = part_mesh->ictemp;
  double *frac_a         = part_mesh->frac_a;
  double *frac_b         = part_mesh->frac_b;
  double *frac_c         = part_mesh->frac_c;
  double *qgrid_tmp_real = part_mesh->qgrid_tmp_real;
  double *qgrid_tmp_imag = part_mesh->qgrid_tmp_imag;
  double **qgrid_now     = part_mesh->qgrid_now;
  double **ua            = part_mesh->ua;
  double **ub            = part_mesh->ub;
  double **uc            = part_mesh->uc;
  double **mn_a          = part_mesh->mn_a;
  double **mn_b          = part_mesh->mn_b;
  double **mn_c          = part_mesh->mn_c;
  double **dmn_a         = part_mesh->dmn_a;
  double **dmn_b         = part_mesh->dmn_b;
  double **dmn_c         = part_mesh->dmn_c;
  int **igrid_a          = part_mesh->igrid_a;
  int **igrid_b          = part_mesh->igrid_b;
  int **igrid_c          = part_mesh->igrid_c;
  int **igrid_now        = part_mesh->igrid_now;

  int np_forc            = class_comm_forc_pkg->num_proc;
  int myid_forc          = class_comm_forc_pkg->myid;
  MPI_Comm comm_forc     = class_comm_forc_pkg->comm;
  int *recv_counts_coef  = para_fft_pkg3d->recv_counts_coef;
  int *recv_dspls_coef   = para_fft_pkg3d->recv_dspls_coef;
  int ncoef,ncoef_use,ncoef_proc,icoef_off,icoef_strt;
  int iperd              = cell->iperd; 

/*==========================================================================*/
/* I) Construct some useful constants                                       */

   nlen_pme       = part_mesh->nlen_pme;
   if(irespa==0){
    ngrid_a        = part_mesh->ngrid_a;
    ngrid_b        = part_mesh->ngrid_b;
    ngrid_c        = part_mesh->ngrid_c;
    n_interp       = part_mesh->n_interp;
    bweight_tot    = part_mesh->bweight_tot;
   }else{
    ngrid_a        = part_mesh->ngrid_a_res;
    ngrid_b        = part_mesh->ngrid_b_res;
    ngrid_c        = part_mesh->ngrid_c_res;
    n_interp       = part_mesh->n_interp_res;
    bweight_tot    = part_mesh->bweight_tot_res;
   }/*endif*/
   nrecip_grid    = 2*ngrid_a*ngrid_b*ngrid_c;
   grid_a         = (double) ngrid_a;
   grid_b         = (double) ngrid_b;
   grid_c         = (double) ngrid_c;
   ngrid_bc       = ngrid_b*ngrid_c;
   ngrid_ab       = ngrid_a*ngrid_b;
   for(j=1;j<=n_interp;j++){
     aj[j] = (double) (j-1);
     rn[j] = (double) (j);
     if(j > 1){rn1[j] = 1.0/((double)(j-1));}
   }/*endfor*/
     rn1[1] = 0.0;

/*==========================================================================*/
/* II) Find charged atoms, get their scaled imaged particle coordinates     */

   idiv        = nchrg / np_forc;
   irem        = nchrg % np_forc;
   natm_proc   = (myid_forc <  irem ? idiv+1 : idiv);
   iatm_str    = (myid_forc <= irem ? myid_forc*(idiv+1)+1 : 
                                      irem*(idiv+1)+1+(myid_forc-irem)*idiv);
   iatm_end    = natm_proc + iatm_str - 1;

#ifdef DEBUG_PME
  for(iproc=1;iproc<=np_forc;iproc++){
    if(iproc==myid_forc+1){
      printf("%d atm_strt %d atm_end %d\n",iproc,iatm_str,iatm_end);
    }/*endif*/
    Barrier(comm_forc);
  }/*endfor*/
  if(myid_forc==0){scanf("%d",&iii);} Barrier(comm_forc);
#endif

   for(iatm = iatm_str;iatm <= iatm_end;++iatm){
      ktemp = clatoms_ichrg[iatm];
      ewd_scr_x[iatm] = clatoms_x[ktemp];
      ewd_scr_y[iatm] = clatoms_y[ktemp];
      ewd_scr_z[iatm] = clatoms_z[ktemp];
      ewd_scr_q[iatm] = clatoms_q[ktemp];
   }/*endfor*/
   for(iatm = iatm_str;iatm <= iatm_end;++iatm){
     atemp = ewd_scr_x[iatm]*cell_hmati[1]
           + ewd_scr_y[iatm]*cell_hmati[4]
           + ewd_scr_z[iatm]*cell_hmati[7];
     btemp = ewd_scr_x[iatm]*cell_hmati[2]
           + ewd_scr_y[iatm]*cell_hmati[5]
           + ewd_scr_z[iatm]*cell_hmati[8];
     ctemp = ewd_scr_x[iatm]*cell_hmati[3]
           + ewd_scr_y[iatm]*cell_hmati[6]
           + ewd_scr_z[iatm]*cell_hmati[9];
     atemp = atemp - NINT((atemp-0.5));
     btemp = btemp - NINT((btemp-0.5));
     ctemp = ctemp - NINT((ctemp-0.5));
     ewd_scr_x[iatm] = atemp*grid_a;
     ewd_scr_y[iatm] = btemp*grid_b;
     ewd_scr_z[iatm] = ctemp*grid_c;
   }/*endfor*/

/*==========================================================================*/
/* III) Calculate the Cardinal B spline interpolation functions of the      */
/*     charge weighted density on the real space grid                       */

   for(i=1;i<=nrecip_grid;i++){
     qgrid[i]=0.0;
   }/*endfor*/

   for(iatm=iatm_str;iatm<=iatm_end;iatm+=nlen_pme){
     iatm1 = iatm-1;
     iend = MIN(iatm_end,iatm1+nlen_pme);
     nnow = iend-iatm1;
/*--------------------------------------------------------------------------*/ 
/* A) Using current fraction, find the grid points on which M_n is non-zero */
     for(i=1;i<=nnow;i++){
      itemp     = i+iatm1;
      iatemp[i] = (int) (ewd_scr_x[itemp]);
      ibtemp[i] = (int) (ewd_scr_y[itemp]);
      ictemp[i] = (int) (ewd_scr_z[itemp]);
      frac_a[i] = ewd_scr_x[itemp] - (double) (iatemp[i]);
      frac_b[i] = ewd_scr_y[itemp] - (double) (ibtemp[i]);
      frac_c[i] = ewd_scr_z[itemp] - (double) (ictemp[i]);
     }/*endfor*/
     for(j=1;j<=n_interp;j++){
      for(i=1;i<=nnow;i++){
       ua[j][i]    = frac_a[i] + aj[j];
       ub[j][i]    = frac_b[i] + aj[j];
       uc[j][i]    = frac_c[i] + aj[j];
       j2       = j-2;
       ia       = iatemp[i] - j2;
       ib       = ibtemp[i] - j2;
       ic       = ictemp[i] - j2;
       ia       = (ia>0 ? ia:ngrid_a+ia);
       ib       = (ib>0 ? ib:ngrid_b+ib);
       ic       = (ic>0 ? ic:ngrid_c+ic);

       igrid_a[j][i] = 2*ia - 1;
       igrid_b[j][i] = 2*(ib - 1)*ngrid_a;
       igrid_c[j][i] = 2*(ic - 1)*ngrid_ab;

      }/*endfor*/
     }/*endfor*/
/*--------------------------------------------------------------------------*/ 
/* B) Initialize M2 and get the Mn's using the recursion relation           */ 
/*    Note: M_n is defined on 0 to n. Since frac between 0 and 1, the       */
/*          calculation is performed in an order that takes advantage of it */
     for(i=1;i<=nnow;i++){
      mn_a[1][i] = 1.0 - fabs(ua[1][i]-1.0);
      mn_b[1][i] = 1.0 - fabs(ub[1][i]-1.0);
      mn_c[1][i] = 1.0 - fabs(uc[1][i]-1.0);
      mn_a[2][i] = 1.0 - fabs(ua[2][i]-1.0);
      mn_b[2][i] = 1.0 - fabs(ub[2][i]-1.0);
      mn_c[2][i] = 1.0 - fabs(uc[2][i]-1.0);
     }/*endfor*/
     for(j=3;j<=n_interp;j++){
      for(i=1;i<=nnow;i++){
       mn_a[j][i]   = 0.0;
       mn_b[j][i]   = 0.0;
       mn_c[j][i]   = 0.0;
      }/*endfor*/
     }/*endfor*/

     for(n=3;n<=n_interp;n++){
       for(j=n;j>=2;j--){
        j1 = j-1;
        for(i=1;i<=nnow;i++){
         mn_a_tmp  = (ua[j][i]*mn_a[j][i]+(rn[n]-ua[j][i])*mn_a[j1][i])*rn1[n];
         mn_b_tmp  = (ub[j][i]*mn_b[j][i]+(rn[n]-ub[j][i])*mn_b[j1][i])*rn1[n];
         mn_c_tmp  = (uc[j][i]*mn_c[j][i]+(rn[n]-uc[j][i])*mn_c[j1][i])*rn1[n];
         mn_a[j][i] = mn_a_tmp;
         mn_b[j][i] = mn_b_tmp;
         mn_c[j][i] = mn_c_tmp;
        }/*end for: i*/
       }/*end for: j*/
       for(i=1;i<=nnow;i++){
        mn_a[1][i] = ua[1][i]*mn_a[1][i]*rn1[n];
        mn_b[1][i] = ub[1][i]*mn_b[1][i]*rn1[n];
        mn_c[1][i] = uc[1][i]*mn_c[1][i]*rn1[n];
      }/*endfor */
     }/*end for: k*/

/*--------------------------------------------------------------------------*/ 
/* C) Stuff Mn's on the grid:Sum is ordered to remove recursive dependencies*/ 
/*                           in qgrid. igrid_now is unique for fixed i      */ 
     for(jc=1;jc<=n_interp;jc++){
     for(jb=1;jb<=n_interp;jb++){
     for(ja=1;ja<=n_interp;ja++){
     for(i=1;i<=nnow;i++){
       itemp     = i+iatm1;
       igrid_now[ja][i] = igrid_a[ja][i]+igrid_b[jb][i]+igrid_c[jc][i];
       qgrid_now[ja][i] = mn_a[ja][i]*mn_b[jb][i]*mn_c[jc][i]*ewd_scr_q[itemp];
     }}/*endfor*/
     for(i=1;i<=nnow;i++){
#ifndef NO_PRAGMA
#pragma IVDEP
#endif
     for(ja=1;ja<=n_interp;ja++){
       qgrid[igrid_now[ja][i]] += qgrid_now[ja][i];
     }}}}/*end for: j*/
   }/*end for: iatm*/

/*==========================================================================*/
/* IV) Fourier Transform qgrid                                              */

   para_fft_gen3d_bck_to_g(qgrid,qgrid_tmp,para_fft_pkg3d); 

   if(np_forc>1){
     sngl_upack_coef(qgrid_tmp,&qgrid_tmp[(nktot+1)],qgrid,para_fft_pkg3d);
     Reduce_scatter(&qgrid_tmp[1],&qgrid_tmp_real[1],
                    &recv_counts_coef[1],MPI_DOUBLE,MPI_SUM,comm_forc);
     Reduce_scatter(&qgrid_tmp[(nktot+2)],&qgrid_tmp_imag[1],
                    &recv_counts_coef[1],MPI_DOUBLE,MPI_SUM,comm_forc);
   }else{
     sngl_upack_coef(qgrid_tmp_real,qgrid_tmp_imag,qgrid,para_fft_pkg3d);
   }/*endif*/

/*==========================================================================*/
/* VI) Compute the potential energy on the spherically cutoff grid.         */

   falp2 = 4.0*alp_ewd*alp_ewd;
   vol = getdeth(cell->hmat);
   tpi = 2.0*M_PI;
   pivol = vol/(4.0*M_PI);
   rvol  = 1.0/vol;  

   for(i=1;i<=9;i++){pvten_tmp[i]=0.0;}


   ncoef       = nktot+1;
   idiv        = ncoef / np_forc;
   irem        = ncoef % np_forc;
   ncoef_proc  = (myid_forc <  irem ? idiv+1 : idiv);
   icoef_strt  = (myid_forc <= irem ? myid_forc*(idiv+1)+1 : 
                     irem*(idiv+1)+1+(myid_forc-irem)*idiv);
   ncoef_use   = (myid_forc+1==np_forc ? ncoef_proc-1 : ncoef_proc);
   icoef_off   = icoef_strt-1;

   vnow = 0.0;
   for(i=1;i <= ncoef_use; ++i) {
     aka = tpi*( (double) kastr[(i+icoef_off)] );
     akb = tpi*( (double) kbstr[(i+icoef_off)] );
     akc = tpi*( (double) kcstr[(i+icoef_off)] );
     xk = (aka*cell_hmati[1] + akb*cell_hmati[2] + akc*cell_hmati[3]);
     yk = (aka*cell_hmati[4] + akb*cell_hmati[5] + akc*cell_hmati[6]);
     zk = (aka*cell_hmati[7] + akb*cell_hmati[8] + akc*cell_hmati[9]);
     g2 = xk*xk + yk*yk + zk*zk;
     preg = exp((-g2/falp2))/(g2*pivol);
     prep = -2.0*preg*((g2/falp2)+1.0)/g2;
     if(iperd != 3){
       preg += clus_corr_r[i]*rvol;
       prep += dclus_corr_r[i]*rvol;
     }/*endif*/
     smag       = (qgrid_tmp_real[i]*qgrid_tmp_real[i]
                  +qgrid_tmp_imag[i]*qgrid_tmp_imag[i])
                  *bweight_tot[(i+icoef_off)];
     vnow += (smag*preg);
     prep    *= smag;
     pvten_tmp[1] += prep*xk*xk;
     pvten_tmp[5] += prep*yk*yk;
     pvten_tmp[9] += prep*zk*zk;
     pvten_tmp[2] += prep*xk*yk;
     pvten_tmp[3] += prep*xk*zk;
     pvten_tmp[6] += prep*yk*zk;
     qgrid_tmp_real[i] *= (preg*bweight_tot[(i+icoef_off)]);
     qgrid_tmp_imag[i] *= (preg*bweight_tot[(i+icoef_off)]);
   }/*endfor*/
   if(ncoef_use<ncoef_proc){
     qgrid_tmp_real[(ncoef_proc)] = 0.0;
     qgrid_tmp_imag[(ncoef_proc)] = 0.0;
   }/*endif*/

   *vrecip += vnow;
   if(iperd!=3 && ncoef_use<ncoef_proc) {
     q_sum1 = dsum1(natm_tot,clatoms_q,1);
     *vrecip += (0.5*q_sum1*q_sum1*clus_corr_r[ncoef_proc]*rvol);
   }/*endif*/

   pvten_tmp[4]  = pvten_tmp[2];
   pvten_tmp[7]  = pvten_tmp[3];
   pvten_tmp[8]  = pvten_tmp[6];
   pvten_tmp[1] += vnow;
   pvten_tmp[5] += vnow;
   pvten_tmp[9] += vnow;


   for(i=1;i<=9;i++){
    pvten[i]     += (pvten_tmp[i]*wght_now);
   }/*endfor*/

   if(iget_pv_real_inter==1&&irespa==0){
    for(i=1;i<=9;i++){   
     pvten_tot[i] += (pvten_tmp[i]);
    }/*endfor*/
   }/*endif*/

/*==========================================================================*/
/* VII) Fourier Transform qgrid                                             */

   if(np_forc>1){
     for(i=1;i<=ncoef_proc;i++){qgrid_tmp[i]=qgrid_tmp_real[i];}
     Allgatherv(&(qgrid_tmp[1]),ncoef_proc,MPI_DOUBLE,&(qgrid_tmp_real[1]),
                &recv_counts_coef[1],&(recv_dspls_coef[1]),
                MPI_DOUBLE,0,comm_forc);
     for(i=1;i<=ncoef_proc;i++){qgrid_tmp[i]=qgrid_tmp_imag[i];}
     Allgatherv(&(qgrid_tmp[1]),ncoef_proc,MPI_DOUBLE,&(qgrid_tmp_imag[1]),
                &recv_counts_coef[1],&(recv_dspls_coef[1]),
                MPI_DOUBLE,0,comm_forc);
   }/*endif*/
   pme_sngl_pack_coef(qgrid_tmp_real,qgrid_tmp_imag,qgrid,para_fft_pkg3d);
   para_fft_gen3d_fwd_to_r(qgrid,qgrid_tmp,para_fft_pkg3d); 

/*==========================================================================*/
/* VIII) Calculate the force                                                */

   for(iatm = iatm_str;iatm <= iatm_end;++iatm){
      ewd_scr_fx[iatm] = 0.0;
      ewd_scr_fy[iatm] = 0.0;
      ewd_scr_fz[iatm] = 0.0;
   }/*endfor*/
   for(iatm=iatm_str;iatm<=iatm_end;iatm+=nlen_pme){
     iatm1 = iatm-1;
     iend = MIN(iatm_end,iatm1+nlen_pme);
     nnow = iend-iatm1;
/*--------------------------------------------------------------------------*/ 
/* A) Using current fraction, find the grid points on which M_n is non-zero */
     for(i=1;i<=nnow;i++){
      itemp     = i+iatm1;
      iatemp[i] = (int) (ewd_scr_x[itemp]);
      ibtemp[i] = (int) (ewd_scr_y[itemp]);
      ictemp[i] = (int) (ewd_scr_z[itemp]);
      frac_a[i] = ewd_scr_x[itemp] - (double) (iatemp[i]);
      frac_b[i] = ewd_scr_y[itemp] - (double) (ibtemp[i]);
      frac_c[i] = ewd_scr_z[itemp] - (double) (ictemp[i]);
     }/*endfor*/
     for(j=1;j<=n_interp;j++){
      for(i=1;i<=nnow;i++){
       ua[j][i]    = frac_a[i] + aj[j];
       ub[j][i]    = frac_b[i] + aj[j];
       uc[j][i]    = frac_c[i] + aj[j];
       j2       = j-2;
       ia       = iatemp[i] - j2;
       ib       = ibtemp[i] - j2;
       ic       = ictemp[i] - j2;
       ia       = (ia>0 ? ia:ngrid_a+ia);
       ib       = (ib>0 ? ib:ngrid_b+ib);
       ic       = (ic>0 ? ic:ngrid_c+ic);

       igrid_a[j][i] = 2*ia - 1;
       igrid_b[j][i] = 2*(ib - 1)*ngrid_a;
       igrid_c[j][i] = 2*(ic - 1)*ngrid_ab;

      }/*endfor*/
     }/*endfor*/
/*--------------------------------------------------------------------------*/ 
/* B) Initialize M2 and get the Mn's using the recursion relation           */ 
/*    Note: M_n is defined on 0 to n. Since frac between 0 and 1, the       */
/*          calculation is performed in an order that takes advantage of it */
     for(i=1;i<=nnow;i++){
      mn_a[1][i] = 1.0 - fabs(ua[1][i]-1.0);
      mn_b[1][i] = 1.0 - fabs(ub[1][i]-1.0);
      mn_c[1][i] = 1.0 - fabs(uc[1][i]-1.0);
      mn_a[2][i] = 1.0 - fabs(ua[2][i]-1.0);
      mn_b[2][i] = 1.0 - fabs(ub[2][i]-1.0);
      mn_c[2][i] = 1.0 - fabs(uc[2][i]-1.0);
     }/*endfor*/
     for(j=3;j<=n_interp;j++){
      for(i=1;i<=nnow;i++){
       mn_a[j][i]   = 0.0;
       mn_b[j][i]   = 0.0;
       mn_c[j][i]   = 0.0;
      }/*endfor*/
     }/*endfor*/

     for(n=3;n<=n_interp;n++){
       for(j=n;j>=2;j--){
        j1 = j-1;
        for(i=1;i<=nnow;i++){
         mn_a_tmp  = (ua[j][i]*mn_a[j][i]+(rn[n]-ua[j][i])*mn_a[j1][i])*rn1[n];
         mn_b_tmp  = (ub[j][i]*mn_b[j][i]+(rn[n]-ub[j][i])*mn_b[j1][i])*rn1[n];
         mn_c_tmp  = (uc[j][i]*mn_c[j][i]+(rn[n]-uc[j][i])*mn_c[j1][i])*rn1[n];
         mn_a[j][i] = mn_a_tmp;
         mn_b[j][i] = mn_b_tmp;
         mn_c[j][i] = mn_c_tmp;
        }/*end for: i*/
       }/*end for: j*/
       for(i=1;i<=nnow;i++){
        mn_a[1][i] = ua[1][i]*mn_a[1][i]*rn1[n];
        mn_b[1][i] = ub[1][i]*mn_b[1][i]*rn1[n];
        mn_c[1][i] = uc[1][i]*mn_c[1][i]*rn1[n];
       }/*endfor */
       if(n==(n_interp-1)){
        for(i=1;i<=nnow;i++){
         dmn_a[1][i] = mn_a[1][i];
         dmn_b[1][i] = mn_b[1][i];
         dmn_c[1][i] = mn_c[1][i];
        }/*endfor*/
        for(j=2;j<=n_interp;j++){    
         j1 = j-1;
         for(i=1;i<=nnow;i++){
          dmn_a[j][i] = mn_a[j][i] - mn_a[j1][i];
          dmn_b[j][i] = mn_b[j][i] - mn_b[j1][i];
          dmn_c[j][i] = mn_c[j][i] - mn_c[j1][i];
         }/*endfor*/
        }/*endfor*/
       }/*endif*/
     }/*end for: n*/
/*--------------------------------------------------------------------------*/ 
/* C) Calculate the force                                                   */ 
     for(jc=1;jc<=n_interp;jc++){
    for(jb=1;jb<=n_interp;jb++){
     for(ja=1;ja<=n_interp;ja++){
     for(i=1;i<=nnow;i++){
       igrid_now[ja][i] = igrid_a[ja][i]+igrid_b[jb][i]+igrid_c[jc][i];
     }}/*endfor*/
     for(ja=1;ja<=n_interp;ja++){
     for(i=1;i<=nnow;i++){
       qgrid_now[ja][i] = qgrid[igrid_now[ja][i]];
     }}/*endfor*/
     for(ja=1;ja<=n_interp;ja++){
     for(i=1;i<=nnow;i++){
        itemp = i+iatm1;
        atemp = dmn_a[ja][i]* mn_b[jb][i]* mn_c[jc][i]*grid_a;
        btemp =  mn_a[ja][i]*dmn_b[jb][i]* mn_c[jc][i]*grid_b;
        ctemp =  mn_a[ja][i]* mn_b[jb][i]*dmn_c[jc][i]*grid_c;
        qgrid_dx = atemp*cell_hmati[1]+btemp*cell_hmati[2]+ctemp*cell_hmati[3];
        qgrid_dy = atemp*cell_hmati[4]+btemp*cell_hmati[5]+ctemp*cell_hmati[6];
        qgrid_dz = atemp*cell_hmati[7]+btemp*cell_hmati[8]+ctemp*cell_hmati[9];
        qgrid_now[ja][i] *= ewd_scr_q[itemp];
        ewd_scr_fx[itemp] -= (qgrid_dx*qgrid_now[ja][i]);
        ewd_scr_fy[itemp] -= (qgrid_dy*qgrid_now[ja][i]);
        ewd_scr_fz[itemp] -= (qgrid_dz*qgrid_now[ja][i]); 
     }}}}/*end for: j*/
   }/*end for: iatm*/

/*==========================================================================*/
/* IX) Scatter the force back and get virial estimator if PIMD on           */

    if(iver_get==1){
#ifndef NO_PRAGMA
#pragma IVDEP
#endif
     for(iatm = iatm_str;iatm <= iatm_end;++iatm){
      ktemp = clatoms_ichrg[iatm];
      clatoms_fxt[ktemp] += ewd_scr_fx[iatm];
      clatoms_fyt[ktemp] += ewd_scr_fy[iatm];
      clatoms_fzt[ktemp] += ewd_scr_fz[iatm];
     }/*endfor*/
    }/*endif*/
    for(iatm = iatm_str;iatm <= iatm_end;++iatm){
      ewd_scr_fx[iatm] *= wght_now;
      ewd_scr_fy[iatm] *= wght_now;
      ewd_scr_fz[iatm] *= wght_now;
    }/*endfor*/
#ifndef NO_PRAGMA
#pragma IVDEP
#endif
    for(iatm = iatm_str;iatm <= iatm_end;++iatm){
      ktemp = clatoms_ichrg[iatm];
      clatoms_fx[ktemp] += ewd_scr_fx[iatm];
      clatoms_fy[ktemp] += ewd_scr_fy[iatm];
      clatoms_fz[ktemp] += ewd_scr_fz[iatm];
    }/*endfor*/

/*-----------------------------------------------------------------------*/ 
  }/*end routine*/
/*=======================================================================*/






/*==========================================================================*/
/* Full g */
/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void ewald3d_recip_pme_full_g(CLATOMS_INFO *clatoms_info,
                       CLATOMS_POS *clatoms_pos,
                       CELL *cell,PTENS *ptens,
                       double alp_ewd,int nktot,
                       int kastr[],int kbstr[],int kcstr[],
                       EWD_SCR *ewd_scr,double *vrecip,
                       double wght_now, int iver_get,PART_MESH *part_mesh,
                       int irespa,CLASS_COMM_FORC_PKG *class_comm_forc_pkg,
                       int iget_pv_real_inter,PARA_FFT_PKG3D *para_fft_pkg3d,
                       FOR_SCR *for_scr, 
                       double *clus_corr_r, double *dclus_corr_r)

/*==========================================================================*/
     {/*begin routine*/ 
/*==========================================================================*/
/*    Local Variables                                                       */

#include "../typ_defs/typ_mask.h"

  int iproc;
  int i,j,iatm,ktemp,n,j2,j1,inow;
  int ja,jb,jc;
  int ia,ib,ic;
  int nk1,nk2,nk3;   
  int skc_interp,j2_skc_interp;
  int ngrid_a,ngrid_b,ngrid_c,nrecip_grid;
  int ngrid_ab;
  int ngrid_bc,n_interp;
  int iatm1,iend,nnow,itemp;
  int n_ovlap,ioff_ovlap,n_ovlap1;
  int iii;
  int iatm_str,iatm_end;


  double mn_a_tmp,mn_b_tmp,mn_c_tmp,prep,vnow;
  double atemp,btemp,ctemp;
  double qgrid_dx,qgrid_dy,qgrid_dz;
  double falp2,vol,tpi,pivol,g2,preg,rvol;
  double aka,akb,akc,q_sum1;
  double xk,yk,zk,smag;
  double grid_a,grid_b,grid_c;

/* Define local pointers                                                  */

  int nchrg               = clatoms_info->nchrg;
  int natm_tot            = clatoms_info->natm_tot;
  int *clatoms_ichrg      = clatoms_info->ichrg;
  double *clatoms_q       = clatoms_info->q;
  double *clatoms_x       = clatoms_pos->x;
  double *clatoms_y       = clatoms_pos->y;
  double *clatoms_z       = clatoms_pos->z;
  double *clatoms_fx      = clatoms_pos->fx;
  double *clatoms_fy      = clatoms_pos->fy;
  double *clatoms_fz      = clatoms_pos->fz;
  double *clatoms_fxt     = clatoms_pos->fxt;
  double *clatoms_fyt     = clatoms_pos->fyt;
  double *clatoms_fzt     = clatoms_pos->fzt;
  double *cell_hmati      = cell->hmati;
  double *cell_hmat       = cell->hmat;
  double *pvten           = ptens->pvten;
  double *pvten_tot       = ptens->pvten_tot;
  double *pvten_tmp       = ptens->pvten_tmp;
  double *ewd_scr_x       = ewd_scr->x;
  double *ewd_scr_y       = ewd_scr->y;
  double *ewd_scr_z       = ewd_scr->z;
  double *ewd_scr_q       = ewd_scr->q;
  double *ewd_scr_fx      = ewd_scr->fx;
  double *ewd_scr_fy      = ewd_scr->fy;
  double *ewd_scr_fz      = ewd_scr->fz; 
  int *map_atm            = for_scr->iexcl;
  list_int *map_atm_tmp1       = for_scr->i_lnk;
  list_int *map_atm_tmp2       = for_scr->j_lnk;

  double *bweight_tot;
  int *nc                 = part_mesh->nc;
  int *ioff_c             = part_mesh->ioff_c;
  double *qgrid           = part_mesh->qgrid;
  double *qgrid_tmp       = part_mesh->qgrid_scr;
  double *aj              = part_mesh->aj;
  double *rn              = part_mesh->rn;
  double *rn1             = part_mesh->rn1;
  int *iatemp             = part_mesh->iatemp;
  int *ibtemp             = part_mesh->ibtemp;
  int *ictemp             = part_mesh->ictemp;
  double *frac_a          = part_mesh->frac_a;
  double *frac_b          = part_mesh->frac_b;
  double *frac_c          = part_mesh->frac_c;
  double *qgrid_tmp_real  = part_mesh->qgrid_tmp_real;
  double *qgrid_tmp_imag  = part_mesh->qgrid_tmp_imag;
  double **qgrid_now      = part_mesh->qgrid_now;
  double **ua             = part_mesh->ua;
  double **ub             = part_mesh->ub;
  double **uc             = part_mesh->uc;
  double **mn_a           = part_mesh->mn_a;
  double **mn_b           = part_mesh->mn_b;
  double **mn_c           = part_mesh->mn_c;
  double **dmn_a          = part_mesh->dmn_a;
  double **dmn_b          = part_mesh->dmn_b;
  double **dmn_c          = part_mesh->dmn_c;
  int **igrid_a           = part_mesh->igrid_a;
  int **igrid_b           = part_mesh->igrid_b;
  int **igrid_c           = part_mesh->igrid_c;
  int **igrid_now         = part_mesh->igrid_now;
  int nlen_pme            = part_mesh->nlen_pme;
  int np_forc             = class_comm_forc_pkg->num_proc;
  int myid_forc           = class_comm_forc_pkg->myid;
  MPI_Comm comm           = class_comm_forc_pkg->comm;
  int nfft_ka_proc        = para_fft_pkg3d->nfft_ka_proc;
  int skb_fft_ka_proc     = para_fft_pkg3d->skb_fft_ka_proc;
  int skc_fft_ka_proc     = para_fft_pkg3d->skc_fft_ka_proc;
  int ekb_fft_ka_proc     = para_fft_pkg3d->ekb_fft_ka_proc;
  int ekc_fft_ka_proc     = para_fft_pkg3d->ekc_fft_ka_proc;
  int ncoef               = para_fft_pkg3d->ncoef_proc;
  int ncoef_use           = para_fft_pkg3d->ncoef_use;
  int icoef_off           = para_fft_pkg3d->icoef_off;

  int *recvcounts_pme_dn  = para_fft_pkg3d->recvcounts_pme_dn;
  int *recvdspls_pme_dn   = para_fft_pkg3d->recvdspls_pme_dn;
  int *sendcounts_pme_dn  = para_fft_pkg3d->sendcounts_pme_dn;
  int *senddspls_pme_dn   = para_fft_pkg3d->senddspls_pme_dn;
  int *map_pme_dn         = para_fft_pkg3d->map_pme_dn;

  int nrecv_up_max        = para_fft_pkg3d->nrecv_up_max;
  int *recvcounts_pme_up  = para_fft_pkg3d->recvcounts_pme_up;
  int *recvdspls_pme_up   = para_fft_pkg3d->recvdspls_pme_up;
  int *sendcounts_pme_up  = para_fft_pkg3d->sendcounts_pme_up;
  int *senddspls_pme_up   = para_fft_pkg3d->senddspls_pme_up;
  int *map_pme_up         = para_fft_pkg3d->map_pme_up;
  int iperd               = cell->iperd; 

/*==========================================================================*/
/* 0) Respa dependent local pointers                                        */

  if(irespa==0){
    ngrid_a           = part_mesh->ngrid_a;
    ngrid_b           = part_mesh->ngrid_b;
    ngrid_c           = part_mesh->ngrid_c;
    n_interp          = part_mesh->n_interp;
    bweight_tot       = part_mesh->bweight_tot;
  }else{
    ngrid_a           = part_mesh->ngrid_a_res;
    ngrid_b           = part_mesh->ngrid_b_res;
    ngrid_c           = part_mesh->ngrid_c_res;
    n_interp          = part_mesh->n_interp_res;
    bweight_tot       = part_mesh->bweight_tot_res;
  }/*endif*/

/*==========================================================================*/
/* I) Construct some useful constants                                       */

  if(np_forc==1){
    n_ovlap     = 2*ngrid_a*(n_interp-1)*ngrid_b;
    ioff_ovlap  = 2*ngrid_a*ngrid_b*ngrid_c;
  }else{
    n_ovlap     = 2*ngrid_a*( (n_interp-1)*ngrid_b + skb_fft_ka_proc-1 );
    ioff_ovlap  = 2*ngrid_a*nfft_ka_proc;
  }/*endif*/
  nrecip_grid = ioff_ovlap + n_ovlap;

  grid_a         = (double) ngrid_a;
  grid_b         = (double) ngrid_b;
  grid_c         = (double) ngrid_c;
  ngrid_bc       = ngrid_b*ngrid_c;
  ngrid_ab       = ngrid_a*ngrid_b;
  skc_interp     = skc_fft_ka_proc - n_interp;

  rn1[1] = 0.0;
  for(j=1;j<=n_interp;j++){
    aj[j] = (double) (j-1);
    rn[j] = (double) (j);
    if(j > 1){rn1[j] = 1.0/((double)(j-1));}
  }/*endfor*/

/*==========================================================================*/
/* II.a) Find charged atoms, get their scaled imaged particle coordinates   */

  for(iatm = 1;iatm <= nchrg;++iatm){
     ktemp = clatoms_ichrg[iatm];
     ewd_scr_x[iatm] = clatoms_x[ktemp];
     ewd_scr_y[iatm] = clatoms_y[ktemp];
     ewd_scr_z[iatm] = clatoms_z[ktemp];
  }/*endfor*/

  for(iatm = 1;iatm <= nchrg;++iatm){
     atemp = ewd_scr_x[iatm]*cell_hmati[1]
           + ewd_scr_y[iatm]*cell_hmati[4]
           + ewd_scr_z[iatm]*cell_hmati[7];
     btemp = ewd_scr_x[iatm]*cell_hmati[2]
           + ewd_scr_y[iatm]*cell_hmati[5]
           + ewd_scr_z[iatm]*cell_hmati[8];
     ctemp = ewd_scr_x[iatm]*cell_hmati[3]
           + ewd_scr_y[iatm]*cell_hmati[6]
           + ewd_scr_z[iatm]*cell_hmati[9];
     atemp = atemp - NINT((atemp-0.5));
     btemp = btemp - NINT((btemp-0.5));
     ctemp = ctemp - NINT((ctemp-0.5));
     ewd_scr_x[iatm] = atemp*grid_a;
     ewd_scr_y[iatm] = btemp*grid_b;
     ewd_scr_z[iatm] = ctemp*grid_c;
  }/*endfor*/

/*==========================================================================*/
/* II.b Assign atoms to processors  based on their scaled coordinates       */

  if(np_forc>1){
     assign_atm_pme(nchrg,clatoms_ichrg,ewd_scr_x,ewd_scr_y,ewd_scr_z,
                    ewd_scr_fx,ewd_scr_fy,ewd_scr_fz,
                    map_atm,map_atm_tmp1,map_atm_tmp2,
                    nc,ioff_c,ngrid_c,skc_fft_ka_proc,ekc_fft_ka_proc,
                    skb_fft_ka_proc,ekb_fft_ka_proc,&iatm_str,&iatm_end);
     clatoms_ichrg = map_atm; /* relocal pointerize to rejiggered map*/
  }else{
     iatm_str=1;
     iatm_end =nchrg;
  }/*endif*/

#ifdef DEBUG_PME
  for(iproc=1;iproc<=np_forc;iproc++){
    if(iproc==myid_forc+1){
      printf("%d atm_strt %d atm_end %d\n",iproc,iatm_str,iatm_end);
    }/*endif*/
    Barrier(comm);
  }/*endfor*/
  if(myid_forc==0){scanf("%d",&iii);} Barrier(comm);
#endif

  for(iatm = iatm_str;iatm <=iatm_end;++iatm){
    ewd_scr_q[iatm] = clatoms_q[clatoms_ichrg[iatm]];
  }/*endfor*/

/*==========================================================================*/
/* III) Calculate the Cardinal B spline interpolation functions of the      */
/*     charge weighted density on the real space grid                       */

  for(i=1;i<=nrecip_grid;i++){
    qgrid[i]=0.0;
  }/*endfor*/

  for(iatm=iatm_str;iatm<=iatm_end;iatm+=nlen_pme){
    iatm1 = iatm-1;
    iend = MIN(iatm_end,iatm1+nlen_pme);
    nnow = iend-iatm1;
/*--------------------------------------------------------------------------*/ 
/* A) Using current fraction, find the grid points on which M_n is non-zero */
    for(i=1;i<=nnow;i++){
       itemp     = i+iatm1;
       iatemp[i] = (int) (ewd_scr_x[itemp]);
       ibtemp[i] = (int) (ewd_scr_y[itemp]);
       ictemp[i] = (int) (ewd_scr_z[itemp]);
       frac_a[i] = ewd_scr_x[itemp] - (double) (iatemp[i]);
       frac_b[i] = ewd_scr_y[itemp] - (double) (ibtemp[i]);
       frac_c[i] = ewd_scr_z[itemp] - (double) (ictemp[i]);
    }/*endfor*/
    for(j=1;j<=n_interp;j++){
       j2            = j-2;
       j2_skc_interp = j2+skc_interp;
       for(i=1;i<=nnow;i++){
         ua[j][i] = frac_a[i] + aj[j];
         ub[j][i] = frac_b[i] + aj[j];
         uc[j][i] = frac_c[i] + aj[j];
         ia       = iatemp[i] - j2;
         ib       = ibtemp[i] - j2;
         ic       = ictemp[i] - j2_skc_interp;
         ia       = (ia>0 ? ia:ngrid_a+ia);
         ib       = (ib>0 ? ib:ngrid_b+ib);
         igrid_a[j][i] = 2*ia - 1;
         igrid_b[j][i] = 2*(ib - 1)*ngrid_a;
         igrid_c[j][i] = 2*(ic - 1)*ngrid_ab;

       }/*endfor*/
    }/*endfor:interpolation pts*/
/*--------------------------------------------------------------------------*/ 
/* B) Initialize M2 and get the Mn's using the recursion relation           */ 
/*    Note: M_n is defined on 0 to n. Since frac between 0 and 1, the       */
/*          calculation is performed in an order that takes advantage of it */
/*  (Number of functions in an order = order ) */
/*  exs:  level 2 interpolation 2 fns, level 3 interpolation 3 fns*/
    for(i=1;i<=nnow;i++){
      mn_a[1][i] = 1.0 - fabs(ua[1][i]-1.0);
      mn_b[1][i] = 1.0 - fabs(ub[1][i]-1.0);
      mn_c[1][i] = 1.0 - fabs(uc[1][i]-1.0);
      mn_a[2][i] = 1.0 - fabs(ua[2][i]-1.0);
      mn_b[2][i] = 1.0 - fabs(ub[2][i]-1.0);
      mn_c[2][i] = 1.0 - fabs(uc[2][i]-1.0);
    }/*endfor*/
    for(j=3;j<=n_interp;j++){
      for(i=1;i<=nnow;i++){
        mn_a[j][i]   = 0.0;
        mn_b[j][i]   = 0.0;
        mn_c[j][i]   = 0.0;
      }/*endfor*/
    }/*endfor*/

    for(n=3;n<=n_interp;n++){/* Interpolation order */

      for(j=n;j>=2;j--){/* Number of non-zero functions in this order */
        j1 = j-1;
        for(i=1;i<=nnow;i++){
         mn_a_tmp  = (ua[j][i]*mn_a[j][i]+(rn[n]-ua[j][i])*mn_a[j1][i])*rn1[n];
         mn_b_tmp  = (ub[j][i]*mn_b[j][i]+(rn[n]-ub[j][i])*mn_b[j1][i])*rn1[n];
         mn_c_tmp  = (uc[j][i]*mn_c[j][i]+(rn[n]-uc[j][i])*mn_c[j1][i])*rn1[n];
         mn_a[j][i] = mn_a_tmp;
         mn_b[j][i] = mn_b_tmp;
         mn_c[j][i] = mn_c_tmp;
        }/*end for: i*/
      }/*end for: : jth interpolation pts*/
      for(i=1;i<=nnow;i++){
         mn_a[1][i] = ua[1][i]*mn_a[1][i]*rn1[n];
         mn_b[1][i] = ub[1][i]*mn_b[1][i]*rn1[n];
         mn_c[1][i] = uc[1][i]*mn_c[1][i]*rn1[n];
      }/*endfor 1st interpolation pt*/

    }/*end for: interpolation order*/

/*--------------------------------------------------------------------------*/ 
/* C) Stuff Mn's on the grid:Sum is ordered to remove recursive dependencies*/ 
/*                           in qgrid. igrid_now is unique for fixed i      */ 
    for(jc=1;jc<=n_interp;jc++){
    for(jb=1;jb<=n_interp;jb++){
    for(ja=1;ja<=n_interp;ja++){
    for(i=1;i<=nnow;i++){
       itemp     = i+iatm1;
       igrid_now[ja][i] = igrid_a[ja][i]+igrid_b[jb][i]+igrid_c[jc][i];
       qgrid_now[ja][i] = mn_a[ja][i]*mn_b[jb][i]*mn_c[jc][i]*ewd_scr_q[itemp];
    }}/*endfor*/
    for(i=1;i<=nnow;i++){
     for(ja=1;ja<=n_interp;ja++){
       qgrid[igrid_now[ja][i]] += qgrid_now[ja][i];
     }}}}/*end for: j*/
  }/*end for: iatm*/

/*--------------------------------------------------------------------------*/ 
/* D) Alltotall to add in particle spread on appropriate processors         */
/*    (even in scalar)                                                      */

#ifdef COMM_GRID
  if(np_forc>1){
    Alltoallv(&qgrid[1],
              &(sendcounts_pme_dn[1]),&(senddspls_pme_dn[1]),MPI_DOUBLE,
              &qgrid_tmp[1],
              &(recvcounts_pme_dn[1]),&(recvdspls_pme_dn[1]),MPI_DOUBLE,comm);
    for(i=1;i<=np_forc;i++){
      for(j=1;j<=recvcounts_pme_dn[i];j++){
        qgrid[(map_pme_dn[i]+j)] += qgrid_tmp[(j+recvdspls_pme_dn[i])];
      }/*endfor*/
    }/*endfor*/
    if(nrecv_up_max>0){
      Alltoallv(&qgrid[1],
               &(sendcounts_pme_up[1]),&(senddspls_pme_up[1]),MPI_DOUBLE,
               &qgrid_tmp[1],
               &(recvcounts_pme_up[1]),&(recvdspls_pme_up[1]),MPI_DOUBLE,comm);
      for(i=1;i<=np_forc;i++){
        for(j=1;j<=recvcounts_pme_up[i];j++){
          qgrid[(map_pme_up[i]+j)] += qgrid_tmp[(j+recvdspls_pme_up[i])];
        }/*endfor*/
      }/*endfor*/
    }/*endif*/
  } else {
    for(i=1;i<=n_ovlap;i++){qgrid[(i+ioff_ovlap)] += qgrid[i];}
  }/* endif */
#endif

/*==========================================================================*/
/* IV) Fourier Transform qgrid                                              */

#ifdef FFT
  para_fft_gen3d_bck_to_g(&qgrid[n_ovlap],qgrid_tmp,para_fft_pkg3d); 
#endif
  sngl_upack_coef(qgrid_tmp_real,qgrid_tmp_imag,&qgrid[n_ovlap],
                  para_fft_pkg3d);

/*==========================================================================*/
/* VI) Compute the potential energy on the spherically cutoff grid.         */

  falp2 = 4.0*alp_ewd*alp_ewd;
  vol = getdeth(cell_hmat);
  tpi = 2.0*M_PI;
  pivol = vol/(4.0*M_PI);
  rvol  = 1.0/vol;  

  for(i=1;i<=9;i++){pvten_tmp[i]=0.0;}

  vnow = 0.0;
  for(i=1;i <= ncoef_use; ++i) {
    aka   = tpi*( (double) kastr[(i+icoef_off)] );
    akb   = tpi*( (double) kbstr[(i+icoef_off)] );
    akc   = tpi*( (double) kcstr[(i+icoef_off)] );
    xk    = (aka*cell_hmati[1] + akb*cell_hmati[2] + akc*cell_hmati[3]);
    yk    = (aka*cell_hmati[4] + akb*cell_hmati[5] + akc*cell_hmati[6]);
    zk    = (aka*cell_hmati[7] + akb*cell_hmati[8] + akc*cell_hmati[9]);
    g2    = xk*xk + yk*yk + zk*zk;
    preg  = exp((-g2/falp2))/(g2*pivol);
    prep  = -2.0*preg*((g2/falp2)+1.0)/g2;
    if(iperd != 3){
      preg += clus_corr_r[i]*rvol;
      prep += dclus_corr_r[i]*rvol;
    }/*endif*/
    smag  = (qgrid_tmp_real[i]*qgrid_tmp_real[i]
           +qgrid_tmp_imag[i]*qgrid_tmp_imag[i])*bweight_tot[(i+icoef_off)];
    vnow += (smag*preg);
    prep *= smag;
    pvten_tmp[1] += prep*xk*xk;
    pvten_tmp[5] += prep*yk*yk;
    pvten_tmp[9] += prep*zk*zk;
    pvten_tmp[2] += prep*xk*yk;
    pvten_tmp[3] += prep*xk*zk;
    pvten_tmp[6] += prep*yk*zk;
    qgrid_tmp_real[i] *= (preg*bweight_tot[(i+icoef_off)]);
    qgrid_tmp_imag[i] *= (preg*bweight_tot[(i+icoef_off)]);
  }/*endfor*/
  if(ncoef_use < ncoef){
    qgrid_tmp_real[ncoef] = 0.0;
    qgrid_tmp_imag[ncoef] = 0.0;
  }/*endif*/

  *vrecip += vnow;
   if(iperd!=3 && ncoef_use < ncoef) {
     q_sum1 = dsum1(natm_tot,clatoms_q,1);
     *vrecip += (0.5*q_sum1*q_sum1*clus_corr_r[ncoef]*rvol);
   }/*endif*/

  pvten_tmp[4]  = pvten_tmp[2];
  pvten_tmp[7]  = pvten_tmp[3];
  pvten_tmp[8]  = pvten_tmp[6];
  pvten_tmp[1] += vnow;
  pvten_tmp[5] += vnow;
  pvten_tmp[9] += vnow;

  for(i=1;i<=9;i++){
    pvten[i]     += (pvten_tmp[i]*wght_now);
  }/*endfor*/

  if(iget_pv_real_inter==1&&irespa==0){
    for(i=1;i<=9;i++){   
      pvten_tot[i] += pvten_tmp[i];
    }/*endfor*/
  }/*endif*/

/*==========================================================================*/
/* VIIA) Fourier Transform qgrid                                             */

  pme_sngl_pack_coef(qgrid_tmp_real,qgrid_tmp_imag,&qgrid[n_ovlap],
                      para_fft_pkg3d);
#ifdef FFT
  para_fft_gen3d_fwd_to_r(&qgrid[n_ovlap],qgrid_tmp,para_fft_pkg3d); 
#endif

/*--------------------------------------------------------------------------*/ 
/* D) Alltotall to put particle spread on appropriate processors            */
/*    (even in scalar)                                                      */

#ifdef COMM_GRID
  if(np_forc > 1){
    if(nrecv_up_max>0){
      for(i=1;i<=np_forc;i++){
        for(j=1;j<=recvcounts_pme_up[i];j++){
          qgrid_tmp[(j+recvdspls_pme_up[i])] = qgrid[(map_pme_up[i]+j)];
        }/*endfor*/
      }/*endfor*/
      Alltoallv(&qgrid_tmp[1],
               &(recvcounts_pme_up[1]),&(recvdspls_pme_up[1]),MPI_DOUBLE,
               &qgrid[1],
               &(sendcounts_pme_up[1]),&(senddspls_pme_up[1]),MPI_DOUBLE,comm);
    }/*endif*/
    for(i=1;i<=np_forc;i++){
      for(j=1;j<=recvcounts_pme_dn[i];j++){
        qgrid_tmp[(j+recvdspls_pme_dn[i])] = qgrid[(map_pme_dn[i]+j)];
      }/*endfor*/
    }/*endfor*/
    Alltoallv(&qgrid_tmp[1],
              &(recvcounts_pme_dn[1]),&(recvdspls_pme_dn[1]),MPI_DOUBLE,
              &qgrid[1],
              &(sendcounts_pme_dn[1]),&(senddspls_pme_dn[1]),MPI_DOUBLE,comm);
  } else {
    for(i=1;i<=n_ovlap;i++){qgrid[i] = qgrid[(i+ioff_ovlap)];}
  }/* endif */
#endif

/*==========================================================================*/
/* VIII) Calculate the force                                                */

  for(iatm = iatm_str;iatm <= iatm_end;++iatm){
    ewd_scr_fx[iatm] = 0.0;
    ewd_scr_fy[iatm] = 0.0;
    ewd_scr_fz[iatm] = 0.0;
  }/*endfor*/
  for(iatm=iatm_str;iatm<=iatm_end;iatm+=nlen_pme){
     iatm1 = iatm-1;
     iend = MIN(iatm_end,iatm1+nlen_pme);
     nnow = iend-iatm1;
/*--------------------------------------------------------------------------*/ 
/* A) Using current fraction, find the grid points on which M_n is non-zero */
     for(i=1;i<=nnow;i++){
       itemp     = i+iatm1;
       iatemp[i] = (int) (ewd_scr_x[itemp]);
       ibtemp[i] = (int) (ewd_scr_y[itemp]);
       ictemp[i] = (int) (ewd_scr_z[itemp]);
       frac_a[i] = ewd_scr_x[itemp] - (double) (iatemp[i]);
       frac_b[i] = ewd_scr_y[itemp] - (double) (ibtemp[i]);
       frac_c[i] = ewd_scr_z[itemp] - (double) (ictemp[i]);
     }/*endfor*/
     for(j=1;j<=n_interp;j++){
       j2       = j-2;
       j2_skc_interp = j2+skc_interp;
       for(i=1;i<=nnow;i++){
         ua[j][i]    = frac_a[i] + aj[j];
         ub[j][i]    = frac_b[i] + aj[j];
         uc[j][i]    = frac_c[i] + aj[j];
         ia       = iatemp[i] - j2;
         ib       = ibtemp[i] - j2;
         ic       = ictemp[i] - j2_skc_interp;
         ia       = (ia>0 ? ia:ngrid_a+ia);
         ib       = (ib>0 ? ib:ngrid_b+ib);
         igrid_a[j][i] = 2*ia - 1;
         igrid_b[j][i] = 2*(ib - 1)*ngrid_a;
         igrid_c[j][i] = 2*(ic - 1)*ngrid_ab;
      }/*endfor*/
    }/* endfor */
/*--------------------------------------------------------------------------*/ 
/* B) Initialize M2 and get the Mn's using the recursion relation           */ 
/*    Note: M_n is defined on 0 to n. Since frac between 0 and 1, the       */
/*          calculation is performed in an order that takes advantage of it */
     for(i=1;i<=nnow;i++){
       mn_a[1][i] = 1.0 - fabs(ua[1][i]-1.0);
       mn_b[1][i] = 1.0 - fabs(ub[1][i]-1.0);
       mn_c[1][i] = 1.0 - fabs(uc[1][i]-1.0);
       mn_a[2][i] = 1.0 - fabs(ua[2][i]-1.0);
       mn_b[2][i] = 1.0 - fabs(ub[2][i]-1.0);
       mn_c[2][i] = 1.0 - fabs(uc[2][i]-1.0);
     }/*endfor*/
     for(j=3;j<=n_interp;j++){
      for(i=1;i<=nnow;i++){
        mn_a[j][i]   = 0.0;
        mn_b[j][i]   = 0.0;
        mn_c[j][i]   = 0.0;
      }/*endfor*/
     }/*endfor*/

     for(n=3;n<=n_interp;n++){
       for(j=n;j>=2;j--){
        j1 = j-1;
        for(i=1;i<=nnow;i++){
         mn_a_tmp  = (ua[j][i]*mn_a[j][i]+(rn[n]-ua[j][i])*mn_a[j1][i])*rn1[n];
         mn_b_tmp  = (ub[j][i]*mn_b[j][i]+(rn[n]-ub[j][i])*mn_b[j1][i])*rn1[n];
         mn_c_tmp  = (uc[j][i]*mn_c[j][i]+(rn[n]-uc[j][i])*mn_c[j1][i])*rn1[n];
         mn_a[j][i] = mn_a_tmp;
         mn_b[j][i] = mn_b_tmp;
         mn_c[j][i] = mn_c_tmp;
        }/*end for: i*/
       }/*end for: j*/
       for(i=1;i<=nnow;i++){
         mn_a[1][i] = ua[1][i]*mn_a[1][i]*rn1[n];
         mn_b[1][i] = ub[1][i]*mn_b[1][i]*rn1[n];
         mn_c[1][i] = uc[1][i]*mn_c[1][i]*rn1[n];
       }/*endfor */
       if(n==(n_interp-1)){
        for(i=1;i<=nnow;i++){
          dmn_a[1][i] = mn_a[1][i];
          dmn_b[1][i] = mn_b[1][i];
          dmn_c[1][i] = mn_c[1][i];
        }/*endfor*/
        for(j=2;j<=n_interp;j++){    
         j1 = j-1;
         for(i=1;i<=nnow;i++){
           dmn_a[j][i] = mn_a[j][i] - mn_a[j1][i];
           dmn_b[j][i] = mn_b[j][i] - mn_b[j1][i];
           dmn_c[j][i] = mn_c[j][i] - mn_c[j1][i];
         }/*endfor*/
        }/*endfor*/
       }/*endif*/
     }/*end for: n*/
/*--------------------------------------------------------------------------*/ 
/* C) Calculate the force                                                   */ 
     for(jc=1;jc<=n_interp;jc++){
     for(jb=1;jb<=n_interp;jb++){
     for(ja=1;ja<=n_interp;ja++){
     for(i=1;i<=nnow;i++){
       igrid_now[ja][i] = igrid_a[ja][i]+igrid_b[jb][i]+igrid_c[jc][i];
     }}/*endfor*/
     for(ja=1;ja<=n_interp;ja++){
     for(i=1;i<=nnow;i++){
       qgrid_now[ja][i] = qgrid[igrid_now[ja][i]];
     }}/*endfor*/
     for(ja=1;ja<=n_interp;ja++){
     for(i=1;i<=nnow;i++){
       itemp = i+iatm1;
       atemp = dmn_a[ja][i]* mn_b[jb][i]* mn_c[jc][i]*grid_a;
       btemp =  mn_a[ja][i]*dmn_b[jb][i]* mn_c[jc][i]*grid_b;
       ctemp =  mn_a[ja][i]* mn_b[jb][i]*dmn_c[jc][i]*grid_c;
       qgrid_dx = atemp*cell_hmati[1]+btemp*cell_hmati[2]+ctemp*cell_hmati[3];
       qgrid_dy = atemp*cell_hmati[4]+btemp*cell_hmati[5]+ctemp*cell_hmati[6];
       qgrid_dz = atemp*cell_hmati[7]+btemp*cell_hmati[8]+ctemp*cell_hmati[9];
       qgrid_now[ja][i] *= ewd_scr_q[itemp];
       ewd_scr_fx[itemp] -= (qgrid_dx*qgrid_now[ja][i]);
       ewd_scr_fy[itemp] -= (qgrid_dy*qgrid_now[ja][i]);
       ewd_scr_fz[itemp] -= (qgrid_dz*qgrid_now[ja][i]); 
     }}}}/*end for: j*/
  }/*end for: iatm*/

/*==========================================================================*/
/* IX) Scatter the force back and get virial estimator if PIMD on           */

  if(iver_get==1){
    for(iatm = iatm_str;iatm <= iatm_end;++iatm){
      ktemp = clatoms_ichrg[iatm];
      clatoms_fxt[ktemp] += ewd_scr_fx[iatm];
      clatoms_fyt[ktemp] += ewd_scr_fy[iatm];
      clatoms_fzt[ktemp] += ewd_scr_fz[iatm];
    }/*endfor*/
  }/*endif iver_get */

  for(iatm = iatm_str;iatm <= iatm_end;++iatm){
    ewd_scr_fx[iatm] *= wght_now;
    ewd_scr_fy[iatm] *= wght_now;
    ewd_scr_fz[iatm] *= wght_now;
  }/*endfor*/
  for(iatm = iatm_str;iatm <= iatm_end;++iatm){
    ktemp = clatoms_ichrg[iatm];
    clatoms_fx[ktemp] += ewd_scr_fx[iatm];
    clatoms_fy[ktemp] += ewd_scr_fy[iatm];
    clatoms_fz[ktemp] += ewd_scr_fz[iatm];
  }/*endfor*/

/*-----------------------------------------------------------------------*/ 
  }/*end routine*/
/*=======================================================================*/






/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void assign_atm_pme(int n,int *index_orig,double *sa,double *sb,double *sc,
                    double *sa_tmp, double *sb_tmp, double *sc_tmp,
                    int *index,list_int *index_tmp,list_int *index_tmp2,
                    int *nc,int *ioff_c, int nkf_c,int skc,int ekc, 
                    int skb,int ekb,int *iatm_str_ret,int *iatm_end_ret) 

/*==========================================================================*/
     {/*begin routine*/ 
/*==========================================================================*/
/*    Local Variables                                                       */

  int i,jjj,isum,iii;
  int iatm_str,iatm_end;

/*==========================================================================*/
/* I) Sort along sc keeping sb and sa commensurate */

  for(i=1;i<=n;i++){index[i]=i;}
  sort_pme(n,sc,index);
  for(i=1;i<=n;i++){
    sa_tmp[i] = sa[index[i]];
    sb_tmp[i] = sb[index[i]];
    sc_tmp[i] = sc[i];
  }/*endif*/

#ifdef DEBUG_SORT1
  printf("===========\n");
  printf("First sort \n");
  printf("===========\n");
  for(i=1;i<=n;i++){
    printf("sc[%d]=%g , sb[%d]=%g \n",i,sc_tmp[i],i,sb_tmp[i]); 
  }/*endfor*/
  scanf("%d",&iii);
#endif

/*==========================================================================*/
/* II) Sort along sb keeping sa commensurate for constant sc                */

  for(i=1;i<=nkf_c;i++){nc[i]=0;}  /* count atm with the same c */
  for(i=1;i<=n;i++){jjj = sc[i]+1;nc[jjj] += 1;}
  for(i=1;i<=n;i++){index_tmp[i]=i;} 

  isum = 0; /* sort along b for atms with the same c */
  for(i=1;i<=nkf_c;i++){
    if(nc[i]> 1){
#ifdef T3E_SCILIB
      sort_pme_short(nc[i],&sb_tmp[isum],&index_tmp[isum]);      
#else
      sort_pme(nc[i],&sb_tmp[isum],&index_tmp[isum]);      
#endif
    }/*endif*/
    isum += nc[i];
  }/*endfor*/

#ifdef DEBUG_SORT2
  printf("===========\n");
  printf("Second sort *\n");
  printf("===========\n");
  for(i=1;i<=n;i++){
    printf("sc[%d]=%g , sb[%d]=%g \n",i,sc_tmp[i],i,sb_tmp[i]); 
  }/*endfor*/
  scanf("%d",&iii);
#endif

  for(i=1;i<=n;i++){
    sa[i]         = sa_tmp[index_tmp[i]];
    sb[i]         = sb_tmp[i];
    sc[i]         = sc_tmp[index_tmp[i]];
    index_tmp2[i] = index[index_tmp[i]];
  }/*endfor*/

/*==========================================================================*/
/* II) Rejigger the indicies to reference original variables                */

  for(i=1;i<=n;i++){index[i] = index_orig[index_tmp2[i]];}

/*==========================================================================*/
/* IV) Assign atoms to procs and account for PME spreading                 */
/*      Find upper limit of leak and lower limit of atoms that leak         */
/*             {n_interp_str_c[i],n_interp_end_c[i]}                        */
/*             {n_interp_str_b[i],n_interp_end_b[i]}                        */
  ioff_c[1] = 0;
  for(i=2;i<=nkf_c;i++){ioff_c[i]=nc[(i-1)]+ioff_c[(i-1)];}

  iatm_str = ioff_c[skc]+1;
  for(i=1;i<=nc[skc];i++){ 
    jjj = (int)(sb[(i+ioff_c[skc])])+1;
    if(jjj>=skb) break;
    iatm_str++;
  }/* endfor */

  iatm_end = ioff_c[ekc]+nc[ekc];
  for(i=nc[ekc];i>=1;i--){
    jjj = (int) (sb[(i+ioff_c[ekc])])+1;
    if(jjj<=ekb) break;
    iatm_end--;
  }/* endfor */

/*==========================================================================*/
/* IV) Set return values */

  *iatm_str_ret = iatm_str;
  *iatm_end_ret = iatm_end;

/*--------------------------------------------------------------------------*/ 
  }/*end routine*/
/*==========================================================================*/




/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void sort_pme(int n, double index[],int jndex[])

/*=======================================================================*/
/*            Begin subprogram:                                          */
   {/*begin routine*/
/*=======================================================================*/
/*          Local variable declarations                                  */

  int m,ir,i,j,rjndex,iii;
  int k,*kndex,*mask,isub,temp;
  double rindex;

/*=======================================================================*/
/* I) Setup                        */

  m  = n/2+1;
  ir = n;

/*=======================================================================*/
/* II) Sort array index keeping jndex commensurrate */

  for(;;){

/*---------------------------------------------------------------------*/
/*  A)hire rindex */
    if(m>1){ 
      m--;
      rindex = index[m];
      rjndex = jndex[m];
/*--------------------------------------------------------------------*/
/*  B)retire/promote index[1] */
    }else{
      rindex = index[ir];
      rjndex = jndex[ir];
      index[ir]=index[1];
      jndex[ir]=jndex[1];
      ir--;
      if(ir==1){
	index[1]=rindex;
	jndex[1]=rjndex;
	break;
      }/*endif*/
    }/*endif*/
/*---------------------------------------------------------------------*/
/*  C)put rindex in appropriate slot */
    i=m;
    j=2*m;
    while(j<=ir){
      /*    a)compare to rindex to underling */
      if((j<ir) && (index[j]< index[(j+1)])) j++;
      /*    b)demote */
      if(rindex<index[j]){
	index[i]=index[j];
	jndex[i]=jndex[j];
	i=j;
	j=2*j;
      }else{
	/*    c)if no demotations exit while */
	j=ir+1;
      }/*endif*/
    } /*endwhile*/
    /*    d)slot rindex */
    index[i] = rindex;
    jndex[i] = rjndex;
  }/*endfor*/

/*-----------------------------------------------------------------------*/
  } /*end routine*/ 
/*==========================================================================*/


/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void sort_pme_short(int n, double index[],list_int jndex[])

/*=======================================================================*/
/*            Begin subprogram:                                          */
   {/*begin routine*/
/*=======================================================================*/
/*          Local variable declarations                                  */

  int m,ir,i,j,rjndex,iii;
  int k,*kndex,*mask,isub,temp;
  double rindex;

/*=======================================================================*/
/* I) Setup                        */

  m  = n/2+1;
  ir = n;

/*=======================================================================*/
/* II) Sort array index keeping jndex commensurrate */

  for(;;){

/*---------------------------------------------------------------------*/
/*  A)hire rindex */
    if(m>1){ 
      m--;
      rindex = index[m];
      rjndex = jndex[m];
/*--------------------------------------------------------------------*/
/*  B)retire/promote index[1] */
    }else{
      rindex = index[ir];
      rjndex = jndex[ir];
      index[ir]=index[1];
      jndex[ir]=jndex[1];
      ir--;
      if(ir==1){
	index[1]=rindex;
	jndex[1]=rjndex;
	break;
      }/*endif*/
    }/*endif*/
/*---------------------------------------------------------------------*/
/*  C)put rindex in appropriate slot */
    i=m;
    j=2*m;
    while(j<=ir){
      /*    a)compare to rindex to underling */
      if((j<ir) && (index[j]< index[(j+1)])) j++;
      /*    b)demote */
      if(rindex<index[j]){
	index[i]=index[j];
	jndex[i]=jndex[j];
	i=j;
	j=2*j;
      }else{
	/*    c)if no demotations exit while */
	j=ir+1;
      }/*endif*/
    } /*endwhile*/
    /*    d)slot rindex */
    index[i] = rindex;
    jndex[i] = rjndex;
  }/*endfor*/

/*-----------------------------------------------------------------------*/
  } /*end routine*/ 
/*==========================================================================*/
