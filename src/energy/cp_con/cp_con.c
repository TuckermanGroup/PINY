/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
/*                                                                          */
/*                         PI_MD:                                           */
/*             The future of simulation technology                          */
/*             ------------------------------------                         */
/*                     Module: cp_con.c                                     */
/*                                                                          */
/* Function performs orthogonalization of the wave functions                */
/*                                                                          */
/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/


#include "standard_include.h"
#include "../typ_defs/typedefs_gen.h"
#include "../typ_defs/typedefs_cp.h"
#include "../proto_defs/proto_energy_cpcon_local.h"
#include "../proto_defs/proto_friend_lib_entry.h"
#include "../proto_defs/proto_math.h"
#include "../proto_defs/proto_communicate_wrappers.h"

#define DEBUG_GS_OFF

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void cp_shake(double *creal,double *cimag,  int icoef_form, int icoef_orth,
              double *vcreal,double *vcimag,int ivcoef_form, 
              double *creal_old,double *cimag_old,int icoef_form_old,
              double c_tolshake,double dt,
              CPSCR_OVMAT *cpscr_ovmat, 
              int *ioff,double *occ,double *rocc_sum,int cp_norb,
              int *iter_shake_cp,CP_COMM_STATE_PKG *cp_comm_state_pkg)

/*========================================================================*/
{/*begin routine*/
/*========================================================================*/
/*             Local variable declarations                                */
   int i,j,is,js,nstate2;
   int ind,iii;
   double alpha,beta,g0_wght;

/*             Local pointer declarations                                */
   int np_states  = cp_comm_state_pkg->num_proc;
   int myid_state = cp_comm_state_pkg->myid;
   int nstate  = cp_comm_state_pkg->nstate;
   double *a   = cpscr_ovmat->ovlap1;
   double *b   = cpscr_ovmat->ovlap2;
   double *p1  = cpscr_ovmat->ovlap3;
   double *p2  = cpscr_ovmat->ovlap4;
   double *pscr = cpscr_ovmat->ovlap5;  /* don't panic */
   double *xl2 = cpscr_ovmat->ovlap5;
   double *xlamb = cpscr_ovmat->ovlap6;
   double *xlold  = cpscr_ovmat->ovlap7;

/*========================================================================*/
/* 0) Checks                                                              */


   if(cp_norb>0){
    if(icoef_orth !=0){
      printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
      printf("The coefficients must be in non-orthogonal form\n");
      printf("on state processor %d in cp_shake \n",myid_state);
      printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
      fflush(stdout);exit(1);
    }/*endif*/
   }/*endif*/

   if(np_states>1){
    if((icoef_form+ivcoef_form+icoef_form_old) !=3){
      printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
      printf("The coefficients must be in transposed form\n");
      printf("on state processor %d in cp_shake \n",myid_state);
      printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
      fflush(stdout);exit(1);
    }/*endif*/
   }/*endif*/

/*========================================================================*/
/* I) Calculate overlap matrices between guess coeffs and old coeffs      */
   /* i) Calculate some useful constants and malloc some useful memory    */

    nstate2 = (nstate)*(nstate);

/*------------------------------------------------------------------------*/
   /* ii) Zero matrices                                             */

   for(is=1;is <= nstate2;is++){
      a[is] = b[is] = 0.0;
   }/* endfor */

/*------------------------------------------------------------------------*/
   /* iii) Assign diagonal elements according to occupation numbers   */

     for(is=1;is <= nstate;is++){
         ind = (is-1)*(nstate) + is;
         a[ind] = b[ind] = occ[is];
     }/* endfor */

/*------------------------------------------------------------------------*/
   /* iv) Calculate overlap integrals and construct matrices in equation  */

  alpha=2.0;
  cp_ovlap_mat_same(creal,cimag,icoef_form,p1,pscr,alpha,ioff,
                    cp_comm_state_pkg);
  cp_ovlap_mat_diff(creal,cimag,icoef_form,creal_old,cimag_old,icoef_form_old,
                    p2,pscr,alpha,ioff,cp_comm_state_pkg);

  for(is=1;is <= nstate2; is++){
   a[is] -= p1[is];
   b[is] -= p2[is];
  }/*endfor*/

/*========================================================================*/
/* II) Iterative solution of matrix equation                              */

  /* i) Initialization of Lagrange multiplier matrix    */

  for(is=1;is <=nstate2; is++){
     xlold[is] = xlamb[is] = rocc_sum[is]*a[is];
  }/* endfor */

/*------------------------------------------------------------------------*/
  /* ii) Set up matrices and iteratively solve      */

  cp_iter_mat_solve_shake(xlamb,xl2,xlold,a,b,p1,p2,occ,rocc_sum,
                          &nstate,iter_shake_cp,c_tolshake);

/*========================================================================*/
/* III) Apply constraint force to coefficients and coefficient velocities */

  alpha = 1.0;  beta = 1.0; g0_wght = 1.0;
  cp_rotation_prim(creal_old,icoef_form_old,creal,icoef_form,xlamb,
                   alpha,beta,g0_wght,ioff,cp_comm_state_pkg);
  cp_rotation_prim(cimag_old,icoef_form_old,cimag,icoef_form,xlamb,
                   alpha,beta,g0_wght,ioff,cp_comm_state_pkg);

  alpha = 1.0/dt;
  cp_rotation_prim(creal_old,icoef_form_old,vcreal,ivcoef_form,xlamb,
                   alpha,beta,g0_wght,ioff,cp_comm_state_pkg);
  cp_rotation_prim(cimag_old,icoef_form_old,vcimag,ivcoef_form,xlamb,
                   alpha,beta,g0_wght,ioff,cp_comm_state_pkg);

/*========================================================================*/
    }/*end routine*/
/*========================================================================*/




/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void cp_shake_mass(double *creal,double *cimag,int icoef_form, int icoef_orth,
                   double *vcreal,double *vcimag,int ivcoef_form,
                   double *creal_old,double *cimag_old,int icoef_form_old,
                   double *cmass,
                   double c_tolshake,double dt,
                   CPSCR_OVMAT *cpscr_ovmat, 
                   int *ioff,double *occ,double *rocc_sum,int cp_norb,
                   int *iter_shake_cp,
                   CP_COMM_STATE_PKG *cp_comm_state_pkg)

/*========================================================================*/
{/*begin routine*/
/*========================================================================*/
/*             Local variable declarations                                */
   int i,is,js,nstate2;
   int icoef,ind;
   int iii;
   double alpha,beta,g0_wght;
   double cm1,scale,rscale;

/*             Local pointer declarations                                */
   int np_states    = cp_comm_state_pkg->num_proc;
   int nstate       = cp_comm_state_pkg->nstate;
   int nstate_proc  = cp_comm_state_pkg->nstate_proc;
   int ncoef_tot      = cp_comm_state_pkg->ncoef;
   int ncoef        = cp_comm_state_pkg->nstate_ncoef_proc;
   int icmoff       = cp_comm_state_pkg->icoef_start-1;
   int myid_state = cp_comm_state_pkg->myid;
   double *a        = cpscr_ovmat->ovlap1;
   double *b        = cpscr_ovmat->ovlap2;
   double *c        = cpscr_ovmat->ovlap3;
   double *p1       = cpscr_ovmat->ovlap4;
   double *p2       = cpscr_ovmat->ovlap5;
   double *pscr     = cpscr_ovmat->ovlap6; /* don't panic*/
   double *xl2      = cpscr_ovmat->ovlap6;
   double *xlamb    = cpscr_ovmat->ovlap7;
   double *xlold    = cpscr_ovmat->ovlap8;

/*========================================================================*/
/* 0) Checks                                                              */

   if(cp_norb>0){
    if(icoef_orth !=0){
      printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
      printf("The coefficients must be in non-orthogonal form\n");
      printf("on state processor %d in cp_shake_mass \n",myid_state);
      printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
      fflush(stdout);exit(1);
    }/*endif*/
   }/*endif*/

   if(np_states>1){
    if((icoef_form+ivcoef_form+icoef_form_old) !=3){
      printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
      printf("The coefficients must be in transposed form\n");
      printf("on state processor %d in cp_shake_mass \n",myid_state);
      printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
      fflush(stdout);exit(1);
    }/*endif*/
   }/*endif*/

/*========================================================================*/
/* I) Calculate overlap matrices between guess coeffs and old coeffs      */
   /* i) Calculate some useful constants                            */
   nstate2 = nstate*nstate;

/*------------------------------------------------------------------------*/
   /* ii) Zero matrices                                             */
  for(is=1;is <=nstate2; is++){
   a[is] = b[is] = c[is] = 0.0;
  }/* endfor */

/*------------------------------------------------------------------------*/
   /* iii) Assign diagonal elements according to occupation numbers   */

     for(is=1;is <= nstate;is++){
         ind = (is-1)*nstate + is;
         a[ind] = b[ind] = occ[is];
     }/* endfor */

/*------------------------------------------------------------------------*/
   /* iv) Scale old coefficients by cmass[ncoef]/cmass[i]   */

   cm1 = cmass[ncoef_tot];
   for(is=1;is<=nstate;is++) {
     for(i=1;i<=ncoef;i++) {
       icoef = i+ioff[is];
       scale = cm1/cmass[(i+icmoff)];       
       creal_old[icoef] *=scale;
       cimag_old[icoef] *=scale;
     }/* endfor i */
   }/* endfor is */

/*------------------------------------------------------------------------*/
   /* v) Get overlap matrices                 */

  alpha=2.0;
  cp_ovlap_mat_same(creal,cimag,icoef_form,p1,pscr,alpha,ioff,
                    cp_comm_state_pkg);
  cp_ovlap_mat_diff(creal,cimag,icoef_form,creal_old,cimag_old,icoef_form_old,
                    p2,pscr,alpha,ioff,cp_comm_state_pkg);
  cp_ovlap_mat_same(creal_old,cimag_old,icoef_form_old,
                    c,pscr,alpha,ioff,cp_comm_state_pkg);

  for(is=1;is <= nstate2; is++){
   a[is] -= p1[is];
   b[is] -= p2[is];
  }/* endfor */

/*========================================================================*/
/* II) Iterative solution of matrix equation                              */

  /* i) Initialization of Lagrange multiplier matrix    */

  for(is=1;is <=nstate2; is++){
     xlold[is] = xlamb[is] = rocc_sum[is]*a[is];
  }/* endfor */

/*------------------------------------------------------------------------*/
  /* ii) Set up matrices and iterative do structure     */

  cp_iter_mat_solve_shake_mass(xlamb,xl2,xlold,a,b,c,p1,p2,rocc_sum,
                               &nstate,iter_shake_cp,c_tolshake);

/*========================================================================*/
/* III) Apply constraint force to coefficients and coefficient velocities */

  alpha = 1.0;  beta = 1.0; g0_wght = 1.0;
  cp_rotation_prim(creal_old,icoef_form_old,creal,icoef_form,xlamb,
                   alpha,beta,g0_wght,ioff,cp_comm_state_pkg);
  cp_rotation_prim(cimag_old,icoef_form_old,cimag,icoef_form,xlamb,
                   alpha,beta,g0_wght,ioff,cp_comm_state_pkg);

  alpha = 1.0/dt;
  cp_rotation_prim(creal_old,icoef_form_old,vcreal,ivcoef_form,xlamb,
                   alpha,beta,g0_wght,ioff,cp_comm_state_pkg);
  cp_rotation_prim(cimag_old,icoef_form_old,vcimag,ivcoef_form,xlamb,
                   alpha,beta,g0_wght,ioff,cp_comm_state_pkg);

/*========================================================================*/
/* IV) Unscale coefficients to return to original state: a nice touch     */

   cm1 = cmass[ncoef_tot];
   for(is=1;is<=nstate;is++) {
     for(i=1;i<=ncoef;i++) {
       icoef = i+ioff[is];
       scale = cm1/cmass[(i+icmoff)];       
       creal_old[icoef] /=scale;
       cimag_old[icoef] /=scale;
     }/* endfor i */
   }/* endfor is */

/*========================================================================*/
}/*end routine*/
/*========================================================================*/




/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void cp_shake_norb(double *creal,double *cimag, int icoef_form,int icoef_orth,
                   double *vcreal,double *vcimag,int ivcoef_form,
                   double *creal_old,double *cimag_old,int icoef_form_old,
                   double dt,int *ioff,double *occ,
                   double *cmass,int icmass_unif,
                   CPSCR_OVMAT *cpscr_ovmat,
                   CP_COMM_STATE_PKG *cp_comm_state_pkg)

/*========================================================================*/
{/*begin routine*/
/*========================================================================*/
/*             Local variable declarations                                */
#include "../typ_defs/typ_mask.h"

  double rdt,wght,descrim;
  double eps=1.e-10;
  int iii,i,icoef,icf,is,i1,i2,indm;
  int np_states   = cp_comm_state_pkg->num_proc;
  int nstate      = cp_comm_state_pkg->nstate;
  int ncoef       = cp_comm_state_pkg->nstate_ncoef_proc;
  int ncoef_tot   = cp_comm_state_pkg->ncoef;
  int myid_state  = cp_comm_state_pkg->myid;
  int icmoff       = cp_comm_state_pkg->icoef_start-1;
  double *a       = cpscr_ovmat->state_vec1;
  double *b       = cpscr_ovmat->state_vec2;
  double *c       = cpscr_ovmat->state_vec3;
  double *xl_norb = cpscr_ovmat->state_vec4;
  double *a_tmp   = cpscr_ovmat->state_vec5;
  double *b_tmp   = cpscr_ovmat->state_vec6;
  double *c_tmp   = cpscr_ovmat->state_vec7;
  MPI_Comm comm   = cp_comm_state_pkg->comm;

/*========================================================================*/
/* 0) Checks                                                              */

   if(icoef_orth !=0){
      printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
      printf("The coefficients must be in non-orthogonal form\n");
      printf("on state processor %d in cp_shake_norb \n",myid_state);
      printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
      fflush(stdout);exit(1);
   }/*endif*/

   if(np_states>1){
    if((icoef_form+ivcoef_form+icoef_form_old) !=3){
      printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
      printf("The coefficients must be in transposed form\n");
      printf("on state processor %d in cp_shake_norb \n",myid_state);
      printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
      fflush(stdout);exit(1);
    }/*endif*/
   }/*endif*/

/*========================================================================*/
/* II) Calculate norms of guess wave functions                            */
/*------------------------------------------------------------------------*/
  /* i) Zero vectors                 */

  for(is=1;is <= nstate;is++){a[is] = b[is] = c[is] = 0.0;}

/*------------------------------------------------------------------------*/
  /* i.v) Scale old coefficients if necessary                             */

  if(icmass_unif != 1) {
   double cm1,scale;
   cm1 = cmass[ncoef_tot];
   for(is=1;is<=nstate;is++) {
     for(i=1;i<=ncoef;i++) {
       icoef = i+ioff[is];
       scale = cm1/cmass[(i+icmoff)];       
       creal_old[icoef] *=scale;
       cimag_old[icoef] *=scale;
     }/* endfor i */
   }/* endfor is */
  }/* endif */

/*------------------------------------------------------------------------*/
  /* ii) Get norms                 */

  wght = 2.0;
  if(myid_state+1==np_states){wght=1.0;}
  for(is=1;is<=nstate;is++) {
   i1 = 1+ioff[is];
   i2 = ncoef-1+ioff[is];
   for(icf=i1;icf<=i2;icf++) {
     a[is] += creal[icf]*creal[icf]         + cimag[icf]*cimag[icf];
     b[is] += creal[icf]*creal_old[icf]     + cimag[icf]*cimag_old[icf];
     c[is] += creal_old[icf]*creal_old[icf] + cimag_old[icf]*cimag_old[icf];
   }/* endfor i */
   icf = ncoef+ioff[is];
   a[is] = 2.0*a[is] + wght*(creal[icf]*creal[icf]+cimag[icf]*cimag[icf]);
   b[is] = 2.0*b[is] + wght*(creal[icf]*creal_old[icf]+
                                           cimag[icf]*cimag_old[icf]);
   c[is] = 2.0*c[is] + wght*(creal_old[icf]*creal_old[icf]+
                                 cimag_old[icf]*cimag_old[icf]);
  }/* endfor is */

  if(np_states>1){
    for(is=1;is <= nstate;is++){
     a_tmp[is] = a[is];
     b_tmp[is] = b[is];
     c_tmp[is] = c[is];
    }/*endfor*/
    Allreduce(&(a_tmp[1]),&(a[1]),nstate,MPI_DOUBLE,MPI_SUM,0,comm); 
    Allreduce(&(b_tmp[1]),&(b[1]),nstate,MPI_DOUBLE,MPI_SUM,0,comm); 
    Allreduce(&(c_tmp[1]),&(c[1]),nstate,MPI_DOUBLE,MPI_SUM,0,comm); 
  }/*endif*/

/*------------------------------------------------------------------------*/
  /* iii) Occupation number conditions                 */

    for(is=1;is <= nstate; is++){a[is] -= occ[is];}

/*========================================================================*/
/* III) Solve quadratic equation for Lagrange multipliers                 */

  for(is=1;is<=nstate;is++) {
   descrim     = b[is]*b[is] - a[is]*c[is];
   if(descrim<-eps){
      printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
      printf("Norb-Shake fails on state processor %d \n",myid_state);
      printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
      fflush(stdout);exit(1);
   }else{
     descrim = fabs(descrim);
   }/*endif*/
   xl_norb[is] = (-b[is] + sqrt(descrim))/c[is];
  }/*endfor*/

/*========================================================================*/
/* IV) Apply constraint force to coefficient and coefficient velocities   */

  for(is=1;is<=nstate;is++) {
   i1 = 1+ioff[is];
   i2 = ncoef+ioff[is];
   for(icf=i1;icf<=i2;icf++) {
     creal[icf] += xl_norb[is]*creal_old[icf];
     cimag[icf] += xl_norb[is]*cimag_old[icf];
   }/* endfor i */
  }/* endfor is */
  
  rdt = 1.0/dt;
  for(is=1;is<=nstate;is++) {
   i1 = 1+ioff[is];
   i2 = ncoef+ioff[is];
   for(icf=i1;icf<=i2;icf++) {
     vcreal[icf] += rdt*xl_norb[is]*creal_old[icf];
     vcimag[icf] += rdt*xl_norb[is]*cimag_old[icf];
   }/* endfor i */
  }/* endfor is */

/*========================================================================*/
  /* V) Unscale old coefficients if necessary                             */

  if(icmass_unif != 1) {
   double cm1,scale;
   cm1 = cmass[ncoef_tot];
   for(is=1;is<=nstate;is++) {
     for(i=1;i<=ncoef;i++) {
       icoef = i+ioff[is];
       scale = cmass[(i+icmoff)]/cm1;       
       creal_old[icoef] *=scale;
       cimag_old[icoef] *=scale;
     }/* endfor i */
   }/* endfor is */
  }/* endif */

/*========================================================================*/
}/*end routine*/
/*========================================================================*/




/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void cp_rattle(double *creal,double *cimag,int icoef_form, int icoef_orth,
               double *vcreal,double *vcimag,int ivcoef_form,
               CPSCR_OVMAT *cpscr_ovmat,
               int *ioff,double *rocc_sum,
               CP_COMM_STATE_PKG *cp_comm_state_pkg,int cp_norb)

/*========================================================================*/
{/*begin routine*/
/*========================================================================*/
/*             Local variable declarations                                */

  int is,js;
  int nstate2;
  int iii;
  int ind1,ind2;
  double alpha,beta,g0_wght;

/*             Local pointer declarations                                */
   int np_states    = cp_comm_state_pkg->num_proc;
   int myid_state  = cp_comm_state_pkg->myid;
   int nstate      = cp_comm_state_pkg->nstate;
   double *a       = cpscr_ovmat->ovlap1;
   double *b       = cpscr_ovmat->ovlap2;
   double *ylamb   = cpscr_ovmat->ovlap3;
   double *a_scr   = cpscr_ovmat->ovlap4;

/*========================================================================*/
/* 0) Checks                                                              */
   if(cp_norb>0){
    if(icoef_orth !=0){
      printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
      printf("The coefficients must be in non-orthogonal form\n");
      printf("on state processor %d in cp_rattle \n",myid_state);
      printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
      fflush(stdout);exit(1);
     }/*endif*/
   }/*endif*/

   if(np_states>1){
    if((icoef_form+ivcoef_form) !=2){
      printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
      printf("The coefficients must be in transposed form\n");
      printf("on state processor %d in cp_rattle \n",myid_state);
      printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
      fflush(stdout);exit(1);
    }/*endif*/
   }/*endif*/

/*========================================================================*/
/* I) Calculate overlap matrices between coeffs and coeff velociites      */
/*------------------------------------------------------------------------*/
   /* i) Calculate some useful constants                            */
    
   nstate2 = nstate*nstate;

/*------------------------------------------------------------------------*/
   /* ii) Zero matrices                                             */

   for(is=1;is <= nstate2;is++){a[is] = ylamb[is] = 0.0;}

/*------------------------------------------------------------------------*/
   /* iii) Get overlap matrix                                             */

  alpha=2.0;
  cp_ovlap_mat_diff(creal,cimag,icoef_form,vcreal,vcimag,ivcoef_form,
                    a,a_scr,alpha,ioff,cp_comm_state_pkg);

/*========================================================================*/
/* II) Calculate Lagrange multiplier matrix             */

  /*--------------------------------*/
  /* create matrix transpose */

   for(is=1;is <=nstate; is++){
     for(js=1;js <=nstate; js++){
       ind1  = (is-1)*nstate + js;
       ind2 =  (js-1)*nstate + is;
       b[ind1] = a[ind2];
     }/* endfor js*/
   }/* endfor is */

   for(is=1;is <= nstate2;is++){
     ylamb[is] = -rocc_sum[is]*(a[is] + b[is]);
   }

/*========================================================================*/
/* III) Apply constraint force to coefficient velocities                  */

  alpha = 1.0;  beta = 1.0; g0_wght = 1.0;
  cp_rotation_prim(creal,icoef_form,vcreal,ivcoef_form,ylamb,
                   alpha,beta,g0_wght,ioff,cp_comm_state_pkg);
  cp_rotation_prim(cimag,icoef_form,vcimag,ivcoef_form,ylamb,
                   alpha,beta,g0_wght,ioff,cp_comm_state_pkg);


/*========================================================================*/
}/*end routine*/
/*========================================================================*/





/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void cp_rattle_mass(double *creal,double *cimag,int icoef_form, int icoef_orth,
                    double *vcreal,double *vcimag,int ivcoef_form,
                    double *cpreal,double *cpimag,
                    double c_tolrattle, double *cmass,
                    CPSCR_OVMAT *cpscr_ovmat, 
                    int *ioff,double *occ, double *rocc_sum,
                    int *iter_rattle_cp,
                    CP_COMM_STATE_PKG *cp_comm_state_pkg,int cp_norb)

/*========================================================================*/
{/*begin routine*/
/*========================================================================*/
/*             Local variable declarations                                */

   int i,is,js,nstate2,iproc;
   int icoef,ind,ind1,ind2,iii;
   int ipcoef_form;
   double alpha,beta,g0_wght;
   double cm1,scale;
   double dymax,dyl;

/*             Local pointer declarations                                */

   int np_states      = cp_comm_state_pkg->num_proc;
   int myid_state     = cp_comm_state_pkg->myid;
   int nstate         = cp_comm_state_pkg->nstate;
   int nstate_proc    = cp_comm_state_pkg->nstate_proc;
   int ncoef_tot      = cp_comm_state_pkg->ncoef;
   int ncoef          = cp_comm_state_pkg->nstate_ncoef_proc;
   int icmoff         = cp_comm_state_pkg->icoef_start-1;

   double *p          = cpscr_ovmat->ovlap1;
   double *q          = cpscr_ovmat->ovlap2;
   double *p1         = cpscr_ovmat->ovlap3;
   double *p2         = cpscr_ovmat->ovlap4;
   double *q_scr      = cpscr_ovmat->ovlap4;  /*don't panic */
   double *ylamb      = cpscr_ovmat->ovlap5;
   double *yold       = cpscr_ovmat->ovlap6;
   double *qt         = cpscr_ovmat->ovlap7;

/*========================================================================*/
/* 0) Checks                                                              */

   if(cp_norb>0){
     if(icoef_orth !=0){
       printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
       printf("The coefficients must be in non-orthogonal form\n");
       printf("on state processor %d in cp_rattle_mass \n",myid_state);
       printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
       fflush(stdout);exit(1);
     }/*endif*/
   }/*endif*/

   if(np_states>1){
    if((icoef_form+ivcoef_form) !=2){
      printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
      printf("The coefficients must be in transposed form\n");
      printf("on state processor %d in cp_rattle_mass \n",myid_state);
      printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
      fflush(stdout);exit(1);
    }/*endif*/
   }/*endif*/

/*========================================================================*/
/* I) Calculate overlap matrices between coeffs and coeff velocities      */
   /* i) Calculate some useful constants                            */

  nstate2 = nstate*nstate;

/*------------------------------------------------------------------------*/
   /* ii) Zero matrices                                             */

   for(is=1;is <= nstate2;is++){
      p[is] = q[is] = ylamb[is] = 0.0;
   }/* endfor */

/*------------------------------------------------------------------------*/
   /* iii) Calculate overlap integrals and construct matrices in equation  */

  alpha=2.0;
  cp_ovlap_mat_diff(creal,cimag,icoef_form,vcreal,vcimag,ivcoef_form,
                    q,q_scr,alpha,ioff,cp_comm_state_pkg);

/*------------------------------------------------------------------------*/
   /* iv) Scale coefficients by cmass[ncoef]/cmass[i]   */

   ipcoef_form = 1;
   cm1 = cmass[ncoef_tot];
   for(is=1;is<=nstate;is++) {
     for(i=1;i<=ncoef;i++) {
       icoef = i+ioff[is];
       scale = cm1/cmass[(i+icmoff)];       
       cpreal[icoef] = creal[icoef]*scale;
       cpimag[icoef] = cimag[icoef]*scale;
     }/* endfor i */
   }/* endfor is */

   
/*------------------------------------------------------------------------*/
   /* v) Calculate overlap integrals   */

  alpha=2.0;
  cp_ovlap_mat_diff(cpreal,cpimag,ipcoef_form,creal,cimag,icoef_form,
                    p,q_scr,alpha,ioff,cp_comm_state_pkg);

/*------------------------------------------------------------------------*/
   /* vi) Assign diagonal elements according matrix equation   */

     for(is=1;is <= nstate;is++){
         ind = (is-1)*nstate + is;
         p[ind] -= occ[is];
     }/* endfor */

  /*--------------------------------*/
  /* matrix transpose */

   for(is=1;is <=nstate; is++){
     for(js=1;js <=nstate; js++){
       ind1  = (is-1)*nstate + js;
       ind2 =  (js-1)*nstate + is;
       qt[ind1] = q[ind2];
     }/* endfor js*/
   }/* endfor is */

/*========================================================================*/
/* II) Iterative solution of matrix equation                              */

  /* i) Initialization of Lagrange multiplier matrix    */

  for(is=1;is <=nstate2; is++){
     yold[is] = ylamb[is] = -rocc_sum[is]*(q[is] + qt[is]);
  }/* endfor */

/*------------------------------------------------------------------------*/
  /* ii) Do iterative solution of matrix equation     */

  cp_iter_mat_solve_rattle_mass(ylamb,yold,p,q,qt,p1,p2,rocc_sum,
                                &nstate,iter_rattle_cp,c_tolrattle);

/*========================================================================*/
/* III) Apply constraint force to coefficient velocities                  */

  alpha = 1.0;  beta = 1.0; g0_wght = 1.0;
  cp_rotation_prim(cpreal,ipcoef_form,vcreal,ivcoef_form,ylamb,
                   alpha,beta,g0_wght,ioff,cp_comm_state_pkg);
  cp_rotation_prim(cpimag,ipcoef_form,vcimag,ivcoef_form,ylamb,
                   alpha,beta,g0_wght,ioff,cp_comm_state_pkg);

/*========================================================================*/
    }/*end routine*/
/*========================================================================*/





/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void cp_rattle_norb(double *creal,double *cimag,int icoef_form, int icoef_orth,
                    double *vcreal,double *vcimag,int ivcoef_form,
                    double *cpreal,double *cpimag,
                    int *ioff,double *occ, 
                    double *cmass,int icmass_unif,
                    CPSCR_OVMAT *cpscr_ovmat,
                    CP_COMM_STATE_PKG *cp_comm_state_pkg)

/*========================================================================*/
{/*begin routine*/
/*========================================================================*/
/*             Local variable declarations                                */
#include "../typ_defs/typ_mask.h"

  int iii,is,i,icoef,icf,i1,i2;
  double wght;
  double cm1,scale;

  int np_states       = cp_comm_state_pkg->num_proc;
  int myid_state      = cp_comm_state_pkg->myid;
  int nstate          = cp_comm_state_pkg->nstate;
  int ncoef           = cp_comm_state_pkg->nstate_ncoef_proc;
  int ncoef_tot       = cp_comm_state_pkg->ncoef;
  int icmoff          = cp_comm_state_pkg->icoef_start-1;
  double *xl_norb     = cpscr_ovmat->state_vec1;
  double *xl_norb_tmp = cpscr_ovmat->state_vec2;

  double *xn_norb     = cpscr_ovmat->state_vec3;
  double *xn_norb_tmp = cpscr_ovmat->state_vec4;
  MPI_Comm comm   = cp_comm_state_pkg->comm;

/*========================================================================*/
/* 0) Checks                                                              */

   if(icoef_orth !=0){
      printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
      printf("The coefficients must be in non-orthogonal form\n");
      printf("on state processor %d in cp_shake_norb \n",myid_state);
      printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
      fflush(stdout);exit(1);
   }/*endif*/

   if(np_states>1){
    if((icoef_form+ivcoef_form) !=2){
      printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
      printf("The coefficients must be in transposed form\n");
      printf("on state processor %d in cp_shake_norb \n",myid_state);
      printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
      fflush(stdout);exit(1);
    }/*endif*/
   }/*endif*/

/*========================================================================*/
/* II) Zero the Lagrange multiplier vector                                */
 
  for(is=1;is <= nstate;is++){
    xl_norb[is] = 0.0;
  }/*endfor*/

/*========================================================================*/
  /* II.V) Scale old coefficients if necessary                             */

  if(icmass_unif != 1) {
   cm1 = cmass[ncoef_tot];
   for(is=1;is<=nstate;is++) {
     xn_norb[is] = 0.0;
     for(i=1;i<=ncoef;i++) {
       icoef = i+ioff[is];
       scale = cm1/cmass[(i+icmoff)];       
       cpreal[icoef] = creal[icoef]*scale;
       cpimag[icoef] = cimag[icoef]*scale;
     }/* endfor i */
   }/* endfor is */
  } else {
   for(is=1;is<=nstate;is++) {
     xn_norb[is] = occ[is];
     for(i=1;i<=ncoef;i++) {
       icoef = i+ioff[is];
       cpreal[icoef] = creal[icoef];
       cpimag[icoef] = cimag[icoef];
     }/* endfor i */
   }/* endfor is */
  }/* endif */

/*========================================================================*/
/* III) Compute diagonal overlap of coefficients with coefficient velociites*/

  wght = 2.0;
  if(myid_state+1==np_states){wght=1.0;}
  for(is=1;is <= nstate;is++){
    i1 = 1+ioff[is];
    i2 = ncoef-1+ioff[is];
    for(icf=i1;icf <= i2; icf++){
      xl_norb[is] += creal[icf]*vcreal[icf] + cimag[icf]*vcimag[icf];
    }/*endfor*/
    if(icmass_unif != 1) {
     for(icf=i1;icf <= i2; icf++){
       xn_norb[is] += cpreal[icf]*creal[icf] + cpimag[icf]*cimag[icf];
     }/*endfor*/
    }/*endif*/
    icf = ncoef + ioff[is];
    xl_norb[is] = 2.0*xl_norb[is] + wght*(cpreal[icf]*vcreal[icf]
                                         +cpimag[icf]*vcimag[icf]);
    if(icmass_unif != 1) {
     xn_norb[is] = 2.0*xn_norb[is] + wght*(cpreal[icf]*creal[icf]
                                          +cpimag[icf]*cimag[icf]);
    }/*endif*/
  }/*endfor*/

  if(np_states>1){
    for(is=1;is <= nstate;is++){xl_norb_tmp[is] = xl_norb[is];}
    Allreduce(&(xl_norb_tmp[1]),&(xl_norb[1]),nstate,MPI_DOUBLE,MPI_SUM,0,
              comm); 
   if(icmass_unif != 1) {
    for(is=1;is <= nstate;is++){xn_norb_tmp[is] = xn_norb[is];}
    Allreduce(&(xn_norb_tmp[1]),&(xn_norb[1]),nstate,MPI_DOUBLE,MPI_SUM,0,
              comm); 
   }/*endfor*/
  }/*endif*/

/*========================================================================*/
/* IV) Occupation number scaling                                          */
  
    for(is=1;is <= nstate; is++){xl_norb[is] /=  xn_norb[is];}

/*========================================================================*/
/* V) Apply constraint force on coefficient velociites                    */

  for(is=1;is <= nstate;is++){
    i1 = 1+ioff[is];
    i2 = ncoef+ioff[is];
    for(icf=i1;icf<= i2; icf++){
      vcreal[icf] -= xl_norb[is]*cpreal[icf];
      vcimag[icf] -= xl_norb[is]*cpimag[icf];
    }/*endfor*/
  }/*endfor*/

/*========================================================================*/
    }/*end routine*/
/*========================================================================*/








