/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
/*                                                                          */
/*                         PI_MD:                                           */
/*             The future of simulation technology                          */
/*             ------------------------------------                         */
/*                     Module: cp_con_utils.c                               */
/*                                                                          */
/* Contains utilities for wave function S matrix evals and constraint solvers*/
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

#define HAND_MULT_OFF
#define MAX_ITER 1000


/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void cp_ovlap_mat_same(double *cr,double *ci, int icoef_form,
                       double *omat,double *omat_tmp,double alpha_in,int *ioff,
                       CP_COMM_STATE_PKG *cp_comm_state_pkg)

/*========================================================================*/
{/*begin routine*/
/*========================================================================*/
/*             Local variable declarations                                */

#include "../typ_defs/typ_mask.h"

  int ncoefm1,iii;
  int nstate2;
  int is,js;
  int ind_ij,ind_nc_i,ind_nc_j;
  int np_states = cp_comm_state_pkg->num_proc;
  int myid_state = cp_comm_state_pkg->myid;
  int nstate = cp_comm_state_pkg->nstate;
  int ncoef = cp_comm_state_pkg->ncoef;
  int itransp = 0;
  int inormal = 1;
  double beta;
  double alpha=alpha_in;
  MPI_Comm comm_states = cp_comm_state_pkg->comm;
  MPI_Comm world       = cp_comm_state_pkg->world;

/*========================================================================*/
/*========================================================================*/
/* I) Serial matrix multiplies                                            */

 if(np_states == 1){

/*------------------------------------------------------------------------*/
   /* A) Zero matrix                                    */

   nstate2=nstate*nstate;
   for(is=1;is <= nstate2;is++){
      omat[is] = 0.0;
   }/* endfor */

/*------------------------------------------------------------------------*/
   /* B) Sum over half of g-space                              */

   beta =   1.0;
   ncoefm1 = ncoef-1;


   GEN_MATMUL(&(cr[1]),&ncoef,&itransp,&(cr[1]),&ncoef,&inormal,
              &(omat[1]),&nstate,&nstate,&nstate,&ncoefm1,
              &alpha,&beta);
   GEN_MATMUL(&(ci[1]),&ncoef,&itransp,&(ci[1]),&ncoef,&inormal,
              &(omat[1]),&nstate,&nstate,&nstate,&ncoefm1,
              &alpha,&beta);

/*------------------------------------------------------------------------*/
   /* C) g=0 contribution                              */

  for(is=1;is <=nstate; is++){
    for(js=1;js <=nstate; js++){
      ind_ij  = (is-1)*nstate + js;
      ind_nc_i = ncoef + ioff[is];
      ind_nc_j = ncoef + ioff[js];
      omat[ind_ij] += cr[ind_nc_j]*cr[ind_nc_i];
    }/* endfor js*/
  }/* endfor is */

 }/*endif*/

/*========================================================================*/
/*========================================================================*/
/* II) Parallel matrix multiplies with g=0 correction                    */

 if(np_states != 1){

   if(icoef_form !=1){
      printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
      printf("The coefficients must be in transposed form\n");
      printf("on state processor %d in cp_ovlmat_mat_same \n",myid_state);
      printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
      fflush(stdout);exit(1);
   }/*endif*/

  nstate2=nstate*nstate;
  for(is=1;is <= nstate2;is++){omat_tmp[is] = 0.0;}

  cp_par_ovlap_same(cr,icoef_form,omat_tmp,alpha,cp_comm_state_pkg);
  cp_par_ovlap_same(ci,icoef_form,omat_tmp,alpha,cp_comm_state_pkg);

  Allreduce(&(omat_tmp[1]),&(omat[1]),nstate2,MPI_DOUBLE,MPI_SUM,0,
                                                      comm_states); 

 }/* endif */


/*========================================================================*/
    }/*end routine*/
/*========================================================================*/



/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void cp_par_ovlap_same(double *c,int icoef_form,
                       double *omat,double alpha,
                       CP_COMM_STATE_PKG *cp_comm_state_pkg)

/*========================================================================*/
{/*begin routine*/
/*========================================================================*/
/*             Local variable declarations                                */

#include "../typ_defs/typ_mask.h"

 double wght_use,alpha_loc,beta;
 int ioff,joff,is,js,ig,iii;
 int ind,jnd;
 int  itransp = 0;
 int  inormal = 1;

 int nstate                = cp_comm_state_pkg->nstate;
 int nstate_max            = cp_comm_state_pkg->nstate_max;
 int ncoef                 = cp_comm_state_pkg->ncoef;
 int nstate_proc           = cp_comm_state_pkg->nstate_proc;
 int nstate_proc_max       = cp_comm_state_pkg->nstate_proc_max;
 int nstate_ncoef_proc_max = cp_comm_state_pkg->nstate_ncoef_proc_max;
 int nstate_ncoef_proc     = cp_comm_state_pkg->nstate_ncoef_proc;
 int num_proc              = cp_comm_state_pkg->num_proc;
 int myid                  = cp_comm_state_pkg->myid;
 MPI_Comm comm  = cp_comm_state_pkg->comm;
 MPI_Comm world = cp_comm_state_pkg->world;

/*========================================================================*/
/* I) Check the forms */

    if(icoef_form !=1){
      printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
      printf("The coefficients must be in transposed form\n");
      printf("on state processor %d in cp_par_ovlap_mat_same \n",myid);
      printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
      fflush(stdout);exit(1);
    }/*endif*/

/*========================================================================*/
/* II) Matrix multiply with correct treatment for g=0                     */


  wght_use = (myid==num_proc-1 ? 1.0 : alpha);
#ifdef HAND_MULT
  ioff = 0;
  for(is=1;is<=nstate;is++){
   joff = 0;
   for(js=1;js<=is;js++){
     ind  = (is-1)*nstate + js;
     for(ig=1;ig<=nstate_ncoef_proc-1;ig++){
       omat[ind] += alpha*(c[(joff+ig)]*c[(ioff+ig)]);
     }/*endfor*/
     ig = nstate_ncoef_proc;
     omat[ind] += wght_use*(c[(joff+ig)]*c[(ioff+ig)]);
     jnd = (js-1)*nstate + is;
     omat[jnd] = omat[ind];
     joff += nstate_ncoef_proc_max;
   }/*endfor*/
   ioff += nstate_ncoef_proc_max;
  }/*endfor*/
#else
  alpha_loc = alpha;
  beta      = 1.0;
  GEN_MATMUL(&(c[1]),&nstate_ncoef_proc_max,&itransp,&(c[1]),
             &nstate_ncoef_proc_max,&inormal,
             &(omat[1]),&nstate,&nstate,&nstate,&nstate_ncoef_proc_max,
             &alpha_loc,&beta);
    
  if(wght_use!=alpha){
    ioff = 0;
    for(is=1;is<=nstate;is++){
      joff = 0;
      for(js=1;js<=is;js++){
        ind  = (is-1)*nstate + js;
        ig = nstate_ncoef_proc;
        omat[ind] += (wght_use-alpha)*(c[(joff+ig)]*c[(ioff+ig)]);
        jnd = (js-1)*nstate + is;
        omat[jnd] = omat[ind];
        joff += nstate_ncoef_proc_max;
      }/*endfor*/
      ioff += nstate_ncoef_proc_max;
    }/*endfor : is*/
  }/*endif:wght_use*/ 
#endif

/*========================================================================*/
 }/* end routine */
/*========================================================================*/




/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void cp_ovlap_mat_diff(double *cr1,double *ci1,int i1coef_form,
                       double *cr2,double *ci2,int i2coef_form,
                       double *omat,double *omat_tmp,double alpha_in,int *ioff,
                       CP_COMM_STATE_PKG *cp_comm_state_pkg)

/*========================================================================*/
{/*begin routine*/
/*========================================================================*/
/*             Local variable declarations                                */

#include "../typ_defs/typ_mask.h"

  int ncoefm1,iproc,iii;
  int nstate2;
  int is,js;
  int ind_ij,ind_nc_i,ind_nc_j;
  int np_states = cp_comm_state_pkg->num_proc;
  int myid_state = cp_comm_state_pkg->myid;
  int nstate = cp_comm_state_pkg->nstate;
  int ncoef = cp_comm_state_pkg->ncoef;
  int itransp = 0;
  int inormal = 1;
  double beta;
  double alpha=alpha_in;
  MPI_Comm comm_states = cp_comm_state_pkg->comm;
  MPI_Comm world       = cp_comm_state_pkg->world;

/*========================================================================*/
/*========================================================================*/
/* I) Serial matrix multiplies                                            */

 if(np_states == 1){

/*------------------------------------------------------------------------*/
   /* A) Zero matrix                                    */

   nstate2=nstate*nstate;
   for(is=1;is <= nstate2;is++){
      omat[is] = 0.0;
   }/* endfor */

/*------------------------------------------------------------------------*/
   /* B) Sum over half of g-space                              */

   beta =   1.0;
   ncoefm1 = ncoef-1;
   GEN_MATMUL(&(cr1[1]),&ncoef,&itransp,&(cr2[1]),&ncoef,&inormal,
              &(omat[1]),&nstate,&nstate,&nstate,&ncoefm1,
              &alpha,&beta);
   GEN_MATMUL(&(ci1[1]),&ncoef,&itransp,&(ci2[1]),&ncoef,&inormal,
              &(omat[1]),&nstate,&nstate,&nstate,&ncoefm1,
              &alpha,&beta);

/*------------------------------------------------------------------------*/
   /* C) g=0 contribution                              */

  for(is=1;is <=nstate; is++){
    for(js=1;js <=nstate; js++){
      ind_ij  = (is-1)*nstate + js;
      ind_nc_i = ncoef + ioff[is];
      ind_nc_j = ncoef + ioff[js];
      omat[ind_ij] += cr1[ind_nc_j]*cr2[ind_nc_i];
    }/* endfor js*/
  }/* endfor is */

 }/*endif*/

/*========================================================================*/
/*========================================================================*/
/* II) Parallel matrix multiplies wirh g=0 correction                    */

 if(np_states != 1){

   if((i1coef_form+i2coef_form) !=2){
      printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
      printf("The coefficients must be in transposed form\n");
      printf("on state processor %d in cp_ovlmat_mat_diff \n",myid_state);
      printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
      fflush(stdout);exit(1);
   }/*endif*/

  nstate2=nstate*nstate;
  for(is=1;is <= nstate2;is++){omat_tmp[is] = 0.0;}

  cp_par_ovlap_diff(cr1,i1coef_form,cr2,i2coef_form,omat_tmp,
                    alpha,cp_comm_state_pkg);
  cp_par_ovlap_diff(ci1,i1coef_form,ci2,i2coef_form,omat_tmp,
                    alpha,cp_comm_state_pkg);
  Allreduce(&(omat_tmp[1]),&(omat[1]),nstate2,MPI_DOUBLE,MPI_SUM,0,
                                                       comm_states); 

 }/* endif */


/*========================================================================*/
}/*end routine*/
/*========================================================================*/




/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void cp_par_ovlap_diff(double *c1,int i1coef_form,
                       double *c2,int i2coef_form,
                       double *omat,double alpha,
                       CP_COMM_STATE_PKG *cp_comm_state_pkg)

/*========================================================================*/
{/*begin routine*/
/*========================================================================*/
/*             Local variable declarations                                */

#include "../typ_defs/typ_mask.h"

 double wght_use,alpha_loc,beta;
 int ioff,joff,is,js,ig,iproc,iii;
 int ind,jnd;
 int itransp = 0;
 int inormal = 1;

 int nstate                = cp_comm_state_pkg->nstate;
 int nstate_max            = cp_comm_state_pkg->nstate_max;
 int ncoef                 = cp_comm_state_pkg->ncoef;
 int nstate_proc           = cp_comm_state_pkg->nstate_proc;
 int nstate_proc_max       = cp_comm_state_pkg->nstate_proc_max;
 int nstate_ncoef_proc_max = cp_comm_state_pkg->nstate_ncoef_proc_max;
 int nstate_ncoef_proc     = cp_comm_state_pkg->nstate_ncoef_proc;
 int num_proc              = cp_comm_state_pkg->num_proc;
 int myid                  = cp_comm_state_pkg->myid;
 MPI_Comm comm   = cp_comm_state_pkg->comm;
 MPI_Comm world  = cp_comm_state_pkg->world;

/*========================================================================*/
/* I) Check the forms */

    if((i1coef_form+i2coef_form) !=2){
      printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
      printf("The coefficients must be in transposed form\n");
      printf("on state processor %d in cp_par_ovlap_mat_diff \n",myid);
      printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
      fflush(stdout);exit(1);
    }/*endif*/

/*========================================================================*/
/* II) Matrix multiply with correct treatment for g=0                     */

  wght_use = (myid==num_proc-1 ? 1.0 : alpha);

#ifdef HAND_MULT
  ioff = 0;
  for(is=1;is<=nstate;is++){
   joff = 0;
   for(js=1;js<=nstate;js++){
     ind  = (is-1)*nstate + js;
     for(ig=1;ig<=nstate_ncoef_proc-1;ig++){
       omat[ind] += alpha*(c1[(joff+ig)]*c2[(ioff+ig)]);
     }/*endfor*/
     ig = nstate_ncoef_proc;
     omat[ind] += wght_use*(c1[(joff+ig)]*c2[(ioff+ig)]);
     joff += nstate_ncoef_proc_max;
   }/*endfor*/
   ioff += nstate_ncoef_proc_max;
  }/*endfor*/
#else

  alpha_loc = alpha;
  beta      = 1.0;
  GEN_MATMUL(&(c1[1]),&nstate_ncoef_proc_max,&itransp,&(c2[1]),
              &nstate_ncoef_proc_max,&inormal,
              &(omat[1]),&nstate,&nstate,&nstate,&nstate_ncoef_proc_max,
              &alpha_loc,&beta);

  if(alpha!=wght_use){
    ioff = 0;
    for(is=1;is<=nstate;is++){
      joff = 0;
      for(js=1;js<=nstate;js++){
        ind  = (is-1)*nstate + js;
        ig = nstate_ncoef_proc;
        omat[ind] += (wght_use-alpha)*(c1[(joff+ig)]*c2[(ioff+ig)]);
        joff += nstate_ncoef_proc_max;
      }/*endfor*/
      ioff += nstate_ncoef_proc_max;
    }/*endfor*/
  }/*endif : alpha*/
#endif

/*========================================================================*/
    }/* end routine */
/*========================================================================*/



/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void cp_iter_mat_solve_shake(double *x,double *x2,double *xold,
                             double *a,double *b,double *p1,double *p2,
                             double *occ,double *rocc_sum,
                             int *nstate,int *iter_shake_cp,
                             double c_tolshake)

/*========================================================================*/
{/*begin routine*/
/*========================================================================*/
/*             Local variable declarations                                */

  double alpha,beta;
  double sqrt_occ;
  int ind1,ind2,is,js;
  int nstate2;
  double xldiff,dxl,dxmax;
  int itransp = 0;
  int inormal = 1;

/*========================================================================*/
  nstate2 = (*nstate)*(*nstate);
  alpha = 1.0;
  beta = 0.0;
  *iter_shake_cp = 0;
  do {
    (*iter_shake_cp)++;

    if(*iter_shake_cp > MAX_ITER){
      printf("Maximum CP SHAKE iterations exceeded\n");
      printf("with # iterations = %d\n",*iter_shake_cp);
      printf("and tolerance = %.12g\n",dxmax);
      printf("Consider increasing your tolerance if dxmax is small\n");
      exit(1);
    }/* endif */

    GEN_MATMUL(&(x[1]),nstate,&inormal,&(b[1]),nstate,&itransp,
               &(p1[1]),nstate,nstate,nstate,nstate,&alpha,&beta);

    for(is=1;is <=*nstate; is++){
      for(js=1;js <=*nstate; js++){
        sqrt_occ = sqrt(occ[js]);
        ind1  = (is-1)*(*nstate) + js;
        x[ind1] *= sqrt_occ;
      }/* endfor js*/
    }/* endfor is */

    GEN_MATMUL(&(x[1]),nstate,&itransp,&(x[1]),nstate,&inormal,
               &(x2[1]),nstate,nstate,nstate,nstate,&alpha,&beta);

  /*--------------------------------*/
  /* matrix transpose */

   for(is=1;is <=*nstate; is++){
     for(js=1;js <=*nstate; js++){
       ind1  = (is-1)*(*nstate) + js;
       ind2 =  (js-1)*(*nstate) + is;
       p2[ind1] = p1[ind2];
     }/* endfor js*/
   }/* endfor is */

   dxmax = 0.0;
   for(is=1;is <= nstate2;is++){
      x[is] = rocc_sum[is]*(a[is] + p1[is] 
                + p2[is] - x2[is]);
      xldiff = x[is]-xold[is];
      dxl = xldiff*xldiff;
      dxmax = MAX(dxl,dxmax);
      xold[is] = x[is];
   }/* endfor is */
    dxmax = sqrt(dxmax);
  } while(dxmax > c_tolshake);

/*========================================================================*/
}/*end routine*/
/*========================================================================*/





/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void cp_iter_mat_solve_shake_mass(double *x,double *x2,double *xold,
                                  double *a,double *b,double *c,
                                  double *p1,double *p2,double *rocc_sum,
                                  int *nstate,int *iter_shake_cp,
                                  double c_tolshake)

/*========================================================================*/
{/*begin routine*/
/*========================================================================*/
/*             Local variable declarations                                */

  double alpha,beta;
  int ind1,ind2,is,js;
  int nstate2;
  double xldiff,dxl,dxmax;
  int itransp = 0;
  int inormal = 1;

/*========================================================================*/
  nstate2 = (*nstate)*(*nstate);
  alpha = 1.0;
  beta = 0.0;
  *iter_shake_cp = 0;
  for(is=1;is <= nstate2; is++){
   p1[is]=p2[is]=0.0;
  }
  do {
    (*iter_shake_cp)++;

    if(*iter_shake_cp > MAX_ITER){
      printf("Maximum CP SHAKE iterations exceeded\n");
      printf("with # iterations = %d\n",*iter_shake_cp);
      printf("and tolerance = %.12g\n",dxmax);
      printf("Consider increasing your tolerance if dxmax is small\n");
      exit(1);
    }/* endif */

    GEN_MATMUL(&(x[1]),nstate,&inormal,&(b[1]),nstate,&itransp,
               &(p1[1]),nstate,nstate,nstate,nstate,&alpha,&beta);

    GEN_MATMUL(&(x[1]),nstate,&inormal,&(c[1]),nstate,&inormal,
               &(p2[1]),nstate,nstate,nstate,nstate,&alpha,&beta);

    GEN_MATMUL(&(p2[1]),nstate,&inormal,&(x[1]),nstate,&inormal,
               &(x2[1]),nstate,nstate,nstate,nstate,&alpha,&beta);

  /*--------------------------------*/
  /* matrix transpose */

   for(is=1;is <=*nstate; is++){
     for(js=1;js <=*nstate; js++){
       ind1  = (is-1)*(*nstate) + js;
       ind2 =  (js-1)*(*nstate) + is;
       p2[ind1] = p1[ind2];
     }/* endfor js*/
   }/* endfor is */

   dxmax = 0.0;
   for(is=1;is <= nstate2;is++){
      x[is] = rocc_sum[is]*(a[is] + p1[is] 
                + p2[is] - x2[is]);
      xldiff = x[is]-xold[is];
      dxl = xldiff*xldiff;
      dxmax = MAX(dxl,dxmax);
      xold[is] = x[is];
   }/* endfor is */
    dxmax = sqrt(dxmax);
  } while(dxmax > c_tolshake);


/*========================================================================*/
}/*end routine*/
/*========================================================================*/





/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void cp_iter_mat_solve_rattle_mass(double *y,double *yold,
                                   double *p,double *q,double *qt,
                                   double *p1,double *p2,double *rocc_sum,
                                   int *nstate,int *iter_rattle_cp,
                                   double c_tolrattle)

/*========================================================================*/
{/*begin routine*/
/*========================================================================*/
/*             Local variable declarations                                */

  double alpha,beta;
  int ind1,ind2,is,js;
  int nstate2;
  double yldiff,dyl,dymax;
  int itransp = 0;
  int inormal = 1;

/*========================================================================*/
  nstate2 = (*nstate)*(*nstate);
  alpha = 1.0;
  beta = 0.0;
  *iter_rattle_cp = 0;
  do {
    (*iter_rattle_cp)++;

    if(*iter_rattle_cp > MAX_ITER){
      printf("Maximum CP SHAKE iterations exceeded\n");
      printf("with # iterations = %d\n",*iter_rattle_cp);
      printf("and tolerance = %.12g\n",dymax);
      printf("Consider increasing your tolerance if dymax is small\n");
      exit(1);
    }/* endif */

    GEN_MATMUL(&(p[1]),nstate,&inormal,&(y[1]),nstate,&itransp,
               &(p1[1]),nstate,nstate,nstate,nstate,&alpha,&beta);

  /*--------------------------------*/
  /* matrix transpose */

   for(is=1;is <=*nstate; is++){
     for(js=1;js <=*nstate; js++){
       ind1  = (is-1)*(*nstate) + js;
       ind2 =  (js-1)*(*nstate) + is;
       p2[ind1] = p1[ind2];
     }/* endfor js*/
   }/* endfor is */

   dymax = 0.0;
   for(is=1;is <= nstate2;is++){
      y[is] = -rocc_sum[is]*(p1[is] + p2[is] 
                + q[is] + qt[is]);
      yldiff = y[is]-yold[is];
      dyl = yldiff*yldiff;
      dymax = MAX(dyl,dymax);
      yold[is] = y[is];
   }/* endfor is */
    dymax = sqrt(dymax);
  } while(dymax > c_tolrattle);

/*========================================================================*/
}/*end routine*/
/*========================================================================*/






























