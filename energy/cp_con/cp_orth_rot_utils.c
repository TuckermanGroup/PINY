/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
/*                                                                          */
/*                         PI_MD:                                           */
/*             The future of simulation technology                          */
/*             ------------------------------------                         */
/*                     Module: cp_rotbasis_utils.c                          */
/*                                                                          */
/* Contains utilities rotating into and out of various basises              */
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
#define HAND_ROT_OFF



/*==========================================================================*/
/* Control gram_schmidt orthog */
/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void control_cp_gram_schmidt(double *creal,double *cimag,int icoef_form,
                             double *occ,double *omat,
                             double *omat_tmp,int *ioff,
                             CP_COMM_STATE_PKG *cp_comm_state_pkg)

/*========================================================================*/
{/*begin routine*/
/*========================================================================*/
/*             Local variable declarations                                */

#include "../typ_defs/typ_mask.h"


  int nstate   = cp_comm_state_pkg->nstate;
  int ncoef    = cp_comm_state_pkg->ncoef;
  int num_proc = cp_comm_state_pkg->num_proc;
  int myid_state= cp_comm_state_pkg->myid;

/*========================================================================*/
/* 0) Do nothing if nstate = 0                                            */

  if(nstate == 0) return;

/*========================================================================*/
/* I) Scalar  GS                                                          */

  if(num_proc ==1){
   cp_gram_schmidt_scalar(creal,cimag,occ,omat,ioff,nstate,ncoef);
  }/*endif*/

/*========================================================================*/
/* II) Parallel GS                                                        */

  if(num_proc > 1){
    if(icoef_form !=1){
      printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
      printf("The coefficients must be in transposed form\n");
      printf("on state processor %d in control_cp_gram_schmidt \n",myid_state);
      printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
      fflush(stdout);exit(1);
    }/*endif*/
   cp_gram_schmidt_par(creal,cimag,icoef_form,occ,ioff,omat,omat_tmp,
                       cp_comm_state_pkg);
  }/*endif*/

/*------------------------------------------------------------------------*/
}/*end routine*/
/*========================================================================*/



/*==========================================================================*/
/* Gram_schmidt orthog in scalar*/
/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void cp_gram_schmidt_scalar(double *cre,double *cim,double *occ,
                            double *ovlap,int *ioff,int nstate,int ncoef)

/*======================================================================*/
/*                Begin Routine */
{   /*begin routine */
/*======================================================================*/
/*               Local variable declarations                            */
  
  double cmix;
  double anorm,scale;
  int len,ig,is,js,joff;
  int i1,i2;
  int ind,jnd,nstate2;

/*==========================================================================*/
/* I) Normalize the states                                                 */

/*--------------------------------------------------------------------------*/
/* i) Find overlap of the each state with iself                            */

  for(is=1;is<=nstate;is++){
    ovlap[is]    = 0.0;
    i1 = ioff[is]+1; i2=ioff[is]+ncoef-1;
    for(ig=i1;ig<=i2;ig++){
      ovlap[is] += (cre[ig]*cre[ig]+cim[ig]*cim[ig]);
    }/*endfor*/
    ovlap[is] *=2.0;
    ig         =  ncoef+ioff[is];
    ovlap[is] += cre[ig]*cre[ig];
  }/*endfor:is*/

/*-------------------------------------------------------------------------*/
/* ii) Normalize each state to 1 for convenience                           */

  for(is=1;is<=nstate;is++){
    scale = sqrt(1.0/ovlap[is]);
    i1 = ioff[is]+1; i2=ioff[is]+ncoef;
    for(ig=i1;ig<=i2;ig++){
      cre[ig] *= scale;
      cim[ig] *= scale;
    }/*endfor:ig*/
  }/*endfor:is*/

/*========================================================================*/
/* III) Gram schmidt orthogonalize                                         */

  for(is=1;is<=nstate-1;is++){

/*-----------------------------------------------------------------------*/
/* i) Find overlap of the is^th state with all higher states             */
/*    Note the joff value                                                */

   len = nstate-is;
   for(js=1;js<=len;js++){
     joff = ioff[(is+js)];
     ovlap[js]     = 0.0;
     for(ig=1;ig<=ncoef-1;ig++){
      ovlap[js]     += (cre[(ig+ioff[is])]*cre[(ig+joff)]
                       +cim[(ig+ioff[is])]*cim[(ig+joff)]);
     }/*endfor:ig*/
     ovlap[js]     *=2.0;
     ig            =  ncoef;
     ovlap[js]     += cre[(ig+ioff[is])]*cre[(ig+joff)];
   }/*endfor:js*/

/*-----------------------------------------------------------*/
/* ii) Orthogonalize all higher states to the is^th state    */
/*     keeping each state properly normalized                */

   for(js=1;js<=len;js++){
     joff = ioff[(is+js)];
     cmix  = -ovlap[js];
     anorm = 1.0 + 2.0*cmix*ovlap[js] + cmix*cmix;
     scale = sqrt(1.0/anorm);
     for(ig=1;ig<=ncoef;ig++){
       cre[(ig+joff)] = (cre[(ig+ioff[is])]*cmix + cre[(ig+joff)])*scale;
       cim[(ig+joff)] = (cim[(ig+ioff[is])]*cmix + cim[(ig+joff)])*scale;
     }/*endfor:ig*/
   }/*endfor:js*/

  }/*endfor:is loop*/

/*========================================================================*/
/* III) Check                                                             */

#ifdef DEBUG_GS

  nstate2 = nstate*nstate;
  for(is=1;is<=nstate;is++){
/*-----------------------------------------------------------------------*/
/* i) Find overlap of the is^th state with all lower states              */
   for(js=1;js<=is;js++){
     ind = (is-1)*nstate + js;
     ovlap[ind]     = 0.0;
     for(ig=1;ig<=ncoef-1;ig++){
      ovlap[ind]     += (cre[(ig+ioff[is])]*cre[(ig+ioff[js])]
                        +cim[(ig+ioff[is])]*cim[(ig+ioff[js])]);
     }/*endfor*/
     ovlap[ind]     *=2.0;
     ig              =  ncoef;
     ovlap[ind]     += wght*(cre[(ig+ioff[is])]*cre[(ig+ioff[js])]
                            +cim[(ig+ioff[is])]*cim[(ig+ioff[js])]);
     jnd = (js-1)*nstate + is;
     ovlap[jnd] = ovlap[ind];
   }/*endfor*/
  }/*endfor:is loop*/
/*-----------------------------------------------------------------------*/
/* ii) Output                                                            */
  for(is=1;is<=nstate;is++){
    for(js=1;js<=is;js++){
      jnd = (js-1)*nstate + is;
      ind = (is-1)*nstate + js;
      printf("is=%d js=%d %g %g \n",is,js,ovlap[ind],ovlap[jnd]);
    }/*endfor*/
   scanf("%d",&iii);
  }/*endfor*/

#endif

/*========================================================================*/
/* III) Adjust norms to occupation numbers                                */

  for(is=1;is<=nstate;is++){
    scale = sqrt(occ[is]);
    i1 = ioff[is]+1; i2=ioff[is]+ncoef;
    for(ig=i1;ig<=i2;ig++){
      cre[ig] *= scale;
      cim[ig] *= scale;
    }/*endfor:ig*/
  }/*endfor:is*/

/*========================================================================*/
}/*end routine*/
/*========================================================================*/



/*==========================================================================*/
/* Gram_schmidt orthog in parallel */
/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void cp_gram_schmidt_par(double *cre,double *cim,int icoef_form,
                         double *occ,int *ioff_vec,
                         double *ovlap, double *ovlap_tmp,
                         CP_COMM_STATE_PKG *cp_comm_state_pkg)

/*======================================================================*/
/*                Begin Routine */
   {/*begin routine */
/*======================================================================*/
/*               Local variable declarations                            */
#include "../typ_defs/typ_mask.h"
  
  double wght;
  double cmix;
  double anorm,scale;
  int len,ig,is,ioff,js,joff;
  int iii,jnd,nstate2,ind,i1,i2;

  int nstate                = cp_comm_state_pkg->nstate;
  int ncoef                 = cp_comm_state_pkg->nstate_ncoef_proc;
  int nstate_ncoef_proc_max = cp_comm_state_pkg->nstate_ncoef_proc_max;
  int nstate_ncoef_proc     = cp_comm_state_pkg->nstate_ncoef_proc;
  int num_proc              = cp_comm_state_pkg->num_proc;
  int myid                  = cp_comm_state_pkg->myid;
  MPI_Comm comm             = cp_comm_state_pkg->comm;

/*========================================================================*/
/* 0) Checks                                                              */

    if(icoef_form !=1){
      printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
      printf("The coefficients must be in transposed form\n");
      printf("on state processor %d in cp_gram_schmidt_par \n",myid);
      printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
      fflush(stdout);exit(1);
    }/*endif*/

/*==========================================================================*/
/* I) Define useful constant                                                */

  wght = (myid != num_proc-1 ? 2.0 : 1.0);

/*==========================================================================*/
/* II) Normalize the states                                                 */

/*--------------------------------------------------------------------------*/
/* i) Find overlap of the each state with iself                            */

  ioff = 0;  
  for(is=1;is<=nstate;is++){
    ovlap_tmp[is]    = 0.0;
    i1 = 1+ioff;i2=ioff+nstate_ncoef_proc-1;
    for(ig=i1;ig<=i2;ig++){
      ovlap_tmp[is] += (cre[ig]*cre[ig]+cim[ig]*cim[ig]);
    }/*endfor:ig*/
    ovlap_tmp[is]   *=2.0;
    ig               =  nstate_ncoef_proc + ioff;
    ovlap_tmp[is]   += wght*(cre[ig]*cre[ig]+cim[ig]*cim[ig]);
    ioff            += nstate_ncoef_proc_max;
  }/*endfor:is*/
  Allreduce(&(ovlap_tmp[1]),&(ovlap[1]),nstate,MPI_DOUBLE,MPI_SUM,0,comm); 

/*-------------------------------------------------------------------------*/
/* ii) Normalize each state to 1                                           */

  ioff = 0;  
  for(is=1;is<=nstate;is++){
    scale = sqrt(1.0/ovlap[is]);
    i1 = 1+ioff;i2=ioff+nstate_ncoef_proc;
    for(ig=i1;ig<=i2;ig++){
      cre[ig] *= scale;
      cim[ig] *= scale;
    }/*endfor:ig*/
    ioff += nstate_ncoef_proc_max;
  }/*endfor:is*/

/*========================================================================*/
/* III) Gram schmidt orthogonalize                                         */

  ioff = 0;
  for(is=1;is<=nstate-1;is++){

/*-----------------------------------------------------------------------*/
/* i) Find overlap of the is^th state with all higher states             */
/*    Note the joff starting value                                       */

   len = nstate-is;
   joff = ioff + nstate_ncoef_proc_max;
   for(js=1;js<=len;js++){
     ovlap_tmp[js]     = 0.0;
     for(ig=1;ig<=nstate_ncoef_proc-1;ig++){
      ovlap_tmp[js]     += (cre[(ig+ioff)]*cre[(ig+joff)]
                           +cim[(ig+ioff)]*cim[(ig+joff)]);
     }/*endfor:ig*/
     ovlap_tmp[js]     *=2.0;
     ig                 =  nstate_ncoef_proc;
     ovlap_tmp[js]     += wght*(cre[(ig+ioff)]*cre[(ig+joff)]
                               +cim[(ig+ioff)]*cim[(ig+joff)]);
     joff += nstate_ncoef_proc_max;
   }/*endfor:js*/
   Allreduce(&(ovlap_tmp[1]),&(ovlap[1]),len,MPI_DOUBLE,MPI_SUM,0,comm); 

/*-----------------------------------------------------------*/
/* ii) Orthogonalize all higher states to the is^th state    */
/*     keeping each state properly normalized                */

   joff = ioff + nstate_ncoef_proc_max;
   for(js=1;js<=len;js++){
     cmix  = -ovlap[js];
     anorm = 1.0 + 2.0*cmix*ovlap[js] + cmix*cmix;
     scale = sqrt(1.0/anorm);
     for(ig=1;ig<=nstate_ncoef_proc;ig++){
       cre[(ig+joff)] = (cre[(ig+ioff)]*cmix+cre[(ig+joff)])*scale;
       cim[(ig+joff)] = (cim[(ig+ioff)]*cmix+cim[(ig+joff)])*scale;
     }/*endfor:ig*/
     joff += nstate_ncoef_proc_max;
   }/*endfor:js*/

/*-----------------------------------------------------------*/
/* iv) increment ioff, the is guys offset                    */

   ioff += nstate_ncoef_proc_max;
  }/*endfor:is*/

/*========================================================================*/
/* III) Check                                                             */
#ifdef DEBUG_GS

  nstate2 = nstate*nstate;
  ioff = 0;
  for(is=1;is<=nstate;is++){
/*-----------------------------------------------------------------------*/
/* i) Find overlap of the is^th state with all lower states              */
   joff = 0;
   for(js=1;js<=is;js++){
     ind = (is-1)*nstate + js;
     ovlap_tmp[ind]     = 0.0;
     for(ig=1;ig<=nstate_ncoef_proc-1;ig++){
      ovlap_tmp[ind]     += (cre[(ig+ioff)]*cre[(ig+joff)]
                            +cim[(ig+ioff)]*cim[(ig+joff)]);
     }/*endfor*/
     ovlap_tmp[ind]     *=2.0;
     ig                  =  nstate_ncoef_proc;
     ovlap_tmp[ind]     += wght*(cre[(ig+ioff)]*cre[(ig+joff)]
                                +cim[(ig+ioff)]*cim[(ig+joff)]);
     jnd = (js-1)*nstate + is;
     ovlap_tmp[jnd] = ovlap_tmp[ind];
     joff += nstate_ncoef_proc_max;
   }/*endfor*/
/*-----------------------------------------------------------*/
/* ii) increment is dudes offset                             */
   ioff += nstate_ncoef_proc_max;
  }/*endfor:is loop*/
/*-----------------------------------------------------------------------*/
/* iii) All reduce                                                       */
  Allreduce(&(ovlap_tmp[1]),&(ovlap[1]),nstate2,MPI_DOUBLE,MPI_SUM,0,comm); 
/*-----------------------------------------------------------------------*/
/* iv) Output                                                            */
  for(is=1;is<=nstate;is++){
   if(myid==0){
    for(js=1;js<=is;js++){
      jnd = (js-1)*nstate + is;
      ind = (is-1)*nstate + js;
      printf("is=%d js=%d %g %g \n",is,js,ovlap[ind],ovlap[jnd]);
    }/*endfor*/
   }/*endif*/
   scanf("%d",&iii);
   Dbx_Barrier(comm);
  }/*endfor*/

#endif
/*========================================================================*/
/* VI) Adjust the norms to the occupation numbers                         */

  for(is=1;is<=nstate;is++){
    scale = sqrt(occ[is]);
    i1 = ioff_vec[is]+1; i2=ioff_vec[is]+ncoef;
    for(ig=i1;ig<=i2;ig++){
      cre[ig] *= scale;
      cim[ig] *= scale;
    }/*endfor:ig*/
  }/*endfor:is*/

/*--------------------------------------------------------------------------*/
  }/*end routine */
/*==========================================================================*/



/*==========================================================================*/
/* Normalize a CP vector */
/*========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void cp_normalize(double *creal,double *cimag,int icoef_form,
                  double *occ,int *ioff,
                  double *norm_save,double *norm_tmp,
                  CP_COMM_STATE_PKG *cp_comm_state_pkg)

/*========================================================================*/
{/*begin routine*/
/*========================================================================*/
/*             Local variable declarations                                */
#include "../typ_defs/typ_mask.h"

  int iii,is,icf,i1,i2;
  double norm,scale,rsnorm,wght;
  int ncoef      = cp_comm_state_pkg->nstate_ncoef_proc;
  int nstate     = cp_comm_state_pkg->nstate;
  int myid_state = cp_comm_state_pkg->myid;
  int np_states  = cp_comm_state_pkg->num_proc;
  MPI_Comm comm   = cp_comm_state_pkg->comm;

/*========================================================================*/
/* 0) Checks                                                              */

    if(icoef_form !=1&&np_states>1){
      printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
      printf("The coefficients must be in transposed form\n");
      printf("on state processor %d in cp_normalize \n",myid_state);
      printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
      fflush(stdout);exit(1);
    }/*endif*/

/*========================================================================*/
/* I) Calculate the norms                                                 */

 wght = 2.0;
 if(myid_state+1==np_states){wght=1.0;}

 for(is=1;is<=nstate;is++){

   norm = 0.0;
   i1 = 1 + ioff[is];
   i2 = ncoef-1 + ioff[is];
   for(icf=i1;icf<=i2; icf++){
     norm += creal[icf]*creal[icf] + cimag[icf]*cimag[icf];    
   }
   icf = ncoef + ioff[is];
   norm_save[is] = 2.0*norm + wght*(creal[icf]*creal[icf]+
                                      cimag[icf]*cimag[icf]);

 }/* endfor is*/

 if(np_states>1){
  for(is=1;is <= nstate;is++){norm_tmp[is] = norm_save[is];}
  Allreduce(&(norm_tmp[1]),&(norm_save[1]),nstate,MPI_DOUBLE,MPI_SUM,0,comm); 
 }/*endif*/

/*========================================================================*/
/* I) Normalize each state                                                */

 for(is=1;is<=nstate;is++){

   rsnorm = sqrt(occ[is]/norm_save[is]);
   i1 = 1 + ioff[is];
   i2 = ncoef + ioff[is];
   for(icf=i1;icf<=i2; icf++){
     creal[icf] *= rsnorm;  
     cimag[icf] *= rsnorm;
   }/*endfor*/

 }/* endfor is*/

/*========================================================================*/
   }/*end routine*/
/*========================================================================*/

/*==========================================================================*/

void cp_normalize_all(double *creal,double *cimag,int icoef_form,
                  double *freal,double *fimag,int ifcoef_form,
                  double *occ,int *ioff,
                  double *norm_save,double *norm_tmp,
                  CP_COMM_STATE_PKG *cp_comm_state_pkg)

/*========================================================================*/
{/*begin routine*/
/*========================================================================*/
/*             Local variable declarations                                */
#include "../typ_defs/typ_mask.h"

  int iii,is,icf,i1,i2;
  double norm,scale,rsnorm,wght;
  int ncoef      = cp_comm_state_pkg->nstate_ncoef_proc;
  int nstate     = cp_comm_state_pkg->nstate;
  int myid_state = cp_comm_state_pkg->myid;
  int np_states  = cp_comm_state_pkg->num_proc;
  MPI_Comm comm   = cp_comm_state_pkg->comm;

/*========================================================================*/
/* 0) Checks                                                              */

    if((icoef_form+ifcoef_form)!=2&&np_states>1){
      printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
      printf("The coefficients must be in transposed form\n");
      printf("on state processor %d in cp_normalize \n",myid_state);
      printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
      fflush(stdout);exit(1);
    }/*endif*/

/*========================================================================*/
/* I) Calculate the norms                                                 */

 wght = 2.0;
 if(myid_state+1==np_states){wght=1.0;}

 for(is=1;is<=nstate;is++){

   norm = 0.0;
   i1 = 1 + ioff[is];
   i2 = ncoef-1 + ioff[is];
   for(icf=i1;icf<=i2; icf++){
     norm += creal[icf]*creal[icf] + cimag[icf]*cimag[icf];    
   }
   icf = ncoef + ioff[is];
   norm_save[is] = 2.0*norm + wght*(creal[icf]*creal[icf]+
                                      cimag[icf]*cimag[icf]);

 }/* endfor is*/

 if(np_states>1){
  for(is=1;is <= nstate;is++){norm_tmp[is] = norm_save[is];}
  Allreduce(&(norm_tmp[1]),&(norm_save[1]),nstate,MPI_DOUBLE,MPI_SUM,0,comm); 
 }/*endif*/

/*========================================================================*/
/* I) Normalize each state                                                */

 for(is=1;is<=nstate;is++){

   rsnorm = sqrt(occ[is]/norm_save[is]);
   i1 = 1 + ioff[is];
   i2 = ncoef + ioff[is];
   for(icf=i1;icf<=i2; icf++){
     creal[icf] *= rsnorm;  
     cimag[icf] *= rsnorm;
     freal[icf] /= rsnorm;
     fimag[icf] /= rsnorm;
   }/*endfor*/

 }/* endfor is*/

/*========================================================================*/
   }/*end routine*/
/*========================================================================*/


/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
/* construct and add the KS contribution to the forces                      */
/*==========================================================================*/

void cp_add_ksmat_force(double *creal,double *cimag,
                        int icoef_form,int icoef_orth,
                        double *fcreal,double *fcimag,
                        int ifcoef_form,int ifcoef_orth,
                        double *ksmat,double *ks_scr,
                        int *ioff,int cp_lsda,int cp_min,
                        double *occ,
                        CP_COMM_STATE_PKG *cp_comm_state_pkg)

/*========================================================================*/
{/*begin routine*/
/*========================================================================*/
/*             Local variable declarations                                */

  int nstate2;
  int iii,i,is,js,icf,jcf;
  int ind,ind1,ind2;
  int ncoefm1,one=1;
  int koff;
  double alpha,beta,g0_wght;

/*             Local pointer declarations                                */
   int np_states    = cp_comm_state_pkg->num_proc;
   int nstate       = cp_comm_state_pkg->nstate;
   int myid_state   = cp_comm_state_pkg->myid;

/*========================================================================*/
/* I) Check the form of coefs and forces                                  */

  if(icoef_orth!=1 || ifcoef_orth!=1){
    printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
    printf("The coefs and coef forces must be in orthogonal form    \n");
    printf("on state processor %d in cp_add_ksmat_force   \n",myid_state);
    printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
    fflush(stdout);exit(1);
  }/*endif*/

  if(np_states>1){
   if(icoef_form!=1 || ifcoef_form!=1){
    printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
    printf("The coefs and coef forces must be in transposed form \n");
    printf("on state processor %d in cp_add_ksmat_force   \n",myid_state);
    printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
    fflush(stdout);exit(1);
   }/*endif*/
  }/*endif*/

/*========================================================================*/
/* II) Compute KS matrix from the coef/orbital overlap  KS=F*C            */

  nstate2 = nstate*nstate;
  for(is=1;is<=nstate2; is++){ksmat[is] = 0.0;}

  alpha=1.0;
  cp_ovlap_mat_diff(creal,cimag,icoef_form,fcreal,fcimag,ifcoef_form,
                    ksmat,ks_scr,alpha,ioff,cp_comm_state_pkg);

/*========================================================================*/
/* III) If minimization: scale the matrix                                 */

  if(cp_min==1){
    alpha = 2.0; 
    koff = 0;
    for(is=1;is<=nstate; is++){
      for(js=1;js<=nstate; js++){
        ksmat[(koff+js)] *= (alpha/occ[is]);    
      }/* endfor */
      koff += nstate;
    }/* endfor */
  }/* endif */

  if(cp_min != 1 && cp_lsda == 1){
     for(is=1;is<=nstate;is++){
	  for(js=1;js<=nstate;js++){
            koff = (is-1)*nstate + js;
            ksmat[koff] *= 2.0/(occ[is]+occ[js]);
	  }/* endfor */
     }/* endfor */
  }/* endif */

  if(cp_min != 1 && cp_lsda == 0){
     for(is=1;is<=nstate;is++){
	  for(js=1;js<=nstate;js++){
            koff = (is-1)*nstate + js;
            ksmat[koff] *= 4.0/(occ[is]+occ[js]);
	  }/* endfor */
     }/* endfor */
  }/* endif */

/*========================================================================*/
/* IV) Add Kohn-Sham matrix contribution to force F -= KS*c               */

  alpha = -1.0;  beta = 1.0; g0_wght = 0.5;
  if(cp_lsda == 1 && cp_min == 0){alpha = -2.0;}
  if(cp_lsda == 0 && nstate == 1 && occ[1] == 1.0 && cp_min == 0)  {alpha = -2.0;}
  
  cp_rotation_prim(creal,icoef_form,fcreal,ifcoef_form,ksmat,
                   alpha,beta,g0_wght,ioff,cp_comm_state_pkg);
  cp_rotation_prim(cimag,icoef_form,fcimag,ifcoef_form,ksmat,
                   alpha,beta,g0_wght,ioff,cp_comm_state_pkg);


/*========================================================================*/
  }/*end routine*/
/*========================================================================*/




/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
/* Rotate back from a basis in which the overlap matrix is diagonal matrix  */
/* of occupation numbers to the original basis                              */
/*==========================================================================*/

void cp_rotate_gen_nonortho(double *creal,double *cimag, 
                            int icoef_form,int *icoef_orth,
                            double *norbmat,int *ioff,
                            double *c1_temp,double *c2_temp,
                            CP_COMM_STATE_PKG *cp_comm_state_pkg)

/*========================================================================*/
   {/*begin routine*/
/*========================================================================*/
/*         Local Variables           */

int myid_state = cp_comm_state_pkg->myid;
int np_states  = cp_comm_state_pkg->num_proc;

/*========================================================================*/
/*  0) Checks                                                             */

  if((*icoef_orth)!=1){
    printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
    printf("The coefs must be in orthogonal form \n");
    printf("on state processor %d in cp_rotate_gen_nonortho   \n",myid_state);
    printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
    fflush(stdout);exit(1);
  }/*endif*/

  if(np_states>1){
   if(icoef_form!=1){
    printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
    printf("The coefficients must be in transposed form \n");
    printf("on state processor %d in cp_rotate_gen_nonortho  \n",myid_state);
    printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
    fflush(stdout);exit(1);
   }/*endif*/
  }/*endif*/

/*========================================================================*/
/*  I) Rotate using the transformation matrix and flip the orth flag     */

  cp_rotate_vector(creal,cimag,icoef_form,norbmat,ioff,
                   c1_temp,c2_temp,cp_comm_state_pkg);
  *icoef_orth = 0;

/*========================================================================*/
  }/*end routine*/
/*========================================================================*/





/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
/* Rotate into a basis in which the overlap matrix is the diagonal matrix   */
/* of occupation numbers                                                    */
/*==========================================================================*/

void cp_rotate_coef_ortho(double *creal,double *cimag,
                          int icoef_form,int *icoef_orth,
                          double *norbmat,double *norbmati,
                          double *ovmat_eigv,
                          double *c1_temp,double *c2_temp,
                          double *occ,int *ioff,
                          double *max_off_diag,double *max_diag,
                          CPSCR_OVMAT *cpscr_ovmat,
                          CP_COMM_STATE_PKG *cp_comm_state_pkg)

/*========================================================================*/
   {/*begin routine*/
/*========================================================================*/
/*         Local Variables           */

  int icoef_orth_now;
  int myid_state = cp_comm_state_pkg->myid;
  int np_states  = cp_comm_state_pkg->num_proc;

/*========================================================================*/
/*  0) Checks                                                             */

  if(*icoef_orth!=0){
    printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
    printf("The coefs must be in nonorthogonal form    \n");
    printf("on state processor %d in cp_rotate_coef_ortho   \n",myid_state);
    printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
    fflush(stdout);exit(1);
  }/*endif*/

  if(np_states>1){
   if(icoef_form!=1){
    printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
    printf("The coefficients must be in transposed form \n");
    printf("on state processor %d in cp_rotate_coef_ortho  \n",myid_state);
    printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
    fflush(stdout);exit(1);
   }/*endif*/
  }/*endif*/

/*========================================================================*/
/*  I) Construct the transformation matrix, its inverse and a             */
/*     transformation matrix to a basis in which S is diagonal but not    */
/*     necessarily equal to the occupation numbers                        */

  icoef_orth_now = 0;
  cp_construct_orth_rotmat(creal,cimag,icoef_form,icoef_orth_now,
                           norbmat,norbmati,ovmat_eigv,
                           occ,ioff,max_off_diag,max_diag,cpscr_ovmat,
                           cp_comm_state_pkg);

/*========================================================================*/
/*  II) Rotate using the transformation matrix and flip the orth flag     */

  cp_rotate_vector(creal,cimag,icoef_form,norbmat,ioff,
                   c1_temp,c2_temp,cp_comm_state_pkg);

  *icoef_orth = 1;

/*========================================================================*/
    }/*end routine*/
/*========================================================================*/




/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
/* Rotate the orthonormal CP vectors to the KS basis                        */
/*==========================================================================*/

void cp_rotate_ks_basis(double *creal,double *cimag,
                        int icoef_form, int icoef_orth,
                        double *vcreal,double *vcimag,
                        int ivcoef_form, int ivcoef_orth, 
                        double *fcreal,double *fcimag,
                        int ifcoef_form, int ifcoef_orth,
                        int *ioff,double *occ,double *rocc_sum,
                        double *ksmat, double *kseig_vals,
                        double *c1_temp,double *c2_temp,
                        CPSCR_OVMAT *cpscr_ovmat,
                        CP_COMM_STATE_PKG *cp_comm_state_pkg,double *kseig_sum)

/*========================================================================*/
   {/*begin routine*/
/*========================================================================*/
/*   Local variables */

  double *kseig_vecs  = cpscr_ovmat->ovlap1;
  double *ks_scr      = cpscr_ovmat->ovlap2;
  double *ovmat       = cpscr_ovmat->ovlap3;
  double *ovmat_scr   = cpscr_ovmat->ovlap4;

  double *rs_scr1    = cpscr_ovmat->state_vec1;
  double *rs_scr2    = cpscr_ovmat->state_vec2;
  int np_states      = cp_comm_state_pkg->num_proc;
  int myid_state     = cp_comm_state_pkg->myid;

/*========================================================================*/
/* 0) Check the form of CP vectors                                        */

   if((icoef_orth+ifcoef_orth+ivcoef_orth)!=3){
    printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
    printf("The coefs and coef forces must be in orthogonal form    \n");
    printf("on state processor %d in cp_rotate_ks_basis   \n",myid_state);
    printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
    fflush(stdout);exit(1);
  }/*endif*/

  if(np_states>1){
   if((icoef_form+ifcoef_form+ivcoef_form)!=3){
    printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
    printf("The CP vectors must be in transposed form \n");
    printf("on state processor %d in cp_rotate_ks_basis   \n",myid_state);
    printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
    fflush(stdout);exit(1);
   }/*endif*/
  }/*endif*/

/*========================================================================*/
/* I) Construct and diagonalize ks matrix                                 */
 
  cp_condiag_ksmat(creal,cimag,icoef_form,icoef_orth,
                   fcreal,fcimag,ifcoef_form,ifcoef_orth,
                   kseig_vals,kseig_vecs,ksmat,ks_scr,rs_scr1,rs_scr2,
                   ioff,cp_comm_state_pkg,kseig_sum);

/*========================================================================*/
/* II) Rotate into the KS-basis (using KS eigenvectors)                   */

  cp_rotate_all(creal,cimag,icoef_form,vcreal,vcimag,ivcoef_form,
                fcreal,fcimag,ifcoef_form,kseig_vecs,ioff,c1_temp,c2_temp,
                cp_comm_state_pkg);


/*========================================================================*/
/* II) Rotate into the KS-basis (using KS eigenvectors)                   */

 rotate_occ_shuffle(creal,cimag,icoef_form,ovmat,ovmat_scr,
                    occ,rocc_sum,ioff,cp_comm_state_pkg);

/*========================================================================*/
   }/*end routine*/
/*========================================================================*/


/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
/*  Construct and diagonalize the KS matrix using orthonormal CP vectors    */
/*==========================================================================*/

void cp_condiag_ksmat(double *creal,double *cimag,
                      int icoef_form, int icoef_orth,
                      double *fcreal,double *fcimag,
                      int ifcoef_form, int ifcoef_orth,
                      double *kseig_vals,double *kseig_vecs,
                      double *ksmat, double *ks_scr,
                      double *rs_scr1,double *rs_scr2,
                      int *ioff,
                      CP_COMM_STATE_PKG *cp_comm_state_pkg,double *kseig_sum)

/*========================================================================*/
   {/*begin routine*/
/*========================================================================*/
/*             Local variable declarations                                */

  int nstate2,i,is,ierr;
  int job=1;
  int ncoef_tot;
  int nline;
  double alpha,beta,g0_wght;
  double kseig_sum_loc;

/*             Local pointer declarations                                */
  int nstate       = cp_comm_state_pkg->nstate;
  int np_states      = cp_comm_state_pkg->num_proc;
  int myid_state     = cp_comm_state_pkg->myid;

/*========================================================================*/
/* 0) Check the form of CP vectors                                        */

  if(icoef_orth!=1 || ifcoef_orth!=1){
    printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
    printf("The coefs and coef forces must be in orthogonal form    \n");
    printf("on state processor %d in cp_condiag_ksmat   \n",myid_state);
    printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
    fflush(stdout);exit(1);
  }/*endif*/

  if(np_states>1){
   if((icoef_form+ifcoef_form)!=2){
    printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
    printf("The CP vectors must be in transposed form \n");
    printf("on state processor %d in cp_condiag_ksmat   \n",myid_state);
    printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
    fflush(stdout);exit(1);
   }/*endif*/
  }/*endif*/

/*========================================================================*/
   /* I) Zero matrix                                                      */

  nstate2 = nstate*nstate;
  for(is=1;is <=nstate2; is++){
    ksmat[is] = 0.0;
  }/*endfor*/

/*========================================================================*/
   /* II) Compute overlap with orbital forces (KS matrix)                */

  alpha=1.0;
  cp_ovlap_mat_diff(creal,cimag,icoef_form,fcreal,fcimag,ifcoef_form,
                    ksmat,ks_scr,alpha,ioff,cp_comm_state_pkg);

  for(is=1;is<=nstate2;is++){
    ksmat[is] *= -0.5;
  }/* endfor */

/*========================================================================*/
   /* III) Diagonalize KS matrix                                          */

  RS(&nstate,&nstate,&(ksmat[1]),&(kseig_vals[1]),&job,&(kseig_vecs[1]),
     &(rs_scr1[1]),&(rs_scr2[1]),&ierr);

/*========================================================================*/
   /* IV) Print out Kohn-Sham eigenvalues                                 */

  if (myid_state==0)
  {
   nline=5;
   PRINT_LINE_DASH
   printf("Kohn-Sham Eigenvalues (Hartrees):\n");
   printf("\n");
   kseig_sum_loc = 0.0;
   for(is=1;is <= nstate;is++){
     kseig_sum_loc += kseig_vals[is];
     printf("%14.8g  ",kseig_vals[is]);
     if((is % nline) == 0) printf("\n");
   } /* endfor */
   printf("\n");
   PRINT_LINE_DASH
   *kseig_sum = kseig_sum_loc;
  }; /* endif */

/*========================================================================*/
   }/*end routine*/
/*========================================================================*/





/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
/*         Construct St^{1/2} and St^{-1/2} such that                       */
/*     St^{-1/2}*S*St^{-1/2} = diag matrix of occ num                       */
/*==========================================================================*/

void cp_construct_orth_rotmat(double *creal,double *cimag,
                              int icoef_form,int icoef_orth,
                              double *norbmat,double *norbmati,
                              double *ovmat_eigv,
                              double *occ,int *ioff,
                              double *max_off_diag,double *max_diag,
                              CPSCR_OVMAT *cpscr_ovmat,
                              CP_COMM_STATE_PKG *cp_comm_state_pkg)

/*========================================================================*/
   {/*begin routine*/
/*========================================================================*/
/*             Local variable declarations                                */

  int ierr,ind,is,js,nstate2;
  int job=1;
  double alpha,beta,occ_rt;
  int itransp = 0;
  int inormal = 1;

/*             Local pointer declarations                                */
  int nstate       = cp_comm_state_pkg->nstate;
  int np_states    = cp_comm_state_pkg->num_proc;
  int myid_state    = cp_comm_state_pkg->myid;
  double *ovmat    = cpscr_ovmat->ovlap1;
  double *ov_scr   = cpscr_ovmat->ovlap2;
  double *oveigs   = cpscr_ovmat->state_vec1;
  double *rs_scr1  = cpscr_ovmat->state_vec2;
  double *rs_scr2  = cpscr_ovmat->state_vec3;
  double *occ_tmp  = cpscr_ovmat->state_vec4;

/*========================================================================*/
/* I) Check the form of the coefs                                         */

  if(icoef_orth!=0){
    printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
    printf("The coefficients must be in nonorthogonal form \n");
    printf("on state processor %d in cp_construct_orth_rotmat  \n",myid_state);
    printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
    fflush(stdout);exit(1);
  }/*endif*/

  if(np_states>1 && icoef_form != 1){
    printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
    printf("The coefficients must be in transposed form \n");
    printf("on state processor %d in cp_construct_orth_rotmat \n",myid_state);
    printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
    fflush(stdout);exit(1);
  }/*endif*/

/*========================================================================*/
/* II) Calculate overlap matrix, S, of nonorthogonal orbitals             */


  alpha=2.0;
  cp_ovlap_mat_same(creal,cimag,icoef_form,ovmat,ov_scr,alpha,ioff,
                    cp_comm_state_pkg);

/*========================================================================*/
/* III) Calculate maximum off-diagonal element of S                       */

  for(is=1;is <= nstate; is++){
    for(js=is+1;js <=nstate; js++){
      ind = (is-1)*nstate + js;
      *max_off_diag = MAX(ovmat[ind],*max_off_diag);
    }/*endfor*/
  }/*endfor*/

/*========================================================================*/
/* IV) Calculate maximum diagonal element of S                            */

  for(is=1;is <= nstate; is++){
      ind = (is-1)*nstate + is;
      *max_diag = MAX(ovmat[ind],*max_diag);
  }/*endfor*/

/*========================================================================*/
/* V) Diagonalize the overlap matrix, S                                   */

  nstate2 = nstate*nstate;
  RS(&nstate,&nstate,&(ovmat[1]),&(oveigs[1]),&job,&(ovmat_eigv[1]),
     &(rs_scr1[1]),&(rs_scr2[1]),&ierr);

/*========================================================================*/
/* Shuffle me                                                             */

 rotate_occ_orth_shuffle(occ,occ_tmp,nstate);


/*========================================================================*/
/* VI) Construct the norb transformation matrix (St^{-1/2} and St^{1/2})  */
     /* i) Get inverse square root of eigenvalues                         */
  for(is=1;is <= nstate; is++){
    oveigs[is]  = 1.0/sqrt(oveigs[is]);
  }/*endfor*/
  
/*------------------------------------------------------------------------*/
    /* ii) Construct the diagonal representation of St^{1/2} and St^{-1/2}*/


  for(is=1;is <= nstate2; is++){
   norbmat[is] = 0.0;
   norbmati[is] = 0.0;
  }/*endfor*/
  for(is=1;is <= nstate; is++){
      ind = (is-1)*nstate + is;
      occ_rt = sqrt(occ_tmp[is]);
      norbmat[ind]  = oveigs[is]*occ_rt;
      norbmati[ind] = 1.0/norbmat[ind];
  }/*endfor*/

/*------------------------------------------------------------------------*/
    /* iii) Rotate diag reps to true reps   */

  alpha = 1.0;  beta = 0.0;
  GEN_MATMUL(&(norbmat[1]),&nstate,&inormal,&(ovmat_eigv[1]),&nstate,&itransp,
             &(ov_scr[1]),&nstate,&nstate,&nstate,&nstate,
             &alpha,&beta);
  GEN_MATMUL(&(ovmat_eigv[1]),&nstate,&inormal,&(ov_scr[1]),&nstate,&inormal,
             &(norbmat[1]),&nstate,&nstate,&nstate,&nstate,
             &alpha,&beta);

  alpha = 1.0;  beta = 0.0;
  GEN_MATMUL(&(norbmati[1]),&nstate,&inormal,&(ovmat_eigv[1]),&nstate,&itransp,
             &(ov_scr[1]),&nstate,&nstate,&nstate,&nstate,
             &alpha,&beta);
  GEN_MATMUL(&(ovmat_eigv[1]),&nstate,&inormal,&(ov_scr[1]),&nstate,&inormal,
             &(norbmati[1]),&nstate,&nstate,&nstate,&nstate,
             &alpha,&beta);

/*========================================================================*/
    }/*end routine*/
/*========================================================================*/


/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
/* Rotate a CP vector (nstate x ncoef) into a different basis               */
/*==========================================================================*/

void cp_rotate_vector(double *creal,double *cimag,int icoef_form, 
                      double *trans_mat,int *ioff,
                      double *c1_temp, double *c2_temp,
                      CP_COMM_STATE_PKG *cp_comm_state_pkg)

/*========================================================================*/
   {/*begin routine*/
/*========================================================================*/
/*           Local variables                      */

  int ntot,i;
  int icoef_form_tmp;
  double alpha,beta,g0_wght;

  int nstate       = cp_comm_state_pkg->nstate;
  int ncoef        = cp_comm_state_pkg->nstate_ncoef_proc_max;
  int np_states    = cp_comm_state_pkg->num_proc;
  int myid_state   = cp_comm_state_pkg->myid;

/*========================================================================*/
/* I) Checks */

  if(np_states>1 && icoef_form !=1){
    printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
    printf("The coefficients must be in transposed form\n");
    printf("on state processor %d in cp_rotate_vector \n",myid_state);
    printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
    fflush(stdout);exit(1);
  }/*endif*/

/*========================================================================*/
/* II) Rotation */

  alpha = 1.0;  beta = 0.0;
  g0_wght = 1.0;
  ntot = ncoef*nstate;
  icoef_form_tmp = icoef_form;

 
  for(i=1;i<=ntot;i++){c1_temp[i]=creal[i];}
  cp_rotation_prim(c1_temp,icoef_form_tmp,creal,icoef_form,trans_mat,
                   alpha,beta,g0_wght,ioff,cp_comm_state_pkg);
  for(i=1;i<=ntot;i++){c2_temp[i]=cimag[i];}
  cp_rotation_prim(c2_temp,icoef_form_tmp,cimag,icoef_form,trans_mat,
                   alpha,beta,g0_wght,ioff,cp_comm_state_pkg);

/*========================================================================*/
    }/*end routine*/
/*========================================================================*/




/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
/* Rotate all CP vectors (nstate x ncoef) onto a different basis             */
/*==========================================================================*/

void cp_rotate_all(double *creal,double *cimag,  int icoef_form,
                   double *vcreal,double *vcimag,int ivcoef_form,
                   double *fcreal,double *fcimag, int ifcoef_form,
                   double *trans_mat,int *ioff,
                   double *c1_temp, double *c2_temp,
                   CP_COMM_STATE_PKG *cp_comm_state_pkg)

/*========================================================================*/
    {/*begin routine*/
/*========================================================================*/
/*           Local variables                      */

  int ntot,i;
  int icoef_form_tmp;
  int ivcoef_form_tmp;
  int ifcoef_form_tmp;
  double alpha,beta,g0_wght;

  int np_states    = cp_comm_state_pkg->num_proc;
  int myid_state   = cp_comm_state_pkg->myid;
  int nstate       = cp_comm_state_pkg->nstate;
  int ncoef        = cp_comm_state_pkg->nstate_ncoef_proc_max;

/*========================================================================*/
/* I) Checks */

  if(np_states>1){
   if((icoef_form+ivcoef_form+ifcoef_form) != 3){
    printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
    printf("The coefficients must be in transposed form \n");
    printf("on state processor %d in cp_rotate_all \n",myid_state);
    printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
    fflush(stdout);exit(1);
   }/*endif*/
  }/*endif*/

/*========================================================================*/
/* II) Rotation */

  alpha = 1.0;  beta = 0.0;
  g0_wght = 1.0;
  ntot = ncoef*nstate;

/*-----------------------------------------------------------------------*/
/*  i) Coefs  */

  icoef_form_tmp  = icoef_form;
  for(i=1;i<=ntot;i++){c1_temp[i]=creal[i];}
  cp_rotation_prim(c1_temp,icoef_form_tmp,creal,icoef_form,trans_mat,
                   alpha,beta,g0_wght,ioff,cp_comm_state_pkg);
  for(i=1;i<=ntot;i++){c2_temp[i]=cimag[i];}
  cp_rotation_prim(c2_temp,icoef_form_tmp,cimag,icoef_form,trans_mat,
                   alpha,beta,g0_wght,ioff,cp_comm_state_pkg);

/*-----------------------------------------------------------------------*/
/*  ii) Coef vels  */

  ivcoef_form_tmp = ivcoef_form;
  for(i=1;i<=ntot;i++){c1_temp[i]=vcreal[i];}
  cp_rotation_prim(c1_temp,ivcoef_form_tmp,vcreal,ivcoef_form,trans_mat,
                   alpha,beta,g0_wght,ioff,cp_comm_state_pkg);
  for(i=1;i<=ntot;i++){c2_temp[i]=vcimag[i];}
  cp_rotation_prim(c2_temp,ivcoef_form_tmp,vcimag,ivcoef_form,trans_mat,
                   alpha,beta,g0_wght,ioff,cp_comm_state_pkg);

/*-----------------------------------------------------------------------*/
/*  iii) Coef forces  */

  ifcoef_form_tmp = ifcoef_form;
  for(i=1;i<=ntot;i++){c1_temp[i]=fcreal[i];}
  cp_rotation_prim(c1_temp,ifcoef_form_tmp,fcreal,ifcoef_form,trans_mat,
                   alpha,beta,g0_wght,ioff,cp_comm_state_pkg);
  for(i=1;i<=ntot;i++){c2_temp[i]=fcimag[i];}
  cp_rotation_prim(c2_temp,ifcoef_form_tmp,fcimag,ifcoef_form,trans_mat,
                   alpha,beta,g0_wght,ioff,cp_comm_state_pkg);

/*========================================================================*/
    }/*end routine*/
/*========================================================================*/



/*==========================================================================*/
/* Rotation primitive: C2 = alpha*MxC1 + beta*C2                            */
/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void cp_rotation_prim(double *c1,int i1coef_form,
                      double *c2,int i2coef_form,
                      double *matrix,double alpha_in,double beta_in,
                      double g0_wgt,int *ioff_st,
                      CP_COMM_STATE_PKG *cp_comm_state_pkg)

/*========================================================================*/
{/*begin routine*/
/*========================================================================*/
/*             Local variable declarations                                */

  int ioff,joff,is,js,iii;
  int ind;
  int ncoefm1;

  int np_states             = cp_comm_state_pkg->num_proc;
  int myid_state            = cp_comm_state_pkg->myid;
  int nstate                = cp_comm_state_pkg->nstate;
  int nstate_max            = cp_comm_state_pkg->nstate_max;
  int ncoef                 = cp_comm_state_pkg->ncoef;
  int itransp = 0;
  int inormal = 1;
  double beta = beta_in;
  double alpha=alpha_in;

/*========================================================================*/
/* I)  Serial */

 if(np_states == 1){

  ncoefm1 = ncoef-1;

  GEN_MATMUL(&(c1[1]),&ncoef,&inormal,&(matrix[1]),&nstate,&inormal,
             &(c2[1]),&ncoef,&ncoefm1,&nstate,&nstate,
             &alpha,&beta);

   for(is=1;is <=nstate; is++){
    ioff = ncoef + ioff_st[is];
     c2[ioff] *= beta_in;
     for(js=1;js <=nstate; js++){
      joff = ncoef + ioff_st[js];
      ind = (is-1)*nstate + js;
      c2[ioff] += g0_wgt*alpha*matrix[ind]*c1[joff];
     }/* endfor js*/
   }/* endfor is */

 }/*endif*/

/*========================================================================*/
/* II) Parallel */

 if(np_states > 1){

   if((i1coef_form+i2coef_form) !=2){
      printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
      printf("The coefficients must be in transposed form\n");
      printf("on state processor %d in cp_rotation_prim \n",myid_state);
      printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
      fflush(stdout);exit(1);
   }/*endif*/

   cp_rotate_prim_par(c1,i1coef_form,c2,i2coef_form,matrix,alpha_in,beta_in,
                      g0_wgt,cp_comm_state_pkg);

 }/* endif */

/*========================================================================*/
   }/*end routine*/
/*========================================================================*/






/*==========================================================================*/
/* Parallel Rotation primitive: C2 = alpha*MxC1 + beta*C2                   */
/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void cp_rotate_prim_par(double *c1,int i1coef_form,
                        double *c2,int i2coef_form,
                        double *matrix,double alpha,double beta,
                        double g0_wgt,CP_COMM_STATE_PKG *cp_comm_state_pkg)

/*========================================================================*/
{/*begin routine*/
/*========================================================================*/
/*             Local variable declarations                                */

#include "../typ_defs/typ_mask.h"

 double wgt_use,sum,alpha_loc,beta_loc;
 int ioff,joff,is,js,ig,ioff_c,i;
 int ind,jnd,iproc,iii,nread,isum;
 int itransp = 0;
 int inormal = 1;

 int ncoef_tot;
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

    if((i1coef_form+i2coef_form) !=2){
      printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
      printf("The coefficients must be in transposed form\n");
      printf("on state processor %d in cp_rotate_prim_par \n",myid);
      printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
      fflush(stdout);exit(1);
    }/*endif*/

/*========================================================================*/
/* 0) Initialize matrix                                                   */


 ncoef_tot=nstate_ncoef_proc_max*nstate_max;
 for(is=1;is <= ncoef_tot;is++){c2[is] *= beta;}

/*========================================================================*/
/* II) Matrix multiply                                                    */

  wgt_use = (myid < num_proc-1 ? 1.0 : g0_wgt);

#ifdef HAND_ROT
  ioff = 0;
  for(is=1;is<=nstate;is++){
   joff = 0;
   for(js=1;js<=nstate;js++){
     ind  = (is-1)*nstate + js;
     for(ig=1;ig<=nstate_ncoef_proc-1;ig++){
       c2[(ioff+ig)] += alpha*matrix[ind]*c1[(joff+ig)];
     }/*endfor*/
     ig = nstate_ncoef_proc;
     c2[(ioff+ig)] += wgt_use*alpha*matrix[ind]*c1[(joff+ig)];
     joff += nstate_ncoef_proc_max;
   }/*endfor*/
   ioff += nstate_ncoef_proc_max;
  }/*endfor*/
#else
  alpha_loc = alpha;
  beta_loc  = beta;
  GEN_MATMUL(&(c1[1]),&nstate_ncoef_proc_max,&inormal,&(matrix[1]),
                &nstate,&inormal,&(c2[1]),&nstate_ncoef_proc_max,
                &nstate_ncoef_proc_max,&nstate,&nstate,&alpha_loc,&beta_loc);

  if(wgt_use!=1.0){
    ioff = 0;
    for(is=1;is<=nstate;is++){
      joff = 0;
      for(js=1;js<=nstate;js++){
        ind  = (is-1)*nstate + js;
        ig = nstate_ncoef_proc;
        c2[(ioff+ig)] += (wgt_use-1.0)*alpha*matrix[ind]*c1[(joff+ig)];
        joff += nstate_ncoef_proc_max;
      }/*endfor*/
      ioff += nstate_ncoef_proc_max;
    }/*endfor*/
  }/*endif : wght_use*/
#endif

/*========================================================================*/
 }/* end routine */
/*========================================================================*/


/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void rotate_occ_shuffle(double *creal,double *cimag,int icoef_form,
                       double *ovmat,double *ovmat_scr,
                       double *occ,double *rocc_sum,
                       int *ioff,CP_COMM_STATE_PKG *cp_comm_state_pkg)
/*========================================================================*/
  {/*begin routine*/
/*========================================================================*/
/*             Local variable declarations                                */

  int i,idiag,iii,j,iocc,ind1;
  int occ_temp;
  double alpha;
  int np_states = cp_comm_state_pkg->num_proc;
  int myid      = cp_comm_state_pkg->myid;
  int nstate    = cp_comm_state_pkg->nstate;

/*========================================================================*/
/* I) Check the forms */

  if(np_states>1){
   if(icoef_form  !=1){
      printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
      printf("The coefficients must be in transposed form\n");
      printf("on state processor %d in rotate_occ_shuffle \n",myid);
      printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
      fflush(stdout);exit(1);
   }/*endif*/
  }/*endif*/

/*========================================================================*/
/* II) Calculate overlap matrix of nonorthogonal orbitals                 */

  alpha=2.0;
  cp_ovlap_mat_same(creal,cimag,icoef_form,ovmat,ovmat_scr,alpha,ioff,
                    cp_comm_state_pkg);

/*========================================================================*/
/* III) Rejigger the occupation numbers                                   */ 

  for(idiag = 1; idiag <= nstate; idiag++){
     ind1 = (idiag -1 )*nstate + idiag;
     occ_temp    = NINT(ovmat[ind1]);
     occ[idiag]  = occ_temp;
  }/*endfor*/

/*  recalculate rocc_sum_up */
  iocc=0;
  for(i=1;i<= nstate ;i++){
   for(j=1;j<= nstate ;j++){
     iocc++;
     rocc_sum[iocc] = 1.0/(occ[i]+occ[j]);
   }/*endfor i*/
  }/* endfor j*/
   
/*========================================================================*/
}/*end routine*/
/*========================================================================*/


/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void  rotate_occ_orth_shuffle(double *occ,double *occ_tmp,int nstate)

/*=======================================================================*/
/*            Begin subprogram:                                          */
{/*begin routine*/
/*=======================================================================*/
/*          Local variable declarations                                  */

  int i;

  for(i=1;i<= nstate; i++){
      occ_tmp[i] = occ[i];
  }/*endfor*/

  if( nstate > 1){
   occ_sort(nstate,occ_tmp);
  }/*endif*/

/*-----------------------------------------------------------------------*/
} /*end routine*/ 
/*==========================================================================*/




/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void occ_sort(int n, double *index)

/*=======================================================================*/
/*            Begin subprogram:                                          */
{/*begin routine*/
/*=======================================================================*/
/*          Local variable declarations                                  */

  int m,ir,i,j,rjndex,iii;
  double rindex;
  int k,*kndex,*mask,isub,temp;

/*=======================================================================*/
/* I) Setup                        */

  m  = n/2+1;
  ir = n;

/*=======================================================================*/
/* II) Sort array index */

  for(;;){

/*---------------------------------------------------------------------*/
/*  A)hire rindex */
    if(m>1){ 
      m--;
      rindex = index[m];
/*--------------------------------------------------------------------*/
/*  B)retire/promote index[1] */
    }else{
      rindex = index[ir];
      index[ir]=index[1];
      ir--;
      if(ir==1){
       index[1]=rindex;
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
       i=j;
       j=2*j;
      }else{
       /*    c)if no demotations exit while */
       j=ir+1;
      }/*endif*/
    } /*endwhile*/
    /*    d)slot rindex */
    index[i] = rindex;
  }/*endfor*/

/*-----------------------------------------------------------------------*/
} /*end routine*/ 
/*==========================================================================*/












