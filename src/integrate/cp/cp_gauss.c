/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
/*                                                                          */
/*                         PI_MD:                                           */
/*             The future of simulation technology                          */
/*             ------------------------------------                         */
/*                     Module: cp_gauss                                     */
/*                                                                          */
/* This subprogram contains functions necessary to integrating              */
/* the CP equations using the Gaussian dynamics method                      */
/*                                                                          */
/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

#include "standard_include.h"
#include "../typ_defs/typedefs_gen.h"
#include "../typ_defs/typedefs_cp.h"
#include "../typ_defs/typedefs_class.h"
#include "../typ_defs/typedefs_bnd.h"
#include "../proto_defs/proto_integrate_cp_local.h"

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void add_gauss_force(CPCOEFFS_INFO *cpcoeffs_info,CPCOEFFS_POS *cpcoeffs_pos,
                     CPSCR_OVMAT *cpscr_ovmat,CPOPTS *cpopts)

/*========================================================================*/
{/*begin routine*/
/*========================================================================*/
/*             Local variable declarations                                */

 int i,is,icoef;
 int ncoef,nstate_up,nstate_dn;
 int nstate_max;
 double rat;
 double penalty;
 double *gforce  = cpscr_ovmat->ovlap1;
 double *dovlap  = cpscr_ovmat->ovlap2;
 double *cre_up  = cpcoeffs_pos->cre_up;
 double *cim_up  = cpcoeffs_pos->cim_up;
 double *cre_dn  = cpcoeffs_pos->cre_dn;
 double *cim_dn  = cpcoeffs_pos->cim_dn;
 double *fcre_up = cpcoeffs_pos->fcre_up;
 double *fcim_up = cpcoeffs_pos->fcim_up;
 double *fcre_dn = cpcoeffs_pos->fcre_dn;
 double *fcim_dn = cpcoeffs_pos->fcim_dn;
 int *ioff_up    = cpcoeffs_info->ioff_up;
 int *ioff_dn    = cpcoeffs_info->ioff_dn;

/*========================================================================*/
/* 0) Define useful constants                                             */

  ncoef = cpcoeffs_info->ncoef;
  nstate_up = cpcoeffs_info->nstate_up;
  nstate_dn = cpcoeffs_info->nstate_dn;
  nstate_max = MAX(nstate_up,nstate_dn);

/*========================================================================*/
/* I) Compute diagonal overlap of coeffs with force and coeffs with coeffs*/
  
 penalty = 0.0;
 for(is=1;is <= nstate_up;is++){
  gforce[is] = 0.0;
  dovlap[is] = 0.0;
 }

 for(is=1;is<=nstate_up;is++) {
   for(i=1;i<=ncoef-1;i++) {
     icoef = i+ioff_up[is];
     dovlap[is] += cre_up[icoef]*cre_up[icoef] +
                   cim_up[icoef]*cim_up[icoef];
   }
  icoef = ncoef+ioff_up[is];
  dovlap[is] = 2.0*dovlap[is] + cre_up[icoef]*cre_up[icoef] +
                                cim_up[icoef]*cim_up[icoef];
  penalty += (dovlap[is]-2.0)*(dovlap[is]-2.0);
  dovlap[is] = 2.0;
 }
 

 for(is=1;is<=nstate_up;is++) {
   for(i=1;i<=ncoef;i++) {
     icoef = i+ioff_up[is];
     gforce[is] += cre_up[icoef]*fcre_up[icoef] +
                   cim_up[icoef]*fcim_up[icoef];
   }
   gforce[is] = 2.0*gforce[is];
 }
 printf("penalty function %.12g\n",penalty);

/*========================================================================*/
/* II) Add Gauss force to cp force                                        */

 for(is=1;is<=nstate_up;is++) {
   rat = gforce[is]/dovlap[is];
   for(i=1;i<=ncoef-1;i++) {
     icoef = i+ioff_up[is];
     fcre_up[icoef] -= cre_up[icoef]*rat;
     fcim_up[icoef] -= cim_up[icoef]*rat;
   }
  icoef = ncoef+ioff_up[is];
  fcre_up[icoef] -= 0.5*cre_up[icoef]*rat;
 }

/*========================================================================*/
/* III) Do same if LSDA                                                   */

  if( (cpopts->cp_lsda) && (nstate_dn != 0) ){
/*------------------------------------------------------------------------*/
/* i) Diagoal overlaps                                                    */

  for(is=1;is <= nstate_dn;is++){
   gforce[is] = 0.0;
   dovlap[is] = 0.0;
  } 

  for(is=1;is<=nstate_dn;is++) {
    for(i=1;i<=ncoef-1;i++) {
      icoef = i+ioff_dn[is];
      dovlap[is] += cre_dn[icoef]*cre_dn[icoef] +
                    cim_dn[icoef]*cim_dn[icoef];
    }
   icoef = ncoef+ioff_dn[is];
   dovlap[is] = 2.0*dovlap[is] + cre_dn[icoef]*cre_dn[icoef] +
                                 cim_dn[icoef]*cim_dn[icoef];
  } 

  for(is=1;is<=nstate_dn;is++) {
    for(i=1;i<=ncoef;i++) {
      icoef = i+ioff_dn[is];
      gforce[is] += cre_dn[icoef]*fcre_dn[icoef] +
                    cim_dn[icoef]*fcim_dn[icoef];
    }
    gforce[is] = 2.0*gforce[is];
  }

/*========================================================================*/
/* ii) Add Gauss force to cp force                                        */

  for(is=1;is<=nstate_dn;is++) {
    rat = gforce[is]/dovlap[is];
    for(i=1;i<=ncoef-1;i++) {
      icoef = i+ioff_dn[is];
      fcre_dn[icoef] -= cre_dn[icoef]*rat;
      fcim_dn[icoef] -= cim_dn[icoef]*rat;
    }
   icoef = ncoef+ioff_dn[is];
   fcre_dn[icoef] -= 0.5*cre_dn[icoef]*rat;
   fcim_dn[icoef] -= 0.5*cim_dn[icoef]*rat;
  }

/*--------------------------------------------------------------------------*/
 }/* endif lsda */
/*--------------------------------------------------------------------------*/

/*--------------------------------------------------------------------------*/
}/*end routine*/
/*==========================================================================*/



/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void apply_lgauss(CPCOEFFS_INFO *cpcoeffs_info,CPCOEFFS_POS *cpcoeffs_pos,
                  CPSCR_WAVE *cpscr_wave,
                  CPSCR_OVMAT *cpscr_ovmat,CPOPTS *cpopts, double dt)

/*========================================================================*/
{/*begin routine*/
/*========================================================================*/
/*             Local variable declarations                                */
 int i,is,icoef,iii;
 int ncoef,nstate_up,nstate_dn;
 int nstate_max;
 double rat;
 double k1,k2,k3,discr,sdelta;
 double cre_0,vcre_0;
 double arg1,arg2;
 double *dovlap   = cpscr_ovmat->ovlap1;
 double *cre_up   = cpcoeffs_pos->cre_up;
 double *cim_up   = cpcoeffs_pos->cim_up;
 double *cre_dn   = cpcoeffs_pos->cre_dn;
 double *cim_dn   = cpcoeffs_pos->cim_dn;
 double *vcre_up  = cpcoeffs_pos->vcre_up;
 double *vcim_up  = cpcoeffs_pos->vcim_up;
 double *vcre_dn  = cpcoeffs_pos->vcre_dn;
 double *vcim_dn  = cpcoeffs_pos->vcim_dn;
 double *a_up     = cpscr_wave->cre_up;
 double *b_up     = cpscr_wave->cim_up;
 double *a_dn     = cpscr_wave->cre_dn;
 double *b_dn     = cpscr_wave->cim_dn;
 int *ioff_up     = cpcoeffs_info->ioff_up;
 int *ioff_dn     = cpcoeffs_info->ioff_dn;

/*========================================================================*/
/* 0) Define useful constants                                             */

  ncoef = cpcoeffs_info->ncoef;
  nstate_up = cpcoeffs_info->nstate_up;
  nstate_dn = cpcoeffs_info->nstate_dn;
  nstate_max = MAX(nstate_up,nstate_dn);

/*========================================================================*/
/* I) Compute diagonal overlap of coeffs with force and coeffs with coeffs*/
  
 for(is=1;is <= nstate_up;is++){
  dovlap[is] = 0.0;
 }

 for(is=1;is<=nstate_up;is++) {
   for(i=1;i<=ncoef-1;i++) {
     icoef = i+ioff_up[is];
     dovlap[is] += cre_up[icoef]*cre_up[icoef] +
                   cim_up[icoef]*cim_up[icoef];
   }
  icoef = ncoef+ioff_up[is];
  dovlap[is] = 2.0*dovlap[is] + cre_up[icoef]*cre_up[icoef] +
                                cim_up[icoef]*cim_up[icoef];
  dovlap[is] = 2.0;
 }

/*========================================================================*/
/* II) Loop over states                                                   */

 for(is=1;is <= nstate_up;is++){
/*--------------------------------------------------------------------------*/
/* i) Compute constants of integration for differential equations           */

   icoef = ncoef+ioff_up[is];
   vcre_0 = vcre_up[icoef];
   cre_0  =  cre_up[icoef];
   for(i=1;i<=ncoef-1;i++){
    icoef = i+ioff_up[is];
    a_up[icoef] = vcre_up[icoef] - cre_up[icoef]*vcre_0/cre_0;
    b_up[icoef] = vcim_up[icoef] - cim_up[icoef]*vcre_0/cre_0;
   }

/*--------------------------------------------------------------------------*/
/* ii) Compute coefficients in quadratic differential equation              */

  rat = 1.0/dovlap[is];
  k1=0.0;
  k2=0.0;
  k3=0.0;
  for(i=1;i<=ncoef-1;i++){
   icoef = i+ioff_up[is];
   k1 -= rat*(cre_up[icoef]*cre_up[icoef] + cim_up[icoef]*cim_up[icoef]);
   k2 -= rat*(cre_up[icoef]*a_up[icoef]   + cim_up[icoef]*b_up[icoef]);
   k3 -= rat*(a_up[icoef]*a_up[icoef] + b_up[icoef]*b_up[icoef]);
  }
  k1 = 2.0*k1 - cre_0*cre_0*rat;
  k1 /= cre_0;
  k2 *= 4.0;
  k3 *= 2.0*cre_0;

 
/*--------------------------------------------------------------------------*/
/* iii) Solve a nasty differential equation in one variable                 */

  discr = 4.0*k1*k3 - k2*k2;
  if(discr < 0.0){
   printf("Discriminant less than 0 -- check for bugs\n");
   exit(1);
  }
  sdelta = sqrt(discr);
  arg1 = (2.0*k1*vcre_0 + k2);
  arg2 = atan2(arg1,sdelta) + 0.25*sdelta*dt;
  vcre_0 = (sdelta*tan(arg2) - k2)/(2.0*k1);


/*--------------------------------------------------------------------------*/
/* iv) Update coefficient velocities using new vcre_0                       */

  for(i=1;i<=ncoef-1;i++){
   icoef = i+ioff_up[is];
   vcre_up[icoef] = cre_up[icoef]*vcre_0/cre_0 + a_up[icoef];
   vcim_up[icoef] = cim_up[icoef]*vcre_0/cre_0 + b_up[icoef];
  }
 

 } /* end loop over states */

/*========================================================================*/
/* II) Same for down states if LSDA                                       */

 if( (cpopts->cp_lsda == 1) && (nstate_dn != 0) ){
/*========================================================================*/
/* I) Compute diagonal overlap of coeffs with force and coeffs with coeffs*/
  
 for(is=1;is <= nstate_dn;is++){
  dovlap[is] = 0.0;
 }

 for(is=1;is<=nstate_dn;is++) {
   for(i=1;i<=ncoef-1;i++) {
     icoef = i+ioff_dn[is];
     dovlap[is] += cre_dn[icoef]*cre_dn[icoef] +
                   cim_dn[icoef]*cim_dn[icoef];
   }
  icoef = ncoef+ioff_dn[is];
  dovlap[is] = 2.0*dovlap[is] + cre_dn[icoef]*cre_dn[icoef] +
                                cim_dn[icoef]*cim_dn[icoef];
 }

/*--------------------------------------------------------------------------*/
/* i) Compute constants of integration for differential equations           */

   icoef = ncoef+ioff_dn[is];
   vcre_0 = vcre_dn[icoef];
   cre_0  =  cre_dn[icoef];
   for(i=1;i<=ncoef-1;i++){
    icoef = i+ioff_dn[is];
    a_dn[icoef] = vcre_dn[icoef] - cre_dn[icoef]*vcre_0/cre_0;
    b_dn[icoef] = vcim_dn[icoef] - cim_dn[icoef]*vcre_0/cre_0;
   }

/*--------------------------------------------------------------------------*/
/* ii) Compute coefficients in quadratic differential equation              */

 for(is=1;is <= nstate_dn;is++){

  rat = 1.0/dovlap[is];
  k1=0.0;
  k2=0.0;
  k3=0.0;
  for(i=1;i<=ncoef-1;i++){
   icoef = i+ioff_dn[is];
   k1 -= rat*(cre_dn[icoef]*cre_dn[icoef] + cim_dn[icoef]*cim_dn[icoef]);
   k2 -= rat*(cre_dn[icoef]*a_dn[icoef]   + cim_dn[icoef]*b_dn[icoef]);
   k3 -= rat*(a_dn[icoef]*a_dn[icoef] + b_dn[icoef]*b_dn[icoef]);
  }
  k1 = 2.0*k1 - cre_0*cre_0*rat;
  k1 /= cre_0;
  k2 *= 4.0;
  k3 *= 2.0;
 
/*--------------------------------------------------------------------------*/
/* iii) Solve a nasty differential equation in one variable                 */

  discr = 4.0*k1*k3 - k2*k2;
  if(discr < 0.0){
   printf("Discriminant less than 0 -- check for bugs\n");
   exit(1);
  }
  sdelta = sqrt(discr);
  arg1 = (2.0*k1*vcre_0 + k2)/sdelta;
  arg2 = atan(arg1) + 0.25*sdelta*dt;
  vcre_0 = (sdelta*tan(arg2) - k2)/(2.0*k1);

/*--------------------------------------------------------------------------*/
/* iv) Update coefficient velocities using new vcre_0                       */

  for(i=1;i<=ncoef-1;i++){
   icoef = i+ioff_dn[is];
   vcre_dn[icoef] = cre_dn[icoef]*vcre_0/cre_0 + a_dn[icoef];
   vcim_dn[icoef] = cim_dn[icoef]*vcre_0/cre_0 + b_dn[icoef];
 }

 } /* end loop over states */
/*--------------------------------------------------------------------------*/
 }/* endif lsda */

/*--------------------------------------------------------------------------*/
}/*end routine*/
/*==========================================================================*/




