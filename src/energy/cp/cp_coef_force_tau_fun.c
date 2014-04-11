/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
/*                                                                          */
/*                         PI_MD:                                           */
/*             The future of simulation technology                          */
/*             ------------------------------------                         */
/*                   Module: control_energy_ee_rho.c                        */
/*                                                                          */
/*                                                                          */
/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/


#include "standard_include.h"
#include "../typ_defs/typedefs_class.h"
#include "../typ_defs/typedefs_gen.h"
#include "../typ_defs/typedefs_cp.h"
#include "../proto_defs/proto_energy_cp_local.h"
#include "../proto_defs/proto_math.h"
#include "../proto_defs/proto_communicate_wrappers.h"


/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
/* calculate the force on the coeff's (up or down, both)  */
/* arising from electron KE density-dependent functionals. */
/*==========================================================================*/

void coef_force_tau_fun_hybrid(CPEWALD *cpewald,int nstate,
                               double *creal,double *cimag, 
                               double *fcreal,double  *fcimag,
                               double *cre_scr,double *cim_scr,
                               double *zfft,double *zfft_tmp,
                               double *v_ks_tau,double *ak2_sm,
                               double *pvten_cp,
                               int cp_ptens_calc,double *hmati,
                               COMMUNICATE *communicate,
                               int icoef_form,int icoef_orth,int ifcoef_form,
                               PARA_FFT_PKG3D *cp_sclr_fft_pkg3d_sm)

/*==========================================================================*/
/*         Begin Routine                                                    */
{/*Begin Routine*/
/*=======================================================================*/
/*            Local variable declarations                               */
#include "../typ_defs/typ_mask.h"

   int is,i,iupper,component;
   double tpi;
   double aka,akb,akc,gk,cfact;
   double eke;
   int ioff,ncoef1,ioff2;
   int iii,iis,iis2,nis;

   int nfft       = cp_sclr_fft_pkg3d_sm->nfft;
   int ncoef      = cp_sclr_fft_pkg3d_sm->ncoef;
   int myid_state = communicate->myid_state;
   int np_states  = communicate->np_states;

/*            Local pointers                                       */

   int  *kastore_sm    =  cpewald->kastr_sm;
   int  *kbstore_sm    =  cpewald->kbstr_sm;
   int  *kcstore_sm    =  cpewald->kcstr_sm; 



#define DEBUG_OFF
#ifdef DEBUG
      int icount;
      double c_g,g2,anorm,sum,vol,cre_now,cim_now;
      double dx,x_pos,y_pos,z_pos,phase_r,phase_i,arg;
      FILE *fp;
#endif

/* ================================================================= */
/*0) Check the form of the coefficients                              */

  if(icoef_orth!=1){
    printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
    printf("The coefficients must be in orthogonal form    \n");
    printf("on state processor %d in coef_force_calc_hybrid\n",myid_state);
    printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
    fflush(stdout);exit(1);
  }/*endif*/

  if(np_states>1) 
   if( ifcoef_form==1 || icoef_form == 1){
    printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
    printf("The coefs/coef forces must be in normal (not transposed) \n");
    printf("form on state processor %d in coef_force_calc_hybrid  \n",
           myid_state);
    printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
    fflush(stdout);exit(1);
  }/*endif*/

/*=================================================================*/
/*  Find the upper state limit                                     */

  tpi = 2.0*M_PI;
  ncoef1 = ncoef - 1;
  iupper = nstate;
  if(nstate % 2 == 1){
     iupper = nstate - 1;
  }

/*==========================================================================*/
/* 1.  Loop over components                                                 */

 for(component = 1; component <= 3; component++){
  get_igcoef_hybrid(creal,cimag,cre_scr,cim_scr,cpewald,hmati,
                    nstate,ncoef,component);

  for(is = 1; is <= iupper; is = is + 2){
     ioff   = (is-1)*ncoef;
     ioff2 = (is)*ncoef;

/*--------------------------------------------------------------------------*/
/*I) pack the complex zfft array with two wavefunctions (real)
    the wavefunctions are reperesented in spherically cuttof
    half g space                                                            */

     dble_pack_coef(&cre_scr[ioff],&cim_scr[ioff],&cre_scr[ioff2],&cim_scr[ioff2],
                    zfft,cp_sclr_fft_pkg3d_sm);

/*--------------------------------------------------------------------------*/
/*II) fourier transform the two wavefunctions to real space
     convention exp(-igr)                                                   */

     para_fft_gen3d_fwd_to_r(zfft,zfft_tmp,cp_sclr_fft_pkg3d_sm);

/*==========================================================================*/
/* 2.  For this component calculate V_TAU |(grad psi)_a > a= x, y or z      */

    cp_vpsi(zfft,v_ks_tau,nfft);

/*--------------------------------------------------------------------------*/
/*  III) fourier transform  to g-space                                       */
/*     convention exp(igr)                                                  */

    para_fft_gen3d_bck_to_g(zfft,zfft_tmp,cp_sclr_fft_pkg3d_sm);

/*==========================================================================*/
/* 3.  For this component, unpack the coefficients and get a contribution   */
/*     to the dot product that is added to the coef forces                  */

/*--------------------------------------------------------------------------*/
/* I) Double unpack */

    dble_upack_coef(&cre_scr[ioff],&cim_scr[ioff],
                    &cre_scr[ioff2],&cim_scr[ioff2],
                    zfft,cp_sclr_fft_pkg3d_sm);

/*--------------------------------------------------------------------------*/
/* II) Calculate part of a dot product                                      */

    for(i=1; i<= ncoef1; i++){

       aka = (double)kastore_sm[i];
       akb = (double)kbstore_sm[i];
       akc = (double)kcstore_sm[i];

       if(component==1) gk = (aka*hmati[1] +akb*hmati[2] +akc*hmati[3])*tpi;
       if(component==2) gk = (aka*hmati[4] +akb*hmati[5] +akc*hmati[6])*tpi;
       if(component==3) gk = (aka*hmati[7] +akb*hmati[8] +akc*hmati[9])*tpi;

       iis = ioff + i;
       iis2 = ioff2 + i;  
       fcreal[iis] += 4.0*gk*cim_scr[iis];
       fcimag[iis] -= 4.0*gk*cre_scr[iis];
       fcreal[iis2] += 4.0*gk*cim_scr[iis2];
       fcimag[iis2] -= 4.0*gk*cre_scr[iis2];

    }/* endfor coef */

/*--------------------------------------------------------------------------*/
/* III) End of state loop                                                   */

  }/*endfor is */

/*==========================================================================*/
/* 4.  If there is an odd number of states, do the same using sngl pack     */

  if((nstate % 2 ) != 0) {

     ioff = (nstate-1)*ncoef;
     sngl_pack_coef(&cre_scr[ioff],&cim_scr[ioff],zfft,cp_sclr_fft_pkg3d_sm);

/*--------------------------------------------------------------------------*/
/* I) fourier transform the last wavefunction to real space
       convention exp(-igr)                                                 */

     para_fft_gen3d_fwd_to_r(zfft,zfft_tmp,cp_sclr_fft_pkg3d_sm);

/*--------------------------------------------------------------------------*/
/* II)  For this component calculate V_TAU |(grad psi)_a > a= x, y or z     */

    cp_vpsi(zfft,v_ks_tau,nfft);

/*--------------------------------------------------------------------------*/
/*   III) fourier transform the result back to g-space                      */
/*         convention exp(igr)                                              */

      para_fft_gen3d_bck_to_g(zfft,zfft_tmp,cp_sclr_fft_pkg3d_sm);

/*==========================================================================*/
/* 5. get forces on coefficients by single unpacking the array zfft         */


      sngl_upack_coef(&cre_scr[ioff],&cim_scr[ioff],zfft,cp_sclr_fft_pkg3d_sm);

/*--------------------------------------------------------------------------*/
/* I) Calculate part of a dot product                                      */

    for(i=1; i<= ncoef1; i++){
       aka = (double)kastore_sm[i];
       akb = (double)kbstore_sm[i];
       akc = (double)kcstore_sm[i];

       if(component==1) gk = (aka*hmati[1] +akb*hmati[2] +akc*hmati[3])*tpi;
       if(component==2) gk = (aka*hmati[4] +akb*hmati[5] +akc*hmati[6])*tpi;
       if(component==3) gk = (aka*hmati[7] +akb*hmati[8] +akc*hmati[9])*tpi;
  
       iis = ioff+i;       
       fcreal[iis] += 4.0*gk*cim_scr[iis];
       fcimag[iis] -= 4.0*gk*cre_scr[iis];

    }/* endfor coef */

/*--------------------------------------------------------------------------*/
/* II) Endif odd number of states                                           */

   }/* endif: odd number of states*/

/*==========================================================================*/

  }/* endfor components */

/*==========================================================================*/
 }/*end routine*/
/*==========================================================================*/




/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
/* calculate the force on the coeff's (up or down, both)  */
/* arising from electron KE density-dependent functionals. */
/*==========================================================================*/

void coef_force_tau_fun_full_g(CPEWALD *cpewald,int nstate,
                               double *creal,double *cimag, 
                               double *fcreal,double  *fcimag,
                               double *cre_scr,double *cim_scr,
                               double *zfft,double *zfft_tmp,
                               double *v_ks_tau,double *ak2_sm,
                               double *pvten_cp,
                               int cp_ptens_calc,double *hmati,
                               COMMUNICATE *communicate,
                               int icoef_form,int icoef_orth,int ifcoef_form,
                               PARA_FFT_PKG3D *cp_para_fft_pkg3d_sm)

/*==========================================================================*/
/*         Begin Routine                                                    */
{/*Begin Routine*/
/*=======================================================================*/
/*            Local variable declarations                               */
#include "../typ_defs/typ_mask.h"

   int is,i,iupper,component;
   double tpi;
   double aka,akb,akc,gk,cfact;
   double eke;
   int ioff,ioff2;
   int iii,iis,iis2,nis;

   int icoef_off  = cp_para_fft_pkg3d_sm->icoef_off;
   int nfft       = cp_para_fft_pkg3d_sm->nfft_proc;
   int ncoef      = cp_para_fft_pkg3d_sm->ncoef_proc;
   int ncoef_use  = cp_para_fft_pkg3d_sm->ncoef_use;
   int myid_state = communicate->myid_state;
   int np_states  = communicate->np_states;

/*            Local pointers                                       */

   int  *kastore_sm    =  cpewald->kastr_sm;
   int  *kbstore_sm    =  cpewald->kbstr_sm;
   int  *kcstore_sm    =  cpewald->kcstr_sm; 



#define DEBUG_OFF
#ifdef DEBUG
      int icount;
      double c_g,g2,anorm,sum,vol,cre_now,cim_now;
      double dx,x_pos,y_pos,z_pos,phase_r,phase_i,arg;
      FILE *fp;
#endif

/* ================================================================= */
/*0) Check the form of the coefficients                              */

  if(icoef_orth!=1){
    printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
    printf("The coefficients must be in orthogonal form    \n");
    printf("on state processor %d in coef_force_tau_calc_full_g\n",myid_state);
    printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
    fflush(stdout);exit(1);
  }/*endif*/

  if(np_states>1) 
   if( ifcoef_form==0 || icoef_form == 0){
    printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
    printf("The coefs/coef forces must be in transposed (not normal) \n");
    printf("form on state processor %d in coef_force_tau_calc_full_g  \n",
           myid_state);
    printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
    fflush(stdout);exit(1);
  }/*endif*/

/*=================================================================*/
/*  Find the upper state limit                                     */


  tpi = 2.0*M_PI;
  iupper = nstate;
  if(nstate % 2 == 1){
     iupper = nstate - 1;
  }

/*==========================================================================*/
/* 1.  Loop over components                                                 */

 for(component = 1; component <= 3; component++){

   get_igcoef_full_g(creal,cimag,cre_scr,cim_scr,cpewald,hmati,
                     icoef_off,nstate,ncoef,ncoef_use,component);

  for(is = 1; is <= iupper; is = is + 2){

     ioff   = (is-1)*ncoef;
     ioff2 = (is)*ncoef;

/*--------------------------------------------------------------------------*/
/*I) pack the complex zfft array with two wavefunctions (real)
    the wavefunctions are reperesented in spherically cuttof
    half g space                                                            */

     dble_pack_coef(&cre_scr[ioff],&cim_scr[ioff],&cre_scr[ioff2],&cim_scr[ioff2],
                    zfft,cp_para_fft_pkg3d_sm);

/*--------------------------------------------------------------------------*/
/*II) fourier transform the two wavefunctions to real space
     convention exp(-igr)                                                   */

     para_fft_gen3d_fwd_to_r(zfft,zfft_tmp,cp_para_fft_pkg3d_sm);

/*==========================================================================*/
/* 2.  For this component calculate V_TAU |(grad psi)_a > a= x, y or z      */

    cp_vpsi(zfft,v_ks_tau,nfft);

/*--------------------------------------------------------------------------*/
/*  III) fourier transform  to g-space                                      */
/*     convention exp(igr)                                                  */

    para_fft_gen3d_bck_to_g(zfft,zfft_tmp,cp_para_fft_pkg3d_sm);

/*==========================================================================*/
/* 3.  For this component, unpack the coefficients and get a contribution   */
/*     to the dot product that is added to the coef forces                  */

/*--------------------------------------------------------------------------*/
/* I) Double unpack */

    dble_upack_coef(&cre_scr[ioff],&cim_scr[ioff],
                    &cre_scr[ioff2],&cim_scr[ioff2],
                    zfft,cp_para_fft_pkg3d_sm);

/*--------------------------------------------------------------------------*/
/* II) Calculate part of a dot product                                      */

    for(i=1; i<= ncoef_use; i++){

       aka = (double)kastore_sm[i];
       akb = (double)kbstore_sm[i];
       akc = (double)kcstore_sm[i];

       if(component==1) gk = (aka*hmati[1] +akb*hmati[2] +akc*hmati[3])*tpi;
       if(component==2) gk = (aka*hmati[4] +akb*hmati[5] +akc*hmati[6])*tpi;
       if(component==3) gk = (aka*hmati[7] +akb*hmati[8] +akc*hmati[9])*tpi;

       iis = ioff+i;
       iis2 = ioff2+i;
       fcreal[iis] += 4.0*gk*cim_scr[iis];
       fcimag[iis] -= 4.0*gk*cre_scr[iis];
       fcreal[iis2] += 4.0*gk*cim_scr[iis2];
       fcimag[iis2] -= 4.0*gk*cre_scr[iis2];

    }/* endfor coef */

/*--------------------------------------------------------------------------*/
/* III) End of state loop                                                   */

  }/*endfor is */

/*==========================================================================*/
/* 4.  If there is an odd number of states, do the same using sngl pack     */

  if((nstate % 2 ) != 0) {

     ioff = (nstate-1)*ncoef;
     sngl_pack_coef(&cre_scr[ioff],&cim_scr[ioff],zfft,cp_para_fft_pkg3d_sm);

/*--------------------------------------------------------------------------*/
/* I) fourier transform the last wavefunction to real space
       convention exp(-igr)                                                 */

     para_fft_gen3d_fwd_to_r(zfft,zfft_tmp,cp_para_fft_pkg3d_sm);

/*--------------------------------------------------------------------------*/
/* II)  For this component calculate V_TAU |(grad psi)_a > a= x, y or z     */

    cp_vpsi(zfft,v_ks_tau,nfft);

/*--------------------------------------------------------------------------*/
/*   III) fourier transform the result back to g-space                      */
/*         convention exp(igr)                                              */

      para_fft_gen3d_bck_to_g(zfft,zfft_tmp,cp_para_fft_pkg3d_sm);

/*==========================================================================*/
/* 5. get forces on coefficients by single unpacking the array zfft         */


      sngl_upack_coef(&cre_scr[ioff],&cim_scr[ioff],zfft,
                          cp_para_fft_pkg3d_sm);

/*--------------------------------------------------------------------------*/
/* I) Calculate part of a dot product                                      */

    for(i=1; i<= ncoef_use; i++){
       aka = (double)kastore_sm[i];
       akb = (double)kbstore_sm[i];
       akc = (double)kcstore_sm[i];

       if(component==1) gk = (aka*hmati[1] +akb*hmati[2] +akc*hmati[3])*tpi;
       if(component==2) gk = (aka*hmati[4] +akb*hmati[5] +akc*hmati[6])*tpi;
       if(component==3) gk = (aka*hmati[7] +akb*hmati[8] +akc*hmati[9])*tpi;
  
       iis = ioff+i;
       fcreal[iis] += 4.0*gk*cim_scr[iis];
       fcimag[iis] -= 4.0*gk*cre_scr[iis];

    }/* endfor coef */

/*--------------------------------------------------------------------------*/
/* II) Endif odd number of states                                           */

   }/* endif: odd number of states*/

/*==========================================================================*/

  }/* endfor components */

/*==========================================================================*/
 }/*end routine*/
/*==========================================================================*/





