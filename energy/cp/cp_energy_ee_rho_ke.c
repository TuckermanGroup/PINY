 /*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
/*                                                                          */
/*                         PI_MD:                                           */
/*             The future of simulation technology                          */
/*             ------------------------------------                         */
/*                Module: control_energy_ee_rho_ke.c                        */
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

void cp_ke_dens_calc_hybrid(CPEWALD *cpewald,CPSCR *cpscr,
                            CELL *cell,double *creal, double *cimag,
                            int icoef_form,int icoef_orth,
                            double *elec_ke_dens,
                            int nstate,int ncoef,int cp_dual_grid_opt,
                            COMMUNICATE *communicate,
                            PARA_FFT_PKG3D *cp_para_fft_pkg3d_lg,
                            PARA_FFT_PKG3D *cp_sclr_fft_pkg3d_lg,
                            PARA_FFT_PKG3D *cp_para_fft_pkg3d_dens_cp_box,
                            PARA_FFT_PKG3D *cp_sclr_fft_pkg3d_dens_cp_box,
                            PARA_FFT_PKG3D *cp_sclr_fft_pkg3d_sm)
      
/*==========================================================================*/
{/*begin routine*/
/*==========================================================================*/
#include "../typ_defs/typ_mask.h"

/* local variables                                                  */

 int iii,ioff,ioff2;
 int is,i,j,iupper,init;
 int component;
 double vol_cp,rvol_cp;
 double rho_tot,rhog_tot,vol,rvol;

/*  Assign local pointers                                           */
 double *zfft           =    cpscr->cpscr_wave.zfft;
 double *zfft_tmp       =    cpscr->cpscr_wave.zfft_tmp;
 double *cre_scr        =    cpscr->cpscr_wave.cre_up;
 double *cim_scr        =    cpscr->cpscr_wave.cim_up;
 double *rhocr_scr       =    cpscr->cpscr_rho.rhocr_scr;
 double *rhoci_scr       =    cpscr->cpscr_rho.rhoci_scr;
 double *rho_scr;
 double *hmati_cp       =    cell->hmati_cp;
 double *hmat_cp        =    cell->hmat_cp;
 int   myid_state       =    communicate->myid_state;
 int   np_states        =    communicate->np_states;

 int  *recv_counts_coef =    cp_para_fft_pkg3d_lg->recv_counts_coef;
 int   nfft_proc        =    cp_para_fft_pkg3d_lg->nfft_proc;
 int   nfft             =    cp_para_fft_pkg3d_lg->nfft;
 int   nfft2            =    nfft/2;
 int   nfft2_proc       =    nfft_proc/2;
 int   box_rat          =    cpewald->box_rat;
 int   ncoef_l          =    cp_para_fft_pkg3d_lg->ncoef;
 int   nfft_dens_cp_box,nfft2_dens_cp_box;
 int   nfft_proc_dens_cp_box,nfft2_proc_dens_cp_box;
 int   ncoef_l_dens_cp_box;
 double *rhocr_dens_cp_box;
 double *rhoci_dens_cp_box;
 int    *recv_counts_coef_dens_cp_box;
 int    *map_dual       = cpscr->cpscr_rho.map_dual;
 double integral;
 double int_tmp;

 MPI_Comm comm_states   =    communicate->comm_states;

 if(cp_dual_grid_opt >= 1){
   nfft_dens_cp_box       =  cp_para_fft_pkg3d_dens_cp_box->nfft;
   nfft_proc_dens_cp_box  =  cp_para_fft_pkg3d_dens_cp_box->nfft_proc;
   nfft2_dens_cp_box      =  nfft_dens_cp_box/2;
   nfft2_proc_dens_cp_box =  nfft_proc_dens_cp_box/2;


   ncoef_l_dens_cp_box  =  cp_para_fft_pkg3d_dens_cp_box->ncoef;
   rhocr_dens_cp_box = cpscr->cpscr_rho.rhocr_up_dens_cp_box;
   rhoci_dens_cp_box = cpscr->cpscr_rho.rhoci_up_dens_cp_box;
   recv_counts_coef_dens_cp_box = 
             cp_para_fft_pkg3d_dens_cp_box->recv_counts_coef;

 }/*endif cp_dual_grid_opt*/

 if(np_states > 1){
   rho_scr = cpscr->cpscr_rho.v_ks_up; 
 } else {
   rho_scr = elec_ke_dens;
 }

  
/* ================================================================= */
/*0) Check the form of the coefficients                              */

  if(icoef_orth!=1){
    printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
    printf("The coefficients must be in orthogonal form    \n");
    printf("on state processor %d in cp_ke_dens_calc_hybrid   \n",myid_state);
    printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
    fflush(stdout);exit(1);
  }/*endif*/

  if(np_states>1 && icoef_form == 1){
    printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
    printf("The coefficients must be in normal (not transposed) \n");
    printf("form on state processor %d in cp_ke_dens_calc_hybrid  \n",myid_state);
    printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
    fflush(stdout);exit(1);
  }/*endif*/


/* ================================================================= */
/*1) zero density and gradients if necessary                         */

 init=0;
 if(cp_dual_grid_opt >= 1){
  for(i=1; i<= nfft2_dens_cp_box ; i++){
      rho_scr[i] = 0.0;
  }/*endfor*/
 }else{
  for(i=1; i<= nfft2 ; i++){
      rho_scr[i] = 0.0;
  }/*endfor*/
 }/*endif cp_dual_grid_opt*/

/*==========================================================================*/
/*==========================================================================*/
/*2) Each component of grad psi_i will be done like the density             */
/*   Two states at a time using the scalar package                          */

 iupper = nstate;
 if((nstate % 2) != 0){ 
     iupper = nstate-1; 
 }/* endif */

/*==========================================================================*/
/* Loop over components                                                     */

 for(component = 1; component <= 3; component++){
  get_igcoef_hybrid(creal,cimag,cre_scr,cim_scr,cpewald,hmati_cp,
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

/*--------------------------------------------------------------------------*/
/* III) add the square of the two wave functions to the density(real space) */

    if(cp_dual_grid_opt >= 1){
      sum_rho(zfft,rho_scr,cp_sclr_fft_pkg3d_dens_cp_box); 
    }else{
      sum_rho(zfft,rho_scr,cp_sclr_fft_pkg3d_lg);
    }/*endif cp_dual_grid_opt*/

  }/*endfor is*/

/*--------------------------------------------------------------------------*/
/*IV) if there is an odd number of states, use single pack
      for the last state (i caen't be chaengin' the laws of physics
      captn, i've got to be usin the single pack!!!!!!)                     */

  if((nstate % 2 ) != 0) {
     ioff = (nstate-1)*ncoef;
     sngl_pack_coef(&cre_scr[ioff],&cim_scr[ioff],zfft,cp_sclr_fft_pkg3d_sm);

/*--------------------------------------------------------------------------*/
/* V) fourier transform the last wavefunction to real space
       convention exp(-igr)                                                 */

     para_fft_gen3d_fwd_to_r(zfft,zfft_tmp,cp_sclr_fft_pkg3d_sm);

/*--------------------------------------------------------------------------*/
/*VI) add the square of the last wave function to the density(real space)   */

    if(cp_dual_grid_opt >= 1){ 
       sum_rho(zfft,rho_scr,cp_sclr_fft_pkg3d_dens_cp_box); 
    }else{
       sum_rho(zfft,rho_scr,cp_sclr_fft_pkg3d_lg);
    }/*endif cp_dual_grid_opt*/
  }/*endif*/

 }/*endfor components */

/*=========================================================================*/
/*=========================================================================*/
/*  3) get density in g space                                              */

 if(cp_dual_grid_opt >= 1){ 
  sngl_pack_rho(zfft,rho_scr,cp_sclr_fft_pkg3d_dens_cp_box); 
 }else{
  sngl_pack_rho(zfft,rho_scr,cp_sclr_fft_pkg3d_lg); 
 }/*endif cp_dual_grid_opt*/


/*--------------------------------------------------------------------------*/
/*  II) back transform to g-space  convention exp(igr)                      */

 if(cp_dual_grid_opt >= 1){ 
  para_fft_gen3d_bck_to_g(zfft,zfft_tmp,cp_sclr_fft_pkg3d_dens_cp_box); 
 }else{
  para_fft_gen3d_bck_to_g(zfft,zfft_tmp,cp_sclr_fft_pkg3d_lg); 
 }/*endif cp_dual_grid_opt*/

/*==========================================================================*/
/*  III) unpack the density                                                 */ 

 if(cp_dual_grid_opt == 0){
  if(np_states == 1){
    sngl_upack_coef(rhocr_scr,rhoci_scr,zfft,cp_sclr_fft_pkg3d_lg);
  }else{
    sngl_upack_coef(zfft_tmp,&zfft_tmp[ncoef_l],zfft,cp_sclr_fft_pkg3d_lg);
  }/*endif*/
 }else{
  if(np_states == 1){
    sngl_upack_coef(rhocr_scr,rhoci_scr,zfft,
                    cp_sclr_fft_pkg3d_dens_cp_box);
  }else{
    sngl_upack_coef(zfft_tmp,&zfft_tmp[ncoef_l_dens_cp_box],zfft,
                    cp_sclr_fft_pkg3d_dens_cp_box);
  }/*endif*/
 }/*endif cp_dual_grid_opt*/


/*==========================================================================*/
/* VII) Reduce rho in g space and get all of it in real space               */
/*      Now in g-level parallel and need parallel packages                  */

 if(np_states > 1){
   if(cp_dual_grid_opt == 0){
      Reduce_scatter(&zfft_tmp[1],&rhocr_scr[1],&recv_counts_coef[1],
                     MPI_DOUBLE,MPI_SUM,comm_states);
      Reduce_scatter(&zfft_tmp[(ncoef_l+1)],&rhoci_scr[1],&recv_counts_coef[1],
                     MPI_DOUBLE,MPI_SUM,comm_states);  

      sngl_pack_coef(rhocr_scr,rhoci_scr,zfft,cp_para_fft_pkg3d_lg);

      para_fft_gen3d_fwd_to_r(zfft,zfft_tmp,cp_para_fft_pkg3d_lg);

      sngl_upack_rho(zfft,elec_ke_dens,cp_para_fft_pkg3d_lg); 

   Barrier(comm_states);
   }else{
     Reduce_scatter(&zfft_tmp[1],&rhocr_scr[1],
                    &recv_counts_coef_dens_cp_box[1],
                    MPI_DOUBLE,MPI_SUM,comm_states);
      Reduce_scatter(&zfft_tmp[(ncoef_l_dens_cp_box+1)],
                     &rhoci_scr[1],&recv_counts_coef_dens_cp_box[1],
                     MPI_DOUBLE,MPI_SUM,comm_states);  

      sngl_pack_coef(rhocr_scr,rhoci_scr,zfft,
                     cp_para_fft_pkg3d_dens_cp_box);

      para_fft_gen3d_fwd_to_r(zfft,zfft_tmp,cp_para_fft_pkg3d_dens_cp_box);

      sngl_upack_rho(zfft,elec_ke_dens,cp_para_fft_pkg3d_dens_cp_box); 
   }/*endif cp_dual_grid_opt*/
 }/*endif np_states*/

/*===========================================================================*/
/* IV) finish the density in real space by dividing by the volume           */
/* DUALED SYSTEMS only keep the real space density on the cp_grid            */

     vol_cp  = getdeth(hmat_cp);
     rvol_cp = 1.0/vol_cp;

     if(cp_dual_grid_opt >= 1){
      for(i=1 ; i<= nfft2_proc_dens_cp_box;i++){
         elec_ke_dens[i] *= rvol_cp;
      }/*endfor*/
     }else{
      for(i=1 ; i<= nfft2_proc;i++){
         elec_ke_dens[i] *= rvol_cp;
      }/*endfor*/
     }/*endif cp_dual_grid_opt*/

/*==============================================================*/
/* Debug code */

#define DEBUG
#ifdef DEBUG
  int_tmp = 0.0;
  for(i=1; i<= nfft2_proc ; i++){
    int_tmp += elec_ke_dens[i];
  }/*endfor*/
  printf("%d %.12g\n",myid_state,int_tmp);
  if(np_states>1) {
    Reduce(&int_tmp,&integral,1,MPI_DOUBLE,MPI_SUM,0,comm_states);
  } else {
    integral=int_tmp;
  }
  if(myid_state==0){
    printf("integral = %.12g\n",integral*vol_cp/((double)nfft2));
  }
#endif

/*==============================================================*/
}/*end routine*/
/*==============================================================*/





/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void cp_ke_dens_calc_full_g(CPEWALD *cpewald,CPSCR *cpscr,
                            CELL *cell,double *creal, double *cimag,
                            int icoef_form,int icoef_orth,
                            double *elec_ke_dens,
                            int nstate,int ncoef,int cp_dual_grid_opt,
                            COMMUNICATE *communicate,
                            PARA_FFT_PKG3D *cp_para_fft_pkg3d_lg,
                            PARA_FFT_PKG3D *cp_para_fft_pkg3d_dens_cp_box,
                            PARA_FFT_PKG3D *cp_para_fft_pkg3d_sm)
      
/*==========================================================================*/
{/*begin routine*/
/*==========================================================================*/
#include "../typ_defs/typ_mask.h"

/* local variables                                                  */

 int iii,ioff,ioff2;
 int is,i,j,iupper,init;
 int component;
 double vol_cp,rvol_cp;
 double rho_tot,rhog_tot,vol,rvol;

/*  Assign local pointers                                           */
 double *zfft           =    cpscr->cpscr_wave.zfft;
 double *zfft_tmp       =    cpscr->cpscr_wave.zfft_tmp;
 double *cre_scr        =    cpscr->cpscr_wave.cre_up;
 double *cim_scr        =    cpscr->cpscr_wave.cim_up;
 double *rho_scr;
 double *hmati_cp       =    cell->hmati_cp;
 double *hmat_cp        =    cell->hmat_cp;
 int   myid_state       =    communicate->myid_state;
 int   np_states        =    communicate->np_states;

 int icoef_off          =    cp_para_fft_pkg3d_sm->icoef_off;
 int ncoef_use          =    cp_para_fft_pkg3d_sm->ncoef_use;
 int  *recv_counts_coef =    cp_para_fft_pkg3d_lg->recv_counts_coef;
 int   nfft_proc        =    cp_para_fft_pkg3d_lg->nfft_proc;
 int   nfft             =    cp_para_fft_pkg3d_lg->nfft;
 int   nfft2            =    nfft/2;
 int   nfft2_proc       =    nfft_proc/2;
 int   box_rat          =    cpewald->box_rat;
 int   nfft_dens_cp_box,nfft2_dens_cp_box;
 int   nfft_proc_dens_cp_box,nfft2_proc_dens_cp_box;
 int   ncoef_l_dens_cp_box;
 double *rhocr_dens_cp_box;
 double *rhoci_dens_cp_box;
 int    *recv_counts_coef_dens_cp_box;
 int    *map_dual       = cpscr->cpscr_rho.map_dual;
 double int_tmp,integral;

 MPI_Comm comm_states   =    communicate->comm_states;

 if(cp_dual_grid_opt >= 1){
   nfft_dens_cp_box       =  cp_para_fft_pkg3d_dens_cp_box->nfft;
   nfft_proc_dens_cp_box  =  cp_para_fft_pkg3d_dens_cp_box->nfft_proc;
   nfft2_dens_cp_box      =  nfft_dens_cp_box/2;
   nfft2_proc_dens_cp_box =  nfft_proc_dens_cp_box/2;


   ncoef_l_dens_cp_box  =  cp_para_fft_pkg3d_dens_cp_box->ncoef;
   rhocr_dens_cp_box = cpscr->cpscr_rho.rhocr_up_dens_cp_box;
   rhoci_dens_cp_box = cpscr->cpscr_rho.rhoci_up_dens_cp_box;
   recv_counts_coef_dens_cp_box = 
             cp_para_fft_pkg3d_dens_cp_box->recv_counts_coef;

 }/*endif cp_dual_grid_opt*/

 rho_scr = elec_ke_dens;

/* ================================================================= */
/*0) Check the form of the coefficients                              */

  if(icoef_orth!=1){
    printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
    printf("The coefficients must be in orthogonal form    \n");
    printf("on state processor %d in cp_ke_dens_full_g   \n",myid_state);
    printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
    fflush(stdout);exit(1);
  }/*endif*/

  if(np_states>1 && icoef_form == 0){
    printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
    printf("The coefficients must be in transposed (not normal) \n");
    printf("form on state processor %d in cp_ke_dens_full_g  \n",myid_state);
    printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
    fflush(stdout);exit(1);
  }/*endif*/
/* ================================================================= */
/*1) zero density and gradients if necessary                         */

 init=0;
 if(cp_dual_grid_opt >= 1){
  for(i=1; i<= nfft2_dens_cp_box ; i++){
      rho_scr[i] = 0.0;
  }/*endfor*/
 }else{
  for(i=1; i<= nfft2_proc ; i++){
      rho_scr[i] = 0.0;
  }/*endfor*/
 }/*endif cp_dual_grid_opt*/

/*==========================================================================*/
/*==========================================================================*/
/*2) Each component of grad psi_i will be done like the density             */
/*   Two states at a time using the scalar package                          */

 iupper = nstate;
 if((nstate % 2) != 0){ 
     iupper = nstate-1; 
 }/* endif */

/*==========================================================================*/
/* Loop over components                                                     */


 for(component = 1; component <= 3; component++){

  get_igcoef_full_g(creal,cimag,cre_scr,cim_scr,cpewald,hmati_cp,
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

/*--------------------------------------------------------------------------*/
/* III) add the square of the two wave functions to the density(real space) */

    if(cp_dual_grid_opt >= 1){
      sum_rho(zfft,rho_scr,cp_para_fft_pkg3d_dens_cp_box); 
    }else{
      sum_rho(zfft,rho_scr,cp_para_fft_pkg3d_lg);
    }/*endif cp_dual_grid_opt*/

  }/*endfor is*/

/*--------------------------------------------------------------------------*/
/*IV) if there is an odd number of states, use single pack
      for the last state (i caen't be chaengin' the laws of physics
      captn, i've got to be usin the single pack!!!!!!)                     */

  if((nstate % 2 ) != 0) {
     ioff = (nstate-1)*ncoef;
     sngl_pack_coef(&cre_scr[ioff],&cim_scr[ioff],zfft,cp_para_fft_pkg3d_sm);

/*--------------------------------------------------------------------------*/
/* V) fourier transform the last wavefunction to real space
       convention exp(-igr)                                                 */

     para_fft_gen3d_fwd_to_r(zfft,zfft_tmp,cp_para_fft_pkg3d_sm);

/*--------------------------------------------------------------------------*/
/*VI) add the square of the last wave function to the density(real space)   */

    if(cp_dual_grid_opt >= 1){ 
       sum_rho(zfft,rho_scr,cp_para_fft_pkg3d_dens_cp_box); 
    }else{
       sum_rho(zfft,rho_scr,cp_para_fft_pkg3d_lg);
    }/*endif cp_dual_grid_opt*/
  }/*endif*/

 }/*endfor components */

/*===========================================================================*/
/* IV) finish the density in real space by dividing by the volume           */
/* DUALED SYSTEMS only keep the real space density on the cp_grid            */

     vol_cp  = getdeth(hmat_cp);
     rvol_cp = 1.0/vol_cp;

     if(cp_dual_grid_opt >= 1){
      for(i=1 ; i<= nfft2_proc_dens_cp_box;i++){
         elec_ke_dens[i] *= rvol_cp;
      }/*endfor*/
     }else{
      for(i=1 ; i<= nfft2_proc;i++){
         elec_ke_dens[i] *= rvol_cp;
      }/*endfor*/
     }/*endif cp_dual_grid_opt*/

/*==============================================================*/
/* Debug code */

#define DEBUG
#ifdef DEBUG
  int_tmp = 0.0;
  for(i=1; i<= nfft2_proc ; i++){
    int_tmp += elec_ke_dens[i];
  }/*endfor*/
  printf("%d %.12g\n",myid_state,int_tmp);
  if(np_states>1) {
    Reduce(&int_tmp,&integral,1,MPI_DOUBLE,MPI_SUM,0,comm_states);
  } else {
    integral=int_tmp;
  }
  if(myid_state==0){
    printf("integral = %.12g\n",integral*vol_cp/((double)nfft2));
  }
#endif

/*==============================================================*/
}/*end routine*/
/*==============================================================*/




/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void get_igcoef_hybrid(double *creal,double *cimag, 
                       double *cre_scr,double *cim_scr,
                       CPEWALD *cpewald,double *hmati,
                       int nstate,int ncoef,int component)
      
/*=========================================================================*/
{/*begin routine*/
/*=========================================================================*/

/*assign local pointers */
     int    *kastore        =    cpewald->kastr_sm;
     int    *kbstore        =    cpewald->kbstr_sm;
     int    *kcstore        =    cpewald->kcstr_sm;

/* local variables */
      
   double  aka,akb,akc,xk,yk,zk,tpi;
   int is,icount;
   int icoef_ind;

/*===================================================================*/
/* I) Get ig times rho (gradient in g-space)                         */

  tpi = 2.0*M_PI;
  switch(component){

   case 1:

    for(is=1;is<=nstate;is++){
      for(icount = 1; icount<= ncoef ; icount++){
         aka = (double)kastore[icount];
         akb = (double)kbstore[icount];
         akc = (double)kcstore[icount];
         xk = (aka*hmati[1] +akb*hmati[2] +akc*hmati[3])*tpi;
         icoef_ind = (is-1)*ncoef + icount;
         cre_scr[icoef_ind] = xk*cimag[icoef_ind];
         cim_scr[icoef_ind] = -xk*creal[icoef_ind];
      }/* endfor icount */
    }/*endfor is*/
    break;

   case 2:

    for(is=1;is<=nstate;is++){
      for(icount = 1; icount<= ncoef ; icount++){
         aka = (double)kastore[icount];
         akb = (double)kbstore[icount];
         akc = (double)kcstore[icount];
         yk = (aka*hmati[4] +akb*hmati[5] +akc*hmati[6])*tpi;
         icoef_ind = (is-1)*ncoef + icount;
         cre_scr[icoef_ind] = yk*cimag[icoef_ind];
         cim_scr[icoef_ind] = -yk*creal[icoef_ind];
      }/* endfor icount */
    }/*endfor is*/
    break;

   case 3:

    for(is=1;is<=nstate;is++){
      for(icount = 1; icount<= ncoef ; icount++){
         aka = (double)kastore[icount];
         akb = (double)kbstore[icount];
         akc = (double)kcstore[icount];
         zk = (aka*hmati[7] +akb*hmati[8] +akc*hmati[9])*tpi;
         icoef_ind = (is-1)*ncoef + icount;
         cre_scr[icoef_ind] = zk*cimag[icoef_ind];
         cim_scr[icoef_ind] = -zk*creal[icoef_ind];
      }/* endfor icount */
    }/*endfor is*/
      break;

  }/* end switch */

  cre_scr[ncoef] = 0.0;
  cim_scr[ncoef] = 0.0;

/*===================================================================*/
}/*end routine*/
/*===================================================================*/






/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void get_igcoef_full_g(double *creal,double *cimag, 
                       double *cre_scr,double *cim_scr,
                       CPEWALD *cpewald,double *hmati,
                       int icoef_off,int nstate,
                       int ncoef,int ncoef_use,int component)
      
/*=========================================================================*/
{/*begin routine*/
/*=========================================================================*/

/*assign local pointers */
     int    *kastore        =    cpewald->kastr_sm;
     int    *kbstore        =    cpewald->kbstr_sm;
     int    *kcstore        =    cpewald->kcstr_sm;

/* local variables */
      
   double  aka,akb,akc,xk,yk,zk,tpi;
   int is,icount;
   int icoef_ind;

/*===================================================================*/
/* I) Get ig times rho (gradient in g-space)                         */

  tpi = 2.0*M_PI;
  switch(component){

   case 1:

    for(is=1;is<=nstate;is++){
      for(icount = 1; icount<= ncoef ; icount++){
         aka = (double)kastore[(icount+icoef_off)];
         akb = (double)kbstore[(icount+icoef_off)];
         akc = (double)kcstore[(icount+icoef_off)];
         xk = (aka*hmati[1] +akb*hmati[2] +akc*hmati[3])*tpi;
         icoef_ind = (is-1)*ncoef + icount;
         cre_scr[icoef_ind] = xk*cimag[icoef_ind];
         cim_scr[icoef_ind] = -xk*creal[icoef_ind];
      }/* endfor icount */
    }/*endfor is*/
    break;

   case 2:

    for(is=1;is<=nstate;is++){
      for(icount = 1; icount<= ncoef ; icount++){
         aka = (double)kastore[(icount+icoef_off)];
         akb = (double)kbstore[(icount+icoef_off)];
         akc = (double)kcstore[(icount+icoef_off)];
         yk = (aka*hmati[4] +akb*hmati[5] +akc*hmati[6])*tpi;
         icoef_ind = (is-1)*ncoef + icount;
         cre_scr[icoef_ind] = yk*cimag[icoef_ind];
         cim_scr[icoef_ind] = -yk*creal[icoef_ind];
      }/* endfor icount */
    }/*endfor is*/
    break;

   case 3:

    for(is=1;is<=nstate;is++){
      for(icount = 1; icount<= ncoef ; icount++){
         aka = (double)kastore[(icount+icoef_off)];
         akb = (double)kbstore[(icount+icoef_off)];
         akc = (double)kcstore[(icount+icoef_off)];
         zk = (aka*hmati[7] +akb*hmati[8] +akc*hmati[9])*tpi;
         icoef_ind = (is-1)*ncoef + icount;
         cre_scr[icoef_ind] = zk*cimag[icoef_ind];
         cim_scr[icoef_ind] = -zk*creal[icoef_ind];
      }/* endfor icount */
    }/*endfor is*/
      break;

  }/* end switch */

  if(ncoef_use != ncoef){
   cre_scr[ncoef] = 0.0;
   cim_scr[ncoef] = 0.0;
  }/* endif */

/*===================================================================*/
}/*end routine*/
/*===================================================================*/






