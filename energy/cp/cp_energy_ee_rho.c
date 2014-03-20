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
/*  cp_rho_calc:                                                            */ 
/*      calculates the density (up,down or total)                           */
/*      using the orbitals. Double packing of ffts and other                */
/*      goodies (most studdly coding)                                       */
/*      Even the option to do full g space (cp_rho_calc_full_g)             */
/*      or the studdly hybrid option in (cp_rho_calc_hybrid)                */
/*==========================================================================*/

void cp_rho_calc_hybrid(CPEWALD *cpewald,CPSCR *cpscr,
                        CPCOEFFS_INFO *cpcoeffs_info,EWALD *ewald,
                        CELL *cell,double *creal, double *cimag,
                        int icoef_form,int icoef_orth,
                        double *rhocr ,double *rhoci,double *rho,
                        double *rhocr_dens_cp_box,double *rhoci_dens_cp_box,
                        double *del_rho_x, double *del_rho_y,
                        double *del_rho_z,  
                        double *del2_rho,int nstate,int ncoef,
                        int cp_gga,int cp_dual_grid_opt,
                        int n_interp_pme_dual,
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
 int is,i,iupper;
 double vol_cp,rvol_cp;
 double temp_r,temp_i;

/*  Assign local pointers                                           */
 double *zfft           =    cpscr->cpscr_wave.zfft;
 double *zfft_tmp       =    cpscr->cpscr_wave.zfft_tmp;
 double *rho_scr;

 double *hmati_cp       =    cell->hmati_cp;
 double *hmat_cp        =    cell->hmat_cp;

 double dbox_rat        =    cpewald->dbox_rat;
 double *bw_r           =    cpscr->cpscr_dual_pme.bw_r;
 double *bw_i           =    cpscr->cpscr_dual_pme.bw_i;
 
 int cp_elf_calc_frq    =    cpcoeffs_info->cp_elf_calc_frq; 

 int   myid_state       =    communicate->myid_state;
 int   np_states        =    communicate->np_states;
 int   laplacian_on     =    cpcoeffs_info->cp_laplacian_on;

 int  *recv_counts_coef =    cp_para_fft_pkg3d_lg->recv_counts_coef;
 int   ncoef_l          =    cp_para_fft_pkg3d_lg->ncoef;
 int   ncoef_l_use      =    cp_para_fft_pkg3d_lg->ncoef_use;
 int   ncoef_l_proc     =    cp_para_fft_pkg3d_lg->ncoef_proc;
 int   nfft_proc        =    cp_para_fft_pkg3d_lg->nfft_proc;
 int   nfft             =    cp_para_fft_pkg3d_lg->nfft;
 int   nfft2            =    nfft/2;
 int   nfft2_proc       =    nfft_proc/2;

 int   nfft_dens_cp_box,nfft2_dens_cp_box;
 int   nfft_proc_dens_cp_box,nfft2_proc_dens_cp_box;
 int   ncoef_l_dens_cp_box;
#ifdef DEBUG_LSDA
 double *rhocr_dens_cp_box;
 double *rhoci_dens_cp_box;
#endif

 double integral,int_tmp;
 int    *recv_counts_coef_dens_cp_box;

 MPI_Comm comm_states   =    communicate->comm_states;

 if(cp_dual_grid_opt >= 1){
   nfft_dens_cp_box       =  cp_para_fft_pkg3d_dens_cp_box->nfft;
   nfft_proc_dens_cp_box  =  cp_para_fft_pkg3d_dens_cp_box->nfft_proc;
   nfft2_dens_cp_box      =  nfft_dens_cp_box/2;
   nfft2_proc_dens_cp_box =  nfft_proc_dens_cp_box/2;


   ncoef_l_dens_cp_box  =  cp_para_fft_pkg3d_dens_cp_box->ncoef;

#ifdef DEBUG_LSDA
   rhocr_dens_cp_box = cpscr->cpscr_rho.rhocr_up_dens_cp_box;
   rhoci_dens_cp_box = cpscr->cpscr_rho.rhoci_up_dens_cp_box;
#endif

   recv_counts_coef_dens_cp_box = 
             cp_para_fft_pkg3d_dens_cp_box->recv_counts_coef;

 }/*endif cp_dual_grid_opt*/

 if(np_states > 1){
   rho_scr = cpscr->cpscr_rho.v_ks_up; 
 } else {
   rho_scr = rho;
 }

/* ================================================================= */
/*0) Check the form of the coefficients                              */

  if(icoef_orth!=1){
    printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
    printf("The coefficients must be in orthogonal form    \n");
    printf("on state processor %d in cp_rho_calc_hybrid   \n",myid_state);
    printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
    fflush(stdout);exit(1);
  }/*endif*/

  if(np_states>1 && icoef_form == 1){
    printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
    printf("The coefficients must be in normal (not transposed) \n");
    printf("form on state processor %d in cp_rho_calc_hybrid  \n",myid_state);
    printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
    fflush(stdout);exit(1);
  }/*endif*/


/* ================================================================= */
/*1) zero density and gradients if necessary                         */

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
/*2)sum the density in real space two states at a time                      */
/*  This is done at state level and uses the scalar packages!               */

 iupper = nstate;
 if((nstate % 2) != 0){ 
     iupper = nstate-1; 
 }/* endif */

 for(is = 1; is <= iupper; is = is + 2){
    ioff   = (is-1)*ncoef;
    ioff2 = (is)*ncoef;

/*--------------------------------------------------------------------------*/
/*I) pack the complex zfft array with two wavefunctions (real) 
    the wavefunctions are reperesented in spherically cuttof 
    half g space                                                            */

    dble_pack_coef(&creal[ioff],&cimag[ioff],&creal[ioff2],&cimag[ioff2],
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
    sngl_pack_coef(&creal[ioff],&cimag[ioff],zfft,cp_sclr_fft_pkg3d_sm);

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

/*=========================================================================*/
/*=========================================================================*/
/*  3) get density in g space                                              */
/*  I)  pack it up                                                         */
/*  rho_scr = rho in scalar                                                */

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
    sngl_upack_coef(rhocr,rhoci,zfft,cp_sclr_fft_pkg3d_lg);
  }else{
    sngl_upack_coef(zfft_tmp,&zfft_tmp[ncoef_l],zfft,cp_sclr_fft_pkg3d_lg);
  }/*endif*/
 }else{
  if(np_states == 1){
    sngl_upack_coef(rhocr_dens_cp_box,rhoci_dens_cp_box,zfft,
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
      Reduce_scatter(&zfft_tmp[1],&rhocr[1],&recv_counts_coef[1],
                     MPI_DOUBLE,MPI_SUM,comm_states);
      Reduce_scatter(&zfft_tmp[(ncoef_l+1)],&rhoci[1],&recv_counts_coef[1],
                     MPI_DOUBLE,MPI_SUM,comm_states);  

      sngl_pack_coef(rhocr,rhoci,zfft,cp_para_fft_pkg3d_lg);

      para_fft_gen3d_fwd_to_r(zfft,zfft_tmp,cp_para_fft_pkg3d_lg);

      sngl_upack_rho(zfft,rho,cp_para_fft_pkg3d_lg); 

   Barrier(comm_states);
   }else{
     Reduce_scatter(&zfft_tmp[1],&rhocr_dens_cp_box[1],
                    &recv_counts_coef_dens_cp_box[1],
                    MPI_DOUBLE,MPI_SUM,comm_states);
     Reduce_scatter(&zfft_tmp[(ncoef_l_dens_cp_box+1)],
                     &rhoci_dens_cp_box[1],&recv_counts_coef_dens_cp_box[1],
                     MPI_DOUBLE,MPI_SUM,comm_states);  

      sngl_pack_coef(rhocr_dens_cp_box,rhoci_dens_cp_box,zfft,
                     cp_para_fft_pkg3d_dens_cp_box);

      para_fft_gen3d_fwd_to_r(zfft,zfft_tmp,cp_para_fft_pkg3d_dens_cp_box);

      sngl_upack_rho(zfft,rho,cp_para_fft_pkg3d_dens_cp_box); 
   }/*endif cp_dual_grid_opt*/
 }/*endif np_states*/

/*===========================================================================*/
/* IF DUALED put rho real space onto the large grid and fft it to g space    */

 if(cp_dual_grid_opt >= 1){
/* sending density*vol_cp on small grid  */
     control_spread_rho(cpscr,rho,cell,dbox_rat,np_states,
                        n_interp_pme_dual,
                        cp_para_fft_pkg3d_dens_cp_box,
                        cp_para_fft_pkg3d_lg,cp_dual_grid_opt);  

  if(np_states == 1){
    para_fft_gen3d_bck_to_g(zfft,zfft_tmp,cp_sclr_fft_pkg3d_lg); 
    sngl_upack_coef(rhocr,rhoci,zfft,cp_sclr_fft_pkg3d_lg);
  }else{
    para_fft_gen3d_bck_to_g(zfft,zfft_tmp,cp_para_fft_pkg3d_lg); 
    sngl_upack_coef(rhocr,rhoci,zfft,cp_para_fft_pkg3d_lg);
  }/*endif np_states*/
 }/*endif cp_dual_grid_opt*/


/*--------------------------------------------------------------------*/
/*  Post-processing for pme grid need to multiply by complex weight factor */
 
  if((n_interp_pme_dual > 1) && (cp_dual_grid_opt == 2)){
    for(i=1; i<= ncoef_l_use; i++){
     temp_r   =  (rhocr[i]*bw_r[i] - rhoci[i]*bw_i[i]); 
     temp_i   =  (rhocr[i]*bw_i[i] + rhoci[i]*bw_r[i]); 
     rhocr[i] =  temp_r;
     rhoci[i] =  temp_i;
    }/*endfor*/

    if((myid_state+1) == np_states){rhocr[ncoef_l_proc]*=bw_r[ncoef_l_proc];}

   }/*endif pme grid */

/*===========================================================================*/
/* IV) finish the density in real space by dividing by the volume           */
/* DUALED SYSTEMS only keep the real space density on the cp_grid            */

     vol_cp  = getdeth(hmat_cp);
     rvol_cp = 1.0/vol_cp;

     if(cp_dual_grid_opt >= 1){
      for(i=1 ; i<= nfft2_proc_dens_cp_box;i++){
         rho[i] *= rvol_cp;
      }/*endfor*/
     }else{
      for(i=1 ; i<= nfft2_proc;i++){
         rho[i] *= rvol_cp;
      }/*endfor*/
     }/*endif cp_dual_grid_opt*/

/*==============================================================*/
/* VII) if doing gradient corrections, get gradient of density    */

  if((cp_gga == 1 || cp_elf_calc_frq > 0)) {
   if(cp_dual_grid_opt >= 1){
    control_grad_rho(cpewald,cpscr,ewald,rhocr_dens_cp_box,rhoci_dens_cp_box,
                     del_rho_x,del_rho_y,del_rho_z,del2_rho,
                     hmati_cp,vol_cp,laplacian_on,cp_dual_grid_opt,
                     cp_para_fft_pkg3d_dens_cp_box);
   }else{
    control_grad_rho(cpewald,cpscr,ewald,rhocr,rhoci,
                     del_rho_x,del_rho_y,del_rho_z,del2_rho,
                     hmati_cp,vol_cp,laplacian_on,cp_dual_grid_opt,
                     cp_para_fft_pkg3d_lg);
   }/*endif cp_dual_grid_opt*/
  }/*endif*/


/*==============================================================*/
}/*end routine*/
/*==============================================================*/




/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
/*==========================================================================*/

void cp_rho_calc_full_g(CPEWALD *cpewald,CPSCR *cpscr,
                        CPCOEFFS_INFO *cpcoeffs_info,EWALD *ewald,
                        CELL *cell,double *creal, double *cimag,
                        int icoef_form,int icoef_orth,
                        double *rhocr ,double *rhoci,double *rho,
                        double *rhocr_dens_cp_box,double *rhoci_dens_cp_box,
                        double *del_rho_x, double *del_rho_y,
                        double *del_rho_z,  
                        double *del2_rho,int nstate,int ncoef,
                        int cp_gga,int cp_dual_grid_opt,
                        int n_interp_pme_dual,
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
 int is,i,iupper;
 double vol_cp,rvol_cp;
 double temp_r,temp_i;

/*  Assign local pointers                                           */
 int cp_elf_calc_frq    =    cpcoeffs_info->cp_elf_calc_frq; 
 double *zfft           =    cpscr->cpscr_wave.zfft;
 double *zfft_tmp       =    cpscr->cpscr_wave.zfft_tmp;
 double *hmati_cp       =    cell->hmati_cp;
 double *hmat_cp        =    cell->hmat_cp;

 double dbox_rat        =    cpewald->dbox_rat;
 double *bw_r           =    cpscr->cpscr_dual_pme.bw_r;
 double *bw_i           =    cpscr->cpscr_dual_pme.bw_i;

 int   myid_state       =    communicate->myid_state;
 int   np_states        =    communicate->np_states;
 int   laplacian_on     =    cpcoeffs_info->cp_laplacian_on;

 int   ncoef_l          =    cp_para_fft_pkg3d_lg->ncoef;
 int   ncoef_l_use      =    cp_para_fft_pkg3d_lg->ncoef_use;
 int   ncoef_l_proc     =    cp_para_fft_pkg3d_lg->ncoef_proc;
 int   nfft_proc        =    cp_para_fft_pkg3d_lg->nfft_proc;
 int   nfft             =    cp_para_fft_pkg3d_lg->nfft;

 int   nfft2_proc       =    nfft_proc/2;
 int   nfft2            =    nfft/2;

 int   nfft_proc_dens_cp_box,nfft2_proc_dens_cp_box;

#ifdef WRITE_DENSITY
#include "../proto_defs/proto_friend_lib_entry.h"
  int nkf1,nkf2,nkf3;  
  int index,ka,kb,kc;
  FILE *fp_rho; 
#endif

#ifdef DEBUG_LSDA
 double *rhocr_dens_cp_box;
 double *rhoci_dens_cp_box;
#endif

 double integral,int_tmp;

 MPI_Comm comm_states   =    communicate->comm_states;

 if(cp_dual_grid_opt >= 1){
   nfft_proc_dens_cp_box   =  cp_para_fft_pkg3d_dens_cp_box->nfft_proc;
   nfft2_proc_dens_cp_box  =  nfft_proc_dens_cp_box/2;

#ifdef DEBUG_LSDA
  if((cp_gga == 1)|| (cp_dual_grid_opt > 1)){
     /*holds the density in g space on cp_box grid*/
   rhocr_dens_cp_box    = cpscr->cpscr_rho.rhocr_up_dens_cp_box;
   rhoci_dens_cp_box    = cpscr->cpscr_rho.rhoci_up_dens_cp_box;
  }/* endif cp_gga */
#endif


 }/*endif cp_dual_grid_opt*/


/* ================================================================= */
/*0) Check the form of the coefficients                              */

  if(icoef_orth!=1){
    printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
    printf("The coefficients must be in orthogonal form    \n");
    printf("on state processor %d in cp_rho_calc_full_g   \n",myid_state);
    printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
    fflush(stdout);exit(1);
  }/*endif*/

  if(np_states>1 && icoef_form == 0){
    printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
    printf("The coefficients must be in transposed (not normal) \n");
    printf("form on state processor %d in cp_rho_calc_full_g  \n",myid_state);
    printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
    fflush(stdout);exit(1);
  }/*endif*/

/* ================================================================= */
/*1) zero density and gradients if necessary                         */

 if(cp_dual_grid_opt >= 1){
  for(i=1; i<= nfft2_proc_dens_cp_box ; i++){
      rho[i] = 0.0;
  }/*endfor*/
 }else{
  for(i=1; i<= nfft2_proc ; i++){
      rho[i] = 0.0;
  }/*endfor*/
 }/*endif cp_dual_grid_opt*/


/*==========================================================================*/
/*==========================================================================*/
/*2)sum the density in real space two states at a time                      */
/*  This is done at state level and uses the scalar packages!               */

 iupper = nstate;
 if((nstate % 2) != 0){ 
     iupper = nstate-1; 
 }/* endif */

 for(is = 1; is <= iupper; is = is + 2){
    ioff   = (is-1)*ncoef;
    ioff2 = (is)*ncoef;

/*--------------------------------------------------------------------------*/
/*I) pack the complex zfft array with two wavefunctions (real) 
    the wavefunctions are reperesented in spherically cuttof 
    half g space                       */

    dble_pack_coef(&creal[ioff],&cimag[ioff],&creal[ioff2],&cimag[ioff2],
                   zfft,cp_para_fft_pkg3d_sm);

/*--------------------------------------------------------------------------*/
/*II) fourier transform the two wavefunctions to real space
     convention exp(-igr)  */
  
    para_fft_gen3d_fwd_to_r(zfft,zfft_tmp,cp_para_fft_pkg3d_sm);

/*--------------------------------------------------------------------------*/
/* III) add the square of the two wave functions to the density  */
    if(cp_dual_grid_opt >= 1){
     sum_rho(zfft,rho,cp_para_fft_pkg3d_dens_cp_box);
    }else{
     sum_rho(zfft,rho,cp_para_fft_pkg3d_lg);
    }/*endif cp_dual_grid_opt */
 }/*endfor is*/


/*--------------------------------------------------------------------------*/
/*IV) if there is an odd number of states, use single pack
      for the last state (i caen't be chaengin' the laws of physics
      captn, i've got to be usin the single pack!!!!!!)  */

 if((nstate % 2 ) != 0) {
    ioff = (nstate-1)*ncoef;
    sngl_pack_coef(&creal[ioff],&cimag[ioff],zfft,cp_para_fft_pkg3d_sm);

/*--------------------------------------------------------------------------*/
/* V) fourier transform the last wavefunction to real space
       convention exp(-igr) */

    para_fft_gen3d_fwd_to_r(zfft,zfft_tmp,cp_para_fft_pkg3d_sm);

/*--------------------------------------------------------------------------*/
/*VI) add the square of the last wave function to the density  */

    if(cp_dual_grid_opt >= 1){
     sum_rho(zfft,rho,cp_para_fft_pkg3d_dens_cp_box);
    }else{ 
     sum_rho(zfft,rho,cp_para_fft_pkg3d_lg);
    }/*endif cp_dual_grid_opt*/
 }/*endif*/

/*--------------------------------------------------------------------------*/
/*rlh Write out the density so that it can be plotted in a graphics program */
/* such as Mathematica. Note: the density includes an extra volume factor   */

#ifdef WRITE_DENSITY
  if(myid_state==0){fp_rho = cfopen("sim_rho_r.dat","o");}

  nkf1    = cp_para_fft_pkg3d_lg->nkf1;
  nkf2    = cp_para_fft_pkg3d_lg->nkf2;
  nkf3    = cp_para_fft_pkg3d_lg->nkf3;

/*--------------------------------------------------------------------------*/
/* write out the density                                                    */
  if(np_states == 1){
    for(kc=1;kc<=nkf3;kc++){
      for(kb=1;kb<=nkf2;kb++){
        for(ka=1;ka<=nkf1;ka++){
          i = (ka-1) + (kb-1)*nkf1 + (kc-1)*nkf1*nkf2 + 1;
          fprintf(fp_rho,"%.5e\n",rho[i]);
        }/* endfor */
      }/* endfor */
    }/* endfor */
  }
  if(myid_state==0){fclose(fp_rho);}
    
#endif /*WRITE_DENSITY*/

/*=========================================================================*/
/*=========================================================================*/
/*  3) get density in g space                                              */
/*  I)  pack it up                                                         */

 if(cp_dual_grid_opt >= 1){
/* sending density*vol_cp on small grid */
  control_spread_rho(cpscr,rho,cell,dbox_rat,np_states,n_interp_pme_dual,
                     cp_para_fft_pkg3d_dens_cp_box,
                     cp_para_fft_pkg3d_lg,cp_dual_grid_opt);  
 }else{
  sngl_pack_rho(zfft,rho,cp_para_fft_pkg3d_lg);
 }/*endif cp_dual_grid_opt*/


/*--------------------------------------------------------------------------*/
/*  II) back transform to g-space  convention exp(igr)   */

  para_fft_gen3d_bck_to_g(zfft,zfft_tmp,cp_para_fft_pkg3d_lg);

/*--------------------------------------------------------------------------*/
/*  III) unpack the density                                    */ 

   sngl_upack_coef(rhocr,rhoci,zfft,cp_para_fft_pkg3d_lg); 

/*--------------------------------------------------------------------*/
/*  Post-processing for pme grid need to multiply by complex weight factor */

  if((n_interp_pme_dual > 1) && (cp_dual_grid_opt == 2)){
    for(i=1; i<= ncoef_l_use; i++){
      temp_r   =  (rhocr[i]*bw_r[i] - rhoci[i]*bw_i[i]);
      temp_i   =  (rhocr[i]*bw_i[i] + rhoci[i]*bw_r[i]);

      rhocr[i] =  temp_r;
      rhoci[i] =  temp_i;
    }/*endfor*/

    if((myid_state+1) == np_states){rhocr[ncoef_l_proc]*=bw_r[ncoef_l_proc];}
  }

/*=========================================================================*/
/* IF CP_DUAL_GRID_OPT and doing GRADIENT CORRECTIONS you need the density in*/
/*   g-space on the mid size grid corresponding to the CP_BOX              */
/*  3.5) get density in G SPACE for CP_BOX                                 */

 if((cp_dual_grid_opt == 1 && cp_gga == 1) || (cp_dual_grid_opt > 1)){
/*  I)  pack it up                                                         */

  sngl_pack_rho(zfft,rho,cp_para_fft_pkg3d_dens_cp_box);

/*--------------------------------------------------------------------------*/
/*  II) back transform to g-space  convention exp(igr)   */

  para_fft_gen3d_bck_to_g(zfft,zfft_tmp,cp_para_fft_pkg3d_dens_cp_box);

/*--------------------------------------------------------------------------*/
/*  III) unpack the density                                                 */ 

  sngl_upack_coef(rhocr_dens_cp_box,rhoci_dens_cp_box,
                  zfft,cp_para_fft_pkg3d_dens_cp_box); 

 }/*endif cp_dual_grid_opt and cp_gga*/

/*--------------------------------------------------------------------------*/
/* IV) finish the density in real space by dividing by the volume */

     vol_cp  = getdeth(hmat_cp);
     rvol_cp = 1.0/vol_cp;
    if(cp_dual_grid_opt >= 1){
     for(i=1 ; i<= nfft2_proc_dens_cp_box;i++){
        rho[i] *= rvol_cp;
     }/*endfor*/
    }else{
     for(i=1 ; i<= nfft2_proc;i++){
        rho[i] *= rvol_cp;
     }/*endfor*/
    }/*endif cp_dual_grid_opt*/

 Barrier(comm_states);
/*==============================================================*/
/* VII) if doing gradient corrections, get gradient of density    */

  if((cp_gga == 1 || cp_elf_calc_frq > 0)) {
   if(cp_dual_grid_opt >= 1){
    control_grad_rho(cpewald,cpscr,ewald,rhocr_dens_cp_box,rhoci_dens_cp_box,
                     del_rho_x,del_rho_y,del_rho_z,del2_rho,
                     hmati_cp,vol_cp,laplacian_on,cp_dual_grid_opt,
                     cp_para_fft_pkg3d_dens_cp_box);
   }else{
    control_grad_rho(cpewald,cpscr,ewald,rhocr,rhoci,
                     del_rho_x,del_rho_y,del_rho_z,del2_rho,
                     hmati_cp,vol_cp,laplacian_on,cp_dual_grid_opt,
                     cp_para_fft_pkg3d_lg);
   }/*endif cp_dual_grid_opt*/
  }/*endif*/




/*==============================================================*/
}/*end routine*/
/*==============================================================*/





/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void coef_force_control(CPOPTS *cpopts,CPCOEFFS_INFO *cpcoeffs_info,
                        CPCOEFFS_POS *cpcoeffs_pos,
                        CPSCR *cpscr,EWALD *ewald,CPEWALD *cpewald,
                        CELL *cell,STAT_AVG *stat_avg,
                        char *vxc_typ,double *pvten_cp,
                        double gc_cut,double alpha_conv_dual,
                        int n_interp_pme_dual,int cp_min_on,
                        COMMUNICATE *communicate,
                        CP_COMM_STATE_PKG *cp_comm_state_pkg_up,
                        CP_COMM_STATE_PKG *cp_comm_state_pkg_dn,
                        PARA_FFT_PKG3D *cp_para_fft_pkg3d_lg,
                        PARA_FFT_PKG3D *cp_sclr_fft_pkg3d_lg,
                        PARA_FFT_PKG3D *cp_para_fft_pkg3d_dens_cp_box,
                        PARA_FFT_PKG3D *cp_sclr_fft_pkg3d_dens_cp_box,
                        PARA_FFT_PKG3D *cp_para_fft_pkg3d_sm,
                        PARA_FFT_PKG3D *cp_sclr_fft_pkg3d_sm,
                        int cp_dual_grid_opt)

/*==========================================================================*/
/*         Begin Routine                                                    */
   {/*Begin Routine*/
/*==========================================================================*/
/* Assign local pointers                                                    */
#include "../typ_defs/typ_mask.h" 
   double cpu1,cpu2;
   int      cp_lda            =  cpopts->cp_lda; 
   int      cp_lsda           =  cpopts->cp_lsda; 
   int      cp_sic            =  cpopts->cp_sic;
   int      cp_nonint         =  cpopts->cp_nonint;
   int      cp_ptens_calc     =  cpopts->cp_ptens_calc;
   int      cp_gga            =  cpopts->cp_gga;
   int      cp_para_opt       =  cpopts->cp_para_opt;

   int      nstate_up         =  cpcoeffs_info->nstate_up_proc;
   int      nstate_dn         =  cpcoeffs_info->nstate_dn_proc;
   int      nstate_up_tot     =  cpcoeffs_info->nstate_up;
   int      nstate_dn_tot     =  cpcoeffs_info->nstate_dn;
   int      laplacian_on      =  cpcoeffs_info->cp_laplacian_on;
   int      cp_tau_functional =  cpcoeffs_info->cp_tau_functional;

   double   *creal_up         =  cpcoeffs_pos->cre_up;
   double   *cimag_up         =  cpcoeffs_pos->cim_up;
   double   *cimag_dn         =  cpcoeffs_pos->cim_dn;
   double   *creal_dn         =  cpcoeffs_pos->cre_dn;
   double   *fcreal_up        =  cpcoeffs_pos->fcre_up;
   double   *fcimag_up        =  cpcoeffs_pos->fcim_up;
   double   *fcimag_dn        =  cpcoeffs_pos->fcim_dn;
   double   *fcreal_dn        =  cpcoeffs_pos->fcre_dn;
   double   *cp_hess_re_up    =  cpcoeffs_pos->cp_hess_re_up;
   double   *cp_hess_im_up    =  cpcoeffs_pos->cp_hess_im_up;
   double   *cp_hess_re_dn    =  cpcoeffs_pos->cp_hess_re_dn;
   double   *cp_hess_im_dn    =  cpcoeffs_pos->cp_hess_im_dn;

   double   *ak2_sm           =  cpewald->ak2_sm;
   double   *v_ks_up          =  cpscr->cpscr_rho.v_ks_up;
   double   *v_ks_dn          =  cpscr->cpscr_rho.v_ks_dn;
   double   *v_ks_tau_up      =  cpscr->cpscr_rho.v_ks_tau_up;
   double   *v_ks_tau_dn      =  cpscr->cpscr_rho.v_ks_tau_dn;
   double   *zfft             =  cpscr->cpscr_wave.zfft;
   double   *zfft_tmp         =  cpscr->cpscr_wave.zfft_tmp;
   double   *cre_scr          =  cpscr->cpscr_wave.cre_up;
   double   *cim_scr          =  cpscr->cpscr_wave.cim_up;

   double   *hmati_cp         =  cell->hmati_cp;

   double   *cp_eke_ret       =   &(stat_avg->cp_eke);
   double   *ks_offset        =   &(cpcoeffs_pos->ks_offset);
   int   myid_state           =  communicate->myid_state;
   int   np_states            =  communicate->np_states;
   int nstate_ncoef_proc_max_up = cpcoeffs_info->nstate_ncoef_proc_max_up;
   int nstate_ncoef_proc_max_dn = cpcoeffs_info->nstate_ncoef_proc_max_dn;
   int nstate_ncoef_proc_up     = cpcoeffs_info->nstate_ncoef_proc_up;
   int nstate_ncoef_proc_dn     = cpcoeffs_info->nstate_ncoef_proc_dn;

   int icoef_orth_up          =  cpcoeffs_pos->icoef_orth_up;
   int icoef_form_up          =  cpcoeffs_pos->icoef_form_up;
   int ifcoef_form_up         =  cpcoeffs_pos->ifcoef_form_up;

   int icoef_orth_dn          =  cpcoeffs_pos->icoef_orth_dn;
   int icoef_form_dn          =  cpcoeffs_pos->icoef_form_dn;
   int ifcoef_form_dn         =  cpcoeffs_pos->ifcoef_form_dn;

/* local variables*/
      int i,iii;
      double cp_eke_dn,cp_eke;

   MPI_Comm comm_states = communicate->comm_states;
   MPI_Comm world       = communicate->world;


/* ================================================================= */
/*0) Check the form of the coefficients                              */

  if(icoef_orth_up!=1){
    printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
    printf("The UP coefficients must be in orthogonal form    \n");
    printf("on state processor %d in coef_force_control\n",myid_state);
    printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
    fflush(stdout);exit(1);
  }/*endif*/
  if((cp_lsda==1) && (nstate_dn > 0)){
   if(icoef_orth_dn!=1){
    printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
    printf("The dn coefficients must be in orthogonal form    \n");
    printf("on state processor %d in coef_force_control \n",myid_state);
    printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
    fflush(stdout);exit(1);
   }/*endif*/
  }/*endif*/

  if(np_states>1 && cp_para_opt == 0) {/* hybrid */
   if( (ifcoef_form_up+icoef_form_up) != 0){
    printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
    printf("The Up coefs/coef forces must be in normal (not transposed) \n");
    printf("form on state processor %d in coef_force_control  \n",myid_state);
    printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
    fflush(stdout);exit(1);
   }/*endif*/
   if((cp_lsda==1) && (nstate_dn > 0) ){
    if( (ifcoef_form_dn+icoef_form_dn) != 0){
     printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
     printf("The Dn coefs/coef forces must be in normal (not transposed) \n");
     printf("form on state processor %d in coef_force_control  \n",myid_state);
     printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
     fflush(stdout);exit(1);
    }/*endif*/
   }/*endif*/
  }/*endif*/

  if(np_states>1 && cp_para_opt == 1) {/* full g */
   if( (ifcoef_form_up+icoef_form_up) != 2){
    printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
    printf("The Up coefs/coef forces must be in transposed (not normal) \n");
    printf("form on state processor %d in coef_force_control  \n",myid_state);
    printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
    fflush(stdout);exit(1);
   }/*endif*/

   if((cp_lsda==1) && (nstate_dn > 0) ){
    if( (ifcoef_form_dn+icoef_form_dn) != 2){
     printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
     printf("The Dn coefs/coef forces must be in transposed (not normal) \n");
     printf("form on state processor %d in coef_force_control  \n",myid_state);
     printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
     fflush(stdout);exit(1);
    }/*endif*/
   }/*endif*/
  }/*endif*/

/*==========================================================================*/
/* II) get the ks potential                                                 */

   cp_get_vks(cpopts,cpscr,cpewald,ewald,communicate,
              cp_comm_state_pkg_up,cp_comm_state_pkg_dn,
              stat_avg,ks_offset,cell,vxc_typ,pvten_cp,
              gc_cut,alpha_conv_dual,n_interp_pme_dual,cp_gga,
              cp_ptens_calc,nstate_dn_tot,
              cp_tau_functional,laplacian_on,cp_sclr_fft_pkg3d_dens_cp_box,
              cp_para_fft_pkg3d_dens_cp_box,cp_sclr_fft_pkg3d_lg,
              cp_para_fft_pkg3d_lg,cp_dual_grid_opt); 

/*==========================================================================*/
/* III) get the force on the states (up and down)                           */

 switch(cp_para_opt){

/*-------------------------------------------------------------------------*/

  case 0: /* hybrid */
 /*-----------------------------------------*/
 /* i)  Up states                           */
   coef_force_calc_hybrid(cpewald,nstate_up,creal_up,cimag_up,
                          fcreal_up,fcimag_up,cre_scr,cim_scr,cp_hess_re_up,cp_hess_im_up,
                          zfft,zfft_tmp,v_ks_up,v_ks_tau_up,ak2_sm,&cp_eke,pvten_cp,
                          cp_ptens_calc,hmati_cp,communicate,icoef_form_up,
                          icoef_orth_up,ifcoef_form_up,cp_tau_functional,cp_min_on,
                          cp_sclr_fft_pkg3d_sm); 
  *cp_eke_ret += cp_eke;

 /*--------------------------------------------*/
 /* ii) down states (if necessary)             */

  if(cp_lsda == 1 && nstate_dn != 0){
   coef_force_calc_hybrid(cpewald,nstate_dn,creal_dn,cimag_dn,
                          fcreal_dn,fcimag_dn,cre_scr,cim_scr,cp_hess_re_dn,cp_hess_im_dn,
                          zfft,zfft_tmp,v_ks_dn,v_ks_tau_dn,ak2_sm,&cp_eke_dn,pvten_cp,
                          cp_ptens_calc,hmati_cp,communicate,icoef_form_dn,
                          icoef_orth_dn,ifcoef_form_dn,cp_tau_functional,cp_min_on,
                          cp_sclr_fft_pkg3d_sm); 
     *cp_eke_ret += cp_eke_dn;
   }/*endif*/

   break;
/*-------------------------------------------------------------------------*/

  case 1: /* full g */
 /*-----------------------------------------*/
 /* i)  Up states                           */

   coef_force_calc_full_g(cpewald,nstate_up_tot,nstate_ncoef_proc_up,
                          nstate_ncoef_proc_max_up,
                          creal_up,cimag_up,fcreal_up,fcimag_up,cre_scr,cim_scr,
                          cp_hess_re_up,cp_hess_im_up,
                          zfft,zfft_tmp,v_ks_up,v_ks_tau_up,ak2_sm,&cp_eke,pvten_cp,
                          cp_ptens_calc,hmati_cp,communicate,icoef_form_up,
                          icoef_orth_up,ifcoef_form_up,cp_tau_functional,cp_min_on,
                          cp_para_fft_pkg3d_sm); 
  *cp_eke_ret += cp_eke;

 /*--------------------------------------------*/
 /* ii) down states (if necessary)             */

  if(cp_lsda == 1 && nstate_dn != 0){
   coef_force_calc_full_g(cpewald,nstate_dn_tot,nstate_ncoef_proc_dn,
                          nstate_ncoef_proc_max_dn,
                          creal_dn,cimag_dn,fcreal_dn,fcimag_dn,cre_scr,cim_scr,
                          cp_hess_re_dn,cp_hess_im_dn,
                          zfft,zfft_tmp,v_ks_dn,v_ks_tau_dn,ak2_sm,&cp_eke_dn,pvten_cp,
                          cp_ptens_calc,hmati_cp,communicate,icoef_form_dn,
                          icoef_orth_dn,ifcoef_form_dn,cp_tau_functional,cp_min_on,
                          cp_para_fft_pkg3d_sm); 
     *cp_eke_ret += cp_eke_dn;
   }/*endif*/

 }/* end switch cp_para_opt */
   
/*======================================================================*/
    }/*end routine*/
/*======================================================================*/




/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void cp_get_vks(CPOPTS *cpopts,CPSCR *cpscr,CPEWALD *cpewald,EWALD *ewald,
                COMMUNICATE *communicate,
                CP_COMM_STATE_PKG *cp_comm_state_pkg_up,
                CP_COMM_STATE_PKG *cp_comm_state_pkg_dn,
                STAT_AVG *stat_avg,double *ks_offset,CELL *cell, char *vxc_typ,
                double *pvten_cp,double gc_cut,double alpha_conv_dual,
                int n_interp_pme_dual,
                int cp_gga,int cp_ptens_calc,int nstate_dn,
                int cp_tau_functional,int laplacian_on,
                PARA_FFT_PKG3D *cp_sclr_fft_pkg3d_dens_cp_box,
                PARA_FFT_PKG3D *cp_para_fft_pkg3d_dens_cp_box,
                PARA_FFT_PKG3D *cp_sclr_fft_pkg3d_lg,
                PARA_FFT_PKG3D *cp_para_fft_pkg3d_lg,
                int cp_dual_grid_opt)

/*==========================================================================*/
/*         Begin Routine                                                    */
{/*Begin Routine*/
/*=======================================================================*/
#include "../typ_defs/typ_mask.h"

/* assign local pointers */
  
   double cpu1,cpu2;

   int       cp_lda          =    cpopts->cp_lda; 
   int       cp_lsda         =    cpopts->cp_lsda; 
   int       cp_nonint       =    cpopts->cp_nonint;
   int       cp_para_opt     =    cpopts->cp_para_opt;
   double    *rhocr          =    cpscr->cpscr_rho.rhocr_up;
   double    *rhoci          =    cpscr->cpscr_rho.rhoci_up;
   double    *rhocr_dens_cp_box =    cpscr->cpscr_rho.rhocr_up_dens_cp_box;
   double    *rhoci_dens_cp_box =    cpscr->cpscr_rho.rhoci_up_dens_cp_box;
   double    *rho_up         =    cpscr->cpscr_rho.rho_up;
   double    *rho_dn         =    cpscr->cpscr_rho.rho_dn;

   double    *vextr          =    cpscr->cpscr_loc.vextr;
   double    *vexti          =    cpscr->cpscr_loc.vexti;
   double    *vextr_dens_cp_box =    cpscr->cpscr_loc.vextr_dens_cp_box;
   double    *vexti_dens_cp_box =    cpscr->cpscr_loc.vexti_dens_cp_box;
   double    *dvextr         =    cpscr->cpscr_loc.dvextr;
   double    *dvexti         =    cpscr->cpscr_loc.dvexti;
   double    *ak2            =    cpewald->ak2;
   double    *ak2_cp_box     =    cpewald->ak2_dens_cp_box;
   double    *v_ks_up        =    cpscr->cpscr_rho.v_ks_up;
   double    *v_ks_dn        =    cpscr->cpscr_rho.v_ks_dn;
   double    *v_ks_tau_up    =    cpscr->cpscr_rho.v_ks_tau_up;
   double    *v_ks_tau_dn    =    cpscr->cpscr_rho.v_ks_tau_dn;
   double    *zfft           =    cpscr->cpscr_wave.zfft;
   double    *zfft_tmp       =    cpscr->cpscr_wave.zfft_tmp;
   double    *bw_r           =    cpscr->cpscr_dual_pme.bw_r;
   double    *bw_i           =    cpscr->cpscr_dual_pme.bw_i;
   double    *clus_corr_r    =    ewald->clus_corr_r;
   int       *kastore        =    ewald->kastr;
   int       *kbstore        =    ewald->kbstr;
   int       *kcstore        =    ewald->kcstr;

   int       myid_state      =    communicate->myid_state;
   int       np_states       =    communicate->np_states;
   int       nstate_up       =    cp_comm_state_pkg_up->nstate;
   int       cp_lyp          =    cpopts->cp_lyp;
   int       cp_lypm1        =    cpopts->cp_lypm1;
   int       nfft            =    cp_para_fft_pkg3d_lg->nfft;
   int       nfft_proc       =    cp_para_fft_pkg3d_lg->nfft_proc;
   int       nfft2           =    nfft/2;
   int       nfft2_proc      =    nfft_proc/2;

   double    *hmat_cp        =    cell->hmat_cp;
   double    *hmati          =    cell->hmati;
   double    *hmat           =    cell->hmat;
   double    *eh_ret         =    &(stat_avg->cp_ehart);
   double    *eext_ret       =    &(stat_avg->cp_eext);
   double    *exc_ret        =    &(stat_avg->cp_exc);
   double    *muxc_ret       =    &(stat_avg->cp_muxc);
   MPI_Comm  comm            =    communicate->comm_states;
   int       ncoef_l         =    cp_para_fft_pkg3d_lg->ncoef_proc;
   int       ncoef_l_dens_cp_box =    
                                  cp_para_fft_pkg3d_dens_cp_box->ncoef_proc;
   int       ncoef_l_use     =    cp_para_fft_pkg3d_lg->ncoef_use;
   int       ncoef_l_use_dens_cp_box =
                                  cp_para_fft_pkg3d_dens_cp_box->ncoef_use;
   int       icoef_off       =    cp_para_fft_pkg3d_lg->icoef_off;
   int       iperd           =    cell->iperd;

   int       *recv_counts_rho;
   int       *displs_rho;
   int       nfft_dens_cp_box;
   int       nfft_proc_dens_cp_box;

   int       nfft2_proc_dens_cp_box;
   int       nfft2_proc_send;

/*         Local Variable declarations                                   */
   double ghfact,pi,fpi;
   double aka,akb,akc,xk,yk,zk,tpi;
   double cfact,eh0,vol,vol_cp,eh,eext,exc,muxc;
   double eext_short_dual,eh_short_dual;
   double pre;
   double temp_r,temp_i;
   int i,j,iii,igo;

   if( cp_dual_grid_opt >= 1){
    nfft_dens_cp_box       = cp_para_fft_pkg3d_dens_cp_box->nfft;
    nfft_proc_dens_cp_box  = cp_para_fft_pkg3d_dens_cp_box->nfft_proc;
    nfft2_proc_dens_cp_box = nfft_proc_dens_cp_box/2;
   }/*endif cp_dual_grid_opt */

/*====================================================================*/
/*  I) calculate hartree potential and add it to the                  */
/*     external potential (sum is over spherically cutoff half space) */

/*--------------------------------------------------------------------*/
/* a) Determine constants and zero energies                           */
   pi     = M_PI; fpi = 4.0*pi; tpi = 2.0*pi;
   pre    = (cp_dual_grid_opt == 2 ? 0.25/(alpha_conv_dual*alpha_conv_dual)
                                     : 0.0);
 
   vol    = getdeth(hmat);
   vol_cp = getdeth(hmat_cp);

   eh   = 0.0;
   eext = 0.0;
   exc  = 0.0;
   muxc = 0.0;


   if(cp_nonint==0){

/*--------------------------------------------------------------------*/
/*  c) Get Hartree + external contributions to VKS -- test for CBCs   */

    if(iperd == 3 ){
     for(i=1 ; i<= ncoef_l_use; i++){
                  /*erf  long range piece on large grid for dual opt */
        ghfact    = fpi*exp(-ak2[i]*pre)/(ak2[i]*vol);  
        eh       +=  0.50*ghfact*(rhocr[i]*rhocr[i] + rhoci[i]*rhoci[i]); 
        eext     +=  vextr[i]*rhocr[i] + vexti[i]*rhoci[i];
        vextr[i] +=  ghfact*rhocr[i];
        vexti[i] +=  ghfact*rhoci[i]; 
     }/*endfor*/
    } else {

     for(i=1 ; i<= ncoef_l_use; i++){            
                  /*erf  long range piece on large grid for dual opt */
        ghfact    = (fpi*exp(-ak2[i]*pre)/ak2[i] + clus_corr_r[i])/vol;
        eh       +=  0.50*ghfact*(rhocr[i]*rhocr[i] + rhoci[i]*rhoci[i]);
        eext     +=  vextr[i]*rhocr[i] + vexti[i]*rhoci[i];

        vextr[i] +=  ghfact*rhocr[i];
        vexti[i] +=  ghfact*rhoci[i]; 


     }/*endfor*/

    }/*endif*/

    eext *= 2.0;

    if((myid_state+1) == np_states){eext +=  vextr[ncoef_l]*rhocr[ncoef_l];}
     eh0 = 2.0*(eh);
     eh  *= 2.0;


/*--------------------------------------------------------------------*/
/*  d) Add in convergent part of long range erf for dualing           */

     if( (cp_dual_grid_opt == 2) && ((myid_state+1) == np_states) ){
       eh -= 0.5*fpi*pre*rhocr[ncoef_l]*rhocr[ncoef_l]/vol; 
       vextr[ncoef_l] -= fpi*pre*rhocr[ncoef_l]/vol;    
     }/*endif*/

/*--------------------------------------------------------------------*/
/*  e) Add in cluster correction if necessary */ 

    if( (iperd != 3 )&& ((myid_state+1) == np_states) ) {
       eh += 0.5*clus_corr_r[ncoef_l]*rhocr[ncoef_l]*rhocr[ncoef_l]/vol;
       vextr[ncoef_l] += clus_corr_r[ncoef_l]*rhocr[ncoef_l]/vol;
     }/*endif*/

/*--------------------------------------------------------------------*/
/*  f) get hartree contribution to pressure tensor                   */

    if(cp_ptens_calc == 1){
      for(i=1; i<= ncoef_l_use; i++){
        aka = (double)kastore[(i+icoef_off)];
        akb = (double)kbstore[(i+icoef_off)];
        akc = (double)kcstore[(i+icoef_off)];

        xk = (aka*hmati[1] + akb*hmati[2] + akc*hmati[3])*tpi;
        yk = (aka*hmati[4] + akb*hmati[5] + akc*hmati[6])*tpi;
        zk = (aka*hmati[7] + akb*hmati[8] + akc*hmati[9])*tpi;
       
        cfact= -2.0*fpi*(rhocr[i]*rhocr[i] + rhoci[i]*rhoci[i])/
               (ak2[i]*ak2[i]*vol);

        pvten_cp[1] +=  xk*xk*cfact;
        pvten_cp[2] +=  xk*yk*cfact;
        pvten_cp[3] +=  xk*zk*cfact;
        pvten_cp[4] +=  xk*yk*cfact;
        pvten_cp[5] +=  yk*yk*cfact;
        pvten_cp[6] +=  yk*zk*cfact;
        pvten_cp[7] +=  xk*zk*cfact;
        pvten_cp[8] +=  yk*zk*cfact;
        pvten_cp[9] +=  zk*zk*cfact;
      }/*endfor*/

      pvten_cp[1] +=  eh0;
      pvten_cp[5] +=  eh0;
      pvten_cp[9] +=  eh0;
    }/*endif ptens_calc*/

   }else{  /* if nonint == 1 (electrons do not interact) */
/*--------------------------------------------------------------------*/
/*  g) Only get external contribution to VKS                          */

    for(i=1; i<= ncoef_l_use; i++){
      eext += (vextr[i]*rhocr[i] + vexti[i]*rhoci[i]);
    }/*endfor*/

    eext *= 2.0;

    if((myid_state+1) == np_states) eext += vextr[ncoef_l]*rhocr[ncoef_l];
   }/*endif nonint*/

/*--------------------------------------------------------------------*/
/*  h) get external contribution to pressure tensor                   */

   if(cp_ptens_calc==1){
     pvten_cp[1] +=  eext;
     pvten_cp[5] +=  eext;
     pvten_cp[9] +=  eext;

     for(i=1; i<= ncoef_l_use; i++){
       aka = (double)kastore[(i+icoef_off)];
       akb = (double)kbstore[(i+icoef_off)];
       akc = (double)kcstore[(i+icoef_off)];

       xk = (aka*hmati[1] +akb*hmati[2] +akc*hmati[3])*tpi;
       yk = (aka*hmati[4] +akb*hmati[5] +akc*hmati[6])*tpi;
       zk = (aka*hmati[7] +akb*hmati[8] +akc*hmati[9])*tpi;

       cfact = 2.0*(dvextr[i]*rhocr[i] + dvexti[i]*rhoci[i]);

       pvten_cp[1] +=  xk*xk*cfact;
       pvten_cp[2] +=  xk*yk*cfact;
       pvten_cp[3] +=  xk*zk*cfact;
       pvten_cp[4] +=  xk*yk*cfact;
       pvten_cp[5] +=  yk*yk*cfact;
       pvten_cp[6] +=  yk*zk*cfact;
       pvten_cp[7] +=  xk*zk*cfact;
       pvten_cp[8] +=  yk*zk*cfact;
       pvten_cp[9] +=  zk*zk*cfact;
     }/*endfor*/
   }/*endif ptens_calc*/

/*--------------------------------------------------------------------*/
/* i) Increment the energies */

   *eext_ret += eext;
   *eh_ret += eh;

/*--------------------------------------------------------------------*/
/* j) Multiply the external potential by the pme g-space weight       */

   if((n_interp_pme_dual > 1) && (cp_dual_grid_opt == 2)){

     for(i=1; i<= ncoef_l_use; i++){
       temp_r   = ( vextr[i]*bw_r[i] + vexti[i]*bw_i[i]);
       temp_i   = (-vextr[i]*bw_i[i] + vexti[i]*bw_r[i]); 
       /*       temp_r   = ( vextr[i]*bw_r[i] - vexti[i]*bw_i[i]);
		temp_i   = (-vextr[i]*bw_i[i] - vexti[i]*bw_r[i]);  */
       vextr[i] =  temp_r;
       vexti[i] =  temp_i;
     }/*endfor*/
     if((myid_state+1) == np_states){vextr[ncoef_l]*=bw_r[ncoef_l];}

   }/*endif pme grid */

/*====================================================================*/
/*  II) single pack the ks potential for fourier transform routine    */

   sngl_pack_coef(vextr,vexti,zfft,cp_para_fft_pkg3d_lg);

/*====================================================================*/
/* III) fourier transform ks potential to real space exp(-igr)        */

   para_fft_gen3d_fwd_to_r(zfft,zfft_tmp,cp_para_fft_pkg3d_lg);

/*====================================================================*/
/* IV) Contract and unpack rho for dual option, otherwise just upack  */

   if(cp_dual_grid_opt >= 1){
      control_contract_rho(cpscr,v_ks_up,cell,np_states,n_interp_pme_dual,
                           cp_para_fft_pkg3d_dens_cp_box,
                           cp_para_fft_pkg3d_lg,cp_dual_grid_opt);   
   }else{
      sngl_upack_rho(zfft,v_ks_up,cp_para_fft_pkg3d_lg);
   }/*endif cp_dual_grid_opt*/

/*=====================================================================*/
/* V) Coulomb and Hartree erfc calculated on small grid for PME dualing*/

  if(cp_dual_grid_opt == 2){
      eh_short_dual = 0.0;
    eext_short_dual = 0.0;
/*--------------------------------------------------------------------*/
/* a) Hartree plus external erfc Coloumb                              */

    if(cp_nonint == 0){
      for(i=1 ; i<= ncoef_l_use_dens_cp_box; i++){
        ghfact    = fpi*(1.0 - exp(-ak2_cp_box[i]*pre))/(ak2_cp_box[i]*vol_cp);  

        eh_short_dual += 0.50*ghfact*
                       (rhocr_dens_cp_box[i]*rhocr_dens_cp_box[i]
                      + rhoci_dens_cp_box[i]*rhoci_dens_cp_box[i]);

        eext_short_dual +=  vextr_dens_cp_box[i]*rhocr_dens_cp_box[i]
                        +   vexti_dens_cp_box[i]*rhoci_dens_cp_box[i];

        vextr_dens_cp_box[i] +=  ghfact*rhocr_dens_cp_box[i];
        vexti_dens_cp_box[i] +=  ghfact*rhoci_dens_cp_box[i]; 
      }/*endfor*/

      eext_short_dual *= 2.0;
      eh_short_dual *= 2.0;

      if((myid_state+1) == np_states){

        eext_short_dual +=  vextr_dens_cp_box[ncoef_l_dens_cp_box]
                           *rhocr_dens_cp_box[ncoef_l_dens_cp_box];

         eh_short_dual  += 0.5*fpi*pre*rhocr_dens_cp_box[ncoef_l_dens_cp_box]
                           *rhocr_dens_cp_box[ncoef_l_dens_cp_box]/vol_cp;

         vextr_dens_cp_box[ncoef_l_dens_cp_box] += fpi*pre
                         *rhocr_dens_cp_box[ncoef_l_dens_cp_box]/vol_cp;   
       }/*endif*/

/*--------------------------------------------------------------------*/
/* b) External erfc Coloumb                                           */

    }else{

       for(i=1; i<= ncoef_l_use_dens_cp_box; i++){
         eext_short_dual += (vextr_dens_cp_box[i]*rhocr_dens_cp_box[i] 
                         +   vexti_dens_cp_box[i]*rhoci_dens_cp_box[i]);
       }/*endfor*/

       eext_short_dual *= 2.0;

       if((myid_state+1) == np_states){
          eext_short_dual += vextr_dens_cp_box[ncoef_l_dens_cp_box]
                          *  rhocr_dens_cp_box[ncoef_l_dens_cp_box];
        }/*endif*/
    }/*endif cp_nonint*/

    eext      += eext_short_dual;
    eh        += eh_short_dual;
    *eext_ret += eext_short_dual;
    *eh_ret   += eh_short_dual;

/*--------------------------------------------------------------------*/
/*  c) add contribution to vks                                        */
    sngl_pack_coef(vextr_dens_cp_box,vexti_dens_cp_box,
                  zfft,cp_para_fft_pkg3d_dens_cp_box);

    para_fft_gen3d_fwd_to_r(zfft,zfft_tmp,cp_para_fft_pkg3d_dens_cp_box);

    sngl_upack_rho_sum(zfft,v_ks_up,cp_para_fft_pkg3d_dens_cp_box);

  }/*endif cp_dual_grid_opt*/

/*====================================================================*/
/* IV) Assign up to down for LSDA because vext is vext        */

   if(cp_lsda==1 && nstate_dn != 0){
     if(cp_dual_grid_opt >= 1){
      for(i=1;i<=nfft2_proc_dens_cp_box;i++){v_ks_dn[i]=v_ks_up[i];}
     }else{
      for(i=1;i<=nfft2_proc;i++){v_ks_dn[i]=v_ks_up[i];}
     }/*endif cp_dual_grid_opt*/
   }/*endif*/

/*====================================================================*/
/* V) calculate the exchange potential and add to v_ks if necessary   */
/*    on the small grid only for dualing                              */

  if(cp_nonint==0){
    igo = 0;
/*--------------------------------------------------------------------*/
/*  a) LDA                                                            */

/*--------------*/
/* Perdew-Zunger*/
/*--------------*/
  
    if(cp_lda==1){
     if(strcasecmp(vxc_typ,"pz_lda")==0){
       if(cp_dual_grid_opt >= 1){
         excpot_pz_lda(v_ks_up,rho_up,&exc,&muxc,nfft_dens_cp_box,
                       nfft_proc_dens_cp_box,vol_cp,
                       cp_lyp,cp_lypm1,pvten_cp,cp_ptens_calc); 
       }else{
         excpot_pz_lda(v_ks_up,rho_up,&exc,&muxc,nfft,nfft_proc,vol_cp,
                      cp_lyp,cp_lypm1,pvten_cp,cp_ptens_calc);

       }/*endif cp_dual_grid_opt*/
      igo = 1;
     }/*endif*/

/*------------------*/
/* Perdew-Wang 1992 */
/*------------------*/

     if(strcasecmp(vxc_typ,"pw_lda")==0 ){
       if(cp_dual_grid_opt >= 1){
         excpot_pw_lda(v_ks_up,rho_up,&exc,&muxc,nfft_dens_cp_box,
                       nfft_proc_dens_cp_box,vol_cp,
                       cp_lyp,cp_lypm1,pvten_cp,cp_ptens_calc);
        }else{
          excpot_pw_lda(v_ks_up,rho_up,&exc,&muxc,nfft,nfft_proc,vol_cp,
                        cp_lyp,cp_lypm1,pvten_cp,cp_ptens_calc);
       }/*endif cp_dual_grid_opt */

      igo = 1;
     }/*endif*/

/*------------------*/
/* Pade */
/*------------------*/

     if(strcasecmp(vxc_typ,"pade_lda")==0 ){
       if(cp_dual_grid_opt >= 1){
         excpot_pade_lda(v_ks_up,rho_up,&exc,&muxc,nfft_dens_cp_box,
                         nfft_proc_dens_cp_box,vol_cp,
                         cp_lyp,cp_lypm1,pvten_cp,cp_ptens_calc);
        }else{
          excpot_pade_lda(v_ks_up,rho_up,&exc,&muxc,nfft,nfft_proc,vol_cp,
                          cp_lyp,cp_lypm1,pvten_cp,cp_ptens_calc);
       }/*endif cp_dual_grid_opt */

      igo = 1;
     }/*endif*/


   }/*endif LDA */

/*--------------------------------------------------------------------*/
/*  b) LSDA                                                           */

   if(cp_lsda == 1) {

/*--------------*/
/* Perdew-Zunger*/
/*--------------*/

     if(strcasecmp(vxc_typ,"pz_lsda")==0){
       if(cp_dual_grid_opt >= 1){
        excpot_pz_lsda(v_ks_up,v_ks_dn,rho_up,rho_dn,&exc,&muxc,
                       nfft_dens_cp_box,nfft_proc_dens_cp_box,
                       vol_cp,cp_lyp,cp_lypm1,pvten_cp,cp_ptens_calc); 
       }else{
       excpot_pz_lsda(v_ks_up,v_ks_dn,rho_up,rho_dn,&exc,&muxc,nfft,nfft_proc,
                      vol_cp,cp_lyp,cp_lypm1,pvten_cp,cp_ptens_calc); 
       }/*endif cp_dual_grid_opt*/
       igo = 1;
     }/*endif*/

/*--------------*/
/* Perdew-Wang  */
/*--------------*/

     if(strcasecmp(vxc_typ,"pw_lsda")==0){
       if(cp_dual_grid_opt >= 1){
        excpot_pw_lsda(v_ks_up,v_ks_dn,rho_up,rho_dn,&exc,&muxc,
                       nfft_dens_cp_box,nfft_proc_dens_cp_box,
                       vol_cp,cp_lyp,cp_lypm1,pvten_cp,cp_ptens_calc); 
       }else{
       excpot_pw_lsda(v_ks_up,v_ks_dn,rho_up,rho_dn,&exc,&muxc,nfft,nfft_proc,
                      vol_cp,cp_lyp,cp_lypm1,pvten_cp,cp_ptens_calc); 
       }/*endif cp_dual_grid_opt*/
       igo = 1;
     }/*endif*/

/*---------*/
/* Pade    */
/*---------*/

     if(strcasecmp(vxc_typ,"pade_lsda")==0){
       if(cp_dual_grid_opt >= 1){
        excpot_pade_lsda(v_ks_up,v_ks_dn,rho_up,rho_dn,&exc,&muxc,
                         nfft_dens_cp_box,nfft_proc_dens_cp_box,
                         vol_cp,cp_lyp,cp_lypm1,pvten_cp,cp_ptens_calc); 
       }else{
       excpot_pade_lsda(v_ks_up,v_ks_dn,rho_up,rho_dn,&exc,&muxc,nfft,nfft_proc,
                        vol_cp,cp_lyp,cp_lypm1,pvten_cp,cp_ptens_calc); 
       }/*endif cp_dual_grid_opt*/
       igo = 1;
     }/*endif*/


/*--------------*/

   }/*endif*/
  

/*--------------------------------------------------------------------*/
/*  c) add gradient corrections if necessary                          */

   if(cp_gga==1){
    if(cp_lsda == 0 ){
     if(cp_dual_grid_opt >= 1){
      grad_corr_lda(cpscr,cpewald,ewald,cell,&exc,&muxc,vol_cp,
                    cpopts,gc_cut,pvten_cp,cp_ptens_calc,
                    cp_tau_functional,laplacian_on,
                    cp_dual_grid_opt,cp_para_fft_pkg3d_dens_cp_box,nstate_up);
     }else{
      grad_corr_lda(cpscr,cpewald,ewald,cell,&exc,&muxc,vol_cp,
                    cpopts,gc_cut,pvten_cp,cp_ptens_calc,
                    cp_tau_functional,laplacian_on,
                    cp_dual_grid_opt,cp_para_fft_pkg3d_lg,nstate_up);
     }/*endif cp_dual_grid_opt*/

    }else {
     if(cp_dual_grid_opt >= 1){
      grad_corr_lsda(cpscr,cpewald,ewald,cell,&exc,&muxc,vol_cp,
                     cpopts,gc_cut,pvten_cp,cp_ptens_calc,cp_tau_functional,
                     nstate_up,nstate_dn,laplacian_on,cp_dual_grid_opt,
                     cp_para_fft_pkg3d_dens_cp_box);
     }else{
      grad_corr_lsda(cpscr,cpewald,ewald,cell,&exc,&muxc,vol_cp,
                     cpopts,gc_cut,pvten_cp,cp_ptens_calc,cp_tau_functional,
                     nstate_up,nstate_dn,laplacian_on,cp_dual_grid_opt,
                     cp_para_fft_pkg3d_lg);
     }/*endif cp_dual_grid_opt*/
    }/* endif cp_lsda */
   }/*endif cp_gga*/

/*--------------------------------------------------------------------*/
/*  d Check for implementation of chosen xc functional (barf if not found)*/

   if(igo == 0){
     fprintf(stderr,"@@@@@@@@@@@@@@@@@@-error-@@@@@@@@@@@@@@@@@@@@@@@ \n\n");
     fprintf(stderr,"user specified correlation functional %s\n",vxc_typ);
     fprintf(stderr,"not found in the guts. contact technical support\n");
     fprintf(stderr,"@@@@@@@@@@@@@@@@@@-error-@@@@@@@@@@@@@@@@@@@@@@@\n");
     exit(1);
   }/*endif*/

  }/*endif nonint==0*/

/*--------------------------------------------------------------------*/
/* e) Increment the energies */

  *exc_ret   += exc;
  *muxc_ret  += muxc;
  *ks_offset  = exc-eh-muxc;

/*====================================================================*/
/* Get the KS matrix in form for state level                          */
/*  that is for hybrid parallel, allgather the puppy                  */

  if(np_states>1 && cp_para_opt == 0){/* hybrid parallel only */

    if(cp_dual_grid_opt >= 1){
      nfft2_proc_send = nfft2_proc_dens_cp_box;
      recv_counts_rho = cp_para_fft_pkg3d_dens_cp_box->recv_counts_rho;
      displs_rho      = cp_para_fft_pkg3d_dens_cp_box->displs_rho;
    }else{
      nfft2_proc_send = nfft2_proc;
      recv_counts_rho = cp_para_fft_pkg3d_lg->recv_counts_rho;
      displs_rho      = cp_para_fft_pkg3d_lg->displs_rho;
    }/*endif cp_dual_grid_opt*/

    for(i=1;i<=nfft2_proc_send;i++){zfft[i] = v_ks_up[i];}
     Allgatherv(&(zfft[1]),nfft2_proc_send,MPI_DOUBLE,&(v_ks_up[1]),
                &recv_counts_rho[1],&(displs_rho[1]),MPI_DOUBLE,0,comm);

     if(cp_tau_functional==1){
      for(i=1;i<=nfft2_proc_send;i++){zfft[i] = v_ks_tau_up[i];}
      Allgatherv(&(zfft[1]),nfft2_proc_send,MPI_DOUBLE,&(v_ks_tau_up[1]),
                 &recv_counts_rho[1],&(displs_rho[1]),MPI_DOUBLE,0,comm);
     }/* endif */
    if(cp_lsda==1){
       for(i=1;i<=nfft2_proc_send;i++){zfft[i]=v_ks_dn[i];}
       Allgatherv(&(zfft[1]),nfft2_proc_send,MPI_DOUBLE,&(v_ks_dn[1]),
                  &recv_counts_rho[1],&(displs_rho[1]),MPI_DOUBLE,0,comm);

       if(cp_tau_functional==1){
         for(i=1;i<=nfft2_proc_send;i++){zfft[i]=v_ks_tau_dn[i];}
         Allgatherv(&(zfft[1]),nfft2_proc_send,MPI_DOUBLE,&(v_ks_tau_dn[1]),
                    &recv_counts_rho[1],&(displs_rho[1]),MPI_DOUBLE,0,comm);
       }/* endif */
     }/*endif*/

  }/*endif*/

    
/*==========================================================================*/
   }/* end routine*/
/*==========================================================================*/


/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
/* calculate the force on the coeff's (up or down, both) */
/* given an appropriate v_ks (up or down, both).         */
/*==========================================================================*/

void coef_force_calc_hybrid(CPEWALD *cpewald,int nstate,
                             double *creal,double *cimag, 
                             double *fcreal,double  *fcimag,
                             double *cre_scr,double *cim_scr,
                             double *cp_hess_re,double *cp_hess_im,
                             double *zfft,double *zfft_tmp,
                             double *v_ks,double *v_ks_tau,double *ak2_sm,
                             double *eke_ret,double *pvten_cp,
                             int cp_ptens_calc,double *hmati,
                             COMMUNICATE *communicate,
                             int icoef_form,int icoef_orth,int ifcoef_form,
                             int cp_tau_functional,int cp_min_on,
                             PARA_FFT_PKG3D *cp_sclr_fft_pkg3d_sm)

/*==========================================================================*/
/*         Begin Routine                                                    */
{/*Begin Routine*/
/*=======================================================================*/
/*            Local variable declarations                               */
#include "../typ_defs/typ_mask.h"

   int is,i,iupper;
   double tpi;
   double aka,akb,akc,xk,yk,zk,cfact;
   double eke;
   int ioff,ncoef1,ioff2;
   int iii,iis,nis;

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

      double sum_check,sum_check_tmp;
      MPI_Comm comm_states = communicate->comm_states;

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

  ncoef1 = ncoef - 1;
  iupper = nstate;
  if(nstate % 2 == 1){
     iupper = nstate - 1;
  }

/*=================================================================*/
/*  get the forces on the coefs of each state                      */

  for(is=1 ; is<= iupper; is+=2 ){

    ioff = (is-1)*ncoef;
    ioff2 = (is)*ncoef;

/*==========================================================================*/
/* 1) get the wave functions in real space two at a time                    */
/*   I) double pack the complex zfft array with two real wavefunctions      */

    dble_pack_coef(&creal[ioff],&cimag[ioff],&creal[ioff2],&cimag[ioff2],
                   zfft,cp_sclr_fft_pkg3d_sm);

/*--------------------------------------------------------------------------*/
/* II) fourier transform the wavefunctions to real space                    */
/*      convention exp(-igr)                                                */

    para_fft_gen3d_fwd_to_r(zfft,zfft_tmp,cp_sclr_fft_pkg3d_sm);

/*==========================================================================*/
/* 2) get v|psi> in g space and store it in zfft                            */
/*   I) get  v|psi> in real space                                           */

    cp_vpsi(zfft,v_ks,nfft);  

/*--------------------------------------------------------------------------*/
/*  II) fourier transform  to g-space                                       */
/*     convention exp(igr)                                                  */

    para_fft_gen3d_bck_to_g(zfft,zfft_tmp,cp_sclr_fft_pkg3d_sm);

/*==========================================================================*/
/* 3) get forces on coefficients by double unpacking the array zfft         */

    dble_upack_coef_sum(&fcreal[ioff],&fcimag[ioff],
                        &fcreal[ioff2],&fcimag[ioff2],
                        zfft,cp_sclr_fft_pkg3d_sm);

  }/*endfor is */

/*==========================================================================*/
/*==========================================================================*/
/* 4) if there is an odd number of states, go through                       */
/*      the same procedure using sng_packs                                  */



  if(nstate % 2 != 0){
     is = nstate;
     ioff = (is -1)*ncoef;

/*--------------------------------------------------------------------------*/
/*   I) sngl pack                                                           */

     sngl_pack_coef(&creal[ioff],&cimag[ioff],zfft,cp_sclr_fft_pkg3d_sm);

/*--------------------------------------------------------------------------*/
/* II) fourier transform the wavefunctions to real space                    */
/*      convention exp(-igr)                                                */
 
      para_fft_gen3d_fwd_to_r(zfft,zfft_tmp,cp_sclr_fft_pkg3d_sm);

/*==========================================================================*/
/* 5) get v|psi> in g space and store it in zfft                            */
/*   I) get  v|psi> in real space                                           */

      cp_vpsi(zfft,v_ks,nfft);

/*--------------------------------------------------------------------------*/
/*   II) fourier transform the result back to g-space */
/*     convention exp(igr)  */

      para_fft_gen3d_bck_to_g(zfft,zfft_tmp,cp_sclr_fft_pkg3d_sm);

/*==========================================================================*/
/* 6) get forces on coefficients by single unpacking the array zfft         */


      sngl_upack_coef_sum(&fcreal[ioff],&fcimag[ioff],zfft,
                          cp_sclr_fft_pkg3d_sm);

  }/* endif: odd number of states*/

/*==========================================================================*/
/* 7) If there is an electron KE density dependent functional, calculate    */
/*    this contribution to the force                                        */

  if(cp_tau_functional==1)
    coef_force_tau_fun_hybrid(cpewald,nstate,creal,cimag,fcreal,fcimag,
                              cre_scr,cim_scr,zfft,zfft_tmp,v_ks_tau,ak2_sm,
                              pvten_cp,cp_ptens_calc,hmati,communicate,
                              icoef_form,icoef_orth,ifcoef_form,
                              cp_sclr_fft_pkg3d_sm);

/*==========================================================================*/
/* 8) If doing minimization, Fourier transform the KS potential to g-space  */
/*    and unpack it into the diagonal Hessian                               */

  if(cp_min_on > 0 ){

    for(i=1;i<=ncoef;i++){ cp_hess_re[i] = cp_hess_im[i] = 0.0;}

    sngl_pack_rho(zfft,v_ks,cp_sclr_fft_pkg3d_sm);

    para_fft_gen3d_bck_to_g(zfft,zfft_tmp,cp_sclr_fft_pkg3d_sm);

    sngl_upack_coef(cp_hess_re,cp_hess_im,zfft,cp_sclr_fft_pkg3d_sm);

  }/* endif cp_min_on */


/*==========================================================================*/
/* 9) calculate the kinetic energy term and add its contribution to the force*/

  tpi = 2.0*M_PI;
  eke = 0.0;
  for(is=1 ; is<= nstate ; is++){
    ioff = (is-1)*ncoef;
    for(i=1; i<= ncoef1 ; i++){
      iis = ioff + i;
      fcreal[iis] -= 2.0*ak2_sm[i]*creal[iis];
      fcimag[iis] -= 2.0*ak2_sm[i]*cimag[iis];
      eke += (2.0*ak2_sm[i]*(creal[iis]*creal[iis] + cimag[iis]*cimag[iis]));
    }/*endfor i*/
   nis = is*ncoef;
   fcimag[nis] = 0.0;
  }/*endfor*/

  eke *= .50;
  *eke_ret = eke;

/*================================================================================*/
/* 10) If doing minimization, calculat kinetic contribution to diagonal Hessian   */

  if(cp_min_on){    
    for(i=1; i<= ncoef ; i++){
      cp_hess_re[i] *= -1;
      cp_hess_im[i] *= -1;
    }/* endfor */
    for(i=1; i<= ncoef1 ; i++){
      cp_hess_re[i] += 2.0*ak2_sm[i];
      cp_hess_im[i] += 2.0*ak2_sm[i];
    }/* endfor */
  }/* endif cp_min_on */

/*===================================================================*/
/* Get a debug tool for this system size                             */

#ifdef DEBUG
  if(myid_state==0){
    printf("Would you like me to write out the model wavefunction: (1/0) \n");
    scanf("%d",&iii);
    if(iii==1){
     vol = 1.0/(hmati[1]*hmati[5]*hmati[9]);
     anorm = 4.0*pow(M_PI,0.75)/sqrt(vol);    
     sum = 0.0;
     fp = fopen("model_lda_gga.coef_in","w");
     fprintf(fp,
   "ncoef_up,ncoef_dn,nstate_up,nstate_dn,dft_typ,restart_typ,time of dump\n");
     fprintf(fp,"%d 0 2 2 lda restart_pos 1000 norb_off\n",ncoef);
     fprintf(fp,"occ up occ dn\n");
     fprintf(fp,"1 1 \n");     
     fprintf(fp,"1 1 \n");     
     fprintf(fp,"Up_State_Real         Up_State_Imag\n");
     dx = 0.2/hmati[1];
     for(is=1;is<=2;is++){
     x_pos = (is-1)*dx;
     y_pos = 0.0;
     z_pos = 0.0;
     for(icount = 1; icount<= (ncoef1) ; icount++){
       aka = (double)kastore_sm[icount];
       akb = (double)kbstore_sm[icount];
       akc = (double)kcstore_sm[icount];
       xk = (aka*hmati[1] +akb*hmati[2] +akc*hmati[3])*tpi;
       yk = (aka*hmati[4] +akb*hmati[5] +akc*hmati[6])*tpi;
       zk = (aka*hmati[7] +akb*hmati[8] +akc*hmati[9])*tpi;
       g2 = xk*xk + yk*yk + zk*zk;
       c_g = exp(-g2/2.0)*anorm;
       arg = xk*x_pos + yk*y_pos + zk*z_pos;
       phase_r = cos(arg);
       phase_i = sin(arg);
       cre_now = c_g*phase_r;
       cim_now = c_g*phase_i;
       fprintf(fp,"%g %g\n",cre_now,cim_now);
       sum = sum + 2.0*c_g*c_g;
     }/*endfor*/
     c_g = anorm;
     sum = sum + c_g*c_g;
     fprintf(fp,"%g 0\n",c_g);
     printf("sum = %g\n",sum);
   }/*endfor*/
     fclose(fp);
     Finalize();
     exit(1);
   }/*endif*/
  }/* endif myid */
  Dbx_Barrier(communicate->comm_states);
#endif

/*==========================================================================*/
/* 9) calculate kinetic contribution to pressure tensor                     */

  if(cp_ptens_calc == 1) {
   tpi = 2.0*M_PI;
   for(is=1; is<= nstate; is++){
     ioff = (is-1)*ncoef;
     for(i=1; i<= ncoef1; i++){
       iis = ioff + i;
       aka = (double)kastore_sm[i];
       akb = (double)kbstore_sm[i];
       akc = (double)kcstore_sm[i];

       xk = (aka*hmati[1] +akb*hmati[2] +akc*hmati[3])*tpi;
       yk = (aka*hmati[4] +akb*hmati[5] +akc*hmati[6])*tpi;
       zk = (aka*hmati[7] +akb*hmati[8] +akc*hmati[9])*tpi;

       cfact = 2.0*(creal[iis]*creal[iis] + cimag[iis]*cimag[iis]);

       pvten_cp[1] += xk*xk*cfact;
       pvten_cp[2] += xk*yk*cfact;
       pvten_cp[3] += xk*zk*cfact;

       pvten_cp[4] += xk*yk*cfact;
       pvten_cp[5] += yk*yk*cfact;
       pvten_cp[6] += yk*zk*cfact;

       pvten_cp[7] += xk*zk*cfact;
       pvten_cp[8] += yk*zk*cfact;
       pvten_cp[9] += zk*zk*cfact;

     }/*endfor*/
   }/*endfor*/
  }/*endif */

/* fine fertig terminado finito */
/*==========================================================================*/
   }/*end routine*/
/*==========================================================================*/





/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
/* calculate the force on the coeff's (up or down, both) */
/* given an appropriate v_ks (up or down, both).         */
/*==========================================================================*/

void coef_force_calc_full_g(CPEWALD *cpewald,int nstate,int ncoef,
                             int ncoef_max,
                             double *creal,double *cimag, 
                             double *fcreal,double  *fcimag,
                             double *cre_scr,double *cim_scr,
                             double *cp_hess_re,double *cp_hess_im,
                             double *zfft,double *zfft_tmp,
                             double *v_ks,double *v_ks_tau,double *ak2_sm,
                             double *eke_ret,double *pvten_cp,
                             int cp_ptens_calc,double *hmati,
                             COMMUNICATE *communicate,
                             int icoef_form,int icoef_orth,int ifcoef_form,
                             int cp_tau_functional,int cp_min_on,
                             PARA_FFT_PKG3D *cp_para_fft_pkg3d_sm)

/*==========================================================================*/
/*         Begin Routine                                                    */
{/*Begin Routine*/
/*=======================================================================*/
/*            Local variable declarations                               */
#include "../typ_defs/typ_mask.h"

   int is,i,iupper;
   int joff;
   double tpi;
   double aka,akb,akc,xk,yk,zk,cfact;
   double eke;
   int ioff,ncoef_use,ioff2;
   int iii,iis;

   int nfft       = cp_para_fft_pkg3d_sm->nfft_proc;
   int icoef_off  = cp_para_fft_pkg3d_sm->icoef_off;
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

      double sum_check,sum_check_tmp;
      MPI_Comm comm_states = communicate->comm_states;

/* ================================================================= */
/*0) Check the form of the coefficients                              */

  if(icoef_orth!=1){
    printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
    printf("The coefficients must be in orthogonal form    \n");
    printf("on state processor %d in coef_force_calc_full_g   \n",myid_state);
    printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
    fflush(stdout);exit(1);
  }/*endif*/

  if(np_states>1) 
   if( ifcoef_form==0 || icoef_form == 0){
    printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
    printf("The coefs/coef forces must be in transposed (not normal) \n");
    printf("form on state processor %d in coef_force_calc_full_g\n",
           myid_state);
    printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
    fflush(stdout);exit(1);
  }/*endif*/

/*=================================================================*/
/*  Find the upper state limit                                     */

  iupper = nstate;
  if(nstate % 2 == 1){
     iupper = nstate - 1;
  }

/*=================================================================*/
/*  get the forces on the coefs of each state                      */

  for(is=1 ; is<= iupper; is+=2 ){

    ioff = (is-1)*ncoef_max;
    ioff2 = (is)*ncoef_max;

/*==========================================================================*/
/* 1) get the wave functions in real space two at a time                    */
/*   I) double pack the complex zfft array with two real wavefunctions      */

    dble_pack_coef(&creal[ioff],&cimag[ioff],&creal[ioff2],&cimag[ioff2],
                   zfft,cp_para_fft_pkg3d_sm);

/*--------------------------------------------------------------------------*/
/* II) fourier transform the wavefunctions to real space                    */
/*      convention exp(-igr)                                                */

    para_fft_gen3d_fwd_to_r(zfft,zfft_tmp,cp_para_fft_pkg3d_sm);

/*==========================================================================*/
/* 2) get v|psi> in g space and store it in zfft                            */
/*   I) get  v|psi> in real space                                           */

    cp_vpsi(zfft,v_ks,nfft);

/*--------------------------------------------------------------------------*/
/*  II) fourier transform  to g-space                                       */
/*     convention exp(igr)                                                  */

    para_fft_gen3d_bck_to_g(zfft,zfft_tmp,cp_para_fft_pkg3d_sm);

/*==========================================================================*/
/* 3) get forces on coefficients by double unpacking the array zfft         */


    dble_upack_coef_sum(&fcreal[ioff],&fcimag[ioff],
                        &fcreal[ioff2],&fcimag[ioff2],
                        zfft,cp_para_fft_pkg3d_sm);

  }/*endfor is */

/*==========================================================================*/
/*==========================================================================*/
/* 4) if there is an odd number of states, go through                       */
/*      the same procedure using sng_packs                                  */


  if(nstate % 2 != 0){
     is = nstate;
     ioff = (is -1)*ncoef_max;

/*--------------------------------------------------------------------------*/
/*   I) sngl pack                                                           */

     sngl_pack_coef(&creal[ioff],&cimag[ioff],zfft,cp_para_fft_pkg3d_sm);

/*--------------------------------------------------------------------------*/
/* II) fourier transform the wavefunctions to real space                    */
/*      convention exp(-igr)                                                */
 
      para_fft_gen3d_fwd_to_r(zfft,zfft_tmp,cp_para_fft_pkg3d_sm);

/*==========================================================================*/
/* 5) get v|psi> in g space and store it in zfft                            */
/*   I) get  v|psi> in real space                                           */

      cp_vpsi(zfft,v_ks,nfft);

/*--------------------------------------------------------------------------*/
/*   II) fourier transform the result back to g-space */
/*     convention exp(igr)  */

      para_fft_gen3d_bck_to_g(zfft,zfft_tmp,cp_para_fft_pkg3d_sm);

/*==========================================================================*/
/* 6) get forces on coefficients by double unpacking the array zfft         */


   sngl_upack_coef_sum(&fcreal[ioff],&fcimag[ioff],zfft,cp_para_fft_pkg3d_sm);

  }/* endif: odd number of states*/

/*==========================================================================*/
/* 7) If there is an electron KE density dependent functional, calculate    */
/*    this contribution to the force                                        */

  if(cp_tau_functional==1)
    coef_force_tau_fun_full_g(cpewald,nstate,creal,cimag,fcreal,fcimag,
                              cre_scr,cim_scr,zfft,zfft_tmp,v_ks_tau,ak2_sm,
                              pvten_cp,cp_ptens_calc,hmati,communicate,
                              icoef_form,icoef_orth,ifcoef_form,
                              cp_para_fft_pkg3d_sm);

/*==========================================================================*/
/* 8) If doing minimization, Fourier transform the KS potential to g-space  */
/*    and unpack it into the diagonal Hessian                               */

  ncoef_use = (myid_state+1 == np_states ? ncoef-1:ncoef);
  if(cp_min_on > 0 ){

    for(i=1;i<=ncoef_use;i++){ cp_hess_re[i] = cp_hess_im[i] = 0.0;}

    sngl_pack_rho(zfft,v_ks,cp_para_fft_pkg3d_sm);

    para_fft_gen3d_bck_to_g(zfft,zfft_tmp,cp_para_fft_pkg3d_sm);

    sngl_upack_coef(cp_hess_re,cp_hess_im,zfft,cp_para_fft_pkg3d_sm);

  }/* endif cp_min_on */

/*==========================================================================*/
/* 9) calculate the kinetic energy term and add its contribution to the force*/

  tpi = 2.0*M_PI;
  eke = 0.0;
  for(is=1 ; is<= nstate ; is++){
    ioff = (is-1)*ncoef_max;
    for(i=1; i<= ncoef_use ; i++){
      joff = i+icoef_off;
      iis = ioff + i;
      fcreal[iis] -= 2.0*ak2_sm[joff]*creal[iis];
      fcimag[iis] -= 2.0*ak2_sm[joff]*cimag[iis];
      eke += (2.0*ak2_sm[joff]
           *(creal[iis]*creal[iis] + cimag[iis]*cimag[iis]));
    }/*endfor i*/
    if(ncoef_use != ncoef){fcimag[(ioff+ncoef)] = 0.0;}
  }/*endfor*/
  eke *= .50;
  *eke_ret = eke;

/*================================================================================*/
/* 10) If doing minimization, calculat kinetic contribution to diagonal Hessian   */

  if(cp_min_on){    
    for(i=1; i<= ncoef_use+1 ; i++){
      cp_hess_re[i] *= -1;
      cp_hess_im[i] *= -1;
    }/* endfor */
    for(i=1; i<= ncoef_use ; i++){
      joff = i+icoef_off;
      cp_hess_re[i] += 2.0*ak2_sm[joff];
      cp_hess_im[i] += 2.0*ak2_sm[joff];
    }/* endfor */
  }/* endif cp_min_on */


/*================================================================================*/
/* 11) calculate kinetic contribution to pressure tensor                     */

  if(cp_ptens_calc == 1) {
   tpi = 2.0*M_PI;
   for(is=1; is<= nstate; is++){
     ioff = (is-1)*ncoef_max;
     for(i=1; i<= ncoef_use; i++){
       joff = icoef_off + i;
       iis = ioff + i;
       aka = (double)kastore_sm[joff];
       akb = (double)kbstore_sm[joff];
       akc = (double)kcstore_sm[joff];

       xk = (aka*hmati[1] +akb*hmati[2] +akc*hmati[3])*tpi;
       yk = (aka*hmati[4] +akb*hmati[5] +akc*hmati[6])*tpi;
       zk = (aka*hmati[7] +akb*hmati[8] +akc*hmati[9])*tpi;

       cfact = 2.0*(creal[iis]*creal[iis] + cimag[iis]*cimag[iis]);

       pvten_cp[1] += xk*xk*cfact;
       pvten_cp[2] += xk*yk*cfact;
       pvten_cp[3] += xk*zk*cfact;

       pvten_cp[4] += xk*yk*cfact;
       pvten_cp[5] += yk*yk*cfact;
       pvten_cp[6] += yk*zk*cfact;

       pvten_cp[7] += xk*zk*cfact;
       pvten_cp[8] += yk*zk*cfact;
       pvten_cp[9] += zk*zk*cfact;

     }/*endfor*/
   }/*endfor*/
  }/*endif */

/*==========================================================================*/
   }/*end routine : fine fertig terminado finito */
/*==========================================================================*/




/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
/* calculate the force on the coeff's (up or down, both) */
/* given an appropriate v_ks (up or down, both).         */
/*==========================================================================*/

void cp_vpsi(double *rifft,double *v_ks,int nfft)

/*==========================================================================*/
/*         Begin Routine                                                    */
{/*Begin Routine*/
/*=======================================================================*/

   int i,m,iii;
/*-----------------------------------------------------------------------*/


   for(i=1,m=1; i<=nfft; i+=2,m++ ){
     rifft[i]   *= v_ks[m];
     rifft[i+1] *= v_ks[m];
    }/*endfor*/

/*==========================================================================*/
 }/*end routine*/
/*==========================================================================*/


/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void cp_pack_vks(double *rifft,double *v_ks,int nfft)

/*==========================================================================*/
/*         Begin Routine                                                    */
{/*Begin Routine*/
/*=======================================================================*/

   int i,m,iii;
/*-----------------------------------------------------------------------*/


   for(i=1,m=1; i<=nfft; i+=2,m++ ){
     rifft[i]   = v_ks[m];
     rifft[i+1] = 0.0;
    }/*endfor*/

/*==========================================================================*/
 }/*end routine*/
/*==========================================================================*/










