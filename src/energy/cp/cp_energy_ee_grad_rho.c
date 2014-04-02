/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
/*                                                                          */
/*                         PI_MD:                                           */
/*             The future of simulation technology                          */
/*             ------------------------------------                         */
/*                   Module: cp_energy_ee_grad_rho.c                        */
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
#include "../proto_defs/proto_friend_lib_entry.h"


/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void control_grad_rho(CPEWALD *cpewald,CPSCR *cpscr,EWALD *ewald,
                      double *rhocr,double *rhoci,
                      double *del_rho_x,double *del_rho_y,double *del_rho_z,
                      double *del2_rho,double *hmati,
                      double vol,int laplacian_on,int cp_dual_grid_opt,
                      PARA_FFT_PKG3D *cp_para_fft_pkg3d)

/*==========================================================================*/
{/*begin routine*/
/*==========================================================================*/
/*assign local pointers */

 double *grhor_x     =    cpscr->cpscr_grho.g_rhor_x;
 double *grhor_y     =    cpscr->cpscr_grho.g_rhor_y;
 double *grhor_z     =    cpscr->cpscr_grho.g_rhor_z;
 double *grhoi_x     =    cpscr->cpscr_grho.g_rhoi_x;
 double *grhoi_y     =    cpscr->cpscr_grho.g_rhoi_y;
 double *grhoi_z     =    cpscr->cpscr_grho.g_rhoi_z;
 double *g2_rhor     =    cpscr->cpscr_grho.g2_rhor;
 double *g2_rhoi     =    cpscr->cpscr_grho.g2_rhoi;
 double *zfft        =    cpscr->cpscr_wave.zfft;
 double *zfft_tmp    =    cpscr->cpscr_wave.zfft_tmp;

 int    *kastore;
 int    *kbstore;
 int    *kcstore;
 int    nktot;   
 int nfft            =    cp_para_fft_pkg3d->nfft_proc;
 int ncoef_l         =    cp_para_fft_pkg3d->ncoef_proc;
 int ncoef_l_use     =    cp_para_fft_pkg3d->ncoef_use;
 int icoef_off       =    cp_para_fft_pkg3d->icoef_off;


/* local variables */
  int i,j,k,jjj;
  int nfft2;
  int iii;

  double rvol = 1.0/vol;

  if(cp_dual_grid_opt >= 1){
   kastore     =    cpewald->kastr_dens_cp_box;
   kbstore     =    cpewald->kbstr_dens_cp_box;
   kcstore     =    cpewald->kcstr_dens_cp_box;
    nktot      =    cpewald->nktot_dens_cp_box;
  }else{
   kastore     =    ewald->kastr;
   kbstore     =    ewald->kbstr;
   kcstore     =    ewald->kcstr;
    nktot      =    ewald->nktot;
  }/*endif cp_dual_grid_opt*/


/*================================================================= */
/*I)  get all three in g space */
/*    Done at g-level and uses parallel packages                    */

 get_igrho(cpscr,kastore,kbstore,kcstore,rhocr,rhoci,
           hmati,ncoef_l_use,ncoef_l,icoef_off);

/*================================================================= */
/*I.V)  Get g^2 *rho if laplacian of density needed */

 if(laplacian_on == 1) {
   get_g2rho(cpscr,kastore,kbstore,kcstore,rhocr,rhoci,
             hmati,ncoef_l_use,ncoef_l,icoef_off);
 }/* endif laplacian */

/*----------------------------------------------------------------- */
/* ii) pack up x and y components of gradient using the double pack */

  dble_pack_coef(grhor_x,grhoi_x,grhor_y,grhoi_y,zfft,cp_para_fft_pkg3d);

/*----------------------------------------------------------------- */
/* iii) transform to real-space convention exp(-igr)  */

  para_fft_gen3d_fwd_to_r(zfft,zfft_tmp,cp_para_fft_pkg3d);

/*----------------------------------------------------------------- */
/*   iv) put x and y components of gradient into del_rho  */

  dble_upack_rho_scal(zfft,del_rho_x,del_rho_y,cp_para_fft_pkg3d,rvol);

/*================================================================= */
/* v) if no laplacian single pack z conponent in real space   */

  if(laplacian_on == 0){

   sngl_pack_coef(grhor_z,grhoi_z,zfft,cp_para_fft_pkg3d);

/*----------------------------------------------------------------- */
/* vi) transform to real-space convention exp(-igr)  */

   para_fft_gen3d_fwd_to_r(zfft,zfft_tmp,cp_para_fft_pkg3d);

/*----------------------------------------------------------------- */
/* vii) single unpack (old get_del_rho)  */

   sngl_upack_rho_scal(zfft,del_rho_z,cp_para_fft_pkg3d,rvol);

/*================================================================= */

  } else {

/*================================================================= */
/* vi) if laplacian double pack z component and laplacian in real space*/

   dble_pack_coef(grhor_z,grhoi_z,g2_rhor,g2_rhoi,zfft,cp_para_fft_pkg3d); 

   para_fft_gen3d_fwd_to_r(zfft,zfft_tmp,cp_para_fft_pkg3d);

   dble_upack_rho_scal(zfft,del_rho_z,del2_rho,cp_para_fft_pkg3d,rvol);

  }/* endif laplacian */



/*================================================================= */
}/*end routine*/
/*================================================================= */




/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void  grad_corr_lda(CPSCR *cpscr,CPEWALD *cpewald,EWALD *ewald,CELL *cell,
                    double *exc,double *muxc,double vol,CPOPTS *cpopts,
                    double gc_cut,double *pvten_cp,int cp_ptens_calc,
                    int cp_tau_functional,int laplacian_on,int cp_dual_grid_opt,
                    PARA_FFT_PKG3D *cp_para_fft_pkg3d_gen,int nstate)

/*==========================================================================*/
/*         Begin Routine                                                    */
   {/*Begin Routine*/
/*=======================================================================*/
/*                 Local variable declarations                           */

/*------------------------------*/
/* density and gradient */
   double rho_p,g_rho2,g_rhoi,g_rho;
   double unit_gx,unit_gy,unit_gz;
   double volscale;
/*------------------------------*/
/*exchange and correlation energy */
   double  ex_gc,ec_gc,muxc_loc;   
/*------------------------------*/
/*  cg functional (f_xc), df_xc/dn, and df_xc/d|grad n|  */
   double fxc_x,dfxc_dn_x,dfxc_dgn_x,dfxc_dtau_x;
   double fxc_c,dfxc_dn_c,dfxc_dgn_c;
   double dfxc_dtau_c;
/*------------------------------*/
/* kinetic energy density varialbes */
   double tau_up;

   double dfxc_dgn_x_c;
   double lap_rho;
   double dfxc_dln_x,dfxc_dln_c;
/*------------------------------*/
/* Switching function */
   double srho,dsrho;        
   double rho_heali,rho_cut,rsw;
/*------------------------------*/
/* loop variables */
   int i,kk,iii,m,jjj,kkk,k,j;
   int iup=1;
/*----------------------------------------------------------*/
/* Becke beta parameter (can vary depending on functional) */
   double beta=0.0042;
/*-----------------------------------------------*/
/* PBE parameters (can vary depending on variety */
   double mu,kappa,beta_pbe,gamma;

/* Assign local pointers                                                     */

    double    *pv_gga;
    double    *hmati_cp        =    cell->hmati_cp;
    double    *rho_up          =    cpscr->cpscr_rho.rho_up;
    double    *del_rho_up_x    =    cpscr->cpscr_grho.d_rhox_up;
    double    *del_rho_up_y    =    cpscr->cpscr_grho.d_rhoy_up;
    double    *del_rho_up_z    =    cpscr->cpscr_grho.d_rhoz_up;
    double    *del2_rho_up     =    cpscr->cpscr_grho.d2_rho_up;
    double    *del2_rho_up_st  =    cpscr->cpscr_grho.d2_rho_up_store;
    double    *rhocr_up        =    cpscr->cpscr_rho.rhocr_up;
    double    *rhoci_up        =    cpscr->cpscr_rho.rhoci_up;
    double    *rhocr_up_st     =    cpscr->cpscr_grho.rhocr_up_st;
    double    *rhoci_up_st     =    cpscr->cpscr_grho.rhoci_up_st;
    double    *elec_ke_dens_up =    cpscr->cpscr_grho.elec_ke_dens_up;
    double    *v_ks_up         =    cpscr->cpscr_rho.v_ks_up;
    double    *v_ks_tau_up     =    cpscr->cpscr_rho.v_ks_tau_up;
    double    *zfft            =    cpscr->cpscr_wave.zfft;
    double    *zfft_tmp        =    cpscr->cpscr_wave.zfft_tmp;
    int cp_becke               = cpopts->cp_becke;
    int cp_fila_1x             = cpopts->cp_fila_1x;
    int cp_fila_2x             = cpopts->cp_fila_2x;
    int cp_pw91x               = cpopts->cp_pw91x;
    int cp_pbe_x               = cpopts->cp_pbe_x;
    int cp_revpbe_x            = cpopts->cp_revpbe_x;
    int cp_rpbe_x              = cpopts->cp_rpbe_x;
    int cp_xpbe_x              = cpopts->cp_rpbe_x;
    int cp_brx89               = cpopts->cp_brx89;
    int cp_brx2k               = cpopts->cp_brx2k;
    int cp_pw91c               = cpopts->cp_pw91c;
    int cp_pbe_c               = cpopts->cp_pbe_c;
    int cp_xpbe_c              = cpopts->cp_xpbe_c;
    int cp_lyp                 = cpopts->cp_lyp;
    int cp_lypm1               = cpopts->cp_lypm1;
    int cp_tau1_c              = cpopts->cp_tau1_c;
    int   ncoef_l_proc         =    cp_para_fft_pkg3d_gen->ncoef_proc;
    int   nfft_proc            =    cp_para_fft_pkg3d_gen->nfft_proc;
    int   nfft                 =    cp_para_fft_pkg3d_gen->nfft;
    int   nfft2                =    nfft/2;
    int   nfft2_proc           =    nfft_proc/2;

    double ec_test_1,ec_test_2;

/*======================================================================*/
/* First set the Becke beta parameter                                   */

    if(cp_tau1_c==1) beta = 0.0045;

/* And the PBE parameters */
    if(cp_pbe_x){
       mu = 0.2195149727645171;
       kappa = 0.804;
    }
    if(cp_revpbe_x){
      mu = 0.2195149727645171;
      kappa = 1.2450;
    }
    if(cp_xpbe_x){
      mu = 0.23214;
      kappa = 0.91954;
    }
    if(cp_pbe_c){
      beta_pbe = 0.06672455060314922;
      gamma = 0.03109069086965489503494086371273;
    }
    if(cp_xpbe_c){
      beta_pbe = 0.089809;
      gamma =  0.020434;
    }

/*======================================================================*/
/* Initialize and malloc variables     rs                               */

   pv_gga = (double *)cmalloc(9*sizeof(double))-1;
   for(kk=1;kk<=9;kk++){pv_gga[kk]=0.0;}

   rho_heali=1.0/(3.9*gc_cut);
   rho_cut = 0.1*gc_cut;

   ex_gc         = 0.0;
   ec_gc         = 0.0;
   muxc_loc      = 0.0;
   fxc_x         = 0.0;
   dfxc_dn_x     = 0.0;
   dfxc_dgn_x    = 0.0;
   dfxc_dln_x    = 0.0;
   dfxc_dtau_x   = 0.0;
   fxc_c         = 0.0;
   dfxc_dn_c     = 0.0;
   dfxc_dgn_c    = 0.0;
   dfxc_dln_c    = 0.0;
   dfxc_dtau_c   = 0.0;

   ec_test_1 = 0.0; ec_test_2 = 0.0;

/*======================================================================*/
/* All important volume scaling factor                                  */

   volscale = vol/((double)(nfft2));

/*==========================================================================*/
/* if ptens and laplacian on store laplacian of density and 
                                                       density coefficients */


   if(laplacian_on == 1 && cp_ptens_calc == 1){
    for(i=1; i<= nfft2_proc; i++){
      del2_rho_up_st[i] = del2_rho_up[i];
    }/* endfor */

    for(i=1;i<=ncoef_l_proc;i++){
     rhocr_up_st[i] = rhocr_up[i];
     rhoci_up_st[i] = rhoci_up[i];
    }
   }/* endfor */

/*=========================================================================*/
/* II) Loop over real space grid                                           */

   
   for(i=1; i<= nfft2_proc; i++){
/*------------------------------------------------------------------------*/
/* i) Strip out the denisty and gradiant and do it up if necessary        */
      rho_p = rho_up[i];
#ifdef VANDERBILT
      if(rho_p < 0.0) rho_p = 0.0;
#endif
      if(rho_p > rho_cut) {
        g_rho2 = del_rho_up_x[i]*del_rho_up_x[i]
               + del_rho_up_y[i]*del_rho_up_y[i]
               + del_rho_up_z[i]*del_rho_up_z[i];

/*-------------------------------------------------------------------------*/
/* ii) laplacian of rho (if laplacian on)                                  */

       if(laplacian_on == 1){lap_rho = del2_rho_up[i];}

/*------------------------------------------------------------------------------*/
/* iii) electron kinetic energy density (if a tau-dependent functional is used) */

       if(cp_tau_functional == 1){tau_up = elec_ke_dens_up[i];}

/*------------------------------------------------------------------------*/
/* iv) Switching functional                                              */

       rsw = (rho_p - rho_cut)*rho_heali;
       rsw = MIN(rsw,1.0);
       rsw = MAX(rsw,0.0);
       srho = rsw*rsw*(3.0-2.0*rsw);
       dsrho = 6.0*rsw*(1.0-rsw)*rho_heali;

/*------------------------------------------------------------------------*/
/* v) Exchange functional                                                */

      if(cp_becke==1)   becke_gcx_lda(rho_p,g_rho2,&fxc_x,&dfxc_dn_x,&dfxc_dgn_x,beta);
      if(cp_fila_1x==1) fila_1x_lda(rho_p,g_rho2,&fxc_x,&dfxc_dn_x,&dfxc_dgn_x);
      if(cp_fila_2x==1) fila_2x_lda(rho_p,g_rho2,&fxc_x,&dfxc_dn_x,&dfxc_dgn_x);
      if(cp_pw91x==1)   pw91_gcx(rho_p,g_rho2,&fxc_x,&dfxc_dn_x,&dfxc_dgn_x);
      if(cp_pbe_x==1)   pbe_gcx(rho_p,g_rho2,&fxc_x,&dfxc_dn_x,&dfxc_dgn_x,
				mu,kappa);
      
      if(cp_revpbe_x==1)pbe_gcx(rho_p,g_rho2,&fxc_x,&dfxc_dn_x,&dfxc_dgn_x,
				mu,kappa);
      if(cp_rpbe_x==1)  rpbe_gcx(rho_p,g_rho2,&fxc_x,&dfxc_dn_x,&dfxc_dgn_x);
      if(cp_xpbe_x==1)  pbe_gcx(rho_p,g_rho2,&fxc_x,&dfxc_dn_x,&dfxc_dgn_x,
				mu,kappa);
      if(cp_brx89==1)   brx89_lda(rho_p,g_rho2,lap_rho,
                                  tau_up,&fxc_x,&dfxc_dn_x,&dfxc_dgn_x,
                                  &dfxc_dln_x,&dfxc_dtau_x);

      if(cp_brx2k==1)   brx2K_lda(rho_p,g_rho2,lap_rho,
                                  tau_up,&fxc_x,&dfxc_dn_x,&dfxc_dgn_x,
                                  &dfxc_dln_x,&dfxc_dtau_x);

      dfxc_dn_x    = (dfxc_dn_x*srho + fxc_x*dsrho);
      dfxc_dgn_x  *= srho;    
      if(laplacian_on == 1) dfxc_dln_x *= srho;
      if(cp_tau_functional == 1) dfxc_dtau_x *= srho;

      fxc_x       *= srho;
       
/*-------------------------------------------------------------------------*/
/* vi) Correlation functional                                               */


     if(cp_pw91c==1)  pw91_gcc(rho_p,g_rho2,&fxc_c,&dfxc_dn_c,&dfxc_dgn_c);
     if(cp_pbe_c==1)  pbe_gcc(rho_p,g_rho2,&fxc_c,&dfxc_dn_c,&dfxc_dgn_c,
			      beta_pbe,gamma);
     if(cp_xpbe_c==1) pbe_gcc(rho_p,g_rho2,&fxc_c,&dfxc_dn_c,&dfxc_dgn_c,
			      beta_pbe,gamma);
     if(cp_lyp==1)    lyp_gcc(rho_p,g_rho2,&fxc_c,&dfxc_dn_c,&dfxc_dgn_c);
     if(cp_lypm1==1)  lypm1(rho_p,g_rho2,lap_rho,&fxc_c,&dfxc_dn_c,&dfxc_dgn_c,
                                   &dfxc_dln_c);
     if(cp_tau1_c==1) tau1_lda(rho_p, g_rho2, lap_rho, tau_up,&fxc_c,
                               &dfxc_dn_c,&dfxc_dgn_c,&dfxc_dln_c,&dfxc_dtau_c,
                               nstate);


      dfxc_dn_c   = (dfxc_dn_c*srho + fxc_c*dsrho);
      dfxc_dgn_c *= srho;    
      if(laplacian_on == 1) dfxc_dln_c *= srho;
      if(cp_tau_functional == 1) dfxc_dtau_c *= srho;

      fxc_c      *= srho;
      
/*-------------------------------------------------------------------------*/
/* vii) Add contributions to exchange and correlation energy                */

      ex_gc +=  fxc_x;
      ec_gc +=  fxc_c;

/*-------------------------------------------------------------------------*/
/* viii) Get non-gradiant contribution to Kohn-Sham potential                 */

      v_ks_up[i] += (dfxc_dn_x + dfxc_dn_c);
      muxc_loc   += (dfxc_dn_x + dfxc_dn_c)*rho_p;
      
/*-------------------------------------------------------------------------*/
/* ix) Store the tau derivative of functional in the tau KS potential    */

      if(cp_tau_functional == 1) v_ks_tau_up[i] = dfxc_dtau_x + dfxc_dtau_c;
      
/*-------------------------------------------------------------------------*/
/* x) Get unit vector del rho/|del rho| and multiply it                   */
/*                                                  times fxn derivative   */

      g_rho = sqrt(g_rho2);
      g_rhoi = 1.0/g_rho;
      unit_gx = del_rho_up_x[i]*g_rhoi;
      unit_gy = del_rho_up_y[i]*g_rhoi;
      unit_gz = del_rho_up_z[i]*g_rhoi;

      dfxc_dgn_x_c = (dfxc_dgn_x + dfxc_dgn_c);
      del_rho_up_x[i] = unit_gx*dfxc_dgn_x_c;
      del_rho_up_y[i] = unit_gy*dfxc_dgn_x_c;
      del_rho_up_z[i] = unit_gz*dfxc_dgn_x_c;

      if(laplacian_on == 1){del2_rho_up[i] = dfxc_dln_x + dfxc_dln_c;}

/*-------------------------------------------------------------------------*/
/* xi) Get contributions to pressure tensor  */

    if(cp_ptens_calc==1){
      ptens_gga(pv_gga,unit_gx,unit_gy,unit_gz,g_rho,fxc_x,
                dfxc_dn_x,dfxc_dgn_x,rho_p);
      ptens_gga(pv_gga,unit_gx,unit_gy,unit_gz,g_rho,fxc_c,
                dfxc_dn_c,dfxc_dgn_c,rho_p);
    }/*endif*/
   } else {
      del_rho_up_x[i] = 0.0;
      del_rho_up_y[i] = 0.0;
      del_rho_up_z[i] = 0.0;
      if(laplacian_on == 1){del2_rho_up[i] = 0.0;}
  }/* endif rho_p>gc_cut*/
 }/*endfor*/


/*===========================================================================*/
/* III) the following routine performs r-->g space fourier                   */
/* transformation with a dot product with ig while                           */
/* in g-space and then a back transformation g-->r                           */
/* to complete the gradient correction term in the                           */
/* kohn-sham potential  and add it in                                       */

   yon_hither(cpscr,cpewald,ewald,del_rho_up_x,del_rho_up_y,del_rho_up_z,
              del2_rho_up,hmati_cp,laplacian_on,cp_dual_grid_opt,
              cp_para_fft_pkg3d_gen); 

   for(i=1,m=1 ; i<= nfft_proc; i+=2,m++){
     v_ks_up[m] -= zfft[i];
     muxc_loc -= zfft[i]*rho_up[m];
   }/* endfor */

/*===========================================================================*/
/*IV) add gradient corrections to energy and kohn-sham potential             */


    ex_gc *= volscale;
    ec_gc *= volscale;
    muxc_loc *= volscale;
    *exc +=  (ex_gc+ec_gc);
    *muxc += muxc_loc;

/*===========================================================================*/
/* V) Pressure tensor contributions including HEINOUS laplacian contribution */

 if(laplacian_on == 1 && cp_ptens_calc == 1){
   ptens_gga_lap(pv_gga,cpscr,cpewald,ewald,hmati_cp,vol,
                 iup,cp_para_fft_pkg3d_gen);
 }

/* everybody must get stoned */
    if(cp_ptens_calc==1){
        pvten_cp[1] +=  pv_gga[1]*volscale;
        pvten_cp[2] +=  pv_gga[2]*volscale;
        pvten_cp[3] +=  pv_gga[3]*volscale;
        pvten_cp[4] +=  pv_gga[4]*volscale;
        pvten_cp[5] +=  pv_gga[5]*volscale;
        pvten_cp[6] +=  pv_gga[6]*volscale;
        pvten_cp[7] +=  pv_gga[7]*volscale;
        pvten_cp[8] +=  pv_gga[8]*volscale;
        pvten_cp[9] +=  pv_gga[9]*volscale;
        pvten_cp[1] -= (ex_gc + ec_gc);
        pvten_cp[5] -= (ex_gc + ec_gc);
        pvten_cp[9] -= (ex_gc + ec_gc);
     }/*endif*/
     cfree(&pv_gga[1]);

/*==========================================================================*/
 }/* end routine */
/*==========================================================================*/




/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/


void  grad_corr_lsda(CPSCR *cpscr,CPEWALD *cpewald,EWALD *ewald,CELL *cell,
                     double *exc,double *muxc,double vol,CPOPTS *cpopts,
                     double gc_cut,
                     double *pvten_cp,int cp_ptens_calc,
                     int cp_tau_functional,int nstate_up,int nstate_dn, 
                     int laplacian_on,
                     int cp_dual_grid_opt,
                     PARA_FFT_PKG3D *cp_para_fft_pkg3d_gen)

/*==========================================================================*/
/*         Begin Routine                                                    */
   {/*Begin Routine*/
/*=======================================================================*/
/*                 Local variable declarations                           */


/*----------------------------------------------------------*/
/* density and gradient */
   double rho_p_up,g_rho2_up,g_rhoi_up,g_rho_up;
   double unit_gx_up,unit_gy_up,unit_gz_up;
   double rho_p_dn,g_rho2_dn,g_rhoi_dn,g_rho_dn;
   double unit_gx_dn,unit_gy_dn,unit_gz_dn;
   double drhox_up,drhoy_up,drhoz_up;
   double drhox_dn,drhoy_dn,drhoz_dn;
   double volscale,rho_p,rho_p_bar;
/*----------------------------------------------------------*/
/*exchange and correlation energy */
   double  ex_gc,ec_gc,muxc_loc;   
/*----------------------------------------------------------*/
/*  cg functional (f_xc), df_xc/dn, and df_xc/d|grad n|  */
   double fxc_x;
   double fxc_x_up,dfxc_dn_x_up,dfxc_dgn_x_up;
   double fxc_x_dn,dfxc_dn_x_dn,dfxc_dgn_x_dn;
   double fxc_c;
   double dfxc_dn_c_up;
   double dfxc_dgxn_c_up,dfxc_dgyn_c_up,dfxc_dgzn_c_up;
   double dfxc_dn_c_dn;
   double dfxc_dgxn_c_dn,dfxc_dgyn_c_dn,dfxc_dgzn_c_dn;
   double dfxc_dgn_x_c;
   double dfxc_dgn_c_up,dfxc_dgn_c_dn;
   double dfxc_dtau_x_up,dfxc_dtau_x_dn;
   double dfxc_dtau_c_up,dfxc_dtau_c_dn;

/*----------------------------------------------------------*/
/* Switching function */
   double srho_up,dsrho_up;        
   double srho_dn,dsrho_dn;        
   double srho_bar,dsrho_bar;        
   double rho_heali,rho_cut,rsw;
/*----------------------------------------------------------*/
/* loop variables */
   int i,kk,iii,m;
   int iup=1,idn=0;
/*----------------------------------------------------------*/
/* Becke beta parameter (can vary depending on functional) */
   double beta=0.0042;
/*-----------------------------------------------*/
/* PBE parameters (can vary depending on variety */
   double mu,kappa,beta_pbe,gamma;

/* Assign local pointers                                                     */

    double    *pv_gga;
    double    *hmati          =    cell->hmati_cp;

    double    *rho_up         =    cpscr->cpscr_rho.rho_up;
    double    *del_rho_up_x   =    cpscr->cpscr_grho.d_rhox_up;
    double    *del_rho_up_y   =    cpscr->cpscr_grho.d_rhoy_up;
    double    *del_rho_up_z   =    cpscr->cpscr_grho.d_rhoz_up;

    double    *rho_dn         =    cpscr->cpscr_rho.rho_dn;
    double    *del_rho_dn_x   =    cpscr->cpscr_grho.d_rhox_dn;
    double    *del_rho_dn_y   =    cpscr->cpscr_grho.d_rhoy_dn;
    double    *del_rho_dn_z   =    cpscr->cpscr_grho.d_rhoz_dn;

    double *del2_rho_up       =    cpscr->cpscr_grho.d2_rho_up;
    double *del2_rho_dn       =    cpscr->cpscr_grho.d2_rho_dn;
    double *del2_rho_up_st    =    cpscr->cpscr_grho.d2_rho_up_store;
    double *del2_rho_dn_st    =    cpscr->cpscr_grho.d2_rho_dn_store;
    double *elec_ke_dens_up   =    cpscr->cpscr_grho.elec_ke_dens_up;
    double *elec_ke_dens_dn   =    cpscr->cpscr_grho.elec_ke_dens_dn;
    double *rhocr_up       =    cpscr->cpscr_rho.rhocr_up;
    double *rhoci_up       =    cpscr->cpscr_rho.rhoci_up;
    double *rhocr_up_st    =    cpscr->cpscr_grho.rhocr_up_st;
    double *rhoci_up_st    =    cpscr->cpscr_grho.rhoci_up_st;
    double *rhocr_dn       =    cpscr->cpscr_rho.rhocr_dn;
    double *rhoci_dn       =    cpscr->cpscr_rho.rhoci_dn;
    double *rhocr_dn_st    =    cpscr->cpscr_grho.rhocr_dn_st;
    double *rhoci_dn_st    =    cpscr->cpscr_grho.rhoci_dn_st;
    double lap_rho_up;                             /* Up laplacian rho*/
    double lap_rho_dn;                             /* Dn laplacian rho */
    double tau_up;                                 /* Up elec. Ke dens */
    double tau_dn;                                 /* Dn elec. Ke dens */
    double dfxc_dln_c_up;
    double dfxc_dln_c_dn;
    double dfxc_dln_x_up;
    double dfxc_dln_x_dn;

    double    *v_ks_up        =    cpscr->cpscr_rho.v_ks_up;
    double    *v_ks_dn        =    cpscr->cpscr_rho.v_ks_dn;
    double    *v_ks_tau_up    =    cpscr->cpscr_rho.v_ks_tau_up;
    double    *v_ks_tau_dn    =    cpscr->cpscr_rho.v_ks_tau_dn;
    double    *zfft           =    cpscr->cpscr_wave.zfft;
    double    *zfft_tmp       =    cpscr->cpscr_wave.zfft_tmp;

    int   ncoef_l             =    cp_para_fft_pkg3d_gen->ncoef_proc;
    int   nfft_proc           =    cp_para_fft_pkg3d_gen->nfft_proc;
    int   nfft                =    cp_para_fft_pkg3d_gen->nfft;
    int   nfft2               =    nfft/2;
    int   nfft2_proc          =    nfft_proc/2;


/* Functional options */

    int cp_becke   = cpopts->cp_becke;
    int cp_brx89   = cpopts->cp_brx89;
    int cp_brx2k   = cpopts->cp_brx2k;
    int cp_fila_1x = cpopts->cp_fila_1x;
    int cp_fila_2x = cpopts->cp_fila_2x;
    int cp_pbe_x   = cpopts->cp_pbe_x;
    int cp_revpbe_x= cpopts->cp_revpbe_x;
    int cp_rpbe_x  = cpopts->cp_rpbe_x;
    int cp_xpbe_x  = cpopts->cp_xpbe_x;
    int cp_lyp     = cpopts->cp_lyp;
    int cp_lypm1   = cpopts->cp_lypm1;
    int cp_tau1_c  = cpopts->cp_tau1_c;
    int cp_pbe_c   = cpopts->cp_pbe_c;
    int cp_xpbe_c  = cpopts->cp_xpbe_c;

/*======================================================================*/
/* First set the Becke beta parameter                                   */

    if(cp_tau1_c==1) beta = 0.0045;

/* And the PBE parameters */
    if(cp_pbe_x){
       mu = 0.2195149727645171;
       kappa = 0.804;
    }
    if(cp_revpbe_x){
      mu = 0.2195149727645171;
      kappa = 1.2450;
    }
    if(cp_xpbe_x){
      mu = 0.23214;
      kappa = 0.91954;
    }
    if(cp_pbe_c){
      beta_pbe = 0.06672455060314922;
      gamma = 0.03109069086965489503494086371273;
    }
    if(cp_xpbe_c){
      beta_pbe = 0.089809;
      gamma =  0.020434;
    }

/*==============================================================================*/
/* Initialize and malloc variables     rs                                       */

   pv_gga = (double *)cmalloc(9*sizeof(double))-1;
   for(kk=1;kk<=9;kk++){pv_gga[kk]=0.0;}

   rho_heali=1.0/(3.9*gc_cut);
   rho_cut = 0.1*gc_cut;

   ex_gc          = 0.0;
   ec_gc          = 0.0;
   muxc_loc       = 0.0;

   fxc_x          = 0.0;

   fxc_x_up       = 0.0;
   dfxc_dn_x_up   = 0.0;
   dfxc_dgn_x_up  = 0.0;
   dfxc_dln_x_up  = 0.0;
   dfxc_dtau_x_up = 0.0;
   fxc_x_dn       = 0.0;
   dfxc_dn_x_dn   = 0.0;
   dfxc_dgn_x_dn  = 0.0;
   dfxc_dln_x_dn  = 0.0;
   dfxc_dtau_x_dn = 0.0;

   fxc_c          = 0.0;

   dfxc_dn_c_up   = 0.0;
   dfxc_dgxn_c_up = 0.0;
   dfxc_dgyn_c_up = 0.0;
   dfxc_dgzn_c_up = 0.0;
   dfxc_dn_c_dn   = 0.0;
   dfxc_dgxn_c_dn = 0.0;
   dfxc_dgyn_c_dn = 0.0;
   dfxc_dgzn_c_dn = 0.0;

   dfxc_dln_c_up = 0.0;
   dfxc_dln_c_dn = 0.0;

   dfxc_dtau_c_up = 0.0;
   dfxc_dtau_c_dn = 0.0;

   nfft2 = nfft/2;
/*==============================================================================*/
/* if ptens and laplacian on store laplacian of density and density coefficients */

   if(laplacian_on == 1 && cp_ptens_calc == 1){
    for(i=1; i<= nfft2_proc; i++){
      del2_rho_up_st[i] = del2_rho_up[i];
    }/* endfor */
    for(i=1;i<=ncoef_l;i++){
     rhocr_up_st[i] = rhocr_up[i];
     rhoci_up_st[i] = rhoci_up[i];
    }
    if(nstate_dn != 0){
     for(i=1; i<= nfft2_proc; i++){
       del2_rho_dn_st[i] = del2_rho_dn[i];
     }/* endfor */
     for(i=1;i<=ncoef_l;i++){
      rhocr_dn_st[i] = rhocr_dn[i];
      rhoci_dn_st[i] = rhoci_dn[i];
     }
    }/* endif */
   }/* endfor */

/*==============================================================================*/
/* If Down State = 0                                                            */

 if(nstate_dn == 0){
   for(i=1; i<= nfft2_proc; i++){
     rho_dn[i]       = 0.0;
     del_rho_dn_x[i] = 0.0;
     del_rho_dn_y[i] = 0.0;
     del_rho_dn_z[i] = 0.0;
   }/*endfor*/
 }/*endif*/

/*===========================================================================*/
/* II) Loop over real space grid                                                */


 for(i=1; i<= nfft2_proc; i++){
/*---------------------------------------------------------------------------*/
/* i) Strip out the density and gradiant and do it up if necessary           */
   rho_p_up = rho_up[i];
   rho_p_dn = rho_dn[i];
#ifdef VANDERBILT
   if(rho_p_up < 0.0) rho_p_up = 0.0;
   if(rho_p_dn < 0.0) rho_p_dn = 0.0;
#endif
   rho_p_bar = 0.5*(rho_p_up+rho_p_dn);
   if((rho_p_up > rho_cut)||(rho_p_dn > rho_cut)) {
     g_rho2_up = del_rho_up_x[i]*del_rho_up_x[i]
       + del_rho_up_y[i]*del_rho_up_y[i]
       + del_rho_up_z[i]*del_rho_up_z[i];
     g_rho2_dn = del_rho_dn_x[i]*del_rho_dn_x[i]
       + del_rho_dn_y[i]*del_rho_dn_y[i]
       + del_rho_dn_z[i]*del_rho_dn_z[i];
     
     drhox_up = del_rho_up_x[i];
     drhoy_up = del_rho_up_y[i];
     drhoz_up = del_rho_up_z[i];
     drhox_dn = del_rho_dn_x[i];
     drhoy_dn = del_rho_dn_y[i];
     drhoz_dn = del_rho_dn_z[i];

/*---------------------------------------------------------------------------*/
/* ii) laplacian of rho (no longer needed)                                   */

     if(laplacian_on == 1){
       lap_rho_up = del2_rho_up[i];
       lap_rho_dn = del2_rho_dn[i];
     }/* endif */
     
/*---------------------------------------------------------------------------*/
/* iii) electron kinetic energy density (if a tau-dependent functional is used) */

     if(cp_tau_functional == 1){
       tau_up = elec_ke_dens_up[i];
       tau_dn = elec_ke_dens_dn[i];
     }

/*---------------------------------------------------------------------------*/
/* iv) Switching functionals                                                 */

     rsw = (rho_p_up - rho_cut)*rho_heali;
     rsw = MIN(rsw,1.0);
     rsw = MAX(rsw,0.0);
     srho_up = rsw*rsw*(3.0-2.0*rsw);
     dsrho_up = 6.0*rsw*(1.0-rsw)*rho_heali;
     
     rsw = (rho_p_dn - rho_cut)*rho_heali;
     rsw = MIN(rsw,1.0);
     rsw = MAX(rsw,0.0);
     srho_dn = rsw*rsw*(3.0-2.0*rsw);
     dsrho_dn = 6.0*rsw*(1.0-rsw)*rho_heali;
     
     rsw = (rho_p_bar - rho_cut)*rho_heali;
     rsw = MIN(rsw,1.0);
     rsw = MAX(rsw,0.0);
     srho_bar = rsw*rsw*(3.0-2.0*rsw);
     dsrho_bar = 3.0*rsw*(1.0-rsw)*rho_heali;
     
/*---------------------------------------------------------------------------*/
/* v) Exchange functional                                                    */

 /* Becke */
 /*--------------*/

     if(cp_becke==1){ 
       becke_gcx_lsda(rho_p_up,g_rho2_up,
		      &fxc_x_up,&dfxc_dn_x_up,&dfxc_dgn_x_up,beta);
       
       dfxc_dn_x_up = (dfxc_dn_x_up*srho_up + fxc_x_up*dsrho_up);
       dfxc_dgn_x_up *= srho_up;    
       
       
       if(nstate_dn > 0){
         becke_gcx_lsda(rho_p_dn,g_rho2_dn,
                        &fxc_x_dn,&dfxc_dn_x_dn,&dfxc_dgn_x_dn,beta);
         dfxc_dn_x_dn = (dfxc_dn_x_dn*srho_dn + fxc_x_dn*dsrho_dn);
         dfxc_dgn_x_dn *= srho_dn;    
       }else{ 
         dfxc_dn_x_dn  = 0.0;
         dfxc_dgn_x_dn = 0.0;
         srho_dn  = 0.0;
         dsrho_dn = 0.0;
         fxc_x_dn = 0.0;
       } /* endif nstate_dn */
       
     }/*endif*/
     
 /* Filatov-Thiel A */
 /*--------------*/
     
     if(cp_fila_1x==1){ 
       fila_1x_lsda(rho_p_up,g_rho2_up,
		    &fxc_x_up,&dfxc_dn_x_up,&dfxc_dgn_x_up);
       
       dfxc_dn_x_up = (dfxc_dn_x_up*srho_up + fxc_x_up*dsrho_up);
       dfxc_dgn_x_up *= srho_up;    
       
       
       if(nstate_dn > 0){
         fila_1x_lsda(rho_p_dn,g_rho2_dn,
                      &fxc_x_dn,&dfxc_dn_x_dn,&dfxc_dgn_x_dn);
         dfxc_dn_x_dn = (dfxc_dn_x_dn*srho_dn + fxc_x_dn*dsrho_dn);
         dfxc_dgn_x_dn *= srho_dn;    
       }else{ 
         dfxc_dn_x_dn  = 0.0;
         dfxc_dgn_x_dn = 0.0;
         srho_dn  = 0.0;
         dsrho_dn = 0.0;
         fxc_x_dn = 0.0;
       } /* endif nstate_dn */
       
     }/*endif*/

 /* Filatov-Thiel B */
 /*--------------*/
     
     if(cp_fila_2x==1){ 
       fila_2x_lsda(rho_p_up,g_rho2_up,
		    &fxc_x_up,&dfxc_dn_x_up,&dfxc_dgn_x_up);
       
       dfxc_dn_x_up = (dfxc_dn_x_up*srho_up + fxc_x_up*dsrho_up);
       dfxc_dgn_x_up *= srho_up;    
       
       
       if(nstate_dn > 0){
         fila_2x_lsda(rho_p_dn,g_rho2_dn,
                      &fxc_x_dn,&dfxc_dn_x_dn,&dfxc_dgn_x_dn);
         dfxc_dn_x_dn = (dfxc_dn_x_dn*srho_dn + fxc_x_dn*dsrho_dn);
         dfxc_dgn_x_dn *= srho_dn;    
       }else{ 
         dfxc_dn_x_dn  = 0.0;
         dfxc_dgn_x_dn = 0.0;
         srho_dn  = 0.0;
         dsrho_dn = 0.0;
         fxc_x_dn = 0.0;
       } /* endif nstate_dn */
       
     }/*endif*/

 /* Becke-Roussel 89 */
 /*------------------*/

     if(cp_brx89==1){
       brx89_lsda(rho_p_up,g_rho2_up,lap_rho_up,
		  tau_up,&fxc_x_up,&dfxc_dn_x_up,&dfxc_dgn_x_up,
		  &dfxc_dln_x_up,&dfxc_dtau_x_up);
       
       dfxc_dn_x_up    = (dfxc_dn_x_up*srho_up + fxc_x_up*dsrho_up);
       dfxc_dgn_x_up  *= srho_up;    
       dfxc_dln_x_up  *= srho_up;
       dfxc_dtau_x_up *= srho_up;
       
       if(nstate_dn > 0){
	 
	 brx89_lsda(rho_p_dn,g_rho2_dn,lap_rho_dn,
		    tau_dn,&fxc_x_dn,&dfxc_dn_x_dn,&dfxc_dgn_x_dn,
		    &dfxc_dln_x_dn,&dfxc_dtau_x_dn);
	 
	 dfxc_dn_x_dn    = (dfxc_dn_x_dn*srho_dn + fxc_x_dn*dsrho_dn);
	 dfxc_dgn_x_dn  *= srho_dn;    
	 dfxc_dln_x_dn  *= srho_dn;
	 dfxc_dtau_x_dn *= srho_dn;
	 
       } else {

	 dfxc_dn_x_dn   = 0.0;
	 dfxc_dgn_x_dn  = 0.0;
	 dfxc_dln_x_dn  = 0.0;
	 dfxc_dtau_x_dn = 0.0;
	 
       }/* endif nstate_dn */
       
     }/* endif brx89 */
     
 /* Becke-Roussel 2K */
 /*------------------*/

     if(cp_brx2k==1){
       brx2K_lsda(rho_p_up,g_rho2_up,lap_rho_up,
		  tau_up,&fxc_x_up,&dfxc_dn_x_up,&dfxc_dgn_x_up,
		  &dfxc_dln_x_up,&dfxc_dtau_x_up);
       
       dfxc_dn_x_up    = (dfxc_dn_x_up*srho_up + fxc_x_up*dsrho_up);
       dfxc_dgn_x_up  *= srho_up;    
       dfxc_dln_x_up  *= srho_up;
       dfxc_dtau_x_up *= srho_up;
       
       if(nstate_dn > 0){
	 
	 brx2K_lsda(rho_p_dn,g_rho2_dn,lap_rho_dn,
		    tau_dn,&fxc_x_dn,&dfxc_dn_x_dn,&dfxc_dgn_x_dn,
		    &dfxc_dln_x_dn,&dfxc_dtau_x_dn);
	 
	 dfxc_dn_x_dn    = (dfxc_dn_x_dn*srho_dn + fxc_x_dn*dsrho_dn);
	 dfxc_dgn_x_dn  *= srho_dn;    
	 dfxc_dln_x_dn  *= srho_dn;
	 dfxc_dtau_x_dn *= srho_dn;
	 
       } else {
	 
	 dfxc_dn_x_dn   = 0.0;
	 dfxc_dgn_x_dn  = 0.0;
	 dfxc_dln_x_dn  = 0.0;
	 dfxc_dtau_x_dn = 0.0;
	 
       }/* endif nstate_dn */
       
     }/* endif brx2k */
     

 /* PBE */
 /*-----*/

     if(cp_pbe_x==1){ 
       pbe_gcx(2.0*rho_p_up,4.0*g_rho2_up,
	       &fxc_x_up,&dfxc_dn_x_up,&dfxc_dgn_x_up,mu,kappa);
       fxc_x_up /= 2.0;
       dfxc_dn_x_up = (dfxc_dn_x_up*srho_up + fxc_x_up*dsrho_up);
       dfxc_dgn_x_up *= srho_up;    
       
       if(nstate_dn > 0){
	 pbe_gcx(2.0*rho_p_dn,4.0*g_rho2_dn,
		 &fxc_x_dn,&dfxc_dn_x_dn,&dfxc_dgn_x_dn,mu,kappa);
	 fxc_x_dn /= 2.0;
	 dfxc_dn_x_dn = (dfxc_dn_x_dn*srho_dn + fxc_x_dn*dsrho_dn);
	 dfxc_dgn_x_dn *= srho_dn;    
       }else{ 
	 dfxc_dn_x_dn  = 0.0;
	 dfxc_dgn_x_dn = 0.0;
	 srho_dn  = 0.0;
	 dsrho_dn = 0.0;
	 fxc_x_dn = 0.0;
       } /* endif nstate_dn */
	 
     }/*endif pbe*/

 /* revPBE */
 /*-----*/

     if(cp_revpbe_x==1){ 
       pbe_gcx(2.0*rho_p_up,4.0*g_rho2_up,
	       &fxc_x_up,&dfxc_dn_x_up,&dfxc_dgn_x_up,mu,kappa);
       fxc_x_up /= 2.0;
       dfxc_dn_x_up = (dfxc_dn_x_up*srho_up + fxc_x_up*dsrho_up);
       dfxc_dgn_x_up *= srho_up;    
       
       if(nstate_dn > 0){
	 pbe_gcx(2.0*rho_p_dn,4.0*g_rho2_dn,
		 &fxc_x_dn,&dfxc_dn_x_dn,&dfxc_dgn_x_dn,mu,kappa);
	 fxc_x_dn /= 2.0;
	 dfxc_dn_x_dn = (dfxc_dn_x_dn*srho_dn + fxc_x_dn*dsrho_dn);
	 dfxc_dgn_x_dn *= srho_dn;    
       }else{ 
	 dfxc_dn_x_dn  = 0.0;
	 dfxc_dgn_x_dn = 0.0;
	 srho_dn  = 0.0;
	 dsrho_dn = 0.0;
	 fxc_x_dn = 0.0;
       } /* endif nstate_dn */
	 
     }/*endif revpbe*/

 /* rPBE */
 /*-----*/

     if(cp_rpbe_x==1){ 
       rpbe_gcx(2.0*rho_p_up,4.0*g_rho2_up,
	       &fxc_x_up,&dfxc_dn_x_up,&dfxc_dgn_x_up);
       fxc_x_up /= 2.0;
       dfxc_dn_x_up = (dfxc_dn_x_up*srho_up + fxc_x_up*dsrho_up);
       dfxc_dgn_x_up *= srho_up;    
       
       if(nstate_dn > 0){
	 rpbe_gcx(2.0*rho_p_dn,4.0*g_rho2_dn,
		 &fxc_x_dn,&dfxc_dn_x_dn,&dfxc_dgn_x_dn);
	 fxc_x_dn /= 2.0;
	 dfxc_dn_x_dn = (dfxc_dn_x_dn*srho_dn + fxc_x_dn*dsrho_dn);
	 dfxc_dgn_x_dn *= srho_dn;    
       }else{ 
	 dfxc_dn_x_dn  = 0.0;
	 dfxc_dgn_x_dn = 0.0;
	 srho_dn  = 0.0;
	 dsrho_dn = 0.0;
	 fxc_x_dn = 0.0;
       } /* endif nstate_dn */
	 
     }/*endif rpbe*/

 /* xPBE */
 /*-----*/

     if(cp_xpbe_x==1){ 
       pbe_gcx(2.0*rho_p_up,4.0*g_rho2_up,
	       &fxc_x_up,&dfxc_dn_x_up,&dfxc_dgn_x_up,mu,kappa);
       fxc_x_up /= 2.0;
       dfxc_dn_x_up = (dfxc_dn_x_up*srho_up + fxc_x_up*dsrho_up);
       dfxc_dgn_x_up *= srho_up;    
       
       if(nstate_dn > 0){
	 pbe_gcx(2.0*rho_p_dn,4.0*g_rho2_dn,
		 &fxc_x_dn,&dfxc_dn_x_dn,&dfxc_dgn_x_dn,mu,kappa);
	 fxc_x_dn /= 2.0;
	 dfxc_dn_x_dn = (dfxc_dn_x_dn*srho_dn + fxc_x_dn*dsrho_dn);
	 dfxc_dgn_x_dn *= srho_dn;    
       }else{ 
	 dfxc_dn_x_dn  = 0.0;
	 dfxc_dgn_x_dn = 0.0;
	 srho_dn  = 0.0;
	 dsrho_dn = 0.0;
	 fxc_x_dn = 0.0;
       } /* endif nstate_dn */
	 
     }/*endif rpbe*/

     fxc_x        = fxc_x_up*srho_up + fxc_x_dn*srho_dn;
     
/*---------------------------------------------------------------------------*/
/* vi) Correlation functional                                                */
     
     if(cp_lyp==1){
       lyp_lsda(rho_p_up,rho_p_dn,drhox_up,drhoy_up,drhoz_up,
		drhox_dn,drhoy_dn,drhoz_dn,
		&fxc_c,&dfxc_dn_c_up,&dfxc_dn_c_dn,
		&dfxc_dgxn_c_up,&dfxc_dgyn_c_up,&dfxc_dgzn_c_up,
		&dfxc_dgxn_c_dn,&dfxc_dgyn_c_dn,&dfxc_dgzn_c_dn);
       
       
       dfxc_dn_c_up   = (dfxc_dn_c_up*srho_bar + fxc_c*dsrho_bar);
       dfxc_dgxn_c_up *= srho_bar;    
       dfxc_dgyn_c_up *= srho_bar;    
       dfxc_dgzn_c_up *= srho_bar;    
       
       if(nstate_dn > 0){
	 dfxc_dn_c_dn   = (dfxc_dn_c_dn*srho_bar + fxc_c*dsrho_bar);
	 dfxc_dgxn_c_dn *= srho_bar;    
	 dfxc_dgyn_c_dn *= srho_bar;    
	 dfxc_dgzn_c_dn *= srho_bar;
       }else{
	 dfxc_dn_c_dn   = 0.0;
	 dfxc_dgxn_c_dn = 0.0;
	 dfxc_dgyn_c_dn = 0.0;
	 dfxc_dgzn_c_dn = 0.0;
       }    
       
       fxc_c  *= srho_bar;
     }/*endif*/
     
     if(cp_tau1_c==1) {
       tau1_lsda(rho_p_up, rho_p_dn,g_rho2_up,g_rho2_dn,
		 lap_rho_up,lap_rho_dn,tau_up,tau_dn,
		 &fxc_c,&dfxc_dn_c_up,&dfxc_dn_c_dn,
		 &dfxc_dgn_c_up,&dfxc_dgn_c_dn,
		 &dfxc_dln_c_up,&dfxc_dln_c_dn,
		 &dfxc_dtau_c_up,&dfxc_dtau_c_dn,nstate_up,nstate_dn);
       
       dfxc_dn_c_up    = (dfxc_dn_c_up*srho_bar + fxc_c*dsrho_bar);
       dfxc_dgn_c_up  *= srho_bar;
       dfxc_dln_c_up  *= srho_bar;
       dfxc_dtau_c_up *= srho_bar;
       if(nstate_dn > 0){
	 dfxc_dn_c_dn    = (dfxc_dn_c_dn*srho_bar + fxc_c*dsrho_bar);
	 dfxc_dgn_c_dn  *= srho_bar;
	 dfxc_dln_c_dn  *= srho_bar;
	 dfxc_dtau_c_dn *= srho_bar;
       } else {
	 dfxc_dn_c_dn   = 0.0;
	 dfxc_dgn_c_dn  = 0.0;
	 dfxc_dln_c_dn  = 0.0;
	 dfxc_dtau_c_dn = 0.0;
       }
       fxc_c  *= srho_bar;
     }/* endif */

/* PBE correlation */
/*-----------------*/
     if(cp_pbe_c==1){
       pbe_gcc_lsda(rho_p_up,rho_p_dn,g_rho2_up,g_rho2_dn,
		    drhox_up,drhoy_up,drhoz_up,
		    drhox_dn,drhoy_dn,drhoz_dn,
		    &fxc_c,
		    &dfxc_dn_c_up,&dfxc_dgn_c_up,
		    &dfxc_dn_c_dn,&dfxc_dgn_c_dn,beta_pbe,gamma);
       
       dfxc_dn_c_up   = (dfxc_dn_c_up*srho_bar + fxc_c*dsrho_bar);
       dfxc_dn_c_dn   = (dfxc_dn_c_dn*srho_bar + fxc_c*dsrho_bar);
       
       dfxc_dgn_c_up *= srho_bar;    
       dfxc_dgn_c_dn *= srho_bar;    
       
       if(nstate_dn <= 0){
	 dfxc_dn_c_dn  = 0.0;
	 dfxc_dgn_c_dn = 0.0;
       }
       
       fxc_c  *= srho_bar;
       
     }/*endif*/
     
/* xPBE correlation */
/*-----------------*/
     if(cp_xpbe_c==1){
       pbe_gcc_lsda(rho_p_up,rho_p_dn,g_rho2_up,g_rho2_dn,
		    drhox_up,drhoy_up,drhoz_up,
		    drhox_dn,drhoy_dn,drhoz_dn,
		    &fxc_c,
		    &dfxc_dn_c_up,&dfxc_dgn_c_up,
		    &dfxc_dn_c_dn,&dfxc_dgn_c_dn,beta_pbe,gamma);
       
       dfxc_dn_c_up   = (dfxc_dn_c_up*srho_bar + fxc_c*dsrho_bar);
       dfxc_dn_c_dn   = (dfxc_dn_c_dn*srho_bar + fxc_c*dsrho_bar);
       
       dfxc_dgn_c_up *= srho_bar;    
       dfxc_dgn_c_dn *= srho_bar;    
       
       if(nstate_dn <= 0){
	 dfxc_dn_c_dn  = 0.0;
	 dfxc_dgn_c_dn = 0.0;
       }
       
       fxc_c  *= srho_bar;
       
     }/*endif*/

/*---------------------------------------------------------------------------*/
/* vii) Add contributions to exchange and correlation energy                 */

     ex_gc +=  fxc_x;
     ec_gc +=  fxc_c;
     muxc_loc += (dfxc_dn_x_up + dfxc_dn_c_up)*rho_p_up +
                 (dfxc_dn_x_dn + dfxc_dn_c_dn)*rho_p_dn;
     
/*-------------------------------------------------------------------------*/
/* viii) Store the tau derivative of functional in the tau KS potential    */

     if(cp_tau_functional == 1) {
       v_ks_tau_up[i] = dfxc_dtau_x_up + dfxc_dtau_c_up;
       v_ks_tau_dn[i] = dfxc_dtau_x_dn + dfxc_dtau_c_dn;
     }
     
/*---------------------------------------------------------------------------*/
/* ix) Get non-gradiant contribution to Kohn-Sham potential                  */

     v_ks_up[i] += (dfxc_dn_x_up + dfxc_dn_c_up);
     v_ks_dn[i] += (dfxc_dn_x_dn + dfxc_dn_c_dn);
     
/*---------------------------------------------------------------------------*/
/* x) Get unit vector del rho/|del rho| and mult by func derivative          */
     
     g_rho_up = sqrt(g_rho2_up);
     g_rhoi_up = 1.0/g_rho_up;
     unit_gx_up = del_rho_up_x[i]*g_rhoi_up;
     unit_gy_up = del_rho_up_y[i]*g_rhoi_up;
     unit_gz_up = del_rho_up_z[i]*g_rhoi_up;
     
     del_rho_up_x[i] = unit_gx_up*dfxc_dgn_x_up;
     del_rho_up_y[i] = unit_gy_up*dfxc_dgn_x_up;
     del_rho_up_z[i] = unit_gz_up*dfxc_dgn_x_up;
     
     if(cp_lyp==1){
       del_rho_up_x[i] += dfxc_dgxn_c_up;
       del_rho_up_y[i] += dfxc_dgyn_c_up;
       del_rho_up_z[i] += dfxc_dgzn_c_up;
     }/* endif */
     
     if(cp_tau1_c==1){
       dfxc_dgn_x_c = dfxc_dgn_x_up + dfxc_dgn_c_up;
       del_rho_up_x[i] += unit_gx_up*dfxc_dgn_c_up;
       del_rho_up_y[i] += unit_gy_up*dfxc_dgn_c_up;
       del_rho_up_z[i] += unit_gz_up*dfxc_dgn_c_up;
     }/* endif */
     
     if(cp_pbe_c==1){
       del_rho_up_x[i] += unit_gx_up*dfxc_dgn_c_up;
       del_rho_up_y[i] += unit_gy_up*dfxc_dgn_c_up;
       del_rho_up_z[i] += unit_gz_up*dfxc_dgn_c_up;
     }/* endif */
     
     g_rho_dn = sqrt(g_rho2_dn);
     g_rhoi_dn = 1.0/g_rho_dn;
     unit_gx_dn = del_rho_dn_x[i]*g_rhoi_dn;
     unit_gy_dn = del_rho_dn_y[i]*g_rhoi_dn;
     unit_gz_dn = del_rho_dn_z[i]*g_rhoi_dn;
     
     del_rho_dn_x[i] = unit_gx_dn*dfxc_dgn_x_dn;
     del_rho_dn_y[i] = unit_gy_dn*dfxc_dgn_x_dn;
     del_rho_dn_z[i] = unit_gz_dn*dfxc_dgn_x_dn;
     
     if(cp_lyp==1){
       del_rho_dn_x[i] += dfxc_dgxn_c_dn;
       del_rho_dn_y[i] += dfxc_dgyn_c_dn;
       del_rho_dn_z[i] += dfxc_dgzn_c_dn;
     }/* endif */
     
     if(cp_tau1_c==1){
       dfxc_dgn_x_c = dfxc_dgn_x_dn + dfxc_dgn_c_dn;
       del_rho_dn_x[i] += unit_gx_dn*dfxc_dgn_c_dn;
       del_rho_dn_y[i] += unit_gy_dn*dfxc_dgn_c_dn;
       del_rho_dn_z[i] += unit_gz_dn*dfxc_dgn_c_dn;
     }/* endif */
     
     if(cp_pbe_c==1){
       del_rho_dn_x[i] += unit_gx_dn*dfxc_dgn_c_dn;
       del_rho_dn_y[i] += unit_gy_dn*dfxc_dgn_c_dn;
       del_rho_dn_z[i] += unit_gz_dn*dfxc_dgn_c_dn;
     }/* endif */
     
/*---------------------------------------------------------------------------*/
/* xi) laplacian for old lyp functional and tau1 functional  */

     if(laplacian_on == 1){
       del2_rho_up[i]  = dfxc_dln_x_up + dfxc_dln_c_up;
       del2_rho_dn[i]  = dfxc_dln_x_dn + dfxc_dln_c_dn;
     }/* endif */


     if(cp_ptens_calc==1){
       ptens_gga(pv_gga,unit_gx_up,unit_gy_up,unit_gz_up,g_rho_up,fxc_x_up,
		 dfxc_dn_x_up,dfxc_dgn_x_up,rho_p_up);
       ptens_gga(pv_gga,unit_gx_dn,unit_gy_dn,unit_gz_dn,g_rho_dn,fxc_x_dn,
		 dfxc_dn_x_dn,dfxc_dgn_x_dn,rho_p_dn);
       if(cp_lyp==1 || cp_lypm1==1){ 
	 ptens_gga_lyp(pv_gga,dfxc_dgxn_c_up,
		       drhox_up,dfxc_dgyn_c_up,drhoy_up,
		       dfxc_dgzn_c_up,drhoz_up,
		       dfxc_dn_c_up,rho_p_up,g_rho_up);
	 ptens_gga_lyp(pv_gga,dfxc_dgxn_c_dn,
		       drhox_dn,dfxc_dgyn_c_dn,drhoy_dn,
		       dfxc_dgzn_c_dn,drhoz_dn,
		       dfxc_dn_c_dn,rho_p_dn,g_rho_dn);
       }/*endif*/
     }/*endif*/
   } else {
     del_rho_up_x[i] = 0.0;
     del_rho_up_y[i] = 0.0;
     del_rho_up_z[i] = 0.0;
     del_rho_dn_x[i] = 0.0;
     del_rho_dn_y[i] = 0.0;
     del_rho_dn_z[i] = 0.0;
     if(laplacian_on == 1){del2_rho_up[i] = del2_rho_dn[i] = 0.0;}
   }/* endif rho_p>gc_cut*/
 }/*endfor*/


/*===========================================================================*/
/* III) the following routine performs r-->g space fourier                   */
/* transformation with a dot product with ig while                           */
/* in g-space and then a back transformation g-->r                           */
/* to complete the gradient correction term in the                           */
/* kohn-sham potential  and add it in                                       */

   yon_hither(cpscr,cpewald,ewald,del_rho_up_x,del_rho_up_y,del_rho_up_z,
              del2_rho_up,hmati,laplacian_on,cp_dual_grid_opt,
              cp_para_fft_pkg3d_gen); 

   /* MAP ? */
   for(kk=1,m=1; kk<= nfft_proc; kk+=2,m++){
     v_ks_up[m] -= zfft[kk];
     muxc_loc   -= zfft[kk]*rho_up[m];
   }/* endfor */

   yon_hither(cpscr,cpewald,ewald,del_rho_dn_x,del_rho_dn_y,del_rho_dn_z,
              del2_rho_dn,hmati,laplacian_on,cp_dual_grid_opt,
              cp_para_fft_pkg3d_gen); 

   for(kk=1,m=1; kk<= nfft_proc; kk+=2,m++){
     v_ks_dn[m] -= zfft[kk];
     muxc_loc   -= zfft[kk]*rho_dn[m];
   }/* endif */

/*===========================================================================*/
/* IV) add gradient corrections to energy and kohn-sham potential             */

    volscale = vol/((double)(nfft2));
    ex_gc *= volscale;
    ec_gc *= volscale;
    muxc_loc *= volscale;
    *exc +=  (ex_gc+ec_gc);
    *muxc += muxc_loc;

/*===========================================================================*/
/* V) Pressure tensor contributions                                          */

 if(laplacian_on == 1 && cp_ptens_calc == 1){
   ptens_gga_lap(pv_gga,cpscr,cpewald,ewald,hmati,vol,iup,cp_para_fft_pkg3d_gen); 
   ptens_gga_lap(pv_gga,cpscr,cpewald,ewald,hmati,vol,idn,cp_para_fft_pkg3d_gen); 
 }/*endfor*/

/* everybody must get stoned */

    if(cp_ptens_calc==1){
        pvten_cp[1] +=  pv_gga[1]*volscale;
        pvten_cp[2] +=  pv_gga[2]*volscale;
        pvten_cp[3] +=  pv_gga[3]*volscale;
        pvten_cp[4] +=  pv_gga[4]*volscale;
        pvten_cp[5] +=  pv_gga[5]*volscale;
        pvten_cp[6] +=  pv_gga[6]*volscale;
        pvten_cp[7] +=  pv_gga[7]*volscale;
        pvten_cp[8] +=  pv_gga[8]*volscale;
        pvten_cp[9] +=  pv_gga[9]*volscale;
        pvten_cp[1] -= (ex_gc + ec_gc);
        pvten_cp[5] -= (ex_gc + ec_gc);
        pvten_cp[9] -= (ex_gc + ec_gc);
     }/*endif*/
     cfree(&pv_gga[1]);


/*==========================================================================*/
 }/*endroutine*/
/*==========================================================================*/


/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void yon_hither(CPSCR *cpscr,CPEWALD *cpewald,EWALD *ewald,
                double *df_dgn_x,double *df_dgn_y,
                double *df_dgn_z,double *df_dln,
                double *hmati,int laplacian_on,int cp_dual_grid_opt,
                PARA_FFT_PKG3D *cp_para_fft_pkg3d_gen)

/*==========================================================================*/
/*         Begin Routine                                                    */
   {/*Begin Routine*/
/*=======================================================================*/
/* assign local pointers  */

    double    *rhocr;
    double    *rhoci;

    double    *grhor_x        =    cpscr->cpscr_grho.g_rhor_x;
    double    *grhor_y        =    cpscr->cpscr_grho.g_rhor_y;
    double    *grhor_z        =    cpscr->cpscr_grho.g_rhor_z;
    double    *grhoi_x        =    cpscr->cpscr_grho.g_rhoi_x;
    double    *grhoi_y        =    cpscr->cpscr_grho.g_rhoi_y;
    double    *grhoi_z        =    cpscr->cpscr_grho.g_rhoi_z;
    double    *g2_rhor        =    cpscr->cpscr_grho.g2_rhor;
    double    *g2_rhoi        =    cpscr->cpscr_grho.g2_rhoi;
    double    *zfft           =    cpscr->cpscr_wave.zfft;
    double    *zfft_tmp       =    cpscr->cpscr_wave.zfft_tmp;
    int       *kastore;
    int       *kbstore;
    int       *kcstore; 
    int        nktot;
    int   ncoef_l_proc        =    cp_para_fft_pkg3d_gen->ncoef_proc;
    int   ncoef_l_use         =    cp_para_fft_pkg3d_gen->ncoef_use;
    int   icoef_l_off         =    cp_para_fft_pkg3d_gen->icoef_off;
    int   nfft_proc           =    cp_para_fft_pkg3d_gen->nfft_proc;
    int   nfft                =    cp_para_fft_pkg3d_gen->nfft;
    int   nfft2               =    nfft/2;
    int   nfft2_proc          =    nfft_proc/2;

/* local variables            */
      int i,iii,j,k,jjj;


   if(cp_dual_grid_opt >= 1){
      kastore        =    cpewald->kastr_dens_cp_box;
      kbstore        =    cpewald->kbstr_dens_cp_box;
      kcstore        =    cpewald->kcstr_dens_cp_box;
      nktot          =    cpewald->nktot_dens_cp_box;
      rhocr          =    cpscr->cpscr_rho.rhocr_up_dens_cp_box;
      rhoci          =    cpscr->cpscr_rho.rhoci_up_dens_cp_box;
   }else{
      kastore        =    ewald->kastr;
      kbstore        =    ewald->kbstr;
      kcstore        =    ewald->kcstr;
      nktot          =    ewald->nktot;
      rhocr          =    cpscr->cpscr_rho.rhocr_up;
      rhoci          =    cpscr->cpscr_rho.rhoci_up;
   }/*endif cp_dual_grid_opt*/


/*=======================================================================*/
/*I)  fourier transform each functional derivative component to g-space  */
/*   i) x and y-components with double pack                              */
/*  a) pack it up                                                        */

   dble_pack_rho(zfft,df_dgn_x,df_dgn_y,cp_para_fft_pkg3d_gen);

/*-----------------------------------------------------------------------*/
/*  b) back transform to g-space exp(igr)                                */

   para_fft_gen3d_bck_to_g(zfft,zfft_tmp,cp_para_fft_pkg3d_gen);

/*-----------------------------------------------------------------------*/
/* c) double unpack                                                      */

   dble_upack_coef(grhor_x,grhoi_x,grhor_y,grhoi_y,zfft,cp_para_fft_pkg3d_gen);

/*-----------------------------------------------------------------------*/
/* if no laplacian z-component by single pack                            */
/*  a) pack it up                                                        */

 if(laplacian_on == 0){
/*-----------------------------------------------------------------------*/

   sngl_pack_rho(zfft,df_dgn_z,cp_para_fft_pkg3d_gen);

/*-----------------------------------------------------------------------*/
/*  b) back transform to g-space exp(igr)                                */
/*     convention exp(igr)                                               */

   para_fft_gen3d_bck_to_g(zfft,zfft_tmp,cp_para_fft_pkg3d_gen);

/*-----------------------------------------------------------------------*/
/*  c) unpack                                                            */

   sngl_upack_coef(grhor_z,grhoi_z,zfft,cp_para_fft_pkg3d_gen);

/*-----------------------------------------------------------------------*/
 } else {
/*-----------------------------------------------------------------------*/
/*   ii) z-component  and laplacian with double pack                     */
/*  a) pack it up                                                        */

   dble_pack_rho(zfft,df_dgn_z,df_dln,cp_para_fft_pkg3d_gen);

/*-----------------------------------------------------------------------*/
/*  b) back transform to g-space exp(igr)                                */
/*     convention exp(igr)                                               */

   para_fft_gen3d_bck_to_g(zfft,zfft_tmp,cp_para_fft_pkg3d_gen);

/*-----------------------------------------------------------------------*/
/* c) double unpack                                                             */

   dble_upack_coef(grhor_z,grhoi_z,g2_rhor,g2_rhoi,zfft,cp_para_fft_pkg3d_gen);

 }/* endif laplacian */

/*=======================================================================*/
/* III) compute dot product with ig                                      */

  get_igdot(cpscr,ewald,cpewald,hmati,ncoef_l_use,ncoef_l_proc,
            icoef_l_off,cp_dual_grid_opt);

  if(laplacian_on == 1){get_g2times(cpscr,ewald,hmati,kastore,kbstore,kcstore,
                                    nktot,ncoef_l_use,ncoef_l_proc,icoef_l_off,
                                    cp_dual_grid_opt);}

/*=======================================================================*/
/*  IV) transform back to real-space                                     */
/*b) pack up (with single pack)and fft                                   */

  sngl_pack_coef(rhocr,rhoci,zfft,cp_para_fft_pkg3d_gen);

/*-----------------------------------------------------------------------*/
/*c) transform to real-space convention exp(-igr)  */

  para_fft_gen3d_fwd_to_r(zfft,zfft_tmp,cp_para_fft_pkg3d_gen);

/*=======================================================================*/
    }/*end routine*/
/*=======================================================================*/


/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
/* subroutine to calculate gradient correction to pressure tensor           */
/*===========================================================================*/

 void ptens_gga(double *pv_gga,double unit_gx,double unit_gy,double unit_gz,
                double g_rho,double fxc,
                double dfxc_dn, double dfxc_dgn,double rho_p)

/*==========================================================================*/
{/* begin routine */
/*==========================================================================*/
/*             Local variables                                              */

  double g_rho_unit_gx_sq, g_rho_unit_gy_sq, g_rho_unit_gz_sq;
  double g_rho_unit_gx_gy, g_rho_unit_gx_gz, g_rho_unit_gy_gz;
  double fxc_x_c;

/*==========================================================================*/
/* I) Assign some useful variables                                          */

    g_rho_unit_gx_sq = unit_gx*unit_gx*g_rho;
    g_rho_unit_gy_sq = unit_gy*unit_gy*g_rho;
    g_rho_unit_gz_sq = unit_gz*unit_gz*g_rho;

    g_rho_unit_gx_gy = unit_gx*unit_gy*g_rho;
    g_rho_unit_gx_gz = unit_gx*unit_gz*g_rho;
    g_rho_unit_gy_gz = unit_gy*unit_gz*g_rho;

/*==========================================================================*/
/* III) Add exchange and correlation contributions to pressure tensor       */

    pv_gga[1] += (dfxc_dn*rho_p  +
                  dfxc_dgn*g_rho +
                  dfxc_dgn*g_rho_unit_gx_sq);

    pv_gga[2] += (dfxc_dgn*g_rho_unit_gx_gy);
    pv_gga[3] += (dfxc_dgn*g_rho_unit_gx_gz);
    pv_gga[4] += (dfxc_dgn*g_rho_unit_gx_gy);
    pv_gga[5] += (dfxc_dn*rho_p  +
                  dfxc_dgn*g_rho +
                  dfxc_dgn*g_rho_unit_gy_sq);

    pv_gga[6] += (dfxc_dgn*g_rho_unit_gy_gz);
    pv_gga[7] += (dfxc_dgn*g_rho_unit_gx_gz);
    pv_gga[8] += (dfxc_dgn*g_rho_unit_gy_gz);
    pv_gga[9] += (dfxc_dn*rho_p  +
                  dfxc_dgn*g_rho +
                  dfxc_dgn*g_rho_unit_gz_sq);

/*======================================================================*/
}/*end routine*/
/*========================================================================*/


/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void ptens_gga_lyp(double *pv_gga,double dfxc_dgxn_c,
                   double drhox,double dfxc_dgyn_c,double drhoy,
                   double dfxc_dgzn_c,double drhoz,
                   double dfxc_dn,double rho_p,
                   double g_rho)

/*==========================================================================*/
   {/*begin routine */
/*==========================================================================*/

 double pv_lyp_dot,pv_diag;
 double pv_lyp_xx,pv_lyp_xy,pv_lyp_xz,pv_lyp_yx,pv_lyp_yy;
 double pv_lyp_yz,pv_lyp_zx,pv_lyp_zy,pv_lyp_zz;

/*==========================================================================*/

      pv_lyp_dot = dfxc_dgxn_c*drhox + dfxc_dgyn_c*drhoy 
                 + dfxc_dgzn_c*drhoz;
      pv_diag    = pv_lyp_dot + dfxc_dn*rho_p;

      pv_lyp_xx = dfxc_dgxn_c*drhox;
      pv_lyp_xy = dfxc_dgxn_c*drhoy;
      pv_lyp_xz = dfxc_dgxn_c*drhoz;
      pv_lyp_yx = dfxc_dgyn_c*drhox;
      pv_lyp_yy = dfxc_dgyn_c*drhoy;
      pv_lyp_yz = dfxc_dgyn_c*drhoz;
      pv_lyp_zx = dfxc_dgzn_c*drhox;
      pv_lyp_zy = dfxc_dgzn_c*drhoy;
      pv_lyp_zz = dfxc_dgzn_c*drhoz;

      pv_gga[1] += (pv_diag + pv_lyp_xx);
      pv_gga[5] += (pv_diag + pv_lyp_yy);
      pv_gga[9] += (pv_diag + pv_lyp_zz);
      pv_gga[2] +=  pv_lyp_xy;
      pv_gga[3] +=  pv_lyp_xz;
      pv_gga[4] +=  pv_lyp_yx;
      pv_gga[6] +=  pv_lyp_yz;
      pv_gga[7] +=  pv_lyp_zx;
      pv_gga[8] +=  pv_lyp_zy;  

/*==========================================================================*/
   }/*end routine*/
/*==========================================================================*/


/*=========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
/* Calculate gradient correction to pressure tensor from laplacian         */
/*=========================================================================*/

void ptens_gga_lap(double *pv_gga,CPSCR *cpscr,CPEWALD *cpewald,
                   EWALD *ewald,double *hmati,double vol,int iup,
                   PARA_FFT_PKG3D *cp_para_fft_pkg3d_lg)

/*==========================================================================*/
/*         Begin Routine                                                    */
{/*Begin Routine*/
/*=======================================================================*/
/* assign local pointers  */

  int i,iopt;  
  double rvol;
  int    nktot       = ewald->nktot;
  double *zfft       = cpscr->cpscr_wave.zfft;
  double *zfft_tmp   = cpscr->cpscr_wave.zfft_tmp;
 double *grhor_x     =    cpscr->cpscr_grho.g_rhor_x;
 double *grhor_y     =    cpscr->cpscr_grho.g_rhor_y;
 double *grhoi_x     =    cpscr->cpscr_grho.g_rhoi_x;
 double *grhoi_y     =    cpscr->cpscr_grho.g_rhoi_y;
 double *rhocr;
 double *rhoci;
 double *del2_rho;
 double *df_dln;
 double *del_rho_x;
 double *del_rho_y;
 double *del_rho_z;
 int   icoef_off           =    cp_para_fft_pkg3d_lg->icoef_off;
 int   ncoef_l_proc        =    cp_para_fft_pkg3d_lg->ncoef_proc;
 int   ncoef_l             =    cp_para_fft_pkg3d_lg->ncoef;
 int   ncoef_l_use         =    cp_para_fft_pkg3d_lg->ncoef_use;
 int   nfft_proc           =    cp_para_fft_pkg3d_lg->nfft_proc;
 int   nfft                =    cp_para_fft_pkg3d_lg->nfft;
 int   nfft2               =    nfft/2;
 int   nfft2_proc          =    nfft_proc/2;

 if(iup == 1){
  rhocr      = cpscr->cpscr_grho.rhocr_up_st;
  rhoci      = cpscr->cpscr_grho.rhoci_up_st;
  del2_rho   = cpscr->cpscr_grho.d2_rho_up_store;
  df_dln     = cpscr->cpscr_grho.d2_rho_up;
  del_rho_x  = cpscr->cpscr_grho.d_rhox_up;
  del_rho_y  = cpscr->cpscr_grho.d_rhoy_up;
  del_rho_z  = cpscr->cpscr_grho.d_rhoz_up;
 } else {
  rhocr      = cpscr->cpscr_grho.rhocr_dn_st;
  rhoci      = cpscr->cpscr_grho.rhoci_dn_st;
  del2_rho   = cpscr->cpscr_grho.d2_rho_dn_store;
  df_dln     = cpscr->cpscr_grho.d2_rho_dn;
  del_rho_x  = cpscr->cpscr_grho.d_rhox_dn;
  del_rho_y  = cpscr->cpscr_grho.d_rhoy_dn;
  del_rho_z  = cpscr->cpscr_grho.d_rhoz_dn;
 }

/*=======================================================================*/
/* 0) Useful constants                                                   */

  rvol = 1.0/vol;

/*=======================================================================*/
/* I) Compute igig tensor two components at a time                       */

 iopt=1;  /* Gets xx and yy components     */

 get_igigrho(cpscr,ewald,rhocr,rhoci,hmati,ncoef_l_use,ncoef_l,iopt,icoef_off);

/*----------------------------------------------------------------- */
/* ii) pack up x and y components of gradient using the double pack */

  dble_pack_coef(grhor_x,grhoi_x,grhor_y,grhoi_y,zfft,cp_para_fft_pkg3d_lg);

/*----------------------------------------------------------------- */
/* iii) transform to real-space convention exp(-igr)  */

  para_fft_gen3d_fwd_to_r(zfft,zfft_tmp,cp_para_fft_pkg3d_lg);

/*----------------------------------------------------------------- */
/*   iv) put x and y components of gradient into del_rho  */

  dble_upack_rho_scal(zfft,del_rho_x,del_rho_y,cp_para_fft_pkg3d_lg,rvol);


/*----------------------------------------------------------------- */
/*   v) Calculate contribution to pressure tensor                   */

   for(i=1;i <=nfft2_proc; i++){
    pv_gga[1] += df_dln[i]*(del2_rho[i] + 2.0*del_rho_x[i]);
    pv_gga[5] += df_dln[i]*(del2_rho[i] + 2.0*del_rho_y[i]);
   }

/*=======================================================================*/
 iopt=2;  /* Gets zz and xy components     */

 get_igigrho(cpscr,ewald,rhocr,rhoci,hmati,ncoef_l_use,ncoef_l,iopt,icoef_off);

/*----------------------------------------------------------------- */
/* ii) pack up x and y components of gradient using the double pack */

  dble_pack_coef(grhor_x,grhoi_x,grhor_y,grhoi_y,zfft,cp_para_fft_pkg3d_lg);

/*----------------------------------------------------------------- */
/* iii) transform to real-space convention exp(-igr)  */

  para_fft_gen3d_fwd_to_r(zfft,zfft_tmp,cp_para_fft_pkg3d_lg);

/*----------------------------------------------------------------- */
/*   iv) put x and y components of gradient into del_rho  */

  dble_upack_rho_scal(zfft,del_rho_x,del_rho_y,cp_para_fft_pkg3d_lg,rvol);


/*----------------------------------------------------------------- */
/*   v) Calculate contribution to pressure tensor                   */

   for(i=1;i <=nfft2_proc; i++){
    pv_gga[9] += df_dln[i]*(del2_rho[i] + 2.0*del_rho_x[i]);
    pv_gga[2] += df_dln[i]*2.0*del_rho_y[i];
    pv_gga[4] += df_dln[i]*2.0*del_rho_y[i];
   }

/*=======================================================================*/
 iopt=3;  /* Gets xz and yz components     */

 get_igigrho(cpscr,ewald,rhocr,rhoci,hmati,ncoef_l_use,ncoef_l,iopt,icoef_off);

/*----------------------------------------------------------------- */
/* ii) pack up x and y components of gradient using the double pack */

  dble_pack_coef(grhor_x,grhoi_x,grhor_y,grhoi_y,zfft,cp_para_fft_pkg3d_lg);

/*----------------------------------------------------------------- */
/* iii) transform to real-space convention exp(-igr)  */

  para_fft_gen3d_fwd_to_r(zfft,zfft_tmp,cp_para_fft_pkg3d_lg);

/*----------------------------------------------------------------- */
/*   iv) put x and y components of gradient into del_rho  */

   dble_upack_rho_scal(zfft,del_rho_x,del_rho_y,cp_para_fft_pkg3d_lg,rvol);

/*----------------------------------------------------------------- */
/*   v) Calculate contribution to pressure tensor                   */

   for(i=1;i <=nfft2_proc; i++){
    pv_gga[3] += df_dln[i]*2.0*del_rho_x[i];
    pv_gga[7] += df_dln[i]*2.0*del_rho_x[i];
    pv_gga[6] += df_dln[i]*2.0*del_rho_y[i];
    pv_gga[8] += df_dln[i]*2.0*del_rho_y[i];
   }

/*=========================================================================*/
  }/*end routine*/
/*========================================================================*/


/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void  get_igdot(CPSCR *cpscr,EWALD *ewald,CPEWALD *cpewald,
                double *hmati,int ncoef_l_use,int ncoef_l, 
                int icoef_off,int cp_dual_grid_opt)

/*==========================================================================*/
/*         Begin Routine                                                    */
   {/*Begin Routine*/
/*=======================================================================*/
/* assign local pointers  */
    double    *rhocr;
    double    *rhoci;
    double    *grhor_x        =    cpscr->cpscr_grho.g_rhor_x;
    double    *grhor_y        =    cpscr->cpscr_grho.g_rhor_y;
    double    *grhor_z        =    cpscr->cpscr_grho.g_rhor_z;
    double    *grhoi_x        =    cpscr->cpscr_grho.g_rhoi_x;
    double    *grhoi_y        =    cpscr->cpscr_grho.g_rhoi_y;
    double    *grhoi_z        =    cpscr->cpscr_grho.g_rhoi_z;
    double    *g2_rhor        =    cpscr->cpscr_grho.g2_rhor;
    double    *g2_rhoi        =    cpscr->cpscr_grho.g2_rhoi;
    int       *kastore;
    int       *kbstore;
    int       *kcstore;
    int        nktot  ;

/* local variables  */
      double  aka,akb,akc,xk,yk,zk,tpi,g2;
      int icount,ncoef_min;

   if(cp_dual_grid_opt >= 1){
      kastore        =    cpewald->kastr_dens_cp_box;
      kbstore        =    cpewald->kbstr_dens_cp_box;
      kcstore        =    cpewald->kcstr_dens_cp_box;
      nktot          =    cpewald->nktot_dens_cp_box;
      rhocr          =    cpscr->cpscr_rho.rhocr_up_dens_cp_box;
      rhoci          =    cpscr->cpscr_rho.rhoci_up_dens_cp_box;
   }else{
      kastore        =    ewald->kastr;
      kbstore        =    ewald->kbstr;
      kcstore        =    ewald->kcstr;
      nktot          =    ewald->nktot;
      rhocr          =    cpscr->cpscr_rho.rhocr_up;
      rhoci          =    cpscr->cpscr_rho.rhoci_up;
   }/*endif cp_dual_grid_opt*/

/*=========================================================================*/
/* I) Compute ig \cdot del rho/|del rho|                                   */

   tpi = 2.0*M_PI;
   ncoef_min = MIN(ncoef_l_use,ncoef_l);
   for(icount=1 ;icount <= ncoef_min ; icount++){
     aka = (double)kastore[(icount+icoef_off)];
     akb = (double)kbstore[(icount+icoef_off)];
     akc = (double)kcstore[(icount+icoef_off)];

     xk = (aka*hmati[1] +akb*hmati[2] +akc*hmati[3])*tpi;
     yk = (aka*hmati[4] +akb*hmati[5] +akc*hmati[6])*tpi;
     zk = (aka*hmati[7] +akb*hmati[8] +akc*hmati[9])*tpi;

     rhocr[icount] = (xk*grhoi_x[icount] + yk*grhoi_y[icount]
                     + zk*grhoi_z[icount]);
     rhoci[icount] = -(xk*grhor_x[icount] + yk*grhor_y[icount]
                    + zk*grhor_z[icount]);
   }/*endfor*/
   if(ncoef_l_use != ncoef_l){
    rhocr[icount] = 0.0;
    rhoci[icount] = 0.0;
   }/* endif */


/*=========================================================================*/
}/*endroutine*/
/*=========================================================================*/



/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void  get_g2times(CPSCR *cpscr,EWALD *ewald,double *hmati,
                  int *kastore,int *kbstore,int *kcstore,int nktot,int ncoef_l_use,
                  int ncoef_l,int icoef_off,int cp_dual_grid_opt)

/*==========================================================================*/
/*         Begin Routine                                                    */
   {/*Begin Routine*/
/*=======================================================================*/
/* assign local pointers  */
    double    *rhocr;
    double    *rhoci;
    double    *g2_rhor        =    cpscr->cpscr_grho.g2_rhor;
    double    *g2_rhoi        =    cpscr->cpscr_grho.g2_rhoi;


/* local variables  */
      double  aka,akb,akc,xk,yk,zk,tpi,g2;
      int icount,ncoef_min;

   if(cp_dual_grid_opt >= 1){
      rhocr          =    cpscr->cpscr_rho.rhocr_up_dens_cp_box;
      rhoci          =    cpscr->cpscr_rho.rhoci_up_dens_cp_box;
   }else{
      rhocr          =    cpscr->cpscr_rho.rhocr_up;
      rhoci          =    cpscr->cpscr_rho.rhoci_up;
   }/*endif cp_dual_grid_opt*/
/*========================================================================*/
/* I) Compute ig \cdot del rho/|del rho|                                  */

   tpi = 2.0*M_PI;
   ncoef_min = MIN(ncoef_l_use,ncoef_l);
   for(icount=1 ;icount <= ncoef_min ; icount++){
     aka = (double)kastore[(icount+icoef_off)];
     akb = (double)kbstore[(icount+icoef_off)];
     akc = (double)kcstore[(icount+icoef_off)];

     xk = (aka*hmati[1] +akb*hmati[2] +akc*hmati[3])*tpi;
     yk = (aka*hmati[4] +akb*hmati[5] +akc*hmati[6])*tpi;
     zk = (aka*hmati[7] +akb*hmati[8] +akc*hmati[9])*tpi;
     g2 = xk*xk + yk*yk + zk*zk;                             
     rhocr[icount] += g2*g2_rhor[icount];     
     rhoci[icount] += g2*g2_rhoi[icount];     

   }/*endfor*/
   if(ncoef_l_use != ncoef_l){
    rhocr[icount] = 0.0;
    rhoci[icount] = 0.0;
   }/* endif */

/*======================================================================*/
}/*endroutine*/
/*======================================================================*/


/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void get_igrho(CPSCR *cpscr,int *kastore,int *kbstore,int *kcstore,
               double *rhocr,double *rhoci,
               double *hmati,int ncoef_l_use,int ncoef_l,int icoef_off)
      
/*=========================================================================*/
{/*begin routine*/
/*=========================================================================*/

/*assign local pointers */
     double *grhor_x        =    cpscr->cpscr_grho.g_rhor_x;
     double *grhor_y        =    cpscr->cpscr_grho.g_rhor_y;
     double *grhor_z        =    cpscr->cpscr_grho.g_rhor_z;
     double *grhoi_x        =    cpscr->cpscr_grho.g_rhoi_x;
     double *grhoi_y        =    cpscr->cpscr_grho.g_rhoi_y;
     double *grhoi_z        =    cpscr->cpscr_grho.g_rhoi_z;

/* local variables */
      
   double  aka,akb,akc,xk,yk,zk,tpi,g2;
   int icount,ncoef_min;

/*===================================================================*/
/* I) Get ig times rho (gradient in g-space)                         */

    tpi = 2.0*M_PI;
    ncoef_min = MIN(ncoef_l_use,ncoef_l);
    for(icount = 1; icount<= ncoef_min ; icount++){
       aka = (double)kastore[(icount+icoef_off)];
       akb = (double)kbstore[(icount+icoef_off)];
       akc = (double)kcstore[(icount+icoef_off)];
       xk = (aka*hmati[1] +akb*hmati[2] +akc*hmati[3])*tpi;
       yk = (aka*hmati[4] +akb*hmati[5] +akc*hmati[6])*tpi;
       zk = (aka*hmati[7] +akb*hmati[8] +akc*hmati[9])*tpi;
       grhor_x[icount] =  xk*rhoci[icount];
       grhor_y[icount] =  yk*rhoci[icount];
       grhor_z[icount] =  zk*rhoci[icount];
       grhoi_x[icount] = -xk*rhocr[icount];
       grhoi_y[icount] = -yk*rhocr[icount];
       grhoi_z[icount] = -zk*rhocr[icount];
    }/*endfor*/
    if(ncoef_l_use != ncoef_l) {
      grhor_x[icount] = 0.0;
      grhor_y[icount] = 0.0;
      grhor_z[icount] = 0.0;
      grhoi_x[icount] = 0.0;
      grhoi_y[icount] = 0.0;
      grhoi_z[icount] = 0.0;
    }/* endif */

/*===================================================================*/
 }/*end routine*/
/*===================================================================*/


/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void get_igigrho(CPSCR *cpscr,EWALD *ewald,double *rhocr,double *rhoci,
                 double *hmati,int ncoef_l_use,int ncoef_l,int iopt,int icoef_off)
      
/*=============================================================================*/
{/*begin routine*/
/*=============================================================================*/

/*assign local pointers */
     double *grhor_x        =    cpscr->cpscr_grho.g_rhor_x;
     double *grhor_y        =    cpscr->cpscr_grho.g_rhor_y;
     double *grhor_z        =    cpscr->cpscr_grho.g_rhor_z;
     double *grhoi_x        =    cpscr->cpscr_grho.g_rhoi_x;
     double *grhoi_y        =    cpscr->cpscr_grho.g_rhoi_y;
     double *grhoi_z        =    cpscr->cpscr_grho.g_rhoi_z;
     double *g2_rhor        =    cpscr->cpscr_grho.g2_rhor;
     double *g2_rhoi        =    cpscr->cpscr_grho.g2_rhoi;
     int *kastore           =    ewald->kastr;
     int *kbstore           =    ewald->kbstr;
     int *kcstore           =    ewald->kcstr;


/* local variables */
      
      double  aka,akb,akc,xk,yk,zk,tpi,g2;
      int icount,ncoef_min;

/*===================================================================*/
/* I) Get ig times rho (gradient in g-space)                         */

  tpi = 2.0*M_PI;
  ncoef_min = MIN(ncoef_l_use,ncoef_l);
  switch(iopt){

   case 1: /* xx and yy components   */

    for(icount = 1; icount<= ncoef_min ; icount++){
       aka = (double)kastore[(icount+icoef_off)];
       akb = (double)kbstore[(icount+icoef_off)];
       akc = (double)kcstore[(icount+icoef_off)];
       xk = (aka*hmati[1] +akb*hmati[2] +akc*hmati[3])*tpi;
       yk = (aka*hmati[4] +akb*hmati[5] +akc*hmati[6])*tpi;
       zk = (aka*hmati[7] +akb*hmati[8] +akc*hmati[9])*tpi;
       grhor_x[icount] = -xk*xk*rhocr[icount];
       grhor_y[icount] = -yk*yk*rhocr[icount];
       grhoi_x[icount] = -xk*xk*rhoci[icount];
       grhoi_y[icount] = -yk*yk*rhoci[icount];
    }/*endfor*/
    if(ncoef_l_use != ncoef_l){
      grhor_x[icount] = 0.0;
      grhor_y[icount] = 0.0;
      grhoi_x[icount] = 0.0;
      grhoi_y[icount] = 0.0;
    }
    break;

   case 2: /* zz and xy components   */

    for(icount = 1; icount<= ncoef_min ; icount++){
       aka = (double)kastore[(icount+icoef_off)];
       akb = (double)kbstore[(icount+icoef_off)];
       akc = (double)kcstore[(icount+icoef_off)];
       xk = (aka*hmati[1] +akb*hmati[2] +akc*hmati[3])*tpi;
       yk = (aka*hmati[4] +akb*hmati[5] +akc*hmati[6])*tpi;
       zk = (aka*hmati[7] +akb*hmati[8] +akc*hmati[9])*tpi;
       grhor_x[icount] = -zk*zk*rhocr[icount];
       grhor_y[icount] = -xk*yk*rhocr[icount];
       grhoi_x[icount] = -zk*zk*rhoci[icount];
       grhoi_y[icount] = -xk*yk*rhoci[icount];
    }/*endfor*/
    if(ncoef_l_use != ncoef_l){
      grhor_x[icount] = 0.0;
      grhor_y[icount] = 0.0;
      grhoi_x[icount] = 0.0;
      grhoi_y[icount] = 0.0;
    }
    break;

   case 3: /* xz and yz components   */

    for(icount = 1; icount<= ncoef_min ; icount++){
       aka = (double)kastore[(icount+icoef_off)];
       akb = (double)kbstore[(icount+icoef_off)];
       akc = (double)kcstore[(icount+icoef_off)];
       xk = (aka*hmati[1] +akb*hmati[2] +akc*hmati[3])*tpi;
       yk = (aka*hmati[4] +akb*hmati[5] +akc*hmati[6])*tpi;
       zk = (aka*hmati[7] +akb*hmati[8] +akc*hmati[9])*tpi;
       grhor_x[icount] = -xk*zk*rhocr[icount];
       grhor_y[icount] = -yk*zk*rhocr[icount];
       grhoi_x[icount] = -xk*zk*rhoci[icount];
       grhoi_y[icount] = -yk*zk*rhoci[icount];
    }/*endfor*/
    if(ncoef_l_use != ncoef_l){
      grhor_x[icount] = 0.0;
      grhor_y[icount] = 0.0;
      grhoi_x[icount] = 0.0;
      grhoi_y[icount] = 0.0;
    }
    break;

  }/* end switch */

/*===================================================================*/
 }/*end routine*/
/*===================================================================*/




/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void get_g2rho(CPSCR *cpscr,int *kastore,int *kbstore,int *kcstore,
               double *rhocr,double *rhoci,double *hmati,
               int ncoef_l_use,int ncoef_l,int icoef_off)

      
/*=============================================================================*/
{/*begin routine*/
/*=============================================================================*/

/*assign local pointers */
     double *g2_rhor        =    cpscr->cpscr_grho.g2_rhor;
     double *g2_rhoi        =    cpscr->cpscr_grho.g2_rhoi;


/* local variables */
      
    double  aka,akb,akc,xk,yk,zk,tpi,g2;
    int icount,ncoef_min;

/*===================================================================*/
/* I) Get g^2 times rho (gradient in g-space)                         */


    tpi = 2.0*M_PI;
    ncoef_min = MIN(ncoef_l_use,ncoef_l);
    for(icount = 1; icount<= ncoef_min ; icount++){
       aka = (double)kastore[(icount+icoef_off)];
       akb = (double)kbstore[(icount+icoef_off)];
       akc = (double)kcstore[(icount+icoef_off)];
       xk = (aka*hmati[1] +akb*hmati[2] +akc*hmati[3])*tpi;
       yk = (aka*hmati[4] +akb*hmati[5] +akc*hmati[6])*tpi;
       zk = (aka*hmati[7] +akb*hmati[8] +akc*hmati[9])*tpi;
       g2 = xk*xk + yk*yk + zk*zk;
       g2_rhor[icount] = -g2*rhocr[icount];
       g2_rhoi[icount] = -g2*rhoci[icount];
    }/*endfor*/
    if(ncoef_l_use != ncoef_l){
      g2_rhor[icount] = 0.0;
      g2_rhoi[icount] = 0.0;
    }/* endif */

/*===================================================================*/
 }/*end routine*/
/*===================================================================*/



