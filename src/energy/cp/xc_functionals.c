/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
/*                                                                          */
/*                         PI_MD:                                           */
/*             The future of simulation technology                          */
/*             ------------------------------------                         */
/*                   Module: xc_functionals.c                               */
/* Contains all the exchange and correlation functionals                    */
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

#define LYP_LAP_OFF


/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
/*  this subroutine calculates the Goedecker-Hutter Pade lda   */
/*  exchange correlation function which approximates the Perdew-Wang 1992   */
/*  functional */

void  excpot_pade_lda(double *v_ks,double *rho,double *exc_ret,double *muxc_ret,
                      int nfft,int nfft_proc,double vol,
                      int cp_lyp,int cp_lypm1,double *pvten_cp,int cp_ptens_calc)

/*==========================================================================*/
/*         Begin Routine                                                    */
 {/*Begin Routine*/
/*=======================================================================*/ 
/* Static parameters */

  static double a0=0.4581652932831429;
  static double a1=2.217058676663745;
  static double a2=0.7405551735357053;
  static double a3=0.01968227878617998;
  static double b1=1.0;
  static double b2=4.504130959426697;
  static double b3=1.110667363742916;
  static double b4=0.02359291751427506;
  static double ot;

/*============================================================================*/
/* Other local variables */

  int i;
  int nfft2 = nfft/2;
  int nfft2_proc = nfft_proc/2;
  double rho_r;
  double num,denom;
  double num_p,denom_p;
  double rs,drs_dn,rs3;
  double exc_now,vxc,pvxc;
  double vscale,dfact;
  
/*============================================================================*/
/* Evaluate the functional */

   ot = 1.0/3.0;
   exc_now = 0.0;
   pvxc   = 0.0;
   vxc    = 0.0;

   for(i=1 ; i<= nfft2_proc; i++){

     rho_r = rho[i]; 
     rs3 = 3.0/(4.0*M_PI*rho_r);
     rs = pow(rs3,ot);

     num = ((a3*rs + a2)*rs + a1)*rs +a0;
     denom = (((b4*rs + b3)*rs + b2)*rs + b1)*rs;
     exc_now -= (num/denom)*rho_r;


/*============================================================================*/
/* Evaluate the derivative of the functional */

    drs_dn = -rs/(3.0*rho_r);
    num_p = (3.0*a3*rs + 2.0*a2)*rs + a1;
    denom_p = ((4.0*b4*rs + 3.0*b3)*rs + 2.0*b2)*rs + b1;
    dfact = ((denom*num_p - denom_p*num)/(denom*denom))*drs_dn;
    v_ks[i] -= dfact*rho_r +  num/denom;
  
  }/* endfor grid loop (i) */

  vscale = vol/((double) nfft2);
  *exc_ret += exc_now*vscale;

/* ADD THE PRESSURE TENSOR CALCULATION TO ME!!!!! */

/*==========================================================================*/
}/*end routine*/
/*==========================================================================*/


/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
/*  this subroutine calculates the Goedecker-Hutter Pade lsda   */
/*  exchange correlation function which approximates the Perdew-Wang 1992   */
/*  functional */

void excpot_pade_lsda(double *v_ks_up,double *v_ks_dn,
                      double *rho_up,double *rho_dn,
                      double *exc_ret,double *muxc_ret,int nfft,int nfft_proc,
                      double vol,int cp_lyp,int cp_lypm1,
                      double *pvten_cp,int cp_ptens_calc)

/*============================================================================*/
{/* Begin function */
/*============================================================================*/
/* Static parameters */

  static double a0  = 0.4581652932831429;
  static double da0 = 0.119086804055547;
  static double a1  = 2.217058676663745;
  static double da1 = 0.6157402568883345;
  static double a2  = 0.7405551735357053;
  static double da2 = 0.1574201515892867;
  static double a3  = 0.01968227878617998;
  static double da3 = 0.003532336663397157;
  static double b1  = 1.0;
  static double db1 = 0.0;
  static double b2  = 4.504130959426697;
  static double db2 = 0.2673612973836267;
  static double b3  = 1.110667363742916;
  static double db3 = 0.2052004607777787;
  static double b4  = 0.02359291751427506;
  static double db4 = 0.004200005045691381;
  static double ot,ft;
  static double fx_denom;

/*============================================================================*/
/* Other local variables */

/*============================================================================*/
/* Other local variables */

  int i;
  int nfft2 = nfft/2;
  int nfft2_proc = nfft_proc/2;
  double exc_now = 0.0;
  double pvxc = 0.0;
  double vxc = 0.0;
  double n,n2;
  double num,denom;
  double num_p,denom_p;
  double dfact,vscale;
  double rs,drs_dn_up,drs_dn_dn,rs3;
  double fx,dfx_dzeta,dfx_dn_up,dfx_dn_dn;
  double zeta,dzeta_dn_up,dzeta_dn_dn;
  double pterm,mterm,dpterm,dmterm;
  double za0,za1,za2,za3,zb1,zb2,zb3,zb4;
  double dza0_up,dza1_up,dza2_up,dza3_up,dzb1_up,dzb2_up,dzb3_up,dzb4_up;
  double dza0_dn,dza1_dn,dza2_dn,dza3_dn,dzb1_dn,dzb2_dn,dzb3_dn,dzb4_dn;
  

/*============================================================================*/
/* Evaluate the functional */

  ot = 1.0/3.0;
  ft = 4.0/3.0;
  fx_denom = 2.0*(pow(2.0,ot)-1.0);

  for(i=1; i<=nfft2_proc; i++){

    n = rho_up[i] + rho_dn[i];
    rs3 = 3.0/(4.0*M_PI*n);
    rs = pow(rs3,ot);
  
    zeta = (rho_up[i]-rho_dn[i])/n;
    pterm = pow((1.0+zeta),ft);
    mterm = pow((1.0-zeta),ft);
    fx = (pterm + mterm - 2.0)/fx_denom;

    za0 = a0 + da0*fx;
    za1 = a1 + da1*fx;
    za2 = a2 + da2*fx;
    za3 = a3 + da3*fx;

    zb1 = b1 + db1*fx;
    zb2 = b2 + db2*fx;
    zb3 = b3 + db3*fx;
    zb4 = b4 + db4*fx;

    num = ((za3*rs + za2)*rs + za1)*rs +za0;
    denom = (((zb4*rs + zb3)*rs + zb2)*rs + zb1)*rs;
    exc_now -= (num/denom)*n;

/*============================================================================*/
/* Evaluate the derivative of the functional with respect to up and down densities*/

    drs_dn_up = drs_dn_dn = -rs/(3.0*n);

    n2 = n*n;
    dzeta_dn_up = 2.0*rho_dn[i]/n2;
    dzeta_dn_dn = -2.0*rho_up[i]/n2;

    dpterm = ft*pow((1.0+zeta),ot);
    dmterm = ft*pow((1.0-zeta),ot);
    dfx_dzeta = (dpterm - dmterm)/fx_denom;
    dfx_dn_up = dfx_dzeta*dzeta_dn_up;
    dfx_dn_dn = dfx_dzeta*dzeta_dn_dn;

    dza0_up = da0*dfx_dn_up;
    dza1_up = da1*dfx_dn_up;
    dza2_up = da2*dfx_dn_up;
    dza3_up = da3*dfx_dn_up;
    dzb1_up = db1*dfx_dn_up;
    dzb2_up = db2*dfx_dn_up;
    dzb3_up = db3*dfx_dn_up;
    dzb4_up = db4*dfx_dn_up;

    dza0_dn = da0*dfx_dn_dn;
    dza1_dn = da1*dfx_dn_dn;
    dza2_dn = da2*dfx_dn_dn;
    dza3_dn = da3*dfx_dn_dn;
    dzb1_dn = db1*dfx_dn_dn;
    dzb2_dn = db2*dfx_dn_dn;
    dzb3_dn = db3*dfx_dn_dn;
    dzb4_dn = db4*dfx_dn_dn;

    num_p = ((3.0*za3*rs + 2.0*za2)*rs + za1)*drs_dn_up + 
            ((dza3_up*rs + dza2_up)*rs + dza1_up)*rs + dza0_up;
    denom_p = (((4.0*zb4*rs + 3.0*zb3)*rs + 2.0*zb2)*rs + zb1)*drs_dn_up + 
              (((dzb4_up*rs + dzb3_up)*rs + dzb2_up)*rs + dzb1_up)*rs;
 
    dfact = ((denom*num_p - denom_p*num)/(denom*denom));

    v_ks_up[i] -= num/denom + dfact*n;

    num_p = ((3.0*za3*rs + 2.0*za2)*rs + za1)*drs_dn_dn + 
            ((dza3_dn*rs + dza2_dn)*rs + dza1_dn)*rs + dza0_dn;
    denom_p = (((4.0*zb4*rs + 3.0*zb3)*rs + 2.0*zb2)*rs + zb1)*drs_dn_dn + 
              (((dzb4_dn*rs + dzb3_dn)*rs + dzb2_dn)*rs + dzb1_dn)*rs;

    dfact = ((denom*num_p - denom_p*num)/(denom*denom));

    v_ks_dn[i] -= num/denom + dfact*n;

  }/* endfor loop over grid */

  vscale = vol/((double) nfft2);
  *exc_ret += exc_now*vscale;

/* ADD THE PRESSURE TENSOR CALCULATION TO ME!!!! */

/*==========================================================================*/
}/*end function*/
/*==========================================================================*/


/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
/*  this subroutine calculates the perdew-zunger lda   */
/*  exchange correlation function form the old prb     */

void  excpot_pz_lda(double *v_ks,double *rho,double *exc_ret,double *muxc_ret,
                    int nfft,int nfft_proc,double vol,
                    int cp_lyp,int cp_lypm1,double *pvten_cp,int cp_ptens_calc)

/*==========================================================================*/
/*         Begin Routine                                                    */
 {/*Begin Routine*/
/*=======================================================================*/


/* Static variables               */
/* lda correlation parameters */
  static double gamma = -0.14230;
  static double beta1 = 1.05290;
  static double beta2 = 0.33340;
  static double au    = 0.03110;
  static double bu    = -0.0480;
  static double cu    = 0.00200;
  static double du    = -0.01160;

/* Local variables */

  double pvx,pvc,vscale;
  double pi,rho_r,power,xfact,cfact,rs,fpi,fpin;
  double ex,ec,mufact,rat1,rat2,cf1,cf2,sqtrs,vxc;
  double lnrs;
  double lyp_fact;
  int    i,kk,iii;
  int nfft2 = nfft/2;
  int nfft2_proc = nfft_proc/2;

/*=================================================================*/
/* I) Assign use constants                                         */

   pi     =   M_PI;
   fpi    = 4.0*pi;
   rat1   = (7.0/6.0)*beta1;
   rat2   = (4.0/3.0)*beta2;
   power    = 1.0/3.0;
   lyp_fact = (double)(1-(cp_lyp || cp_lypm1));
   ex     = 0.0;
   ec     = 0.0;
   pvx    = 0.0;
   pvc    = 0.0;
   vxc    = 0.0;

   for(i=1 ; i<= nfft2_proc; i++){
     rho_r = rho[i]; 
     fpin = fpi*rho_r;
     rs   = pow((3.0/fpin),power);
     sqtrs = sqrt(rs);

      xfact = pow(((3.0*rho_r/pi)),power);
      xfact = -xfact;
/* rs > 1  */
   if(rs >= 1.0){
     cf1 = 1.0 + beta1*sqtrs + beta2*rs;
     cf2 = 1.0 + rat1*sqtrs + rat2*rs;
     cfact = gamma/cf1;
      mufact = cfact*(cf2/cf1);
/*rs < 1   */
   }else{
    lnrs = log(rs);
    cfact = au*lnrs + bu + cu*rs*lnrs + du*rs;
    mufact = au*lnrs + (bu - power*au) + 2.0*power*cu*rs*lnrs
           + power*(2.0*du - cu)*rs;
   }/*endif*/
       
   cfact = cfact*lyp_fact;
   ex = ex + xfact*rho_r;
   ec = ec + cfact*rho_r;
   mufact = mufact*lyp_fact;

   v_ks[i] += (xfact + mufact);

   vxc      += (xfact + mufact)*rho_r;

   pvx +=  xfact*rho_r;
   pvc += (mufact - cfact)*rho_r;

  }/*endfor*/

   ex *= 0.750;
   pvx *= 0.250;

   vscale = vol/((double)(nfft2));
   *exc_ret += (ex+ec)*vscale;
   *muxc_ret += vxc*vscale;
   pvx *= vscale;
   pvc *= vscale;

   if(cp_ptens_calc==1){
      pvten_cp[1] += (pvx + pvc);
      pvten_cp[5] += (pvx + pvc);
      pvten_cp[9] += (pvx + pvc);
   }/*endif*/

/*==========================================================================*/
}/*end routine*/
/*==========================================================================*/

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void  excpot_pw_lda(double *v_ks,double *rho,double *exc_ret,double *muxc_ret,
                    int nfft,int nfft_proc,double vol,
                    int cp_lyp,int cp_lypm1,double *pvten_cp,int cp_ptens_calc)

/*==========================================================================*/
/*         Begin Routine                                                    */
   {/*Begin Routine*/
/*=======================================================================*/


/*this subroutine calculates the perdew-wang lda   */
/*exchange correlation function from the 1992 prb  */

/*local variables  */

 double pi,rho_r,power,xfact,cfact,rs,fpi,fpin;
 double ex,ec,vxc,mufact,cf1,cf2,sqtrs,rat1,rat3,rat4,logterm;
 double tc1,tc2,aalph1,sixa,mtwoa,opars;
 double pvx,pvc,vscale;
 double lyp_fact;
 int    i,kk,iii;
 int nfft2 = nfft/2;
 int nfft2_proc = nfft_proc/2;


/*lda correlation parameters */
/* static variables          */
 static double a = 0.03109070;
 static double alpha1 = 0.213700;
 static double beta1  = 7.59570;
 static double beta2  = 3.58760;
 static double beta3  = 1.63820;
 static double beta4  = 0.492940;

/*=======================================================================*/
/* I) Assign useful constants                                            */

  pi =  M_PI;
  fpi = 4.0*pi;
  power = 1.0/3.0;
  rat1 = 0.50*beta1;
  rat3 = 1.50*beta3;
  rat4 = 2.0*beta4;
  aalph1 = a*alpha1;
  sixa = 6.0*a;
  mtwoa = -2.0*a;
  lyp_fact = (double)(1-(cp_lyp || cp_lypm1));
  ex = 0.0;
  ec = 0.0;
  vxc = 0.0;
  pvx = 0.0;
  pvc = 0.0;

/*=======================================================================*/
/* II) Loop over real space grid                                         */

  for(i=1; i<= nfft2_proc; i++){
    rho_r = rho[i];
    fpin = fpi*rho_r;
    rs = pow((3.0/fpin),power);
    sqtrs = sqrt(rs);
      
    xfact = pow((3.0*rho_r/pi),power);
    xfact = -xfact;

    cf1 = 2.0*a*(beta1*sqtrs + beta2*rs + 
                 beta3*pow(rs,1.50) + beta4*rs*rs);
    cf2 = 2.0*a*(rat1/sqtrs + beta2 + 
                 rat3*sqtrs + rat4*rs);
    logterm = log(1.0 + 1.0/cf1);
    opars = 1.0 + alpha1*rs;
    cfact = mtwoa*opars*logterm;
    tc1 = mtwoa*(1.0 + 2.0*alpha1*rs/3.0);
    tc2 = mtwoa*rs*opars/(3.0*cf1*(1.0 + cf1));
    mufact = tc1*logterm + tc2*cf2;

    cfact = cfact*lyp_fact;
    ex +=  xfact*rho_r;
    ec +=  cfact*rho_r;
    mufact *= lyp_fact;

    pvx +=  xfact*rho_r;
    pvc +=  ((mufact - cfact)*rho_r);
    v_ks[i] += (xfact + mufact);
    vxc += (xfact + mufact)*rho_r;

  }/*endfor*/

/*=======================================================================*/
/* III) Pressure volume stuff                                            */

   ex *= 0.750;
   pvx *= 0.250;
   vscale = vol/((double)(nfft2));
   *exc_ret += (ex+ec)*vscale;
   *muxc_ret += vxc*vscale;
   pvx *= vscale;
   pvc *= vscale;

   if(cp_ptens_calc==1){
     pvten_cp[1] += (pvx + pvc);
     pvten_cp[5] += (pvx + pvc);
     pvten_cp[9] += (pvx + pvc);
   }/*endif*/

/*==========================================================================*/
 }/*end routine*/
/*==========================================================================*/

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void excpot_pz_lsda(double *v_ks_up,double *v_ks_dn,
                     double *rho_up,double *rho_dn,
                     double *exc_ret,double *muxc_ret,int nfft,int nfft_proc,
                     double vol,int cp_lyp,int cp_lypm1,
                     double *pvten_cp,int cp_ptens_calc)

/*==========================================================================*/
/*         Begin Routine                                                    */
   {/*Begin Routine*/
/*=======================================================================*/
/*this subroutine calculates the perdew-zunger lsda   */
/*exchange correlation function form the old prb      */

/*local variables */
   double pi,rho_r,rho_ru,rho_rd,power,exc;
   double xfact,rs,fpi,fpin,zeta,xff,xfp,muxf,muxp,rat1,rat2;
   double ex,vxc,muxfact_up,muxfact_dn,ec,mucf,mucp,mucfact_up;
   double mucfact_dn,cfact,cff,cfp;
   double cff1,cfp1,cff2,cfp2,ratbf1,ratbf2;
   double ratbp1,ratbp2;
   double fzeta,fpzeta,opzeta,omzeta,rden,sqtrs;
   double pvx,pvc,vscale;
   double lyp_fact;
   int i,kk;
   int nfft2 = nfft/2;
   int nfft2_proc = nfft_proc/2;

/*lsda correlation parameters   */
/* Static variables             */
   static double cf = 0.57730;
   static double cp = 0.45820;
   static double gammaf = -0.08430;
   static double gammap = -0.14230;
   static double beta1f = 1.05290;
   static double beta1p = 1.39810;
   static double beta2f = 0.33340;
   static double beta2p = 0.26110;

/*=======================================================================*/
/* I) Assign useful constants                                            */

   pi = M_PI;
   fpi = 4.0*pi;
   rat1 = 2.0/3.0;
   rat2 = 4.0/3.0;
   ratbf1 = (7.0/6.0)*beta1f;
   ratbf2 = (4.0/3.0)*beta2f;
   ratbp1 = (7.0/6.0)*beta1p;
   ratbp2 = (4.0/3.0)*beta2p;
   power = 1.0/3.0;
   rden = 1.0/( pow(2.0,rat2) - 2.0);
   lyp_fact = (double)(1-(cp_lyp || cp_lypm1));

   ex = 0.0;
   ec = 0.0;
   vxc = 0.0;
   pvx = 0.0;
   pvc = 0.0;

/*=======================================================================*/
/* II) Loop over real space grid                                         */

   for(i=1; i<= nfft2_proc; i++){
     rho_ru = rho_up[i];
     rho_rd = rho_dn[i];
     rho_r = rho_ru + rho_rd;
     zeta = (rho_ru - rho_rd)/rho_r;
     opzeta = 1.0 + zeta;
     omzeta = 1.0 - zeta;
     fpin = fpi*rho_r;
     rs = pow((3.0/fpin),power);
     sqtrs = sqrt(rs);
     xff = -cf/rs;
     xfp = -cp/rs;
     muxf = rat2*xff;
     muxp = rat2*xfp;
     cff1 = 1.0 + beta1f*sqtrs + beta2f*rs;
     cff2 = 1.0 + ratbf1*sqtrs + ratbf2*rs;
     cfp1 = 1.0 + beta1p*sqtrs + beta2p*rs;
     cfp2 = 1.0 + ratbp1*sqtrs + ratbp2*rs;
     cff = gammaf/cff1;
     cfp = gammap/cfp1;
     mucf = cff*(cff2/cff1);
     mucp = cfp*(cfp2/cfp1);
     fzeta = rden*(pow(opzeta,rat2) + pow(omzeta,rat2) - 2.0);
     fpzeta = rat2*rden*(pow(opzeta,power) - pow(omzeta,power));
     xfact = xfp + fzeta*(xff - xfp);
     cfact = cfp + fzeta*(cff - cfp);

     muxfact_up = muxp + fzeta*(muxf-muxp) + 
                     omzeta*fpzeta*(xff - xfp);
     muxfact_dn = muxp + fzeta*(muxf-muxp) - 
                     opzeta*fpzeta*(xff - xfp);

     mucfact_up = mucp + fzeta*(mucf-mucp) + 
                     omzeta*fpzeta*(cff - cfp);
     mucfact_dn = mucp + fzeta*(mucf-mucp) - 
                     opzeta*fpzeta*(cff - cfp);

     cfact = cfact*lyp_fact;
     ex += xfact*rho_r;
     ec += cfact*rho_r;
     mucfact_up *= lyp_fact;
     mucfact_dn *= lyp_fact;

     pvx +=  (muxfact_up*rho_ru + muxfact_dn*rho_rd - xfact*rho_r);
     pvc +=  (mucfact_up*rho_ru + mucfact_dn*rho_rd - cfact*rho_r);

     v_ks_up[i] += (muxfact_up + mucfact_up);
     v_ks_dn[i] += (muxfact_dn + mucfact_dn);
     vxc += (muxfact_up + mucfact_up)*rho_ru + (muxfact_dn + mucfact_dn)*rho_rd;
   }/*endfor*/

/*=======================================================================*/
/* III) Pressure/volume stuff                                           */
    vscale = vol/((double)(nfft2));
    *exc_ret += (ex+ec)*vscale;
    *muxc_ret += vxc*vscale;
    pvx *= vscale;
    pvc *= vscale;

    if(cp_ptens_calc==1){
     pvten_cp[1] += (pvx + pvc);
     pvten_cp[5] += (pvx + pvc);
     pvten_cp[9] += (pvx + pvc);
    }/*endif*/

/*=======================================================================*/
}/*end routine*/
/*=======================================================================*/

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void excpot_pw_lsda(double *v_ks_up,double *v_ks_dn,
                     double *rho_up,double *rho_dn,
                     double *exc_ret,double *muxc_ret,int nfft,int nfft_proc,
                     double vol,int cp_lyp,int cp_lypm1,
                     double *pvten_cp,int cp_ptens_calc)

/*==========================================================================*/
/*         Begin Routine                                                    */
   {/*Begin Routine*/
/*=======================================================================*/
/*this subroutine calculates the perdew-wang lsda   */
/*exchange correlation function from prb 1992 and prl 1996  */

/*local variables */
   double pi,rho_r,rho_ru,rho_rd;
   double aa,power,fth,kf,xfact_up,xfact_dn;
   double ex,vxc,muxfact_up,muxfact_dn,ec,mucfact_up,mucfact_dn;
   double cfact;
   double pvx,pvc,vscale;
   double lyp_fact;
   int i,kk;
   int nfft2 = nfft/2;
   int nfft2_proc = nfft_proc/2;

/*=======================================================================*/
/* I) Assign useful constants                                            */

   pi = M_PI;
   aa = 3.0*pi*pi;
   power = 1.0/3.0;
   fth = 4.0/3.0;
   lyp_fact = (double)(1-(cp_lyp || cp_lypm1));

   ex = 0.0;
   ec = 0.0;
   vxc = 0.0;
   pvx = 0.0;
   pvc = 0.0;

/*=======================================================================*/
/* II) Loop over real space grid                                         */

   for(i=1; i<= nfft2_proc; i++){
     rho_ru = rho_up[i];
     rho_rd = rho_dn[i];
     rho_r = rho_ru + rho_rd;

     /* Do the exchange part */
     /* up density => 2*rho_up in spin unpolarized exchange */
     kf = pow(aa*2.0*rho_ru,power);
     xfact_up = -(0.75/pi)*kf; 
     muxfact_up = fth*xfact_up;

     /* dn density => 2*rho_dn in spin unpolarized exchange */
     kf = pow(aa*2.0*rho_rd,power);
     xfact_dn = -(0.75/pi)*kf; 
     muxfact_dn = fth*xfact_dn;    

     ex += (rho_ru*xfact_up + rho_rd*xfact_dn);
     
     /* Do the correlation part */

     pw_c_lsda(rho_ru,rho_rd,&cfact,&mucfact_up,&mucfact_dn);

     cfact = cfact*lyp_fact;
     ec += cfact*rho_r;
     mucfact_up *= lyp_fact;
     mucfact_dn *= lyp_fact;

     /* I'm not sure the formula for the pressure is correct */
     pvx +=  ((muxfact_up-xfact_up)*rho_ru + (muxfact_dn-xfact_dn)*rho_rd);
     pvc +=  (mucfact_up*rho_ru + mucfact_dn*rho_rd - cfact*rho_r);

     v_ks_up[i] += (muxfact_up + mucfact_up);
     v_ks_dn[i] += (muxfact_dn + mucfact_dn);
     vxc += (muxfact_up + mucfact_up)*rho_ru 
          + (muxfact_dn + mucfact_dn)*rho_rd;
   }/*endfor*/

/*=======================================================================*/
/* III) Pressure/volume stuff                                           */
    vscale = vol/((double)(nfft2));
    *exc_ret += (ex+ec)*vscale;
    *muxc_ret += vxc*vscale;
    pvx *= vscale;
    pvc *= vscale;

    if(cp_ptens_calc==1){
     pvten_cp[1] += (pvx + pvc);
     pvten_cp[5] += (pvx + pvc);
     pvten_cp[9] += (pvx + pvc);
    }/*endif*/

/*=======================================================================*/
}/*end routine*/
/*=======================================================================*/

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void pw_c_lsda(double rho_up,double rho_dn,double *ec_ret,double *dfxc_dn_c_up,
	       double *dfxc_dn_c_dn)

/*==========================================================================*/
/*         Begin Routine                                                    */
   {/*Begin Routine*/
/*=======================================================================*/
/*this subroutine calculates the perdew-wang lsda   */
/*exchange correlation function from prb 1992 and prl 1996  */

/*local variables */
     double rho,rs,zeta,zeta3,zeta4,F,dF_dz;
     double ec,eu,ep,ac,ediff;
     double deu_drs,dep_drs,dac_drs,dec_drs,dec_dz,base;
     
     /* gam = 2^(4/3) - 2 */
     /* d2F_dz0 = d^2(F)/d_zeta [0] = 8/(9*gam) */ 
     static double pow1 = 1.0/3.0;
     static double fth  = 4.0/3.0;
     static double gam  = 0.5198420997897460;
     static double d2F_dz0 = 1.70992093416137;

     eu = 0.0;
     ep = 0.0;
     ac = 0.0;
     deu_drs=0.0;
     dep_drs=0.0;
     dac_drs=0.0;

     rho = rho_up+rho_dn;
     rs = pow((3.0/(4.0*M_PI*rho)),pow1);
     zeta = (rho_up-rho_dn)/rho;
     zeta3 = zeta*zeta*zeta;
     zeta4 = zeta3*zeta;
     F = (pow(1+zeta,fth)+pow(1-zeta,fth)-2.0)/gam;
     dF_dz = fth * (pow(1+zeta,pow1) - pow(1-zeta,pow1))/gam;

     /* spin unpolarized */
     pw_ec(0.0310907,0.21370,7.5957,3.5876,1.6382,0.49294, 
	   rs, &eu, &deu_drs);

     /* spin polarized */
     pw_ec(0.01554535,0.20548,14.1189,6.1977,3.3662,0.62517, 
	   rs, &ep, &dep_drs);

     /* spin stiffness */
     pw_ec(0.0168869,0.11125,10.357,3.6231,0.88026,0.49671,
	   rs, &ac, &dac_drs);
     /*returns minus of the spin stiffness */
     ac *= -1.0;
     dac_drs *= -1.0;

     /* Calculate the functional */     
     ediff = ep - eu;
     ec = eu + ac*F/d2F_dz0*(1.0-zeta4) + ediff*F*zeta4;

     /* Calculate the derivative of the functional wrt rs */     
     dec_drs = deu_drs*(1.0 - F*zeta4) + dep_drs*F*zeta4 
       + dac_drs*F*(1.0-zeta4)/d2F_dz0;

     /* Calculate the derivative of the functional wrt zeta */
     dec_dz = 4.0*zeta3*F*(ediff - ac/d2F_dz0)  
       + dF_dz*(zeta4*ediff + (1.0 - zeta4)*ac/d2F_dz0); 
  
     /* Calculate the derivative of the functional wrt density */
     base = ec - rs/3.0*dec_drs - zeta*dec_dz;
     *dfxc_dn_c_up = base + dec_dz;
     *dfxc_dn_c_dn = base - dec_dz;
     *ec_ret = ec;
     
/*=======================================================================*/
}/*end routine*/
/*=======================================================================*/

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void pw_ec(double a, double alpha1, double beta1, double beta2, double beta3, 
	   double beta4, double rs, double *ec, double *dec_drs)

/*==========================================================================*/
/*         Begin Routine                                                    */
{/*Begin Routine*/
/*=======================================================================*/
/*local variables */

  double rr,drr;
  double srs,rstt,rs2,rs4;
  double logterm;
  double preterm;

/*===========================================================================*/
/* calculate ec (spin neutral, spin polarized, or spin stiffness depending   */
/* on the contants that are passed in) within Perdue-Wang LDA (prb 1992).    */

  rs2 = rs*rs;
  rs4 = rs2*rs2;
  rstt = pow(rs,1.5);
  srs = sqrt(rs);
  rr = 2.0*a*(beta1*srs + beta2*rs + beta3*rstt + beta4*rs2);
  logterm = log(1.0 + 1.0/rr);
  preterm = 1.0 + alpha1*rs;
  *ec = -2.0*a*preterm*logterm;

  /* calculate the d(ec)/d(rs) */
  drr = 2.0*a*(0.5*beta1/srs + beta2 + 1.5*beta3*srs + 2.0*beta4*rs);
  *dec_drs = -2.0*a*(alpha1*logterm - preterm*drr/(rr*(1.0 + rr)));

/*=======================================================================*/
}/*end routine*/
/*=======================================================================*/


/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

 void becke_gcx_lsda(double rho,double g_rho2,
                     double *fn,double *df_dn,double *df_dgn,double beta)

/*==========================================================================*/
/*         Begin Routine                                                    */
   {/*Begin Routine*/
/*=======================================================================*/
/*local variables */
   double gmag;
   double rho43,rho13,x,root,asinch;
   double bterm,bx;
   double bpterm,bpx,bx_by_x;

/*   static double beta; */
   static double sixb;
   static double twelveb;
   static double fth;

   /*   beta = 0.0042; */
   sixb = 6.0*beta;
   twelveb = 12.0*beta;
   fth = 4.0/3.0;


/*=======================================================================*/
/* calculate b(x)   the Becke function */

   rho43 = pow(rho,(4.0/3.0));
   rho13 = pow(rho,(1.0/3.0));
   gmag = sqrt(g_rho2);
   x = gmag/rho43;
   root = sqrt(1.0 + x*x);
   asinch = log(x + root);
   bterm = 1.0 + 6.0*beta*x*asinch;
   bx = x*x/bterm;
   bx_by_x = x/bterm;

/*=======================================================================*/
/* calculate b'(x)  */

   bpterm = (asinch + x/root);
   bpx = 2.0*bx_by_x - sixb*(bx_by_x)*(bx_by_x)*bpterm;

/*=======================================================================*/
/* calculate function and its derivatives wrt */
/* density and the gradient of the density    */

   *fn = -beta*rho43*bx;
   *df_dn = -beta*(fth*rho13*(bx - x*bpx));
   *df_dgn = -beta*bpx;

/*==============================================================================*/
}/*end routine*/
/*==============================================================================*/



/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

 void becke_gcx_lda(double rho,double g_rho2,
                    double *fn,double *df_dn,double *df_dgn,double beta)

/*==========================================================================*/
/*         Begin Routine                                                    */
   {/*Begin Routine*/
/*=======================================================================*/
/*local variables */
   double gmag;
   double rho43,rho13,x,root,asinch;
   double bterm,bx;
   double bpterm,bpx,bx_by_x;

/*   static double beta; */
   static double sixb;
   static double twelveb;
   static double fth;

/*   beta = 0.0042; */
   sixb = 6.0*beta;
   twelveb = 12.0*beta;
   fth = 4.0/3.0;


/*=======================================================================*/
/* calculate b(x)   the Becke function */

   rho    *= 0.5;
   g_rho2 *= 0.25;
   
   rho43 = pow(rho,(4.0/3.0));
   rho13 = pow(rho,(1.0/3.0));
   gmag = sqrt(g_rho2);
   x = gmag/rho43;
   root = sqrt(1.0 + x*x);
   asinch = log(x + root);
   bterm = 1.0 + 6.0*beta*x*asinch;
   bx = x*x/bterm;
   bx_by_x = x/bterm;

/*=======================================================================*/
/* calculate b'(x)  */

   bpterm = (asinch + x/root);
   bpx = 2.0*bx_by_x - sixb*(bx_by_x)*(bx_by_x)*bpterm;

/*=======================================================================*/
/* calculate function and its derivatives wrt */
/* density and the gradient of the density    */

   *fn = -2.0*beta*rho43*bx;
   *df_dn = -beta*(fth*rho13*(bx - x*bpx));
   *df_dgn = -beta*bpx;

/*==============================================================================*/
}/*end routine*/
/*==============================================================================*/


/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

 void fila_1x_lsda(double rho,double g_rho2,
                   double *fn,double *df_dn,double *df_dgn)

/*==========================================================================*/

/*         Begin Routine                                                    */
   {/*Begin Routine*/
/*=======================================================================*/
/*local variables */

   double rho43,rho13,rho83,y,root,asinch,asinch2;
   double T1,T1_by_y,T1p,sterm,spterm;
   double brack;

   static double beta;
   static double beta2,nineb2;
   static double fth,nhlf;

   beta = 0.00293;
   beta2 = beta*beta;
   nineb2 = 9.0*beta2;
   fth = 4.0/3.0;
   nhlf = 9.0/2.0;


/*=======================================================================*/
/* calculate T_1(y)   the modified Becke function */

   rho43 = pow(rho,(4.0/3.0));
   rho13 = rho43/rho;
   rho83 = rho43*rho43;
   y = g_rho2/rho83;
   root = sqrt(1.0 + y*y);
   asinch = log(y + root);
   asinch2 = asinch*asinch;
   sterm = 1.0 + 9.0*beta2*y*asinch2;
   T1 = beta*y/sqrt(sterm);
   T1_by_y = T1/y;

/*=======================================================================*/
/* calculate T_1'  */

   brack = asinch + 2.0*y/root;
   spterm = 1.0 - nhlf*T1_by_y*T1*asinch*brack;
   T1p = T1_by_y*spterm;   

/*=======================================================================*/
/* calculate function and its derivatives wrt */
/* density and the gradient of the density    */

   *fn = -rho43*T1;
   *df_dn = -fth*rho13*(T1 - 2.0*y*T1p);
   *df_dgn = -2.0*T1p*sqrt(y);

/*==============================================================================*/
}/*end routine*/
/*==============================================================================*/


/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

 void fila_2x_lsda(double rho,double g_rho2,
                   double *fn,double *df_dn,double *df_dgn)

/*==========================================================================*/

/*         Begin Routine                                                    */
   {/*Begin Routine*/
/*=======================================================================*/
/*local variables */

   double g_rho;
   double rho43,rho13,rho83,y,root,asinch,asinch2;
   double T1,T1_by_y,T1p,T1b,sterm,spterm;
   double brack;
   double betaf,betaf2,nineb2;
   double betafp;
   double denom;

   static double beta0,beta1,beta2;
   static double fth,nhlf;

   beta0 = 0.002913644;
   beta1 = 0.000947417;
   beta2 = 2501.149;
   fth = 4.0/3.0;
   nhlf = 9.0/2.0;


/*=======================================================================*/
/* calculate T_1(y)   the modified Becke function */

   rho43 = pow(rho,(4.0/3.0));
   rho13 = rho43/rho;
   rho83 = rho43*rho43;
   y = g_rho2/rho83;
   g_rho = sqrt(g_rho2);
   root = sqrt(1.0 + y*y);
   asinch = log(y + root);
   asinch2 = asinch*asinch;
   denom = beta2 + g_rho2;
   betaf = beta0 + beta1*(g_rho2/denom);
   betaf2 = betaf*betaf;
   sterm = 1.0 + 9.0*betaf2*y*asinch2;
   T1 = betaf*y/sqrt(sterm);
   T1_by_y = T1/y;

/*=======================================================================*/
/* calculate dT_1/dy, dT_1/d\beta  */

   betafp = 2.0*beta1*(g_rho/denom - g_rho*g_rho2/(denom*denom));
   brack = asinch + 2.0*y/root;
   spterm = 1.0 - nhlf*T1_by_y*T1*asinch*brack;
   T1p = T1_by_y*spterm;   
   T1b = (T1/betaf)*(1.0 - 9.0*T1*T1*asinch2/y);

/*=======================================================================*/
/* calculate function and its derivatives wrt */
/* density and the gradient of the density    */

   *fn = -rho43*T1;
   *df_dn = -fth*rho13*(T1 - 2.0*y*T1p);
   *df_dgn = -2.0*T1p*sqrt(y) - rho43*T1b*betafp;

/*==============================================================================*/
}/*end routine*/
/*==============================================================================*/



/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

 void fila_1x_lda(double rho,double g_rho2,
                  double *fn,double *df_dn,double *df_dgn)

/*==========================================================================*/

/*         Begin Routine                                                    */
   {/*Begin Routine*/
/*=======================================================================*/
/*local variables */

   double rho43,rho13,rho83,y,root,asinch,asinch2;
   double T1,T1_by_y,T1p,sterm,spterm;
   double brack;

   static double beta;
   static double beta2,nineb2;
   static double fth,nhlf;

   beta = 0.00293;
   beta2 = beta*beta;
   nineb2 = 9.0*beta2;
   fth = 4.0/3.0;
   nhlf = 9.0/2.0;


/*=======================================================================*/
/* calculate T_1(y)   the modified Becke function */

   rho *= 0.5;
   g_rho2 *= 0.25;

   rho43 = pow(rho,(4.0/3.0));
   rho13 = rho43/rho;
   rho83 = rho43*rho43;
   y = g_rho2/rho83;
   root = sqrt(1.0 + y*y);
   asinch = log(y + root);
   asinch2 = asinch*asinch;
   sterm = 1.0 + 9.0*beta2*y*asinch2;
   T1 = beta*y/sqrt(sterm);
   T1_by_y = T1/y;

/*=======================================================================*/
/* calculate T_1'  */

   brack = asinch + 2.0*y/root;
   spterm = 1.0 - nhlf*T1_by_y*T1*asinch*brack;
   T1p = T1_by_y*spterm;   

/*=======================================================================*/
/* calculate function and its derivatives wrt */
/* density and the gradient of the density    */

   *fn = -2.0*rho43*T1;
   *df_dn = -fth*rho13*(T1 - 2.0*y*T1p);
   *df_dgn = -2.0*T1p*sqrt(y);

/*==============================================================================*/
}/*end routine*/
/*==============================================================================*/



/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

 void fila_2x_lda(double rho,double g_rho2,
                  double *fn,double *df_dn,double *df_dgn)

/*==========================================================================*/

/*         Begin Routine                                                    */
   {/*Begin Routine*/
/*=======================================================================*/
/*local variables */

   double g_rho;
   double rho43,rho13,rho83,y,root,asinch,asinch2;
   double T1,T1_by_y,T1p,T1b,sterm,spterm;
   double brack;
   double betaf,betaf2,nineb2;
   double betafp;
   double denom;

   static double beta0,beta1,beta2;
   static double fth,nhlf;

   beta0 = 0.002913644;
   beta1 = 0.000947417;
   beta2 = 2501.149;
   fth = 4.0/3.0;
   nhlf = 9.0/2.0;


/*=======================================================================*/
/* calculate T_1(y)   the modified Becke function */

   rho *= 0.5;
   g_rho2 *= 0.25; 

   rho43 = pow(rho,(4.0/3.0));
   rho13 = rho43/rho;
   rho83 = rho43*rho43;
   y = g_rho2/rho83;
   g_rho = sqrt(g_rho2);
   root = sqrt(1.0 + y*y);
   asinch = log(y + root);
   asinch2 = asinch*asinch;
   denom = beta2 + g_rho2;
   betaf = beta0 + beta1*(g_rho2/denom);
   betaf2 = betaf*betaf;
   sterm = 1.0 + 9.0*betaf2*y*asinch2;
   T1 = betaf*y/sqrt(sterm);
   T1_by_y = T1/y;

/*=======================================================================*/
/* calculate dT_1/dy, dT_1/d\beta  */

   betafp = 2.0*beta1*(g_rho/denom - g_rho*g_rho2/(denom*denom));
   brack = asinch + 2.0*y/root;
   spterm = 1.0 - nhlf*T1_by_y*T1*asinch*brack;
   T1p = T1_by_y*spterm;   
   T1b = (T1/betaf)*(1.0 - 9.0*T1*T1*asinch2/y);

/*=======================================================================*/
/* calculate function and its derivatives wrt */
/* density and the gradient of the density    */

   *fn = -2.0*rho43*T1;
   *df_dn = -fth*rho13*(T1 - 2.0*y*T1p);
   *df_dgn = -2.0*T1p*sqrt(y) - rho43*T1b*betafp;

/*==============================================================================*/
}/*end routine*/
/*==============================================================================*/



/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void  pw91_gcx(double rho,double g_rho2,double *fn,
               double *df_dn,double *df_dgn)

/*==========================================================================*/
/*         Begin Routine                                                    */
{/*Begin Routine*/
/*=======================================================================*/
/*local variables */
  double drho;
  double pi,ax,xkf,s,ds_dn,ds_dgn,root;
  double sinch,ex_cent,exp_term,xnum,denom,rs;
  double xnum_p,denom_p,rps;
  double rhop4,a2s,s2;
  int iii;

/* static variables */
  static double a1 = 0.196450;
  static double a2 = 7.79560;
  static double a3 = 0.27430;
  static double a4 = -0.15080;
  static double a5 = 0.0040;
  static double power;
  static double pow4;

/*=======================================================================*/
/* I) Use parts of functional                                            */

   power = 1.0/3.0;
   pow4 = 4.0/3.0;
   drho = sqrt(g_rho2);
   pi = M_PI;
   ax = -0.750*pow((3.0/pi),power);

   xkf = pow((3.0*pi*pi*rho),power);
   s = drho/(2.0*xkf*rho);
   ds_dn = -2.0*drho/(3.0*xkf*rho*rho);
   ds_dgn = 1.0/(2.0*xkf*rho);

   a2s = a2*s;
   s2 = s*s;
   root = sqrt(1.0 + a2s*a2s);
   sinch = a1*log(a2s + root);
   ex_cent = exp(-100.0*s2);
   exp_term = s2*(a3 + a4*ex_cent);
   xnum = 1.0 + sinch*s + exp_term;
   denom = 1.0 + sinch*s + a5*s2*s2;
   rs = xnum/denom;

/*=======================================================================*/
/*  II) calculate the PW91x function                                     */

   rhop4 = ax*pow(rho,pow4);
   *fn = rhop4*(rs - 1.0);
   xnum_p = sinch + a1*a2*s/root + 2.0*s*a3 
          + 2.0*s*a4*(1.0 - 100.0*s2)*ex_cent;
   denom_p = sinch + a1*a2*s/root + 4.0*a5*s2*s;

/*=======================================================================*/
/* III) calculate the derivative of the function wrt density             */

   rps = xnum_p/denom - xnum*denom_p/(denom*denom);
   *df_dn = pow4*ax*pow(rho,power)*(rs - 1.0)
          + rhop4*rps*ds_dn;
 
/*=======================================================================*/
/* IV) calculate the derivative of the function wrt grad density         */

   *df_dgn = rhop4*rps*ds_dgn;

/*==============================================================================*/
}/*endroutine*/
/*==========================================================================*/


/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void pbe_gcx(double rho,double g_rho2,double *fn,double *df_dn,double *df_dgn,
	     double mu,double kappa)

/*=======================================================================*/
/* Perdew-Burke-Ernzerhof generalized gradient approximation to the      */
/* density functional for exchange-correlation energy of a many-electron */
/* system. See PRL 77 (1996) 3865 and White & Bird PRB 50 (1994) 4954.   */
/* Change of parameters -> different functionals:                        */
/* revPBE see Zhang and Wang, PRL 80 (1998) 890.  (kappa differs)        */
/* xPBE   see Xu and Goddard, JCP 121 (2004) 4068.(kappa and mu differ)  */
{/*Begin Routine*/
/*=======================================================================*/
/*local variables */

 double s,fx,dfx_dn,dfx_dgn;
 double denom,denom2;
 double gmag;
 double kf,cfact;

/* Constants */

 static double pow1 = 1.0/3.0;
 static double aa = 3.0*M_PI*M_PI;
 /*static double mu = 0.21951;
   static double kappa = 0.804;*/
 static double fth = 4.0/3.0;
 
/*===========================================================================*/
/* Calculate the F_X function */
/* Note: The lda xchange part is already included in excpot_pw_lda => drop 1 */
/* from fx. */ 

  gmag = sqrt(g_rho2);
  kf = pow(aa*rho,pow1);
  cfact = -(0.75/M_PI)*kf; 
  s = gmag/(2.0*kf*rho);
  denom = 1.0 + mu*s*s/kappa;
  fx = kappa - kappa/denom;

/*==========================================================================*/
/* Calculate the functional */

  *fn = rho*cfact*fx;

/*==========================================================================*/
/* Calculate the derivative of the F_X function wrt density */

  denom2 = denom*denom;
  dfx_dn = -8.0*mu*s*s/(3.0*denom2*rho);

/*==========================================================================*/
/* Calculate the derivative of the functional wrt density */

  *df_dn = fth*cfact*fx + cfact*rho*dfx_dn;

/*==========================================================================*/
/* Calculate the derivative of the F_X function wrt |grad n| */

  dfx_dgn = mu*s/(denom2*kf*rho);

/*==========================================================================*/
/* Calculate the derivative of the functional wrt |grad n| */

 *df_dgn = cfact*rho*dfx_dgn;

/*==========================================================================*/
}/*end routine*/
/*==========================================================================*/

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void rpbe_gcx(double rho,double g_rho2,double *fn,double *df_dn,double *df_dgn)

/*=======================================================================*/
/* Same as pbe, with a different form for F                              */
/*  Hammer, Hansen and Norskov PRB 59 (1999) 7413.                       */
{/*Begin Routine*/
/*=======================================================================*/
/*local variables */

 double s,fx,dfx_dn,dfx_dgn;
 double P0,dfx_ds,ds_dn,ds_dgn;
 double gmag;
 double kf,cfact;

/* Constants */

 static double pow1 = 1.0/3.0;
 static double aa = 3.0*M_PI*M_PI;
 static double mu = 0.21951;
 static double kappa = 0.8040;
 static double fth = 4.0/3.0;

/*===========================================================================*/
/* Calculate the F_X function */
/* Note: The lda xchange part is already included in excpot_pw_lda => drop 1 */
/* from fx. */ 

  gmag = sqrt(g_rho2);
  kf = pow(aa*rho,pow1);
  cfact = -(0.75/M_PI)*kf; 
  s = gmag/(2.0*kf*rho);
  P0 = exp(-mu*s*s/kappa);
  fx = kappa*(1.0-P0);

/*==========================================================================*/
/* Calculate the functional */

  *fn = rho*cfact*fx;

/*==========================================================================*/
/* Calculate the derivative of the F_X function wrt density */

  dfx_ds = 2.0*mu*s*P0;
  ds_dn  = -fth*s/rho;
  dfx_dn = dfx_ds*ds_dn;

/*==========================================================================*/
/* Calculate the derivative of the functional wrt density */

  *df_dn = fth*cfact*fx + cfact*rho*dfx_dn;

/*==========================================================================*/
/* Calculate the derivative of the F_X function wrt |grad n| */

  ds_dgn  = 1.0/(2.0*kf*rho); 
  dfx_dgn = dfx_ds*ds_dgn;

/*==========================================================================*/
/* Calculate the derivative of the functional wrt |grad n| */

 *df_dgn = cfact*rho*dfx_dgn;

/*==========================================================================*/
}/*end routine*/
/*==========================================================================*/


/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void pw91_gcc(double rho,double g_rho2,double *fn,
              double *df_dn,double *df_dgn)

/*==========================================================================*/
/*         Begin Routine                                                    */
   {/*Begin Routine*/
/*=======================================================================*/
/* local variables */
  double rho2,drho;
  double fk,sk,s,t,ds_dn,dt_dn,rs,rs2,rs3;
  double cnum,cdenom,cn;
  double dcnum1,dcnum2,dcnum3,dc_dn;
  double denom_86,dln_term,ec_86,dr_num,dr_den;
  double deriv,dr_term,r_term,dec_dn_86;
  double ec_a,dec_a_dn,arg_a,ex_a,a_per;
  double pre,da_dn,t2,t4,a2,onum,oden;
  double omega,t3,ot1,ot2,pre_ot2,brack,do_dn;
  double hund,hln_term,c_exp,h_per,om_term;
  double poly,ex_term;
  double dc_term,dh_dn,do_term;
  double c_term,p_fun;

/* static variables */
  static double alpha = 0.090;
  static double beta =0.06672632120;
  static double c0 = 15.75590;
  static double c1 = 0.0035210;
  static double c1p = 0.0016670;
  static double c2p = 0.0025680;
  static double c3p = 0.0232660;
  static double c4p = 7.389e-6;
  static double c5p = 8.7230;
  static double c6p = 0.4720;
  static double c7p = 7.389e-2;
  static double a_86 = 0.03109070;
  static double alpha1_86 = 0.213700;
  static double beta1_86 = 7.59570;
  static double beta2_86 = 3.58760;
  static double beta3_86 = 1.63820;
  static double beta4_86 = 0.492940;
  static double pi = M_PI;
  static double pow3;
  static double pi2;
/*=======================================================================*/
/* I) Get square root of g_rho2                                          */

  pow3 = 1.0/3.0;
  pi2 = M_PI*M_PI;
  drho = sqrt(g_rho2);

/*=======================================================================*/
/* II) calculate s and t and their derivatives                           */

  fk = pow((3.0*pi2*rho),pow3);
  sk = sqrt(4.0*fk/pi);
  s = drho/(2.0*fk*rho);
  t = drho/(2.0*sk*rho);
  rho2 = rho*rho;
  ds_dn = -2.0*drho/(3.0*fk*rho2);
  dt_dn = -7.0*drho/(12.0*sk*rho2);

/*=======================================================================*/
/* III) calculate c(n) and its derivative                                */

  rs = pow((3.0/(4.0*pi*rho)),pow3);
  rs2 = rs*rs;
  rs3 = rs*rs2;
  cnum = c2p + c3p*rs + c4p*rs2;
  cdenom = 1.0 + c5p*rs + c6p*rs2 + c7p*rs3;
  cn = c1p + cnum/cdenom;
  dcnum1 = c3p + 2.0*c4p*rs;
  dcnum2 = c2p + c3p*rs + c4p*rs2;
  dcnum3 = c5p + 2.0*c6p*rs + 3.0*c7p*rs2 ;
  dc_dn = (dcnum1 - dcnum2*dcnum3/cdenom)/cdenom;
  dc_dn = -rs*dc_dn/(3.0*rho);

/*=======================================================================*/
/* IV) get the old pewdew-wang 1986 lda correlation                      */
/*     functional and its derivative                                     */

   denom_86 = beta1_86*sqrt(rs)  + beta2_86*rs +
              beta3_86*pow(rs,1.50) + beta4_86*rs*rs;
   denom_86 = denom_86*2.0*a_86;
   dln_term = log(1.0 + 1.0/denom_86);
   ec_86 = -2.0*rho*a_86*(1.0 + alpha1_86*rs)*dln_term;
   dr_num = a_86*rs*(1.0 + alpha1_86*rs);
   dr_den = denom_86*(1.0 + denom_86);
   deriv = 0.50*beta1_86/sqrt(rs) + beta2_86 +
           1.50*beta3_86*sqrt(rs) + 2.0*beta4_86*rs;
   deriv = deriv*2.0*a_86;
   dr_term = -2.0*(dr_num/dr_den)*(deriv/3.0);
   r_term = -2.0*a_86*(1.0 + 2.0*alpha1_86*rs/3.0)*dln_term;
   dec_dn_86 = r_term + dr_term;

/*=======================================================================*/
/* V) use ec_86 to compute exponential factor                            */

   ec_a = ec_86/rho;
   dec_a_dn = (dec_dn_86 - ec_86/rho)/rho;

/*=======================================================================*/
/*  VI) compute the 'a' functional and its derivative                    */

   arg_a = -2.0*alpha*ec_a/(beta*beta);
   ex_a = exp(arg_a) - 1.0;
   a_per = (2.0*alpha/beta)/ex_a;

   pre = 4.0*alpha*alpha/(beta*beta*beta);
   da_dn = (pre/(ex_a*ex_a))*dec_a_dn*exp(arg_a);

/*=======================================================================*/
/*  V) compute the omega functional                                      */
/*  (argument of the log term in perdew-wang language)                    */

   t2 = t*t;
   t4 = t2*t2;
   a2 = a_per*a_per;
   onum = t2 + a_per*t4;
   oden = 1.0 + a_per*t2 + a2*t4;
   omega = 1.0 + (2.0*alpha/beta)*(onum/oden);

   t3 = t2*t;
   ot1 = (2.0*t + 4.0*a_per*t3)*dt_dn + t4*da_dn;
   ot2 = (t2 + 2.0*a_per*t4)*da_dn +
         (2.0*a_per*t + 4.0*a2*t3)*dt_dn;
   pre_ot2 = onum/oden;
   ot2 = ot2*pre_ot2;
   brack = ot1 - ot2;
   do_dn = (2.0*alpha/beta)*brack/oden;

/*=======================================================================*/
/*  VI) compute the 'h' functional and its derivative                    */

   hund = -100.0*s*s;
   hln_term = (beta*beta/(2.0*alpha))*log(omega);
   c_exp = t2*exp(hund);
   h_per = c0*(cn - c1)*c_exp + hln_term;

   om_term = (beta*beta/(2.0*alpha))*do_dn/omega;
   poly = 2.0*t*dt_dn - 200.0*s*t2*ds_dn;
   ex_term = c0*(cn - c1)*poly*exp(hund);
   dc_term = c0*dc_dn*c_exp;
   dh_dn = om_term + ex_term + dc_term;

/*=======================================================================*/
/*  VII) calculate the function and its derivative wrt rho               */

  *fn = rho*h_per;
  *df_dn = h_per + rho*dh_dn;

/*=======================================================================*/
/* VIII) get derivative of f wrt density gradient                        */

   pre = 0.50/(sk*sk);
   do_term = (1.0 + 2.0*a_per*t2)/(oden*oden);
   do_term = do_term*(beta/omega);
   c_term = c0*(cn - c1)*(1.0 + hund)*exp(hund);
   p_fun = pre*(do_term + c_term);
   *df_dgn = p_fun*(drho/rho);

/*==============================================================================*/
}/*end routine*/
/*==============================================================================*/

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void pbe_gcc(double rho,double g_rho2,double *fn,double *df_dn,double *df_dgn,
	     double beta,double gamma)

/*==========================================================================*/
/*         Begin Routine                                                    */
/*=======================================================================*/
/* Perdew-Burke-Ernzerhof generalized gradient approximation to the      */
/* density functional for exchange-correlation energy of a many-electron */
/* system. See PRL 77 (1996) 3865 and White & Bird PRB 50 (1994) 4954.   */
/* Change of parameters for extended PBE only   :                        */
/* revPBE see Zhang and Wang, PRL 80 (1998) 890.                         */
/* rPBE see Hammer, Hansen and Norskov PRB 59 (1999) 7413.               */
/* xPBE   see Xu and Goddard, JCP 121 (2004) 4068.(beta and gamma differ)*/
{/*Begin Routine*/
/*=======================================================================*/
/*local variables */

 double ec,dec_dn,rr;
 double rs,drs;
 double srs,rstt,rs2,rs4;
 double gmag;
 double logterm;
 double preterm;
 double kf,ks,t,t2,t4;
 double expf,Afun,A2;
 double num,denom;
 double Hfun,Hlogterm;
 double rat;
 double Bfun;

 double dt_dn;
 double dA_dn;
 double dP_dn;
 double dB_dn;
 double dH_dn;

 double dt_dgn;
 double dP_dgn;
 double dB_dgn;
 double dH_dgn;

/* static variables for lda functional         */
 static double a = 0.03109070;
 static double aa = 3.0*M_PI*M_PI;
 static double alpha1 = 0.213700;
 static double beta1  = 7.59570;
 static double beta2  = 3.58760;
 static double beta3  = 1.63820;
 static double beta4  = 0.492940;

/* static variables for gc part */

 /*static double gamma = 0.031091;*/
 /* static double gamma = 0.03109069086965489503494086371273;*/
 /*static double beta = 0.066725;*/
 /* static double beta = 0.06672455060314922;*/
 static double pow1 = 1.0/3.0;
 static double bet_by_gam;

/*===========================================================================*/
/* First do e_c and its derivative  */

  rs = pow((3.0/(4.0*M_PI*rho)),pow1);
  rs2 = rs*rs;
  rs4 = rs2*rs2;
  rstt = pow(rs,1.5);
  srs = sqrt(rs);
  rr = 2.0*a*(beta1*srs + beta2*rs + beta3*rstt + beta4*rs2);
  logterm = log(1.0 + 1.0/rr);
  preterm = 1.0 + alpha1*rs;
  ec = -2.0*a*preterm*logterm;


  drs = 2.0*a*(0.5*beta1/srs + beta2 + 1.5*beta3*srs + 2.0*beta4*rs);
  dec_dn = (8.0*M_PI*a*rs4/9.0)*(alpha1*logterm - preterm*drs/(rr*(1.0 + rr)));

/*==============================================================================*/
/* Calculate the H function  */

  bet_by_gam = beta/gamma;
  gmag = sqrt(g_rho2);
  kf = pow(aa*rho,pow1);
  ks = sqrt(4.0*kf/M_PI);
  t = gmag/(2.0*ks*rho);

/* Afun */

  expf = exp(-ec/gamma) - 1.0;
  Afun = bet_by_gam/expf;
  
/* rest of H */

 t2 = t*t;
 t4 = t2*t2;
 A2 = Afun*Afun;
 num = 1.0 + Afun*t2;
 denom = 1.0 + Afun*t2 + A2*t4;
 rat = num/denom;
 Bfun = t2*bet_by_gam*rat;
 Hlogterm = log(1.0 + Bfun);
 Hfun = gamma*Hlogterm;

/*==============================================================================*/
/* Calculate the functional  */
/* Note: the lda correlation is already calculated in excpot_pw_lda, so it is */
/* omitted here. */

 *fn = rho*Hfun;

/*==============================================================================*/
/* Calculate the derivative of Hfun wrt density  */

 dt_dn = -7.0*t/(6.0*rho);
 dA_dn = bet_by_gam*(expf + 1.0)*dec_dn/(expf*expf*gamma);
 dP_dn = (t2*dA_dn + 2.0*Afun*t*dt_dn)/denom 
       - (num/(denom*denom))*(t2*(1.0 + 2.0*Afun*t2)*dA_dn 
                            + 2.0*Afun*t*(1.0 + 2.0*Afun*t2)*dt_dn);

 dB_dn = bet_by_gam*(2.0*t*dt_dn*rat + t2*dP_dn);
 dH_dn = gamma*dB_dn/(1.0 + Bfun);

/*==============================================================================*/
/* Calculate the derivative of the functional wrt density  */

 *df_dn = Hfun + rho*dH_dn; 
 
/*==============================================================================*/
/*  Calculate the derivative of the Hfun wrt |grad n|  */

  dt_dgn = 1.0/(2.0*ks*rho);
  dP_dgn = 2.0*Afun*t*dt_dgn/denom
         - num*2.0*Afun*t*(1.0 + 2.0*Afun*t2)*dt_dgn/(denom*denom);
  dB_dgn = bet_by_gam*(2.0*t*dt_dgn*rat + t2*dP_dgn);
  dH_dgn = gamma*dB_dgn/(1.0 + Bfun);

/*==============================================================================*/
/* Calculate the derivative of the functional wrt density  */

 *df_dgn = rho*dH_dgn;

/*==============================================================================*/
}/*end routine*/
/*==============================================================================*/



/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void lyp_gcc(double rho,double g_rho2,double *fn,double *df_dn,double *df_dgn)

/*==========================================================================*/
/*         Begin Routine                                                    */
{/*Begin Routine*/
/*=======================================================================*/
/*  local variables   */
   double rho2,drho;
   double denom,denom2,cexp,omega,domega,delta,ddelta;
   double arho,darho;
   double dprho;
   double rho13,rho113,rho43;
   double term1,brack1,brack2,brack3,term2;
   double dterm1,dterm2;

/* static variables */
   static double a = 0.049180;
   static double b = 0.1320;
   static double c = 0.25330;
   static double d = 0.3490;
   static double pi = M_PI;
   static double pow13;
   static double pow113;
   static double pow43;
   static double pow53;
   static double pow83;
   static double cf;

/*=======================================================================*/
/* I) Get square root density coefficient                                */

   pow13 = -1.0/3.0;
   pow113 = 11.0*pow13;
   pow43 = 4.0*pow13;
   pow53 = -5.0*pow13;
   pow83 = -8.0*pow13;
   cf = 0.30*pow((3.0*pi*pi),(2.0/3.0));
   drho = sqrt(g_rho2);

/*=======================================================================*/
/*  II) evaluate some grouped quantities                                 */

   rho13 = pow(rho,pow13);
   rho113 = pow(rho,pow113);
   rho43 = rho13/rho;
   denom = 1.0 + d*rho13;
   denom2 = denom*denom;
   cexp = exp(-c*rho13);

   omega = (cexp*rho113)/denom;
   brack1 = (d*rho43 + denom*(c*rho43 - 11.0/rho))/3.0;
   domega = omega*brack1/denom;

   delta = c*rho13 + d*rho13/denom;
   brack1 = (d*rho13 - denom)*d*rho43/denom2;
   ddelta = (-c*rho43 + brack1)/3.0;

/*=======================================================================*/
/* III) get all the brackets                                             */

   rho2 = rho*rho;
   term1 = -a*rho/denom;
   brack1 = ((5.0 - 7.0/6.0*delta)/3.0)*g_rho2;
   brack2 = 4.0*cf*pow(rho,pow83) + brack1;
   brack3 = (0.250*brack2 - 11.0*g_rho2/24.0)*rho2;
   term2 = -a*b*omega*brack3;

   arho = 0.250*brack2 - 11.0*g_rho2/24.0;
   brack1 = -7.0/18.0*ddelta*g_rho2;
   darho = 0.250*(32.0*cf*(pow(rho,pow53))/3.0 + brack1);

   dprho = 2.0*arho*rho + darho*rho2;
   dterm1 = -a*(1.0 - pow43*d*rho13)/denom2;
   dterm2 = -a*b*(domega*brack3 + omega*dprho);

/*=======================================================================*/
/* IV) evaluate the function                                             */

   *fn = term1 + term2;
/*=======================================================================*/
/*  V) evaluate the derivatives of the funcion                           */

    *df_dn = dterm1 + dterm2;
    brack1 = 0.50*((5.0 - 7.0/6.0*delta)/3.0) - 11.0/12.0;
    *df_dgn = -a*b*omega*brack1*rho2*drho;

/*==============================================================================*/
}/*end routine*/
/*==============================================================================*/


/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

 void lypm1(double rho,double g_rho2,double lap_rho,
            double *fn,double *df_dn,double *df_dgn,double *df_dln)

/*==========================================================================*/
/*         Begin Routine                                                    */
 {/*Begin Routine*/
/*=======================================================================*/
/*local variables */

 double fn_ret;
 double drho,ddrho;
 double denom,f1,f2,df1_dn,xnum,df2_dn;
 double f2fact,ddf2_dn,tw,p;
 double dpdn,ex;
 double epsi,deps_dn,deps_dgn;
 double rhop,rhop2,rhop8;

 static double tpi2,tt;
 static double a,b,c,d,k;
 static double  npow,npow2,npow4,npow5,npow6,npow7,npowp5,npowp2,npowp8;
 static double cf;

/*=======================================================================*/
/* Assign static variables                                              */

 tpi2 = 3.0*M_PI*M_PI;
 tt = 2.0/3.0;
 cf = 0.3*pow(tpi2,tt);
 npow  = -1.0/3.0;
 npow2 = 2.0*npow;
 npow4 = 4.0*npow;
 npow5 = 5.0*npow;
 npow7 = 7.0*npow;
 npowp5 = -npow5;
 npowp2 = -npow2;
 npowp8 = -8.0*npow;

#ifdef LYP_LAP
 a=0.04918;
 k=0.0;
#else
 a=0.052;
 k=0.002;
#endif
 b=0.132;
 c=0.2533;
 d=0.349;

/*=======================================================================*/
/* Local gradient and laplacian of density                              */

 drho = sqrt(g_rho2);
 ddrho = lap_rho;
 rhop = pow(rho,npow);
 rhop2 = pow(rho,npow2);
 rhop8 = pow(rho,npowp8);

/*=======================================================================*/
/* f1,f2 and their derivatives (two for f2) */

/* Ivanov pieces */
 epsi = exp(-k*g_rho2/rhop8);
 deps_dn = npowp8*k*g_rho2*epsi/(rho*rhop8);
 deps_dgn = -2.0*k*drho*epsi/rhop8;

 denom = 1.0 + d*rhop;
 f1 = rho/denom;
 f2 = b*rhop2/denom;

 xnum = 1.0 + 4.0*d*rhop/3.0;
 df1_dn = xnum/(denom*denom);

 xnum = 1.0 + denom;
 df2_dn = -(b*pow(rho,npow5)/3.0)*xnum/(denom*denom);

 f2fact = b/(9.0*rhop8);
 xnum = 10.0 + 12.0*d*rhop + 4.0*d*d*rhop2;
 ddf2_dn = f2fact*xnum/(denom*denom*denom);

/*=======================================================================*/
/* p and its derivative */

 tw = (drho*drho/rho - ddrho)/8.0;
 p = cf*pow(rho,npowp5) - 2.0*tw + (tw + 0.5*ddrho)/9.0;

 dpdn = 5.0*cf*pow(rho,npowp2)/3.0 + 17.0*(drho/rho)*(drho/rho)/72.0;

/*=======================================================================*/
/* calculate the function */

  ex = exp(-c*rhop);
  fn_ret = -a*(f1 + f2*p*ex);
  *fn = fn_ret*epsi;

/*=======================================================================*/
/* derivative of function wrt density, grad density and del^2 density */

 *df_dn = -a*(df1_dn + (df2_dn*p + f2*dpdn + c*pow(rho,npow4)*f2*p/3.0)*ex)*epsi
        + fn_ret*deps_dn;
 *df_dgn = 17.0*a*f2*ex*(drho/rho)*epsi/36.0 + fn_ret*deps_dgn;
 *df_dln = -21.0*a*f2*ex*epsi/72.0;

/*=======================================================================*/
}/*end routine*/
/*=======================================================================*/

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void lyp_lsda(double rho_up,double rho_dn,
              double drhox_up,double drhoy_up,double drhoz_up,
              double drhox_dn,double drhoy_dn,double drhoz_dn,
              double *fn,double *df_dn_up,
              double *df_dn_dn,double *dfx_up,double *dfy_up,
              double *dfz_up,double *dfx_dn,double *dfy_dn,double *dfz_dn)


/*==========================================================================*/
/*         Begin Routine                                                    */
{/*Begin Routine*/
/*=======================================================================*/
/*  local variables   */

  int iii;
  double g_rho2_up;
  double g_rho2_dn;
  double dot_up_dn;
  double drho_up,drho_dn;

  double denom,denom2,cexp,omega,delta;
  double domega,ddelta;
  double rho,rho2,rho13,rho113,rho43;
  double drho,drho2;
  double term1,term2;
  double dterm1_up,dterm2_up;
  double dterm1_dn,dterm2_dn;
  double ibrack1,ibrack2,ibrack3,ibrack4;
  double dibrack1_up,dibrack2_up,dibrack3_up,dibrack4_up;
  double dibrack1_dn,dibrack2_dn,dibrack3_dn,dibrack4_dn;
  double del_fact;
  double dgi2_up,dgi3_up,dgi4_up;
  double dgi2_dn,dgi3_dn,dgi4_dn;
  double dgi2x_up,dgi2y_up,dgi2z_up;
  double dgi3x_up,dgi3y_up,dgi3z_up;
  double dgi4x_up,dgi4y_up,dgi4z_up;
  double dgi2x_dn,dgi2y_dn,dgi2z_dn;
  double dgi3x_dn,dgi3y_dn,dgi3z_dn;
  double dgi4x_dn,dgi4y_dn,dgi4z_dn;
  double mbrack1,mbrack2,mbrack3;
  double dmbrack1_up,dmbrack2_up,dmbrack3_up;
  double dmbrack1_dn,dmbrack2_dn,dmbrack3_dn;
  double dgm1_up,dgm2_up,dgm3_up;
  double dgm1_dn,dgm2_dn,dgm3_dn;
  double dgm1x_up,dgm1y_up,dgm1z_up;
  double dgm2x_up,dgm2y_up,dgm2z_up;
  double dgm3x_up,dgm3y_up,dgm3z_up;
  double dgm1x_dn,dgm1y_dn,dgm1z_dn;
  double dgm2x_dn,dgm2y_dn,dgm2z_dn;
  double dgm3x_dn,dgm3y_dn,dgm3z_dn;
  double ux_up,uy_up,uz_up;
  double ux_dn,uy_dn,uz_dn;
  double da_up,da_dn,db_up,db_dn;
  double dp_up,dp_dn;
  double obrack,dobrack_up,dobrack_dn;
  double drho_upi,drho_dni;
  static double pi,cf;

/* Static variables */
   static double a = 0.049180;
   static double b = 0.1320;
   static double c = 0.25330;
   static double d = 0.3490;
   static double pow13;
   static double pow43;
   static double pow53;
   static double pow113;
   static double pow83;
   static double t113;
   static double eps = 1.0e-8;

/*=======================================================================*/
/* I) Assign useful constants                                            */

   pow13 = -1.0/3.0;
   pow43 = -4.0*pow13;
   pow53 = -5.0*pow13;
   pow113 = 11.0*pow13;
   pow83 = -8.0*pow13;
   t113  = pow(2.0,(-pow113));
   pi = M_PI;
   cf = 0.30*pow((3.0*pi*pi),(2.0/3.0));

/*=======================================================================*/
/* II) Calculate useful parts of functional                              */

   g_rho2_up = drhox_up*drhox_up + drhoy_up*drhoy_up
             + drhoz_up*drhoz_up;
   g_rho2_dn = drhox_dn*drhox_dn + drhoy_dn*drhoy_dn
             + drhoz_dn*drhoz_dn;
   dot_up_dn = drhox_up*drhox_dn 
             + drhoy_up*drhoy_dn
             + drhoz_up*drhoz_dn;
   drho_up = sqrt(g_rho2_up);
   drho_dn = sqrt(g_rho2_dn);

   drho_up = MAX(eps,drho_up);
   drho_dn = MAX(eps,drho_dn);

   drho_upi = 1.0/drho_up;
   drho_dni = 1.0/drho_dn;

   ux_up = drhox_up*drho_upi;
   uy_up = drhoy_up*drho_upi;
   uz_up = drhoz_up*drho_upi;
   ux_dn = drhox_dn*drho_dni;
   uy_dn = drhoy_dn*drho_dni;
   uz_dn = drhoz_dn*drho_dni;

/*=======================================================================*/
/* III) evaluate some grouped quantities                                 */

   rho = rho_up + rho_dn;
   rho2 = rho*rho;
   drho2 = g_rho2_up + g_rho2_dn + 2.0*dot_up_dn;
   drho = sqrt(drho2);
   rho13 = pow(rho,pow13);
   rho113 = pow(rho,pow113);
   rho43 = pow(rho,(-pow43));
   denom = 1.0 + d*rho13;
   denom2 = denom*denom;
   cexp = exp(-c*rho13);

   omega = (cexp*rho113)/denom;
   ibrack1 = (d*rho43 + denom*(c*rho43 - 11.0/rho))/3.0;
   domega = omega*ibrack1/denom;

   delta = c*rho13 + d*rho13/denom;
   ibrack1 = (d*rho13 - denom)*d*rho43/denom2;
   ddelta = (-c*rho43 + ibrack1)/3.0;
      
/*=======================================================================*/
/* IV) get all the brackets                                              */

   term1 = -4.0*a*rho_up*rho_dn/(denom*rho);
   dterm1_up = pow43*d*rho13*rho_up*rho_dn/(rho2*denom2)
             + 4.0*rho_dn/(rho*denom)*(1.0 - rho_up/rho);

   dterm1_dn = pow43*d*rho13*rho_up*rho_dn/(rho2*denom2)
             + 4.0*rho_up/(rho*denom)*(1.0 - rho_dn/rho);
   dterm1_up = -a*dterm1_up;
   dterm1_dn = -a*dterm1_dn;

   ibrack1 = t113*cf*(pow(rho_up,pow83) + pow(rho_dn,pow83));
   dibrack1_up = 2.0*pow43*t113*cf*pow(rho_up,pow53);
   dibrack1_dn = 2.0*pow43*t113*cf*pow(rho_dn,pow53);

   del_fact = (47.0 - 7.0*delta)/18.0;
   ibrack2 = del_fact*drho*drho;
   dibrack2_up = -7.0*ddelta*drho*drho/18.0;
   dibrack2_dn = dibrack2_up;
   dgi2_up = del_fact*2.0*drho_up;
   dgi2x_up = dgi2_up*ux_up + del_fact*2.0*drhox_dn;
   dgi2y_up = dgi2_up*uy_up + del_fact*2.0*drhoy_dn;
   dgi2z_up = dgi2_up*uz_up + del_fact*2.0*drhoz_dn;
   dgi2_dn = del_fact*2.0*drho_dn;
   dgi2x_dn = dgi2_dn*ux_dn + del_fact*2.0*drhox_up;
   dgi2y_dn = dgi2_dn*uy_dn + del_fact*2.0*drhoy_up;
   dgi2z_dn = dgi2_dn*uz_dn + del_fact*2.0*drhoz_up;

   del_fact = (5.0 - 1.0*delta/9.0)/2.0;
   ibrack3 = -del_fact*
             (g_rho2_up + g_rho2_dn);

   dibrack3_up = ddelta*(g_rho2_up + g_rho2_dn)/18.0;
   dibrack3_dn = dibrack3_up;
   dgi3_up = -del_fact*2.0*drho_up;
   dgi3x_up = dgi3_up*ux_up;
   dgi3y_up = dgi3_up*uy_up;
   dgi3z_up = dgi3_up*uz_up;
   dgi3_dn = -del_fact*2.0*drho_dn;
   dgi3x_dn = dgi3_dn*ux_dn;
   dgi3y_dn = dgi3_dn*uy_dn;
   dgi3z_dn = dgi3_dn*uz_dn;

   del_fact = (delta - 11.0)/9.0;
   ibrack4 = -del_fact*(rho_up*g_rho2_up/rho +
                        rho_dn*g_rho2_dn/rho);
   dibrack4_up = -(ddelta/(9.0*rho))
                 *(rho_up*g_rho2_up + rho_dn*g_rho2_dn)
                - (del_fact/rho)*((1.0 - rho_up/rho)*g_rho2_up
                - rho_dn/rho*g_rho2_dn);
   dibrack4_dn = -(ddelta/(9.0*rho))*
                  (rho_up*g_rho2_up + rho_dn*g_rho2_dn)
               -  (del_fact/rho)*(-rho_up/rho*g_rho2_up
               +  (1.0 - rho_dn/rho)*g_rho2_dn);
   dgi4_up = -del_fact*2.0*rho_up*drho_up/rho;
   dgi4x_up = dgi4_up*ux_up;
   dgi4y_up = dgi4_up*uy_up;
   dgi4z_up = dgi4_up*uz_up;
   dgi4_dn = -del_fact*2.0*rho_dn*drho_dn/rho;
   dgi4x_dn = dgi4_dn*ux_dn;
   dgi4y_dn = dgi4_dn*uy_dn;
   dgi4z_dn = dgi4_dn*uz_dn;

   mbrack1 = -2.0*rho*rho*drho*drho/3.0;
   dmbrack1_up = 2.0*mbrack1/rho;
   dmbrack1_dn = dmbrack1_up;
   dgm1x_up = -pow43*rho2*(drho_up*ux_up + drhox_dn);
   dgm1y_up = -pow43*rho2*(drho_up*uy_up + drhoy_dn);
   dgm1z_up = -pow43*rho2*(drho_up*uz_up + drhoz_dn);
   dgm1x_dn = -pow43*rho2*(drho_dn*ux_dn + drhox_up);
   dgm1y_dn = -pow43*rho2*(drho_dn*uy_dn + drhoy_up);
   dgm1z_dn = -pow43*rho2*(drho_dn*uz_dn + drhoz_up);

   del_fact = 2.0*rho*rho/3.0 - rho_up*rho_up;
   mbrack2 = del_fact*g_rho2_dn;
   dmbrack2_up = (pow43*rho - 2.0*rho_up)*g_rho2_dn;
   dmbrack2_dn = pow43*rho*g_rho2_dn;
   dgm2x_up = 0.0;
   dgm2y_up = 0.0;
   dgm2z_up = 0.0;
   dgm2x_dn = 2.0*del_fact*drho_dn*ux_dn;
   dgm2y_dn = 2.0*del_fact*drho_dn*uy_dn;
   dgm2z_dn = 2.0*del_fact*drho_dn*uz_dn;

   del_fact = 2.0*rho*rho/3.0 - rho_dn*rho_dn;
   mbrack3 = (2.0*rho*rho/3.0 - rho_dn*rho_dn)*g_rho2_up;
   dmbrack3_up = pow43*rho*g_rho2_up;
   dmbrack3_dn = (pow43*rho - 2.0*rho_dn)*g_rho2_up;
   dgm3x_up = 2.0*del_fact*drho_up*ux_up;
   dgm3y_up = 2.0*del_fact*drho_up*uy_up;
   dgm3z_up = 2.0*del_fact*drho_up*uz_up;
   dgm3x_dn = 0.0;
   dgm3y_dn = 0.0;
   dgm3z_dn = 0.0;

   obrack = rho_up*rho_dn*(ibrack1+ibrack2+ibrack3+ibrack4) +
            mbrack1+mbrack2+mbrack3;
   da_up = dibrack1_up + dibrack2_up + dibrack3_up + dibrack4_up;
   da_dn = dibrack1_dn + dibrack2_dn + dibrack3_dn + dibrack4_dn;
   db_up = dmbrack1_up + dmbrack2_up + dmbrack3_up;
   db_dn = dmbrack1_dn + dmbrack2_dn + dmbrack3_dn;
   dp_up = rho_dn*(ibrack1+ibrack2+ibrack3+ibrack4) +
           rho_up*rho_dn*da_up + db_up;
   dp_dn = rho_up*(ibrack1+ibrack2+ibrack3+ibrack4) +
          rho_up*rho_dn*da_dn + db_dn;
   term2 = -a*b*omega*obrack;

   dterm2_up = -a*b*(domega*obrack + omega*dp_up);
   dterm2_dn = -a*b*(domega*obrack + omega*dp_dn);

/*=======================================================================*/
/*  V) evaluate the function                                             */

  (*fn) = term1 + term2;

/*=======================================================================*/
/*  VI) evaluate the derivatives with respect to density                 */

   (*df_dn_up) = dterm1_up + dterm2_up;
   (*df_dn_dn) = dterm1_dn + dterm2_dn;
        

/*=======================================================================*/
/* VII) evaluate the derivatives with respect to gradients               */

   (*dfx_up) = -a*b*omega*
               (rho_up*rho_dn*(dgi2x_up+dgi3x_up+dgi4x_up)
             +             dgm1x_up+dgm2x_up+dgm3x_up);
   (*dfy_up) = -a*b*omega*
                (rho_up*rho_dn*(dgi2y_up+dgi3y_up+dgi4y_up)
             +                 dgm1y_up+dgm2y_up+dgm3y_up);
   (*dfz_up) = -a*b*omega*
              (rho_up*rho_dn*(dgi2z_up+dgi3z_up+dgi4z_up)
            +                 dgm1z_up+dgm2z_up+dgm3z_up);

   (*dfx_dn) = -a*b*omega*
              (rho_up*rho_dn*(dgi2x_dn+dgi3x_dn+dgi4x_dn)
            +                 dgm1x_dn+dgm2x_dn+dgm3x_dn);
   (*dfy_dn) = -a*b*omega*
              (rho_up*rho_dn*(dgi2y_dn+dgi3y_dn+dgi4y_dn)
            +                 dgm1y_dn+dgm2y_dn+dgm3y_dn);
   (*dfz_dn) = -a*b*omega*
              (rho_up*rho_dn*(dgi2z_dn+dgi3z_dn+dgi4z_dn)
            +                 dgm1z_dn+dgm2z_dn+dgm3z_dn);

/*==============================================================================*/
}/*end routine*/
/*==============================================================================*/



/*==============================================================*/
/*==============================================================*/

void brx89_lda(double n, double g_rho2, double ln, double tau,double *f,
                double *df_dn,double *df_dgn,double *df_dln,double *df_dtau)
/*==============================================================*/
{/* begin function */

  double x,a;
  double D,Q,rnq2,gamma=0.8;
  double eee1,eee1x,eee,b,db_dn,dx_dn,dmu_dn,db_dgn,dx_dgn;
  double dx_dln,db_dln,dx_dtau,db_dtau;
  double arg,xm2,den,brack,xfact;
  double fact,bp;
  double nft,ntt;
  double dD_dn;
  double gn;
  static double ft = 5.0/3.0;
  static double tt = 2.0/3.0;
  static double ot = 1.0/3.0;
  static double pi = M_PI;
  static double pi_tt;
  int ierr,iii;
  
/*==============================================================*/
/* I) Calculate right side of trancedental equation             */

   n *= 0.5;
   g_rho2 *= 0.25;
   ln *= 0.5;
   tau *= 0.5;

   pi_tt = pow(M_PI,(2.0/3.0));
   nft = pow(n,ft);
   ntt = pow(n,tt);
   gn = sqrt(g_rho2);
   D = tau-0.25*gn*gn/n;
   Q = (ln-2.0*gamma*D)/6.0;
   a = tt*pow(pi,tt)*nft/Q;
   rnq2 = 1.0/(9.0*Q*Q);

/*==============================================================*/
/* II) Solve the equation for x                                 */

   x = newt_raph(a,&ierr); 
   if(ierr==0 || isnan(x)==1) x = bisect(a);

/*==============================================================*/
/* III) Calculate the functional                                 */

   eee1 = exp(-x);
   eee1x = eee1*(1.0+x);
   bp   = x*x*x*eee1/(8.0*M_PI*n);
   b = pow(fabs(bp),ot)*bp/fabs(bp);
   xfact = (1.0-eee1-0.5*x*eee1)/b;
   *f  = -n*xfact;

/*==============================================================*/
/* IV) Calculate df/dn                                          */

   arg = 2.0*x/3.0;
   xm2 = x-2.0;
   den = 1.0 - arg - x/xm2;
   eee = exp(arg);
   fact = xm2*eee/den;
   brack = (tt*pi_tt*ntt/Q)*(ft + gamma*gn*gn/(12.0*Q*n));
   dx_dn = brack*fact;
   db_dn = (x*x/(24.0*M_PI*b*b*n))*((3.0-x)*dx_dn - x/n)*eee1;
   *df_dn = 0.5*(xfact*(n*db_dn/b-1.0) 
          - 0.5*(n*dx_dn/b)*eee1x);
   
/*==============================================================*/
/* IV) Calculate df/dgn                                         */

   dx_dgn = -(gamma*pi_tt*nft*rnq2)*(gn/n)*fact;
   db_dgn = (x*x/(24.0*M_PI*b*b*n))*dx_dgn*eee1*(3.0-x);
   *df_dgn = 0.5*((n*xfact/b)*db_dgn 
          - 0.5*(n*dx_dgn/b)*eee1x);
   
/*==============================================================*/
/* V) Calculate df/dln                                         */

   dx_dln = -(pi_tt*nft*rnq2)*fact;
   db_dln = (x*x/(24.0*M_PI*b*b*n))*dx_dln*eee1*(3.0-x);
   *df_dln = 0.5*((n*xfact/b)*db_dln 
          - 0.5*(n*dx_dln/b)*eee1x);
   
/*==============================================================*/
/* V) Calculate df/dtau                                         */

   dx_dtau = 2.0*(gamma*pi_tt*nft*rnq2)*fact;
   db_dtau = (x*x/(24.0*M_PI*b*b*n))*dx_dtau*eee1*(3.0-x);
   *df_dtau = 0.5*((n*xfact/b)*db_dtau 
          - 0.5*(n*dx_dtau/b)*eee1x);

/*==============================================================*/
}/* end function */
/*==============================================================*/



/*==============================================================*/
/*==============================================================*/

void brx2K_lda(double n, double g_rho2, double ln, double tau,double *f,
                double *df_dn,double *df_dgn,double *df_dln,double *df_dtau)
/*==============================================================*/
{/* begin function */

  double f_now;
  double x,a;
  double D,Q,rnq2;
  double eee1,eee1x,eee,b,db_dn,dx_dn,dmu_dn,db_dgn,dx_dgn;
  double dt_dn,dt_dtau,dw_dn,dw_dtau,dg_dn,dg_dtau;
  double dx_dln,db_dln,dx_dtau,db_dtau;
  double arg,xm2,den,brack,xfact;
  double fact,bp;
  double nft,ntt;
  double dD_dn;
  double gn;
  double gfact,ogfact,tterm,wterm,w2;
  static double ft = 5.0/3.0;
  static double tt = 2.0/3.0;
  static double ot = 1.0/3.0;
  static double pi = M_PI;
  static double pi_tt;
  static double ck,spi2;
  int ierr,iii;

/*==============================================================*/
/* Adjustable parameters (Should be gamma = 1, at = 1)          */

  double gamma = 1.0;
  double at = 1.0;

/*==============================================================*/
/* I) Calculate right side of trancedental equation             */

   n *= 0.5;
   g_rho2 *= 0.25;
   ln *= 0.5;
   tau *= 0.5;

   pi_tt = pow(M_PI,(2.0/3.0));
   nft = pow(n,ft);
   ntt = pow(n,tt);
   gn = sqrt(g_rho2);
   D = tau-0.25*gn*gn/n;
   Q = (ln-2.0*gamma*D)/6.0;
   a = tt*pow(pi,tt)*nft/Q;
   rnq2 = 1.0/(9.0*Q*Q);

/*--------------------------------------------------------------*/
/* These are needed for BRx2K                                   */

   spi2 = 6.0*M_PI*M_PI;
   ck = pow(spi2,tt)/ft;

/*==============================================================*/
/* II) Solve the equation for x                                 */

   x = newt_raph(a,&ierr); 
   if(ierr==0 || isnan(x)==1) x = bisect(a);

/*==============================================================*/
/* III) Calculate the functional                                 */

   eee1 = exp(-x);
   eee1x = eee1*(1.0+x);
   bp   = x*x*x*eee1/(8.0*M_PI*n);
   b = pow(fabs(bp),ot)*bp/fabs(bp);
   xfact = (1.0-eee1-0.5*x*eee1)/b;

/*--------------------------------------------------------------*/
/* The following is a modification for BRx2K                    */

   tterm = ck*nft/tau;
   wterm = (tterm-1.0)/(tterm+1.0);
   w2 = wterm*wterm;
   gfact = at*wterm*((w2 - 2.0)*w2 + 1);
   ogfact = 1.0 + gfact;

   f_now = -n*xfact;
   *f = f_now*ogfact;

/*==============================================================*/
/* IV) Calculate df/dn                                          */

   arg = 2.0*x/3.0;
   xm2 = x-2.0;
   den = 1.0 - arg - x/xm2;
   eee = exp(arg);
   fact = xm2*eee/den;
   brack = (tt*pi_tt*ntt/Q)*(ft + gamma*gn*gn/(12.0*Q*n));
   dx_dn = brack*fact;
   db_dn = (x*x/(24.0*M_PI*b*b*n))*((3.0-x)*dx_dn - x/n)*eee1;

/*--------------------------------------------------------------*/
/* The following is a modification for BRx2K                    */
 
   dt_dn = ft*ck*ntt/tau;
   dw_dn = (1.0 - wterm)*dt_dn/(tterm + 1.0);
   dg_dn = at*((5.0*w2 - 6.0)*w2 + 1.0)*dw_dn;

   *df_dn = 0.5*(xfact*(n*db_dn/b-1.0)
          - 0.5*(n*dx_dn/b)*eee1x)*ogfact + 0.5*f_now*dg_dn;

/*==============================================================*/
/* IV) Calculate df/dgn                                         */

   dx_dgn = -(gamma*pi_tt*nft*rnq2)*(gn/n)*fact;
   db_dgn = (x*x/(24.0*M_PI*b*b*n))*dx_dgn*eee1*(3.0-x);
   *df_dgn = 0.5*((n*xfact/b)*db_dgn 
          - 0.5*(n*dx_dgn/b)*eee1x)*ogfact;
   
/*==============================================================*/
/* V) Calculate df/dln                                         */

   dx_dln = -(pi_tt*nft*rnq2)*fact;
   db_dln = (x*x/(24.0*M_PI*b*b*n))*dx_dln*eee1*(3.0-x);
   *df_dln = 0.5*((n*xfact/b)*db_dln 
          - 0.5*(n*dx_dln/b)*eee1x)*ogfact;
   
/*==============================================================*/
/* V) Calculate df/dtau                                         */

   dx_dtau = 2.0*(gamma*pi_tt*nft*rnq2)*fact;
   db_dtau = (x*x/(24.0*M_PI*b*b*n))*dx_dtau*eee1*(3.0-x);

/*--------------------------------------------------------------*/
/* The following is a modification for BRx2K                    */
 
   dt_dtau = -ck*nft/(tau*tau);
   dw_dtau = (1.0 - wterm)*dt_dtau/(tterm + 1.0);
   dg_dtau = at*((5.0*w2 - 6.0)*w2 + 1.0)*dw_dtau;

   *df_dtau = 0.5*((n*xfact/b)*db_dtau
          - 0.5*(n*dx_dtau/b)*eee1x)*ogfact + 0.5*f_now*dg_dtau;

/*==============================================================*/
}/* end function */
/*==============================================================*/



/*==============================================================*/
/*==============================================================*/

void brx89_lsda(double n, double g_rho2, double ln, double tau,double *f,
                double *df_dn,double *df_dgn,double *df_dln,double *df_dtau)
/*==============================================================*/
{/* begin function */

  double x,a;
  double D,Q,rnq2,gamma=0.8;
  double eee1,eee1x,eee,b,db_dn,dx_dn,dmu_dn,db_dgn,dx_dgn;
  double dx_dln,db_dln,dx_dtau,db_dtau;
  double arg,xm2,den,brack,xfact;
  double fact,bp;
  double nft,ntt;
  double dD_dn;
  double gn;
  static double ft = 5.0/3.0;
  static double tt = 2.0/3.0;
  static double ot = 1.0/3.0;
  static double pi = M_PI;
  static double pi_tt;
  int ierr,iii;
  
/*==============================================================*/
/* I) Calculate right side of trancedental equation             */

   pi_tt = pow(M_PI,(2.0/3.0));
   nft = pow(n,ft);
   ntt = pow(n,tt);
   gn = sqrt(g_rho2);
   D = tau-0.25*gn*gn/n;
   Q = (ln-2.0*gamma*D)/6.0;
   a = tt*pow(pi,tt)*nft/Q;
   rnq2 = 1.0/(9.0*Q*Q);

/*==============================================================*/
/* II) Solve the equation for x                                 */

   x = newt_raph(a,&ierr); 
   if(ierr==0 || isnan(x)==1) x = bisect(a);

/*==============================================================*/
/* III) Calculate the functional                                 */

   eee1 = exp(-x);
   eee1x = eee1*(1.0+x);
   bp   = x*x*x*eee1/(8.0*M_PI*n);
   b = pow(fabs(bp),ot)*bp/fabs(bp);
   xfact = (1.0-eee1-0.5*x*eee1)/b;
   *f  = -0.5*n*xfact;

/*==============================================================*/
/* IV) Calculate df/dn                                          */

   arg = 2.0*x/3.0;
   xm2 = x-2.0;
   den = 1.0 - arg - x/xm2;
   eee = exp(arg);
   fact = xm2*eee/den;
   brack = (tt*pi_tt*ntt/Q)*(ft + gamma*gn*gn/(12.0*Q*n));
   dx_dn = brack*fact;
   db_dn = (x*x/(24.0*M_PI*b*b*n))*((3.0-x)*dx_dn - x/n)*eee1;
   *df_dn = 0.5*(xfact*(n*db_dn/b-1.0) 
          - 0.5*(n*dx_dn/b)*eee1x);
   
/*==============================================================*/
/* IV) Calculate df/dgn                                         */

   dx_dgn = -(gamma*pi_tt*nft*rnq2)*(gn/n)*fact;
   db_dgn = (x*x/(24.0*M_PI*b*b*n))*dx_dgn*eee1*(3.0-x);
   *df_dgn = 0.5*((n*xfact/b)*db_dgn 
          - 0.5*(n*dx_dgn/b)*eee1x);
   
/*==============================================================*/
/* V) Calculate df/dln                                         */

   dx_dln = -(pi_tt*nft*rnq2)*fact;
   db_dln = (x*x/(24.0*M_PI*b*b*n))*dx_dln*eee1*(3.0-x);
   *df_dln = 0.5*((n*xfact/b)*db_dln 
          - 0.5*(n*dx_dln/b)*eee1x);
   
/*==============================================================*/
/* V) Calculate df/dtau                                         */

   dx_dtau = 2.0*(gamma*pi_tt*nft*rnq2)*fact;
   db_dtau = (x*x/(24.0*M_PI*b*b*n))*dx_dtau*eee1*(3.0-x);
   *df_dtau = 0.5*((n*xfact/b)*db_dtau 
          - 0.5*(n*dx_dtau/b)*eee1x);

/*==============================================================*/
}/* end function */
/*==============================================================*/



/*==============================================================*/
/*==============================================================*/

void brx2K_lsda(double n, double g_rho2, double ln, double tau,double *f,
                double *df_dn,double *df_dgn,double *df_dln,double *df_dtau)
/*==============================================================*/
{/* begin function */

  double f_now;
  double x,a;
  double D,Q,rnq2;
  double eee1,eee1x,eee,b,db_dn,dx_dn,dmu_dn,db_dgn,dx_dgn;
  double dt_dn,dt_dtau,dw_dn,dw_dtau,dg_dn,dg_dtau;
  double dx_dln,db_dln,dx_dtau,db_dtau;
  double arg,xm2,den,brack,xfact;
  double fact,bp;
  double nft,ntt;
  double dD_dn;
  double gn;
  double gfact,ogfact,tterm,wterm,w2;
  static double ft = 5.0/3.0;
  static double tt = 2.0/3.0;
  static double ot = 1.0/3.0;
  static double pi = M_PI;
  static double pi_tt;
  static double ck,spi2;
  int ierr,iii;

/*==============================================================*/
/* Adjustable parameters (Should be gamma = 1, at = 1)          */

  double gamma = 0.8;
  double at = 0.928;

/*==============================================================*/
/* I) Calculate right side of trancedental equation             */

   pi_tt = pow(M_PI,(2.0/3.0));
   nft = pow(n,ft);
   ntt = pow(n,tt);
   gn = sqrt(g_rho2);
   D = tau-0.25*gn*gn/n;
   Q = (ln-2.0*gamma*D)/6.0;
   a = tt*pow(pi,tt)*nft/Q;
   rnq2 = 1.0/(9.0*Q*Q);

/*--------------------------------------------------------------*/
/* These are needed for BRx2K                                   */

   spi2 = 6.0*M_PI*M_PI;
   ck = pow(spi2,tt)/ft;

/*==============================================================*/
/* II) Solve the equation for x                                 */

   x = newt_raph(a,&ierr); 
   if(ierr==0 || isnan(x)==1) x = bisect(a);

/*==============================================================*/
/* III) Calculate the functional                                 */

   eee1 = exp(-x);
   eee1x = eee1*(1.0+x);
   bp   = x*x*x*eee1/(8.0*M_PI*n);
   b = pow(fabs(bp),ot)*bp/fabs(bp);
   xfact = (1.0-eee1-0.5*x*eee1)/b;

/*--------------------------------------------------------------*/
/* The following is a modification for BRx2K                    */

   tterm = ck*nft/tau;
   wterm = (tterm-1.0)/(tterm+1.0);
   w2 = wterm*wterm;
   gfact = at*wterm*((w2 - 2.0)*w2 + 1);
   ogfact = 1.0 + gfact;

   f_now = -0.5*n*xfact;
   *f = f_now*ogfact;

/*==============================================================*/
/* IV) Calculate df/dn                                          */

   arg = 2.0*x/3.0;
   xm2 = x-2.0;
   den = 1.0 - arg - x/xm2;
   eee = exp(arg);
   fact = xm2*eee/den;
   brack = (tt*pi_tt*ntt/Q)*(ft + gamma*gn*gn/(12.0*Q*n));
   dx_dn = brack*fact;
   db_dn = (x*x/(24.0*M_PI*b*b*n))*((3.0-x)*dx_dn - x/n)*eee1;

/*--------------------------------------------------------------*/
/* The following is a modification for BRx2K                    */
 
   dt_dn = ft*ck*ntt/tau;
   dw_dn = (1.0 - wterm)*dt_dn/(tterm + 1.0);
   dg_dn = at*((5.0*w2 - 6.0)*w2 + 1.0)*dw_dn;

   *df_dn = 0.5*(xfact*(n*db_dn/b-1.0)
          - 0.5*(n*dx_dn/b)*eee1x)*ogfact + f_now*dg_dn;

/*==============================================================*/
/* IV) Calculate df/dgn                                         */

   dx_dgn = -(gamma*pi_tt*nft*rnq2)*(gn/n)*fact;
   db_dgn = (x*x/(24.0*M_PI*b*b*n))*dx_dgn*eee1*(3.0-x);
   *df_dgn = 0.5*((n*xfact/b)*db_dgn 
          - 0.5*(n*dx_dgn/b)*eee1x)*ogfact;
   
/*==============================================================*/
/* V) Calculate df/dln                                         */

   dx_dln = -(pi_tt*nft*rnq2)*fact;
   db_dln = (x*x/(24.0*M_PI*b*b*n))*dx_dln*eee1*(3.0-x);
   *df_dln = 0.5*((n*xfact/b)*db_dln 
          - 0.5*(n*dx_dln/b)*eee1x)*ogfact;
   
/*==============================================================*/
/* V) Calculate df/dtau                                         */

   dx_dtau = 2.0*(gamma*pi_tt*nft*rnq2)*fact;
   db_dtau = (x*x/(24.0*M_PI*b*b*n))*dx_dtau*eee1*(3.0-x);

/*--------------------------------------------------------------*/
/* The following is a modification for BRx2K                    */
 
   dt_dtau = -ck*nft/(tau*tau);
   dw_dtau = (1.0 - wterm)*dt_dtau/(tterm + 1.0);
   dg_dtau = at*((5.0*w2 - 6.0)*w2 + 1.0)*dw_dtau;

   *df_dtau = 0.5*((n*xfact/b)*db_dtau
          - 0.5*(n*dx_dtau/b)*eee1x)*ogfact + f_now*dg_dtau;

/*==============================================================*/
}/* end function */
/*==============================================================*/





/*===============================================================*/

double newt_raph(double a, int *ierr)

/*===============================================================*/
{/* begin function */

   int count=0,iii;
   double x = -0.375*(log(a)-5.0);
   double f=100.0,arg,eee,xm2,fp;
   static double tol = 1.0e-12;
   double eps=1.0e-12;
   static int max_iter=50;

   *ierr=1;
   if(a > 0.0 && x < 2.0) x = 2.0+eps;
   if(a < 0.0) x = 2.0-eps;
   while(fabs(f) > tol) {
     xm2 = x-2.0;
     arg = 2.0*x/3.0;
     eee = exp(-arg);
     f = x*eee/xm2 - a;
     fp = (eee/xm2)*(1.0 - arg - x/xm2);
     x -= f/fp;
     x = fabs(x);
     count++;
     if(count>max_iter) {*ierr=0; break;}
   } /* endwhile */
   
   return x;

/*==============================================================*/
}/* end function */
/*==============================================================*/


/*===============================================================*/

double bisect(double a)

/*===============================================================*/
{/* begin function */

   int count=0;
   double f=100.0,arg,eee,xm2;
   static double tol = 1.0e-10;
   double x;
   double x1,x2;
   double f1,f2;
   double eps=1.0e-12;
   int iii;
   static int max_iter=500;

   if(a > 0.0) {x1 = 2.0+eps; x2 = 10.0;}
   if(a < 0.0) {x2 = 2.0-eps; x1 = 0.0;}
   while(fabs(f) > tol){
     x = 0.5*(x1+x2);
     xm2 = x-2.0;
     arg = 2.0*x/3.0;
     eee = exp(-arg);
     f = x*eee/xm2 - a;
     if(f>0.0) x1=x;
     if(f<0.0) x2=x;
     count++;
     if(count>max_iter) {
      printf("*********************WARNING********************************\n");
      printf("Convergence not fully reached\n");
      printf("in BRx2K functional: bisect, f = %.12g\n",f);
      if(fabs(f) < 10e-6){
        printf("This is probably not catastrophic, but a quick\n");
        printf("check of forces using the debug_cp option would be wise\n");
      } else {
        printf("This could cause errors or drifts or indicate a problem\n");
        printf("with your system setup.  A thorough check would be wise\n");
      }
      printf("*********************WARNING********************************\n");
      break;
     }/* endif */
   }/* end while */
   
   return x;

/*==============================================================*/
}/* end function */
/*==============================================================*/



/*==============================================================*/
/*==============================================================*/

void tau1_lda(double n, double g_rho2, double ln, double tau,double *f,
                double *df_dn,double *df_dgn,double *df_dln,double *df_dtau,
                int nstate)
/*==============================================================*/
{/* begin function */

 static double b1,b2,b3,b4,b5,b6;
 static double c1,c2,c3,c4,c5,c6;
 static double ae,ae2,beta2,c_exp;
 static double gamma2,cp;
 static double tt,ot,ft;
 double fpd,q,q2,k,k2,k3,k4,bsq,not,nft;
 double dk_dn,dk_dtau,dq_dk,dq2_dk;
 double dbsq_dn,dbsq_dln,dbsq_dgn;
 double dfpd_dn,dfpd_dgn,dfpd_dln,dfpd_dtau;
 double r1,r2,r3,r4,r5,r6;
 double dr1,dr2,dr3,dr4,dr5,dr6;
 double s1,s2,s3,s4,s5,s6;
 double ds1,ds2,ds3,ds4,ds5,ds6;
 double eee,n2,n3;
 double fpp,dns,bpp2,pre;
 double dfpp_dn,dfpp_dgn,dfpp_dln,dfpp_dtau;
 double dbpp2_dn,dbpp2_dgn,dbpp2_dln,dbpp2_dtau;
 double bfact;
 double gn;

/*==============================================================*/
/* Assign constants                                             */

 b1 = 2.763169;
 b2 = 1.757515;
 b3 = 1.470000;
 b4 = 0.568985;
 b5 = 1.572202;
 b6 = 1.885389;

 c1 = 0.583090;
 c2 = 1.757515;
 c3 = 4.377624;
 c4 = 0.568985;
 c5 = 0.331770;
 c6 = 2.302031;

 c_exp = 0.253;
 ae    = 1.299;
 ae2   = ae*ae;
 beta2 = 0.0;

 tt = 2.0/3.0;
 ot = 1.0/3.0;
 ft = 4.0/3.0;

 gamma2 = 0.175;
 cp     = 0.045;

/*==============================================================*/
/* Calculate anti-parallel part                                 */

 gn  = sqrt(g_rho2);
 not = pow(n,ot);
 nft = pow(n,ft);
 k   = tt*ae2*tau/n;
 k2  = k*k;
 k3  = k2*k;
 k4  = k2*k2;
 eee = exp(-c_exp/not);
 bsq = (beta2/24.0)*(n*ln - gn*gn)*eee;

 r1 = b1/(1.0+b2*k);
 r2 = b3/k;
 r3 = (b4+k)/k;
 r4 = b5/k;
 r5 = b6/k2;
 q = -r1 + r2*log(r3) + r4 - r5;
 
 s1 = c1/(k2*(1.0+c2*k));
 s2 = c3/k3;
 s3 = (c4+k)/k;
 s4 = c5/k3;
 s5 = c6/k4;
 q2 = -s1 - s2*log(s3) + s4 + s5;

 n2 = n*n;
 n3 = n*n*n;
 fpd = 0.25*n2*q + bsq*q2;

/*==============================================================*/
/* Calculate density derivative of anti-parallel part           */

 dk_dn = -ae2*tt*tau/(n*n);

 dr1 = r1*b2/(1.0+b2*k);
 dr2 = b3/k2;
 dr3 = (b3/k)*(1.0/(b4+k));
 dr4 = 1.0 - (b4+k)/k;
 dr5 = b5/k2;
 dr6 = 2.0*b6/k3;
 dq_dk = dr1 - dr2*log(r3) + dr3*dr4 - dr5 + dr6;

 ds1 = (s1*(2.0 + 3.0*c2*k)/k)*1.0/(1.0+c2*k);
 ds2 = 3.0*c3/k4;
 ds3 = (c3/k3)*1.0/(c4+k);
 ds4 = 1.0 - (c4+k)/k;
 ds5 = 3.0*c5/k4;
 ds6 = 4.0*c6/(k*k4);
 dq2_dk = ds1 + ds2*log(s3) - ds3*ds4 - ds5 - ds6;

 dbsq_dn = bsq*c_exp/(3.0*nft) + beta2*ln*eee/24.0;

 dfpd_dn = 0.5*n*q + 0.25*n2*dq_dk*dk_dn + dbsq_dn*q2 
         + bsq*dq2_dk*dk_dn;

 
/*==============================================================*/
/* Calculate df(anti-parallel)/d_gn                             */

  dbsq_dgn = -(beta2/12.0)*gn*eee;
  dfpd_dgn = q2*dbsq_dgn;

/*==============================================================*/
/* Calculate df(anti-parallel)/d_ln                             */

  dbsq_dln = beta2*n*eee/24.0;
  dfpd_dln = q2*dbsq_dln;

/*==============================================================*/
/* Calculate df(anti-parallel)/d_tau                            */

 dk_dtau = ae2*2.0/(3.0*n);
 dfpd_dtau = 0.25*n2*dq_dk*dk_dtau + bsq*dq2_dk*dk_dtau;

/*==============================================================*/
/* Calculate parallel part                                      */

   dns = (double)nstate;
   pre = 1.0 - 1.0/dns;

   bfact = 0.125*ln + tau - 0.25*gn*gn/n;
   bpp2 = (gamma2/n)*bfact*eee;
   fpp = 2.0*pre*0.125*n2*cp*(q + bpp2*q2);

/*==============================================================*/
/* Calculate df(parallel)/dn                                    */

   dbpp2_dn = (gamma2/n)*bfact*(c_exp/(3.0*nft) - 1.0/n)*eee
            + (gamma2/n3)*0.25*gn*gn*eee;
   dfpp_dn = 2.0*pre*0.25*cp*n*(q + bpp2*q2)
           + 2.0*pre*0.125*n2*cp*
                (dq_dk*dk_dn + bpp2*dq2_dk*dk_dn + q2*dbpp2_dn);

/*==============================================================*/
/* Calculate df(parallel)/dgn                                   */

   dbpp2_dgn = -0.5*gamma2*gn*eee/n2;
   dfpp_dgn = 2.0*pre*0.125*n2*cp*q2*dbpp2_dgn;

/*==============================================================*/
/* Calculate df(parallel)/dln                                   */

  dbpp2_dln = 0.125*gamma2*eee/n;
  dfpp_dln = 2.0*pre*0.125*n2*cp*q2*dbpp2_dln;

/*==============================================================*/
/* Calculate df(parallel)/dtau                                  */

  dbpp2_dtau = gamma2*eee/n;
  dfpp_dtau = 2.0*pre*0.125*n2*cp*
              (dq_dk*dk_dtau + bpp2*dq2_dk*dk_dtau + q2*dbpp2_dtau);

/*==============================================================*/
/* Total functional and its derivatives                         */

  *f       = fpd + fpp;
  *df_dn   = dfpd_dn + dfpp_dn;
  *df_dgn  = dfpd_dgn + dfpp_dgn;
  *df_dln  = dfpd_dln + dfpp_dln;
  *df_dtau = dfpd_dtau + dfpp_dtau;

/*==============================================================================*/
}/*end routine*/
/*==============================================================================*/



/*==============================================================*/
/*==============================================================*/

void tau1_lsda(double n_up, double n_dn,double g_rho2_up,double g_rho2_dn,
               double ln_up,double ln_dn, double tau_up,double tau_dn,
               double *f,double *df_dn_up,double *df_dn_dn,
               double *df_dgn_up,double *df_dgn_dn,
               double *df_dln_up,double *df_dln_dn,
               double *df_dtau_up,double *df_dtau_dn,int nstate_up,int nstate_dn)

/*==============================================================*/
{/* begin function */

 static double b1,b2,b3,b4,b5,b6;
 static double c1,c2,c3,c4,c5,c6;
 static double ae,ae2,beta2,c_exp;
 static double gamma2,cp;
 static double tt,ot,ft;
 double fpd,qpd,k_pd,k2_pd,k3_pd,k4_pd,n,not,nft;
 double fpp,qpp,q2pp,k2_up,k3_up,k4_up;
 double fdd,qdd,q2dd,k2_dn,k3_dn,k4_dn;
 double k_up,k_dn;

 double dk_pd_dn_up,dk_pd_dtau_up,dk_pd_dn_dn,dk_pd_dtau_dn;
 double dk_up_dn_up,dk_up_dtau_up,dk_dn_dn_dn,dk_dn_dtau_dn;
 double dq_dkpd,dq_dk_up,dq_dk_dn,dq2_dk_up,dq2_dk_dn;
 double dbsq_dn,dbsq_dln,dbsq_dgn;
 double dfpd_dn_up,dfpd_dn_dn,dfpd_dgn_up,dfpd_dgn_dn;
 double dfpd_dln_up,dfpd_dln_dn,dfpd_dtau_up,dfpd_dtau_dn;
 double r1,r2,r3,r4,r5,r6;
 double rr1,rr2,rr3,rr4,rr5,rr6;
 double dr1,dr2,dr3,dr4,dr5,dr6;
 double drr1,drr2,drr3,drr4,drr5,drr6;
 double s1,s2,s3,s4,s5,s6;
 double ss1,ss2,ss3,ss4,ss5,ss6;
 double ds1,ds2,ds3,ds4,ds5,ds6;
 double dss1,dss2,dss3,dss4,dss5,dss6;
 double eee,n2_up,n3_up,n2_dn,n3_dn;
 double dns_up,bpp2,pre_up;
 double dns_dn,bdd2,pre_dn;
 double dfpp_dn_up,dfpp_dn_dn,dfdd_dn_up,dfdd_dn_dn;
 double dfpp_dgn_up,dfpp_dln_up,dfpp_dtau_up;
 double dfdd_dgn_dn,dfdd_dln_dn,dfdd_dtau_dn;
 double dbpp2_dn_up,dbdd2_dn_dn,dbpp2_dn_dn,dbdd2_dn_up;
 double dbpp2_dgn_up,dbpp2_dln_up,dbpp2_dtau_up;
 double dbdd2_dgn_dn,dbdd2_dln_dn,dbdd2_dtau_dn;
 double bfact_up,bfact_dn;
 double gn_up,gn_dn;

/*==============================================================*/
/* Assign constants                                             */

 b1 = 2.763169;
 b2 = 1.757515;
 b3 = 1.470000;
 b4 = 0.568985;
 b5 = 1.572202;
 b6 = 1.885389;

 c1 = 0.583090;
 c2 = 1.757515;
 c3 = 4.377624;
 c4 = 0.568985;
 c5 = 0.331770;
 c6 = 2.302031;

 c_exp = 0.253;
 ae    = 1.299;
 ae2   = ae*ae;
 beta2 = 0.0;

 tt = 2.0/3.0;
 ot = 1.0/3.0;
 ft = 4.0/3.0;

 gamma2 = 0.175;
 cp     = 0.045;

/*==============================================================*/
/* Calculate anti-parallel part                                 */

 gn_up  = sqrt(g_rho2_up);
 gn_dn  = sqrt(g_rho2_dn);
 n = n_up + n_dn;
 not = pow(n,ot);
 nft = pow(n,ft);
 n2_up = n_up*n_up;
 n2_dn = n_dn*n_dn;
 n3_up = n2_up*n_up;
 n3_dn = n2_dn*n_dn;

 k_up   = tt*ae2*tau_up/n_up;
 k_dn   = tt*ae2*tau_dn/n_dn;
 k2_up  = k_up*k_up;
 k3_up  = k2_up*k_up;
 k4_up  = k2_up*k2_up;
 k2_dn  = k_dn*k_dn;
 k3_dn  = k2_dn*k_dn;
 k4_dn  = k2_dn*k2_dn;
 
 k_pd = 2.0*k_up*k_dn/(k_up+k_dn);
 k2_pd  = k_pd*k_pd;
 k3_pd  = k2_pd*k_pd;
 k4_pd  = k2_pd*k2_pd;
 eee = exp(-c_exp/not);

 r1 = b1/(1.0+b2*k_pd);
 r2 = b3/k_pd;
 r3 = (b4+k_pd)/k_pd;
 r4 = b5/k_pd;
 r5 = b6/k2_pd;
 qpd = -r1 + r2*log(r3) + r4 - r5;
 
 fpd = n_up*n_dn*qpd;


/*==============================================================*/
/* Calculate density derivatives of anti-parallel part          */

 dk_up_dn_up = -ae2*tt*tau_up/(n_up*n_up);
 dk_dn_dn_dn = -ae2*tt*tau_dn/(n_dn*n_dn);

 dk_pd_dn_up = (1.0/(k_up+k_dn))*(2.0*k_dn - k_pd);
 dk_pd_dn_dn = (1.0/(k_up+k_dn))*(2.0*k_up - k_pd);

 dr1 = r1*b2/(1.0+b2*k_pd);
 dr2 = b3/k2_pd;
 dr3 = (b3/k_pd)*(1.0/(b4+k_pd));
 dr4 = 1.0 - (b4+k_pd)/k_pd;
 dr5 = b5/k2_pd;
 dr6 = 2.0*b6/k3_pd;
 dq_dkpd = dr1 - dr2*log(r3) + dr3*dr4 - dr5 + dr6;

 dfpd_dn_up = n_dn*qpd + n_up*n_dn*dq_dkpd*dk_pd_dn_up*dk_up_dn_up;
 dfpd_dn_dn = n_up*qpd + n_up*n_dn*dq_dkpd*dk_pd_dn_dn*dk_dn_dn_dn;


/*==============================================================*/
/* Calculate df(anti-parallel)/d_tau                             */

  dk_up_dtau_up = ae2*2.0/(3.0*n_up);
  dk_dn_dtau_dn = ae2*2.0/(3.0*n_dn);
  dfpd_dtau_up = n_up*n_dn*dq_dkpd*dk_pd_dn_up*dk_up_dtau_up;
  dfpd_dtau_dn = n_up*n_dn*dq_dkpd*dk_pd_dn_dn*dk_dn_dtau_dn;

  *df_dtau_up = dfpd_dtau_up;
  *df_dtau_dn = dfpd_dtau_dn;

/*==============================================================*/
/* Calculate parallel part                                      */

   dns_up = (double)nstate_up;
   dns_dn = (double)nstate_dn;
   pre_up = 1.0 - 1.0/dns_up;
   pre_dn = 1.0 - 1.0/dns_dn;

   r1 = b1/(1.0+b2*k_up);
   r2 = b3/k_up;
   r3 = (b4+k_up)/k_up;
   r4 = b5/k_up;
   r5 = b6/k2_up;
   qpp = -r1 + r2*log(r3) + r4 - r5;
 
   rr1 = b1/(1.0+b2*k_dn);
   rr2 = b3/k_dn;
   rr3 = (b4+k_dn)/k_dn;
   rr4 = b5/k_dn;
   rr5 = b6/k2_dn;
   qdd = -rr1 + rr2*log(rr3) + rr4 - rr5;
 
   s1 = c1/(k2_up*(1.0+c2*k_up));
   s2 = c3/k3_up;
   s3 = (c4+k_up)/k_up;
   s4 = c5/k3_up;
   s5 = c6/k4_up;
   q2pp = -s1 - s2*log(s3) + s4 + s5;

   ss1 = c1/(k2_dn*(1.0+c2*k_dn));
   ss2 = c3/k3_dn;
   ss3 = (c4+k_dn)/k_dn;
   ss4 = c5/k3_dn;
   ss5 = c6/k4_dn;
   q2dd = -ss1 - ss2*log(ss3) + ss4 + ss5;

   bfact_up = 0.125*ln_up + tau_up - 0.25*gn_up*gn_up/n_up;
   bfact_dn = 0.125*ln_dn + tau_dn - 0.25*gn_dn*gn_dn/n_dn;

   bpp2 = (gamma2/n_up)*bfact_up*eee;
   bdd2 = (gamma2/n_dn)*bfact_dn*eee;

   fpp = pre_up*0.5*n2_up*cp*(qpp + bpp2*q2pp);
   fdd = pre_dn*0.5*n2_dn*cp*(qdd + bdd2*q2dd);


/*==============================================================*/
/* Calculate df(parallel)/dn                                    */

   dbpp2_dn_up = (gamma2/n_up)*bfact_up*(c_exp/(3.0*nft) - 1.0/n_up)*eee
               + 0.25*(gamma2/n3_up)*gn_up*gn_up*eee;

   dbdd2_dn_dn = (gamma2/n_dn)*bfact_dn*(c_exp/(3.0*nft) - 1.0/n_dn)*eee
               + 0.25*(gamma2/n3_dn)*gn_dn*gn_dn*eee;

   dbpp2_dn_dn = (gamma2/n_up)*bfact_up*c_exp/(3.0*nft)*eee;
   dbdd2_dn_up = (gamma2/n_dn)*bfact_dn*c_exp/(3.0*nft)*eee;

   dr1 = r1*b2/(1.0+b2*k_up);
   dr2 = b3/k2_up;
   dr3 = (b3/k_up)*(1.0/(b4+k_up));
   dr4 = 1.0 - (b4+k_up)/k_up;
   dr5 = b5/k2_up;
   dr6 = 2.0*b6/k3_up;
   dq_dk_up = dr1 - dr2*log(r3) + dr3*dr4 - dr5 + dr6;

   drr1 = rr1*b2/(1.0+b2*k_dn);
   drr2 = b3/k2_dn;
   drr3 = (b3/k_dn)*(1.0/(b4+k_dn));
   drr4 = 1.0 - (b4+k_dn)/k_dn;
   drr5 = b5/k2_dn;
   drr6 = 2.0*b6/k3_dn;
   dq_dk_dn = drr1 - drr2*log(rr3) + drr3*drr4 - drr5 + drr6;

   ds1 = (s1*(2.0 + 3.0*c2*k_up)/k_up)*1.0/(1.0+c2*k_up);
   ds2 = 3.0*c3/k4_up;
   ds3 = (c3/k3_up)*1.0/(c4+k_up);
   ds4 = 1.0 - (c4+k_up)/k_up;
   ds5 = 3.0*c5/k4_up;
   ds6 = 4.0*c6/(k_up*k4_up);
   dq2_dk_up = ds1 + ds2*log(s3) - ds3*ds4 - ds5 - ds6;

   dss1 = (ss1*(2.0 + 3.0*c2*k_dn)/k_dn)*1.0/(1.0+c2*k_dn);
   dss2 = 3.0*c3/k4_dn;
   dss3 = (c3/k3_dn)*1.0/(c4+k_dn);
   dss4 = 1.0 - (c4+k_dn)/k_dn;
   dss5 = 3.0*c5/k4_dn;
   dss6 = 4.0*c6/(k_dn*k4_dn);
   dq2_dk_dn = dss1 + dss2*log(ss3) - dss3*dss4 - dss5 - dss6;

   dfpp_dn_up = pre_up*cp*n_up*(qpp + bpp2*q2pp)
           + pre_up*0.5*n2_up*cp*
             (dq_dk_up*dk_up_dn_up + bpp2*dq2_dk_up*dk_up_dn_up + q2pp*dbpp2_dn_up);

   dfdd_dn_dn = pre_dn*cp*n_dn*(qdd + bdd2*q2dd)
           + pre_dn*0.5*n2_dn*cp*
             (dq_dk_dn*dk_dn_dn_dn + bdd2*dq2_dk_dn*dk_dn_dn_dn + q2dd*dbdd2_dn_dn);

   dfpp_dn_dn = pre_up*0.5*n2_up*cp*q2pp*dbpp2_dn_dn;
   dfdd_dn_up = pre_dn*0.5*n2_dn*cp*q2dd*dbdd2_dn_up;

/*==============================================================*/
/* Calculate df(parallel)/dgn                                   */

   dbpp2_dgn_up = -0.5*gamma2*gn_up*eee/n2_up;
   dbdd2_dgn_dn = -0.5*gamma2*gn_dn*eee/n2_dn;

   dfpp_dgn_up = pre_up*0.5*n2_up*cp*q2pp*dbpp2_dgn_up;
   dfdd_dgn_dn = pre_dn*0.5*n2_dn*cp*q2dd*dbdd2_dgn_dn;

/*==============================================================*/
/* Calculate df(parallel)/dln                                   */

   dbpp2_dln_up = 0.125*gamma2*eee/n_up;
   dbdd2_dln_dn = 0.125*gamma2*eee/n_dn;

   dfpp_dln_up = pre_up*0.5*n2_up*cp*q2pp*dbpp2_dln_up;
   dfdd_dln_dn = pre_dn*0.5*n2_dn*cp*q2dd*dbdd2_dln_dn;


/*==============================================================*/
/* Calculate df(parallel)/dtau                                  */

  dbpp2_dtau_up = gamma2*eee/n_up;
  dbdd2_dtau_dn = gamma2*eee/n_dn;

  dfpp_dtau_up = pre_up*0.5*n2_up*cp*
         (dq_dk_up*dk_up_dtau_up + bpp2*dq2_dk_up*dk_up_dtau_up + q2pp*dbpp2_dtau_up);

  dfdd_dtau_dn = pre_dn*0.5*n2_dn*cp*
         (dq_dk_dn*dk_dn_dtau_dn + bdd2*dq2_dk_dn*dk_dn_dtau_dn + q2dd*dbdd2_dtau_dn);

/*==============================================================*/
/* Total functional and its derivatives                         */

  *f          = fpd + fpp + fdd;
  *df_dn_up   = dfpd_dn_up + dfpp_dn_up + dfdd_dn_up;
  *df_dn_dn   = dfpd_dn_dn + dfpp_dn_dn + dfdd_dn_dn;
  *df_dgn_up  = dfpp_dgn_up;
  *df_dgn_dn  = dfdd_dgn_dn;
  *df_dln_up  = dfpp_dln_up;
  *df_dln_dn  = dfdd_dln_dn;
  *df_dtau_up = dfpd_dtau_up + dfpp_dtau_up;
  *df_dtau_dn = dfpd_dtau_dn + dfdd_dtau_dn;

/*==============================================================================*/
}/*end routine*/
/*==============================================================================*/

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void pbe_gcc_lsda(double rho_up, double rho_dn, double g_rho2_up, 
		  double g_rho2_dn,double drhox_up,double drhoy_up,
		  double drhoz_up,double drhox_dn,double drhoy_dn,
		  double drhoz_dn,double *fn,
		  double *df_dn_up, double *df_dgn_up, double *df_dn_dn,
		  double *df_dgn_dn,double beta,double gamma)

/*=======================================================================*/
/*         Begin Routine                                                 */
/*=======================================================================*/
/* Perdew-Burke-Ernzerhof generalized gradient approximation to the      */
/* density functional for exchange-correlation energy of a many-electron */
/* system. See PRL 77 (1996) 3865 and White & Bird PRB 50 (1994) 4954.   */
/* Change of parameters for extended PBE only   :                        */
/* revPBE see Zhang and Wang, PRL 80 (1998) 890.                         */
/* rPBE see Hammer, Hansen and Norskov PRB 59 (1999) 7413.               */
/* xPBE   see Xu and Goddard, JCP 121 (2004) 4068.(beta and gamma differ)*/
{/*Begin Routine*/
/*=======================================================================*/
/*local variables */

 double ec,dec_dn_up,dec_dn_dn;
 double rho,g_rho2,g_rho_up,g_rho_dn,g_rho,grup_dot_grdn;
 double zeta,phi,phi2,phi3;
 double kf,ks,t,t2,t3,t4;
 double expf,Afun,A2;
 double num,denom;
 double Hfun,Hlogterm;
 double rat;
 double Bfun;

 double dphi_dz,dphi_dn_up,dphi_dn_dn;
 double dz_dn_up,dz_dn_dn;
 double dt_dn,dt_dn_up,dt_dn_dn;
 double preA,dA_dn_up,dA_dn_dn;
 double dP_dt,dP_dA;
 double preH,dH_dn_up,dH_dn_dn;

 double dt_dgn_up,dt_dgn_dn;
 double dH_dgn_up,dH_dgn_dn;

/* static variables for lda functional         */
 static double aa = 3.0*M_PI*M_PI;
 static double twothird = 2.0/3.0;
/* static variables for gc part */

/* static double gamma = 0.031091;*/
/* static double gamma = 0.03109069086965489503494086371273;*/
/* static double beta = 0.066725;*/
/* static double beta = 0.06672455060314922;*/
 static double pow1 = 1.0/3.0;
 static double msixth = -1.0/6.0;
 static double eta = 1.0e-12;
 static double bet_by_gam;

/*==========================================================================*/
/* Set some variables to zero  */
  ec = 0.0;
  dec_dn_up = 0.0;
  dec_dn_dn = 0.0;

/*==========================================================================*/
/* Calculate some spin and density properties  */

  grup_dot_grdn = drhox_up*drhox_dn+drhoy_up*drhoy_dn+drhoz_up*drhoz_dn;
  g_rho2 = g_rho2_up + g_rho2_dn + 2.0*grup_dot_grdn;
  g_rho = sqrt(g_rho2);
  g_rho_up = sqrt(g_rho2_up);
  g_rho_dn = sqrt(g_rho2_dn);
  rho = rho_up + rho_dn;  
  zeta = (rho_up-rho_dn)/rho;
  phi = (pow(1.0+zeta,twothird)+pow(1.0-zeta,twothird))/2.0;
  phi2 = phi*phi;
  phi3 = phi*phi2;

/*==========================================================================*/
/* First do e_c and its derivative  */

  pw_c_lsda(rho_up,rho_dn,&ec,&dec_dn_up,&dec_dn_dn);
  dec_dn_up = (dec_dn_up -ec)/rho;
  dec_dn_dn = (dec_dn_dn - ec)/rho;

/*=========================================================================*/
/* Calculate the H function  */

  bet_by_gam = beta/gamma;
  kf = pow(aa*rho,pow1);
  ks = sqrt(4.0*kf/M_PI);
  t = g_rho/(2.0*phi*ks*rho);

/* Afun */

  expf = exp(-ec/(phi3*gamma)) - 1.0;
  Afun = bet_by_gam/expf;
  
/* rest of H */

  t2 = t*t;
  t3 = t2*t;
  t4 = t2*t2;
  A2 = Afun*Afun;
  num = 1.0 + Afun*t2;
  denom = 1.0 + Afun*t2 + A2*t4;
  rat = num/denom;
  Bfun = t2*bet_by_gam*rat;
  Hlogterm = log(1.0 + Bfun);
  Hfun = gamma*phi3*Hlogterm;
 
/*==========================================================================*/
/* Calculate the functional  */

  *fn = rho*Hfun;

/*==========================================================================*/
/* Calculate the derivative of Hfun wrt density  */

/* terms that are independent of up/down density */

/* prevents a blow-up when zeta -> 1.0 */
  dphi_dz = (pow((1.0+zeta)*(1.0+zeta)+eta, msixth) 
	    - pow((1.0-zeta)*(1.0-zeta)+eta, msixth))/3.0;
  dt_dn = -7.0*t/(6.0*rho);
  preH  =  beta*phi3/(1.0 + Bfun);
  dP_dt =  (2.0*t+4.0*Afun*t3)/(denom*denom);
  dP_dA = -t4*(2.0*t2*Afun+t4*A2)/(denom*denom);
  preA  = -bet_by_gam/(gamma*phi3)*(expf + 1.0)/(expf*expf);

/* terms that depend on up/down density */

  dz_dn_up   =  2.0*rho_dn/(rho*rho);
  dz_dn_dn   = -2.0*rho_up/(rho*rho);
  dphi_dn_up =  dphi_dz*dz_dn_up;
  dphi_dn_dn =  dphi_dz*dz_dn_dn;
  dt_dn_up   =  dt_dn - t*dphi_dn_up/phi;
  dt_dn_dn   =  dt_dn - t*dphi_dn_dn/phi;
  dA_dn_up   =  preA*(3.0*ec*dphi_dn_up/phi - dec_dn_up);
  dA_dn_dn   =  preA*(3.0*ec*dphi_dn_dn/phi - dec_dn_dn);

  dH_dn_up = preH*(dP_dt*dt_dn_up + dP_dA*dA_dn_up) + 3.0*Hfun/phi*dphi_dn_up;
  dH_dn_dn = preH*(dP_dt*dt_dn_dn + dP_dA*dA_dn_dn) + 3.0*Hfun/phi*dphi_dn_dn;
  
/*=========================================================================*/
/* Calculate the derivative of the functional wrt density  */

 *df_dn_up = Hfun + rho*dH_dn_up;
 *df_dn_dn = Hfun + rho*dH_dn_dn;

/*=========================================================================*/
/*  Calculate the derivative of the Hfun wrt |grad n|  */

  dt_dgn_up = t*(g_rho_up + grup_dot_grdn/g_rho_up)/g_rho2;
  dt_dgn_dn = t*(g_rho_dn + grup_dot_grdn/g_rho_dn)/g_rho2;
  dH_dgn_up = preH*dP_dt*dt_dgn_up;
  dH_dgn_dn = preH*dP_dt*dt_dgn_dn;

/*=========================================================================*/
/* Calculate the derivative of the functional wrt density  */
  
  *df_dgn_up = rho*dH_dgn_up;
  *df_dgn_dn = rho*dH_dgn_dn;

/*==========================================================================*/
}/*end routine*/
/*==========================================================================*/

/*=======================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

 void construct_elf(double *cp_elf,double *elec_ke_dens,
                    double *del_rho_x,double *del_rho_y,double *del_rho_z,
                    double *rho,double gc_cut,
                    PARA_FFT_PKG3D *cp_para_fft_pkg3d_lg)

/*==========================================================================*/
/*         Begin Routine                                                    */
   {/*Begin Routine*/
/*=======================================================================*/
/*local variables */

   int nfft            =    cp_para_fft_pkg3d_lg->nfft;
   int nfft_proc       =    cp_para_fft_pkg3d_lg->nfft_proc;
   int nfft2           =    nfft/2;
   int nfft2_proc      =    nfft_proc/2;
   int i,iii;
   double rho_cut;
   double g_rho_2;
   double d,d0;
   double chi;
   double integral=0;
   static double pow53,pow23,pre_fact,spi2;
   static double tol = 1.0e-30;

/*=======================================================================*/
/* I) Assign useful constants                                            */

   pow53 = 5.0/3.0;
   pow23 = 2.0/3.0;
   spi2 = 6.0*M_PI*M_PI;
   pre_fact = pow(spi2,pow23)/pow53;
   rho_cut = gc_cut;

/*=======================================================================*/
/* II) Construct the ELF                                                */


   for(i=1; i<= nfft2_proc; i++){
    if(rho[i] > rho_cut){
       g_rho_2 = del_rho_x[i]*del_rho_x[i]
               + del_rho_y[i]*del_rho_y[i]
               + del_rho_z[i]*del_rho_z[i];
       d  = elec_ke_dens[i] - 0.25*(g_rho_2/rho[i]);
       d0 = pre_fact*pow(rho[i],pow53);
       chi = d/d0;
       cp_elf[i] = 1.0/(1.0 + chi*chi);
    } else {
       cp_elf[i] = 0.5;
    }
   }/* endfor */

/*=======================================================================*/
}/*end routine*/
/*=======================================================================*/
    

/*=======================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

 void debug97x(double rho,double g_rho2,
               double *fn,double *df_dn,double *df_dgn)

/*==========================================================================*/
/*         Begin Routine                                                    */
   {/*Begin Routine*/
/*=======================================================================*/
/*local variables */
 
 double g_rho,x;

 g_rho = sqrt(g_rho2);
 x = g_rho/rho;
 *fn = rho*exp(-(x*x));
 *df_dn  =  (2.0*x*x + 1.0)*exp(-(x*x));
 *df_dgn = -2.0*x*exp(-(x*x));

/*=======================================================================*/
}/*end routine*/
/*=======================================================================*/


/*=======================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

 void debug_lap1(double rho,double g_rho2,double lap_rho,
                 double *fn,double *df_dn,double *df_dgn,double *df_dln)


/*==========================================================================*/
/*         Begin Routine                                                    */
   {/*Begin Routine*/
/*=======================================================================*/
/*local variables */

  *fn = lap_rho*exp(-rho);
  *df_dn = -lap_rho*exp(-rho);
  *df_dgn = 0.0;
  *df_dln = exp(-rho);


/*=======================================================================*/
 }/*end routine*/
/*=======================================================================*/




/*=======================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

 void debug_lap2(double rho,double g_rho2,
               double *fn,double *df_dn,double *df_dgn)

/*==========================================================================*/
/*         Begin Routine                                                    */
   {/*Begin Routine*/
/*=======================================================================*/
/*local variables */

  *fn = g_rho2*exp(-rho);
  *df_dn = -g_rho2*exp(-rho);
  *df_dgn = 2.0*sqrt(g_rho2)*exp(-rho);


/*=======================================================================*/
 }/*end routine*/
/*=======================================================================*/


/*=======================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void debug_tau(double rho,double g_rho2,double tau,
               double *fn,double *df_dn,double *df_dgn,double *df_dtau)

/*==========================================================================*/
/*         Begin Routine                                                    */
   {/*Begin Routine*/
/*=======================================================================*/
/*local variables */

 int iii;

 *fn = exp(-tau*tau); 
 *df_dn  =  0.0;
 *df_dgn = 0.0;
 *df_dtau = -2.0*tau*exp(-tau*tau); 

/*=======================================================================*/
}/*end routine*/
/*=======================================================================*/

