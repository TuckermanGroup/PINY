/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
/*                                                                          */
/*                         PI_MD:                                           */
/*             The future of simulation technology                          */
/*             ------------------------------------                         */
/*                     Module: spline_fit.c                                 */
/*                                                                          */
/* Bin the interactions on a 1D grid.                                       */
/*                                                                          */
/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/


#include "standard_include.h"
#include "../typ_defs/typedefs_par.h"
#include "../typ_defs/typedefs_class.h"
#include "../proto_defs/proto_inter_params_local.h"
#include "../proto_defs/proto_friend_lib_entry.h"


/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void spline_vdv(double rmin, double rmax, 
		double *cv0, double *cdv0, 
		int nsplin, double *dr, double *dri,double sig, double eps,
		double a, double b, double c, double rm, 
                double c6, double c8, double c9, double c10, 
		double alp, double qij, int iperd,int itype, int ishift,
                double rheal, int diele_opt, double diele_rheal, double 
                diele_cut,double diele_eps)

/*==========================================================================*/ 
/*           Begin Routine                                                  */
{ /*begin routine*/
/*==========================================================================*/
/*           Local Variables                                                */

  int i,iii;
  double *r,*sw,*dsw;
  double rp,shift;

/*==========================================================================*/
/* 0) Malloc some temporary memory                                          */

  r = (double *)cmalloc(nsplin*sizeof(double))-1;
  sw = (double *)cmalloc(nsplin*sizeof(double))-1;
  dsw = (double *)cmalloc(nsplin*sizeof(double))-1;

/*==========================================================================*/
/* II) get positions at which to evaluate potential and force               */
  
  *dr = (rmax-rmin) /(double)(nsplin - 5);
  *dri = 1.0/(*dr);

  for (i = 1; i <= nsplin; ++i) {
    r[i] = (*dr) * ((double)(i-3)) + rmin;
    sw[i] = 1.0;
    dsw[i] = 0.0;  
  }/*endfor*/
  if(ishift==2){
    for (i = 1; i <= nsplin; ++i) {
      rp    = (r[i]-rmax+rheal)/rheal;
      rp    = MIN(rp,1.0);
      rp    = MAX(rp,0.0);
      sw[i] = 1.0 + rp*rp*(2.0*rp-3.0);
     dsw[i] = -6.0*(rp/rheal)*(rp-1.0)/r[i];
    }/*endfor*/
  }/*endif*/

/*==========================================================================*/
/* III) Bin the potential energy                                            */
/*       Recommended range is 1.0 to 20a0 (or cutoff)                       */

  switch(itype){
   /* Coulomb */
    case 1: vcoul_bin(nsplin,r,cv0,cdv0,sw,dsw,qij,iperd,alp,diele_opt,
                      diele_rheal,diele_cut,diele_eps); break; 
   /* LJ      */
    case 2: vlj_bin(nsplin,r,cv0,cdv0,sw,dsw,eps,sig);  break; 
   /* Exp-6 */
    case 3: vwill_bin(nsplin,r,cv0,cdv0,sw,dsw,a,b,c6,c8,c10); break;
   /* Null potential */
    case 4: vnull_bin(nsplin,r,cv0,cdv0); break;  
   /* Exp-6 plus LJ */
    case 5: vwill_lj_bin(nsplin,r,cv0,cdv0,sw,dsw,a,b,c6,c8,c10,eps,sig);break;
   /* Aziz-Chen potential */
    case 6: vaziz_bin(nsplin,r,cv0,cdv0,sw,dsw,a,b,c,rm,c6,c8,c9,c10); break;
  }/*end switch*/
/*==========================================================================*/
/* IV) Shift potential energy                                            */

  if (ishift == 1) {
    shift  = cv0[(nsplin-2)];
    for (i = 1; i <= nsplin; ++i) {
      cv0[i] -= shift;
    }/*endfor*/
  }/*endif*/

/*==========================================================================*/
/* V) Free the memory                                                       */

  cfree(&r[1]);
  cfree(&sw[1]);
  cfree(&dsw[1]);

/*--------------------------------------------------------------------------*/
 }/* end routine */ 
/*==========================================================================*/




/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void vnull_bin(int n, double r[], double v[],double dv[])

/*==========================================================================*/
{ /*begin routine */
/*==========================================================================*/
/*      Local variables       */

  int i;

/*==========================================================================*/
/* Null Potential */

  for(i=1;i<=n;i++){
    v[i]  = 0.0;
    dv[i] = 0.0;
  }/*endfor*/

/*--------------------------------------------------------------------------*/
} /* end routine */
/*==========================================================================*/



/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void vlj_bin(int n, double r[], double v[],double dv[], double sw[],
             double dsw[],double eps, double sig)

/*==========================================================================*/
{ /*begin routine */
/*==========================================================================*/
/*      Local variables       */

  int i;
  double vnow,dvnow;
  double sig6,sig6r;
  double r2;

/*==========================================================================*/
/* Lennard Jones Potential */

  sig6 = sig*sig*sig*sig*sig*sig;
  for(i=1;i<=n;i++){
    r2     = r[i] * r[i];
    sig6r  = sig6 / (r2 * r2 * r2);
    vnow   = (eps * 4. * sig6r * (sig6r - 1.));
    dvnow  = (eps * 24. * sig6r * (sig6r * 2. - 1.) / r2);
    v[i]   = vnow*sw[i];
    dv[i]  = dvnow*sw[i] + vnow*dsw[i];
  }/*endfor*/

/*--------------------------------------------------------------------------*/
} /* end routine */
/*==========================================================================*/



/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void vwill_bin(int n,double r[],double v[],double dv[],
               double sw[],double dsw[],double a,double b,
               double c6,double c8,double c10)

/*==========================================================================*/
{/*begin routine */
/*==========================================================================*/
/*      Local variables       */

  int i;
  double vnow,dvnow;
  double r2, r6, r8, r10;

/*==========================================================================*/
/* Williams Potential */

  for(i=1;i<=n;i++){
    r2    = r[i] * r[i];
    r6    = r2 * r2 * r2;
    r8    = r2 * r2 * r2 * r2;
    r10   = r8 * r2;
    vnow  = a*exp(-b*r[i]) - (c6/r6 + c8/r8 + c10/r10);
    dvnow = a*b*exp(-b*r[i])/r[i] - (c6 * 6./r6 + c8*8./r8 + c10*10./r10)/r2;
    v[i]  = vnow*sw[i];
    dv[i] = dvnow*sw[i] + vnow*dsw[i];
  }/*endfor*/

/*--------------------------------------------------------------------------*/
} /* end routine */
/*==========================================================================*/



/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void vwill_lj_bin(int n,double r[],double v[],double dv[],
               double sw[],double dsw[],double a,double b,
               double c6,double c8,double c10,double eps,double sig)

/*==========================================================================*/
{/*begin routine */
/*==========================================================================*/
/*      Local variables       */

  int i;
  double vnow,dvnow;
  double sig6,sig6r;
  double r2, r6, r8, r10;

/*==========================================================================*/
/* LJ+Williams Potential */

  sig6 = sig*sig*sig*sig*sig*sig;
  for(i=1;i<=n;i++){
    r2    = r[i] * r[i];
    r6    = r2 * r2 * r2;
    r8    = r2 * r2 * r2 * r2;
    r10   = r8 * r2;
    sig6r  = sig6 / (r2 * r2 * r2);
    vnow  = (a*exp(-b*r[i]) - (c6/r6 + c8/r8 + c10/r10))
          + (eps * 4. * sig6r * (sig6r - 1.));
    dvnow = (a*b*exp(-b*r[i])/r[i] - (c6 * 6./r6 + c8*8./r8 + c10*10./r10)/r2)
          + (eps * 24. * sig6r * (sig6r * 2. - 1.) / r2);
    v[i]  = vnow*sw[i];
    dv[i] = dvnow*sw[i] + vnow*dsw[i];
  }/*endfor*/

/*--------------------------------------------------------------------------*/
} /* end routine */
/*==========================================================================*/




/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void vaziz_bin(int n,double r[],double v[],double dv[],
               double sw[],double dsw[],double a,double b,double c,
               double rm,double c6,double c8,double c9,double c10)

/*==========================================================================*/
{/*begin routine */
/*==========================================================================*/
/*      Local variables       */

  int i;
  double vnow,dvnow;
  double dvtmp,r2,r6,r8,r9,r10,rs;
  double arg1,arg2,temp1,temp2,temp3;

/*==========================================================================*/
/* Aziz Chen Potential */

  for(i=1;i<=n;i++){
   rs    = MIN(r[i],rm);
   r2    = r[i]*r[i];
   r6    = r2*r2*r2;
   r8    = r6*r2;
   r9    = r8*r[i];
   r10   = r8*r2;
   arg1  = -b*r[i]-c*r2;
   temp1 = a*exp(arg1);
   arg2  = -((rm/rs-1.0)*(rm/rs-1.0));
   temp2 = exp(arg2);
   temp3 = (c6/r6+c8/r8-c9/r9+c10/r10);
   vnow  =  (temp1-temp2*temp3);
   dvtmp = (b+2.0*c*r[i])*temp1
         -(6.0*c6/r6+8.0*c8/r8-9.0*c9/r9+10.0*c10/r10)*temp2/r[i]
         +temp3*temp2*2.0*(rm/rs-1.0)*rm/(rs*rs);
   dvnow = dvtmp/r[i];
   v[i]  = vnow*sw[i];
   dv[i] = dvnow*sw[i] + vnow*dsw[i];
  }/*endfor*/  

/*--------------------------------------------------------------------------*/
} /* end routine */
/*==========================================================================*/




/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void vcoul_bin(int n,double r[],double v[],double dv[],
               double sw[],double dsw[],
               double qij,int iperd, double alp, int diele_opt, 
               double diele_rheal, double diele_cut,double diele_eps)

/*==========================================================================*/
{/*begin routine */
/*==========================================================================*/
/*         Local Variables */

  int i,iii;
  double vnow,dvnow;
  double palp, ralp, talp2, gerfc, r2, tt, eee,dgerfc;
  double p=0.3614;
  double e1 = 0.2041422096422003, e2 = 0.1997535956961481;
  double e3 = 0.2213176596405576, e4 = 0.03360430734640255;
  double e5 = 0.4732592578721755, e6 =-0.509078520069735;
  double e7 = 0.6772631491947646, e8 =-0.369912979092217;
  double e9 = 0.06965131976970335;
  double de1,de2,de3;
  double de4,de5,de6;
  double de7,de8,de9;
  double rt,s_eps,eps_r,deps_r;
  de1 = 1.0*e1; de2 = 2.0*e2; de3 = 3.0*e3;
  de4 = 4.0*e4; de5 = 5.0*e5; de6 = 6.0*e6;
  de7 = 7.0*e7; de8 = 8.0*e8; de9 = 9.0*e9;

/*==========================================================================*/
/* Regular old Mr. c */

  if (iperd == 0) {

   if(diele_opt==1){
    for(i=1;i<=n;i++){
      rt = (r[i]-diele_cut+diele_rheal)/diele_rheal;
      rt = MAX(rt,0.0);
      rt = MIN(rt,1.0);
      s_eps = rt*rt*(3.0-2.0*rt);
      eps_r  = 1.0 + s_eps*(diele_eps-1.0);
      vnow  = qij/(r[i]*eps_r);
      deps_r = (diele_eps - 1.0)*(6.0*rt*(1.0-rt))/diele_rheal;
      dvnow = qij*(eps_r+r[i]*deps_r)/(r[i]*r[i]*r[i]*eps_r*eps_r);
      v[i]  = vnow*sw[i];
      dv[i] = dvnow*sw[i] + vnow*dsw[i];
    }/*endfor*/
   }/*endif*/

   if(diele_opt==0){
    for(i=1;i<=n;i++){
      vnow  = qij/r[i];
      dvnow = qij/(r[i]*r[i]*r[i]);
      v[i]  = vnow*sw[i];
      dv[i] = dvnow*sw[i] + vnow*dsw[i];

    }/*endfor*/
  }/*endif*/

 }/*endif*/

/*==========================================================================*/
/* Real space ewald */

  if (iperd > 0) {
    talp2 = 2.0*alp*alp;
    palp  = p*alp;
    for(i=1;i<=n;i++){
      r2     = r[i]*r[i];
      ralp   = r[i] * alp;
      eee    = exp(-ralp*ralp);
      tt     = 1.0/(1.0+p*ralp);
      gerfc  = ((((((((e9*tt+e8)*tt+e7)*tt+e6)*tt+e5)*tt
                      +e4)*tt+e3)*tt+e2)*tt+e1)*tt*eee;
      dgerfc = ((((((((de9*tt+de8)*tt+de7)*tt+de6)*tt+de5)*tt
                           +de4)*tt+de3)*tt+de2)*tt+de1)*tt*tt*eee*palp
                           +talp2*gerfc*r[i];
      vnow   = qij * gerfc/r[i];
      dvnow  = (gerfc/r2 + dgerfc/r[i])*qij/r[i];
      v[i]   = vnow*sw[i];
      dv[i]  = dvnow*sw[i] + vnow*dsw[i];
    }/*endfor*/
  }/*endif*/

/*--------------------------------------------------------------------------*/
} /* end routine */
/*==========================================================================*/
