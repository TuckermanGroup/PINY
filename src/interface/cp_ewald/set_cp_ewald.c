/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
/*                                                                          */
/*                         PI_MD:                                           */
/*             The future of simulation technology                          */
/*             ------------------------------------                         */
/*                     Module: set_ewald                                    */
/*                                                                          */
/* This reads in and sets up the electron-atom interaction pseudopotential  */
/*                                                                          */
/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

#include "standard_include.h"
#include "../typ_defs/typedefs_gen.h"
#include "../typ_defs/typedefs_class.h"
#include "../typ_defs/typedefs_cp.h"
#include "../proto_defs/proto_cp_ewald_local.h"
#include "../proto_defs/proto_math.h"
#include "../proto_defs/proto_friend_lib_entry.h"

#define CP_EMAGIC 0.00050


/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void set_cpmass(int ncoef,int *kastore,int *kbstore,int *kcstore,
                double *cmass,double *hmati,
                double *cmass_tau_def,double cmass_cut_def,
                int *icmass_unif)

/*==========================================================================*/
/*            Begin subprogram:                                          */
    {/*begin routine*/
/*=======================================================================*/
/*            Local variable declarations:                               */

  double aka,akb,akc,xk,yk,zk;              /* Num: k-vector components */
  double tpi,enow,cmass0;                /* Num: Some useful consts  */
  int i,iii;                              /* Num: For loop counter    */

/*=======================================================================*/
/*  A) Set up CP masses                                                  */

  *icmass_unif = 1;
  tpi = 2.0*M_PI;
  *cmass_tau_def /= TIME_CONV;
  cmass0 = 4.0*(*cmass_tau_def)*(*cmass_tau_def)*CP_EMAGIC;
  for(i=1;i<=ncoef-1;i++) {
    aka = (double)(kastore[i]);
    akb = (double)(kbstore[i]);
    akc = (double)(kcstore[i]);
    xk = (aka*hmati[1]+akb*hmati[2]+akc*hmati[3])*tpi;
    yk = (aka*hmati[4]+akb*hmati[5]+akc*hmati[6])*tpi;
    zk = (aka*hmati[7]+akb*hmati[8]+akc*hmati[9])*tpi;
    enow = (xk*xk+yk*yk+zk*zk)*0.5;
    if(enow > (0.5*cmass_cut_def)) {
      cmass[i] =  cmass0*enow/(0.5*cmass_cut_def);
      *icmass_unif = 0;
    } else {
      cmass[i] =  cmass0;
    } /* endif */
  } /* endfor */

  cmass[ncoef] = 0.5*cmass0;


/*-------------------------------------------------------------------------*/
} /* end routine */
/*=======================================================================*/






/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

 void check_kvec(int nktot,int kastore[],int kbstore[],int kcstore[],
                 int nktot_sm,
                 int kastore_sm[],int kbstore_sm[],int kcstore_sm[])

/*==========================================================================*/
/*               Begin subprogram:                                          */
      {/*begin routine*/
/*==========================================================================*/
/*               Local variable declarations:                               */
      int i,icount;                           /* Num: Counters              */
/*==========================================================================*/
      icount = 0;

      for(i=1;i<=nktot;i++) {
       if((kastore[i] == kastore_sm[icount+1]) &&
          (kbstore[i] == kbstore_sm[icount+1]) &&
          (kcstore[i] == kcstore_sm[icount+1])) {
        icount++;
       } /* endif */
      } /* endfor */
      if(icount != nktot_sm) {
        printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
        printf("Mismatch in large and small kvectors\n");
        printf("%d vs %d\n",icount,nktot_sm);
        printf("Contact technical support\n");
        printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
        fflush(stdout);
        exit(1);
      } /* endif */

/*-------------------------------------------------------------------------*/
     } /* end routine */
/*==========================================================================*/



/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void calc_cutoff(int kmax_ewd, double *ecut,double *ecut_cp,int cp_on,
                 int *kmax_cp, int *kmaxv, double *hmatik, double deth)

/*==========================================================================*/
/*               Begin subprogram:                                          */
      {/*begin routine*/
/*==========================================================================*/
/*               Local variable declarations:                               */

  int iii;
  double rtwoth,tpi,rvol23;
  double d1,d2,d3;
  double try1,try2,try3;
  double temp1,temp2,temp3;

/*==========================================================================*/
/* II) Calculate Ewald Cutoff */

   rtwoth = -(2./3.);
   rvol23 = pow(deth,rtwoth);
   tpi = M_PI * 2.0;
   *ecut = M_PI * .5 * M_PI * (double) (kmax_ewd * kmax_ewd)*rvol23;

/*==========================================================================*/
/* III) Compare Ewald to CP cutoff  */

   if (cp_on == 1) {
     if (*ecut_cp * .5 < *ecut) {
        printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
	printf("Ewald cutoff greater than cp cutoff %g vs %g \n",
                                     (*ecut)*2.,*ecut_cp);
        printf("Therefore, you must lower the maximum k-vector\n");
        printf("required by your Ewald sum, perhaps change alp_ewd\n");
        printf("and try again\n");
        printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
        fflush(stdout);
        exit(1);
     }/*endif*/
     d1 = *ecut, d2 = *ecut_cp * .5;
     *ecut = MAX(d1,d2);
   }/*endif*/
   *ecut_cp = *ecut;

/*==========================================================================*/
/* III) Adjust shape of reciprocal space                                    */

   d1 = hmatik[1];  d2 = hmatik[4];   d3 = hmatik[7];
   try1 = sqrt(d1 * d1 + d2 * d2 + d3 * d3);
   d1 = hmatik[2];  d2 = hmatik[5];   d3 = hmatik[8];
   try2 = sqrt(d1 * d1 + d2 * d2 + d3 * d3);
   d1 = hmatik[3];  d2 = hmatik[6];   d3 = hmatik[9];
   try3 = sqrt(d1 * d1 + d2 * d2 + d3 * d3);
   temp1 = sqrt(*ecut * .5) / (M_PI * try1);
   temp2 = sqrt(*ecut * .5) / (M_PI * try2);
   temp3 = sqrt(*ecut * .5) / (M_PI * try3);
   if(cp_on==1){
      kmax_cp[1] = NINT(temp1);
      kmax_cp[2] = NINT(temp2);
      kmax_cp[3] = NINT(temp3);
      radixme(&(kmax_cp[1]),&(kmax_cp[2]),&(kmax_cp[3]));
      kmaxv[1] = kmax_cp[1] << 1;  /* <<1 multiply by 2 */
      kmaxv[2] = kmax_cp[2] << 1;
      kmaxv[3] = kmax_cp[3] << 1;
   }else{
      kmaxv[1] = (int) ( 2.0*temp1);
      kmaxv[2] = (int) ( 2.0*temp2);
      kmaxv[3] = (int) ( 2.0*temp3);
    }/*endif*/

/*-------------------------------------------------------------------------*/
     } /* end routine */
/*==========================================================================*/



/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void set_pme_grid(double ecut_now,double deth,double *hmatik,int *kmaxv,
                  int *ngrid_a,int *ngrid_b,int *ngrid_c,int n_interp,
                  int kmax_pme)

/*==========================================================================*/
/*               Begin subprogram:                                          */
      {/*begin routine*/
/*==========================================================================*/
/*               Local variable declarations:                               */

  int iii,igo;
  double rtwoth,tpi,rvol23;
  double d1,d2,d3;
  double try1,try2,try3;
  double temp1,temp2,temp3,ecut_pme,ecut;
  int ktemp1,ktemp2,ktemp3;


/*==========================================================================*/
/* IV) Calculate PME cutoff                                                 */

   rtwoth = -(2./3.);
   rvol23 = pow(deth,rtwoth);
   tpi = M_PI * 2.0;
   ecut_pme = M_PI * .5 * M_PI * (double) (kmax_pme * kmax_pme)*rvol23;
   if(ecut_pme < ecut_now){
    printf("$$$$$$$$$$$$$$$$$$$$_Warning_$$$$$$$$$$$$$$$$$$$$\n");
    printf("Warning PME cutoff taken too small             \n");
    printf("Using an appropriate larger value.             \n");
    printf("$$$$$$$$$$$$$$$$$$$$_Warning_$$$$$$$$$$$$$$$$$$$$\n");
    fflush(stdout);
   }/*endif*/
   ecut = MAX(ecut_pme,ecut_now);
  

/*==========================================================================*/
/* III) Adjust shape of reciprocal space  */

   d1 = hmatik[1];  d2 = hmatik[4];   d3 = hmatik[7];
   try1 = sqrt(d1 * d1 + d2 * d2 + d3 * d3);
   d1 = hmatik[2];  d2 = hmatik[5];   d3 = hmatik[8];
   try2 = sqrt(d1 * d1 + d2 * d2 + d3 * d3);
   d1 = hmatik[3];  d2 = hmatik[6];   d3 = hmatik[9];
   try3 = sqrt(d1 * d1 + d2 * d2 + d3 * d3);
   temp1 = sqrt(ecut * .5) / (M_PI * try1);
   temp2 = sqrt(ecut * .5) / (M_PI * try2);
   temp3 = sqrt(ecut * .5) / (M_PI * try3);
   ktemp1 = NINT(temp1);
   ktemp2 = NINT(temp2);
   ktemp3 = NINT(temp3);
   radixme(&ktemp1,&ktemp2,&ktemp3);
   *ngrid_a = 4*(ktemp1 + 1);
   *ngrid_b = 4*(ktemp2 + 1);
   *ngrid_c = 4*(ktemp3 + 1);
   igo=0;
   if(*ngrid_a < 2*(kmaxv[1]+2)){igo=1;}
   if(*ngrid_b < 2*(kmaxv[2]+2)){igo=1;}
   if(*ngrid_c < 2*(kmaxv[3]+2)){igo=1;}
   while(igo==1){
       if(*ngrid_a < 2*(kmaxv[1]+2)){temp1+=1.0;}
       if(*ngrid_b < 2*(kmaxv[2]+2)){temp2+=1.0;}
       if(*ngrid_c < 2*(kmaxv[3]+2)){temp3+=1.0;}
       ktemp1 = NINT(temp1);
       ktemp2 = NINT(temp2);
       ktemp3 = NINT(temp3);
       radixme(&ktemp1,&ktemp2,&ktemp3);
       *ngrid_a = 4*(ktemp1 + 1);
       *ngrid_b = 4*(ktemp2 + 1);
       *ngrid_c = 4*(ktemp3 + 1);
       igo = 0;
       if(*ngrid_a < 2*(kmaxv[1]+2)){igo=1;}
       if(*ngrid_b < 2*(kmaxv[2]+2)){igo=1;}
       if(*ngrid_c < 2*(kmaxv[3]+2)){igo=1;}
   }/*endwhile*/
   if((n_interp > *ngrid_a) || 
      (n_interp > *ngrid_b) || 
      (n_interp > *ngrid_c) ){
       printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
       printf("PME parameter n_interp too large for pme cutoff\n");
       printf("%d > %d or %d or %d \n",
            n_interp, *ngrid_a, *ngrid_b, *ngrid_c);
       printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
       fflush(stdout);
       exit(1);
   }/*endif*/
   

/*-------------------------------------------------------------------------*/
     } /* end routine */
/*==========================================================================*/




/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
/* Count the number of k vectors on large grid */
/*==========================================================================*/

void countkvec3d(int *nktot,double ecut,int *kmaxv,double *hmatik)

/*==========================================================================*/
/*       Begin routine */
{/*begin routine */
/*==========================================================================*/
/* Local variables */

  int iii,icount;
  int i, kbmin, kcmin, kbmax, kcmax, kamax, ka, kb, kc;
  double xk, yk, zk;
  double aka, akb, akc;
  double tpi,try;

/*==========================================================================*/
/* Count the kvectors */

  tpi = 2.0*M_PI;
  icount = 0;
  kamax = kmaxv[1];

/*********************************/

   for (i = 1; i <= kamax; ++i) {
    aka = (double) i;
    xk = aka * hmatik[1] * tpi;
    yk = aka * hmatik[4] * tpi;
    zk = aka * hmatik[7] * tpi;
    try = (xk * xk + yk * yk + zk * zk) * .5;
    if (try > ecut * 4.) {
      break;
    }
  }

  kamax = i - 1;

/***********************************/

  for (ka = 0; ka <= kamax; ++ka) {
    aka = (double) ka;
    kbmin = -kmaxv[2];
    if (ka == 0) {
      kbmin = 0;
    }
    for (i = kbmin; i <= 0; ++i) {
      akb = (double) i;
      xk = (aka * hmatik[1] + akb * hmatik[2]) * tpi;
      yk = (aka * hmatik[4] + akb * hmatik[5]) * tpi;
      zk = (aka * hmatik[7] + akb * hmatik[8]) * tpi;
      try = (xk * xk + yk * yk + zk * zk) * .5;
      if (try <= ecut * 4.) {
	break;
      }
    }

/*********************************/

    kbmin = i;
    for (i = 1; i <= kmaxv[2]; ++i) {
      akb = (double) i;
      xk = (aka * hmatik[1] + akb * hmatik[2]) * tpi;
      yk = (aka * hmatik[4] + akb * hmatik[5]) * tpi;
      zk = (aka * hmatik[7] + akb * hmatik[8]) * tpi;
      try = (xk * xk + yk * yk + zk * zk) * .5;
      if (try > ecut * 4.) {
	break;
      }
    }
    
    kbmax = i - 1;
    for (kb = kbmin; kb <= kbmax; ++kb) {
      
      akb = (double) kb;
      kcmin = -kmaxv[3];
      if (ka == 0 && kb == 0) {
	kcmin = 1;
      }
      for (i = kcmin; i <= 0; ++i) {
	akc = (double) i;
	xk = (aka * hmatik[1] + akb * hmatik[2] + akc * hmatik[3]) * tpi;
	yk = (aka * hmatik[4] + akb * hmatik[5] + akc * hmatik[6]) * tpi;
	zk = (aka * hmatik[7] + akb * hmatik[8] + akc * hmatik[9]) * tpi;
	try = (xk * xk + yk * yk + zk * zk) * .5;
	if (try <= ecut * 4.) {
	  break;
	}
      }
/*********************************/

      kcmin = i;
      for (i = 1; i <= kmaxv[3]; ++i) {
	akc = (double) i;
	xk = (aka * hmatik[1] + akb * hmatik[2] + akc * hmatik[3]) * tpi;
	yk = (aka * hmatik[4] + akb * hmatik[5] + akc * hmatik[6]) * tpi;
	zk = (aka * hmatik[7] + akb * hmatik[8] + akc * hmatik[9]) * tpi;
	try = (xk * xk + yk * yk + zk * zk) * .5;
	if (try > ecut * 4.) {
	  break;
	}
      }

      kcmax = i - 1;
      akc = (double) kcmin;
      for (kc = kcmin; kc <= kcmax; ++kc) {
	++icount;
      }
    }
  }
  *nktot = icount;
/*--------------------------------------------------------------------------*/
} /* countkvec3d */
/*==========================================================================*/



/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
/* This subroutine determines the allowed spherically truncated             */
/* half space k (.i.e. g) vectors given a cutoff. it is used by             */
/* both the ewald and cp modules.                                           */
/* Sets up the k vectors for a given cutoff and shape                       */
/*==========================================================================*/

void setkvec3d(int nktot,double ecut,int *kmaxv,double *hmatik,
		int *kastore, int *kbstore, int *kcstore, 
                int *ibreak1, int *ibreak2, int cp_on, 
                double *gmin_spl, double *gmin_true,double *gmax_spl)

/*==========================================================================*/
/*       Begin routine */
{/*begin routine */
/*==========================================================================*/
/* Local variables */

    int iii;
    int i1, i2, i3;

    int i, kbmin, kcmin, kbmax, kcmax, kamax, ka, kb, kc;
    double xk, yk, zk;
    int icount;
    double aka, akb, akc;
    double tpi;
    double try;

/*==========================================================================*/
/* Count the K-vectors */

    tpi = 2.0*M_PI;
    for(i=1;i<=nktot;i++){
      ibreak1[i] = 0;
      ibreak2[i] = 0;
    }/*endfor*/
    icount = 0;
    (*gmin_spl) = 1.0e10;
    (*gmax_spl) = 0.0;

/*=============================*/

    kamax = kmaxv[1];
    i1 = kamax;
    for (ka = 0; ka <= i1; ++ka) {
      aka = (double) ka;
      kbmin = -kmaxv[2];
      if (ka == 0) {
	kbmin = 0;
      }
      for (i = kbmin; i <= 0; ++i) {
	akb = (double) i;
	xk = (aka * hmatik[1] + akb * hmatik[2]) * tpi;
	yk = (aka * hmatik[4] + akb * hmatik[5]) * tpi;
	zk = (aka * hmatik[7] + akb * hmatik[8]) * tpi;
	try = (xk * xk + yk * yk + zk * zk) * .5;
	if (try <= ecut * 4.) {
	  break;
	}
      }
      kbmin = i;
      i2 = kmaxv[2];
      for (i = 1; i <= i2; ++i) {
	akb = (double) i;
	xk = (aka * hmatik[1] + akb * hmatik[2]) * tpi;
	yk = (aka * hmatik[4] + akb * hmatik[5]) * tpi;
	zk = (aka * hmatik[7] + akb * hmatik[8]) * tpi;
	try = (xk * xk + yk * yk + zk * zk) * .5;
	if (try > ecut * 4.) {
	  break;
	}
      }

/*=============================*/

      kbmax = i - 1;
      i2 = kbmax;
      for (kb = kbmin; kb <= i2; ++kb) {
	akb = (double) kb;
	kcmin = -kmaxv[3];
	if (ka == 0 && kb == 0) {
	  kcmin = 1;
	}
	for (i = kcmin; i <= 0; ++i) {
	  akc = (double) i;
	  xk = (aka * hmatik[1] + akb * hmatik[2] + akc * hmatik[3]) * tpi;
	  yk = (aka * hmatik[4] + akb * hmatik[5] + akc * hmatik[6]) * tpi;
	  zk = (aka * hmatik[7] + akb * hmatik[8] + akc * hmatik[9]) * tpi;
	  try = (xk * xk + yk * yk + zk * zk) * .5;
	  if (try <= ecut * 4.) {
	    break;
	  }
	}

/*=============================*/

	kcmin = i;
	i3 = kmaxv[3];
	for (i = 1; i <= i3; ++i) {
	  akc = (double) i;
	  xk = (aka * hmatik[1] + akb * hmatik[2] + akc * hmatik[3]) * tpi;
	  yk = (aka * hmatik[4] + akb * hmatik[5] + akc * hmatik[6]) * tpi;
	  zk = (aka * hmatik[7] + akb * hmatik[8] + akc * hmatik[9]) * tpi;
	  try = (xk * xk + yk * yk + zk * zk) * .5;
	  if (try > ecut * 4.) {
	    break;
	  }
	}
	kcmax = i - 1;
	i3 = kcmax;
	for (kc = kcmin; kc <= i3; ++kc) {
	  ++icount;
	  aka = (double) ka;
	  akb = (double) kb;
	  akc = (double) kc;
	  xk = (aka * hmatik[1] + akb * hmatik[2] + akc * hmatik[3]) * tpi;
	  yk = (aka * hmatik[4] + akb * hmatik[5] + akc * hmatik[6]) * tpi;
	  zk = (aka * hmatik[7] + akb * hmatik[8] + akc * hmatik[9]) * tpi;
	  kastore[icount] = ka;
	  kbstore[icount] = kb;
	  kcstore[icount] = kc;
	  try = sqrt(xk * xk + yk * yk + zk * zk);
	  (*gmin_spl) = MIN(try,(*gmin_spl));
	  (*gmax_spl) = MAX(try,(*gmax_spl));
	  if (kc == kcmin) {
	    ibreak1[icount] = 1;
	  }
	  if (kc < kcmax) {
	    ibreak2[icount] = 1;
	  }
	}
      }
    }
    if(nktot!=icount){
        printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
        printf("Mismatch number of kvectors\n");
        printf("%d vs %d\n",icount,nktot);
        printf("Contact technical support\n");
        printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
        fflush(stdout);
        exit(1);
    }
    if (cp_on == 1) {
      ++icount;
      kastore[icount] = 0;
      kbstore[icount] = 0;
      kcstore[icount] = 0;
      (*gmin_true) = (*gmin_spl);
      (*gmin_spl) *= 0.75;
      (*gmax_spl) *= (4.0/3.0);
    }
/*--------------------------------------------------------------------------*/
} /* setkvec3d */
/*==========================================================================*/




/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
/* This subroutine determines the allowed spherically truncated             */
/* half space k (.i.e. g) vectors given a cutoff. It is used by             */
/* cp modules. Sets up the k vectors for a given cutoff and shape           */
/*==========================================================================*/

void setkvec3d_sm(int nktot,double ecut,int *kmax_cp,double *hmatik,
                  int *kastore, int *kbstore, int *kcstore, 
		  int *ibreak1, int *ibreak2, double *gmin, double *gmax)

/*==========================================================================*/
/*       Begin routine */
{/*begin routine */
/*==========================================================================*/
/* Local variables */

  int i1, i2, i3;

  int i, kbmin, kcmin, kbmax, kcmax, kamax, ka, kb, kc;
  double xk, yk, zk,g;
  int icount;
  double aka, akb, akc;
  double tpi, try;
  tpi = M_PI * 2.;

/*==========================================================================*/
/* SETUP THE KVECTORS */
  
  for(i=1;i<=nktot;i++){
      ibreak1[i] = 0;
      ibreak2[i] = 0;
  }

  *gmin = 10000.0;
  *gmax = .0;
  icount = 0;

  i1 = kmax_cp[1];
  for (i = 1; i <= i1; ++i) {
    aka = (double) i;
    xk = aka * hmatik[1] * tpi;
    yk = aka * hmatik[4] * tpi;
    zk = aka * hmatik[7] * tpi;
    try = (xk * xk + yk * yk + zk * zk) * .5;
    if (try > ecut) {
      break;
    }
  }

  kamax = i - 1;
  i1 = kamax;
  for (ka = 0; ka <= i1; ++ka) {
    aka = (double) ka;
    kbmin = -kmax_cp[2];
    if (ka == 0) {
      kbmin = 0;
    }
    for (i = kbmin; i <= 0; ++i) {
      akb = (double) i;
      xk = (aka * hmatik[1] + akb * hmatik[2]) * tpi;
      yk = (aka * hmatik[4] + akb * hmatik[5]) * tpi;
      zk = (aka * hmatik[7] + akb * hmatik[8]) * tpi;
      try = (xk * xk + yk * yk + zk * zk) * .5;
      if (try <= ecut) {
	break;
      }
    }
    kbmin = i;
    i2 = kmax_cp[2];
    for (i = 1; i <= i2; ++i) {
      akb = (double) i;
      xk = (aka * hmatik[1] + akb * hmatik[2]) * tpi;
      yk = (aka * hmatik[4] + akb * hmatik[5]) * tpi;
      zk = (aka * hmatik[7] + akb * hmatik[8]) * tpi;
      try = (xk * xk + yk * yk + zk * zk) * .5;
      if (try > ecut) {
	break;
      }
    }

    kbmax = i - 1;
    i2 = kbmax;
    for (kb = kbmin; kb <= i2; ++kb) {
      akb = (double) kb;
      kcmin = -kmax_cp[3];
      if (ka == 0 && kb == 0) {
	kcmin = 1;
      }
      for (i = kcmin; i <= 0; ++i) {
	akc = (double) i;
	xk = (aka * hmatik[1] + akb * hmatik[2] + akc * hmatik[3]) * tpi;
	yk = (aka * hmatik[4] + akb * hmatik[5] + akc * hmatik[6]) * tpi;
	zk = (aka * hmatik[7] + akb * hmatik[8] + akc * hmatik[9]) * tpi;
	try = (xk * xk + yk * yk + zk * zk) * .5;
	if (try <= ecut) {
	  break;
	}
      }
      
      kcmin = i;
      i3 = kmax_cp[3];
      for (i = 1; i <= i3; ++i) {
	akc = (double) i;
	xk = (aka * hmatik[1] + akb * hmatik[2] + akc * hmatik[3]) * tpi;
	yk = (aka * hmatik[4] + akb * hmatik[5] + akc * hmatik[6]) * tpi;
	zk = (aka * hmatik[7] + akb * hmatik[8] + akc * hmatik[9]) * tpi;
	try = (xk * xk + yk * yk + zk * zk) * .5;
	if (try > ecut) {
	  break;
	}
      }
      
      kcmax = i - 1;
      i3 = kcmax;
      for (kc = kcmin; kc <= i3; ++kc) {
	++icount;
	aka = (double) ka;
	akb = (double) kb;
	akc = (double) kc;
	xk = (aka * hmatik[1] + akb * hmatik[2] + akc * hmatik[3]) * tpi;
	yk = (aka * hmatik[4] + akb * hmatik[5] + akc * hmatik[6]) * tpi;
	zk = (aka * hmatik[7] + akb * hmatik[8] + akc * hmatik[9]) * tpi;
        g  = sqrt(xk*xk+yk*yk+zk*zk);
        *gmin = MIN(*gmin,g);
        *gmax = MAX(*gmax,g);
	kastore[icount] = ka;
	kbstore[icount] = kb;
	kcstore[icount] = kc;
	if (kc == kcmin) {
	  ibreak1[icount] = 1;
	}
	if (kc < kcmax) {
	  ibreak2[icount] = 1;
	}
      }
    }
  }
  if(nktot!=icount){
       printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
       printf("Mismatch number of small kvectors\n");
       printf("%d vs %d\n",icount,nktot);
       printf("Contact technical support\n");
       printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
       fflush(stdout);
       exit(1);
  }

  ++icount;
  kastore[icount] = 0;
  kbstore[icount] = 0;
  kcstore[icount] = 0;

/*--------------------------------------------------------------------------*/
    }/* setkvec3d_sm */
/*==========================================================================*/




/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
/* Determines which k vectors of the full set belong to respa               */
/*==========================================================================*/

void setkvec3d_res(int kmax_res, double *hmatik, 
		  int *kastore, int *kbstore, int *kcstore, int *ibreak3, 
		  int nktot, int nktot_res)

/*==========================================================================*/
/*       Begin routine */
{/*begin routine */
/*==========================================================================*/
/* Local variables */

  int i1;
  
  double rvol23, akmax2_res;
  int ka, kb, kc;
  double xk, yk, zk;
  int icount,jcount;
  double aka, akb, akc;
  double try,voli,temp;

/*==========================================================================*/
  /* START THE ROUTINE */

  voli =  getdeth(hmatik);
  temp =  (2./3.);
  rvol23 = pow(voli,temp);
  akmax2_res = (double) ((kmax_res) * (kmax_res));
  /* SETUP THE KVECTORS */
  i1 = nktot;
  jcount = 0;
  for (icount = 1; icount <= i1; ++icount) {
    ka = kastore[icount];
    kb = kbstore[icount];
    kc = kcstore[icount];
    aka = (double) ka;
    akb = (double) kb;
    akc = (double) kc;
    xk = aka * hmatik[1] + akb * hmatik[2] + akc * hmatik[3];
    yk = aka * hmatik[4] + akb * hmatik[5] + akc * hmatik[6];
    zk = aka * hmatik[7] + akb * hmatik[8] + akc * hmatik[9];
    try = (xk * xk + yk * yk + zk * zk)/rvol23;
    ibreak3[icount] = 0;
    if (try <= akmax2_res) {
      ibreak3[icount] = 1;
      jcount +=1;
    }
  }
  if(jcount!=nktot_res){
    printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
    printf("Error setting break3 in set_recip.c\n");
    printf("%d k-vectors found versus %d expected\n",jcount,nktot_res);
    printf("Your respa kmax is %d\n",kmax_res);
    printf("*Contact technical support*\n");
    printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
    fflush(stdout);
    exit(1);
  }
/*--------------------------------------------------------------------------*/
} /* setkvec3d_res */
/*==========================================================================*/



/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
/* This routine makes sure that the k1,k2,k3 for the fft                    */
/*  satisfy the radix condition                                             */
/*==========================================================================*/

void radixme(int *kmax1, int *kmax2, int *kmax3)

/*==========================================================================*/
/* Calculate the quantity to be radicized */
/* it written with the plus one to get rid */
/* of the stupid annoying edge vectors. */
/* the factor of 4 appears because the kmax's are the */
/* maximum k vector along a direction. the normal fft */
/* grid is therefore (2*(kmax1+1))(2*(kmax2+1))(2*(kmax3+1)). */
/* the additional factor of two comes from the fact the density */
/* needs to be defined on twice as fine an fft grid */
/* (4*(kmax1+1))(4*(kmax2+1))(4*(kmax3+1)). */
/*==========================================================================*/
/*       Begin routine */
{/*begin routine */
/*==========================================================================*/
/* Local variables */

  int i1, i2;
  int krad[50], nrad, i, k1, k2, k3, iii, jjj, kkk;

/*==========================================================================*/
  k1 = (*kmax1 + 1) << 2;
  k2 = (*kmax2 + 1) << 2;
  k3 = (*kmax3 + 1) << 2;

  nrad = 50;


/* radix condition IBM_ESSL (2^h)(3^i)(5^j)(7^k)(11^m)                */
/* radix condition HP_VECLIB SGI_COMPLIB   (2^h)(3^i)(5^j)             */
#ifdef IBM_ESSL
  krad[0] = 4;
  krad[1] = 8;
  krad[2] = 12;
  krad[3] = 16;
  krad[4] = 20;
  krad[5] = 24;
  krad[6] = 32;
  krad[7] = 36;
  krad[8] = 40;
  krad[9] = 48;
  krad[10] = 60;
  krad[11] = 64;
  krad[12] = 72;
  krad[13] = 80;
  krad[14] = 96;
  krad[15] = 112;
  krad[16] = 112;
  krad[17] = 120;
  krad[18] = 128;
  krad[19] = 144;
  krad[20] = 160;
  krad[21] = 180;
  krad[22] = 192;
  krad[23] = 220;
  krad[24] = 224;
  krad[25] = 240;
  krad[26] = 252;
  krad[27] = 288;
  krad[28] = 308;
  krad[29] = 320;
  krad[30] = 336;
  krad[31] = 360;
  krad[32] = 384;
  krad[33] = 420;
  krad[34] = 440;
  krad[35] = 480;
  krad[36] = 504;
  krad[37] = 512;
  krad[38] = 560;
  krad[39] = 576;
  krad[40] = 616;
  krad[41] = 640;
  krad[42] = 660;
  krad[43] = 720;
  krad[44] = 768;
  krad[45] = 792;
  krad[46] = 880;
  krad[47] = 896;
  krad[48] = 960;
  krad[49] = 1008;
#else
  krad[0] = 4;
  krad[1] = 8;
  krad[2] = 12;
  krad[3] = 16;
  krad[4] = 20;
  krad[5] = 24;
  krad[6] = 32;
  krad[7] = 36;
  krad[8] = 40;
  krad[9] = 48;
  krad[10] = 60;
  krad[11] = 64;
  krad[12] = 72;
  krad[13] = 80;
  krad[14] = 96;
  krad[15] = 100;
  krad[16] = 108;
  krad[17] = 120;
  krad[18] = 128;
  krad[19] = 144;
  krad[20] = 160;
  krad[21] = 180;
  krad[22] = 192;
  krad[23] = 200;
  krad[24] = 216;
  krad[25] = 240;
  krad[26] = 256;
  krad[27] = 288;
  krad[28] = 300;
  krad[29] = 320;
  krad[30] = 324;
  krad[31] = 360;
  krad[32] = 384;
  krad[33] = 400;
  krad[34] = 432;
  krad[35] = 480;
  krad[36] = 500;
  krad[37] = 512;
  krad[38] = 540;
  krad[39] = 576;
  krad[40] = 600;
  krad[41] = 640;
  krad[42] = 648;
  krad[43] = 720;
  krad[44] = 768;
  krad[45] = 800;
  krad[46] = 864;
  krad[47] = 900;
  krad[48] = 960;
  krad[49] = 972;
#endif

  i1 = nrad;
  for (i = 1; i <= i1; ++i) {
    if (krad[i - 1] > k1) {
      i2 = i - 1;
      iii = MAX(i2,1);
      jjj = krad[(i - 1)];
      kkk = krad[(iii - 1)];
      if (k1 > kkk && k1 < jjj) {
	jjj = (i2 = k1 - krad[(i - 1)], abs(i2));
	kkk = (i2 = k1 - krad[(iii - 1)], abs(i2));
	if (jjj < kkk) {
	  k1 = krad[(i - 1)];
	} else {
	  k1 = krad[(iii - 1)];
	}
      }
    }
    if (krad[(i - 1)] > k2) {
      i2 = i - 1;
      iii = MAX(i2,1);
      jjj = krad[(i - 1)];
      kkk = krad[(iii - 1)];
      if (k2 > kkk && k2 < jjj) {
	jjj = (i2 = k2 - krad[(i - 1)], abs(i2));
	kkk = (i2 = k2 - krad[(iii - 1)], abs(i2));
	if (jjj < kkk) {
	  k2 = krad[(i - 1)];
	} else {
	  k2 = krad[(iii - 1)];
	}
      }
    }
    if (krad[(i - 1)] > k3) {
      i2 = i - 1;
      iii = MAX(i2,1);
      jjj = krad[(i - 1)];
      kkk = krad[(iii - 1)];
      if (k3 > kkk && k3 < jjj) {
	jjj = (i2 = k3 - krad[(i - 1)], abs(i2));
	kkk = (i2 = k3 - krad[(iii - 1)], abs(i2));
	if (jjj < kkk) {
	  k3 = krad[(i - 1)];
	} else {
	  k3 = krad[(iii - 1)];
	}
      }
    }
  }
  *kmax1 = k1 / 4 - 1;
  *kmax2 = k2 / 4 - 1;
  *kmax3 = k3 / 4 - 1;
#ifdef DEBUG
  printf("k1 %d kmax1 %d \n",k1,*kmax1);
  printf("k2 %d kmax1 %d \n",k2,*kmax2);
  printf("k3 %d kmax1 %d \n",k3,*kmax3);
#endif

/*--------------------------------------------------------------------------*/
} /* radixme */
/*==========================================================================*/




/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
/* Count number of k vectors on small grid */
/*==========================================================================*/

void countkvec3d_sm(int *nktot, double ecut, int *kmax_cp, double *hmatik )

/*==========================================================================*/
/*       Begin routine */
{/*begin routine */
/*==========================================================================*/
/* Local variables */

  int i1, i2, i3;
  int i, kbmin, kcmin, kbmax, kcmax, kamax, ka, kb, kc;
  double xk, yk, zk;
  int icount;
  double aka, akb, akc;
  double tpi, try;

/*==========================================================================*/

  tpi = M_PI * 2.;

  icount = 0;
  i1 = kmax_cp[1];
  for (i = 1; i <= i1; ++i) {
    aka = (double) i;
    xk = aka * hmatik[1] * tpi;
    yk = aka * hmatik[4] * tpi;
    zk = aka * hmatik[7] * tpi;
    try = (xk * xk + yk * yk + zk * zk) * .5;
    if (try > ecut) {
      break;
    }
  }
  kamax = i - 1;
  i1 = kamax;
  for (ka = 0; ka <= i1; ++ka) {
    aka = (double) ka;
    kbmin = -kmax_cp[2];
    if (ka == 0) {
      kbmin = 0;
    }
    for (i = kbmin; i <= 0; ++i) {
      akb = (double) i;
      xk = (aka * hmatik[1] + akb * hmatik[2]) * tpi;
      yk = (aka * hmatik[4] + akb * hmatik[5]) * tpi;
      zk = (aka * hmatik[7] + akb * hmatik[8]) * tpi;
      try = (xk * xk + yk * yk + zk * zk) * .5;
      if (try <= ecut) {
	break;
      }
    }

    kbmin = i;
    i2 = kmax_cp[2];
    for (i = 1; i <= i2; ++i) {
      akb = (double) i;
      xk = (aka * hmatik[1] + akb * hmatik[2]) * tpi;
      yk = (aka * hmatik[4] + akb * hmatik[5]) * tpi;
      zk = (aka * hmatik[7] + akb * hmatik[8]) * tpi;
      try = (xk * xk + yk * yk + zk * zk) * .5;
      if (try > ecut) {
	break;
      }
    }

    kbmax = i - 1;
    i2 = kbmax;
    for (kb = kbmin; kb <= i2; ++kb) {
      akb = (double) kb;
      kcmin = -kmax_cp[3];
      if (ka == 0 && kb == 0) {
	kcmin = 1;
      }
      for (i = kcmin; i <= 0; ++i) {
	akc = (double) i;
	xk = (aka * hmatik[1] + akb * hmatik[2] + akc * hmatik[3]) * tpi;
	yk = (aka * hmatik[4] + akb * hmatik[5] + akc * hmatik[6]) * tpi;
	zk = (aka * hmatik[7] + akb * hmatik[8] + akc * hmatik[9]) * tpi;
	try = (xk * xk + yk * yk + zk * zk) * .5;
	if (try <= ecut) {
	  break;
	}
      }

      kcmin = i;
      i3 = kmax_cp[3];
      for (i = 1; i <= i3; ++i) {
	akc = (double) i;
	xk = (aka * hmatik[1] + akb * hmatik[2] + akc * hmatik[3]) * tpi;
	yk = (aka * hmatik[4] + akb * hmatik[5] + akc * hmatik[6]) * tpi;
	zk = (aka * hmatik[7] + akb * hmatik[8] + akc * hmatik[9]) * tpi;
	try = (xk * xk + yk * yk + zk * zk) * .5;
	if (try > ecut) {
	  break;
	}
      }

      kcmax = i - 1;
      i3 = kcmax;
      for (kc = kcmin; kc <= i3; ++kc) {
	++icount;
      }
    }
  }
  *nktot = icount;

/*--------------------------------------------------------------------------*/
} /* countkvec3d_sm */
/*==========================================================================*/



/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
/* Set the PME maps   */
/*==========================================================================*/

void set_pme_wght(int nktot,int *kastore,int *kbstore,int *kcstore,
                  int nkf1,int nkf2,int nkf3,
                  int ncoef_proc,int ncoef_use,int icoef_off,
                  int pme_b_opt,double *bfact_r,double *bfact_i,
                  double *bweight_tot, int n_interp,
                  double *aj,double *rn,double *rn1)

/*==========================================================================*/
/*       Begin routine */
{/*begin routine */
/*==========================================================================*/
/* Local variables */

  size_t dim_k;
  int i,ka,kb,kc,kap,kbp,kcp;
  int ngrid_a,ngrid_b,ngrid_c;
  double tmp_a_r,tmp_a_i;
  double tmp_b_r,tmp_b_i;
  double tmp_ab_r,tmp_ab_i;
  double tmp_c_r,tmp_c_i;
  double *bden_a_r,*bden_a_i,*bweight_a;
  double *bden_b_r,*bden_b_i,*bweight_b;
  double *bden_c_r,*bden_c_i,*bweight_c;
  double *uk,*mn_k;
  int *map_a,*map_b,*map_c;
  
/*==========================================================================*/
/* I) Spherical Map */

  dim_k = (size_t) nktot;
  map_a         = (int *) cmalloc(dim_k*sizeof(int))-1;
  map_b         = (int *) cmalloc(dim_k*sizeof(int))-1;
  map_c         = (int *) cmalloc(dim_k*sizeof(int))-1;
  for(i=1;i<=nktot;i++){
    ka = kastore[i];
    kb = kbstore[i];
    kc = kcstore[i];
    if (kc < 0) {
         kcp = kc + nkf3 + 1;
    } else {
         kcp = kc + 1;
    }/*endif*/
    if (kb < 0) {
         kbp = kb + nkf2 + 1;
    } else {
         kbp = kb + 1;
    }/*endif*/
    if (ka < 0) {
         kap = ka + nkf1 + 1;
    } else {
         kap = ka + 1;
    }/*endif*/
    map_a[i]   = kap;
    map_b[i]   = kbp;
    map_c[i]   = kcp;
  }/*endfor*/

/*==========================================================================*/
/* V) Calculate bweight on the spherically cutoff grid                      */

/*--------------------------------------------------------------------------*/
/*     A) Malloc memory and define constants                                */
   ngrid_a = nkf1;
   ngrid_b = nkf2;
   ngrid_c = nkf3;

   dim_k = (size_t) (ngrid_a);
   bden_a_r  = (double *)cmalloc(dim_k*sizeof(double))-1;
   bden_a_i  = (double *)cmalloc(dim_k*sizeof(double))-1;
   bweight_a = (double *)cmalloc(dim_k*sizeof(double))-1;

   dim_k = (size_t) (ngrid_b);
   bden_b_r  = (double *)cmalloc(dim_k*sizeof(double))-1;
   bden_b_i  = (double *)cmalloc(dim_k*sizeof(double))-1;
   bweight_b = (double *)cmalloc(dim_k*sizeof(double))-1;

   dim_k = (size_t) (ngrid_c);
   bden_c_r  = (double *)cmalloc(dim_k*sizeof(double))-1;
   bden_c_i  = (double *)cmalloc(dim_k*sizeof(double))-1;
   bweight_c = (double *)cmalloc(dim_k*sizeof(double))-1;

   dim_k = (size_t) (n_interp);
   uk   = (double *)cmalloc(dim_k*sizeof(double))-1;
   mn_k = (double *)cmalloc(dim_k*sizeof(double))-1;

/*--------------------------------------------------------------------------*/
/*     B) Construct the weighting Function                                  */

   get_bspline_wght1d(n_interp,ngrid_a,aj,rn,rn1,mn_k,uk,
                      bden_a_r,bden_a_i,bweight_a);
   get_bspline_wght1d(n_interp,ngrid_b,aj,rn,rn1,mn_k,uk,
                      bden_b_r,bden_b_i,bweight_b);
   get_bspline_wght1d(n_interp,ngrid_c,aj,rn,rn1,mn_k,uk,
                      bden_c_r,bden_c_i,bweight_c);
   if(pme_b_opt > 0){
    for(i=1;i <= nktot; ++i) {
     bweight_tot[i] = bweight_a[map_a[i]]
                     *bweight_b[map_b[i]]
                     *bweight_c[map_c[i]];
    }/*endfor*/
   }/*endif*/

   if(pme_b_opt == 0 || pme_b_opt ==2){
    for(i=1;i <=ncoef_use; ++i){
      tmp_a_r = bden_a_r[map_a[i+icoef_off]];
      tmp_a_i = bden_a_i[map_a[i+icoef_off]];
      tmp_b_r = bden_b_r[map_b[i+icoef_off]];
      tmp_b_i = bden_b_i[map_b[i+icoef_off]];
      tmp_c_r = bden_c_r[map_c[i+icoef_off]];
      tmp_c_i = bden_c_i[map_c[i+icoef_off]];
      tmp_ab_r = tmp_a_r*tmp_b_r - tmp_a_i*tmp_b_i;
      tmp_ab_i = tmp_a_i*tmp_b_r + tmp_a_r*tmp_b_i;
      bfact_r[i] = tmp_ab_r*tmp_c_r - tmp_ab_i*tmp_c_i;
      bfact_i[i] = tmp_ab_i*tmp_c_r + tmp_ab_r*tmp_c_i;
    }/*endfor*/
    if(ncoef_proc > ncoef_use){
      bfact_r[(ncoef_proc)] = bden_a_r[1]*bden_b_r[1]*bden_c_r[1];
      bfact_i[(ncoef_proc)] = 0.0;
    }
  }/*endif*/

/*========================================================================*/
/* Free memory */

  cfree(&map_a[1]);
  cfree(&map_b[1]);
  cfree(&map_c[1]);
  cfree(&bden_a_r[1]); 
  cfree(&bden_a_i[1]); 
  cfree(&bweight_a[1]);
  cfree(&bden_b_r[1]); 
  cfree(&bden_b_i[1]); 
  cfree(&bweight_b[1]);
  cfree(&bden_c_r[1]); 
  cfree(&bden_c_i[1]); 
  cfree(&bweight_c[1]);
  cfree(&uk[1]); 
  cfree(&mn_k[1]); 
/*--------------------------------------------------------------------------*/
} /* set_pme_map */
/*==========================================================================*/





/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
/* Calculate the B spline weighting function   */
/*==========================================================================*/

void get_bspline_wght1d(int n_interp,int ngrid,double *aj,double *rn,
                        double *rn1,double *mn_k,double *uk,
                        double *bden_r,double *bden_i,double *bweight)

/*==========================================================================*/
/*       Begin routine */
{/*begin routine */
/*==========================================================================*/
/* Local variables */

  int k,k1,n,m,j;
  int nk,iopt,ierr,incl,nn,incn;
  double arg;
  double tpi_n,mn_k_tmp;
  double bnum_real,bnum_imag;
  double bden_real,bden_imag,denom;
  double tmp_real,tmp_imag;
  double grid;

/*==========================================================================*/
/* I) Get B spline coefficients                                         */

   grid  = (double) ngrid;
   for(j=1;j<=n_interp;j++){
     aj[j] = (double) (j-1);
     rn[j] = (double) (j);
     if(j > 1){rn1[j] = 1.0/((double)(j-1));}
   }/*endfor*/
   rn1[1] = 0.0;
   mn_k[1] = 1.0; 
   uk[1]   = 1.0;
   for(k=2;k<=n_interp;k++){
     uk[k] = (double) (k);
     mn_k[k] = 0.0; 
   }/*endfor*/
   for(n=3;n<=n_interp;n++){
     for(k=n;k>=2;k--){
       k1 = k-1;
       mn_k_tmp  = (uk[k]*mn_k[k]+(rn[n]-uk[k])*mn_k[k1])*rn1[n];
       mn_k[k] = mn_k_tmp;
     }/*endfor*/
     mn_k[1] = uk[1]*mn_k[1]*rn1[n];
   }/*endfor*/

/*==========================================================================*/
/* II) Transform coeffs                                                     */

/*----------------------------------------------*/
/*       i)On HP perform transform using 1D FFT */
#ifdef HP_VECLIB
   for(k=1;k<=ngrid;k++){
     bden_r[k] = 0.0;
     bden_i[k] = 0.0;
   }/*endfor*/
   for(k=1;k<=n_interp-1;k++){
     bden_r[k] = mn_k[k];
   }/*endfor*/
    iopt = 1;ierr = 0;incl = 1;nn = 1;incn = 1;
    nk = ngrid;
    DFFTS(&bden_r[1],&bden_i[1],&nk,&incl,&nn,&incn,&iopt,&ierr);
   for(k=1;k<=ngrid;k++){
     bden_i[k] *= -1.0;
   }/*endfor*/
/*-----------------------------------------------*/
/*     ii) Not HP perform transform using 1D SFT */
#else
   tpi_n = 2.0*M_PI/grid;
   for(m=1;m<=ngrid;m++){
     bden_r[m] = 0.0;
     bden_i[m] = 0.0;
     for(k=1;k<=n_interp-1;k++){
       arg = tpi_n*((double)((k-1)*(m-1)));
       bden_r[m] += cos(arg)*mn_k[k];
       bden_i[m] += sin(arg)*mn_k[k];
     }/*endfor*/
   }/*endfor*/
#endif

/*==========================================================================*/
/* III) Make separable bweights                                             */

   tpi_n = 2.0*M_PI*((double) (n_interp-1))/grid;
   for(m=1;m<=ngrid;m++){

     arg = tpi_n*((double)(m-1));
     bnum_real  = cos(arg);
     bnum_imag  = sin(arg);
     bden_real  = bden_r[m];
     bden_imag  = bden_i[m];

     denom      = bden_real*bden_real + bden_imag*bden_imag;
     tmp_real   = (bnum_real*bden_real + bnum_imag*bden_imag)/denom;
     tmp_imag   = (bnum_imag*bden_real - bnum_real*bden_imag)/denom;
     bweight[m] =  tmp_real*tmp_real + tmp_imag*tmp_imag;

     bden_r[m] = tmp_real;
     bden_i[m] = tmp_imag;

   }/*endfor*/

/*--------------------------------------------------------------------------*/
     }/*end routine */
/*==========================================================================*/












