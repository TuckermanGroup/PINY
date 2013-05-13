/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
/*                                                                          */
/*                         PI_MD:                                           */
/*             The future of simulation technology                          */
/*             ------------------------------------                         */
/*                     Module: weight_node_gauss_hermite                    */
/*                                                                          */
/* Sets up Gauss-Hermite integration stuff                                  */ 
/*                                                                          */
/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/


#include "standard_include.h"
#include "../typ_defs/typedefs_gen.h"
#include "../typ_defs/typedefs_par.h"
#include "../proto_defs/proto_vps_params_local.h"
#include "../proto_defs/proto_friend_lib_entry.h"


/*==================================================================*/
void weight_node_gauss_hermite(int ngh,double *rgh,double *wgh)
/*==================================================================*/
{/*begin routine*/
/*==================================================================*/
 int nmax = 180;
 int index,i,m,n,nup,np,niter;
 int nmall;

 double *a,*ap,*da;
 double *root0,*root;
 double fact,r,rapp;
 double pi,sqpi;
 double w,x,y,y0,yp,yp1,dy,zero;
   
/*--------------------------------------------*/  
/* Local mallocs                              */

   nmall = nmax + 1;

 /* start at zero */

   a = (double *) calloc(nmall,sizeof(double));
  da = (double *) calloc(nmall,sizeof(double));
  ap = (double *) calloc(nmall,sizeof(double));

 /* start at one */

  root0 = (double *) calloc(nmall,sizeof(double))-1;
  root  = (double *) calloc(nmall,sizeof(double))-1;

/*----------------------------------------------*/  
/* double number of gauss-hermite points        */
/*   note: roots are symmetric about the orgin  */
    n=2*ngh;

/* Compute the coefficients wgh  */
   if(n>nmax){
     fprintf(stderr,"Too many Gauss-Hermite integration points \n");
     exit(1); 
   }

   nup = ngh;

/* Zero the arrays */

   for(i=0; i<= n; i++){
     a[i]=0.0000000000000000000000000000000;
    da[i]=0.0000000000000000000000000000000;
    ap[i]=0.0000000000000000000000000000000;
   }/*endfor*/


   fact=(double)(nup+1);

  for(i=(nup+2); i<= n; i++){
     fact *= (double)i;
  }/*endfor*/

  a[0]=fact*pow((-1.0000000000000000000000),nup);
  fact=1.00000000000000000000000000000000;
     r=0.00000000000000000000000000;

  for(i=1; i<=n; i++){
   r += 1.00000000000000000000000;
   fact *= sqrt(r);
  }/*endfor*/

   pi= M_PI;
  sqpi=sqrt(pi);
  sqpi=sqrt(sqpi);

  a[0]= a[0]/(pow(2.000000000000,nup))/fact/sqpi;

  for(i=(nup-1); i>= 0;i--){
   m=i+1;
   index=n-2*i;
   a[index] = -a[index-2]*4.0000000000000000000000000
            *  (double)(m)/(double)(n-m-m+2)/(double)(n-m-m+1);
  }/*endfor*/

/*Derivative da */
 for(i=0; i<=(n-1);i++){
  da[i]=((double)(i) + 1.0000000000000000000000)*a[i+1];
 }/*endfor*/

    np=n+1;
    nup=np/2;
    fact=(double)(nup+1);

  for(i=(nup+2);i<=np; i++){
     fact *= (double)i;
  }/*endfor*/

  ap[1]= 2.000000000000000000000*(pow(-1.000000000000000000000,nup))*fact;
  fact = 1.0000000000000000000000;
     r = 0.00000000000000000000000000;

  for(i=1; i<=np; i++){
    r += 1.00000000000000000000000;
    fact *= sqrt(r);
  }/*endfor*/

  ap[1]=ap[1]/(pow(2.000000000000,nup))/fact/sqpi/sqrt(2.000000000000);

  for(i=(nup-1);i>=0; i--){
    m=i+1;
    index=np-2*i;
    ap[index] = -ap[index-2]*4.000000000000000000
              * (double)(m)/(double)(np-m-m+2)/(double)(np-m-m+1);
  }/*endfor*/


/* First approximation to roots */
   zero = 0.000000000000000000;
   newt_eval(a,n,&zero,&y0); 


   index = 0;

   for(i=1; i<= 5000; i++){
     x = 0.00333*(double)i;
     newt_eval(a,n,&x,&y);
      if(y*y0 < 0.0){
     index++;
     root0[index] = x;
     y0=y;
    }/*endif*/
   }/*endfor*/

/* Newton search for the roots */
   for(i=1; i<=index; i++){
    niter=0;
    newt_eval(a,n,&(root0[i]),&y);

    while(fabs(y) > 1.e-16){
      newt_eval(a,n,&(root0[i]),&y);
      niter++;
      if(niter > 20)break;

      newt_eval(da,n-1,&(root0[i]),&yp);

      root0[i]=root0[i]-y/yp;

      newt_eval(a,n,&(root0[i]),&y);

    }/*endif*/

     root[i] = root0[i];
     rgh[i]  = root[i];

   }/*endfor*/


/* Weights */

   rapp = -ap[n+1]/a[n];

   for(i=1; i<=index; i++){
     newt_eval(ap,n+1,&(root[i]),&yp1);
     newt_eval(da,n-1,&(root[i]),&dy);

     w=rapp/yp1/dy;
     wgh[i]=w;
   }/*endfor*/ 


/*===========================================*/
/*Free locally assigned memory */

   free(a);
   free(da);
   free(ap);

   free(&root0[1]);
   free(&root[1]);

/*==================================================================*/
}/*end routine*/
/*==================================================================*/

/*==================================================================*/
void newt_eval(double *coef,int n,double *px,double *py) 
{
  int i;
  double x,y;

  x = *px;
  y = coef[n];

  for(i=(n-1); i>=0; i--){
   y = coef[i]+x*y;
  }/*endfor*/

  *px = x;
  *py = y;

}
/*==================================================================*/


/*==================================================================*/
/*==================================================================*/

void limit_gauss_hermite(int *pngh,double *rgh,double *wgh)

/*==================================================================*/
{/*begin routine*/

  int i,iii;
  int ngh_count;
  int ngh = *pngh;

    ngh_count = 0;

    for(i=1; i<=ngh;i++){
      if(wgh[i] > 1.0e-10){
        ngh_count++;
        wgh[ngh_count] = wgh[i];
        rgh[ngh_count] = rgh[i];
      }/*endif*/
    }/*endfor*/

    *pngh = ngh_count;

}/*end routine*/
/*==================================================================*/
