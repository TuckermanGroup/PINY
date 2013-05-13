/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
/*                                                                          */
/*                         PI_MD:                                           */
/*             The future of simulation technology                          */
/*             ------------------------------------                         */
/*                   Module: group_bond_con_util.c                          */
/*                                                                          */
/* This file contains utilities used by the atom based group constraint     */
/* routines                                                                 */
/*                                                                          */
/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

#include "standard_include.h"
#include "../typ_defs/typedefs_class.h"
#include "../typ_defs/typedefs_gen.h"
#include "../typ_defs/typedefs_bnd.h"
#include "../proto_defs/proto_intra_con_local.h"
#include "../proto_defs/proto_friend_lib_entry.h"
#include "../proto_defs/proto_math.h"

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
void nrerror(char error_text[])
/* Numerical Recipes standard error handler */
{/*begin routine*/
   fprintf(stderr,"Numerical Recipes run-time error...\n");
   fprintf(stderr,"%s\n",error_text);
   fprintf(stderr,"...now exiting to system...\n");
   exit(1);
}/* end routine */
/*==========================================================================*/

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
double *dvector(int nl, int nh)
/* allocate a double vector with subscript range v[nl..nh] */
{/*begin routine*/
   double *v;

   v=(double *)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(double)));
   if (!v) nrerror("allocation failure in dvector()");
   return v-nl+NR_END;
}/* end routine */
/*==========================================================================*/

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
double **dmatrix(int nrl, int nrh, int ncl, int nch)
/* allocate a double matrix with subscript range m[nrl..nrh][ncl..nch] */
{/*begin routine*/
   int i, nrow=nrh-nrl+1,ncol=nch-ncl+1;
   double **m;

   /* allocate pointers to rows */
   m=(double **) malloc((size_t)((nrow+NR_END)*sizeof(double*)));
   if (!m) nrerror("allocation failure 1 in dmatrix()");
   m += NR_END;
   m -= nrl;

   /* allocate rows and set pointers to them */
   m[nrl]=(double *) malloc((size_t)((nrow*ncol+NR_END)*sizeof(double)));
   if (!m[nrl]) nrerror("allocation failure 2 in dmatrix()");
   m[nrl] += NR_END;
   m[nrl] -= ncl;

   for(i=nrl+1;i<=nrh;i++) m[i]=m[i-1]+ncol;

   /* return pointer to array of pointers to rows */
   return m;
}/* end routine */
/*==========================================================================*/

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
double ***d3tensor(int nrl, int nrh, int ncl, int nch,
                          int ndl, int ndh)
/* allocate a double 3tensor with range t[nrl..nrh][ncl..nch][ndl..ndh] */
{/*begin routine*/
   int i,j,nrow=nrh-nrl+1,ncol=nch-ncl+1,ndep=ndh-ndl+1;
   double ***t;

   /* allocate pointers to pointers to rows */
   t=(double ***) malloc((size_t)((nrow+NR_END)*sizeof(double**)));
   if (!t) nrerror("allocation failure 1 in d3tensor()");
   t += NR_END;
   t -= nrl;

   /* allocate pointers to rows and set pointers to them */
   t[nrl]=(double **) malloc((size_t)((nrow*ncol+NR_END)*sizeof(double*)));
   if (!t[nrl]) nrerror("allocation failure 2 in d3tensor()");
   t[nrl] += NR_END;
   t[nrl] -= ncl;

   /* allocate rows and set pointers to them */
   t[nrl][ncl]=(double *)
   malloc((size_t)((nrow*ncol*ndep+NR_END)*sizeof(double)));
   if (!t[nrl][ncl]) nrerror("allocation failure 3 in d3tensor()");
   t[nrl][ncl] += NR_END;
   t[nrl][ncl] -= ndl;

   for(j=ncl+1;j<=nch;j++) t[nrl][j]=t[nrl][j-1]+ndep;
   for(i=nrl+1;i<=nrh;i++) {
      t[i]=t[i-1]+ncol;
      t[i][ncl]=t[i-1][ncl]+ncol*ndep;
      for(j=ncl+1;j<=nch;j++) t[i][j]=t[i][j-1]+ndep;
   }

   /* return pointer to array of pointers to rows */
   return t;
}/* end routine */
/*==========================================================================*/

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
void free_dvector(double *v, int nl, int nh)
/* free a double vector allocated with dvector() */
{/*begin routine*/
   free((v+nl-NR_END));
}/* end routine */
/*==========================================================================*/

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
void free_dmatrix(double **m, int nrl, int nrh, int ncl, int nch)
/* free a double matrix allocated by dmatrix() */
{/*begin routine*/
   free((char *)(m[nrl]+ncl-NR_END));
   free((char *)(m+nrl-NR_END));
}/* end routine */
/*==========================================================================*/

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
void free_d3tensor(double ***t, int nrl, int nrh, int ncl, int nch,
   int ndl, int ndh)
/* free a double d3tensor allocated by d3tensor() */
{/*begin routine*/
   free((char *)(t[nrl][ncl]+ndl-NR_END));
   free((char *)(t[nrl]+ncl-NR_END));
   free((char *)(t+nrl-NR_END));
}/* end routine */
/*==========================================================================*/

