/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
/*                                                                          */
/*                         PI_MD:                                           */
/*             The future of simulation technology                          */
/*             ------------------------------------                         */
/*                     Module: get_clong.c                                  */
/*                                                                          */
/* Get the long range correction                                            */
/*                                                                          */
/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
/*               Includes:                                                  */

#include "standard_include.h"
#include "../typ_defs/typedefs_par.h"
#include "../typ_defs/typedefs_class.h"
#include "../proto_defs/proto_inter_params_local.h"
#include "../proto_defs/proto_friend_lib_entry.h"

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void get_clong(int natm_tot,int natm_typ,int ninter,
	       int iatm_typ[],double c6m[],
               int inter_map_index[],
	       double *clong,double *clong_res,double cutoff[],
	       double cutoff_res[],int iswit, double rheal)

/*=======================================================================*/
    {/*begin routine*/
/*=======================================================================*/
/*            Local variable declarations                                */

  int iii,i,j,ntype1;
  int inter_map_index_i;
  int *index,*num_typ;
  double fpi3,rmin,fpi,fact,fact_res;
  double *clongm;

/*======================================================================*/
/* 0) Allocate temporary arrays                                         */

  num_typ = (int *) cmalloc(natm_typ*sizeof(int))-1;
  index   = (int *) cmalloc(natm_typ*sizeof(int))-1;
  clongm  = (double *) cmalloc(ninter*sizeof(double))-1;

/*=======================================================================*/
/* I) Determine how many of each atom type there are                     */
  
  for(i=1;i<=natm_typ;i++) num_typ[i] = 0;  
  for(i=1;i<=natm_tot;i++) {num_typ[(iatm_typ[i])] ++;}
  
/*=======================================================================*/
/* II) Get the long range correction                                     */
  
  ntype1=natm_typ+1;
  for(i=1;i<=natm_typ;i++){
    for(j=1;j<=i;j++){
      index[j] = ((j-1)*(2*ntype1-j))/2 + i - j + 1;
    }
    for(j=1;j<=i-1;j++){
      clongm[(index[j])] = num_typ[i]*num_typ[j];
    }
    clongm[(index[i])]= (num_typ[i]*(num_typ[i]-1))/2;
/*    clongm[(index[i])]= (num_typ[i]*(num_typ[i]))/2; */
  } 
  
  *clong=0.0;
  *clong_res= 0.0;
  fpi3=4.0*M_PI/3.0;
  fpi=4.0*M_PI;
  if(iswit==0){
   for(i=1;i<=ninter;i++){
     inter_map_index_i = inter_map_index[i];
     rmin = (cutoff[inter_map_index_i] < cutoff_res[inter_map_index_i] 
           ? cutoff[inter_map_index_i]:cutoff_res[inter_map_index_i]);
     *clong -= fpi3*c6m[inter_map_index_i]*clongm[i]/pow(cutoff[inter_map_index_i],3.0);
     *clong_res -= fpi*c6m[inter_map_index_i]*clongm[i]*(
                   (1.0/(rheal*rheal))*(1.0/rmin+1.0/(rmin-rheal))
                  +(2.0/(rheal*rheal*rheal))*log((rmin-rheal)/rmin));
   }/*endfor*/
  }else{
   for(i=1;i<=ninter;i++){
     inter_map_index_i = inter_map_index[i];
     if((c6m[inter_map_index_i]>0)&&(cutoff[inter_map_index_i]-rheal>0)){
      rmin        = (cutoff[inter_map_index_i] < cutoff_res[inter_map_index_i] 
                  ?  cutoff[inter_map_index_i] : cutoff_res[inter_map_index_i]);
      fact        = (1.0/(rheal*rheal))*(1.0/cutoff[inter_map_index_i]+1.0
                  / (cutoff[inter_map_index_i]-rheal))
                  + (2.0/(rheal*rheal*rheal))*log((cutoff[inter_map_index_i]-rheal)
                  /  cutoff[inter_map_index_i]);
      *clong     -= (fpi*c6m[inter_map_index_i]*clongm[i]*fact);
      fact_res    = (1.0/(rheal*rheal))*(1.0/rmin+1.0/(rmin-rheal))
                   +(2.0/(rheal*rheal*rheal))*log((rmin-rheal)/rmin);
      *clong_res -= fpi*c6m[inter_map_index_i]*clongm[i]*fact_res;
     }/*endif*/
   }/*endfor*/
  }/*endif*/

/*======================================================================*/
/* III) Free temporaries                                                */
  
  cfree(&num_typ[1]);
  cfree(&index[1]);
  cfree(&clongm[1]);

/*----------------------------------------------------------------------*/
  }/*end routine */
/*======================================================================*/

