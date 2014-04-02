/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
/*                                                                          */
/*                         PI_MD:                                           */
/*             The future of simulation technology                          */
/*             ------------------------------------                         */
/*                     Module: exl_sort                                     */
/*                                                                          */
/* This subprogram sorts the exclusion list eliminating repeaters           */
/*                                                                          */
/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

#include "standard_include.h"
#include "../typ_defs/typedefs_class.h"
#include "../typ_defs/typedefs_gen.h"
#include "../typ_defs/typedefs_bnd.h"
#include "../typ_defs/typedefs_par.h"
#include "../proto_defs/proto_lists_local.h"
#include "../proto_defs/proto_friend_lib_entry.h"

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void exl_sort(int *n, int index[],int jndex[],int ind_max)

/*=======================================================================*/
/*            Begin subprogram:                                          */
{/*begin routine*/
/*=======================================================================*/
/*          Local variable declarations                                  */

  int m,ir,i,j,rindex,rjndex,iii;
  int k,*kndex,*mask,isub,temp;

/*=======================================================================*/
/* I) Setup                        */

  m  = *n/2+1;
  ir = *n;

/*=======================================================================*/
/* II) Sort array index keeping jndex commensurrate */

  for(;;){

/*---------------------------------------------------------------------*/
/*  A)hire rindex */
    if(m>1){ 
      m--;
      rindex = index[m];
      rjndex = jndex[m];
/*--------------------------------------------------------------------*/
/*  B)retire/promote index[1] */
    }else{
      rindex = index[ir];
      rjndex = jndex[ir];
      index[ir]=index[1];
      jndex[ir]=jndex[1];
      ir--;
      if(ir==1){
	index[1]=rindex;
	jndex[1]=rjndex;
	break;
      }/*endif*/
    }/*endif*/
/*---------------------------------------------------------------------*/
/*  C)put rindex in appropriate slot */
    i=m;
    j=2*m;
    while(j<=ir){
      /*    a)compare to rindex to underling */
      if((j<ir) && (index[j]< index[(j+1)])) j++;
      /*    b)demote */
      if(rindex<index[j]){
	index[i]=index[j];
	jndex[i]=jndex[j];
	i=j;
	j=2*j;
      }else{
	/*    c)if no demotations exit while */
	j=ir+1;
      }/*endif*/
    } /*endwhile*/
    /*    d)slot rindex */
    index[i] = rindex;
    jndex[i] = rjndex;
  }/*endfor*/

/*========================================================================*/
/* III) Eliminate repeaters          */

  kndex    = (int *)cmalloc((ind_max)*sizeof(int))-1;
  mask     = (int *)cmalloc((*n)*sizeof(int))-1;

/*----------------------------------------------------------------------*/
/*   A)  Define sub blocks based on value of index */

  for(i=1;i<=ind_max;i++) kndex[i]=0;
  for(i=1;i<=*n;i++) mask[i]=1;
  for(i=1;i<=*n;i++) {
    kndex[(index[i])]++;
  }

/*----------------------------------------------------------------------*/
/*   B)  Find repeaters in each sub block */
  isub = 0;
  temp=0;
  for(i=1;i<=ind_max;i++){
    for(j=isub+1;j<=isub+kndex[i];j++){
      for(k=j+1;k<=isub+kndex[i];k++){
	if(jndex[k]==jndex[j]){mask[k]=0;}
      }/*endfor*/
    }/*endfor*/
    isub = isub + kndex[i];
  } /*endfor*/
/*-----------------------------------------------------------------------*/
/*   C)  Eliminate repeaters in each sub block */
  isub = 1;
  for(i=1;i<=*n;i++){
    index[isub]=index[i];
    jndex[isub]=jndex[i];
    if(mask[i]==1){isub += 1;}
  }/*endfor*/
  *n = isub-1;
  cfree(&kndex[1]);
  cfree(&mask[1]);

/*-----------------------------------------------------------------------*/
} /*end routine*/ 
/*==========================================================================*/
