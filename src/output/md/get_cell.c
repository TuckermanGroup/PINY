/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
/*                                                                          */
/*                         PI_MD:                                           */
/*             The future of simulation technology                          */
/*             ------------------------------------                         */
/*                     Module: get_cell.c                                   */
/*                                                                          */
/* This subprogram gets cell lenths and angles                              */
/*                                                                          */
/*==========================================================================*/


#include "standard_include.h"
#include "../typ_defs/typedefs_bnd.h"
#include "../typ_defs/typedefs_gen.h"
#include "../typ_defs/typedefs_class.h"
#include "../proto_defs/proto_output_local.h"


/*==========================================================================*/
/*               Header:                                                    */
void get_cell(double *hmat,double *a,  double *b,  double *c,
                         double *tab,double *tbc,double *tac)
/*==========================================================================*/
{/*begin routine*/
/*==========================================================================*/

  double tpi;

/*==========================================================================*/

  tpi = 2.0*M_PI;
  (*a) = sqrt(hmat[1]*hmat[1]+hmat[2]*hmat[2]+hmat[3]*hmat[3]);
  (*b) = sqrt(hmat[1]*hmat[1]+hmat[4]*hmat[4]+hmat[7]*hmat[7]);
/*  if((*a)<(*b)){ */
    (*a)    = sqrt(hmat[1]*hmat[1]+hmat[2]*hmat[2]+hmat[3]*hmat[3]);
    (*b)    = sqrt(hmat[4]*hmat[4]+hmat[5]*hmat[5]+hmat[6]*hmat[6]);
    (*c)    = sqrt(hmat[7]*hmat[7]+hmat[8]*hmat[8]+hmat[9]*hmat[9]);
    (*tab)  = (hmat[1]*hmat[4]+hmat[2]*hmat[5]+hmat[3]*hmat[6]);
    (*tab)  = acos((*tab)/((*a)*(*b)))*360.0/tpi;
    (*tac)  = (hmat[1]*hmat[7]+hmat[2]*hmat[8]+hmat[3]*hmat[9]);
    (*tac)  = acos((*tac)/((*a)*(*c)))*360.0/tpi;
    (*tbc)  = (hmat[4]*hmat[7]+hmat[5]*hmat[8]+hmat[6]*hmat[9]);
    (*tbc)  = acos((*tbc)/((*b)*(*c)))*360.0/tpi;
/* Kari's thing which I don't like */
/*  }else{
    (*a)    = sqrt(hmat[1]*hmat[1]+hmat[4]*hmat[4]+hmat[7]*hmat[7]);
    (*b)    = sqrt(hmat[2]*hmat[2]+hmat[5]*hmat[5]+hmat[8]*hmat[8]);
    (*c)    = sqrt(hmat[3]*hmat[3]+hmat[6]*hmat[6]+hmat[9]*hmat[9]);
    (*tab)  = (hmat[1]*hmat[2]+hmat[4]*hmat[5]+hmat[7]*hmat[8]);
    (*tab)  = acos((*tab)/((*a)*(*b)))*360.0/tpi;
    (*tac)  = (hmat[1]*hmat[3]+hmat[4]*hmat[6]+hmat[7]*hmat[9]);
    (*tac)  = acos((*tac)/((*a)*(*c)))*360.0/tpi;
    (*tbc)  = (hmat[2]*hmat[3]+hmat[5]*hmat[6]+hmat[8]*hmat[9]);
    (*tbc)  = acos((*tbc)/((*b)*(*c)))*360.0/tpi;
  } */

/*--------------------------------------------------------------------------*/
} /*end routine*/ 
/*==========================================================================*/





/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
void getdeth_avg(double *hmat,double *vol)
/*==========================================================================*/
{/*begin routine*/
/*==========================================================================*/

  (*vol) = (hmat[1]*(hmat[5]*hmat[9]-hmat[6]*hmat[8]) +
	    hmat[2]*(hmat[6]*hmat[7]-hmat[4]*hmat[9]) +
	    hmat[3]*(hmat[4]*hmat[8]-hmat[5]*hmat[7]));

/*--------------------------------------------------------------------------*/
} /*end routine*/ 
/*==========================================================================*/





