/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
/*                                                                          */
/* This subprogram calculates the intermediate scattering functions         */
/* - assuming orthogonal primitives vectors e.g. an isotropic system        */
/* - assuming a 3D system                                                   */
/*                                                                          */
/*               I_c(Q,t)= <rho(-k,0)* rho(k,t)>                            */
/*                rho(k,t)=sum_j (exp[+iQR_j(t)])                           */
/*                                                                          */
/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

/*==========================  INCLUDE FILES ================================*/
#include <stddef.h>
#include "standard_include.h"
#include "../typ_defs/typedefs_stat.h"
#include "../typ_defs/typedefs_class.h"
#include "../typ_defs/typedefs_gen.h"
#include "../typ_defs/typedefs_bnd.h"
#include "../proto_defs/proto_analysis_local_entry.h"
#include "../proto_defs/proto_friend_lib_entry.h"
/*==========================================================================*/

/*==========================  MISCELLANEOUS ================================*/
#define DEBUG_OFF
#define DEBUG_POS_ZERO_OFF
#define DEBUG_KVEC_OFF
#define nptk_max (50)     /* Maximun number ok k-vectors                    */
#define kmax_max (30.0)   /* Maximun value that can take kmax in Angstrom-1 */
/*==========================================================================*/

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void calcul_ickt_iso_md(CLASS *class,GENERAL_DATA *general_data,
                     ANALYSIS *analysis)
/*==========================================================================*/
/*  Begin subprogram                                                        */
/*==========================================================================*/
{ 
 /*=========================================================================*/

 /*-------------------------------------------------------------------------*/
 /* Local variable declarations                                             */
 /*-------------------------------------------------------------------------*/

  unsigned int nstep,nrun;
  unsigned int njump_ickt_iso;
  unsigned int nbpas_ickt_iso;
  unsigned int nnn_ickt_iso;

 /*-------------------------------------------------------------------------*/
 /* Begin                                                                   */
 /*-------------------------------------------------------------------------*/

/*==========================================================================*/
} /* end calcul_ickt_iso */
/*==========================================================================*/

