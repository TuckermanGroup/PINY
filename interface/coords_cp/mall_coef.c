/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
/*                                                                          */
/*                         PI_MD:                                           */
/*             The future of simulation technology                          */
/*             ------------------------------------                         */
/*                     Module: mall_coef                                    */
/*                                                                          */
/* This subprogram reads atm-atm_NHC vol-vol_NHC input for a MD on a        */ 
/* LD-classical potential energy surface (LD-PES)                           */
/*                                                                          */
/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/


#include "standard_include.h"
#include "../typ_defs/typedefs_class.h"
#include "../typ_defs/typedefs_par.h"
#include "../typ_defs/typedefs_gen.h"
#include "../typ_defs/typedefs_cp.h"
#include "../proto_defs/proto_coords_cp_entry.h"
#include "../proto_defs/proto_coords_cp_local.h"
#include "../proto_defs/proto_friend_lib_entry.h"

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void mall_coef(CP *cp,SIMOPTS *simopts,int pi_beads_proc)

/*======================================================================*/
/*                Begin Routine */
{   /*begin routine */
/*======================================================================*/
/*               Local variable declarations                            */

  int cp_on,iii;

/*========================================================================*/
/* I) Malloc the vector structures                                        */


   cp->cpcoeffs_pos = (CPCOEFFS_POS *)cmalloc(pi_beads_proc *
                        sizeof(CPCOEFFS_POS))-1;
   cp->cptherm_pos = (CPTHERM_POS *)cmalloc(pi_beads_proc*
                                   sizeof(CPTHERM_POS))-1;


/*========================================================================*/
} /* end routine */
/*==========================================================================*/




