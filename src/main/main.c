/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
/*                                                                          */
/*                         PI_MD:                                           */
/*             The future of simulation technology                          */
/*             ------------------------------------                         */
/*                      Module: main                                        */
/*                                                                          */
/* This program performs MD on a classical potential energy surface (PES),  */
/* minimization on a classical PES,                                         */
/* MD on a mixed classical-density functional PES,                          */
/* minimization on a mixed classical-density functional PES,                */
/* PIMD on a classical PES,                                                 */
/* centroid minimization on a classical PES,                                */
/* PIMD on a mixed classical-density functional PES,                        */
/* centroid minimization on a mixed classical-density functional PES.       */
/*                                                                          */
/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

#include "standard_include.h"
#include "../typ_defs/typedefs_class.h"
#include "../typ_defs/typedefs_bnd.h"
#include "../typ_defs/typedefs_gen.h"
#include "../typ_defs/typedefs_cp.h"
#include "../typ_defs/typedefs_par.h"
#include "../typ_defs/typedefs_stat.h"
#include "../proto_defs/proto_main_entry.h"
#include "../proto_defs/proto_main_local.h"
#include "../proto_defs/proto_main_cp_local.h"
#include "../proto_defs/proto_parse_entry.h"
#include "../proto_defs/proto_communicate_wrappers.h"
#include "../proto_defs/proto_communicate_entry.h"
#if defined PLUMED
#include "../proto_defs/proto_plumed.h"
#endif


/*==========================================================================*/
void Init_PINY(int argc, char *argv[],
               CLASS* class, GENERAL_DATA* general_data) {
/*==========================================================================*/

  Init(&argc, &argv, &((*class).communicate.world));
  Comm_size((*class).communicate.world, &((*class).communicate.np));
  Comm_rank((*class).communicate.world, &((*class).communicate.myid));
  (*general_data).error_check_on = ((*class).communicate.myid==0?1:0);

}
/*--------------------------------------------------------------------------*/


/*==========================================================================*/
int main (int argc, char *argv[]) {
/*==========================================================================*/

  /* local variables */
  CLASS class;
  BONDED bonded;
  GENERAL_DATA general_data;
  CP cp;
  ANALYSIS analysis;

  /*======================*/
  /* Check for input file */

  if(argc < 2) {
    printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
    printf("No input file specified.\n");
    printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
    fflush(stdout);
    exit(1);
  }

  /*================*/
  /* initialize MPI */

  Init_PINY(argc, argv, &class, &general_data);

  /*=======================*/
  /* invoke user interface */

  parse(&class, &bonded, &general_data, &cp, &analysis, argv[1]);

  #if defined PLUMED
  plumed_piny_init(&general_data, &class);
  #endif

  /*====================*/
  /* perform simulation */

  /*-----------------------*/
  /* run MD/PIMD/CP/CPPIMD */

  if (general_data.simopts.md ==1 ) {
    control_md(&class, &bonded, &general_data, &analysis);
  }
  if (general_data.simopts.pimd == 1) {
    control_pimd(&class, &bonded, &general_data, &analysis);
  }
  if (general_data.simopts.cp_wave + general_data.simopts.cp == 1) {
    control_cp(&class, &bonded, &general_data, &cp, &analysis);
  }
  if (general_data.simopts.cp_pimd + general_data.simopts.cp_wave_pimd == 1) {
      control_cp_pimd(&class, &bonded, &general_data, &cp, &analysis);
  }

  /*-------------------------*/
  /* debug MD/PIMD/CP/CPPIMD */

  if (general_data.simopts.debug == 1) {
    control_debug(&class, &bonded, &general_data);
  }
  if (general_data.simopts.debug_pimd == 1) {
    control_debug_pimd(&class, &bonded, &general_data);
  }
  if (general_data.simopts.debug_cp == 1) {
    control_debug_cp(&class, &bonded, &general_data, &cp);
  }
  if (general_data.simopts.debug_cp_pimd == 1) {
    control_debug_cp_pimd(&class, &bonded, &general_data, &cp);
  }

  /*----------------------------*/
  /* minimize MD/PIMD/CP/CPPIMD */

  if (general_data.simopts.minimize == 1) {
    control_min(&class, &bonded, &general_data, &analysis);
  }
  if (general_data.simopts.cp_wave_min + general_data.simopts.cp_min == 1) {
    control_cp_min(&class,&bonded,&general_data,&cp,&analysis);
  }
  if (general_data.simopts.cp_wave_min_pimd == 1 ) {
    control_cp_pimd_min(&class, &bonded, &general_data, &cp, &analysis);
  }

  /*===========================*/
  /* finalize and exit program */

  #if defined PLUMED
  plumed_piny_finalize();
  #endif

  if (class.communicate.np > 1) {
    Barrier(class.communicate.world);
    Finalize();
  }

  fflush(stdout);
  return 0;

}
/*--------------------------------------------------------------------------*/

