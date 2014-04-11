/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
/*                                                                          */
/*                         PI_MD:                                           */
/*             The future of simulation technology                          */
/*             ------------------------------------                         */
/*                     Module: auto_exit                                    */
/*                                                                          */
/* Allows user to stop a run before number of steps has been reached        */
/*                                                                          */
/*                                                                          */
/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/


#include "standard_include.h"
#include "../typ_defs/typedefs_class.h"
#include "../typ_defs/typedefs_bnd.h"
#include "../typ_defs/typedefs_gen.h"
#include "../typ_defs/typedefs_par.h"
#include "../typ_defs/typedefs_stat.h"
#include "../proto_defs/proto_main_local.h"

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void check_auto_exit(int *exit_flag)

/*=======================================================================*/
{ /*begin function */
/*=======================================================================*/
/*            Local variable declarations                                */

  const char *fn="EXIT";
  FILE *ef;

/*=======================================================================*/
/* Check to see if EXIT file exists.  If so, print out and exit         */

  if((ef=fopen(fn,"r"))==NULL) return;
  *exit_flag=1; 

  printf("@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n");
  printf("You have issued the automatic EXIT command.\n");
  printf("I will now perform one last step and\n");
  printf("then exit.  Thanks for using PINY_MD!\n");
  printf("We now return you to your regularly scheduled programming\n");
  printf("@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n");
  system("/bin/rm EXIT");

  fclose(ef);

/*-----------------------------------------------------------------------*/
  }/*end function */
/*========================================================================*/



