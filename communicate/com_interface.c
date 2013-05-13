/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
/*                                                                          */
/*                         PI_MD:                                           */
/*             The future of simulation technology                          */
/*             ------------------------------------                         */
/*                     Module: comm_interface.c                             */
/*                                                                          */
/* Subprogram contains MPI utils and communication routines for interface   */
/*                                                                          */
/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/


#include "standard_include.h"
#include "../typ_defs/typedefs_gen.h"
#include "../typ_defs/typedefs_par.h"
#include "../typ_defs/typedefs_bnd.h"
#include "../typ_defs/typedefs_class.h"
#include "../typ_defs/typedefs_cp.h"
#include "../proto_defs/proto_communicate_entry.h"
#include "../proto_defs/proto_communicate_wrappers.h"
#include "../proto_defs/proto_communicate_local.h"

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void communicate_interface(CLASS *class,BONDED *bonded,CP *cp,
                GENERAL_DATA *general_data,NULL_INTER_PARSE *null_inter_parse,
                CLASS_PARSE *class_parse,CP_PARSE *cp_parse)

/*=======================================================================*/
/*             Begin routine                                              */
{/*begin routine */

#include "../typ_defs/typ_mask.h"


/*=======================================================================*/
/*             Local variable declarations                                */

  int pimd_on,cp_on,natm_mall,natm_typ_mall,iii;
  int myid = class->communicate.myid;
  MPI_Comm world = class->communicate.world;

/*=======================================================================*/
/*            Output                                                     */

  if(myid==0){
      PRINT_LINE_STAR;
      printf("Communicating the interface\n");
      PRINT_LINE_DASH;printf("\n");
  }/*endif*/

/*=======================================================================*/
/*            General_Data                                               */

  Barrier(world);
  communicate_general_data(general_data,cp,world);
  Barrier(world);

  cp_on = general_data->simopts.cp_min+general_data->simopts.cp_wave_min
        +general_data->simopts.cp+general_data->simopts.cp_wave
        +general_data->simopts.debug_cp+general_data->simopts.cp_pimd
        +general_data->simopts.debug_cp_pimd+general_data->simopts.cp_wave_pimd
        +general_data->simopts.cp_wave_min_pimd;

  pimd_on = general_data->simopts.pimd + general_data->simopts.cp_pimd 
          + general_data->simopts.cp_wave_pimd 
          + general_data->simopts.cp_wave_min_pimd 
          + general_data->simopts.debug_pimd 
          + general_data->simopts.debug_cp_pimd;

/*=======================================================================*/
/*           Class                                                     */

   Barrier(world);
   communicate_class_info(class,general_data,class_parse);
   Barrier(world);
   if(myid!=0){
        mall_class(class,general_data,class_parse,pimd_on);    
   }/*endif*/
   Barrier(world);
   communicate_class_data(class,general_data,class_parse,pimd_on);
   Barrier(world);
   communicate_class_list(class,general_data,class_parse,pimd_on);

/*=======================================================================*/
/*            Bonded                                                     */

   Barrier(world);
   communicate_bond_info(bonded,null_inter_parse,world); 
   Barrier(world);
   if(myid!=0){
       mall_bond(bonded,null_inter_parse); 
   }/*endif*/
   Barrier(world);
   communicate_bond_data(bonded,null_inter_parse,world); 
   Barrier(world);

/*=======================================================================*/
/*            CP                                                         */

   communicate_cp_info(cp,cp_parse,world,myid); 
   Barrier(world);
   natm_mall     = class->clatoms_info.natm_mall;
   natm_typ_mall = class->atommaps.natm_typ_mall;

   if(cp_on==1){
     communicate_cp_data(cp,natm_mall,natm_typ_mall,world,myid);
   }/*endif cp_on==1*/
   Barrier(world);

/*=======================================================================*/
/*         More Output                                                   */

   if(myid==0){
      printf("\n");PRINT_LINE_DASH;
      printf("Completed communicating the interface\n");
      PRINT_LINE_STAR;printf("\n");
   }/*endif*/

/*------------------------------------------------------------------------*/
} /*end routine*/ 
/*==========================================================================*/











