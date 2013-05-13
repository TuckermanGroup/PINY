/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
/*                                                                          */
/*                         PI_MD:                                           */
/*             The future of simulation technology                          */
/*             ------------------------------------                         */
/*                     Module: set_inter_dict.c                             */
/*                                                                          */
/* These subprograms set the intermolecular dictionaries                    */
/*                                                                          */
/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
/*               Includes:                                                  */

#include "standard_include.h"
#include "../typ_defs/typedefs_par.h"
#include "../typ_defs/typedefs_class.h"
#include "../proto_defs/proto_surf_params_local.h"
#include "../proto_defs/proto_friend_lib_entry.h"

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
/*  set_potinter_dict:set up the atom  key word dictionary                  */ 
/*==========================================================================*/

void set_potsurf_dict(DICT_WORD *surf_dict[],int *num_surf_dict, int ifirst)

/*=======================================================================*/
/* Begin Routine */
{/*begin routine*/
/*=======================================================================*/
/* Local Variables */

  int i;

/*=======================================================================*/
/* 0) Malloc the dictionary */

  if(ifirst==1){
    *num_surf_dict=7;
    *surf_dict = (DICT_WORD *)cmalloc(*num_surf_dict*sizeof(DICT_WORD))-1;  
  }/*endif*/

/*=======================================================================*/
/* I) Assign the users set flags 0 */

  for(i=1;i<=*num_surf_dict;i++){(*surf_dict)[i].iuset    = 0;}
  for(i=1;i<=*num_surf_dict;i++){(*surf_dict)[i].key_type = 0;}

/*=======================================================================*/
/* II) Fill the dictionary with words */
  /*-----------------------------------------------------------------------*/ 
  /*  1) /surface_type{} */
  strcpy((*surf_dict)[1].keyword,"surface_type");
  strcpy((*surf_dict)[1].keyarg,"");
  strcpy((*surf_dict)[1].error_mes,"");
  (*surf_dict)[1].key_type = 1;  /* must spec*/
  /*-----------------------------------------------------------------------*/ 
  /*  2) /ad_atom{} */
  strcpy((*surf_dict)[2].keyword,"ad_atom");
  strcpy((*surf_dict)[2].keyarg,"");
  strcpy((*surf_dict)[2].error_mes,"");
  (*surf_dict)[2].key_type = 1;  /* must spec*/
  /*-----------------------------------------------------------------------*/ 
  /*  3) /pot_type{} */
  strcpy((*surf_dict)[3].keyword,"pot_type");
  strcpy((*surf_dict)[3].keyarg,"");
  strcpy((*surf_dict)[3].error_mes,
       "lennard-jones_12-3,lennard-jones_9-3,null");
  (*surf_dict)[3].key_type = 1;  /* must spec*/
  /*-----------------------------------------------------------------------*/ 
  /*  4) /min_dist{} */
  strcpy((*surf_dict)[4].keyword,"min_dist");
  strcpy((*surf_dict)[4].keyarg,"");
  strcpy((*surf_dict)[4].keyarg,"0.1");
  strcpy((*surf_dict)[4].error_mes,"a number >0");
  /*-----------------------------------------------------------------------*/ 
  /*  5) /max_dist{} */
  strcpy((*surf_dict)[5].keyword,"max_dist");
  strcpy((*surf_dict)[5].keyarg,"10.0");
  strcpy((*surf_dict)[5].error_mes,"a number >0");
  (*surf_dict)[5].key_type = 1;  /* must spec the cutoff in case coulomb on*/
  /*-----------------------------------------------------------------------*/ 
  /*  6) /sig{} */
  strcpy((*surf_dict)[6].keyword,"sig");
  strcpy((*surf_dict)[6].keyarg,"0.0");
  strcpy((*surf_dict)[6].error_mes,"a number >=0");
  /*-----------------------------------------------------------------------*/ 
  /*  7) /eps{} */
  strcpy((*surf_dict)[7].keyword,"eps");
  strcpy((*surf_dict)[7].keyarg,"0.0");
  strcpy((*surf_dict)[7].error_mes,"a number >=0");
  /*-----------------------------------------------------------------------*/ 
}  /*end routine*/
/*==========================================================================*/

