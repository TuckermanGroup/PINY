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
#include "../proto_defs/proto_inter_params_local.h"
#include "../proto_defs/proto_friend_lib_entry.h"

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
/*  set_potinter_dict:set up the atom  key word dictionary                  */ 
/*==========================================================================*/

void set_potinter_dict(DICT_WORD *inter_dict[],int *num_inter_dict, int ifirst)

/*=======================================================================*/
/* Begin Routine */
{/*begin routine*/
/*=======================================================================*/
/* Local Variables */

  int i;

/*=======================================================================*/
/* 0) Malloc the dictionary */

  if(ifirst==1){
    *num_inter_dict=16;
    *inter_dict = (DICT_WORD *)cmalloc(*num_inter_dict*sizeof(DICT_WORD))-1;  
  }/*endif*/

/*=======================================================================*/
/* I) Assign the users set flags 0 */

  for(i=1;i<=*num_inter_dict;i++){(*inter_dict)[i].iuset    = 0;}
  for(i=1;i<=*num_inter_dict;i++){(*inter_dict)[i].key_type = 0;}

/*=======================================================================*/
/* II) Fill the dictionary with words */
  /*-----------------------------------------------------------------------*/ 
  /*  1) /atom1{} */
  strcpy((*inter_dict)[1].keyword,"atom1");
  strcpy((*inter_dict)[1].keyarg,"");
  strcpy((*inter_dict)[1].error_mes,"");
  (*inter_dict)[1].key_type = 1;  /* must spec*/
  /*-----------------------------------------------------------------------*/ 
  /*  2) /atom2{} */
  strcpy((*inter_dict)[2].keyword,"atom2");
  strcpy((*inter_dict)[2].keyarg,"");
  strcpy((*inter_dict)[2].error_mes,"");
  (*inter_dict)[2].key_type = 1;  /* must spec*/
  /*-----------------------------------------------------------------------*/ 
  /*  3) /pot_type{} */
  strcpy((*inter_dict)[3].keyword,"pot_type");
  strcpy((*inter_dict)[3].keyarg,"");
  strcpy((*inter_dict)[3].error_mes,"lennard-jones,williams,aziz-chen,null");
  (*inter_dict)[3].key_type = 1;  /* must spec*/
  /*-----------------------------------------------------------------------*/ 
  /*  4) /min_dist{} */
  strcpy((*inter_dict)[4].keyword,"min_dist");
  strcpy((*inter_dict)[4].keyarg,"");
  strcpy((*inter_dict)[4].keyarg,"0.1");
  strcpy((*inter_dict)[4].error_mes,"a number >0");
  /*-----------------------------------------------------------------------*/ 
  /*  5) /max_dist{} */
  strcpy((*inter_dict)[5].keyword,"max_dist");
  strcpy((*inter_dict)[5].keyarg,"10.0");
  strcpy((*inter_dict)[5].error_mes,"a number >0");
  (*inter_dict)[5].key_type = 1;  /* must spec the cutoff in case coulomb on*/
  /*-----------------------------------------------------------------------*/ 
  /*  6) /res_dist{} */
  strcpy((*inter_dict)[6].keyword,"res_dist");
  strcpy((*inter_dict)[6].keyarg,"10.0");
  strcpy((*inter_dict)[6].error_mes,"a number >0");
  /*-----------------------------------------------------------------------*/ 
  /*  7) /sig{} */
  strcpy((*inter_dict)[7].keyword,"sig");
  strcpy((*inter_dict)[7].keyarg,"0.0");
  strcpy((*inter_dict)[7].error_mes,"a number >=0");
  /*-----------------------------------------------------------------------*/ 
  /*  8) /eps{} */
  strcpy((*inter_dict)[8].keyword,"eps");
  strcpy((*inter_dict)[8].keyarg,"0.0");
  strcpy((*inter_dict)[8].error_mes,"a number >=0");
  /*-----------------------------------------------------------------------*/ 
  /*  9) /c6{} */
  strcpy((*inter_dict)[9].keyword,"c6");
  strcpy((*inter_dict)[9].keyarg,"0.0");
  strcpy((*inter_dict)[9].error_mes,"a number >=0");
  /*-----------------------------------------------------------------------*/ 
  /* 10) /c8{} */
  strcpy((*inter_dict)[10].keyword,"c8");
  strcpy((*inter_dict)[10].keyarg,"0.0");
  strcpy((*inter_dict)[10].error_mes,"a number 0");
  /*-----------------------------------------------------------------------*/ 
  /* 11) /c10{} */
  strcpy((*inter_dict)[11].keyword,"c10");
  strcpy((*inter_dict)[11].keyarg,"0.0");
  strcpy((*inter_dict)[11].error_mes,"a number 0");
  /*-----------------------------------------------------------------------*/ 
  /* 12) /Awill{} */
  strcpy((*inter_dict)[12].keyword,"Awill");
  strcpy((*inter_dict)[12].keyarg,"0.0");
  strcpy((*inter_dict)[12].error_mes,"a number >=0");
  /*-----------------------------------------------------------------------*/ 
  /* 13) /Bwill{} */
  strcpy((*inter_dict)[13].keyword,"Bwill");
  strcpy((*inter_dict)[13].keyarg,"0.0");
  strcpy((*inter_dict)[13].error_mes,"a number >=0");
  /*-----------------------------------------------------------------------*/ 
  /* 14) /Cwill{} */
  strcpy((*inter_dict)[14].keyword,"Cwill");
  strcpy((*inter_dict)[14].keyarg,"0.0");
  strcpy((*inter_dict)[14].error_mes,"a number >=0");
  /*-----------------------------------------------------------------------*/ 
  /* 15) /rm_swit{} */
  strcpy((*inter_dict)[15].keyword,"rm_swit");
  strcpy((*inter_dict)[15].keyarg,"0.0");
  strcpy((*inter_dict)[15].error_mes,"a number >=0");
  /*-----------------------------------------------------------------------*/ 
  /* 16) /c9{} */
  strcpy((*inter_dict)[16].keyword,"c9");
  strcpy((*inter_dict)[16].keyarg,"0.0");
  strcpy((*inter_dict)[16].error_mes,"a number 0");
  /*-----------------------------------------------------------------------*/ 
}  /*end routine*/
/*==========================================================================*/

