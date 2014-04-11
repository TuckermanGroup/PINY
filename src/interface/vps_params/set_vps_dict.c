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
#include "../proto_defs/proto_vps_params_local.h"
#include "../proto_defs/proto_friend_lib_entry.h"

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
/*  set_potvps_dict:set up the atom  key word dictionary                    */ 
/*==========================================================================*/
void set_potvps_dict(DICT_WORD *vps_dict[],int *num_vps_dict, int ifirst)
/*==========================================================================*/
{/*begin routine*/
int i;
/*==========================================================================*/
/* 0) Malloc the dictionary */
 if(ifirst==1){
  *num_vps_dict=11;
  *vps_dict    = (DICT_WORD *)cmalloc(*num_vps_dict*sizeof(DICT_WORD))-1;  
 /*endif*/}
/*==========================================================================*/
/* I) Assign the users set flags 0 */
  for(i=1;i<=*num_vps_dict;i++){(*vps_dict)[i].iuset    = 0;}
  for(i=1;i<=*num_vps_dict;i++){(*vps_dict)[i].key_type = 0;}
/*==========================================================================*/
/* II) Fill the dictionary with words */
/*-------------------------------------------------------------------------*/ 
/*  1) /atom1{} */
  strcpy((*vps_dict)[1].keyword,"atom1");
  strcpy((*vps_dict)[1].keyarg,"");
  strcpy((*vps_dict)[1].error_mes,"");
  (*vps_dict)[1].key_type = 1;  /* must spec*/
/*--------------------------------------------------------------------------*/ 
/*  2) /vps_typ{} */
  strcpy((*vps_dict)[2].keyword,"vps_typ");
  strcpy((*vps_dict)[2].keyarg,"");
  strcpy((*vps_dict)[2].error_mes,"local,kb,gauss_hermite,vdb,null,goedecker");
  (*vps_dict)[2].key_type = 1;  /* must spec*/
/*--------------------------------------------------------------------------*/ 
/*  3) /vps_file{} */
  strcpy((*vps_dict)[3].keyword,"vps_file");
  strcpy((*vps_dict)[3].keyarg,"");
  strcpy((*vps_dict)[3].error_mes,"file not found");
  (*vps_dict)[3].key_type = 1;  /* must spec*/
/*--------------------------------------------------------------------------*/ 
/*  4) /n_ang{} */
  strcpy((*vps_dict)[4].keyword,"n_ang");
  strcpy((*vps_dict)[4].keyarg,"0");
  strcpy((*vps_dict)[4].error_mes,"a number >= 0 and <=3");
/*--------------------------------------------------------------------------*/ 
/*  5) /loc_opt{} */
  strcpy((*vps_dict)[5].keyword,"loc_opt");
  strcpy((*vps_dict)[5].keyarg,"0");
  strcpy((*vps_dict)[5].error_mes,"a number <= n_ang >=0");
/*--------------------------------------------------------------------------*/ 
/*  6) /rcut_nl{} */
  strcpy((*vps_dict)[6].keyword,"rcut_nl");
  strcpy((*vps_dict)[6].keyarg,"10000.0");
  strcpy((*vps_dict)[6].error_mes,"a number >= 0");
/*--------------------------------------------------------------------------*/
/*  7) /nrad_0{} */
  strcpy((*vps_dict)[7].keyword,"nrad_0");
  strcpy((*vps_dict)[7].keyarg,"1");
  strcpy((*vps_dict)[7].error_mes,"a number >= 0");
/*--------------------------------------------------------------------------*/
/*  8) /nrad_1{} */
  strcpy((*vps_dict)[8].keyword,"nrad_1");
  strcpy((*vps_dict)[8].keyarg,"1");
  strcpy((*vps_dict)[8].error_mes,"a number >= 0");
/*--------------------------------------------------------------------------*/
/*  9) /nrad_2{} */
  strcpy((*vps_dict)[9].keyword,"nrad_2");
  strcpy((*vps_dict)[9].keyarg,"1");
  strcpy((*vps_dict)[9].error_mes,"a number >= 0");
/*--------------------------------------------------------------------------*/
/*  10) /nrad_3{} */
  strcpy((*vps_dict)[10].keyword,"nrad_3");
  strcpy((*vps_dict)[10].keyarg,"1");
  strcpy((*vps_dict)[10].error_mes,"a number >= 0");
/*--------------------------------------------------------------------------*/
/*  11) /num_gauss_hermite{} */
  strcpy((*vps_dict)[11].keyword,"num_gauss_hermite");
  strcpy((*vps_dict)[11].keyarg,"0");
  strcpy((*vps_dict)[11].error_mes,"An even number > 0 and <= 180");
/*--------------------------------------------------------------------------*/ 
/*end routine*/}
/*==========================================================================*/










