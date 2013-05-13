/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
/*                                                                          */
/*                         PI_MD:                                           */
/*             The future of simulation technology                          */
/*             ------------------------------------                         */
/*                     Module: set_res_def_parms.c                          */
/*                                                                          */
/* This subprogram sets molecule/cp setup parameters                        */
/*                                                                          */
/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

#include "standard_include.h"
#include "../typ_defs/typedefs_par.h"
#include "../typ_defs/typedefs_class.h"
#include "../typ_defs/typedefs_bnd.h"
#include "../proto_defs/proto_intra_params_local.h"
#include "../proto_defs/proto_friend_lib_entry.h"
#include "../proto_defs/proto_handle_entry.h"

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void set_res_def_params(char filename[],char fun_key[],
			DICT_WORD res_def_dict[], int num_res_def_dict,
			ATOMMAPS *atommaps,NAME res_param_name[],
			int nres_now,
                        BUILD_INTRA *build_intra,int ires_off,
                        FILENAME_PARSE *filename_parse)

/*==========================================================================*/
/*      Begin Routine                                                       */
{ /*begin routine */

/*==========================================================================*/
/*      Local Variables */
  int i,num,index,ires_ind,ifound;

/*=======================================================================*/
/* I) Check for missing key words  */
  
  for(i=1;i<num_res_def_dict;i++){
    if(res_def_dict[i].iuset==0 && res_def_dict[i].key_type==1){
      keyword_miss(res_def_dict,filename,fun_key,i);
    }
  }    /*endfor*/

/*======================================================================*/
/* II) Set the params */
/*-----------------------------------------------------------------------*/ 
/*  1)\residue_index{} */

  sscanf(res_def_dict[1].keyarg,"%d",&ires_ind);
  index = 1;
  if((ires_ind <=0)|| (ires_ind > nres_now)){
    keyarg_barf(res_def_dict,filename,fun_key,index);
  }
  build_intra->ires_ind_chk[ires_ind]++;

  /*-----------------------------------------------------------------------*/ 
  /*  2)\natom{} */

  sscanf(res_def_dict[2].keyarg,"%d",&num);
  atommaps->natm_jres_jmol_typ[ires_ind+ires_off] = num;
  index = 2;
  if((num<=0)){
    keyarg_barf(res_def_dict,filename,fun_key,index);
  }
  
  /*-----------------------------------------------------------------------*/ 
  /*  3)\residue_def_parm_file{} */

  strcpy(filename_parse->res_param_name[ires_ind],res_def_dict[3].keyarg);

  /*-----------------------------------------------------------------------*/ 
  /*  4)\residue_def_name{} */

  ifound = 0;
  for(i=1;i<=atommaps->nres_typ;i++){
    if(strcasecmp(res_def_dict[4].keyarg,atommaps->res_typ[i])==0){
      ifound = i;
      break;
    }
  } /*endfor */
  if(ifound == 0){
    atommaps->nres_typ++;
    if(atommaps->nres_typ > build_intra->nres_typ_max){
      build_intra->nres_typ_max+=NMEM_MIN;
      atommaps->res_typ = (NAME *)
        crealloc(&(atommaps->res_typ[1]),
                  (build_intra->nres_typ_max)*sizeof(NAME))-1;
    }/*endif*/
    strcpy(atommaps->res_typ[atommaps->nres_typ],res_def_dict[4].keyarg);
    ifound = atommaps->nres_typ;
  }
  atommaps->ires_typ_jres_jmol_typ[ires_ind+ires_off] = ifound;
  /*-----------------------------------------------------------------------*/ 
  /* 5 )\residue_def_fix_file{} */
  strcpy(filename_parse->res_fix_name[ires_ind],res_def_dict[5].keyarg);

/*-----------------------------------------------------------------------*/ 
}  /*end routine*/
/*==========================================================================*/

