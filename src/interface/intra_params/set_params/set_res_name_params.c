/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
/*  set_res_name_params:set up the resname params                           */
/*==========================================================================*/

#include "standard_include.h"
#include "../typ_defs/typedefs_class.h"
#include "../typ_defs/typedefs_par.h"
#include "../typ_defs/typedefs_bnd.h"
#include "../proto_defs/proto_intra_params_local.h"
#include "../proto_defs/proto_handle_entry.h"

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void set_res_name_params(DICT_WORD res_name_dict[],int num_res_name_dict,
                         char fun_key[],char file_name[],NAME res_typ[],
                         int natm_jres_jmol_typ[],int jmol_typ,
                         int jres,int jres_off,
                         int ires_typ_jres_jmol_typ[])

/*==========================================================================*/
{
  int num,i,index,index2;
  /*=======================================================================*/
  /* I) Check for missing key words*/
  for(i=1;i<=num_res_name_dict;i++){
    if(res_name_dict[i].iuset==0 && res_name_dict[i].key_type==1){
      keyword_miss(res_name_dict,file_name,fun_key,i);}
  }    /*endfor*/
  /*======================================================================*/
  /* II) Set the keywords */
  /*-------------------------------------------------------------------------*/
  /* 1) \residue_name{} */

  index = 1;  
  index2 = ires_typ_jres_jmol_typ[(jres+jres_off)];
  if(strcasecmp(res_typ[index2],res_name_dict[1].keyarg)!=0){
    keyarg_barf(res_name_dict,file_name,fun_key,index);}

  /*-------------------------------------------------------------------------*/
  /* 3) \natom{} */

  sscanf(res_name_dict[3].keyarg,"%d",&num);
  index = 3;
  if(num<0||num!=natm_jres_jmol_typ[(jres+jres_off)]){
    keyarg_barf(res_name_dict,file_name,fun_key,index);
  }
  natm_jres_jmol_typ[jres+jres_off] = num;

  /*---------------------------------------------------------------------*/
}  /*end routine*/
/*==========================================================================*/





