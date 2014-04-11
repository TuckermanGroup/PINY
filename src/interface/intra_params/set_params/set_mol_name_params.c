/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
/*  set_mol_name_params:set up the molname params                           */
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

void set_mol_name_params(DICT_WORD mol_name_dict[],int num_mol_name_dict,
                         char fun_key[],char filename[],NAME mol_typ[],
                         int natm_1mol_jmol_typ[],int jmol_typ,int nres_expt)

/*==========================================================================*/
{
  int num,i,index;

/*=======================================================================*/
/* I) Check for missing key words*/

  for(i=1;i<=num_mol_name_dict;i++){
    if(mol_name_dict[i].iuset==0 && mol_name_dict[i].key_type==1){
      keyword_miss(mol_name_dict,filename,fun_key,i);}
  }    /*endfor*/

/*======================================================================*/
/* II) Set the keywords */
/*-------------------------------------------------------------------------*/
  /* 1) \mol_name{} */
  index = 1;  
  if(strcasecmp(mol_typ[jmol_typ],mol_name_dict[1].keyarg)!=0){
    keyarg_barf(mol_name_dict,filename,fun_key,index);}
  /*-------------------------------------------------------------------------*/
  /* 2) \nresidue{} */
  sscanf(mol_name_dict[2].keyarg,"%d",&num);
  index = 2;
  if(num < 0){
    keyarg_barf(mol_name_dict,filename,fun_key,index);
  }
  if(num != nres_expt){
   printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");    
   printf("Number of residue's specified %d not equal to number expected %d\n",
          num,nres_expt);
   printf("in file %s\n",filename);
   printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");    
   fflush(stdout);
   exit(1);    
  }
  /*-------------------------------------------------------------------------*/
  /* 3) \natom{} */
  sscanf(mol_name_dict[3].keyarg,"%d",&num);
  index = 3;
  
  if(num<0){
    keyarg_barf(mol_name_dict,filename,fun_key,index);
  }
  natm_1mol_jmol_typ[jmol_typ] = num;
  /*---------------------------------------------------------------------*/
}  /*end routine*/
/*==========================================================================*/




