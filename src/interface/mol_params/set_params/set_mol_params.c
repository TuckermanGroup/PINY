/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
/*                                                                          */
/*                         PI_MD:                                           */
/*             The future of simulation technology                          */
/*             ------------------------------------                         */
/*                     Module: set_mol_parms.c                              */
/*                                                                          */
/* This subprogram sets molecule/cp setup parameters                        */
/*                                                                          */
/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

#include "standard_include.h"
#include "../typ_defs/typedefs_gen.h"
#include "../typ_defs/typedefs_par.h"
#include "../typ_defs/typedefs_bnd.h"
#include "../typ_defs/typedefs_cp.h"
#include "../typ_defs/typedefs_class.h"
#include "../proto_defs/proto_mol_params_local.h"
#include "../proto_defs/proto_friend_lib_entry.h"
#include "../proto_defs/proto_handle_entry.h"

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
void set_mol_params(FILENAME_PARSE*filename_parse,char fun_key[],
		    DICT_WORD mol_dict[], int num_mol_dict,
		    CLASS_PARSE *class_parse,ATOMMAPS *atommaps,int mol_ind_chk[],
                    int pi_beads)
{
  int i,num,index,mol_ind,ifound,inum;
  char *strip,*strip2;
  double dnum;
  /*=======================================================================*/
  /* I) Check for missing key words*/
  
  for(i=1;i<num_mol_dict;i++){
    if(mol_dict[i].iuset==0 && mol_dict[i].key_type==1){
      keyword_miss(mol_dict,filename_parse->molsetname,fun_key,i);}
  }    /*endfor*/

  /*======================================================================*/
  /* II) Set the params */
  /*-----------------------------------------------------------------------*/ 
  /*  1)\mol_index{} */

  sscanf(mol_dict[1].keyarg,"%d",&mol_ind);
  index = 1;
  if((mol_ind <=0)|| (mol_ind> atommaps->nmol_typ)){
    keyarg_barf(mol_dict,filename_parse->molsetname,fun_key,index);
  }
  mol_ind_chk[mol_ind] = mol_ind_chk[mol_ind] + 1;

  /*-----------------------------------------------------------------------*/ 
  /*  2)\num_mol{} */

  sscanf(mol_dict[2].keyarg,"%d",&num);
  atommaps->nmol_jmol_typ[mol_ind] = num;
  index = 2;
  if((num<=0)){
    keyarg_barf(mol_dict,filename_parse->molsetname,fun_key,index);
  }
  
  /*-----------------------------------------------------------------------*/ 
  /*  3)\mol_parm_file{} */

  strcpy(filename_parse->mol_param_name[mol_ind],mol_dict[3].keyarg);

  /*-----------------------------------------------------------------------*/ 
  /*  4)\mol_nhc_opt{} */

  ifound = 0;
  if(strcasecmp(mol_dict[4].keyarg,"none")==0){
    ifound=1; class_parse->imol_nhc_opt[mol_ind]=0;
  }
  if(strcasecmp(mol_dict[4].keyarg,"global")==0){
    ifound=1; class_parse->imol_nhc_opt[mol_ind]=1;
  }
  if(strcasecmp(mol_dict[4].keyarg,"glob_mol")==0){
    ifound=1; class_parse->imol_nhc_opt[mol_ind]=2;
  }
  if(strcasecmp(mol_dict[4].keyarg,"ind_mol")==0){
    ifound=1;class_parse->imol_nhc_opt[mol_ind]=3;
  }
  if(strcasecmp(mol_dict[4].keyarg,"res_mol")==0){
    ifound=1;class_parse->imol_nhc_opt[mol_ind]=4;
  }
  if(strcasecmp(mol_dict[4].keyarg,"atm_mol")==0){
    ifound=1;class_parse->imol_nhc_opt[mol_ind]=5;
  }
  if(strcasecmp(mol_dict[4].keyarg,"mass_mol")==0){
    ifound=1;class_parse->imol_nhc_opt[mol_ind]=6;
  }
  index = 4;
  if(ifound==0){
    keyarg_barf(mol_dict,filename_parse->molsetname,fun_key,index);
  }

/*-----------------------------------------------------------------------*/ 
/*  5)\mol_tau_nhc{} */

  sscanf(mol_dict[5].keyarg,"%lg",&dnum);
  class_parse->tau_nhc_mol[mol_ind]=dnum;
  index = 5;
  if((dnum <=0)){
    keyarg_barf(mol_dict,filename_parse->molsetname,fun_key,index);
  }

/*-----------------------------------------------------------------------*/ 
/*  6)\mol_text_nhc{} */

  sscanf(mol_dict[6].keyarg,"%lg",&dnum);
  class_parse->text_nhc_mol[mol_ind]=dnum;
  index = 6;
  if((dnum <=0)){
    keyarg_barf(mol_dict,filename_parse->molsetname,fun_key,index);}

/*-----------------------------------------------------------------------*/ 
/*  7)\mol_name{} */

  strcpy(atommaps->mol_typ[mol_ind],mol_dict[7].keyarg);
  
/*-----------------------------------------------------------------------*/ 
/*  8)\num_residue{} */

  sscanf(mol_dict[8].keyarg,"%d",&inum);
  index = 8;
  if(inum <0){
    keyarg_barf(mol_dict,filename_parse->molsetname,fun_key,index);   
  }/*endif*/
  atommaps->nres_1mol_jmol_typ[mol_ind] = inum;

/*-----------------------------------------------------------------------*/ 
/*  9)\onefour_opt{} */

  index  = 9;
  ifound = 0;
  if(strcasecmp(mol_dict[9].keyarg,"on")==0){
    ifound=1;class_parse->ionfo_opt[mol_ind] = 1;
  }
  if(strcasecmp(mol_dict[9].keyarg,"off")==0){
    ifound=1;class_parse->ionfo_opt[mol_ind] = 0;
  }
  if(ifound==0){
    keyarg_barf(mol_dict,filename_parse->molsetname,fun_key,index);    
  }/*endif*/

/*-----------------------------------------------------------------------*/ 
/*  10)\res_bond_convention{} */

  index  = 10;
  ifound = 0;
  if(strcasecmp(mol_dict[10].keyarg,"pi_md")==0){
    ifound=1;class_parse->ires_bond_conv[mol_ind] = 1;
  }
  if(strcasecmp(mol_dict[10].keyarg,"charmm")==0){
    ifound=1;class_parse->ires_bond_conv[mol_ind] = 2;
  }
  if(strcasecmp(mol_dict[10].keyarg,"sybyl")==0){
    ifound=1;class_parse->ires_bond_conv[mol_ind] = 3;
  }
  if(ifound==0){
    keyarg_barf(mol_dict,filename_parse->molsetname,fun_key,index);    
  }/*endif*/

/*-----------------------------------------------------------------------*/ 
/*  11)\mol_freeze_opt{} */

  index  = 11;
  ifound = 0;
  if(strcasecmp(mol_dict[11].keyarg,"none")==0){
    ifound=1;class_parse->mol_freeze_opt[mol_ind] = 0;
  }
  if(strcasecmp(mol_dict[11].keyarg,"all")==0){
    ifound=1;class_parse->mol_freeze_opt[mol_ind] = 1;
  }
  if(strcasecmp(mol_dict[11].keyarg,"backbone")==0){
    ifound=1;class_parse->mol_freeze_opt[mol_ind] = 2;
  }
  if(strcasecmp(mol_dict[11].keyarg,"heavy_atoms")==0){
    ifound=1;class_parse->mol_freeze_opt[mol_ind] = 3;
  }
  if(ifound==0){
    keyarg_barf(mol_dict,filename_parse->molsetname,fun_key,index);    
  }/*endif*/

/*-----------------------------------------------------------------------*/
/*  12)\hydrog_mass_opt{} */
  index  = 12;
  strip  = (char *)cmalloc(MAXWORD*sizeof(char));
  strip2 = (char *)cmalloc(MAXWORD*sizeof(char));
  strcpy(strip,mol_dict[12].keyarg);
  parse_hydrog_mass(strip,strip2,&(class_parse->mol_hydrog_mass_opt[mol_ind]),
                 &(class_parse->mol_hydrog_mass_val[mol_ind]),&ifound);
  if(ifound==0){
    keyarg_barf(mol_dict,filename_parse->molsetname,fun_key,index);
  }/*endif*/
  cfree(strip);
  cfree(strip2);

/*-----------------------------------------------------------------------*/
/*  13)\hydrog_con_opt{} */
  index  = 13;
  ifound = 0;
  if(strcasecmp(mol_dict[13].keyarg,"off")==0){
    ifound=1;class_parse->mol_hydrog_con_opt[mol_ind] = 0;
  }
  if(strcasecmp(mol_dict[13].keyarg,"all")==0){
    ifound=1;class_parse->mol_hydrog_con_opt[mol_ind] = 1;
  }
  if(strcasecmp(mol_dict[13].keyarg,"polar")==0){
    ifound=1;class_parse->mol_hydrog_con_opt[mol_ind] = 2;
  }
  if(ifound==0){
    keyarg_barf(mol_dict,filename_parse->molsetname,fun_key,index);
  }/*endif*/
  if(class_parse->mol_hydrog_con_opt[mol_ind] != 0 && pi_beads>1){
      printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");        
      printf("The code, while most excellent, cannot as yet perform  \n");
      printf("path integral simulations with constrained bonds, which\n");
      printf("are somewhat nasty.  Contact sales rep.       \n");
      printf("Turn off the \\hydrog_con option \n");
      printf("in file %s molecule %s\n",
              filename_parse->molsetname,atommaps->mol_typ[mol_ind]);
      printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
      fflush(stdout);
      exit(1);
  }


/*==========================================================================*/
/* III) Error checking */

   if((atommaps->nres_1mol_jmol_typ[mol_ind]==0)&&
      (class_parse->imol_nhc_opt[mol_ind]==4)){
     printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");        
     printf("The molecule type specified %s in file %s      \n",
             atommaps->mol_typ[mol_ind],filename_parse->molsetname);
     printf("does not contain residues. Therefore, the res_mol \n");
     printf("thermostat option is not an appropriate choice    \n");
     printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
     fflush(stdout);
     exit(1);
   }/*endif*/

/*-----------------------------------------------------------------------*/ 
}  /*end routine*/
/*==========================================================================*/






