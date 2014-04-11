/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
/*                                                                          */
/*                         PI_MD:                                           */
/*             The future of simulation technology                          */
/*             ------------------------------------                         */
/*                     Module: set_wave_params.c                            */
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
#include "../typ_defs/typedefs_class.h"
#include "../typ_defs/typedefs_cp.h"
#include "../proto_defs/proto_mol_params_local.h"
#include "../proto_defs/proto_handle_entry.h"

/*==========================================================================*/

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
void set_wave_params(char *molsetname,char fun_key[],
		     DICT_WORD wave_dict[], int num_wave_dict,
		     CPOPTS *cpopts,CPCOEFFS_INFO *cpcoeffs_info,
                     CP_PARSE *cp_parse)
{
  int i,num,index,ifound;
  double dnum;
  /*=======================================================================*/
  /* I) Check for missing key words*/
  for(i=1;i<num_wave_dict;i++){
    if(wave_dict[i].iuset==0 && wave_dict[i].key_type==1){
      keyword_miss(wave_dict,molsetname,fun_key,i);}
  }    /*endfor*/
  /*=======================================================================*/
  /* II) Set the params */
  /*-----------------------------------------------------------------------*/ 
  /* 1) /nstate_up{} */
  sscanf(wave_dict[1].keyarg,"%d",&num);
  cpcoeffs_info->nstate_up = num;
  index = 1;
  if((num <=0)){
    keyarg_barf(wave_dict,molsetname,fun_key,index);
  }
  /*-----------------------------------------------------------------------*/ 
  /* 2) /nstate_dn{} */
  sscanf(wave_dict[2].keyarg,"%d",&num);
  cpcoeffs_info->nstate_dn = num;

  if(cpopts->cp_lda == 1){ 
      cpcoeffs_info->nstate_dn = cpcoeffs_info->nstate_up;
  }

  index = 2;
  if(num < 0){
    keyarg_barf(wave_dict,molsetname,fun_key,index);
  }

  if( (cpopts->cp_lda==1) && 
      (cpcoeffs_info->nstate_up != num)){
    keyarg_barf(wave_dict,molsetname,fun_key,index);
  }
    
  /*-----------------------------------------------------------------------*/ 
  /* 3) /cp_nhc_opt{} */
  ifound = 0;
  if(strcasecmp(wave_dict[3].keyarg,"none")==0)
    {ifound=1;cp_parse->istate_nhc_opt=0;}
  if(strcasecmp(wave_dict[3].keyarg,"global")==0)
    {ifound=1;cp_parse->istate_nhc_opt=1;}
  if(strcasecmp(wave_dict[3].keyarg,"glob_st")==0)
    {ifound=1;cp_parse->istate_nhc_opt=2;}
  if(strcasecmp(wave_dict[3].keyarg,"ind_st")==0)
    {ifound=1;cp_parse->istate_nhc_opt=3;}
  if(strcasecmp(wave_dict[3].keyarg,"mass_coef")==0)
    {ifound=1;cp_parse->istate_nhc_opt=4;}
  index = 3;
  if(ifound==0){
    keyarg_barf(wave_dict,molsetname,fun_key,index);
  }
  if(cp_parse->istate_nhc_opt==3 && cpopts->cp_norb < 2){
      printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
      printf("The \"ind_st\" thermostatting option for the\n");
      printf("coefficients, specified in the set up file %s, \n",
              molsetname);
      printf("while most excellent, only works with \"cp_norb\", \n");
      printf("using either the \"norm_only\" or \"none\"\n"); 
      printf("constraint options\n");
      printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
      fflush(stdout);
      exit(1);
  }/*endif*/
  if(cp_parse->istate_nhc_opt==4 && cpopts->cp_norb !=3){ 
      printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
      printf("The \"mass_coef\" thermostatting option for the\n");
      printf("coefficients, specified in the set up file %s, \n",
              molsetname);
      printf("while most excellent, only works with \"cp_norb\", \n");
      printf("using the \"none\" constraint option\n");
      printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
      fflush(stdout);
      exit(1);
 }/*endif*/
 if((cpopts->cp_isok_opt == 1) && (cp_parse->istate_nhc_opt != 0)){
      printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
      printf("The CP isokinetic option, while most excellent,\n");
      printf("will not work in tandem with the CP Nose-Hoover chain\n");
      printf("thermostatting scheme.  Please set the CP_NHC_OPT\n");
      printf("option in your .set file to NONE\n");
      printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
      fflush(stdout);
      exit(1);
 }/* endif */
  /*-----------------------------------------------------------------------*/ 
  /* 4) /cp_tau_nhc{} */
  sscanf(wave_dict[4].keyarg,"%lg",&dnum);
  cp_parse->cp_tau_nhc=dnum;
  index = 4;
  if((dnum <=0)){
    keyarg_barf(wave_dict,molsetname,fun_key,index);
  }
  /*-----------------------------------------------------------------------*/ 
  /* 5) /cp_energy_cut{} */
  sscanf(wave_dict[5].keyarg,"%lg",&dnum);
  cp_parse->cp_ecut=dnum;
  index = 5;
  if((dnum <=0)){
    keyarg_barf(wave_dict,molsetname,fun_key,index);
  }
  /*-----------------------------------------------------------------------*/ 
  /* 6) /cp_mass_cut{} */
  sscanf(wave_dict[6].keyarg,"%lg",&dnum);
  cp_parse->cp_mass_cut=dnum;
  index = 6;
  if((dnum <=0)){
    keyarg_barf(wave_dict,molsetname,fun_key,index);
  }
  /*-----------------------------------------------------------------------*/ 
  /* 7) /cp_mass_tau{} */
  sscanf(wave_dict[7].keyarg,"%lg",&dnum);
  cp_parse->cp_mass_tau=dnum;
  index = 7;
  if((dnum <=0)){
    keyarg_barf(wave_dict,molsetname,fun_key,index);
  }
  /*-----------------------------------------------------------------------*/ 
  /* 8) /cp_energy_cut_dual_grid{} */
  sscanf(wave_dict[8].keyarg,"%lg",&dnum);
  cp_parse->cp_ecut_dual_grid=dnum;
  index = 8;
  if((dnum <=0)){
    keyarg_barf(wave_dict,molsetname,fun_key,index);
  }
  /*-----------------------------------------------------------------------*/ 
}/*end routine*/
/*==========================================================================*/




