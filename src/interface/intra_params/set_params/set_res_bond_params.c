/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
/*                                                                          */
/*                         PI_MD:                                           */
/*             The future of simulation technology                          */
/*             ------------------------------------                         */
/*                     Module: set_mol_bond_parms.c                         */
/*                                                                          */
/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

#include "standard_include.h"
#include "../typ_defs/typedefs_par.h"
#include "../typ_defs/typedefs_class.h"
#include "../typ_defs/typedefs_bnd.h"
#include "../proto_defs/proto_intra_params_local.h"
#include "../proto_defs/proto_handle_entry.h"

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void set_res_bond_params(char *filename,char fun_key[],int ires_bond,
			 DICT_WORD res_bond_dict[],int num_res_bond_dict,
			 RESBOND_PRM resbond_prm[],
                         int nres_now,ATOMMAPS *atommaps,int ires_off,
                         int pi_beads)

/*==========================================================================*/
{
  int num,i,index,ires1,ires2,ires_typ;
  
  /*========================================================================*/
  /* I) Check for missing key words*/
  for(i=1;i<=num_res_bond_dict;i++){
    if(res_bond_dict[i].iuset==0 && res_bond_dict[i].key_type==1){
      keyword_miss(res_bond_dict,filename,fun_key,i);}
  } /*endfor*/

  /*=======================================================================*/
  /* II) Set the params */
  /*-----------------------------------------------------------------------*/ 
  /*  1)\res1_typ{} */

  strcpy(resbond_prm[ires_bond].res1_typ,res_bond_dict[1].keyarg);

  /*-----------------------------------------------------------------------*/ 
  /*  2)\res2_typ{} */

  strcpy(resbond_prm[ires_bond].res2_typ,res_bond_dict[2].keyarg);

  /*-----------------------------------------------------------------------*/ 
  /*  3)\res1_index{} */

  sscanf(res_bond_dict[3].keyarg,"%d",&ires1);
  resbond_prm[ires_bond].res1_index=ires1;
  index = 3;
  if((ires1 <=0)|| (ires1 > nres_now)){
    keyarg_barf(res_bond_dict,filename,fun_key,index);
  }

  /*-----------------------------------------------------------------------*/ 
  /*  4)\res2_index{} */

  sscanf(res_bond_dict[4].keyarg,"%d",&ires2);
  resbond_prm[ires_bond].res2_index=ires2;
  index = 4;
  if((ires2 <=0 || (ires2>nres_now))){
    keyarg_barf(res_bond_dict,filename,fun_key,index);
  }
  /*-----------------------------------------------------------------------*/ 
  /*  5)\res1_bond_site{} */

  strcpy(resbond_prm[ires_bond].res1_site,res_bond_dict[5].keyarg);
  if(strcasecmp(res_bond_dict[5].keyarg,"null")==0){
        printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
        printf("You have defined the null bond site!!!     \n");
        printf("Atoms are defaulted to the null bond site. \n");
        printf("Choose another name for your bond site and \n");
        printf("be excellent to your residue.              \n");
        printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
        fflush(stdout);
        exit(1);
  }/*endif*/

  /*-----------------------------------------------------------------------*/ 
  /*  6)\res2_bond_site{} */

  strcpy(resbond_prm[ires_bond].res2_site,res_bond_dict[6].keyarg);
  if(strcasecmp(res_bond_dict[6].keyarg,"null")==0){
        printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
        printf("You have defined the null bond site!!!    \n");
        printf("Atoms are defaulted to the null bond site. \n");
        printf("Choose another name for your bond site and \n");
        printf("be excellent to your residue.              \n");
        printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
        fflush(stdout);
        exit(1);
  }/*endif*/

  /*-----------------------------------------------------------------------*/ 
  /*  7)\res1_bondfile{} */

  strcpy(resbond_prm[ires_bond].res1_file,res_bond_dict[7].keyarg);

  /*-----------------------------------------------------------------------*/ 
  /*  8)\res2_bondfile{} */

  strcpy(resbond_prm[ires_bond].res2_file,res_bond_dict[8].keyarg);

  /*-----------------------------------------------------------------------*/ 
  /*  9)\modifier{} */

  index = 9;
  num=-1;
  if(strcasecmp(res_bond_dict[9].keyarg,"on")==0){num= 0;}
  if(strcasecmp(res_bond_dict[9].keyarg,"con")==0){num= 1;}
  if(strcasecmp(res_bond_dict[9].keyarg,"null")==0){num= 2;}
  if((num <0)){keyarg_barf(res_bond_dict,filename,fun_key,index);}
  resbond_prm[ires_bond].opt = num;
  if(num == 1 && pi_beads>1){ 
        printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
        printf("Constraints not implemented under  \n");
        printf(" path integral dynamics.            \n");
        printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
        fflush(stdout);
        exit(1);
  }/*endif*/

  /*-----------------------------------------------------------------------*/ 
  /*  10)\label{} */

  strcpy(resbond_prm[ires_bond].label,res_bond_dict[10].keyarg);

  /*=======================================================================*/
  /* III) Check the res types versus the res indicies given                */
    
  ires_typ = atommaps->ires_typ_jres_jmol_typ[ires1+ires_off];
  if(strcasecmp(resbond_prm[ires_bond].res1_typ,
                atommaps->res_typ[ires_typ])!=0){
     printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");     
     printf("Input error in ~res_bond_def\n");
     printf("1st Residue type, 1st residue index mismatch in file %s\n",
             filename);
     printf("Residue index %d is of type %s not of type %s\n",
           ires1,resbond_prm[ires_bond].res1_typ,atommaps->res_typ[ires_typ]);
     printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
     fflush(stdout);
     exit(1);
  }/*endif*/
  ires_typ = atommaps->ires_typ_jres_jmol_typ[ires2+ires_off];
  if(strcasecmp(resbond_prm[ires_bond].res2_typ,
		atommaps->res_typ[ires_typ])!=0){
     printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");     
     printf("Input error in ~res_bond_def\n");
     printf("2nd Residue type, 2nd residue index mismatch in file %s\n",
             filename);
     printf("Residue index %d is of type %s not of type %s\n",
            ires2,resbond_prm[ires_bond].res2_typ,
	    atommaps->res_typ[ires_typ]);
     printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
     fflush(stdout);
     exit(1);
  }/*endif*/

  /*-----------------------------------------------------------------------*/ 
}/*end routine*/
/*==========================================================================*/
