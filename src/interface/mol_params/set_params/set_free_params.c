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
#include "../typ_defs/typedefs_class.h"
#include "../typ_defs/typedefs_cp.h"
#include "../typ_defs/typedefs_par.h"
#include "../typ_defs/typedefs_bnd.h"
#include "../proto_defs/proto_mol_params_local.h"
#include "../proto_defs/proto_handle_entry.h"
#include "../proto_defs/proto_friend_lib_entry.h"

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void set_bond_free_params(char *molsetname,char fun_key[],
                          DICT_WORD bond_free_dict[], int num_bond_free_dict,
                          BOND_FREE *bond_free,FREE_PARSE *free_parse,
                          int nmol_typ)     

/*========================================================================*/
  {/*begin routine */
/*========================================================================*/
/* Local variables */
  int num,i,index;
  double dnum;
/*========================================================================*/
/* I) Check for missing key words*/

  for(i=1;i<num_bond_free_dict;i++){
    if(bond_free_dict[i].iuset==0 && bond_free_dict[i].key_type==1){
      keyword_miss(bond_free_dict,molsetname,fun_key,i);}
  }/*endfor*/

/*=========================================================================*/
/* II) Set the params */
/*-----------------------------------------------------------------------*/ 
  /*  1)\atom1_moltyp_ind{} */
  sscanf(bond_free_dict[1].keyarg,"%d",&num);
  free_parse->imoltyp_bond_free[1] = num; 
  index = 1;
  if((num <=0)|| (num>nmol_typ)){
    keyarg_barf(bond_free_dict,molsetname,fun_key,index);
  }
  /*-----------------------------------------------------------------------*/ 
  /*  2)\atom2_moltyp_ind{} */
  sscanf(bond_free_dict[2].keyarg,"%d",&num);
  free_parse->imoltyp_bond_free[2] = num; 
  index = 2;
  if((num <=0)|| (num>nmol_typ)){
    keyarg_barf(bond_free_dict,molsetname,fun_key,index);
  }
  /*-----------------------------------------------------------------------*/ 
  /*  3)\atom1_mol_ind{}    */
  sscanf(bond_free_dict[3].keyarg,"%d",&num);
  free_parse->imol_bond_free[1] = num; 
  index = 3;
  if((num <=0)){keyarg_barf(bond_free_dict,molsetname,fun_key,index);}

  /*-----------------------------------------------------------------------*/ 
  /*  4)\atom2_mol_ind{}    */
  sscanf(bond_free_dict[4].keyarg,"%d",&num);
  free_parse->imol_bond_free[2] = num; 
  index = 4;
  if((num <=0)){keyarg_barf(bond_free_dict,molsetname,fun_key,index);}

  /*-----------------------------------------------------------------------*/ 
  /*  5)\atom1_atm_ind{}    */
  sscanf(bond_free_dict[5].keyarg,"%d",&num);
  free_parse->iatm_bond_free[1] = num; 
  index = 5;
  if((num <=0)){keyarg_barf(bond_free_dict,molsetname,fun_key,index);}

  /*-----------------------------------------------------------------------*/ 
  /*  6)\atom2_atm_ind{}    */
  sscanf(bond_free_dict[6].keyarg,"%d",&num);
  free_parse->iatm_bond_free[2] = num; 
  index = 6;
  if((num <=0)){keyarg_barf(bond_free_dict,molsetname,fun_key,index);}

  /*-----------------------------------------------------------------------*/ 
  /* 7)\eq{}               */
  sscanf(bond_free_dict[7].keyarg,"%lg",&dnum);
  bond_free->eq = dnum/BOHR;
  index = 7;
  /*-----------------------------------------------------------------------*/ 
  /* 8)\fk{}               */
  sscanf(bond_free_dict[8].keyarg,"%lg",&dnum);
  bond_free->fk = dnum*BOHR*BOHR/BOLTZ;
  index = 8;
  if((dnum < 0)){keyarg_barf(bond_free_dict,molsetname,fun_key,index);}
  /*-----------------------------------------------------------------------*/ 
  /* 9)\rmin_hist{}         */
  sscanf(bond_free_dict[9].keyarg,"%lg",&dnum);
  bond_free->rmin = dnum/BOHR;
  /*-----------------------------------------------------------------------*/ 
  /* 10)\rmax_hist{}         */
  sscanf(bond_free_dict[10].keyarg,"%lg",&dnum);
  bond_free->rmax = dnum/BOHR;
  /*-----------------------------------------------------------------------*/ 
  /* 11)\num_hist{}         */
  sscanf(bond_free_dict[11].keyarg,"%d",&num);
  bond_free->nhist = num;
  index = 11;
  if((num <=0)){keyarg_barf(bond_free_dict,molsetname,fun_key,index);}
  /*-----------------------------------------------------------------------*/ 
  /* 12)\power{}            */
  sscanf(bond_free_dict[12].keyarg,"%d",&num);
  bond_free->npow = num;
  index = 12;
  if((num <=0)){keyarg_barf(bond_free_dict,molsetname,fun_key,index);}
  /*-----------------------------------------------------------------------*/ 
  /* 13)\hist_file{}        */
  strcpy(bond_free->file,bond_free_dict[13].keyarg);
  /*-----------------------------------------------------------------------*/
  /*  14)\atom1_residue_ind{} */
  sscanf(bond_free_dict[14].keyarg,"%d",&num);
  free_parse->ires_bond_free[1] = num; 
  index = 14;
  if((num <=0)){
    keyarg_barf(bond_free_dict,molsetname,fun_key,index);
  }
  /*-----------------------------------------------------------------------*/ 
  /*  15)\atom2_residue_ind{} */
  sscanf(bond_free_dict[15].keyarg,"%d",&num);
  free_parse->ires_bond_free[2] = num; 
  index = 15;
  if((num <=0)){
    keyarg_barf(bond_free_dict,molsetname,fun_key,index);
  }
  /*-----------------------------------------------------------------------*/

  bond_free->del =((bond_free->rmax)- (bond_free->rmin))/
                       ((double)(bond_free->nhist));
/*-----------------------------------------------------------------------*/
/*end routine*/}
/*==========================================================================*/





/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void set_bend_free_params(char *molsetname,char fun_key[],
                          DICT_WORD bend_free_dict[], int num_bend_free_dict,
                          BEND_FREE *bend_free,FREE_PARSE *free_parse,
                          int nmol_typ)

/*==========================================================================*/
{/*begin routine*/
/*==========================================================================*/
/* Local variables*/
  int i,num,index;
  double dnum;
  /*=======================================================================*/
  /* I) Check for missing key words*/
  for(i=1;i<num_bend_free_dict;i++){
    if(bend_free_dict[i].iuset==0 && bend_free_dict[i].key_type==1){
      keyword_miss(bend_free_dict,molsetname,fun_key,i);}
  } /*endfor*/

  /*=======================================================================*/
  /* II) Set the params */
  /*-----------------------------------------------------------------------*/ 
  /*  1)\atom1_moltyp_ind{} */
  sscanf(bend_free_dict[1].keyarg,"%d",&num);
  free_parse->imoltyp_bend_free[1] = num; 
  index = 1;
  if((num <=0)|| (num>nmol_typ)){
    keyarg_barf(bend_free_dict,molsetname,fun_key,index);
  }

  /*-----------------------------------------------------------------------*/ 
  /*  2)\atom2_moltyp_ind{} */
  sscanf(bend_free_dict[2].keyarg,"%d",&num);
  free_parse->imoltyp_bend_free[2] = num; 
  index = 2;
  if((num <=0)|| (num>nmol_typ)){
    keyarg_barf(bend_free_dict,molsetname,fun_key,index);
  }
  /*-----------------------------------------------------------------------*/ 
  /*  3)\atom3_moltyp_ind{} */
  sscanf(bend_free_dict[3].keyarg,"%d",&num);
  free_parse->imoltyp_bend_free[3] = num; 
  index = 3;
  if((num <=0)|| (num>nmol_typ)){
    keyarg_barf(bend_free_dict,molsetname,fun_key,index);
  }
  /*-----------------------------------------------------------------------*/ 
  /*  4)\atom1_mol_ind{}    */
  sscanf(bend_free_dict[4].keyarg,"%d",&num);
  free_parse->imol_bend_free[1] = num; 
  index = 4;
  if((num <=0)){
    keyarg_barf(bend_free_dict,molsetname,fun_key,index);
  }
  /*-----------------------------------------------------------------------*/ 
  /*  5)\atom2_mol_ind{}    */
  sscanf(bend_free_dict[5].keyarg,"%d",&num);
  free_parse->imol_bend_free[2] = num; 
  index = 5;
  if((num <=0)){
    keyarg_barf(bend_free_dict,molsetname,fun_key,index);
  }
  /*-----------------------------------------------------------------------*/ 
  /*  6)\atom3_mol_ind{}    */
  sscanf(bend_free_dict[6].keyarg,"%d",&num);
  free_parse->imol_bend_free[3] = num; 
  index = 6;
  if((num <=0)){
    keyarg_barf(bend_free_dict,molsetname,fun_key,index);
  }
  /*-----------------------------------------------------------------------*/ 
  /*  7)\atom1_atm_ind{}    */
  sscanf(bend_free_dict[7].keyarg,"%d",&num);
  free_parse->iatm_bend_free[1] = num; 
  index = 7;
  if((num <=0)){
    keyarg_barf(bend_free_dict,molsetname,fun_key,index);
  }
  /*-----------------------------------------------------------------------*/ 
  /*  8)\atom2_atm_ind{}    */
  sscanf(bend_free_dict[8].keyarg,"%d",&num);
  free_parse->iatm_bend_free[2] = num; 
  index = 8;
  if((num <=0)){
    keyarg_barf(bend_free_dict,molsetname,fun_key,index);
  }
  /*-----------------------------------------------------------------------*/ 
  /*  9)\atom3_atm_ind{}    */
  sscanf(bend_free_dict[9].keyarg,"%d",&num);
  free_parse->iatm_bend_free[3] = num; 
  index = 9;
  if((num <=0)){
    keyarg_barf(bend_free_dict,molsetname,fun_key,index);
  }
  /*-----------------------------------------------------------------------*/ 
  /* 10)\eq{}               */
  sscanf(bend_free_dict[10].keyarg,"%lg",&dnum);
  bend_free->eq = dnum*M_PI/180.0;
  index = 10;
  if((dnum <0.0)||(dnum >180.0)){
    keyarg_barf(bend_free_dict,molsetname,fun_key,index);
  }
  /*-----------------------------------------------------------------------*/ 
  /* 11)\fk{}               */
  sscanf(bend_free_dict[11].keyarg,"%lg",&dnum);
  bend_free->fk = dnum/BOLTZ;
  index = 11;
  if((dnum < 0.)){
    keyarg_barf(bend_free_dict,molsetname,fun_key,index);}
  /*-----------------------------------------------------------------------*/ 
  /* 12)\num_hist{}         */
  sscanf(bend_free_dict[12].keyarg,"%d",&num);
  bend_free->nhist = num;
  index = 12;
  if(num <=0){
    keyarg_barf(bend_free_dict,molsetname,fun_key,index);}
  /*-----------------------------------------------------------------------*/ 
  /* 13)\power{}            */
  sscanf(bend_free_dict[13].keyarg,"%d",&num);
  bend_free->npow = num;
  if(num <=0){
    keyarg_barf(bend_free_dict,molsetname,fun_key,index);}
  /*-----------------------------------------------------------------------*/ 
  /* 14)\hist_file{}        */
  strcpy(bend_free->file,bend_free_dict[14].keyarg);
  /*-----------------------------------------------------------------------*/
  /*  15)\atom1_residue_ind{} */
  sscanf(bend_free_dict[15].keyarg,"%d",&num);
  free_parse->ires_bend_free[1] = num; 
  index = 15;
  if(num <=0){
    keyarg_barf(bend_free_dict,molsetname,fun_key,index);
  }
  /*-----------------------------------------------------------------------*/ 
  /*  16)\atom2_residue_ind{} */
  sscanf(bend_free_dict[16].keyarg,"%d",&num);
  free_parse->ires_bend_free[2] = num; 
  index = 16;
  if(num <=0){
    keyarg_barf(bend_free_dict,molsetname,fun_key,index);
  }
  /*-----------------------------------------------------------------------*/ 
  /*  17)\atom3_residue_ind{} */
  sscanf(bend_free_dict[17].keyarg,"%d",&num);
  free_parse->ires_bend_free[3] = num; 
  index = 17;
  if(num <=0){
    keyarg_barf(bend_free_dict,molsetname,fun_key,index);
  }
  /*-----------------------------------------------------------------------*/
  bend_free->del = M_PI/((double)(bend_free->nhist));
  /*=======================================================================*/
}  /*end routine*/
/*==========================================================================*/





/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void set_tors_free_params(char molsetname[],char fun_key[],
                          DICT_WORD tors_free_dict[], int num_tors_free_dict,
                          TORS_FREE *tors_free,FREE_PARSE *free_parse,
                          int nmol_typ)

/*==========================================================================*/
{/*begin routine*/
/*==========================================================================*/
/* Local variables */
  int i,num,index,ntors;
  double dnum;
/*==========================================================================*/
  /* I) Check for missing key words*/
  for(i=1;i<=num_tors_free_dict;i++){
    if(tors_free_dict[i].iuset==0 && tors_free_dict[i].key_type==1){
      keyword_miss(tors_free_dict,molsetname,fun_key,i);
    }
  }    /*endfor*/

 /*==========================================================================*/
  /* II) Set the params */

  /*-----------------------------------------------------------------------*/ 
  /* 39)\ntors{} */
  index = 39;
  sscanf(tors_free_dict[index].keyarg,"%d",&num);
  ntors = num;
  tors_free->num = num; 
  if(num <=0 || num > 2){
    keyarg_barf(tors_free_dict,molsetname,fun_key,index);
  }

  if(ntors>1){

    for(i=22;i<=33;i++){
      if(tors_free_dict[i].iuset==0){
        keyword_miss(tors_free_dict,molsetname,fun_key,i);
      }
    }/*endfor*/
    if(tors_free_dict[38].iuset==0){
      keyword_miss(tors_free_dict,molsetname,fun_key,i);
    }/*endif*/

  }/*endif*/


  /*-----------------------------------------------------------------------*/ 
  /*  1)\atom1.1_moltyp_ind{} */
  sscanf(tors_free_dict[1].keyarg,"%d",&num);
  free_parse->imoltyp_tors_free[1] = num; 
  index = 1;
  if((num <=0)|| (num>nmol_typ)){
    keyarg_barf(tors_free_dict,molsetname,fun_key,index);
  }
  /*-----------------------------------------------------------------------*/ 
  /*  2)\atom2.1_moltyp_ind{} */
  sscanf(tors_free_dict[2].keyarg,"%d",&num);
  free_parse->imoltyp_tors_free[2] = num; 
  index = 2;
  if((num <=0)|| (num>nmol_typ)){
    keyarg_barf(tors_free_dict,molsetname,fun_key,index);
  }
  /*-----------------------------------------------------------------------*/ 
  /*  3)\atom3.1_moltyp_ind{} */
  sscanf(tors_free_dict[3].keyarg,"%d",&num);
  free_parse->imoltyp_tors_free[3] = num; 
  index = 3;
  if((num <=0)|| (num>nmol_typ)){
    keyarg_barf(tors_free_dict,molsetname,fun_key,index);
  }
  /*-----------------------------------------------------------------------*/ 
  /*  4)\atom4.1_moltyp_ind{} */
  sscanf(tors_free_dict[4].keyarg,"%d",&num);
  free_parse->imoltyp_tors_free[4] = num; 
  index = 4;
  if((num <=0)|| (num>nmol_typ)){
    keyarg_barf(tors_free_dict,molsetname,fun_key,index);
  }
  /*-----------------------------------------------------------------------*/ 
  /*  5)\atom1.1_mol_ind{}    */
  sscanf(tors_free_dict[5].keyarg,"%d",&num);
  free_parse->imol_tors_free[1] = num; 
  index = 5;
  if((num <=0)){
    keyarg_barf(tors_free_dict,molsetname,fun_key,index);
  }
  /*-----------------------------------------------------------------------*/ 
  /*  6)\atom2.1_mol_ind{}    */
  sscanf(tors_free_dict[6].keyarg,"%d",&num);
  free_parse->imol_tors_free[2] = num; 
  index = 6;
  if((num <=0)){
    keyarg_barf(tors_free_dict,molsetname,fun_key,index);
  }
  /*-----------------------------------------------------------------------*/ 
  /*  7)\atom3.1_mol_ind{}    */
  sscanf(tors_free_dict[7].keyarg,"%d",&num);
  free_parse->imol_tors_free[3] = num; 
  index = 7;
  if((num <=0)){
    keyarg_barf(tors_free_dict,molsetname,fun_key,index);
  }
  /*-----------------------------------------------------------------------*/ 
  /*  8)\atom4.1_mol_ind{}    */
  sscanf(tors_free_dict[8].keyarg,"%d",&num);
  free_parse->imol_tors_free[4] = num; 
  index = 8;
  if((num <=0)){
    keyarg_barf(tors_free_dict,molsetname,fun_key,index);
  }
  /*-----------------------------------------------------------------------*/ 
  /*  9)\atom1.1_atm_ind{}    */
  sscanf(tors_free_dict[9].keyarg,"%d",&num);
  free_parse->iatm_tors_free[1] = num; 
  index = 9;
  if((num <=0)){
    keyarg_barf(tors_free_dict,molsetname,fun_key,index);
  }
  /*-----------------------------------------------------------------------*/ 
  /* 10)\atom2.1_atm_ind{}    */
  sscanf(tors_free_dict[10].keyarg,"%d",&num);
  free_parse->iatm_tors_free[2] = num; 
  index = 10;
  if((num <=0)){
    keyarg_barf(tors_free_dict,molsetname,fun_key,index);
  }
  /*-----------------------------------------------------------------------*/ 
  /* 11)\atom3.1_atm_ind{}    */
  sscanf(tors_free_dict[11].keyarg,"%d",&num);
  free_parse->iatm_tors_free[3] = num; 
  index = 11;
  if((num <=0)){
    keyarg_barf(tors_free_dict,molsetname,fun_key,index);
  }
  /*-----------------------------------------------------------------------*/ 
  /* 12)\atom4.1_atm_ind{}    */
  sscanf(tors_free_dict[12].keyarg,"%d",&num);
  free_parse->iatm_tors_free[4] = num; 
  index = 12;
  if((num <=0)){
    keyarg_barf(tors_free_dict,molsetname,fun_key,index);
  }
  /*-----------------------------------------------------------------------*/ 
  /* 13)\eq_1{}               */
  sscanf(tors_free_dict[13].keyarg,"%lg",&dnum);
  tors_free->eq[1] = dnum*M_PI/180.0;
  index = 13;
  if((dnum <-180.0)||(dnum >180.0)){
    keyarg_barf(tors_free_dict,molsetname,fun_key,index);
  }
  /*-----------------------------------------------------------------------*/ 
  /* 14)\fk{}               */
  sscanf(tors_free_dict[14].keyarg,"%lg",&dnum);
  tors_free->fk = dnum/BOLTZ;
  index = 14;
  if((dnum < 0.)){
    keyarg_barf(tors_free_dict,molsetname,fun_key,index);
  }
  /*-----------------------------------------------------------------------*/ 
  /* 15)\num_hist{}         */
  sscanf(tors_free_dict[15].keyarg,"%d",&num);
  tors_free->nhist = num;
  index = 15;
  if((num <=0)){
    keyarg_barf(tors_free_dict,molsetname,fun_key,index);
  }
  /*-----------------------------------------------------------------------*/ 
  /* 16)\power{}            */
  sscanf(tors_free_dict[16].keyarg,"%d",&num);
  tors_free->npow = num;
  index = 16;
  if((num <=0)){
    keyarg_barf(tors_free_dict,molsetname,fun_key,index);
  }
  /*-----------------------------------------------------------------------*/ 
  /* 17)\hist_file{}        */
  strcpy(tors_free->file,tors_free_dict[17].keyarg);

  /*-----------------------------------------------------------------------*/
  /*  18)\atom1.1_residue_ind{} */
  sscanf(tors_free_dict[18].keyarg,"%d",&num);
  free_parse->ires_tors_free[1] = num; 
  index = 18;
  if((num <=0)){
    keyarg_barf(tors_free_dict,molsetname,fun_key,index);
  }
  /*-----------------------------------------------------------------------*/ 
  /*  19)\atom2.1_iresidue_ind{} */
  sscanf(tors_free_dict[19].keyarg,"%d",&num);
  free_parse->ires_tors_free[2] = num; 
  index = 19;
  if((num <=0)){
    keyarg_barf(tors_free_dict,molsetname,fun_key,index);
  }
  /*-----------------------------------------------------------------------*/ 
  /*  20)\atom3.1_iresidue_ind{} */
  sscanf(tors_free_dict[20].keyarg,"%d",&num);
  free_parse->ires_tors_free[3] = num; 
  index = 20;
  if((num <=0)){
    keyarg_barf(tors_free_dict,molsetname,fun_key,index);
  }
  /*-----------------------------------------------------------------------*/ 
  /*  21)\atom4.1_iresidue_ind{} */
  sscanf(tors_free_dict[21].keyarg,"%d",&num);
  free_parse->ires_tors_free[4] = num; 
  index = 21;
  if(num <=0){
    keyarg_barf(tors_free_dict,molsetname,fun_key,index);
  }
  /*-----------------------------------------------------------------------*/
  /*-----------------------------------------------------------------------*/ 
  /*  22)\atom1.2_moltyp_ind{} */
  index = 22;
  sscanf(tors_free_dict[index].keyarg,"%d",&num);
  free_parse->imoltyp_tors_free[5] = num; 
  if((num <=0)|| (num>nmol_typ)){
    keyarg_barf(tors_free_dict,molsetname,fun_key,index);
  }
  /*-----------------------------------------------------------------------*/ 
  /* 23)\atom2.2_moltyp_ind{} */
  index = 23;
  sscanf(tors_free_dict[index].keyarg,"%d",&num);
  free_parse->imoltyp_tors_free[6] = num; 
  if((num <=0)|| (num>nmol_typ)){
    keyarg_barf(tors_free_dict,molsetname,fun_key,index);
  }
  /*-----------------------------------------------------------------------*/ 
  /* 24)\atom3.2_moltyp_ind{} */
  index = 24;
  sscanf(tors_free_dict[index].keyarg,"%d",&num);
  free_parse->imoltyp_tors_free[7] = num; 
  if((num <=0)|| (num>nmol_typ)){
    keyarg_barf(tors_free_dict,molsetname,fun_key,index);
  }
  /*-----------------------------------------------------------------------*/ 
  /* 25)\atom4.2_moltyp_ind{} */
  index = 25;
  sscanf(tors_free_dict[index].keyarg,"%d",&num);
  free_parse->imoltyp_tors_free[8] = num; 
  if((num <=0)|| (num>nmol_typ)){
    keyarg_barf(tors_free_dict,molsetname,fun_key,index);
  }
  /*-----------------------------------------------------------------------*/ 
  /* 26)\atom1.2_mol_ind{}    */
  index = 26;
  sscanf(tors_free_dict[index].keyarg,"%d",&num);
  free_parse->imol_tors_free[5] = num; 
  if((num <=0)){
    keyarg_barf(tors_free_dict,molsetname,fun_key,index);
  }
  /*-----------------------------------------------------------------------*/ 
  /* 27)\atom2.2_mol_ind{}    */
  index = 27;
  sscanf(tors_free_dict[index].keyarg,"%d",&num);
  free_parse->imol_tors_free[6] = num; 
  if((num <=0)){
    keyarg_barf(tors_free_dict,molsetname,fun_key,index);
  }
  /*-----------------------------------------------------------------------*/ 
  /* 28)\atom3.2_mol_ind{}    */
  index = 28;
  sscanf(tors_free_dict[index].keyarg,"%d",&num);
  free_parse->imol_tors_free[7] = num; 
  if((num <=0)){
    keyarg_barf(tors_free_dict,molsetname,fun_key,index);
  }
  /*-----------------------------------------------------------------------*/ 
  /* 29)\atom4.2_mol_ind{}    */
  index = 29;
  sscanf(tors_free_dict[index].keyarg,"%d",&num);
  free_parse->imol_tors_free[8] = num; 
  if((num <=0)){
    keyarg_barf(tors_free_dict,molsetname,fun_key,index);
  }
  /*-----------------------------------------------------------------------*/ 
  /* 30)\atom1.2_atm_ind{}    */
  index = 30;
  sscanf(tors_free_dict[index].keyarg,"%d",&num);
  free_parse->iatm_tors_free[5] = num; 
  if((num <=0)){
    keyarg_barf(tors_free_dict,molsetname,fun_key,index);
  }
  /*-----------------------------------------------------------------------*/ 
  /* 31)\atom2.2_atm_ind{}    */
  index = 31;
  sscanf(tors_free_dict[index].keyarg,"%d",&num);
  free_parse->iatm_tors_free[6] = num; 
  if((num <=0)){
    keyarg_barf(tors_free_dict,molsetname,fun_key,index);
  }
  /*-----------------------------------------------------------------------*/ 
  /* 32)\atom3.2_atm_ind{}    */
  index = 32;
  sscanf(tors_free_dict[index].keyarg,"%d",&num);
  free_parse->iatm_tors_free[7] = num; 
  if((num <=0)){
    keyarg_barf(tors_free_dict,molsetname,fun_key,index);
  }
  /*-----------------------------------------------------------------------*/ 
  /* 33)\atom4.2_atm_ind{}    */
  index = 33;
  sscanf(tors_free_dict[index].keyarg,"%d",&num);
  free_parse->iatm_tors_free[8] = num; 
  if((num <=0)){
    keyarg_barf(tors_free_dict,molsetname,fun_key,index);
  }
  /*-----------------------------------------------------------------------*/
  /*  34)\atom1.2_residue_ind{} */
  index = 34;
  sscanf(tors_free_dict[index].keyarg,"%d",&num);
  free_parse->ires_tors_free[5] = num; 
  if((num <=0)){
    keyarg_barf(tors_free_dict,molsetname,fun_key,index);
  }
  /*-----------------------------------------------------------------------*/ 
  /* 35)\atom2.2_iresidue_ind{} */
  index = 35;
  sscanf(tors_free_dict[index].keyarg,"%d",&num);
  free_parse->ires_tors_free[6] = num; 
  if((num <=0)){
    keyarg_barf(tors_free_dict,molsetname,fun_key,index);
  }
  /*-----------------------------------------------------------------------*/ 
  /* 36)\atom3.2_iresidue_ind{} */
  index = 36;
  sscanf(tors_free_dict[index].keyarg,"%d",&num);
  free_parse->ires_tors_free[7] = num; 
  if((num <=0)){
    keyarg_barf(tors_free_dict,molsetname,fun_key,index);
  }
  /*-----------------------------------------------------------------------*/ 
  /* 37)\atom4.2_iresidue_ind{} */
  index = 37;
  sscanf(tors_free_dict[index].keyarg,"%d",&num);
  free_parse->ires_tors_free[8] = num; 
  if(num <=0){
    keyarg_barf(tors_free_dict,molsetname,fun_key,index);
  }
  /*-----------------------------------------------------------------------*/ 
  /* 38)\eq_2{}               */
  index = 38;
  sscanf(tors_free_dict[index].keyarg,"%lg",&dnum);
  tors_free->eq[2] = dnum*M_PI/180.0;
  if((dnum <-180.0)||(dnum >180.0)){
    keyarg_barf(tors_free_dict,molsetname,fun_key,index);
  }
  /*-----------------------------------------------------------------------*/ 

  tors_free->del = (2.0*M_PI)/((double)(tors_free->nhist));

  /*-----------------------------------------------------------------------*/
   }  /*end routine*/
 /*========================================================================*/



/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void set_rbar_free_params(char *molsetname,char fun_key[],
                          DICT_WORD rbar_free_dict[], int num_rbar_free_dict,
                          RBAR_SIG_FREE *rbar_sig_free,FREE_PARSE *free_parse,
                          int nmol_typ)     

/*========================================================================*/
  {/*begin routine */
/*========================================================================*/
/* Local variables */

  int num,i,index;
  double dnum;
  FILE *fp;

/*========================================================================*/
/* I) Check for missing key words*/

  for(i=1;i<num_rbar_free_dict;i++){
    if((rbar_free_dict[i].iuset==0) && (rbar_free_dict[i].key_type==1)){
      keyword_miss(rbar_free_dict,molsetname,fun_key,i);
    }/*endif*/
  }/*endfor*/

/*=========================================================================*/
/* II) Set the params */

  /*-----------------------------------------------------------------------*/ 
  /* 1)\num_bond{}     */
   sscanf(rbar_free_dict[1].keyarg,"%d",&num);
   rbar_sig_free->nfree  = num;
   free_parse->nbar_bond = num;
   index = 1;
   if((num <=0)){keyarg_barf(rbar_free_dict,molsetname,fun_key,index);}
  /*-----------------------------------------------------------------------*/ 
  /* 3)\eq_bar{}               */
   sscanf(rbar_free_dict[3].keyarg,"%lg",&dnum);
   rbar_sig_free->eq_bar = dnum/BOHR;
   index = 3;
   if((dnum < 0)){keyarg_barf(rbar_free_dict,molsetname,fun_key,index);}
  /*-----------------------------------------------------------------------*/ 
  /* 4)\fk_bar{}               */
   sscanf(rbar_free_dict[4].keyarg,"%lg",&dnum);
   rbar_sig_free->fk_bar = dnum*BOHR*BOHR/BOLTZ;
   index = 4;
   if((dnum < 0)){keyarg_barf(rbar_free_dict,molsetname,fun_key,index);}
  /*-----------------------------------------------------------------------*/ 
  /* 5)\eq_sig{}               */
   sscanf(rbar_free_dict[5].keyarg,"%lg",&dnum);
   rbar_sig_free->eq_sigma = dnum/BOHR;
   index = 5;
   if((dnum < 0)){keyarg_barf(rbar_free_dict,molsetname,fun_key,index);}
  /*-----------------------------------------------------------------------*/ 
  /* 6)\fk_sig{}               */
   sscanf(rbar_free_dict[6].keyarg,"%lg",&dnum);
   rbar_sig_free->fk_sigma = dnum*BOHR*BOHR/BOLTZ;
   index = 6;
   if((dnum < 0)){keyarg_barf(rbar_free_dict,molsetname,fun_key,index);}
  /*-----------------------------------------------------------------------*/ 
  /* 7)\rmin_hist{}         */
   sscanf(rbar_free_dict[7].keyarg,"%lg",&dnum);
   rbar_sig_free->rmin = dnum/BOHR;
   index = 7;
   if((dnum < 0)){keyarg_barf(rbar_free_dict,molsetname,fun_key,index);}
  /*-----------------------------------------------------------------------*/ 
  /* 8)\rmax_hist{}         */
   sscanf(rbar_free_dict[8].keyarg,"%lg",&dnum);
   rbar_sig_free->rmax = dnum/BOHR;
   index = 8;
   if((dnum < 0)|| (dnum< rbar_sig_free->rmin)){
                       keyarg_barf(rbar_free_dict,molsetname,fun_key,index);}
  /*-----------------------------------------------------------------------*/ 
  /* 9)\num_rbar_hist{}         */
   sscanf(rbar_free_dict[9].keyarg,"%d",&num);
   rbar_sig_free->nhist_bar = num;
   index = 9;
   if((num <=0)){keyarg_barf(rbar_free_dict,molsetname,fun_key,index);}
  /*-----------------------------------------------------------------------*/ 
  /* 10)\smin_hist{}         */
   sscanf(rbar_free_dict[10].keyarg,"%lg",&dnum);
   rbar_sig_free->smin = dnum/BOHR;
   index = 10;
   if((dnum < 0)){keyarg_barf(rbar_free_dict,molsetname,fun_key,index);}
  /*-----------------------------------------------------------------------*/ 
  /* 11)\smax_hist{}         */
   sscanf(rbar_free_dict[11].keyarg,"%lg",&dnum);
   rbar_sig_free->smax = dnum/BOHR;
   index = 11;
   if((dnum < 0)|| (dnum< rbar_sig_free->smin)){
                       keyarg_barf(rbar_free_dict,molsetname,fun_key,index);}
  /*-----------------------------------------------------------------------*/ 
  /* 12)\num_sig_hist{}         */
   sscanf(rbar_free_dict[12].keyarg,"%d",&num);
   rbar_sig_free->nhist_sig = num;
   index = 12;
   if((num <=0)){keyarg_barf(rbar_free_dict,molsetname,fun_key,index);}
  /*-----------------------------------------------------------------------*/ 
  /* 13)\hist_file{}        */
   strcpy(rbar_sig_free->file,rbar_free_dict[13].keyarg);
  /*-----------------------------------------------------------------------*/

/*=========================================================================*/
/* III) Use the data and malloc                                            */

  rbar_sig_free->del_bar = ((rbar_sig_free->rmax)- (rbar_sig_free->rmin))/
                           ((double)(rbar_sig_free->nhist_bar));
  rbar_sig_free->del_sig = ((rbar_sig_free->smax)- (rbar_sig_free->smin))/
                           ((double)(rbar_sig_free->nhist_sig));

  num                    = rbar_sig_free->nfree;
  rbar_sig_free->rnfree  = (double) num;

  free_parse->imoltyp_rbar1_free = (int *)cmalloc(num*sizeof(int))-1;
  free_parse->imoltyp_rbar2_free = (int *)cmalloc(num*sizeof(int))-1;
  free_parse->imol_rbar1_free    = (int *)cmalloc(num*sizeof(int))-1;
  free_parse->imol_rbar2_free    = (int *)cmalloc(num*sizeof(int))-1;
  free_parse->ires_rbar1_free    = (int *)cmalloc(num*sizeof(int))-1;
  free_parse->ires_rbar2_free    = (int *)cmalloc(num*sizeof(int))-1;
  free_parse->iatm_rbar1_free    = (int *)cmalloc(num*sizeof(int))-1;
  free_parse->iatm_rbar2_free    = (int *)cmalloc(num*sizeof(int))-1;

/*=========================================================================*/
/* IV) Read in the indicies                                                */

  fp = cfopen(rbar_free_dict[2].keyarg,"r");

  fscanf(fp,"%d",&num);readtoendofline(fp);
  if(num!=rbar_sig_free->nfree){
     printf("To have no errors,\n");
     printf("Would be life without meaning,\n");
     printf("No struggle, no joy.\n");
     printf("Number of free energy bonds %d in file %s,\n",num,
                                             rbar_free_dict[2].keyarg);
     printf("does not match the request value %d \n",rbar_sig_free->nfree);
     fflush(stdout);
     exit(1);
  }/*endif*/

  for(i=1;i<=num;i++){
    fscanf(fp,"%d %d %d %d",
               &(free_parse->imoltyp_rbar1_free[i]),
               &(free_parse->imol_rbar1_free[i]),
               &(free_parse->ires_rbar1_free[i]),
               &(free_parse->iatm_rbar1_free[i]));readtoendofline(fp);
    fscanf(fp,"%d %d %d %d",
               &(free_parse->imoltyp_rbar2_free[i]),
               &(free_parse->imol_rbar2_free[i]),
               &(free_parse->ires_rbar2_free[i]),
               &(free_parse->iatm_rbar2_free[i]));readtoendofline(fp);
  }/*endfor*/

  fclose(fp);

/*-----------------------------------------------------------------------*/
   }/*end routine*/
/*==========================================================================*/
