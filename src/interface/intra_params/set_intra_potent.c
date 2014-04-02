/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
/*                                                                          */
/*                         PI_MD:                                           */
/*             The future of simulation technology                          */
/*             ------------------------------------                         */
/*                     Module: set_intra_potent.c                           */
/*                                                                          */
/* This subprogram set molecular and intramolecular data sets               */
/*                                                                          */
/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

#include "standard_include.h"
#include "../typ_defs/typedefs_class.h"
#include "../typ_defs/typedefs_par.h"
#include "../typ_defs/typedefs_bnd.h"
#include "../proto_defs/proto_intra_params_local.h"
#include "../proto_defs/proto_search_entry.h"
#include "../proto_defs/proto_friend_lib_entry.h"

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
/*  control_set_intra_potent:                                               */
/*==========================================================================*/

void set_intra_potent(BONDED *bonded,BUILD_INTRA *build_intra,
                      NAME def_intra_name[],NAME user_intra_name[])

/*========================================================================*/
    { /* begin routine */
/*========================================================================*/
  /*            Local Variables                                             */
  int i,iii,ifirst,isum,n;
  int num_fun_dict;
  DICT_WORD *fun_dict;
  /*=======================================================================*/
  /*-----------------------------------------------------------------------*/
  /*=======================================================================*/
  /*     Control set_intra_potent:                                         */
  /*=======================================================================*/
  /* 0) Set up                                                             */
  PRINT_LINE_STAR;
  printf("Searching the data bases\n");
  printf("(both user defined and default)\n");
  printf("for the intramolecular potential parameters\n");
  PRINT_LINE_DASH;printf("\n");

  ifirst =1;
  set_potfun_dict(&fun_dict,&num_fun_dict,ifirst);

  /*=======================================================================*/
  /* I) set the bond parameters                                            */
  printf("Setting the bonds\n");
  isum = bonded->bond.ntyp_pow
       + bonded->bond.ntyp_con
       + bonded->grp_bond_con.ntyp_21
       + bonded->grp_bond_con.ntyp_23
       + bonded->grp_bond_con.ntyp_33
       + bonded->grp_bond_watts.ntyp_33
       + bonded->grp_bond_con.ntyp_43
       + bonded->grp_bond_con.ntyp_46;
  if(isum > 0){
    set_bond_potent(&(bonded->bond),&(bonded->grp_bond_con),
                        &(bonded->grp_bond_watts),build_intra,
                        def_intra_name[1],user_intra_name[1],
                        fun_dict,num_fun_dict);
  }/*endif*/

  /*=======================================================================*/
  /* II) set the bend parameters                                           */
  printf("Setting the bends\n");
  if((bonded->bend.ntyp_con) > 0){
    set_bend_potent(&(bonded->bend),build_intra,
                        def_intra_name[2],user_intra_name[2],
                        fun_dict,num_fun_dict);
  }/*endif*/

  /*=======================================================================*/
  /* IV) set the onfo parameters                                            */
  printf("Setting the onefours\n");
  if((bonded->onfo.ntyp) > 0){
    set_onfo_potent(&(bonded->onfo),build_intra,
                        def_intra_name[4],user_intra_name[4],
                        fun_dict,num_fun_dict);
  }/*endif*/

  /*=======================================================================*/
  /* V) set the bend_bnd parameters                                        */
  printf("Setting the Urey-Bradleys\n");
  if((bonded->bend_bnd.ntyp) > 0){
    n = bonded->bend_bnd.ntyp;
    build_intra->ibend_bnd_typ_pure = (int *)cmalloc(n*sizeof(int))-1;
    build_intra->ibend_bnd_typ_map = (int *)cmalloc(n*sizeof(int))-1;
    set_bend_bnd_potent(&(bonded->bend_bnd),build_intra,
                        def_intra_name[2],user_intra_name[2],
                        fun_dict,num_fun_dict);
    printf("Extracting the pure bends from the Urey-Bradley list\n");
    extract_pure_bends(&(bonded->bend_bnd),&(bonded->bend),build_intra);
    cfree(&(build_intra->ibend_bnd_typ_pure[1]));
    cfree(&(build_intra->ibend_bnd_typ_map[1]));
  }/*endif*/

  /*=======================================================================*/
  /* III) set the tors parameters                                          */
  printf("Setting the torsions\n");
  if((bonded->tors.ntyp_pow+bonded->tors.ntyp_con) > 0){
    set_tors_potent(&(bonded->tors),build_intra,
                        def_intra_name[3],user_intra_name[3],
                        fun_dict,num_fun_dict);
  }/*endif*/

/*=======================================================================*/
  cfree(&fun_dict[1]);

  PRINT_LINE_DASH;
  printf("Completed the search of the data bases ");
  printf("(both user defined and default)\n");
  printf("for the intramolecular potential parameters\n");
/*  printf("enter an integer: "); scanf("%d",&iii);*/
  PRINT_LINE_STAR;printf("\n");

/*-----------------------------------------------------------------------*/
}  /*end routine*/
/*==========================================================================*/



/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
/*  set_bond_potent:                                                        */
/*==========================================================================*/

void set_bond_potent(BOND *bond,GRP_BOND_CON *grp_bond_con,
                          GRP_BOND_WATTS *grp_bond_watts,
                          BUILD_INTRA *build_intra,
                          NAME def_intra_name, NAME user_intra_name,
                          DICT_WORD fun_dict[],int num_fun_dict)

/*=======================================================================*/
{ /*begin routine */ 
  /*=======================================================================*/
  /*           Local Variables                                             */
  
  int i;
  int *ifound_pow,*isearch_pow,*igood_pow;
  int *ifound_con,*isearch_con,*igood_con;
  int *ifound_21,*isearch_21,*igood_21;
  int *ifound_23,*isearch_23,*igood_23;
  int *ifound_33,*isearch_33,*igood_33;
  int *ifound_watt_33,*isearch_watt_33,*igood_watt_33;
  int *ifound_43,*isearch_43,*igood_43;
  int *ifound_46,*isearch_46,*igood_46;
  DATA_BASE_ENTRIES *bond_base;
  CATM_LAB *cbond_base;
  CATM_LAB *cbond_pow,*cbond_con;
  CATM_LAB *cbond_21,*cbond_23,*cbond_33,*cbond_watt_33,*cbond_43,*cbond_46;
  int nbase,nbase2,ibase_want;
  int nsearch,natm_srch,icon_flag;
  int ntyp_pow,ntyp_con;
  int ntyp_21,ntyp_23,ntyp_33,ntyp_watt_33,ntyp_43,ntyp_46;
  int ntyp_21_tot,ntyp_23_tot,ntyp_33_tot,ntyp_watt_33_tot;
  int ntyp_43_tot,ntyp_46_tot;
  int j1,j2,j3,j4,j5,j6;
  char typ[5];

/*=======================================================================*/
/* 0) Search initialization. */
  natm_srch = 2;
  strcpy(typ,"bond");

/*-----------------------------------------------------------------------*/
/* A) Power series                                                       */

  ntyp_pow = bond->ntyp_pow;
  if(ntyp_pow > 0){
    ifound_pow = (int *)cmalloc(ntyp_pow*sizeof(int))-1;
    isearch_pow = (int *)cmalloc(ntyp_pow*sizeof(int))-1;
    igood_pow = (int *)cmalloc(ntyp_pow*sizeof(int))-1;
    cbond_pow = (CATM_LAB *)cmalloc(ntyp_pow*sizeof(CATM_LAB))-1;
    for(i=1;i<=ntyp_pow;i++){
      strcpy(cbond_pow[i].atm1,build_intra->cbond_typ_pow[i].atm1);
      strcpy(cbond_pow[i].atm2,build_intra->cbond_typ_pow[i].atm2);
      strcpy(cbond_pow[i].label,build_intra->cbond_typ_pow[i].label);
    }/*endfor*/
    for(i=1;i<=ntyp_pow;i++){ifound_pow[i]=0;}
    for(i=1;i<=ntyp_pow;i++){igood_pow[i] =6;}
  }/*endif*/
/*-----------------------------------------------------------------------*/
/* B.1) Normal cons                                                      */

  ntyp_con = bond->ntyp_con;
  if(ntyp_con > 0){
    ifound_con = (int *)cmalloc(ntyp_con*sizeof(int))-1;
    isearch_con = (int *)cmalloc(ntyp_con*sizeof(int))-1;
    igood_con = (int *)cmalloc(ntyp_con*sizeof(int))-1;
    cbond_con = (CATM_LAB *)cmalloc(ntyp_con*sizeof(CATM_LAB))-1;
    for(i=1;i<=ntyp_con;i++){
      strcpy(cbond_con[i].atm1,build_intra->cbond_typ_con[i].atm1);
      strcpy(cbond_con[i].atm2,build_intra->cbond_typ_con[i].atm2);
      strcpy(cbond_con[i].label,build_intra->cbond_typ_con[i].label);
    }/*endfor*/
    for(i=1;i<=ntyp_con;i++){ifound_con[i]=0;}
    for(i=1;i<=ntyp_con;i++){igood_con[i] =6;}
  }/*endif*/
/*-----------------------------------------------------------------------*/
/* B.2) Group cons 21 [1-2]                                              */

  ntyp_21 = grp_bond_con->ntyp_21;
  if(ntyp_21 > 0){
    ifound_21 = (int *)cmalloc(ntyp_21*sizeof(int))-1;
    isearch_21 = (int *)cmalloc(ntyp_21*sizeof(int))-1;
    igood_21 = (int *)cmalloc(ntyp_21*sizeof(int))-1;
    cbond_21 = (CATM_LAB *)cmalloc(ntyp_21*sizeof(CATM_LAB))-1;
    for(i=1;i<=ntyp_21;i++){
      strcpy(cbond_21[i].atm1,build_intra->cgrp_typ_21[i].atm1);
      strcpy(cbond_21[i].atm2,build_intra->cgrp_typ_21[i].atm2);
      strcpy(cbond_21[i].label,"");
    }/*endfor*/
    for(i=1;i<=ntyp_21;i++){ifound_21[i]=0;}
    for(i=1;i<=ntyp_21;i++){igood_21[i] =6;}
  }/*endif*/
/*-----------------------------------------------------------------------*/
/* B.3) Group cons 23  [1-2,1-3]                                        */

  ntyp_23 = grp_bond_con->ntyp_23;
  if(ntyp_23 > 0){
    ntyp_23_tot = 2*grp_bond_con->ntyp_23;
    ifound_23 = (int *)cmalloc(ntyp_23_tot*sizeof(int))-1;
    isearch_23 = (int *)cmalloc(ntyp_23_tot*sizeof(int))-1;
    igood_23 = (int *)cmalloc(ntyp_23_tot*sizeof(int))-1;
    cbond_23 = (CATM_LAB *)cmalloc(ntyp_23_tot*sizeof(CATM_LAB))-1;
    for(i=1;i<=ntyp_23;i++){
      j1 = 2*i-1;
      j2 = 2*i;
      strcpy(cbond_23[j1].atm1,build_intra->cgrp_typ_23[i].atm1);
      strcpy(cbond_23[j1].atm2,build_intra->cgrp_typ_23[i].atm2);
      strcpy(cbond_23[j1].label,"");
      strcpy(cbond_23[j2].atm1,build_intra->cgrp_typ_23[i].atm1);
      strcpy(cbond_23[j2].atm2,build_intra->cgrp_typ_23[i].atm3);
      strcpy(cbond_23[j2].label,"");
    }/*endfor*/
    for(i=1;i<=ntyp_23_tot;i++){ifound_23[i]=0;}
    for(i=1;i<=ntyp_23_tot;i++){igood_23[i] =6;}
  }/*endif*/

/*-----------------------------------------------------------------------*/
/* B.4) Group cons 33  [1-2,1-3,2-3]                                     */

  ntyp_33 = grp_bond_con->ntyp_33;
  if(ntyp_33 > 0){
    ntyp_33_tot = 3*grp_bond_con->ntyp_33;
    ifound_33 = (int *)cmalloc(ntyp_33_tot*sizeof(int))-1;
    isearch_33 = (int *)cmalloc(ntyp_33_tot*sizeof(int))-1;
    igood_33 = (int *)cmalloc(ntyp_33_tot*sizeof(int))-1;
    cbond_33 = (CATM_LAB *)cmalloc(ntyp_33_tot*sizeof(CATM_LAB))-1;
    for(i=1;i<=ntyp_33;i++){
      j1 = 3*i-2;
      j2 = 3*i-1;
      j3 = 3*i;
      strcpy(cbond_33[j1].atm1,build_intra->cgrp_typ_33[i].atm1);
      strcpy(cbond_33[j1].atm2,build_intra->cgrp_typ_33[i].atm2);
      strcpy(cbond_33[j1].label,"");
      strcpy(cbond_33[j2].atm1,build_intra->cgrp_typ_33[i].atm1);
      strcpy(cbond_33[j2].atm2,build_intra->cgrp_typ_33[i].atm3);
      strcpy(cbond_33[j2].label,"");
      strcpy(cbond_33[j3].atm1,build_intra->cgrp_typ_33[i].atm2);
      strcpy(cbond_33[j3].atm2,build_intra->cgrp_typ_33[i].atm3);
      strcpy(cbond_33[j3].label,"");
    }/*endfor*/
    for(i=1;i<=ntyp_33_tot;i++){ifound_33[i]=0;}
    for(i=1;i<=ntyp_33_tot;i++){igood_33[i] =6;}
  }/*endif*/

/*-----------------------------------------------------------------------*/
/* B.4) Group watts 33  [1-2,1-3,2-3]                                    */


  ntyp_watt_33 = grp_bond_watts->ntyp_33;
  if(ntyp_watt_33 > 0){
    ntyp_watt_33_tot = 3*grp_bond_watts->ntyp_33;
    ifound_watt_33 = (int *)cmalloc(ntyp_watt_33_tot*sizeof(int))-1;
    isearch_watt_33 = (int *)cmalloc(ntyp_watt_33_tot*sizeof(int))-1;
    igood_watt_33 = (int *)cmalloc(ntyp_watt_33_tot*sizeof(int))-1;
    cbond_watt_33 = (CATM_LAB *)cmalloc(ntyp_watt_33_tot*sizeof(CATM_LAB))-1;
    for(i=1;i<=ntyp_watt_33;i++){
      j1 = 3*i-2;
      j2 = 3*i-1;
      j3 = 3*i;
      strcpy(cbond_watt_33[j1].atm1,build_intra->cgrp_typ_watt_33[i].atm1);
      strcpy(cbond_watt_33[j1].atm2,build_intra->cgrp_typ_watt_33[i].atm2);
      strcpy(cbond_watt_33[j1].label,"");
      strcpy(cbond_watt_33[j2].atm1,build_intra->cgrp_typ_watt_33[i].atm1);
      strcpy(cbond_watt_33[j2].atm2,build_intra->cgrp_typ_watt_33[i].atm3);
      strcpy(cbond_watt_33[j2].label,"");
      strcpy(cbond_watt_33[j3].atm1,build_intra->cgrp_typ_watt_33[i].atm2);
      strcpy(cbond_watt_33[j3].atm2,build_intra->cgrp_typ_watt_33[i].atm3);
      strcpy(cbond_watt_33[j3].label,"");
    }/*endfor*/
    for(i=1;i<=ntyp_watt_33_tot;i++){ifound_watt_33[i]=0;}
    for(i=1;i<=ntyp_watt_33_tot;i++){igood_watt_33[i] =6;}
  }/*endif*/

/*-----------------------------------------------------------------------*/
/* B.5) Group cons 43    [1-2,1-3,1-4]                                   */

  ntyp_43 = grp_bond_con->ntyp_43;
  if(ntyp_43 > 0){
    ntyp_43_tot = 3*grp_bond_con->ntyp_43;
    ifound_43 = (int *)cmalloc(ntyp_43_tot*sizeof(int))-1;
    isearch_43 = (int *)cmalloc(ntyp_43_tot*sizeof(int))-1;
    igood_43 = (int *)cmalloc(ntyp_43_tot*sizeof(int))-1;
    cbond_43 = (CATM_LAB *)cmalloc(ntyp_43_tot*sizeof(CATM_LAB))-1;
    for(i=1;i<=ntyp_43;i++){
      j1 = 3*i-2;
      j2 = 3*i-1;
      j3 = 3*i;
      strcpy(cbond_43[j1].atm1,build_intra->cgrp_typ_43[i].atm1);
      strcpy(cbond_43[j1].atm2,build_intra->cgrp_typ_43[i].atm2);
      strcpy(cbond_43[j1].label,"");
      strcpy(cbond_43[j2].atm1,build_intra->cgrp_typ_43[i].atm1);
      strcpy(cbond_43[j2].atm2,build_intra->cgrp_typ_43[i].atm3);
      strcpy(cbond_43[j2].label,"");
      strcpy(cbond_43[j3].atm1,build_intra->cgrp_typ_43[i].atm1);
      strcpy(cbond_43[j3].atm2,build_intra->cgrp_typ_43[i].atm4);
      strcpy(cbond_43[j3].label,"");
    }/*endfor*/
    for(i=1;i<=ntyp_43_tot;i++){ifound_43[i]=0;}
    for(i=1;i<=ntyp_43_tot;i++){igood_43[i] =6;}
  }/*endif*/

/*-----------------------------------------------------------------------*/
/* B.6) Group cons 46 [1-2,1-3,1-4,2-3,2-4,3-4]                          */

  ntyp_46 = grp_bond_con->ntyp_46;
  if(ntyp_46 > 0){
    ntyp_46_tot = 6*grp_bond_con->ntyp_46;
    ifound_46 = (int *)cmalloc(ntyp_46_tot*sizeof(int))-1;
    isearch_46 = (int *)cmalloc(ntyp_46_tot*sizeof(int))-1;
    igood_46 = (int *)cmalloc(ntyp_46_tot*sizeof(int))-1;
    cbond_46 = (CATM_LAB *)cmalloc(ntyp_46_tot*sizeof(CATM_LAB))-1;
    for(i=1;i<=ntyp_46;i++){
      j1 = 6*i-5;
      j2 = 6*i-4;
      j3 = 6*i-3;
      j4 = 6*i-2;
      j5 = 6*i-1;
      j6 = 6*i;
      strcpy(cbond_46[j1].atm1,build_intra->cgrp_typ_46[i].atm1);
      strcpy(cbond_46[j1].atm2,build_intra->cgrp_typ_46[i].atm2);
      strcpy(cbond_46[j1].label,"");
      strcpy(cbond_46[j2].atm1,build_intra->cgrp_typ_46[i].atm1);
      strcpy(cbond_46[j2].atm2,build_intra->cgrp_typ_46[i].atm3);
      strcpy(cbond_46[j2].label,"");
      strcpy(cbond_46[j3].atm1,build_intra->cgrp_typ_46[i].atm1);
      strcpy(cbond_46[j3].atm2,build_intra->cgrp_typ_46[i].atm4);
      strcpy(cbond_46[j3].label,"");
      strcpy(cbond_46[j4].atm1,build_intra->cgrp_typ_46[i].atm2);
      strcpy(cbond_46[j4].atm2,build_intra->cgrp_typ_46[i].atm3);
      strcpy(cbond_46[j4].label,"");
      strcpy(cbond_46[j5].atm1,build_intra->cgrp_typ_46[i].atm2);
      strcpy(cbond_46[j5].atm2,build_intra->cgrp_typ_46[i].atm4);
      strcpy(cbond_46[j5].label,"");
      strcpy(cbond_46[j6].atm1,build_intra->cgrp_typ_46[i].atm3);
      strcpy(cbond_46[j6].atm2,build_intra->cgrp_typ_46[i].atm4);
      strcpy(cbond_46[j6].label,"");
    }/*endfor*/
    for(i=1;i<=ntyp_46_tot;i++){ifound_46[i]=0;}
    for(i=1;i<=ntyp_46_tot;i++){igood_46[i] =6;}
  }/*endif*/

/*=======================================================================*/
/* I) Search the user base. */
  if(strlen(user_intra_name) != 0){
    nsearch = 1;
/*-----------------------------------------------------------------------*/
/* A) Count the base, malloc and read the base.                          */
    ibase_want = 2;
    count_data_base(user_intra_name,fun_dict,num_fun_dict,&nbase,ibase_want);
    if(nbase>0){
      nbase2 = 2*nbase;
      bond_base = (DATA_BASE_ENTRIES *)
                  cmalloc(nbase*sizeof(DATA_BASE_ENTRIES))-1;
      cbond_base = (CATM_LAB *)cmalloc(nbase2*sizeof(CATM_LAB))-1;
      read_data_base(user_intra_name,fun_dict,num_fun_dict,bond_base,
                   cbond_base,ibase_want,nbase);
    }/*endif*/

/*-----------------------------------------------------------------------*/
/* B.1) Find the power series bonds.                                     */

    if(ntyp_pow > 0 &&(nbase > 0)){
      search_base(nbase,nbase2,cbond_base,ntyp_pow,
                  cbond_pow,igood_pow,ifound_pow,
                  isearch_pow,nsearch,natm_srch,user_intra_name);
      icon_flag = 0;
      assign_base_bond(bond_base,nbase,ifound_pow,bond,isearch_pow,
                       nsearch,icon_flag);
    }/*endif*/

/*-----------------------------------------------------------------------*/
/* B.2) Find the constrained bonds.                                      */

    if(ntyp_con > 0 &&(nbase > 0)){
      search_base(nbase,nbase2,cbond_base,ntyp_con,
                  cbond_con,igood_con,ifound_con,
                  isearch_con,nsearch,natm_srch,user_intra_name);
      icon_flag = 1;
      assign_base_bond(bond_base,nbase,ifound_con,bond,isearch_con,
                       nsearch,icon_flag);
    }/*endif*/

/*-----------------------------------------------------------------------*/
/* B.3) Find the group 21                                                */

  if((ntyp_21 > 0)&&(nbase > 0)){
      search_base(nbase,nbase2,cbond_base,ntyp_21,
                  cbond_21,igood_21,ifound_21,
                  isearch_21,nsearch,natm_srch,user_intra_name);
      icon_flag = 1;
      assign_base_grpbnd(bond_base,nbase,ifound_21,grp_bond_con,
                         grp_bond_watts,
                         isearch_21,nsearch,icon_flag);
  }/*endif*/

/*-----------------------------------------------------------------------*/
/* B.4) Find the group 23                                                */

  if((ntyp_23 > 0)&&(nbase > 0)){
      search_base(nbase,nbase2,cbond_base,ntyp_23_tot,
                  cbond_23,igood_23,ifound_23,
                  isearch_23,nsearch,natm_srch,user_intra_name);
      icon_flag = 2;
      assign_base_grpbnd(bond_base,nbase,ifound_23,grp_bond_con,
                         grp_bond_watts,
                         isearch_23,nsearch,icon_flag);
  }/*endif*/

/*-----------------------------------------------------------------------*/
/* B.5) Find the group 33                                                */

  if((ntyp_33 > 0)&&(nbase > 0)){
      search_base(nbase,nbase2,cbond_base,ntyp_33_tot,
                  cbond_33,igood_33,ifound_33,
                  isearch_33,nsearch,natm_srch,user_intra_name);
      icon_flag = 3;
      assign_base_grpbnd(bond_base,nbase,ifound_33,grp_bond_con,
                         grp_bond_watts,
                         isearch_33,nsearch,icon_flag);
  }/*endif*/

/*-----------------------------------------------------------------------*/
/* B.5) Find the group Watts 33                                          */

  if((ntyp_watt_33 > 0)&&(nbase > 0)){
      search_base(nbase,nbase2,cbond_base,ntyp_watt_33_tot,
                  cbond_watt_33,igood_watt_33,ifound_watt_33,
                  isearch_watt_33,nsearch,natm_srch,user_intra_name);
      icon_flag = 6;
      assign_base_grpbnd(bond_base,nbase,ifound_watt_33,grp_bond_con,
                         grp_bond_watts,
                         isearch_watt_33,nsearch,icon_flag);
  }/*endif*/

/*-----------------------------------------------------------------------*/
/* B.6) Find the group 43                                                */

  if((ntyp_43 > 0)&& (nbase > 0)){
      search_base(nbase,nbase2,cbond_base,ntyp_43_tot,
                  cbond_43,igood_43,ifound_43,
                  isearch_43,nsearch,natm_srch,user_intra_name);
      icon_flag = 4;
      assign_base_grpbnd(bond_base,nbase,ifound_43,grp_bond_con,
                         grp_bond_watts,
                         isearch_43,nsearch,icon_flag);
  }/*endif*/

/*-----------------------------------------------------------------------*/
/* B.7) Find the group 46                                                */

  if((ntyp_46 > 0)&& (nbase > 0)){
      search_base(nbase,nbase2,cbond_base,ntyp_46_tot,
                  cbond_46,igood_46,ifound_46,
                  isearch_46,nsearch,natm_srch,user_intra_name);
      icon_flag = 5;
      assign_base_grpbnd(bond_base,nbase,ifound_46,grp_bond_con,
                         grp_bond_watts,
                         isearch_46,nsearch,icon_flag);
  }/*endif*/

/*-----------------------------------------------------------------------*/
/* D) Free the bases. */
    if(nbase>0){
      cfree(&bond_base[1]);
      cfree(&cbond_base[1]);
    }/*endif*/
  }/*endif*/

/*=======================================================================*/
/* II) Search the default base. */
  nsearch = 2;
/*-----------------------------------------------------------------------*/
/* A) Count the base, malloc and read the base.                          */
  ibase_want = 2;
  count_data_base(def_intra_name,fun_dict,num_fun_dict,&nbase,ibase_want);
  if(nbase>0){
    nbase2 = 2*nbase;
    bond_base = (DATA_BASE_ENTRIES *)
                      cmalloc(nbase*sizeof(DATA_BASE_ENTRIES))-1;
    cbond_base = (CATM_LAB *)cmalloc(nbase2*sizeof(CATM_LAB))-1;
    read_data_base(def_intra_name,fun_dict,num_fun_dict,bond_base,
                   cbond_base,ibase_want,nbase);
  }/*endif*/

/*-----------------------------------------------------------------------*/
/* B.1) Find the power series bonds.                                     */

  if((ntyp_pow > 0)&&(nbase > 0)){
    search_base(nbase,nbase2,cbond_base,ntyp_pow,
                cbond_pow,igood_pow,ifound_pow,
                isearch_pow,nsearch,natm_srch,def_intra_name);
    icon_flag = 0;
    assign_base_bond(bond_base,nbase,ifound_pow,bond,isearch_pow,
                     nsearch,icon_flag);
  }/*endif*/

/*-----------------------------------------------------------------------*/
/* B.2) Find the constrained bonds.                                      */

  if((ntyp_con > 0)&&( nbase > 0)){

    search_base(nbase,nbase2,cbond_base,ntyp_con,
                cbond_con,igood_con,ifound_con,
                isearch_con,nsearch,natm_srch,def_intra_name);
    icon_flag = 1;
    assign_base_bond(bond_base,nbase,ifound_con,bond,isearch_con,
                     nsearch,icon_flag);
  }/*endif*/

/*-----------------------------------------------------------------------*/
/* B.3) Find the group 21                                                */

  if((ntyp_21 > 0)&&(nbase > 0)){
      search_base(nbase,nbase2,cbond_base,ntyp_21,
                  cbond_21,igood_21,ifound_21,
                  isearch_21,nsearch,natm_srch,def_intra_name);
      icon_flag = 1;
      assign_base_grpbnd(bond_base,nbase,ifound_21,grp_bond_con,
                         grp_bond_watts,
                         isearch_21,nsearch,icon_flag);
  }/*endif*/

/*-----------------------------------------------------------------------*/
/* B.4) Find the group 23                                                */

  if((ntyp_23 > 0)&&(nbase > 0)){
      search_base(nbase,nbase2,cbond_base,ntyp_23_tot,
                  cbond_23,igood_23,ifound_23,
                  isearch_23,nsearch,natm_srch,def_intra_name);
      icon_flag = 2;
      assign_base_grpbnd(bond_base,nbase,ifound_23,grp_bond_con,
                         grp_bond_watts,
                         isearch_23,nsearch,icon_flag);
  }/*endif*/

/*-----------------------------------------------------------------------*/
/* B.5) Find the group 33                                                */

  if((ntyp_33 > 0)&&(nbase > 0)){
      search_base(nbase,nbase2,cbond_base,ntyp_33_tot,
                  cbond_33,igood_33,ifound_33,
                  isearch_33,nsearch,natm_srch,def_intra_name);
      icon_flag = 3;
      assign_base_grpbnd(bond_base,nbase,ifound_33,grp_bond_con,
                         grp_bond_watts,
                         isearch_33,nsearch,icon_flag);
  }/*endif*/

/*-----------------------------------------------------------------------*/
/* B.5) Find the group Watts 33                                          */

  if((ntyp_watt_33 > 0)&&(nbase > 0)){
      search_base(nbase,nbase2,cbond_base,ntyp_watt_33_tot,
                  cbond_watt_33,igood_watt_33,ifound_watt_33,
                  isearch_watt_33,nsearch,natm_srch,def_intra_name);
      icon_flag = 6;
      assign_base_grpbnd(bond_base,nbase,ifound_watt_33,grp_bond_con,
                         grp_bond_watts,
                         isearch_watt_33,nsearch,icon_flag);
  }/*endif*/

/*-----------------------------------------------------------------------*/
/* B.6) Find the group 43                                                */

  if((ntyp_43 > 0)&&(nbase > 0)){
      search_base(nbase,nbase2,cbond_base,ntyp_43_tot,
                  cbond_43,igood_43,ifound_43,
                  isearch_43,nsearch,natm_srch,def_intra_name);
      icon_flag = 4;
      assign_base_grpbnd(bond_base,nbase,ifound_43,grp_bond_con,
                         grp_bond_watts,
                         isearch_43,nsearch,icon_flag);
  }/*endif*/

/*-----------------------------------------------------------------------*/
/* B.7) Find the group 46                                                */

  if((ntyp_46 > 0)&&(nbase > 0)){
      search_base(nbase,nbase2,cbond_base,ntyp_46_tot,
                  cbond_46,igood_46,ifound_46,
                  isearch_46,nsearch,natm_srch,def_intra_name);
      icon_flag = 5;
      assign_base_grpbnd(bond_base,nbase,ifound_46,grp_bond_con,
                         grp_bond_watts,
                         isearch_46,nsearch,icon_flag);
  }/*endif*/

/*-----------------------------------------------------------------------*/
/* D) Free the bases. */
  if(nbase>0){
    cfree(&bond_base[1]);
    cfree(&cbond_base[1]);
  }/*endif*/

/*=======================================================================*/
/* III) Die if bonds are missing. */
  if(ntyp_pow>0){atmlst_not_found(ntyp_pow,cbond_pow,ifound_pow,
                                  natm_srch,typ);}
  if(ntyp_con>0){atmlst_not_found(ntyp_con,cbond_con,ifound_con,
                                  natm_srch,typ);}
  if(ntyp_21>0){atmlst_not_found(ntyp_21,cbond_21,ifound_21,
                                 natm_srch,typ);}
  if(ntyp_23>0){atmlst_not_found(ntyp_23,cbond_23,ifound_23,
                                 natm_srch,typ);}
  if(ntyp_33>0){atmlst_not_found(ntyp_33,cbond_33,ifound_33,
                                 natm_srch,typ);}
  if(ntyp_watt_33>0){atmlst_not_found(ntyp_watt_33,cbond_watt_33,
                            ifound_watt_33,natm_srch,typ);}
  if(ntyp_43>0){atmlst_not_found(ntyp_43,cbond_43,ifound_43,
                                 natm_srch,typ);}
  if(ntyp_46>0){atmlst_not_found(ntyp_46,cbond_46,ifound_46,
                                 natm_srch,typ);}
/*=======================================================================*/
/* IV) Free the memory. */
  if(ntyp_pow>0){
    cfree(&ifound_pow[1]);
    cfree(&isearch_pow[1]);
    cfree(&igood_pow[1]);
    cfree(&cbond_pow[1]);
  }/*endif*/
  if(ntyp_con>0){
    cfree(&ifound_con[1]);
    cfree(&isearch_con[1]);
    cfree(&igood_con[1]);
    cfree(&cbond_con[1]);
  }/*endif*/
  if(ntyp_21>0){
    cfree(&ifound_21[1]);
    cfree(&isearch_21[1]);
    cfree(&igood_21[1]);
    cfree(&cbond_21[1]);
  }/*endif*/
  if(ntyp_23>0){
    cfree(&ifound_23[1]);
    cfree(&isearch_23[1]);
    cfree(&igood_23[1]);
    cfree(&cbond_23[1]);
  }/*endif*/
  if(ntyp_33>0){
    cfree(&ifound_33[1]);
    cfree(&isearch_33[1]);
    cfree(&igood_33[1]);
    cfree(&cbond_33[1]);
  }/*endif*/
  if(ntyp_watt_33>0){
    cfree(&ifound_watt_33[1]);
    cfree(&isearch_watt_33[1]);
    cfree(&igood_watt_33[1]);
    cfree(&cbond_watt_33[1]);
  }/*endif*/
  if(ntyp_43>0){
    cfree(&ifound_43[1]);
    cfree(&isearch_43[1]);
    cfree(&igood_43[1]);
    cfree(&cbond_43[1]);
  }/*endif*/
/*-----------------------------------------------------------------------*/
}  /*end routine*/
/*==========================================================================*/




/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
/*  set_bend_potent:                                                        */
/*==========================================================================*/

void set_bend_potent(BEND *bend,BUILD_INTRA *build_intra,
                          NAME def_intra_name, NAME user_intra_name,
                          DICT_WORD fun_dict[],int num_fun_dict)

/*=======================================================================*/
{ /*begin routine */ 
  /*=======================================================================*/
  /*           Local Variables                                             */
  
  int i;
  int *ifound_con,*isearch_con,*igood_con;
  DATA_BASE_ENTRIES *bend_base;
  CATM_LAB *cbend_base;
  CATM_LAB *cbend_con;
  int nbase,nbase2,ibase_want;
  int nsearch,natm_srch,icon_flag;
  int ntyp_con;
  char typ[5];

/*=======================================================================*/
/* 0) Search initialization. */

  natm_srch = 3;
  strcpy(typ,"bend");
  ntyp_con = bend->ntyp_con;
  if(ntyp_con > 0){
    ifound_con = (int *)cmalloc(ntyp_con*sizeof(int))-1;
    isearch_con = (int *)cmalloc(ntyp_con*sizeof(int))-1;
    igood_con = (int *)cmalloc(ntyp_con*sizeof(int))-1;
    cbend_con = (CATM_LAB *)cmalloc(ntyp_con*sizeof(CATM_LAB))-1;
    for(i=1;i<=ntyp_con;i++){
      strcpy(cbend_con[i].atm1,build_intra->cbend_typ_con[i].atm1);
      strcpy(cbend_con[i].atm2,build_intra->cbend_typ_con[i].atm2);
      strcpy(cbend_con[i].atm3,build_intra->cbend_typ_con[i].atm3);
      strcpy(cbend_con[i].label,build_intra->cbend_typ_con[i].label);
    }/*endfor*/
    for(i=1;i<=ntyp_con;i++){ifound_con[i]=0;}
    for(i=1;i<=ntyp_con;i++){igood_con[i]=6;}
  }/*endif*/
  
/*=======================================================================*/
/* I) Search the user base. */
  if(strlen(user_intra_name) != 0){
    nsearch = 1;
/*-----------------------------------------------------------------------*/
/* A) Count the base, malloc and read the base.                          */
    ibase_want = 3;
    count_data_base(user_intra_name,fun_dict,num_fun_dict,&nbase,ibase_want);
    if(nbase>0){
      nbase2 = 2*nbase;
      bend_base = (DATA_BASE_ENTRIES *)
                     cmalloc(nbase*sizeof(DATA_BASE_ENTRIES))-1;
      cbend_base = (CATM_LAB *)cmalloc(nbase2*sizeof(CATM_LAB))-1;
      read_data_base(user_intra_name,fun_dict,num_fun_dict,bend_base,
                   cbend_base,ibase_want,nbase);
    }/*endif*/

/*-----------------------------------------------------------------------*/
/* B) Find the constrained bends.                                        */
    if(ntyp_con > 0 && nbase > 0){
      search_base(nbase,nbase2,cbend_base,ntyp_con,
                       cbend_con,igood_con,ifound_con,
                       isearch_con,nsearch,natm_srch,def_intra_name);
      icon_flag = 1;
      assign_base_bend(bend_base,nbase,ifound_con,bend,isearch_con,
                       nsearch,icon_flag);
    }/*endif*/

/*-----------------------------------------------------------------------*/
/* D) Free the bases. */
    if(nbase>0){
      cfree(&bend_base[1]);
      cfree(&cbend_base[1]);
    }/*endif*/
  }/*endif*/

/*=======================================================================*/
/* II) Search the default base. */
  nsearch = 2;
/*-----------------------------------------------------------------------*/
/* A) Count the base, malloc and read the base.                          */
  ibase_want = 3;
  count_data_base(user_intra_name,fun_dict,num_fun_dict,&nbase,ibase_want);
  if(nbase>0){
    nbase2 = 2*nbase;
    bend_base = (DATA_BASE_ENTRIES *)cmalloc(nbase*sizeof(DATA_BASE_ENTRIES))-1;
    cbend_base = (CATM_LAB *)cmalloc(nbase2*sizeof(CATM_LAB))-1;
    read_data_base(def_intra_name,fun_dict,num_fun_dict,bend_base,
                 cbend_base,ibase_want,nbase);
  }/*endif*/

/*-----------------------------------------------------------------------*/
/* B) Find the constrained bends.                                        */
  if(ntyp_con > 0 && nbase > 0){

    search_base(nbase,nbase2,cbend_base,ntyp_con,
                     cbend_con,igood_con,ifound_con,
                     isearch_con,nsearch,natm_srch,def_intra_name);
    icon_flag = 1;
    assign_base_bend(bend_base,nbase,ifound_con,bend,isearch_con,
                     nsearch,icon_flag);
  }/*endif*/
/*-----------------------------------------------------------------------*/
/* D) Free the bases. */
  if(nbase>0){
    cfree(&bend_base[1]);
    cfree(&cbend_base[1]);
  }/*endif*/

/*=======================================================================*/
/* III) Die if bends are missing. */
  if(ntyp_con>0){atmlst_not_found(ntyp_con,cbend_con,ifound_con,
                 natm_srch,typ);}
/*=======================================================================*/
/* IV) Free the memory. */
  if(ntyp_con>0){
    cfree(&ifound_con[1]);
    cfree(&isearch_con[1]);
    cfree(&igood_con[1]);
    cfree(&cbend_con[1]);
  }/*endif*/
/*-----------------------------------------------------------------------*/
}  /*end routine*/
/*==========================================================================*/



/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
/*  set_tors_potent:                                                        */
/*==========================================================================*/

void set_tors_potent(TORS *tors,BUILD_INTRA *build_intra,
                          NAME def_intra_name, NAME user_intra_name,
                          DICT_WORD fun_dict[],int num_fun_dict)

/*=======================================================================*/
{ /*begin routine */ 
  /*=======================================================================*/
  /*           Local Variables                                             */
  
  int i;
  int *ifound_pow,*isearch_pow,*igood_pow;
  int *ifound_con,*isearch_con,*igood_con;
  DATA_BASE_ENTRIES *tors_base;
  CATM_LAB *ctors_base;
  CATM_LAB *ctors_pow,*ctors_con;
  int nbase,nbase2,ibase_want;
  int nsearch,natm_srch,icon_flag;
  int ntyp_pow,ntyp_con;
  char typ[5];

/*=======================================================================*/
/* 0) Search initialization. */
  natm_srch = 4;
  strcpy(typ,"tors");
  ntyp_pow = tors->ntyp_pow;
  ntyp_con = tors->ntyp_con;
  if(ntyp_pow > 0){
    ifound_pow = (int *)cmalloc(ntyp_pow*sizeof(int))-1;
    isearch_pow = (int *)cmalloc(ntyp_pow*sizeof(int))-1;
    igood_pow = (int *)cmalloc(ntyp_pow*sizeof(int))-1;
    ctors_pow = (CATM_LAB *)cmalloc(ntyp_pow*sizeof(CATM_LAB))-1;
    for(i=1;i<=ntyp_pow;i++){
      strcpy(ctors_pow[i].atm1,build_intra->ctors_typ_pow[i].atm1);
      strcpy(ctors_pow[i].atm2,build_intra->ctors_typ_pow[i].atm2);
      strcpy(ctors_pow[i].atm3,build_intra->ctors_typ_pow[i].atm3);
      strcpy(ctors_pow[i].atm4,build_intra->ctors_typ_pow[i].atm4);
      strcpy(ctors_pow[i].label,build_intra->ctors_typ_pow[i].label);
    }/*endfor*/
    for(i=1;i<=ntyp_pow;i++){ifound_pow[i]=0;}
    for(i=1;i<=ntyp_pow;i++){igood_pow[i]=6;}
  }/*endif*/
  if(ntyp_con > 0){
    ifound_con = (int *)cmalloc(ntyp_con*sizeof(int))-1;
    isearch_con = (int *)cmalloc(ntyp_con*sizeof(int))-1;
    igood_con = (int *)cmalloc(ntyp_con*sizeof(int))-1;
    ctors_con = (CATM_LAB *)cmalloc(ntyp_pow*sizeof(CATM_LAB))-1;
    for(i=1;i<=ntyp_con;i++){
      strcpy(ctors_con[i].atm1,build_intra->ctors_typ_con[i].atm1);
      strcpy(ctors_con[i].atm2,build_intra->ctors_typ_con[i].atm2);
      strcpy(ctors_con[i].atm3,build_intra->ctors_typ_con[i].atm3);
      strcpy(ctors_con[i].atm4,build_intra->ctors_typ_con[i].atm4);
      strcpy(ctors_con[i].label,build_intra->ctors_typ_con[i].label);
    }/*endfor*/
    for(i=1;i<=ntyp_con;i++){ifound_con[i]=0;}
    for(i=1;i<=ntyp_con;i++){igood_con[i]=6;}
  }/*endif*/
  
/*=======================================================================*/
/* I) Search the user base. */
  if(strlen(user_intra_name) != 0){
    nsearch = 1;
/*-----------------------------------------------------------------------*/
/* A) Count the base, malloc and read the base.                          */
    ibase_want = 4;
    count_data_base(user_intra_name,fun_dict,num_fun_dict,&nbase,ibase_want);
    if(nbase>0){
      nbase2 = 2*nbase;
      tors_base = (DATA_BASE_ENTRIES *)
                          cmalloc(nbase*sizeof(DATA_BASE_ENTRIES))-1;
      ctors_base = (CATM_LAB *)cmalloc(nbase2*sizeof(CATM_LAB))-1;
      read_data_base(user_intra_name,fun_dict,num_fun_dict,tors_base,
                   ctors_base,ibase_want,nbase);
    }/*endif*/
/*-----------------------------------------------------------------------*/
/* B) Find the power series torss.                                       */
    if(ntyp_pow > 0 && nbase > 0){
      search_base(nbase,nbase2,ctors_base,ntyp_pow,
                       ctors_pow,igood_pow,ifound_pow,
                       isearch_pow,nsearch,natm_srch,user_intra_name);
      icon_flag = 0;
      assign_base_tors(tors_base,nbase,ifound_pow,tors,isearch_pow,
                       nsearch,icon_flag);
    }/*endif*/

/*-----------------------------------------------------------------------*/
/* C) Find the constrained torss.                                        */
    if(ntyp_con > 0 && nbase > 0){
      search_base(nbase,nbase2,ctors_base,ntyp_con,
                       ctors_con,igood_con,ifound_con,
                       isearch_con,nsearch,natm_srch,user_intra_name);
      icon_flag = 1;
      assign_base_tors(tors_base,nbase,ifound_con,tors,isearch_con,
                       nsearch,icon_flag);
    }/*endif*/

/*-----------------------------------------------------------------------*/
/* D) Free the bases. */
    if(nbase>0){
      cfree(&tors_base[1]);
      cfree(&ctors_base[1]);
    }/*endif*/
  }/*endif*/

/*=======================================================================*/
/* II) Search the default base. */
  nsearch = 2;
/*-----------------------------------------------------------------------*/
/* A) Count the base, malloc and read the base.                          */
  ibase_want = 4;
  count_data_base(def_intra_name,fun_dict,num_fun_dict,&nbase,ibase_want);
  if(nbase>0){
    nbase2 = 2*nbase;
    tors_base = (DATA_BASE_ENTRIES *)
                        cmalloc(nbase*sizeof(DATA_BASE_ENTRIES))-1;
    ctors_base = (CATM_LAB *)cmalloc(nbase2*sizeof(CATM_LAB))-1;
    read_data_base(def_intra_name,fun_dict,num_fun_dict,tors_base,
                 ctors_base,ibase_want,nbase);
  }/*endif*/

/*-----------------------------------------------------------------------*/
/* B) Find the power series torss.                                       */
  if(ntyp_pow > 0 && nbase > 0){
    search_base(nbase,nbase2,ctors_base,ntyp_pow,
                     ctors_pow,igood_pow,ifound_pow,
                     isearch_pow,nsearch,natm_srch,def_intra_name);
    icon_flag = 0;
    assign_base_tors(tors_base,nbase,ifound_pow,tors,isearch_pow,
                       nsearch,icon_flag);
  }/*endif*/
/*-----------------------------------------------------------------------*/
/* C) Find the constrained torss.                                        */
  if(ntyp_con > 0 && nbase > 0){

    search_base(nbase,nbase2,ctors_base,ntyp_con,
                     ctors_con,igood_con,ifound_con,
                     isearch_con,nsearch,natm_srch,def_intra_name);
    icon_flag = 1;
    assign_base_tors(tors_base,nbase,ifound_con,tors,isearch_con,
                     nsearch,icon_flag);
  }/*endif*/
/*-----------------------------------------------------------------------*/
/* D) Free the bases. */
  if(nbase>0){
    cfree(&tors_base[1]);
    cfree(&ctors_base[1]);
  }/*endif*/

/*=======================================================================*/
/* III) Die if torss are missing. */
  if(ntyp_pow>0){atmlst_not_found(ntyp_pow,ctors_pow,ifound_pow,
                 natm_srch,typ);}
  if(ntyp_con>0){atmlst_not_found(ntyp_con,ctors_con,ifound_con,
                 natm_srch,typ);}
/*=======================================================================*/
/* IV) Free the memory. */
  if(ntyp_pow>0){
    cfree(&ifound_pow[1]);
    cfree(&isearch_pow[1]);
    cfree(&igood_pow[1]);
    cfree(&ctors_pow[1]);
  }/*endif*/
  if(ntyp_con>0){
    cfree(&ifound_con[1]);
    cfree(&isearch_con[1]);
    cfree(&igood_con[1]);
    cfree(&ctors_con[1]);
  }/*endif*/
/*-----------------------------------------------------------------------*/
}  /*end routine*/
/*==========================================================================*/



/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
/*  set_onfo_potent:                                                        */
/*==========================================================================*/

void set_onfo_potent(ONFO *onfo,BUILD_INTRA *build_intra,
                          NAME def_intra_name, NAME user_intra_name,
                          DICT_WORD fun_dict[],int num_fun_dict)

/*=======================================================================*/
{ /*begin routine */ 
  /*=======================================================================*/
  /*           Local Variables                                             */
  
  int i;
  int *ifound,*isearch,*igood;
  DATA_BASE_ENTRIES *onfo_base;
  CATM_LAB *confo_base;
  CATM_LAB *confo;
  int nbase,nbase2,ibase_want;
  int nsearch,natm_srch,icon_flag;
  int ntyp;
  char typ[5];

/*=======================================================================*/
/* 0) Search initialization. */
  natm_srch = 2;
  strcpy(typ,"onfo");
  ntyp= onfo->ntyp;
  if(ntyp> 0){
    ifound= (int *)cmalloc(ntyp*sizeof(int))-1;
    isearch= (int *)cmalloc(ntyp*sizeof(int))-1;
    igood= (int *)cmalloc(ntyp*sizeof(int))-1;
    confo = (CATM_LAB *)cmalloc(ntyp*sizeof(CATM_LAB))-1;
    for(i=1;i<=ntyp;i++){
      strcpy(confo[i].atm1,build_intra->confo_typ[i].atm1);
      strcpy(confo[i].atm2,build_intra->confo_typ[i].atm2);
      strcpy(confo[i].label,build_intra->confo_typ[i].label);
    }/*endfor*/
    for(i=1;i<=ntyp;i++){ifound[i]=0;}
    for(i=1;i<=ntyp;i++){igood[i]=6;}
  }/*endif*/
  
/*=======================================================================*/
/* I) Search the user base. */
  if(strlen(user_intra_name) != 0){
    nsearch = 1;
/*-----------------------------------------------------------------------*/
/* A) Count the base, malloc and read the base.                          */
    ibase_want = 5;
    count_data_base(user_intra_name,fun_dict,num_fun_dict,&nbase,ibase_want);
    if(nbase>0){
      nbase2 = 2*nbase;
      onfo_base = (DATA_BASE_ENTRIES *)
                           cmalloc(nbase*sizeof(DATA_BASE_ENTRIES))-1;
      confo_base = (CATM_LAB *)cmalloc(nbase2*sizeof(CATM_LAB))-1;
      read_data_base(user_intra_name,fun_dict,num_fun_dict,onfo_base,
                   confo_base,ibase_want,nbase);
    }/*endif*/
/*-----------------------------------------------------------------------*/
/* B) Find the onfos.                                       */
    if(ntyp> 0 && nbase > 0){
      search_base(nbase,nbase2,confo_base,ntyp,
                       confo,igood,ifound,
                       isearch,nsearch,natm_srch,user_intra_name);
      icon_flag = 0;
      assign_base_onfo(onfo_base,nbase,ifound,onfo,isearch,
                       nsearch,icon_flag);
    }/*endif*/

/*-----------------------------------------------------------------------*/
/* D) Free the bases. */
    if(nbase>0){
      cfree(&onfo_base[1]);
      cfree(&confo_base[1]);
    }/*endif*/
  }/*endif*/

/*=======================================================================*/
/* II) Search the default base. */
  nsearch = 2;
/*-----------------------------------------------------------------------*/
/* A) Count the base, malloc and read the base.                          */
  ibase_want = 5;
  count_data_base(def_intra_name,fun_dict,num_fun_dict,&nbase,ibase_want);
  if(nbase>0){
    nbase2 = 2*nbase;
    onfo_base = (DATA_BASE_ENTRIES *)
                      cmalloc(nbase*sizeof(DATA_BASE_ENTRIES))-1;
    confo_base = (CATM_LAB *)cmalloc(nbase2*sizeof(CATM_LAB))-1;
    read_data_base(def_intra_name,fun_dict,num_fun_dict,onfo_base,
                 confo_base,ibase_want,nbase);
  }/*endif*/

/*-----------------------------------------------------------------------*/
/* B) Find the onfos.                                       */
  if(ntyp> 0 && nbase > 0){
    search_base(nbase,nbase2,confo_base,ntyp,
                     confo,igood,ifound,
                     isearch,nsearch,natm_srch,def_intra_name);
    icon_flag = 0;
    assign_base_onfo(onfo_base,nbase,ifound,onfo,isearch,
                       nsearch,icon_flag);
  }/*endif*/
/*-----------------------------------------------------------------------*/
/* D) Free the bases. */
  if(nbase>0){
    cfree(&onfo_base[1]);
    cfree(&confo_base[1]);
  }/*endif*/

/*=======================================================================*/
/* III) Die if onfos are missing. */
  if(ntyp>0){atmlst_not_found(ntyp,confo,ifound,
                 natm_srch,typ);}
/*=======================================================================*/
/* IV) Free the memory. */
  if(ntyp>0){
    cfree(&ifound[1]);
    cfree(&isearch[1]);
    cfree(&igood[1]);
    cfree(&confo[1]);
  }/*endif*/
/*-----------------------------------------------------------------------*/
}  /*end routine*/
/*==========================================================================*/



/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
/*  set_bend_bnd_potent:                                                    */
/*==========================================================================*/

void set_bend_bnd_potent(BEND_BND *bend_bnd,BUILD_INTRA *build_intra,
                              NAME def_intra_name, NAME user_intra_name,
                              DICT_WORD fun_dict[],int num_fun_dict)

/*=======================================================================*/
{ /*begin routine */ 
  /*=======================================================================*/
  /*           Local Variables                                             */
  
  int i;
  int *ifound,*isearch,*igood;
  DATA_BASE_ENTRIES *bend_bnd_base;
  CATM_LAB *cbend_bnd_base;
  CATM_LAB *cbend_bnd;
  int nbase,nbase2,ibase_want;
  int nsearch,natm_srch,icon_flag;
  int ntyp;
  char typ[5];

/*=======================================================================*/
/* 0) Search initialization. */
  natm_srch = 3;
  strcpy(typ,"bend");
  ntyp= bend_bnd->ntyp;
  if(ntyp> 0){
    ifound= (int *)cmalloc(ntyp*sizeof(int))-1;
    isearch= (int *)cmalloc(ntyp*sizeof(int))-1;
    igood= (int *)cmalloc(ntyp*sizeof(int))-1;
    cbend_bnd = (CATM_LAB *)cmalloc(ntyp*sizeof(CATM_LAB))-1;
    for(i=1;i<=ntyp;i++){
      strcpy(cbend_bnd[i].atm1,build_intra->cbend_bnd_typ[i].atm1);
      strcpy(cbend_bnd[i].atm2,build_intra->cbend_bnd_typ[i].atm2);
      strcpy(cbend_bnd[i].atm3,build_intra->cbend_bnd_typ[i].atm3);
      strcpy(cbend_bnd[i].label,build_intra->cbend_bnd_typ[i].label);
    }/*endfor*/
    for(i=1;i<=ntyp;i++){ifound[i]=0;}
    for(i=1;i<=ntyp;i++){igood[i]=6;}
  }/*endif*/
  
/*=======================================================================*/
/* I) Search the user base. */
  if(strlen(user_intra_name) != 0){
    nsearch = 1;
/*-----------------------------------------------------------------------*/
/* A) Count the base, malloc and read the base.                          */
    ibase_want = 3;
    count_data_base(user_intra_name,fun_dict,num_fun_dict,&nbase,ibase_want);
    if(nbase>0){
      nbase2 = 2*nbase;
      bend_bnd_base=(DATA_BASE_ENTRIES *)
                    cmalloc(nbase*sizeof(DATA_BASE_ENTRIES))-1;
      cbend_bnd_base = (CATM_LAB *)cmalloc(nbase2*sizeof(CATM_LAB))-1;
      read_data_base(user_intra_name,fun_dict,num_fun_dict,bend_bnd_base,
                   cbend_bnd_base,ibase_want,nbase);
    }/*endif*/
/*-----------------------------------------------------------------------*/
/* B) Find the bend_bnds.                                       */
    if(ntyp> 0 && nbase > 0){
      search_base(nbase,nbase2,cbend_bnd_base,ntyp,
                       cbend_bnd,igood,ifound,
                       isearch,nsearch,natm_srch,user_intra_name);
      icon_flag = 0;
      assign_base_bend_bnd(bend_bnd_base,nbase,ifound,bend_bnd,isearch,
                       nsearch,icon_flag,build_intra);
    }/*endif*/

/*-----------------------------------------------------------------------*/
/* D) Free the bases. */
    if(nbase>0){
      cfree(&bend_bnd_base[1]);
      cfree(&cbend_bnd_base[1]);
    }/*endif*/
  }/*endif*/

/*=======================================================================*/
/* II) Search the default base. */
  nsearch = 2;
/*-----------------------------------------------------------------------*/
/* A) Count the base, malloc and read the base.                          */
  ibase_want = 3;
  count_data_base(def_intra_name,fun_dict,num_fun_dict,&nbase,ibase_want);
  if(nbase>0){
    nbase2 = 2*nbase;
    bend_bnd_base = (DATA_BASE_ENTRIES *)
                    cmalloc(nbase*sizeof(DATA_BASE_ENTRIES))-1;
    cbend_bnd_base = (CATM_LAB *)cmalloc(nbase2*sizeof(CATM_LAB))-1;
    read_data_base(def_intra_name,fun_dict,num_fun_dict,bend_bnd_base,
                 cbend_bnd_base,ibase_want,nbase);
  }/*endif*/

/*-----------------------------------------------------------------------*/
/* B) Find the bend_bnds.                                       */
  if(ntyp> 0 && nbase > 0){
    search_base(nbase,nbase2,cbend_bnd_base,ntyp,
                     cbend_bnd,igood,ifound,
                     isearch,nsearch,natm_srch,def_intra_name);
    icon_flag = 0;
    assign_base_bend_bnd(bend_bnd_base,nbase,ifound,bend_bnd,isearch,
                       nsearch,icon_flag,build_intra);
  }/*endif*/
/*-----------------------------------------------------------------------*/
/* D) Free the bases. */
  if(nbase>0){
    cfree(&bend_bnd_base[1]);
    cfree(&cbend_bnd_base[1]);
  }/*endif*/

/*=======================================================================*/
/* III) Die if bend_bnds are missing. */
  if(ntyp>0){atmlst_not_found(ntyp,cbend_bnd,ifound,
                 natm_srch,typ);}
/*=======================================================================*/
/* IV) Free the memory. */
  if(ntyp>0){
    cfree(&ifound[1]);
    cfree(&isearch[1]);
    cfree(&igood[1]);
    cfree(&cbend_bnd[1]);
  }/*endif*/
/*-----------------------------------------------------------------------*/
}  /*end routine*/
/*=======================================================================*/



/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
/*  extract_pure_bends                                                      */
/*==========================================================================*/

void extract_pure_bends(BEND_BND *bend_bnd,BEND *bend,
                        BUILD_INTRA *build_intra)

/*=======================================================================*/
{ /*begin routine */ 
/*=======================================================================*/
/*           Local Variables                                             */
  
  int i,nbend,nbend_bnd,ityp,ntyp,ntyp_pure;
  int nbend_typ,nbend_bnd_typ,ipure;

/*=======================================================================*/
/* 0) Construct the type maps of bends and bend_bnds                     */

   ntyp_pure = 0;   ntyp  = 0;   
   for(i=1;i<=bend_bnd->ntyp;i++){
      ipure = build_intra->ibend_bnd_typ_pure[i];
      if(ipure==0){ntyp++;build_intra->ibend_bnd_typ_map[i]=ntyp;}
      if(ipure==1){ntyp_pure++;
                   build_intra->ibend_bnd_typ_map[i]=ntyp_pure;}
   }/*endfor*/

/*=======================================================================*/
/* I) Malloc the bend memory                                            */

   bend->j1_pow   = (int *) cmalloc(bend_bnd->num*sizeof(int))-1;
   bend->j2_pow   = (int *) cmalloc(bend_bnd->num*sizeof(int))-1;
   bend->j3_pow   = (int *) cmalloc(bend_bnd->num*sizeof(int))-1;
   bend->jtyp_pow = (int *) cmalloc(bend_bnd->num*sizeof(int))-1;

/*=======================================================================*/
/* II) Loop over the bend_bnds picking out the pure bends                 */

   nbend = 0;   nbend_bnd = 0;
   for(i=1;i<=bend_bnd->num;i++){
    ityp = bend_bnd->jtyp[i];
    if(build_intra->ibend_bnd_typ_pure[ityp]==1){
     nbend++;
     bend->j1_pow[nbend]   = bend_bnd->j1[i];
     bend->j2_pow[nbend]   = bend_bnd->j2[i];
     bend->j3_pow[nbend]   = bend_bnd->j3[i];
     bend->jtyp_pow[nbend] = build_intra->ibend_bnd_typ_map[ityp];
    }/*endif*/
    if(build_intra->ibend_bnd_typ_pure[ityp]==0){
     nbend_bnd++;
     bend_bnd->j1[nbend_bnd] = bend_bnd->j1[i];
     bend_bnd->j2[nbend_bnd] = bend_bnd->j2[i];
     bend_bnd->j3[nbend_bnd] = bend_bnd->j3[i];
     bend_bnd->jtyp[nbend_bnd]   = build_intra->ibend_bnd_typ_map[ityp];
    }/*endif*/
   }/*endfor*/
   bend_bnd->num = nbend_bnd;
   bend->npow    = nbend;

/*=======================================================================*/
/* II) Realloc the  memory                                               */

   bend_bnd->nbend_bnd_mall = bend_bnd->num;
   bend_bnd->j1  = (int *)   crealloc(&(bend_bnd->j1)[1],
				      bend_bnd->num*sizeof(int))-1;
   bend_bnd->j2  = (int *)   crealloc(&(bend_bnd->j2)[1],
				      bend_bnd->num*sizeof(int))-1;
   bend_bnd->j3  = (int *)   crealloc(&(bend_bnd->j3)[1],
				      bend_bnd->num*sizeof(int))-1;
   bend_bnd->jtyp = (int *)  crealloc(&(bend_bnd->jtyp)[1],
				      bend_bnd->num*sizeof(int))-1;
   bend->j1_pow  = (int *)   crealloc(&(bend->j1_pow)[1],
                                      bend->npow*sizeof(int))-1;
   bend->j2_pow  = (int *)   crealloc(&(bend->j2_pow)[1],
                                      bend->npow*sizeof(int))-1;
   bend->j3_pow  = (int *)   crealloc(&(bend->j3_pow)[1],
                                      bend->npow*sizeof(int))-1;
   bend->jtyp_pow = (int *)  crealloc(&(bend->jtyp_pow)[1],
                                      bend->npow*sizeof(int))-1;
   bend->nbend_pow_mall = bend->npow;
   bend->nbend_typ_pow_mall = bend_bnd->ntyp;

/*=======================================================================*/
/* III) Malloc the bend typ memory                                      */

    bend->eq_pow = (double *)cmalloc(bend_bnd->ntyp*sizeof(double))-1;
    bend->c_0    = (double *)cmalloc(bend_bnd->ntyp*sizeof(double))-1;
    bend->c_1    = (double *)cmalloc(bend_bnd->ntyp*sizeof(double))-1;
    bend->c_2    = (double *)cmalloc(bend_bnd->ntyp*sizeof(double))-1;
    bend->c_3    = (double *)cmalloc(bend_bnd->ntyp*sizeof(double))-1;
    bend->c_4    = (double *)cmalloc(bend_bnd->ntyp*sizeof(double))-1;
    bend->c_5    = (double *)cmalloc(bend_bnd->ntyp*sizeof(double))-1;
    bend->c_6    = (double *)cmalloc(bend_bnd->ntyp*sizeof(double))-1;
    bend->s_0    = (double *)cmalloc(bend_bnd->ntyp*sizeof(double))-1;
    bend->s_1    = (double *)cmalloc(bend_bnd->ntyp*sizeof(double))-1;
    bend->s_2    = (double *)cmalloc(bend_bnd->ntyp*sizeof(double))-1;
    bend->s_3    = (double *)cmalloc(bend_bnd->ntyp*sizeof(double))-1;
    bend->s_4    = (double *)cmalloc(bend_bnd->ntyp*sizeof(double))-1;
    bend->s_5    = (double *)cmalloc(bend_bnd->ntyp*sizeof(double))-1;
    bend->s_6    = (double *)cmalloc(bend_bnd->ntyp*sizeof(double))-1;
    bend->dc_0   = (double *)cmalloc(bend_bnd->ntyp*sizeof(double))-1;
    bend->dc_1   = (double *)cmalloc(bend_bnd->ntyp*sizeof(double))-1;
    bend->dc_2   = (double *)cmalloc(bend_bnd->ntyp*sizeof(double))-1;
    bend->dc_3   = (double *)cmalloc(bend_bnd->ntyp*sizeof(double))-1;
    bend->dc_4   = (double *)cmalloc(bend_bnd->ntyp*sizeof(double))-1;
    bend->dc_5   = (double *)cmalloc(bend_bnd->ntyp*sizeof(double))-1;
    bend->dc_6   = (double *)cmalloc(bend_bnd->ntyp*sizeof(double))-1;
    bend->ds_0   = (double *)cmalloc(bend_bnd->ntyp*sizeof(double))-1;
    bend->ds_1   = (double *)cmalloc(bend_bnd->ntyp*sizeof(double))-1;
    bend->ds_2   = (double *)cmalloc(bend_bnd->ntyp*sizeof(double))-1;
    bend->ds_3   = (double *)cmalloc(bend_bnd->ntyp*sizeof(double))-1;
    bend->ds_4   = (double *)cmalloc(bend_bnd->ntyp*sizeof(double))-1;
    bend->ds_5   = (double *)cmalloc(bend_bnd->ntyp*sizeof(double))-1;
    bend->ds_6   = (double *)cmalloc(bend_bnd->ntyp*sizeof(double))-1;

/*=======================================================================*/
/* III) Loop over the bend_bnds picking out the pure bend interactions    */

   nbend_typ = 0;   nbend_bnd_typ = 0;
   for(i=1;i<=bend_bnd->ntyp;i++){
    if(build_intra->ibend_bnd_typ_pure[i]==1){
     nbend_typ++;
     bend->eq_pow[nbend_typ] = bend_bnd->eq_bend[i];
     bend->c_0[nbend_typ]    = bend_bnd->cbend_0[i];
     bend->c_1[nbend_typ]    = bend_bnd->cbend_1[i];
     bend->c_2[nbend_typ]    = bend_bnd->cbend_2[i];
     bend->c_3[nbend_typ]    = bend_bnd->cbend_3[i];
     bend->c_4[nbend_typ]    = bend_bnd->cbend_4[i];
     bend->c_5[nbend_typ]    = bend_bnd->cbend_5[i];
     bend->c_6[nbend_typ]    = bend_bnd->cbend_6[i];
     bend->s_0[nbend_typ]    = bend_bnd->sbend_0[i];
     bend->s_1[nbend_typ]    = bend_bnd->sbend_1[i];
     bend->s_2[nbend_typ]    = bend_bnd->sbend_2[i];
     bend->s_3[nbend_typ]    = bend_bnd->sbend_3[i];
     bend->s_4[nbend_typ]    = bend_bnd->sbend_4[i];
     bend->s_5[nbend_typ]    = bend_bnd->sbend_5[i];
     bend->s_6[nbend_typ]    = bend_bnd->sbend_6[i]; 
     bend->dc_0[nbend_typ]   = bend_bnd->dcbend_0[i];
     bend->dc_1[nbend_typ]   = bend_bnd->dcbend_1[i];
     bend->dc_2[nbend_typ]   = bend_bnd->dcbend_2[i]; 
     bend->dc_3[nbend_typ]   = bend_bnd->dcbend_3[i]; 
     bend->dc_4[nbend_typ]   = bend_bnd->dcbend_4[i]; 
     bend->dc_5[nbend_typ]   = bend_bnd->dcbend_5[i]; 
     bend->dc_6[nbend_typ]   = bend_bnd->dcbend_6[i];
     bend->ds_0[nbend_typ]   = bend_bnd->dsbend_0[i]; 
     bend->ds_1[nbend_typ]   = bend_bnd->dsbend_1[i]; 
     bend->ds_2[nbend_typ]   = bend_bnd->dsbend_2[i];  
     bend->ds_3[nbend_typ]   = bend_bnd->dsbend_3[i];  
     bend->ds_4[nbend_typ]   = bend_bnd->dsbend_4[i];  
     bend->ds_5[nbend_typ]   = bend_bnd->dsbend_5[i];  
     bend->ds_6[nbend_typ]   = bend_bnd->dsbend_6[i]; 
    }/*endif*/
    if(build_intra->ibend_bnd_typ_pure[i]==0){
     nbend_bnd_typ++;
     bend_bnd->eq_bend[nbend_bnd_typ]  = bend_bnd->eq_bend[i];
     bend_bnd->cbend_0[nbend_bnd_typ]  = bend_bnd->cbend_0[i];
     bend_bnd->cbend_1[nbend_bnd_typ]  = bend_bnd->cbend_1[i];
     bend_bnd->cbend_2[nbend_bnd_typ]  = bend_bnd->cbend_2[i];
     bend_bnd->cbend_3[nbend_bnd_typ]  = bend_bnd->cbend_3[i];
     bend_bnd->cbend_4[nbend_bnd_typ]  = bend_bnd->cbend_4[i];
     bend_bnd->cbend_5[nbend_bnd_typ]  = bend_bnd->cbend_5[i];
     bend_bnd->cbend_6[nbend_bnd_typ]  = bend_bnd->cbend_6[i];
     bend_bnd->sbend_0[nbend_bnd_typ]  = bend_bnd->sbend_0[i];
     bend_bnd->sbend_1[nbend_bnd_typ]  = bend_bnd->sbend_1[i];
     bend_bnd->sbend_2[nbend_bnd_typ]  = bend_bnd->sbend_2[i];
     bend_bnd->sbend_3[nbend_bnd_typ]  = bend_bnd->sbend_3[i];
     bend_bnd->sbend_4[nbend_bnd_typ]  = bend_bnd->sbend_4[i];
     bend_bnd->sbend_5[nbend_bnd_typ]  = bend_bnd->sbend_5[i];
     bend_bnd->sbend_6[nbend_bnd_typ]  = bend_bnd->sbend_6[i]; 
     bend_bnd->dcbend_0[nbend_bnd_typ] = bend_bnd->dcbend_0[i];
     bend_bnd->dcbend_1[nbend_bnd_typ] = bend_bnd->dcbend_1[i];
     bend_bnd->dcbend_2[nbend_bnd_typ] = bend_bnd->dcbend_2[i]; 
     bend_bnd->dcbend_3[nbend_bnd_typ] = bend_bnd->dcbend_3[i]; 
     bend_bnd->dcbend_4[nbend_bnd_typ] = bend_bnd->dcbend_4[i]; 
     bend_bnd->dcbend_5[nbend_bnd_typ] = bend_bnd->dcbend_5[i]; 
     bend_bnd->dcbend_6[nbend_bnd_typ] = bend_bnd->dcbend_6[i];
     bend_bnd->dsbend_0[nbend_bnd_typ] = bend_bnd->dsbend_0[i]; 
     bend_bnd->dsbend_1[nbend_bnd_typ] = bend_bnd->dsbend_1[i]; 
     bend_bnd->dsbend_2[nbend_bnd_typ] = bend_bnd->dsbend_2[i];  
     bend_bnd->dsbend_3[nbend_bnd_typ] = bend_bnd->dsbend_3[i];  
     bend_bnd->dsbend_4[nbend_bnd_typ] = bend_bnd->dsbend_4[i];  
     bend_bnd->dsbend_5[nbend_bnd_typ] = bend_bnd->dsbend_5[i];  
     bend_bnd->dsbend_6[nbend_bnd_typ] = bend_bnd->dsbend_6[i];

     bend_bnd->eq_bond[nbend_bnd_typ]  = bend_bnd->eq_bond[i];
     bend_bnd->cbond_0[nbend_bnd_typ]  = bend_bnd->cbond_0[i];
     bend_bnd->cbond_1[nbend_bnd_typ]  = bend_bnd->cbond_1[i];
     bend_bnd->cbond_2[nbend_bnd_typ]  = bend_bnd->cbond_2[i];
     bend_bnd->cbond_3[nbend_bnd_typ]  = bend_bnd->cbond_3[i];
     bend_bnd->cbond_4[nbend_bnd_typ]  = bend_bnd->cbond_4[i];
     bend_bnd->cbond_5[nbend_bnd_typ]  = bend_bnd->cbond_5[i];
     bend_bnd->cbond_6[nbend_bnd_typ]  = bend_bnd->cbond_6[i];
     bend_bnd->dcbond_0[nbend_bnd_typ] = bend_bnd->dcbond_0[i];
     bend_bnd->dcbond_1[nbend_bnd_typ] = bend_bnd->dcbond_1[i];
     bend_bnd->dcbond_2[nbend_bnd_typ] = bend_bnd->dcbond_2[i]; 
     bend_bnd->dcbond_3[nbend_bnd_typ] = bend_bnd->dcbond_3[i]; 
     bend_bnd->dcbond_4[nbend_bnd_typ] = bend_bnd->dcbond_4[i]; 
     bend_bnd->dcbond_5[nbend_bnd_typ] = bend_bnd->dcbond_5[i]; 
     bend_bnd->dcbond_6[nbend_bnd_typ] = bend_bnd->dcbond_6[i];
    }/*endif*/
   }/*endfor*/
   bend_bnd->ntyp = nbend_bnd_typ;
   bend->ntyp_pow = nbend_typ;

/*=======================================================================*/
/* III) Realloc the memory                                               */

    bend->eq_pow = (double *)crealloc(&(bend->eq_pow)[1],
                    bend->ntyp_pow*sizeof(double))-1;
    bend->c_0    = (double *)crealloc(&(bend->c_0)[1],
                    bend->ntyp_pow*sizeof(double))-1;
    bend->c_1    = (double *)crealloc(&(bend->c_1)[1],
                    bend->ntyp_pow*sizeof(double))-1;
    bend->c_2    = (double *)crealloc(&(bend->c_2)[1],
                    bend->ntyp_pow*sizeof(double))-1;
    bend->c_3    = (double *)crealloc(&(bend->c_3)[1],
                    bend->ntyp_pow*sizeof(double))-1;
    bend->c_4    = (double *)crealloc(&(bend->c_4)[1],
                    bend->ntyp_pow*sizeof(double))-1;
    bend->c_5    = (double *)crealloc(&(bend->c_5)[1],
                    bend->ntyp_pow*sizeof(double))-1;
    bend->c_6    = (double *)crealloc(&(bend->c_6)[1],
                    bend->ntyp_pow*sizeof(double))-1;
    bend->s_0    = (double *)crealloc(&(bend->s_0)[1],
                    bend->ntyp_pow*sizeof(double))-1;
    bend->s_1    = (double *)crealloc(&(bend->s_1)[1],
                    bend->ntyp_pow*sizeof(double))-1;
    bend->s_2    = (double *)crealloc(&(bend->s_2)[1],
                    bend->ntyp_pow*sizeof(double))-1;
    bend->s_3    = (double *)crealloc(&(bend->s_3)[1],
                    bend->ntyp_pow*sizeof(double))-1;
    bend->s_4    = (double *)crealloc(&(bend->s_4)[1],
                    bend->ntyp_pow*sizeof(double))-1;
    bend->s_5    = (double *)crealloc(&(bend->s_5)[1],
                    bend->ntyp_pow*sizeof(double))-1;
    bend->s_6    = (double *)crealloc(&(bend->s_6)[1],
                    bend->ntyp_pow*sizeof(double))-1;
    bend->dc_0   = (double *)crealloc(&(bend->dc_0)[1],
                    bend->ntyp_pow*sizeof(double))-1;
    bend->dc_1   = (double *)crealloc(&(bend->dc_1)[1],
                    bend->ntyp_pow*sizeof(double))-1;
    bend->dc_2   = (double *)crealloc(&(bend->dc_2)[1],
                    bend->ntyp_pow*sizeof(double))-1;
    bend->dc_3   = (double *)crealloc(&(bend->dc_3)[1],
                    bend->ntyp_pow*sizeof(double))-1;
    bend->dc_4   = (double *)crealloc(&(bend->dc_4)[1],
                    bend->ntyp_pow*sizeof(double))-1;
    bend->dc_5   = (double *)crealloc(&(bend->dc_5)[1],
                    bend->ntyp_pow*sizeof(double))-1;
    bend->dc_6   = (double *)crealloc(&(bend->dc_6)[1],
                    bend->ntyp_pow*sizeof(double))-1;
    bend->ds_0   = (double *)crealloc(&(bend->ds_0)[1],
                    bend->ntyp_pow*sizeof(double))-1;
    bend->ds_1   = (double *)crealloc(&(bend->ds_1)[1],
                    bend->ntyp_pow*sizeof(double))-1;
    bend->ds_2   = (double *)crealloc(&(bend->ds_2)[1],
                    bend->ntyp_pow*sizeof(double))-1;
    bend->ds_3   = (double *)crealloc(&(bend->ds_3)[1],
                    bend->ntyp_pow*sizeof(double))-1;
    bend->ds_4   = (double *)crealloc(&(bend->ds_4)[1],
                    bend->ntyp_pow*sizeof(double))-1;
    bend->ds_5   = (double *)crealloc(&(bend->ds_5)[1],
                    bend->ntyp_pow*sizeof(double))-1;
    bend->ds_6   = (double *)crealloc(&(bend->ds_6)[1],
                    bend->ntyp_pow*sizeof(double))-1;
/*=======================================================================*/
/* III) Realloc more memory                                               */

    bend_bnd->eq_bend = (double *)crealloc(&(bend_bnd->eq_bend)[1],
                    bend_bnd->ntyp*sizeof(double))-1;
    bend_bnd->cbend_0    = (double *)crealloc(&(bend_bnd->cbend_0)[1],
                    bend_bnd->ntyp*sizeof(double))-1;
    bend_bnd->cbend_1    = (double *)crealloc(&(bend_bnd->cbend_1)[1],
                    bend_bnd->ntyp*sizeof(double))-1;
    bend_bnd->cbend_2    = (double *)crealloc(&(bend_bnd->cbend_2)[1],
                    bend_bnd->ntyp*sizeof(double))-1;
    bend_bnd->cbend_3    = (double *)crealloc(&(bend_bnd->cbend_3)[1],
                    bend_bnd->ntyp*sizeof(double))-1;
    bend_bnd->cbend_4    = (double *)crealloc(&(bend_bnd->cbend_4)[1],
                    bend_bnd->ntyp*sizeof(double))-1;
    bend_bnd->cbend_5    = (double *)crealloc(&(bend_bnd->cbend_5)[1],
                    bend_bnd->ntyp*sizeof(double))-1;
    bend_bnd->cbend_6    = (double *)crealloc(&(bend_bnd->cbend_6)[1],
                    bend_bnd->ntyp*sizeof(double))-1;
    bend_bnd->sbend_0    = (double *)crealloc(&(bend_bnd->sbend_0)[1],
                    bend_bnd->ntyp*sizeof(double))-1;
    bend_bnd->sbend_1    = (double *)crealloc(&(bend_bnd->sbend_1)[1],
                    bend_bnd->ntyp*sizeof(double))-1;
    bend_bnd->sbend_2    = (double *)crealloc(&(bend_bnd->sbend_2)[1],
                    bend_bnd->ntyp*sizeof(double))-1;
    bend_bnd->sbend_3    = (double *)crealloc(&(bend_bnd->sbend_3)[1],
                    bend_bnd->ntyp*sizeof(double))-1;
    bend_bnd->sbend_4    = (double *)crealloc(&(bend_bnd->sbend_4)[1],
                    bend_bnd->ntyp*sizeof(double))-1;
    bend_bnd->sbend_5    = (double *)crealloc(&(bend_bnd->sbend_5)[1],
                    bend_bnd->ntyp*sizeof(double))-1;
    bend_bnd->sbend_6    = (double *)crealloc(&(bend_bnd->sbend_6)[1],
                    bend_bnd->ntyp*sizeof(double))-1;
    bend_bnd->dcbend_0   = (double *)crealloc(&(bend_bnd->dcbend_0)[1],
                    bend_bnd->ntyp*sizeof(double))-1;
    bend_bnd->dcbend_1   = (double *)crealloc(&(bend_bnd->dcbend_1)[1],
                    bend_bnd->ntyp*sizeof(double))-1;
    bend_bnd->dcbend_2   = (double *)crealloc(&(bend_bnd->dcbend_2)[1],
                    bend_bnd->ntyp*sizeof(double))-1;
    bend_bnd->dcbend_3   = (double *)crealloc(&(bend_bnd->dcbend_3)[1],
                    bend_bnd->ntyp*sizeof(double))-1;
    bend_bnd->dcbend_4   = (double *)crealloc(&(bend_bnd->dcbend_4)[1],
                    bend_bnd->ntyp*sizeof(double))-1;
    bend_bnd->dcbend_5   = (double *)crealloc(&(bend_bnd->dcbend_5)[1],
                    bend_bnd->ntyp*sizeof(double))-1;
    bend_bnd->dcbend_6   = (double *)crealloc(&(bend_bnd->dcbend_6)[1],
                    bend_bnd->ntyp*sizeof(double))-1;
    bend_bnd->dsbend_0   = (double *)crealloc(&(bend_bnd->dsbend_0)[1],
                    bend_bnd->ntyp*sizeof(double))-1;
    bend_bnd->dsbend_1   = (double *)crealloc(&(bend_bnd->dsbend_1)[1],
                    bend_bnd->ntyp*sizeof(double))-1;
    bend_bnd->dsbend_2   = (double *)crealloc(&(bend_bnd->dsbend_2)[1],
                    bend_bnd->ntyp*sizeof(double))-1;
    bend_bnd->dsbend_3   = (double *)crealloc(&(bend_bnd->dsbend_3)[1],
                    bend_bnd->ntyp*sizeof(double))-1;
    bend_bnd->dsbend_4   = (double *)crealloc(&(bend_bnd->dsbend_4)[1],
                    bend_bnd->ntyp*sizeof(double))-1;
    bend_bnd->dsbend_5   = (double *)crealloc(&(bend_bnd->dsbend_5)[1],
                    bend_bnd->ntyp*sizeof(double))-1;
    bend_bnd->dsbend_6   = (double *)crealloc(&(bend_bnd->dsbend_6)[1],
                    bend_bnd->ntyp*sizeof(double))-1;

    bend_bnd->nbend_bnd_typ_mall = bend_bnd->ntyp;
    bend_bnd->eq_bond = (double *)crealloc(&(bend_bnd->eq_bond)[1],
                    bend_bnd->ntyp*sizeof(double))-1;
    bend_bnd->cbond_0    = (double *)crealloc(&(bend_bnd->cbond_0)[1],
                    bend_bnd->ntyp*sizeof(double))-1;
    bend_bnd->cbond_1    = (double *)crealloc(&(bend_bnd->cbond_1)[1],
                    bend_bnd->ntyp*sizeof(double))-1;
    bend_bnd->cbond_2    = (double *)crealloc(&(bend_bnd->cbond_2)[1],
                    bend_bnd->ntyp*sizeof(double))-1;
    bend_bnd->cbond_3    = (double *)crealloc(&(bend_bnd->cbond_3)[1],
                    bend_bnd->ntyp*sizeof(double))-1;
    bend_bnd->cbond_4    = (double *)crealloc(&(bend_bnd->cbond_4)[1],
                    bend_bnd->ntyp*sizeof(double))-1;
    bend_bnd->cbond_5    = (double *)crealloc(&(bend_bnd->cbond_5)[1],
                    bend_bnd->ntyp*sizeof(double))-1;
    bend_bnd->cbond_6    = (double *)crealloc(&(bend_bnd->cbond_6)[1],
                    bend_bnd->ntyp*sizeof(double))-1;
    bend_bnd->dcbond_0   = (double *)crealloc(&(bend_bnd->dcbond_0)[1],
                    bend_bnd->ntyp*sizeof(double))-1;
    bend_bnd->dcbond_1   = (double *)crealloc(&(bend_bnd->dcbond_1)[1],
                    bend_bnd->ntyp*sizeof(double))-1;
    bend_bnd->dcbond_2   = (double *)crealloc(&(bend_bnd->dcbond_2)[1],
                    bend_bnd->ntyp*sizeof(double))-1;
    bend_bnd->dcbond_3   = (double *)crealloc(&(bend_bnd->dcbond_3)[1],
                    bend_bnd->ntyp*sizeof(double))-1;
    bend_bnd->dcbond_4   = (double *)crealloc(&(bend_bnd->dcbond_4)[1],
                    bend_bnd->ntyp*sizeof(double))-1;
    bend_bnd->dcbond_5   = (double *)crealloc(&(bend_bnd->dcbond_5)[1],
                    bend_bnd->ntyp*sizeof(double))-1;
    bend_bnd->dcbond_6   = (double *)crealloc(&(bend_bnd->dcbond_6)[1],
                    bend_bnd->ntyp*sizeof(double))-1;
/*-----------------------------------------------------------------------*/
}  /*end routine*/
/*==========================================================================*/






