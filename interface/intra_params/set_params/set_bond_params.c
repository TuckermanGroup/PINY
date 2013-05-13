/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
/*  set_bond_params:set up intramolecular interaction                       */
/*                    key word dictionary                                   */
/*==========================================================================*/

#include "standard_include.h"
#include "../typ_defs/typedefs_par.h"
#include "../typ_defs/typedefs_gen.h"
#include "../typ_defs/typedefs_class.h"
#include "../typ_defs/typedefs_bnd.h"
#include "../proto_defs/proto_intra_params_local.h"
#include "../proto_defs/proto_friend_lib_entry.h"
#include "../proto_defs/proto_handle_entry.h"

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

 void set_bond_params(DICT_WORD *intra_dict,int num_intra_dict,
		      char *fun_key,char *file_name,int jmol_typ,
		      CLATOMS_INFO *clatoms_info,ATOMMAPS *atommaps,
		      BOND *bond,NULL_INTER_PARSE *null_inter_parse,
		      BUILD_INTRA *build_intra, int iresidue, int ires_off,
                      int mol_hydrog_con_opt)

/*==========================================================================*/
{ /*begin routine*/

  int num,index,ifound,igo;
  int itype1,itype2;
  int iatm_ind1,iatm_ind2;
  int imask1,imask2;
  int i,itype;
  int mass_now1,mass_now2;

/*=======================================================================*/
/* I) Check for missing key words*/

  for(i=1;i<=2;i++){
    if(intra_dict[i].iuset==0 && intra_dict[i].key_type==1){
      keyword_miss(intra_dict,file_name,fun_key,i);
    }
  }/*endfor*/

/*=======================================================================*/
/* II) Fill the dictionary with words */
/*-----------------------------------------------------------------------*/
  /*  1) \atom1{}    */
  index = 1;
  sscanf(intra_dict[1].keyarg,"%d",&num);
  iatm_ind1 = num;
  if(iatm_ind1>build_intra->natmind_1res_now||iatm_ind1<0){
    keyarg_barf(intra_dict,file_name,fun_key,index);
  }
  imask1 = build_intra->mask_atm[iatm_ind1];  
  if(imask1>0)iatm_ind1 = build_intra->index_atm[iatm_ind1];
  /*-----------------------------------------------------------------------*/
  /*  2) \atom2{}    */
  index = 2;
  sscanf(intra_dict[2].keyarg,"%d",&num);
  iatm_ind2 = num;
  if(iatm_ind2>build_intra->natmind_1res_now||iatm_ind2<0){
            keyarg_barf(intra_dict,file_name,fun_key,index);}
  imask2 = build_intra->mask_atm[iatm_ind2];  
  if(imask2>0)iatm_ind2 = build_intra->index_atm[iatm_ind2];
  /*----------------------------------------------------------------------*/
  /*  5) \modifier{} */
  index  = 5;
  ifound = 0;
  if(strcasecmp(intra_dict[5].keyarg,"on")==0) {ifound = 1;}
  if(strcasecmp(intra_dict[5].keyarg,"con")==0){ifound = 2;}
  if(strcasecmp(intra_dict[5].keyarg,"off")==0){ifound = 3;}
  if(ifound==0){
    keyarg_barf(intra_dict,file_name,fun_key,index);
  }
  if(ifound == 2 && clatoms_info->pi_beads>1){ 
        printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
        printf("Constraints not implemented under  \n");
        printf(" path integral molecular dynamics. \n");
        printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
        fflush(stdout);
        exit(1);
  }/*endif*/
  /*-----------------------------------------------------------------------*/
  /*  6) bond type and correct modifier if hydrog_con_opt is on  */
  igo = imask1*imask2;
  if(igo==1){
    itype1 = atommaps->iatm_atm_typ[(clatoms_info->natm_tot+iatm_ind1)];
    itype2 = atommaps->iatm_atm_typ[(clatoms_info->natm_tot+iatm_ind2)];
    strcpy(build_intra->cbond_typ_now->atm1,atommaps->atm_typ[itype1]);
    strcpy(build_intra->cbond_typ_now->atm2,atommaps->atm_typ[itype2]);
    strcpy(build_intra->cbond_typ_now->label,intra_dict[6].keyarg);
    mass_now1 =
            (int)NINT(clatoms_info->mass[(clatoms_info->natm_tot+iatm_ind1)]);
    mass_now2 =
            (int)NINT(clatoms_info->mass[(clatoms_info->natm_tot+iatm_ind2)]);
  /* check to see if this is any atom-H bond                      */
    if((mass_now1<=2) || (mass_now2<=2)){
      if(mol_hydrog_con_opt==1){ifound=2;}
    }/*endif*/
  /* check to see if this is a polar atom-H bond (polar=N,O,S)    */
    if((mass_now1<=2) && ((mass_now2==14) || (mass_now2==16)
                                          || (mass_now2==32))){
      if(mol_hydrog_con_opt==2){ifound=2;}
    }/*endif*/
    if((mass_now2<=2) && ((mass_now1==14) || (mass_now1==16)
                                          || (mass_now1==32))){
      if(mol_hydrog_con_opt==2){ifound=2;}
    }/*endif*/
  }/*endif*/
  /*======================================================================*/
  /* III) Spread Power bonds */
  if(ifound==1&&igo==1){
    /*---------------------------------------------------------------------*/
    /* A) Add more space */
    if(bond->npow+1 > build_intra->nbond_pow_max){
      build_intra->nbond_pow_max += NMEM_MIN;
      bond->j1_pow = (int *)crealloc(&(bond->j1_pow)[1],
				     build_intra->nbond_pow_max*sizeof(int))-1;
      bond->j2_pow = (int *)crealloc(&(bond->j2_pow)[1],
				     build_intra->nbond_pow_max*sizeof(int))-1;
      bond->jtyp_pow=(int *)crealloc(&(bond->jtyp_pow)[1],
				     build_intra->nbond_pow_max*sizeof(int))-1;
    }/*endif*/
    /*----------------------------------------------------------------------*/
    /* B) Check type */
    itype = (bond->ntyp_pow)+1;
    for(i=1;i<=bond->ntyp_pow;i++){
      if((strcasecmp(build_intra->cbond_typ_pow[i].atm1,
		     build_intra->cbond_typ_now->atm1)==0)
         &&(strcasecmp(build_intra->cbond_typ_pow[i].atm2,
		       build_intra->cbond_typ_now->atm2)==0)
	 &&(strcasecmp(build_intra->cbond_typ_pow[i].label,
		       build_intra->cbond_typ_now->label)==0)) {itype=i;}
      if((strcasecmp(build_intra->cbond_typ_pow[i].atm1,
		     build_intra->cbond_typ_now->atm2)==0)
         &&(strcasecmp(build_intra->cbond_typ_pow[i].atm2,
		       build_intra->cbond_typ_now->atm1)==0)
	 &&(strcasecmp(build_intra->cbond_typ_pow[i].label,
		       build_intra->cbond_typ_now->label)==0)) {itype=i;}
    }/*endfor*/
    /*----------------------------------------------------------------------*/
    /* C) Add space */
    if(itype>build_intra->nbond_typ_pow_max){
      build_intra->nbond_typ_pow_max += NMEM_MIN;
      build_intra->cbond_typ_pow      = (CBOND *) 
	crealloc(&build_intra->cbond_typ_pow[1],
		 build_intra->nbond_typ_pow_max*sizeof(CBOND))-1;
    } /*endif*/
    /*----------------------------------------------------------------------*/
    /* D) Add a type */
    if(itype==(bond->ntyp_pow)+1){
      bond->ntyp_pow+=1;
      strcpy(build_intra->cbond_typ_pow[itype].atm1,
	     build_intra->cbond_typ_now->atm1);
      strcpy(build_intra->cbond_typ_pow[itype].atm2,
	     build_intra->cbond_typ_now->atm2);
      strcpy(build_intra->cbond_typ_pow[itype].label,
	     build_intra->cbond_typ_now->label);
    }/*endif*/
    /*----------------------------------------------------------------------*/
    /* E) Assign */

    bond->npow += 1;
    bond->j1_pow[bond->npow] = iatm_ind1 + clatoms_info->natm_tot;
    bond->j2_pow[bond->npow] = iatm_ind2 + clatoms_info->natm_tot;
    bond->jtyp_pow[bond->npow] = itype;
  }/*endif*/
  /*=======================================================================*/
  /* IV) Assign Con bonds */
  if(ifound==2&&igo==1){
    ires_off = atommaps->jres_jmol_typ_strt[jmol_typ]-1;
    atommaps->nfree_jres_jmol_typ[ires_off+iresidue] -= 1;
    atommaps->icons_jres_jmol_typ[ires_off+iresidue] = 1;
    /*---------------------------------------------------------------------*/
    /* A) Add space */
    if(bond->ncon+1 > build_intra->nbond_con_max){
      build_intra->nbond_con_max += NMEM_MIN;
      bond->j1_con = (int *)crealloc(&(bond->j1_con)[1],
				     build_intra->nbond_con_max*sizeof(int))-1;
      bond->j2_con = (int *)crealloc(&(bond->j2_con)[1],
				     build_intra->nbond_con_max*sizeof(int))-1;
      bond->jtyp_con=(int *)crealloc(&(bond->jtyp_con)[1],
				     build_intra->nbond_con_max*sizeof(int))-1;
    }/*endif*/
    /*---------------------------------------------------------------------*/
    /* B) Check type */
    itype = (bond->ntyp_con)+1;
    for(i=1;i<=bond->ntyp_con;i++){
      if((strcasecmp(build_intra->cbond_typ_con[i].atm1,
		     build_intra->cbond_typ_now->atm1)==0)
         &&(strcasecmp(build_intra->cbond_typ_con[i].atm2,
		       build_intra->cbond_typ_now->atm2)==0)
	 &&(strcasecmp(build_intra->cbond_typ_con[i].label,
		       build_intra->cbond_typ_now->label)==0)) {itype=i;}
      if((strcasecmp(build_intra->cbond_typ_con[i].atm1,
		     build_intra->cbond_typ_now->atm2)==0)
         &&(strcasecmp(build_intra->cbond_typ_con[i].atm2,
		       build_intra->cbond_typ_now->atm1)==0)
	 &&(strcasecmp(build_intra->cbond_typ_con[i].label,
		       build_intra->cbond_typ_now->label)==0)) {itype=i;}
    }/*endfor*/
    /*---------------------------------------------------------------------*/
    /* C) Add space */
    if(itype>build_intra->nbond_typ_con_max){
      build_intra->nbond_typ_con_max += NMEM_MIN;
      build_intra->cbond_typ_con      = (CBOND *) 
	crealloc(&(build_intra->cbond_typ_con)[1],
		 build_intra->nbond_typ_con_max*sizeof(CBOND))-1;
    }/*endif*/
    /*---------------------------------------------------------------------*/
    /* D) Add a type */
    if(itype==(bond->ntyp_con)+1){
      bond->ntyp_con+=1;
      strcpy(build_intra->cbond_typ_con[itype].atm1,
	     build_intra->cbond_typ_now->atm1);
      strcpy(build_intra->cbond_typ_con[itype].atm2,
	     build_intra->cbond_typ_now->atm2);
      strcpy(build_intra->cbond_typ_con[itype].label,
	     build_intra->cbond_typ_now->label);
    }/*endif*/
    /*---------------------------------------------------------------------*/
    /* E) Assign */

    bond->ncon += 1;
    bond->j1_con[bond->ncon] = iatm_ind1 + clatoms_info->natm_tot;
    bond->j2_con[bond->ncon] = iatm_ind2 + clatoms_info->natm_tot;
    bond->jtyp_con[bond->ncon] = itype;

    /*---------------------------------------------------------------------*/
    /* F) Error */
    if((atommaps->ighost_flag[bond->j1_con[bond->ncon]]!=0)||
       (atommaps->ighost_flag[bond->j2_con[bond->ncon]]!=0)){
        printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
        printf("Ghost atoms are not permitted in constrained bonds\n");
        printf("in molecule index %d in residue index %d in file %s\n",
                jmol_typ,iresidue,file_name);
        printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
        fflush(stdout);
        exit(1);
    }/*endif*/      

  }/*endif*/
  /*======================================================================*/
  /* V) Spread nul bonds */
  if(ifound==3&&igo==1){
    /*---------------------------------------------------------------------*/
    /* A) Add more space */
    if((null_inter_parse->nbond_nul)+1> build_intra->nbond_nul_max) {
      build_intra->nbond_nul_max += NMEM_MIN;
      null_inter_parse->jbond1_nul     = 
	(int *) crealloc(&(null_inter_parse->jbond1_nul)[1],
			 build_intra->nbond_nul_max*sizeof(int))-1;
      null_inter_parse->jbond2_nul     = 
	(int *) crealloc(&(null_inter_parse->jbond2_nul)[1],
			 build_intra->nbond_nul_max*sizeof(int))-1;
    }  /*endif*/
    /*---------------------------------------------------------------------*/
    /* B) Spread */
    null_inter_parse->nbond_nul += 1;
    null_inter_parse->jbond1_nul[null_inter_parse->nbond_nul] 
      = iatm_ind1 + clatoms_info->natm_tot;
    null_inter_parse->jbond2_nul[null_inter_parse->nbond_nul] 
      = iatm_ind2 + clatoms_info->natm_tot;
  }/*endfor*/
  /*-----------------------------------------------------------------------*/
}/*end routine*/
/*==========================================================================*/






