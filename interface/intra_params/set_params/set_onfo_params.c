/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
/*  set_onfo_params:set up intramolecular interaction                       */
/*                    key word dictionary                                   */
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

void set_onfo_params(DICT_WORD *intra_dict,int num_intra_dict,
		     char *fun_key,char *file_name,int jmol_typ,
		     CLATOMS_INFO *clatoms_info,ATOMMAPS *atommaps,
		     ONFO *onfo,NULL_INTER_PARSE *null_inter_parse,
		     BUILD_INTRA *build_intra, int iresidue, int ires_off)

/*==========================================================================*/
{/*begin routine */
  int num,index,ifound,igo;
  int itype1,itype2;
  int iatm_ind1,iatm_ind2;
  int imask1,imask2;
  int i,itype;

/*=======================================================================*/
/* I) Check for missing key words*/
  for(i=1;i<=2;i++){
    if(intra_dict[i].iuset==0 && intra_dict[i].key_type==1){
      keyword_miss(intra_dict,file_name,fun_key,i);}
  }   /*endfor*/
  /*=======================================================================*/
  /* II) Fill the dictionary with words */
  /*-----------------------------------------------------------------------*/
  /*  1) \atom1{}    */
  index = 1;
  sscanf(intra_dict[1].keyarg,"%d",&num);
  iatm_ind1 = num;
  if(iatm_ind1>build_intra->natmind_1res_now||iatm_ind1<0){
    keyarg_barf(intra_dict,file_name,fun_key,index);}
  imask1 = build_intra->mask_atm[iatm_ind1];  
  if(imask1>0)iatm_ind1 = build_intra->index_atm[iatm_ind1];
  /*------------------------------------------------------------------------*/
  /*  2) \atom2{}    */
  index = 2;
  sscanf(intra_dict[2].keyarg,"%d",&num);
  iatm_ind2 = num;
  if(iatm_ind2>build_intra->natmind_1res_now||iatm_ind2<0){
    keyarg_barf(intra_dict,file_name,fun_key,index);}
  imask2 = build_intra->mask_atm[iatm_ind2];  
  if(imask2>0)iatm_ind2 = build_intra->index_atm[iatm_ind2];
  /*-----------------------------------------------------------------------*/
  /*  5) \modifier{} */
  index  = 5;
  ifound = 0;
  if(strcasecmp(intra_dict[5].keyarg,"on")==0) {ifound = 1;}
  if(strcasecmp(intra_dict[5].keyarg,"off")==0){ifound = 2;}
  if(ifound==0){
    keyarg_barf(intra_dict,file_name,fun_key,index);
  }
  /*-----------------------------------------------------------------------*/
  /*  6) bond type    */
  igo = imask1*imask2;
  if(igo==1){
    itype1 = atommaps->iatm_atm_typ[(clatoms_info->natm_tot+iatm_ind1)];
    itype2 = atommaps->iatm_atm_typ[(clatoms_info->natm_tot+iatm_ind2)];
    strcpy(build_intra->confo_typ_now->atm1,atommaps->atm_typ[itype1]);
    strcpy(build_intra->confo_typ_now->atm2,atommaps->atm_typ[itype2]);
    strcpy(build_intra->confo_typ_now->label,intra_dict[6].keyarg);
  }/*endif*/
  /*======================================================================*/
  /* III) Spread onfos */
  if(ifound==1&&igo==1){
    /*---------------------------------------------------------------------*/
    /* A) Add more space */
    if(onfo->num+1 > build_intra->nonfo_max){
      build_intra->nonfo_max += NMEM_MIN;
      onfo->j1 =(int *) crealloc(&(onfo->j1)[1],
				 build_intra->nonfo_max*sizeof(int))-1;
      onfo->j2 =(int *) crealloc(&(onfo->j2)[1],
				 build_intra->nonfo_max*sizeof(int))-1;
      onfo->jtyp =(int *) crealloc(&(onfo->jtyp)[1],
				   build_intra->nonfo_max*sizeof(int))-1;
    }/*endif*/
    /*--------------------------------------------------------------------*/
    /* C) Check type */
    itype = (onfo->ntyp)+1;
    for(i=1;i<=onfo->ntyp;i++){
      if((strcasecmp(build_intra->confo_typ[i].atm1,
		     build_intra->confo_typ_now->atm1)==0)
         &&(strcasecmp(build_intra->confo_typ[i].atm2,
		       build_intra->confo_typ_now->atm2)==0)
	 &&(strcasecmp(build_intra->confo_typ[i].label,
		       build_intra->confo_typ_now->label)==0)) {itype=i;}
      if((strcasecmp(build_intra->confo_typ[i].atm1,
		     build_intra->confo_typ_now->atm2)==0)
         &&(strcasecmp(build_intra->confo_typ[i].atm2,
		       build_intra->confo_typ_now->atm1)==0)
	 &&(strcasecmp(build_intra->confo_typ[i].label,
		       build_intra->confo_typ_now->label)==0)) {itype=i;}
    }/*endfor*/
    /*---------------------------------------------------------------------*/
    /* D) Add space */
    if(itype>build_intra->nonfo_typ_max){
      build_intra->nonfo_typ_max += NMEM_MIN;
      build_intra->confo_typ      = (CBOND *) 
	crealloc(&(build_intra->confo_typ)[1],
		 build_intra->nonfo_typ_max*sizeof(CBOND))-1;
    }/*endif*/
    /*---------------------------------------------------------------------*/
    /* E) Add a type */
    if(itype==(onfo->ntyp)+1){
      onfo->ntyp+=1;
      strcpy(build_intra->confo_typ[itype].atm1,
	     build_intra->confo_typ_now->atm1);
      strcpy(build_intra->confo_typ[itype].atm2,
	     build_intra->confo_typ_now->atm2);
      strcpy(build_intra->confo_typ[itype].label,
	     build_intra->confo_typ_now->label);
    }/*endif*/
    /*---------------------------------------------------------------------*/
    /* B) Spread */

    onfo->num += 1;
    onfo->j1[onfo->num] = iatm_ind1 + clatoms_info->natm_tot;
    onfo->j2[onfo->num] = iatm_ind2 + clatoms_info->natm_tot;
    onfo->jtyp[onfo->num] = itype;
  }/*endif*/
  /*======================================================================*/
  /* V) Spread nul onfos */
  if(ifound==3&&igo==1){
    /*--------------------------------------------------------------------*/
    /* A) Add more space */
    if(null_inter_parse->nonfo_nul+1 > build_intra->nonfo_nul_max){
      build_intra->nonfo_nul_max += NMEM_MIN;
      null_inter_parse->jonfo1_nul     = 
	(int *) crealloc(&(null_inter_parse->jonfo1_nul)[1],
			 build_intra->nonfo_nul_max*sizeof(int))-1;
      null_inter_parse->jonfo2_nul     = 
	(int *) crealloc(&(null_inter_parse->jonfo2_nul)[1],
			 build_intra->nonfo_nul_max*sizeof(int))-1;
    }/*endif*/
    /*--------------------------------------------------------------------*/
    /* B) Spread */
    null_inter_parse->nonfo_nul += 1;
    null_inter_parse->jonfo1_nul[null_inter_parse->nonfo_nul]
                                         = iatm_ind1 + clatoms_info->natm_tot;
    null_inter_parse->jonfo2_nul[null_inter_parse->nonfo_nul] 
                                         = iatm_ind2 + clatoms_info->natm_tot;
  }/*endif*/
  /*----------------------------------------------------------------------*/
}/*end routine*/
/*==========================================================================*/





