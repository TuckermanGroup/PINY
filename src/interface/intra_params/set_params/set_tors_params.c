/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
/*  set_tors_params:set up intramolecular interaction                       */
/*                    key word dictionary                                   */
/*==========================================================================*/

#include "standard_include.h"
#include "../typ_defs/typedefs_class.h"
#include "../typ_defs/typedefs_par.h"
#include "../typ_defs/typedefs_bnd.h"
#include "../proto_defs/proto_intra_params_local.h"
#include "../proto_defs/proto_handle_entry.h"
#include "../proto_defs/proto_friend_lib_entry.h"

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void set_tors_params(DICT_WORD intra_dict[],int num_intra_dict,
		     char *fun_key,char *file_name,int jmol_typ,
		     CLATOMS_INFO *clatoms_info,ATOMMAPS *atommaps,
		     TORS *tors,NULL_INTER_PARSE *null_inter_parse,
		     BUILD_INTRA *build_intra, int iresidue, int ires_off)

/*==========================================================================*/
{
  int num,index,ifound,igo;
  int itype1,itype2,itype3,itype4;
  int iatm_ind1,iatm_ind2,iatm_ind3,iatm_ind4;
  int imask1,imask2,imask3,imask4;
  int i,itype;

  /*=======================================================================*/
  /* I) Check for missing key words*/
  for(i=1;i<=num_intra_dict;i++){
    if(intra_dict[i].iuset==0 && intra_dict[i].key_type==1){
      keyword_miss(intra_dict,file_name,fun_key,i);
    }
  }
  /*======================================================================*/
  /* II) Fill the dictionary with words */
  /*----------------------------------------------------------------------*/
  /*  1) \atom1{}    */
  index = 1;
  sscanf(intra_dict[1].keyarg,"%d",&num);
  iatm_ind1 = num;
  if(iatm_ind1>build_intra->natmind_1res_now||iatm_ind1<0){
            keyarg_barf(intra_dict,file_name,fun_key,index);}
  imask1 = build_intra->mask_atm[iatm_ind1];  
  if(imask1>0)iatm_ind1 = build_intra->index_atm[iatm_ind1];

  /*----------------------------------------------------------------------*/
  /*  2) \atom2{}    */
  index = 2;
  sscanf(intra_dict[2].keyarg,"%d",&num);
  iatm_ind2 = num;
  if(iatm_ind2>build_intra->natmind_1res_now||iatm_ind2<0){
            keyarg_barf(intra_dict,file_name,fun_key,index);}
  imask2 = build_intra->mask_atm[iatm_ind2];  
  if(imask2>0)iatm_ind2 = build_intra->index_atm[iatm_ind2];
  /*------------------------------------------------------------------------*/
  /*  3) \atom3{}    */
  index = 3;
  sscanf(intra_dict[3].keyarg,"%d",&num);
  iatm_ind3 = num;
  if(iatm_ind3>build_intra->natmind_1res_now||iatm_ind3<0){
            keyarg_barf(intra_dict,file_name,fun_key,index);}
  imask3 = build_intra->mask_atm[iatm_ind3];  
  if(imask3>0)iatm_ind3 = build_intra->index_atm[iatm_ind3];

  /*----------------------------------------------------------------------*/
  /*  4) \atom4{}    */
  index = 4;
  sscanf(intra_dict[4].keyarg,"%d",&num);
  iatm_ind4 = num;
  if(iatm_ind4>build_intra->natmind_1res_now||iatm_ind4<0){
            keyarg_barf(intra_dict,file_name,fun_key,index);}
  imask4 = build_intra->mask_atm[iatm_ind4];  
  if(imask4>0)iatm_ind4 = build_intra->index_atm[iatm_ind4];

  /*---------------------------------------------------------------------*/
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
        printf(" path integral dynamics.            \n");
        printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
        fflush(stdout);
        exit(1);
  }/*endif*/
  /*-----------------------------------------------------------------------*/
  /*  6) tors type    */
  igo = imask1*imask2*imask3*imask4;
  if(igo==1){
    itype1 = atommaps->iatm_atm_typ[(clatoms_info->natm_tot+iatm_ind1)];
    itype2 = atommaps->iatm_atm_typ[(clatoms_info->natm_tot+iatm_ind2)];
    itype3 = atommaps->iatm_atm_typ[(clatoms_info->natm_tot+iatm_ind3)];
    itype4 = atommaps->iatm_atm_typ[(clatoms_info->natm_tot+iatm_ind4)];
    strcpy(build_intra->ctors_typ_now->atm1,atommaps->atm_typ[itype1]);
    strcpy(build_intra->ctors_typ_now->atm2,atommaps->atm_typ[itype2]);
    strcpy(build_intra->ctors_typ_now->atm3,atommaps->atm_typ[itype3]);
    strcpy(build_intra->ctors_typ_now->atm4,atommaps->atm_typ[itype4]);
    strcpy(build_intra->ctors_typ_now->label,intra_dict[6].keyarg);  
    if(strcasecmp(build_intra->ctors_typ_now->label,"improper")==0){
      tors->nimpr++;
    }
  }/*endif*/

  /*=======================================================================*/
  /* III) Spread Power torss */
  /*---------------------------------------------------------------------*/
  if(ifound==1&&igo==1){
    /* A) Add more space */
    if((tors->npow)+1 > build_intra->ntors_pow_max){
      build_intra->ntors_pow_max += NMEM_MIN;
      tors->j1_pow =(int*)crealloc(&(tors->j1_pow)[1],
				   build_intra->ntors_pow_max*sizeof(int))-1;
      tors->j2_pow =(int*)crealloc(&(tors->j2_pow)[1],
				   build_intra->ntors_pow_max*sizeof(int))-1;
      tors->j3_pow =(int*)crealloc(&(tors->j3_pow)[1],
				   build_intra->ntors_pow_max*sizeof(int))-1;
      tors->j4_pow =(int*)crealloc(&(tors->j4_pow)[1],
				   build_intra->ntors_pow_max*sizeof(int))-1;
      tors->jtyp_pow=(int*)crealloc(&(tors->jtyp_pow)[1],
				    build_intra->ntors_pow_max*sizeof(int))-1;
    }/*endif*/
    
    /*----------------------------------------------------------------------*/
    /* C) Check type */
    
    itype = tors->ntyp_pow+1;
    for(i=1;i<=tors->ntyp_pow;i++){
      if((strcasecmp(build_intra->ctors_typ_pow[i].atm1,
		     build_intra->ctors_typ_now->atm1)==0)
         &&(strcasecmp(build_intra->ctors_typ_pow[i].atm2,
		       build_intra->ctors_typ_now->atm2)==0)
         &&(strcasecmp(build_intra->ctors_typ_pow[i].atm3,
		       build_intra->ctors_typ_now->atm3)==0)
         &&(strcasecmp(build_intra->ctors_typ_pow[i].atm4,
		       build_intra->ctors_typ_now->atm4)==0)
	 &&(strcasecmp(build_intra->ctors_typ_pow[i].label,
		       build_intra->ctors_typ_now->label)==0)) {itype=i;}
      
      if((strcasecmp(build_intra->ctors_typ_pow[i].atm1,
		     build_intra->ctors_typ_now->atm4)==0)
         &&(strcasecmp(build_intra->ctors_typ_pow[i].atm2,
		       build_intra->ctors_typ_now->atm3)==0)
         &&(strcasecmp(build_intra->ctors_typ_pow[i].atm3,
		       build_intra->ctors_typ_now->atm2)==0)
         &&(strcasecmp(build_intra->ctors_typ_pow[i].atm4,
		       build_intra->ctors_typ_now->atm1)==0)
	 &&(strcasecmp(build_intra->ctors_typ_pow[i].label,
		       build_intra->ctors_typ_now->label)==0)) {itype=i;}
    }     /*endfor*/
    /*--------------------------------------------------------------------*/
    /* D) Add space */
    if(itype>build_intra->ntors_typ_pow_max){
      build_intra->ntors_typ_pow_max += NMEM_MIN;
      build_intra->ctors_typ_pow = 
	(CTORS *) crealloc(&(build_intra->ctors_typ_pow)[1],
			   build_intra->ntors_typ_pow_max*sizeof(CTORS))-1;
    }
    /*--------------------------------------------------------------------*/
    /* E) Add a type */
    if(itype==tors->ntyp_pow+1){
      tors->ntyp_pow++;
      strcpy(build_intra->ctors_typ_pow[itype].atm1,
	     build_intra->ctors_typ_now->atm1);
       strcpy(build_intra->ctors_typ_pow[itype].atm2,
	      build_intra->ctors_typ_now->atm2);
       strcpy(build_intra->ctors_typ_pow[itype].atm3,
	      build_intra->ctors_typ_now->atm3);
       strcpy(build_intra->ctors_typ_pow[itype].atm4,
	      build_intra->ctors_typ_now->atm4);
       strcpy(build_intra->ctors_typ_pow[itype].label,
	      build_intra->ctors_typ_now->label);
    }/*endif*/
    /*---------------------------------------------------------------------*/
    /* B) Spread */

    tors->npow++;
    tors->j1_pow[tors->npow] = iatm_ind1 + clatoms_info->natm_tot;
    tors->j2_pow[tors->npow] = iatm_ind2 + clatoms_info->natm_tot;
    tors->j3_pow[tors->npow] = iatm_ind3 + clatoms_info->natm_tot;
    tors->j4_pow[tors->npow] = iatm_ind4 + clatoms_info->natm_tot;
    tors->jtyp_pow[tors->npow] = itype;
  }  /*endif*/

  /*======================================================================*/
  /* IV) Spread Con torss */
  if(ifound==2&&igo==1){
    ires_off = atommaps->jres_jmol_typ_strt[jmol_typ]-1;
    atommaps->nfree_jres_jmol_typ[ires_off+iresidue] -= 1;
    atommaps->icons_jres_jmol_typ[ires_off+iresidue] = 1;
      /*---------------------------------------------------------------------*/
     /* A) Add more space */
     if(tors->ncon+1 > build_intra->ntors_con_max){
       build_intra->ntors_con_max += NMEM_MIN;
       tors->j1_con     = (int *) crealloc(&(tors->j1_con)[1],
			  build_intra->ntors_con_max*sizeof(int))-1;
       tors->j2_con     = (int *) crealloc(&(tors->j2_con)[1],
			  build_intra->ntors_con_max*sizeof(int))-1;
       tors->j3_con     = (int *) crealloc(&(tors->j3_con)[1],
			  build_intra->ntors_con_max*sizeof(int))-1;
       tors->j4_con     =  (int *) crealloc(&(tors->j4_con)[1],
			  build_intra->ntors_con_max*sizeof(int))-1;
       tors->jtyp_con     = (int *) crealloc(&(tors->jtyp_con)[1],
			  build_intra->ntors_con_max*sizeof(int))-1;
     }     /*endif*/
     /*--------------------------------------------------------------------*/
     /* C) Check type */
     itype = tors->ntyp_con+1;
     for(i=1;i<= tors->ntyp_con;i++){
       if((strcasecmp(build_intra->ctors_typ_con[i].atm1,
		      build_intra->ctors_typ_now->atm1)==0)
	  &&(strcasecmp(build_intra->ctors_typ_con[i].atm2,
			build_intra->ctors_typ_now->atm2)==0)
	  &&(strcasecmp(build_intra->ctors_typ_con[i].atm3,
			build_intra->ctors_typ_now->atm3)==0)
	  &&(strcasecmp(build_intra->ctors_typ_con[i].atm4,
			build_intra->ctors_typ_now->atm4)==0)
	  &&(strcasecmp(build_intra->ctors_typ_con[i].label,
			build_intra->ctors_typ_now->label)==0)) {itype=i;}
       
       if((strcasecmp(build_intra->ctors_typ_con[i].atm1,
		      build_intra->ctors_typ_now->atm4)==0)
	  &&(strcasecmp(build_intra->ctors_typ_con[i].atm2,
			build_intra->ctors_typ_now->atm3)==0)
	  &&(strcasecmp(build_intra->ctors_typ_con[i].atm3,
			build_intra->ctors_typ_now->atm2)==0)
	  &&(strcasecmp(build_intra->ctors_typ_con[i].atm4,
			build_intra->ctors_typ_now->atm1)==0)
	  &&(strcasecmp(build_intra->ctors_typ_con[i].label,
			build_intra->ctors_typ_now->label)==0)) {itype=i;}
     }       /*endfor*/
     /*--------------------------------------------------------------------*/
     /* D) Add space */
     if(itype>build_intra->ntors_typ_con_max){
       build_intra->ntors_typ_con_max += NMEM_MIN;
       build_intra->ctors_typ_con      = (CTORS *) 
	 crealloc(&(build_intra->ctors_typ_con)[1],
		  build_intra->ntors_typ_con_max*sizeof(CTORS))-1;
     }/*endif*/

     /*-------------------------------------------------------------------*/
     /* E) Add a type */
     if(itype==tors->ntyp_con+1){
       tors->ntyp_con+=1;
       strcpy(build_intra->ctors_typ_con[itype].atm1,
	      build_intra->ctors_typ_now->atm1);
       strcpy(build_intra->ctors_typ_con[itype].atm2,
	      build_intra->ctors_typ_now->atm2);
       strcpy(build_intra->ctors_typ_con[itype].atm3,
	      build_intra->ctors_typ_now->atm3);
       strcpy(build_intra->ctors_typ_con[itype].atm4,
	      build_intra->ctors_typ_now->atm4);
       strcpy(build_intra->ctors_typ_con[itype].label,
	      build_intra->ctors_typ_now->label);
     }/*endif*/

     /*--------------------------------------------------------------------*/
     /* F) Assign */
     tors->ncon += 1;
     tors->j1_con[tors->ncon] = iatm_ind1 + clatoms_info->natm_tot;
     tors->j2_con[tors->ncon] = iatm_ind2 + clatoms_info->natm_tot;
     tors->j3_con[tors->ncon] = iatm_ind3 + clatoms_info->natm_tot;
     tors->j4_con[tors->ncon] = iatm_ind4 + clatoms_info->natm_tot;
     tors->jtyp_con[tors->ncon] = itype;

    /*---------------------------------------------------------------------*/
    /* G) Error */
    if((atommaps->ighost_flag[tors->j1_con[tors->ncon]]!=0)||
       (atommaps->ighost_flag[tors->j2_con[tors->ncon]]!=0)||
       (atommaps->ighost_flag[tors->j3_con[tors->ncon]]!=0)||
       (atommaps->ighost_flag[tors->j4_con[tors->ncon]]!=0)){
        printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
        printf("Ghost atoms are not permitted in constrained torsion\n");
        printf("in molecule index %d in residue index %d in file %s\n",
                jmol_typ,iresidue,file_name);
        printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
        fflush(stdout);
        exit(1);
    }/*endif*/      
   }  /*endif*/
  /*=====================================================================*/
  /* V) Spread nul torss */
  if(ifound==3&&igo==1){
    /*---------------------------------------------------------------------*/
    /* A) Add more space */
    if(null_inter_parse->ntors_nul+1> build_intra->ntors_nul_max){
      build_intra->ntors_nul_max += NMEM_MIN;
      null_inter_parse->jtors1_nul     = 
	(int *) crealloc(&(null_inter_parse->jtors1_nul[1]),
			 build_intra->ntors_nul_max*sizeof(int))-1;
      null_inter_parse->jtors2_nul     = 
         (int *) crealloc(&(null_inter_parse->jtors2_nul[1]),
			  build_intra->ntors_nul_max*sizeof(int))-1;
      null_inter_parse->jtors3_nul     = 
	(int *) crealloc(&(null_inter_parse->jtors3_nul[1]),
			 build_intra->ntors_nul_max*sizeof(int))-1;
      null_inter_parse->jtors4_nul     = 
	(int *) crealloc(&(null_inter_parse->jtors4_nul[1]),
			 build_intra->ntors_nul_max*sizeof(int))-1;
    }      /*endif*/

    /*---------------------------------------------------------------------*/

    null_inter_parse->ntors_nul += 1;

    null_inter_parse->jtors1_nul[null_inter_parse->ntors_nul] = 
                                      iatm_ind1 + clatoms_info->natm_tot ;
    null_inter_parse->jtors2_nul[null_inter_parse->ntors_nul] = 
                                        iatm_ind2 + clatoms_info->natm_tot;
    null_inter_parse->jtors3_nul[null_inter_parse->ntors_nul] = 
	                                iatm_ind3 + clatoms_info->natm_tot;
    null_inter_parse->jtors4_nul[null_inter_parse->ntors_nul] = 
	                                iatm_ind4 + clatoms_info->natm_tot;
  }/*endif*/
  /*----------------------------------------------------------------------*/
}/*end routine*/
/*=======================================================================*/
