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

void set_atm_morph(DICT_WORD *atm_dict,int num_atm_dict,
                   char *fun_key,char *filename,
                   int jmol_typ,ATOMMAPS *atommaps,
                   BUILD_INTRA *build_intra,CLATOMS_INFO *clatoms_info,
                   GHOST_ATOMS *ghost_atoms)

/*=======================================================================*/
{/*begin routine*/
/*=======================================================================*/

  NAME site;
  int i,num,num2,num3,imask,index,isum,num1,num4; 
  double dnum;
  int itype,iatm_ind,iatm_ind_sm,iii;

/*=======================================================================*/
/* I) Get atm index                   */

  sscanf(atm_dict[2].keyarg,"%d",&iatm_ind);
  index = 2;
  if((iatm_ind>build_intra->natm_1res_pure_now) || (iatm_ind <1)){
      keyarg_barf(atm_dict,filename,fun_key,index);
  }/*endif*/
  imask = build_intra->mask_atm[iatm_ind];

/*=======================================================================*/
/*=======================================================================*/

  if(imask>0){

/*=======================================================================*/
/* II) Get rejiggered atm index                                         */
/*-----------------------------------------------------------------------*/

    iatm_ind = build_intra->index_atm[iatm_ind];
    iatm_ind_sm = iatm_ind;
    iatm_ind = iatm_ind + clatoms_info->natm_tot;

/*=======================================================================*/
/* III) Fill the dictionary with words */
/*-----------------------------------------------------------------------*/
  /*  3) \mass{} */
    if(atm_dict[3].iuset==1){
      sscanf(atm_dict[3].keyarg,"%lg",&dnum);
      index = 3;
      if(dnum<=0.0){
        keyarg_barf(atm_dict,filename,fun_key,index);}
        clatoms_info->mass[iatm_ind] = dnum; 
    }/*endif*/
  /*---------------------------------------------------------------------*/
  /*  4) \charge{} */
    if(atm_dict[4].iuset==1){
      sscanf(atm_dict[4].keyarg,"%lg",&dnum);
      clatoms_info->q[iatm_ind] = dnum; 
    }/*endif*/
  /*---------------------------------------------------------------------*/
  /*  5) \alpha_pol{} */
    if(atm_dict[5].iuset==1){
      sscanf(atm_dict[5].keyarg,"%lg",&dnum);
      if(dnum<0.0){
        keyarg_barf(atm_dict,filename,fun_key,index);}
      clatoms_info->alp_pol[iatm_ind] = dnum; 
    }/*endif*/
  /*----------------------------------------------------------------------*/
  /*  6) \b_neut{} */
    if(atm_dict[6].iuset==1){
      sscanf(atm_dict[6].keyarg,"%lg",&dnum);
      clatoms_info->b_neut[iatm_ind] = dnum; 
    }/*endif*/
  /*---------------------------------------------------------------------*/
  /*  7) \valence{} */
    if(atm_dict[7].iuset==1){
      sscanf(atm_dict[7].keyarg,"%d",&num);
      build_intra->bond_site[iatm_ind_sm].valence = num;
      index = 7;
      if(num<0||num>MAX_VALENCE){
        keyarg_barf(atm_dict,filename,fun_key,index);
      }/*endif*/
    }/*endif*/
  /*--------------------------------------------------------------------*/
/*  8) \improper_def{} */
    if(atm_dict[8].iuset==1){
     strcpy(build_intra->strip1,atm_dict[8].keyarg);
     parse_improp(build_intra->strip1,build_intra->strip2,&num1,&num2,&num3,
                                                                      &num4);
     build_intra->bond_site[iatm_ind_sm].improper_ind[1]=num1+1;
     build_intra->bond_site[iatm_ind_sm].improper_ind[2]=num2+1;
     build_intra->bond_site[iatm_ind_sm].improper_ind[3]=num3+1;
     build_intra->bond_site[iatm_ind_sm].improper_ind[4]=num4+1;
     index = 8;
     if((num1+num2+num3+num4)!=0&&(num1+num2+num3+num4)!=6){
       keyarg_barf(atm_dict,filename,fun_key,index);}
     if((num1+num2+num3+num4)==6){
      if(num1<=0||num1>3||num2<=0||num2>3||num3<=0||num3>3||num4<=0||num4>3||
         num1==num2||num1==num3||num2==num3||num1==num4||num2==num4||
         num3==num4){
       keyarg_barf(atm_dict,filename,fun_key,index);}}
    }/*endif*/
  /*---------------------------------------------------------------------*/
  /*  9) \bond_site_2{} */
    if(atm_dict[9].iuset==1){
      strcpy(build_intra->strip1,atm_dict[9].keyarg);
      parse_bond_site(build_intra->strip1,build_intra->strip2,
                      site,&num2,&num3);
      strcpy(build_intra->bond_site[iatm_ind_sm].bond_site_name[1],site);
      build_intra->bond_site[iatm_ind_sm].branch_1[1]     = num2;
      build_intra->bond_site[iatm_ind_sm].branch_2[1] = num3;
      index = 8;
     if(strlen(site)==0||num2<-1||num3<-1||num2>MAX_VALENCE||num3>MAX_VALENCE){
        keyarg_barf(atm_dict,filename,fun_key,index);}
    }/*endif*/
  /*---------------------------------------------------------------------*/
  /*  10) \bond_site_3{} */
    if(atm_dict[10].iuset==1){
      strcpy(build_intra->strip1,atm_dict[10].keyarg);
      parse_bond_site(build_intra->strip1,build_intra->strip2,
                      site,&num2,&num3);
      strcpy(build_intra->bond_site[iatm_ind_sm].bond_site_name[2],site);
      build_intra->bond_site[iatm_ind_sm].branch_1[2]     = num2;
      build_intra->bond_site[iatm_ind_sm].branch_2[2] = num3;
      index = 10;
     if(strlen(site)==0||num2<-1||num3<-1||num2>MAX_VALENCE||num3>MAX_VALENCE){
        keyarg_barf(atm_dict,filename,fun_key,index);}
    }/*endif*/
  /*---------------------------------------------------------------------*/
  /*  11) \bond_site_4{} */
    if(atm_dict[11].iuset==1){
      strcpy(build_intra->strip1,atm_dict[11].keyarg);
      parse_bond_site(build_intra->strip1,build_intra->strip2,
                      site,&num2,&num3);
      strcpy(build_intra->bond_site[iatm_ind_sm].bond_site_name[3],site);
      build_intra->bond_site[iatm_ind_sm].branch_1[3]     = num2;
      build_intra->bond_site[iatm_ind_sm].branch_2[3] = num3;
      index = 11;
     if(strlen(site)==0||num2<-1||num3<-1||num2>MAX_VALENCE||num3>MAX_VALENCE){
        keyarg_barf(atm_dict,filename,fun_key,index);}
    }/*endif*/
  /*---------------------------------------------------------------------*/
  /*  12) \bond_site_1{} */
    if(atm_dict[12].iuset==1){
      strcpy(build_intra->strip1,atm_dict[12].keyarg);
      parse_bond_site(build_intra->strip1,build_intra->strip2,
                      site,&num2,&num3);
      strcpy(build_intra->bond_site[iatm_ind_sm].bond_site_name[4],site);
      build_intra->bond_site[iatm_ind_sm].branch_1[4]     =  num2;
      build_intra->bond_site[iatm_ind_sm].branch_2[4] = num3;
      index = 12;
     if(strlen(site)<=0||num2<-1||num3<-1||num2>MAX_VALENCE||num3>MAX_VALENCE){
        keyarg_barf(atm_dict,filename,fun_key,index);}
    }/*endif*/
  /*---------------------------------------------------------------------*/
  /*  13) \cp_vlnc_up{} */
    if(atm_dict[13].iuset==1){
      sscanf(atm_dict[13].keyarg,"%d",&num);
      clatoms_info->cp_vlnc_up[iatm_ind] = num; 
    }/*endif*/

  /*---------------------------------------------------------------------*/
  /*  14) \cp_vlnc_dn{} */
    if(atm_dict[14].iuset==1){
      sscanf(atm_dict[14].keyarg,"%d",&num);
      clatoms_info->cp_vlnc_dn[iatm_ind] = num; 
    }/*endif*/

  /*---------------------------------------------------------------------*/
  /*  15) \cp_atom{yes,no} */

    if(atm_dict[15].iuset==1){
     if(strcasecmp(atm_dict[15].keyarg,"yes")==0){
       clatoms_info->cp_atm_flag[iatm_ind] = 1;
     }
     if(strcasecmp(atm_dict[15].keyarg,"no")==0){
       clatoms_info->cp_atm_flag[iatm_ind] = 0;
     }
    }/*endif*/

  /*---------------------------------------------------------------------*/
  /*  16) \def ghost1{} */
    isum = 0;
    for(i=1;i<=NCOEF_GHOST_MAX;i++){isum+=atm_dict[(i+11)].iuset;}
    if(isum>0){
     set_ghost(ghost_atoms,clatoms_info,atommaps,build_intra,atm_dict,
               num_atm_dict,
               filename,fun_key,iatm_ind);
    }/*endif*/

  /*=====================================================================*/
  /* III) Atm type */
    itype = atommaps->iatm_atm_typ[iatm_ind + (clatoms_info->natm_tot)];
    if(atm_dict[1].iuset==1){
      itype = atommaps->natm_typ + 1;
      for(i=1;i<=atommaps->natm_typ;i++){
        if(strcasecmp(atm_dict[1].keyarg,atommaps->atm_typ[i])==0)itype=i;
      }/*endfor*/
      if(itype>build_intra->natm_typ_max){
        build_intra->natm_typ_max+=NMEM_MIN;
        atommaps->atm_typ = (NAME *) 
          crealloc(&(atommaps->atm_typ[1]),
                   (build_intra->natm_typ_max)*sizeof(NAME));
      }/*endif*/
      if(itype==atommaps->natm_typ+1){
        atommaps->natm_typ+=1;
        strcpy(atommaps->atm_typ[itype],atm_dict[1].keyarg);
      }/*endif*/
    }/*endif*/
  /*======================================================================*/
  /*======================================================================*/

    }/*endif*/

  /*======================================================================*/
}/*end routine*/
/*==========================================================================*/




