/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
/*  set_grp_bond_params:set up intramolecular interaction                   */
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

void set_grp_bond_params(DICT_WORD intra_dict[],int num_intra_dict,
                     char *fun_key,char *file_name,int jmol_typ,
                     CLATOMS_INFO *clatoms_info,ATOMMAPS *atommaps,
                     GRP_BOND_CON *grp_bond_con,
                     GRP_BOND_WATTS *grp_bond_watts,
                     NULL_INTER_PARSE *null_inter_parse,
                     BUILD_INTRA *build_intra, int iresidue, int ires_off)

/*==========================================================================*/
{/*begin routine*/
/*==========================================================================*/

  int num,index,ifound,igo;
  int itype1,itype2,itype3,itype4,iupper;
  int iatm_ind1,iatm_ind2,iatm_ind3,iatm_ind4;
  int imask1,imask2,imask3,imask4;
  int i,itype,igrp_type,iii;

  /*======================================================================*/
  /* II) Fill the dictionary with words */
#ifdef DEBUG_DAWN
   printf("filling dictionary with words \n");
#endif
  /*----------------------------------------------------------------------*/
  /* 12) Type of grp bond                                                  */
  i = 12;
  if(intra_dict[i].iuset==0 && intra_dict[i].key_type==1){
      keyword_miss(intra_dict,file_name,fun_key,i);
  }
  index  = 12;
  igrp_type = 0;
  if(strcasecmp(intra_dict[12].keyarg,"2x1")==0) {igrp_type = 1;iupper=2;}
  if(strcasecmp(intra_dict[12].keyarg,"2x3")==0) {igrp_type = 2;iupper=3;}
  if(strcasecmp(intra_dict[12].keyarg,"3x3")==0) {igrp_type = 3;iupper=3;}
  if(strcasecmp(intra_dict[12].keyarg,"watts_3x3")==0) 
                                               {igrp_type = 6;iupper=3;}
  if(strcasecmp(intra_dict[12].keyarg,"4x6")==0) {igrp_type = 4;iupper=4;}
  if(strcasecmp(intra_dict[12].keyarg,"4x3")==0) {igrp_type = 5;iupper=4;}
  if(igrp_type==0){
    keyarg_barf(intra_dict,file_name,fun_key,index);
  }/*endif*/
  /*---------------------------------------------------------------------*/
  /*  5) \modifier{} */
  index  = 5;
  ifound = 0;
  if(strcasecmp(intra_dict[5].keyarg,"on")==0){ifound = 3;}
  if(strcasecmp(intra_dict[5].keyarg,"con")==0){ifound = 2;}
  if(strcasecmp(intra_dict[5].keyarg,"off")==0){ifound = 1;}
  if(ifound==0){
    keyarg_barf(intra_dict,file_name,fun_key,index);
  }/*endif*/
  if(ifound != 3 && igrp_type==6){ 
        printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
        printf("Watts bonds must be on.                        \n");
        printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
        fflush(stdout);
        exit(1);
  }
  if(ifound == 3 && igrp_type!=6){ 
        printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
        printf("Only Watts bonds can be on.                    \n");
        printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
        fflush(stdout);
        exit(1);
  }
  if(ifound == 2 && clatoms_info->pi_beads>1){ 
        printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
        printf("Constraints not implemented under  \n");
        printf(" path integral dynamics.            \n");
        printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
        fflush(stdout);
        exit(1);
  }/*endif*/
  if(ifound == 1 && igrp_type!=6){ 
        printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
        printf("Only Watts 3x3 implemented.  \n");
        printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
        fflush(stdout);
        exit(1);
  }/*endif*/
  /*---------------------------------------------------------------------*/
  /* 1) Check for missing atom specifiers                                 */
  for(i=1;i<=iupper;i++){
    if(intra_dict[i].iuset==0 && intra_dict[i].key_type==1){
      keyword_miss(intra_dict,file_name,fun_key,i);
    }
  }
  /*---------------------------------------------------------------------*/
  /*  1) \atom1{}    */
  index = 1;
  sscanf(intra_dict[1].keyarg,"%d",&num);
  iatm_ind1 = num;
  if(iatm_ind1>build_intra->natmind_1res_now||iatm_ind1<0){
            keyarg_barf(intra_dict,file_name,fun_key,index);}
  imask1 = build_intra->mask_atm[iatm_ind1];  
  if(imask1>0){iatm_ind1 = build_intra->index_atm[iatm_ind1];}

  /*----------------------------------------------------------------------*/
  /*  2) \atom2{}    */
  index = 2;
  sscanf(intra_dict[2].keyarg,"%d",&num);
  iatm_ind2 = num;
  if(iatm_ind2>build_intra->natmind_1res_now||iatm_ind2<0){
            keyarg_barf(intra_dict,file_name,fun_key,index);}
  imask2 = build_intra->mask_atm[iatm_ind2];  
  if(imask2>0){iatm_ind2 = build_intra->index_atm[iatm_ind2];}
  /*------------------------------------------------------------------------*/
  /*  3) \atom3{}    */
  if(iupper > 2){
   index = 3;
   sscanf(intra_dict[3].keyarg,"%d",&num);
   iatm_ind3 = num;
   if(iatm_ind3>build_intra->natmind_1res_now||iatm_ind3<0){
             keyarg_barf(intra_dict,file_name,fun_key,index);}
   imask3 = build_intra->mask_atm[iatm_ind3];  
   if(imask3>0){iatm_ind3 = build_intra->index_atm[iatm_ind3];}

  /*----------------------------------------------------------------------*/
  /*  4) \atom4{}    */
   if(iupper > 3){
    index = 4;
    if(intra_dict[index].iuset==0 && intra_dict[index].key_type==1){
      keyword_miss(intra_dict,file_name,fun_key,index);
    }/*endif*/
    sscanf(intra_dict[4].keyarg,"%d",&num);
    iatm_ind4 = num;
    if(iatm_ind4>build_intra->natmind_1res_now||iatm_ind4<0){
            keyarg_barf(intra_dict,file_name,fun_key,index);}
    imask4 = build_intra->mask_atm[iatm_ind4];  
    if(imask4>0){iatm_ind4 = build_intra->index_atm[iatm_ind4];}
   }/*endif*/
  }/* endif */
  /*-----------------------------------------------------------------------*/
  /*  Construct the present grp constraint type*/
  igo = imask1*imask2;
  if(iupper > 2) igo = imask1*imask2*imask3;
  if(iupper > 3){igo = imask1*imask2*imask3*imask4;}
  if(igo==1){
     itype1 = atommaps->iatm_atm_typ[(clatoms_info->natm_tot+iatm_ind1)];
     itype2 = atommaps->iatm_atm_typ[(clatoms_info->natm_tot+iatm_ind2)];
     strcpy(build_intra->cgrp_typ_now->atm1,atommaps->atm_typ[itype1]);
     strcpy(build_intra->cgrp_typ_now->atm2,atommaps->atm_typ[itype2]);
     if(iupper > 2){
      itype3 = atommaps->iatm_atm_typ[(clatoms_info->natm_tot+iatm_ind3)];
      strcpy(build_intra->cgrp_typ_now->atm3,atommaps->atm_typ[itype3]);
      if(iupper>3){
        itype4 = atommaps->iatm_atm_typ[(clatoms_info->natm_tot+iatm_ind4)];
        strcpy(build_intra->cgrp_typ_now->atm4,atommaps->atm_typ[itype4]);
      }else{
        strcpy(build_intra->cgrp_typ_now->atm4,"");
      }/* endif */
     }else{
       strcpy(build_intra->cgrp_typ_now->atm3,"");
     }/* endif */
      strcpy(build_intra->cgrp_typ_now->label[1],intra_dict[6].keyarg);  
      strcpy(build_intra->cgrp_typ_now->label[2],intra_dict[7].keyarg);  
      strcpy(build_intra->cgrp_typ_now->label[3],intra_dict[8].keyarg);  
      strcpy(build_intra->cgrp_typ_now->label[4],intra_dict[9].keyarg);  
      strcpy(build_intra->cgrp_typ_now->label[5],intra_dict[10].keyarg);  
      strcpy(build_intra->cgrp_typ_now->label[6],intra_dict[11].keyarg);  
  }/*endif*/

/*======================================================================*/
/*======================================================================*/
/*                     Spread GPRCONS */


  if(igo==1){

    ires_off = atommaps->jres_jmol_typ_strt[jmol_typ]-1;
    if(igrp_type!=6){
      atommaps->icons_jres_jmol_typ[ires_off+iresidue] = 1;
    }/*endif*/

/*======================================================================*/
/* )) 21 */
    if(igrp_type==1){
#ifdef DEBUG_DAWN
     printf("in igrp 21 \n");
#endif
     atommaps->nfree_jres_jmol_typ[ires_off+iresidue] -= 1;
/*---------------------------------------------------------------------*/
/* A) Add more space */

     if(grp_bond_con->num_21+1 > build_intra->ngrp_21_max){
       build_intra->ngrp_21_max += NMEM_MIN;
       grp_bond_con->j1_21     = (int *) crealloc(&(grp_bond_con->j1_21)[1],
                          build_intra->ngrp_21_max*sizeof(int))-1;
       grp_bond_con->j2_21     = (int *) crealloc(&(grp_bond_con->j2_21)[1],
                          build_intra->ngrp_21_max*sizeof(int))-1;
       grp_bond_con->jtyp_21   = (int *) crealloc(&(grp_bond_con->jtyp_21)[1],
                          build_intra->ngrp_21_max*sizeof(int))-1;
     }     /*endif*/
/*--------------------------------------------------------------------*/
/* C) Check type */
     itype = grp_bond_con->ntyp_21+1;
     for(i=1;i<= grp_bond_con->ntyp_21;i++){
       if((strcasecmp(build_intra->cgrp_typ_21[i].atm1,
                      build_intra->cgrp_typ_now->atm1)==0)
          &&(strcasecmp(build_intra->cgrp_typ_21[i].atm2,
                        build_intra->cgrp_typ_now->atm2)==0)
          &&(strcasecmp(build_intra->cgrp_typ_21[i].atm3,
                        build_intra->cgrp_typ_now->atm3)==0)
          &&(strcasecmp(build_intra->cgrp_typ_21[i].atm4,
                        build_intra->cgrp_typ_now->atm4)==0)
          &&(strcasecmp(build_intra->cgrp_typ_21[i].label[1],
                  build_intra->cgrp_typ_now->label[1])==0)
          &&(strcasecmp(build_intra->cgrp_typ_21[i].label[2],
                  build_intra->cgrp_typ_now->label[2])==0)
                                                          ) {itype=i;}
       
       if((strcasecmp(build_intra->cgrp_typ_21[i].atm1,
                      build_intra->cgrp_typ_now->atm4)==0)
          &&(strcasecmp(build_intra->cgrp_typ_21[i].atm2,
                        build_intra->cgrp_typ_now->atm3)==0)
          &&(strcasecmp(build_intra->cgrp_typ_21[i].atm3,
                        build_intra->cgrp_typ_now->atm2)==0)
          &&(strcasecmp(build_intra->cgrp_typ_21[i].atm4,
                        build_intra->cgrp_typ_now->atm1)==0)
          &&(strcasecmp(build_intra->cgrp_typ_21[i].label[1],
                  build_intra->cgrp_typ_now->label[2])==0)
          &&(strcasecmp(build_intra->cgrp_typ_21[i].label[2],
                  build_intra->cgrp_typ_now->label[1])==0)
                                                          ) {itype=i;}
     }       /*endfor*/
/*--------------------------------------------------------------------*/
/* D) Add space */
     if(itype>build_intra->ngrp_typ_21_max){
       build_intra->ngrp_typ_21_max += NMEM_MIN;
       build_intra->cgrp_typ_21      = (CGRP_CONS *) 
                         crealloc(&(build_intra->cgrp_typ_21)[1],
              build_intra->ngrp_typ_21_max*sizeof(CGRP_CONS))-1;
     }/*endif*/

/*-------------------------------------------------------------------*/
/* E) Add a type */
     if(itype==grp_bond_con->ntyp_21+1){
       grp_bond_con->ntyp_21+=1;
       strcpy(build_intra->cgrp_typ_21[itype].atm1,
              build_intra->cgrp_typ_now->atm1);
       strcpy(build_intra->cgrp_typ_21[itype].atm2,
              build_intra->cgrp_typ_now->atm2);
       strcpy(build_intra->cgrp_typ_21[itype].atm3,
              build_intra->cgrp_typ_now->atm3);
       strcpy(build_intra->cgrp_typ_21[itype].atm4,
              build_intra->cgrp_typ_now->atm4);
       strcpy(build_intra->cgrp_typ_21[itype].label[1],
              build_intra->cgrp_typ_now->label[1]);
       strcpy(build_intra->cgrp_typ_21[itype].label[2],
              build_intra->cgrp_typ_now->label[2]);
     }/*endif*/

/*--------------------------------------------------------------------*/
/* B) Spread */
     grp_bond_con->num_21 += 1;
     grp_bond_con->j1_21[grp_bond_con->num_21] = iatm_ind1 + 
                                                    clatoms_info->natm_tot;
     grp_bond_con->j2_21[grp_bond_con->num_21] = iatm_ind2 + 
                                                    clatoms_info->natm_tot;
     grp_bond_con->jtyp_21[grp_bond_con->num_21] = itype;
     if((atommaps->ighost_flag[grp_bond_con->j1_21[grp_bond_con->num_21]]!=0)||
        (atommaps->ighost_flag[grp_bond_con->j2_21[grp_bond_con->num_21]]!=0)){
        printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
        printf("Ghost atoms are not permitted in group constraints\n");
        printf("in molecule index %d in residue index %d in file %s\n",
                jmol_typ,iresidue,file_name);
        printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
        fflush(stdout);
        exit(1);
     }/*endif*/      

    }/*endif:2x1*/

/*======================================================================*/
/* I) 23 */
    if(igrp_type==2){
     atommaps->nfree_jres_jmol_typ[ires_off+iresidue] -= 2;
/*---------------------------------------------------------------------*/
/* A) Add more space */
     if(grp_bond_con->num_23+1 > build_intra->ngrp_23_max){
       build_intra->ngrp_23_max += NMEM_MIN;
       grp_bond_con->j1_23     = (int *) crealloc(&(grp_bond_con->j1_23)[1],
                          build_intra->ngrp_23_max*sizeof(int))-1;
       grp_bond_con->j2_23     = (int *) crealloc(&(grp_bond_con->j2_23)[1],
                          build_intra->ngrp_23_max*sizeof(int))-1;
       grp_bond_con->j3_23     = (int *) crealloc(&(grp_bond_con->j3_23)[1],
                          build_intra->ngrp_23_max*sizeof(int))-1;
       grp_bond_con->jtyp_23   = (int *) crealloc(&(grp_bond_con->jtyp_23)[1],
                          build_intra->ngrp_23_max*sizeof(int))-1;
     }     /*endif*/
/*--------------------------------------------------------------------*/
/* C) Check type */
     itype = grp_bond_con->ntyp_23+1;
     for(i=1;i<= grp_bond_con->ntyp_23;i++){
       if((strcasecmp(build_intra->cgrp_typ_23[i].atm1,
                      build_intra->cgrp_typ_now->atm1)==0)
          &&(strcasecmp(build_intra->cgrp_typ_23[i].atm2,
                        build_intra->cgrp_typ_now->atm2)==0)
          &&(strcasecmp(build_intra->cgrp_typ_23[i].atm3,
                        build_intra->cgrp_typ_now->atm3)==0)
          &&(strcasecmp(build_intra->cgrp_typ_23[i].atm4,
                        build_intra->cgrp_typ_now->atm4)==0)
          &&(strcasecmp(build_intra->cgrp_typ_23[i].label[1],
                  build_intra->cgrp_typ_now->label[1])==0)
          &&(strcasecmp(build_intra->cgrp_typ_23[i].label[2],
                  build_intra->cgrp_typ_now->label[2])==0)
                                                          ) {itype=i;}
       
       if((strcasecmp(build_intra->cgrp_typ_23[i].atm1,
                      build_intra->cgrp_typ_now->atm4)==0)
          &&(strcasecmp(build_intra->cgrp_typ_23[i].atm2,
                        build_intra->cgrp_typ_now->atm3)==0)
          &&(strcasecmp(build_intra->cgrp_typ_23[i].atm3,
                        build_intra->cgrp_typ_now->atm2)==0)
          &&(strcasecmp(build_intra->cgrp_typ_23[i].atm4,
                        build_intra->cgrp_typ_now->atm1)==0)
          &&(strcasecmp(build_intra->cgrp_typ_23[i].label[1],
                  build_intra->cgrp_typ_now->label[2])==0)
          &&(strcasecmp(build_intra->cgrp_typ_23[i].label[2],
                  build_intra->cgrp_typ_now->label[1])==0)
                                                          ) {itype=i;}
     }       /*endfor*/
/*--------------------------------------------------------------------*/
/* D) Add space */
     if(itype>build_intra->ngrp_typ_23_max){
       build_intra->ngrp_typ_23_max += NMEM_MIN;
       build_intra->cgrp_typ_23      = (CGRP_CONS *) 
                         crealloc(&(build_intra->cgrp_typ_23)[1],
              build_intra->ngrp_typ_23_max*sizeof(CGRP_CONS))-1;
     }/*endif*/

/*-------------------------------------------------------------------*/
/* E) Add a type */
     if(itype==grp_bond_con->ntyp_23+1){
       grp_bond_con->ntyp_23+=1;
       strcpy(build_intra->cgrp_typ_23[itype].atm1,
              build_intra->cgrp_typ_now->atm1);
       strcpy(build_intra->cgrp_typ_23[itype].atm2,
              build_intra->cgrp_typ_now->atm2);
       strcpy(build_intra->cgrp_typ_23[itype].atm3,
              build_intra->cgrp_typ_now->atm3);
       strcpy(build_intra->cgrp_typ_23[itype].atm4,
              build_intra->cgrp_typ_now->atm4);
       strcpy(build_intra->cgrp_typ_23[itype].label[1],
              build_intra->cgrp_typ_now->label[1]);
       strcpy(build_intra->cgrp_typ_23[itype].label[2],
              build_intra->cgrp_typ_now->label[2]);
     }/*endif*/

/*--------------------------------------------------------------------*/
/* B) Spread */
     grp_bond_con->num_23 += 1;
     grp_bond_con->j1_23[grp_bond_con->num_23] = iatm_ind1 + 
                                                    clatoms_info->natm_tot;
     grp_bond_con->j2_23[grp_bond_con->num_23] = iatm_ind2 + 
                                                    clatoms_info->natm_tot;
     grp_bond_con->j3_23[grp_bond_con->num_23] = iatm_ind3 + 
                                                    clatoms_info->natm_tot;
     grp_bond_con->jtyp_23[grp_bond_con->num_23] = itype;
     if((atommaps->ighost_flag[grp_bond_con->j1_23[grp_bond_con->num_23]]!=0)||
        (atommaps->ighost_flag[grp_bond_con->j2_23[grp_bond_con->num_23]]!=0)||
        (atommaps->ighost_flag[grp_bond_con->j3_23[grp_bond_con->num_23]]!=0)){
        printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
        printf("Ghost atoms are not permitted in group constraints\n");
        printf("in molecule index %d in residue index %d in file %s\n",
                jmol_typ,iresidue,file_name);
        printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
        fflush(stdout);
        exit(1);
     }/*endif*/      

    }/*endif:2x3*/

/*======================================================================*/
/* I) 33 */

    if(igrp_type==3){
     atommaps->nfree_jres_jmol_typ[ires_off+iresidue] -= 3;
/*---------------------------------------------------------------------*/
/* A) Add more space */
     if(grp_bond_con->num_33+1 > build_intra->ngrp_33_max){
       build_intra->ngrp_33_max += NMEM_MIN;
       grp_bond_con->j1_33     = (int *) crealloc(&(grp_bond_con->j1_33)[1],
                          build_intra->ngrp_33_max*sizeof(int))-1;
       grp_bond_con->j2_33     = (int *) crealloc(&(grp_bond_con->j2_33)[1],
                          build_intra->ngrp_33_max*sizeof(int))-1;
       grp_bond_con->j3_33     = (int *) crealloc(&(grp_bond_con->j3_33)[1],
                          build_intra->ngrp_33_max*sizeof(int))-1;
       grp_bond_con->jtyp_33   = (int *) crealloc(&(grp_bond_con->jtyp_33)[1],
                          build_intra->ngrp_33_max*sizeof(int))-1;
     }     /*endif*/
/*--------------------------------------------------------------------*/
/* C) Check type */
     itype = grp_bond_con->ntyp_33+1;
     for(i=1;i<= grp_bond_con->ntyp_33;i++){
       if((strcasecmp(build_intra->cgrp_typ_33[i].atm1,
                      build_intra->cgrp_typ_now->atm1)==0)
          &&(strcasecmp(build_intra->cgrp_typ_33[i].atm2,
                        build_intra->cgrp_typ_now->atm2)==0)
          &&(strcasecmp(build_intra->cgrp_typ_33[i].atm3,
                        build_intra->cgrp_typ_now->atm3)==0)
          &&(strcasecmp(build_intra->cgrp_typ_33[i].atm4,
                        build_intra->cgrp_typ_now->atm4)==0)
          &&(strcasecmp(build_intra->cgrp_typ_33[i].label[1],
                  build_intra->cgrp_typ_now->label[1])==0)
          &&(strcasecmp(build_intra->cgrp_typ_33[i].label[2],
                  build_intra->cgrp_typ_now->label[2])==0)
          &&(strcasecmp(build_intra->cgrp_typ_33[i].label[3],
                  build_intra->cgrp_typ_now->label[3])==0)
                                                          ) {itype=i;}
       
       if((strcasecmp(build_intra->cgrp_typ_33[i].atm1,
                      build_intra->cgrp_typ_now->atm4)==0)
          &&(strcasecmp(build_intra->cgrp_typ_33[i].atm2,
                        build_intra->cgrp_typ_now->atm3)==0)
          &&(strcasecmp(build_intra->cgrp_typ_33[i].atm3,
                        build_intra->cgrp_typ_now->atm2)==0)
          &&(strcasecmp(build_intra->cgrp_typ_33[i].atm4,
                        build_intra->cgrp_typ_now->atm1)==0)
          &&(strcasecmp(build_intra->cgrp_typ_33[i].label[1],
                  build_intra->cgrp_typ_now->label[3])==0)
          &&(strcasecmp(build_intra->cgrp_typ_33[i].label[2],
                  build_intra->cgrp_typ_now->label[2])==0)
          &&(strcasecmp(build_intra->cgrp_typ_33[i].label[3],
                  build_intra->cgrp_typ_now->label[1])==0)
                                                          ) {itype=i;}
     }       /*endfor*/
/*--------------------------------------------------------------------*/
/* D) Add space */
     if(itype>build_intra->ngrp_typ_33_max){
       build_intra->ngrp_typ_33_max += NMEM_MIN;
       build_intra->cgrp_typ_33      = (CGRP_CONS *) 
         crealloc(&(build_intra->cgrp_typ_33)[1],
              build_intra->ngrp_typ_33_max*sizeof(CGRP_CONS))-1;
     }/*endif*/

/*-------------------------------------------------------------------*/
/* E) Add a type */
     if(itype==grp_bond_con->ntyp_33+1){
       grp_bond_con->ntyp_33+=1;
       strcpy(build_intra->cgrp_typ_33[itype].atm1,
              build_intra->cgrp_typ_now->atm1);
       strcpy(build_intra->cgrp_typ_33[itype].atm2,
              build_intra->cgrp_typ_now->atm2);
       strcpy(build_intra->cgrp_typ_33[itype].atm3,
              build_intra->cgrp_typ_now->atm3);
       strcpy(build_intra->cgrp_typ_33[itype].atm4,
              build_intra->cgrp_typ_now->atm4);
       strcpy(build_intra->cgrp_typ_33[itype].label[1],
              build_intra->cgrp_typ_now->label[1]);
       strcpy(build_intra->cgrp_typ_33[itype].label[2],
              build_intra->cgrp_typ_now->label[2]);
       strcpy(build_intra->cgrp_typ_33[itype].label[3],
              build_intra->cgrp_typ_now->label[3]);
     }/*endif*/

/*--------------------------------------------------------------------*/
/* B) Spread */
     grp_bond_con->num_33 += 1;
     grp_bond_con->j1_33[grp_bond_con->num_33] = iatm_ind1 + 
                                                    clatoms_info->natm_tot;
     grp_bond_con->j2_33[grp_bond_con->num_33] = iatm_ind2 + 
                                                    clatoms_info->natm_tot;
     grp_bond_con->j3_33[grp_bond_con->num_33] = iatm_ind3 + 
                                                    clatoms_info->natm_tot;
     grp_bond_con->jtyp_33[grp_bond_con->num_33] = itype;
     if((atommaps->ighost_flag[grp_bond_con->j1_33[grp_bond_con->num_33]]!=0)||
        (atommaps->ighost_flag[grp_bond_con->j2_33[grp_bond_con->num_33]]!=0)||
        (atommaps->ighost_flag[grp_bond_con->j3_33[grp_bond_con->num_33]]!=0)){
        printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
        printf("Ghost atoms are not permitted in group constraints\n");
        printf("in molecule index %d in residue index %d in file %s\n",
                jmol_typ,iresidue,file_name);
        printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
        fflush(stdout);
        exit(1);
     }/*endif*/      
    }/*endif:3x3*/

/*======================================================================*/
/* I) Watts 33 */

    if(igrp_type==6){
/*---------------------------------------------------------------------*/
/* A) Add more space */
     if(grp_bond_watts->num_33+1 > build_intra->ngrp_watt_33_max){
       build_intra->ngrp_watt_33_max += NMEM_MIN;
       grp_bond_watts->j1_33     = 
                          (int *) crealloc(&(grp_bond_watts->j1_33)[1],
                          build_intra->ngrp_watt_33_max*sizeof(int))-1;
       grp_bond_watts->j2_33     = 
                          (int *) crealloc(&(grp_bond_watts->j2_33)[1],
                          build_intra->ngrp_watt_33_max*sizeof(int))-1;
       grp_bond_watts->j3_33     = 
                          (int *) crealloc(&(grp_bond_watts->j3_33)[1],
                          build_intra->ngrp_watt_33_max*sizeof(int))-1;
       grp_bond_watts->jtyp_33   = 
                          (int *) crealloc(&(grp_bond_watts->jtyp_33)[1],
                          build_intra->ngrp_watt_33_max*sizeof(int))-1;
     }     /*endif*/
/*--------------------------------------------------------------------*/
/* C) Check type */
     itype = grp_bond_watts->ntyp_33+1;
     for(i=1;i<= grp_bond_watts->ntyp_33;i++){
       if((strcasecmp(build_intra->cgrp_typ_watt_33[i].atm1,
                      build_intra->cgrp_typ_now->atm1)==0)
          &&(strcasecmp(build_intra->cgrp_typ_watt_33[i].atm2,
                        build_intra->cgrp_typ_now->atm2)==0)
          &&(strcasecmp(build_intra->cgrp_typ_watt_33[i].atm3,
                        build_intra->cgrp_typ_now->atm3)==0)
          &&(strcasecmp(build_intra->cgrp_typ_watt_33[i].atm4,
                        build_intra->cgrp_typ_now->atm4)==0)
          &&(strcasecmp(build_intra->cgrp_typ_watt_33[i].label[1],
                  build_intra->cgrp_typ_now->label[1])==0)
          &&(strcasecmp(build_intra->cgrp_typ_watt_33[i].label[2],
                  build_intra->cgrp_typ_now->label[2])==0)
          &&(strcasecmp(build_intra->cgrp_typ_watt_33[i].label[3],
                  build_intra->cgrp_typ_now->label[3])==0)
                                                          ) {itype=i;}
       
       if((strcasecmp(build_intra->cgrp_typ_watt_33[i].atm1,
                      build_intra->cgrp_typ_now->atm4)==0)
          &&(strcasecmp(build_intra->cgrp_typ_watt_33[i].atm2,
                        build_intra->cgrp_typ_now->atm3)==0)
          &&(strcasecmp(build_intra->cgrp_typ_watt_33[i].atm3,
                        build_intra->cgrp_typ_now->atm2)==0)
          &&(strcasecmp(build_intra->cgrp_typ_watt_33[i].atm4,
                        build_intra->cgrp_typ_now->atm1)==0)
          &&(strcasecmp(build_intra->cgrp_typ_watt_33[i].label[1],
                  build_intra->cgrp_typ_now->label[3])==0)
          &&(strcasecmp(build_intra->cgrp_typ_watt_33[i].label[2],
                  build_intra->cgrp_typ_now->label[2])==0)
          &&(strcasecmp(build_intra->cgrp_typ_watt_33[i].label[3],
                  build_intra->cgrp_typ_now->label[1])==0)
                                                          ) {itype=i;}
     }       /*endfor*/
/*--------------------------------------------------------------------*/
/* D) Add space */
     if(itype>build_intra->ngrp_typ_watt_33_max){
       build_intra->ngrp_typ_watt_33_max += NMEM_MIN;
       build_intra->cgrp_typ_watt_33      = (CGRP_CONS *) 
         crealloc(&(build_intra->cgrp_typ_watt_33)[1],
              build_intra->ngrp_typ_watt_33_max*sizeof(CGRP_CONS))-1;
     }/*endif*/

/*-------------------------------------------------------------------*/
/* E) Add a type */
     if(itype==grp_bond_watts->ntyp_33+1){
       grp_bond_watts->ntyp_33+=1;
       strcpy(build_intra->cgrp_typ_watt_33[itype].atm1,
              build_intra->cgrp_typ_now->atm1);
       strcpy(build_intra->cgrp_typ_watt_33[itype].atm2,
              build_intra->cgrp_typ_now->atm2);
       strcpy(build_intra->cgrp_typ_watt_33[itype].atm3,
              build_intra->cgrp_typ_now->atm3);
       strcpy(build_intra->cgrp_typ_watt_33[itype].atm4,
              build_intra->cgrp_typ_now->atm4);
       strcpy(build_intra->cgrp_typ_watt_33[itype].label[1],
              build_intra->cgrp_typ_now->label[1]);
       strcpy(build_intra->cgrp_typ_watt_33[itype].label[2],
              build_intra->cgrp_typ_now->label[2]);
       strcpy(build_intra->cgrp_typ_watt_33[itype].label[3],
              build_intra->cgrp_typ_now->label[3]);
     }/*endif*/

/*--------------------------------------------------------------------*/
/* B) Spread */
     grp_bond_watts->num_33 += 1;
     grp_bond_watts->j1_33[grp_bond_watts->num_33] = iatm_ind1 + 
                                                    clatoms_info->natm_tot;
     grp_bond_watts->j2_33[grp_bond_watts->num_33] = iatm_ind2 + 
                                                    clatoms_info->natm_tot;
     grp_bond_watts->j3_33[grp_bond_watts->num_33] = iatm_ind3 + 
                                                    clatoms_info->natm_tot;
     grp_bond_watts->jtyp_33[grp_bond_watts->num_33] = itype;
    }/*endif:watts_3x3*/

/*======================================================================*/
/* I) 43 */

    if(igrp_type==5){
     atommaps->nfree_jres_jmol_typ[ires_off+iresidue] -= 3;
/*---------------------------------------------------------------------*/
/* A) Add more space */
     if(grp_bond_con->num_43+1 > build_intra->ngrp_43_max){
       build_intra->ngrp_43_max += NMEM_MIN;
       grp_bond_con->j1_43     = (int *) crealloc(&(grp_bond_con->j1_43)[1],
                          build_intra->ngrp_43_max*sizeof(int))-1;
       grp_bond_con->j2_43     = (int *) crealloc(&(grp_bond_con->j2_43)[1],
                          build_intra->ngrp_43_max*sizeof(int))-1;
       grp_bond_con->j3_43     = (int *) crealloc(&(grp_bond_con->j3_43)[1],
                          build_intra->ngrp_43_max*sizeof(int))-1;
       grp_bond_con->j4_43     =  (int *) crealloc(&(grp_bond_con->j4_43)[1],
                          build_intra->ngrp_43_max*sizeof(int))-1;
       grp_bond_con->jtyp_43   = (int *) crealloc(&(grp_bond_con->jtyp_43)[1],
                          build_intra->ngrp_43_max*sizeof(int))-1;
     }     /*endif*/
/*--------------------------------------------------------------------*/
/* C) Check type */
     itype = grp_bond_con->ntyp_43+1;
     for(i=1;i<= grp_bond_con->ntyp_43;i++){
       if((strcasecmp(build_intra->cgrp_typ_43[i].atm1,
                      build_intra->cgrp_typ_now->atm1)==0)
          &&(strcasecmp(build_intra->cgrp_typ_43[i].atm2,
                        build_intra->cgrp_typ_now->atm2)==0)
          &&(strcasecmp(build_intra->cgrp_typ_43[i].atm3,
                        build_intra->cgrp_typ_now->atm3)==0)
          &&(strcasecmp(build_intra->cgrp_typ_43[i].atm4,
                        build_intra->cgrp_typ_now->atm4)==0)
          &&(strcasecmp(build_intra->cgrp_typ_43[i].label[1],
                  build_intra->cgrp_typ_now->label[1])==0) 
          &&(strcasecmp(build_intra->cgrp_typ_43[i].label[2],
                  build_intra->cgrp_typ_now->label[2])==0) 
          &&(strcasecmp(build_intra->cgrp_typ_43[i].label[3],
                  build_intra->cgrp_typ_now->label[3])==0) 
          &&(strcasecmp(build_intra->cgrp_typ_43[i].label[4],
                  build_intra->cgrp_typ_now->label[4])==0) 
          &&(strcasecmp(build_intra->cgrp_typ_43[i].label[5],
                  build_intra->cgrp_typ_now->label[5])==0) 
          &&(strcasecmp(build_intra->cgrp_typ_43[i].label[6],
                  build_intra->cgrp_typ_now->label[6])==0) 
                                                         ){itype=i;}
       
       if((strcasecmp(build_intra->cgrp_typ_43[i].atm1,
                      build_intra->cgrp_typ_now->atm4)==0)
          &&(strcasecmp(build_intra->cgrp_typ_43[i].atm2,
                        build_intra->cgrp_typ_now->atm3)==0)
          &&(strcasecmp(build_intra->cgrp_typ_43[i].atm3,
                        build_intra->cgrp_typ_now->atm2)==0)
          &&(strcasecmp(build_intra->cgrp_typ_43[i].atm4,
                        build_intra->cgrp_typ_now->atm1)==0)
          &&(strcasecmp(build_intra->cgrp_typ_43[i].label[1],
                  build_intra->cgrp_typ_now->label[6])==0) 
          &&(strcasecmp(build_intra->cgrp_typ_43[i].label[2],
                  build_intra->cgrp_typ_now->label[5])==0) 
          &&(strcasecmp(build_intra->cgrp_typ_43[i].label[3],
                  build_intra->cgrp_typ_now->label[3])==0) 
          &&(strcasecmp(build_intra->cgrp_typ_43[i].label[4],
                  build_intra->cgrp_typ_now->label[4])==0) 
          &&(strcasecmp(build_intra->cgrp_typ_43[i].label[5],
                  build_intra->cgrp_typ_now->label[2])==0) 
          &&(strcasecmp(build_intra->cgrp_typ_43[i].label[6],
                  build_intra->cgrp_typ_now->label[1])==0) 
                                                         ){itype=i;}
     }       /*endfor*/
/*--------------------------------------------------------------------*/
/* D) Add space */
     if(itype>build_intra->ngrp_typ_43_max){
       build_intra->ngrp_typ_43_max += NMEM_MIN;
       build_intra->cgrp_typ_43      = (CGRP_CONS *) 
         crealloc(&(build_intra->cgrp_typ_43)[1],
              build_intra->ngrp_typ_43_max*sizeof(CGRP_CONS))-1;
     }/*endif*/

/*-------------------------------------------------------------------*/
/* E) Add a type */
     if(itype==grp_bond_con->ntyp_43+1){
       grp_bond_con->ntyp_43+=1;
       strcpy(build_intra->cgrp_typ_43[itype].atm1,
              build_intra->cgrp_typ_now->atm1);
       strcpy(build_intra->cgrp_typ_43[itype].atm2,
              build_intra->cgrp_typ_now->atm2);
       strcpy(build_intra->cgrp_typ_43[itype].atm3,
              build_intra->cgrp_typ_now->atm3);
       strcpy(build_intra->cgrp_typ_43[itype].atm4,
              build_intra->cgrp_typ_now->atm4);
       strcpy(build_intra->cgrp_typ_43[itype].label[1],
              build_intra->cgrp_typ_now->label[1]);
       strcpy(build_intra->cgrp_typ_43[itype].label[2],
              build_intra->cgrp_typ_now->label[2]);
       strcpy(build_intra->cgrp_typ_43[itype].label[3],
              build_intra->cgrp_typ_now->label[3]);
       strcpy(build_intra->cgrp_typ_43[itype].label[4],
              build_intra->cgrp_typ_now->label[4]);
       strcpy(build_intra->cgrp_typ_43[itype].label[5],
              build_intra->cgrp_typ_now->label[5]);
       strcpy(build_intra->cgrp_typ_43[itype].label[6],
              build_intra->cgrp_typ_now->label[6]);
     }/*endif*/

/*--------------------------------------------------------------------*/
/* B) Spread */
     grp_bond_con->num_43 += 1;
     grp_bond_con->j1_43[grp_bond_con->num_43] = iatm_ind1 + 
                                                     clatoms_info->natm_tot;
     grp_bond_con->j2_43[grp_bond_con->num_43] = iatm_ind2 + 
                                                     clatoms_info->natm_tot;
     grp_bond_con->j3_43[grp_bond_con->num_43] = iatm_ind3 + 
                                                     clatoms_info->natm_tot;
     grp_bond_con->j4_43[grp_bond_con->num_43] = iatm_ind4 + 
                                                     clatoms_info->natm_tot;
     grp_bond_con->jtyp_43[grp_bond_con->num_43] = itype;
     if((atommaps->ighost_flag[grp_bond_con->j1_43[grp_bond_con->num_43]]!=0)||
        (atommaps->ighost_flag[grp_bond_con->j2_43[grp_bond_con->num_43]]!=0)||
        (atommaps->ighost_flag[grp_bond_con->j3_43[grp_bond_con->num_43]]!=0)||
        (atommaps->ighost_flag[grp_bond_con->j4_43[grp_bond_con->num_43]]!=0)){
        printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
        printf("Ghost atoms are not permitted in group constraints\n");
        printf("in molecule index %d in residue index %d in file %s\n",
                jmol_typ,iresidue,file_name);
        printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
        fflush(stdout);
        exit(1);
     }/*endif*/      

    }/*endif:4x6*/
/*======================================================================*/
/* I) 46 */

    if(igrp_type==4){
     atommaps->nfree_jres_jmol_typ[ires_off+iresidue] -= 6;
/*---------------------------------------------------------------------*/
/* A) Add more space */
     if(grp_bond_con->num_46+1 > build_intra->ngrp_46_max){
       build_intra->ngrp_46_max += NMEM_MIN;
       grp_bond_con->j1_46     = (int *) crealloc(&(grp_bond_con->j1_46)[1],
                          build_intra->ngrp_46_max*sizeof(int))-1;
       grp_bond_con->j2_46     = (int *) crealloc(&(grp_bond_con->j2_46)[1],
                          build_intra->ngrp_46_max*sizeof(int))-1;
       grp_bond_con->j3_46     = (int *) crealloc(&(grp_bond_con->j3_46)[1],
                          build_intra->ngrp_46_max*sizeof(int))-1;
       grp_bond_con->j4_46     =  (int *) crealloc(&(grp_bond_con->j4_46)[1],
                          build_intra->ngrp_46_max*sizeof(int))-1;
       grp_bond_con->jtyp_46   = (int *) crealloc(&(grp_bond_con->jtyp_46)[1],
                          build_intra->ngrp_46_max*sizeof(int))-1;
     }     /*endif*/
/*--------------------------------------------------------------------*/
/* C) Check type */
     itype = grp_bond_con->ntyp_46+1;
     for(i=1;i<= grp_bond_con->ntyp_46;i++){
       if((strcasecmp(build_intra->cgrp_typ_46[i].atm1,
                      build_intra->cgrp_typ_now->atm1)==0)
          &&(strcasecmp(build_intra->cgrp_typ_46[i].atm2,
                        build_intra->cgrp_typ_now->atm2)==0)
          &&(strcasecmp(build_intra->cgrp_typ_46[i].atm3,
                        build_intra->cgrp_typ_now->atm3)==0)
          &&(strcasecmp(build_intra->cgrp_typ_46[i].atm4,
                        build_intra->cgrp_typ_now->atm4)==0)
          &&(strcasecmp(build_intra->cgrp_typ_46[i].label[1],
                  build_intra->cgrp_typ_now->label[1])==0) 
          &&(strcasecmp(build_intra->cgrp_typ_46[i].label[2],
                  build_intra->cgrp_typ_now->label[2])==0) 
          &&(strcasecmp(build_intra->cgrp_typ_46[i].label[3],
                  build_intra->cgrp_typ_now->label[3])==0) 
          &&(strcasecmp(build_intra->cgrp_typ_46[i].label[4],
                  build_intra->cgrp_typ_now->label[4])==0) 
          &&(strcasecmp(build_intra->cgrp_typ_46[i].label[5],
                  build_intra->cgrp_typ_now->label[5])==0) 
          &&(strcasecmp(build_intra->cgrp_typ_46[i].label[6],
                  build_intra->cgrp_typ_now->label[6])==0) 
                                                         ){itype=i;}
       
       if((strcasecmp(build_intra->cgrp_typ_46[i].atm1,
                      build_intra->cgrp_typ_now->atm4)==0)
          &&(strcasecmp(build_intra->cgrp_typ_46[i].atm2,
                        build_intra->cgrp_typ_now->atm3)==0)
          &&(strcasecmp(build_intra->cgrp_typ_46[i].atm3,
                        build_intra->cgrp_typ_now->atm2)==0)
          &&(strcasecmp(build_intra->cgrp_typ_46[i].atm4,
                        build_intra->cgrp_typ_now->atm1)==0)
          &&(strcasecmp(build_intra->cgrp_typ_46[i].label[1],
                  build_intra->cgrp_typ_now->label[6])==0) 
          &&(strcasecmp(build_intra->cgrp_typ_46[i].label[2],
                  build_intra->cgrp_typ_now->label[5])==0) 
          &&(strcasecmp(build_intra->cgrp_typ_46[i].label[3],
                  build_intra->cgrp_typ_now->label[3])==0) 
          &&(strcasecmp(build_intra->cgrp_typ_46[i].label[4],
                  build_intra->cgrp_typ_now->label[4])==0) 
          &&(strcasecmp(build_intra->cgrp_typ_46[i].label[5],
                  build_intra->cgrp_typ_now->label[2])==0) 
          &&(strcasecmp(build_intra->cgrp_typ_46[i].label[6],
                  build_intra->cgrp_typ_now->label[1])==0) 
                                                         ){itype=i;}
     }       /*endfor*/
/*--------------------------------------------------------------------*/
/* D) Add space */
     if(itype>build_intra->ngrp_typ_46_max){
       build_intra->ngrp_typ_46_max += NMEM_MIN;
       build_intra->cgrp_typ_46      = (CGRP_CONS *) 
         crealloc(&(build_intra->cgrp_typ_46)[1],
              build_intra->ngrp_typ_46_max*sizeof(CGRP_CONS))-1;
     }/*endif*/

/*-------------------------------------------------------------------*/
/* E) Add a type */
     if(itype==grp_bond_con->ntyp_46+1){
       grp_bond_con->ntyp_46+=1;
       strcpy(build_intra->cgrp_typ_46[itype].atm1,
              build_intra->cgrp_typ_now->atm1);
       strcpy(build_intra->cgrp_typ_46[itype].atm2,
              build_intra->cgrp_typ_now->atm2);
       strcpy(build_intra->cgrp_typ_46[itype].atm3,
              build_intra->cgrp_typ_now->atm3);
       strcpy(build_intra->cgrp_typ_46[itype].atm4,
              build_intra->cgrp_typ_now->atm4);
       strcpy(build_intra->cgrp_typ_46[itype].label[1],
              build_intra->cgrp_typ_now->label[1]);
       strcpy(build_intra->cgrp_typ_46[itype].label[2],
              build_intra->cgrp_typ_now->label[2]);
       strcpy(build_intra->cgrp_typ_46[itype].label[3],
              build_intra->cgrp_typ_now->label[3]);
       strcpy(build_intra->cgrp_typ_46[itype].label[4],
              build_intra->cgrp_typ_now->label[4]);
       strcpy(build_intra->cgrp_typ_46[itype].label[5],
              build_intra->cgrp_typ_now->label[5]);
       strcpy(build_intra->cgrp_typ_46[itype].label[6],
              build_intra->cgrp_typ_now->label[6]);
     }/*endif*/

/*--------------------------------------------------------------------*/
/* B) Spread */
     grp_bond_con->num_46 += 1;
     grp_bond_con->j1_46[grp_bond_con->num_46] = iatm_ind1 + 
                                                     clatoms_info->natm_tot;
     grp_bond_con->j2_46[grp_bond_con->num_46] = iatm_ind2 + 
                                                     clatoms_info->natm_tot;
     grp_bond_con->j3_46[grp_bond_con->num_46] = iatm_ind3 + 
                                                     clatoms_info->natm_tot;
     grp_bond_con->j4_46[grp_bond_con->num_46] = iatm_ind4 + 
                                                     clatoms_info->natm_tot;
     grp_bond_con->jtyp_46[grp_bond_con->num_46] = itype;
     if((atommaps->ighost_flag[grp_bond_con->j1_46[grp_bond_con->num_46]]!=0)||
        (atommaps->ighost_flag[grp_bond_con->j2_46[grp_bond_con->num_46]]!=0)||
        (atommaps->ighost_flag[grp_bond_con->j3_46[grp_bond_con->num_46]]!=0)||
        (atommaps->ighost_flag[grp_bond_con->j4_46[grp_bond_con->num_46]]!=0)){
        printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
        printf("Ghost atoms are not permitted in group constraints\n");
        printf("in molecule index %d in residue index %d in file %s\n",
                jmol_typ,iresidue,file_name);
        printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
        fflush(stdout);
        exit(1);
     }/*endif*/      

    }/*endif:4x6*/
  }/*endif:igo*/

/*----------------------------------------------------------------------*/
   }/*end routine*/
/*=======================================================================*/






