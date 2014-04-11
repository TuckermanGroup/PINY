/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
/*                                                                          */
/*                         PI_MD:                                           */
/*             The future of simulation technology                          */
/*             ------------------------------------                         */
/*                     Module: set_intra_parms.c                            */
/*                                                                          */
/* This subprogram set molecular and intramolecular data sets               */
/*                                                                          */
/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
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
/*  set_atm_params:set up the atom params                                   */
/*==========================================================================*/

void set_atm_params(DICT_WORD atm_dict[],int num_atm_dict,
                    char fun_key[],char filename[],
                    int jmol_typ,int jres, int jres_off,ATOMMAPS *atommaps,
                    BUILD_INTRA *build_intra,CLATOMS_INFO *clatoms_info,
                    GHOST_ATOMS *ghost_atoms)

/*=======================================================================*/
{/*begin routine*/
/*=======================================================================*/
/*  Local Variables */

  NAME site;
  int i,iii,num,num2,num3,num4,iatm_ind,imask,index,num1; 
  int itype,iatm_ind_sm,ifound;
  double dnum;

/*=======================================================================*/
/* I) Check for missing key words      */

  for(i=1;i<=num_atm_dict;i++){
    if(atm_dict[i].iuset==0 && atm_dict[i].key_type==1){
      keyword_miss(atm_dict,filename,fun_key,i);}
  }/*endfor*/

/*=======================================================================*/
/* II) Get atm index                   */

  sscanf(atm_dict[2].keyarg,"%d",&iatm_ind);
  index = 2;
  if(iatm_ind>build_intra->natmind_1res_now||iatm_ind<0){
             keyarg_barf(atm_dict,filename,fun_key,index);}
  imask = build_intra->mask_atm[iatm_ind];

/*=======================================================================*/
/*=======================================================================*/

  if(imask>0){

/*=======================================================================*/
/* III) Get rejiggered atm index                                         */
/*-----------------------------------------------------------------------*/

    iatm_ind = (build_intra->index_atm)[iatm_ind];
    if(iatm_ind>build_intra->natm_1res_now||iatm_ind<0){
             keyarg_barf(atm_dict,filename,fun_key,index);}

    build_intra->iatm_ind_chk[iatm_ind]++;
    iatm_ind_sm = iatm_ind;
    iatm_ind = iatm_ind + clatoms_info->natm_tot;

/*=======================================================================*/
/* IV) Fill the dictionary with words */
/*-----------------------------------------------------------------------*/
     /*  3) \mass{} */
     sscanf(atm_dict[3].keyarg,"%lg",&dnum);
     index = 3;
     if(dnum<=0.0){
       keyarg_barf(atm_dict,filename,fun_key,index);
     }/*endif*/
     clatoms_info->mass[iatm_ind] = dnum; 
/*--------------------------------------------------------------------*/
/*  4) \charge{} */
     sscanf(atm_dict[4].keyarg,"%lg",&dnum);
     clatoms_info->q[iatm_ind] = dnum; 
/*--------------------------------------------------------------------*/
/*  5) \alpha_pol{} */
     sscanf(atm_dict[5].keyarg,"%lg",&dnum);
     index = 5;
     if(dnum<0.0){
       keyarg_barf(atm_dict,filename,fun_key,index);
     }/*endif*/
     clatoms_info->alp_pol[iatm_ind] = dnum; 
/*--------------------------------------------------------------------*/
/*  6) \b_neut{} */
     sscanf(atm_dict[6].keyarg,"%lg",&dnum);
     clatoms_info->b_neut[iatm_ind] = dnum; 
/*--------------------------------------------------------------------*/
/*  7) \valence{} */
     sscanf(atm_dict[7].keyarg,"%d",&num);
     build_intra->bond_site[iatm_ind_sm].valence = num;
     index = 7;
     if(num<0||num>MAX_VALENCE){
       keyarg_barf(atm_dict,filename,fun_key,index);}
/*--------------------------------------------------------------------*/
/*  8) \improper_def{} */
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
      if(num1<0||num1>3||num2<0||num2>3||num3<0||num3>3||num4<0||num4>3||
         num1==num2||num1==num3||num2==num3||num1==num4||num2==num4||
         num3==num4){
       keyarg_barf(atm_dict,filename,fun_key,index);}}
/*--------------------------------------------------------------------*/
/*  9) \bond_site_1{} */
     strcpy(build_intra->strip1,atm_dict[9].keyarg);
     parse_bond_site(build_intra->strip1,build_intra->strip2,site,&num2,&num3);
     strcpy(build_intra->bond_site[iatm_ind_sm].bond_site_name[1],site);
     build_intra->bond_site[iatm_ind_sm].branch_1[1]     = num2;
     build_intra->bond_site[iatm_ind_sm].branch_2[1] = num3;
     index = 9;
     if(strlen(site)==0||num2<-1||num3<-1||num2>MAX_VALENCE||num3>MAX_VALENCE){
       keyarg_barf(atm_dict,filename,fun_key,index);
     }/*endif*/
/*--------------------------------------------------------------------*/
/*  10) \bond_site_2{} */
     strcpy(build_intra->strip1,atm_dict[10].keyarg);
     parse_bond_site(build_intra->strip1,build_intra->strip2,site,&num2,&num3);
     strcpy(build_intra->bond_site[iatm_ind_sm].bond_site_name[2],site);
     build_intra->bond_site[iatm_ind_sm].branch_1[2]     = num2;
     build_intra->bond_site[iatm_ind_sm].branch_2[2] = num3;
     index = 10;
     if(strlen(site)==0||num2<-1||num3<-1||num2>MAX_VALENCE||num3>MAX_VALENCE){
       keyarg_barf(atm_dict,filename,fun_key,index);
     }
/*---------------------------------------------------------------------*/
/*  11) \bond_site_3{} */
     strcpy(build_intra->strip1,atm_dict[11].keyarg);
     parse_bond_site(build_intra->strip1,build_intra->strip2,site,&num2,&num3);
     strcpy(build_intra->bond_site[iatm_ind_sm].bond_site_name[3],site);
     build_intra->bond_site[iatm_ind_sm].branch_1[3]     = num2;
     build_intra->bond_site[iatm_ind_sm].branch_2[3] = num3;
     index = 11;
     if(strlen(site)==0||num2<-1||num3<-1||num2>MAX_VALENCE||num3>MAX_VALENCE){
       keyarg_barf(atm_dict,filename,fun_key,index);}
/*---------------------------------------------------------------------*/
/*  12) \bond_site_4{} */
     strcpy(build_intra->strip1,atm_dict[12].keyarg);
     parse_bond_site(build_intra->strip1,build_intra->strip2,site,&num2,&num3);
     strcpy(build_intra->bond_site[iatm_ind_sm].bond_site_name[4],site);
     build_intra->bond_site[iatm_ind_sm].branch_1[4]     = num2;
     build_intra->bond_site[iatm_ind_sm].branch_2[4]     = num3;
     index = 12;
     if(strlen(site)==0||num2<-1||num3<-1||num2>MAX_VALENCE||num3>MAX_VALENCE){
       keyarg_barf(atm_dict,filename,fun_key,index);}
/*--------------------------------------------------------------------*/
/*  13) \cp_valence_up{} */
     sscanf(atm_dict[13].keyarg,"%d",&num);
     clatoms_info->cp_vlnc_up[iatm_ind] = num; 

  
/*--------------------------------------------------------------------*/
/*  14) \cp_valence_dn{} */
     sscanf(atm_dict[14].keyarg,"%d",&num);
     clatoms_info->cp_vlnc_dn[iatm_ind] = num; 

    /* if cp_valence_up is assigned but down is not assign dn = to up */
     if(atm_dict[13].iuset == 1 && atm_dict[14].iuset == 0){
      clatoms_info->cp_vlnc_dn[iatm_ind] = clatoms_info->cp_vlnc_up[iatm_ind];
     }/*endif*/

/*------------------------------------------------------------------------*/ 
/*  15) \cp_atom{} */

     ifound = 0;
     if(strcasecmp(atm_dict[15].keyarg,"yes")==0){
       ifound=1;clatoms_info->cp_atm_flag[iatm_ind] = 1;
     }
     if(strcasecmp(atm_dict[15].keyarg,"no")==0){
       ifound=1;clatoms_info->cp_atm_flag[iatm_ind] = 0;
     }

/*------------------------------------------------------------------------*/ 
/*  16) /def_ghost1{} */

     atommaps->ighost_flag[iatm_ind]  = 0;
     if(atm_dict[16].iuset!=0){
      set_ghost(ghost_atoms,clatoms_info,atommaps,build_intra,
               atm_dict,num_atm_dict,
               filename,fun_key,iatm_ind);
     }/*endif*/

/*------------------------------------------------------------------------*/ 
/*  17) \label{} */
     index  = 16+NCOEF_GHOST_MAX;
     ifound = 0;
     if(strcasecmp(atm_dict[index].keyarg,"standard")==0){
       ifound=1;atommaps->atom_label[iatm_ind] = 0;
     }
     if(strcasecmp(atm_dict[index].keyarg,"backbone")==0){
       ifound=1;atommaps->atom_label[iatm_ind] = 1;
     }
     if(strcasecmp(atm_dict[index].keyarg,"sidechain")==0){
       ifound=1;atommaps->atom_label[iatm_ind] = 2;
     }


    if(ifound==0){
     keyarg_barf(atm_dict,filename,fun_key,index);    
    }/*endif*/
/*------------------------------------------------------------------------*/ 
/*====================================================================*/
/* IV) Do up atom types  */

     itype = atommaps->natm_typ + 1;
     for(i=1;i<=atommaps->natm_typ;i++){
       if(strcasecmp(atm_dict[1].keyarg,atommaps->atm_typ[i])==0)
         itype=i;
     } /*endfor*/
     if(itype>build_intra->natm_typ_max){
       build_intra->natm_typ_max+=NMEM_MIN;
       atommaps->atm_typ = (NAME *) 
         crealloc(&(atommaps->atm_typ[1]),
                  (build_intra->natm_typ_max)*sizeof(NAME))-1;
     }/*endif*/
     if(itype==atommaps->natm_typ+1){
       atommaps->natm_typ+=1;
       strcpy(atommaps->atm_typ[itype],atm_dict[1].keyarg);
     }/*endif*/

/*====================================================================*/
/* V) Do up types                                                     */

     atommaps->iatm_mol_typ[iatm_ind] = jmol_typ;
     atommaps->iatm_mol_num[iatm_ind] = 1;
     atommaps->iatm_res_num[iatm_ind] = jres;
     atommaps->iatm_atm_typ[iatm_ind] = itype;
     atommaps->iatm_res_typ[iatm_ind] = 
                        atommaps->ires_typ_jres_jmol_typ[jres+jres_off];

/*====================================================================*/
/*====================================================================*/

   } /*endif (imask>0) */

/*======================================================================*/
}  /*end routine*/
/*==========================================================================*/




/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
/*  set_ghost_params:set up the ghost params                                */
/*==========================================================================*/

void set_ghost(GHOST_ATOMS *ghost_atoms,CLATOMS_INFO *clatoms_info,
               ATOMMAPS *atommaps,
               BUILD_INTRA *build_intra,DICT_WORD *atm_dict,
               int num_atm_dict,
               char *filename,char *fun_key,int iatm_ind)

/*==========================================================================*/
{ /*begin routine*/
/*==========================================================================*/
/*       Local Variables */

   int nstrip=NCOEF_GHOST_MAX,imask;
   int numg[NCOEF_GHOST_MAX1],index,i;
   double anumg[NCOEF_GHOST_MAX1];
   int nghost_old,ncomp_old,nghost;
   int nghost_new,ncomp_new,ncomp,iii;

/*==========================================================================*/
/* I) Strip out indices of Ghost atoms */

   ncomp = 0;index=16; /*INDEX for GHOST*/
   for(i=1;i<=nstrip;i++){
     if(atm_dict[index].iuset==0){break;}
      ncomp = i;
      strcpy(build_intra->strip1,atm_dict[index].keyarg);
      parse_ghost(build_intra->strip1,build_intra->strip2,&numg[i],&anumg[i]);
      index++;
   }/*endfor*/

/*==========================================================================*/
/* II) Check Mask of indices making up Ghost atoms and rejigger them*/

   index = 16;
   for(i=1;i<=ncomp;i++){
      if((numg[i]<=0)||(numg[i]>build_intra->natmind_1res_now)){
                     keyarg_barf(atm_dict,filename,fun_key,index);}
      imask = build_intra->mask_atm[numg[i]];
      if(imask==0){
        printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
        printf("You have managed to nuke one of the atoms making \n");
        printf("up a ghost atom without nuking the ghost atom, \n");
        printf("itself, in residue file  %s\n",filename);
        printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
        fflush(stdout);
        exit(1);
      }/*endif*/
      numg[i] = build_intra->index_atm[numg[i]];
      if(numg[i]>build_intra->natm_1res_now||numg[i]<0){
             keyarg_barf(atm_dict,filename,fun_key,index);}
      numg[i] += clatoms_info->natm_tot;
      index++;
   }/*endfor*/

/*==========================================================================*/
/* III) Memory allocation                                                   */

   if(ncomp>0){
       ghost_atoms->nghost_tot++;
       build_intra->nghost_now++;
       nghost = ghost_atoms->nghost_tot;
       if(nghost > build_intra->nghost_tot_max){
           nghost_old  = build_intra->nghost_tot_max;
           ncomp_old   = NCOEF_GHOST_MAX;
           build_intra->nghost_tot_max+=NMEM_MIN;
           nghost_new  = build_intra->nghost_tot_max;
           ncomp_new   = NCOEF_GHOST_MAX;

           ghost_atoms->ighost_map = (int *) crealloc(
                                     &(ghost_atoms->ighost_map)[1],
                                       nghost_new*sizeof(int))-1;
           ghost_atoms->natm_comp  = (int *) crealloc(  
                                     &(ghost_atoms->natm_comp)[1],
                                       nghost_new*sizeof(int))-1;
           ghost_atoms->ighost_map = (int *) crealloc(
                                   &(ghost_atoms->ighost_map)[1],
                                     nghost_new*sizeof(int))-1;
           ghost_atoms->iatm_comp  = creall_int_mat(ghost_atoms->iatm_comp,
                                                    1,ncomp_old,1,nghost_old,
                                                    1,ncomp_new,1,nghost_new);
           ghost_atoms->coef       = creall_mat(ghost_atoms->coef,
                                                1,ncomp_old,1,nghost_old,
                                                1,ncomp_new,1,nghost_new);
        }/*endif*/
   }/*endif*/

   if(ncomp>NCOEF_GHOST_MAX){
        printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
        printf("You have managed to use too many ghost_defs!!!\n");
        printf("Technical support was too tired to do the  \n");
        printf("heinous work necessary to support an arbitray number.\n");
        printf("Wait for the next upgrade you power tool.      \n");
        printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
        fflush(stdout);
        exit(1);
   }/*endif*/

/*==========================================================================*/
/* IV) Assign Ghosts parameters                                             */

   if(ncomp>0){         
     atommaps->ighost_flag[iatm_ind]  = nghost;
     ghost_atoms->natm_comp_max = MAX(ghost_atoms->natm_comp_max,ncomp);
     ghost_atoms->ighost_map[nghost] = iatm_ind;
     ghost_atoms->natm_comp[nghost]  = ncomp;
     for(i=1;i<=ncomp;i++){
       ghost_atoms->iatm_comp[i][nghost] = numg[i];
       ghost_atoms->coef[i][nghost]      = anumg[i];
     }/*endfor*/
   }/*endif*/

/*--------------------------------------------------------------------------*/
}/*end routine */
/*==========================================================================*/








