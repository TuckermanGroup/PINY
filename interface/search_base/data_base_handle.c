/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
/*                                                                          */
/*                         PI_MD:                                           */
/*             The future of simulation technology                          */
/*             ------------------------------------                         */
/*                     Module: interface_hand                               */
/*                                                                          */
/* Subprograms that search the data bases                                   */
/*                                                                          */
/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

#include "standard_include.h"
#include "../typ_defs/typedefs_class.h"
#include "../typ_defs/typedefs_par.h"
#include "../typ_defs/typedefs_bnd.h"
#include "../proto_defs/proto_search_entry.h"
#include "../proto_defs/proto_inter_params_entry.h"
#include "../proto_defs/proto_inter_params_local.h"
#include "../proto_defs/proto_intra_params_local.h"
#include "../proto_defs/proto_surf_params_local.h"
#include "../proto_defs/proto_handle_entry.h"
#include "../proto_defs/proto_friend_lib_entry.h"

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
/* count_data_base: Finds length of an intra data base                      */
/*==========================================================================*/
void count_data_base(char filename[],DICT_WORD fun_dict[], int num_fun_dict,
                      int *nbase, int ibase_want)
{/*begin routine*/
  /*=======================================================================*/
  /*          Local variable declarations                                  */
  int nline,nkey,nfun_key,nbase_now,num;
  NAME fun_key;
  DICT_WORD word;
  FILE *fp;
  /*=======================================================================*/
  fp       = cfopen(filename,"r");
  nline    = 0;
  nkey     = 0;
  nfun_key = 0;
  nbase_now = 0;
  while(get_fun_key_cnt(fp,fun_key,&nline,&nfun_key,filename)){
    get_fun_key_index(fun_key,num_fun_dict,fun_dict,nline,nfun_key,
                      filename,&num);
    if(num==ibase_want){nbase_now++;}
  }/*endwhile*/
  fclose(fp);
  *nbase = nbase_now;
/*-----------------------------------------------------------------------*/
} /*end routine*/
/*==========================================================================*/

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
/* read_data_base: Reads data bases                                         */
/*==========================================================================*/
void read_data_base(char filename[],DICT_WORD fun_dict[], int num_fun_dict,
                    DATA_BASE_ENTRIES *data_base,CATM_LAB *cbase,
                    int ibase_want, int nbase)
{/*begin routine*/
  /*=======================================================================*/
  /*          Local variable declarations                                  */
  int nline,nkey,i,num,nfun_key;  
  int ifirst,ibase;
  NAME fun_key;
  DICT_WORD word;
  DICT_WORD *base_dict;
  int num_base_dict;
  FILE *fp;
  /*=======================================================================*/
  /* I) Set up the dictionaries.  */

  ifirst   = 1;
  switch(ibase_want){
    case 1: set_potinter_dict(&base_dict,&num_base_dict,ifirst); break;
    case 2: set_potbond_dict(&base_dict,&num_base_dict,ifirst); break;
    case 3: set_potbend_bnd_dict(&base_dict,&num_base_dict,ifirst); break;
    case 4: set_pottors_dict(&base_dict,&num_base_dict,ifirst); break;
    case 5: set_potonfo_dict(&base_dict,&num_base_dict,ifirst); break;
    case 7: set_potbend_bnd_dict(&base_dict,&num_base_dict,ifirst); break;
    case 8: set_potsurf_dict(&base_dict,&num_base_dict,ifirst); break;
  }/*endswitch*/

  /*=======================================================================*/
  /* II) Pack the data base into data_base and cbase.  */

  ifirst   = 0;
  fp       = cfopen(filename,"r");
  nline    = 0;
  nkey     = 0;
  nfun_key = 0;
  ibase = 0;
  while(get_fun_key(fp,fun_key,&nline,&nfun_key,filename)){
    get_fun_key_index(fun_key,num_fun_dict,fun_dict,nline,
                      nfun_key,filename,&num);
    if(num==ibase_want){
      ibase++;
      switch(ibase_want){
        case 1: set_potinter_dict(&base_dict,&num_base_dict,ifirst); break;
        case 2: set_potbond_dict(&base_dict,&num_base_dict,ifirst); break;
        case 3: set_potbend_bnd_dict(&base_dict,&num_base_dict,ifirst);
                break;
        case 4: set_pottors_dict(&base_dict,&num_base_dict,ifirst); break;
        case 5: set_potonfo_dict(&base_dict,&num_base_dict,ifirst); break;
        case 7: set_potbend_bnd_dict(&base_dict,&num_base_dict,ifirst);break;
        case 8: set_potsurf_dict(&base_dict,&num_base_dict,ifirst);break;
      }/*endswitch*/
      while(get_word(fp,&word,&nline,&nkey,nfun_key,filename)){
        put_word_dict(&word,base_dict,num_base_dict,
                      fun_key,nline,nkey,nfun_key,filename);
      }/*endwhile*/
      switch(ibase_want){
        case 1: inter_coef(base_dict,filename,fun_key,data_base,
                           cbase,ibase);
                break;
        case 2: bond_coef(base_dict,filename,fun_key,data_base,
                              cbase,ibase);
                break;
        case 3: bend_bnd_coef(base_dict,filename,fun_key,data_base,
                              cbase,ibase);
                break;
        case 4: tors_coef(base_dict,filename,fun_key,data_base,
                          cbase,ibase);
                break;
        case 5: onfo_coef(base_dict,filename,fun_key,data_base,
                          cbase,ibase);
                break;
        case 7: bend_bnd_coef(base_dict,filename,fun_key,data_base,
                              cbase,ibase);
        case 8: surf_coef(base_dict,filename,fun_key,data_base,
                          cbase,ibase);
                break;
      }/*endswitch*/
    }else{
      close_fun_key_cnt(fp,fun_key,&nline,nfun_key,filename);      
    }/*endelse*/
 }/*end while*/
  fclose(fp);

  /*=======================================================================*/
  /* III) Double and flip the base.                                        */

  if((ibase_want==1) || (ibase_want==2) || (ibase_want==5)){
    for(i=1;i<=nbase;i++){
      strcpy(cbase[(i+nbase)].atm2,cbase[i].atm1);
      strcpy(cbase[(i+nbase)].atm1,cbase[i].atm2);
      if(strcasecmp(cbase[i].atm1,cbase[i].atm2)!=0){
        strcpy(cbase[(i+nbase)].label,cbase[i].label);
      }else{
        strcpy(cbase[(i+nbase)].label,"XSYMX");
      }/*endif*/
    }/*endfor*/
  }/*endif*/

  if( (ibase_want==8) ){
    for(i=1;i<=nbase;i++){
      strcpy(cbase[(i+nbase)].atm2,cbase[i].atm1);
      strcpy(cbase[(i+nbase)].atm1,cbase[i].atm2);
      strcpy(cbase[(i+nbase)].label,"XNOSYMX"); /* no atm switches */
    }/*endfor*/
  }/*endif*/

  if((ibase_want==3) || (ibase_want==7)){
    for(i=1;i<=nbase;i++){
      strcpy(cbase[(i+nbase)].atm3,cbase[i].atm1);
      strcpy(cbase[(i+nbase)].atm2,cbase[i].atm2);
      strcpy(cbase[(i+nbase)].atm1,cbase[i].atm3);
      if(strcasecmp(cbase[i].atm1,cbase[i].atm3)!=0){
        strcpy(cbase[(i+nbase)].label,cbase[i].label);
      }else{
        strcpy(cbase[(i+nbase)].label,"XSYMX");
      }/*endif*/
    }/*endfor*/
  }/*endif*/

  if(ibase_want==4){
    for(i=1;i<=nbase;i++){
      strcpy(cbase[(i+nbase)].atm4,cbase[i].atm1);
      strcpy(cbase[(i+nbase)].atm3,cbase[i].atm2);
      strcpy(cbase[(i+nbase)].atm2,cbase[i].atm3);
      strcpy(cbase[(i+nbase)].atm1,cbase[i].atm4);
      if((strcasecmp(cbase[i].atm1,cbase[i].atm4)==0)
      && (strcasecmp(cbase[i].atm2,cbase[i].atm3)==0)){
        strcpy(cbase[(i+nbase)].label,"XSYMX");
      }else{
        strcpy(cbase[(i+nbase)].label,cbase[i].label);
      }/*endif*/
    }/*endfor*/
  }/*endif*/

  /*=======================================================================*/
  /* IV) Free the memory.                                                  */
  cfree(&base_dict[1]);
/*-------------------------------------------------------------------------*/
}  /*end routine*/
/*=========================================================================*/
