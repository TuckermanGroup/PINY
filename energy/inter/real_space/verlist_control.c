/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
/*                         PI_MD:                                           */
/*             The future of simulation technology                          */
/*             ------------------------------------                         */
/*                     Module: verlist_control.c                            */
/*                                                                          */
/* These routines generate a verlet neighbor list using either a            */
/* link_lst or no list at all to provide possible interaction sets          */
/* A root-branch update option is provided for each generation procedure    */
/*==========================================================================*/



#include "standard_include.h"
#include "../typ_defs/typedefs_class.h"
#include "../typ_defs/typedefs_gen.h"
#include "../typ_defs/typedefs_bnd.h"
#include "../proto_defs/proto_real_space_local.h"
#include "../proto_defs/proto_friend_lib_entry.h"
#include "../proto_defs/proto_communicate_wrappers.h"



/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void verlist_control(CLATOMS_INFO *clatoms_info,CLATOMS_POS *clatoms_pos,
                     NBR_LIST *nbr_list,EXCL *excl,
                     ATOMMAPS *atommaps,CELL *cell,INTRA_SCR *intra_scr,
                     FOR_SCR *for_scr,TIMEINFO *timeinfo,
                     INTERACT *interact,int error_check_on,
                     CLASS_COMM_FORC_PKG *class_comm_forc_pkg)
/*======================================================================*/
/*      Begin Routine                                                   */
 {/*Begin Routine*/
/*======================================================================*/
/*      Local variable declarations                                     */

#include "../typ_defs/typ_mask.h"

  int i,isave,iii;
  int nver,nver_res,npairs_tmp=0,npairs_res_tmp=0;
  int isuccess;
  int brnch_root_list_opt   = nbr_list->brnch_root_list_opt; 
  int np_forc               = class_comm_forc_pkg->num_proc;
  int myid_forc             = class_comm_forc_pkg->myid;
  MPI_Comm comm_forc        = class_comm_forc_pkg->comm;

/*=======================================================================*/
/* I) Save the present positions                                         */
  
  for(i=1;i<=clatoms_info->natm_tot;i++){
    nbr_list->x0[i]= clatoms_pos->x[i];
    nbr_list->y0[i]= clatoms_pos->y[i];
    nbr_list->z0[i]= clatoms_pos->z[i];
  }/*endfor*/

/*=====================================================================*/
/* II) Generate Verlet list using no list                              */
  
  if(nbr_list->verlist.nolst_ver_update == 1){

    nbr_list->verlist.iver_fill=1;nbr_list->verlist.iver_count=1;   
                                                /* both fill and count */
    if(brnch_root_list_opt==0){
      nolst_ver_gen(clatoms_info,clatoms_pos,
                    for_scr,atommaps,cell,interact,
                    timeinfo,nbr_list,excl,intra_scr,&nver,&nver_res,
                    error_check_on,class_comm_forc_pkg);
    }else{

      nolst_ver_gen_root(clatoms_info,clatoms_pos,
                         for_scr,atommaps,cell,interact,
                         timeinfo,nbr_list,excl,intra_scr,&nver,&nver_res,
                         error_check_on,class_comm_forc_pkg);
    }/*endif*/

  }/*endif:no_lst_ver*/

/*=====================================================================*/
/* II) Generate Verlet list using lnk lst: 1st generates the lnk_list  */
/*                     then updates verlist using a padding scheme     */

  if(nbr_list->verlist.lnk_ver_update == 1){

    isave = timeinfo->int_res_ter;timeinfo->int_res_ter = 0; 
    make_lnk_lst(clatoms_info,clatoms_pos,
                 nbr_list,cell,for_scr,intra_scr,timeinfo,excl,interact,
                 error_check_on); 
    timeinfo->int_res_ter=isave;
    isuccess = 0;
    while(isuccess==0){
      if(brnch_root_list_opt==0){
        lnk_ver_gen_pad(clatoms_info,clatoms_pos,
                        for_scr,atommaps,cell,interact,
                        timeinfo,nbr_list,excl,intra_scr,&nver,&nver_res,
                        &isuccess,error_check_on,class_comm_forc_pkg);
      }else{
        lnk_ver_gen_pad_root(clatoms_info,clatoms_pos,
                        for_scr,atommaps,cell,interact,
                        timeinfo,nbr_list,excl,intra_scr,&nver,&nver_res,
                        &isuccess,error_check_on,class_comm_forc_pkg);
      }/*endif*/
      if(isuccess==0){
        printf("\n====================================\n"); 
        printf("New Verlet list padding = %d myid = %d\n",
                nbr_list->verlist.jver_pad,myid_forc);
        printf("====================================\n\n"); 
      }/*endif*/
    }/*endwhile: padding*/

  }/*endif:lnkver*/

/*=====================================================================*/
/* III) Assign the number of pairs                                     */

  nbr_list->verlist.nver_lst_now = nver;
  nbr_list->verlist.nver_lst_now_res = nver_res;

/*----------------------------------------------------------------------*/
  }/*end routine */
/*==========================================================================*/




/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void nolst_ver_gen(CLATOMS_INFO *clatoms_info,CLATOMS_POS *clatoms_pos,
                   FOR_SCR *for_scr,ATOMMAPS *atommaps,
                   CELL *cell,INTERACT *interact,
                   TIMEINFO *timeinfo, NBR_LIST *nbr_list, EXCL *excl,
                   INTRA_SCR *intra_scr,int *nver,int *nver_res,
                   int error_check_on,
                   CLASS_COMM_FORC_PKG *class_comm_forc_pkg)

/*==========================================================================*/
/*                      Brief description                                   */
/*                                                                          */
/* All possible interactions are sent to verlist create routine in bundles  */
/* length of nlen. Break points are tabluated                               */
/*                                                                          */
/*==========================================================================*/
{/*Begin Routine*/
/*==========================================================================*/
/*            Local Variable declarations               */

    int ipart, jpart,i1;           /* ith particle, jth particle          */
    int intact_tot;             /* total number interactions           */
    int intact_tot_tmp;             /* total number interactions           */
    int intact_save;            /* total number interactions           */
    int num_call;               /*  number of force calls              */
    int num_call_tmp;               /*  number of force calls              */
    int num_call_save;          /*  number of force calls              */
    int lower,upper;            /* lower and upper limits on for loop  */
    int ifor_call;              /* force call flag                     */
    int intact_now;             /* interactions now                    */
    int jind_off;
    int iswitch;
    int nver_now,nver_res_now;
    int iii,mtemp;
    int jstart,jend,joff,k; 
    int lst_typ,irealloc;
    int ipart_st,nskip;
    int inc1,inc2;
/* Local pointers */

    int excl_nlst              = excl->nlst;
    int *for_scr_intact        = &(for_scr->intact);
    int *for_scr_num_brk       = &(for_scr->num_brk);
    int *for_scr_num_brk_i     = for_scr->num_brk_i;
    int *for_scr_iexcl         = for_scr->iexcl;
    int *for_scr_i_index       = for_scr->i_index;  
    int *for_scr_j_index       = for_scr->j_index;  
    int *excl_j_off            = excl->j_off; 
    int *excl_num              = excl->num;
    int *excl_j                = excl->j;
    int *nbr_list_nter         = nbr_list->verlist.nter;
    int *nbr_list_jver_off     = nbr_list->verlist.jver_off;
    int *nbr_list_nter_res     = nbr_list->verlist.nter_res;
    int *nbr_list_jver_off_res = nbr_list->verlist.jver_off_res;

    int int_res_ter            = timeinfo->int_res_ter;
    int iver_count             = nbr_list->verlist.iver_count;
    int iver_fill              = nbr_list->verlist.iver_fill;
    int natm_tot               = clatoms_info->natm_tot;  
    int nlen                   = for_scr->nlen;
    int nmem_min_lst           = nbr_list->verlist.nmem_min_lst;
    int iver_init              = nbr_list->verlist.iver_init;
    double mem_safe            = nbr_list->verlist.mem_safe;
    int np_forc                = class_comm_forc_pkg->num_proc;
    int myid_forc              = class_comm_forc_pkg->myid;
    MPI_Comm comm_forc         = class_comm_forc_pkg->comm;

/*==========================================================================*/
/* I) Count the interactions                                               */

   intact_save      = (((natm_tot)*(natm_tot-1))/2-(excl_nlst));
   num_call_save    = intact_save/nlen;
   if( (intact_save % nlen)!=0 ){num_call_save += 1;}

/*==========================================================================*/
/* II) Initialize                                                           */

   nver_now            = 0;
   nver_res_now        = 0;
   (*for_scr_intact)   = 0;
   (*for_scr_num_brk)  = 0;
   num_call            = 0;
   intact_tot          = 0;
   iswitch             = 0;
   lst_typ             = 1;

   for(jpart=1;jpart<=natm_tot;jpart++){
      nbr_list_nter[jpart]     = 0;
    }/*endfor*/
   if(int_res_ter==1){
     for(jpart=1;jpart<=natm_tot;jpart++){
       nbr_list_nter_res[jpart]     = 0;
     }/*endfor*/
   }/*endif*/

   if(iver_count==1){
     for(jpart=1;jpart<=natm_tot;jpart++){
       nbr_list_jver_off[jpart] = 0;
     }/*endfor*/
     if(int_res_ter==1){
       for(jpart=1;jpart<=natm_tot;jpart++){
         nbr_list_jver_off_res[jpart] = 0;
       }/*endfor*/
    }/*endif*/    
   }/*endif*/    

/*==========================================================================*/
/* III) Loop over the possible interactions:                                */

   ipart_st = 2+myid_forc;
   nskip    = np_forc;
   inc1=0;
   inc2=0;
   for(ipart=ipart_st;ipart<=natm_tot;ipart+=nskip){

/*==========================================================================*/
/* A) Check for exclusions and create an interaction list                   */

      if(excl_num[ipart]!=0){ 
        lower = excl_j_off[ipart]+1;
        upper = excl_num[ipart]+excl_j_off[ipart];
        jstart = 1;
        intact_now = 0;
        for(k=lower;k<=upper;k++){
          jend = excl_j[k]-1;
          joff = intact_now-jstart+1;
          for(jpart=jstart;jpart<=jend;jpart++){
            for_scr_iexcl[(jpart+joff)]  = jpart;          
          }/*endfor*/
          intact_now += MAX((jend-jstart+1),0);
          jstart = excl_j[k]+1;
        }/*endfor*/
          jend = ipart-1;
          joff = intact_now-jstart+1;
          for(jpart=jstart;jpart<=jend;jpart++){
            for_scr_iexcl[(jpart+joff)]  = jpart;          
          }        
          intact_now += MAX((jend-jstart+1),0);
      }else{
        for(jpart=1;jpart<=ipart-1;jpart++){
           for_scr_iexcl[jpart] = jpart;
        }/*endfor*/
        intact_now = (ipart-1);
      }/*endif*/ 

/*==========================================================================*/
/*  B) Add interactions to the force routine interaction list               */
/*     in chunks no greater than nlen. Keep track of the break points       */

      lower  = 1;
      upper  = 0;
      while(upper!=intact_now){
        upper = intact_now;
        if(upper-lower+1+(*for_scr_intact)>(nlen)){
           upper = nlen - (*for_scr_intact) + lower-1;
        }/*endif*/
        jind_off = (*for_scr_intact)-lower+1;
        for(jpart=lower;jpart<=upper;jpart++){
           mtemp = jpart+jind_off;
           for_scr_j_index[mtemp]=for_scr_iexcl[jpart];
           for_scr_i_index[mtemp]=ipart;
        }/*endfor*/
        (*for_scr_intact) += upper-lower+1;
        (*for_scr_num_brk)++;
        for_scr_num_brk_i[(*for_scr_num_brk)] = upper-lower+1;
        lower = upper+1;
      
/*==========================================================================*/
/*  C) If enough interactions have accumulated (or you are done)            */
/*       make a verlist create call and then update and reinitialize        */
/*       appropriate counters                                               */

/*-------------------------------------------------------------------------*/
/*  i)  List creat conditions                                              */

        ifor_call = 0;
        if(((*for_scr_intact)==(nlen))){ifor_call=1;}
        if((ipart==(natm_tot))&&(upper==intact_now)&&
          ((*for_scr_intact)!=0)   ){
            ifor_call=1;
        }/*endif*/

/*-------------------------------------------------------------------------*/
/*  ii)  List create call and memory checks                                */

        if(ifor_call==1){
          inc1++;
          verlist_create(clatoms_info,clatoms_pos,
                         intra_scr,for_scr,atommaps,cell,
                         nbr_list,timeinfo,interact,iswitch,
                         &nver_now,&nver_res_now,lst_typ);
          num_call        += 1;
          intact_tot      += (*for_scr_intact);
          (*for_scr_intact)  = 0;
          (*for_scr_num_brk) = 0;

          ver_list_mem_check(nver_now,&(nbr_list->verlist.nver_lst),
                             iver_fill,nmem_min_lst,iver_init,mem_safe,
                             &irealloc,myid_forc);
          if(irealloc==1){
            nbr_list->verlist.jter = (list_int *)
                   crealloc(&(nbr_list->verlist.jter)[1],
                             (nbr_list->verlist.nver_lst)*sizeof(list_int))-1;
          }/*endif*/

          if(int_res_ter==1){
            ver_list_mem_check(nver_res_now,&(nbr_list->verlist.nver_lst_res),
                               iver_fill,nmem_min_lst,iver_init,mem_safe,
                               &irealloc,myid_forc);
            if(irealloc==1){
              nbr_list->verlist.jter_res = (list_int *)
                    crealloc(&(nbr_list->verlist.jter_res)[1],
                      (nbr_list->verlist.nver_lst_res)*sizeof(list_int))-1;
            }/*endif*/
          }/*endif:RESPA*/

        }/*endif:for_call*/

      }/*endwhile:getting interactions*/

   }/*endfor:each particles interactions*/

/*==========================================================================*/
/* IV) Final list create call and memory check */

  if((*for_scr_intact)>0){
          inc2++;
    verlist_create(clatoms_info,clatoms_pos,
                   intra_scr,for_scr,atommaps,cell,
                   nbr_list,timeinfo,interact,iswitch,
                   &nver_now,&nver_res_now,lst_typ);
    num_call          += 1;
    intact_tot        += (*for_scr_intact);

    ver_list_mem_check(nver_now,&(nbr_list->verlist.nver_lst),
                       iver_fill,nmem_min_lst,iver_init,mem_safe,
                       &irealloc,myid_forc);
    if(irealloc==1){
         nbr_list->verlist.jter = (list_int *)
                crealloc(&(nbr_list->verlist.jter)[1],
                       (nbr_list->verlist.nver_lst)*sizeof(list_int))-1;
    }/*endif*/

    if(int_res_ter==1){
      ver_list_mem_check(nver_res_now,&(nbr_list->verlist.nver_lst_res),
                         iver_fill,nmem_min_lst,iver_init,mem_safe,
                         &irealloc,myid_forc);
      if(irealloc==1){
         nbr_list->verlist.jter_res = (list_int *)
                crealloc(&(nbr_list->verlist.jter_res)[1],
                       (nbr_list->verlist.nver_lst_res)*sizeof(list_int))-1;
      }/*endif*/
    }/*endif:RESPA*/

  }/*endif:last force call*/

/*==========================================================================*/
/* V) Assign present size of list to appropriate variable                   */

   (*nver) = nver_now;
   (*nver_res) = nver_res_now;

/*==========================================================================*/
/* VI) Calculate jver off                                                   */

   for(jpart=2;jpart<=natm_tot;jpart++){
      mtemp = jpart-1;
      nbr_list_jver_off[jpart] = nbr_list_jver_off[mtemp]
                               + nbr_list_nter[mtemp];
   }/*endfor*/
   if((int_res_ter==1)){
     for(jpart=2;jpart<=natm_tot;jpart++){
        mtemp = jpart-1;
        nbr_list_jver_off_res[jpart] = nbr_list_jver_off_res[mtemp]
                                     + nbr_list_nter_res[mtemp];
     }/*endfor*/
   }/*endif*/

/*==========================================================================*/
/* VI) Check                                                                */

    
 if(np_forc==1){
    if(intact_tot!=intact_save){
       printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
       printf("Internal Error:\n ");
       printf("Incorrect number of interactions calculated\n ");
       printf("   in routine nolst_ver_gen \n ");
       printf("%d vs %d \n",intact_tot,intact_save);
       printf("Contact product support \n ");
       printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
       fflush(stdout);
       exit(1);
   }/*endif*/
   if(num_call!=num_call_save){
       printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
       printf("Internal Error:\n ");
       printf("Incorrect number of force calls\n ");
       printf("   in routine nolst_ver_gen \n ");
       printf("%d vs %d \n",num_call,num_call_save);
       printf("Contact product support \n ");
       printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
       fflush(stdout);
       exit(1);
   }/*endif*/
 }/*endif*/

/*--------------------------------------------------------------------------*/
    }/*end routine */
/*==========================================================================*/





/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void nolst_ver_gen_root(CLATOMS_INFO *clatoms_info,CLATOMS_POS *clatoms_pos,
                   FOR_SCR *for_scr,ATOMMAPS *atommaps,
                   CELL *cell,INTERACT *interact,
                   TIMEINFO *timeinfo, NBR_LIST *nbr_list, EXCL *excl,
                   INTRA_SCR *intra_scr,int *nver,int *nver_res,
                   int error_check_on,
                   CLASS_COMM_FORC_PKG *class_comm_forc_pkg)

/*==========================================================================*/
/*                      Brief description                                   */
/*                                                                          */
/* All possible interactions are sent to verlist create routine in bundles  */
/* length of nlen. Break points are tabluated                               */
/*                                                                          */
/*==========================================================================*/
{/*Begin Routine*/
/*==========================================================================*/
/*            Local Variable declarations               */

#include "../typ_defs/typ_mask.h"

    int ipart, jpart,i1;           /* ith particle, jth particle          */
    int intact_tot;             /* total number interactions           */
    int intact_tot_tmp;             /* total number interactions           */
    int intact_save;            /* total number interactions           */
    int num_call;               /*  number of force calls              */
    int num_call_tmp;               /*  number of force calls              */
    int num_call_save;          /*  number of force calls              */
    int lower,upper;            /* lower and upper limits on for loop  */
    int ifor_call;              /* force call flag                     */
    int intact_now;             /* interactions now                    */
    int jind_off;
    int iswitch;
    int nver_now,nver_res_now;
    int iii,mtemp;
    int jstart,jend,joff,k; 
    int lst_typ,i,j,ngo,ires_flag;
    int my_root,jind,ind,my_ind,nlist_start,nlist_start_res;
    int irealloc;
    int ipart_st,nskip;    

/* Local pointer */
    int excl_nlst              = excl->nlst_root;
    int *for_scr_intact        = &(for_scr->intact);
    int *for_scr_num_brk       = &(for_scr->num_brk);
    int *for_scr_num_brk_i     = for_scr->num_brk_i;
    int *for_scr_iexcl         = for_scr->iexcl;
    int *for_scr_i_index       = for_scr->i_index;  
    int *for_scr_j_index       = for_scr->j_index;  
    int *excl_j_off            = excl->j_off_root; 
    int *excl_num              = excl->num_root;
    int *excl_j                = excl->j_root;
    int *nbr_list_nter         = nbr_list->verlist.nter;
    int *nbr_list_jver_off     = nbr_list->verlist.jver_off;
    int *nbr_list_nter_res     = nbr_list->verlist.nter_res;
    int *nbr_list_jver_off_res = nbr_list->verlist.jver_off_res;
    list_int *nbr_list_jter;
    list_int *nbr_list_jter_res;

    int *brnch_atm_root        = nbr_list->brnch_root.brnch_atm_root;
    int *brnch_atm_list        = nbr_list->brnch_root.brnch_atm_list;
    int *natm_add              = nbr_list->brnch_root.natm_add;
    int *iatm_add              = nbr_list->brnch_root.iatm_add;
    int *iatm_add_off          = nbr_list->brnch_root.iatm_add_off;
    int *root_atm_list         = nbr_list->brnch_root.root_atm_list;
    int *root_atm_map          = nbr_list->brnch_root.root_atm_map;
    int *brnch_atm_map         = nbr_list->brnch_root.brnch_atm_map;
    int **ibrnch_of_root       = nbr_list->brnch_root.ibrnch_of_root;
    int *nbrnch_of_root        = nbr_list->brnch_root.nbrnch_of_root;
    int nbrnch_tot             = nbr_list->brnch_root.nbrnch_tot;
    int nroot_tot              = nbr_list->brnch_root.nroot_tot;

    int int_res_ter            = timeinfo->int_res_ter;
    int iver_count             = nbr_list->verlist.iver_count;
    int iver_fill              = nbr_list->verlist.iver_fill;
    int natm_tot               = clatoms_info->natm_tot;  
    int nlen                   = for_scr->nlen;
    int nmem_min_lst           = nbr_list->verlist.nmem_min_lst;
    int iver_init              = nbr_list->verlist.iver_init;
    double mem_safe            = nbr_list->verlist.mem_safe;
    int brnch_root_list_opt   = nbr_list->brnch_root_list_opt; 
    int np_forc                = class_comm_forc_pkg->num_proc;
    int myid_forc              = class_comm_forc_pkg->myid;
    MPI_Comm comm_forc         = class_comm_forc_pkg->comm;
    int nolst_ver_update       = nbr_list->verlist.nolst_ver_update;  

/*==========================================================================*/
/* I) Count the interactions                                               */

   intact_save      = (((nroot_tot)*(nroot_tot-1))/2-(excl_nlst));
   num_call_save = intact_save/nlen;
   if( (intact_save % nlen)!=0 ){num_call_save += 1;}
  
/*==========================================================================*/
/* II) Initialize                                                           */

   nver_now         = 0;
   nver_res_now     = 0;
   (*for_scr_intact)   = 0;
   (*for_scr_num_brk)  = 0;
   num_call         = 0;
   intact_tot       = 0;
   iswitch          = 0;
   lst_typ          = 1;

/* Zero interactions of each particle   */
   for(jpart=1;jpart<=natm_tot;jpart++){
      nbr_list_nter[jpart]     = 0;
   }/*endfor*/
   if(int_res_ter==1){
     for(jpart=1;jpart<=natm_tot;jpart++){
       nbr_list_nter_res[jpart]     = 0;
     }/*endfor*/
   }/*endif*/

/* Zero particle offsets into interaction list */
   if(iver_count==1){
     for(jpart=1;jpart<=natm_tot;jpart++){
       nbr_list_jver_off[jpart] = 0;
     }/*endfor*/
     if(int_res_ter==1){
       for(jpart=1;jpart<=natm_tot;jpart++){
         nbr_list_jver_off_res[jpart] = 0;
       }/*endfor*/
    }/*endif*/    
   }/*endif*/    

/*==========================================================================*/
/* III) Loop over the possible root interactions:                           */

   ipart_st = 1+myid_forc;
   nskip    = np_forc;
   for(ipart=ipart_st;ipart<=nroot_tot;ipart+=nskip){

/*==========================================================================*/
/* A) Check for exclusions and create an interaction list                   */

      if(excl_num[ipart]!=0){ 
        lower = excl_j_off[ipart]+1;
        upper = excl_num[ipart]+excl_j_off[ipart];
        jstart = 1;
        intact_now = 0;
        for(k=lower;k<=upper;k++){
          jend = excl_j[k]-1;
          joff = intact_now-jstart+1;
          for(jpart=jstart;jpart<=jend;jpart++){
            for_scr_iexcl[(jpart+joff)]  = root_atm_list[jpart];          
          }/*endfor*/
          intact_now += MAX((jend-jstart+1),0);
          jstart = excl_j[k]+1;
        }/*endfor*/
          jend = ipart-1;
          joff = intact_now-jstart+1;
          for(jpart=jstart;jpart<=jend;jpart++){
            for_scr_iexcl[(jpart+joff)]  = root_atm_list[jpart];          
          }        
          intact_now += MAX((jend-jstart+1),0);
      }else{
        for(jpart=1;jpart<=ipart-1;jpart++){
           for_scr_iexcl[jpart] = root_atm_list[jpart];
        }/*endfor*/
        intact_now = (ipart-1);
      }/*endif*/ 

/*==========================================================================*/
/*  B) Add interactions to the force routine interaction list               */
/*     in chunks no greater than nlen. Keep track of the break points       */

      lower  = 1;
      upper  = 0;
      while(upper!=intact_now){
        upper = intact_now;
        if(upper-lower+1+(*for_scr_intact)>(nlen)){
            upper = nlen - (*for_scr_intact) + lower-1;
         }/*endif*/
        jind_off = (*for_scr_intact)-lower+1;
        for(jpart=lower;jpart<=upper;jpart++){
           mtemp = jpart+jind_off;
           for_scr_j_index[mtemp]=for_scr_iexcl[jpart];
           for_scr_i_index[mtemp]=root_atm_list[ipart];
        }/*endfor*/
        (*for_scr_intact) += upper-lower+1;
        (*for_scr_num_brk)++;
        for_scr_num_brk_i[(*for_scr_num_brk)] = upper-lower+1;
        lower = upper+1;
/*==========================================================================*/
/*  C) If enough interactions have accumulated (or you are done)            */
/*       make a verlist create call and then update and reinitialize        */
/*       appropriate counters                                               */

/*------------------------------------------------------------------------*/
/*   i) Force call conditions                                             */
        ifor_call = 0;
        if(((*for_scr_intact)==(nlen))){ifor_call=1;}
        if((ipart==(natm_tot))&&(upper==intact_now)&&
          ((*for_scr_intact)!=0)   ){
            ifor_call=1;
        }/*endif*/

/*------------------------------------------------------------------------*/
/*  ii) List add call, add interaction, zero counters, check memory       */
        if(ifor_call==1){

          verlist_create(clatoms_info,clatoms_pos,
                          intra_scr,for_scr,atommaps,cell,
                          nbr_list,timeinfo,interact,iswitch,
                          &nver_now,&nver_res_now,lst_typ);
          num_call        += 1;
          intact_tot      += (*for_scr_intact);
          (*for_scr_intact)  = 0;
          (*for_scr_num_brk) = 0;

          ver_list_mem_check(nver_now,&(nbr_list->verlist.nver_lst),
                             iver_fill,nmem_min_lst,iver_init,mem_safe,
                             &irealloc,myid_forc);
          if(irealloc==1){
            nbr_list->verlist.jter = (list_int *)
                  crealloc(&(nbr_list->verlist.jter)[1],
                            (nbr_list->verlist.nver_lst)*sizeof(list_int))-1;
          }/*endif*/

          if(int_res_ter==1){
            ver_list_mem_check(nver_res_now,&(nbr_list->verlist.nver_lst_res),
                               iver_fill,nmem_min_lst,iver_init,mem_safe,
                               &irealloc,myid_forc);
            if(irealloc==1){
              nbr_list->verlist.jter_res = (list_int *)
                    crealloc(&(nbr_list->verlist.jter_res)[1],
                         (nbr_list->verlist.nver_lst_res)*sizeof(list_int))-1;
            }/*endif*/
          }/*endif:RESPA*/

        }/*endif:force_call*/

      }/*endwhile:accumulating interaction sets*/

   }/*endfor:each particles interactions       */

/*==========================================================================*/
/* IV) Final Force call and memory checks                                  */

  if((*for_scr_intact)>0){

     verlist_create(clatoms_info,clatoms_pos,
                    intra_scr,for_scr,atommaps,cell,
                    nbr_list,timeinfo,interact,iswitch,
                    &nver_now,&nver_res_now,lst_typ);
     num_call          += 1;
     intact_tot        += (*for_scr_intact);

     ver_list_mem_check(nver_now,&(nbr_list->verlist.nver_lst),
                        iver_fill,nmem_min_lst,iver_init,mem_safe,
                        &irealloc,myid_forc);
     if(irealloc==1){
        nbr_list->verlist.jter = (list_int *)
                  crealloc(&(nbr_list->verlist.jter)[1],
                       (nbr_list->verlist.nver_lst)*sizeof(list_int))-1;
     }/*endif*/

     if(int_res_ter==1){
       ver_list_mem_check(nver_res_now,&(nbr_list->verlist.nver_lst_res),
                          iver_fill,nmem_min_lst,iver_init,mem_safe,
                          &irealloc,myid_forc);
       if(irealloc==1){
         nbr_list->verlist.jter_res = (list_int *)
                 crealloc(&(nbr_list->verlist.jter_res)[1],
                       (nbr_list->verlist.nver_lst_res)*sizeof(list_int))-1;
       }/*endif*/
     }/*endif:RESPA*/

  }/*endif: last force call*/

/*==========================================================================*/
/* VI) Count only list expansion         */

  if(iver_fill==0 && iver_count ==1){
    expand_brnch_root_cnt(&nver_now,nbr_list_nter,nbr_list_jver_off,
                        nbrnch_tot,brnch_atm_root,brnch_atm_list,
                        nroot_tot,root_atm_list,natm_tot,natm_add,
                        ibrnch_of_root,nbrnch_of_root,
                        class_comm_forc_pkg,nolst_ver_update);
    if(int_res_ter==1){
     expand_brnch_root_cnt(&nver_res_now,nbr_list_nter_res,
                           nbr_list_jver_off_res,
                           nbrnch_tot,brnch_atm_root,brnch_atm_list,
                           nroot_tot,root_atm_list,natm_tot,natm_add,
                           ibrnch_of_root,nbrnch_of_root,
                           class_comm_forc_pkg,nolst_ver_update);
    }/*endif*/

  }/*endif*/


/*==========================================================================*/
/* VI) Count and fill list expansion    */

  if(iver_fill==1 && iver_count ==1){

    ires_flag = 0;
    expand_brnch_root_cntfll(nbr_list,ires_flag,&nver_now,
                             nbrnch_tot,brnch_atm_root,brnch_atm_list,
                             brnch_atm_map,
                             nroot_tot,root_atm_list,nbrnch_of_root,
                             ibrnch_of_root,natm_tot,natm_add,
                             iatm_add_off,iatm_add,root_atm_map,
                             iver_fill,
                             nmem_min_lst,iver_init,mem_safe,
                             for_scr_iexcl,
                             class_comm_forc_pkg);
    if(int_res_ter==1){
      ires_flag = 1;
      expand_brnch_root_cntfll(nbr_list,ires_flag,&nver_res_now,
                               nbrnch_tot,brnch_atm_root,brnch_atm_list,
                               brnch_atm_map,
                               nroot_tot,root_atm_list,nbrnch_of_root,
                               ibrnch_of_root,natm_tot,natm_add,
                               iatm_add_off,iatm_add,root_atm_map,
                               iver_fill,
                               nmem_min_lst,iver_init,mem_safe,
                               for_scr_iexcl,
                               class_comm_forc_pkg);

    }/*endif:RESPA*/

  }/*endif*/

/*==========================================================================*/
/* V) Check                                                                */

 if(np_forc==1){
    if(intact_tot!=intact_save){
       printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
       printf("Internal Error:\n ");
       printf("Incorrect number of interactions calculated\n ");
       printf("   in routine nolst_ver_gen_root \n ");
       printf("%d vs %d \n",intact_tot,intact_save);
       printf("Contact product support \n ");
       printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
       fflush(stdout);
       exit(1);
   }/*endif*/
   if(num_call!=num_call_save){
       printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
       printf("Internal Error:\n ");
       printf("Incorrect number of force calls\n ");
       printf("   in routine nolst_ver_gen_root \n ");
       printf("%d vs %d \n",num_call,num_call_save);
       printf("Contact product support \n ");
       printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
       fflush(stdout);
       exit(1);
   }/*endif*/
 }
/*==========================================================================*/
/* VI) Pare down the list if filling                                        */
   
   if(iver_fill==1 && brnch_root_list_opt==2){
      ver_pare_down_root(clatoms_info,clatoms_pos,for_scr,atommaps,
                        cell,interact,timeinfo,nbr_list,excl,
                        intra_scr,(&nver_now),&(nver_res_now),myid_forc);
   }/*endif*/

/*==========================================================================*/
/* V) Assign present size of list to appropriate variable                   */

   (*nver)     = nver_now;
   (*nver_res) = nver_res_now;

/*--------------------------------------------------------------------------*/
    }/*end routine */
/*==========================================================================*/




/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void lnk_ver_gen(CLATOMS_INFO *clatoms_info,CLATOMS_POS *clatoms_pos,
                 FOR_SCR *for_scr,ATOMMAPS *atommaps,
                 CELL *cell,INTERACT *interact,
                 TIMEINFO *timeinfo,NBR_LIST *nbr_list,
                 EXCL *excl,INTRA_SCR *intra_scr, int *nver, int *nver_res,
                 int error_check_on,
                 CLASS_COMM_FORC_PKG *class_comm_forc_pkg)

/*==========================================================================*/
/*                      Brief description                                   */

/* lnk_list[ia][ib][ic][iatm] stores the indices of iatmth atm in           */
/*                           a cell labled ia,ib,ic                         */
/*                           (a pipette cut out of the full system cell).   */
/*                           It is packed with zeros for                    */
/*                           iatm greater than the number of atoms          */
/*                           in a given cell for computational              */
/*                           efficiency. Zero index atoms give              */
/*                           rise to ``null'' interactions that are         */
/*                           removed with logic                             */
/* ishft_a[i],ishft_b[i],ishft_c[i]                                      */
/*                           are vectors containing the allowed sets of     */
/*                           shifts along a,b,c                             */
/* iexl_check[i]             indicates whether or not exclusions needed     */
/*                           to be removed in a particular shift set        */
/* When the lnk_list is shifted interactions are generated between cells    */
/* of the shifted list and old unshifted list.                              */
/* These interactions are sent to the force routine in bundles of nlen.     */
/* Break points where memory conflicts can occur when summing the forces    */
/* are tablulated.                                                          */

/*==========================================================================*/
{/*Begin Routine*/
/*==========================================================================*/
/*            Local Variable declarations                                   */

   int itemp,nlen_use,i1;              /* value of nlen to be used           */
   int ncell_tot, nlist_tot;        /* tot number of cells, size of list  */
   int natm_shft;                   /* number of atom shifts              */
   int iadd_shft;                   /* atm shift flag                     */
   int ishft, iatm_shft;           /* cell and atom shift loops indicies  */
   int ja,jb,jc,jatm;               /* the present shift                  */
   int icheck;                      /* excl check on                      */
   int intact_now, intact_tmp;      /* interaction temporaries            */
   int intact_tot, num_call;        /* tot interactions,num of force calls*/
   int i,j,k;                           /* loop counter                       */
   int imin,jmax;                   /* min/max atm ind in a interact pair */
   int nact;                        /* interaction is good flag           */
   int iact;                        /* excl interaction loop              */
   int lower,upper;                 /* limits on excl interaction loop    */
   int ifor_call;                   /* force call flag                    */
   int jind_off;                    /* integer offset                     */
   int iswitch;
   int jpart;
   int nver_now,nver_res_now,irealloc;
   int iii;
   int intact_init,intact_half,nlist_tot2,mtemp;
   int ntemp,ioff;
   int lst_typ;
   int ishft_st,nskip;

/* Local Pointer */
    int ncell_a   = nbr_list->lnklist.ncell_a;
    int ncell_b   = nbr_list->lnklist.ncell_b;
    int ncell_c   = nbr_list->lnklist.ncell_c;
    int natm_cell_max = nbr_list->lnklist.natm_cell_max;
    int nshft_lnk  = nbr_list->lnklist.nshft_lnk;

    int excl_nlst              = excl->nlst;
    int *for_scr_intact        = &(for_scr->intact);
    int *for_scr_num_brk       = &(for_scr->num_brk);
    int *for_scr_num_brk_i     = for_scr->num_brk_i;
    int *for_scr_i_index       = for_scr->i_index;  
    int *for_scr_j_index       = for_scr->j_index;  
    int *excl_j_off            = excl->j_off; 
    int *excl_num              = excl->num;
    int *excl_j                = excl->j;
    list_int *for_scr_i_lnk         = for_scr->i_lnk;
    list_int *for_scr_j_lnk         = for_scr->j_lnk;
    int *nbr_list_nter         = nbr_list->verlist.nter;
    int *nbr_list_jver_off     = nbr_list->verlist.jver_off;
    int *nbr_list_nter_res     = nbr_list->verlist.nter_res;
    int *nbr_list_jver_off_res = nbr_list->verlist.jver_off_res;
    list_int *nbr_list_lnk_list     = nbr_list->lnklist.lnk_list;
    int *nbr_list_ishft_a      = nbr_list->lnklist.ishft_a;
    int *nbr_list_ishft_b      = nbr_list->lnklist.ishft_b;
    int *nbr_list_ishft_c      = nbr_list->lnklist.ishft_c;
    int *nbr_list_iexl_chk     = nbr_list->lnklist.iexl_chk;
    double *nbr_list_shft_wght = nbr_list->lnklist.shft_wght;
    double *for_scr_wght_lnk   = &(for_scr->wght_lnk);
    int *for_scr_iexcl         = for_scr->iexcl;

    int int_res_ter            = timeinfo->int_res_ter;
    int iver_count             = nbr_list->verlist.iver_count;
    int iver_fill              = nbr_list->verlist.iver_fill;
    int natm_tot               = clatoms_info->natm_tot;  
    int nlen                   = for_scr->nlen;
    int nmem_min_lst           = nbr_list->verlist.nmem_min_lst;
    int iver_init              = nbr_list->verlist.iver_init;
    double mem_safe            = nbr_list->verlist.mem_safe;
    int np_forc                = class_comm_forc_pkg->num_proc;
    int myid_forc              = class_comm_forc_pkg->myid;

/*==========================================================================*/
/* Error check */

  if(iver_count==1 && iver_fill==1){
    printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
    printf("Internal error: No ver_count plus ver_fill with lnk lst \n");    
    printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
    fflush(stdout);
    exit(1);
  }/*endif*/

/*==========================================================================*/
/* I) Get a nice value of nlen: a value divisible by natm_tot or vice versa */

  if(nlen>natm_tot){
    nlen_use =  nlen -  (nlen % natm_tot) ;
  }else{
    itemp    = natm_tot/nlen + 1;
    nlen_use = natm_tot/itemp;
  }/*endif*/

/*==========================================================================*/
/* II) Get some useful constants*/

  ncell_tot = ncell_a*ncell_b*ncell_c;
  nlist_tot = ncell_tot*natm_cell_max;

/*==========================================================================*/
/* III) Inititialize           */
 
  nver_now           = 0;
  nver_res_now       = 0;
  (*for_scr_intact)  = 0;   
  (*for_scr_num_brk) = 0;
  num_call           = 0;
  intact_tot         = 0;
  iswitch            = 1;
  lst_typ            = 3;

/*-------------------------*/
/* zero the off set vector */
  if(iver_count==1){

     for(jpart=1;jpart<=natm_tot;jpart++){
       nbr_list_jver_off[jpart] = 0;
     }/*endfor*/

     if(int_res_ter==1){
       for(jpart=1;jpart<=natm_tot;jpart++){
         nbr_list_jver_off_res[jpart] = 0;
       }/*endfor*/
     }/*endif*/    

  }/*endif*/    

/*---------------------------------------------------*/
/* set the interactions/particle equal to the offset */
  for(jpart=1;jpart<=natm_tot;jpart++){
     nbr_list_nter[jpart]  =  nbr_list_jver_off[jpart];
  }/*endfor*/

  if(int_res_ter==1){
     for(jpart=1;jpart<=natm_tot;jpart++){
       nbr_list_nter_res[jpart]  = nbr_list_jver_off_res[jpart];
     }/*endfor*/
  }/*endif*/

/*--------------------------------------------------*/
/* Find half a shift                                */
  intact_init = natm_tot;
  intact_half = 0; 
  nlist_tot2  = 4*natm_cell_max*ncell_a*ncell_b*ncell_c
              - 4*ncell_a*ncell_b*ncell_c
              - 2*ncell_a*ncell_b
              - ncell_a;
  for(i=1;i<=natm_tot;i++){
    if(for_scr_iexcl[i]<=nlist_tot2){intact_half=i;}
  }/*endfor*/

/*==========================================================================*/
/* IV) Loop over the allowed cell shifts: ja,jb,jc                          */

   ishft_st = 1 + myid_forc;
   nskip    = np_forc;
   for(ishft=ishft_st;ishft<=nshft_lnk;ishft+=nskip){
      ja                  = nbr_list_ishft_a[ishft];
      jb                  = nbr_list_ishft_b[ishft];
      jc                  = nbr_list_ishft_c[ishft];
      icheck              = nbr_list_iexl_chk[ishft];
      (*for_scr_wght_lnk) = nbr_list_shft_wght[ishft];

/*==========================================================================*/
/* V) Loop over all the atm shifts: jatm                                    */
/*     Bookeeping for the first cell shift to take care                     */
/*           of self interactions                                           */

      natm_shft = natm_cell_max;
      iadd_shft = 0;
      if(ishft==1){
         iadd_shft  = ((natm_cell_max-1) % 2);
         natm_shft  = (natm_cell_max-1)/2 + iadd_shft;
      }/*endif*/
      for(iatm_shft=1;iatm_shft<=natm_shft;iatm_shft++){
         jatm = iatm_shft;

/*==========================================================================*/
/*  B) Book keeping for last atm shift of 1st cell shift:                  */
/*         On the last iatm_shft of the 1st cell shift get rid of           */
/*         repeated interactions that occur if natm_cell_max is even.       */
/*         This simply means taking 1/2 of the interactions found           */

         intact_now = intact_init;
         if((ishft==1)&&(iadd_shft==1)&&(iatm_shft==natm_shft))
                 {intact_now = intact_half;}

/*==========================================================================*/
/*   C) Eliminate null interactions caused by packing the lnk_list          */
/*       with zero index atoms                                              */

         intact_tmp = intact_now;
         intact_now = 0;
         ioff =  ja + jb*ncell_a*2 + jc*ncell_a*ncell_b*4 
               + jatm*ncell_a*ncell_b*ncell_c*8;
         for(i=1;i<=intact_tmp;i++){
           ntemp = for_scr_iexcl[i]+ioff;
           if(nbr_list_lnk_list[ntemp]!=0){
               intact_now += 1;
               mtemp = for_scr_iexcl[i];
               for_scr_i_lnk[intact_now] = nbr_list_lnk_list[mtemp];
               for_scr_j_lnk[intact_now] = nbr_list_lnk_list[ntemp];
           }/*endif*/
         }/*endfor*/
/*==========================================================================*/
/*   E) Eliminate standard exclusions in shifts that might have them.       */
/*       Cells that are far enough away in space won't have any excls.      */

         if((icheck==1)&&(excl_nlst>0)){
            intact_tmp = intact_now;
            intact_now  = 0;
            for(i=1;i<=intact_tmp;i++){
               imin  = MIN( (for_scr_i_lnk)[i],(for_scr_j_lnk)[i] );
               jmax  = MAX( (for_scr_i_lnk)[i],(for_scr_j_lnk)[i] );
               lower = (excl_j_off)[jmax]+1;
               upper = (excl_num)[jmax]+(excl_j_off)[jmax];
               nact  = 1;
               for(iact=lower;iact<=upper;iact++){
                 if(imin==(excl_j)[iact]){nact = 0;}
               }/*endfor*/
               if(nact==1){
                 intact_now += 1;
                 (for_scr_i_lnk)[intact_now] = (for_scr_i_lnk)[i];
                 (for_scr_j_lnk)[intact_now] = (for_scr_j_lnk)[i];
               }/*endif*/
            }/*endfor*/
         }/*endif*/

/*==========================================================================*/
/*  G) Add interactions to the temporary interaction list                   */
/*     in chunks no greater than nlen_use: Keep track of the break points   */

         lower = 1;
         upper = 0;
         while(upper!=intact_now){
           upper = intact_now;
           if(upper-lower+1+(*for_scr_intact)>nlen_use){
               upper = nlen_use-(*for_scr_intact)+lower-1;
           }/*endif*/
           jind_off = (*for_scr_intact)-lower+1;
           for(i=lower;i<=upper;i++){
              mtemp =  (i+jind_off);
              (for_scr_j_index)[mtemp]=(for_scr_j_lnk)[i];
              (for_scr_i_index)[mtemp]=(for_scr_i_lnk)[i];
           }/*endfor*/
           (*for_scr_intact)                       += upper-lower+1;
           (*for_scr_num_brk)                      += 1;
           (for_scr_num_brk_i)[(*for_scr_num_brk)] = upper-lower+1;
           lower                                    = upper+1;

/*==========================================================================*/
/*  H) If enough interactions have accumulated (or you are done)            */
/*       make a verlist create call and then update and reinitialize        */
/*       appropriate counters. Realloc if the tot num of interactions       */
/*       is too large                                                       */

           ifor_call = 0;
           if( (*for_scr_intact)>=nlen_use){ifor_call=1;}
           if(  (ishft==(nshft_lnk))&&
                 (iatm_shft==natm_shft)&&
                 (upper==intact_now)&&
                 ((*for_scr_intact)!=0) )  {ifor_call=1;}

           if(ifor_call==1){
             verlist_create(clatoms_info,clatoms_pos,
                            intra_scr,for_scr,atommaps,cell,
                            nbr_list,timeinfo,interact,iswitch,
                            &nver_now,&nver_res_now,lst_typ);
             num_call          += 1;
             intact_tot        += (*for_scr_intact);
             (*for_scr_intact)  = 0;   
             (*for_scr_num_brk) = 0;

             ver_list_mem_check(nver_now,&(nbr_list->verlist.nver_lst),
                                iver_fill,nmem_min_lst,iver_init,mem_safe,
                                &irealloc,myid_forc);
             if(irealloc==1){
               nbr_list->verlist.jter = (list_int *)
                  crealloc(&(nbr_list->verlist.jter)[1],
                            (nbr_list->verlist.nver_lst)*sizeof(list_int))-1;
             }/*endif*/

             if(int_res_ter==1){
               ver_list_mem_check(nver_res_now,
                                  &(nbr_list->verlist.nver_lst_res),
                                  iver_fill,nmem_min_lst,iver_init,mem_safe,
                                  &irealloc,myid_forc);
               if(irealloc==1){
                nbr_list->verlist.jter_res = (list_int *)
                    crealloc(&(nbr_list->verlist.jter_res)[1],
                         (nbr_list->verlist.nver_lst_res)*sizeof(list_int))-1;
               }/*endif*/
             }/*endif:RESPA*/


           }/*endif:list create call*/         

         }/*endwhile: accumulating interactions*/

/*==========================================================================*/

      }/*endfor:iatm_shft*/
  }/*endfor:ishft*/

/*==========================================================================*/
/* VI) Final list create call */

  if((*for_scr_intact)>0){

     verlist_create(clatoms_info,clatoms_pos,
                    intra_scr,for_scr,atommaps,cell,
                    nbr_list,timeinfo,interact,iswitch,
                    &nver_now,&nver_res_now,lst_typ);
     num_call          += 1;
     intact_tot        += (*for_scr_intact);

     ver_list_mem_check(nver_now,&(nbr_list->verlist.nver_lst),
                        iver_fill,nmem_min_lst,iver_init,mem_safe,
                        &irealloc,myid_forc);
     if(irealloc==1){
        nbr_list->verlist.jter = (list_int *)
               crealloc(&(nbr_list->verlist.jter)[1],
                      (nbr_list->verlist.nver_lst)*sizeof(list_int))-1;
     }/*endif*/

     if(int_res_ter==1){
        ver_list_mem_check(nver_res_now,&(nbr_list->verlist.nver_lst_res),
                           iver_fill,nmem_min_lst,iver_init,mem_safe,
                           &irealloc,myid_forc);
        if(irealloc==1){
           nbr_list->verlist.jter_res = (list_int *)
                 crealloc(&(nbr_list->verlist.jter_res)[1],
                      (nbr_list->verlist.nver_lst_res)*sizeof(list_int))-1;
        }/*endif*/
     }/*endif:RESPA*/

  }/*endif:final list create call*/

/*==========================================================================*/
/* VII) Assign present size of list to appropriate variable                 */

   *nver = nver_now;
   *nver_res = nver_res_now;

/*==========================================================================*/
/* VIII) If you are just filling depad the interactions                     */

   if((iver_count==0)&&(iver_fill==1)){

     for(jpart=2;jpart<=natm_tot;jpart++){
       nbr_list_nter[jpart]  -=  (nbr_list_jver_off[jpart]);
     }/*endfor*/

     if(int_res_ter==1){
       for(jpart=2;jpart<=natm_tot;jpart++){
         nbr_list_nter_res[jpart]  -= nbr_list_jver_off_res[jpart];
       }/*endfor*/
     }/*endif*/

  }/*endif*/

/*==========================================================================*/
/* IX) If you are just counting fill the offset list                      */

   if((iver_count==1)&&(iver_fill==0)){

     nbr_list_jver_off[1] = 0;
     nbr_list_nter[1]     = 0;
     for(jpart=2;jpart<=natm_tot;jpart++){
        mtemp = jpart-1; 
        nbr_list_jver_off[jpart] = nbr_list_jver_off[mtemp]
                                 + nbr_list_nter[mtemp];
     }/*endfor*/

     if(int_res_ter==1){
        nbr_list_jver_off_res[1] = 0;
        nbr_list_nter_res[1]     = 0;
        for(jpart=2;jpart<=natm_tot;jpart++){
           mtemp = jpart-1; 
           nbr_list_jver_off_res[jpart] = nbr_list_jver_off_res[mtemp]
                                        + nbr_list_nter_res[mtemp];
        }/*endfor*/
     }/*endif*/




   }/*endif*/

/*--------------------------------------------------------------------------*/
         }/*end routine */
/*==========================================================================*/






/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void lnk_ver_gen_root(CLATOMS_INFO *clatoms_info,CLATOMS_POS *clatoms_pos,
                      FOR_SCR *for_scr,ATOMMAPS *atommaps,
                      CELL *cell,INTERACT *interact,
                      TIMEINFO *timeinfo,NBR_LIST *nbr_list,
                      EXCL *excl,INTRA_SCR *intra_scr, int *nver, 
                      int *nver_res,int error_check_on,
                      CLASS_COMM_FORC_PKG *class_comm_forc_pkg)

/*==========================================================================*/
/*                      Brief description                                   */

/* lnk_list[ia][ib][ic][iatm] stores the indices of iatmth atm in           */
/*                           a cell labled ia,ib,ic                         */
/*                           (a pipette cut out of the full system cell).   */
/*                           It is packed with zeros for                    */
/*                           iatm greater than the number of atoms          */
/*                           in a given cell for computational              */
/*                           efficiency. Zero index atoms give              */
/*                           rise to ``null'' interactions that are         */
/*                           removed with logic                             */
/* ishft_a[i],ishft_b[i],ishft_c[i]                                      */
/*                           are vectors containing the allowed sets of     */
/*                           shifts along a,b,c                             */
/* iexl_check[i]             indicates whether or not exclusions needed     */
/*                           to be removed in a particular shift set        */
/* When the lnk_list is shifted interactions are generated between cells    */
/* of the shifted list and old unshifted list.                              */
/* These interactions are sent to the force routine in bundles of nlen.     */
/* Break points where memory conflicts can occur when summing the forces    */
/* are tablulated.                                                          */

/*==========================================================================*/
{/*Begin Routine*/
/*==========================================================================*/
/*            Local Variable declarations                                   */

   int itemp,jtemp,nlen_use,i1;     /* value of nlen to be used           */
   int ncell_tot, nlist_tot;        /* tot number of cells, size of list  */
   int natm_shft;                   /* number of atom shifts              */
   int iadd_shft;                   /* atm shift flag                     */
   int ishft, iatm_shft;           /* cell and atom shift loops indicies  */
   int ja,jb,jc,jatm;               /* the present shift                  */
   int icheck;                      /* excl check on                      */
   int intact_now, intact_tmp;      /* interaction temporaries            */
   int intact_tot, num_call;        /* tot interactions,num of force calls*/
   int i,k,ires_flag;                 /* loop counter                       */
   int imin,jmax;                   /* min/max atm ind in a interact pair */
   int nact;                        /* interaction is good flag           */
   int iact;                        /* excl interaction loop              */
   int lower,upper;                 /* limits on excl interaction loop    */
   int ifor_call;                   /* force call flag                    */
   int jind_off;                    /* integer offset                     */
   int iswitch;
   int jpart;
   int nver_now,nver_res_now;
   int iii;
   int intact_init,intact_half,nlist_tot2,mtemp,irealloc;
   int ntemp,ioff;
   int lst_typ;
   int my_ind,jind,ind,ngo,my_root,j;
   int ishft_st,nskip;

/* Local Pointer */
    int ncell_a   = nbr_list->lnklist.ncell_a;
    int ncell_b   = nbr_list->lnklist.ncell_b;
    int ncell_c   = nbr_list->lnklist.ncell_c;
    int natm_cell_max = nbr_list->lnklist.natm_cell_max;
    int nshft_lnk  = nbr_list->lnklist.nshft_lnk;

    int excl_nlst              = excl->nlst_root;
    int *for_scr_intact        = &(for_scr->intact);
    int *for_scr_num_brk       = &(for_scr->num_brk);
    int *for_scr_num_brk_i     = for_scr->num_brk_i;
    int *for_scr_i_index       = for_scr->i_index;  
    int *for_scr_j_index       = for_scr->j_index;  
    int *excl_j_off            = excl->j_off_root; 
    int *excl_num              = excl->num_root;
    int *excl_j                = excl->j_root;
    list_int *for_scr_i_lnk         = for_scr->i_lnk;
    list_int *for_scr_j_lnk         = for_scr->j_lnk;
    int *nbr_list_nter         = nbr_list->verlist.nter;
    int *nbr_list_jver_off     = nbr_list->verlist.jver_off;
    int *nbr_list_nter_res     = nbr_list->verlist.nter_res;
    int *nbr_list_jver_off_res = nbr_list->verlist.jver_off_res;
    list_int *nbr_list_lnk_list     = nbr_list->lnklist.lnk_list;
    int *nbr_list_ishft_a      = nbr_list->lnklist.ishft_a;
    int *nbr_list_ishft_b      = nbr_list->lnklist.ishft_b;
    int *nbr_list_ishft_c      = nbr_list->lnklist.ishft_c;
    int *nbr_list_iexl_chk     = nbr_list->lnklist.iexl_chk;
    double *nbr_list_shft_wght = nbr_list->lnklist.shft_wght;
    double *for_scr_wght_lnk   = &(for_scr->wght_lnk);
    int *for_scr_iexcl         = for_scr->iexcl;
    list_int *nbr_list_jter;        /* assigned below in case of reallocs*/
    list_int *nbr_list_jter_res;    /* assigned below in case of reallocs*/

    int *brnch_atm_root        = nbr_list->brnch_root.brnch_atm_root;
    int *brnch_atm_list        = nbr_list->brnch_root.brnch_atm_list;
    int *natm_add              = nbr_list->brnch_root.natm_add;
    int *iatm_add              = nbr_list->brnch_root.iatm_add;
    int *iatm_add_off          = nbr_list->brnch_root.iatm_add_off;
    int *root_atm_list         = nbr_list->brnch_root.root_atm_list;
    int *root_atm_map          = nbr_list->brnch_root.root_atm_map;
    int *brnch_atm_map         = nbr_list->brnch_root.brnch_atm_map;
    int **ibrnch_of_root       = nbr_list->brnch_root.ibrnch_of_root;
    int *nbrnch_of_root        = nbr_list->brnch_root.nbrnch_of_root;
    int nbrnch_tot             = nbr_list->brnch_root.nbrnch_tot;
    int nroot_tot              = nbr_list->brnch_root.nroot_tot;

    int int_res_ter            = timeinfo->int_res_ter;
    int iver_count             = nbr_list->verlist.iver_count;
    int iver_fill              = nbr_list->verlist.iver_fill;
    int natm_tot               = clatoms_info->natm_tot;  
    int nlen                   = for_scr->nlen;
    int nmem_min_lst           = nbr_list->verlist.nmem_min_lst;
    int iver_init              = nbr_list->verlist.iver_init;
    double mem_safe            = nbr_list->verlist.mem_safe;
    int brnch_root_list_opt   = nbr_list->brnch_root_list_opt; 
    int np_forc                = class_comm_forc_pkg->num_proc;
    int myid_forc              = class_comm_forc_pkg->myid;
    int nolst_ver_update  = nbr_list->verlist.nolst_ver_update;

/*==========================================================================*/
/*  Error check */

  if(iver_count==1 && iver_fill==1){
    printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
    printf("Internal error: No ver_count plus ver_fill with lnk lst \n");    
    printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
    fflush(stdout);
    exit(1);
  }/*endif*/

/*==========================================================================*/
/* I) Get a nice value of nlen: a value divisible by nroot_totor vice versa */

  if(nlen>nroot_tot){
    nlen_use =  nlen -  (nlen % nroot_tot) ;
  }else{
    itemp    = nroot_tot/nlen + 1;
    nlen_use = nroot_tot/itemp;
  }/*endif*/

/*==========================================================================*/
/* II) Get some useful constants*/

  ncell_tot = ncell_a*ncell_b*ncell_c;
  nlist_tot = ncell_tot*natm_cell_max;

/*==========================================================================*/
/* III) Inititialize           */
 
  nver_now           = 0;
  nver_res_now       = 0;
  (*for_scr_intact)  = 0;   
  (*for_scr_num_brk) = 0;
  num_call           = 0;
  intact_tot         = 0;
  iswitch            = 1;
  lst_typ            = 3;

/*-------------------------*/
/* zero the off set vector */
  if(iver_count==1){

     for(jpart=1;jpart<=natm_tot;jpart++){
       nbr_list_jver_off[jpart] = 0;
     }/*endfor*/

     if(int_res_ter==1){
       for(jpart=1;jpart<=natm_tot;jpart++){
         nbr_list_jver_off_res[jpart] = 0;
       }/*endfor*/
     }/*endif*/    

  }/*endif*/    

/*---------------------------------------------------*/
/* set the interactions/particle equal to the offset */
  for(jpart=1;jpart<=natm_tot;jpart++){
     nbr_list_nter[jpart]  =  nbr_list_jver_off[jpart];
  }/*endfor*/

  if(int_res_ter==1){
     for(jpart=1;jpart<=natm_tot;jpart++){
       nbr_list_nter_res[jpart]  = nbr_list_jver_off_res[jpart];
     }/*endfor*/
  }/*endif*/

/*--------------------------------------------------*/
/* Find half a shift                                */
  intact_init = nroot_tot;
  intact_half = 0; 
  nlist_tot2  = 4*natm_cell_max*ncell_a*ncell_b*ncell_c
              - 4*ncell_a*ncell_b*ncell_c
              - 2*ncell_a*ncell_b
              - ncell_a;
  for(i=1;i<=nroot_tot;i++){
    if(for_scr_iexcl[i]<=nlist_tot2){intact_half=i;}
  }/*endfor*/

/*==========================================================================*/
/* IV) Loop over the allowed cell shifts: ja,jb,jc                          */

   ishft_st = 1 + myid_forc;
   nskip    = np_forc;
   for(ishft=ishft_st;ishft<=nshft_lnk;ishft+=nskip){
      ja                  = nbr_list_ishft_a[ishft];
      jb                  = nbr_list_ishft_b[ishft];
      jc                  = nbr_list_ishft_c[ishft];
      icheck              = nbr_list_iexl_chk[ishft];
      (*for_scr_wght_lnk) = nbr_list_shft_wght[ishft];

/*==========================================================================*/
/* V) Loop over all the atm shifts: jatm                                    */
/*     Bookeeping for the first cell shift to take care                     */
/*           of self interactions                                           */

      natm_shft = natm_cell_max;
      iadd_shft = 0;
      if(ishft==1){
         iadd_shft  = ((natm_cell_max-1) % 2);
         natm_shft  = (natm_cell_max-1)/2 + iadd_shft;
      }/*endif*/
      for(iatm_shft=1;iatm_shft<=natm_shft;iatm_shft++){
         jatm = iatm_shft;

/*==========================================================================*/
/*  B) Book keeping for last atm shift of 1st cell shift:                  */
/*         On the last iatm_shft of the 1st cell shift get rid of           */
/*         repeated interactions that occur if natm_cell_max is even.       */
/*         This simply means taking 1/2 of the interactions found           */

         intact_now = intact_init;
         if((ishft==1)&&(iadd_shft==1)&&(iatm_shft==natm_shft))
                 {intact_now = intact_half;}

/*==========================================================================*/
/*   C) Eliminate null interactions caused by packing the lnk_list          */
/*       with zero index atoms                                              */

         intact_tmp = intact_now;
         intact_now = 0;
         ioff =  ja + jb*ncell_a*2 + jc*ncell_a*ncell_b*4 
               + jatm*ncell_a*ncell_b*ncell_c*8;
         for(i=1;i<=intact_tmp;i++){
           ntemp = for_scr_iexcl[i]+ioff;
           if(nbr_list_lnk_list[ntemp]!=0){
               intact_now += 1;
               mtemp = for_scr_iexcl[i];
               for_scr_i_lnk[intact_now] = nbr_list_lnk_list[mtemp];
               for_scr_j_lnk[intact_now] = nbr_list_lnk_list[ntemp];
           }/*endif*/
         }/*endfor*/

/*==========================================================================*/
/*   E) Eliminate standard exclusions in shifts that might have them.       */
/*       Cells that are far enough away in space won't have any excls.      */

         if((icheck==1)&&(excl_nlst>0)){

            intact_tmp = intact_now;
            intact_now  = 0;
            for(i=1;i<=intact_tmp;i++){
               itemp = root_atm_map[(for_scr_i_lnk)[i]];
               jtemp = root_atm_map[(for_scr_j_lnk)[i]];
               imin  = MIN(itemp,jtemp);
               jmax  = MAX(itemp,jtemp);
               lower = (excl_j_off)[jmax]+1;
               upper = (excl_num)[jmax]+(excl_j_off)[jmax];
               nact  = 1;
               for(iact=lower;iact<=upper;iact++){
                 if(imin==(excl_j)[iact]){nact = 0;}
               }/*endfor*/
               if(nact==1){
                 intact_now += 1;
                 (for_scr_i_lnk)[intact_now] = (for_scr_i_lnk)[i];
                 (for_scr_j_lnk)[intact_now] = (for_scr_j_lnk)[i];
               }/*endif*/
            }/*endfor*/

         }/*endif*/

/*==========================================================================*/
/*  G) Add interactions to the temporary interaction list                   */
/*     in chunks no greater than nlen_use: Keep track of the break points   */

         lower = 1;
         upper = 0;
         while(upper!=intact_now){
           upper = intact_now;
           if(upper-lower+1+(*for_scr_intact)>nlen_use){
               upper = nlen_use-(*for_scr_intact)+lower-1;
           }/*endif*/
           jind_off = (*for_scr_intact)-lower+1;
           for(i=lower;i<=upper;i++){
              mtemp =  (i+jind_off);
              (for_scr_j_index)[mtemp]=(for_scr_j_lnk)[i];
              (for_scr_i_index)[mtemp]=(for_scr_i_lnk)[i];
           }/*endfor*/
           (*for_scr_intact)                       += upper-lower+1;
           (*for_scr_num_brk)                      += 1;
           (for_scr_num_brk_i)[(*for_scr_num_brk)] = upper-lower+1;
           lower                                    = upper+1;

/*==========================================================================*/
/*  H) If enough interactions have accumulated (or you are done)            */
/*       make a verlist create call and then update and reinitialize        */
/*       appropriate counters. Realloc if the tot num of interactions       */
/*       is too large                                                       */


           ifor_call = 0;
           if( (*for_scr_intact)>=nlen_use){ifor_call=1;}
           if(  (ishft==(nshft_lnk))&&
                 (iatm_shft==natm_shft)&&
                 (upper==intact_now)&&
                 ((*for_scr_intact)!=0) )  {ifor_call=1;}

           if(ifor_call==1){
             verlist_create(clatoms_info,clatoms_pos,
                            intra_scr,for_scr,atommaps,cell,
                            nbr_list,timeinfo,interact,iswitch,
                            &nver_now,&nver_res_now,lst_typ);
             num_call          += 1;
             intact_tot        += (*for_scr_intact);
             (*for_scr_intact)  = 0;   
             (*for_scr_num_brk) = 0;

             ver_list_mem_check(nver_now,&(nbr_list->verlist.nver_lst),
                                iver_fill,nmem_min_lst,iver_init,mem_safe,
                                &irealloc,myid_forc);
             if(irealloc==1){
               nbr_list->verlist.jter = (list_int *)
                  crealloc(&(nbr_list->verlist.jter)[1],
                        (nbr_list->verlist.nver_lst)*sizeof(list_int))-1;
             }/*endif*/

             if(int_res_ter==1){
               ver_list_mem_check(nver_res_now,
                                  &(nbr_list->verlist.nver_lst_res),
                                  iver_fill,nmem_min_lst,iver_init,mem_safe,
                                  &irealloc,myid_forc);
               if(irealloc==1){
                nbr_list->verlist.jter_res = (list_int *)
                    crealloc(&(nbr_list->verlist.jter_res)[1],
                        (nbr_list->verlist.nver_lst_res)*sizeof(list_int))-1;
               }/*endif*/
             }/*endif:RESPA*/


           }/*endif:list create call*/         

         }/*endwhile: accumulating interactions*/

/*==========================================================================*/

      }/*endfor:iatm_shft*/
  }/*endfor:ishft*/

/*==========================================================================*/
/* VI) Final list create call */

  if((*for_scr_intact)>0){

     verlist_create(clatoms_info,clatoms_pos,
                    intra_scr,for_scr,atommaps,cell,
                    nbr_list,timeinfo,interact,iswitch,
                    &nver_now,&nver_res_now,lst_typ);
     num_call          += 1;
     intact_tot        += (*for_scr_intact);

     ver_list_mem_check(nver_now,&(nbr_list->verlist.nver_lst),
                        iver_fill,nmem_min_lst,iver_init,mem_safe,
                        &irealloc,myid_forc);
     if(irealloc==1){
        nbr_list->verlist.jter = (list_int *)
               crealloc(&(nbr_list->verlist.jter)[1],
                      (nbr_list->verlist.nver_lst)*sizeof(list_int))-1;
     }/*endif*/

     if(int_res_ter==1){
        ver_list_mem_check(nver_res_now,&(nbr_list->verlist.nver_lst_res),
                           iver_fill,nmem_min_lst,iver_init,mem_safe,
                           &irealloc,myid_forc);
        if(irealloc==1){
           nbr_list->verlist.jter_res = (list_int *)
                 crealloc(&(nbr_list->verlist.jter_res)[1],
                     (nbr_list->verlist.nver_lst_res)*sizeof(list_int))-1;
        }/*endif*/
     }/*endif:RESPA*/

  }/*endif:final list create call*/

/*==========================================================================*/
/* VIII) Fill only                                                          */

  if(iver_fill==1 && iver_count ==0){
    ires_flag = 0;
    expand_brnch_root_fll(nbr_list,ires_flag,&nver_now,
                          nbrnch_tot,brnch_atm_root,brnch_atm_list,
                          brnch_atm_map,
                          nroot_tot,root_atm_list,nbrnch_of_root,
                          ibrnch_of_root,natm_tot,natm_add,
                          iatm_add_off,iatm_add,root_atm_map,
                          irealloc,iver_fill,
                          nmem_min_lst,iver_init,mem_safe,
                          class_comm_forc_pkg);
    if(int_res_ter==1){
      ires_flag = 1;
      expand_brnch_root_fll(nbr_list,ires_flag,&nver_res_now,
                            nbrnch_tot,brnch_atm_root,brnch_atm_list,
                            brnch_atm_map,
                            nroot_tot,root_atm_list,nbrnch_of_root,
                            ibrnch_of_root,natm_tot,natm_add,
                            iatm_add_off,iatm_add,root_atm_map,
                            irealloc,iver_fill,
                            nmem_min_lst,iver_init,mem_safe,
                            class_comm_forc_pkg);

    }/*endif:RESPA*/

  }/*endif*/

/*==========================================================================*/
/* VI) Count only list expansion         */

  if(iver_fill==0 && iver_count ==1){
     expand_brnch_root_cnt(&nver_now,nbr_list_nter,nbr_list_jver_off,
                           nbrnch_tot,brnch_atm_root,brnch_atm_list,
                           nroot_tot,root_atm_list,natm_tot,natm_add,
                           ibrnch_of_root,nbrnch_of_root,
                           class_comm_forc_pkg,nolst_ver_update);
    if(int_res_ter==1){
     expand_brnch_root_cnt(&nver_res_now,nbr_list_nter_res,
                           nbr_list_jver_off_res,
                           nbrnch_tot,brnch_atm_root,brnch_atm_list,
                           nroot_tot,root_atm_list,natm_tot,natm_add,
                           ibrnch_of_root,nbrnch_of_root,
                           class_comm_forc_pkg,nolst_ver_update);
    }/*endif*/

  }/*endif*/

/*==========================================================================*/
/* VI) Pare down the list if filling                                        */
   
   if(iver_fill==1&& brnch_root_list_opt==2){
      ver_pare_down_root(clatoms_info,clatoms_pos,for_scr,atommaps,
                        cell,interact,timeinfo,nbr_list,excl,
                        intra_scr,&nver_now,&nver_res_now,myid_forc);

   }/*endif*/




/*==========================================================================*/
/* VII) Assign present size of list to appropriate variable                 */

   *nver     = nver_now;
   *nver_res = nver_res_now;

/*--------------------------------------------------------------------------*/
     }/*end routine */
/*==========================================================================*/





/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void lnk_ver_gen_pad(CLATOMS_INFO *clatoms_info,CLATOMS_POS *clatoms_pos,
                     FOR_SCR *for_scr,ATOMMAPS *atommaps,
                     CELL *cell,INTERACT *interact,
                     TIMEINFO *timeinfo,NBR_LIST *nbr_list,
                     EXCL *excl,INTRA_SCR *intra_scr, int *nver, 
                     int *nver_res, int *isuccess,int error_check_on,
                     CLASS_COMM_FORC_PKG *class_comm_forc_pkg)

/*==========================================================================*/
/*                      Brief description                                   */

/* lnk_list[ia][ib][ic][iatm] stores the indices of iatmth atm in           */
/*                           a cell labled ia,ib,ic                         */
/*                           (a pipette cut out of the full system cell).   */
/*                           It is packed with zeros for                    */
/*                           iatm greater than the number of atoms          */
/*                           in a given cell for computational              */
/*                           efficiency. Zero index atoms give              */
/*                           rise to ``null'' interactions that are         */
/*                           removed with logic                             */
/* ishft_a[i],ishft_b[i],ishft_c[i]                                      */
/*                           are vectors containing the allowed sets of     */
/*                           shifts along a,b,c                             */
/* iexl_check[i]             indicates whether or not exclusions needed     */
/*                           to be removed in a particular shift set        */
/* When the lnk_list is shifted interactions are generated between cells    */
/* of the shifted list and old unshifted list.                              */
/* These interactions are sent to the force routine in bundles of nlen.     */
/* Break points where memory conflicts can occur when summing the forces    */
/* are tablulated.                                                          */

/*==========================================================================*/
{/*Begin Routine*/
/*==========================================================================*/
/*            Local Variable declarations                                   */

   int itemp,nlen_use;              /* value of nlen to be used           */
   int ncell_tot, nlist_tot;        /* tot number of cells, size of list  */
   int natm_shft;                   /* number of atom shifts              */
   int iadd_shft;                   /* atm shift flag                     */
   int ishft, iatm_shft;           /* cell and atom shift loops indicies  */
   int ja,jb,jc,jatm;               /* the present shift                  */
   int icheck;                      /* excl check on                      */
   int intact_now, intact_tmp;      /* interaction temporaries            */
   int intact_tot, num_call;        /* tot interactions,num of force calls*/
   int i;                           /* loop counter                       */
   int imin,jmax;                   /* min/max atm ind in a interact pair */
   int nact;                        /* interaction is good flag           */
   int iact;                        /* excl interaction loop              */
   int lower,upper;                 /* limits on excl interaction loop    */
   int ifor_call;                   /* force call flag                    */
   int jind_off;                    /* integer offset                     */
   int iswitch;
   int jpart;
   int nver_now,nver_res_now;
   int iii;
   int jver_pad_new;
   int ioff_tot,iter,irealloc;
   int jver_tmp,iver_tmp;
   int i1,i2;
   int npad_tot;
   int intact_init,intact_half,nlist_tot2,mtemp;
   int ioff,ntemp;
   int lst_typ;
   int iver_count;
   int iver_fill;
   int ishft_st,nskip;

/* Local Pointers */

   int ncell_a       = nbr_list->lnklist.ncell_a;
   int ncell_b       = nbr_list->lnklist.ncell_b;
   int ncell_c       = nbr_list->lnklist.ncell_c;
   int natm_cell_max = nbr_list->lnklist.natm_cell_max;
   int nshft_lnk     = nbr_list->lnklist.nshft_lnk;

   int excl_nlst              = excl->nlst;
   int *for_scr_intact        = &(for_scr->intact);
   int *for_scr_num_brk       = &(for_scr->num_brk);
   double *for_scr_wght_lnk   = &(for_scr->wght_lnk);
   int *for_scr_num_brk_i     = for_scr->num_brk_i;
   int *for_scr_i_index       = for_scr->i_index;  
   int *for_scr_j_index       = for_scr->j_index;  
   int *excl_j_off            = excl->j_off; 
   int *excl_num              = excl->num;
   int *excl_j                = excl->j;
   list_int *for_scr_i_lnk         = for_scr->i_lnk;
   list_int *for_scr_j_lnk         = for_scr->j_lnk;
   int *nbr_list_nter         = nbr_list->verlist.nter;
   int *nbr_list_jver_off     = nbr_list->verlist.jver_off;
   int *nbr_list_nter_res     = nbr_list->verlist.nter_res;
   int *nbr_list_jver_off_res = nbr_list->verlist.jver_off_res;
   list_int *nbr_list_lnk_list     = nbr_list->lnklist.lnk_list;
   int *nbr_list_ishft_a      = nbr_list->lnklist.ishft_a;
   int *nbr_list_ishft_b      = nbr_list->lnklist.ishft_b;
   int *nbr_list_ishft_c      = nbr_list->lnklist.ishft_c;
   int *nbr_list_iexl_chk     = nbr_list->lnklist.iexl_chk;
   double *nbr_list_shft_wght = nbr_list->lnklist.shft_wght;
   int *for_scr_iexcl         = for_scr->iexcl;
   list_int *nbr_list_jter;        /* Assigned below because of reallocs */
   list_int *nbr_list_jter_res;


   int int_res_ter            = timeinfo->int_res_ter;
   int natm_tot               = clatoms_info->natm_tot;  
   int nlen                   = for_scr->nlen;
   int nmem_min_lst           = nbr_list->verlist.nmem_min_lst;
   int iver_init              = nbr_list->verlist.iver_init;
   double mem_safe            = nbr_list->verlist.mem_safe;
   int jver_pad               = nbr_list->verlist.jver_pad;
   int nver_lst               = nbr_list->verlist.nver_lst;
   int nver_lst_res           = nbr_list->verlist.nver_lst_res;
   int np_forc                = class_comm_forc_pkg->num_proc;
   int myid_forc              = class_comm_forc_pkg->myid;

/*==========================================================================*/
/* I) Get a nice value of nlen: a value divisible by natm_tot or vice versa */

  if(nlen>natm_tot){
    nlen_use =  nlen - ( nlen % natm_tot );
  }else{
    itemp    = natm_tot/nlen + 1;
    nlen_use = natm_tot/itemp;
  }/*endif*/

/*==========================================================================*/
/* II) Get some useful constants*/

  ncell_tot = ncell_a*ncell_b*ncell_c;
  nlist_tot = ncell_tot*natm_cell_max;

  iver_fill  = 1;
  iver_count = 0;
  nbr_list->verlist.iver_fill=1;    /* This routine fills only !!!!!!! */  
  nbr_list->verlist.iver_count=0;   

/*==========================================================================*/
/* III) Inititialize scalars          */
 
  nver_now           = 0;
  nver_res_now       = 0;
  (*for_scr_intact)  = 0;   
  (*for_scr_num_brk) = 0;
  num_call           = 0;
  intact_tot         = 0;
  iswitch            = 1;
  lst_typ            = 3;

/*==========================================================================*/
/* III) Inititialize lists          */

/*----------------------------------*/
/* Pad the offset lists             */

  nbr_list_jver_off[1] = 0;
  for(jpart=2;jpart<=natm_tot;jpart++){
     nbr_list_jver_off[jpart] += jver_pad*(jpart-1);
  }/*endfor*/

  if(int_res_ter==1){
    nbr_list_jver_off_res[1] = 0;
    for(jpart=2;jpart<=natm_tot;jpart++){
      nbr_list_jver_off_res[jpart] += jver_pad*(jpart-1);
    }/*endfor*/
  }/*endif*/    
  npad_tot = jver_pad*(natm_tot-1);

/*----------------------------------------------------*/
/* Set the interactions/particle equal to the offsets */

  for(jpart=1;jpart<=natm_tot;jpart++){
     nbr_list_nter[jpart]     = nbr_list_jver_off[jpart];
  }/*endfor*/
  if( (nbr_list_nter[natm_tot]+nmem_min_lst) > nver_lst){
    printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
    printf("Neighor list incorrectly allocated\n");
    printf("Contact Technical Support \n");
    printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
    fflush(stdout);
    exit(1);
  }/*endif*/

  if(int_res_ter==1){
     for(jpart=1;jpart<=natm_tot;jpart++){
       nbr_list_nter_res[jpart] = nbr_list_jver_off_res[jpart];
     }/*endfor*/
     if((nbr_list_nter_res[natm_tot]+nmem_min_lst)>nver_lst_res){
        printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
        printf("RESPA Neighor list incorrectly allocated\n");
        printf("Contact Technical Support \n");
        printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
        fflush(stdout);
        exit(1);
     }/*endif*/
  }/*endif*/

/*------------------------------------*/
/* Find half a shift                  */

  intact_init = natm_tot;
  nlist_tot2  = 4*natm_cell_max*ncell_a*ncell_b*ncell_c
              - 4*ncell_a*ncell_b*ncell_c
              - 2*ncell_a*ncell_b
              - ncell_a;
  intact_half = 0; 
  for(i=1;i<=natm_tot;i++){
    if(for_scr_iexcl[i]<=nlist_tot2){intact_half=i;}
  }/*endif*/

/*==========================================================================*/
/* IV) Loop over the allowed cell shifts: ja,jb,jc                          */

   ishft_st = 1 + myid_forc;
   nskip    = np_forc;
   for(ishft=ishft_st;ishft<=nshft_lnk;ishft+=nskip){
      ja                  = nbr_list_ishft_a[ishft];
      jb                  = nbr_list_ishft_b[ishft];
      jc                  = nbr_list_ishft_c[ishft];
      icheck              = nbr_list_iexl_chk[ishft];
      (*for_scr_wght_lnk) = nbr_list_shft_wght[ishft];

/*==========================================================================*/
/* V) Loop over all the atm shifts: jatm                                    */
/*     Bookeeping for the first cell shift to take care                     */
/*           of self interactions                                           */

      natm_shft = natm_cell_max;
      iadd_shft = 0;
      if(ishft==1){
         iadd_shft  = ((natm_cell_max-1) % 2);
         natm_shft  = (natm_cell_max-1)/2 + iadd_shft;
      }/*endif*/
      for(iatm_shft=1;iatm_shft<=natm_shft;iatm_shft++){
         jatm = iatm_shft;

/*==========================================================================*/
/*  B) Book keeping for last atm shift of 1st cell shift:                  */
/*         On the last iatm_shft of the 1st cell shift get rid of           */
/*         repeated interactions that occur if natm_cell_max is even.       */
/*         This simply means taking 1/2 of the interactions found           */

         intact_now = intact_init;
         if((ishft==1)&&(iadd_shft==1)&&(iatm_shft==natm_shft))
                 {intact_now = intact_half;}

/*==========================================================================*/
/*   C) Eliminate null interactions caused by packing the lnk_list          */
/*       with zero index atoms                                              */

         intact_tmp = intact_now;
         intact_now = 0;
         ioff =  ja + jb*ncell_a*2 + jc*ncell_a*ncell_b*4 
               + jatm*ncell_a*ncell_b*ncell_c*8;
         for(i=1;i<=intact_tmp;i++){
           itemp = for_scr_iexcl[i];
           ntemp = itemp+ioff;
           mtemp = nbr_list_lnk_list[ntemp];
           if(mtemp!=0){
            intact_now += 1;
               for_scr_i_lnk[intact_now] = nbr_list_lnk_list[itemp];
               for_scr_j_lnk[intact_now] = mtemp;
           }/*endif*/
         }/*endfor*/

/*==========================================================================*/
/*   E) Eliminate standard exclusions in shifts that might have them.       */
/*       Cells that are far enough away in space won't have any excls.      */

         if((icheck==1)&&(excl_nlst>0)){
            intact_tmp = intact_now;
            intact_now  = 0;
            for(i=1;i<=intact_tmp;i++){
               imin  = MIN( (for_scr_i_lnk)[i],(for_scr_j_lnk)[i] );
               jmax  = MAX( (for_scr_i_lnk)[i],(for_scr_j_lnk)[i] );
               lower = (excl_j_off)[jmax]+1;
               upper = (excl_num)[jmax]+(excl_j_off)[jmax];
               nact  = 1;
               for(iact=lower;iact<=upper;iact++){
                 if(imin==(excl_j)[iact]){nact = 0;}
               }/*endfor*/
               if(nact==1){
                 intact_now += 1;
                 (for_scr_i_lnk)[intact_now] = (for_scr_i_lnk)[i];
                 (for_scr_j_lnk)[intact_now] = (for_scr_j_lnk)[i];
               }/*endif*/
            }/*endfor*/
         }/*endif*/

/*==========================================================================*/
/*  G) Add interactions to the temporary interaction list                   */
/*     in chunks no greater than nlen_use: Keep track of the break points   */

        lower = 1;
        upper = 0;
        while(upper!=intact_now){
           upper = intact_now;
           if(upper-lower+1+(*for_scr_intact)>nlen_use){
               upper = nlen_use-(*for_scr_intact)+lower-1;
           }/*endif*/
           jind_off = (*for_scr_intact)-lower+1;
           for(i=lower;i<=upper;i++){
              mtemp = (i+jind_off);
              (for_scr_j_index)[mtemp]=(for_scr_j_lnk)[i];
              (for_scr_i_index)[mtemp]=(for_scr_i_lnk)[i];
           }/*endfor*/
           (*for_scr_intact)                       += upper-lower+1;
           (*for_scr_num_brk)                      += 1;
           (for_scr_num_brk_i)[(*for_scr_num_brk)] = upper-lower+1;
           lower                                   = upper+1;

/*==========================================================================*/
/*  H) If enough interactions have accumulated (or you are done)            */
/*       make a verlist create call and then update and reinitialize        */
/*       appropriate counters. Realloc if the tot num of interactions       */
/*       is too large                                                       */

           ifor_call = 0;
           if( (*for_scr_intact)>=nlen_use){ifor_call=1;}
           if(  (ishft==(nshft_lnk))&&
                 (iatm_shft==natm_shft)&&
                 (upper==intact_now)&&
                 ((*for_scr_intact)!=0) ) {ifor_call=1;}

           if(ifor_call==1){
             verlist_create(clatoms_info,clatoms_pos,
                            intra_scr,for_scr,atommaps,cell,
                            nbr_list,timeinfo,interact,iswitch,
                            &nver_now,&nver_res_now,lst_typ);
              num_call          += 1;
              intact_tot        += (*for_scr_intact);
              (*for_scr_intact)  = 0;   
              (*for_scr_num_brk) = 0;

              verpad_list_mem_check(nver_now,&(nbr_list->verlist.nver_lst),
                                     nbr_list_nter[natm_tot],nmem_min_lst,
                                     npad_tot,&irealloc,myid_forc);
              if(irealloc==1){
                  nbr_list->verlist.jter = (list_int *)
                      crealloc(&(nbr_list->verlist.jter)[1],
                    (nbr_list->verlist.nver_lst)*sizeof(list_int))-1;
              }/*endif*/
            
              if(int_res_ter==1){
                 verpad_list_mem_check(nver_res_now,
                                      &(nbr_list->verlist.nver_lst_res),
                                     nbr_list_nter_res[natm_tot],nmem_min_lst,
                                     npad_tot,&irealloc,myid_forc);
                 if(irealloc==1){
                   nbr_list->verlist.jter_res = (list_int *)
                          crealloc(&(nbr_list->verlist.jter_res)[1],
                      (nbr_list->verlist.nver_lst_res)*sizeof(list_int))-1;
                 }/*endif*/
	      }/*endif:RESPA*/

           }/*endif:for call*/         

        }/*endwhile: accumulating interactions*/

/*==========================================================================*/

      }/*endfor:iatm_shft*/
   }/*endfor:ishft*/

/*==========================================================================*/
/* VI) Final Dump */

  if((*for_scr_intact)>0){

     verlist_create(clatoms_info,clatoms_pos,
                    intra_scr,for_scr,atommaps,cell,
                    nbr_list,timeinfo,interact,iswitch,
                    &nver_now,&nver_res_now,lst_typ);
     num_call          += 1;
     intact_tot        += (*for_scr_intact);

     verpad_list_mem_check(nver_now,&(nbr_list->verlist.nver_lst),
                           nbr_list_nter[natm_tot],nmem_min_lst,
                           npad_tot,&irealloc,myid_forc);
     if(irealloc==1){
        nbr_list->verlist.jter = (list_int *)
                 crealloc(&(nbr_list->verlist.jter)[1],
                  (nbr_list->verlist.nver_lst)*sizeof(list_int))-1;
     }/*endif*/
      
     if(int_res_ter==1){
         verpad_list_mem_check(nver_res_now,&(nbr_list->verlist.nver_lst_res),
                               nbr_list_nter_res[natm_tot],nmem_min_lst,
                               npad_tot,&irealloc,myid_forc);
         if(irealloc==1){
             nbr_list->verlist.jter_res = (list_int *)
                          crealloc(&(nbr_list->verlist.jter_res)[1],
                       (nbr_list->verlist.nver_lst_res)*sizeof(list_int))-1;
         }/*endif*/
     }/*endif:RESPA*/

  }/*endif*/

/*==========================================================================*/
/* VI) Unpad offsets                                                        */

  for(jpart=1;jpart<=natm_tot;jpart++){
     nbr_list_nter[jpart]  -=  nbr_list_jver_off[jpart];
  }/*endfor*/

  if(int_res_ter==1){
    for(jpart=1;jpart<=natm_tot;jpart++){
      nbr_list_nter_res[jpart]  -= nbr_list_jver_off_res[jpart];
    }/*endfor*/
  }/*endif*/

/*==========================================================================*/
/* VII) Check for successful compleation and unpad the list                 */

/*--------------------------------------------------------------------------*/
/* A) Regular list */

   nbr_list_jter     = nbr_list->verlist.jter;   /* Assigned after realoc*/
   (*isuccess)       = 1;
   jver_pad_new      = jver_pad;
   ioff_tot          = 0;
   nver_now          = 0;
   lower             = 1;
   for(jpart=2;jpart<=natm_tot;jpart++){
      mtemp    = jpart-1;
      nver_now+= nbr_list_nter[mtemp];    
      jver_tmp = lower+nbr_list_nter[mtemp]-1;  
      iver_tmp = jver_tmp-nbr_list_jver_off[jpart];
      if(iver_tmp>0){(*isuccess)=0;}
      i1=jver_pad_new;i2=jver_pad+iver_tmp;jver_pad_new=MAX(i1,i2);
      lower         = nbr_list_jver_off[jpart]+1;
      upper         = lower+nbr_list_nter[jpart]-1;
      ioff_tot     += iver_tmp;
      nbr_list_jver_off[jpart] = nver_now;
      if((*isuccess==1)){
          for(iter=lower;iter<=upper;iter++){
             mtemp =  iter+ioff_tot;
             nbr_list_jter[mtemp] = nbr_list_jter[iter];  
          }/*endfor*/
      }/*endif*/
   }/*endfor*/
   nver_now += nbr_list_nter[natm_tot];

/*--------------------------------------------------------------------------*/
/* B) Respa list */
   if(int_res_ter==1){
      nbr_list_jter_res  = nbr_list->verlist.jter_res; 
                                           /* Assigned after realoc*/
      ioff_tot          = 0;
      nver_res_now      = 0;
      lower             = 1;
      for(jpart=2;jpart<=natm_tot;jpart++){
         mtemp    = jpart-1;
         nver_res_now+= nbr_list_nter_res[mtemp];    
         jver_tmp = lower+nbr_list_nter_res[mtemp]-1;  
         iver_tmp = jver_tmp-nbr_list_jver_off_res[jpart];
         if(iver_tmp>0){(*isuccess)=0;}
         i1=jver_pad_new;i2=jver_pad+iver_tmp;jver_pad_new=MAX(i1,i2);
         lower         = nbr_list_jver_off_res[jpart]+1;
         upper         = lower+nbr_list_nter_res[jpart]-1;
         ioff_tot     += iver_tmp;
         nbr_list_jver_off_res[jpart] = nver_res_now;
         if((*isuccess==1)){
            for(iter=lower;iter<=upper;iter++){
               mtemp =  iter+ioff_tot;
               nbr_list_jter_res[mtemp] = nbr_list_jter_res[iter];  
            }/*endfor*/
         }/*endif*/
      }/*endfor*/
      nver_res_now += nbr_list_nter_res[natm_tot];
   }/*endif:RESPA*/

/*--------------------------------------------------------------------------*/
/* C) Assign new padding  */
   nbr_list->verlist.jver_pad = jver_pad_new; 
   npad_tot                   = jver_pad_new*(natm_tot-1);

/*==========================================================================*/
/* VIII)Assign present size of list to appropriate variable                 */

   (*nver)     = nver_now;
   (*nver_res) = nver_res_now;

/*==========================================================================*/
/* IX)Check list sizes one more time using new pad value                    */

   verpad_list_mem_check(nver_now,&(nbr_list->verlist.nver_lst),
                         nbr_list_nter[natm_tot],nmem_min_lst,
                         npad_tot,&irealloc,myid_forc);
   if(irealloc==1){
        nbr_list->verlist.jter = (list_int *)
                 crealloc(&(nbr_list->verlist.jter)[1],
                     (nbr_list->verlist.nver_lst)*sizeof(list_int))-1;
   }/*endif*/
      
   if(int_res_ter==1){
     verpad_list_mem_check(nver_res_now,&(nbr_list->verlist.nver_lst_res),
                           nbr_list_nter_res[natm_tot],nmem_min_lst,
                           npad_tot,&irealloc,myid_forc);
     if(irealloc==1){
        nbr_list->verlist.jter_res = (list_int *)
                 crealloc(&(nbr_list->verlist.jter_res)[1],
                      (nbr_list->verlist.nver_lst_res)*sizeof(list_int))-1;
     }/*endif*/
   }/*endif:RESPA*/

/*--------------------------------------------------------------------------*/
   }/*end routine */
/*==========================================================================*/





/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void lnk_ver_gen_pad_root(CLATOMS_INFO *clatoms_info,CLATOMS_POS *clatoms_pos,
                          FOR_SCR *for_scr,ATOMMAPS *atommaps,
                          CELL *cell,INTERACT *interact,
                          TIMEINFO *timeinfo,NBR_LIST *nbr_list,
                          EXCL *excl,INTRA_SCR *intra_scr, int *nver, 
                          int *nver_res, int *isuccess,int error_check_on,
                          CLASS_COMM_FORC_PKG *class_comm_forc_pkg)

/*==========================================================================*/
/*                      Brief description                                   */

/* lnk_list[ia][ib][ic][iatm] stores the indices of iatmth atm in           */
/*                           a cell labled ia,ib,ic                         */
/*                           (a pipette cut out of the full system cell).   */
/*                           It is packed with zeros for                    */
/*                           iatm greater than the number of atoms          */
/*                           in a given cell for computational              */
/*                           efficiency. Zero index atoms give              */
/*                           rise to ``null'' interactions that are         */
/*                           removed with logic                             */
/* ishft_a[i],ishft_b[i],ishft_c[i]                                      */
/*                           are vectors containing the allowed sets of     */
/*                           shifts along a,b,c                             */
/* iexl_check[i]             indicates whether or not exclusions needed     */
/*                           to be removed in a particular shift set        */
/* When the lnk_list is shifted interactions are generated between cells    */
/* of the shifted list and old unshifted list.                              */
/* These interactions are sent to the force routine in bundles of nlen.     */
/* Break points where memory conflicts can occur when summing the forces    */
/* are tablulated.                                                          */

/*==========================================================================*/
{/*Begin Routine*/
/*==========================================================================*/
/*            Local Variable declarations                                   */

   int itemp,nlen_use;              /* value of nlen to be used           */
   int ncell_tot, nlist_tot;        /* tot number of cells, size of list  */
   int natm_shft;                   /* number of atom shifts              */
   int iadd_shft;                   /* atm shift flag                     */
   int ishft, iatm_shft;            /* cell and atom shift loops indicies  */
   int ja,jb,jc,jatm;               /* the present shift                  */
   int icheck,ires_flag;            /* excl check on                      */
   int intact_now, intact_tmp;      /* interaction temporaries            */
   int intact_tot, num_call;        /* tot interactions,num of force calls*/
   int i;                           /* loop counter                       */
   int imin,jmax;                   /* min/max atm ind in a interact pair */
   int nact;                        /* interaction is good flag           */
   int iact;                        /* excl interaction loop              */
   int lower,upper;                 /* limits on excl interaction loop    */
   int ifor_call;                   /* force call flag                    */
   int jind_off;                    /* integer offset                     */
   int iswitch;
   int jpart;
   int nver_now,nver_res_now;
   int iii;
   int jver_pad_new;
   int ioff_tot,iter,irealloc;
   int jver_tmp,iver_tmp;
   int i1,i2;
   int npad_tot;
   int intact_init,intact_half,nlist_tot2,mtemp;
   int ioff,ntemp;
   int lst_typ;
   int iver_count;
   int iver_fill;
   int my_ind,jind,ind,ngo,my_root,j,k;
   int ishft_st,nskip;

/* Local Pointers */

   int ncell_a       = nbr_list->lnklist.ncell_a;
   int ncell_b       = nbr_list->lnklist.ncell_b;
   int ncell_c       = nbr_list->lnklist.ncell_c;
   int natm_cell_max = nbr_list->lnklist.natm_cell_max;
   int nshft_lnk     = nbr_list->lnklist.nshft_lnk;

   int excl_nlst              = excl->nlst;
   int *for_scr_intact        = &(for_scr->intact);
   int *for_scr_num_brk       = &(for_scr->num_brk);
   double *for_scr_wght_lnk   = &(for_scr->wght_lnk);
   int *for_scr_num_brk_i     = for_scr->num_brk_i;
   int *for_scr_i_index       = for_scr->i_index;  
   int *for_scr_j_index       = for_scr->j_index;  
   int *excl_j_off            = excl->j_off; 
   int *excl_num              = excl->num;
   int *excl_j                = excl->j;
   list_int *for_scr_i_lnk         = for_scr->i_lnk;
   list_int *for_scr_j_lnk         = for_scr->j_lnk;
   int *nbr_list_nter         = nbr_list->verlist.nter;
   int *nbr_list_jver_off     = nbr_list->verlist.jver_off;
   int *nbr_list_nter_res     = nbr_list->verlist.nter_res;
   int *nbr_list_jver_off_res = nbr_list->verlist.jver_off_res;
   list_int *nbr_list_lnk_list     = nbr_list->lnklist.lnk_list;
   int *nbr_list_ishft_a      = nbr_list->lnklist.ishft_a;
   int *nbr_list_ishft_b      = nbr_list->lnklist.ishft_b;
   int *nbr_list_ishft_c      = nbr_list->lnklist.ishft_c;
   int *nbr_list_iexl_chk     = nbr_list->lnklist.iexl_chk;
   double *nbr_list_shft_wght = nbr_list->lnklist.shft_wght;
   int *for_scr_iexcl         = for_scr->iexcl;
   list_int *nbr_list_jter;        /* Assigned below because of reallocs */
   list_int *nbr_list_jter_res;


    int *brnch_atm_root        = nbr_list->brnch_root.brnch_atm_root;
    int *brnch_atm_list        = nbr_list->brnch_root.brnch_atm_list;
    int *natm_add              = nbr_list->brnch_root.natm_add;
    int *iatm_add              = nbr_list->brnch_root.iatm_add;
    int *iatm_add_off          = nbr_list->brnch_root.iatm_add_off;
    int *root_atm_list         = nbr_list->brnch_root.root_atm_list;
    int *root_atm_map          = nbr_list->brnch_root.root_atm_map;
    int *brnch_atm_map         = nbr_list->brnch_root.brnch_atm_map;
    int **ibrnch_of_root       = nbr_list->brnch_root.ibrnch_of_root;
    int *nbrnch_of_root        = nbr_list->brnch_root.nbrnch_of_root;
    int nbrnch_tot             = nbr_list->brnch_root.nbrnch_tot;
    int nroot_tot              = nbr_list->brnch_root.nroot_tot;

   int int_res_ter            = timeinfo->int_res_ter;
   int natm_tot               = clatoms_info->natm_tot;  
   int nlen                   = for_scr->nlen;
   int nmem_min_lst           = nbr_list->verlist.nmem_min_lst;
   int iver_init              = nbr_list->verlist.iver_init;
   double mem_safe            = nbr_list->verlist.mem_safe;
   int jver_pad               = nbr_list->verlist.jver_pad;
   int nver_lst               = nbr_list->verlist.nver_lst;
   int nver_lst_res           = nbr_list->verlist.nver_lst_res;
   int brnch_root_list_opt    = nbr_list->brnch_root_list_opt; 
   int np_forc                = class_comm_forc_pkg->num_proc;
   int myid_forc              = class_comm_forc_pkg->myid;
   MPI_Comm comm_forc         = class_comm_forc_pkg->comm;
   int nolst_ver_update  = nbr_list->verlist.nolst_ver_update;

/*==========================================================================*/
/* I) Get a nice value of nlen: a value divisible by natm_tot or vice versa */

  if(nlen>nroot_tot){
    nlen_use =  nlen - ( nlen % nroot_tot );
  }else{
    itemp    = nroot_tot/nlen + 1;
    nlen_use = nroot_tot/itemp;
  }/*endif*/

/*==========================================================================*/
/* II) Get some useful constants*/

  ncell_tot = ncell_a*ncell_b*ncell_c;
  nlist_tot = ncell_tot*natm_cell_max;

  iver_fill  = 1;
  iver_count = 0;
  nbr_list->verlist.iver_fill=1;    /* This routine fills only !!!!!!! */  
  nbr_list->verlist.iver_count=0;   

/*==========================================================================*/
/* III) Inititialize scalars          */
 
  nver_now           = 0;
  nver_res_now       = 0;
  (*for_scr_intact)  = 0;   
  (*for_scr_num_brk) = 0;
  num_call           = 0;
  intact_tot         = 0;
  iswitch            = 1;
  lst_typ            = 3;

/*==========================================================================*/
/* III) Inititialize lists          */

/*----------------------------------*/
/* Pad the offset lists             */

  nbr_list_jver_off[1] = 0;
  for(jpart=2;jpart<=natm_tot;jpart++){
     nbr_list_jver_off[jpart] += jver_pad*(jpart-1);
  }/*endfor*/

  if(int_res_ter==1){
    nbr_list_jver_off_res[1] = 0;
    for(jpart=2;jpart<=natm_tot;jpart++){
      nbr_list_jver_off_res[jpart] += jver_pad*(jpart-1);
    }/*endfor*/
  }/*endif*/    
  npad_tot = jver_pad*(natm_tot-1);

/*----------------------------------------------------*/
/* Set the interactions/particle equal to the offsets */

  for(jpart=1;jpart<=natm_tot;jpart++){
     nbr_list_nter[jpart]     = nbr_list_jver_off[jpart];
  }/*endfor*/
  if( (nbr_list_nter[natm_tot]+nmem_min_lst) > nver_lst){
    printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
    printf("Neighor list incorrectly allocated\n");
    printf("Contact Technical Support \n");
    printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
    fflush(stdout);
    exit(1);
  }/*endif*/

  if(int_res_ter==1){
     for(jpart=1;jpart<=natm_tot;jpart++){
       nbr_list_nter_res[jpart] = nbr_list_jver_off_res[jpart];
     }/*endfor*/
     if((nbr_list_nter_res[natm_tot]+nmem_min_lst)>nver_lst_res){
        printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
        printf("RESPA Neighor list incorrectly allocated\n");
        printf("Contact Technical Support \n");
        printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
        fflush(stdout);
        exit(1);
     }/*endif*/
  }/*endif*/

/*------------------------------------*/
/* Find half a shift                  */

  intact_init = nroot_tot;
  nlist_tot2  = 4*natm_cell_max*ncell_a*ncell_b*ncell_c
              - 4*ncell_a*ncell_b*ncell_c
              - 2*ncell_a*ncell_b
              - ncell_a;
  intact_half = 0; 
  for(i=1;i<=nroot_tot;i++){
    if(for_scr_iexcl[i]<=nlist_tot2){intact_half=i;}
  }/*endif*/

/*==========================================================================*/
/* IV) Loop over the allowed cell shifts: ja,jb,jc                          */

   ishft_st = 1 + myid_forc;
   nskip    = np_forc;
   for(ishft=ishft_st;ishft<=nshft_lnk;ishft+=nskip){

      ja                  = nbr_list_ishft_a[ishft];
      jb                  = nbr_list_ishft_b[ishft];
      jc                  = nbr_list_ishft_c[ishft];
      icheck              = nbr_list_iexl_chk[ishft];
      (*for_scr_wght_lnk) = nbr_list_shft_wght[ishft];

/*==========================================================================*/
/* V) Loop over all the atm shifts: jatm                                    */
/*     Bookeeping for the first cell shift to take care                     */
/*           of self interactions                                           */

      natm_shft = natm_cell_max;
      iadd_shft = 0;
      if(ishft==1){
         iadd_shft  = ((natm_cell_max-1) % 2);
         natm_shft  = (natm_cell_max-1)/2 + iadd_shft;
      }/*endif*/
      for(iatm_shft=1;iatm_shft<=natm_shft;iatm_shft++){
         jatm = iatm_shft;

/*==========================================================================*/
/*  B) Book keeping for last atm shift of 1st cell shift:                  */
/*         On the last iatm_shft of the 1st cell shift get rid of           */
/*         repeated interactions that occur if natm_cell_max is even.       */
/*         This simply means taking 1/2 of the interactions found           */

         intact_now = intact_init;
         if((ishft==1)&&(iadd_shft==1)&&(iatm_shft==natm_shft))
                 {intact_now = intact_half;}

/*==========================================================================*/
/*   C) Eliminate null interactions caused by packing the lnk_list          */
/*       with zero index atoms                                              */

         intact_tmp = intact_now;
         intact_now = 0;
         ioff =  ja + jb*ncell_a*2 + jc*ncell_a*ncell_b*4 
               + jatm*ncell_a*ncell_b*ncell_c*8;
         for(i=1;i<=intact_tmp;i++){
           itemp = for_scr_iexcl[i];
           ntemp = itemp+ioff;
           mtemp = nbr_list_lnk_list[ntemp];
           if(mtemp!=0){
            intact_now += 1;
               for_scr_i_lnk[intact_now] = nbr_list_lnk_list[itemp];
               for_scr_j_lnk[intact_now] = mtemp;
           }/*endif*/
         }/*endfor*/

/*==========================================================================*/
/*   E) Eliminate standard exclusions in shifts that might have them.       */
/*       Cells that are far enough away in space won't have any excls.      */

         if((icheck==1)&&(excl_nlst>0)){
            intact_tmp = intact_now;
            intact_now  = 0;
            for(i=1;i<=intact_tmp;i++){
               imin  = MIN( (for_scr_i_lnk)[i],(for_scr_j_lnk)[i] );
               jmax  = MAX( (for_scr_i_lnk)[i],(for_scr_j_lnk)[i] );
               lower = (excl_j_off)[jmax]+1;
               upper = (excl_num)[jmax]+(excl_j_off)[jmax];
               nact  = 1;
               for(iact=lower;iact<=upper;iact++){
                 if(imin==(excl_j)[iact]){nact = 0;}
               }/*endfor*/
               if(nact==1){
                 intact_now += 1;
                 (for_scr_i_lnk)[intact_now] = (for_scr_i_lnk)[i];
                 (for_scr_j_lnk)[intact_now] = (for_scr_j_lnk)[i];
               }/*endif*/
            }/*endfor*/
         }/*endif*/
/*==========================================================================*/
/*  G) Add interactions to the temporary interaction list                   */
/*     in chunks no greater than nlen_use: Keep track of the break points   */

        lower = 1;
        upper = 0;
        while(upper!=intact_now){
           upper = intact_now;
           if(upper-lower+1+(*for_scr_intact)>nlen_use){
               upper = nlen_use-(*for_scr_intact)+lower-1;
           }/*endif*/
           jind_off = (*for_scr_intact)-lower+1;
           for(i=lower;i<=upper;i++){
              mtemp = (i+jind_off);
              (for_scr_j_index)[mtemp]=(for_scr_j_lnk)[i];
              (for_scr_i_index)[mtemp]=(for_scr_i_lnk)[i];
           }/*endfor*/
           (*for_scr_intact)                       += upper-lower+1;
           (*for_scr_num_brk)                      += 1;
           (for_scr_num_brk_i)[(*for_scr_num_brk)] = upper-lower+1;
           lower                                   = upper+1;

/*==========================================================================*/
/*  H) If enough interactions have accumulated (or you are done)            */
/*       make a verlist create call and then update and reinitialize        */
/*       appropriate counters. Realloc if the tot num of interactions       */
/*       is too large                                                       */

           ifor_call = 0;
           if( (*for_scr_intact)>=nlen_use){ifor_call=1;}
           if(  (ishft==(nshft_lnk))&&
                 (iatm_shft==natm_shft)&&
                 (upper==intact_now)&&
                 ((*for_scr_intact)!=0) ) {ifor_call=1;}

           if(ifor_call==1){
             verlist_create(clatoms_info,clatoms_pos,
                            intra_scr,for_scr,atommaps,cell,
                            nbr_list,timeinfo,interact,iswitch,
                            &nver_now,&nver_res_now,lst_typ);
              num_call          += 1;
              intact_tot        += (*for_scr_intact);
              (*for_scr_intact)  = 0;   
              (*for_scr_num_brk) = 0;

              verpad_list_mem_check(nver_now,&(nbr_list->verlist.nver_lst),
                                     nbr_list_nter[natm_tot],nmem_min_lst,
                                     npad_tot,&irealloc,myid_forc);
              if(irealloc==1){
                  nbr_list->verlist.jter = (list_int *)
                      crealloc(&(nbr_list->verlist.jter)[1],
                     (nbr_list->verlist.nver_lst)*sizeof(list_int))-1;
              }/*endif*/
            
              if(int_res_ter==1){
                 verpad_list_mem_check(nver_res_now,
                                      &(nbr_list->verlist.nver_lst_res),
                                     nbr_list_nter_res[natm_tot],nmem_min_lst,
                                     npad_tot,&irealloc,myid_forc);
                 if(irealloc==1){
                   nbr_list->verlist.jter_res = (list_int *)
                          crealloc(&(nbr_list->verlist.jter_res)[1],
                     (nbr_list->verlist.nver_lst_res)*sizeof(list_int))-1;
                 }/*endif*/
	      }/*endif:RESPA*/

           }/*endif:for call*/         

        }/*endwhile: accumulating interactions*/

/*==========================================================================*/

      }/*endfor:iatm_shft*/
   }/*endfor:ishft*/

/*==========================================================================*/
/* VI) Final Dump */

  if((*for_scr_intact)>0){

     verlist_create(clatoms_info,clatoms_pos,
                    intra_scr,for_scr,atommaps,cell,
                    nbr_list,timeinfo,interact,iswitch,
                    &nver_now,&nver_res_now,lst_typ);
     num_call          += 1;
     intact_tot        += (*for_scr_intact);

     verpad_list_mem_check(nver_now,&(nbr_list->verlist.nver_lst),
                           nbr_list_nter[natm_tot],nmem_min_lst,
                           npad_tot,&irealloc,myid_forc);
     if(irealloc==1){
        nbr_list->verlist.jter = (list_int *)
                 crealloc(&(nbr_list->verlist.jter)[1],
                     (nbr_list->verlist.nver_lst)*sizeof(list_int))-1;
     }/*endif*/
      
     if(int_res_ter==1){
         verpad_list_mem_check(nver_res_now,&(nbr_list->verlist.nver_lst_res),
                               nbr_list_nter_res[natm_tot],nmem_min_lst,
                               npad_tot,&irealloc,myid_forc);
         if(irealloc==1){
             nbr_list->verlist.jter_res = (list_int *)
                          crealloc(&(nbr_list->verlist.jter_res)[1],
                      (nbr_list->verlist.nver_lst_res)*sizeof(list_int))-1;
         }/*endif*/
     }/*endif:RESPA*/

  }/*endif*/

/*==========================================================================*/
/*Fill only */

  ires_flag = 0;
  expand_brnch_root_fll(nbr_list,ires_flag,&nver_now,
                          nbrnch_tot,brnch_atm_root,brnch_atm_list,
                          brnch_atm_map,
                          nroot_tot,root_atm_list,nbrnch_of_root,
                          ibrnch_of_root,natm_tot,natm_add,
                          iatm_add_off,iatm_add,root_atm_map,
                          irealloc,iver_fill,
                          nmem_min_lst,iver_init,mem_safe,
                          class_comm_forc_pkg);
  if(int_res_ter==1){
    ires_flag = 1;
    expand_brnch_root_fll(nbr_list,ires_flag,&nver_res_now,
                            nbrnch_tot,brnch_atm_root,brnch_atm_list,
                            brnch_atm_map,
                            nroot_tot,root_atm_list,nbrnch_of_root,
                            ibrnch_of_root,natm_tot,natm_add,
                            iatm_add_off,iatm_add,root_atm_map,
                            irealloc,iver_fill,
                            nmem_min_lst,iver_init,mem_safe,
                            class_comm_forc_pkg);

  }/*endif:RESPA*/

/*==========================================================================*/
/* VII) Check for successful compleation and unpad the list                 */

/*--------------------------------------------------------------------------*/
/* A) Regular list */

   nbr_list_jter     = nbr_list->verlist.jter;   /* Assigned after realoc*/
   (*isuccess)       = 1;
   jver_pad_new      = jver_pad;
   ioff_tot          = 0;
   nver_now          = 0;
   lower             = 1;
   for(jpart=2;jpart<=natm_tot;jpart++){
      mtemp    = jpart-1;
      nver_now+= nbr_list_nter[mtemp];    
      jver_tmp = lower+nbr_list_nter[mtemp]-1;  
      iver_tmp = jver_tmp-nbr_list_jver_off[jpart];
      if(iver_tmp>0){(*isuccess)=0;}
      i1=jver_pad_new;i2=jver_pad+iver_tmp;jver_pad_new=MAX(i1,i2);
      lower         = nbr_list_jver_off[jpart]+1;
      upper         = lower+nbr_list_nter[jpart]-1;
      ioff_tot     += iver_tmp;
      nbr_list_jver_off[jpart] = nver_now;
      if((*isuccess==1)){
          for(iter=lower;iter<=upper;iter++){
             mtemp =  iter+ioff_tot;
             nbr_list_jter[mtemp] = nbr_list_jter[iter];  
          }/*endfor*/
      }/*endif*/
   }/*endfor*/
   nver_now += nbr_list_nter[natm_tot];

/*--------------------------------------------------------------------------*/
/* B) Respa list */
   if(int_res_ter==1){
      nbr_list_jter_res  = nbr_list->verlist.jter_res; 
                                           /* Assigned after realoc*/
      ioff_tot          = 0;
      nver_res_now      = 0;
      lower             = 1;
      for(jpart=2;jpart<=natm_tot;jpart++){
         mtemp    = jpart-1;
         nver_res_now+= nbr_list_nter_res[mtemp];    
         jver_tmp = lower+nbr_list_nter_res[mtemp]-1;  
         iver_tmp = jver_tmp-nbr_list_jver_off_res[jpart];
         if(iver_tmp>0){(*isuccess)=0;}
         i1=jver_pad_new;i2=jver_pad+iver_tmp;jver_pad_new=MAX(i1,i2);
         lower         = nbr_list_jver_off_res[jpart]+1;
         upper         = lower+nbr_list_nter_res[jpart]-1;
         ioff_tot     += iver_tmp;
         nbr_list_jver_off_res[jpart] = nver_res_now;
         if((*isuccess==1)){
            for(iter=lower;iter<=upper;iter++){
               mtemp =  iter+ioff_tot;
               nbr_list_jter_res[mtemp] = nbr_list_jter_res[iter];  
            }/*endfor*/
         }/*endif*/
      }/*endfor*/
      nver_res_now += nbr_list_nter_res[natm_tot];
   }/*endif:RESPA*/

/*--------------------------------------------------------------------------*/
/* C) Assign new padding  */
   nbr_list->verlist.jver_pad = jver_pad_new; 
   npad_tot                   = jver_pad_new*(natm_tot-1);

/*==========================================================================*/
/* VIII)Assign present size of list to appropriate variable                 */

   (*nver)     = nver_now;
   (*nver_res) = nver_res_now;

/*==========================================================================*/
/* IX)Check list sizes one more time using new pad value                    */


   verpad_list_mem_check(nver_now,&(nbr_list->verlist.nver_lst),
                         nbr_list_nter[natm_tot],nmem_min_lst,
                         npad_tot,&irealloc,myid_forc);
   if(irealloc==1){
        nbr_list->verlist.jter = (list_int *)
                 crealloc(&(nbr_list->verlist.jter)[1],
                 (nbr_list->verlist.nver_lst)*sizeof(list_int))-1;
   }/*endif*/
      
   if(int_res_ter==1){
     verpad_list_mem_check(nver_res_now,&(nbr_list->verlist.nver_lst_res),
                           nbr_list_nter_res[natm_tot],nmem_min_lst,
                           npad_tot,&irealloc,myid_forc);
     if(irealloc==1){
        nbr_list->verlist.jter_res = (list_int *)
                 crealloc(&(nbr_list->verlist.jter_res)[1],
                    (nbr_list->verlist.nver_lst_res)*sizeof(list_int))-1;
     }/*endif*/
   }/*endif:RESPA*/
  
/*==========================================================================*/
/* VI) Pare down the list if filling                                        */
   
   if(iver_fill==1 && brnch_root_list_opt==2 && (*isuccess==1)){
      ver_pare_down_root(clatoms_info,clatoms_pos,for_scr,atommaps,
                        cell,interact,timeinfo,nbr_list,excl,
                        intra_scr,nver,nver_res,myid_forc); 
   }/*endif*/


/*--------------------------------------------------------------------------*/
   }/*end routine */
/*==========================================================================*/




/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void ver_list_mem_check(int nver_now,int *nver_lst,int iver_fill,
                        int nmem_min_lst,int iver_init,
                        double mem_safe,int *irealloc,int myid)

/*==========================================================================*/
{/*Begin Routine*/

  int i1;

/*==========================================================================*/
/* I) Over flow check:            */

 if( ( nver_now>(*nver_lst) )&&( iver_fill==1 ) ){
    printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
    printf("Over flow in neighbor list routine on proc %d\n",myid);
    printf("Contact Technical Support \n");
    printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
    fflush(stdout);
    exit(1);
 }/*endif*/

/*==========================================================================*/
/* II) Safe memory realloc:                                                 */

 (*irealloc) = 0;
 if( ( (nver_now+nmem_min_lst)>(*nver_lst) )&&( iver_fill==1 ) ){
    i1          = (*nver_lst);
    (*nver_lst) = (nver_now+2*nmem_min_lst)*mem_safe;
    if(iver_init==0){
      printf("\n====================================\n");
      printf("Allocating more verlist memory on proc %d: %d %d\n",
              myid,i1,(*nver_lst) );
      printf("====================================\n\n"); 
    }/*endif*/
    (*irealloc) = 1;
 }/*endif:ver memory*/        

/*--------------------------------------------------------------------------*/
   }/*end routine */
/*==========================================================================*/



/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void verpad_list_mem_check(int nver_now,int *nver_lst,int nter_last,
                           int nmem_min_lst,int npad_tot, int *irealloc,
                           int myid)

/*==========================================================================*/
{/*Begin Routine*/

  int i1,nmax,i2;

/*==========================================================================*/
/* I) Over flow check:            */

  if( ( nver_now>(*nver_lst) ) || (nter_last >(*nver_lst) ) ){
    printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
    printf("Over flow in neighbor list routine on proc %d\n",myid);
    printf("Contact Technical Support \n");
    printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
    fflush(stdout);
    exit(1);
  }/*endif*/

/*==========================================================================*/
/* II) Safe memory realloc:                                                 */

  (*irealloc) = 0;
  if( ( (nver_now+nmem_min_lst+npad_tot)>(*nver_lst) )||
      ( (nter_last+nmem_min_lst)        >(*nver_lst) )   ){
      i2          = (nver_now+npad_tot);
      nmax        = MAX(i2,nter_last);
      i1          = (*nver_lst);
      (*nver_lst) =  (nmax+2*nmem_min_lst);
      printf("\n====================================\n"); 
      printf("Allocating more verlist memory on proc %d:  %d %d\n",
                           myid,i1,(*nver_lst));
      printf("====================================\n\n"); 
      fflush(stdout);
     (*irealloc)  = 1;
   }/*endif:ver memory*/        

/*--------------------------------------------------------------------------*/
   }/*end routine */
/*==========================================================================*/








/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void expand_brnch_root_cnt(int *nver_now,int *nbr_list_nter,
                           int *nbr_list_jver_off,
                           int nbrnch_tot,int *brnch_atm_root,
                           int *brnch_atm_list,
                           int nroot_tot ,int *root_atm_list,
                           int natm_tot,int *natm_add,
                           int **ibrnch_of_root,int *nbrnch_of_root,
                           CLASS_COMM_FORC_PKG *class_comm_forc_pkg,
                           int nolst_ver_update)

/*==========================================================================*/
   {/*Begin Routine*/
/*==========================================================================*/
/*              Local variables */

  int my_root,my_ind,i,jpart,mtemp,iroot;
  int iroot_st,nskip,iii;

/*==========================================================================*/
/*              Local pointers                                              */

  int np_forc    = class_comm_forc_pkg->num_proc;
  int myid_forc  = class_comm_forc_pkg->myid;

/*==========================================================================*/
/* I) Set the # of interactions of each particle.                           */
/*    Give each branch same # of interactions as its root.                  */
/*    Add the number of additions (natm_add) as well                        */

/* First,assign the # of branch neighbors equal to the # of its neighbors*/
/* of its root plus its own additions */

  if(nolst_ver_update==1){
    iroot_st = 1+myid_forc;
    nskip    = np_forc;
  }else{
    iroot_st = 1;
    nskip    = 1;
    if(myid_forc>0){for(i=1;i<=natm_tot;i++){natm_add[i]=0;}}
  }/*endif*/

  for(iroot=iroot_st;iroot<=nroot_tot;iroot+=nskip){
      my_root = root_atm_list[iroot];
      for(i=1;i<=nbrnch_of_root[iroot];i++){
        my_ind    = ibrnch_of_root[iroot][i];
        nbr_list_nter[my_ind] = nbr_list_nter[my_root] + natm_add[my_ind];
        (*nver_now) += (nbr_list_nter[my_ind]);
     }/*endfor*/
  }/*endfor*/

/* Next, the # of Root atom neighbors can be updated */
  for(iroot=iroot_st;iroot<=nroot_tot;iroot+=nskip){
      my_ind                = root_atm_list[iroot];
      nbr_list_nter[my_ind] += natm_add[my_ind];
      (*nver_now)           += natm_add[my_ind];
  }/*endfor*/

/*==========================================================================*/
/* II) Calculate the off set into the root-root interaction list            */

  nbr_list_jver_off[1] = 0;
  for(jpart=2;jpart<=natm_tot;jpart++){
      mtemp = jpart-1;
      nbr_list_jver_off[jpart] = nbr_list_jver_off[mtemp]
                               + nbr_list_nter[mtemp];
  }/*endfor*/

/*--------------------------------------------------------------------------*/
   }/*end routine */
/*==========================================================================*/




/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void expand_brnch_root_cntfll(NBR_LIST *nbr_list,int ires_flag,
                              int *nver_now_ret,
                              int nbrnch_tot,int *brnch_atm_root,
                              int *brnch_atm_list,int *brnch_atm_map,
                              int nroot_tot,int *root_atm_list,
                              int *nbrnch_of_root,int **ibrnch_of_root,
                              int natm_tot,int *natm_add,
                              int *iatm_add_off,int *iatm_add,
                              int *root_atm_map,
                              int iver_fill,int nmem_min_lst,
                              int iver_init,double mem_safe,int *iamgood,
                              CLASS_COMM_FORC_PKG *class_comm_forc_pkg)

/*==========================================================================*/
   {/*Begin Routine*/
/*==========================================================================*/
/*              Local variables */

   int jpart,mtemp,ind,ngo,jind,my_ind,i,my_root,nlist_start;
   int natm_add_now,nskip;
   int j,joff,iii,iroot_st,iroot;
   int *nbr_list_jver_off,*nbr_list_nter,*nver_lst;
   list_int *nbr_list_jter;
   int *nbr_list_jver_ioff;
   int nver_now;
   int myid_forc = class_comm_forc_pkg->myid;
   int np_forc   = class_comm_forc_pkg->num_proc;
   int nolst_ver_update  = nbr_list->verlist.nolst_ver_update;
   int irealloc;

/*==========================================================================*/
/*  0) Error check and assign the pointers */

  if(nolst_ver_update!=1){
     printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
     printf("Internal Error: Root count fill called by lnklist update\n ");
     printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
     fflush(stdout);
     exit(1);
  }/*endif*/

  nver_now  = (*nver_now_ret);
  if(ires_flag==0){
   nbr_list_jter      = nbr_list->verlist.jter;
   nbr_list_jver_off  = nbr_list->verlist.jver_off;
   nbr_list_nter      = nbr_list->verlist.nter;
  }/*endif*/
  if(ires_flag==1){
   nbr_list_jter      = nbr_list->verlist.jter_res;
   nbr_list_jver_off  = nbr_list->verlist.jver_off_res;
   nbr_list_nter      = nbr_list->verlist.nter_res;
  }/*endif*/

/*==========================================================================*/
/* I) Calculate the off set into the root-root interaction list            */

  nbr_list_jver_off[1] = 0;
  for(jpart=2;jpart<=natm_tot;jpart++){
      mtemp = jpart-1;
      nbr_list_jver_off[jpart] = nbr_list_jver_off[mtemp]
                               + nbr_list_nter[mtemp];
  }/*endfor*/

/*==========================================================================*/
/* II) Expand root-root list to include the branches off roots              */

  jind = 1;   
  while(jind<=nver_now){
    jpart = nbr_list_jter[jind];    /* index of root neigh in atm list  */
    ind   = root_atm_map[jpart];    /* index of root neigh in root list */ 
    ngo   = nbrnch_of_root[ind];    /* # of branches of rt neigh        */
    for(i=1;i<=ngo;i++){            /* add branches of root neighbor   */
      nbr_list_jter[(jind+i)] = ibrnch_of_root[ind][i];
    }/*endfor*/     
    jind += (ngo + 1);
  }/*endwhile*/     

/*==========================================================================*/
/* III) Calculate the size of the expanded list                             */
/*      Add number of neighbors and ``additions'' for each particle.        */

  for(i=1;i<=natm_tot;i++){iamgood[i]=0;} 
  nver_now = 0;
  iroot_st = 1+myid_forc;
  nskip    = np_forc;
/* Brnch atoms */
  for(iroot=iroot_st;iroot<=nroot_tot;iroot+=nskip){
      my_root = root_atm_list[iroot];
    for(i=1;i<=nbrnch_of_root[iroot];i++){
      my_ind                = ibrnch_of_root[iroot][i];
      iamgood[my_ind]        = 1; 
      nver_now             += (nbr_list_nter[my_root] + natm_add[my_ind]);
      nbr_list_nter[my_ind] = nbr_list_nter[my_root];
    }/*endfor*/
  }/*endfor*/


/* Root atoms */
   for(iroot=iroot_st;iroot<=nroot_tot;iroot+=nskip){
     my_ind         = root_atm_list[iroot];
     iamgood[my_ind] = 1;   
     nver_now      += (nbr_list_nter[my_ind] + natm_add[my_ind]);
  }/*endfor*/


/*==========================================================================*/
/* IV) Chek the memory and reallocate if necessary                          */

  if(ires_flag == 0 ){
    ver_list_mem_check(nver_now,&(nbr_list->verlist.nver_lst),
                       iver_fill,nmem_min_lst,iver_init,mem_safe,
                       &irealloc,myid_forc);
    if(irealloc==1){
     nbr_list->verlist.jter = (list_int *)
            crealloc(&(nbr_list->verlist.jter)[1],
                      (nbr_list->verlist.nver_lst)*sizeof(list_int))-1;
     nbr_list_jter     = nbr_list->verlist.jter;
    }/*endif*/
  }/*endif*/

  if(ires_flag == 1){
    ver_list_mem_check(nver_now,&(nbr_list->verlist.nver_lst_res),
                       iver_fill,nmem_min_lst,iver_init,mem_safe,
                       &irealloc,myid_forc);
    if(irealloc==1){
     nbr_list->verlist.jter_res = (list_int *)
            crealloc(&(nbr_list->verlist.jter_res)[1],
                   (nbr_list->verlist.nver_lst_res)*sizeof(list_int))-1;
     nbr_list_jter     = nbr_list->verlist.jter_res;
    }/*endif*/
  }/*endif*/

/*==========================================================================*/
/* V)    Copy the root interactions to the correct place in list and        */
/*       fill the offset list.                                              */
/*       Start at the back of the list to prevent overwriting.              */
/*       Do the root interactions, first, to prevent overwriting.           */
/*       The roots are ordered but a branch can have an index < its roots.  */

  nlist_start = nver_now;
  for(i=natm_tot;i>=1;i--){
/* root atom fill; compute offset*/
   natm_add_now = (iamgood[i]>0 ? natm_add[i]:0);
   if(root_atm_map[i] > 0){      
      nlist_start -= (nbr_list_nter[i] + natm_add_now);
      for(j=nbr_list_nter[i];j>=1;j--){
       nbr_list_jter[(j+nlist_start)] = 
                  nbr_list_jter[(j+nbr_list_jver_off[i])];
      }/*endfor*/
      nbr_list_jver_off[i] = nlist_start;
/* skip brnch atoms; compute offset*/
    }else{                                   
      nlist_start -= (nbr_list_nter[i] + natm_add_now);
      nbr_list_jver_off[i] = nlist_start;
    }/*endelse*/
  }/*endfor*/

/*==========================================================================*/
/* VI) Now copy the branch interactions to the correct spot in the list   */
/*     The branch has the same list as the roots up till now              */

  for(iroot=iroot_st;iroot<=nroot_tot;iroot+=nskip){
     my_root = root_atm_list[iroot];
    for(i=1;i<=nbrnch_of_root[iroot];i++){
      my_ind  = ibrnch_of_root[iroot][i];
      for(j=1;j<=nbr_list_nter[my_ind];j++){
       nbr_list_jter[(j+nbr_list_jver_off[my_ind])] = 
                     nbr_list_jter[(j+nbr_list_jver_off[my_root])];
      }/*endfor*/
    }/*endfor*/
  }/*endfor*/

/*==========================================================================*/
/* VII) Add the addition-list to the ver_list and adjust the               */
/*       number of neighbors  per atom  (nter or nter_res)                  */

  for(i=1;i<=natm_tot;i++){
    joff = nbr_list_nter[i]+nbr_list_jver_off[i];
    natm_add_now = (iamgood[i]>0 ? natm_add[i]:0);
    for(j=1;j<=natm_add_now;j++){
       nbr_list_jter[(j+joff)] = iatm_add[(j+iatm_add_off[i])];
    }/*endfor*/
    nbr_list_nter[i] += natm_add_now;
  }/*endfor*/

/*==========================================================================*/
/* VIII) Set return value */

 (*nver_now_ret) = nver_now;

/*--------------------------------------------------------------------------*/
   }/*end routine */
/*==========================================================================*/




/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void expand_brnch_root_fll(NBR_LIST *nbr_list,int ires_flag,
                           int *nver_now_ret,
                           int nbrnch_tot,int *brnch_atm_root,
                           int *brnch_atm_list,int *brnch_atm_map,
                           int nroot_tot,int *root_atm_list,
                           int *nbrnch_of_root,int **ibrnch_of_root,
                           int natm_tot,int *natm_add,
                           int *iatm_add_off,int *iatm_add,
                           int *root_atm_map,
                           int irealloc,int iver_fill,int nmem_min_lst,
                           int iver_init,double mem_safe,
                           CLASS_COMM_FORC_PKG *class_comm_forc_pkg)

/*==========================================================================*/
   {/*Begin Routine*/
/*==========================================================================*/
/*              Local variables */

   int jpart,mtemp,ind,ngo,jind,my_ind,my_root,nlist_start;
   int j,joff,k,iii,iroot;
   int *nbr_list_jver_off,*nbr_list_nter,*nver_lst;
   list_int  *nbr_list_jter;
   int nver_now;
   int ipad_ok,iend;

   int myid_forc = class_comm_forc_pkg->myid;
   int np_forc   = class_comm_forc_pkg->num_proc;
   int nolst_ver_update  = nbr_list->verlist.nolst_ver_update;

/*==========================================================================*/
/*  0) Assign the pointers */

  if(nolst_ver_update!=0){
     printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
     printf("Internal Error: Root fill called by nolist update\n ");
     printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
     fflush(stdout);
     exit(1);
  }/*endif*/

  nver_now  = (*nver_now_ret);
  if(ires_flag==0){
   nbr_list_jter      = nbr_list->verlist.jter;
   nbr_list_jver_off  = nbr_list->verlist.jver_off;
   nbr_list_nter      = nbr_list->verlist.nter;
  }/*endif*/
  if(ires_flag==1){
   nbr_list_jter      = nbr_list->verlist.jter_res;
   nbr_list_jver_off  = nbr_list->verlist.jver_off_res;
   nbr_list_nter      = nbr_list->verlist.nter_res;
  }/*endif*/

/*==========================================================================*/
/* I) Check for padding errors */

  ipad_ok  = 1;
  for(jpart=1;jpart<=natm_tot-1;jpart++){
    if(brnch_atm_map[jpart]==0){
     iend = nbr_list_nter[jpart];
    }else{
     my_root = brnch_atm_root[brnch_atm_map[jpart]];
     iend  = (nbr_list_nter[my_root]-nbr_list_jver_off[my_root]
             +nbr_list_jver_off[jpart]);
    }/*endif*/
    if(myid_forc==0){iend+=natm_add[jpart];}
    if(iend>nbr_list_jver_off[(jpart+1)]){ipad_ok=0;}
  }/*endfor*/


/*==========================================================================*/
/* II) Add branches of each root to the list                                 */

 if(ipad_ok==1){
  for(iroot=1;iroot<=nroot_tot;iroot++){
     my_ind = root_atm_list[iroot];       /*my atm index in atm list */
     jind   = nbr_list_jver_off[my_ind]+1;/*starting index of my neighbors*/
     while(jind<=nbr_list_nter[my_ind]){
       jpart = nbr_list_jter[jind];        /*rt neigh atm index in atm list*/
       ind   = root_atm_map[jpart];        /*rt neigh atm index in rt list */ 
       ngo   = nbrnch_of_root[ind];        /* # of branches of rt neigh    */
       for(k=1;k<=ngo;k++){                /* add branches of root neigh   */
         nbr_list_jter[(jind+k)] = ibrnch_of_root[ind][k];
       }/*endfor*/     
       jind += (ngo + 1);
     }/*endwhile*/
  }/*endfor*/
 }/*endfor*/

/*==========================================================================*/
/* III) Assign each branch atm the neighbors of its root                     */

  for(iroot=1;iroot<=nroot_tot;iroot++){
    my_root = root_atm_list[iroot];
    for(k=1;k<=nbrnch_of_root[iroot];k++){
      my_ind    = ibrnch_of_root[iroot][k];
      ngo       = nbr_list_nter[my_root]-nbr_list_jver_off[my_root];
      if(ipad_ok==1){
       for(j=1;j<=ngo;j++){
        nbr_list_jter[(j+nbr_list_jver_off[my_ind])] = 
             nbr_list_jter[(j+nbr_list_jver_off[my_root])];
       }/*endfor*/    
      }/*endif*/
      nbr_list_nter[my_ind]+=ngo;
    }/*endfor*/
  }/*endfor*/    

/*==========================================================================*/
/* IV)  Add the additions  */

 if(myid_forc==0){

  for(iroot=1;iroot<=nroot_tot;iroot++){
   my_root = root_atm_list[iroot];
   joff = nbr_list_nter[my_root];
   if(ipad_ok==1){
    for(j=1;j<=natm_add[my_root];j++){
      nbr_list_jter[(j+joff)] = iatm_add[(j+iatm_add_off[my_root])];
    }/*endfor*/
   }/*endif*/
   nbr_list_nter[my_root] += natm_add[my_root];
  }/*endfor*/

  for(iroot=1;iroot<=nroot_tot;iroot++){
    for(k=1;k<=nbrnch_of_root[iroot];k++){
      ind = ibrnch_of_root[iroot][k];
      joff = nbr_list_nter[ind];
      if(ipad_ok==1){
       for(j=1;j<=natm_add[ind];j++){
         nbr_list_jter[(j+joff)] = iatm_add[(j+iatm_add_off[ind])];
       }/*endfor*/
      }/*endif*/
      nbr_list_nter[ind] += natm_add[ind];
    }/*endfor*/
  }/*endfor*/

 }/*endif*/


/*==========================================================================*/
/* V)  Depad the number of interactions per particle                        */

  nver_now = 0;
  for(jpart=1;jpart<=natm_tot;jpart++){
     nbr_list_nter[jpart]  -=  (nbr_list_jver_off[jpart]);
     nver_now              += nbr_list_nter[jpart];
  }/*endfor*/

/*==========================================================================*/
/* VI) Set return value                                                      */

 (*nver_now_ret) = nver_now;

/*--------------------------------------------------------------------------*/
   }/*end routine */
/*==========================================================================*/






/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void ver_pare_down_root(CLATOMS_INFO *clatoms_info,CLATOMS_POS *clatoms_pos,
                        FOR_SCR *for_scr,ATOMMAPS *atommaps,
                        CELL *cell,INTERACT *interact,
                        TIMEINFO *timeinfo, NBR_LIST *nbr_list, EXCL *excl,
                        INTRA_SCR *intra_scr,int *nver,int *nver_res,
                        int myid_forc)

/*==========================================================================*/
/*                      Brief description                                   */
/*                                                                          */
/* All possible interactions are sent to verlist create routine in bundles  */
/* length of nlen. Break points are tabluated                               */
/*                                                                          */
/*==========================================================================*/
{/*Begin Routine*/
/*==========================================================================*/
/*            Local Variable declarations               */


    int ipart, jpart,i1;           /* ith particle, jth particle          */
    int intact_tot;             /* total number interactions           */
    int intact_save;            /* total number interactions           */
    int num_call;               /*  number of force calls              */
    int num_call_save;          /*  number of force calls              */
    int lower,upper;            /* lower and upper limits on for loop  */
    int ifor_call;              /* force call flag                     */
    int intact_now;             /* interactions now                    */
    int jind_off;
    int iswitch;
    int nver_now,nver_res_now;
    int iii,mtemp,num;
    int jstart,jend,joff,k; 
    int lst_typ,irealloc;
    int iver_count,iver_fill;

/* Local pointers */

    int excl_nlst              = excl->nlst;
    int *for_scr_intact        = &(for_scr->intact);
    int *for_scr_num_brk       = &(for_scr->num_brk);
    int *for_scr_num_brk_i     = for_scr->num_brk_i;
    int *for_scr_iexcl         = for_scr->iexcl;
    int *for_scr_i_index       = for_scr->i_index;  
    int *for_scr_j_index       = for_scr->j_index;  
    int *excl_j_off            = excl->j_off; 
    int *excl_num              = excl->num;
    int *excl_j                = excl->j;
    int *nbr_list_nter         = nbr_list->verlist.nter;
    list_int *nbr_list_jter         = nbr_list->verlist.jter;
    list_int *nbr_list_jter_res     = nbr_list->verlist.jter_res;
    int *nbr_list_jver_off     = nbr_list->verlist.jver_off;
    int *nbr_list_nter_res     = nbr_list->verlist.nter_res;
    int *nbr_list_jver_off_res = nbr_list->verlist.jver_off_res;

    int int_res_ter            = timeinfo->int_res_ter;
    int natm_tot               = clatoms_info->natm_tot;  
    int nlen                   = for_scr->nlen;
    int nmem_min_lst           = nbr_list->verlist.nmem_min_lst;
    int iver_init              = nbr_list->verlist.iver_init;
    double mem_safe            = nbr_list->verlist.mem_safe;

/*==========================================================================*/
/* I) Count the interactions                                               */

   intact_save      = 0;
   for(ipart=1;ipart<=natm_tot;ipart++){
    intact_save+=nbr_list_nter[ipart];
   }
   num_call_save    = intact_save/nlen;
   if( (intact_save % nlen)!=0 ){num_call_save += 1;}

/*==========================================================================*/
/* II) Initialize                                                           */

   iver_fill       = 1;
   iver_count      = 0;
   iver_init       = 0;
   nbr_list->verlist.iver_init  = 0;
   nbr_list->verlist.iver_fill  = 1;
   nbr_list->verlist.iver_count = 0;
   nver_now            = 0;
   nver_res_now        = 0;
   (*for_scr_intact)   = 0;
   (*for_scr_num_brk)  = 0;
   num_call            = 0;
   intact_tot          = 0;
   iswitch             = 0;
   lst_typ             = 1;
   nbr_list->brnch_root_list_opt = 0; 

/*==========================================================================*/
/* III) Loop over the possible interactions:                                */

    for(ipart=1;ipart<= natm_tot;ipart++){

/*------------------------------*/
/*  0) Get an  interaction sets */
      for(jpart=1;jpart<=nbr_list_nter[ipart];jpart++){
        for_scr_iexcl[jpart] = nbr_list_jter[jpart+nbr_list_jver_off[ipart]];
      }/*endfor*/
      intact_now = nbr_list_nter[ipart];
/*---------------------*/
/*  B) Set the padding */
      nbr_list_nter[ipart] = nbr_list_jver_off[ipart];
      if(int_res_ter==1){
        nbr_list_nter_res[ipart]=nbr_list_jver_off_res[ipart];
      }

/*==========================================================================*/
/*  B) Add interactions to the force routine interaction list               */
/*     in chunks no greater than nlen. Keep track of the break points       */

      lower  = 1;
      upper  = 0;
      while(upper!=intact_now){
        upper = intact_now;
        if(upper-lower+1+(*for_scr_intact)>(nlen)){
           upper = nlen - (*for_scr_intact) + lower-1;
        }/*endif*/
        jind_off = (*for_scr_intact)-lower+1;
        for(jpart=lower;jpart<=upper;jpart++){
           mtemp = jpart+jind_off;
           for_scr_j_index[mtemp]=for_scr_iexcl[jpart];
           for_scr_i_index[mtemp]=ipart;
        }/*endfor*/
        (*for_scr_intact) += upper-lower+1;
        (*for_scr_num_brk)++;
        for_scr_num_brk_i[(*for_scr_num_brk)] = upper-lower+1;
        lower = upper+1;
      
/*==========================================================================*/
/*  C) If enough interactions have accumulated (or you are done)            */
/*       make a verlist create call and then update and reinitialize        */
/*       appropriate counters                                               */

/*-------------------------------------------------------------------------*/
/*  i)  List creat conditions                                              */

        ifor_call = 0;
        if(((*for_scr_intact)==(nlen))){ifor_call=1;}
        if((ipart==(natm_tot))&&(upper==intact_now)&&
          ((*for_scr_intact)!=0)   ){
            ifor_call=1;
        }/*endif*/

/*-------------------------------------------------------------------------*/
/*  ii)  List create call and memory checks                                */

        if(ifor_call==1){

          verlist_create(clatoms_info,clatoms_pos,
                         intra_scr,for_scr,atommaps,cell,
                         nbr_list,timeinfo,interact,iswitch,
                         &nver_now,&nver_res_now,lst_typ);
          num_call        += 1;
          intact_tot      += (*for_scr_intact);
          (*for_scr_intact)  = 0;
          (*for_scr_num_brk) = 0;

          ver_list_mem_check(nver_now,&(nbr_list->verlist.nver_lst),
                             iver_fill,nmem_min_lst,iver_init,mem_safe,
                             &irealloc,myid_forc);
          if(irealloc==1){
            nbr_list->verlist.jter = (list_int *)
                   crealloc(&(nbr_list->verlist.jter)[1],
                         (nbr_list->verlist.nver_lst)*sizeof(list_int))-1;
          }/*endif*/

          if(int_res_ter==1){
            ver_list_mem_check(nver_res_now,&(nbr_list->verlist.nver_lst_res),
                               iver_fill,nmem_min_lst,iver_init,mem_safe,
                               &irealloc,myid_forc);
            if(irealloc==1){
              nbr_list->verlist.jter_res = (list_int *)
                    crealloc(&(nbr_list->verlist.jter_res)[1],
                      (nbr_list->verlist.nver_lst_res)*sizeof(list_int))-1;
            }/*endif*/
          }/*endif:RESPA*/

        }/*endif:for_call*/

      }/*endwhile:getting interactions*/

   }/*endfor:each particles interactions*/

/*==========================================================================*/
/* IV) Final list create call and memory check */

  if((*for_scr_intact)>0){
    verlist_create(clatoms_info,clatoms_pos,
                   intra_scr,for_scr,atommaps,cell,
                   nbr_list,timeinfo,interact,iswitch,
                   &nver_now,&nver_res_now,lst_typ);
    num_call          += 1;
    intact_tot        += (*for_scr_intact);

    ver_list_mem_check(nver_now,&(nbr_list->verlist.nver_lst),
                       iver_fill,nmem_min_lst,iver_init,mem_safe,
                       &irealloc,myid_forc);
    if(irealloc==1){
         nbr_list->verlist.jter = (list_int *)
                crealloc(&(nbr_list->verlist.jter)[1],
                     (nbr_list->verlist.nver_lst)*sizeof(list_int))-1;
    }/*endif*/

    if(int_res_ter==1){
      ver_list_mem_check(nver_res_now,&(nbr_list->verlist.nver_lst_res),
                         iver_fill,nmem_min_lst,iver_init,mem_safe,
                         &irealloc,myid_forc);
      if(irealloc==1){
         nbr_list->verlist.jter_res = (list_int *)
                crealloc(&(nbr_list->verlist.jter_res)[1],
                   (nbr_list->verlist.nver_lst_res)*sizeof(list_int))-1;
      }/*endif*/
    }/*endif:RESPA*/

  }/*endif:last force call*/

/*==========================================================================*/
/* V) Depad nter                                                           */

   for(jpart=1;jpart<=natm_tot;jpart++){
      nbr_list_nter[jpart] -= nbr_list_jver_off[jpart];
   }/*endfor*/

   if((int_res_ter==1)){
     for(jpart=1;jpart<=natm_tot;jpart++){
        nbr_list_nter_res[jpart] -= nbr_list_jver_off_res[jpart];
     }/*endfor*/
   }/*endif*/


/*==========================================================================*/
/* VI) Assign present size of list to appropriate variable                   */

   (*nver) = nver_now;
   (*nver_res) = nver_res_now;

/*==========================================================================*/
/* VI) Check                                                                */

    if(intact_tot!=intact_save){
       printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
       printf("Internal Error:\n ");
       printf("Incorrect number of interactions calculated\n ");
       printf("   in routine ver_pare_down_root \n ");
       printf("%d vs %d \n",intact_tot,intact_save);
       printf("Contact product support \n ");
       printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
       fflush(stdout);
       exit(1);
   }/*endif*/
   if(num_call!=num_call_save){
       printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
       printf("Internal Error:\n ");
       printf("Incorrect number of force calls\n ");
       printf("   in routine ver_pare_down_root \n ");
       printf("%d vs %d \n",num_call,num_call_save);
       printf("Contact product support \n ");
       printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
       fflush(stdout);
       exit(1);
   }/*endif*/

/*==========================================================================*/
/* VI) Reset list option                                                    */

   nbr_list->brnch_root_list_opt = 2; 

/*--------------------------------------------------------------------------*/
    }/*end routine */
/*==========================================================================*/



