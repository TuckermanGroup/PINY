/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
/*                                                                          */
/*                         PI_MD:                                           */
/*             The future of simulation technology                          */
/*             ------------------------------------                         */
/*                     Module: mall_make_lists.c                            */
/*                                                                          */
/* This subprogram allocates the memory for the neighbor lists and then     */
/* creates them                                                             */
/*                                                                          */
/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

#include "standard_include.h"
#include "../typ_defs/typedefs_par.h"
#include "../typ_defs/typedefs_class.h"
#include "../typ_defs/typedefs_gen.h"
#include "../typ_defs/typedefs_bnd.h"
#include "../proto_defs/proto_lists_entry.h"
#include "../proto_defs/proto_lists_local.h"
#include "../proto_defs/proto_real_space_local.h"
#include "../proto_defs/proto_friend_lib_entry.h"
#include "../proto_defs/proto_communicate_wrappers.h"

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void mall_make_lists(CLASS *class,GENERAL_DATA *general_data,
                                    BONDED *bonded,int error_check_on)

/*==========================================================================*/
{/*begin routine */

int isave,iii;
int np         = class->communicate.np;
int rank       = class->communicate.myid;
MPI_Comm world = class->communicate.world;

/*==========================================================================*/
/* 0) Initialize */

if(rank==0){
  printf("\n");PRINT_LINE_STAR;
  printf("Setting up the neighbor lists\n");
  PRINT_LINE_DASH;printf("\n");
}/*endif*/
if(np>1){Barrier(world);}

/*==========================================================================*/
/* I) Malloc and create the lnk lists: Must be done first.                  */
/*                                     RESPA lnk lst shifts not             */
/*                                     needed by lnk_ver_update option      */

 if( (class->nbr_list.ilnk==1)||
    ((class->nbr_list.verlist.lnk_ver_update+class->nbr_list.iver)==2)){
   isave = general_data->timeinfo.int_res_ter;    
   if(class->nbr_list.iver==1){general_data->timeinfo.int_res_ter=0;}
   if(rank==0){
    printf("Setting up a link list\n"); 
   }/*endif*/
   if(np>1){Barrier(world);}

   lnk_mall_make(&(class->clatoms_info),&(class->clatoms_pos[1]),
                 &(class->nbr_list),&(bonded->excl),
                 &(class->atommaps),&(general_data->cell),
                 &(bonded->intra_scr),
                 &(class->for_scr),&(general_data->timeinfo),
                 &(class->interact),&(class->tot_memory),rank,
                   error_check_on,world,np);

   if(class->nbr_list.iver==1&&rank==0){printf("\n");}
   if(class->nbr_list.iver==1){general_data->timeinfo.int_res_ter=isave;}
 /*endif*/}

/*==========================================================================*/
/* II) Malloc and create the ver lists */

 if(class->nbr_list.iver==1){
   (general_data->stat_avg.itime_update) = 0;
  if(rank==0){
    printf("Setting up a Verlet list\n"); 
  }/*endif*/
  if(np>1){Barrier(world);}
   ver_mall_make(&(class->clatoms_info),&(class->clatoms_pos[1]),
                 &(class->nbr_list),&(bonded->excl),
                 &(class->atommaps),&(general_data->cell),
                 &(bonded->intra_scr),
                 &(class->for_scr),&(general_data->timeinfo),
                 &(class->interact),&(class->tot_memory),rank,
                 error_check_on,world,&(class->class_comm_forc_pkg),np);
 /*endif*/}

/*==========================================================================*/
/* III) Done */

if(rank==0){
  printf("\n");
  PRINT_LINE_DASH;
  printf("Completed neighbor list set up\n");
  PRINT_LINE_STAR;printf("\n");
}/*endif*/
if(np>1){Barrier(world);}

/*--------------------------------------------------------------------------*/
/*end routine*/}
/*==========================================================================*/


/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void ver_mall_make(CLATOMS_INFO *clatoms_info,CLATOMS_POS *clatoms_pos,
                   NBR_LIST *nbr_list,EXCL *excl,
                   ATOMMAPS *atommaps,CELL *cell,INTRA_SCR *intra_scr,
                   FOR_SCR *for_scr,TIMEINFO *timeinfo,
                   INTERACT *interact, double *tot_memory,int rank,
                   int error_check_on,MPI_Comm world,
                   CLASS_COMM_FORC_PKG *class_comm_forc_pkg,int np)

/*==========================================================================*/
{/*begin routine */
/*==========================================================================*/
/*                     Local variable definitions                           */

#include "../typ_defs/typ_mask.h"

int i,natm_mall,nmem_int,nmem_dbl,iii,iflag;
int nver,nver_res,natm_tot;
int itest;
int npairs,npairs_tmp,np_tot;
double now_memory;
int brnch_root_list_opt   = nbr_list->brnch_root_list_opt; 
MPI_Comm comm_forc = class_comm_forc_pkg->comm;

/*==========================================================================*/
/* 0) Output to Screen                                                      */

/*==========================================================================*/
/* I) Malloc the stuff of length natm_tot: Must be done first              */

    natm_tot  = clatoms_info->natm_tot;
    natm_mall = clatoms_info->natm_tot;
    if((natm_mall           % 2)==0){natm_mall      += 1;}
    nbr_list->verlist.nter     = (int *)cmalloc(natm_mall*sizeof(int))-1;
    nbr_list->verlist.jver_off = (int *)cmalloc(natm_mall*sizeof(int))-1;
    nbr_list->x0       = (double *)cmalloc(natm_mall*sizeof(double))-1;
    nbr_list->y0       = (double *)cmalloc(natm_mall*sizeof(double))-1;
    nbr_list->z0       = (double *)cmalloc(natm_mall*sizeof(double))-1;
    nmem_int = 2*natm_mall;
    nmem_dbl = 3*natm_mall;
    if(timeinfo->int_res_ter==1){
       nmem_int += 2*natm_mall;
       nbr_list->verlist.nter_res  = (int *)cmalloc(natm_mall*sizeof(int))-1; 
       nbr_list->verlist.jver_off_res=(int *)cmalloc(natm_mall*sizeof(int))-1;
    /*endif*/}
    for(i=1;i<=natm_tot;i++){
      nbr_list->x0[i] = clatoms_pos->x[i];
      nbr_list->y0[i] = clatoms_pos->y[i];
      nbr_list->z0[i] = clatoms_pos->z[i];
    }/*endfor*/


/*==========================================================================*/
/* II) Construct Verlet list  */

   nbr_list->verlist.iver_init = 1;

/*--------------------------------------------------------------------------*/
/*  A) No list option */

 if((nbr_list->verlist.nolst_ver_update) == 1){

/* i) Count      */

    nbr_list->verlist.iver_fill=0;nbr_list->verlist.iver_count=1;       
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

/* ii) Malloc    */                                                        

    nbr_list->verlist.nver_lst  = (nver)*nbr_list->verlist.mem_safe;
    itest = nbr_list->verlist.nver_lst
          - (nver+2*nbr_list->verlist.nmem_min_lst);
    if(itest<0){nbr_list->verlist.nver_lst-=itest;}
    if((nbr_list->verlist.nver_lst % 2)==0){nbr_list->verlist.nver_lst+=1;}
    nbr_list->verlist.jter = 
        (list_int *) cmalloc(nbr_list->verlist.nver_lst*sizeof(list_int))-1;

    if(timeinfo->int_res_ter==1){
      nbr_list->verlist.nver_lst_res = (nver_res)*nbr_list->verlist.mem_safe;
      itest = nbr_list->verlist.nver_lst_res
            - (nver+2*nbr_list->verlist.nmem_min_lst);
      if(itest<0){nbr_list->verlist.nver_lst_res-=itest;}
      if((nbr_list->verlist.nver_lst_res % 2)==0){ 
                         nbr_list->verlist.nver_lst_res+=1;}
      nbr_list->verlist.jter_res=(list_int *) cmalloc(
                 nbr_list->verlist.nver_lst_res*sizeof(list_int))-1;
    }/*endif*/

/* iii) Fill     */

    nbr_list->verlist.iver_fill=1;nbr_list->verlist.iver_count=1;       
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
    nbr_list->verlist.nver_lst_now = nver;
    nbr_list->verlist.nver_lst_now_res = nver_res;


  }/*endif*/



/*--------------------------------------------------------------------------*/
/*  B) Lnk list option */

 if(nbr_list->verlist.lnk_ver_update == 1){

/* i) Count      */

    nbr_list->verlist.iver_fill=0;nbr_list->verlist.iver_count=1;       
    if(brnch_root_list_opt==0){
      lnk_ver_gen(clatoms_info,clatoms_pos,
                  for_scr,atommaps,cell,interact,
                  timeinfo,nbr_list,excl,intra_scr,&nver,&nver_res,
                  error_check_on,class_comm_forc_pkg);
    }else{

      lnk_ver_gen_root(clatoms_info,clatoms_pos,
                       for_scr,atommaps,cell,interact,
                       timeinfo,nbr_list,excl,intra_scr,&nver,&nver_res,
                       error_check_on,class_comm_forc_pkg);
    }/*endif*/

/* ii) Malloc    */

    nbr_list->verlist.nver_lst  = nver*nbr_list->verlist.mem_safe;
    itest = nbr_list->verlist.nver_lst
          - (nbr_list->verlist.jver_pad*(natm_tot-1)
          + nver+2*nbr_list->verlist.nmem_min_lst);
    if(itest<0){nbr_list->verlist.nver_lst-=itest;}
    if((nbr_list->verlist.nver_lst % 2)==0){nbr_list->verlist.nver_lst+=1;}
    nbr_list->verlist.jter =
          (list_int *) cmalloc(nbr_list->verlist.nver_lst*sizeof(list_int))-1;

    if(timeinfo->int_res_ter==1){
      nbr_list->verlist.nver_lst_res  = nver_res*nbr_list->verlist.mem_safe;
      itest = nbr_list->verlist.nver_lst_res
        - (nbr_list->verlist.jver_pad*(natm_tot-1)
        + nver_res+2*nbr_list->verlist.nmem_min_lst);
      if(itest<0){nbr_list->verlist.nver_lst_res-=itest;}
      if((nbr_list->verlist.nver_lst_res % 2)==0){
                                       nbr_list->verlist.nver_lst_res+=1;}
      nbr_list->verlist.jter_res=
     (list_int *) cmalloc(nbr_list->verlist.nver_lst_res*sizeof(list_int))-1;
    }/*endif*/

/* iii) Fill     */
    nbr_list->verlist.iver_fill=1;nbr_list->verlist.iver_count=0;   
    if(brnch_root_list_opt==0){
      lnk_ver_gen(clatoms_info,clatoms_pos,
                  for_scr,atommaps,cell,interact,
                  timeinfo,nbr_list,excl,intra_scr,&nver,&nver_res,
                  error_check_on,class_comm_forc_pkg);
    }else{
      lnk_ver_gen_root(clatoms_info,clatoms_pos,
                       for_scr,atommaps,cell,interact,
                       timeinfo,nbr_list,excl,intra_scr,&nver,&nver_res,
                       error_check_on,class_comm_forc_pkg);
    }/*endif*/
    nbr_list->verlist.nver_lst_now = nver;
    nbr_list->verlist.nver_lst_now_res = nver_res;
 }/*endif*/

/*--------------------------------------------------------------------------*/
/* C) increment memory and set intialization flag to 0                      */

  nbr_list->verlist.iver_init = 0;
  nmem_int += nbr_list->verlist.nver_lst;
  if(timeinfo->int_res_ter==1){nmem_int += nbr_list->verlist.nver_lst_res;}
/*==========================================================================*/
/* VI)  More output to screen                                               */

  now_memory   = (sizeof(list_int)*nmem_int +sizeof(double)*nmem_dbl)*1.e-06;
  *tot_memory += now_memory;

  if(rank==0){
    printf("Verlet list allocation: %g Mbytes; Total memory: %g Mbytes\n",
            now_memory,*tot_memory);
  }/*endif*/
  if(np>1){Barrier(world);}

  npairs = nbr_list->verlist.nver_lst_now;
  np_tot = (clatoms_info->natm_tot)*(clatoms_info->natm_tot-1)/2 
         - (excl->nlst);


  if(class_comm_forc_pkg->num_proc>1){
    Reduce(&npairs,&npairs_tmp,1,MPI_INT,MPI_SUM,0,comm_forc);
  }else{
    npairs_tmp = npairs;
  }/*endif*/

  if(rank==0){
    printf("Number of pairs    = %d out of %d \n",npairs_tmp,np_tot);
  }/*endif*/
  if(np>1){Barrier(world);}

  if(timeinfo->int_res_ter==1){
    npairs = nbr_list->verlist.nver_lst_now_res;
    if(class_comm_forc_pkg->num_proc>1){
      Reduce(&npairs,&npairs_tmp,1,MPI_INT,MPI_SUM,0,comm_forc);
    }else{
      npairs_tmp = npairs;
    }/*endif*/
    if(rank==0){
      printf("Number RESPA pairs = %d out of %d \n",npairs_tmp,np_tot);
    }/*endif*/
    if(np>1){Barrier(world);}
  }/*endif*/

/*--------------------------------------------------------------------------*/
/*end routine*/}
/*==========================================================================*/


/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void lnk_mall_make(CLATOMS_INFO *clatoms_info,CLATOMS_POS *clatoms_pos,
                   NBR_LIST *nbr_list,EXCL *excl,
                   ATOMMAPS *atommaps,CELL *cell,INTRA_SCR *intra_scr,
                   FOR_SCR *for_scr,TIMEINFO *timeinfo,
                   INTERACT *interact, double *tot_memory,int rank,
                   int error_check_on,MPI_Comm world, int np)
/*==========================================================================*/
{/*begin routine */
/*==========================================================================*/
/*                     Local variable definitions                           */

int nmem_int,nmem_dbl,iii;
double now_memory;
double rexl_max_now;
int ncell_tot,nsh_tot;

/*==========================================================================*/
/* I) Get maximum cutoff distance between particles                         */
/*    and if we are using lnkcell to update the verlist add the skin        */

     max_cut_part(interact,nbr_list);
     if(nbr_list->iver==1){
        nbr_list->lnklist.rcut_max_res += (interact->skin);
        nbr_list->lnklist.rcut_max     += (interact->skin);
     /*endif*/}
     if(rank==0){
     printf("Using maximum cutoff %gA to construct the lnk shifts\n",
             nbr_list->lnklist.rcut_max*BOHR);
     }/*endif*/
     if(np>1){Barrier(world);}
     if(timeinfo->int_res_ter==1){
      if(rank==0){
        printf("Using maximum cutoff %gA to construct the RESPA lnk shifts\n",
               nbr_list->lnklist.rcut_max_res*BOHR);
      }/*endif*/
      if(np>1){Barrier(world);}
     }/*endif*/

/*==========================================================================*/
/* II) Get maximum separation between exclusions, add 1au for saftey        */

    nbr_list->lnklist.rexl_max = 0.0;
    rexl_max_now = 0.0;
    if((excl->nlst)!=0){
      max_cut_excl(clatoms_info,clatoms_pos,cell,nbr_list,excl);
     rexl_max_now = nbr_list->lnklist.rexl_max*BOHR;
     nbr_list->lnklist.rexl_max += nbr_list->lnklist.lnk_excl_safe;
    }/*endif*/
    if(rank==0){
     printf("Using maximum exclusion distance %g A",rexl_max_now);
     printf("to establish excl checking\n");
    }/*endif*/
    if(np>1){Barrier(world);}
  
/*==========================================================================*/
/* III) Make the list                                                       */

    nbr_list->lnklist.hmat_lnk  = 
                           (double *)cmalloc((size_t)9*sizeof(double))-1;
    nbr_list->lnklist.hmati_lnk = 
                           (double *)cmalloc((size_t)9*sizeof(double))-1;
    nbr_list->lnklist.ilnk_init = 1;  

    make_lnk_lst(clatoms_info,clatoms_pos,nbr_list,cell,for_scr,intra_scr,
                 timeinfo,excl,interact,error_check_on);

    nbr_list->lnklist.ilnk_init = 0;  

/*==========================================================================*/
/* VI)  More output to screen                                               */

  ncell_tot = (nbr_list->lnklist.ncell_a)*(nbr_list->lnklist.ncell_b)*
                                (nbr_list->lnklist.ncell_c);
  nmem_int  = ncell_tot + nbr_list->lnklist.lnk_list_dim 
                           +  nbr_list->lnklist.nshft_lnk*4;
  nmem_dbl  = nbr_list->lnklist.nshft_lnk;
  if(timeinfo->int_res_ter==1){
    nmem_int  += nbr_list->lnklist.nshft_lnk_res*4;
    nmem_dbl  += nbr_list->lnklist.nshft_lnk_res;
  }/*endif*/
  now_memory   = (sizeof(int)*0 + sizeof(double)*2*9 + 
                  sizeof(list_int)*nmem_int+sizeof(double )*nmem_dbl )*1.e-06;
  *tot_memory  += now_memory;
  if(rank==0){
   printf("Link list allocation: %g Mbytes; Total memory: %g Mbytes\n",
          now_memory,*tot_memory);
  }/*endif*/
  if(np>1){Barrier(world);}

  nsh_tot = ((nbr_list->lnklist.ncell_a)*
             (nbr_list->lnklist.ncell_b)*
             (nbr_list->lnklist.ncell_c)-1)/2 +
           (((nbr_list->lnklist.ncell_a)*
             (nbr_list->lnklist.ncell_b)*
             (nbr_list->lnklist.ncell_c)-1) % 2) + 1;
  if(rank==0){
   printf("Shape: ncell_a = %d, ncell_b = %d, ncell_c = %d \n",
                 nbr_list->lnklist.ncell_a,
                 nbr_list->lnklist.ncell_b,
                 nbr_list->lnklist.ncell_c);
   printf("# of shifts = %d out of %d\n",
           nbr_list->lnklist.nshft_lnk,nsh_tot);
  }/*endif*/
  if(np>1){Barrier(world);}


/*--------------------------------------------------------------------------*/
}/*end routine*/
/*==========================================================================*/





