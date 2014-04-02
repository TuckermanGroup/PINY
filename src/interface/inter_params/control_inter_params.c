/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
/*                                                                          */
/*                         PI_MD:                                           */
/*             The future of simulation technology                          */
/*             ------------------------------------                         */
/*                     Module: control_inter_params                         */
/*                                                                          */
/* This reads in and sets up the electron-atom interaction pseudopotential  */ 
/*                                                                          */
/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

#include "standard_include.h"
#include "../typ_defs/typedefs_class.h"
#include "../typ_defs/typedefs_par.h"
#include "../typ_defs/typedefs_bnd.h"
#include "../proto_defs/proto_inter_params_entry.h"
#include "../proto_defs/proto_inter_params_local.h"
#include "../proto_defs/proto_intra_params_local.h"
#include "../proto_defs/proto_search_entry.h"
#include "../proto_defs/proto_friend_lib_entry.h"
#include "../proto_defs/proto_handle_entry.h"
#include "../proto_defs/proto_communicate_wrappers.h"

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void control_inter_params(INTERACT *interact,SPLINE_PARSE *spline_parse,
                          FILENAME_PARSE *filename_parse,
                          double alp_ewd,int ncharge,
                          int natm_tot,int natm_typ,
                          NAME atm_typ[],int iatm_typ[],
                          int iperd,int ishift_pot,double *tot_memory,
                          int int_res_ter,int myid, MPI_Comm comm,
                          int num_proc )

/*======================================================================*/
/*  Begin routine */
     {/*begin routine*/
/*======================================================================*/
/*          Local variable declarations                                */
#include "../typ_defs/typ_mask.h"
  
  DATA_BASE_ENTRIES *inter_base;          /* Lst: Database parameters    */
  double *eps,*sig;                      /* Lst: Lennard-Jones params   */
  double *awill,*bwill,*c6m,*c8m,*c10m;  /* Lst: Williams params        */
  double *cwill ,*rm_swit, *c9m;         /* Lst: Aziz-chen params       */
  double *temp_cutoff,*temp_cutoff_res,*temp_cutti;
  int *inter_label;
  int *ifound,*isearch,*igood;           /* Lst: found,search goodness flags*/
  int i,j,iii;                           /* Num: Counters               */
  int ninter;                            /* Num: Number of interactions */
  double now_mem;                        /* Num: Memory allocated here  */
  char typ[6];
  int nbase,nbase2,ibase_want;
  CATM_LAB *cinter,*cinter_base;
  char *fun_key;
  DICT_WORD *fun_dict;
  int nsearch,natm_srch;
  int num_fun_dict;
  int ifirst,ityp;
  int nsplin_mall;
  int nsplin_mall_tot;
  int ninter_mall;

  int ninter_unique;                     /* Num: number of interactions with 
                                             unqiue paramters */
  int ninter_unique_mall;
 
/*======================================================================*/
/* 0) Write to screen                                                   */
  
  if(myid==0){
    ninter = natm_typ*(natm_typ + 1)/2;
    putchar('\n');
    PRINT_LINE_STAR;
    printf("Searching the data bases (both user defined and default)\n");
    printf("for the %d intermolecular interaction sets\n",ninter);
    PRINT_LINE_DASH;printf("\n");
  }/*endif*/

/*======================================================================*/
/*  I) Malloc the memory                                                 */

  ninter = natm_typ*(natm_typ + 1)/2;
  ninter_mall = ninter;
  if((ninter_mall!=0)&&((ninter_mall %2)==0)){ninter_mall +=1;}
  inter_label = (int *) cmalloc(ninter_mall*sizeof(int))-1;
  eps        = (double *) cmalloc(ninter_mall*sizeof(double))-1;
  sig        = (double *) cmalloc(ninter_mall*sizeof(double))-1;
  awill      = (double *) cmalloc(ninter_mall*sizeof(double))-1;
  bwill      = (double *) cmalloc(ninter_mall*sizeof(double))-1;
  cwill      = (double *) cmalloc(ninter_mall*sizeof(double))-1;
  rm_swit    = (double *) cmalloc(ninter_mall*sizeof(double))-1;
  c6m        = (double *) cmalloc(ninter_mall*sizeof(double))-1;
  c8m        = (double *) cmalloc(ninter_mall*sizeof(double))-1;
  c9m        = (double *) cmalloc(ninter_mall*sizeof(double))-1;
  c10m       = (double *) cmalloc(ninter_mall*sizeof(double))-1;
  fun_key    = (char *)cmalloc(MAXWORD*sizeof(char));  
  cinter     = (CATM_LAB *)cmalloc(ninter*sizeof(CATM_LAB))-1;  
  ifound     = (int *)cmalloc(ninter*sizeof(int))-1;
  isearch    = (int *)cmalloc(ninter*sizeof(int))-1;
  igood      = (int *)cmalloc(ninter*sizeof(int))-1;

  temp_cutoff     = (double *) cmalloc(ninter_mall*sizeof(double))-1;
  temp_cutoff_res = (double *) cmalloc(ninter_mall*sizeof(double))-1;
  temp_cutti      = (double *) cmalloc(ninter_mall*sizeof(double))-1;

  interact->cutoff     = (double *) cmalloc(ninter_mall*sizeof(double))-1;
  interact->cutoff_res = (double *) cmalloc(ninter_mall*sizeof(double))-1;
  interact->cutti      = (double *) cmalloc(ninter_mall*sizeof(double))-1;

  interact->cutskin    = (double *) cmalloc(ninter_mall*sizeof(double))-1;
  interact->cutskin_res= (double *) cmalloc(ninter_mall*sizeof(double))-1;
  interact->cutskin_root   = (double *) cmalloc(ninter_mall*sizeof(double))-1;
  interact->cutskin_root_res= (double *) cmalloc(ninter_mall*sizeof(double))-1;

  interact->inter_map_index = (int *) cmalloc(ninter_mall*sizeof(int))-1;

/*======================================================================*/
/*  II) Set up the data structures                                      */

if(myid==0){
  ifirst =1;
  set_potfun_dict(&fun_dict,&num_fun_dict,ifirst);
  ityp = 0;
  for(i=1;i <= natm_typ; i++) {
    for(j=i;j <= natm_typ; j++) {
      ityp++;
      strcpy(cinter[ityp].atm1,atm_typ[i]);
      strcpy(cinter[ityp].atm2,atm_typ[j]);
    }/*endfor*/
  }/*endfor*/
  for(i=1;i<=ninter;i++){ifound[i]=0;}
  for(i=1;i<=ninter;i++){igood[i]=6;}
}/*endif*/

/*======================================================================*/
/*  III) Search the user defined data base                              */

if(myid==0){
  natm_srch = 2;
  if(strlen(filename_parse->user_inter_name) != 0){
    nsearch = 1;
    ibase_want = 1;
    count_data_base(filename_parse->user_inter_name,fun_dict,num_fun_dict,
                    &nbase,ibase_want);
    if(nbase>0){
      nbase2 = 2*nbase;
      inter_base  = (DATA_BASE_ENTRIES *)
                       cmalloc(nbase2*sizeof(DATA_BASE_ENTRIES))-1;
      cinter_base = (CATM_LAB *)cmalloc(nbase2*sizeof(CATM_LAB))-1;
      read_data_base(filename_parse->user_inter_name,fun_dict,num_fun_dict,
                     inter_base,cinter_base,ibase_want,nbase);
      search_base(nbase,nbase2,cinter_base,ninter,cinter,igood,ifound,
                  isearch,nsearch,natm_srch,filename_parse->user_inter_name);

      assign_base_inter(inter_base,nbase,ifound,ninter,
                        sig,eps,awill,bwill,cwill,rm_swit,c6m,c8m,c9m,c10m,
                        inter_label,interact->cutoff,interact->cutoff_res,
                        interact->cutti,
                        isearch,nsearch,cinter,cinter_base);
      cfree(&inter_base[1]);
      cfree(&cinter_base[1]);
    }/*endif*/
  }/*endif*/
}/*endif*/
/*======================================================================*/
/*  IV) Search the default defined data base                            */

if(myid==0){
  if(strlen(filename_parse->def_inter_name) != 0){
    nsearch = 2;
    ibase_want = 1;
    count_data_base(filename_parse->def_inter_name,fun_dict,num_fun_dict,
                    &nbase,ibase_want);
    if(nbase>0){
      nbase2 = 2*nbase;
      inter_base = (DATA_BASE_ENTRIES *)
                    cmalloc(nbase*sizeof(DATA_BASE_ENTRIES))-1;
      cinter_base = (CATM_LAB *)cmalloc(nbase2*sizeof(CATM_LAB))-1;
      read_data_base(filename_parse->def_inter_name,fun_dict,num_fun_dict,
                     inter_base,cinter_base,ibase_want,nbase);
      search_base(nbase,nbase2,cinter_base,ninter,cinter,igood,ifound,
                  isearch,nsearch,natm_srch,filename_parse->def_inter_name);

      assign_base_inter(inter_base,nbase,ifound,ninter,
                        sig,eps,awill,bwill,cwill,rm_swit,c6m,c8m,c9m,c10m,
                        inter_label,interact->cutoff,interact->cutoff_res,
                        interact->cutti,
                        isearch,nsearch,cinter,cinter_base);
      cfree(&inter_base[1]);
      cfree(&cinter_base[1]);
    }/*endif*/
  }/*endif*/
}/*endif*/

/*======================================================================*/
/* V) Check list for missing entries                                    */

if(myid==0){
  strcpy(typ,"inter");
  atmlst_not_found(ninter,cinter,ifound,natm_srch,typ);
}/*endif*/

/*======================================================================*/
/*Find unique values of epsilon,sigma,rcut for LJ and null interactions */
/* This involves rearranging all the intermolecular interactions.       */

if(myid==0){

 for(i=1; i<= ninter; i++){
   temp_cutoff[i]     = interact->cutoff[i];
   temp_cutoff_res[i] = interact->cutoff_res[i];
   temp_cutti[i]      = interact->cutti[i];
 }/*endfor*/

#define BYPASS_OFF
#ifdef BYPASS

 ninter_unique = ninter;
 for(i=1;i<=ninter;i++){
   interact->inter_map_index[i] = i;
 }/*endfor*/

#else
 sort_inter_params(eps,sig,
                   awill,bwill,cwill,
                   rm_swit,
                   c6m,c8m,c9m,c10m,
                   temp_cutoff,temp_cutoff_res,temp_cutti,
                   inter_label,interact->inter_map_index,
                   &ninter_unique,ninter);
#endif

}/*endif myid*/

if(num_proc > 1){ Bcast(&ninter_unique,1,MPI_INT,0,comm);}

 ninter_unique_mall = ninter_unique;
 if((ninter_unique_mall!=0)&&((ninter_unique_mall %2)==0))
    {ninter_unique_mall++;}

/*======================================================================*/
/* V) Broadcast the parameters                                          */

 if(num_proc>1){
   Bcast(&(sig[1]),ninter_unique,MPI_DOUBLE,0,comm);
   Bcast(&(eps[1]),ninter_unique,MPI_DOUBLE,0,comm);
   Bcast(&(awill[1]),ninter_unique,MPI_DOUBLE,0,comm);
   Bcast(&(bwill[1]),ninter_unique,MPI_DOUBLE,0,comm);
   Bcast(&(cwill[1]),ninter_unique,MPI_DOUBLE,0,comm);
   Bcast(&(rm_swit[1]),ninter_unique,MPI_DOUBLE,0,comm);
   Bcast(&(c6m[1]),ninter_unique,MPI_DOUBLE,0,comm);
   Bcast(&(c8m[1]),ninter_unique,MPI_DOUBLE,0,comm);
   Bcast(&(c9m[1]),ninter_unique,MPI_DOUBLE,0,comm);
   Bcast(&(c10m[1]),ninter_unique,MPI_DOUBLE,0,comm);

   Bcast(&(inter_label[1]),ninter_unique,MPI_INT,0,comm);

   Bcast(&(temp_cutoff[1]),ninter_unique,MPI_DOUBLE,0,comm);
   Bcast(&(temp_cutoff_res[1]),ninter_unique,MPI_DOUBLE,0,comm);
   Bcast(&(temp_cutti[1]),ninter_unique,MPI_DOUBLE,0,comm);

   Bcast(&(interact->cutoff[1]),ninter,MPI_DOUBLE,0,comm);
   Bcast(&(interact->cutoff_res[1]),ninter,MPI_DOUBLE,0,comm);
   Bcast(&(interact->cutti[1]),ninter,MPI_DOUBLE,0,comm);
   Bcast(&(interact->inter_map_index[1]),ninter,MPI_INT,0,comm);
 }/*endif*/

/*======================================================================*/
/* VI) Allocate spline arrays                                           */

  interact->nter_typ = ninter;
  interact->nsplin_tot = interact->nsplin*ninter;
  now_mem = ((double)(interact->nsplin*ninter_unique*(sizeof(double)*4
                    +sizeof(int)*0)+
                    + sizeof(int)*0)  +
                  ninter_unique*(sizeof(double)*4 + sizeof(int)*0 )
           )*1.e-06;

  *tot_memory += now_mem;
  
  if(myid==0){
   printf("There are %d unique interactions \n",ninter_unique);
   printf("Intermolecular allocation: %g Mbytes; Total memory: %g Mbytes\n",
           now_mem,*tot_memory);
  }/*endif*/

  nsplin_mall_tot = interact->nsplin*ninter_unique;
  if((nsplin_mall_tot!=0)&&((nsplin_mall_tot % 2)==0)){nsplin_mall_tot += 1;}

  nsplin_mall = interact->nsplin;
  if((nsplin_mall!=0)&&((nsplin_mall % 2)==0)){nsplin_mall += 1;}

  interact->cv0    = (double *) cmalloc(nsplin_mall_tot*sizeof(double))-1;
  interact->cdv0   = (double *) cmalloc(nsplin_mall_tot*sizeof(double))-1;
  interact->cv0_c  = (double *) cmalloc(nsplin_mall_tot*sizeof(double))-1;
  interact->cdv0_c = (double *) cmalloc(nsplin_mall_tot*sizeof(double))-1;

  interact->vcut_coul  = (double *) cmalloc(ninter_unique_mall*sizeof(double))-1;
  interact->rmin_spl   = (double *) cmalloc(ninter_unique_mall*sizeof(double))-1;
  interact->dr_spl     = (double *) cmalloc(ninter_unique_mall*sizeof(double))-1;
  interact->dri_spl    = (double *) cmalloc(ninter_unique_mall*sizeof(double))-1;


/*======================================================================*/
/*          Assign mall variables                                 */

  interact->ninter_mall        = ninter_mall;
  interact->nsplin_mall_tot    = nsplin_mall_tot;

/*=======================================================================*/
/* VII) Set up the splines for the real-space intermolecular potential   */
/*     energy and forces                                                 */

/* Spline the intermolecular interaction */

 set_inter_splin(sig,eps,awill,bwill,cwill,rm_swit,c6m,c8m,c9m,c10m,
                  alp_ewd,ishift_pot,inter_label,
                  spline_parse,interact,
                  temp_cutti,temp_cutoff,
                  ninter_unique,ncharge,iperd,myid);


/*=======================================================================*/
/*  VIII) Get the long range correction                                   */

  interact->clong = 0.0;
  interact->clong_res = 0.0;      
  if(iperd == 3) {

    get_clong(natm_tot,natm_typ,ninter,iatm_typ,c6m,
              interact->inter_map_index,
             &(interact->clong),
             &(interact->clong_res),temp_cutoff,temp_cutoff_res,
             interact->iswit_vdw,interact->rheal_res);

  } /* endif */
  
  if(myid==0){
   printf("Dispersion long range parameter %.15g \n",interact->clong);
   if(int_res_ter==1){
    printf("Dispersion long range parameter(RESPA) %g \n",
         interact->clong_res);
   }/*endif*/
  }/*endif*/

/*=======================================================================*/
/*  IX) Free temporary memory                                            */

  cfree(&inter_label[1]);
  cfree(&eps[1]);
  cfree(&sig[1]);
  cfree(&awill[1]);
  cfree(&bwill[1]);
  cfree(&cwill[1]);
  cfree(&rm_swit[1]);
  cfree(&c6m[1]);
  cfree(&c8m[1]);
  cfree(&c9m[1]);
  cfree(&c10m[1]);
  cfree(&temp_cutoff[1]);
  cfree(&temp_cutoff_res[1]);
  cfree(&temp_cutti[1]);

  if(myid==0){
    cfree(&fun_dict[1]);
  }/*endif*/
  cfree(fun_key);
  cfree(&cinter[1]);
  cfree(&ifound[1]);
  cfree(&isearch[1]);
  cfree(&igood[1]);

/*=======================================================================*/
/* X) Write to screen                                                    */

if(myid==0){
  printf("\n");
  PRINT_LINE_DASH;
  printf("Completed the data bases searches\n");
  PRINT_LINE_STAR;
  putchar('\n');
}/*endif*/

/*----------------------------------------------------------------------*/
} /* end routine */
/*==========================================================================*/



/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void inter_coef(DICT_WORD *dict,char filename[],char fun_key[],
               DATA_BASE_ENTRIES *inter_base,CATM_LAB *cinter_base,int ibase)

/*==========================================================================*/
/*               Begin subprogram:                                          */
      {/*begin routine*/
/*==========================================================================*/
/*               Local variable declarations:                               */

  int index,iii;
  double bohr2,bohr6,bohr8,bohr9,bohr10;
  double cutti,cutoff,cutoff_res;
  double sig,eps;
  double c6m,c8m,c9m,c10m,awill,bwill,cwill,rm_swit;

/*==========================================================================*/
/* I) Fill atom types and label part of the data base     */

  strcpy(cinter_base[ibase].atm1,dict[1].keyarg);
  strcpy(cinter_base[ibase].atm2,dict[2].keyarg);
  strcpy(cinter_base[ibase].label,"");

/*=======================================================================*/
/* II) Set up */

  inter_base[ibase].inter_label = -1;
  bohr2  = BOHR*BOHR;
  bohr6  = bohr2*bohr2*bohr2;
  bohr8  = bohr6*bohr2;
  bohr9  = bohr8*BOHR;
  bohr10 = bohr8*bohr2;

/*=======================================================================*/
/* III) Convert and assign cutoffs */
  
  sscanf(dict[4].keyarg,"%lf",&cutti);
  sscanf(dict[5].keyarg,"%lf",&cutoff);
  sscanf(dict[6].keyarg,"%lf",&cutoff_res);
  if(cutti<0){
    index=4;
    keyarg_barf(dict,filename,fun_key,index);
  }/*endif*/
  if(cutoff<0){
    index=5;
    keyarg_barf(dict,filename,fun_key,index);
  }/*endif*/
  if(cutoff_res<0){
    index=6;
    keyarg_barf(dict,filename,fun_key,index);
  }/*endif*/
  inter_base[ibase].cutti      = cutti/BOHR;
  inter_base[ibase].cutoff     = cutoff/BOHR;
  inter_base[ibase].cutoff_res = cutoff_res/BOHR;

/*=======================================================================*/
/* IV) Convert and assign Lennard-Jones                                  */
  
  if(strcasecmp(dict[3].keyarg,"lennard-Jones") == 0) {
    inter_base[ibase].inter_label = 2;
    sscanf(dict[7].keyarg,"%lg",&sig);
    sscanf(dict[8].keyarg,"%lg",&eps);
    if(sig<0){
      index=7;
      keyarg_barf(dict,filename,fun_key,index);
    }/*endif*/
    if(eps<0){
      index=8;
      keyarg_barf(dict,filename,fun_key,index);
    }/*endif*/
    eps /= BOLTZ;
    sig /= BOHR;
    inter_base[ibase].eps        = eps;
    inter_base[ibase].sig        = sig;
    inter_base[ibase].c6m        = (4.0*eps)*(sig*sig*sig*sig*sig*sig);
  }/*endif*/

/*=======================================================================*/
/*  3) Convert and assign Williams                                        */
  
  if(strcasecmp(dict[3].keyarg,"williams") == 0) {
    inter_base[ibase].inter_label = 3;
    sscanf(dict[9].keyarg,"%lg",&c6m);
    sscanf(dict[10].keyarg,"%lg",&c8m);
    sscanf(dict[11].keyarg,"%lg",&c10m);
    sscanf(dict[12].keyarg,"%lg",&awill);
    sscanf(dict[13].keyarg,"%lg",&bwill);
    if(awill<0){
      index=12;
      keyarg_barf(dict,filename,fun_key,index);
    }/*endif*/
    if(bwill<0){
      index=13;
      keyarg_barf(dict,filename,fun_key,index);
    }/*endif*/
    inter_base[ibase].awill      = awill/BOLTZ;
    inter_base[ibase].bwill      = bwill*BOHR;
    inter_base[ibase].c6m        = c6m/(BOLTZ*bohr6);
    inter_base[ibase].c8m        = c8m/(BOLTZ*bohr8);
    inter_base[ibase].c10m       = c10m/(BOLTZ*bohr10);
  }/*endif*/

/*=======================================================================*/
/*  5) Convert and assign Null                                           */
  
  if(strcasecmp(dict[3].keyarg,"null") == 0) {
    inter_base[ibase].inter_label = 4;
    inter_base[ibase].c6m = 0.0;
  }/*endif*/

/*=======================================================================*/
/*  4) Convert and assign Aziz-Chen                                      */
  
  if(strcasecmp(dict[3].keyarg,"aziz-chen") == 0) {
    inter_base[ibase].inter_label = 6;
    sscanf(dict[9].keyarg,"%lg",&c6m);
    sscanf(dict[10].keyarg,"%lg",&c8m);
    sscanf(dict[16].keyarg,"%lg",&c9m);
    sscanf(dict[11].keyarg,"%lg",&c10m);
    sscanf(dict[12].keyarg,"%lg",&awill);
    sscanf(dict[13].keyarg,"%lg",&bwill);
    sscanf(dict[14].keyarg,"%lg",&cwill);
    sscanf(dict[15].keyarg,"%lg",&rm_swit);
    if(awill<0){
      index=12;
      keyarg_barf(dict,filename,fun_key,index);
    }/*endif*/
    if(bwill<0){
      index=13;
      keyarg_barf(dict,filename,fun_key,index);
    }/*endif*/
    if(cwill<0){
      index=14;
      keyarg_barf(dict,filename,fun_key,index);
    }/*endif*/
    if(rm_swit<0){
      index=15;
      keyarg_barf(dict,filename,fun_key,index);
    }/*endif*/
    inter_base[ibase].awill      = awill/BOLTZ;
    inter_base[ibase].bwill      = bwill*BOHR;
    inter_base[ibase].cwill      = cwill*bohr2;
    inter_base[ibase].rm_swit    = rm_swit/BOHR; 
    inter_base[ibase].c6m        = c6m/(BOLTZ*bohr6);
    inter_base[ibase].c8m        = c8m/(BOLTZ*bohr8);
    inter_base[ibase].c9m        = c9m/(BOLTZ*bohr9);
    inter_base[ibase].c10m       = c10m/(BOLTZ*bohr10);
  }/*endif*/

/*=======================================================================*/
/*  6) Check Potential type                                              */

  if(inter_base[ibase].inter_label==-1){
    index=3;
    keyarg_barf(dict,filename,fun_key,index);
  }/*endif*/


/*--------------------------------------------------------------------------*/
   }/*end routine*/
/*==========================================================================*/



/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
/* set_inter_splin:This subroutine fits the necessary terms to a spline     */
/*--------------------------------------------------------------------------*/

void assign_base_inter(DATA_BASE_ENTRIES *inter_base,int nbase,int *ifound,
                      int ninter,
                       double *sig, double *eps,double *awill,double *bwill,
                       double *cwill,double *rm_swit,double *c6m,
                       double *c8m,double *c9m,double *c10m,
                       int *inter_label,double *cutoff,
                       double *cutoff_res,
                       double *cutti,
                       int *isearch,int nsearch,CATM_LAB *cinter,
                       CATM_LAB *cinter_base)

/*=======================================================================*/
/*            Begin subprogram:                                          */
  {/*begin routine*/
/*=======================================================================*/
/*             Local variable declarations                                */
  int ibase,i,iii;
/*=======================================================================*/

  for(i=1;i<=ninter;i++){
    if(ifound[i] > 0 && isearch[i]==nsearch){
      ibase = ifound[i];

      eps[i]         = 0.0;
      sig[i]         = 0.0;
      awill[i]       = 0.0;
      bwill[i]       = 0.0;
      cwill[i]       = 0.0;
      rm_swit[i]     = 0.0;
      c8m[i]         = 0.0;  
      c9m[i]         = 0.0;  
      c10m[i]        = 0.0;  

      inter_label[i] = inter_base[ibase].inter_label;
      cutti[i]       = inter_base[ibase].cutti;
      cutoff[i]      = inter_base[ibase].cutoff;     
      cutoff_res[i]  = inter_base[ibase].cutoff_res;

      switch(inter_label[i]){
        case 2:  eps[i]     = inter_base[ibase].eps;
                 sig[i]     = inter_base[ibase].sig;
                 c6m[i]     = inter_base[ibase].c6m;  
               break;
        case 3:  awill[i]   = inter_base[ibase].awill;
                 bwill[i]   = inter_base[ibase].bwill;
                 c6m[i]     = inter_base[ibase].c6m;
                 c8m[i]     = inter_base[ibase].c8m;
                 c10m[i]    = inter_base[ibase].c10m; 
               break;
        case 4:  c6m[i]     =  0.0;  
               break;
        case 6:  awill[i]   = inter_base[ibase].awill;
                 bwill[i]   = inter_base[ibase].bwill;
                 cwill[i]   = inter_base[ibase].cwill;
                 rm_swit[i] = inter_base[ibase].rm_swit;
                 c6m[i]     = inter_base[ibase].c6m;
                 c8m[i]     = inter_base[ibase].c8m;
                 c9m[i]     = inter_base[ibase].c9m;
                 c10m[i]    = inter_base[ibase].c10m; 
               break;
      }/*endswitch*/
    }/*endif*/
  }/*endfor*/

/*--------------------------------------------------------------------------*/
   }/*end routine*/
/*==========================================================================*/
/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
/* set_inter_splin:This subroutine fits the necessary terms to a spline     */
/*--------------------------------------------------------------------------*/

void set_inter_splin(double sig[],double eps[],double a[],double b[],
                     double c[],double rm[],
                     double c6[],double c8[],double c9[],double c10[],
                     double alp_ewd,int ishift,int inter_label[],
                     SPLINE_PARSE *spline_parse,
                     INTERACT *interact,
                     double cutti[],double cutoff[],
                     int ninter,int ncharge,int iperd,
                     int myid)

/*=======================================================================*/
/*            Begin subprogram:                                          */
  {/*begin routine*/
/*=======================================================================*/
/*             Local variable declarations                                */

  double sigt,epst;                      /* Num: LJ potential params    */
  double at,bt,c6t,c8t,c9t,c10t;             /* Num: Williams pot params    */
  double ct,rmt;                         /* Num: Aziz-Chen pot params   */
  int ioff;                              /* Num: Offset for spline      */
  int i,iii;                             /* Num: For loop counter       */
  int itype;                             /* Num: Potential type label   */
  int iiperd;                            /* Num: Periodicity            */
  int ishift_now;
  double rmax_c,rmax;                    /* Num: Potential cutoff radii */
  double qijp;                           /* Num: product of charges     */
  double dr_tmp,rmin;
  double dri_tmp;

  double *rmin_spl  = interact->rmin_spl;
  double *vcut_coul = interact->vcut_coul;

/*======================================================================*/
/*  I) Spline each interaction                                          */

  ioff = 0;
  qijp = 0.0;
  rmax_c = 0.0;

  for(i=1;i <= ninter;i++) {

    itype       = inter_label[i];
    rmin_spl[i] = cutti[i];
    rmax        = cutoff[i];
    rmax_c      = MAX(rmax,rmax_c);
    iiperd      = iperd;

    spline_vdv(interact->rmin_spl[i],rmax,&(interact->cv0)[ioff],
             &(interact->cdv0)[ioff],interact->nsplin,&(interact->dr_spl)[i],
             &(interact->dri_spl)[i],sig[i],eps[i],a[i],b[i],c[i],rm[i],
             c6[i],c8[i],c9[i],c10[i],alp_ewd,qijp,iiperd,itype,ishift,
             interact->rheal_res,interact->dielectric_opt,
             interact->dielectric_rheal,interact->dielectric_cut,
             interact->dielectric_eps);

    ioff += interact->nsplin;
  } /* endfor:spline each interaction */

  interact->cutoff_max = rmax_c;

/*======================================================================*/
/*  II) Zero the coulomb shift                                          */

  if(ishift == 0) {
    for(i=1;i <= ninter;i++) {
      vcut_coul[i] = 0.0;
    }/* endfor */
  }/* endif: shift off */
  
/*======================================================================*/
/*  III) Spline the ewald coulomb stuff                                 */

  if (ncharge > 0 ){
   if( iperd >0 || interact->dielectric_opt==1 ){
    itype = 1;
    qijp = 1.0;
    iiperd = iperd;
    sigt = 1.0;
    epst = 0.0;
    at = 0.0;
    bt = 0.0;
    ct = 0.0;
    rmt = 0.0;
    c6t = 0.0;
    c8t = 0.0;
    c9t = 0.0;
    c10t = 0.0;
    if(iperd >0){
     if((alp_ewd*rmax_c)<3.4){
      if(myid==0){
       printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
       printf("Ewald convergence parameter too small !\n");
       printf("%g*rmax requested, at least 3.4*rmax required\n", 
             (alp_ewd*rmax_c));
       printf("If this criteria is not met, energy will drift.\n");
       printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
       fflush(stdout);
      }/*endif*/
      /* exit(1);*/
     }/*endif*/    
    }/*endif for iperd>0*/    
    ioff = 0;
    ishift_now = ishift;
    for(i=1;i <= ninter;i++) {
       rmin = cutti[i];
       rmax = cutoff[i];
       iiperd = iperd;
       spline_vdv(rmin,rmax,&(interact->cv0_c)[ioff],&(interact->cdv0_c)[ioff],
             interact->nsplin,&dr_tmp,&dri_tmp,sigt,epst,at,bt,ct,rmt,c6t,
             c8t,c9t,c10t,alp_ewd,qijp,iiperd,itype,ishift_now,
             interact->rheal_res,interact->dielectric_opt,
             interact->dielectric_rheal,interact->dielectric_cut,
             interact->dielectric_eps);
       ioff += interact->nsplin;
    }/*endfor*/
   } /* endif:spline special coulomb */
  }/*endif for ncharge*/

/*======================================================================*/
/*  V) Shift straight coulomb                                           */

 if(ncharge > 0 && iperd == 0 &&  ishift == 1) {
     for(i=1;i <= ninter;i++) {
       vcut_coul[i] = 1.0/cutoff[i];
     }/* endfor */
 }/* endif */

/*--------------------------------------------------------------------*/
  }/* end routine */
/*==========================================================================*/


/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
/* sort_inter_params:This subroutine sorts intermolecular parameters        */
/*  and finds the unique eps,sigma,and rcut for lj/null parameters          */
/*--------------------------------------------------------------------------*/

void sort_inter_params(double *eps,double *sig,
                       double *awill,double *bwill,double *cwill,
                       double *rm_swit,
                       double *c6m,double *c8m,double *c9m,double *c10m,
                       double *cutoff,double *cutoff_res,double *cutti,
                       int *inter_label,int *inter_map_index,
                       int *pninter_unique,int ninter)

/*=======================================================================*/
/*            Begin subprogram:                                          */
{/*begin routine*/
/*=======================================================================*/
/*             Local variable declarations                                */

  double *epst,*sigt;  
  double *awillt,*bwillt,*cwillt,*rm_switt;
  double *c6mt,*c8mt,*c9mt,*c10mt;
  double *cutofft,*cutoff_rest,*cuttit;
  int *ilabelt;

  double temp;
  int nunique,istart,ioff,icount,index;
  int nunique_mall;
  int nlj,nnull,n3,n5,n6,ntot,nother;
  int i,j,k,iii;
  int jlow,jup,ilow,iup;

/*=======================================================================*/
/*-------------------------------------------------*/
/*I) Malloc local inter parameters and assign them */
/*-------------------------------------------------*/

  epst     = (double *) calloc(ninter,sizeof(double))-1; 
  sigt     = (double *) calloc(ninter,sizeof(double))-1; 
  awillt   = (double *) calloc(ninter,sizeof(double))-1; 
  bwillt   = (double *) calloc(ninter,sizeof(double))-1; 
  cwillt   = (double *) calloc(ninter,sizeof(double))-1; 
  rm_switt = (double *) calloc(ninter,sizeof(double))-1; 

  c6mt     = (double *) calloc(ninter,sizeof(double))-1; 
  c8mt     = (double *) calloc(ninter,sizeof(double))-1; 
  c9mt     = (double *) calloc(ninter,sizeof(double))-1; 
  c10mt    = (double *) calloc(ninter,sizeof(double))-1; 

  cutofft     = (double *) calloc(ninter,sizeof(double))-1; 
  cutoff_rest = (double *) calloc(ninter,sizeof(double))-1; 
  cuttit      = (double *) calloc(ninter,sizeof(double))-1; 

  ilabelt = (int    *) calloc(ninter,sizeof(int))-1; 

  for(i=1; i<= ninter; i++){
    epst[i]     = eps[i];
    sigt[i]     = sig[i];
    awillt[i]   = awill[i];
    bwillt[i]   = bwill[i];
    cwillt[i]   = cwill[i];
    rm_switt[i] = rm_swit[i];
    c6mt[i]     = c6m[i];
    c8mt[i]     = c8m[i];
    c9mt[i]     = c9m[i];
    c10mt[i]    = c10m[i];

    cutofft[i]     = cutoff[i];
    cutoff_rest[i] = cutoff_res[i];
    cuttit[i]      = cutti[i];
    ilabelt[i]     = inter_label[i];
  }

  for(i=1; i<= ninter; i++){
   inter_map_index[i] = i;
  }

/*=======================================================================*/
/*-----------------------------------------------------------------*/
/*II) Count up the number of different intermolecular interactions */ 
/*-----------------------------------------------------------------*/

  nlj   = 0; 
  n3    = 0;
  nnull = 0; 
  n5    = 0;
  n6    = 0;
 for(i=1; i<= ninter; i++){
  switch(inter_label[i]){
   case 2: nlj++; break; 
   case 3: n3++;  break;
   case 4: nnull++;  break;
   case 5: n5++;  break;
   case 6: n6++;  break;
  }/*endswitch*/
 }/*endfor*/

   ntot = nlj + n3 + nnull + n5 + n6;

  if( ntot != ninter){
   printf("@@@@@@@@@@@@@@@@@@@@@ ERROR @@@@@@@@@@@@@@@@@@@@@@ \n");
   printf("TOTAL NUMBER OF INTERACTIONS NOT EQUAL TO NINTER\n");
   printf("CONTACT TECHNICAL SUPPORT \n");
   printf("@@@@@@@@@@@@@@@@@@@@@ ERROR @@@@@@@@@@@@@@@@@@@@@@ \n");
   exit(1);
  }/*endif*/

/*=======================================================================*/
/*-------------------------------------------------*/
/* III) SORT  intermolecular  interactions by type */
/*-------------------------------------------------*/

/*---------------------------------*/
/* 1. inter_label=6 Aziz-Chen last */
/*---------------------------------*/

 for(i = ninter; i > 1; i--){
  if(ilabelt[i] != 6){
   for(j=i-1; j >=  1; j--){
     if(ilabelt[j] == 6){
       switchij(&epst[i],&epst[j]);
       switchij(&sigt[i],&sigt[j]);
       switchij(&awillt[i],&awillt[j]);
       switchij(&bwillt[i],&bwillt[j]);
       switchij(&cwillt[i],&cwillt[j]);
       switchij(&rm_switt[i],&rm_switt[j]);
       switchij(&c6mt[i],&c6mt[j]);
       switchij(&c8mt[i],&c8mt[j]);
       switchij(&c9mt[i],&c9mt[j]);
       switchij(&c10mt[i],&c10mt[j]);
       switchij(&cutofft[i],&cutofft[j]);
       switchij(&cutoff_rest[i],&cutoff_rest[j]);
       switchij(&cuttit[i],&cuttit[j]);
       iswitchij(&ilabelt[i],&ilabelt[j]);
       break;
     }/*endif*/
   }/*endfor*/
  }/*endif*/
 }/*endfor*/

/*---------------------------------------------*/
/* 2. inter_label=5 Williams-LJ second to last */
/*---------------------------------------------*/

   ioff = n6;
   iup  = ninter - ioff; 
 for(i = iup; i > 1; i--){
  if(ilabelt[i] != 5){
   for(j=i-1; j >=  1; j--){
     if(ilabelt[j] == 5){
       switchij(&epst[i],&epst[j]);
       switchij(&sigt[i],&sigt[j]);
       switchij(&awillt[i],&awillt[j]);
       switchij(&bwillt[i],&bwillt[j]);
       switchij(&cwillt[i],&cwillt[j]);
       switchij(&rm_switt[i],&rm_switt[j]);
       switchij(&c6mt[i],&c6mt[j]);
       switchij(&c8mt[i],&c8mt[j]);
       switchij(&c9mt[i],&c9mt[j]);
       switchij(&c10mt[i],&c10mt[j]);
       switchij(&cutofft[i],&cutofft[j]);
       switchij(&cutoff_rest[i],&cutoff_rest[j]);
       switchij(&cuttit[i],&cuttit[j]);
       iswitchij(&ilabelt[i],&ilabelt[j]);
       break;
     }/*endif*/
   }/*endfor*/
  }/*endif*/
 }/*endfor*/

/*------------------------------------------*/
/* 3. inter_label= 3 Williams third to last */
/*------------------------------------------*/

   ioff = n6 + n5;
   iup  = ninter - ioff; 
 for(i = iup; i > 1; i--){
  if(ilabelt[i] != 3){
   for(j=i-1; j >=  1; j--){
     if(ilabelt[j] == 3){
       switchij(&epst[i],&epst[j]);
       switchij(&sigt[i],&sigt[j]);
       switchij(&awillt[i],&awillt[j]);
       switchij(&bwillt[i],&bwillt[j]);
       switchij(&cwillt[i],&cwillt[j]);
       switchij(&rm_switt[i],&rm_switt[j]);
       switchij(&c6mt[i],&c6mt[j]);
       switchij(&c8mt[i],&c8mt[j]);
       switchij(&c9mt[i],&c9mt[j]);
       switchij(&c10mt[i],&c10mt[j]);
       switchij(&cutofft[i],&cutofft[j]);
       switchij(&cutoff_rest[i],&cutoff_rest[j]);
       switchij(&cuttit[i],&cuttit[j]);
       iswitchij(&ilabelt[i],&ilabelt[j]);
       break;
     }/*endif*/
   }/*endfor*/
  }/*endif*/
 }/*endfor*/

/*--------------------------------------------*/
/* 2. PLACE inter_label=4 null fourth to last */
/*--------------------------------------------*/

   ioff = n6 + n5 + n3;
   iup  = ninter - ioff; 
 for(i = iup; i > 1; i--){
  if(ilabelt[i] != 4){
   for(j=i-1; j >=  1; j--){
     if(ilabelt[j] == 4){
       switchij(&epst[i],&epst[j]);
       switchij(&sigt[i],&sigt[j]);
       switchij(&awillt[i],&awillt[j]);
       switchij(&bwillt[i],&bwillt[j]);
       switchij(&cwillt[i],&cwillt[j]);
       switchij(&rm_switt[i],&rm_switt[j]);
       switchij(&c6mt[i],&c6mt[j]);
       switchij(&c8mt[i],&c8mt[j]);
       switchij(&c9mt[i],&c9mt[j]);
       switchij(&c10mt[i],&c10mt[j]);
       switchij(&cutofft[i],&cutofft[j]);
       switchij(&cutoff_rest[i],&cutoff_rest[j]);
       switchij(&cuttit[i],&cuttit[j]);
       iswitchij(&ilabelt[i],&ilabelt[j]);
       break;
     }/*endif*/
   }/*endfor*/
  }/*endif*/
 }/*endfor*/

 
/*=======================================================================*/
/*--------------------------------------------*/
/*IV) SORT LJ and NULL Parameters             */
/*    same eps sig rcut in order              */
/*--------------------------------------------*/

 ioff   = n3+n5+n6; 
 i = 1;
 while(i< ninter-ioff){
  k = i;
  for(j=i+1; j <=  ninter-ioff; j++){
    /* Bubble down matching interactions*/
    if( (epst[i]        == epst[j]) && 
        (sigt[i]        == sigt[j]) &&
        (cutofft[i]     == cutofft[j]) &&
        (cutti[i]       == cutti[j]) && 
        (cutoff_rest[i] == cutoff_rest[j])){
      k++;
      switchij(&epst[k],&epst[j]);
      switchij(&sigt[k],&sigt[j]);
      switchij(&awillt[k],&awillt[j]);
      switchij(&bwillt[k],&bwillt[j]);
      switchij(&cwillt[k],&cwillt[j]);
      switchij(&rm_switt[k],&rm_switt[j]);
      switchij(&c6mt[k],&c6mt[j]);
      switchij(&c8mt[k],&c8mt[j]);
      switchij(&c9mt[k],&c9mt[j]);
      switchij(&c10mt[k],&c10mt[j]);
      switchij(&cutofft[k],&cutofft[j]);
      switchij(&cutoff_rest[k],&cutoff_rest[j]);
      switchij(&cuttit[k],&cuttit[j]);
      iswitchij(&ilabelt[k],&ilabelt[j]);
    }/*endif*/
  }/*endfor*/
  i = k+1;      /* The next unique guy is at k+1 */
 }/*end while */

/*=======================================================================*/
/*--------------------------------------------*/
/*V) Determine number of unique LJ parameters */
/*    and reassign in new order               */
/*--------------------------------------------*/


/*Assign unique LJ/Null parameters */

 nother = n3+n5+n6;
 i = 1;
 nunique = 0;
 while(i<=ninter-ioff){
   nunique++;
   j = i;
   while( ((epst[i]        == epst[j])        &&
           (sigt[i]        == sigt[j])        &&
           (cutofft[i]     == cutofft[j])     &&
           (cutoff_rest[i] == cutoff_rest[j]) &&
           (cuttit[i]      == cuttit[j])     )){
      epst[nunique]        = epst[j];
      sigt[nunique]        = sigt[j];
      awillt[nunique]      = awillt[j];
      bwillt[nunique]      = bwillt[j];
      cwillt[nunique]      = cwillt[j];
      rm_switt[nunique]    = rm_switt[j];
      c6mt[nunique]        = c6mt[j];
      c8mt[nunique]        = c8mt[j];
      c9mt[nunique]        = c9mt[j];
      c10mt[nunique]       = c10mt[j];

      cutofft[nunique]     = cutofft[j];
      cutoff_rest[nunique] = cutoff_rest[j];
      cuttit[nunique]      = cuttit[j];
      ilabelt[nunique]     = ilabelt[j];
      j++;
      if(j==ninter-ioff+1){break;}
    }/*while*/
    i = j;
  }/*endfor*/


/*Reassign all other inter parameters*/
     icount = 0;

 for(i=ninter-nother+1; i<= ninter; i++){
   icount++; 
   index = nunique+icount;
     epst[index]      = epst[i];
     sigt[index]      = sigt[i];
   awillt[index]      = awillt[i];
   bwillt[index]      = bwillt[i];
   cwillt[index]      = cwillt[i];
   rm_switt[index]    = rm_switt[i];
     c6mt[index]      = c6mt[i];
     c8mt[index]      = c8mt[i];
     c9mt[index]      = c9mt[i];
    c10mt[index]      = c10mt[i];
   cutofft[index]     = cutofft[i];
   cutoff_rest[index] = cutoff_rest[i];
   cuttit[index]      = cuttit[i];
   ilabelt[index]     = ilabelt[i];

 }/*endfor*/


   nunique += nother;
  *pninter_unique = nunique;

  if(nunique>ninter){
   printf("@@@@@@@@@@@@@@@@@@@@@ ERROR @@@@@@@@@@@@@@@@@@@@@@ \n");
   printf("Total number of unique interactions > then ninter!!\n");
   printf("Contact technical support : %d vs %d\n",ninter,nunique);
   printf("@@@@@@@@@@@@@@@@@@@@@ ERROR @@@@@@@@@@@@@@@@@@@@@@ \n");
   exit(1);
  }

  if(nunique==0 && ninter!=0){
   printf("@@@@@@@@@@@@@@@@@@@@@ ERROR @@@@@@@@@@@@@@@@@@@@@@ \n");
   printf("Total number of unique interactions  = 0!!\n");
   printf("Contact technical support");
   printf("@@@@@@@@@@@@@@@@@@@@@ ERROR @@@@@@@@@@@@@@@@@@@@@@ \n");
   exit(1);
  }


/*=======================================================================*/
/*----------------------------------------------*/
/*VI) Assign lj_ind base on unique eps sig rcut */
/*----------------------------------------------*/

  for(i=1; i<= ninter; i++){
   inter_map_index[i] = 0;
  }

  for(i=1; i<= ninter; i++){
   if(inter_label[i] == 2 || 
      inter_label[i] == 4){

    for(j=1; j<= (nunique-nother); j++){
     if( (eps[i] == epst[j]) && 
         (sig[i] == sigt[j]) && 
         (cutoff[i] == cutofft[j]) &&
         (cutti[i] == cuttit[j]) &&
         (cutoff_res[i] == cutoff_rest[j]) ){
       inter_map_index[i] = j; 
     }/*endif*/ 
    }/*endfor*/

   }/*endif*/
  }/*endfor*/

/*-----------------------------------------------------------------*/
/* Assign  Williams  compare awill bwill c6m c8m c10m and cutoffs  */
/*-----------------------------------------------------------------*/


   jlow = nunique - nother + 1;
   jup  = jlow + n3;

  for(i=1; i<= ninter; i++){
   if(inter_label[i] == 3 ){
    for(j=jlow; j<= jup; j++){ 
     if((awill[i] == awillt[j]) && 
        (bwill[i] == bwillt[j]) &&
        (c6m[i]==c6mt[j]) && 
        (c8m[i]==c8mt[j]) && 
        (c10m[i] == c10mt[j]) &&
        (cutoff[i] == cutofft[j]) && 
        (cutti[i] == cuttit[j]) &&
        (cutoff_res[i] == cutoff_res[j]) ){
       inter_map_index[i] = j; 
     }/*endif*/ 
    }/*endfor*/
   }/*endif*/
  }/*endfor*/


/*---------------------------------------------------------------------------*/
/* Assign  Williams-LJ compare awill bwill c6m c8m c10m eps sig and cutoffs  */
/*---------------------------------------------------------------------------*/
   jlow = nunique - nother + n3 + 1; 
   jup  = jlow + n5;

  for(i=1; i<= ninter; i++){
   if(inter_label[i] == 5 ){
    for(j=jlow; j<= jup; j++){
     if( (eps[i] == epst[j]) && 
         (sig[i] == sigt[j]) && 
         (awill[i] == awillt[j]) && 
         (bwill[i] == bwillt[j]) &&
         (c6m[i] == c6mt[j]) && 
         (c8m[i] == c8mt[j]) && 
         (c10m[i] == c10mt[j]) &&
         (cutoff[i] == cutofft[j]) && 
         (cutti[i] == cuttit[j]) &&
         (cutoff_res[i] == cutoff_rest[j])  ){
       inter_map_index[i] = j; 
     }/*endif*/ 
    }/*endfor*/
   }/*endif*/
  }/*endfor*/

/*--------------------------------------------------------------------------*/
/* Assign Aziz-Chen compare awill bwill cwill rm_swit c6m c8m c9m c10m  */
/* and cutoffs  */
/*--------------------------------------------------------------------------*/

   jlow = nunique - nother + n3 + n5 + 1; 
   jup  = jlow + n6;


  for(i=1; i<= ninter; i++){
   if(inter_label[i] == 6 ){

    for(j=jlow; j<= jup; j++){
     if((awill[i] == awillt[j]) &&
        (bwill[i] == bwillt[j]) &&
        (cwill[i] == cwillt[j]) && 
        (rm_swit[i] == rm_switt[j]) && 
        (c6m[i] == c6mt[j]) && 
        (c8m[i] == c8mt[j]) && 
        (c9m[i] == c9mt[j]) && 
        (c10m[i] == c10mt[j]) &&
        (cutoff[i] == cutofft[j]) && 
        (cutti[i] == cuttit[j]) &&
        (cutoff_res[i] == cutoff_rest[j])  ){
          inter_map_index[i] = j; 
      }/*endif*/ 
    }/*endfor*/

   }/*endif*/
  }/*endfor*/


/*-----------------------------------------------*/
/*  Reassign the unique atom types               */
/*-----------------------------------------------*/

   for(i=1; i<= nunique; i++){
     eps[i]        = epst[i];
     sig[i]        = sigt[i];
     awill[i]      = awillt[i];
     bwill[i]      = bwillt[i];
     cwill[i]      = cwillt[i];
     rm_swit[i]    = rm_switt[i];
     c6m[i]        = c6mt[i];
     c8m[i]        = c8mt[i];
     c9m[i]        = c9mt[i];
     c10m[i]       = c10mt[i];
     cutoff[i]     = cutofft[i];
     cutti[i]      = cuttit[i];
     cutoff_res[i] = cutoff_rest[i];
    inter_label[i] = ilabelt[i];
   }/*endfor*/

/*-----------------------------------------------*/
/* free locally assigned memory                  */
/*-----------------------------------------------*/

  cfree(&epst[1]);
  cfree(&sigt[1]);
  cfree(&awillt[1]);
  cfree(&bwillt[1]);
  cfree(&cwillt[1]);
  cfree(&rm_switt[1]);

  cfree(&c6mt[1]);
  cfree(&c8mt[1]);
  cfree(&c9mt[1]);
  cfree(&c10mt[1]);

  cfree(&cutofft[1]);
  cfree(&cutoff_rest[1]);
  cfree(&cuttit[1]);

  cfree(&ilabelt[1]);

/*==========================================================================*/
}/*end routine*/
/*==========================================================================*/

/*==========================================================================*/
/*==========================================================================*/
void switchij(double *valuei,double *valuej)
/*==========================================================================*/
{/*begin routine*/
/*==========================================================================*/
 double temp;

  temp     = *valuei;
  *valuei  = *valuej;
  *valuej  = temp;
 

/*==========================================================================*/
}/*end routine*/
/*==========================================================================*/
/*==========================================================================*/
/*==========================================================================*/
void iswitchij(int *valuei,int *valuej)
/*==========================================================================*/
{/*begin routine*/
/*==========================================================================*/
 int temp;

  temp     = *valuei;
  *valuei  = *valuej;
  *valuej  = temp;
 

/*==========================================================================*/
}/*end routine*/
/*==========================================================================*/
