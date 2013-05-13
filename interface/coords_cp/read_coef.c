/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
/*                                                                          */
/*                         PI_MD:                                           */
/*             The future of simulation technology                          */
/*             ------------------------------------                         */
/*                     Module: read_coef.c                                  */
/*                                                                          */
/* This subprogram provides output for a MD on a                            */
/* LD-classical potential energy surface (LD-PES)                           */
/*                                                                          */
/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

#include "standard_include.h"
#include "../typ_defs/typedefs_class.h"
#include "../typ_defs/typedefs_par.h"
#include "../typ_defs/typedefs_gen.h"
#include "../typ_defs/typedefs_cp.h"
#include "../proto_defs/proto_energy_cpcon_entry.h"
#include "../proto_defs/proto_coords_cp_entry.h"
#include "../proto_defs/proto_handle_entry.h"
#include "../proto_defs/proto_friend_lib_entry.h"
#include "../proto_defs/proto_communicate_wrappers.h"

void read_coef_sys_info(CP *, FILE *,int , char *);
void read_coef_alloc_init(CP *, int ,double *);
void read_coef_fetch_occs(CP *, FILE *, char *);
void read_coef_fetch_coefs(CP *, FILE *, char *,int);
void read_coef_fetch_vcoefs(CP *, FILE *, char *);
void read_coef_fetch_vnhc(CP *, FILE *, char *);
void read_coef_init_nhc(CP *);
void read_coef_transpose(CP *,int ,int, int);

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void read_coef(CP *cp,GENERAL_DATA *general_data,CLASS *class,
               FILENAME_PARSE *filename_parse,
               CP_PARSE *cp_parse,double *tot_memory)

/*==========================================================================*/
/*               Begin subprogram:                                          */
  {/*begin routine*/
/*==========================================================================*/
/*               Local variable declarations                                */
#include "../typ_defs/typ_mask.h"

  FILE *fp_dnameci;   

/* Local pointers */
  SIMOPTS *simopts       = &(general_data->simopts);
  int istart             = cp_parse->istart_cp;
  int myid               = cp->communicate.myid;
  MPI_Comm world         = cp->communicate.world;
  int cp_wave_min        = simopts->cp_wave_min;
  int cp_wave_min_pimd   = simopts->cp_wave_min_pimd;
  int cp_min             = simopts->cp_min;
  int cp_min_on;
  int initial_spread_opt = simopts->initial_spread_opt;
  int np_states          = cp->communicate.np_states;
  int num_proc           = cp->communicate.np;
  char *dnameci          = filename_parse->dnameci;

  int ibinary            = cp->cpopts.iread_coef_binary;
  int iii;

  cp_min_on = cp_wave_min + cp_wave_min_pimd + cp_min;

/*========================================================================*/
/*  II)Write to screen:                                                   */

  if( myid==0 ){
    printf("\n");PRINT_LINE_STAR;
   if(istart >= 1){ printf("Reading user specified CP coordinate file %s\n",dnameci);}
    if(istart==1){printf("using the `initial' restart option\n");}
    if(istart==2){printf("using the `restart_pos' restart option\n");}
    if(istart==3){printf("using the `restart_posvel' restart option\n");}
    if(istart==4){printf("using the `restart_all' restart option\n");}
    if(istart==0){printf("using the `gen_wave' restart option\n");}
    PRINT_LINE_DASH;printf("\n");
  }/*endif */   

/*========================================================================*/
/* III) Open coefficient dump file                                        */

  if((myid==0) && (istart >= 1)){
  if(ibinary == 0){
   fp_dnameci = cfopen(dnameci,"r");
  }else{
   fp_dnameci = cfopen(dnameci,"rb");
  }/*endif ibinary*/

  }/*endif*/

/*========================================================================*/
/* IV) Read System information                                            */

  if(myid==0 && istart >= 1){
   read_coef_sys_info(cp,fp_dnameci,istart,dnameci);
  }/*endif*/
  if(num_proc>1){Barrier(world);}


/*========================================================================*/
/*  V) Allocate and initialize coefficient arrays                         */

  read_coef_alloc_init(cp,cp_min_on,tot_memory);
  if(num_proc>1){Barrier(world);}

/*========================================================================*/
/*  VI) Read and communicate the occupation numbers                       */

  if(istart >= 1){
   read_coef_fetch_occs(cp,fp_dnameci,dnameci);
  }/*endif*/
  if(num_proc>1){Barrier(world);}

/*========================================================================*/
/*  VII) Read/Spread the Plane wave coefficients                          */

  if(istart >= 1){
   read_coef_fetch_coefs(cp,fp_dnameci,dnameci,initial_spread_opt);
  }else{
   gen_wave(class,general_data,cp,cp_parse,filename_parse->vps_name);
  }
  if(num_proc>1){Barrier(world);}

/*========================================================================*/
/* VIII) Read/Communicate Coeffs Vels                                     */

  if((istart >= 3)&&((cp_wave_min+cp_min)==0)) {
    read_coef_fetch_vcoefs(cp,fp_dnameci,dnameci);
  } /* endif for istart >=3 */
  if(num_proc>1){Barrier(world);}

/*========================================================================*/
/* IX) Read/communicate Extended system variables                         */
  
  if((istart == 4)&&((cp_wave_min+cp_min)==0)) {
    read_coef_fetch_vnhc(cp,fp_dnameci,dnameci);
  } /* endif:istart=4 */
  if(num_proc>1){Barrier(world);}
  read_coef_init_nhc(cp);

/*=======================================================================*/
/* X) Close the coef dump file. Reading time is done                     */

  if( (myid==0) && (istart > 0) ){
    fclose(fp_dnameci);
  }/*endif*/

/*=======================================================================*/
/* XI) Transpose the coefs, coef vels and thermostats (massiv only)      */

  if(np_states>1){
    read_coef_transpose(cp,istart,cp_wave_min,cp_min);
  }/*endif np_states > 1*/
  if(num_proc>1){Barrier(world);}

/*=======================================================================*/
/* XII) Write to the screen                                              */

  if(myid==0){

    printf("\n");
    PRINT_LINE_DASH;
    printf("Done reading user specified CP coordinate file %s\n",dnameci);
    PRINT_LINE_STAR;printf("\n");

  }/*endif myid==0*/


/*-----------------------------------------------------------------------*/
  }/* end routine */
/*==========================================================================*/





/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void read_coef_sys_info(CP *cp, FILE *fp_dnameci,int istart, char *dnameci)
              
/*==========================================================================*/
 {/*begin routine */
/*==========================================================================*/
/*         Local Variables */
#include "../typ_defs/typ_mask.h"
  int n,ibinary;
  int iii;

  char *c_array;
  char *dft_type_now,    *norb_opt_now;
  char *restart_type_now,*restart_type_spec;
  int ncoef_up_now,  ncoef_dn_now;
  int nstate_up_now, nstate_dn_now;
  int itime_dump,    istart_now;
  int cp_lda_now,    cp_lsda_now;

/*         Local Pointers */
  int cp_norb        = cp->cpopts.cp_norb;
  int cp_init_orthog = cp->cpopts.cp_init_orthog;
  int cp_lda         = cp->cpopts.cp_lda;
  int cp_lsda        = cp->cpopts.cp_lsda;
  int ncoef_up       = cp->cpcoeffs_info.ncoef;
  int ncoef_dn       = cp->cpcoeffs_info.ncoef;
  int nstate_up      = cp->cpcoeffs_info.nstate_up;
  int nstate_dn      = cp->cpcoeffs_info.nstate_dn;
  if(cp_lda==1){ncoef_dn=0;}

  ibinary = cp->cpopts.iread_coef_binary;


/*==========================================================================*/
/* 0) Malloc up some memory */

  restart_type_now  = (char *) cmalloc(MAXWORD*sizeof(char));
  restart_type_spec = (char *) cmalloc(MAXWORD*sizeof(char));
  dft_type_now      = (char *) cmalloc(MAXWORD*sizeof(char));
  norb_opt_now      = (char *) cmalloc(MAXWORD*sizeof(char));

  if(ibinary == 1){c_array = (char *) cmalloc (MAXWORD*sizeof(char)); }


/*==========================================================================*/
/* I) Read the header  */

  if(ibinary == 0){
    readtoendofline(fp_dnameci);
    fscanf(fp_dnameci,"%d %d %d %d %s %s %d %s",
           &ncoef_up_now,&ncoef_dn_now,
           &nstate_up_now,&nstate_dn_now,
           dft_type_now,restart_type_now,&itime_dump,norb_opt_now);
    readtoendofline(fp_dnameci);
  }else{

    fread(c_array,sizeof(char),MAXWORD,fp_dnameci);
    fread(c_array,sizeof(char),MAXWORD,fp_dnameci);
    fread(c_array,sizeof(char),MAXWORD,fp_dnameci);
    fread(c_array,sizeof(char),MAXWORD,fp_dnameci);
    fread(c_array,sizeof(char),MAXWORD,fp_dnameci);
    fread(c_array,sizeof(char),MAXWORD,fp_dnameci);
    fread(c_array,sizeof(char),MAXWORD,fp_dnameci);

    n = 1;
    fread(&ncoef_up_now,sizeof(int),n,fp_dnameci);
    fread(&ncoef_dn_now ,sizeof(int),n,fp_dnameci);
    fread(&nstate_up_now,sizeof(int),n,fp_dnameci);
    fread(&nstate_dn_now,sizeof(int),n,fp_dnameci);

    fread(dft_type_now,sizeof(char),MAXWORD,fp_dnameci);
    fread(restart_type_now ,sizeof(char),MAXWORD,fp_dnameci);
    fread(&itime_dump ,sizeof(int),n,fp_dnameci);
    fread(norb_opt_now,sizeof(char),MAXWORD,fp_dnameci);

  }/*endif*/

/*==========================================================================*/
/* II) Check the start  */

  istart_now = -1;
  if(strcasecmp(restart_type_now,"initial") == 0)       istart_now = 1;
  if(strcasecmp(restart_type_now,"restart_pos") == 0)   istart_now = 2;
  if(strcasecmp(restart_type_now,"restart_posvel") == 0)istart_now = 3;
  if(strcasecmp(restart_type_now,"restart_all") == 0)   istart_now = 4;

  if(istart==1){strcpy(restart_type_spec,"initial");}
  if(istart==2){strcpy(restart_type_spec,"restart_pos");}
  if(istart==3){strcpy(restart_type_spec,"restart_posvel");}
  if(istart==4){strcpy(restart_type_spec,"restart_all");}

  if(istart_now < istart) {
     printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
     printf("Start up option, %s, in\n",restart_type_now);
     printf("user specified coordinate file %s \n",dnameci);
     printf("incompatible with system setup,%s\n",restart_type_spec);
     printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
     fflush(stdout);exit(1);
  }/* endif */

/*==========================================================================*/
/* III) Check the norb opt  */

  if((strcasecmp(norb_opt_now,"norb_on")==0)&&(cp_norb==0)
     &&(cp_init_orthog==0)){
     printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
     printf("A nonorthogonal wave function cannot be used\n");
     printf("to start a standard CP calculation without \n");
     printf("an initial orthogonalization %s\n",norb_opt_now);
     printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
     fflush(stdout);exit(1);
  }/*endif*/

/*==========================================================================*/
/* IV) Check lsda/lda stuff  */

  cp_lda_now=-1;cp_lsda_now=-1;
  if(strcasecmp(dft_type_now,"lda") == 0){cp_lda_now=1;cp_lsda_now=0;}
  if(strcasecmp(dft_type_now,"lsda") == 0){cp_lda_now=0;cp_lsda_now=1;}

  if(cp_lda_now != cp_lda) {
     printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
     printf("CP LDA option in \n");
     printf("user specified plane wave coefficient file %s\n",dnameci);
     printf("incompatible with system set up\n");
     printf("%d vs %d\n",cp_lda,cp_lda_now);
     printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
     fflush(stdout);exit(1);
  } /* endif */
  if(cp_lsda_now != cp_lsda) {
     printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
     printf("CP LSDA option in \n");
     printf("user specified plane wave coefficient file %s\n",dnameci);
     printf("incompatible with system set up\n");
     printf("%d vs %d\n",cp_lsda,cp_lsda_now);
     printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
     fflush(stdout);exit(1);
  } /* endif */

/*==========================================================================*/
/* V) Check coef/state sizes  */

  if(ncoef_up_now != ncoef_up) {
     printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
     printf("Number of spin-up plane wave coefficients in\n");
     printf("user specified plane wave coefficient file %s\n",dnameci);
     printf("incompatible with system set up\n");
     printf("%d vs %d\n",ncoef_up,ncoef_up_now);
     printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
     fflush(stdout);exit(1);
  } /* endif */

  if(ncoef_dn_now != ncoef_dn) {
     printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
     printf("Number of spin-dn plane wave coefficients in\n");
     printf("user specified plane wave coefficient file %s\n",dnameci);
     printf("incompatible with system set up\n");
     printf("%d vs %d\n",ncoef_dn,ncoef_dn_now);
     printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
     fflush(stdout);exit(1);
  } /* endif */

  if(nstate_up_now != nstate_up) {
     printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
     printf("Number of spin-up states in\n");
     printf("user specified plane wave coefficient file %s\n",dnameci);
     printf("incompatible with system set up\n");
     printf("%d vs %d\n",nstate_up,nstate_up_now);
     printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
     fflush(stdout);exit(1);
  } /* endif */

  if(nstate_dn_now != nstate_dn) {
     printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
     printf("Number of spin-dn states in\n");
     printf("user specified plane wave coefficient file %s\n",dnameci);
     printf("incompatible with system set up\n");
     printf("%d vs %d\n",nstate_dn,nstate_dn_now);
     printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
     fflush(stdout);exit(1);
  } /* endif */

/*==========================================================================*/
/* VI) Free up some memory */

  cfree(restart_type_now);
  cfree(restart_type_spec);
  cfree(dft_type_now);
  cfree(norb_opt_now);

  if (ibinary == 1){
   cfree(c_array);
  }

/*-----------------------------------------------------------------------*/
    } /* end routine */
/*==========================================================================*/




/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void read_coef_alloc_init(CP *cp,int cp_min_on,double *tot_memory)
              
/*==========================================================================*/
  {/*begin routine */
/*==========================================================================*/
/*    Local Variables   */
#include "../typ_defs/typ_mask.h"
  int myid = cp->communicate.myid;
  int i,ip,nread,is;
  int par_size_up,par_size_dn,ncoef_up_tot;
  int nstate,nstate2,ncoef_dn_tot;
  double *cre,*cim,*vcre,*vcim,*fcre,*fcim;
  int *ioff_up,*ioff_upt,*ioff_dn,*ioff_dnt;
  double mem_test;

/*  Local Pointers */
  int pi_beads       = cp->cpcoeffs_info.pi_beads;
  int pi_beads_proc  = cp->cpcoeffs_info.pi_beads_proc;
  int nstate_up_proc = cp->cpcoeffs_info.nstate_up_proc;
  int nstate_dn_proc = cp->cpcoeffs_info.nstate_dn_proc;
  int nstate_up      = cp->cpcoeffs_info.nstate_up;
  int nstate_dn      = cp->cpcoeffs_info.nstate_dn;
  int cp_lda         = cp->cpopts.cp_lda;
  int cp_lsda        = cp->cpopts.cp_lsda;
  int ncoef_up       = cp->cpcoeffs_info.ncoef;
  int ncoef_dn       = cp->cpcoeffs_info.ncoef;
  int cp_norb        = cp->cpopts.cp_norb;

  int np_states               = cp->communicate.np_states;
  int ncoef_up_proc           = cp->cp_comm_state_pkg_up.nstate_ncoef_proc_max;
  int ncoef_dn_proc           = cp->cp_comm_state_pkg_dn.nstate_ncoef_proc_max;
  int nstate_max_up           =cp->cp_comm_state_pkg_up.nstate_max;
  int nstate_ncoef_proc_max_up=cp->cp_comm_state_pkg_up.nstate_ncoef_proc_max;
  int nstate_max_dn           =cp->cp_comm_state_pkg_dn.nstate_max;
  int nstate_ncoef_proc_max_dn=cp->cp_comm_state_pkg_dn.nstate_ncoef_proc_max;

  if(cp_lda == 1){ncoef_dn = 0;}



/*==========================================================================*/
/* 0) Calculate the sizes */

  par_size_up = nstate_max_up*(nstate_ncoef_proc_max_up);
  ncoef_up_tot = nstate_up_proc*ncoef_up;
  ncoef_up_tot = MAX(ncoef_up_tot,par_size_up);

  par_size_dn = nstate_max_dn*(nstate_ncoef_proc_max_dn);
  if(cp_lda ==1){par_size_dn = 0;}
  ncoef_dn_tot = MAX(nstate_dn_proc*ncoef_dn,1);
  ncoef_dn_tot = MAX(ncoef_dn_tot,par_size_dn);

  nstate  = MAX(nstate_up,nstate_dn);
  nstate2 = nstate*nstate;

/*==========================================================================*/
/* I) Malloc the variables */

 for(i=1;i<=pi_beads_proc;i++){
  cp->cpcoeffs_pos[i].cre_up =(double *)cmalloc(ncoef_up_tot*sizeof(double))-1;
  cp->cpcoeffs_pos[i].cim_up =(double *)cmalloc(ncoef_up_tot*sizeof(double))-1;
  cp->cpcoeffs_pos[i].cre_dn =(double *)cmalloc(ncoef_dn_tot*sizeof(double))-1;
  cp->cpcoeffs_pos[i].cim_dn =(double *)cmalloc(ncoef_dn_tot*sizeof(double))-1;
  cp->cpcoeffs_pos[i].vcre_up=(double *)cmalloc(ncoef_up_tot*sizeof(double))-1;
  cp->cpcoeffs_pos[i].vcim_up=(double *)cmalloc(ncoef_up_tot*sizeof(double))-1;
  cp->cpcoeffs_pos[i].vcre_dn=(double *)cmalloc(ncoef_dn_tot*sizeof(double))-1;
  cp->cpcoeffs_pos[i].vcim_dn=(double *)cmalloc(ncoef_dn_tot*sizeof(double))-1;
  cp->cpcoeffs_pos[i].fcre_up=(double *)cmalloc(ncoef_up_tot*sizeof(double))-1;
  cp->cpcoeffs_pos[i].fcim_up=(double *)cmalloc(ncoef_up_tot*sizeof(double))-1;
  cp->cpcoeffs_pos[i].fcre_dn=(double *)cmalloc(ncoef_dn_tot*sizeof(double))-1;
  cp->cpcoeffs_pos[i].fcim_dn=(double *)cmalloc(ncoef_dn_tot*sizeof(double))-1;
  cp->cpcoeffs_pos[i].ksmat_up    =(double *)cmalloc(nstate2*sizeof(double))-1;
  cp->cpcoeffs_pos[i].ksmat_dn    =(double *)cmalloc(nstate2*sizeof(double))-1;
  cp->cpcoeffs_pos[i].ksmat_eig_up=(double *)cmalloc(nstate*sizeof(double))-1;
  cp->cpcoeffs_pos[i].ksmat_eig_dn=(double *)cmalloc(nstate*sizeof(double))-1;
  cp->cpcoeffs_pos[i].norbmat_up  =(double *)cmalloc(nstate2*sizeof(double))-1;
  cp->cpcoeffs_pos[i].norbmat_dn  =(double *)cmalloc(nstate2*sizeof(double))-1;
  cp->cpcoeffs_pos[i].norbmati_up =(double *)cmalloc(nstate2*sizeof(double))-1;
  cp->cpcoeffs_pos[i].norbmati_dn =(double *)cmalloc(nstate2*sizeof(double))-1;
  cp->cpcoeffs_pos[i].ovmat_eigv_up=(double *)
                                  cmalloc(nstate2*sizeof(double))-1;
  cp->cpcoeffs_pos[i].ovmat_eigv_dn=(double *)
                                  cmalloc(nstate2*sizeof(double))-1;
  if(cp_min_on > 0){  /* Then we need to allocate the diagonal Hessian */
    cp->cpcoeffs_pos[i].cp_hess_re_up = (double *) cmalloc(ncoef_up*sizeof(double))-1;
    cp->cpcoeffs_pos[i].cp_hess_im_up = (double *) cmalloc(ncoef_up*sizeof(double))-1;
    cp->cpcoeffs_pos[i].cp_hess_re_dn = (double *) cmalloc(ncoef_dn*sizeof(double))-1;
    cp->cpcoeffs_pos[i].cp_hess_im_dn = (double *) cmalloc(ncoef_dn*sizeof(double))-1;
  }/* endif cp_min_on */
 }/*endfor*/

  cp->cpcoeffs_info.ioff_up  = (int *)cmalloc(nstate*sizeof(int))-1;
  cp->cpcoeffs_info.ioff_dn  = (int *)cmalloc(nstate*sizeof(int))-1;
  cp->cpcoeffs_info.ioff_upt = (int *)cmalloc(nstate*sizeof(int))-1;
  cp->cpcoeffs_info.ioff_dnt = (int *)cmalloc(nstate*sizeof(int))-1;
  cp->cpopts.occ_up          = (double *) cmalloc(nstate*sizeof(double))-1;
  cp->cpopts.occ_dn          = (double *) cmalloc(nstate*sizeof(double))-1;
  cp->cpopts.rocc_sum_up = (double *) cmalloc(nstate*nstate*sizeof(double))-1;
  cp->cpopts.rocc_sum_dn = (double *) cmalloc(nstate*nstate*sizeof(double))-1;

 if(myid==0){
   mem_test = (pi_beads_proc*(ncoef_up_tot + ncoef_dn_tot)*6*sizeof(double)
             + nstate*2*sizeof(int)) *1.0e-6;
   if(cp_min_on > 0) mem_test += (pi_beads_proc*(2*ncoef_up + 2*ncoef_dn))*1.0e-6;
   *tot_memory += mem_test;
   printf("CP wave func allocation: %g Mbytes; Total memory: %g Mbytes\n",
          mem_test,*tot_memory);

 }/*endif for myid==0*/

/*==========================================================================*/
/* II) Assign the offsets */

  ioff_up  = cp->cpcoeffs_info.ioff_up;
  ioff_upt = cp->cpcoeffs_info.ioff_upt;
  ioff_dn  = cp->cpcoeffs_info.ioff_dn;
  ioff_dnt = cp->cpcoeffs_info.ioff_dnt;

  for(is=1;is<=nstate;is++){ioff_up[is]=(is-1)*ncoef_up;}
  for(is=1;is<=nstate;is++){ioff_dn[is]=(is-1)*ncoef_dn;}

  if(np_states==1){
   for(is=1;is<=nstate;is++){ioff_upt[is]=(is-1)*ncoef_up;}
   for(is=1;is<=nstate;is++){ioff_dnt[is]=(is-1)*ncoef_dn;}
  }else{
   for(is=1;is<=nstate;is++){ioff_upt[is]=(is-1)*ncoef_up_proc;}
   for(is=1;is<=nstate;is++){ioff_dnt[is]=(is-1)*ncoef_dn_proc;}
  }/*endif*/

/*========================================================================*/
/* III) Initialize coeficient arrays and form flags                       */

 for(ip=1;ip<=pi_beads_proc;ip++){

   cp->cpcoeffs_pos[ip].icoef_form_up  = 0;
   cp->cpcoeffs_pos[ip].ivcoef_form_up = 0;
   cp->cpcoeffs_pos[ip].ifcoef_form_up = 1;
   cp->cpcoeffs_pos[ip].icoef_orth_up  = 1;
   if(cp_norb>0){cp->cpcoeffs_pos[ip].icoef_orth_up = 0;}
   cp->cpcoeffs_pos[ip].ivcoef_orth_up  = 1;
   if(cp_norb>0){cp->cpcoeffs_pos[ip].ivcoef_orth_up = 0;}

   cre = cp->cpcoeffs_pos[ip].cre_up;
   cim = cp->cpcoeffs_pos[ip].cim_up;
   vcre = cp->cpcoeffs_pos[ip].vcre_up;
   vcim = cp->cpcoeffs_pos[ip].vcim_up;
   fcre = cp->cpcoeffs_pos[ip].fcre_up;
   fcim = cp->cpcoeffs_pos[ip].fcim_up;
   nread = ncoef_up_tot;
   for(i=1;i<=nread;i++){
    cre[i] = 0.0;
    cim[i] = 0.0;
    vcre[i] = 0.0;
    vcim[i] = 0.0;
    fcre[i] = 0.0;
    fcim[i] = 0.0;
   }/*endfor*/

   if(cp_lsda==1){

    cp->cpcoeffs_pos[ip].icoef_form_dn  = 0; 
    cp->cpcoeffs_pos[ip].ivcoef_form_dn = 0;
    cp->cpcoeffs_pos[ip].ifcoef_form_dn = 1;
    cp->cpcoeffs_pos[ip].icoef_orth_dn  = 1;
    if(cp_norb>0){cp->cpcoeffs_pos[ip].icoef_orth_dn = 0;}
    cp->cpcoeffs_pos[ip].ivcoef_orth_dn  = 1;
    if(cp_norb>0){cp->cpcoeffs_pos[ip].ivcoef_orth_dn = 0;}

    cre = cp->cpcoeffs_pos[ip].cre_dn;
    cim = cp->cpcoeffs_pos[ip].cim_dn;
    vcre = cp->cpcoeffs_pos[ip].vcre_dn;
    vcim = cp->cpcoeffs_pos[ip].vcim_dn;
    fcre = cp->cpcoeffs_pos[ip].fcre_dn;
    fcim = cp->cpcoeffs_pos[ip].fcim_dn;
    nread = ncoef_dn_tot;
    for(i=1;i<=nread;i++){
     cre[i] = 0.0;
     cim[i] = 0.0;
     vcre[i] = 0.0;
     vcim[i] = 0.0;
     fcre[i] = 0.0;
     fcim[i] = 0.0;
    }/*endfor*/ 

   }/*endif:lsda*/

 }/*endfor:pi_bead*/

/*-----------------------------------------------------------------------*/
  }/* end routine */
/*==========================================================================*/




/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void read_coef_fetch_occs(CP *cp, FILE *fp_dnameci, char *dnameci)
              
/*==========================================================================*/
  {/*begin routine */
/*==========================================================================*/
/*               Local variable declarations                                */
#include "../typ_defs/typ_mask.h"

  int nstate,i,j,iocc,is,js;

/* Local pointers */
  int myid            = cp->communicate.myid;
  int nproc           = cp->communicate.np;
  MPI_Comm world      = cp->communicate.world;
  int cp_norb         = cp->cpopts.cp_norb;
  int cp_lsda         = cp->cpopts.cp_lsda;
  int cp_lda          = cp->cpopts.cp_lda;
  double *occ_dn      = cp->cpopts.occ_dn;
  double *occ_up      = cp->cpopts.occ_up;
  double *rocc_sum_up = cp->cpopts.rocc_sum_up;
  double *rocc_sum_dn = cp->cpopts.rocc_sum_dn;
  int nstate_up       = cp->cpcoeffs_info.nstate_up;
  int nstate_dn       = cp->cpcoeffs_info.nstate_dn;

  char *c_array1,*c_array2;
  int ibinary,n;

  ibinary = cp->cpopts.iread_coef_binary;
  
  nstate  = MAX(nstate_up,nstate_dn);

/*==========================================================================*/
/*  I) Reading in occupation numbers                                        */

  if( (ibinary == 1) && (myid == 0) ){
     c_array1 = (char *) cmalloc(MAXWORD*sizeof(char ));
     c_array2 = (char *) cmalloc(MAXWORD*sizeof(char ));
  }/*endif*/

 if(myid==0){

   if(ibinary == 0){
     printf("Reading in the occupation numbers\n");
     readtoendofline(fp_dnameci);

     for(i=1;i<=nstate;i++){
      fscanf(fp_dnameci,"%lf %lf",&occ_up[i],&occ_dn[i]);
      readtoendofline(fp_dnameci);
      if(occ_up[i]>1 || occ_dn[i] > 1){
         printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
         printf("The occupation number of each state must be <= 1\n");
         printf("The %dth state is given as %g %g\n",i,occ_up[i],occ_dn[i]);
         printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
         fflush(stdout);exit(1);
      }/*endif*/
      if(cp_lda==1){occ_up[i] += occ_dn[i];}
     }/*endfor*/
   }else{

     printf("Reading in the occupation numbers\n");

     fread(c_array1,sizeof(char),MAXWORD,fp_dnameci);
     fread(c_array2,sizeof(char),MAXWORD,fp_dnameci);

 
     for(i=1;i<=nstate;i++){
       n = 1;  
       fread(&occ_up[i],sizeof(double),n,fp_dnameci);
       fread(&occ_dn[i],sizeof(double),n,fp_dnameci);
       if(occ_up[i]>1 || occ_dn[i] > 1){
         printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
         printf("The occupation number of each state must be <= 1\n");
         printf("The %dth state is given as %g %g\n",i,occ_up[i],occ_dn[i]);
         printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
         fflush(stdout);exit(1);
       }/*endif*/
       if(cp_lda==1){occ_up[i] += occ_dn[i];}
     }/*endfor*/

   }/*endif ibinary*/

   for(i=1;i<=nstate_up;i++){
    if(occ_up[i]==0){
      printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
      printf("The occupation number of each state must be > 0\n");
      printf("That is unoccupied states must be at the end\n");
      printf("of the state list.                              \n");
      printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
      fflush(stdout);exit(1);
    }/*endif*/
   }/*endfor*/
   if(cp_lsda==1){
    for(i=1;i<=nstate_dn;i++){
      if(occ_dn[i]==0){
       printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
       printf("The occupation number of each state must be > 0\n");
       printf("That is unoccupied states must be at the end\n");
       printf("of the state list.                              \n");
       printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
       fflush(stdout);exit(1);
      }/*endif*/
    }/*endfor*/
   }/*endif*/

   if(cp_norb >= 2){
     for(i=2; i<= nstate_up ; i++){
      if(occ_up[i] != occ_up[1]){
       printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
       printf("The occupation number of each state must be identical\n");
       printf("if you are using norb with no constraints or norm_only. \n");
       printf("Hasta la vista baby!!!!! \n");
       printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
       fflush(stdout);exit(1);
       }/*endif*/
    }/*endfor*/
    if(cp_lsda==1){
     for(i=2; i<= nstate_dn; i++){
       if(occ_dn[i] != occ_dn[1]){
        printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
        printf("The occupation number of each state must be identical\n");
        printf("if you are using no constraints. Hasta la vista baby!!!!! \n");
        printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
        fflush(stdout);exit(1);
       }/*endif*/
     }/*endfor*/
    }/*endif*/
   }/*endif*/

   iocc=0;
   for(i=1;i<=nstate_up;i++){
    for(j=1;j<=nstate_up;j++){
     iocc++;
     rocc_sum_up[iocc] = 1.0/(occ_up[i]+occ_up[j]);
    }/*endfor i*/
   }/* endfor j*/

   if(cp_lsda==1){
    iocc=0;
    for(i=1;i<=nstate_dn;i++){
     for(j=1;j<=nstate_dn;j++){
      iocc++;
      rocc_sum_dn[iocc] = 1.0/(occ_dn[i]+occ_dn[j]);
     }/*endfor i*/
    }/* endfor j*/
   }/*endif*/

 }/*endif for myid==0*/

 if( (ibinary == 1) && (myid == 0)){
     cfree(c_array1);
     cfree(c_array2);
 }/*endif*/


   if(nproc>1){Barrier(world);}

/*==========================================================================*/
/*  II) Communicate occupation variables                                    */

 if(nproc>1){

  Bcast(&(occ_up[1]),nstate,MPI_DOUBLE,0,world);
  if(cp_lsda==1){
    Bcast(&(occ_dn[1]),nstate,MPI_DOUBLE,0,world);
  }/*endif*/

  Bcast(&(rocc_sum_up[1]),nstate*nstate,MPI_DOUBLE,0,world);
  if(cp_lsda==1){
    Bcast(&(rocc_sum_dn[1]),nstate*nstate,MPI_DOUBLE,0,world);
  }/*endif*/

 }/* endif */


/*-----------------------------------------------------------------------*/
  }/* end routine */
/*==========================================================================*/





/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void read_coef_fetch_coefs(CP *cp, FILE *fp_dnameci, char *dnameci,
                          int initial_spread_opt)
              
/*==========================================================================*/
  {/*begin routine */
/*==========================================================================*/
/*  Local Variables */
#include "../typ_defs/typ_mask.h"

  int upper,pflag,ip,is,igo,i;
  double *cre_up,*cim_up,*cre_dn,*cim_dn;
  int  isoff,ipoff,iii;
  int ihaveit,iwantit;
  int myid_send,myid_recv;

/* Local Pointers */
  int myid               = cp->communicate.myid;
  int nproc              = cp->communicate.np;
  MPI_Comm world         = cp->communicate.world;
  int pi_beads           = cp->cpcoeffs_info.pi_beads;
  int nstate_up          = cp->cpcoeffs_info.nstate_up;
  int nstate_dn          = cp->cpcoeffs_info.nstate_dn;
  int istate_up_st       = cp->cpcoeffs_info.istate_up_st;
  int istate_up_end      = cp->cpcoeffs_info.istate_up_end;
  int istate_dn_st       = cp->cpcoeffs_info.istate_dn_st;
  int istate_dn_end      = cp->cpcoeffs_info.istate_dn_end;
  int pi_beads_st        = cp->cpcoeffs_info.pi_beads_proc_st;
  int pi_beads_end       = cp->cpcoeffs_info.pi_beads_proc_end;
  double *cre_up_tmp     = cp->cpscr.cpscr_wave.cre_up;
  double *cim_up_tmp     = cp->cpscr.cpscr_wave.cim_up;
  double *cre_dn_tmp     = cp->cpscr.cpscr_wave.cre_dn;
  double *cim_dn_tmp     = cp->cpscr.cpscr_wave.cim_dn;
  int ncoef              = cp->cpcoeffs_info.ncoef;
  int cp_lsda            = cp->cpopts.cp_lsda;
  int cp_lda             = cp->cpopts.cp_lda;

  int ibinary,n;
  char *c_array1,*c_array2;
  float cre_dum,cim_dum;

  ibinary = cp->cpopts.iread_coef_binary;

  if( (myid == 0) && (ibinary == 1)){
    c_array1 = (char *) cmalloc(MAXWORD*sizeof(char ));
    c_array2 = (char *) cmalloc(MAXWORD*sizeof(char ));
  }/*endif*/

/*==========================================================================*/
/* 0) Start reading */

  if(myid==0){
    printf("Reading in plane wave coefficients\n");
   if(ibinary == 0){
     readtoendofline(fp_dnameci);
   }else{
     fread(c_array1,sizeof(char),MAXWORD,fp_dnameci);
     fread(c_array2,sizeof(char),MAXWORD,fp_dnameci);
   }/*endif ibinary */
  }/*endif for myid==0*/
  if(nproc>1){Barrier(world);}

/*==========================================================================*/
/*  I) Read Up states                                                    */

  upper = pi_beads;
  pflag = 1;
  if(initial_spread_opt == 1){upper = 1;pflag=2;}

  fflush(stdout);  if(nproc>1){Barrier(world);}
  for(ip=1;ip<=upper;ip++){
    for(is=1;is<=nstate_up;is++){
       igo=0;
        if((is>=istate_up_st) && (is<=istate_up_end) && 
           (ip>=pi_beads_st)  && (ip<=pi_beads_end)){
          if(myid != 0){
            Ssend(&myid,1,MPI_INT,0,0,world);
          }/*endif*/
          igo=1;
       }/* endif */
       if(igo==0 && myid==0){
         Recv(&myid_recv,1,MPI_INT,MPI_ANY_SOURCE,MPI_ANY_TAG,world);
       }
       if(myid==0){
        for(i=1;i<=ncoef;i++){
         if(ibinary == 0){
          fscanf(fp_dnameci,"%lf %lf",&(cre_up_tmp[i]),&(cim_up_tmp[i]));
          readtoendofline(fp_dnameci);
         }else{
          n = 1;
          fread(&(cre_dum),sizeof(float),n,fp_dnameci);
          fread(&(cim_dum),sizeof(float),n,fp_dnameci);
            cre_up_tmp[i] = (double) cre_dum;
            cim_up_tmp[i] = (double) cim_dum;
          
         }/*endif ibinary*/
        }/*endfor: write*/
        if(igo==0){
          Ssend(&(cre_up_tmp[1]),ncoef,MPI_DOUBLE,myid_recv,myid_recv,world);
          Ssend(&(cim_up_tmp[1]),ncoef,MPI_DOUBLE,myid_recv,myid_recv,world);
        }/* endif igo*/
       } /* endif */
       if(igo==1){
         if(myid != 0){
          Recv(&(cre_up_tmp[1]),ncoef,MPI_DOUBLE,0,MPI_ANY_TAG,world);
          Recv(&(cim_up_tmp[1]),ncoef,MPI_DOUBLE,0,MPI_ANY_TAG,world);
         }/* endif */
         isoff=(is-istate_up_st)*ncoef;
         ipoff=(ip-pi_beads_st+1);
         cre_up = cp->cpcoeffs_pos[ipoff].cre_up;
         cim_up = cp->cpcoeffs_pos[ipoff].cim_up;
         for(i=1;i<=ncoef;i++){
           cre_up[(i+isoff)] = cre_up_tmp[i];
           cim_up[(i+isoff)] = cim_up_tmp[i];
         }/* endfor */
       }/* endif igo */
       if(nproc>1){Barrier(world);}
    }/*endfor:states*/
  }/* endfor ip*/


/*==========================================================================*/
/*  II) Read Dn states                                                    */

 if(cp_lsda == 1){

   if(myid==0){
     if(ibinary == 0){
     readtoendofline(fp_dnameci);
     }else{
     fread(c_array1,sizeof(char),MAXWORD,fp_dnameci);
     fread(c_array2,sizeof(char),MAXWORD,fp_dnameci);

     }/*endif ibinary*/
   }/*endif myid*/
   for(ip=1;ip<=upper;ip++){
     for(is=1;is<=nstate_dn;is++){
       igo=0;
       if((is>=istate_dn_st) && (is<=istate_dn_end) && 
          (ip>=pi_beads_st) && (ip<=pi_beads_end)){
        if(myid != 0){Ssend(&myid,1,MPI_INT,0,0,world);}
        igo=1;
       }/* endif */
       if(igo==0 && myid==0){
          Recv(&myid_recv,1,MPI_INT,MPI_ANY_SOURCE,MPI_ANY_TAG,world);
        }
       if(myid==0){
        for(i=1;i<=ncoef;i++){
       if(ibinary == 0){
         fscanf(fp_dnameci,"%lf %lf",&(cre_dn_tmp[i]),&(cim_dn_tmp[i]));
         readtoendofline(fp_dnameci);
       }else{
         n = 1;
         fread(&(cre_dum),sizeof(float),n,fp_dnameci);
         fread(&(cim_dum),sizeof(float),n,fp_dnameci);
         cre_dn_tmp[i] = (double) cre_dum;
         cim_dn_tmp[i] = (double) cim_dum;
       }/*endif ibinary*/
        }/*endfor: write*/
        if(igo==0){
          Ssend(&(cre_dn_tmp[1]),ncoef,MPI_DOUBLE,myid_recv,myid_recv,world);
          Ssend(&(cim_dn_tmp[1]),ncoef,MPI_DOUBLE,myid_recv,myid_recv,world);
        }/* endif igo*/
       } /* endif */
       if(igo==1){
         if(myid != 0){
          Recv(&(cre_dn_tmp[1]),ncoef,MPI_DOUBLE,0,MPI_ANY_TAG,world);
          Recv(&(cim_dn_tmp[1]),ncoef,MPI_DOUBLE,0,MPI_ANY_TAG,world);
         }/* endif */
         isoff=(is-istate_dn_st)*ncoef;
         ipoff=(ip-pi_beads_st+1);
         cre_dn = cp->cpcoeffs_pos[ipoff].cre_dn;
         cim_dn = cp->cpcoeffs_pos[ipoff].cim_dn;
         for(i=1;i<=ncoef;i++){
           cre_dn[(i+isoff)] = cre_dn_tmp[i];
           cim_dn[(i+isoff)] = cim_dn_tmp[i];
         }/* endfor */
       }/* endif igo */
       if(nproc>1){Barrier(world);}
     }/*endfor:states*/
    }/* endfor ip*/


 }/* endif cp_lsda */

/*==========================================================================*/
/* III) Spread Up states */

   if(upper < pi_beads){

    for(ip=2;ip<=pi_beads;ip++){
     for(is=1;is<=nstate_up;is++){
      ihaveit = 0; iwantit = 0;
      if(pi_beads_st==1 && ((is >= istate_up_st) && (is <= istate_up_end))){
        ihaveit = myid+1;
        isoff=(is-istate_up_st)*ncoef;
        cre_up = cp->cpcoeffs_pos[1].cre_up;
        cim_up = cp->cpcoeffs_pos[1].cim_up;
        for(i=1;i<=ncoef;i++){
          cre_up_tmp[i] = cre_up[(i+isoff)];
          cim_up_tmp[i] = cim_up[(i+isoff)];
        }/* endfor */
      }/* endif */
      if((is>=istate_up_st) && (is<=istate_up_end) && 
          (ip>=pi_beads_st) && (ip<=pi_beads_end)){
        iwantit = myid+1;
      }/* endif */
      myid_send=0; myid_recv=0;
      if(nproc>1){
       Allreduce(&ihaveit,&myid_send,1,MPI_INT,MPI_SUM,0,world);
       Allreduce(&iwantit,&myid_recv,1,MPI_INT,MPI_SUM,0,world);
       myid_send--; myid_recv--;
      } /* endif */
      if(myid_send != myid_recv){
       if(myid == myid_send){
        Ssend(&(cre_up_tmp[1]),ncoef,MPI_DOUBLE,myid_recv,myid_recv,world);
        Ssend(&(cim_up_tmp[1]),ncoef,MPI_DOUBLE,myid_recv,myid_recv,world);
       }
       if(myid == myid_recv){
         Recv(&(cre_up_tmp[1]),ncoef,MPI_DOUBLE,myid_send,MPI_ANY_TAG,world);
         Recv(&(cim_up_tmp[1]),ncoef,MPI_DOUBLE,myid_send,MPI_ANY_TAG,world);
       }
      } /* endif */
      if(iwantit==myid+1){
        isoff=(is-istate_up_st)*ncoef;
        ipoff=(ip-pi_beads_st+1);
        cre_up = cp->cpcoeffs_pos[ipoff].cre_up;
        cim_up = cp->cpcoeffs_pos[ipoff].cim_up;
        for(i=1;i<=ncoef;i++){
          cre_up[(i+isoff)] = cre_up_tmp[i];
          cim_up[(i+isoff)] = cim_up_tmp[i];
        }/* endfor */
      }/* endif iwantit*/
      if(nproc>1){Barrier(world);}
     }/* endfor is */
    }/* endfor ip */

   }/* endif upper<pi_beads*/

/*==========================================================================*/
/* III) Spread Dn states    */

   if(upper < pi_beads && cp_lsda == 1 && nstate_dn != 0){

     for(ip=2;ip<=pi_beads;ip++){
      for(is=1;is<=nstate_dn;is++){
       ihaveit = 0; iwantit = 0;
       if(pi_beads_st==1 && ((is >= istate_dn_st) && (is <= istate_dn_end))){
         ihaveit = myid+1;
         isoff=(is-istate_dn_st)*ncoef;
         cre_dn = cp->cpcoeffs_pos[1].cre_dn;
         cim_dn = cp->cpcoeffs_pos[1].cim_dn;
         for(i=1;i<=ncoef;i++){
           cre_dn_tmp[i] = cre_dn[(i+isoff)];
           cim_dn_tmp[i] = cim_dn[(i+isoff)];
         }/* endfor */
       }/* endif */
       if((is>=istate_dn_st) && (is<=istate_dn_end) && 
           (ip>=pi_beads_st) && (ip<=pi_beads_end)){
         iwantit = myid+1;
       }/* endif */
       myid_send=0; myid_recv=0;
       if(nproc>1){
        Allreduce(&myid_send,&ihaveit,1,MPI_INT,MPI_SUM,0,world);
        Allreduce(&myid_recv,&iwantit,1,MPI_INT,MPI_SUM,0,world);
        myid_send--; myid_recv--;
       } /* endif */
       if(myid_send != myid_recv){
        if(myid == myid_send){
          Ssend(&(cre_dn_tmp[1]),ncoef,MPI_DOUBLE,myid_recv,myid_recv,world);
          Ssend(&(cim_dn_tmp[1]),ncoef,MPI_DOUBLE,myid_recv,myid_recv,world);
        }
        if(myid == myid_recv){
          Recv(&(cre_dn_tmp[1]),ncoef,MPI_DOUBLE,myid_send,MPI_ANY_TAG,world);
          Recv(&(cim_dn_tmp[1]),ncoef,MPI_DOUBLE,myid_send,MPI_ANY_TAG,world);
        }
       } /* endif */
       if(iwantit==myid+1){
         isoff=(is-istate_dn_st)*ncoef;
         ipoff=(ip-pi_beads_st+1);
         cre_dn = cp->cpcoeffs_pos[ipoff].cre_dn;
         cim_dn = cp->cpcoeffs_pos[ipoff].cim_dn;
         for(i=1;i<=ncoef;i++){
           cre_dn[(i+isoff)] = cre_dn_tmp[i];
           cim_dn[(i+isoff)] = cim_dn_tmp[i];
         }/* endfor */
       }/* endif iwantit*/
       if(nproc>1){Barrier(world);}
      }/* endfor is */
     }/* endfor ip */
   
   }/* endif upper<pi_beads and lsda*/

/* free locally assigned memory */
  if( (myid == 0) && (ibinary == 1)){
    cfree(c_array1);
    cfree(c_array2);
  }/*endif*/

/*-----------------------------------------------------------------------*/
  } /* end routine */
/*==========================================================================*/






/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void read_coef_fetch_vcoefs(CP *cp, FILE *fp_dnameci, char *dnameci)
              
/*==========================================================================*/
  {/*begin routine */
/*==========================================================================*/
/*               Local variable declarations                                */
#include "../typ_defs/typ_mask.h"

  int upper,pflag,ip,is,igo,i;
  double *vcre_up,*vcim_up,*vcre_dn,*vcim_dn;
  int  isoff,ipoff;
  int ihaveit,iwantit;
  int myid_send,myid_recv;

/* Local Pointers */
  int myid               = cp->communicate.myid;
  int nproc              = cp->communicate.np;
  MPI_Comm world         = cp->communicate.world;
  int pi_beads           = cp->cpcoeffs_info.pi_beads;
  int nstate_up          = cp->cpcoeffs_info.nstate_up;
  int nstate_dn          = cp->cpcoeffs_info.nstate_dn;
  int istate_up_st       = cp->cpcoeffs_info.istate_up_st;
  int istate_up_end      = cp->cpcoeffs_info.istate_up_end;
  int istate_dn_st       = cp->cpcoeffs_info.istate_dn_st;
  int istate_dn_end      = cp->cpcoeffs_info.istate_dn_end;
  int pi_beads_st        = cp->cpcoeffs_info.pi_beads_proc_st;
  int pi_beads_end       = cp->cpcoeffs_info.pi_beads_proc_end;
  double *cre_up_tmp     = cp->cpscr.cpscr_wave.cre_up;
  double *cim_up_tmp     = cp->cpscr.cpscr_wave.cim_up;
  double *cre_dn_tmp     = cp->cpscr.cpscr_wave.cre_dn;
  double *cim_dn_tmp     = cp->cpscr.cpscr_wave.cim_dn;
  int ncoef              = cp->cpcoeffs_info.ncoef;
  int cp_lsda            = cp->cpopts.cp_lsda;
  int cp_lda             = cp->cpopts.cp_lda;

  char *c_array1,*c_array2;
  int n;  /* for binary read */
  int ibinary;

  float cre_dum,cim_dum;

  ibinary = cp->cpopts.iread_coef_binary;

  if( (myid == 0) && (ibinary == 1)){
    c_array1 = (char *) cmalloc(MAXWORD*sizeof(char ));
    c_array2 = (char *) cmalloc(MAXWORD*sizeof(char ));
  }/*endif*/

/*==========================================================================*/
/* 0) Announce to the screen */

   if(myid==0){
     printf("Reading in plane wave velocities\n");
     if(ibinary == 0){
      readtoendofline(fp_dnameci);
     }else{
       fread(c_array1,sizeof(char),MAXWORD,fp_dnameci);
       fread(c_array2,sizeof(char),MAXWORD,fp_dnameci);
     }/*endif ibinary */
   }/*endif for myid==0*/

/*==========================================================================*/
/*  I) Read up states velocities                                            */

   for(ip=1;ip<=pi_beads;ip++){
     if(nproc>1){Barrier(world);}
     for(is=1;is<=nstate_up;is++){
       igo=0;
       if((is>=istate_up_st) && (is<=istate_up_end) && 
          (ip>=pi_beads_st) && (ip<=pi_beads_end)){
        if(myid != 0){Ssend(&myid,1,MPI_INT,0,0,world);}
        igo=1;
       }/* endif */
       if(igo==0 && myid==0){
          Recv(&myid_recv,1,MPI_INT,MPI_ANY_SOURCE,MPI_ANY_TAG,world);
        }
       if(myid==0){
        for(i=1;i<=ncoef;i++){
       if(ibinary == 0){
         fscanf(fp_dnameci,"%lf %lf",&(cre_up_tmp[i]),&(cim_up_tmp[i]));
         readtoendofline(fp_dnameci);
       }else{
         n = 1;
         fread(&(cre_dum),sizeof(float),n,fp_dnameci);
         fread(&(cim_dum),sizeof(float),n,fp_dnameci);
          cre_up_tmp[i] = (double) cre_dum;
          cim_up_tmp[i] = (double) cim_dum;
       }/*endif ibinary */
        }/*endfor: write*/
        if(igo==0){
          Ssend(&(cre_up_tmp[1]),ncoef,MPI_DOUBLE,myid_recv,myid_recv,world);
          Ssend(&(cim_up_tmp[1]),ncoef,MPI_DOUBLE,myid_recv,myid_recv,world);
        }/* endif igo*/
       } /* endif */
       if(igo==1){
         if(myid != 0){
          Recv(&(cre_up_tmp[1]),ncoef,MPI_DOUBLE,0,MPI_ANY_TAG,world);
          Recv(&(cim_up_tmp[1]),ncoef,MPI_DOUBLE,0,MPI_ANY_TAG,world);
         }/* endif */
         isoff=(is-istate_up_st)*ncoef;
         ipoff=(ip-pi_beads_st+1);
         vcre_up =   cp->cpcoeffs_pos[ipoff].vcre_up;
         vcim_up =   cp->cpcoeffs_pos[ipoff].vcim_up;
         for(i=1;i<=ncoef;i++){
           vcre_up[(i+isoff)] = cre_up_tmp[i];
           vcim_up[(i+isoff)] = cim_up_tmp[i];
         }/* endfor */
       }/* endif igo */
       if(nproc>1){Barrier(world);}
     }/*endfor:states*/
    }/* endfor ip*/


/*-----------------------------------------------------------------------*/
/*  A) Read dn states                                                    */

 if(cp_lsda == 1){

   if(myid==0){
     if(ibinary == 0){
     readtoendofline(fp_dnameci);
     }else{
      fread(c_array1,sizeof(char),MAXWORD,fp_dnameci);
      fread(c_array2,sizeof(char),MAXWORD,fp_dnameci);
     }/*endif ibinary*/
   }/*endif myid*/
   for(ip=1;ip<=pi_beads;ip++){
     for(is=1;is<=nstate_dn;is++){
       igo=0;
       if((is>=istate_dn_st) && (is<=istate_dn_end) && 
          (ip>=pi_beads_st) && (ip<=pi_beads_end)){
        if(myid != 0){Ssend(&myid,1,MPI_INT,0,0,world);}
        igo=1;
       }/* endif */
       if(igo==0 && myid==0){
          Recv(&myid_recv,1,MPI_INT,MPI_ANY_SOURCE,MPI_ANY_TAG,world);
        }
       if(myid==0){
        for(i=1;i<=ncoef;i++){
        if(ibinary == 0){
         fscanf(fp_dnameci,"%lf %lf",&(cre_dn_tmp[i]),&(cim_dn_tmp[i]));
         readtoendofline(fp_dnameci);
        }else{
          n = 1;
          fread(&(cre_dum),sizeof(float),n,fp_dnameci);
          fread(&(cim_dum),sizeof(float),n,fp_dnameci);
          cre_dn_tmp[i] = (double) cre_dum;
          cim_dn_tmp[i] = (double) cim_dum;
        }/*endif ibinary*/
        }/*endfor: write*/
        if(igo==0){
          Ssend(&(cre_dn_tmp[1]),ncoef,MPI_DOUBLE,myid_recv,myid_recv,world);
          Ssend(&(cim_dn_tmp[1]),ncoef,MPI_DOUBLE,myid_recv,myid_recv,world);
        }/* endif igo*/
       } /* endif */
       if(igo==1){
         if(myid != 0){
          Recv(&(cre_dn_tmp[1]),ncoef,MPI_DOUBLE,0,MPI_ANY_TAG,world);
          Recv(&(cim_dn_tmp[1]),ncoef,MPI_DOUBLE,0,MPI_ANY_TAG,world);
         }/* endif */
         isoff=(is-istate_dn_st)*ncoef;
         ipoff=(ip-pi_beads_st+1);
         vcre_dn =   cp->cpcoeffs_pos[ipoff].vcre_dn;
         vcim_dn =   cp->cpcoeffs_pos[ipoff].vcim_dn;
         for(i=1;i<=ncoef;i++){
           vcre_dn[(i+isoff)] = cre_dn_tmp[i];
           vcim_dn[(i+isoff)] = cim_dn_tmp[i];
         }/* endfor */
       }/* endif igo */
       if(nproc>1){Barrier(world);}
     }/*endfor:states*/
    }/* endfor ip*/

  }/* endif lsda */

/* free locally assigned memory */
  if( (myid == 0) && (ibinary == 1)){
    c_array1 = (char *) cmalloc(MAXWORD*sizeof(char ));
    c_array2 = (char *) cmalloc(MAXWORD*sizeof(char ));
  }/*endif*/

/*-----------------------------------------------------------------------*/
  }/* end routine */
/*==========================================================================*/





/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void read_coef_fetch_vnhc(CP *cp, FILE *fp_dnameci, char *dnameci)
              
/*==========================================================================*/
  {/*begin routine */
/*==========================================================================*/
/*               Local variable declarations                                */
#include "../typ_defs/typ_mask.h"

  int ichain,inhc,ip,iproc,is;
  int nrecv,myid_recv;
  int num_c_nhc_now,len_c_nhc_now;
  int iwantit;
  int nstate_proc_now;
  int nstate_proc_lat;
  int ipoff,ioff_re,ioff_im;
  double **vc_nhc;

/* Local Pointers */
  int np_states        = cp->communicate.np_states;
  int myid             = cp->communicate.myid;
  int myid_st          = cp->communicate.myid_state;
  int nproc            = cp->communicate.np;
  MPI_Comm world       = cp->communicate.world;
  int pi_beads         = cp->cpcoeffs_info.pi_beads;
  int pi_beads_proc    = cp->cpcoeffs_info.pi_beads_proc;
  int ncoef            = cp->cpcoeffs_info.ncoef;
  int num_c_nhc        = cp->cptherm_info.num_c_nhc;
  int num_c_nhc_proc   = cp->cptherm_info.num_c_nhc_proc;
  int massiv_flag      = cp->cptherm_info.massiv_flag;
  int len_c_nhc        = cp->cptherm_info.len_c_nhc;
  int istate_nhc_opt   = cp->cptherm_info.istate_nhc_opt;
  int nstate_up        = cp->cpcoeffs_info.nstate_up;
  int nstate_dn        = cp->cpcoeffs_info.nstate_dn;
  int nstate_up_proc   = cp->cpcoeffs_info.nstate_up_proc;
  int nstate_dn_proc   = cp->cpcoeffs_info.nstate_dn_proc;
  int pi_beads_st      = cp->cpcoeffs_info.pi_beads_proc_st;
  int pi_beads_end     = cp->cpcoeffs_info.pi_beads_proc_end;
  int cp_lsda          = cp->cpopts.cp_lsda;
  int cp_lda           = cp->cpopts.cp_lda;
  double *cre_up_tmp   = cp->cpscr.cpscr_wave.cre_up;

  char *c_array1,*c_array2;
  int ibinary,n;

  float cre_dum,cim_dum;

  ibinary = cp->cpopts.iread_coef_binary;

  if( (myid == 0) && (ibinary == 1)){
     c_array1 = (char *) cmalloc(MAXWORD*sizeof(char));
     c_array2 = (char *) cmalloc(MAXWORD*sizeof(char));
  }/*endif*/

/*==========================================================================*/
/* 0) Announce to the screen and check the data */

  if(myid==0){

    printf("Reading in plane wave coefficients NHCs\n");    
    if(ibinary == 0){
     readtoendofline(fp_dnameci);
     fscanf(fp_dnameci,"%d %d\n",&num_c_nhc_now,&len_c_nhc_now);
     readtoendofline(fp_dnameci);
    }else{
     fread(c_array1,sizeof(char),MAXWORD,fp_dnameci);
     fread(c_array2,sizeof(char),MAXWORD,fp_dnameci);
      n = 1;
     fread(&num_c_nhc_now,sizeof(int),n,fp_dnameci);
     fread(&len_c_nhc_now,sizeof(int),n,fp_dnameci);
     fread(c_array1,sizeof(char),MAXWORD,fp_dnameci);
    }/*endif ibinary*/

    if(num_c_nhc_now != num_c_nhc) {
      printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
      printf("Mismatched number of coefficient Nose-Hoover chains\n");
      printf("%d vs %d\n",num_c_nhc,num_c_nhc_now);
      printf("in user specified plane wave coefficient file %s\n",dnameci);
      printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
      fflush(stdout);exit(1);
    } /* endif */

    if(len_c_nhc_now != len_c_nhc) {
      printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
      printf("Mismatched length of coefficient Nose-Hoover chains\n");
      printf("%d vs %d\n",len_c_nhc,len_c_nhc_now);
      printf("in user specified plane wave coefficient file %s\n",dnameci);
      printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
      fflush(stdout);exit(1);
    } /* endif */

  }/*endif for myid==0*/

/*==========================================================================*/
/* I) Initialize form flags and read stuff in                               */
   
  for(ip=1;ip<=pi_beads_proc;ip++){
    cp->cptherm_pos[ip].itherm_form_up = 0;
    cp->cptherm_pos[ip].itherm_form_dn = 0;
  }/*endfor*/
   
  for(ichain=1;ichain<=len_c_nhc;ichain++){
    for(ip=1;ip<=pi_beads;ip++){

  /*-----------------------------------------------------------------------*/
  /* i) Massive */
      if(istate_nhc_opt==4){

         nrecv = 2*ncoef;
         for(iproc=0;iproc<np_states;iproc++){
           iwantit=0;
           nstate_proc_now=0;
           nstate_proc_lat=0;
           if(ip>=pi_beads_st && ip <= pi_beads_end && myid_st == iproc){
             nstate_proc_now = nstate_up_proc;
             nstate_proc_lat = nstate_dn_proc;
             if(myid != 0){
               Ssend(&(myid),1,MPI_INT,0,0,world);
               Ssend(&(nstate_up_proc),1,MPI_INT,0,0,world);
               Ssend(&(nstate_dn_proc),1,MPI_INT,0,0,world);
             }/* endif */
             iwantit=1;
           }/* endif */
           if(myid==0){
             if(iwantit != 1){
              Recv(&(myid_recv),1,MPI_INT,MPI_ANY_SOURCE,MPI_ANY_TAG,world);
              Recv(&(nstate_proc_now),1,MPI_INT,MPI_ANY_SOURCE,MPI_ANY_TAG,
                     world);
              Recv(&(nstate_proc_lat),1,MPI_INT,MPI_ANY_SOURCE,MPI_ANY_TAG,
                     world);
             }/* endif */
           }/*endif*/
           for(is=1;is<=nstate_proc_now;is++){
             if(myid==0){
               for(inhc=1;inhc<=nrecv;inhc++){
                if(ibinary == 0){
                 fscanf(fp_dnameci,"%lf",&(cre_up_tmp[inhc]));
                 readtoendofline(fp_dnameci);
                }else{
                  n = 1;
                  fread(&(cre_dum),sizeof(float),n,fp_dnameci);
                  cre_up_tmp[inhc] = (double) cre_dum;

                }/*endif ibinary*/
               }/*endfor*/
               if(iwantit != 1){
                 Ssend(&(cre_up_tmp[1]),nrecv,MPI_DOUBLE,
                                         myid_recv,myid_recv,world);
               }/*endif*/
             }/*endif : myid==0*/
             if(iwantit==1){
               if(myid !=0){
                 Recv(&(cre_up_tmp[1]),nrecv,MPI_DOUBLE,
                                     MPI_ANY_SOURCE,MPI_ANY_TAG,world);
               }/* endif myid */
               ipoff=(ip-pi_beads_st+1);
               ioff_re = (is-1)*ncoef;
               ioff_im = (nstate_proc_now+is-1)*ncoef;
               vc_nhc  = cp->cptherm_pos[ipoff].vc_nhc;
               for(inhc=1;inhc<=ncoef;inhc++){
                 vc_nhc[ichain][inhc+ioff_re] = cre_up_tmp[inhc];
               }/* endfor */
               for(inhc=1;inhc<=ncoef;inhc++){
                 vc_nhc[ichain][inhc+ioff_im] =  cre_up_tmp[inhc+ncoef];
               }/* endfor */
             }/* endif iwantit*/
           }/*endfor : is*/

           if(cp_lsda==1){
             for(is=1;is<=nstate_proc_lat;is++){
               if(myid==0){
                 for(inhc=1;inhc<=nrecv;inhc++){
       	          if(ibinary == 0){
                   fscanf(fp_dnameci,"%lf",&(cre_up_tmp[inhc]));
                   readtoendofline(fp_dnameci);
		  }else{
                   n = 1;
                   fwrite(&(cre_up_tmp[inhc]),sizeof(double),n,fp_dnameci);
		  }/*endif ibinary*/
                 }/*endfor*/
                 if(iwantit != 1){
                   Ssend(&(cre_up_tmp[1]),nrecv,MPI_DOUBLE,
                                         myid_recv,myid_recv,world);
                 }/*endif*/
               }/*endif : myid==0*/
               if(iwantit==1){
                 if(myid !=0){
                   Recv(&(cre_up_tmp[1]),nrecv,MPI_DOUBLE,
                                       MPI_ANY_SOURCE,MPI_ANY_TAG,world);
                 }/* endif myid */
                 ipoff=(ip-pi_beads_st+1);
                 ioff_re = (is-1+nstate_up_proc*2)*ncoef;
                 ioff_im = (nstate_proc_now+is-1+nstate_up_proc*2)*ncoef;
                 vc_nhc  = cp->cptherm_pos[ipoff].vc_nhc;
                 for(inhc=1;inhc<=ncoef;inhc++){
                   vc_nhc[ichain][inhc+ioff_re] = cre_up_tmp[inhc];
                 }/* endfor */
                 for(inhc=1;inhc<=ncoef;inhc++){
                   vc_nhc[ichain][inhc+ioff_im] = cre_up_tmp[inhc+ncoef];
                 }/* endfor */
               }/* endif iwantit*/
              }/*endfor : is*/
	   }/*endif : lsda*/
           if(nproc>1){Barrier(world);}
         }/* endfor : iproc*/

      }/*endif : istate_nhc_opt==4*/
  /*-----------------------------------------------------------------------*/
  /* ii) Normal */
      if(istate_nhc_opt<=3){

         if(myid==0){
           for(inhc=1;inhc<=num_c_nhc;inhc++){
	     if(ibinary == 0){
              fscanf(fp_dnameci,"%lf",&(cre_up_tmp[inhc]));
              readtoendofline(fp_dnameci);
	     }else{
              n = 1;
              fread(&(cre_dum),sizeof(float),n,fp_dnameci);
               cre_up_tmp[inhc] = (double) cre_dum;
             }/*endif ibinary*/
           }/*endfor*/
         }/*endif:myid==0*/         
         if(nproc>1){Bcast(&cre_up_tmp[1],num_c_nhc,MPI_DOUBLE,0,world);}
         if(ip>=pi_beads_st && ip <= pi_beads_end && myid_st < np_states){
           ipoff=(ip-pi_beads_st+1);
           vc_nhc  = cp->cptherm_pos[ipoff].vc_nhc;
           for(inhc=1;inhc<=num_c_nhc;inhc++){
             vc_nhc[ichain][inhc] = cre_up_tmp[inhc];
           }/* endfor */
         }/*endif*/

      }/*endif : istate_nhc_opt==3*/
  /*-----------------------------------------------------------------------*/

    }/* endfor ip*/
  }/*endfor ichain*/

/*free locally assigned memory */
  if( (myid == 0) && (ibinary == 1)){
     cfree(c_array1);
     cfree(c_array2);
  }/*endif*/

/*-----------------------------------------------------------------------*/
  } /* end routine */
/*==========================================================================*/






/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void read_coef_init_nhc(CP *cp)
              
/*==========================================================================*/
  {/*begin routine */
/*==========================================================================*/
/*               Local variable declarations                                */
#include "../typ_defs/typ_mask.h"

 int ichain,inhc,ip;
 double **c_nhc;

/* Local pointers */
 int pi_beads         = cp->cpcoeffs_info.pi_beads;
 int pi_beads_proc    = cp->cpcoeffs_info.pi_beads_proc;
 int num_c_nhc        = cp->cptherm_info.num_c_nhc;
 int num_c_nhc_proc   = cp->cptherm_info.num_c_nhc_proc;
 int massiv_flag      = cp->cptherm_info.massiv_flag;
 int len_c_nhc        = cp->cptherm_info.len_c_nhc;

/*==========================================================================*/
/* I) standard thermos */

    if(massiv_flag == 0 && num_c_nhc >0){

     for(ip=1;ip<=pi_beads_proc;ip++){
      c_nhc = cp->cptherm_pos[ip].c_nhc;
      for(ichain=1;ichain<=len_c_nhc;ichain++){
       for(inhc=1;inhc<=num_c_nhc_proc;inhc++){
         c_nhc[ichain][inhc] = 0.0;
       } /*endfor*/
      } /*endfor*/
     } /*endfor*/

    }/*endif*/

/*==========================================================================*/
/* II) Massive thermos */

   if(massiv_flag != 0 && num_c_nhc >0){

     for(ip=1;ip<=pi_beads_proc;ip++){
        cp->cptherm_pos[ip].c_nhc_massiv = 0.0;
     }/*endfor*/

    }/*endif*/

/*-----------------------------------------------------------------------*/
  } /* end routine */
/*==========================================================================*/





/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void read_coef_transpose(CP *cp,int istart,int cp_wave_min,int cp_min)
              
/*==========================================================================*/
  {/*begin routine */
/*==========================================================================*/
/*               Local variable declarations                                */
#include "../typ_defs/typ_mask.h"

  int ip;

/* Local pointers */
  int pi_beads_proc    = cp->cpcoeffs_info.pi_beads_proc;

/*==========================================================================*/


  if((istart == 4) && ((cp_wave_min+cp_min)==0)) {
     control_coef_transpose_fwd(cp,3);
  }/* endif */

  if((istart == 3) &&((cp_wave_min+cp_min)==0)) {
     control_coef_transpose_fwd(cp,2);
  }/* endif */

  if((istart <= 2) || ((cp_wave_min+cp_min)==1)){
     control_coef_transpose_fwd(cp,1);
     for(ip=1;ip<=pi_beads_proc;ip++){
         cp->cpcoeffs_pos[ip].ivcoef_form_up = 1;
         cp->cpcoeffs_pos[ip].ivcoef_form_dn = 1;
     }/*endfor:ip*/
  }/*endif*/

  if(istart <= 3){
     for(ip=1;ip<=pi_beads_proc;ip++){
        cp->cptherm_pos[ip].itherm_form_up = 0;
        cp->cptherm_pos[ip].itherm_form_dn = 0;
     }/*endfor*/
  }/*endif*/

/*-----------------------------------------------------------------------*/
   }/* end routine */
/*==========================================================================*/






