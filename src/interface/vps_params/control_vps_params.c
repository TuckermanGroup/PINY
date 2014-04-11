/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
/*                                                                          */
/*                         PI_MD:                                           */
/*             The future of simulation technology                          */
/*             ------------------------------------                         */
/*                     Module: control_vps_params                           */
/*                                                                          */
/* This reads in and sets up the electron-atom interaction pseudopotential  */ 
/*                                                                          */
/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/


#include "standard_include.h"
#include "../typ_defs/typedefs_gen.h"
#include "../typ_defs/typedefs_cp.h"
#include "../typ_defs/typedefs_par.h"
#include "../typ_defs/typedefs_class.h"
#include "../typ_defs/typedefs_bnd.h"
#include "../proto_defs/proto_vps_params_entry.h"
#include "../proto_defs/proto_search_entry.h"
#include "../proto_defs/proto_intra_params_local.h"
#include "../proto_defs/proto_vps_params_local.h"
#include "../proto_defs/proto_handle_entry.h"
#include "../proto_defs/proto_friend_lib_entry.h"
#include "../proto_defs/proto_communicate_wrappers.h"
#include "../proto_defs/proto_math.h"

#define DEBUG_DKNY_OFF

#define JUERG_FACTOR_ON
#ifdef  JUERG_FACTOR_ON
#define JUERG_FACTOR 0.72
#else
#define JUERG_FACTOR 1.0
#endif

typedef struct vps_file{
 char name[MAXWORD];
}VPS_FILE;

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void control_vps_params(PSEUDO *pseudo,CELL *cell,
                        FILENAME_PARSE *filename_parse,
                        SPLINE_PARSE *spline_parse,int natm_typ,NAME *atm_typ,
                        double *tot_memory,int natm_tot,int natm_ab_init,
                        int cp_ptens_calc,
                        int cp_dual_grid_opt,
                        COMMUNICATE *communicate,double ecut_cp)

/*==========================================================================*/
/*               Begin subprogram:                                          */
      {/*begin routine*/
/*==========================================================================*/
/*               Local variable declarations                                */
#include "../typ_defs/typ_mask.h"

   int i;                       /* Num:  For loop counters             */
   int ifound;                  /* Num:  Data base match flag          */
   int ishift,ishift2,ishift3;  /* Num:  Angular momentum shifts       */
   int ngh_now;                 /* Num:  Number of ngh points read in  */
   int ngh_max;                 /* Num:  Max number of ngh points      */
   VPS_FILE *vps_file;          /* Fle:  Pseudopotential file          */
   double now_mem;              /* Num:  Current memory usage          */

  /* Dictionary memory */
   char *filename;              /* Char: temp file name                */
   CVPS *cvps_typ;
   char *fun_key;
   DICT_WORD *word;
   DICT_WORD *fun_dict;
   int num_fun_dict;
   DICT_WORD *vps_dict,*vps_dict_tmp;
   int num_vps_dict,ifirst,iii;
   int natm_typ_mall,natm_mall,nsplin_mall,norm_mall,nlist_mall;
   double dummy0,dummy1,dummy2,dummy3;
   int nmall_gh,natm_typ_gh = 0;   /* starting malloc value for gauss-hermite */
  /* Volume and grid variables */
   double alpha_conv_dual = pseudo->alpha_conv_dual;
   double vol_cp   = cell->vol_cp;

  /* Parallel processing */
   MPI_Comm comm = communicate->world;
   int myid      = communicate->myid;
   int num_proc  = communicate->np;

/*==========================================================================*/
/* 0) Output                                                                */

  if(myid==0){
   putchar('\n');
   PRINT_LINE_STAR
   printf("Searching the data bases both user defined and default\n");
   printf("for the electron-atom pseudopotentials\n");
   PRINT_LINE_DASH;printf("\n");
  }/*endif*/


/*==========================================================================*/
/* 0.5) Toast bad spline points : (note conversion back to Ry)              */

  if( (pseudo->nsplin_g< 4000) && (2.0*ecut_cp > 60.0)){
     if(myid==0){
       printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
       printf("Now, dude, lets clean up those files and use a \n");
       printf("reasonable number of psuedo spline points at large \n");
       printf("cutoffs. Its clear, %d points at %g Ry, won't do!\n",
               pseudo->nsplin_g,2.0*ecut_cp);
       printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
       fflush(stdout);
     }/*endif*/
     exit(1);
  }/*endif*/

/*==========================================================================*/
/* I) Convert alpha_conv_dual and Allocate some temporary character arrays  */

   alpha_conv_dual        /= (pow(vol_cp,1.0/3.0));
   pseudo->alpha_conv_dual = alpha_conv_dual;

  if(myid==0){
   filename_parse->vps_name = (NAME *) cmalloc(natm_typ*sizeof(NAME))-1;
   vps_file  = (VPS_FILE *) cmalloc(natm_typ*sizeof(VPS_FILE))-1;
   fun_key   = (char *)cmalloc(MAXWORD*sizeof(char));  
   filename  = (char *)cmalloc(MAXWORD*sizeof(char));  
   word      = (DICT_WORD *)cmalloc(sizeof(DICT_WORD))-1;  
   cvps_typ  = (CVPS *)cmalloc(sizeof(CVPS));  
   ifirst    = 1;
   set_potfun_dict(&fun_dict,&num_fun_dict,ifirst);
   set_potvps_dict(&vps_dict,&num_vps_dict,ifirst);      
   set_potvps_dict(&vps_dict_tmp,&num_vps_dict,ifirst);      
 }/*endif*/


/*==========================================================================*/
/* III) Malloc up the vps stuff                                             */ 

   natm_typ_mall      = natm_typ;

   if( (natm_typ_mall % 2)==0){natm_typ_mall++;}

   pseudo->n_ang      = (int *) cmalloc(natm_typ_mall*sizeof(int))-1;
   pseudo->loc_opt    = (int *) cmalloc(natm_typ_mall*sizeof(int))-1;
   pseudo->ivps_label = (int *) cmalloc(natm_typ_mall*sizeof(int))-1;
   pseudo->rcut_nl    = (double *) cmalloc(natm_typ_mall*sizeof(double))-1;
   pseudo->q_pseud    = (double *) cmalloc(natm_typ_mall*sizeof(double))-1;

   pseudo->nrad_0 = (int *) cmalloc(natm_typ_mall*sizeof(double))-1;
   pseudo->nrad_1 = (int *) cmalloc(natm_typ_mall*sizeof(double))-1;
   pseudo->nrad_2 = (int *) cmalloc(natm_typ_mall*sizeof(double))-1;
   pseudo->nrad_3 = (int *) cmalloc(natm_typ_mall*sizeof(double))-1;
  
   now_mem = (natm_typ_mall*(sizeof(double)*2 + sizeof(int))*7)*1.e-06;
  *tot_memory += now_mem;

/*==========================================================================*/
/* III) Loop over all unique atom types and get the vps stuff               */ 

 if(myid==0){
   pseudo->n_rad_max    = 0;
   pseudo->n_ang_max    = 0;
   pseudo->n_ang_max_kb = 0;
   pseudo->n_ang_max_gh = 0;
   ngh_max = 0;
   for(i=1;i<=natm_typ;i++) {
/*--------------------------------------------------------------------------*/
/*     A) First search the user defined data base                           */

      ifound = 0;
      strcpy(cvps_typ->atm1,atm_typ[i]);
      if(strcasecmp(filename_parse->user_vps_name,"")!=0) {
         search_base_vps(filename_parse->user_vps_name,
                         cvps_typ,fun_dict,num_fun_dict,
                         &vps_dict_tmp,vps_dict,num_vps_dict,&ifound);
         if(ifound==1){
            set_vps_params(vps_dict,
                           filename_parse->user_vps_name,fun_key,
                           &(pseudo->ivps_label[i]),filename,
                           &(pseudo->loc_opt[i]),&(pseudo->n_ang[i]),
                           &(pseudo->rcut_nl[i]),&ngh_now,
                           &(pseudo->nrad_0[i]),
                           &(pseudo->nrad_1[i]),&(pseudo->nrad_2[i]),
                           &(pseudo->nrad_3[i]));
         }/*endif*/
     }/*endif*/
/*--------------------------------------------------------------------------*/
/*     B) If you haven't found it search the default data base              */

      if(ifound == 0) {
         search_base_vps(filename_parse->def_vps_name,
                         cvps_typ,fun_dict,num_fun_dict,
                         &vps_dict_tmp,vps_dict,num_vps_dict,&ifound);
         if(ifound==1){
            set_vps_params(vps_dict,
                           filename_parse->def_vps_name,fun_key,
                           &(pseudo->ivps_label[i]),filename,
                           &(pseudo->loc_opt[i]),&(pseudo->n_ang[i]),
                           &(pseudo->rcut_nl[i]),&ngh_now,
                           &(pseudo->nrad_0[i]),
                           &(pseudo->nrad_1[i]),&(pseudo->nrad_2[i]),
                           &(pseudo->nrad_3[i]));
         }/*endif*/
      }/*endif*/
/*--------------------------------------------------------------------------*/
/*     C) Make sure you have now found this puppy, if not exit              */

      if(ifound == 0) {
         printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
         printf("Electron pseudopotential interaction with\n"); 
         printf("%s\n",atm_typ[i]);
         printf("not found in default interaction data base\n");
         printf("pi_md.vps\n");
         if(strlen(filename_parse->user_vps_name) > 0)  {
             printf("or in user defined pseudopot data base\n");
             printf("%s\n",filename_parse->user_vps_name);
         /*endif*/}
        putchar('\n');
        printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
        fflush(stdout);
        exit(1);
     }/*endif*/

/*--------------------------------------------------------------------------*/
/*     D) Find maximum angular momentum component                           */
      pseudo->n_ang_max = (pseudo->n_ang_max > pseudo->n_ang[i] ? 
                           pseudo->n_ang_max : pseudo->n_ang[i]);

     if(pseudo->ivps_label[i] != 2){
      pseudo->n_ang_max_kb = (pseudo->n_ang_max_kb > pseudo->n_ang[i] ? 
                              pseudo->n_ang_max_kb : pseudo->n_ang[i]);
     }else{
      pseudo->n_ang_max_gh = (pseudo->n_ang_max_gh > pseudo->n_ang[i] ? 
                              pseudo->n_ang_max_gh : pseudo->n_ang[i]);
      natm_typ_gh++;
     }/*endif*/

      pseudo->n_rad_max = (pseudo->n_rad_max > pseudo->nrad_0[i] ? 
                           pseudo->n_rad_max : pseudo->nrad_0[i]);
      pseudo->n_rad_max = (pseudo->n_rad_max > pseudo->nrad_1[i] ?
                           pseudo->n_rad_max : pseudo->nrad_1[i]);
      pseudo->n_rad_max = (pseudo->n_rad_max > pseudo->nrad_2[i] ?
                           pseudo->n_rad_max : pseudo->nrad_2[i]);
      pseudo->n_rad_max = (pseudo->n_rad_max > pseudo->nrad_3[i] ?
                           pseudo->n_rad_max : pseudo->nrad_3[i]);
      ngh_max           = MAX(ngh_max,ngh_now);

      strcpy(vps_file[i].name,filename);
      strcpy(filename_parse->vps_name[i],filename);
   }/*endfor natm_typ*/
      pseudo->ngh = ngh_max;
      pseudo->natm_typ_gh  = natm_typ_gh;
 }/*endif : myid==0*/

 if(num_proc>1){
    Bcast(&(pseudo->ivps_label[1]),natm_typ,MPI_INT,0,comm);
    Bcast(&(pseudo->loc_opt[1]),natm_typ,MPI_INT,0,comm);
    Bcast(&(pseudo->n_ang[1]),natm_typ,MPI_INT,0,comm);
    Bcast(&(pseudo->rcut_nl[1]),natm_typ,MPI_DOUBLE,0,comm);
    Bcast(&(pseudo->nrad_0[1]),natm_typ,MPI_INT,0,comm);
    Bcast(&(pseudo->nrad_1[1]),natm_typ,MPI_INT,0,comm);
    Bcast(&(pseudo->nrad_2[1]),natm_typ,MPI_INT,0,comm);
    Bcast(&(pseudo->nrad_3[1]),natm_typ,MPI_INT,0,comm);

    Bcast(&(pseudo->n_rad_max),1,MPI_INT,0,comm);
    Bcast(&(pseudo->n_ang_max),1,MPI_INT,0,comm);
    Bcast(&(pseudo->n_ang_max_kb),1,MPI_INT,0,comm);
    Bcast(&(pseudo->n_ang_max_gh),1,MPI_INT,0,comm);

    Bcast(&(pseudo->ngh),1,MPI_INT,0,comm);
    Bcast(&(pseudo->natm_typ_gh),1,MPI_INT,0,comm);
    Bcast(&(ngh_max),1,MPI_INT,0,comm);

    natm_typ_gh = pseudo->natm_typ_gh;

 }/*endif*/

/*==========================================================================*/
/*  III) Allocate more memory for pseudopotentials                          */

/*--------------------------------------------------------------------------*/
/* i) Malloc Pseudo spline and other stuff */
   pseudo->nsplin_g_tot = (pseudo->n_ang_max+1)*(pseudo->n_rad_max)
                         *(pseudo->nsplin_g)*natm_typ;
   nsplin_mall = pseudo->nsplin_g_tot;
   norm_mall   = (pseudo->n_ang_max+1)*natm_typ*
                 (pseudo->n_rad_max)*(pseudo->n_rad_max);
   if((nsplin_mall % 2)==0){nsplin_mall++;}
   if((norm_mall % 2)==0){norm_mall++;}
   pseudo->nsplin_g_mall = nsplin_mall;
   pseudo->nvpsnorm_mall = norm_mall;

   pseudo->vps0    = (double *) cmalloc(nsplin_mall*sizeof(double))-1;
   pseudo->vps1    = (double *) cmalloc(nsplin_mall*sizeof(double))-1;
   pseudo->vps2    = (double *) cmalloc(nsplin_mall*sizeof(double))-1;
   pseudo->vps3    = (double *) cmalloc(nsplin_mall*sizeof(double))-1;
 
   now_mem   = ( nsplin_mall*4 *sizeof(double))*1.e-06;
  *tot_memory += now_mem;

   if(cp_ptens_calc == 1){
     pseudo->dvps0    = (double *) cmalloc(nsplin_mall*sizeof(double))-1;
     pseudo->dvps1    = (double *) cmalloc(nsplin_mall*sizeof(double))-1;
     pseudo->dvps2    = (double *) cmalloc(nsplin_mall*sizeof(double))-1;
     pseudo->dvps3    = (double *) cmalloc(nsplin_mall*sizeof(double))-1;
     now_mem   = ( nsplin_mall*4 *sizeof(double))*1.e-06;
    *tot_memory += now_mem;
   }/* endif */
   pseudo->vpsnorm = (double *) cmalloc(norm_mall*sizeof(double))-1;
   pseudo->gzvps   = (double *) cmalloc(natm_typ_mall*sizeof(double))-1;

   pseudo->gzvps0  = (double *) cmalloc(natm_typ_mall*pseudo->n_rad_max*
                                        sizeof(double))-1;
   pseudo->nrad_max_l = (int *)cmalloc(natm_typ_mall*5*sizeof(int))-1;

   now_mem = ((norm_mall+natm_typ_mall+natm_typ_mall*pseudo->n_rad_max)*sizeof(double)
             +(natm_typ_mall*5)*sizeof(int))*1.e-06;

   *tot_memory += now_mem;

/*--------------------------------------------------------------------------*/
/* iii) Pseudo list : MAJOR HACKET JOB */

   natm_mall = natm_tot;
   if((natm_mall % 2)==0){natm_mall++;}
   nlist_mall = (pseudo->n_ang_max+1)*natm_tot;
   if((nlist_mall % 2)==0){nlist_mall++;}
   pseudo->n_ang_mall    = pseudo->n_ang_max;
   pseudo->nlist_mall    = nlist_mall;

   pseudo->x0w   = (double *) cmalloc(natm_mall*sizeof(double))-1;
   pseudo->y0w   = (double *) cmalloc(natm_mall*sizeof(double))-1;
   pseudo->z0w   = (double *) cmalloc(natm_mall*sizeof(double))-1;

   pseudo->np_nl        = (int *) cmalloc((pseudo->n_ang_max+1)*sizeof(int))-1;
   pseudo->np_nl_gh     = (int *) cmalloc((pseudo->n_ang_max+1)*sizeof(int))-1;
   pseudo->ip_nl        = (int *) cmalloc(nlist_mall*sizeof(int))-1;
   pseudo->ip_nl_gh     = (int *) cmalloc(nlist_mall*sizeof(int))-1;
   pseudo->ip_nl_rev    = (int *) cmalloc(nlist_mall*sizeof(int))-1;
   pseudo->ip_nl_rev_gh = (int *) cmalloc(nlist_mall*sizeof(int))-1;

   pseudo->map_nl = (int *) cmalloc(natm_mall*sizeof(int))-1;
   pseudo->ip_loc_cp_box = (int *) cmalloc(natm_mall *sizeof(int))-1;

   pseudo->np_nl_rad_str  = cmall_int_mat(1,pseudo->n_ang_max+1,
                                          1,pseudo->n_rad_max);
   pseudo->np_nl_rad_end  = cmall_int_mat(1,pseudo->n_ang_max+1,
                                          1,pseudo->n_rad_max);

   pseudo->rgh = (double *)cmalloc(ngh_max*sizeof(double))-1;

   nmall_gh =  ngh_max*(pseudo->natm_typ_gh)*(pseudo->n_ang_max+1);

   pseudo->wgh = (double *)cmalloc((nmall_gh)*sizeof(double))-1;

  for(i=1; i<= nmall_gh; i++){
   pseudo->wgh[i] = 0.0;
  }/*endfor*/ 

   now_mem = ((natm_mall*3 + 2*ngh_max)*sizeof(double)
            + nmall_gh*sizeof(double)
            +((pseudo->n_ang_max+1)+nlist_mall+natm_mall*2
            + 2*(pseudo->n_ang_max+1)*pseudo->n_rad_max )
             *sizeof(int))*1.e-06;

   *tot_memory += now_mem;

/*--------------------------------------------------------------------------*/
/* iv) Output */

   if(myid==0){
    printf("Pseudopotential allocation: %g Mbytes; Total memory: %g Mbytes\n",
           now_mem,*tot_memory);
   }/*endif*/

/*==========================================================================*/
/*  IV) Spline up the stuff                                                 */


   for(i=1;i<=natm_typ;i++){
     pseudo->nrad_max_l[1] = MAX(pseudo->nrad_max_l[1],pseudo->nrad_0[i]);
     pseudo->nrad_max_l[2] = MAX(pseudo->nrad_max_l[2],pseudo->nrad_1[i]);
     pseudo->nrad_max_l[3] = MAX(pseudo->nrad_max_l[3],pseudo->nrad_2[i]);
     pseudo->nrad_max_l[4] = MAX(pseudo->nrad_max_l[4],pseudo->nrad_3[i]);
   }/*endfor*/

   pseudo->dg_spl = ((pseudo->gmax_spl)-(pseudo->gmin_spl))
                    /((double)(pseudo->nsplin_g));

   for(i=1;i<=natm_typ*(pseudo->n_rad_max);i++){
     (pseudo->gzvps0)[i] = 0;
   }/*endfor*/

   for(i=1;i<=natm_typ;i++) {
      ishift  = (i-1)*(pseudo->n_ang_max+1)*(pseudo->nsplin_g)
                     *(pseudo->n_rad_max);
      ishift2 = (i-1)*(pseudo->n_ang_max+1)*(pseudo->n_rad_max)
                     *(pseudo->n_rad_max);
      ishift3 = (i-1)*(pseudo->n_rad_max);
      if(myid==0){strcpy(filename,vps_file[i].name);}
      if(cp_ptens_calc == 1){
          make_vps_splin(filename,pseudo->loc_opt[i],pseudo->n_ang[i],
                         pseudo->ivps_label[i],
                         pseudo->nsplin_g,pseudo->dg_spl,
                         pseudo->gmin_spl,
                         pseudo->gmax_spl,pseudo->gmin_true,
                         &(pseudo->vps0)[ishift],&(pseudo->vps1)[ishift],
                         &(pseudo->vps2)[ishift],&(pseudo->vps3)[ishift],
                         &(pseudo->dvps0)[ishift],&(pseudo->dvps1)[ishift],
                         &(pseudo->dvps2)[ishift],&(pseudo->dvps3)[ishift],
                         &(pseudo->gzvps)[i],&(pseudo->gzvps0)[ishift3],
                         &(pseudo->q_pseud[i]),
                         &(pseudo->vpsnorm)[ishift2],
                         (pseudo->nrad_0[i]),
                         (pseudo->nrad_1[i]),(pseudo->nrad_2[i]),
                         (pseudo->nrad_3[i]),
                         cp_ptens_calc,myid,comm,num_proc,
                         cp_dual_grid_opt,alpha_conv_dual,pseudo->n_rad_max,
                         pseudo->n_ang_max_gh,
                         &(pseudo->ngh),pseudo->rgh,pseudo->wgh);
      } else {
          make_vps_splin(filename,pseudo->loc_opt[i],pseudo->n_ang[i],
                         pseudo->ivps_label[i],
                         pseudo->nsplin_g,pseudo->dg_spl,
                         pseudo->gmin_spl,
                         pseudo->gmax_spl,pseudo->gmin_true,
                         &(pseudo->vps0)[ishift],&(pseudo->vps1)[ishift],
                         &(pseudo->vps2)[ishift],&(pseudo->vps3)[ishift],
                         &dummy0,&dummy1,
                         &dummy2,&dummy3,
                         &(pseudo->gzvps)[i],&(pseudo->gzvps0)[ishift3],
                         &(pseudo->q_pseud[i]),
                         &(pseudo->vpsnorm)[ishift2],
                         (pseudo->nrad_0[i]),
                         (pseudo->nrad_1[i]),(pseudo->nrad_2[i]),
                         (pseudo->nrad_3[i]),
                         cp_ptens_calc,myid,comm,num_proc,
                         cp_dual_grid_opt,alpha_conv_dual,pseudo->n_rad_max,
                         pseudo->n_ang_max_gh,
                         &(pseudo->ngh),pseudo->rgh,pseudo->wgh);
      }/* endif */
   }/*endfor*/


   if(myid==0){
      printf("Electron-atom interactions succesfully assigned\n");
   }/*endif*/

/*==========================================================================*/
/*  V) Free                                                 */

  if(myid==0){
    cfree(&vps_file[1]);
    cfree(fun_key);
    cfree(filename);
    cfree(&word[1]);
    cfree(&fun_dict[1]);
    cfree(&vps_dict[1]);
    cfree(&vps_dict_tmp[1]);
    cfree(cvps_typ);
  }/*endif*/

/*==========================================================================*/
/*  VI) Output                                                              */

  if(myid==0){
    printf("\n");PRINT_LINE_DASH
    printf("Completed the pseudopotential data bases searches\n");
    PRINT_LINE_STAR;putchar('\n');
 }/*endif*/

/*--------------------------------------------------------------------------*/
  }/*end routine*/
/*==========================================================================*/






/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void set_vps_params(DICT_WORD vps_dict[],char *filename,char *fun_key,
                    int *ivps_label,char *vps_file,
                    int *loc_opt,int *n_ang,double *rcut_nl,int *pngh,
                    int *nrad_0, int *nrad_1, int *nrad_2, int *nrad_3)

/*==========================================================================*/
/*               Begin subprogram:                                          */
{/*begin routine*/
int index,ierr,ngh;
double tmp;
/*--------------------------------------------------------------------------*/
/*  0) Set up                                                               */
         strcpy(fun_key,"vps_parm");
/*--------------------------------------------------------------------------*/
/*  1) Assign vps file                                                      */
      strcpy(vps_file,vps_dict[3].keyarg);
/*--------------------------------------------------------------------------*/
/*  2) Assign radial channel info                                           */

       sscanf(vps_dict[7].keyarg,"%lg",&tmp);
       (*nrad_0)  = (int )(tmp);
       sscanf(vps_dict[8].keyarg,"%lg",&tmp);
       (*nrad_1)  = (int )(tmp);
       sscanf(vps_dict[9].keyarg,"%lg",&tmp);
       (*nrad_2)  = (int )(tmp);
       sscanf(vps_dict[10].keyarg,"%lg",&tmp);
       (*nrad_3)  = (int )(tmp);

       if((*nrad_0)<0){index=7;
                 keyarg_barf(vps_dict,filename,fun_key,index);}
       if((*nrad_1)<0){index=7;
                 keyarg_barf(vps_dict,filename,fun_key,index);}
       if((*nrad_2)<0){index=7;
                 keyarg_barf(vps_dict,filename,fun_key,index);}
       if((*nrad_3)<0){index=7;
                 keyarg_barf(vps_dict,filename,fun_key,index);}

/*--------------------------------------------------------------------------*/
/*  4) Assign vps type                                                      */

      (*ivps_label) = -1;
      if(strcasecmp(vps_dict[2].keyarg,"loc") == 0)            {(*ivps_label) = 0;}
      if(strcasecmp(vps_dict[2].keyarg,"kb")  == 0)            {(*ivps_label) = 1;}
      if(strcasecmp(vps_dict[2].keyarg,"gauss_hermite")  == 0) {(*ivps_label) = 2;}
      if(strcasecmp(vps_dict[2].keyarg,"vdb") == 0)            {(*ivps_label) = 3;}
      if(strcasecmp(vps_dict[2].keyarg,"null")== 0)            {(*ivps_label) = 4;}
      if(strcasecmp(vps_dict[2].keyarg,"goedecker")== 0)       {(*ivps_label) = 5;}
      if((*ivps_label)==-1){index=2;
                 keyarg_barf(vps_dict,filename,fun_key,index);}

/*--------------------------------------------------------------------------*/
/*  2) Assign non-local info                                                */

       sscanf(vps_dict[4].keyarg,"%lg",&tmp);
       (*n_ang)  = (int )(tmp);
       if((*n_ang)<0){index=4;
                      keyarg_barf(vps_dict,filename,fun_key,index);}
       if((*n_ang)>4){index=4;
                      keyarg_barf(vps_dict,filename,fun_key,index);}

       sscanf(vps_dict[6].keyarg,"%lg",&tmp);
       (*rcut_nl) = tmp;
       if((*rcut_nl)<0){index=6;
                 keyarg_barf(vps_dict,filename,fun_key,index);}

/*--------------------------------------------------------------------------*/
/*  4) Local option                                                        */

    sscanf(vps_dict[5].keyarg,"%lg",&tmp);
    (*loc_opt) = (int )(tmp);

    if((*ivps_label)!=5){
      if(((*loc_opt)<0)||((*loc_opt)>(*n_ang))){index=5;
                 keyarg_barf(vps_dict,filename,fun_key,index);}
    }/*endif*/

    if(*ivps_label == 5){
     if( (*loc_opt)!=(*n_ang+1) ){
      if( ((*loc_opt)==0) && ((*nrad_0)==0) ){
       index=1;
      }else{
       printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
       printf("Dude, the local option must be set to \n"); 
       printf("the number of available angular momentum channels+1\n");
       printf("%d vs %d \n",(*loc_opt),(*n_ang+1));
       printf("for a Goedecker-type pseudopotential\n");
       printf("unless the lmax=0 and n_rad_0=0,\n");
       printf("that is, unless Stephan has made a local pseudo\n");
       printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
       fflush(stdout);
       exit(1);
      }/*endif*/
     }/*endif*/
     ierr = 0;
     if( ((*loc_opt)==1)&&((*nrad_0)==0) ){ierr++;}
     if( ((*loc_opt)==2)&&((*nrad_1)==0) ){ierr++;}
     if( ((*loc_opt)==3)&&((*nrad_2)==0) ){ierr++;}
     if( ((*loc_opt)==4)&&((*nrad_3)==0) ){ierr++;} 
     if(ierr>0){
       printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
       printf("You have set up the last angular channel with \n");
       printf("no radial channels for a Goedecker-type pseudopotential\n");
       printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
       fflush(stdout);
       exit(1);
     }/*endif*/
     if( (*loc_opt)!=0 ){*n_ang = (*n_ang + 1);}
    }/*endif*/

/*--------------------------------------------------------------------------*/
/*  5) Gauss hermite pseudopotentials                                      */

    index = 11;
    sscanf(vps_dict[11].keyarg,"%lg",&tmp);
    ngh = (int) tmp;
    if(*ivps_label == 2 && ngh > 180) {
        keyarg_barf(vps_dict,filename,fun_key,index);
    }
   if(*ivps_label == 2 && ngh == 0){
       printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
       printf("You have set a gauss-hermite nonlocal atom\n");
       printf("but have zero integration points: ngh =  %d.\n",ngh);
       printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
       fflush(stdout);
       exit(1);
    }/*endif*/

    *pngh = ngh;
/*--------------------------------------------------------------------------*/
/* Set Radial channels if not Goedecker */

    if(*ivps_label != 5){
      *nrad_0=1;
      *nrad_1=1;
      *nrad_2=1;
      *nrad_3=1;
    }/*endif*/

/*-----------------------------------------------------------------------*/
 }/*end routine*/
/*==========================================================================*/







/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void make_vps_splin(char *vps_file,int loc_opt,int n_ang,
                    int ivps_label,int nsplin_g,double dg_spl,
                    double gmin_spl,double gmax_spl,double gmin_true,
                    double *vps0,double *vps1,double *vps2,double *vps3,
                    double *dvps0,double *dvps1,double *dvps2,double *dvps3,
                    double *gzvps,double *gzvps0,double *q_pseud,
                    double *vpsnorm, int nrad_0,int nrad_1,int nrad_2,
                    int nrad_3,
                    int cp_ptens_calc,int myid, MPI_Comm comm,int num_proc,
                    int cp_dual_grid_opt,double alpha_conv_dual,
                    int n_rad_max,int n_ang_max_gh,
                    int *pngh,double *rgh,double *wgh)

/*==========================================================================*/
/*               Begin subprogram:                                          */
      {/*begin routine*/
/*==========================================================================*/
/*               Local variable declarations                                */
#include "../typ_defs/typ_mask.h"

      double rmax;                   /* Num: Maximum spline radius       */
      double z_1,z_2;                /* Num: Charges for long range piece */
      double alpha_1,alpha_2;        /* Num: Ewald alpha's for long range */ 
      double dr;                     /* Num: delta r for spline           */
      double amat;                   /* Num: Useful temporary for spline  */
      double v_now;                  /* Num: current value of pp          */
      double rphi_now;               /* Num: current value of pseudo wf   */
      double zpol;                   /* Num: Polarization charge          */
      double gamma;                  /* Num: Useful temporary for long range */
      double dummy;                  /* Num: Really dum for pressure tensor*/
      int i;                         /* Num: Generic counter             */
      int ilong;                     /* Num: Flag for long range piece   */
      int n_ang_now;                 /* Num: Number of angular momentum  
                                               components                  */
      int nchan_tot;                 /* Num: Total number of radial channels*/
      int n_rad0_now,n_rad1_now;     /* Num: Number of radial channels for */
      int n_rad2_now,n_rad3_now;     /*   for angular momentum channel     */
      int iii;
      int iang;                      /* Num: Angular momentum counter    */
      int ishift_now;                /* Num: Angular momentum shift      */
      int iang_now;                  /* Num: Angular momentum counter    */
      int irad,jrad;                 /* Num: Radial channel counter      */
      int ioff_rad;                  /* Num: Radial channel shift        */
      int *nrad_l_tmp;               /* Temp array  # radial channels for
                                         given angular momentum channel */
      double *g;                     /* Lst: temporary spline array      */
      double *v_rphi,*v_loc,*r;      /* Lst: temporary spline arrays     */
      int ir;                        /* Num: counter for above arrays    */
      int nr;                        /* Num: Length of above arrays      */ 
      FILE *fp_vps_file;             /* File pointer to pseudo pot file  */
 
      int ierr;
      double tmp;
      int n_rad_max_sq = n_rad_max*n_rad_max;
      int ngh = *pngh;

 /* Temp gauss-hermite spline arrays */
      static int natm_typ_gh = 0;
      double *wgh_tmp;
      double *c0,*c1,*c2,*c3;
      double *dvl;
      double mach_eps = 1.0e-16;

/*==========================================================================*/
/*  I) Allocate g-space vector for spline                                   */

   g = (double *) cmalloc((nsplin_g)*sizeof(double))-1;

/*==========================================================================*/
/*  II) Open the electron-atom pseudopotential file for reading             */

   if(myid==0){
     fp_vps_file = cfopen(vps_file,"r");
   }/*endif*/

/*==========================================================================*/
/*  III) Set up g vectors                                                   */

   for(i=1;i <= nsplin_g;i++){ 
     g[i] = dg_spl*((double) (i-1)) + gmin_spl;
   }/*endfor*/

/*==========================================================================*/
/*  IV) Set up a local pseudo potential                                     */

   if(ivps_label <= 4) {

     if(myid==0){
       if(fscanf(fp_vps_file,"%d %lf %d\n",&nr,&rmax,&n_ang_now) != 3) 
                      {vps_read_error(vps_file);}
       if(fscanf(fp_vps_file,"%lf %lf %lf %lf\n",
                              &z_1,&alpha_1,&z_2,&alpha_2) != 4) 
                      {vps_read_error(vps_file);}
       if(fscanf(fp_vps_file,"%lf %lf\n",&zpol,&gamma) != 2)      
                 {vps_read_error(vps_file);}
     }/*endif*/

     if(num_proc>1){
      Bcast(&(nr),1,MPI_INT,0,comm);
      Bcast(&(rmax),1,MPI_DOUBLE,0,comm);
      Bcast(&(n_ang_now),1,MPI_INT,0,comm);
      Bcast(&(z_1),1,MPI_DOUBLE,0,comm);
      Bcast(&(alpha_1),1,MPI_DOUBLE,0,comm);
      Bcast(&(z_2),1,MPI_DOUBLE,0,comm);
      Bcast(&(alpha_2),1,MPI_DOUBLE,0,comm);
      Bcast(&(zpol),1,MPI_DOUBLE,0,comm);
      Bcast(&(gamma),1,MPI_DOUBLE,0,comm);
     }/*endif*/

   }/*endif ivps_label */
/*--------------------------------------*/
/* ivps_label == 5  GOEDECKER Potential */

   if(ivps_label == 5){
     if(myid==0){
       if(fscanf(fp_vps_file,"%d %lf %d\n",&nr,&rmax,&n_ang_now) != 3) 
                      {vps_read_error(vps_file);}

       if(fscanf(fp_vps_file,"%d %d %d %d\n",&n_rad0_now,
                                             &n_rad1_now,
                                             &n_rad2_now, 
                                             &n_rad3_now) != 4) 
                      {vps_read_error(vps_file);}
       if(fscanf(fp_vps_file,"%lf %lf %lf %lf\n",
                              &z_1,&alpha_1,&z_2,&alpha_2) != 4) 
                      {vps_read_error(vps_file);}
       if(fscanf(fp_vps_file,"%lf %lf\n",&zpol,&gamma) != 2) 
                      {vps_read_error(vps_file);}

      /*-------------------------------------*/
      /* Read in Radial Channel Matrices     */
      /*-------------------------------------*/

      /*-----------  S CHANNEL --------------*/
       ioff_rad = 0;
       for(irad=1; irad <= n_rad0_now; irad++){
        for(jrad=1; jrad <= n_rad0_now; jrad++){
         fscanf(fp_vps_file,"%lf ",&(vpsnorm[jrad+ioff_rad]));
        }/*endfor*/
        ioff_rad += n_rad_max;
       }/*endfor*/

      /*-----------  P CHANNEL --------------*/
       ioff_rad = n_rad_max*n_rad_max;
       for(irad=1; irad <= n_rad1_now; irad++){
        for(jrad=1; jrad <= n_rad1_now; jrad++){
         fscanf(fp_vps_file,"%lf ",&(vpsnorm[jrad+ioff_rad]));
        }/*endfor*/
        ioff_rad += n_rad_max;
       }/*endfor*/

      /*-----------  D CHANNEL --------------*/
       ioff_rad = 2.0*n_rad_max*n_rad_max;
       for(irad=1; irad <= n_rad2_now; irad++){
        for(jrad=1; jrad <= n_rad2_now; jrad++){
         fscanf(fp_vps_file,"%lf ",&(vpsnorm[jrad+ioff_rad]));
        }/*endfor*/
        ioff_rad += n_rad_max;
       }/*endfor*/

      /*-----------  F CHANNEL --------------*/
       ioff_rad = 3.0*n_rad_max*n_rad_max;
       for(irad=1; irad <= n_rad3_now; irad++){
        for(jrad=1; jrad <= n_rad3_now; jrad++){
         fscanf(fp_vps_file,"%lf ",&(vpsnorm[jrad+ioff_rad]));
        }/*endfor*/
        ioff_rad += n_rad_max;
       }/*endfor*/

     }/*endif myid*/

     if(num_proc>1){
       Bcast(&(nr),1,MPI_INT,0,comm);
       Bcast(&(rmax),1,MPI_DOUBLE,0,comm);
       Bcast(&(n_rad0_now),1,MPI_INT,0,comm);
       Bcast(&(n_rad1_now),1,MPI_INT,0,comm);
       Bcast(&(n_rad2_now),1,MPI_INT,0,comm);
       Bcast(&(n_rad3_now),1,MPI_INT,0,comm);
       Bcast(&(n_ang_now),1,MPI_INT,0,comm);
       Bcast(&(z_1),1,MPI_DOUBLE,0,comm);
       Bcast(&(alpha_1),1,MPI_DOUBLE,0,comm);
       Bcast(&(z_2),1,MPI_DOUBLE,0,comm);
       Bcast(&(alpha_2),1,MPI_DOUBLE,0,comm);
       Bcast(&(zpol),1,MPI_DOUBLE,0,comm);
       Bcast(&(gamma),1,MPI_DOUBLE,0,comm);
       Bcast(&(vpsnorm[1]),((n_ang_now+1)*n_rad_max*n_rad_max),
                  MPI_DOUBLE,0,comm);
     }/*endif*/

   }/*endif ivps_label goedecker*/

   v_rphi = (double *) cmalloc(nr*sizeof(double))-1;
   v_loc  = (double *) cmalloc(nr*sizeof(double))-1;
   r      = (double *) cmalloc(nr*sizeof(double))-1;

   *q_pseud = z_1 + z_2;

/*--------------------------------------------------------------------------*/
/*   A) Allocate real space arrays                                          */

   dr = rmax/((double ) nr);
   for(i=1;i <= nr;i++){
     r[i] = ((double ) (i-1))*dr;
   }/*endfor*/

   if( (n_ang_now < n_ang) && (myid==0) && (ivps_label!=5) ) {
       printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
       printf("Dude, like given the number of angular momentum\n"); 
       printf("components you've specified in the\n");
       printf("pseudopotential file %s\n",vps_file);
       printf("proceeding with the simulation would be a\n");
       printf("pointless exercise leading to most bogus results\n");
       printf("%d vs %d\n",n_ang_now,n_ang);
       putchar('\n');
       printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
       fflush(stdout);
       exit(1);
   }/*endif*/

   if( (n_ang_now != n_ang-1) && (ivps_label==5) && (myid==0)) {
     if( (n_ang_now==0)&&(loc_opt==0) && (n_rad0_now==0)){
      iii = 1;
     }else{
       printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
       printf("Dude, the number of angular momentum\n"); 
       printf("components you've specified in the\n");
       printf("pseudopotential file %s\n",vps_file);
       printf("does not match does the value in the controller,\n");
       printf("%d vs %d\n,",n_ang_now,n_ang-1);
       printf("for a Goedecker-type pseudopotential. Stephan would\n");
       printf("not approve and neither do we!\n");
       putchar('\n');
       printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
       fflush(stdout);
       exit(1);
     }/*endif*/
   }/*endif*/

   if( (ivps_label==5) && (myid==0)) {   
      ierr = 0;
      if( (nrad_0 != n_rad0_now) && (n_ang-1>=0)){ierr++;}
      if( (nrad_1 != n_rad1_now) && (n_ang-1>=1)){ierr++;}
      if( (nrad_2 != n_rad2_now) && (n_ang-1>=2)){ierr++;}
      if( (nrad_3 != n_rad3_now) && (n_ang-1>=3)){ierr++;}
      if(ierr>0){
       printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
       printf("Dude, the number of radial channel \n"); 
       printf("components you've specified in the\n");
       printf("pseudopotential file %s\n",vps_file);
       printf("does not match does the value in the controller,\n");
       printf("for a Goedecker-type pseudopotential. Stephan would\n");
       printf("not approve and neither do we!\n");
       if( (nrad_0 != n_rad0_now) && (n_ang-1>0)){
         printf("l=0 : %d vs %d\n",nrad_0,n_rad0_now);
       }/*endif*/
       if( (nrad_1 != n_rad1_now) && (n_ang-1>1)){
         printf("l=1 : %d vs %d\n",nrad_1,n_rad1_now);
       }/*endif*/
       if( (nrad_2 != n_rad2_now) && (n_ang-1>2)){
         printf("l=2 : %d vs %d\n",nrad_2,n_rad2_now);
       }/*endif*/
       if( (nrad_3 != n_rad3_now) && (n_ang-1>3)){
         printf("l=3 : %d vs %d\n",nrad_3,n_rad3_now);
       }/*endif*/
       putchar('\n');
       printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
       fflush(stdout);
       exit(1);
     }/*endif*/
   }/*endif*/


/*--------------------------------------------------------------------------*/
/*   B) Find v_loc                                                         */

   if(ivps_label < 5){ /* NOT GOEDECKER TYPE */

     if(myid==0){

       for(iang=1;iang <= (loc_opt + 1);iang++) {
         for(ir=1;ir <= nr;ir++) {          
           if(fscanf(fp_vps_file,"%lf %lf\n",&v_now,&rphi_now) != 2) 
                           {vps_read_error(vps_file);}
           v_loc[ir] = v_now;
         }/* endfor */ 
       } /* endfor */
     }/*endif*/

   }else{ /*GOEDECKER TYPE */

    if(myid==0){
      nchan_tot = n_rad0_now+n_rad1_now+n_rad2_now+n_rad3_now;
      for(iang=1; iang<=(nchan_tot+1); iang++){
        for(ir=1;ir <= nr;ir++) {          
           if(fscanf(fp_vps_file,"%lf %lf\n",&v_now,&rphi_now) != 2) 
                           {vps_read_error(vps_file);}
           v_loc[ir] = v_now;
        }/* endfor */ 
       }/*endfor*/

    }/*endif myid*/

   }/*endif ivps_label*/

   if(num_proc>1){Bcast(&(v_loc[1]),nr,MPI_DOUBLE,0,comm);}

   for(ir=1;ir <= nr;ir++) {v_rphi[ir] = v_loc[ir]*r[ir];}

     ilong      = 1;
     iang_now   = 0;
     ishift_now = (loc_opt)*n_rad_max*nsplin_g;
     if(cp_ptens_calc == 1){
        slow_bess_vps(v_rphi,nr,dr,r,&(vps0)[ishift_now],
                     &(dvps0)[ishift_now],
                      nsplin_g,g,gmin_true,
                      z_1,alpha_1,z_2,alpha_2,zpol,gamma,ilong,iang_now,
                      gzvps,&gzvps0[1],cp_ptens_calc,cp_dual_grid_opt,
                      alpha_conv_dual);
     } else {
        slow_bess_vps(v_rphi,nr,dr,r,&(vps0)[ishift_now],
                      &dummy,
                      nsplin_g,g,gmin_true,
                      z_1,alpha_1,z_2,alpha_2,zpol,gamma,ilong,iang_now,
                      gzvps,&gzvps0[1],cp_ptens_calc,cp_dual_grid_opt,
                      alpha_conv_dual);
     }/* endif */
     spline_fit(&(vps0)[ishift_now],&(vps1)[ishift_now],
                &(vps2)[ishift_now],&(vps3)[ishift_now],g,nsplin_g);
     if(cp_ptens_calc == 1){
       spline_fit(&(dvps0)[ishift_now],&(dvps1)[ishift_now],
                  &(dvps2)[ishift_now],&(dvps3)[ishift_now],g,nsplin_g);
     }/* endif */

/*==========================================================================*/
/*  V) Set up a null local pseudo potential                                 */

     if(ivps_label == 4) {

       loc_opt  = 0;
       n_ang    = 0;
       *q_pseud = 0.0;
       for(i=1;i <= nsplin_g; i++) {
         vps0[i] = 0.0;
         vps1[i] = 0.0;
         vps2[i] = 0.0;
         vps3[i] = 0.0;
       }/*endfor*/
       if(cp_ptens_calc == 1){
         for(i=1;i <= nsplin_g; i++) {
           dvps0[i] = 0.0;
           dvps1[i] = 0.0;
           dvps2[i] = 0.0;
           dvps3[i] = 0.0;
         } /* endfor */
       }/* endif */

     }/* endif:null potential */

/*==========================================================================*/
/* VI) KB nonlocal potential                                                */

     if(ivps_label == 1) {
/*--------------------------------------------------------------------------*/
/*    A) rewind                                                      */

       if(myid==0){
         rewind(fp_vps_file);      
         if(fscanf(fp_vps_file,"%d %lf %d\n",&nr,&rmax,&n_ang_now) != 3) 
                     {vps_read_error(vps_file);}
         if(fscanf(fp_vps_file,"%lf %lf %lf %lf\n",
                               &z_1,&alpha_1,&z_2,&alpha_2)!= 4) 
                    {vps_read_error(vps_file);}
         if(fscanf(fp_vps_file,"%lf %lf\n",&zpol,&gamma) != 2) 
                    {vps_read_error(vps_file);}
       }/*endif*/
/*--------------------------------------------------------------------------*/
/*    B) Spline projection operator                                         */
 
       vpsnorm[(loc_opt*n_rad_max_sq + 1)] = 0.0;
       for(iang = 1; iang <= n_ang + 1; iang++) {
         amat = 0.0;
         if(myid==0){
          for(ir=1; ir <= nr; ir++) {
           if(fscanf(fp_vps_file,"%lf %lf\n",&v_now,&rphi_now) != 2) 
                       {vps_read_error(vps_file); }
            v_rphi[ir] = (v_now-v_loc[ir])*rphi_now;
            amat += rphi_now*v_rphi[ir]*dr;
          } /* endfor */
         }/*endif*/
         if(num_proc>1){
          Bcast(&(v_rphi[1]),nr,MPI_DOUBLE,0,comm);
          Bcast(&(amat),1,MPI_DOUBLE,0,comm);
         }/*endif*/
         if(iang != loc_opt+1) {
           if(fabs(amat) > mach_eps)
              vpsnorm[((iang-1)*n_rad_max_sq+1)] = (1.0/amat); 
           else 
	     vpsnorm[((iang-1)*n_rad_max_sq+1)] = 1.0;
           ishift_now = (iang-1)*n_rad_max*nsplin_g;
           ilong = 0;
           iang_now = iang-1;
           if(cp_ptens_calc == 1){
            slow_bess_vps(v_rphi,nr,dr,r,&(vps0[ishift_now]),
                          &(dvps0[ishift_now]),
                          nsplin_g,g,gmin_true,
                          z_1,alpha_1,z_2,alpha_2,zpol,gamma,ilong,iang_now,
                          gzvps,&gzvps0[1],cp_ptens_calc,cp_dual_grid_opt,
                          alpha_conv_dual);
           } else {
            slow_bess_vps(v_rphi,nr,dr,r,&(vps0[ishift_now]),
                          &dummy,
                          nsplin_g,g,gmin_true,
                          z_1,alpha_1,z_2,alpha_2,zpol,gamma,ilong,iang_now,
                          gzvps,&gzvps0[1],cp_ptens_calc,cp_dual_grid_opt,
                          alpha_conv_dual);
           }/* endif */
           spline_fit(&(vps0[ishift_now]),&(vps1[ishift_now]),
                      &(vps2[ishift_now]),&(vps3[ishift_now]),g,nsplin_g);

          if(cp_ptens_calc == 1){
           spline_fit(&(dvps0[ishift_now]),&(dvps1[ishift_now]),
                      &(dvps2[ishift_now]),&(dvps3[ishift_now]),g,nsplin_g);
          }/* endif */

         } /* endif:loc_opt */
       } /* endfor:channels */
     } /* endif:KB nonlocal */

/*==========================================================================*/
/* VII) Gauss-Hermite nonlocal pseudopotentials                             */

     if(ivps_label == 2){
       natm_typ_gh++;
       wgh_tmp = (double *) cmalloc(ngh*sizeof(double))-1;

       weight_node_gauss_hermite(ngh,rgh,wgh_tmp);

 /*Only use Gauss-Hermite integration points with a weight > 1.0e-10 */   
 /* CHECK WHETHER THIS LIMIT STATEMENT IS NEEDED */

   limit_gauss_hermite(&(ngh),rgh,wgh_tmp);  

   *pngh = ngh;


/*--------------------------------------------------------------------------*/
/* Read in the pseudo potential                                             */

/*    A) rewind                                                             */

       if(myid==0){
         rewind(fp_vps_file);      
         if(fscanf(fp_vps_file,"%d %lf %d\n",&nr,&rmax,&n_ang_now) != 3) 
                     {vps_read_error(vps_file);}
         if(fscanf(fp_vps_file,"%lf %lf %lf %lf\n",
                               &z_1,&alpha_1,&z_2,&alpha_2)!= 4) 
                    {vps_read_error(vps_file);}
         if(fscanf(fp_vps_file,"%lf %lf\n",&zpol,&gamma) != 2) 
                    {vps_read_error(vps_file);}
       }/*endif*/
/*--------------------------------------------------------------------------*/
/* This only works when vloc is the highest angular momentum channel        */

       for(iang = 1; iang < (loc_opt + 1); iang++) {
         if(myid==0){
          for(ir=1; ir <= nr; ir++) {
           if(fscanf(fp_vps_file,"%lf %lf\n",&v_now,&rphi_now) != 2) 
                       {vps_read_error(vps_file); }
            v_rphi[ir] = (v_now-v_loc[ir]); 
          } /* endfor */
         }/*endif myid*/

         if(num_proc>1){
          Bcast(&(v_rphi[1]),nr,MPI_DOUBLE,0,comm);
         }/*endif*/

/*--------------------------------------------------------------------------*/
/* i) make r                                                                */
   printf("gh 6\n");
       c0  = (double *) cmalloc(nr*sizeof(double))-1;  
       c1  = (double *) cmalloc(nr*sizeof(double))-1;  
       c2  = (double *) cmalloc(nr*sizeof(double))-1;  
       c3  = (double *) cmalloc(nr*sizeof(double))-1;  
       r   = (double *) cmalloc(nr*sizeof(double))-1;  
       dvl = (double *) cmalloc(ngh*sizeof(double))-1;

      dr = rmax/((double) nr);

      for(i=1; i<= nr; i++){
        r[i]  = (double)(i-1)*dr;
        c0[i] = v_rphi[i];
      }/*endfor*/

/*--------------------------------------------------------------------------*/
/* ii) Fit to spline                                                         */

       spline_fit(c0,c1,c2,c3,r,nr);

/*--------------------------------------------------------------------------*/
/* ii) Fetch dvl                                                           */

       get_dvl(ngh,rgh,dvl,c0,c1,c2,c3,dr,nr);

/*--------------------------------------------------------------------------*/
/* iii) Make the general weight                                             */

     make_weight_gen(wgh_tmp,wgh,rgh,dvl,ngh,
                     iang,natm_typ_gh,n_ang_max_gh,myid);
     }/* endfor:channels */

     cfree(&(wgh_tmp[1]));
     cfree(&(c0[1]));
     cfree(&(c1[1]));
     cfree(&(c2[1]));
     cfree(&(c3[1]));
     cfree(&(r[1]));
     cfree(&(dvl[1]));

    }/*endif GAUSS-HERMITE*/

/*==========================================================================*/
/* VIIII) GOEDECKER nonlocal potential                                         */

    if(ivps_label == 5) {
/*--------------------------------------------------------------------------*/
/*    A) rewind                                                      */
      if(myid==0){
       rewind(fp_vps_file);      
       if(fscanf(fp_vps_file,"%d %lf %d\n",&nr,&rmax,&n_ang_now) != 3) 
                      {vps_read_error(vps_file);}
       if(fscanf(fp_vps_file,"%d %d %d %d\n",&n_rad0_now,
                                             &n_rad1_now,
                                             &n_rad2_now,
                                             &n_rad3_now) != 4) 
                      {vps_read_error(vps_file);}
       if(fscanf(fp_vps_file,"%lf %lf %lf %lf\n",
                              &z_1,&alpha_1,&z_2,&alpha_2) != 4) 
                      {vps_read_error(vps_file);}
       if(fscanf(fp_vps_file,"%lf %lf\n",&zpol,&gamma) != 2) 
                      {vps_read_error(vps_file);}
      /*-----------  Skip the matrix  --------------*/
       for(irad=1; irad <= (n_rad0_now*n_rad0_now); irad++){
         fscanf(fp_vps_file,"%lf ",&tmp);
       }/*endfor*/
       for(irad=1; irad <= (n_rad1_now*n_rad1_now); irad++){
         fscanf(fp_vps_file,"%lf ",&tmp);
       }/*endfor*/
       for(irad=1; irad <= (n_rad2_now*n_rad2_now); irad++){
         fscanf(fp_vps_file,"%lf ",&tmp);
       }/*endfor*/
       for(irad=1; irad <= (n_rad3_now*n_rad3_now); irad++){
         fscanf(fp_vps_file,"%lf ",&tmp);
       }/*endfor*/
      }/*endif myid*/

/*--------------------------------------------------------------------------*/
/*    B) Spline projection operator                                         */

      nrad_l_tmp     = (int *) cmalloc(4*sizeof(int));
      nrad_l_tmp[0] = n_rad0_now;
      nrad_l_tmp[1] = n_rad1_now;
      nrad_l_tmp[2] = n_rad2_now;
      nrad_l_tmp[3] = n_rad3_now;

      /* loc_opt=n_ang so loop stops at n_ang*/
      for(iang = 1; iang <= n_ang ; iang++) {
       for(irad = 1; irad <= nrad_l_tmp[(iang-1)]; irad++){ 
        if(myid==0){
          for(ir=1; ir <= nr; ir++) {
           if(fscanf(fp_vps_file,"%lf %lf\n",&rphi_now,&tmp) != 2) 
                       {vps_read_error(vps_file); }
            v_rphi[ir] = rphi_now;
          } /* endfor */
        }/*endif*/
        if(num_proc>1){
          Bcast(&(v_rphi[1]),nr,MPI_DOUBLE,0,comm);
        }/*endif*/
        ishift_now = (iang-1)*n_rad_max*nsplin_g + (irad-1)*nsplin_g;
        ilong = 0;
        iang_now = iang-1;
        if(cp_ptens_calc == 1){
          slow_bess_vps(v_rphi,nr,dr,r,&(vps0[ishift_now]),
                        &(dvps0[ishift_now]),
                        nsplin_g,g,gmin_true,
                        z_1,alpha_1,z_2,alpha_2,zpol,gamma,ilong,iang_now,
                        gzvps,&gzvps0[irad],cp_ptens_calc,cp_dual_grid_opt,
                        alpha_conv_dual);
        } else {
          slow_bess_vps(v_rphi,nr,dr,r,&(vps0[ishift_now]),
                        &dummy,
                        nsplin_g,g,gmin_true,
                        z_1,alpha_1,z_2,alpha_2,zpol,gamma,ilong,iang_now,
                        gzvps,&gzvps0[irad],cp_ptens_calc,cp_dual_grid_opt,
                        alpha_conv_dual);
        }/* endif */
        spline_fit(&(vps0[ishift_now]),&(vps1[ishift_now]),
                   &(vps2[ishift_now]),&(vps3[ishift_now]),g,nsplin_g);

        if(cp_ptens_calc == 1){
          spline_fit(&(dvps0[ishift_now]),&(dvps1[ishift_now]),
                     &(dvps2[ishift_now]),&(dvps3[ishift_now]),g,nsplin_g);
        }/* endif */

       }/*endfor: radial channels */
      }/* endfor:angular channels */

      cfree(nrad_l_tmp);

    } /* endif:GOEDECKER nonlocal */

/*==========================================================================*/
/*   VII) Free memory                                                       */

    if(ivps_label != 4) {
      cfree(&(v_rphi[1]));
      cfree(&(v_loc[1]));
      if(ivps_label != 2) cfree(&(r[1]));
    }/*endif*/
    cfree(&(g[1]));


/*==========================================================================*/
/*   VIII) Close the file and done                                          */

   if(myid==0){
     fclose(fp_vps_file);
   }/*endif*/

/*-----------------------------------------------------------------------*/
  }/*end routine*/
/*==========================================================================*/







/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

 void slow_bess_vps(double v_rphi[],int nr,double dr,double r[],
                double fv_rphi[],double fdv_rphi[],int nsplin_g,double g[],
                double gmin_true,
                double z_1,double alpha_1,double z_2,double alpha_2,
                double zpol,double gamma,int ilong,int iang,
                double *gzvps,double *gzvps0,int cp_ptens_calc,
                int cp_dual_grid_opt, double alpha_conv_dual)

/*==========================================================================*/
/*               Begin subprogram:                                          */
      {/*begin routine*/
/*==========================================================================*/
/*               Local variable declarations                                */

   double c[8],angamma[8],angamma2[8];  /* Lst: Useful temporary arrays */
   double phi;                    /* Num: Angle used in 
                                          polarization corrections*/
   double g2,rj0,rj1,rj2,rj3,gzero,arg,fpi,pi,tpi,fpidr;
                                  /* Num: Useful constants and
                                                    temporaries             */
   double r2dj0,r2dj1,r2dj2,r2dj3,j4,j3,j2,j1,j0;
   int ig,ir,i;              /* Num: Counters                */
   int iii;
   double falpha_12,falpha_22;

  double dg = 0.00001;
  double fvp,fvm,gt,df;
  double gm1,gamm_plu,gamm_min,pre_c;
  double ztot,falpha_conv_dual;
 

/*==========================================================================*/
/* II) Get some constants                                                    */

   pi  = M_PI;
   tpi = 2.0*pi;
   fpi = 4.0*pi;
   fpidr = fpi*dr;
   for(ig=1;ig <= nsplin_g;ig++) {
     fv_rphi[ig] = 0.0;
   }
   if(cp_ptens_calc ==1 ){
    for(ig=1;ig <= nsplin_g;ig++) {
      fdv_rphi[ig] = 0.0;
    }/* endfor */
   }/* endif */

/*==========================================================================*/
/*==========================================================================*/

   switch(iang) {

/*==========================================================================*/
/* I) L=0 Term */

    case 0:

/*--------------------------------------------------------------------------*/
/*    i) g = 0 : fdv = 0                                                    */

    gzero = 0.0;
    for(ir=2;ir <= nr;ir++) {gzero += fpidr*r[ir]*v_rphi[ir];}
    if(ilong == 1) {*gzvps  = gzero;}
    if(ilong != 1) {*gzvps0 = gzero;}

/*--------------------------------------------------------------------------*/
/*    ii) g ne 0                                                            */

     for(ig=1;ig <= nsplin_g; ig++) {
       for(ir=2;ir <= nr; ir++) {
        arg = r[ir]*g[ig];
        rj0 = (sin(arg)/arg)*r[ir];
        fv_rphi[ig]  += fpidr*rj0*v_rphi[ir];
       } /* endfor */
      if(cp_ptens_calc ==1 ){
        for(ir=2;ir <= nr; ir++) {
         arg = r[ir]*g[ig];
         j1 = (sin(arg)/arg - cos(arg))/arg;
         r2dj0 = -(j1)*r[ir]*r[ir];
         fdv_rphi[ig] += fpidr*r2dj0*v_rphi[ir]/g[ig];
        } /* endfor */  
       }/* endif */
      } /* endfor */
    

    break;

/*==========================================================================*/
/* II) L=1 Term */

    case 1:

/*--------------------------------------------------------------------------*/
/*    i) g ne 0                                                             */

      for(ig=1;ig <= nsplin_g; ig++) {
        for(ir=2;ir <= nr; ir++) {
          arg = r[ir]*g[ig];
          rj1 = ((sin(arg)/arg - cos(arg))/arg)*r[ir];
          fv_rphi[ig] += fpidr*rj1*v_rphi[ir];
        } /* endfor */
       if(cp_ptens_calc ==1 ){
        for(ir=2;ir <= nr; ir++) {
          arg = r[ir]*g[ig];
          j0 = 1.0*(sin(arg)/arg);
          j2 = ((3.0/(arg*arg)-1.0)*sin(arg)-3.0*cos(arg)/arg)/arg;
          r2dj1 = (j0-2.0*j2 )*r[ir]*r[ir]/3.0;
          fdv_rphi[ig] += fpidr*r2dj1*v_rphi[ir]/g[ig];
        } /* endfor */
       }/* endif */
      } /* endfor */

    break;

/*==========================================================================*/
/* III) L=2 Term */

    case 2:

/*--------------------------------------------------------------------------*/
/*    i) g ne 0                                                             */

      for(ig=1;ig <= nsplin_g; ig++) {
        for(ir=2;ir <= nr; ir++) {
          arg = r[ir]*g[ig];
          rj2 = (((3.0/(arg*arg)-1.0)*sin(arg)-3.0*cos(arg)/arg)/arg)*r[ir];
          fv_rphi[ig] += fpidr*rj2*v_rphi[ir];
        } /* endfor */
       if(cp_ptens_calc ==1 ){
        for(ir=2;ir <= nr; ir++) {
          arg = r[ir]*g[ig];
          j1  = (sin(arg)/arg - cos(arg))/arg;
          j3  = ((15.0/(arg*arg) - 6.0)*sin(arg)/arg + 
                          (1.0 - 15.0/(arg*arg))*cos(arg))/arg;
          r2dj2 = ( 2.0*j1 - 3.0*j3 )*r[ir]*r[ir]/5.0;
          fdv_rphi[ig] += fpidr*r2dj2*v_rphi[ir]/g[ig];
        } /* endfor */
       }/* endif */
      } /* endfor */

    break;

/*==========================================================================*/
/* IV) L=3 Term */

    case 3:

/*--------------------------------------------------------------------------*/
/*    i) g ne 0                                                             */

      for(ig=1;ig <= nsplin_g; ig++) {
        for(ir=2;ir <= nr; ir++) {
          arg = r[ir]*g[ig];
          rj3 = (((15.0/(arg*arg) - 6.0)*sin(arg)/arg + 
                 (1.0 - 15.0/(arg*arg))*cos(arg))/arg)*r[ir];
          fv_rphi[ig] += fpidr*rj3*v_rphi[ir];
        } /* endfor */
       if(cp_ptens_calc ==1 ){
        for(ir=2;ir <= nr; ir++) {
          arg = r[ir]*g[ig];
          j3 = ((15.0/(arg*arg) - 6.0)*sin(arg)/arg + 
                 (1.0 - 15.0/(arg*arg))*cos(arg))/arg;
          j2 = ((3.0/(arg*arg)-1.0)*sin(arg)-3.0*cos(arg)/arg)/arg;
          j4 = (7.0/arg)*j3 - j2;
          r2dj3 = ( 3.0*j2-4.0*j4 )*r[ir]*r[ir]/7.0;
          fdv_rphi[ig] += fpidr*r2dj3*v_rphi[ir]/g[ig];
        } /* endfor */
       }/* endif */
      } /* endfor */

    break;

   } /* end switch */

/*==========================================================================*/
/*==========================================================================*/


/*==========================================================================*/
/* V) Add in long range part if necessary                                  */

   if(ilong == 1) {
/*--------------------------------------------------------------------------*/
/*   A) Coulomb                                                             */
     if((z_1 != 0.0) || (z_2 != 0.0)) {
       ztot = 0.0;
       if(cp_dual_grid_opt==2){ztot = z_1 + z_2;}
       for(ig=1;ig <= nsplin_g;ig++) {
         g2 = g[ig]*g[ig];
         fv_rphi[ig] -= 
               (z_1*fpi*exp((-0.25*g2/(alpha_1*alpha_1)))/g2
             +  z_2*fpi*exp((-0.25*g2/(alpha_2*alpha_2)))/g2
             -  ztot*fpi*(
                        exp((-0.25*g2/(alpha_conv_dual*alpha_conv_dual))))/g2
               );
    } /* endfor */
       falpha_12        = 4.0*alpha_1*alpha_1;
       falpha_22        = 4.0*alpha_2*alpha_2; 
       falpha_conv_dual = 4.0*alpha_conv_dual*alpha_conv_dual; 
       (*gzvps) += fpi*(z_1/falpha_12 + z_2/falpha_22
                       -ztot/falpha_conv_dual);
       if(cp_ptens_calc ==1 ){
        for(ig=1;ig <= nsplin_g;ig++) {
         g2 = g[ig]*g[ig];
         fdv_rphi[ig] += 
               ( (2.0*z_1*fpi*exp((-g2/falpha_12))/g2)
                *(1.0/(falpha_12)+1.0/g2)
              +  (2.0*z_2*fpi*exp((-g2/falpha_22))/g2)
                *(1.0/(falpha_22)+1.0/g2)
               );
        } /* endfor */
       }/* endif */
     } /* endif */
/*--------------------------------------------------------------------------*/
/*   B) Polarization corrections: */
     if(zpol != 0) {
        for(i=1;i<=7;i++) {
          angamma[i]  = ((double) (i-1))*gamma;
          angamma2[i] = angamma[i]*angamma[i];
        } /* endfor */
        c[1] =  1.0;  c[2] = -6.0;  c[3] = 15.0;
        c[4] =-20.0;  c[5] = 15.0;  c[6] = -6.0;
        c[7] =  1.0;
        for(i=2;i<=7;i++) {
          (*gzvps)  -= 0.50*zpol*fpi*c[i]*angamma[i]*log(angamma[i]);
        } /* endfor */
        if(cp_ptens_calc ==1 ){
         for(i=1;i<=7;i++) {
          for(ig=1;ig<=nsplin_g;ig++) {
            g2            = g[ig]*g[ig];  
            gm1           = 1.0/g[ig];
            gamm_plu      = angamma2[i]+g2;
            gamm_min      = angamma2[i]-g2;
            phi           = atan2(g[ig],angamma[i]);
            pre_c         = 0.50*zpol*tpi*c[i];
            fv_rphi[ig]  -= pre_c*(gm1*gamm_min*phi+angamma[i]*log(gamm_plu));
            fdv_rphi[ig] -= pre_c*(gm1*gamm_min*(-phi*gm1+angamma[i]/gamm_plu)
                                  -2.0*phi+2.0*g[ig]*angamma[i]/gamm_plu)*gm1;
          } /* endfor ig */
         } /* endfor i */
        }else{
         for(i=1;i<=7;i++) {
          for(ig=1;ig<=nsplin_g;ig++) {
            g2            = g[ig]*g[ig];  
            gm1           = 1.0/g[ig];
            gamm_plu      = angamma2[i]+g2;
            gamm_min      = angamma2[i]-g2;
            phi           = atan2(g[ig],angamma[i]);
            pre_c         = 0.50*zpol*tpi*c[i];
            fv_rphi[ig]  -= pre_c*(gm1*gamm_min*phi+angamma[i]*log(gamm_plu));
          } /* endfor ig */
         } /* endfor i */
        }/*endif*/
    } /* endif zpol */

   } /* endif ilong */

/*-----------------------------------------------------------------------*/
   }/*end routine*/
/*==========================================================================*/




/*==========================================================================*/

void vps_read_error(char *vps_file)

/*==========================================================================*/
{/*Begin routine*/
/*-----------------------------------------------------------------------*/

   printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
   printf("Error reading input from electron-atom\n"); 
   printf("pseudopotential file %s\n",vps_file);
   printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
   fflush(stdout);
   exit(1);

/*-----------------------------------------------------------------------*/
}/*end routine*/
/*==========================================================================*/



/*==========================================================================*/
/*==========================================================================*/

void get_dvl(int ngh,double *rgh,double *dvl,double *c0,
              double *c1, double *c2, double *c3,double dr,int nr)

/*==========================================================================*/
{/*Begin routine*/
  
  int igh,iii;
  double r,h,h0;  
/*-----------------------------------------------------------------------*/

  for(igh=1;igh<=ngh;igh++){
    r = rgh[igh]*JUERG_FACTOR;
    iii = r/dr + 1;
    iii = MIN(iii,nr);
    iii = MAX(iii,1);
    h0  = (double)(iii-1)*dr;
    h = r-h0;
    dvl[igh] = ((c3[iii]*h+c2[iii])*h+c1[iii])*h+c0[iii];
  }/* endfor */

/*-----------------------------------------------------------------------*/
}/*end routine*/
/*==========================================================================*/



/*==========================================================================*/
/*==========================================================================*/

void make_weight_gen(double *wgh_tmp,double *wgh,double *rgh, double *dvl,
                     int ngh,int iang,int natm_typ_gh,
                     int n_ang_max_gh,int myid)

/*==========================================================================*/
{/*Begin routine*/
  
  int igh,iii,ioff;
  double juerg_wght,arg;
  double r;
  double fpisq;

/*-----------------------------------------------------------------------*/

   fpisq = (4.0*M_PI)*(4.0*M_PI);

   ioff = (natm_typ_gh-1)*n_ang_max_gh*ngh + (iang-1)*ngh;

  for(igh=1;igh<=ngh;igh++){
    r   = rgh[igh]*JUERG_FACTOR;
    arg = rgh[igh]*rgh[igh];
    juerg_wght = JUERG_FACTOR*exp(arg);

/* FIX UP FOR DUAL GRID OPTION */

#ifdef  JUERG_FACTOR_ON
   wgh[igh+ioff] = (wgh_tmp[igh]*fpisq*dvl[igh]*r*r*juerg_wght);  

#define CFPI_SQ_OFF
#ifdef  CFPI_SQ
   if(iang == 1){
     wgh[igh+ioff] = fpisq;
   }else{
     wgh[igh+ioff] = fpisq;
   }
#endif

#else
    wgh[igh+ioff] = (wgh_tmp[igh]*fpisq*dvl[igh]*r*r);
#endif


  }/* endfor */

/*-----------------------------------------------------------------------*/
}/*end routine*/
/*==========================================================================*/


