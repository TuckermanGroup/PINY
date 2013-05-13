/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
/*                                                                          */
/*                         PI_MD:                                           */
/*             The future of simulation technology                          */
/*             ------------------------------------                         */
/*                     Module: set_ewald                                    */
/*                                                                          */
/* This reads in and sets up the k-space                                    */
/*                                                                          */
/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

#include "standard_include.h"
#include "../typ_defs/typedefs_class.h"
#include "../typ_defs/typedefs_gen.h"
#include "../typ_defs/typedefs_cp.h"
#include "../typ_defs/typedefs_par.h"
#include "../typ_defs/typedefs_bnd.h"
#include "../proto_defs/proto_cp_ewald_entry.h"
#include "../proto_defs/proto_cp_ewald_local.h"
#include "../proto_defs/proto_energy_cp_local.h"
#include "../proto_defs/proto_math.h"
#include "../proto_defs/proto_friend_lib_entry.h"
#include "../proto_defs/proto_communicate_wrappers.h"

#define DEBUG_CLUS_CORR_OFF
#define CHECK_CLUS_OFF

#define MAX_INT 12.0
#define MIN3(A,B,C) (MIN(MIN((A),(B)),(C)))

#define ORIG_OFF
#define PME

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
/* Control ewald/cp g-space initialization */
/*==========================================================================*/

void control_set_cp_ewald(SIMOPTS *simopts,CELL *cell,
                          CPCOEFFS_INFO *cpcoeffs_info,
                          EWALD *ewald, CPEWALD *cpewald,CP_PARSE *cp_parse,
                          double *gmin_true,double *gmin_spl,double *gmax_spl,
                          EWD_SCR *ewd_scr,int kmax_ewd,int kmax_res,
                          double *tot_memory,int int_res_ter,
                          PART_MESH *part_mesh,ECOR *ecor, int myid,
                          int cp_lsda,int cp_min_diis,int cp_dual_grid_opt_on) 

/*=======================================================================*/
/*            Begin subprogram:                                          */
   {/*begin routine*/
/*=======================================================================*/
/*            Local variable declarations:                               */

   double ecut_now;                    /* Num: Energy cutoff                 */
   double ecut_dens_cp_box_now;        /* Num: Energy cutoff                 */
   double ecut_res,ecut_tmp;
   double ecut_lg;                     /* Num: Energy cutoff for dens        */
   double ecut_sm;                     /* Num: Energy cutoff for wf          */
   double *hmati_ewd;                  /* Num: Inverse ewald h matrix        */
   double *hmati_ewd_cp;               /* Num: Inverse cp ewald h matrix     */
   double deth,deth_cp,side;           /* Num: Volumes and sizes             */
   double now_mem;                     /* Num: Memory used here              */
   int idum1=1,idum2=1,idum3=1;
   int iii;                            /* Num: Debug tool                    */
   int nmall;
   int cp_on;                          /* Num: CP on flag                    */
   int nk1,nk2,nk3;                    /* Num: Number of k vec on small grid */
   int nkf1_dens_cp_box,
       nkf2_dens_cp_box,
       nkf3_dens_cp_box;               /* Num: Number of k vec on large  grid*/
   int nkf1,nkf2,nkf3;                 /* Num: Number of k vec on mega  grid */

   int nktot,nktot_res,nktot_sm;       /* Num: Spherically cutoff k vecs     */
   int nktot_dens_cp_box;
   int ncoef,ncoef_dens_cp_box,ncoef_l;

   int ngrid_tot,nlen_pme;             /* Num: PME sizes                     */
   int ngrid_a_res,ngrid_b_res,ngrid_c_res;
   int ngrid_a, ngrid_b, ngrid_c;

   int *kmaxv;                         /* Lst: K-vector ranges               */
   int *kmax_cp_tmp,*kmaxv_res,cp_on_tmp;
   int *kmax_cp;
   int *kmaxv_dens_cp_box;
   int *kmax_cp_dens_cp_box;

   int pme_b_opt;
   double *bfact_r, *bfact_i;

   double gmin_spl_tmp,gmin_true_tmp;  /* Num : Min/Max g-vectors            */
   double gmax_spl_tmp;

   /*--------------------------------------*/
   /* Local pointers                       */

   int pme_on          = part_mesh->pme_on;
   int pme_res_on      = part_mesh->pme_res_on;
   int n_interp        = part_mesh->n_interp;
   int n_interp_res    = part_mesh->n_interp_res;
   int iperd           = cell->iperd;

   double *hmat_ewd    = cell->hmat_ewd;
   double *hmat_ewd_cp = cell->hmat_ewd_cp;

   int box_rat      =   cpewald->box_rat;
   double dbox_rat  =   cpewald->dbox_rat;

/*=======================================================================*/
/* 0) Output to screen                                                   */

   if(myid==0){
     printf("\n");PRINT_LINE_STAR;
     printf("Setting up reciprocal space\n");
     PRINT_LINE_DASH;printf("\n");
   }/*endif*/

/*=======================================================================*/
/* I) Set cp switch and initialize respa kvectors                        */

   cp_on = simopts->cp_min 
          +simopts->cp_wave_min
          +simopts->cp
          +simopts->cp_wave
          +simopts->debug_cp
          +simopts->cp_pimd
          +simopts->debug_cp_pimd+simopts->cp_wave_pimd
          +simopts->cp_wave_min_pimd;

   if(int_res_ter==0){ewald->nktot_res=0;ecor->nktot_res=0;}

/*=======================================================================*/
/* II) Allocate simple memory                                            */

   hmati_ewd      = (double *) cmalloc((size_t)9*sizeof(double))-1;
   hmati_ewd_cp   = (double *) cmalloc((size_t)9*sizeof(double))-1;
   kmaxv          =    (int *) cmalloc((size_t)3*sizeof(int))-1;
   kmaxv_res      =    (int *) cmalloc((size_t)3*sizeof(int))-1;
   kmax_cp_tmp    =    (int *) cmalloc((size_t)3*sizeof(int))-1;

   if(cp_on == 1){
     cpewald->kmax_cp = (int *) cmalloc((size_t)3*sizeof(int))-1;
     kmax_cp          = cpewald->kmax_cp;
     cpewald->kmax_cp_dens_cp_box = (int *) cmalloc((size_t)3*sizeof(int))-1;
     kmax_cp_dens_cp_box          = cpewald->kmax_cp_dens_cp_box;
     if(cp_dual_grid_opt_on >= 1){
       kmaxv_dens_cp_box = (int *) cmalloc((size_t)3*sizeof(int))-1;
     }/*endif*/
   }/* endif cp_on */


/*==========================================================================*/
/* III) Get inverse cell matrix and convert ewald_alpha                     */

   gethinv(hmat_ewd,hmati_ewd,&deth,iperd);    
   gethinv(hmat_ewd_cp,hmati_ewd_cp,&deth_cp,iperd);    

   side  = pow(deth,(1.0/3.0));  
   (ewald->alp_ewd) /= side;
  
/*==========================================================================*/
/* IV) Calculate cutoff, count number k vectors, Malloc and Fill            */

/*----------------------------------------------------------------------*/
/* A) Calculate cutoff, count number k vectors, malloc and fill        */
/*    Without the dual box this is standard grid for cp/ewald          */
/*    With the dual box this is the small box calculation              */

   calc_cutoff(kmax_ewd,&ecut_now,&(cp_parse->cp_ecut),cp_on,
                kmax_cp,kmaxv,hmati_ewd_cp,deth_cp);  

   countkvec3d(&(ewald->nktot),ecut_now,kmaxv,hmati_ewd_cp); 

   nktot                   = ewald->nktot;
   cpcoeffs_info->ecut     = ecut_now;
   cpcoeffs_info->ncoef_l  = nktot+1;
   ncoef_l                 = nktot+1;
   ecor->ecut              = 4.0*ecut_now;
   ewald->ecut             = 4.0*ecut_now;
   ewald->nkc_max          = kmaxv[3];

/*----------------------------------------------------------------------*/
/* A.1) For dualing : Calculate cutoff and count kvectors for the large */
/*      box and save the small box.                                     */
if(cp_on == 1){
   switch(cp_dual_grid_opt_on){
    case 0:  
      ecut_dens_cp_box_now   = ecut_now;
      kmax_cp_dens_cp_box[1] = kmax_cp[1];
      kmax_cp_dens_cp_box[2] = kmax_cp[2];
      kmax_cp_dens_cp_box[3] = kmax_cp[3];
    break;
    case 1:  
      ecut_dens_cp_box_now               = ecut_now;
      nktot_dens_cp_box                  = nktot;
      cpewald->nktot_dens_cp_box         = nktot;
      cpcoeffs_info->ncoef_l_dens_cp_box = nktot + 1;
      ncoef_dens_cp_box                  = nktot + 1;
      cpcoeffs_info->ecut_dens_cp_box    = ecut_now;

      kmaxv_dens_cp_box[1]   = kmaxv[1];
      kmaxv_dens_cp_box[2]   = kmaxv[2];
      kmaxv_dens_cp_box[3]   = kmaxv[3];

      kmax_cp_dens_cp_box[1] = kmax_cp[1];
      kmax_cp_dens_cp_box[2] = kmax_cp[2];
      kmax_cp_dens_cp_box[3] = kmax_cp[3];

      nk1 = 4*(kmax_cp[1]+1);
      nk2 = 4*(kmax_cp[2]+1);
      nk3 = 4*(kmax_cp[3]+1);
/*!!!!  DERIVE THIS EXPRESSION !!!! */
/* How mamy g-vectors do you need to get the density on the large grid */
/* correctly if you fft to g-space with a spherical cut-off and fft back */
/* Glenn's analysis indicates thou shalt not truncate g space  at all   */
/*    for this option , great since this makes the pme a more important */
/*    method  */

      kmaxv[1] = 2*(box_rat*(kmax_cp[1]+1)-1);
      kmaxv[2] = 2*(box_rat*(kmax_cp[2]+1)-1);
      kmaxv[3] = 2*(box_rat*(kmax_cp[3]+1)-1);

      countkvec3d(&(ewald->nktot),ecut_now,kmaxv,hmati_ewd);

      nktot                   = ewald->nktot;
      cpcoeffs_info->ecut     = ecut_now;
      cpcoeffs_info->ncoef_l  = nktot+1;
      ncoef_l                 = nktot+1;
      ecor->ecut              = 4.0*ecut_now;
      ewald->ecut             = 4.0*ecut_now;
      ewald->nkc_max          = kmaxv[3];
    break;
#ifdef  ORIG
    case 2:
      ecut_dens_cp_box_now               = ecut_now;
      nktot_dens_cp_box                  = nktot;
      cpewald->nktot_dens_cp_box         = nktot;
      cpcoeffs_info->ncoef_l_dens_cp_box = nktot + 1;
      ncoef_dens_cp_box                  = nktot + 1;
      cpcoeffs_info->ecut_dens_cp_box    = ecut_now;

      kmaxv_dens_cp_box[1]   = kmaxv[1];
      kmaxv_dens_cp_box[2]   = kmaxv[2];
      kmaxv_dens_cp_box[3]   = kmaxv[3];

      kmax_cp_dens_cp_box[1] = kmax_cp[1];
      kmax_cp_dens_cp_box[2] = kmax_cp[2];
      kmax_cp_dens_cp_box[3] = kmax_cp[3];

      nk1 = 4*(kmax_cp[1]+1);
      nk2 = 4*(kmax_cp[2]+1);
      nk3 = 4*(kmax_cp[3]+1);
      kmaxv[1] = 2*(box_rat*(kmax_cp[1]+1)-1);
      kmaxv[2] = 2*(box_rat*(kmax_cp[2]+1)-1);
      kmaxv[3] = 2*(box_rat*(kmax_cp[3]+1)-1);

      countkvec3d(&(ewald->nktot),ecut_now,kmaxv,hmati_ewd);

      nktot                   = ewald->nktot;
      cpcoeffs_info->ecut     = ecut_now;
      cpcoeffs_info->ncoef_l  = nktot+1;
      ncoef_l                 = nktot+1;
      ecor->ecut              = 4.0*ecut_now;
      ewald->ecut             = 4.0*ecut_now;
      ewald->nkc_max          = kmaxv[3];
    break;
#endif
#ifdef  PME
    case 2:
      ecut_dens_cp_box_now               = ecut_now;
      nktot_dens_cp_box                  = nktot;
      cpewald->nktot_dens_cp_box         = nktot;
      cpcoeffs_info->ncoef_l_dens_cp_box = nktot + 1;
      ncoef_dens_cp_box                  = nktot + 1;
      cpcoeffs_info->ecut_dens_cp_box    = ecut_now;

      kmaxv_dens_cp_box[1]   = kmaxv[1];
      kmaxv_dens_cp_box[2]   = kmaxv[2];
      kmaxv_dens_cp_box[3]   = kmaxv[3];

      kmax_cp_dens_cp_box[1] = kmax_cp[1];
      kmax_cp_dens_cp_box[2] = kmax_cp[2];
      kmax_cp_dens_cp_box[3] = kmax_cp[3];

      calc_cutoff(kmax_ewd,&ecut_now,&(cp_parse->cp_ecut_dual_grid),cp_on,
                  kmax_cp,kmaxv,hmati_ewd,deth);  
      countkvec3d(&(ewald->nktot),ecut_now,kmaxv,hmati_ewd); 

      nktot                   = ewald->nktot;
      cpcoeffs_info->ecut     = ecut_now;
      cpcoeffs_info->ncoef_l  = nktot+1;
      ncoef_l                 = nktot+1;
      ecor->ecut              = 4.0*ecut_now;
      ewald->ecut             = 4.0*ecut_now;
      ewald->nkc_max          = kmaxv[3];
    break;
#endif
   }/*end switch*/
}/*endif cp_on*/
 
/*----------------------------------------------------------------------*/
/* B) Malloc                                                            */

   nmall =  nktot+1;   if((nmall % 2)==0){nmall++;}
   ewald->nktot_mall = nmall;
   now_mem    = (nmall*(sizeof(double)*0 + sizeof(int)*5))*1.e-06;
   if(int_res_ter!=0){
      now_mem  += (nmall*(sizeof(double)*0 + sizeof(int)*6 ))*1.e-06;
   }/*endif*/
   (*tot_memory) += now_mem;

   ewald->kastr = (int *) cmalloc(nmall*sizeof(int))-1;
   ewald->kbstr = (int *) cmalloc(nmall*sizeof(int))-1;
   ewald->kcstr = (int *) cmalloc(nmall*sizeof(int))-1;
   ewald->ibrk1 = (int *) cmalloc(nmall*sizeof(int))-1;
   ewald->ibrk2 = (int *) cmalloc(nmall*sizeof(int))-1;
   if(int_res_ter != 0) {
      ewald->nktot_res_mall = nmall;
      ewald->kastr_res      = (int *) cmalloc(nmall*sizeof(int))-1;
      ewald->kbstr_res      = (int *) cmalloc(nmall*sizeof(int))-1;
      ewald->kcstr_res      = (int *) cmalloc(nmall*sizeof(int))-1;
      ewald->ibrk1_res      = (int *) cmalloc(nmall*sizeof(int))-1;
      ewald->ibrk2_res      = (int *) cmalloc(nmall*sizeof(int))-1;
      ewald->ibrk3          = (int *) cmalloc(nmall*sizeof(int))-1;
   }/*endif*/

   if(myid==0){
      printf("Ewald allocation: %g Mbytes; Total memory %g Mbytes\n",
              now_mem,*tot_memory);
   }/*endif*/

/*------------------------------------------------------------------------*/

   if( cp_dual_grid_opt_on >= 1 && cp_on == 1){
     nmall      =  nktot_dens_cp_box+1;   if((nmall % 2)==0){nmall++;}
     now_mem    = (nmall*(sizeof(double)*0 + sizeof(int)*5))*1.e-06;
     (*tot_memory) += now_mem;

     cpewald->kastr_dens_cp_box    = (int *) cmalloc(nmall*sizeof(int))-1;
     cpewald->kbstr_dens_cp_box    = (int *) cmalloc(nmall*sizeof(int))-1;
     cpewald->kcstr_dens_cp_box    = (int *) cmalloc(nmall*sizeof(int))-1;
     cpewald->ibrk1_dens_cp_box    = (int *) cmalloc(nmall*sizeof(int))-1;
     cpewald->ibrk2_dens_cp_box    = (int *) cmalloc(nmall*sizeof(int))-1;

     if(myid==0){
       printf("cp grid allocation: %g Mbytes; Total memory %g Mbytes\n",
               now_mem,*tot_memory);
     }/*endif*/
   }/*endif cp_dual_grid_opt_on*/

/*------------------------------------------------------------------------*/
/* C) Fill                                                                */
   setkvec3d(nktot,ecut_now,kmaxv,hmati_ewd,
             ewald->kastr,ewald->kbstr,ewald->kcstr,
             ewald->ibrk1,ewald->ibrk2,cp_on,
             gmin_spl,gmin_true,gmax_spl);
/*------------------------------------------------------------------------*/
/* C) Fill DENS_CP_BOX                                                    */

   if( cp_dual_grid_opt_on >= 1 && cp_on == 1){
     setkvec3d(nktot_dens_cp_box,ecut_dens_cp_box_now,
               kmaxv_dens_cp_box,hmati_ewd_cp,
               cpewald->kastr_dens_cp_box,
               cpewald->kbstr_dens_cp_box,cpewald->kcstr_dens_cp_box,
               cpewald->ibrk1_dens_cp_box,cpewald->ibrk2_dens_cp_box,cp_on,
               &gmin_spl_tmp,&gmin_true_tmp,&gmax_spl_tmp);
   }/*endif cp_dual_grid_opt_on*/

   if( cp_dual_grid_opt_on == 2 && cp_on == 1){
     *gmin_spl  = gmin_spl_tmp;
     *gmin_true = gmin_true_tmp;
     *gmax_spl  = gmax_spl_tmp;
   }/*endif*/

/*=======================================================================*/
/* V) Setup PME                                                          */

   if(pme_on==1){
/*-----------------------------------------------------------------------*/
/*  A) Get grid size */

      part_mesh->nktot_pme = nktot;
      set_pme_grid(ecut_now,deth,hmati_ewd,kmaxv,
                  &(part_mesh->ngrid_a),&(part_mesh->ngrid_b),
                  &(part_mesh->ngrid_c),n_interp,
                  part_mesh->kmax_pme);
      part_mesh->ecut     = 4.0*ecut_now;
      ngrid_a = part_mesh->ngrid_a;
      ngrid_b = part_mesh->ngrid_b;
      ngrid_c = part_mesh->ngrid_c;
      if((ngrid_a<n_interp) || (ngrid_b<n_interp) || (ngrid_c<n_interp)){
       if(myid==0){
        printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n");
        printf("The PME n_interp parameter > number of grid points \n");
        printf("This is not allowed\n");      
        printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n");
       }/*endif*/
       fflush(stdout);
       exit(1);      
      }/*endif*/
/*-----------------------------------------------------------------------*/
/*  B) Malloc */
      now_mem = ( (nktot)*(sizeof(double))
                +3*(n_interp) *(sizeof(double))  )*1.e-06;
      if(pme_res_on == 1) {
        now_mem += ((nktot) *(sizeof(double)) )*1.e-06;
      }/*endif*/
     *tot_memory += now_mem;

   /*------------------------*/
   /*  Malloc the bweights   */
     nmall                     =  nktot;
     part_mesh->bweight_tot    = (double *) cmalloc(nmall*sizeof(double))-1;
     bfact_r =  part_mesh->bweight_tot;
     bfact_i =  part_mesh->bweight_tot;
     if(pme_res_on == 1) {
      part_mesh->bweight_tot_res = (double *) cmalloc(nmall*sizeof(double))-1;
     }/*endif*/

   /*------------------------*/
   /*  Malloc some scratch  */
      nmall                          = n_interp;
      part_mesh->ninterp_tot_mall    = n_interp;

      part_mesh->aj    = (double *)cmalloc(nmall*sizeof(double))-1;
      part_mesh->rn    = (double *)cmalloc(nmall*sizeof(double))-1;
      part_mesh->rn1   = (double *)cmalloc(nmall*sizeof(double))-1;

   /*-----------------------*/
   /*  Output               */
     if(myid==0){
      printf("PME allocation: %g Mbytes; Total memory %g Mbytes\n",
                                           now_mem,*tot_memory);
     }/*endif*/
/*-----------------------------------------------------------------------*/
/*  C) Create PME Bweight */
      pme_b_opt = 1;
      set_pme_wght(nktot,ewald->kastr,ewald->kbstr,ewald->kcstr,
                   ngrid_a,ngrid_b,ngrid_c,
                   idum1,idum2,idum3,
                   pme_b_opt,bfact_r,bfact_i,
                   part_mesh->bweight_tot,n_interp,
                   part_mesh->aj,part_mesh->rn,part_mesh->rn1);
   }/*endif*/

/*=======================================================================*/
/* VI) Set up RESPA k vectors and RESPA PME                              */

   if(int_res_ter != 0) {
/*-----------------------------------------------------------------------*/
/*  A) Repsa K vectors */
      ecut_tmp = 0.0;cp_on_tmp=0;

      calc_cutoff(kmax_res,&ecut_res,&ecut_tmp,cp_on_tmp,
                  kmax_cp_tmp,kmaxv_res,hmati_ewd,deth);

      ecor->ecut_res   = 4.0*ecut_res;
      ewald->ecut_res  = 4.0*ecut_res;
      countkvec3d(&(ewald->nktot_res),ecut_res,kmaxv_res,hmati_ewd);
      nktot_res         = ewald->nktot_res;
      ecor->nktot_res   = nktot_res;
      setkvec3d(nktot_res,ecut_res,kmaxv_res,hmati_ewd,
                ewald->kastr_res,ewald->kbstr_res,ewald->kcstr_res,
                ewald->ibrk1_res,ewald->ibrk2_res,cp_on_tmp,
                &gmin_spl_tmp,&gmin_true_tmp,&gmax_spl_tmp);
      setkvec3d_res(kmax_res,hmati_ewd,
                    ewald->kastr,ewald->kbstr,ewald->kcstr,ewald->ibrk3,
                    nktot,nktot_res);
      if(pme_res_on==1){
/*-----------------------------------------------------------------------*/
/*  B) Set the PME GRID */
       part_mesh->nktot_pme_res = nktot_res;
       set_pme_grid(ecut_res,deth,hmati_ewd,kmaxv_res,
                  &(part_mesh->ngrid_a_res),&(part_mesh->ngrid_b_res),
                  &(part_mesh->ngrid_c_res),n_interp_res,
                  part_mesh->kmax_pme_res);
       part_mesh->ecut_res     = 4.0*ecut_res;
       ngrid_a_res = part_mesh->ngrid_a_res;
       ngrid_b_res = part_mesh->ngrid_b_res;
       ngrid_c_res = part_mesh->ngrid_c_res;
       if((ngrid_a_res<n_interp_res) || (ngrid_b_res<n_interp_res) || 
          (ngrid_c_res<n_interp_res)){
        if(myid==0){
         printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n");
         printf("The RESPA PME n_interp parameter > number of grid points \n");
         printf("This is not allowed\n");      
         printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n");
        }/*endif*/
        fflush(stdout);
        exit(1);      
       }/*endif*/
/*-----------------------------------------------------------------------*/
/*  C) Create PME Bweight */
       pme_b_opt = 1;
       set_pme_wght(nktot_res,ewald->kastr_res,ewald->kbstr_res,
                  ewald->kcstr_res,ngrid_a_res,ngrid_b_res,ngrid_c_res,
                  idum1,idum2,idum3,
                  pme_b_opt,bfact_r,bfact_i,
                  part_mesh->bweight_tot_res,n_interp_res,
                  part_mesh->aj,part_mesh->rn,part_mesh->rn1);
      }/*endif*/
    }/*endif*/

/*=======================================================================*/
/* VII) Setup CP coefs and k vectors                                     */

   if(cp_on == 1) {

/*--------------------------------------------------------------------*/
/*  A)  Count the k-vectors                                           */
      ecut_sm = ecut_dens_cp_box_now;
      countkvec3d_sm(&(cpewald->nktot_sm),ecut_sm,
                     kmax_cp_dens_cp_box,hmati_ewd_cp);
      nktot_sm = cpewald->nktot_sm;
      cpcoeffs_info->ncoef   = nktot_sm+1;
      ncoef                  = nktot_sm+1;
/*--------------------------------------------------------------------*/
/*  B)  Malloc                                                       */
      nmall =  (nktot_sm+1);if((nmall % 2)==0){nmall++;}
      cpewald->nktot_cp_sm_mall = nmall;
      now_mem = (nmall*(sizeof(double)*0 + sizeof(int)*5 ))*1.e-06;
      *tot_memory += now_mem;
      cpewald->kastr_sm = (int *) cmalloc(nmall*sizeof(int))-1;
      cpewald->kbstr_sm = (int *) cmalloc(nmall*sizeof(int))-1;
      cpewald->kcstr_sm = (int *) cmalloc(nmall*sizeof(int))-1;
      cpewald->ibrk1_sm = (int *) cmalloc(nmall*sizeof(int))-1;
      cpewald->ibrk2_sm = (int *) cmalloc(nmall*sizeof(int))-1;

      nmall =  ncoef; if((nmall % 2)==0){nmall++;}
      cpcoeffs_info->cmass = (double *)cmalloc(nmall*sizeof(double))-1;

      if(myid==0){
        printf("CP allocation: %g Mbytes; Total memory %g Mbytes\n",
                 now_mem,*tot_memory);
      }/*endif*/
/*--------------------------------------------------------------------*/
/*  C)  Fill and check                                                */

      setkvec3d_sm(nktot_sm,ecut_sm,kmax_cp_dens_cp_box,hmati_ewd_cp,
                   cpewald->kastr_sm,cpewald->kbstr_sm,cpewald->kcstr_sm,
                   cpewald->ibrk1_sm,cpewald->ibrk2_sm,
                   &(cpewald->gw_gmin),&(cpewald->gw_gmax));

    if(cp_dual_grid_opt_on == 0 && cp_on == 1){
      check_kvec(ewald->nktot,ewald->kastr,ewald->kbstr,ewald->kcstr,nktot_sm,
                  cpewald->kastr_sm,cpewald->kbstr_sm,cpewald->kcstr_sm);
    }

      if(cp_dual_grid_opt_on >= 1 && cp_on == 1){
        check_kvec(cpewald->nktot_dens_cp_box,cpewald->kastr_dens_cp_box,
                   cpewald->kbstr_dens_cp_box,
                   cpewald->kcstr_dens_cp_box,nktot_sm,
                   cpewald->kastr_sm,cpewald->kbstr_sm,cpewald->kcstr_sm);
      }/*endif cp_dual_grid_opt_on */   
/*--------------------------------------------------------------------*/
/*  D) Set up the cp masses                                           */

      set_cpmass(ncoef,cpewald->kastr_sm,
                 cpewald->kbstr_sm,cpewald->kcstr_sm,
                 cpcoeffs_info->cmass,hmati_ewd_cp,
                 &(cp_parse->cp_mass_tau_def),cp_parse->cp_mass_cut_def,
                 &(cpcoeffs_info->icmass_unif));

   }/*endif:cpon*/

/*=======================================================================*/
/* VIII) Output time has arrived                                         */

   if(myid ==0){
 /*-----------------------------------------------------------------------*/
 /* A) CP output */
     if(cp_on == 1) {
       switch(cp_dual_grid_opt_on){
         case 0:
           nkf1 = 4*(kmax_cp[1]+1);
           nkf2 = 4*(kmax_cp[2]+1);
           nkf3 = 4*(kmax_cp[3]+1);
         break;
         case 1:
           nkf1 = 4*box_rat*(kmax_cp_dens_cp_box[1] + 1);
           nkf2 = 4*box_rat*(kmax_cp_dens_cp_box[2] + 1);
           nkf3 = 4*box_rat*(kmax_cp_dens_cp_box[3] + 1);

           nkf1_dens_cp_box = 4*(kmax_cp_dens_cp_box[1] + 1);
           nkf2_dens_cp_box = 4*(kmax_cp_dens_cp_box[2] + 1);
           nkf3_dens_cp_box = 4*(kmax_cp_dens_cp_box[3] + 1);
          break;
#ifdef ORIG
         case 2:
           nkf1 = 4*box_rat*(kmax_cp_dens_cp_box[1] + 1);
           nkf2 = 4*box_rat*(kmax_cp_dens_cp_box[2] + 1);
           nkf3 = 4*box_rat*(kmax_cp_dens_cp_box[3] + 1);

           nkf1_dens_cp_box = 4*(kmax_cp_dens_cp_box[1] + 1);
           nkf2_dens_cp_box = 4*(kmax_cp_dens_cp_box[2] + 1);
           nkf3_dens_cp_box = 4*(kmax_cp_dens_cp_box[3] + 1);
         break;
#endif
#ifdef  PME
         case 2:
           nkf1 = 4*(kmax_cp[1] + 1);
           nkf2 = 4*(kmax_cp[2] + 1);
           nkf3 = 4*(kmax_cp[3] + 1);

           nkf1_dens_cp_box = 4*(kmax_cp_dens_cp_box[1] + 1);
           nkf2_dens_cp_box = 4*(kmax_cp_dens_cp_box[2] + 1);
           nkf3_dens_cp_box = 4*(kmax_cp_dens_cp_box[3] + 1);
         break;
#endif
       }/*end switch*/

       nk1  = 2*kmax_cp_dens_cp_box[1];
       nk2  = 2*kmax_cp_dens_cp_box[2];
       nk3  = 2*kmax_cp_dens_cp_box[3];

       ecut_lg = 4.0*ecut_now;
       printf("Your large cp-fft grid is  %d by %d by %d\n",nkf1,nkf2,nkf3);
       printf("There are  %d total k-vectors ",ncoef_l); 
       printf("upon spherical truncation. \n");
       printf("The large energy cutoff is 4*Ecut= %f Ryd\n",2.0*ecut_lg);
       putchar('\n');

       if(cp_dual_grid_opt_on >= 1){
         printf("Your large density grid for the cp-fft box is %d by %d by %d",
                 nkf1_dens_cp_box,nkf2_dens_cp_box,nkf3_dens_cp_box);
         printf("\nThere are  %d total k-vectors ",ncoef_dens_cp_box); 
         printf("upon spherical truncation. \n");
         printf("The large energy cutoff is 4*Ecut= %f Ryd\n",8.0*ecut_sm);
         putchar('\n');
       }/*endif cp_dual_grid_opt_on */

       printf("Your small cp-fft grid is  %d by %d by %d\n",nk1,nk2,nk3);
       printf("There are %d total k-vectors ",ncoef); 
       printf("upon spherical truncation. \n");
       printf("The small energy cutoff is Ecut= %f Ryd\n",2.0*ecut_sm);

     }else {
 /*-----------------------------------------------------------------------*/
 /* B) Non-CP output */
       printf("You are using %d k-vectors in your ewald sum\n",nktot);
       printf("Your reciprocal space shape: (-%d,%d) by (-%d,%d) by (-%d,%d)",
              kmaxv[1],kmaxv[1],kmaxv[2],kmaxv[2],kmaxv[3],kmaxv[3]);
       printf("\n");
     } /* endif:cp_on */

 /*-----------------------------------------------------------------------*/
 /* C) PME output  */
     if(pme_on==1){
       printf("Your particle mesh grid shape: %d by %d by %d\n",
                                                   ngrid_a,ngrid_b,ngrid_c);
       if((pme_res_on==1)&&(int_res_ter != 0)&&(kmax_res > 0)){
         printf("Your respa particle mesh grid shape: %d by %d by %d\n",
                                       ngrid_a_res,ngrid_b_res,ngrid_c_res);
       }/*endif*/
     }/*endif*/

     putchar('\n');PRINT_LINE_DASH;
     printf("Completed reciprocal space set up\n");
     PRINT_LINE_STAR;printf("\n");

   }/*endif:myid==0*/

/*=======================================================================*/
/* IX) Free excess memory                                                */

   cfree(&(hmati_ewd)[1]);
   cfree(&(hmati_ewd_cp)[1]);
   cfree(&(kmaxv)[1]);
   cfree(&(kmaxv_res)[1]);
   cfree(&(kmax_cp_tmp)[1]);

/*-----------------------------------------------------------------------*/ 
  }/* end routine */
/*=======================================================================*/




/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
/* Control the setup of the FFT packages                                    */
/*==========================================================================*/

void control_fft_pkg(PARA_FFT_PKG3D *cp_sclr_fft_pkg_sm,
                     PARA_FFT_PKG3D *cp_para_fft_pkg_sm,
                     PARA_FFT_PKG3D *cp_sclr_fft_pkg_dens_cp_box,
                     PARA_FFT_PKG3D *cp_para_fft_pkg_dens_cp_box,
                     PARA_FFT_PKG3D *cp_sclr_fft_pkg_lg,
                     PARA_FFT_PKG3D *cp_para_fft_pkg_lg,
                     PARA_FFT_PKG3D *pme_fft_pkg,
                     PARA_FFT_PKG3D *pme_res_fft_pkg,
                     EWALD* ewald,CPEWALD *cpewald,
                     PART_MESH *part_mesh, CPCOEFFS_INFO *cpcoeffs_info,
                     COMMUNICATE *communicate,int cp_on,int cp_lsda,
                     double *tot_memory,int int_res_ter,
                     int cp_para_opt,int cp_dual_grid_opt_on)

/*=========================================================================*/
     {/*Begin Routine*/
/*=========================================================================*/
/*    Local Variables */

  int nkf1,nkf2,nkf3,nfft_ext,iii;
  int nkf1_dens_cp_box,nkf2_dens_cp_box,nkf3_dens_cp_box;
          /* used for denisty on cp grid when have 2 boxes*/

  int nstate_up          = cpcoeffs_info->nstate_up;
  int nstate_dn          = cpcoeffs_info->nstate_dn;
  int ncoef              = cpcoeffs_info->ncoef;
  int ncoef_dens_cp_box  = cpcoeffs_info->ncoef_l_dens_cp_box;
  int ncoef_l            = cpcoeffs_info->ncoef_l;
  int nktot              = ewald->nktot;
  int nktot_dens_cp_box  = cpewald->nktot_dens_cp_box;
  int box_rat            = cpewald->box_rat;
  int nktot_res          = ewald->nktot_res;
  int *kmax_cp           = cpewald->kmax_cp;
  int *kmax_cp_dens_cp_box = cpewald->kmax_cp_dens_cp_box;
  int myid               = communicate->myid;
  int np_states          = communicate->np_states;
  int myid_state         = communicate->myid_state;
  int pme_on             = part_mesh->pme_on;
  int pme_res_on         = part_mesh->pme_res_on;
  int myid_forc          = communicate->myid_forc;
  int np_forc            = communicate->np_forc;
  int ngrid_a            = part_mesh->ngrid_a;
  int ngrid_b            = part_mesh->ngrid_b;
  int n_interp           = part_mesh->n_interp;
  int n_interp_res       = part_mesh->n_interp_res;
  int pme_para_opt       = part_mesh->pme_para_opt;

/*=========================================================================*/
/* 0) Print to screen and check for nproc > nstate error */

  if(myid == 0){
    printf("\n");PRINT_LINE_STAR;
    printf("Setting up FFTs\n");
    PRINT_LINE_DASH;printf("\n");
  }/* endif myid */

  if(cp_on == 1){
    if(cp_lsda == 0){ 
     if(nstate_up < np_states){
      printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@@@@@@\n");
       printf("Number of states less than number of processors\n");
       printf("If possible, reduce number of processors to be\n");
       printf("less than the number of states or run a bigger system.\n");
       printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@@@@@@\n");
       fflush(stdout);
       exit(1);
     }/* endif */
    } else {
     if(nstate_up + nstate_dn < np_states){
      printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@@@@@@\n");
       printf("Number of states less than number of processors\n");
       printf("If possible, reduce number of processors to be\n");
       printf("less than the number of states or run a bigger system.\n");
       printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@@@@@@\n");
       fflush(stdout);
       exit(1);
     }/* endif */
    }/* endif lsda */
  }/* endif */

/*=========================================================================*/
/* 0.1) Set CP FFT Size  */

  if(cp_on==1){
    nkf1 = 4*(kmax_cp_dens_cp_box[1]+1);
    nkf2 = 4*(kmax_cp_dens_cp_box[2]+1);
    nkf3 = 4*(kmax_cp_dens_cp_box[3]+1);
  }/*endif*/

/*=========================================================================*/
/* I) DENS_CP_BOX CP scalar package                                        */

 if(cp_dual_grid_opt_on >= 1 && cp_para_opt == 0){/* hybrid option */

    cp_sclr_fft_pkg_dens_cp_box->nkf1       = nkf1;
    cp_sclr_fft_pkg_dens_cp_box->nkf2       = nkf2;
    cp_sclr_fft_pkg_dens_cp_box->nkf3       = nkf3;
      
    cp_sclr_fft_pkg_dens_cp_box->nktot      = nktot_dens_cp_box;
    cp_sclr_fft_pkg_dens_cp_box->ncoef      = ncoef_dens_cp_box;

    cp_sclr_fft_pkg_dens_cp_box->myid       = 0;
    cp_sclr_fft_pkg_dens_cp_box->myidp1     = 1;
    cp_sclr_fft_pkg_dens_cp_box->num_proc   = 1;
    cp_sclr_fft_pkg_dens_cp_box->comm       = communicate->comm_faux;


    create_para_fft_pkg3d(cp_sclr_fft_pkg_dens_cp_box,
                          cpewald->kastr_dens_cp_box,
                          cpewald->kbstr_dens_cp_box,
                          cpewald->kcstr_dens_cp_box,cp_dual_grid_opt_on);
  }/*endif*/


/*=========================================================================*/
/* II) DENSITY_CP_BOX  parallel package                                    */
/*       This package must be made for both hybrid and full_g options      */

  if(cp_dual_grid_opt_on >= 1){

    cp_para_fft_pkg_dens_cp_box->nkf1       = nkf1;
    cp_para_fft_pkg_dens_cp_box->nkf2       = nkf2;
    cp_para_fft_pkg_dens_cp_box->nkf3       = nkf3;
      
    cp_para_fft_pkg_dens_cp_box->nktot      = nktot_dens_cp_box;
    cp_para_fft_pkg_dens_cp_box->ncoef      = ncoef_dens_cp_box;
   
    cp_para_fft_pkg_dens_cp_box->myid       = myid_state;
    cp_para_fft_pkg_dens_cp_box->myidp1     = myid_state+1;
    cp_para_fft_pkg_dens_cp_box->num_proc   = np_states;
    cp_para_fft_pkg_dens_cp_box->comm       = communicate->comm_states;

    create_para_fft_pkg3d(cp_para_fft_pkg_dens_cp_box,
                          cpewald->kastr_dens_cp_box,
                          cpewald->kbstr_dens_cp_box,
                          cpewald->kcstr_dens_cp_box,cp_dual_grid_opt_on);


  }/*endif*/

/*=========================================================================*/
/* I) Large CP scalar package                                              */


 if(cp_on == 1 && cp_para_opt == 0){/* hybrid option */

   switch(cp_dual_grid_opt_on){
    case 0:
     nkf1 = 4*(kmax_cp[1]+1);
     nkf2 = 4*(kmax_cp[2]+1);
     nkf3 = 4*(kmax_cp[3]+1);
    break;
    case 1:
     nkf1 = 4*box_rat*(kmax_cp_dens_cp_box[1]+1);
     nkf2 = 4*box_rat*(kmax_cp_dens_cp_box[2]+1);
     nkf3 = 4*box_rat*(kmax_cp_dens_cp_box[3]+1);
    break;

#ifdef ORIG
    case 2:
     nkf1 = 4*box_rat*(kmax_cp_dens_cp_box[1]+1);
     nkf2 = 4*box_rat*(kmax_cp_dens_cp_box[2]+1);
     nkf3 = 4*box_rat*(kmax_cp_dens_cp_box[3]+1);
    break;
#endif
#ifdef  PME
    case 2:
     nkf1 = 4*(kmax_cp[1]+1);
     nkf2 = 4*(kmax_cp[2]+1);
     nkf3 = 4*(kmax_cp[3]+1);
    break;
#endif
   }/*end switch */

    cp_sclr_fft_pkg_lg->nkf1       = nkf1;
    cp_sclr_fft_pkg_lg->nkf2       = nkf2;
    cp_sclr_fft_pkg_lg->nkf3       = nkf3;
      
    cp_sclr_fft_pkg_lg->nktot      = nktot;
    cp_sclr_fft_pkg_lg->ncoef      = ncoef_l;

    cp_sclr_fft_pkg_lg->myid       = 0;
    cp_sclr_fft_pkg_lg->myidp1     = 1;
    cp_sclr_fft_pkg_lg->num_proc   = 1;
    cp_sclr_fft_pkg_lg->comm       = communicate->comm_faux;

    create_para_fft_pkg3d(cp_sclr_fft_pkg_lg,
                          ewald->kastr,ewald->kbstr,
                          ewald->kcstr,cp_dual_grid_opt_on);
  }/*endif*/


/*=========================================================================*/
/* II) Large CP parallel package                                           */

  if(cp_on == 1){

   switch(cp_dual_grid_opt_on){
    case 0:
     nkf1 = 4*(kmax_cp[1]+1);
     nkf2 = 4*(kmax_cp[2]+1);
     nkf3 = 4*(kmax_cp[3]+1);
    break;
    case 1:
     nkf1 = 4*box_rat*(kmax_cp_dens_cp_box[1]+1);
     nkf2 = 4*box_rat*(kmax_cp_dens_cp_box[2]+1);
     nkf3 = 4*box_rat*(kmax_cp_dens_cp_box[3]+1);
    break;
#ifdef ORIG
    case 2:
     nkf1 = 4*box_rat*(kmax_cp_dens_cp_box[1]+1);
     nkf2 = 4*box_rat*(kmax_cp_dens_cp_box[2]+1);
     nkf3 = 4*box_rat*(kmax_cp_dens_cp_box[3]+1);
    break;
#endif
#ifdef  PME
    case 2:
     nkf1 = 4*(kmax_cp[1]+1);
     nkf2 = 4*(kmax_cp[2]+1);
     nkf3 = 4*(kmax_cp[3]+1);
    break;
#endif
   }/*end switch */

    cp_para_fft_pkg_lg->nkf1       = nkf1;
    cp_para_fft_pkg_lg->nkf2       = nkf2;
    cp_para_fft_pkg_lg->nkf3       = nkf3;
      
    cp_para_fft_pkg_lg->nktot      = nktot;
    cp_para_fft_pkg_lg->ncoef      = ncoef_l;
   
    cp_para_fft_pkg_lg->myid       = myid_state;
    cp_para_fft_pkg_lg->myidp1     = myid_state+1;
    cp_para_fft_pkg_lg->num_proc   = np_states;
    cp_para_fft_pkg_lg->comm       = communicate->comm_states;

    create_para_fft_pkg3d(cp_para_fft_pkg_lg,
                          ewald->kastr,ewald->kbstr,
                          ewald->kcstr,cp_dual_grid_opt_on);
  }/*endif*/

/*=========================================================================*/
/* III) Small  CP scalar package                                           */

  if(cp_on == 1 && cp_para_opt == 0){/* hybrid option */

    nkf1 = 4*(kmax_cp_dens_cp_box[1]+1);
    nkf2 = 4*(kmax_cp_dens_cp_box[2]+1);
    nkf3 = 4*(kmax_cp_dens_cp_box[3]+1);

    cp_sclr_fft_pkg_sm->nkf1       = nkf1;
    cp_sclr_fft_pkg_sm->nkf2       = nkf2;
    cp_sclr_fft_pkg_sm->nkf3       = nkf3;
      
    cp_sclr_fft_pkg_sm->nktot      = ncoef-1;
    cp_sclr_fft_pkg_sm->ncoef      = ncoef;

    cp_sclr_fft_pkg_sm->myid       = 0;
    cp_sclr_fft_pkg_sm->myidp1     = 1;
    cp_sclr_fft_pkg_sm->num_proc   = 1;
    cp_sclr_fft_pkg_sm->comm       = communicate->comm_faux;

    create_para_fft_pkg3d(cp_sclr_fft_pkg_sm,
                          cpewald->kastr_sm,cpewald->kbstr_sm,
                          cpewald->kcstr_sm,cp_dual_grid_opt_on);

  }/*endif*/

/*=========================================================================*/
/* IV) Small  CP parallel package                                         */

  if(cp_on == 1 && cp_para_opt == 1){/* full g option */

    nkf1 = 4*(kmax_cp_dens_cp_box[1]+1);
    nkf2 = 4*(kmax_cp_dens_cp_box[2]+1);
    nkf3 = 4*(kmax_cp_dens_cp_box[3]+1);

    cp_para_fft_pkg_sm->nkf1       = nkf1;
    cp_para_fft_pkg_sm->nkf2       = nkf2;
    cp_para_fft_pkg_sm->nkf3       = nkf3;
      
    cp_para_fft_pkg_sm->nktot      = ncoef-1;
    cp_para_fft_pkg_sm->ncoef      = ncoef;

    cp_para_fft_pkg_sm->myid       = myid_state;
    cp_para_fft_pkg_sm->myidp1     = myid_state+1;
    cp_para_fft_pkg_sm->num_proc   = np_states;
    cp_para_fft_pkg_sm->comm       = communicate->comm_states;

    create_para_fft_pkg3d(cp_para_fft_pkg_sm,
                         cpewald->kastr_sm,cpewald->kbstr_sm,
                          cpewald->kcstr_sm,cp_dual_grid_opt_on);

  }/*endif*/

/*=========================================================================*/
/* V) PME package                                                         */

  if(cp_on==0 && pme_on ==1 ){

    pme_fft_pkg->nkf1       = ngrid_a;
    pme_fft_pkg->nkf2       = ngrid_b;
    pme_fft_pkg->nkf3       = part_mesh->ngrid_c;
      
    pme_fft_pkg->nktot      = ewald->nktot;
    pme_fft_pkg->ncoef      = ewald->nktot+1;

    if(pme_para_opt==2){
     pme_fft_pkg->myid       = myid_forc;
     pme_fft_pkg->myidp1     = myid_forc+1;
     pme_fft_pkg->num_proc   = np_forc;
     pme_fft_pkg->comm       = communicate->comm_forc;
   }else{
     pme_fft_pkg->myid       = 0;
     pme_fft_pkg->myidp1     = 1;
     pme_fft_pkg->num_proc   = 1;
     pme_fft_pkg->comm       = communicate->comm_faux;
   }/*endif*/

    create_para_fft_pkg3d(pme_fft_pkg,
                          ewald->kastr,ewald->kbstr,
                          ewald->kcstr,cp_dual_grid_opt_on);

    pme_fft_pkg->scale_opt = 0;
#ifdef HP_VECLIB
    if(pme_fft_pkg->igeneric_opt==0){pme_fft_pkg->scale_opt = -1;}
#endif

    if(pme_para_opt==2){
      nfft_ext  = 2*ngrid_a*(  (pme_fft_pkg->nfft_ka_proc)
                             + (pme_fft_pkg->skb_fft_ka_proc) - 1
                             + (pme_fft_pkg->ekb_fft_ka_proc) - ngrid_b
                             + (n_interp-1)*ngrid_b );
      pme_fft_pkg->nfft_size = MAX(pme_fft_pkg->nfft_size,nfft_ext);

      if(np_forc>1){create_pme_comm_full_g(n_interp,pme_fft_pkg);}
    }/*endif*/

    if((pme_para_opt==1)&&(np_forc>1)){
       create_pme_comm_hybr(np_forc,myid_forc,communicate->comm_forc,
                             pme_fft_pkg);
    }/*endif*/

  }/*endif*/


/*=========================================================================*/
/* VI) PME_RES package                                                      */

  if(cp_on==0 && int_res_ter == 1 && pme_res_on==1 && nktot_res > 0){

    pme_res_fft_pkg->nkf1       = part_mesh->ngrid_a_res;
    pme_res_fft_pkg->nkf2       = part_mesh->ngrid_b_res;
    pme_res_fft_pkg->nkf3       = part_mesh->ngrid_c_res;
      
    pme_res_fft_pkg->nktot      = nktot_res;
    pme_res_fft_pkg->ncoef      = nktot_res+1;

    if(pme_para_opt==2){
      pme_res_fft_pkg->myid       = myid_forc;
      pme_res_fft_pkg->myidp1     = myid_forc+1;
      pme_res_fft_pkg->num_proc   = np_forc;
      pme_res_fft_pkg->comm       = communicate->comm_forc;
    }else{
      pme_res_fft_pkg->myid       = 0;
      pme_res_fft_pkg->myidp1     = 1;
      pme_res_fft_pkg->num_proc   = 1;
      pme_res_fft_pkg->comm       = communicate->comm_faux;
    }/*endif*/

    create_para_fft_pkg3d(pme_res_fft_pkg,
                          ewald->kastr_res,ewald->kbstr_res,
                          ewald->kcstr_res,cp_dual_grid_opt_on);

    pme_res_fft_pkg->scale_opt = 0;
#ifdef HP_VECLIB
    if(pme_res_fft_pkg->igeneric_opt==0){pme_res_fft_pkg->scale_opt = -1;}
#endif

    if(np_forc>1 && pme_para_opt==2){
      create_pme_comm_full_g(n_interp_res,pme_res_fft_pkg);
    }/*endif*/

    pme_fft_pkg->nfft_size = MAX(pme_fft_pkg->nfft_size,
                                 pme_res_fft_pkg->nfft_size);

    if((pme_para_opt==1)&&(np_forc>1)){ 
      create_pme_comm_hybr(np_forc,myid_forc,communicate->comm_forc,
                           pme_res_fft_pkg);
    }/*endif*/

  }/*endif*/

/*=========================================================================*/
/* VI) Output */

  if(myid == 0){
    printf("\n");PRINT_LINE_DASH;
    printf("Finished setting up FFTs\n");
    PRINT_LINE_STAR;printf("\n");
  }/* endif myid */

/*-------------------------------------------------------------------------*/
     }/*end routine */
/*=========================================================================*/



/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
/* For cluster boundary conditions, get the long range correction term      */
/*==========================================================================*/

void get_coul_clus_corr(EWALD *ewald,
                        COMMUNICATE *communicate,
                        CELL *cell,int kmax_ewd,
                        int cp_on,double ecut_small_cp,
                        double ecut_small_dual,
                        int igeneric_opt,int pme_para_opt,
                        double *tot_memory,int cp_dual_grid_opt_on,
                        int nsplin_g)

/*=========================================================================*/
     {/*Begin Routine*/
/*=========================================================================*/
/*    Local Variables */
#include "../typ_defs/typ_mask.h"

  int nkf1,nkf2,nkf3;
  int *kmax_cp_clus,*kmaxv_clus;
  int nfft;
  int i,j,k,iii,ig,ir;
  int count;
  int index;  
  int skc_fft_ka_proc;
  int ekc_fft_ka_proc;
  int skb_fft_ka_proc;
  int ekb_fft_ka_proc;
  int ka,kb,kc;
  int ka_str,kb_str,kb_end;
  int dual_off = 0;

  double now_memory,phase;
  double ecut_tmp,ecut_small;
  double ecut_clus_good,edge;
  double *hmati,deth,rvol23,rtwoth;
  double gerf;
  double x,y,z;
  double alp,r,ralp,eee,tt,temp;
  double da,db,dc;
  double sa,sb,sc;
  double aka,akb,akc;
  double xk,yk,zk;
  double tpi,al2,g2,g;
  double ft_erf;
  double erf_rinv_lim;

  double clus_corr_short;
  double rcut,rp,rj0,arg,g_diff;
  double swit,swit1,gmax,dg_spl,dr,fpidr,gmax_true;
  double *g_spl,*f_spl;
  double *spl_corr_short0;
  double *spl_corr_short1;
  double *spl_corr_short2;
  double *spl_corr_short3;

  double bx,by,bz,gerf_true;
  double swx,swy,swz;
  double xcut,ycut,zcut;

  PARA_FFT_PKG3D *cp_para_fft_pkg3d_clus;
  double *zfft;
  double *zfft_tmp;
  double *dclus_corr_r;
  double *clus_corr_r;
  double *clus_corr_i;
  int ncoef_l_use,ncoef_l_proc;
  int icoef_off;

#ifdef DEBUG_CLUS_CORR
  double eext,ehart,ghfact;
  double eext_tmp,ehart_tmp;
  double ans_ext,alp_ch;
  double ans_ehart,rhog;
#endif

  double p  = 0.3614;
  double e1 = 0.2041422096422003, e2 = 0.1997535956961481;
  double e3 = 0.2213176596405576, e4 = 0.03360430734640255;
  double e5 = 0.4732592578721755, e6 =-0.509078520069735;
  double e7 = 0.6772631491947646, e8 =-0.369912979092217;
  double e9 = 0.06965131976970335;
  double eps = 1.0e-30;

  double *hmat;
  double *hmat_ewd;
  int *kastore         = ewald->kastr;
  int *kbstore         = ewald->kbstr;
  int *kcstore         = ewald->kcstr;
  int nktot            = ewald->nktot;
  int ncoef_l          = nktot + 1;

  double  alp_clus     = ewald->alp_clus; 
  double  alp_ewd      = ewald->alp_ewd; 
  double  ecut_clus    = ewald->ecut_clus;

  int myid_now;
  int myid = communicate->myid;
  int np_now;
  MPI_Comm comm_now;
/*
  int npts_rad  = ewald->npts_rad;
*/
  int npts_rad = 3001;
  double rheal;

  hmat         = cell->hmat;
  hmat_ewd     = cell->hmat_ewd;

  tpi = M_PI * 2.0;
  if(cp_on==1){
    myid_now  = communicate->myid_state;
    np_now    = communicate->np_states;
    comm_now  = communicate->comm_states;

    if(cp_dual_grid_opt_on == 2){
      ecut_small = ecut_small_dual;
    }else{
      ecut_small = ecut_small_cp;
    }/*endif cp_dual_grid_opt*/
  }else{
    deth       = getdeth(hmat);
    rtwoth     = -(2./3.);
    rvol23     = pow(deth,rtwoth);
    ecut_small = M_PI * .5 * M_PI * ((double) (kmax_ewd * kmax_ewd))*rvol23;
    if(pme_para_opt>0){
      myid_now  = communicate->myid_forc;
      np_now    = communicate->np_forc;
    }else{
      myid_now  = 0;
      np_now    = 1;
    }/*endif*/
    comm_now  = communicate->comm_forc;
  }/*endif*/

/*=========================================================================*/
/* I) Output to the screen */

  if(myid == 0){
    printf("\n");PRINT_LINE_STAR;
    printf("Setting up long range cluster correction\n");
    PRINT_LINE_DASH;printf("\n");
  }/* endif myid_now */

/*=========================================================================*/
/* II)  Malloc some temporary memory                                      */

  kmax_cp_clus  = (int *) cmalloc(3*sizeof(int))-1;
  kmaxv_clus    = (int *) cmalloc(3*sizeof(int))-1;
  hmati         = (double *)cmalloc(9*sizeof(double))-1;
  cp_para_fft_pkg3d_clus = (PARA_FFT_PKG3D *) cmalloc(sizeof(PARA_FFT_PKG3D));

/*=========================================================================*/
/* III) Get size of fine real-space grid for evaluation of Coulomb corr.   */

  gethinv(hmat_ewd,hmati,&deth,3);    
  edge      = pow(deth,(1.0/3.0));
  alp_clus /= edge;

  if(1.1*ecut_small > 0.5*ecut_clus){
    if(myid == 0){
      printf("$$$$$$$$$$$$$$$$$$$$_warning_$$$$$$$$$$$$$$$$$$$$\n");    
      printf("The cluster descreening cutoff is less than 1.1 times \n");
      printf("the true cutoff. Resetting to 1.1 times the true cutoff, %g.\n",
                         (2.0*ecut_small)*1.1);
      printf("$$$$$$$$$$$$$$$$$$$$_warning_$$$$$$$$$$$$$$$$$$$$\n\n");    
    }/*endif*/
    ecut_clus = (2.0*ecut_small)*1.1;
  }/*endif*/

  ecut_clus_good = 0.5*(2.0*M_PI*alp_clus)*(2.0*M_PI*alp_clus);
  if(ecut_clus_good > 0.5*ecut_clus){
    if(myid == 0){
      printf("$$$$$$$$$$$$$$$$$$$$_warning_$$$$$$$$$$$$$$$$$$$$\n");    
      printf("The cluster descreening cutoff less than the \n");
      printf("the recommended value, 0.5*(2*pi*alp/V^{1/3})^2. \n");
      printf("Resetting to  %g \n",ecut_clus_good*2.0);
      printf("$$$$$$$$$$$$$$$$$$$$_warning_$$$$$$$$$$$$$$$$$$$$\n\n");    
    }/*endif*/
    ecut_clus = (2.0*ecut_clus_good);
  }/*endif*/

  calc_cutoff(kmax_ewd,&ecut_tmp,&ecut_clus,1,
              kmax_cp_clus,kmaxv_clus,hmati,deth);
  
  nkf1 = 4*(kmax_cp_clus[1]+1);
  nkf2 = 4*(kmax_cp_clus[2]+1);
  nkf3 = 4*(kmax_cp_clus[3]+1);

  if(myid == 0){
   printf("Using an FFT grid size of %d by %d by %d\n",nkf1,nkf2,nkf3);
  }/* endif */

/*=========================================================================*/
/* IV) Create a temp parallel FFT package and strip out the useful info    */

  cp_para_fft_pkg3d_clus->nkf1         = nkf1;
  cp_para_fft_pkg3d_clus->nkf2         = nkf2;
  cp_para_fft_pkg3d_clus->nkf3         = nkf3;
      
  cp_para_fft_pkg3d_clus->nktot        = nktot;
  cp_para_fft_pkg3d_clus->ncoef        = ncoef_l;
   
  cp_para_fft_pkg3d_clus->myid         = myid_now;
  cp_para_fft_pkg3d_clus->myidp1       = myid_now+1;
  cp_para_fft_pkg3d_clus->num_proc     = np_now;
  cp_para_fft_pkg3d_clus->comm         = comm_now;
  cp_para_fft_pkg3d_clus->igeneric_opt = igeneric_opt;

  create_para_fft_pkg3d(cp_para_fft_pkg3d_clus,kastore,kbstore,kcstore,
                        dual_off);
  nfft            = cp_para_fft_pkg3d_clus->nfft_size;
  ncoef_l_proc    = cp_para_fft_pkg3d_clus->ncoef_proc;
  ncoef_l_use     = cp_para_fft_pkg3d_clus->ncoef_use;
  icoef_off       = cp_para_fft_pkg3d_clus->icoef_off;
  skc_fft_ka_proc = cp_para_fft_pkg3d_clus->skc_fft_ka_proc;
  ekc_fft_ka_proc = cp_para_fft_pkg3d_clus->ekc_fft_ka_proc;
  skb_fft_ka_proc = cp_para_fft_pkg3d_clus->skb_fft_ka_proc;
  ekb_fft_ka_proc = cp_para_fft_pkg3d_clus->ekb_fft_ka_proc;

/*=========================================================================*/
/* V) Malloc the memory                    */

  ewald->clus_corr_r  = (double *) cmalloc(ncoef_l_proc*sizeof(double))-1;
  ewald->dclus_corr_r = (double *) cmalloc(ncoef_l_proc*sizeof(double))-1;
  clus_corr_i         = (double *) cmalloc(ncoef_l_proc*sizeof(double))-1;

  clus_corr_r         = ewald->clus_corr_r;
  dclus_corr_r        = ewald->dclus_corr_r;

  g_spl             = (double *) cmalloc(nsplin_g*sizeof(double))-1;
  spl_corr_short3   = (double *) cmalloc(nsplin_g*sizeof(double))-1;
  spl_corr_short2   = (double *) cmalloc(nsplin_g*sizeof(double))-1;
  spl_corr_short1   = (double *) cmalloc(nsplin_g*sizeof(double))-1;
  spl_corr_short0   = (double *) cmalloc(nsplin_g*sizeof(double))-1;
  f_spl             = (double *) cmalloc(npts_rad*sizeof(double))-1;

  zfft                = (double *) cmalloc(nfft*sizeof(double))-1;
  zfft_tmp            = (double *) cmalloc(nfft*sizeof(double))-1;

/*=========================================================================*/
/* VI) Fill ZFFT with the function to be transformed on 3 space            */
/*     (1-S(r))*erf(alpha_clus*r)/r  and spline the short range part       */
/*     S(r)*erf(alpha_clus*r)/r using a radial grid                        */

   da = 1.0/((double) nkf1);
   db = 1.0/((double) nkf2);
   dc = 1.0/((double) nkf3);

   xcut = 0.5*hmat_ewd[1];
   ycut = 0.5*hmat_ewd[5];
   zcut = 0.5*hmat_ewd[9];
  
   rcut  = MIN3(xcut,ycut,zcut);  /* radial distance that fits in box    */
   rheal = rcut/10.0;             /* healing length is 0.1 of the radius */

 /*------------------------------------------------------------------------*/
 /* i) Long range part */

   index=0;
   erf_rinv_lim = alp_clus*2.0/sqrt(M_PI);
   for(kc=skc_fft_ka_proc;kc<=ekc_fft_ka_proc;kc++){
     kb_str = (kc==skc_fft_ka_proc ? skb_fft_ka_proc : 1);
     kb_end = (kc==ekc_fft_ka_proc ? ekb_fft_ka_proc : nkf2);
     for(kb=kb_str;kb<=kb_end;kb++){
       for(ka=1;ka<=nkf1;ka++){
         sa = da*((double)(ka-1)) - 0.5;
         sb = db*((double)(kb-1)) - 0.5;
         sc = dc*((double)(kc-1)) - 0.5;
         x  = sa*hmat[1]+sb*hmat[4]+sc*hmat[7];
         y  = sa*hmat[2]+sb*hmat[5]+sc*hmat[8];
         z  = sa*hmat[3]+sb*hmat[6]+sc*hmat[9];
         r  = sqrt(x*x+y*y+z*z);

         ralp    = r*(alp_clus);
         eee     = exp(-ralp*ralp);
         tt      = 1.0/(1.0+p*ralp);
         temp    = ((((((((e9*tt+e8)*tt+e7)*tt+e6)*tt+e5)*tt
                           +e4)*tt+e3)*tt+e2)*tt+e1)*tt*eee;
         gerf    = (1.0-temp);

         rp    = (r-rcut+rheal)/rheal;
         rp    = MIN(rp,1.0);
         rp    = MAX(rp,0.0);
         swit1 = -rp*rp*(2.0*rp-3.0);
#ifndef CHECK_CLUS
         zfft[(index+2*ka-1)] = (r > eps ? gerf*swit1/r : erf_rinv_lim*swit1);
         zfft[(index+2*ka)]   = 0.0;
#else
         zfft[(index+2*ka-1)] = swit1*exp(-ralp*ralp);
         zfft[(index+2*ka)]   = 0.0;
#endif
       }/*endfor:ka*/
       index += 2*nkf1;

     }/*endfor:kb*/
   }/*endfor:kc*/

 /*------------------------------------------------------------------------*/
 /* ii) Short range part via slow bessel transformation                    */

   gmax     = 0.0;
   for(i = 1; i<= nktot ; i++){
     aka = (double)kastore[i];
     akb = (double)kbstore[i];
     akc = (double)kcstore[i];
     xk = (aka*hmati[1] +akb*hmati[2] +akc*hmati[3])*tpi;
     yk = (aka*hmati[4] +akb*hmati[5] +akc*hmati[6])*tpi;
     zk = (aka*hmati[7] +akb*hmati[8] +akc*hmati[9])*tpi;
     g2 = xk*xk + yk*yk + zk*zk;
     gmax = MAX(gmax,g2);
   }/*endfor*/
   gmax     = 1.2*sqrt(gmax);

   dg_spl   = gmax/((double)nsplin_g-1);
   dr       = rcut/((double)npts_rad-1);
   fpidr    = 4.0*M_PI*dr;

   for(ig=1;ig <= nsplin_g;ig++){ 
     g_spl[ig] = dg_spl*((double) (ig-1));
   }/*endfor*/

   for(ir=1;ir<=npts_rad;ir++){
     r     = ((double)(ir-1))*dr;
     ralp  = r*(alp_clus);
     eee   = exp(-ralp*ralp);
     tt    = 1.0/(1.0+p*ralp);
     temp  = ((((((((e9*tt+e8)*tt+e7)*tt+e6)*tt+e5)*tt
                          +e4)*tt+e3)*tt+e2)*tt+e1)*tt*eee;
     gerf  = (1.0-temp);

     rp    = (r-rcut+rheal)/rheal;
     rp    = MIN(rp,1.0);
     rp    = MAX(rp,0.0);
     swit  = 1.0+rp*rp*(2.0*rp-3.0);
#ifndef CHECK_CLUS
     f_spl[ir] = gerf*swit;
#else
     f_spl[ir] = r*swit*exp(-ralp*ralp);
#endif
   }/*endfor*/

   for(ig=1;ig<=nsplin_g;ig++){spl_corr_short0[ig] = 0.0;}

   for(ir=1;ir<=npts_rad;ir++){
     r   = ((double)(ir-1))*dr;
     spl_corr_short0[1]  += fpidr*f_spl[ir]*r;
   }/*endfor*/
   for(ig=2;ig<=nsplin_g;ig++){
     for(ir=2;ir<=npts_rad;ir++){
       r   = ((double)(ir-1))*dr;
       arg = r*g_spl[ig];
       rj0 = (sin(arg)/arg)*r;
       spl_corr_short0[ig]  += fpidr*rj0*f_spl[ir];
     }/*endfor*/
   }/*endfor*/

   spline_fit(spl_corr_short0,spl_corr_short1,spl_corr_short2,
              spl_corr_short3,g_spl,nsplin_g);

/*=========================================================================*/
/* VII) FFT the long range function to g-space                             */

  para_fft_gen3d_bck_to_g(zfft,zfft_tmp,cp_para_fft_pkg3d_clus);
  sngl_upack_coef(clus_corr_r,clus_corr_i,zfft,cp_para_fft_pkg3d_clus);

/*=========================================================================*/
/* VIII) Subtract off the analytic transform of erf(alpha_clus*r)/r        */
/*       and add together the long and short range bits                    */

  tpi = 2.0*M_PI;
  al2 = alp_clus*alp_clus;
  for(count = 1; count<= ncoef_l_use ; count++){
     aka = (double)kastore[(count+icoef_off)];
     akb = (double)kbstore[(count+icoef_off)];
     akc = (double)kcstore[(count+icoef_off)];
     phase = cos(M_PI*(aka+akb+akc));
     xk = (aka*hmati[1] +akb*hmati[2] +akc*hmati[3])*tpi;
     yk = (aka*hmati[4] +akb*hmati[5] +akc*hmati[6])*tpi;
     zk = (aka*hmati[7] +akb*hmati[8] +akc*hmati[9])*tpi;
     g2 = xk*xk + yk*yk + zk*zk;
     g  = sqrt(g2);
    /* Spline of short range part */
     index       = (int)(g/dg_spl + 1.0);
     g_diff      = g - g_spl[index];
     clus_corr_short = ((spl_corr_short3[index]*g_diff
                        +spl_corr_short2[index])*g_diff
                        +spl_corr_short1[index])*g_diff
                        +spl_corr_short0[index];
    /* Infinite system reference */
     ft_erf = 2.0*tpi*exp(-0.25*g2/al2)/g2;
#ifdef CHECK_CLUS
     if(myid==0){
       temp = pow((M_PI/al2),1.5)*exp(-0.25*g2/al2);
       printf("g %g exact %g approx %g \n",g,temp,
                    (phase*clus_corr_r[count]*deth + clus_corr_short));
       if((count % 20)==0){scanf("%d",&ir);}
     }/*endif*/
#endif
    /* Construct cluster correction */
     clus_corr_r[count]  = phase*clus_corr_r[count]*deth - ft_erf
                         + clus_corr_short;
     dclus_corr_r[count] = 0.0;

  }/* endfor */

  /* g=0 term */
  if(ncoef_l_use< ncoef_l_proc){
     clus_corr_r[ncoef_l_proc] = clus_corr_r[ncoef_l_proc]*deth + M_PI/al2
                               + spl_corr_short0[1];
     dclus_corr_r[ncoef_l_proc] = 0.0;
  }/*endif*/

/*=========================================================================*/
#ifdef DEBUG_CLUS_CORR

  eext = 0.0;
  ehart = 0.0;
  alp_ch = 12.0/hmat[1];
  al2 = alp_ch*alp_ch;
  for(count=1;count<=ncoef_l_use;count++){
     aka = (double)kastore[(count+icoef_off)];
     akb = (double)kbstore[(count+icoef_off)];
     akc = (double)kcstore[(count+icoef_off)];
     xk = (aka*hmati[1] +akb*hmati[2] +akc*hmati[3])*tpi;
     yk = (aka*hmati[4] +akb*hmati[5] +akc*hmati[6])*tpi;
     zk = (aka*hmati[7] +akb*hmati[8] +akc*hmati[9])*tpi;
     g2 = xk*xk + yk*yk + zk*zk;
     rhog   = exp(-0.25*g2/al2);
     ghfact = (2.0*tpi/g2 + clus_corr_r[count])/deth;
     eext  += ghfact*rhog;
     ehart += (ghfact*rhog*rhog);
  }/*endfor*/
  eext     *= 2.0;
  ehart    *= 2.0;
  if(ncoef_l_use < ncoef_l_proc){
    eext      += (clus_corr_r[ncoef_l_proc]/deth);
    ehart     += (clus_corr_r[ncoef_l_proc]/deth); 
  }/*endif*/

  ans_ext   = 2.0*alp_ch/sqrt(M_PI);
  ans_ehart = alp_ch*sqrt(2.0/M_PI);

  if(np_now>1){
    Allreduce(&eext, &eext_tmp, 1,MPI_DOUBLE,MPI_SUM,0,comm_now);
    Allreduce(&ehart,&ehart_tmp,1,MPI_DOUBLE,MPI_SUM,0,comm_now);
    eext = eext_tmp;ehart=ehart_tmp;
  }/*endif*/
  if(myid==0){
     printf("eext  = %.12g  %.12g  %.12g\n",eext,ans_ext,
                               fabs((eext-ans_ext)/ans_ext));
     printf("ehart = %.12g  %.12g  %.12g\n",ehart,ans_ehart,
                              fabs((ehart-ans_ehart)/ans_ehart));
     printf("Enter an integer\n");
     scanf("%d",&iii);
  }/*endif*/
  Dbx_Barrier(comm_now);

#endif

/*=========================================================================*/
/* IX) Destroy the FFT package and free temporary memory                   */

  destroy_para_fft_pkg3d(cp_para_fft_pkg3d_clus); 

  cfree(&(kmax_cp_clus[1]));
  cfree(&(kmaxv_clus[1]));
  cfree(&(hmati[1]));
  cfree(&(zfft[1]));
  cfree(&(zfft_tmp[1]));
  cfree(&clus_corr_i[1]);
  cfree(&spl_corr_short3[1]);
  cfree(&spl_corr_short2[1]);
  cfree(&spl_corr_short1[1]);
  cfree(&spl_corr_short0[1]);
  cfree(&f_spl[1]);
  cfree(&g_spl[1]);
  cfree(cp_para_fft_pkg3d_clus);

/*=========================================================================*/
/* X) Output to screen */

  now_memory   = (ncoef_l_proc*sizeof(double))*1.0e-06;
  *tot_memory += now_memory;

  if(myid == 0){
    printf("Cluster corr. allocation: %g Mbytes; Total memory: %g Mbytes\n",
            now_memory,*tot_memory);
 
    printf("\n");PRINT_LINE_DASH;
    printf("Completed long range cluster correction");
    printf("\n");PRINT_LINE_STAR;
  }/* endif myid */

/*-------------------------------------------------------------------------*/ 
  }/* end routine */
/*=======================================================================*/

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
/* For 1D periodic boundary conditions, get the long range correction term  */
/*==========================================================================*/

void get_coul_1D_corr(EWALD *ewald,COMMUNICATE *communicate,
                      CELL *cell,int nsplin_g,int cp_on,int pme_on,
                      PARA_FFT_PKG3D *cp_para_fft_pkg_lg,
                      PARA_FFT_PKG3D *pme_fft_pkg,
                      double *tot_memory,int ifirst)

/*=========================================================================*/
     {/*Begin Routine*/
/*=========================================================================*/
/*    Local Variables */
#include "../typ_defs/typ_mask.h"

  int i,iw,ig,igr,index,kb_abs,kc_abs,ir;
  int iii,kndex;
  double w2;
  double hmati[10],deth;
  double aka,akb,akc,xk,yk,zk,g2,gs,alp,tpi;
  double fpi,falp2i,preg,gs2,wmin,wmax,wmax_y,wmax_z,func,dgs_now;
  double pregsw,sum_w_y,sum_w_z,one,two,three,phase,phase_y,phase_z,t_radpi;
  double gc,gc2,dgc,prey,prez,w,dy,dz,yalp,zalp,delta_gs,preg_ewd;  
  double edge,preg_clus;
  double sum_g0_y,sum_g0_z;
  double eee,eee_y,eee_z,gerf_y,gerf_z,tt,gerfc,screen_g0;
  double w_y,w_z;
  int ngo,irem,iend,istart;
  double now_memory;

  double p  = 0.3614;
  double e1 = 0.2041422096422003, e2 = 0.1997535956961481;
  double e3 = 0.2213176596405576, e4 = 0.03360430734640255;
  double e5 = 0.4732592578721755, e6 =-0.509078520069735;
  double e7 = 0.6772631491947646, e8 =-0.369912979092217;
  double e9 = 0.06965131976970335;

  int ngauss;
  double *anode,*weight,*w2_y,*w2_z,*prew_y,*prew_z,*gs_vec;

  int nfft_y,nfft_z,incl=1,num=1,iscale_opt,iopt=1,ier,incn_y,incn_z;
  int nwork_1=60000,nwork_2=60000,nfax=13;
  int *ifax_y,*ifax_z;
  double *work_1_y,*work_2_y,*work_1_z,*work_2_z,*zfft_y,*zfft_z;
  double **iy,**iz;

  int nktot        = ewald->nktot;
  int *ka          = ewald->kastr;
  int *kb          = ewald->kbstr;
  int *kc          = ewald->kcstr;
  int nkb_max      = ewald->nkb_max;
  int nkc_max      = ewald->nkc_max;
  double alp_clus  = ewald->alp_clus;
  double *hmat_ewd = cell->hmat_ewd;
  int myid         = communicate->myid;
  int myid_forc    = communicate->myid_forc;
  int np_forc      = communicate->np_forc;


  int igeneric_opt = cp_para_fft_pkg_lg->igeneric_opt;
  int ncoef_l_use;
  int ncoef_l_proc;
  int icoef_off;


  double *cvgs_2d0,*cvgs_2d1,*cvgs_2d2,*cvgs_2d3;/* after malloc */
  double *vgc_1d_y,*vgc_1d_z;        /* after malloc */
  double *clus_corr_r,*dclus_corr_r; /* after malloc */
  double gs_max;                     /* after assign */


/*=========================================================================*/
/* Output */

  if(myid == 0){
    printf("\n");PRINT_LINE_STAR;
    printf("Setting up 1D finite system coulomb correction \n");
    PRINT_LINE_DASH;printf("\n");
  }/* endif myid */

/*=========================================================================*/
/* Find the parallelization of g-space */

  switch(cp_on){
    case 0: 
       if(pme_on==1){
          ncoef_l_use  = pme_fft_pkg->ncoef_use;
          ncoef_l_proc = pme_fft_pkg->ncoef_proc;
          icoef_off    = pme_fft_pkg->icoef_off;
        }else{
          ngo      = nktot/np_forc;
          irem     = nktot % np_forc;
          if(myid_forc<=irem){
            istart = (ngo + 1)*myid_forc + 1;
          }else{
            istart = (ngo + 1)*irem + ngo*(myid_forc - irem) + 1;
          }/*endif*/
          if(myid_forc<irem){ngo++;}
          iend         = istart + ngo - 1;
          ncoef_l_proc = ngo;
          ncoef_l_use  = ngo;
          if(myid_forc+1==np_forc){ncoef_l_proc++;}
          icoef_off = istart-1;
        }/*endif*/
        break;
    case 1: 
        ncoef_l_use  = cp_para_fft_pkg_lg->ncoef_use;
        ncoef_l_proc = cp_para_fft_pkg_lg->ncoef_proc;
        icoef_off    = cp_para_fft_pkg_lg->icoef_off;
       break;
  }/*endif: switch*/

/*=========================================================================*/
/* I) Setup and Malloc the stuff  */

  gethinv(hmat_ewd,hmati,&deth,3);    
  edge      = pow(deth,(1.0/3.0));
  if(ifirst==1){alp_clus /= edge;}

  if(ifirst==1){
    ewald->nsplin_g     = nsplin_g; 
    ewald->cvgs_2d0     = (double *) cmalloc(nsplin_g*sizeof(double))-1;
    ewald->cvgs_2d1     = (double *) cmalloc(nsplin_g*sizeof(double))-1;
    ewald->cvgs_2d2     = (double *) cmalloc(nsplin_g*sizeof(double))-1;
    ewald->cvgs_2d3     = (double *) cmalloc(nsplin_g*sizeof(double))-1;
    ewald->vgc_1d_y     = (double *) cmalloc((nkb_max+1)*sizeof(double))-1;
    ewald->vgc_1d_z     = (double *) cmalloc((nkc_max+1)*sizeof(double))-1;
    ewald->clus_corr_r  = (double *) cmalloc(ncoef_l_proc*sizeof(double))-1;
    ewald->dclus_corr_r = (double *) cmalloc(ncoef_l_proc*sizeof(double))-1;
  }/*endif*/
  vgc_1d_y     = ewald->vgc_1d_y;
  vgc_1d_z     = ewald->vgc_1d_z;
  clus_corr_r  = ewald->clus_corr_r;
  dclus_corr_r = ewald->dclus_corr_r;

/*=========================================================================*/
/* II) Calculate the g_z dependent part of the descreening function        */

  if(ifirst==1){

 /*----------------------------------------------------*/
 /*  i) malloc and set up 1d FFTs                      */

    nfft_y   = MAX(2048,2*nkb_max);
    nfft_z   = MAX(2048,2*nkc_max);
    incn_y   = nfft_y;
    incn_z   = nfft_z;
    ifax_y   = (int *) cmalloc(13*sizeof(int))-1;
    ifax_z   = (int *) cmalloc(13*sizeof(int))-1;
    zfft_y   = (double *) cmalloc(2*nfft_y*sizeof(double))-1;
    zfft_z   = (double *) cmalloc(2*nfft_z*sizeof(double))-1;
    work_1_y = (double *) cmalloc(nwork_1*sizeof(double))-1;
    work_2_y = (double *) cmalloc(nwork_2*sizeof(double))-1;
    work_1_z = (double *) cmalloc(nwork_1*sizeof(double))-1;
    work_2_z = (double *) cmalloc(nwork_2*sizeof(double))-1;
    fft_gen1d_init(nfft_y,incl,num,incn_y,iopt,&ier,work_1_y,
                   nwork_1,work_2_y,nwork_2,
                   ifax_y,&iscale_opt,igeneric_opt);
    fft_gen1d_init(nfft_z,incl,num,incn_z,iopt,&ier,work_1_z,
                   nwork_1,work_2_z,nwork_2,
                   ifax_z,&iscale_opt,igeneric_opt);
 /*----------------------------------------------------*/
 /*  ii)  setup data and perform 1d FFTs               */

    dy       = hmat_ewd[5]/((double)(nfft_y));
    dz       = hmat_ewd[9]/((double)(nfft_z));
    prey     = dy*alp_clus/sqrt(M_PI);
    prez     = dz*alp_clus/sqrt(M_PI);
    for(i=1,ir=1;ir<=2*nfft_y;i++,ir+=2){
      yalp         = ( ((double)(i-1))*dy - 0.5*hmat_ewd[5] )*alp_clus;
      zfft_y[ir]     = prey*exp(-yalp*yalp);
      zfft_y[(ir+1)] = 0.0;
    }/*endfor*/

    fft_gen1d(zfft_y,nfft_y,incl,num,incn_y,iopt,&ier,
              work_1_y,nwork_1,work_2_y,nwork_2,
              ifax_y,igeneric_opt);

    for(i=1,ir=1;ir<=2*nfft_z;i++,ir+=2){
      zalp         = ( ((double)(i-1))*dz - 0.5*hmat_ewd[9] )*alp_clus;
      zfft_z[ir]     = prez*exp(-zalp*zalp);
      zfft_z[(ir+1)] = 0.0;
    }/*endfor*/

    fft_gen1d(zfft_z,nfft_z,incl,num,incn_z,iopt,&ier,
              work_1_z,nwork_1,work_2_z,nwork_2,
              ifax_z,igeneric_opt);
 /*----------------------------------------------------*/
 /* iii) modify the results appropriately              */

    for(ig=1,igr=1;ig<=(nkb_max+1);ig++,igr+=2){
      phase      = ( ((ig-1) % 2 == 0) ? 1.0 : -1.0);
      vgc_1d_y[ig] = zfft_y[igr]*phase;
    }/*endfor*/

    for(ig=1,igr=1;ig<=(nkc_max+1);ig++,igr+=2){
      phase      = ( ((ig-1) % 2 == 0) ? 1.0 : -1.0);
      vgc_1d_z[ig] = zfft_z[igr]*phase;
    }/*endfor*/

  }/*endif ifirst*/

/*=========================================================================*/
/* III) Spline the g_x dependent part of the descreening function           */

  if(ifirst==1){

 /*----------------------------------------------------*/
 /*  ii)  Set up the integration                       */

     ngauss  = 128;
     anode   = (double *) cmalloc(ngauss*sizeof(double))-1;
     weight  = (double *) cmalloc(ngauss*sizeof(double))-1;
     w2_y    = (double *) cmalloc(ngauss*sizeof(double))-1;
     w2_z    = (double *) cmalloc(ngauss*sizeof(double))-1;
     prew_y  = (double *) cmalloc(ngauss*sizeof(double))-1;
     prew_z  = (double *) cmalloc(ngauss*sizeof(double))-1;
     iy      = cmall_mat(1,ngauss,1,(nkb_max+1));
     iz      = cmall_mat(1,ngauss,1,(nkc_max+1));
#include "../proto_defs/weights_nodes_128.h"
     prey     = dy*4.0/M_PI;
     prez     = dz*4.0/M_PI;
#ifdef SAVE
     wmax_y    = MIN(0.5*alp_clus*hmat_ewd[5],12.0);
     wmax_z    = MIN(0.5*alp_clus*hmat_ewd[9],12.0);
#endif
     wmax_y    = 0.5*alp_clus*hmat_ewd[5];
     wmax_z    = 0.5*alp_clus*hmat_ewd[9];
     for(iw=1;iw<=ngauss;iw++){     
       w_y        = 0.5*wmax_y*(anode[iw]+1.0);
       w_z        = 0.5*wmax_z*(anode[iw]+1.0);
       w2_y[iw]   = w_y*w_y;
       w2_z[iw]   = w_z*w_z;
       prew_y[iw] = (0.5*wmax_y*weight[iw])*exp(-w2_y[iw]);
       prew_z[iw] = (0.5*wmax_z*weight[iw])*exp(-w2_z[iw]);

/*-------------------------*/
/* z internal integral */

       for(i=1,ir=1;ir<=2*nfft_z;i++,ir+=2){
         zalp  = ( ((double)(i-1))*dz - 0.5*hmat_ewd[9] )*2.0*w_y/hmat_ewd[5];
         zfft_z[ir]     = prez*exp(-zalp*zalp);
         zfft_z[(ir+1)] = 0.0;
       }/*endfor*/

      fft_gen1d(zfft_z,nfft_z,incl,num,incn_z,iopt,&ier,
                work_1_z,nwork_1,work_2_z,nwork_2,ifax_z,igeneric_opt);

      for(ig=1,igr=1;ig<=(nkc_max+1);ig++,igr+=2){
        phase      = ( ((ig-1) % 2 == 0) ? 1.0 : -1.0);
        iz[iw][ig] = zfft_z[igr]*phase;
      }/*endfor*/


/*-------------------------*/
/* y internal integral */

       for(i=1,ir=1;ir<=2*nfft_y;i++,ir+=2){
         yalp   = ( ((double)(i-1))*dy - 0.5*hmat_ewd[5] )*2.0*w_z/hmat_ewd[9];
         zfft_y[ir]     = prey*exp(-yalp*yalp);
         zfft_y[(ir+1)] = 0.0;
       }/*endfor*/

      fft_gen1d(zfft_y,nfft_y,incl,num,incn_y,iopt,&ier,
                work_1_y,nwork_1,work_2_y,nwork_2,ifax_y,igeneric_opt);

      for(ig=1,igr=1;ig<=(nkb_max+1);ig++,igr+=2){
        phase      = ( ((ig-1) % 2 == 0) ? 1.0 : -1.0);
        iy[iw][ig] = zfft_y[igr]*phase;
      }/*endfor*/

     }/* endfor */


  }/*endif*/
/*=========================================================================*/
/* III) Construct the descreening function                                 */

  tpi      = 2.0*M_PI;
  fpi      = 4.0*M_PI;
  falp2i   = 1.0/(4.0*alp_clus*alp_clus);

  for(i=1;i<=ncoef_l_use ;i++){

     kndex       = i+icoef_off;
     aka         = (double)(ka[kndex]);
     akb         = (double)(kb[kndex]);
     akc         = (double)(kc[kndex]);
     xk          = aka*hmati[1]*tpi;
     yk          = akb*hmati[5]*tpi;
     zk          = akc*hmati[9]*tpi;
     g2          = xk*xk+yk*yk+zk*zk;
     gs2         = xk*xk;

 /*----------------------------------------------------*/
 /*  iii)  Evaluate the Gauss quadrature integrals     */

      kb_abs      = ( (kb[kndex] >= 0) ? kb[kndex] : -kb[kndex]);
      kc_abs      = ( (kc[kndex] >= 0) ? kc[kndex] : -kc[kndex]);
      pregsw     = gs2*hmat_ewd[5]*hmat_ewd[5]/16.0;
      sum_w_y      = 0.0;
      for(iw=1;iw<=ngauss;iw++){
        w_y    = 0.5*wmax_y*(anode[iw]+1.0);
        func   = w_y*exp(-pregsw/w2_y[iw]);
        sum_w_y += func*prew_y[iw]*iz[iw][(kc_abs+1)];
      }/*endfor*/
       pregsw     = gs2*hmat_ewd[9]*hmat_ewd[9]/16.0;
      sum_w_z      = 0.0;
      for(iw=1;iw<=ngauss;iw++){
        w_z    = 0.5*wmax_z*(anode[iw]+1.0);
        func   = w_z*exp(-pregsw/w2_z[iw]);
        sum_w_z += func*prew_z[iw]*iy[iw][(kb_abs+1)];
      }/*endfor*/
 /*----------------------------------------------------*/

    preg        = fpi/g2;    
    phase_y     = ( (kb[kndex] % 2 == 0) ? -1 : 1);
    phase_z     = ( (kc[kndex] % 2 == 0) ? -1 : 1);
    one         = exp(-gs2*falp2i)*vgc_1d_y[(kb_abs+1)]*vgc_1d_z[(kc_abs+1)] 
                - exp(-g2*falp2i);
    two         = phase_y*sum_w_y/hmat_ewd[5];
    three       = phase_z*sum_w_z/hmat_ewd[9];

    clus_corr_r[i]  = preg*(one + two + three);
    dclus_corr_r[i] = 0.0;


#ifdef DEBUG_1DCORR
    if(ka[kndex]==8 && kb[kndex]==1 && kc[kndex]==1){
      for(iw=1;iw<=ngauss;iw++){
       printf("%.12g %.12g\n",
              iz[iw][(kc_abs+1)]*(dz/prez),iy[iw][(kb_abs+1)]*(dy/prey));
      }
       printf("sums y %.12g z %.12g\n",sum_w_y,sum_w_z);
       printf("g2 %.12g\n",g2);
       printf("preg %.12g Ly %.12g Ly %.12g\n",preg,hmat_ewd[5],hmat_ewd[9]);
       printf("one %.12g two %.12g three %.12g\n",
               preg*one,preg*two,preg*three);
       printf("alp %.12g code answer %.12g\n",alp_clus,clus_corr_r[i]);
       scanf("%d",&iii);
    }
#endif



  }/*endfor*/

 /*----------------------------------------------------*/
 /*  iv)  Evaluate the g=0 term                        */

     wmin     = 0.5*alp_clus*hmat_ewd[5];
     wmax     = MAX_INT;
      sum_g0_y   = 0.0;   
     if(wmax > wmin){
       for(iw=1;iw<=ngauss;iw++){     
          w       = 0.5*( (wmax-wmin)*anode[iw]+wmax+wmin );
         w2      = w*w;
         eee_y   = exp(-w2);
         yalp    = w*hmat_ewd[9]/hmat_ewd[5];
         eee     = exp(-yalp*yalp);
        tt      = 1.0/(1.0+p*wmin);
        gerf_y  = 1.0-((((((((e9*tt+e8)*tt+e7)*tt+e6)*tt+e5)*tt
                               +e4)*tt+e3)*tt+e2)*tt+e1)*tt*eee;
        sum_g0_y += (weight[iw]*eee_y*gerf_y/w2);
       }/*endfor*/
     }/* endif */
       sum_g0_y *= alp_clus*alp_clus*hmat_ewd[5]*hmat_ewd[5]*
                   ( (wmax-wmin)/sqrt(M_PI) );
     yalp    = wmin;
     eee     = exp(-yalp*yalp);
     tt      = 1.0/(1.0+p*wmin);
     gerf_y  = 1.0-((((((((e9*tt+e8)*tt+e7)*tt+e6)*tt+e5)*tt
                              +e4)*tt+e3)*tt+e2)*tt+e1)*tt*eee;
   
     wmin     = 0.5*alp_clus*hmat_ewd[9];
     sum_g0_z   = 0.0;   
     if(wmax > wmin){
       for(iw=1;iw<=ngauss;iw++){     
         w       = 0.5*( (wmax-wmin)*anode[iw]+wmax+wmin );
         w2      = w*w;
         eee_z   = exp(-w2);
         zalp    = w*hmat_ewd[5]/hmat_ewd[9];
         eee     = exp(-zalp*zalp);
         tt      = 1.0/(1.0+p*wmin);
         gerf_z  = 1.0-((((((((e9*tt+e8)*tt+e7)*tt+e6)*tt+e5)*tt
                               +e4)*tt+e3)*tt+e2)*tt+e1)*tt*eee;
   
         sum_g0_z += (weight[iw]*eee_z*gerf_z/w2);
       }/*endfor*/
     }/* endif */
     sum_g0_z *= alp_clus*alp_clus*hmat_ewd[9]*hmat_ewd[9]*
                 ( (wmax-wmin)/sqrt(M_PI) );
     zalp    = wmin;
     eee     = exp(-zalp*zalp);
     tt      = 1.0/(1.0+p*wmin);
     gerf_z  = 1.0-((((((((e9*tt+e8)*tt+e7)*tt+e6)*tt+e5)*tt
                              +e4)*tt+e3)*tt+e2)*tt+e1)*tt*eee;
     screen_g0 = -M_PI/(alp_clus*alp_clus)*(gerf_y*gerf_z 
               +  0.25*(sum_g0_y + sum_g0_z));
     screen_g0 += M_PI/(alp_clus*alp_clus);

 /*----------------------------------------------------*/
  if(ncoef_l_use<ncoef_l_proc){
    clus_corr_r[ncoef_l_proc]   = screen_g0;
    dclus_corr_r[ncoef_l_proc]  = 0.0;
  }/*endif*/  

#ifdef DEBUG_1DCORR
  printf("screen_g0 %.12g\n",screen_g0);
#endif

/*=========================================================================*/
/* IV) Free memory */

     cfree(&zfft_y[1]);
     cfree(&zfft_z[1]);
     cfree(&work_1_y[1]);
     cfree(&work_2_y[1]);
     cfree(&work_1_z[1]);
     cfree(&work_2_z[1]);
     cfree(&ifax_y[1]);
     cfree(&ifax_z[1]);

     cfree(&anode[1]);
     cfree(&weight[1]);
     cfree(&w2_y[1]);
     cfree(&w2_z[1]);
     cfree(&prew_y[1]);
     cfree(&prew_z[1]);

     cfree_mat(iy,1,ngauss,1,(nkb_max+1));
     cfree_mat(iz,1,ngauss,1,(nkc_max+1));


/*=========================================================================*/
/* IV) Output */

  now_memory   = ((2*ncoef_l_proc+4*nsplin_g + nkc_max+1)*sizeof(double))
                 *1.0e-06;
  *tot_memory += now_memory;

  if(myid==0){
    printf("1D ewald corr. allocation: %g Mbytes; Total memory: %g Mbytes\n",
            now_memory,*tot_memory);
 
    printf("\n");PRINT_LINE_DASH;
    printf("Completed 1D finite system coulomb correction\n");
    printf("\n");PRINT_LINE_STAR;
  }/*endif*/


/*-------------------------------------------------------------------------*/ 
    }/* end routine */
/*=========================================================================*/


/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
/* For 2D periodic boundary conditions, get the long range correction term  */
/*==========================================================================*/

void get_coul_2D_corr(EWALD *ewald,COMMUNICATE *communicate,
                      CELL *cell,int nsplin_g,int cp_on,int pme_on,
                      PARA_FFT_PKG3D *cp_para_fft_pkg_lg,
                      PARA_FFT_PKG3D *pme_fft_pkg,
                      double *tot_memory,int ifirst)

/*=========================================================================*/
     {/*Begin Routine*/
/*=========================================================================*/
/*    Local Variables */
#include "../typ_defs/typ_mask.h"

  int i,iw,ig,igr,index,kc_abs,ir;
  int iii,kndex;
  double hmati[10],deth;
  double aka,akb,akc,xk,yk,zk,g2,gs,alp,tpi;
  double fpi,falp2i,preg,gs2,wmin,wmax,func,dgs_now;
  double pregsw,sum_w,one,two,phase,t_radpi;
  double gc,gc2,dgc,prez,w,dz,zalp,delta_gs,preg_ewd;  
  double edge,preg_clus;
  double sum_g0;
  double eee,tt,screen_g0;
  double alp_rpi,eglm,eglp;
  double arg,alp2,ralp,cerf,nerf;
  double t1,t2,t3;
  int ngo,irem,iend,istart;
  double now_memory;

  int nfft,incl=1,num=1,iscale_opt,iopt=1,ier,incn;
  int nwork_1=60000,nwork_2=60000,nfax=13;
  int *ifax;
  double *work_1,*work_2,*zfft;

  int nktot        = ewald->nktot;
  int *ka          = ewald->kastr;
  int *kb          = ewald->kbstr;
  int *kc          = ewald->kcstr;
  int nkc_max      = ewald->nkc_max;
  double alp_clus  = ewald->alp_clus;
  double *hmat_ewd = cell->hmat_ewd;
  int myid         = communicate->myid;
  int myid_forc    = communicate->myid_forc;
  int np_forc      = communicate->np_forc;


  int igeneric_opt = cp_para_fft_pkg_lg->igeneric_opt;
  int ncoef_l_use;
  int ncoef_l_proc;
  int icoef_off;


  double *cvgs_2d0,*cvgs_2d1,*cvgs_2d2,*cvgs_2d3,*vgc_2d; /* after malloc */
  double *clus_corr_r,*dclus_corr_r;                      /* after malloc */
  double gs_max;                                          /* after assign */

/*=========================================================================*/
/* Output */

  if(myid == 0){
    printf("\n");PRINT_LINE_STAR;
    printf("Setting up 2D finite system coulomb correction \n");
    PRINT_LINE_DASH;printf("\n");
  }/* endif myid */

/*=========================================================================*/
/* Find the parallelization of g-space */

  switch(cp_on){
    case 0: 
       if(pme_on==1){
          ncoef_l_use  = pme_fft_pkg->ncoef_use;
          ncoef_l_proc = pme_fft_pkg->ncoef_proc;
          icoef_off    = pme_fft_pkg->icoef_off;
        }else{
          ngo      = nktot/np_forc;
          irem     = nktot % np_forc;
          if(myid_forc<=irem){
            istart = (ngo + 1)*myid_forc + 1;
          }else{
            istart = (ngo + 1)*irem + ngo*(myid_forc - irem) + 1;
          }/*endif*/
          if(myid_forc<irem){ngo++;}
          iend         = istart + ngo - 1;
          ncoef_l_proc = ngo;
          ncoef_l_use  = ngo;
          if(myid_forc+1==np_forc){ncoef_l_proc++;}
          icoef_off = istart-1;
        }/*endif*/
        break;
    case 1: 
        ncoef_l_use  = cp_para_fft_pkg_lg->ncoef_use;
        ncoef_l_proc = cp_para_fft_pkg_lg->ncoef_proc;
        icoef_off    = cp_para_fft_pkg_lg->icoef_off;
       break;
  }/*endif: switch*/

/*=========================================================================*/
/* I) Setup and Malloc the stuff  */

  gethinv(hmat_ewd,hmati,&deth,3);    
  edge      = pow(deth,(1.0/3.0));
  if(ifirst==1){alp_clus /= edge;ewald->alp_clus = alp_clus;}
  alp2 = alp_clus*alp_clus;

  if(ifirst==1){
    ewald->nsplin_g     = nsplin_g; 
    ewald->cvgs_2d0     = (double *) cmalloc(nsplin_g*sizeof(double))-1;
    ewald->cvgs_2d1     = (double *) cmalloc(nsplin_g*sizeof(double))-1;
    ewald->cvgs_2d2     = (double *) cmalloc(nsplin_g*sizeof(double))-1;
    ewald->cvgs_2d3     = (double *) cmalloc(nsplin_g*sizeof(double))-1;
    ewald->vgc_2d       = (double *) cmalloc((nkc_max+1)*sizeof(double))-1;
    ewald->clus_corr_r  = (double *) cmalloc(ncoef_l_proc*sizeof(double))-1;
    ewald->dclus_corr_r = (double *) cmalloc(ncoef_l_proc*sizeof(double))-1;
  }/*endif*/
  cvgs_2d0     = ewald->cvgs_2d0;     /* must be outside if statement */
  cvgs_2d1     = ewald->cvgs_2d1;     /* after the malloc */
  cvgs_2d2     = ewald->cvgs_2d2;
  cvgs_2d3     = ewald->cvgs_2d3;
  vgc_2d       = ewald->vgc_2d;
  clus_corr_r  = ewald->clus_corr_r;
  dclus_corr_r = ewald->dclus_corr_r;

/*=========================================================================*/
/* II) Calculate the g_z dependent part of the descreening function        */

  if(ifirst==1){

 /*----------------------------------------------------*/
 /*  i) malloc and set up 1d FFT                       */

    nfft   = MAX(8192,2*nkc_max);
    incn   = nfft;
    ifax   = (int *) cmalloc(13*sizeof(int))-1;
    zfft   = (double *) cmalloc(2*nfft*sizeof(double))-1;
    work_1 = (double *) cmalloc(nwork_1*sizeof(double))-1;
    work_2 = (double *) cmalloc(nwork_2*sizeof(double))-1;
    fft_gen1d_init(nfft,incl,num,incn,iopt,&ier,work_1,nwork_1,work_2,nwork_2,
                   ifax,&iscale_opt,igeneric_opt);
 /*-----------------------------------------------------------*/
 /*  ii)  setup data and perform 1d FFT to get complex erfc   */

    dz       = hmat_ewd[9]/((double)(nfft));
    prez     = dz*alp_clus/sqrt(M_PI);
    for(i=1,ir=1;ir<=2*nfft;i++,ir+=2){
      zalp         = ( ((double)(i-1))*dz - 0.5*hmat_ewd[9] )*alp_clus;
      zfft[ir]     = prez*exp(-zalp*zalp);
      zfft[(ir+1)] = 0.0;
    }/*endfor*/

    fft_gen1d(zfft,nfft,incl,num,incn,iopt,&ier,work_1,nwork_1,work_2,nwork_2,
              ifax,igeneric_opt);

    for(ig=1,igr=1;ig<=(nkc_max+1);ig++,igr+=2){
      vgc_2d[ig] = ((ig % 2) == 0 ? -zfft[igr]:zfft[igr]);
    }

 /*----------------------------------------------------*/
 /* iv) free the memory                                */

    cfree(&zfft[1]);
    cfree(&work_1[1]);
    cfree(&work_2[1]);
    cfree(&ifax[1]);

  }/*endif*/

/*=========================================================================*/
/* III) Calculate the remaining part of screening function                 */

  if(ifirst==1){

 /*----------------------------------------------------*/
 /* i)    Find g_s^{max} and add some saftey           */

    tpi    = 2.0*M_PI;
    gs_max = 0.0;
    for(i=1;i<=nktot;i++){
      aka    = (double)(ka[i]);
      akb    = (double)(kb[i]);
      xk     = (aka*hmati[1]+akb*hmati[2])*tpi;
      yk     = (aka*hmati[4]+akb*hmati[5])*tpi;
      gs     = sqrt(xk*xk+yk*yk);
      gs_max =  MAX(gs,gs_max);
    }/*endfor*/
    gs_max       *=1.2;
    ewald->gs_max = gs_max; 


 /*----------------------------------------------------*/
 /*  iv)  Evaluate the g=0 term                        */

    arg = 0.5*hmat_ewd[9]*alp_clus;
    ralp   = 0.5*alp_clus*hmat_ewd[9];
    cerf   = gerfc(ralp);
#define MARK
#ifdef MARK
    screen_g0 = M_PI*cerf*(1.0/alp2 + 0.5*hmat_ewd[9]*hmat_ewd[9])
              - (sqrt(M_PI)*hmat_ewd[9]/alp_clus)*exp(-arg*arg);
#endif
#ifdef GLENN_1
    nerf = gerf(ralp);
    screen_g0 = -M_PI*nerf/alp2 - 0.5*nerf*hmat_ewd[9]*hmat_ewd[9]
              - (sqrt(M_PI)*hmat_ewd[9]/alp_clus)*exp(-arg*arg);
#endif
#ifdef GLENN_2
    nerf = gerf(ralp);
    screen_g0 = M_PI*cerf/alp2 - 0.5*nerf*hmat_ewd[9]*hmat_ewd[9]
              - (sqrt(M_PI)*hmat_ewd[9]/alp_clus)*exp(-arg*arg);
#endif

  }/*endif*/
  gs_max   = ewald->gs_max;              /* Outside of if statement*/
  delta_gs = gs_max/((double) nsplin_g); /* for ifirst != 1 */

/*=========================================================================*/
/* III) Construct the descreening function                                 */

  tpi      = 2.0*M_PI;
  fpi      = 4.0*M_PI;
  falp2i   = 1.0/(4.0*alp_clus*alp_clus);
  alp_rpi  = alp_clus/sqrt(M_PI);

  for(i=1;i<=ncoef_l_use ;i++){

    kndex       = i+icoef_off;
    aka         = (double)(ka[kndex]);
    akb         = (double)(kb[kndex]);
    akc         = (double)(kc[kndex]);
    xk          = (aka*hmati[1]+akb*hmati[2])*tpi;
    yk          = (aka*hmati[4]+akb*hmati[5])*tpi;
    zk          =  akc*hmati[9]*tpi;
    g2          = xk*xk+yk*yk+zk*zk;
    gs          = sqrt(xk*xk+yk*yk);
    gs2         = gs*gs;
    arg         = 0.5*gs*hmat_ewd[9];
    eglm        = exp(-arg);
    if(arg > 100.0) arg = 100.0; 
    eglp        = exp(arg); 

    preg        = fpi/g2;    
    kc_abs      = ( (kc[kndex] >= 0) ? kc[kndex] : -kc[kndex]);
    phase       = ( (kc[kndex] % 2 == 0) ? -1.0 : 1.0);
    one         = exp(-g2*falp2i) - exp(-gs2*falp2i)*vgc_2d[(kc_abs+1)];

    ralp   = 0.5*(alp2*hmat_ewd[9]-gs)/alp_clus;
    if(ralp > 0.0){
     cerf  = gerfc(ralp);
    } else {
     cerf = 1.0 + gerf(ralp);
    }
    t2 = 0.5*eglm*cerf;

    ralp   = 0.5*(alp2*hmat_ewd[9]+gs)/alp_clus;
    cerf  = gerfc(ralp);
    t3 = 0.5*eglp*cerf;

    two         = phase*(eglm - t2 - t3) - one;

    clus_corr_r[i]  = preg*two;
    dclus_corr_r[i] = 0.0;

  }/*endfor*/

  if(ncoef_l_use<ncoef_l_proc){
    clus_corr_r[ncoef_l_proc]   = screen_g0;
    dclus_corr_r[ncoef_l_proc]  = 0.0;
  }/*endif*/  

/*=========================================================================*/
/* IV) Output */

  now_memory   = ((2*ncoef_l_proc+4*nsplin_g + nkc_max+1)*sizeof(double))
                 *1.0e-06;
  *tot_memory += now_memory;

  if(myid==0){
    printf("2D ewald corr. allocation: %g Mbytes; Total memory: %g Mbytes\n",
            now_memory,*tot_memory);
 
    printf("\n");PRINT_LINE_DASH;
    printf("Completed 2D finite system coulomb correction\n");
    printf("\n");PRINT_LINE_STAR;
  }/*endif*/

/*-------------------------------------------------------------------------*/ 
    }/* end routine */
/*=========================================================================*/



/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
/* Initalize the pme routine get the weights                                */
/*==========================================================================*/

void init_cp_dual_pme_bw(CPSCR_DUAL_PME *cpscr_dual_pme,EWALD *ewald,
                         CPEWALD *cpewald,int n_interp_pme_dual,
                         PARA_FFT_PKG3D *para_fft_pkg_lg)

/*=========================================================================*/
     {/*Begin Routine*/
/*=========================================================================*/
/*    Local Variables */
#include "../typ_defs/typ_mask.h"

/* local pointers */
 int nkf1,nkf2,nkf3;

 int nktot       = ewald->nktot;
 int *kastore    = ewald->kastr;
 int *kbstore    = ewald->kbstr;
 int *kcstore    = ewald->kcstr;
 int *kmax_cp    = cpewald->kmax_cp;
 int *kmax_cp_dens_cp_box = cpewald->kmax_cp_dens_cp_box;
 int box_rat     = cpewald->box_rat;
 double *bw_r    = cpscr_dual_pme->bw_r;
 double *bw_i    = cpscr_dual_pme->bw_i;

 double *aj      = cpscr_dual_pme->aj;
 double *rn      = cpscr_dual_pme->rn;
 double *rn1     = cpscr_dual_pme->rn1;
 double *bw_mag;
 int pme_b_opt = 0;

/* local pointers */
  int ncoef_proc   = para_fft_pkg_lg->ncoef_proc;
  int ncoef_use    = para_fft_pkg_lg->ncoef_use;
  int icoef_off    = para_fft_pkg_lg->icoef_off;

#ifdef PME
  nkf1 = 4*(kmax_cp[1] + 1);
  nkf2 = 4*(kmax_cp[2] + 1);
  nkf3 = 4*(kmax_cp[3] + 1);
#endif

#ifdef ORIG
  nkf1 = 4*box_rat*(kmax_cp_dens_cp_box[1] + 1);
  nkf2 = 4*box_rat*(kmax_cp_dens_cp_box[2] + 1);
  nkf3 = 4*box_rat*(kmax_cp_dens_cp_box[3] + 1);
#endif

 bw_mag = (double *) cmalloc(sizeof(double))-1;

 set_pme_wght(nktot,kastore,kbstore,kcstore,
              nkf1,nkf2,nkf3,ncoef_proc,ncoef_use,
              icoef_off,pme_b_opt,
              bw_r,bw_i,bw_mag,n_interp_pme_dual,
              aj,rn,rn1);


/*-------------------------------------------------------------------------*/ 
    }/* end routine */
/*=========================================================================*/

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
/* Initalize the pme routine get the weights                                */
/*==========================================================================*/

void init_cp_atom_pme_bw(CPSCR_ATOM_PME *cpscr_atom_pme,EWALD *ewald,
                         CPEWALD *cpewald,PARA_FFT_PKG3D *para_fft_pkg_lg)

/*=========================================================================*/
     {/*Begin Routine*/
/*=========================================================================*/
/*    Local Variables */
#include "../typ_defs/typ_mask.h"

/* local pointers */
 int nkf1,nkf2,nkf3;

 int nktot       = ewald->nktot;
 int *kastore    = ewald->kastr;
 int *kbstore    = ewald->kbstr;
 int *kcstore    = ewald->kcstr;
 int *kmax_cp    = cpewald->kmax_cp;
 int *kmax_cp_dens_cp_box = cpewald->kmax_cp_dens_cp_box;
 int box_rat     = cpewald->box_rat;
 int n_interp_pme = cpscr_atom_pme->n_interp;
 double *bw_r    = cpscr_atom_pme->bw_r;
 double *bw_i    = cpscr_atom_pme->bw_i;
 double *bweight_tot = cpscr_atom_pme->bweight_tot;

 double *aj      = cpscr_atom_pme->aj;
 double *rn      = cpscr_atom_pme->rn;
 double *rn1     = cpscr_atom_pme->rn1;
 double *bw_mag;
 int pme_b_opt = 2;

/* local pointers */
  int ncoef_proc   = para_fft_pkg_lg->ncoef_proc;
  int ncoef_use    = para_fft_pkg_lg->ncoef_use;
  int icoef_off    = para_fft_pkg_lg->icoef_off;

#ifdef PME
  nkf1 = 4*(kmax_cp[1] + 1);
  nkf2 = 4*(kmax_cp[2] + 1);
  nkf3 = 4*(kmax_cp[3] + 1);
#endif

 set_pme_wght(nktot,kastore,kbstore,kcstore,
              nkf1,nkf2,nkf3,ncoef_proc,ncoef_use,
              icoef_off,pme_b_opt,
              bw_r,bw_i,bweight_tot,n_interp_pme,
              aj,rn,rn1);

/*-------------------------------------------------------------------------*/ 
    }/* end routine */
/*=========================================================================*/
