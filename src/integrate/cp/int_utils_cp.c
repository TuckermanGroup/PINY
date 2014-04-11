/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
/*                                                                          */
/*                         PI_MD:                                           */
/*             The future of simulation technology                          */
/*             ------------------------------------                         */
/*                     Module: int_utilities                                */
/*                                                                          */
/* This subprogram provides some integrator utility routines                */
/*                                                                          */
/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

#include "standard_include.h"
#include "../typ_defs/typedefs_gen.h"
#include "../typ_defs/typedefs_cp.h"
#include "../typ_defs/typedefs_class.h"
#include "../typ_defs/typedefs_bnd.h"
#include "../proto_defs/proto_integrate_cp_local.h"
#include "../proto_defs/proto_communicate_wrappers.h"

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
 
void get_cpke(CPCOEFFS_INFO *cpcoeffs_info,CPCOEFFS_POS *cpcoeffs_pos,
              STAT_AVG *stat_avg, int cp_lsda,int np_states)

/*========================================================================*/
{/*begin routine*/
/*========================================================================*/
/*             Local variable declarations                                */

#include "../typ_defs/typ_mask.h"

    int is,icoef,i,j,iii;

    double kinet_cp_up,kinet_cp_up_tmp,kinet_cp_dn,kinet_cp_dn_tmp;
    double kinet_old,istart,iend;
    double *cpcoeffs_vcre_up  = cpcoeffs_pos->vcre_up;
    double *cpcoeffs_vcim_up  = cpcoeffs_pos->vcim_up;
    double *cpcoeffs_vcre_dn  = cpcoeffs_pos->vcre_dn;
    double *cpcoeffs_vcim_dn  = cpcoeffs_pos->vcim_dn;
    int ivcoef_form_up        = cpcoeffs_pos->ivcoef_form_up;
    int ivcoef_form_dn        = cpcoeffs_pos->ivcoef_form_dn;
    double *cpcoeffs_cmass    = cpcoeffs_info->cmass;
    int *cpcoeffs_ioff_up     = cpcoeffs_info->ioff_upt;
    int *cpcoeffs_ioff_dn     = cpcoeffs_info->ioff_dnt;
    int ncoef     = cpcoeffs_info->ncoef;
    int nstate_up = cpcoeffs_info->nstate_up;
    int nstate_dn = cpcoeffs_info->nstate_dn;
    int icmoff_up        = cpcoeffs_info->icoef_start_up-1;
    int icmoff_dn        = cpcoeffs_info->icoef_start_dn-1;
    int ncoef_up,ncoef_dn;

/*========================================================================*/
/* 0) Checks and Assigns */

  if(np_states>1){
    if((ivcoef_form_up)!=1){
     printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
     printf("Up CP vectors are not in transposed form \n");
     printf("in get_cpke \n");
     printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
     fflush(stdout);
     exit(1);
    }/*endif*/
    if(cp_lsda==1){
     if((ivcoef_form_dn)!=1){
      printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
      printf("Up CP vectors are not in transposed form \n");
      printf("in get_cpke \n");
      printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
      fflush(stdout);
      exit(1);
     }/*endif*/
    }/*endif*/
  }/*endif*/

  if(np_states==1){
     ncoef_up     = cpcoeffs_info->ncoef;
     ncoef_dn     = cpcoeffs_info->ncoef;
  }else{
     ncoef_up     = cpcoeffs_info->nstate_ncoef_proc_up;
     ncoef_dn     = cpcoeffs_info->nstate_ncoef_proc_dn;
  }/*endif*/

/*========================================================================*/
  /* I) Get the CP kinetic energy (Although up and down are separate here,*/ 
  /*    ultimately the total fict KE is the sum of the two) */

   kinet_cp_up   = 0.0;
   for(is=1;is<=nstate_up;is++) {
     for(i=1;i<=ncoef_up;i++) {
       icoef = i+cpcoeffs_ioff_up[is];
       kinet_cp_up += (cpcoeffs_vcre_up[icoef]*cpcoeffs_vcre_up[icoef]
                      *cpcoeffs_cmass[(i+icmoff_up)]);
       kinet_cp_up += (cpcoeffs_vcim_up[icoef]*cpcoeffs_vcim_up[icoef]
                      *cpcoeffs_cmass[(i+icmoff_up)]);
     }/*endfor*/
   }/*endfor*/

   kinet_cp_dn = 0.0;
   if( (cp_lsda == 1) && (nstate_dn != 0) ){
     for(is=1;is<=nstate_dn;is++) {
       for(i=1;i<=ncoef_dn;i++) {
         icoef = i+cpcoeffs_ioff_dn[is];
         kinet_cp_dn += (cpcoeffs_vcre_dn[icoef]*cpcoeffs_vcre_dn[icoef]
                        *cpcoeffs_cmass[(i+icmoff_dn)]);
         kinet_cp_dn += (cpcoeffs_vcim_dn[icoef]*cpcoeffs_vcim_dn[icoef]
                        *cpcoeffs_cmass[(i+icmoff_dn)]);
       }/*endfor*/
     }/*endfor*/
   }/* endif */

   kinet_cp_up /= 2.0;
   kinet_cp_dn /= 2.0;
   stat_avg->kinet_cp_up = kinet_cp_up;
   stat_avg->kinet_cp_dn = kinet_cp_dn;
   stat_avg->kinet_cp    = kinet_cp_up + kinet_cp_dn;

/*------------------------------------------------------------------------*/
   /*end routine*/}
/*========================================================================*/




/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
 
void nhc_cp_potkin(CPTHERM_INFO *cptherm_info,CPTHERM_POS *cptherm_pos,
                   STAT_AVG *stat_avg,int myid_state,MPI_Comm comm_states)

/*========================================================================*/
{/*begin routine*/
/*========================================================================*/
/*             Local variable declarations                                */
#include "../typ_defs/typ_mask.h"

int inhc,ichain;
int len_c_nhc,num_c_nhc;
double kinet_nhc,vpotnhc;
double **cptherm_c_nhc     = cptherm_pos->c_nhc;
double **cptherm_vc_nhc    = cptherm_pos->vc_nhc;
double **cptherm_cmass_nhc = cptherm_info->cmass_nhc;
double **cptherm_c_gkt     = cptherm_info->c_gkt;

/*==========================================================================*/
/* I) nhc stuff */

    vpotnhc     = 0.0;
    kinet_nhc   = 0.0;
    if(myid_state==0){
      len_c_nhc   = cptherm_info->len_c_nhc;
      num_c_nhc   = cptherm_info->num_c_nhc_proc;
      for(ichain=1;ichain<=len_c_nhc;ichain++){
        for(inhc=1;inhc<=num_c_nhc;inhc++){
          kinet_nhc += (cptherm_cmass_nhc[ichain][inhc]
                 *cptherm_vc_nhc[ichain][inhc]*cptherm_vc_nhc[ichain][inhc]);
          vpotnhc += (cptherm_c_gkt[ichain][inhc]*
                                cptherm_c_nhc[ichain][inhc]);   
        }/*endfor*/
      }/*endfor*/
    }/*endif*/
    kinet_nhc /= 2.0;
    stat_avg->vpotnhc_cp   = vpotnhc;
    stat_avg->kinet_nhc_cp = kinet_nhc;

/*------------------------------------------------------------------------*/
/*end routine*/}
/*========================================================================*/



/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
 
void nhc_cp_potkin_massiv(CPTHERM_INFO *cptherm_info,CPTHERM_POS *cptherm_pos,
                          STAT_AVG *stat_avg)

/*========================================================================*/
{/*begin routine*/
/*========================================================================*/
/*             Local variable declarations                                */
#include "../typ_defs/typ_mask.h"

int inhc,ichain;
int len_c_nhc,num_c_nhc;
double kinet_nhc,vpotnhc;
double c_nhc     = cptherm_pos->c_nhc_massiv;
double **vc_nhc  = cptherm_pos->vc_nhc;
double cmass_nhc = cptherm_info->cmass_nhc_massiv;
double c_gkt     = cptherm_info->c_gkt_massiv;

/*==========================================================================*/
/* I) nhc stuff */

    vpotnhc     = 0.0;
    kinet_nhc   = 0.0;
    len_c_nhc   = cptherm_info->len_c_nhc;
    num_c_nhc   = cptherm_info->num_c_nhc_proc;
    for(ichain=1;ichain<=len_c_nhc;ichain++){
      for(inhc=1;inhc<=num_c_nhc;inhc++){
        kinet_nhc += (cmass_nhc*vc_nhc[ichain][inhc]*vc_nhc[ichain][inhc]);
      }/*endfor*/
    }/*endfor*/
    vpotnhc = c_gkt*c_nhc;   
    kinet_nhc /= 2.0;
    stat_avg->vpotnhc_cp   = vpotnhc;
    stat_avg->kinet_nhc_cp = kinet_nhc;

/*------------------------------------------------------------------------*/
/*end routine*/}
/*========================================================================*/




/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void init_cp_NHC(CPTHERM_INFO *cptherm_info,CPTHERM_POS *cptherm_pos,
                 CPCOEFFS_INFO *cpcoeffs_info,CPCOEFFS_POS *cpcoeffs_pos,
                 CPSCR *cpscr,int cp_lsda,MPI_Comm comm_states,int np_states,
                 int myid_state)

/*========================================================================*/
{/*begin routine*/
/*========================================================================*/
/*             Local variable declarations                                */
#include "../typ_defs/typ_mask.h"


    int icoef,inhc,ichain;              /* Num: for loop counters */
    int len_c_nhc,len_c_nhcm1,len_c_nhcp1;    /* Num: length of chains  */
    int len_c_nhcm2;
    int iii;
    int is, i,ichain1;
    int myid=myid_state;

/* Define local pointers                                          */
      double *cpscr_coef_kin     = cpscr->cpscr_therm.coef_kin;
      double *cpscr_coef_kin_temp= cpscr->cpscr_therm.sc_cp;
      int *icmapup_nhc  = cptherm_info->icmapup_nhc;
      int *icmapdn_nhc  = cptherm_info->icmapdn_nhc;
      double *cpcoeffs_cmass   = cpcoeffs_info->cmass;
      double *cpcoeffs_vcre_up = cpcoeffs_pos->vcre_up;
      double *cpcoeffs_vcim_up = cpcoeffs_pos->vcim_up;
      double *cpcoeffs_vcre_dn = cpcoeffs_pos->vcre_dn;
      double *cpcoeffs_vcim_dn = cpcoeffs_pos->vcim_dn;
      int ivcoef_form_up        = cpcoeffs_pos->ivcoef_form_up;
      int ivcoef_form_dn        = cpcoeffs_pos->ivcoef_form_dn;
      double **fc_nhc  = cptherm_pos->fc_nhc;
      double **vc_nhc  = cptherm_pos->vc_nhc;
      double **c_gkt   = cptherm_info->c_gkt;
      double **cmass_nhc = cptherm_info->cmass_nhc;
      int *cpcoeffs_ioff_up     = cpcoeffs_info->ioff_upt;
      int *cpcoeffs_ioff_dn     = cpcoeffs_info->ioff_dnt;
      int nstate_up = cpcoeffs_info->nstate_up;
      int nstate_dn = cpcoeffs_info->nstate_dn;
      int ncoef = cpcoeffs_info->ncoef;
      int num_c_nhc = cptherm_info->num_c_nhc_proc;
      int icmoff_up        = cpcoeffs_info->icoef_start_up-1;
      int icmoff_dn        = cpcoeffs_info->icoef_start_dn-1;
      int ncoef_up,ncoef_dn;

/*========================================================================*/
/* 0) Checks and Assigns */

  if(np_states>1){
    if((ivcoef_form_up)!=1){
     printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
     printf("Up CP vectors are not in transposed form \n");
     printf("in int_cp_NHC \n");
     printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
     fflush(stdout);
     exit(1);
    }/*endif*/
    if(cp_lsda==1){
     if((ivcoef_form_dn)!=1){
      printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
      printf("Up CP vectors are not in transposed form \n");
      printf("in int_NVE_cp \n");
      printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
      fflush(stdout);
      exit(1);
     }/*endif*/
    }/*endif*/
  }/*endif*/

  if(np_states==1){
     ncoef_up     = cpcoeffs_info->ncoef;
     ncoef_dn     = cpcoeffs_info->ncoef;
  }else{
     ncoef_up     = cpcoeffs_info->nstate_ncoef_proc_up;
     ncoef_dn     = cpcoeffs_info->nstate_ncoef_proc_dn;
  }/*endif*/

/*==========================================================================*/
/* I) Map the ke of each particle to the appropriate nose-hoover chain      */
/*     and assemble the total particle ke associated                        */
/*     with each chain. The num_c_nhc+1 cpthermo is the null cptherm        */

    for(inhc=1;inhc<=num_c_nhc+1;inhc++){
      cpscr_coef_kin[inhc] = 0.0;
    /*endfor*/}

    for(is=1;is<=nstate_up;is++) {
     for(i=1;i<=ncoef_up;i++) {
       icoef = i+cpcoeffs_ioff_up[is];
       cpscr_coef_kin[icmapup_nhc[is]] += 
          cpcoeffs_cmass[(i+icmoff_up)]
         *cpcoeffs_vcre_up[icoef]*cpcoeffs_vcre_up[icoef];
     }/*endfor*/
     for(i=1;i<=ncoef_up;i++) {
      icoef = i+cpcoeffs_ioff_up[is];
      cpscr_coef_kin[icmapup_nhc[is]] += 
        cpcoeffs_cmass[(i+icmoff_up)]
       *cpcoeffs_vcim_up[icoef]*cpcoeffs_vcim_up[icoef];
     }/*endfor*/
    }/* endfor */

    if( (cp_lsda == 1) && (nstate_dn != 0) ){
      for(is=1;is<=nstate_dn;is++) {
        for(i=1;i<=ncoef_dn;i++) {
          icoef = i+cpcoeffs_ioff_dn[is];
          cpscr_coef_kin[icmapdn_nhc[is]] += 
           cpcoeffs_cmass[(i+icmoff_dn)]
          *cpcoeffs_vcre_dn[icoef]*cpcoeffs_vcre_dn[icoef];
        }/*endfor*/
        for(i=1;i<=ncoef_dn;i++) {
          icoef = i+cpcoeffs_ioff_dn[is];
          cpscr_coef_kin[icmapdn_nhc[is]] += 
          cpcoeffs_cmass[(i+icmoff_dn)]
         *cpcoeffs_vcim_dn[icoef]*cpcoeffs_vcim_dn[icoef];
        }/*endfor*/
      }/* endfor */
    }/* endif cp_lsda */


 printf("in init_nhc %d 0\n",myid);


    if(np_states>1){
      Allreduce(&(cpscr_coef_kin[1]), &(cpscr_coef_kin_temp[1]),num_c_nhc,
                   MPI_DOUBLE,MPI_SUM,0,comm_states);
      for(inhc=1;inhc<=num_c_nhc;inhc++){
        cpscr_coef_kin[inhc] = cpscr_coef_kin_temp[inhc];
      }/*endfor*/

    }/*endif*/

 printf("in init_nhc %d 1\n",myid);

/*==========================================================================*/

    len_c_nhc   = (cptherm_info->len_c_nhc);
    len_c_nhcm1 = (cptherm_info->len_c_nhc)-1;
    len_c_nhcm2 = (cptherm_info->len_c_nhc)-2;
    len_c_nhcp1 = (cptherm_info->len_c_nhc)+1;
    for(inhc=1;inhc<=num_c_nhc;inhc++){
      fc_nhc[1][inhc] =  (cpscr_coef_kin[inhc]-c_gkt[1][inhc])
                         /cmass_nhc[1][inhc];
    }/*endfor*/
    for(ichain=1;ichain<=len_c_nhcm1;ichain++){
     ichain1 = ichain+1;
     if(ichain1!=len_c_nhcm1){
       for(inhc=1;inhc<=num_c_nhc;inhc++){
         fc_nhc[ichain1][inhc] = (cmass_nhc[ichain][inhc]*
               vc_nhc[ichain][inhc]*vc_nhc[ichain][inhc]
                     -c_gkt[ichain1][inhc])
                      /cmass_nhc[ichain1][inhc];
       }/*endfor*/
     }else{
       for(inhc=1;inhc<=num_c_nhc;inhc++){
         fc_nhc[len_c_nhcm1][inhc] = (cmass_nhc[len_c_nhc][inhc]*
            vc_nhc[len_c_nhc][inhc]*vc_nhc[len_c_nhc][inhc]
                      +cmass_nhc[len_c_nhcm2][inhc]*
            vc_nhc[len_c_nhcm2][inhc]*vc_nhc[len_c_nhcm2][inhc]
            - c_gkt[len_c_nhcm1][inhc])/cmass_nhc[len_c_nhcm1][inhc];
       }/*endfor*/
     }/*endif*/
    }/*endfor*/

 printf("in init_nhc %d 2\n",myid);

/*--------------------------------------------------------------------------*/
   /*end routine*/}
/*==========================================================================*/






/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void apply_c_nhc(CPTHERM_INFO *cptherm_info,CPTHERM_POS *cptherm_pos,
                 CPCOEFFS_INFO *cpcoeffs_info,CPCOEFFS_POS *cpcoeffs_pos,
                 CPSCR *cpscr,int cp_lsda,COMMUNICATE *communicate)

/*==========================================================================*/
/*             Begin Routine                                                */
{/*begin routine*/
/*========================================================================*/
/*             Local variable declarations                                */
#include "../typ_defs/typ_mask.h"

    int icoef,inhc,ichain;              /* Num: for loop counters */
    int iresn,iyosh;                    /* Num: for loop counters */
    int len_c_nhc,len_c_nhcm1,len_c_nhcp1;    /* Num: length of chains  */
    int len_c_nhcm2;
    double arg,aa;                      /* Num: scalar temps      */
    int iii;
    int nstate_up,nstate_dn;
    int ncoef;
    int ncoef_up_tot,ncoef_dn_tot;
    int num_c_nhc;
    int is, i,ichain1,itemp,jtemp;

/* Define local pointers                                          */
      double *cpscr_coef_kin   = cpscr->cpscr_therm.coef_kin;
      double *cpscr_sc         = cpscr->cpscr_therm.sc_cp;
      double *vc_nhc_tmp       = cpscr->cpscr_therm.coef_kin;
      double *c_nhc_tmp        = cpscr->cpscr_therm.sc_cp;
      int *icmapup_nhc  = cptherm_info->icmapup_nhc;
      int *icmapdn_nhc  = cptherm_info->icmapdn_nhc;
      double *cpcoeffs_cmass   = cpcoeffs_info->cmass;
      double *cpcoeffs_vcre_up = cpcoeffs_pos->vcre_up;
      double *cpcoeffs_vcim_up = cpcoeffs_pos->vcim_up;
      double *cpcoeffs_vcre_dn = cpcoeffs_pos->vcre_dn;
      double *cpcoeffs_vcim_dn = cpcoeffs_pos->vcim_dn;
      double **fc_nhc  = cptherm_pos->fc_nhc;
      double **vc_nhc  = cptherm_pos->vc_nhc;
      double **c_nhc   = cptherm_pos->c_nhc;
      double **c_gkt   = cptherm_info->c_gkt;
      double **cmass_nhc = cptherm_info->cmass_nhc;
      double *wdti2     = cptherm_info->wdti2;
      double *wdti4     = cptherm_info->wdti4;
      double *wdti8     = cptherm_info->wdti8;
      int *cpcoeffs_ioff_up     = cpcoeffs_info->ioff_upt;
      int *cpcoeffs_ioff_dn     = cpcoeffs_info->ioff_dnt;
      int myid_state            = communicate->myid_state;
      MPI_Comm comm_states      = communicate->comm_states;
      int icmoff_up        = cpcoeffs_info->icoef_start_up-1;
      int icmoff_dn        = cpcoeffs_info->icoef_start_dn-1;
      int np_states             = communicate->np_states;
      int ncoef_up,ncoef_dn;

     int ivcoef_form_up        = cpcoeffs_pos->ivcoef_form_up;
     int ivcoef_form_dn        = cpcoeffs_pos->ivcoef_form_dn;

/*==========================================================================*/
/* 0) Checks                                                                */

  if(np_states>1){
    if(ivcoef_form_up!=1){
     printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
     printf("Up CP velocities are not in transposed form \n");
     printf("on state processor %d in apply_c_nhc \n",myid_state);
     printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
     fflush(stdout);
     exit(1);
    }/*endif*/

    if(cp_lsda==1){
    if(ivcoef_form_dn!=1){
      printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
      printf("Dn CP velocities are not in transposed form \n");
      printf("on state processor %d in int_NVE_cp \n",myid_state);
      printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
      fflush(stdout);
      exit(1);
     }/*endif*/
    }/*endif*/
  }/*endif*/

/*==========================================================================*/
/* 0) Simple definitions */

    nstate_up    = cpcoeffs_info->nstate_up;
    nstate_dn    = cpcoeffs_info->nstate_dn;
    if(np_states==1){
     ncoef_up       = cpcoeffs_info->ncoef;
     ncoef_dn       = cpcoeffs_info->ncoef;
    }else{
     ncoef_up     = cpcoeffs_info->nstate_ncoef_proc_up;
     ncoef_dn     = cpcoeffs_info->nstate_ncoef_proc_dn;
    }/*endif*/
    num_c_nhc    = cptherm_info->num_c_nhc_proc;

/*==========================================================================*/
/* I) Map the ke of each particle to the appropriate nose-hoover chain      */
/*     and assemble the total particle ke associated                        */
/*     with each chain. The num_c_nhc+1 cpthermo is the null cptherm        */

    for(inhc=1;inhc<=num_c_nhc+1;inhc++){
      cpscr_coef_kin[inhc] = 0.0;
    /*endfor*/}

    for(is=1;is<=nstate_up;is++) {
     for(i=1;i<=ncoef_up;i++) {
       icoef = i+cpcoeffs_ioff_up[is];
       cpscr_coef_kin[icmapup_nhc[is]] += 
          cpcoeffs_cmass[(i+icmoff_up)]
          *cpcoeffs_vcre_up[icoef]*cpcoeffs_vcre_up[icoef];
     }/*endfor*/
     for(i=1;i<=ncoef_up;i++) {
      icoef = i+cpcoeffs_ioff_up[is];
      cpscr_coef_kin[icmapup_nhc[is]] += 
        cpcoeffs_cmass[(i+icmoff_up)]
        *cpcoeffs_vcim_up[icoef]*cpcoeffs_vcim_up[icoef];
     }/*endfor*/
    }/* endfor */

    if( (cp_lsda == 1) && (nstate_dn != 0) ){
      for(is=1;is<=nstate_dn;is++) {
        for(i=1;i<=ncoef_dn;i++) {
          icoef = i+cpcoeffs_ioff_dn[is];
          cpscr_coef_kin[icmapdn_nhc[is]] += 
           cpcoeffs_cmass[(i+icmoff_dn)]
          *cpcoeffs_vcre_dn[icoef]*cpcoeffs_vcre_dn[icoef];
        }/*endfor*/
        for(i=1;i<=ncoef_dn;i++) {
          icoef = i+cpcoeffs_ioff_dn[is];
          cpscr_coef_kin[icmapdn_nhc[is]] += 
          cpcoeffs_cmass[(i+icmoff_dn)]
          *cpcoeffs_vcim_dn[icoef]*cpcoeffs_vcim_dn[icoef];
        }/*endfor*/
      }/* endfor */
    }/* endif cp_lsda */

    if(np_states>1){

       Allreduce(&(cpscr_coef_kin[1]), &(cpscr_sc[1]),num_c_nhc,
                 MPI_DOUBLE,MPI_SUM,0,comm_states);
       for(inhc=1;inhc<=num_c_nhc;inhc++){
          cpscr_coef_kin[inhc] = cpscr_sc[inhc];
       }/*endfor*/

    }/*endif*/

/*==========================================================================*/
/* III) Get the force on the first NHC in each chain                        */

    for(inhc=1;inhc<=num_c_nhc;inhc++){
      fc_nhc[1][inhc] = (cpscr_coef_kin[inhc]-c_gkt[1][inhc])
                                /cmass_nhc[1][inhc];
    }/*endfor*/
  
/*==========================================================================*/
/* IV) Apply the nhc evolution operator using RESPA                         */

    len_c_nhc   = (cptherm_info->len_c_nhc);
    len_c_nhcm1 = (cptherm_info->len_c_nhc)-1;
    len_c_nhcm2 = (cptherm_info->len_c_nhc)-2;
    len_c_nhcp1 = (cptherm_info->len_c_nhc)+1;
    for(inhc=1;inhc<=num_c_nhc+1;inhc++){
      cpscr_sc[inhc]     = 1.0;
    }/*endfor*/
    for(iresn=1;iresn<=cptherm_info->nres_c_nhc;iresn++){
      for(iyosh=1;iyosh<=cptherm_info->nyosh_c_nhc;iyosh++){
/*--------------------------------------------------------------------------*/
/*  1) Evolve the last cptherm velocity in each chain                       */

       if(len_c_nhc>2){
        for(inhc=1;inhc<=num_c_nhc;inhc++){
          arg = -wdti8[iyosh]*vc_nhc[len_c_nhcm1][inhc];
          aa = exp(arg);
          vc_nhc[len_c_nhc][inhc] =  vc_nhc[len_c_nhc][inhc]*aa*aa
                                   + wdti4[iyosh]*fc_nhc[len_c_nhc][inhc]*aa;
        }/*endfor*/
        for(inhc=1;inhc<=num_c_nhc;inhc++){
          fc_nhc[len_c_nhcm1][inhc] = 
                (cmass_nhc[len_c_nhc][inhc]*
                 vc_nhc[len_c_nhc][inhc]*vc_nhc[len_c_nhc][inhc]
                +cmass_nhc[len_c_nhcm2][inhc]*
                 vc_nhc[len_c_nhcm2][inhc]*vc_nhc[len_c_nhcm2][inhc]
                -c_gkt[len_c_nhcm1][inhc])/cmass_nhc[len_c_nhcm1][inhc];
        }/*endfor*/
      }else{
         for(inhc=1;inhc<=num_c_nhc;inhc++){
           vc_nhc[len_c_nhc][inhc] =  vc_nhc[len_c_nhc][inhc]
                                   + wdti4[iyosh]*fc_nhc[len_c_nhc][inhc];
         }/*endfor*/
      }/*endif*/

/*--------------------------------------------------------------------------*/
/*  2) Evolve the last-1 to the first cpthermo velocitiy in each chain      */
        for(ichain=1;ichain<=len_c_nhcm1;ichain++){
          itemp = (len_c_nhc-ichain);
          jtemp = (len_c_nhcp1-ichain);
          for(inhc=1;inhc<=num_c_nhc;inhc++){
            arg = -wdti8[iyosh]*vc_nhc[jtemp][inhc];
            aa = exp(arg);
            vc_nhc[itemp][inhc] = vc_nhc[itemp][inhc]*aa*aa
                                + wdti4[iyosh]*fc_nhc[itemp][inhc]*aa;
          }/*endfor*/
        }/*endfor*/
/*--------------------------------------------------------------------------*/
/*  3) Evolve the particle velocities (by adding to the scaling factor)     */
        for(inhc=1;inhc<=num_c_nhc;inhc++){
          arg = -wdti2[iyosh]*vc_nhc[1][inhc];
          aa = exp(arg);
          cpscr_sc[inhc]       *= aa;
          cpscr_coef_kin[inhc] *= (aa*aa);
        }/*endfor*/

/*--------------------------------------------------------------------------*/
/*  4) Evolve the cptherm positions                                         */
        for(ichain=1;ichain<=len_c_nhc;ichain++){
          for(inhc=1;inhc<=num_c_nhc;inhc++){
            c_nhc[ichain][inhc] += 
                    vc_nhc[ichain][inhc]*wdti2[iyosh];
          }/*endfor*/
        }/*endfor*/
/*--------------------------------------------------------------------------*/
/*  5) Evolve the 1 to last-1 cptherm velocity in each chain                */
/*     calculting cptherm forces as you go along                            */
/*     careful with the len-1st puppy.                                      */
        for(inhc=1;inhc<=num_c_nhc;inhc++){
           fc_nhc[1][inhc] =  (cpscr_coef_kin[inhc]-c_gkt[1][inhc])
                             /cmass_nhc[1][inhc];
        }/*endfor*/
        for(ichain=1;ichain<=len_c_nhcm1;ichain++){
          ichain1 = ichain+1;
          for(inhc=1;inhc<=num_c_nhc;inhc++){
            arg = -wdti8[iyosh]*vc_nhc[ichain1][inhc];
            aa = exp(arg);
            vc_nhc[ichain][inhc] = vc_nhc[ichain][inhc]*aa*aa
                                 + wdti4[iyosh]*fc_nhc[ichain][inhc]*aa;
          }/*endfor*/
          if(ichain1!=len_c_nhcm1 || len_c_nhc<=2){
            for(inhc=1;inhc<=num_c_nhc;inhc++){
              fc_nhc[ichain1][inhc] = 
                   (cmass_nhc[ichain][inhc]*
                    vc_nhc[ichain][inhc]*vc_nhc[ichain][inhc]
                   -c_gkt[ichain1][inhc])/cmass_nhc[ichain1][inhc];
            }/*endfor*/
          }else{
            for(inhc=1;inhc<=num_c_nhc;inhc++){
              fc_nhc[len_c_nhcm1][inhc] = 
                     (cmass_nhc[len_c_nhc][inhc]*
                      vc_nhc[len_c_nhc][inhc]*vc_nhc[len_c_nhc][inhc]
                     +cmass_nhc[len_c_nhcm2][inhc]*
                      vc_nhc[len_c_nhcm2][inhc]*vc_nhc[len_c_nhcm2][inhc]
                     -c_gkt[len_c_nhcm1][inhc])/cmass_nhc[len_c_nhcm1][inhc];
            }/*endfor*/
          }/*endif*/
        }/*endfor*/
/*--------------------------------------------------------------------------*/
/*  6) Evolve the last cptherm velocotiy in each chain                      */
       if(len_c_nhc>2){
        for(inhc=1;inhc<=num_c_nhc;inhc++){
          arg = -wdti8[iyosh]*vc_nhc[len_c_nhcm1][inhc];
          aa = exp(arg);
          vc_nhc[len_c_nhc][inhc] = vc_nhc[len_c_nhc][inhc]*aa*aa
                                   + wdti4[iyosh]*fc_nhc[len_c_nhc][inhc]*aa;
        }/*endfor*/
       }else{
         for(inhc=1;inhc<=num_c_nhc;inhc++){
          arg = -wdti8[iyosh]*vc_nhc[len_c_nhcm1][inhc];
          aa = exp(arg);
          vc_nhc[len_c_nhc][inhc] = vc_nhc[len_c_nhc][inhc]*aa*aa
                                   + wdti4[iyosh]*fc_nhc[len_c_nhc][inhc]*aa;
         }/*endfor*/
       }/*endif*/
/*--------------------------------------------------------------------------*/
/* 7) End Respa                                                             */
      /*endfor: iyosh */}
    /*endfor: iresn*/}

/*==========================================================================*/
/* V) Apply the accumulated scaling factor to the velocities                */

    for(is=1;is<=nstate_up;is++) {
     for(i=1;i<=ncoef_up;i++) {
       icoef = i+cpcoeffs_ioff_up[is];
       cpcoeffs_vcre_up[icoef] *= cpscr_sc[icmapup_nhc[is]];
     }/*endfor*/
     for(i=1;i<=ncoef_up;i++) {
      icoef = i+cpcoeffs_ioff_up[is];
      cpcoeffs_vcim_up[icoef] *= cpscr_sc[icmapup_nhc[is]];
     } /*endfor*/
    } /*endfor*/

   if( (cp_lsda == 1) && (nstate_dn != 0)){
    for(is=1;is<=nstate_dn;is++) {
     for(i=1;i<=ncoef_dn;i++) {
       icoef = i+cpcoeffs_ioff_dn[is];
       cpcoeffs_vcre_dn[icoef] *= cpscr_sc[icmapdn_nhc[is]];
     }/*endfor*/
     for(i=1;i<=ncoef_dn;i++) {
      cpcoeffs_vcim_dn[icoef] *= cpscr_sc[icmapdn_nhc[is]];
     }/*endfor*/
    }/*endfor*/
   }/* endif cp_lsda */

/*==========================================================================*/
/* V) Broadcast the results to keep consistency                             */

  if(np_states>1){

      for(ichain=1;ichain<=len_c_nhc;ichain++){
    if(myid_state==0){
        for(inhc=1;inhc<=num_c_nhc;inhc++){
         vc_nhc_tmp[inhc] = vc_nhc[ichain][inhc];
         c_nhc_tmp[inhc]  = c_nhc[ichain][inhc];
      }/*endfor : inhc*/
    }/*endif*/
    Bcast(&vc_nhc_tmp[1],num_c_nhc,MPI_DOUBLE,0,comm_states);
    Bcast(&c_nhc_tmp[1],num_c_nhc,MPI_DOUBLE,0,comm_states);
    for(inhc=1;inhc<=num_c_nhc;inhc++){
      vc_nhc[ichain][inhc] = vc_nhc_tmp[inhc];
      c_nhc[ichain][inhc]  = c_nhc_tmp[inhc];
    }/*endfor : inhc*/
   }/*endfor : ichain*/

 }/*endif*/

/*--------------------------------------------------------------------------*/
    }/*end routine*/
/*==========================================================================*/






/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void apply_c_nhc_massiv(CPTHERM_INFO *cptherm_info,CPTHERM_POS *cptherm_pos,
                      CPCOEFFS_INFO *cpcoeffs_info,CPCOEFFS_POS *cpcoeffs_pos,
                      CPSCR *cpscr,int cp_lsda,int myid_state, int np_states)

/*==========================================================================*/
/*             Begin Routine                                                */
{/*begin routine*/
/*========================================================================*/
/*             Local variable declarations                                */

    int icoef,inhc,ichain;              /* Num: for loop counters */
    int iresn,iyosh;                    /* Num: for loop counters */
    int len_c_nhc,len_c_nhcm1,len_c_nhcp1;    /* Num: length of chains  */
    int len_c_nhcm2;
    double arg,aa;                      /* Num: scalar temps      */
    int iii;
    int ncoef;
    int ncoef_up_tot,ncoef_dn_tot;
    int num_c_nhc,num_c_nhc_tot,ichainm1;
    int is, i,ichain1,itemp,jtemp,ktemp;
    int ioff;
    double fc_nhc;
    int nstate_dn,nstate_up;
/* Define local pointers                                          */
      double *cmass   = cpcoeffs_info->cmass;
      double *vcre_up = cpcoeffs_pos->vcre_up;
      double *vcim_up = cpcoeffs_pos->vcim_up;
      double *vcre_dn = cpcoeffs_pos->vcre_dn;
      double *vcim_dn = cpcoeffs_pos->vcim_dn;
      double **vc_nhc   = cptherm_pos->vc_nhc;
      double c_nhc      = cptherm_pos->c_nhc_massiv;
      double c_gkt      = cptherm_info->c_gkt_massiv;
      double c_gkt_div,cmass_nhci;
      double cmass_nhc  = cptherm_info->cmass_nhc_massiv;
      double *wdti2     = cptherm_info->wdti2;
      double *wdti4     = cptherm_info->wdti4;
      double *wdti8     = cptherm_info->wdti8;
      int *cpcoeffs_ioff_up     = cpcoeffs_info->ioff_upt;
      int *cpcoeffs_ioff_dn     = cpcoeffs_info->ioff_dnt;
      int icmoff_up        = cpcoeffs_info->icoef_start_up-1;
      int icmoff_dn        = cpcoeffs_info->icoef_start_dn-1;
      int ncoef_up,ncoef_dn;
      int ivcoef_form_up        = cpcoeffs_pos->ivcoef_form_up;
      int ivcoef_form_dn        = cpcoeffs_pos->ivcoef_form_dn;

/*==========================================================================*/
/* 0) Checks                                                                */

  if(np_states>1){
    if(ivcoef_form_up!=1){
     printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
     printf("Up CP velocities are not in transposed form \n");
     printf("on state processor %d in apply_c_nhc \n",myid_state);
     printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
     fflush(stdout);
     exit(1);
    }/*endif*/

    if(cp_lsda==1){
    if(ivcoef_form_dn!=1){
      printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
      printf("Dn CP velocities are not in transposed form \n");
      printf("on state processor %d in int_NVE_cp \n",myid_state);
      printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
      fflush(stdout);
      exit(1);
     }/*endif*/
    }/*endif*/
  }/*endif*/

/*==========================================================================*/
/* 0) Simple definitions */

    c_gkt_div  = c_gkt/cmass_nhc;
    cmass_nhci = 1.0/cmass_nhc;
    nstate_up    = cpcoeffs_info->nstate_up;
    nstate_dn    = cpcoeffs_info->nstate_dn;
    if(np_states==1){
     ncoef_up     = cpcoeffs_info->ncoef;
     ncoef_dn     = cpcoeffs_info->ncoef;
    }else{
     ncoef_up     = cpcoeffs_info->nstate_ncoef_proc_up;
     ncoef_dn     = cpcoeffs_info->nstate_ncoef_proc_dn;
    }
    ncoef_up_tot = ncoef_up*nstate_up;
    ncoef_dn_tot = ncoef_dn*nstate_dn;
    
    if(cp_lsda==0){ncoef_dn_tot=0;}
    num_c_nhc     = cptherm_info->num_c_nhc_proc;
    if(num_c_nhc!=2*(ncoef_up_tot+ncoef_dn_tot)){
       printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n");    
       printf("Internal error in massive coefficient thermostatting\n");
       printf("occured in routine apply_c_nhc_massiv'\n");
       printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n");    
       fflush(stdout);
       exit(1);
    }/*endif*/
    len_c_nhc   = (cptherm_info->len_c_nhc);
    len_c_nhcm1 = (cptherm_info->len_c_nhc)-1;
    len_c_nhcm2 = (cptherm_info->len_c_nhc)-2;
    len_c_nhcp1 = (cptherm_info->len_c_nhc)+1;

/*==========================================================================*/
/* IV) Apply the nhc evolution operator using RESPA                         */

    for(iresn=1;iresn<=cptherm_info->nres_c_nhc;iresn++){
      for(iyosh=1;iyosh<=cptherm_info->nyosh_c_nhc;iyosh++){
/*--------------------------------------------------------------------------*/
/*  1) Evolve the last cptherm velocity in each chain                       */
        if(num_c_nhc>1){
         for(inhc=1;inhc<=num_c_nhc;inhc++){
           fc_nhc  = vc_nhc[len_c_nhcm1][inhc]*vc_nhc[len_c_nhcm1][inhc]
                     -c_gkt_div;
           vc_nhc[len_c_nhc][inhc] = vc_nhc[len_c_nhc][inhc]
                                   + wdti4[iyosh]*fc_nhc;
         }/*endfor*/
        }/*endif*/
/*--------------------------------------------------------------------------*/
/*  2) Evolve the last-1 to the 2nd cpthermo velocitiy in each chain      */
        for(ichain=1;ichain<=len_c_nhcm2;ichain++){
          itemp = len_c_nhc-ichain;
          jtemp = itemp+1;
          ktemp = itemp-1;
          for(inhc=1;inhc<=num_c_nhc;inhc++){
            arg    = -wdti8[iyosh]*vc_nhc[jtemp][inhc];
            aa     = exp(arg);
            fc_nhc = vc_nhc[ktemp][inhc]*vc_nhc[ktemp][inhc]-c_gkt_div;
            vc_nhc[itemp][inhc] = vc_nhc[itemp][inhc]*aa*aa
                                + wdti4[iyosh]*fc_nhc*aa;
          }/*endfor*/
        }/*endfor*/
/*--------------------------------------------------------------------------*/
/*  2.5) Evolve the 1st cpthermo velocitiy in each chain                    */

        ioff = 0;
        for(is=1;is<=nstate_up;is++) {
         for(i=1;i<=ncoef_up;i++) {
            inhc = i+ioff;
            itemp = i+cpcoeffs_ioff_up[is];
            arg    = -wdti8[iyosh]*vc_nhc[2][inhc];
            aa     = exp(arg);
            fc_nhc = (vcre_up[itemp]*vcre_up[itemp]*cmass[(i+icmoff_up)]
                      -c_gkt)*cmass_nhci;
            vc_nhc[1][inhc] = vc_nhc[1][inhc]*aa*aa
                                + wdti4[iyosh]*fc_nhc*aa;
         }/*endfor*/
         ioff += ncoef_up;
        }/*endfor*/
        for(is=1;is<=nstate_up;is++) {
         for(i=1;i<=ncoef_up;i++) {
            inhc = i+ioff;
            itemp = i+cpcoeffs_ioff_up[is];
            arg    = -wdti8[iyosh]*vc_nhc[2][inhc];
            aa     = exp(arg);
            fc_nhc = (vcim_up[itemp]*vcim_up[itemp]*cmass[(i+icmoff_up)]
                     -c_gkt)*cmass_nhci;
            vc_nhc[1][inhc] = vc_nhc[1][inhc]*aa*aa
                                + wdti4[iyosh]*fc_nhc*aa;
         }/*endfor*/
         ioff += ncoef_up;
        }/*endfor*/

        if((cp_lsda==1)&&(nstate_dn>0)){
         for(is=1;is<=nstate_dn;is++) {
          for(i=1;i<=ncoef_dn;i++) {
            inhc = i+ioff;
            itemp  = i + cpcoeffs_ioff_dn[is];
            arg    = -wdti8[iyosh]*vc_nhc[2][inhc];
            aa     = exp(arg);
            fc_nhc = (vcre_dn[itemp]*vcre_dn[itemp]*cmass[(i+icmoff_dn)]
                      -c_gkt)*cmass_nhci;
            vc_nhc[1][inhc] = vc_nhc[1][inhc]*aa*aa
                                + wdti4[iyosh]*fc_nhc*aa;
           }/*endfor*/
           ioff += ncoef_dn;
          }/*endfor*/
          for(is=1;is<=nstate_dn;is++) {
           for(i=1;i<=ncoef_dn;i++) {
            inhc = i+ioff;
            itemp  = i + cpcoeffs_ioff_dn[is];
            arg    = -wdti8[iyosh]*vc_nhc[2][inhc];
            aa     = exp(arg);
            fc_nhc = (vcim_dn[itemp]*vcim_dn[itemp]*cmass[(i+icmoff_dn)]
                      -c_gkt)*cmass_nhc;
            vc_nhc[1][inhc] = vc_nhc[1][inhc]*aa*aa
                                + wdti4[iyosh]*fc_nhc*aa;
           }/*endfor*/
           ioff += ncoef_dn;
          }/*endfor*/
	}/*endif:lsda*/
/*--------------------------------------------------------------------------*/
/*  3) Evolve the particle velocities (by adding to the scaling factor)     */

        ioff=0;
        for(is=1;is<=nstate_up;is++) {
         for(i=1;i<=ncoef_up;i++) {
            inhc = i+ioff;
            itemp = i+cpcoeffs_ioff_up[is];
            arg    = -wdti2[iyosh]*vc_nhc[1][inhc];
            aa     = exp(arg);
            vcre_up[itemp] *=aa;
         }/*endfor*/
         ioff+=ncoef_up;
        }/*endfor*/

        for(is=1;is<=nstate_up;is++) {
         for(i=1;i<=ncoef_up;i++) {
            inhc = i+ioff;
            itemp  = i+cpcoeffs_ioff_up[is];
            arg    = -wdti2[iyosh]*vc_nhc[1][inhc];
            aa     = exp(arg);
            vcim_up[itemp] *=aa;
         }/*endfor*/
         ioff += ncoef_up;
        }/*endfor*/

       if(cp_lsda == 1){
        for(is=1;is<=nstate_dn;is++) {
         for(i=1;i<=ncoef_dn;i++) {
            inhc = i+ioff;
            itemp  = i+cpcoeffs_ioff_dn[is];
            arg    = -wdti2[iyosh]*vc_nhc[1][inhc];
            aa     = exp(arg);
            vcre_dn[itemp] *=aa;
         }/*endfor*/
         ioff += ncoef_dn;
        }/*endfor*/

        for(is=1;is<=nstate_dn;is++) {
         for(i=1;i<=ncoef_dn;i++) {
            inhc = i+ioff;
            itemp  = i+cpcoeffs_ioff_dn[is];
            arg    = -wdti2[iyosh]*vc_nhc[1][inhc];
            aa     = exp(arg);
            vcim_dn[itemp] *=aa;
         }/*endfor*/
         ioff += ncoef_dn;
        }/*endfor*/
       }/* endif lsda */

/*--------------------------------------------------------------------------*/
/*  4) Evolve the cptherm positions                                         */
        for(ichain=1;ichain<=len_c_nhc;ichain++){
          for(inhc=1;inhc<=num_c_nhc;inhc++){
            c_nhc += vc_nhc[ichain][inhc]*wdti2[iyosh];
          }/*endfor*/
        }/*endfor*/
/*--------------------------------------------------------------------------*/
/*  4.5) Evolve the 1st cpthermo velocitiy in each chain                    */

        ioff = 0;
        for(is=1;is<=nstate_up;is++) {
         for(i=1;i<=ncoef_up;i++) {
            inhc = i+ioff;
            itemp = i+cpcoeffs_ioff_up[is];
            arg    = -wdti8[iyosh]*vc_nhc[2][inhc];
            aa     = exp(arg);
            fc_nhc = (vcre_up[itemp]*vcre_up[itemp]*cmass[(i+icmoff_up)]
                      -c_gkt)*cmass_nhci;
            vc_nhc[1][inhc] = vc_nhc[1][inhc]*aa*aa
                                + wdti4[iyosh]*fc_nhc*aa;
         }/*endfor*/
         ioff += ncoef_up;
        }/*endfor*/
        for(is=1;is<=nstate_up;is++) {
         for(i=1;i<=ncoef_up;i++) {
            inhc = i+ioff;
            itemp = i+cpcoeffs_ioff_up[is];
            arg    = -wdti8[iyosh]*vc_nhc[2][inhc];
            aa     = exp(arg);
            fc_nhc = (vcim_up[itemp]*vcim_up[itemp]*cmass[(i+icmoff_up)]
                     -c_gkt)*cmass_nhci;
            vc_nhc[1][inhc] = vc_nhc[1][inhc]*aa*aa
                                + wdti4[iyosh]*fc_nhc*aa;
         }/*endfor*/
         ioff += ncoef_up;
        }/*endfor*/

        if((cp_lsda==1)&&(nstate_dn>0)){
         for(is=1;is<=nstate_dn;is++) {
          for(i=1;i<=ncoef_dn;i++) {
            inhc = i+ioff;
            itemp  = i + cpcoeffs_ioff_dn[is];
            arg    = -wdti8[iyosh]*vc_nhc[2][inhc];
            aa     = exp(arg);
            fc_nhc = (vcre_dn[itemp]*vcre_dn[itemp]*cmass[(i+icmoff_dn)]
                      -c_gkt)*cmass_nhci;
            vc_nhc[1][inhc] = vc_nhc[1][inhc]*aa*aa
                                + wdti4[iyosh]*fc_nhc*aa;
           }/*endfor*/
           ioff += ncoef_dn;
          }/*endfor*/
          for(is=1;is<=nstate_dn;is++) {
           for(i=1;i<=ncoef_dn;i++) {
            inhc = i+ioff;
            itemp  = i + cpcoeffs_ioff_dn[is];
            arg    = -wdti8[iyosh]*vc_nhc[2][inhc];
            aa     = exp(arg);
            fc_nhc = (vcim_dn[itemp]*vcim_dn[itemp]*cmass[(i+icmoff_dn)]
                      -c_gkt)*cmass_nhc;
            vc_nhc[1][inhc] = vc_nhc[1][inhc]*aa*aa
                                + wdti4[iyosh]*fc_nhc*aa;
           }/*endfor*/
           ioff += ncoef_dn;
          }/*endfor*/
	}/*endif:lsda*/

/*--------------------------------------------------------------------------*/
/*  5) Evolve the 2 to last-1 cptherm velocity in each chain                */
/*     calculting cptherm forces as you go along                            */
/*     careful with the len-1st puppy.                                      */

        for(ichain=2;ichain<=len_c_nhcm1;ichain++){
          ichain1 = ichain+1;
          ichainm1 = ichain-1;
          for(inhc=1;inhc<=num_c_nhc;inhc++){
            arg    = -wdti8[iyosh]*vc_nhc[ichain1][inhc];
            aa     = exp(arg);
            fc_nhc = vc_nhc[ichainm1][inhc]*vc_nhc[ichainm1][inhc] - c_gkt_div;
            vc_nhc[ichain][inhc] = vc_nhc[ichain][inhc]*aa*aa
                                 + wdti4[iyosh]*fc_nhc*aa;
          }/*endfor*/
        }/*endfor*/
/*--------------------------------------------------------------------------*/
/*  6) Evolve the last cptherm velocotiy in each chain                      */
       if(num_c_nhc>1){
        for(inhc=1;inhc<=num_c_nhc;inhc++){
          fc_nhc   = vc_nhc[len_c_nhcm1][inhc]*vc_nhc[len_c_nhcm1][inhc]
                     -c_gkt_div;
          vc_nhc[len_c_nhc][inhc] =  vc_nhc[len_c_nhc][inhc]
                                   + wdti4[iyosh]*fc_nhc;
        }/*endfor*/
       }/*endif*/
/*--------------------------------------------------------------------------*/
/* 7) End Respa                                                             */
     }/*endfor: iyosh */
   }/*endfor: iresn*/
   cptherm_pos->c_nhc_massiv = c_nhc;


/*--------------------------------------------------------------------------*/
}/*end routine*/
/*==========================================================================*/






















