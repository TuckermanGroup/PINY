/*===================================================================*/
/*ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*===================================================================*/
/*                                                                   */
/*                         PI_MD:                                    */
/*             The future of simulation technology                   */
/*             ------------------------------------                  */
/*                     Module: samp_vel.c                            */
/*                                                                   */
/* These subprograms sample the velocities                           */
/*                                                                   */
/*===================================================================*/
/*ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*===================================================================*/

#include "standard_include.h"
#include "../typ_defs/typedefs_gen.h"
#include "../typ_defs/typedefs_cp.h"
#include "../typ_defs/typedefs_class.h"
#include "../proto_defs/proto_vel_sampl_cp_local.h"
#include "../proto_defs/proto_math.h"

/*===================================================================*/
/*ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*===================================================================*/
/* Coefficient Velocities */
/*===================================================================*/

void sampl_vc(CPOPTS *cpopts,
              CPCOEFFS_INFO *cpcoeffs_info,CPCOEFFS_POS *cpcoeffs_pos,
              STATEPOINT *statepoint,int *iseed,int *iseed2,double *qseed,
              COMMUNICATE *communicate,CP_COMM_STATE_PKG *cp_comm_state_pkg_up,
              CP_COMM_STATE_PKG *cp_comm_state_pkg_dn)

/*===================================================================*/
{/*begin routine*/
/* =================================================================*/
/*               Local variable declarations                        */

   double width;
   int i,j,ip,nvelc,iii,ioff,icoef;

/*  Local pointers */

   int cp_lsda    = cpopts->cp_lsda;
   int pi_beads   = cpcoeffs_info->pi_beads_proc; 
   int ncoef      = cpcoeffs_info->ncoef;
   int ncoef_use  = cpcoeffs_info->ncoef;
   int *ioff_upt  = cpcoeffs_info->ioff_upt;
   int *ioff_dnt  = cpcoeffs_info->ioff_dnt;
   int nstate_up  = cpcoeffs_info->nstate_up;
   int nstate_dn  = cpcoeffs_info->nstate_dn;
   int myid       = cp_comm_state_pkg_up->myid;
   int np_states  = cp_comm_state_pkg_up->num_proc;
   double te_ext  = cpopts->te_ext;
   double *cmass  = cpcoeffs_info->cmass;
   int icmoff_up  = cpcoeffs_info->icoef_start_up-1;
   int icmoff_dn  = cpcoeffs_info->icoef_start_dn-1;
   double *vcre_up,*vcim_up;
   double *vcre_dn,*vcim_dn;

/*==================================================================*/
/* I) Sample the up state when ready, gridley                       */


   if(np_states>1){
    for(ip=1;ip<=pi_beads;ip++){
      if(cpcoeffs_pos[ip].ivcoef_form_up!=1){
        printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
        printf("Up Coef velocities are not in transposed form \n");
        printf("on state processor %d in samp_vel_cp \n",myid);
        printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
        fflush(stdout);
        exit(1);
      }/*endif*/
     }/*endfor*/
   }/*endif*/

   if(np_states>1){
     ncoef    = cp_comm_state_pkg_up->nstate_ncoef_proc_max;
     ncoef_use= cp_comm_state_pkg_up->nstate_ncoef_proc;
   }/*endif*/

   nvelc = ncoef*nstate_up;
   if(cpopts->zero_cp_vel < 1){
     for(ip=1;ip<=pi_beads;ip++){
       gaussran(nvelc,iseed,iseed2,qseed,cpcoeffs_pos[ip].vcre_up);
       gaussran(nvelc,iseed,iseed2,qseed,cpcoeffs_pos[ip].vcim_up);
     }/*endfor*/

     for(ip=1;ip<=pi_beads;ip++){

       vcre_up = cpcoeffs_pos[ip].vcre_up;
       vcim_up = cpcoeffs_pos[ip].vcim_up;
       for(i=1;i<=nstate_up;i++){
         for(j=1;j<=ncoef_use;j++){
           icoef = ioff_upt[i]+j;
           width = sqrt(te_ext/(cmass[(j+icmoff_up)]*BOLTZ));
           vcre_up[icoef] *= width;
           vcim_up[icoef] *= width;
         }/*endfor:ncoef*/
         if(myid+1==np_states){
          icoef = ncoef_use+ioff_upt[i];
          vcim_up[icoef] = 0.0;
         }/*endif*/
         for(j=ncoef_use+1+ioff_upt[i];j<=ncoef+ioff_upt[i];j++){
           vcre_up[j] = 0.0;
           vcim_up[j] = 0.0;
         }/*endfor*/
       }/*endfor:nstate_up*/

     }/*endfor:pi_beads*/

   } else {

     for(ip=1;ip<=pi_beads;ip++){

       vcre_up = cpcoeffs_pos[ip].vcre_up;
       vcim_up = cpcoeffs_pos[ip].vcim_up;
       for(i=1;i<=nstate_up;i++){
         for(j=1;j<=ncoef_use;j++){
           icoef = ioff_upt[i]+j;
           vcre_up[icoef] = 0.0;
           vcim_up[icoef] = 0.0;
         }/*endfor:ncoef*/
         if(myid+1==np_states){
          icoef = ncoef_use+ioff_upt[i];
          vcim_up[icoef] = 0.0;
         }/*endif*/
         for(j=ncoef_use+1+ioff_upt[i];j<=ncoef+ioff_upt[i];j++){
           vcre_up[j] = 0.0;
         vcim_up[j] = 0.0;
         }/*endfor*/
       }/*endfor:nstate_up*/

     }/*endfor:pi_beads*/
   }/* endif zero_cp_vel */
/*==================================================================*/
/* II) Sample the dn state when ready, gridley                       */

 if(cp_lsda==1){
   if(np_states>1){
    for(ip=1;ip<=pi_beads;ip++){
      if(cpcoeffs_pos[ip].ivcoef_form_dn!=1){
        printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
        printf("Dn Coef velocities are not in transposed form \n");
        printf("on state processor %d in samp_vel_cp \n",myid);
        printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
        fflush(stdout);
        exit(1);
      }/*endif*/
     }/*endfor*/
   }/*endif*/

    if(np_states>1){
      ncoef    = cp_comm_state_pkg_dn->nstate_ncoef_proc_max;
      ncoef_use= cp_comm_state_pkg_dn->nstate_ncoef_proc;
    }/*endif*/

    nvelc = ncoef*nstate_dn;
    if(cpopts->zero_cp_vel < 1){
      for(ip=1;ip<=pi_beads;ip++){
        gaussran(nvelc,iseed,iseed2,qseed,cpcoeffs_pos[ip].vcre_dn);
        gaussran(nvelc,iseed,iseed2,qseed,cpcoeffs_pos[ip].vcim_dn);
      }/*endfor*/

     for(ip=1;ip<=pi_beads;ip++){
  
       vcre_dn = cpcoeffs_pos[ip].vcre_dn;
       vcim_dn = cpcoeffs_pos[ip].vcim_dn;
       for(i=1;i<=nstate_dn;i++){
         for(j=1;j<=ncoef_use;j++){
           icoef = ioff_dnt[i]+j;
           width = sqrt(te_ext/(cmass[(j+icmoff_dn)]*BOLTZ));
           vcre_dn[icoef] *= width;
           vcim_dn[icoef] *= width;
         }/*endfor:ncoef*/
         if(myid+1==np_states){
          icoef = ncoef_use+ioff_dnt[i];
          vcim_dn[icoef] = 0.0;
         }/*endif*/
         for(j=ncoef_use+1+ioff_dnt[i];j<=ncoef+ioff_dnt[i];j++){
           vcre_dn[j] = 0.0;
           vcim_dn[j] = 0.0;
         }/*endfor*/
       }/*endfor:nstate_dn*/

    }/*endfor:pi_beads*/

  } else {


     for(ip=1;ip<=pi_beads;ip++){
  
       vcre_dn = cpcoeffs_pos[ip].vcre_dn;
       vcim_dn = cpcoeffs_pos[ip].vcim_dn;
       for(i=1;i<=nstate_dn;i++){
         for(j=1;j<=ncoef_use;j++){
           icoef = ioff_dnt[i]+j;
           vcre_dn[icoef] = 0.0;
           vcim_dn[icoef] = 0.0;
         }/*endfor:ncoef*/
         if(myid+1==np_states){
          icoef = ncoef_use+ioff_dnt[i];
          vcim_dn[icoef] = 0.0;
         }/*endif*/
         for(j=ncoef_use+1+ioff_dnt[i];j<=ncoef+ioff_dnt[i];j++){
           vcre_dn[j] = 0.0;
           vcim_dn[j] = 0.0;
         }/*endfor*/
       }/*endfor:nstate_dn*/

    }/*endfor:pi_beads*/
  }/* endif zero_cp_vel */

 }/*endif lsda*/
  

/*-------------------------------------------------------------------*/
  }/*end routine*/
/*===================================================================*/



/*===================================================================*/
/*ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*===================================================================*/
/* Coefficient NHC Velocities */
/*===================================================================*/

void sampl_vcnhc(CPTHERM_INFO *cptherm_info,CPTHERM_POS *cptherm_pos,
                 CPSCR *cpscr,CPOPTS *cpopts, int pi_beads,
                 int *iseed,int *iseed2,double *qseed, 
                 CPCOEFFS_INFO *cpcoeffs_info,int np_states)

/*===================================================================*/
{/*begin routine*/
/* ===================================================================*/
/*             Local variable declarations                            */

 double width;
 int i,j,ip,iii;
 int ioff_im_up,ioff_re_dn,ioff_im_dn;
 int ncoef_up_tot,ncoef_dn_tot;
 int ncoef_up,ncoef_dn;

 double **vc_nhc;
 double **cmass_nhc = cptherm_info->cmass_nhc;
 double *coef_kin = cpscr->cpscr_therm.coef_kin;
 double *cre_up   = cpscr->cpscr_wave.cre_up;
 double *cre_dn   = cpscr->cpscr_wave.cre_dn;

 double cmass_nhc_massiv = cptherm_info->cmass_nhc_massiv;
 int cp_lsda     = cpopts->cp_lsda;
 int massiv_flag = cptherm_info->massiv_flag;
 double te_ext   = cpopts->te_ext;
 int ncoef       = cpcoeffs_info->ncoef;
 int nstate_up   = cpcoeffs_info->nstate_up;
 int nstate_dn   = cpcoeffs_info->nstate_dn;
 int num_c_nhc   = cptherm_info->num_c_nhc_proc;
 int len_c_nhc   = cptherm_info->len_c_nhc;
 int nstate_ncoef_proc_up = cpcoeffs_info->nstate_ncoef_proc_up;
 int nstate_ncoef_proc_dn = cpcoeffs_info->nstate_ncoef_proc_dn;
 double heat_fact   = cptherm_info->cp_therm_heat_fact;

/*====================================================================*/
/* I) Sample away                                                       */

 if(massiv_flag==0){

  for(j=1;j<=len_c_nhc;j++){
   for(ip=1;ip<=pi_beads;ip++){
    gaussran(num_c_nhc,iseed,iseed2,qseed,coef_kin);
    vc_nhc = cptherm_pos[ip].vc_nhc;
    for(i=1;i<=num_c_nhc;i++){
     width = sqrt(te_ext*heat_fact/(cmass_nhc[j][i]*BOLTZ));   
     vc_nhc[j][i] = coef_kin[i]*width;
    }/*endfor*/
   }/*endfor*/
  }/*endfor*/

 }/*endif:massiv*/

/*====================================================================*/
/* Sample carefully                                                  */

 if(massiv_flag==1){

  if(np_states==1) {
    ncoef_up_tot = ncoef*nstate_up;
    ncoef_dn_tot = ncoef*nstate_dn;
  }else{
    ncoef_up_tot = nstate_ncoef_proc_up*nstate_up;
    ncoef_dn_tot = nstate_ncoef_proc_dn*ncoef*nstate_dn;
  }
  ioff_im_up   = ncoef_up_tot;
  ioff_re_dn   = 2*ncoef_up_tot;
  ioff_im_dn   = 2*ncoef_up_tot+ncoef_dn_tot;

  for(j=1;j<=len_c_nhc;j++){
   for(ip=1;ip<=pi_beads;ip++){
    gaussran(ncoef_up_tot,iseed,iseed2,qseed,cre_up);
    vc_nhc = cptherm_pos[ip].vc_nhc;
    for(i=1;i<=ncoef_up_tot;i++){
     width = sqrt(te_ext/(cmass_nhc_massiv*BOLTZ));   
     vc_nhc[j][i] = cre_up[i]*width;
    }/*endfor*/
    gaussran(ncoef_up_tot,iseed,iseed2,qseed,cre_up);
    for(i=1;i<=ncoef_up_tot;i++){
     width = sqrt(te_ext/(cmass_nhc_massiv*BOLTZ));   
     vc_nhc[j][(i+ioff_im_up)] = cre_up[i]*width;
    }/*endfor*/

    if((cp_lsda==1) && (nstate_dn > 0) ){
     gaussran(ncoef_dn_tot,iseed,iseed2,qseed,cre_dn);
     for(i=1;i<=ncoef_up_tot;i++){
      width = sqrt(te_ext/(cmass_nhc_massiv*BOLTZ));   
      vc_nhc[j][(i+ioff_re_dn)] = cre_dn[i]*width;
     }/*endfor*/
     gaussran(ncoef_dn_tot,iseed,iseed2,qseed,cre_dn);
     for(i=1;i<=ncoef_dn_tot;i++){
      width = sqrt(te_ext/(cmass_nhc_massiv*BOLTZ));   
      vc_nhc[j][(i+ioff_im_dn)] = cre_dn[i]*width;
     }/*endfor*/
    }/*endif*/

    cptherm_pos[ip].itherm_form_up = 1;
    cptherm_pos[ip].itherm_form_dn = 1;

   }/*endfor:ip*/
  }/*endfor:len*/

 }/*endif:massiv*/

/*-------------------------------------------------------------------*/
}/*end routine*/
/*-------------------------------------------------------------------*/






