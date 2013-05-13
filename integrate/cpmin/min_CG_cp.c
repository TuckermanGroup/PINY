/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
/*                                                                          */
/*                         PI_MD:                                           */
/*             The future of simulation technology                          */
/*             ------------------------------------                         */
/*                     Module: int_NVE                                      */
/*                                                                          */
/* This subprogram integrates the system using Vel Verlet                   */
/*                                                                          */
/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/


#include "standard_include.h"
#include "../typ_defs/typedefs_gen.h"
#include "../typ_defs/typedefs_class.h"
#include "../typ_defs/typedefs_bnd.h"
#include "../typ_defs/typedefs_cp.h"
#include "../proto_defs/proto_integrate_cpmin_entry.h"
#include "../proto_defs/proto_integrate_cpmin_local.h"
#include "../proto_defs/proto_energy_cpcon_entry.h"
#include "../proto_defs/proto_intra_con_entry.h"
#include "../proto_defs/proto_energy_ctrl_cp_entry.h"
#include "../proto_defs/proto_energy_ctrl_local.h"
#include "../proto_defs/proto_energy_cp_local.h"
#include "../proto_defs/proto_friend_lib_entry.h"
#include "../proto_defs/proto_communicate_wrappers.h"
#include "../proto_defs/proto_math.h"

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
void min_CG_cp(CLASS *class,BONDED *bonded,GENERAL_DATA *general_data,CP *cp,
               int ip_now)
/*========================================================================*/
{/*begin routine*/
/*========================================================================*/
/*             Local variable declarations                                */
#include "../typ_defs/typ_mask.h"

    int i,ipart,icoef,is;
    double dt;
    int iii,iv,ioff_hyb;
    int ncoef_up_tot,ncoef_dn_tot;
    double gamma_up,gamma_dn,wght1,wght2;
    double fc_mag_up,fc_mag_dn;
    static double fovlap_up,fovlap_dn;
    static double fovlap_tmp;
    double fovlap_up_old,fovlap_dn_old;
    double *zeta_up,*zeta_dn;

    int natm_tot  = class->clatoms_info.natm_tot;
    double *class_clatoms_fx   = class->clatoms_pos[ip_now].fx;
    double *class_clatoms_fy   = class->clatoms_pos[ip_now].fy;
    double *class_clatoms_fz   = class->clatoms_pos[ip_now].fz;
    double *class_clatoms_mass = class->clatoms_info.mass;
    double *ptens_pvten     = general_data->ptens.pvten;
    double *ptens_pvten_tot =  general_data->ptens.pvten_tot;

    double *cp_hess_re_up = cp->cpcoeffs_pos[ip_now].cp_hess_re_up;
    double *cp_hess_im_up = cp->cpcoeffs_pos[ip_now].cp_hess_im_up;
    double *cp_hess_re_dn = cp->cpcoeffs_pos[ip_now].cp_hess_re_dn;
    double *cp_hess_im_dn = cp->cpcoeffs_pos[ip_now].cp_hess_im_dn;

    double *cre_up  = cp->cpcoeffs_pos[ip_now].cre_up;
    double *cim_up  = cp->cpcoeffs_pos[ip_now].cim_up;
    double *cre_dn  = cp->cpcoeffs_pos[ip_now].cre_dn;
    double *cim_dn  = cp->cpcoeffs_pos[ip_now].cim_dn;
    double *fcre_up  = cp->cpcoeffs_pos[ip_now].fcre_up;
    double *fcim_up  = cp->cpcoeffs_pos[ip_now].fcim_up;
    double *fcre_dn  = cp->cpcoeffs_pos[ip_now].fcre_dn;
    double *fcim_dn  = cp->cpcoeffs_pos[ip_now].fcim_dn;
    double *hcre_up  = cp->cpcoeffs_pos[ip_now].vcre_up;
    double *hcim_up  = cp->cpcoeffs_pos[ip_now].vcim_up;
    double *hcre_dn  = cp->cpcoeffs_pos[ip_now].vcre_dn;
    double *hcim_dn  = cp->cpcoeffs_pos[ip_now].vcim_dn;
    double *cmass    = cp->cpcoeffs_info.cmass;
    int icmoff_up        = cp->cpcoeffs_info.icoef_start_up-1;
    int icmoff_dn        = cp->cpcoeffs_info.icoef_start_dn-1;
    int icoef_form_up        = cp->cpcoeffs_pos[ip_now].icoef_form_up;
    int icoef_orth_up        = cp->cpcoeffs_pos[ip_now].icoef_orth_up;
    int ifcoef_form_up        = cp->cpcoeffs_pos[ip_now].ifcoef_form_up;
    int ifcoef_orth_up        = cp->cpcoeffs_pos[ip_now].ifcoef_orth_up;
    int icoef_form_dn        = cp->cpcoeffs_pos[ip_now].icoef_form_dn;
    int icoef_orth_dn        = cp->cpcoeffs_pos[ip_now].icoef_orth_dn;
    int ifcoef_form_dn        = cp->cpcoeffs_pos[ip_now].ifcoef_form_dn;
    int ifcoef_orth_dn        = cp->cpcoeffs_pos[ip_now].ifcoef_orth_dn;
    int *ioff_up       = cp->cpcoeffs_info.ioff_upt;
    int *ioff_dn       = cp->cpcoeffs_info.ioff_dnt;

    int nstate_up = cp->cpcoeffs_info.nstate_up;
    int nstate_dn = cp->cpcoeffs_info.nstate_dn;
    int ncoef     = cp->cpcoeffs_info.ncoef;
    int cp_lsda   = cp->cpopts.cp_lsda;
    int cg_reset_flag = cp->cpcoeffs_info.cg_reset_flag;
    int cp_norb =  cp->cpopts.cp_norb;

    int np_states        = class->communicate.np_states;
    MPI_Comm comm_states = class->communicate.comm_states;
    int myid_state           = cp->communicate.myid_state;
    int ncoef_up,ncoef_dn,ncoef_up_max,ncoef_dn_max;
    int icoef_off_up,icoef_off_dn;
 
    int cp_cg_line_min_len = general_data->minopts.cp_cg_line_min_len;
    double eenergy,eenergy_temp;
    double gamma=1.0;
    double sum_check = 0.0,sum_check_tmp=0.0;
    int cp_para_opt = cp->cpopts.cp_para_opt;

/*==========================================================================*/
/* 0) Checks                                                                */

  if(cp_norb>0){
    if((icoef_orth_up+ifcoef_orth_up)!=0){
      printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
      printf("Up CP vectors are in orthonormal form \n");
      printf("on state processor %d in min_CG_cp \n",myid_state);
      printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
      fflush(stdout);
      exit(1);
    }/*endif*/
    if(cp_lsda==1){
     if((icoef_orth_dn+ifcoef_orth_dn)!=0){
       printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
       printf("Dn CP vectors are in orthonormal form \n");
       printf("on state processor %d in min_CG_cp \n",myid_state);
       printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
       fflush(stdout);
       exit(1);
     }/*endif*/
    }/*endif*/
  }/*endif*/

  if(np_states>1){
    if((icoef_form_up+ifcoef_form_up)!=2){
     printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
     printf("Up CP vectors are not in transposed form \n");
     printf("on state processor %d in min_CG_cp \n",myid_state);
     printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
     fflush(stdout);
     exit(1);
    }/*endif*/
    if(cp_lsda==1){
     if((icoef_form_dn+ifcoef_form_dn)!=2){
      printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
      printf("Up CP vectors are not in transposed form \n");
      printf("on state processor %d in min_CG_cp \n",myid_state);
      printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
      fflush(stdout);
      exit(1);
     }/*endif*/
    }/*endif*/
  }/*endif*/


/*==========================================================================*/
/* 0) Useful constants                                                      */
    general_data->timeinfo.int_res_tra    = 0;
    general_data->timeinfo.int_res_ter    = 0;
    dt                              = general_data->timeinfo.dt;
    general_data->stat_avg.iter_shake     = 0; 
    general_data->stat_avg.iter_ratl      = 0; 

    if(np_states==1){
     ncoef_up      = cp->cpcoeffs_info.ncoef;
     ncoef_dn      = cp->cpcoeffs_info.ncoef;
     ncoef_up_max  = cp->cpcoeffs_info.ncoef;
     ncoef_dn_max  = cp->cpcoeffs_info.ncoef;
     /*RLH Add icoef_off_up, icoef_off_dn */
     icoef_off_up = cp->cpcoeffs_info.icoef_start_up-1;
     icoef_off_dn = cp->cpcoeffs_info.icoef_start_dn-1;
    }else{
     ncoef_up     = cp->cp_comm_state_pkg_up.nstate_ncoef_proc;
     ncoef_dn     = cp->cp_comm_state_pkg_dn.nstate_ncoef_proc;
     ncoef_up_max = cp->cp_comm_state_pkg_up.nstate_ncoef_proc_max;
     ncoef_dn_max = cp->cp_comm_state_pkg_dn.nstate_ncoef_proc_max;
     icoef_off_up=cp->cp_comm_state_pkg_up.icoef_start-1;
     icoef_off_dn=cp->cp_comm_state_pkg_dn.icoef_start-1;
    }/*endif*/

    ncoef_up_tot = ncoef_up_max*nstate_up;
    ncoef_dn_tot = ncoef_dn_max*nstate_dn; 

    zeta_up = (double *) cmalloc(ncoef*sizeof(double))-1;
    if(cp_lsda==1 && nstate_dn != 0){
      zeta_dn = (double *) cmalloc(ncoef*sizeof(double))-1;
    }/*endif*/


/*==========================================================================*/
/* 0.1) Zero conjugate gradients                                            */

   if(cg_reset_flag == 1){
     for(i=1;i<=ncoef_up_tot; i++){
      hcre_up[i] = 0.0;
      hcim_up[i] = 0.0;
     }
     if( (cp_lsda== 1) && (nstate_dn != 0) ){
      for(i=1;i<=ncoef_dn_tot; i++){
       hcre_dn[i] = 0.0;
       hcim_dn[i] = 0.0;
      }
     }/* endif */
     gamma_up = 0.0;
     gamma_dn = 0.0;
     fovlap_up = 1.0;
     fovlap_dn = 1.0;
   }/* endif */

/*==========================================================================*/
/* I) Get forces                                                            */

   for(i=1;i<= natm_tot;i++){
     class_clatoms_fx[i]  = 0.0;
     class_clatoms_fy[i]  = 0.0;
     class_clatoms_fz[i]  = 0.0;
   }/*endfor*/

   for(i=1;i<=ncoef_up_tot;i++){
     fcre_up[i] = 0.0;
     fcim_up[i] = 0.0;
   }/*endfor*/
   if( (cp_lsda== 1) && (nstate_dn != 0) ){
     for(i=1;i<=ncoef_dn_tot; i++){
      fcre_dn[i] = 0.0;
      fcim_dn[i] = 0.0;
     }/*endfor*/
   }/*endif*/
   general_data->stat_avg.cp_ehart = 0.0;
   general_data->stat_avg.cp_exc   = 0.0;
   general_data->stat_avg.cp_muxc  = 0.0;
   general_data->stat_avg.cp_eext  = 0.0;
   general_data->stat_avg.cp_enl   = 0.0;
   general_data->stat_avg.cp_eke   = 0.0;
   general_data->stat_avg.vrecip   = 0.0;

   for(i=1;i<=9;i++){
     ptens_pvten[i]     = 0.0;
     ptens_pvten_tot[i] = 0.0;
   }/*endfor*/

   cp_ks_energy_ctrl(cp,ip_now,&(general_data->ewald),&(class->ewd_scr),
                       &(general_data->cell),
                       &(class->clatoms_info),
                       &(class->clatoms_pos[ip_now]),
                       &(class->atommaps),&(general_data->stat_avg),
                       &(general_data->ptens),
                       &(general_data->simopts),
                       &(class->for_scr));

   get_diag_cp_hess(cp,ip_now,&(general_data->cell),gamma);

    
   eenergy_temp = general_data->stat_avg.cp_ehart
               + general_data->stat_avg.cp_exc + general_data->stat_avg.cp_eext
               + general_data->stat_avg.cp_enl + general_data->stat_avg.cp_eke  ;



   if(np_states>1){
        Allreduce(&(eenergy_temp),&(eenergy),1,MPI_DOUBLE,MPI_SUM,0,comm_states);
   }else{
         eenergy = eenergy_temp;
   }


/*==========================================================================*/
/* II) Calculate the gamma's                                                */

   fovlap_up_old = fovlap_up;
   fovlap_up = 0.0;
   wght1 = 2.0;
   wght2 = 2.0;
   if(myid_state+1==np_states){wght1=1.0;wght2=0.0;}
   for(i=1;i<=nstate_up;i++){
    fovlap_up += diag_ovlap(ncoef_up,wght1,wght2,fcre_up,
                            fcim_up,fcre_up,
                            fcim_up,ioff_up[i]);
   }/*endfor*/
   if(np_states>1){
    fovlap_tmp = fovlap_up; 
    Allreduce(&(fovlap_tmp),&(fovlap_up),1,MPI_DOUBLE,MPI_SUM,0,comm_states);
   }/*endif*/
   if(cg_reset_flag != 1) {gamma_up = fovlap_up/fovlap_up_old;}

   if( (cp_lsda== 1) && (nstate_dn != 0) ){
     wght1 = 2.0;
     wght2 = 2.0;
     if(myid_state+1==np_states){wght1=1.0;wght2=0.0;}
     fovlap_dn_old = fovlap_dn;
     fovlap_dn = 0.0;
     for(i=1;i<=nstate_dn;i++){
      fovlap_dn += diag_ovlap(ncoef_dn,wght1,wght2,fcre_dn,
                             fcim_dn,fcre_dn,
                             fcim_dn,ioff_dn[i]);
     }/* endfor */
     if(np_states>1){
      fovlap_tmp = fovlap_dn; 
      Allreduce(&(fovlap_tmp),&(fovlap_dn),1,MPI_DOUBLE,MPI_SUM,0,comm_states);
     }/*endif*/
     if(cg_reset_flag != 1) {gamma_dn = fovlap_dn/fovlap_dn_old;}
   }/* endif */

/*==========================================================================*/
/* II.V) Evolve gradients                                                   */

   for(i=1;i<=ncoef_up_tot; i++){
     hcre_up[i] = fcre_up[i] + gamma_up*hcre_up[i];
     hcim_up[i] = fcim_up[i] + gamma_up*hcim_up[i];
   }/*endfor*/

   if( (cp_lsda== 1) && (nstate_dn != 0) ){
    for(i=1;i<=ncoef_dn_tot; i++){
      hcre_dn[i] = fcre_dn[i] + gamma_dn*hcre_dn[i];
      hcim_dn[i] = fcim_dn[i] + gamma_dn*hcim_dn[i];
    }/*endfor*/
   }/* endif */

/*==========================================================================*/
/* III) Calculate the step length                                           */
 
  /* RLH Fix offset for hybrid parallelization */
  /* ioff_hyb = (cp_para_opt == 0 ? myid_state*ncoef_up_max : 0);*/
  ioff_hyb = (cp_para_opt == 0 ? icoef_off_up : 0);
  if(cp_cg_line_min_len == 0){
    for(i=1;i<=ncoef_up;i++) {
      zeta_up[i] = 1.0/cp_hess_re_up[ioff_hyb + i];
    }/* endfor */
    if( (cp_lsda== 1) && (nstate_dn != 0) ){
      /*RLH Fix offset for hybrid parallelization*/
      /*ioff_hyb = (cp_para_opt == 0 ? myid_state*ncoef_dn_max : 0);*/
      ioff_hyb = (cp_para_opt == 0 ? icoef_off_dn : 0);
      for(i=1;i<=ncoef_dn;i++) {
        zeta_dn[i] = 1.0/cp_hess_re_dn[ioff_hyb + i];
      }/* endfor */
    }/* endif */
 }else{
     line_min_cp(class,bonded,general_data,cp,zeta_up,zeta_dn,ip_now,eenergy);
 }/*endif*/


/*==========================================================================*/


/*==========================================================================*/
/* IV) Evolve positions and coefficients                                   */ 

   for(is=1;is<=nstate_up;is++) {
     for(i=1;i<=ncoef_up;i++) {
       icoef = i+ioff_up[is];
       cre_up[icoef] +=zeta_up[i]*hcre_up[icoef];
       cim_up[icoef] +=zeta_up[i]*hcim_up[icoef];
     }/*endfor*/
   }/*endfor*/

   if( (cp_lsda== 1) && (nstate_dn != 0) ){
     for(is=1;is<=nstate_dn;is++) {
       for(i=1;i<=ncoef_dn;i++) {
         icoef = i+ioff_dn[is];
         cre_dn[icoef] += zeta_dn[i]*hcre_dn[icoef];
         cim_dn[icoef] += zeta_dn[i]*hcim_dn[icoef];
       }/*endfor*/
     }/*endfor*/
   }/* endif */

/*==========================================================================*/
/* V) Orthogonalize wave functions                                          */
  
   orthog_control_cp(cp,ip_now);
 
/*==========================================================================*/
/* ii) Free the memory                                                      */

   cfree(&(zeta_up[1]));
   if( (cp_lsda== 1) && (nstate_dn != 0) ){
    cfree(&(zeta_dn[1]));
   }

/*-------------------------------------------------------------------------*/
    }/*end routine*/
/*==========================================================================*/



/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
 
double diag_ovlap(int ncomp,double wght1,double wght2,double *v1_re,double *v1_im,
                  double *v2_re,double *v2_im,int ioff)
/*========================================================================*/
{/*begin routine */
/*========================================================================*/

 int ig,i;
 double ovlap;

 ovlap = 0.0;
 for(ig=1;ig<=ncomp-1;ig++){
  i = ig+ioff;
  ovlap += 2.0*(v1_re[i]*v2_re[i] + v1_im[i]*v2_im[i]);
      }
  i = ncomp+ioff;
  ovlap += wght1*v1_re[i]*v2_re[i] + wght2*v1_im[i]*v2_im[i];
  return ovlap;

/*========================================================================*/
}/* end routine */
/*========================================================================*/




/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void line_min_cp(CLASS *class,BONDED *bonded,GENERAL_DATA *general_data,CP *cp,
                 double *zeta_up,double *zeta_dn,int ip_now,double eenergy)

  
/*========================================================================*/
{/*begin routine */
/*========================================================================*/
#include "../typ_defs/typ_mask.h"
    int k,j,flag;
    double *guess,*etot;
    double energy_temp,energy,dum;

    int i,ip=1,ipart,icoef,is;
    double dt;
    int iii,iv;
    int ncoef_up,ncoef_up_tot,ncoef_up_max;
    int ncoef_dn,ncoef_dn_tot,ncoef_dn_max;
    int maxcoef;

/*========================================================================*/
/* Assign local pointers                                                  */
    double *class_clatoms_fx   = class->clatoms_pos[ip_now].fx;
    double *class_clatoms_fy   = class->clatoms_pos[ip_now].fy;
    double *class_clatoms_fz   = class->clatoms_pos[ip_now].fz;
    double *class_clatoms_mass = class->clatoms_info.mass;
    double *ptens_pvten        = general_data->ptens.pvten;
    double *ptens_pvten_tot    =  general_data->ptens.pvten_tot;


    double *cre_up   = cp->cpcoeffs_pos[ip_now].cre_up;
    double *cim_up   = cp->cpcoeffs_pos[ip_now].cim_up;
    double *cre_dn   = cp->cpcoeffs_pos[ip_now].cre_dn;
    double *cim_dn   = cp->cpcoeffs_pos[ip_now].cim_dn;
    double *fcre_up  = cp->cpcoeffs_pos[ip_now].fcre_up;
    double *fcim_up  = cp->cpcoeffs_pos[ip_now].fcim_up;
    double *fcre_dn  = cp->cpcoeffs_pos[ip_now].fcre_dn;
    double *fcim_dn  = cp->cpcoeffs_pos[ip_now].fcim_dn;

    double *hcre_up  = cp->cpcoeffs_pos[ip_now].vcre_up;
    double *hcim_up  = cp->cpcoeffs_pos[ip_now].vcim_up;
    double *hcre_dn  = cp->cpcoeffs_pos[ip_now].vcre_dn;
    double *hcim_dn  = cp->cpcoeffs_pos[ip_now].vcim_dn;
    double *cmass    = cp->cpcoeffs_info.cmass;

    int icmoff_up        = cp->cpcoeffs_info.icoef_start_up-1;
    int icmoff_dn        = cp->cpcoeffs_info.icoef_start_dn-1;
    int icoef_form_up        = cp->cpcoeffs_pos[ip_now].icoef_form_up;
    int icoef_orth_up        = cp->cpcoeffs_pos[ip_now].icoef_orth_up;
    int ifcoef_form_up        = cp->cpcoeffs_pos[ip_now].ifcoef_form_up;
    int ifcoef_orth_up        = cp->cpcoeffs_pos[ip_now].ifcoef_orth_up;
    int icoef_form_dn        = cp->cpcoeffs_pos[ip_now].icoef_form_dn;
    int icoef_orth_dn        = cp->cpcoeffs_pos[ip_now].icoef_orth_dn;
    int ifcoef_form_dn        = cp->cpcoeffs_pos[ip_now].ifcoef_form_dn;
    int ifcoef_orth_dn        = cp->cpcoeffs_pos[ip_now].ifcoef_orth_dn; 

    int *ioff_up       = cp->cpcoeffs_info.ioff_upt;
    int *ioff_dn       = cp->cpcoeffs_info.ioff_dnt;

    int np_states        = class->communicate.np_states;
    MPI_Comm comm_states = class->communicate.comm_states;
    int myid_state           = cp->communicate.myid_state;

    int natm_tot  = class->clatoms_info.natm_tot;

    int cp_lsda        = cp->cpopts.cp_lsda;

    int nstate_up = cp->cpcoeffs_info.nstate_up;
    int nstate_dn = cp->cpcoeffs_info.nstate_dn;
    int ncoef     = cp->cpcoeffs_info.ncoef;


/*========================================================================*/
/* Assign local variables and arrays                                      */

    double *cre_up_old;
    double *cim_up_old;
    double *cre_dn_old;
    double *cim_dn_old;


    double *A,*Ainv; 
    double *C,*B;    /*vector of coefficients */
    double dtmin;
    double min;
    int    imin;
    double deth;
    int iperd = 3;
    int  nguess= general_data->minopts.cp_cg_line_min_len ;
                      /* number of guesses for dt used */

/*==========================================================================*/
  if(np_states>1){
    if((icoef_form_up+ifcoef_form_up)!=2){
     printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
     printf("Up CP vectors are not in transposed form \n");
     printf("on state processor %d in min_CG_cp \n",myid_state);
     printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
     fflush(stdout);
     exit(1);
    }/*endif*/
    if(cp_lsda==1){
     if((icoef_form_dn+ifcoef_form_dn)!=2){
      printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
      printf("Up CP vectors are not in transposed form \n");
      printf("on state processor %d in min_CG_cp \n",myid_state);
      printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
      fflush(stdout);
      exit(1);
     }/*endif*/
    }/*endif*/
  }/*endif*/
/*==========================================================================*/
    if(np_states==1){
     ncoef_up      = cp->cpcoeffs_info.ncoef;
     ncoef_dn      = cp->cpcoeffs_info.ncoef;
     ncoef_up_max  = cp->cpcoeffs_info.ncoef;
     ncoef_dn_max  = cp->cpcoeffs_info.ncoef;
    }else{
     ncoef_up     = cp->cp_comm_state_pkg_up.nstate_ncoef_proc;
     ncoef_dn     = cp->cp_comm_state_pkg_dn.nstate_ncoef_proc;
     ncoef_up_max = cp->cp_comm_state_pkg_up.nstate_ncoef_proc_max;
     ncoef_dn_max = cp->cp_comm_state_pkg_dn.nstate_ncoef_proc_max;
    }/*endif*/

    ncoef_up_tot = ncoef_up_max*nstate_up;
    ncoef_dn_tot = ncoef_dn_max*nstate_dn; 

    dt           = general_data->timeinfo.dt;

/*==========================================================================*/
/* malloc up local arrays                                                   */

    guess  = (double *) cmalloc(nguess*sizeof(double ))-1;
    etot   = (double *) cmalloc((nguess+1)*sizeof(double ))-1;

    A      =  (double *) cmalloc(9*sizeof(double )) -1 ;
    Ainv   =  (double *) cmalloc(9*sizeof(double )) -1 ;


    C      = (double *) cmalloc(3*sizeof(double)) -1;
    B      = (double *) cmalloc(3*sizeof(double)) -1;

   cre_up_old = (double *) cmalloc(ncoef_up_tot*sizeof(double))-1;
   cim_up_old = (double *) cmalloc(ncoef_up_tot*sizeof(double))-1;

 if(cp_lsda == 1){
    cre_dn_old = (double *) cmalloc(ncoef_dn_tot*sizeof(double))-1;
    cim_dn_old = (double *) cmalloc(ncoef_dn_tot*sizeof(double))-1;
 }/*end if*/


/*==========================================================================*/
/* save  coefficients                                       */

   for(is=1;is<=nstate_up;is++) {
     for(i=1;i<=ncoef_up;i++) {
       icoef = i+ioff_up[is];
         cre_up_old[icoef] = cre_up[icoef];
         cim_up_old[icoef] = cim_up[icoef];
     }/*endfor*/
   }/*endfor*/

 if( (cp_lsda == 1) && (nstate_dn != 0) ){
   for(is=1;is<=nstate_dn;is++) {
     for(i=1;i<=ncoef_dn;i++) {
       icoef = i+ioff_dn[is];
         cre_dn_old[icoef] = cre_dn[icoef];
         cim_dn_old[icoef] = cim_dn[icoef];
     }/*endfor*/
   }/*endfor*/
 }/*endif*/


/*==========================================================================*/
/* III) Calculate the step length                                           */

 if(nguess == 3 ){
  guess[1] = dt;
  guess[2] = dt*2.0;
  guess[3] = dt*4.0;
/*  guess[2] = dt*4.0;
  guess[3] = dt*8.0; */
 }

 if(nguess > 3){
    guess[1] = dt;
   for(i=2; i<= nguess; i++){
     guess[i] = dt*1.5*((double)i);
   }/*endfor*/
 }/*endif*/

 for(k=1; k<=nguess ;k++){

   for(i=1;i<=ncoef_up;i++) {
    zeta_up[i] = guess[k]/cmass[(i+icmoff_up)];  /*guess[k] replaced dt */
   }/* endfor */
   if( (cp_lsda == 1) && (nstate_dn != 0) ){
    for(i=1;i<=ncoef_dn;i++) {
     zeta_dn[i] = guess[k]/cmass[(i+icmoff_dn)];
    }/* endfor */
   }/* endif */


/*==========================================================================*/


/*==========================================================================*/
/* IV) Evolve positions and coefficients                                   */ 

   for(is=1;is<=nstate_up;is++) {
     for(i=1;i<=ncoef_up;i++) {
       icoef = i+ioff_up[is];
         cre_up[icoef] = cre_up_old[icoef] + zeta_up[i]*hcre_up[icoef];
         cim_up[icoef] = cim_up_old[icoef] + zeta_up[i]*hcim_up[icoef]; 
     }/*endfor*/
   }/*endfor*/


 if( (cp_lsda == 1) && (nstate_dn != 0) ){
   for(is=1;is<=nstate_dn;is++) {
     for(i=1;i<=ncoef_dn;i++) {
       icoef = i+ioff_dn[is];
         cre_dn[icoef] = cre_dn_old[icoef] + zeta_dn[i]*hcre_dn[icoef];
         cim_dn[icoef] = cim_dn_old[icoef] + zeta_dn[i]*hcim_dn[icoef]; 
     }/*endfor*/
   }/*endfor*/
  }/* endif */


/*========================================================================*/
/* V) Orthogonalize wave functions                                          */
 orthog_control_cp(cp,ip_now);

/*========================================================================*/
/*  for line min */
/*========================================================================*/
/* Calculate the energy for each step made using guess[k]                 */
/*==========================================================================*/
/* I) Get forces                                                            */

   for(i=1;i<= natm_tot;i++){
     class_clatoms_fx[i]  = 0.0;
     class_clatoms_fy[i]  = 0.0;
     class_clatoms_fz[i]  = 0.0;
   }/*endfor*/

   for(i=1;i<=ncoef_up_tot;i++){
     fcre_up[i] = 0.0;
     fcim_up[i] = 0.0;
   }/*endfor*/
   if( (cp_lsda== 1) && (nstate_dn != 0) ){
     for(i=1;i<=ncoef_dn_tot; i++){
      fcre_dn[i] = 0.0;
      fcim_dn[i] = 0.0;
     }/*endfor*/
   }/*endif*/
   general_data->stat_avg.cp_ehart = 0.0;
   general_data->stat_avg.cp_exc   = 0.0;
   general_data->stat_avg.cp_muxc  = 0.0;
   general_data->stat_avg.cp_eext  = 0.0;
   general_data->stat_avg.cp_enl   = 0.0;
   general_data->stat_avg.cp_eke   = 0.0;
   general_data->stat_avg.vrecip   = 0.0;

   for(i=1;i<=9;i++){
     ptens_pvten[i]     = 0.0;
     ptens_pvten_tot[i] = 0.0;
   }/*endfor*/

   cp_ks_energy_ctrl(cp,ip_now,&(general_data->ewald),&(class->ewd_scr),
                       &(general_data->cell),
                       &(class->clatoms_info),
                       &(class->clatoms_pos[ip_now]),
                       &(class->atommaps),&(general_data->stat_avg),
                       &(general_data->ptens),
                       &(general_data->simopts),
                       &(class->for_scr));
    
   energy_temp =  general_data->stat_avg.cp_ehart
               + general_data->stat_avg.cp_exc + general_data->stat_avg.cp_eext
               + general_data->stat_avg.cp_enl + general_data->stat_avg.cp_eke ;

   if(np_states>1){
        Allreduce(&(energy_temp),&(energy),1,MPI_DOUBLE,MPI_SUM,0,comm_states);
   }else{
         energy = energy_temp;
   }

       etot[k] = energy;

}/*endfor nguess*/


/*==========================================================================*/
/*==========================================================================*/

    imin = 0;
    min = 1000.0;

/* Determine if the end point move is the lowest energy                   */
/* If so, set dtmin = guess[1] or guess[nguess]                           */

  for(i=1; i<= nguess; i++){
    if( etot[i] <= min ){
        min = etot[i];
        imin = i;
    }/*endif*/
  }/*endfor*/


 if(imin==1){
   general_data->timeinfo.dt /= 1.20;
/*  printf("scaling the dt /1.20 \n"); */
 }
 if(imin==nguess){
   general_data->timeinfo.dt *= 1.20;
/*    printf("scaling the dt *1.20 \n"); */
 }

 if((imin == 1) || (imin == nguess)  ){
       dtmin = guess[imin];
   for(i=1;i<=ncoef_up;i++) {
    zeta_up[i] = dtmin/cmass[(i+icmoff_up)];  /*guess[k] replaced dt */

   }/* endfor */
   if( (cp->cpopts.cp_lsda == 1) && (nstate_dn != 0) ){
     for(i=1;i<=ncoef_dn;i++) {
      zeta_dn[i] = dtmin/cmass[(i+icmoff_dn)];
     }/* endfor */
 }/* endif lsda */
 }else{
/*------------------------------------------------------------*/
/* nguess is the number of guess (dt) values there are        */
/* compute matrix A(m,k) = sum over nx (guess^(k+m))          */

/*printf("!!!MATRIX PART !!!!! \n");  */
        iii=1;
   for(i=0; i< 3; i++){
       for(j=0; j< 3; j++){
             A[iii] = 0.0;
           for(k=0; k< nguess; k++){
              dum = pow(guess[k+1],(i+j));
             A[iii] += dum;
           }
              iii++;

       }
   }


/*------------------------------------------------------------*/
/* add in the energy for dt=0  */
/* note that this contributes only to the A[1] element  */
           A[1] += eenergy;     

  
    gethinv(A,Ainv,&deth,iperd);

/*------------------------------------------------------------*/
/* compute B lhs of equation */

 
   for(i=1 ; i<=3; i++){
       B[i] = 0.0;
     for(j=1; j<=  nguess ; j++){
        B[i] += (etot[j]*pow(guess[j],(i-1)));
     }/*endfor*/
   }/*endfor*/


/*------------------------------------------------------------*/
/* compute coefficients in C for polynomial c0 + c1*z + c2*z^2  */
/* where z used as short hand */

     C[1] = Ainv[1]*B[1] + Ainv[2]*B[2] + Ainv[3]*B[3]; 
     C[2] = Ainv[4]*B[1] + Ainv[5]*B[2] + Ainv[6]*B[3]; 
     C[3] = Ainv[7]*B[1] + Ainv[8]*B[2] + Ainv[9]*B[3]; 

    dtmin = -C[2]/(2.0*C[3]);

/*------------------------------------------------------------*/

   for(i=1;i<=ncoef_up;i++) {
    zeta_up[i] = dtmin/cmass[(i+icmoff_up)]; 
   }/* endfor */
   if( (cp_lsda == 1) && (nstate_dn != 0) ){
    for(i=1;i<=ncoef_dn;i++) {
     zeta_dn[i] = dtmin/cmass[(i+icmoff_dn)];
    }/* endfor */
 }/* endif lsda */

 }/*endif*/

/*========================================================================*/
/* reassign coefficients                                                  */


   for(is=1;is<=nstate_up;is++) {
     for(i=1;i<=ncoef_up;i++) {
       icoef = i+ioff_up[is];
         cre_up[icoef] = cre_up_old[icoef];
         cim_up[icoef] = cim_up_old[icoef];
     }/*endfor*/
   }/*endfor*/

 if( (cp_lsda == 1) && (nstate_dn != 0) ){
   for(is=1;is<=nstate_dn;is++) {
     for(i=1;i<=ncoef_dn;i++) {
       icoef = i+ioff_dn[is];
         cre_dn[icoef] = cre_dn_old[icoef];
         cim_dn[icoef] = cim_dn_old[icoef];
     }/*endfor*/
   }/*endfor*/
 }/*endif*/


  
/*========================================================================*/
/* free locally assigned memory                                           */

   cfree(&(cre_up_old[1]));
   cfree(&(cim_up_old[1]));
 if(cp_lsda == 1){
   cfree(&(cre_dn_old[1]));
   cfree(&(cim_dn_old[1]));
 }/*end if*/

   cfree(&(guess[1]));
   cfree(&(etot[1]));
   cfree(&(C[1]));
   cfree(&(B[1]));
   cfree(&(A[1]));
   cfree(&(Ainv[1]));

/*========================================================================*/
}/* end routine */
/*========================================================================*/










