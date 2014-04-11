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
#include "../proto_defs/proto_energy_cp_entry.h"
#include "../proto_defs/proto_energy_cp_local.h"
#include "../proto_defs/proto_communicate_wrappers.h"

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void min_STD_cp(CLASS *class,BONDED *bonded,GENERAL_DATA *general_data,
                CP *cp,int ip_now)

/*========================================================================*/
{/*begin routine*/
/*========================================================================*/
/*             Local variable declarations                                */
#include "../typ_defs/typ_mask.h"


    int i,ipart,icoef,is;
    double dt;
    int iii;
    int ncoef_up_tot,ncoef_dn_tot;


    double *cp_cpcoeffs_cre_up  = cp->cpcoeffs_pos[ip_now].cre_up;
    double *cp_cpcoeffs_cim_up  = cp->cpcoeffs_pos[ip_now].cim_up;
    double *cp_cpcoeffs_cre_dn  = cp->cpcoeffs_pos[ip_now].cre_dn;
    double *cp_cpcoeffs_cim_dn  = cp->cpcoeffs_pos[ip_now].cim_dn;
    double *cp_cpcoeffs_fcre_up  = cp->cpcoeffs_pos[ip_now].fcre_up;
    double *cp_cpcoeffs_fcim_up  = cp->cpcoeffs_pos[ip_now].fcim_up;
    double *cp_cpcoeffs_fcre_dn  = cp->cpcoeffs_pos[ip_now].fcre_dn;
    double *cp_cpcoeffs_fcim_dn  = cp->cpcoeffs_pos[ip_now].fcim_dn;
    double *cp_cpcoeffs_cmass    = cp->cpcoeffs_info.cmass;
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
    int *cpcoeffs_ioff_up       = cp->cpcoeffs_info.ioff_upt;
    int *cpcoeffs_ioff_dn       = cp->cpcoeffs_info.ioff_dnt;

    double *class_clatoms_fx  = class->clatoms_pos[ip_now].fx;
    double *class_clatoms_fy  = class->clatoms_pos[ip_now].fy;
    double *class_clatoms_fz  = class->clatoms_pos[ip_now].fz;

    int nstate_up = cp->cpcoeffs_info.nstate_up;
    int nstate_dn = cp->cpcoeffs_info.nstate_dn;
    int ncoef     = cp->cpcoeffs_info.ncoef;
    int natm_tot  = class->clatoms_info.natm_tot;

    int np_states            = cp->communicate.np_states;
    int myid_state           = cp->communicate.myid_state;
    int cp_norb               = cp->cpopts.cp_norb;
    int cp_lsda               =  cp->cpopts.cp_lsda;
    int ncoef_up,ncoef_dn,ncoef_up_max,ncoef_dn_max;

/*==========================================================================*/
/* 0) Checks                                                                */

  if(cp_norb>0){
    if((icoef_orth_up+ifcoef_orth_up)!=0){
      printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
      printf("Up CP vectors are in orthonormal form \n");
      printf("on state processor %d in min_STD_cp \n",myid_state);
      printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
      fflush(stdout);
      exit(1);
    }/*endif*/
    if(cp_lsda==1){
     if((icoef_orth_dn+ifcoef_orth_dn)!=0){
       printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
       printf("Dn CP vectors are in orthonormal form \n");
       printf("on state processor %d in min_STD_cp \n",myid_state);
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
     printf("on state processor %d in min_STD_cp \n",myid_state);
     printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
     fflush(stdout);
     exit(1);
    }/*endif*/
    if(cp_lsda==1){
     if((icoef_form_dn+ifcoef_form_dn)!=2){
      printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
      printf("Up CP vectors are not in transposed form \n");
      printf("on state processor %d in min_STD_cp \n",myid_state);
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
    dt  = general_data->timeinfo.dt;
    general_data->stat_avg.iter_shake     = 0; 
    general_data->stat_avg.iter_ratl     = 0; 

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
                   

/*==========================================================================*/
/* I) Get forces                                                            */

   for(i=1;i<= natm_tot;i++){
     class_clatoms_fx[i]  = 0.0;
     class_clatoms_fy[i]  = 0.0;
     class_clatoms_fz[i]  = 0.0;
   }/*endfor*/

   for(i=1;i<=ncoef_up_tot;i++){
     cp_cpcoeffs_fcre_up[i] = 0.0;
     cp_cpcoeffs_fcim_up[i] = 0.0;
   }/*endfor*/
   if( (cp->cpopts.cp_lsda == 1) && (nstate_dn != 0) ){
     for(i=1;i<=ncoef_dn_tot; i++){
      cp_cpcoeffs_fcre_dn[i] = 0.0;
      cp_cpcoeffs_fcim_dn[i] = 0.0;
     }/*endfor*/
   }/*endif*/

   for(i=1;i<=9;i++){
     general_data->ptens.pvten[i]     = 0.0;
     general_data->ptens.pvten_tot[i] = 0.0;
   }/*endfor*/

   general_data->stat_avg.cp_ehart = 0.0;
   general_data->stat_avg.cp_exc   = 0.0;
   general_data->stat_avg.cp_muxc  = 0.0;
   general_data->stat_avg.cp_eext  = 0.0;
   general_data->stat_avg.cp_enl   = 0.0;
   general_data->stat_avg.cp_eke   = 0.0;
   general_data->stat_avg.vrecip   = 0.0;
   cp_ks_energy_ctrl(cp,ip_now,&(general_data->ewald),&(class->ewd_scr),
                     &(general_data->cell),
                     &(class->clatoms_info),
                     &(class->clatoms_pos[ip_now]),
                     &(class->atommaps),&(general_data->stat_avg),
                     &(general_data->ptens),
                     &(general_data->simopts),
                     &(class->for_scr));

/*==========================================================================*/
/* IV) Evolve positions and coefficients                                   */ 

   for(is=1;is<=nstate_up;is++) {
     for(i=1;i<=ncoef_up;i++) {
       icoef = i+cpcoeffs_ioff_up[is];
       cp_cpcoeffs_cre_up[icoef] += 
          dt*cp_cpcoeffs_fcre_up[icoef]/cp_cpcoeffs_cmass[(i+icmoff_up)];
       cp_cpcoeffs_cim_up[icoef] += 
          dt*cp_cpcoeffs_fcim_up[icoef]/cp_cpcoeffs_cmass[(i+icmoff_up)];
     }/*endfor*/
   }/*endfor*/

   if( (cp->cpopts.cp_lsda == 1) && (nstate_dn != 0) ){
    for(is=1;is<=nstate_dn;is++) {
      for(i=1;i<=ncoef_dn;i++) {
       icoef = i+cpcoeffs_ioff_dn[is];
       cp_cpcoeffs_cre_dn[icoef] += 
          dt*cp_cpcoeffs_fcre_dn[icoef]/cp_cpcoeffs_cmass[(i+icmoff_dn)];
       cp_cpcoeffs_cim_dn[icoef] += 
          dt*cp_cpcoeffs_fcim_dn[icoef]/cp_cpcoeffs_cmass[(i+icmoff_dn)];
      }/*endfor*/
    }/*endfor*/
   }/* endif */

/*==========================================================================*/
/* IV.V) Orthogonalize Wave functions                                        */

  cp_shuffle_states(cp,ip_now);
  orthog_control_cp(cp,ip_now);

/*-------------------------------------------------------------------------*/
    }/*end routine*/
/*==========================================================================*/


