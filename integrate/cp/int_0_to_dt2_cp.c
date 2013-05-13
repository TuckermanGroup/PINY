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
#include "../typ_defs/typedefs_cp.h"
#include "../typ_defs/typedefs_class.h"
#include "../typ_defs/typedefs_bnd.h"
#include "../proto_defs/proto_integrate_cp_entry.h"
#include "../proto_defs/proto_integrate_cp_local.h"
#include "../proto_defs/proto_integrate_md_local.h"
#include "../proto_defs/proto_energy_cpcon_entry.h"
#include "../proto_defs/proto_energy_ctrl_cp_entry.h"
#include "../proto_defs/proto_intra_con_entry.h"
#include "../proto_defs/proto_friend_lib_entry.h"
#include "../proto_defs/proto_communicate_wrappers.h"
#include "../proto_defs/proto_coords_local.h"


#define DEBUG_OFF
#define DEBUG_NHC_OFF

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void int_0_to_dt2_cp(CLASS *class,BONDED *bonded,
                     GENERAL_DATA *general_data,CP *cp)

/*========================================================================*/
{/*begin routine*/
/*========================================================================*/
/*             Local variable declarations                                */
#include "../typ_defs/typ_mask.h"

    int i,ip,ipart,icoef,is,js,ifirst=1,iproc;
    double dt,dt2,tol_glob;
    int natm_tot,iii,jjj;
/* PM */ 
/* int nstate_up,nstate_dn; */  
    int ncoef,ncoef_up_tot,ncoef_dn_tot;
    int ncoef_up,ncoef_dn,ncoef_up_max,ncoef_dn_max;

    double *cpscr_cre_up  = cp->cpscr.cpscr_wave.cre_up;
    double *cpscr_cim_up  = cp->cpscr.cpscr_wave.cim_up;
    double *cpscr_cre_dn  = cp->cpscr.cpscr_wave.cre_dn;
    double *cpscr_cim_dn  = cp->cpscr.cpscr_wave.cim_dn;
    int *cpcoeffs_ioff_up    = cp->cpcoeffs_info.ioff_upt;
    int *cpcoeffs_ioff_dn    = cp->cpcoeffs_info.ioff_dnt;
    int massiv_flag          = cp->cptherm_info.massiv_flag;
    int myid_state           = class->communicate.myid_state;
    int np_states            = class->communicate.np_states;
    MPI_Comm comm_states     = class->communicate.comm_states;
    int icmoff_up        = cp->cpcoeffs_info.icoef_start_up-1;
    int icmoff_dn        = cp->cpcoeffs_info.icoef_start_dn-1;
    int icoef_form_up        = cp->cpcoeffs_pos[1].icoef_form_up;
    int icoef_orth_up        = cp->cpcoeffs_pos[1].icoef_orth_up;
    int ivcoef_form_up        = cp->cpcoeffs_pos[1].ivcoef_form_up;
    int ivcoef_orth_up        = cp->cpcoeffs_pos[1].ivcoef_orth_up;
    int ifcoef_form_up        = cp->cpcoeffs_pos[1].ifcoef_form_up;
    int ifcoef_orth_up        = cp->cpcoeffs_pos[1].ifcoef_orth_up;
    int icoef_form_dn        = cp->cpcoeffs_pos[1].icoef_form_dn;
    int icoef_orth_dn        = cp->cpcoeffs_pos[1].icoef_orth_dn;
    int ivcoef_form_dn        = cp->cpcoeffs_pos[1].ivcoef_form_dn;
    int ivcoef_orth_dn        = cp->cpcoeffs_pos[1].ivcoef_orth_dn;
    int ifcoef_form_dn        = cp->cpcoeffs_pos[1].ifcoef_form_dn;
    int ifcoef_orth_dn        = cp->cpcoeffs_pos[1].ifcoef_orth_dn;
    int cp_norb               = cp->cpopts.cp_norb;
    int cp_lsda               =  cp->cpopts.cp_lsda;
    int cp_isok_opt           = cp->cpopts.cp_isok_opt;
    int pi_beads_proc         = class->clatoms_info.pi_beads_proc;
    double *cpcoeffs_cre_up;
    double *cpcoeffs_cim_up;
    double *cpcoeffs_cre_dn;
    double *cpcoeffs_cim_dn;
    double *cpcoeffs_vcre_up;
    double *cpcoeffs_vcim_up;
    double *cpcoeffs_vcre_dn;
    double *cpcoeffs_vcim_dn;
    double *cpcoeffs_fcre_up;
    double *cpcoeffs_fcim_up;
    double *cpcoeffs_fcre_dn;
    double *cpcoeffs_fcim_dn;
    double *cpcoeffs_cmass    = cp->cpcoeffs_info.cmass;

    int nscale_up; 
    int nscale_dn;
    int nstate_up        = cp->cpcoeffs_info.nstate_up;
    int nstate_dn        = cp->cpcoeffs_info.nstate_dn;
    int ncoef_tot        = cp->cpcoeffs_info.ncoef;

    double sum_check;
     
    double temp_ext = cp->cpopts.te_ext; 
    double ac_up = 0.0;
    double bc_up = 0.0;
    double ac_dn = 0.0;
    double bc_dn = 0.0;
    double kinet_cp = 0.0; 
    double k0_up    = cp->cpcoeffs_info.k0_up;
    double k0_dn    = cp->cpcoeffs_info.k0_dn;
    nscale_up = ncoef_tot*nstate_up; 
    nscale_dn = ncoef_tot*nstate_dn;  
 
/*==========================================================================*/
/* 0) Checks                                                                */


  if(cp_norb>0){
    if((icoef_orth_up+ivcoef_orth_up+ifcoef_orth_up)!=0){
      printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
      printf("Up CP vectors are in orthonormal form \n");
      printf("on state processor %d in int_NVE_cp \n",myid_state);
      printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
      fflush(stdout);
      exit(1);
    }/*endif*/
    if(cp_lsda==1){
     if((icoef_orth_dn+ivcoef_orth_dn+ifcoef_orth_dn)!=0){
       printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
       printf("Dn CP vectors are in orthonormal form \n");
       printf("on state processor %d in int_NVE_cp \n",myid_state);
       printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
       fflush(stdout);
       exit(1);
     }/*endif*/
    }/*endif*/
  }/*endif*/

  if(np_states>1){
    if((icoef_form_up+ivcoef_form_up+ifcoef_form_up)!=3){
     printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
     printf("Up CP vectors are not in transposed form \n");
     printf("on state processor %d in int_NVE_cp \n",myid_state);
     printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
     fflush(stdout);
     exit(1);
    }/*endif*/
    if(cp_lsda==1){
     if((icoef_form_dn+ivcoef_form_dn+ifcoef_form_dn)!=3){
      printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
      printf("Up CP vectors are not in transposed form \n");
      printf("on state processor %d in int_NVE_cp \n",myid_state);
      printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
      fflush(stdout);
      exit(1);
     }/*endif*/
    }/*endif*/
  }/*endif*/

/*==========================================================================*/
/* 0) Useful constants                                                      */

    natm_tot = class->clatoms_info.natm_tot;
    general_data->timeinfo.int_res_tra    = 0;
    general_data->timeinfo.int_res_ter    = 0;
    dt  = general_data->timeinfo.dt;
    dt2 = general_data->timeinfo.dt/2.0;

    zero_constrt_iters(&(general_data->stat_avg));

    if(np_states==1){
     ncoef_up     = cp->cpcoeffs_info.ncoef;
     ncoef_dn     = cp->cpcoeffs_info.ncoef;
     ncoef_up_max = cp->cpcoeffs_info.ncoef;
     ncoef_dn_max  = cp->cpcoeffs_info.ncoef;
    }else{
     ncoef_up     = cp->cp_comm_state_pkg_up.nstate_ncoef_proc;
     ncoef_dn     = cp->cp_comm_state_pkg_dn.nstate_ncoef_proc;
     ncoef_up_max = cp->cp_comm_state_pkg_up.nstate_ncoef_proc_max;
     ncoef_dn_max = cp->cp_comm_state_pkg_dn.nstate_ncoef_proc_max;
    }/*endif*/

    ncoef_up_tot = ncoef_up_max*nstate_up;
    ncoef_dn_tot = ncoef_dn_max*nstate_dn; 

    if((cp->cptherm_info.num_c_nhc > 0) && (cp->cptherm_info.len_c_nhc > 0)) {
      cp->cptherm_info.dt_nhc  = dt;
      cp->cptherm_info.dti_nhc = dt/( (double)(cp->cptherm_info.nres_c_nhc) );
      set_yosh(cp->cptherm_info.nyosh_c_nhc,cp->cptherm_info.dti_nhc ,
               cp->cptherm_info.wdti,cp->cptherm_info.wdti2,
               cp->cptherm_info.wdti4,cp->cptherm_info.wdti8,
               cp->cptherm_info.wdti16);
    }/*endif:therms on*/


/*==========================================================================*/
/* 0.V) If gaussian dynamics apply the velocity dependent part              */
/*      of the Liouville operator                                           */

  if(cp->cpopts.cp_gauss == 1){
   apply_lgauss(&(cp->cpcoeffs_info),&(cp->cpcoeffs_pos[1]),
                &(cp->cpscr.cpscr_wave),&(cp->cpscr.cpscr_ovmat),
                &(cp->cpopts),dt);
  }

/*==========================================================================*/
/* I) First coefficient thermostat application                               */

  if((cp->cptherm_info.num_c_nhc > 0) && (cp->cptherm_info.len_c_nhc > 0)) {
      for(ip=1;ip<=pi_beads_proc;ip++) {
        if(massiv_flag==0){
          apply_c_nhc(&(cp->cptherm_info),&(cp->cptherm_pos[ip]),
                 &(cp->cpcoeffs_info),&(cp->cpcoeffs_pos[ip]),
                 &(cp->cpscr),cp->cpopts.cp_lsda,&(class->communicate));
        }else{
          apply_c_nhc_massiv(&(cp->cptherm_info),&(cp->cptherm_pos[ip]),
                          &(cp->cpcoeffs_info),&(cp->cpcoeffs_pos[ip]),
                          &(cp->cpscr),cp->cpopts.cp_lsda,
                           class->communicate.myid_state,
                           class->communicate.np_states);
        }/*endif*/
      }/*endfor*/
    }/*endif*/

/*==========================================================================*/
/* I.V)Save positions and coefficients                                        */

   for(ip=1;ip<=pi_beads_proc;ip++) {
/*-------------------------------------------------------------------------*/
     cpcoeffs_cre_up  = cp->cpcoeffs_pos[ip].cre_up;
     cpcoeffs_cim_up  = cp->cpcoeffs_pos[ip].cim_up;
     cpcoeffs_vcre_up  = cp->cpcoeffs_pos[ip].vcre_up;
     cpcoeffs_vcim_up  = cp->cpcoeffs_pos[ip].vcim_up;
     cpcoeffs_fcre_up  = cp->cpcoeffs_pos[ip].fcre_up;
     cpcoeffs_fcim_up  = cp->cpcoeffs_pos[ip].fcim_up;
     for(i=1;i<=ncoef_up_tot;i++){ 
       cpscr_cre_up[i] = cpcoeffs_cre_up[i];
       cpscr_cim_up[i] = cpcoeffs_cim_up[i];
     }/*endfor*/

     if( (cp->cpopts.cp_lsda==1) && (nstate_dn != 0) ){
       cpcoeffs_cre_dn  = cp->cpcoeffs_pos[ip].cre_dn;
       cpcoeffs_cim_dn  = cp->cpcoeffs_pos[ip].cim_dn;
       cpcoeffs_vcre_dn  = cp->cpcoeffs_pos[ip].vcre_dn;
       cpcoeffs_vcim_dn  = cp->cpcoeffs_pos[ip].vcim_dn;
       cpcoeffs_fcre_dn  = cp->cpcoeffs_pos[ip].fcre_dn;
       cpcoeffs_fcim_dn  = cp->cpcoeffs_pos[ip].fcim_dn;
       for(i=1;i<=ncoef_dn_tot;i++){ 
         cpscr_cre_dn[i] = cpcoeffs_cre_dn[i];
         cpscr_cim_dn[i] = cpcoeffs_cim_dn[i];
       }/*endfor*/
     }/*endif*/

/*==========================================================================*/
/* Calculate some factors for CP isokinetic method */

      if(cp_isok_opt == 1){
        get_isok_ab_facts(cpcoeffs_vcre_up,cpcoeffs_vcim_up,cpcoeffs_fcre_up,
                          cpcoeffs_fcim_up,cpcoeffs_cmass,icmoff_up,
                          cpcoeffs_ioff_up,nstate_up,ncoef_up,k0_up,
                          np_states,comm_states,&ac_up,&bc_up);

        if(cp_lsda == 1 && nstate_dn != 0)
           get_isok_ab_facts(cpcoeffs_vcre_dn,cpcoeffs_vcim_dn,cpcoeffs_fcre_dn,
                             cpcoeffs_fcim_dn,cpcoeffs_cmass,icmoff_dn,
                             cpcoeffs_ioff_dn,nstate_dn,ncoef_dn,k0_dn,
                             np_states,comm_states,&ac_dn,&bc_dn);
     }/* endif isok opt */

/*==========================================================================*/
/* II) Evolve Coef velocities                                               */

/*-------------------------------------------------------------------------*/
/* Up states, including isokinetic option */

     if(cp_isok_opt == 1){

       for(is=1;is<=nstate_up;is++) {
         for(i=1;i<=ncoef_up;i++) {
           icoef = i+cpcoeffs_ioff_up[is];
           cpcoeffs_vcre_up[icoef] = (cpcoeffs_vcre_up[icoef]+HF(dt2,ac_up,bc_up) 
                                   *cpcoeffs_fcre_up[icoef]/cpcoeffs_cmass[(i+icmoff_up)])
                                                         /HFDOT(dt2,ac_up,bc_up); 
           cpcoeffs_vcim_up[icoef] = (cpcoeffs_vcim_up[icoef]+HF(dt2,ac_up,bc_up) 
                     *cpcoeffs_fcim_up[icoef]/cpcoeffs_cmass[(i+icmoff_up)])
                                                         /HFDOT(dt2,ac_up,bc_up); 
       }/*endfor*/
      }/*endfor*/

     } else {/* Not isokinetic */

       for(is=1;is<=nstate_up;is++) {
         for(i=1;i<=ncoef_up;i++) {
           icoef = i+cpcoeffs_ioff_up[is];
           cpcoeffs_vcre_up[icoef] +=
             dt2*cpcoeffs_fcre_up[icoef]/cpcoeffs_cmass[(i+icmoff_up)];
           cpcoeffs_vcim_up[icoef] +=
             dt2*cpcoeffs_fcim_up[icoef]/cpcoeffs_cmass[(i+icmoff_up)];

       }/*endfor*/
      }/*endfor*/

     }/* endif isokinetic option */

/*-------------------------------------------------------------------------*/
/* Down states, including isokinetic option */


     if( (cp->cpopts.cp_lsda == 1) && (nstate_dn != 0) ){

       if(cp_isok_opt == 1){

         for(is=1;is<=nstate_dn;is++) {
            for(i=1;i<=ncoef_dn;i++) {
             icoef = i+cpcoeffs_ioff_dn[is];
             cpcoeffs_vcre_dn[icoef] = (cpcoeffs_vcre_dn[icoef]+HF(dt2,ac_dn,bc_dn) 
                                     *cpcoeffs_fcre_dn[icoef]/cpcoeffs_cmass[(i+icmoff_dn)])
                                                                    /HFDOT(dt2,ac_dn,bc_dn); 
             cpcoeffs_vcim_dn[icoef] = (cpcoeffs_vcim_dn[icoef]+HF(dt2,ac_dn,bc_dn) 
                                    *cpcoeffs_fcim_dn[icoef]/cpcoeffs_cmass[(i+icmoff_dn)])
                                                                   /HFDOT(dt2,ac_dn,bc_dn); 
          }/*endfor*/
        }/*endfor*/

       } else {/* Not isokinetic */

        for(is=1;is<=nstate_dn;is++) {
          for(i=1;i<=ncoef_dn;i++) {
            icoef = i+cpcoeffs_ioff_dn[is];
            cpcoeffs_vcre_dn[icoef] +=
                dt2*cpcoeffs_fcre_dn[icoef]/cpcoeffs_cmass[(i+icmoff_dn)];
            cpcoeffs_vcim_dn[icoef] +=
                dt2*cpcoeffs_fcim_dn[icoef]/cpcoeffs_cmass[(i+icmoff_dn)];
          }/*endfor*/
        }/*endfor*/

       }/* endif isok_opt */

     }/* endif lsda */


/*==========================================================================*/
/* IV) Evolve coefficients                                                  */ 

     for(i=1;i<=ncoef_up_tot;i++){ 
       cpcoeffs_cre_up[i] += dt*cpcoeffs_vcre_up[i];
       cpcoeffs_cim_up[i] += dt*cpcoeffs_vcim_up[i];
     }/*endfor*/

     if( (cp->cpopts.cp_lsda == 1) && (nstate_dn != 0) ){
       for(i=1;i<=ncoef_dn_tot;i++){ 
         cpcoeffs_cre_dn[i] += dt*cpcoeffs_vcre_dn[i];
         cpcoeffs_cim_dn[i] += dt*cpcoeffs_vcim_dn[i];
       }/*endfor*/
     }/* endif */

/*==========================================================================*/
/* VI) CP-Shake                                                              */

     if( (cp->cpopts.cp_gauss == 0) && (cp->cpopts.cp_norb < 3) ){
       shake_control_cp(cp,&(general_data->stat_avg.iter_shake_cp),
                      general_data->timeinfo.dt,ip);
     }/*endif*/


/*-------------------------------------------------------------------------*/
  }/*endfor : ip*/
/*-------------------------------------------------------------------------*/
     }/*end routine*/
/*==========================================================================*/



/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

double HF(double TIME, double ac, double bc)

/*==========================================================================*/
{/* Begin */
/*==========================================================================*/
/* Local Variables */

 double temp;

/*--------------------------------------------------------------------------*/
/* Calculate integration factor */

 temp =  (ac/bc)*(cosh(sqrt(bc)*TIME)-1.0);
 temp += (1.0/sqrt(bc))*sinh(sqrt(bc)*TIME);

/*-------------------------------------------------------------------------*/
 return temp;

}/*end*/ 

/*==========================================================================*/
      



/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

double HFDOT(double TIME, double ac, double bc) 

/*==========================================================================*/
{/* Begin */
/*==========================================================================*/
/* Local Variables */

  double temp;

/*--------------------------------------------------------------------------*/
/* Calculate integration factor */

  temp =  ac/sqrt(bc)*sinh(sqrt(bc)*TIME);
  temp += cosh(sqrt(bc)*TIME);

/*-------------------------------------------------------------------------*/
  return temp; 

}/*end*/ 

/*==========================================================================*/


/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void get_isok_ab_facts(double *vcre,double *vcim,double *fcre,double *fcim,
                       double *cmass,int icmoff,
                       int *ioff,int nstate,int ncoef,double k0,
                       int np_states,MPI_Comm comm_states,
                       double *ac, double *bc)

/*==========================================================================*/
{/* Begin */
/*==========================================================================*/
/* Local Variables */
#include "../typ_defs/typ_mask.h"

  double ac_ret=0.0,bc_ret=0.0;
  double ac_tmp,bc_tmp;
  int is,i,icoef;

/*-------------------------------------------------------------------------*/
/* Calculate the a and b factors */


     for(is=1;is<=nstate;is++) { 
       for(i=1;i<=ncoef;i++) {
          icoef = i+ioff[is];
          ac_ret += vcre[icoef]*fcre[icoef] + vcim[icoef]*fcim[icoef];
          bc_ret += (fcre[icoef]*fcre[icoef] +fcim[icoef]*fcim[icoef])
                   /cmass[(i+icmoff)];
       }/*endfor*/
     }/*endfor*/ 

     if(np_states > 1) {
       ac_tmp = ac_ret;
       bc_tmp = bc_ret;
       Allreduce(&ac_tmp,&ac_ret,1,MPI_DOUBLE,MPI_SUM,0,comm_states);
       Allreduce(&bc_tmp,&bc_ret,1,MPI_DOUBLE,MPI_SUM,0,comm_states);
     }/*endif */

    ac_ret /= (2.0*k0);
    bc_ret /= (2.0*k0);   

    *ac = ac_ret;
    *bc = bc_ret;

/*-------------------------------------------------------------------------*/
}/*end*/ 
/*==========================================================================*/









