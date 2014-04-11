/*===================================================================*/
/*ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*===================================================================*/
/*                                                                   */
/*                         PI_MD:                                    */
/*             The future of simulation technology                   */
/*             ------------------------------------                  */
/*                     Module: proj_vel.c                            */
/*                                                                   */
/* These subprograms projects the velocities onto surface of         */
/* constraint                                                        */
/*                                                                   */
/*===================================================================*/
/*ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*===================================================================*/

#include "standard_include.h"
#include "../typ_defs/typedefs_gen.h"
#include "../typ_defs/typedefs_class.h"
#include "../typ_defs/typedefs_cp.h"
#include "../typ_defs/typedefs_bnd.h"
#include "../proto_defs/proto_vel_sampl_cp_local.h"
#include "../proto_defs/proto_vel_sampl_cp_entry.h"
#include "../proto_defs/proto_energy_cpcon_entry.h"

/*===================================================================*/
/*ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*===================================================================*/
/* Sample the coefficient velocities                                 */
/*===================================================================*/

void proj_velc(CP *cp,int ip)

/*===================================================================*/
{/*begin routine*/
/*===================================================================*/
/*             Local variable declarations                            */
/*===================================================================*/

  int i,iii;
  int ncoef_tot_up,ncoef_tot_dn;
  double dt           = 0.0001;
  int iter_shake_cp   = 0;

/* Local pointers */

  int ncoef_up        = cp->cpcoeffs_info.ncoef;
  int ncoef_dn        = cp->cpcoeffs_info.ncoef;
  int nstate_up       = cp->cpcoeffs_info.nstate_up;
  int nstate_dn       = cp->cpcoeffs_info.nstate_dn;
  int cp_lsda         = cp->cpopts.cp_lsda;
  double *cre_up      = cp->cpcoeffs_pos[ip].cre_up; 
  double *cim_up      = cp->cpcoeffs_pos[ip].cim_up; 
  double *cre_dn      = cp->cpcoeffs_pos[ip].cre_dn; 
  double *cim_dn      = cp->cpcoeffs_pos[ip].cim_dn; 
  double *vcre_up     = cp->cpcoeffs_pos[ip].vcre_up; 
  double *vcim_up     = cp->cpcoeffs_pos[ip].vcim_up; 
  double *vcre_dn     = cp->cpcoeffs_pos[ip].vcre_dn; 
  double *vcim_dn     = cp->cpcoeffs_pos[ip].vcim_dn; 
  int ivcoef_form_up  = cp->cpcoeffs_pos[ip].ivcoef_form_up;
  int icoef_form_up   = cp->cpcoeffs_pos[ip].icoef_form_up;
  int icoef_orth_up   = cp->cpcoeffs_pos[ip].icoef_orth_up;
  int ivcoef_form_dn  = cp->cpcoeffs_pos[ip].ivcoef_form_dn;
  int icoef_form_dn   = cp->cpcoeffs_pos[ip].icoef_form_dn;
  int icoef_orth_dn   = cp->cpcoeffs_pos[ip].icoef_orth_dn;

  double *scr_vcre_up = cp->cpcoeffs_pos[ip].fcre_up; 
  double *scr_vcim_up = cp->cpcoeffs_pos[ip].fcim_up; 
  double *scr_vcre_dn = cp->cpcoeffs_pos[ip].fcre_dn; 
  double *scr_vcim_dn = cp->cpcoeffs_pos[ip].fcim_dn; 
  double *scr_cre_up  = cp->cpscr.cpscr_wave.cre_up; 
  double *scr_cim_up  = cp->cpscr.cpscr_wave.cim_up; 
  double *scr_cre_dn  = cp->cpscr.cpscr_wave.cre_dn; 
  double *scr_cim_dn  = cp->cpscr.cpscr_wave.cim_dn; 
  int cp_norb         = cp->cpopts.cp_norb;
  int myid            = cp->cp_comm_state_pkg_up.myid;
  int np_states       = cp->cp_comm_state_pkg_up.num_proc;

/*===================================================================*/
/* 0) Checks                                                         */

   if(np_states>1){
      ncoef_up = cp->cp_comm_state_pkg_up.nstate_ncoef_proc_max;
      if(ivcoef_form_up!=1&&icoef_form_up!=1){
        printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
        printf("Up Coef and/or coef velocities are not in transposed form \n");
        printf("on state processor %d in proj_vel_cp \n",myid);
        printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
        fflush(stdout);
        exit(1);
      }/*endif*/
      if(cp_lsda==1){
       ncoef_dn = cp->cp_comm_state_pkg_dn.nstate_ncoef_proc_max;
       if(ivcoef_form_dn!=1&&icoef_form_dn!=1){
        printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
        printf("Dn Coef and/or coef velocities are not in transposed form \n");
        printf("on state processor %d in proj_vel_cp \n",myid);
        printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
        fflush(stdout);
        exit(1);
       }/*endif*/
      }/*endif*/
   }/*endif*/

   if(icoef_orth_up!=0 && cp_norb>0){
      printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
      printf("Up Coefs must be in nonorthonormal form under norb \n");
      printf("in proj_vel_cp \n");
      printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
      fflush(stdout);
      exit(1);
   }/*endif*/
   if(cp_lsda==1){
    if(icoef_orth_dn!=0 && cp_norb>0){
      printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
      printf("Dn Coefs must be in nonorthonormal form under norb \n");
      printf("in proj_vel_cp \n");
      printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
      fflush(stdout);
      exit(1);
    }/*endif*/
   }/*endif*/

/*===================================================================*/
/* I) Copy system                                                    */

  ncoef_tot_up = ncoef_up*nstate_up;
  for(i=1;i<=ncoef_tot_up;i++){
    scr_cre_up[i]  =  cre_up[i];   
    scr_cim_up[i]  =  cim_up[i];   
    scr_vcre_up[i] = vcre_up[i];
    scr_vcim_up[i] = vcim_up[i];
  }/*endfor*/

  if(cp_lsda == 1){
    ncoef_tot_dn = ncoef_dn*nstate_dn;
    for(i=1;i<=ncoef_tot_dn;i++){
      scr_cre_dn[i]  = cre_dn[i];   
      scr_cim_dn[i]  = cim_dn[i];   
      scr_vcre_dn[i] = vcre_dn[i];
      scr_vcim_dn[i] = vcim_dn[i];
    }/*endfor*/
  }/*endif*/

/*===================================================================*/
/* II) Evolve system and shake once to get on surface of constraint  */

  for(i=1;i<=ncoef_tot_up;i++){
    vcre_up[i] = 0.0;
    vcim_up[i] = 0.0;
    cre_up[i]  = scr_cre_up[i] + scr_vcre_up[i]*dt;  
    cim_up[i]  = scr_cim_up[i] + scr_vcim_up[i]*dt;  
  }/*endfor*/

  if(cp_lsda==1) {
    for(i=1;i<=ncoef_tot_dn;i++){
      vcre_dn[i] = 0.0;
      vcim_dn[i] = 0.0;
      cre_dn[i] = scr_cre_dn[i] + scr_vcre_dn[i]*dt;  
      cim_dn[i] = scr_cim_dn[i] + scr_vcim_dn[i]*dt;  
    }/*endfor*/
  }/*endif*/

  shake_control_cp(cp,&iter_shake_cp,dt,ip);


  for(i=1;i<=ncoef_tot_up;i++){
    scr_cre_up[i] = cre_up[i];   
    scr_cim_up[i] = cim_up[i];   
  }/*endfor*/

  if(cp_lsda == 1){
    for(i=1;i<=ncoef_tot_dn;i++){
      scr_cre_dn[i] = cre_dn[i];   
      scr_cim_dn[i] = cim_dn[i];   
    }/*endfor*/
  }/*endif*/

/*===================================================================*/
/* III) Evolve system and shake again to get velocities              */

  for(i=1;i<=ncoef_tot_up;i++){
    vcre_up[i] = 0.0;
    vcim_up[i] = 0.0;
    cre_up[i] = scr_cre_up[i] + scr_vcre_up[i]*dt;  
    cim_up[i] = scr_cim_up[i] + scr_vcim_up[i]*dt;  
  }/*endfor*/

  if(cp_lsda==1) {
   for(i=1;i<=ncoef_tot_dn;i++){
     vcre_dn[i] = 0.0;
     vcim_dn[i] = 0.0;
     cre_dn[i] = scr_cre_dn[i] + scr_vcre_dn[i]*dt;  
     cim_dn[i] = scr_cim_dn[i] + scr_vcim_dn[i]*dt;  
   }/*endfor*/
  }/*endif*/

  shake_control_cp(cp,&iter_shake_cp,dt,ip);

  for(i=1;i<=ncoef_tot_up;i++){
    vcre_up[i] = (cre_up[i]-scr_cre_up[i])/dt;
    vcim_up[i] = (cim_up[i]-scr_cim_up[i])/dt;
  }/*endfor*/

  if(cp_lsda==1) {
    for(i=1;i<=ncoef_tot_dn;i++){
      vcre_dn[i] = (cre_dn[i]-scr_cre_dn[i])/dt;
      vcim_dn[i] = (cim_dn[i]-scr_cim_dn[i])/dt;
    }/*endfor*/
  }/*endif*/

/*-------------------------------------------------------------------*/
   }/*end routine*/
/*===================================================================*/









