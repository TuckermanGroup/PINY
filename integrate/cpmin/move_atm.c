/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
/*                                                                          */
/*                         PI_MD:                                           */
/*             The future of simulation technology                          */
/*             ------------------------------------                         */
/*                     Module: move_atm                                     */
/*                                                                          */
/* This subprogram minimizes atomic positions using either                  */
/* steepest descent, conjugate gradient or DIIS                             */
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
#include "../proto_defs/proto_friend_lib_entry.h"
#include "../proto_defs/proto_communicate_wrappers.h"
#include "../proto_defs/proto_math.h"

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void move_atm_cg(CLASS *class,BONDED *bonded,GENERAL_DATA *general_data,
                 CP *cp,int ifirst)

/*========================================================================*/
{/*begin routine*/
/*========================================================================*/
/*             Local variable declarations                                */

#include "../typ_defs/typ_mask.h"

    int i,ip=1,ipart,icoef,is;
    double dt,dts;
    int natm_tot,iii;
    int nstate_up,nstate_dn;
    int maxcoef;
    double gamma;
    static double fovlap;
    double mass_sc_fact  = class->clatoms_info.mass_sc_fact;
    double fovlap_old;
    double *clatoms_x    = class->clatoms_pos[1].x;
    double *clatoms_y    = class->clatoms_pos[1].y;
    double *clatoms_z    = class->clatoms_pos[1].z;
    double *clatoms_fx   = class->clatoms_pos[1].fx;
    double *clatoms_fy   = class->clatoms_pos[1].fy;
    double *clatoms_fz   = class->clatoms_pos[1].fz;
    double *grad_x       = class->ewd_scr.fx;
    double *grad_y       = class->ewd_scr.fy;
    double *grad_z       = class->ewd_scr.fz;
    double *hess_xx      = class->clatoms_pos[1].hess_xx;
    double *hess_xy      = class->clatoms_pos[1].hess_xy;
    double *hess_xz      = class->clatoms_pos[1].hess_xz;
    double *hess_yy      = class->clatoms_pos[1].hess_yy;
    double *hess_yz      = class->clatoms_pos[1].hess_yz;
    double *hess_zz      = class->clatoms_pos[1].hess_zz;
    double *chx          = class->clatoms_pos[1].vx;
    double *chy          = class->clatoms_pos[1].vy;
    double *chz          = class->clatoms_pos[1].vz;
    double *scr_x        = class->clatoms_info.xold;
    double *scr_y        = class->clatoms_info.yold;
    double *scr_z        = class->clatoms_info.zold;
    double *clatoms_mass = class->clatoms_info.mass;
    double *cre_up       = cp->cpcoeffs_pos[1].cre_up;
    double *cim_up       = cp->cpcoeffs_pos[1].cim_up;
    double *cre_dn       = cp->cpcoeffs_pos[1].cre_dn;
    double *cim_dn       = cp->cpcoeffs_pos[1].cim_dn;
    double *fcre_up      = cp->cpcoeffs_pos[1].fcre_up;
    double *fcim_up      = cp->cpcoeffs_pos[1].fcim_up;
    double *fcre_dn      = cp->cpcoeffs_pos[1].fcre_dn;
    double *fcim_dn      = cp->cpcoeffs_pos[1].fcim_dn;
    double *zeta_up,*zeta_dn;
    double *cmass        = cp->cpcoeffs_info.cmass;
    int *ioff_up         = cp->cpcoeffs_info.ioff_upt;
    int *ioff_dn         = cp->cpcoeffs_info.ioff_dnt;
    double fc_mag_up, fc_mag_dn;
    int np_states        = class->communicate.np_states;
    int myid             = class->communicate.myid;
    int myid_state = class->communicate.myid_state;
    MPI_Comm comm_states = class->communicate.comm_states;

    int hess_calc        = class->clatoms_info.hess_calc;
    int cp_norb          = cp->cpopts.cp_norb;
    int cp_lsda          = cp->cpopts.cp_lsda;
    int icmoff_up        = cp->cpcoeffs_info.icoef_start_up-1;
    int icmoff_dn        = cp->cpcoeffs_info.icoef_start_dn-1;
    int icoef_form_up    = cp->cpcoeffs_pos[1].icoef_form_up;
    int icoef_orth_up    = cp->cpcoeffs_pos[1].icoef_orth_up;
    int ifcoef_form_up   = cp->cpcoeffs_pos[1].ifcoef_form_up;
    int ifcoef_orth_up   = cp->cpcoeffs_pos[1].ifcoef_orth_up;
    int icoef_form_dn    = cp->cpcoeffs_pos[1].icoef_form_dn;
    int icoef_orth_dn    = cp->cpcoeffs_pos[1].icoef_orth_dn;
    int ifcoef_form_dn   = cp->cpcoeffs_pos[1].ifcoef_form_dn;
    int ifcoef_orth_dn   = cp->cpcoeffs_pos[1].ifcoef_orth_dn;
    int min_atm_com_fix_opt = general_data->minopts.min_atm_com_fix_opt;
    int ncoef_up,ncoef_dn;

/*==========================================================================*/
/* 0) Checks                                                                */


  if(cp_norb>0){
    if((icoef_orth_up+ifcoef_orth_up)!=0){
      printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
      printf("Up CP vectors are in orthonormal form \n");
      printf("on state processor %d in move_atm_cg \n",myid_state);
      printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
      fflush(stdout);
      exit(1);
    }/*endif*/
    if(cp_lsda==1){
     if((icoef_orth_dn+ifcoef_orth_dn)!=0){
       printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
       printf("Dn CP vectors are in orthonormal form \n");
       printf("on state processor %d in move_atm_cg \n",myid_state);
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
     printf("on state processor %d in move_atm_cg \n",myid_state);
     printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
     fflush(stdout);
     exit(1);
    }/*endif*/
    if(cp_lsda==1){
     if((icoef_form_dn+ifcoef_form_dn)!=2){
      printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
      printf("Up CP vectors are not in transposed form \n");
      printf("on state processor %d in move_atm_cg \n",myid_state);
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
    dts = dt*mass_sc_fact;
    zero_constrt_iters(&(general_data->stat_avg));

    nstate_up    = cp->cpcoeffs_info.nstate_up;
    nstate_dn    = cp->cpcoeffs_info.nstate_dn;
    if(np_states==1){
     ncoef_up      = cp->cpcoeffs_info.ncoef;
     ncoef_dn      = cp->cpcoeffs_info.ncoef;
    }else{
     ncoef_up     = cp->cp_comm_state_pkg_up.nstate_ncoef_proc;
     ncoef_dn     = cp->cp_comm_state_pkg_dn.nstate_ncoef_proc;
    }/*endif*/

/*==========================================================================*/
/* 0.1) Zero conjugate gradients                                            */

   if(ifirst == 1){
     for(i=1;i<=natm_tot; i++){
      chx[i] = 0.0;
      chy[i] = 0.0;
      chz[i] = 0.0;
     }
     gamma = 0.0;
     fovlap = 1.0;
    }/* endif */
    
/*==========================================================================*/
/* I) Get forces                                                            */

    (class->energy_ctrl.iget_full_inter)= 1;
    (class->energy_ctrl.iget_res_inter) = 0;
    (class->energy_ctrl.iget_full_intra)= 1;
    (class->energy_ctrl.iget_res_intra) = 0;

    cp_energy_control(class,bonded,general_data,cp);
/*==========================================================================*/
/* II) Calculate the gamma's                                                */
 

 if(myid_state==0){
   if(min_atm_com_fix_opt==1){
     proj_com_out(natm_tot,clatoms_fx,clatoms_fy,clatoms_fz);
   }/*endif*/
#ifdef WILL_PRINT
   printf("================\n");
   printf("Forces on Proc 0\n");
   for(i=1;i<=natm_tot; i++){
     printf("atm %d : %g %g %g\n",i,clatoms_fx[i],clatoms_fy[i],clatoms_fz[i]);
   }/*endfor*/
   printf("================\n");
#endif

   fovlap_old = fovlap;
   fovlap = 0.0;
   for(i=1;i<=natm_tot;i++){
    fovlap += clatoms_fx[i]*clatoms_fx[i]
            + clatoms_fy[i]*clatoms_fy[i]
            + clatoms_fz[i]*clatoms_fz[i];
   }
   if(ifirst != 1) {gamma = fovlap/fovlap_old;}

/*==========================================================================*/
/* II.V) Evolve gradients                                                   */

     for(i=1;i<=natm_tot; i++){
      chx[i] = clatoms_fx[i] + gamma*chx[i];
      chy[i] = clatoms_fy[i] + gamma*chy[i];
      chz[i] = clatoms_fz[i] + gamma*chz[i];
     }/* endfor */

/*==========================================================================*/
/* III) Calculate the step length                                           */

     if(hess_calc >= 2)
       act_hess_inv_on_grad(hess_xx,hess_xy,hess_xz,hess_yy,hess_yz,hess_zz,
                                 grad_x,grad_y,grad_z,chx,chy,chz,
                                 natm_tot);
     if(hess_calc == 1)
       act_hess_inv_on_grad_diag(hess_xx,hess_xy,hess_xz,hess_yy,hess_yz,hess_zz,
                            grad_x,grad_y,grad_z,chx,chy,chz,
                            natm_tot);
     if(hess_calc == 0){
      for(i=1;i<=natm_tot;i++){
        grad_x[i] = chx[i];
        grad_y[i] = chy[i];
        grad_z[i] = chz[i];
      }/* endfor */
    }/* end switch */


/*==========================================================================*/
/* IV) Evolve positions                                                     */ 

  for(i=1;i<=natm_tot;i++){
    clatoms_x[i] += dt*grad_x[i];
    clatoms_y[i] += dt*grad_y[i];
    clatoms_z[i] += dt*grad_z[i];
   }/*endfor*/

   get_ghost_pos(&(class->clatoms_info),&(class->clatoms_pos[1]),
                 &(class->ghost_atoms));
 }/* endif myid_state */

 if(np_states > 1){
   Bcast(&(clatoms_x[1]),natm_tot,MPI_DOUBLE,0,comm_states);
   Bcast(&(clatoms_y[1]),natm_tot,MPI_DOUBLE,0,comm_states);
   Bcast(&(clatoms_z[1]),natm_tot,MPI_DOUBLE,0,comm_states);
 }/* endif */

/*==========================================================================*/
/* IV.V) Do a steepest descent step for coefficients                        */ 

   for(is=1;is<=nstate_up;is++) {
     for(i=1;i<=ncoef_up;i++) {
       icoef = i+ioff_up[is];
       cre_up[icoef] += dt*fcre_up[icoef]/cmass[(i+icmoff_up)];
       cim_up[icoef] += dt*fcim_up[icoef]/cmass[(i+icmoff_up)];
     }/*endfor*/
   }/*endfor*/


 if( (cp->cpopts.cp_lsda == 1) && (nstate_dn != 0) ){
   for(is=1;is<=nstate_dn;is++) {
     for(i=1;i<=ncoef_dn;i++) {
       icoef = i+ioff_dn[is];
       cre_dn[icoef] += dt*fcre_dn[icoef]/cmass[(i+icmoff_dn)];
       cim_dn[icoef] += dt*fcim_dn[icoef]/cmass[(i+icmoff_dn)];
     }/*endfor*/
   }/*endfor*/
  }/* endif */

/*==========================================================================*/
/* V) Orthogonalize wave functions                                          */

  orthog_control_cp(cp,ip);

/*==========================================================================*/
/* ii) Free memory and shuffle states                                       */

   cp_shuffle_states(cp,ip);

/*========================================================================*/
}/*end routine*/
/*========================================================================*/




/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
 
void move_atm_std(CLASS *class,BONDED *bonded,GENERAL_DATA *general_data,
                  CP *cp,int ifirst)

/*========================================================================*/
{/*begin routine*/
/*========================================================================*/
/*             Local variable declarations                                */

#include "../typ_defs/typ_mask.h"
    int i,ip=1,ipart,icoef,is;
    double dt,dts;
    int natm_tot,iii;
    int nstate_up,nstate_dn;
    int hess_calc        = class->clatoms_info.hess_calc;
    double mass_sc_fact  = class->clatoms_info.mass_sc_fact;
    double *clatoms_x    = class->clatoms_pos[1].x;
    double *clatoms_y    = class->clatoms_pos[1].y;
    double *clatoms_z    = class->clatoms_pos[1].z;
    double *clatoms_fx   = class->clatoms_pos[1].fx;
    double *clatoms_fy   = class->clatoms_pos[1].fy;
    double *clatoms_fz   = class->clatoms_pos[1].fz;
    double *grad_x       = class->ewd_scr.fx;
    double *grad_y       = class->ewd_scr.fy;
    double *grad_z       = class->ewd_scr.fz;
    double *hess_xx      = class->clatoms_pos[1].hess_xx;
    double *hess_xy      = class->clatoms_pos[1].hess_xy;
    double *hess_xz      = class->clatoms_pos[1].hess_xz;
    double *hess_yy      = class->clatoms_pos[1].hess_yy;
    double *hess_yz      = class->clatoms_pos[1].hess_yz;
    double *hess_zz      = class->clatoms_pos[1].hess_zz;
    double *chx          = class->clatoms_pos[1].vx;
    double *chy          = class->clatoms_pos[1].vy;
    double *chz          = class->clatoms_pos[1].vz;
    double *scr_x        = class->clatoms_info.xold;
    double *scr_y        = class->clatoms_info.yold;
    double *scr_z        = class->clatoms_info.zold;
    double *clatoms_mass = class->clatoms_info.mass;
    double *cre_up       = cp->cpcoeffs_pos[1].cre_up;
    double *cim_up       = cp->cpcoeffs_pos[1].cim_up;
    double *cre_dn       = cp->cpcoeffs_pos[1].cre_dn;
    double *cim_dn       = cp->cpcoeffs_pos[1].cim_dn;
    double *fcre_up      = cp->cpcoeffs_pos[1].fcre_up;
    double *fcim_up      = cp->cpcoeffs_pos[1].fcim_up;
    double *fcre_dn      = cp->cpcoeffs_pos[1].fcre_dn;
    double *fcim_dn      = cp->cpcoeffs_pos[1].fcim_dn;
    double *zeta_up,*zeta_dn;
    double *cmass        = cp->cpcoeffs_info.cmass;
    int *ioff_up         = cp->cpcoeffs_info.ioff_upt;
    int *ioff_dn         = cp->cpcoeffs_info.ioff_dnt;
    double fc_mag_up, fc_mag_dn;
    int ncoef_tot;
    int np_states        = class->communicate.np_states;
    int myid_state       = class->communicate.myid_state;
    MPI_Comm comm_states = class->communicate.comm_states;

    int cp_norb               = cp->cpopts.cp_norb;
    int cp_lsda               = cp->cpopts.cp_lsda;
    int icmoff_up             = cp->cpcoeffs_info.icoef_start_up-1;
    int icmoff_dn             = cp->cpcoeffs_info.icoef_start_dn-1;
    int icoef_form_up         = cp->cpcoeffs_pos[1].icoef_form_up;
    int icoef_orth_up         = cp->cpcoeffs_pos[1].icoef_orth_up;
    int ifcoef_form_up        = cp->cpcoeffs_pos[1].ifcoef_form_up;
    int ifcoef_orth_up        = cp->cpcoeffs_pos[1].ifcoef_orth_up;
    int icoef_form_dn         = cp->cpcoeffs_pos[1].icoef_form_dn;
    int icoef_orth_dn         = cp->cpcoeffs_pos[1].icoef_orth_dn;
    int ifcoef_form_dn        = cp->cpcoeffs_pos[1].ifcoef_form_dn;
    int ifcoef_orth_dn        = cp->cpcoeffs_pos[1].ifcoef_orth_dn;
    int min_atm_com_fix_opt = general_data->minopts.min_atm_com_fix_opt;
    int ncoef_up,ncoef_dn;

/*==========================================================================*/
/* 0) Checks                                                                */

  if(cp_norb>0){
    if((icoef_orth_up+ifcoef_orth_up)!=0){
      printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
      printf("Up CP vectors are in orthonormal form \n");
      printf("on state processor %d in move_atm_std \n",myid_state);
      printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
      fflush(stdout);
      exit(1);
    }/*endif*/
    if(cp_lsda==1){
     if((icoef_orth_dn+ifcoef_orth_dn)!=0){
       printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
       printf("Dn CP vectors are in orthonormal form \n");
       printf("on state processor %d in move_atm_std \n",myid_state);
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
     printf("on state processor %d in move_atm_std \n",myid_state);
     printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
     fflush(stdout);
     exit(1);
    }/*endif*/
    if(cp_lsda==1){
     if((icoef_form_dn+ifcoef_form_dn)!=2){
      printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
      printf("Up CP vectors are not in transposed form \n");
      printf("on state processor %d in move_atm_std \n",myid_state);
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
    dts = dt*mass_sc_fact;
    zero_constrt_iters(&(general_data->stat_avg));

    nstate_up    = cp->cpcoeffs_info.nstate_up;
    nstate_dn    = cp->cpcoeffs_info.nstate_dn;
    if(np_states==1){
     ncoef_up      = cp->cpcoeffs_info.ncoef;
     ncoef_dn      = cp->cpcoeffs_info.ncoef;
    }else{
     ncoef_up     = cp->cp_comm_state_pkg_up.nstate_ncoef_proc;
     ncoef_dn     = cp->cp_comm_state_pkg_dn.nstate_ncoef_proc;
    }/*endif*/

/*==========================================================================*/
/* I) Get forces                                                            */

    (class->energy_ctrl.iget_full_inter)= 1;
    (class->energy_ctrl.iget_res_inter) = 0;
    (class->energy_ctrl.iget_full_intra)= 1;
    (class->energy_ctrl.iget_res_intra) = 0;

    cp_energy_control(class,bonded,general_data,cp);
/*==========================================================================*/
/* II) Get forces                                                            */

 if(myid_state == 0){

    if(min_atm_com_fix_opt==1){
      proj_com_out(natm_tot,clatoms_fx,clatoms_fy,clatoms_fz);
    }/*endif*/

     if(hess_calc >= 2)
      act_hess_inv_on_grad(hess_xx,hess_xy,hess_xz,hess_yy,hess_yz,hess_zz,
                                grad_x,grad_y,grad_z,clatoms_fx,clatoms_fy,clatoms_fz,
                                natm_tot);
     if(hess_calc == 1)
      act_hess_inv_on_grad_diag(hess_xx,hess_xy,hess_xz,hess_yy,hess_yz,hess_zz,
                           grad_x,grad_y,grad_z,clatoms_fx,clatoms_fy,clatoms_fz,
                           natm_tot);
     if(hess_calc == 0){
      for(i=1;i<=natm_tot;i++){
        grad_x[i] = clatoms_fx[i];
        grad_y[i] = clatoms_fy[i];
        grad_z[i] = clatoms_fz[i];
      }/* endfor */
    }/* endif */

/*==========================================================================*/
/* III) Evolve positions                                                     */ 

  for(i=1;i<=natm_tot;i++){
      clatoms_x[i] += dts*grad_x[i]/clatoms_mass[i];
      clatoms_y[i] += dts*grad_y[i]/clatoms_mass[i];
      clatoms_z[i] += dts*grad_z[i]/clatoms_mass[i];
   }/*endfor*/

   get_ghost_pos(&(class->clatoms_info),&(class->clatoms_pos[1]),
                 &(class->ghost_atoms));
}/* endif myid_state */

 if(np_states > 1){
   Bcast(&(clatoms_x[1]),natm_tot,MPI_DOUBLE,0,comm_states);
   Bcast(&(clatoms_y[1]),natm_tot,MPI_DOUBLE,0,comm_states);
   Bcast(&(clatoms_z[1]),natm_tot,MPI_DOUBLE,0,comm_states);
 }/* endif */

/*==========================================================================*/
/* IV.V) Do a steepest descent step for coefficients                        */ 

   for(is=1;is<=nstate_up;is++) {
     for(i=1;i<=ncoef_up;i++) {
       icoef = i+ioff_up[is];
       cre_up[icoef] += dt*fcre_up[icoef]/cmass[(i+icmoff_up)];
       cim_up[icoef] += dt*fcim_up[icoef]/cmass[(i+icmoff_up)];
     }/*endfor*/
   }/*endfor*/

 if( (cp->cpopts.cp_lsda == 1) && (nstate_dn != 0) ){
   for(is=1;is<=nstate_dn;is++) {
     for(i=1;i<=ncoef_dn;i++) {
       icoef = i+ioff_dn[is];
       cre_dn[icoef] += dt*fcre_dn[icoef]/cmass[(i+icmoff_dn)];
       cim_dn[icoef] += dt*fcim_dn[icoef]/cmass[(i+icmoff_dn)];
     }/*endfor*/
   }/*endfor*/
  }/* endif */

/*==========================================================================*/
/* V) Orthogonalize wave functions                                          */

  orthog_control_cp(cp,ip);

/*==========================================================================*/
/* ii) Free memory and shuffle states                                       */

  cp_shuffle_states(cp,ip);

/*========================================================================*/
}/*end routine*/
/*========================================================================*/
 



/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
 
void move_atm_diis(CLASS *class,BONDED *bonded,GENERAL_DATA *general_data,
                   CP *cp, int ifirst)

/*========================================================================*/
{/*begin routine*/
/*========================================================================*/
/* Local Variable declarations                                          */
#include "../typ_defs/typ_mask.h"
    int i,ip=1,ipart,icoef,is,ist,iist,job;
    int info=0;
    double dt,dts;
    int natm_tot,iii;
    int nstate_up,nstate_dn;
    static int itcount=1;
    int mvec,mp1;
    double mass_sc_fact  = class->clatoms_info.mass_sc_fact;
    double *clatoms_x    = class->clatoms_pos[1].x;
    double *clatoms_y    = class->clatoms_pos[1].y;
    double *clatoms_z    = class->clatoms_pos[1].z;
    double *clatoms_fx   = class->clatoms_pos[1].fx;
    double *clatoms_fy   = class->clatoms_pos[1].fy;
    double *clatoms_fz   = class->clatoms_pos[1].fz;
    double *grad_x       = class->ewd_scr.fx;
    double *grad_y       = class->ewd_scr.fy;
    double *grad_z       = class->ewd_scr.fz;
    double *hess_xx      = class->clatoms_pos[1].hess_xx;
    double *hess_xy      = class->clatoms_pos[1].hess_xy;
    double *hess_xz      = class->clatoms_pos[1].hess_xz;
    double *hess_yy      = class->clatoms_pos[1].hess_yy;
    double *hess_yz      = class->clatoms_pos[1].hess_yz;
    double *hess_zz      = class->clatoms_pos[1].hess_zz;
    double *chx          = class->clatoms_pos[1].vx;
    double *chy          = class->clatoms_pos[1].vy;
    double *chz          = class->clatoms_pos[1].vz;
    double *scr_x        = class->clatoms_info.xold;
    double *scr_y        = class->clatoms_info.yold;
    double *scr_z        = class->clatoms_info.zold;
    double *clatoms_mass = class->clatoms_info.mass;
    double *cre_up       = cp->cpcoeffs_pos[1].cre_up;
    double *cim_up       = cp->cpcoeffs_pos[1].cim_up;
    double *cre_dn       = cp->cpcoeffs_pos[1].cre_dn;
    double *cim_dn       = cp->cpcoeffs_pos[1].cim_dn;
    double *fcre_up      = cp->cpcoeffs_pos[1].fcre_up;
    double *fcim_up      = cp->cpcoeffs_pos[1].fcim_up;
    double *fcre_dn      = cp->cpcoeffs_pos[1].fcre_dn;
    double *fcim_dn      = cp->cpcoeffs_pos[1].fcim_dn;
    double *zeta_up,*zeta_dn;
    double *cmass        = cp->cpcoeffs_info.cmass;
    int *ioff_up         = cp->cpcoeffs_info.ioff_upt;
    int *ioff_dn         = cp->cpcoeffs_info.ioff_dnt;
    double fc_mag_up, fc_mag_dn;
    int ncoef_tot;
    int np_states        = class->communicate.np_states;
    int myid_state       = class->communicate.myid_state;
    MPI_Comm comm_states = class->communicate.comm_states;
    int MVEC_MAX_ATM     = general_data->minopts.diis_hist_len;
    int hess_calc        = class->clatoms_info.hess_calc;

    int cp_norb          = cp->cpopts.cp_norb;
    int cp_lsda          = cp->cpopts.cp_lsda;
    int icmoff_up        = cp->cpcoeffs_info.icoef_start_up-1;
    int icmoff_dn        = cp->cpcoeffs_info.icoef_start_dn-1;
    int icoef_form_up    = cp->cpcoeffs_pos[1].icoef_form_up;
    int icoef_orth_up    = cp->cpcoeffs_pos[1].icoef_orth_up;
    int ifcoef_form_up   = cp->cpcoeffs_pos[1].ifcoef_form_up;
    int ifcoef_orth_up   = cp->cpcoeffs_pos[1].ifcoef_orth_up;
    int icoef_form_dn    = cp->cpcoeffs_pos[1].icoef_form_dn;
    int icoef_orth_dn    = cp->cpcoeffs_pos[1].icoef_orth_dn;
    int ifcoef_form_dn   = cp->cpcoeffs_pos[1].ifcoef_form_dn;
    int ifcoef_orth_dn   = cp->cpcoeffs_pos[1].ifcoef_orth_dn;
    int min_atm_com_fix_opt = general_data->minopts.min_atm_com_fix_opt;
    int ncoef_up,ncoef_dn;

/* DIIS local pointers */
    static int *ipvt;
    static double *pmat;
    static double *vect;
    static double *xst,*yst,*zst;
    static double *fxst,*fyst,*fzst;

/*==========================================================================*/
/* 0) Checks                                                                */

  if(cp_norb>0){
    if((icoef_orth_up+ifcoef_orth_up)!=0){
      printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
      printf("Up CP vectors are in orthonormal form \n");
      printf("on state processor %d in move_atm_std \n",myid_state);
      printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
      fflush(stdout);
      exit(1);
    }/*endif*/
    if(cp_lsda==1){
     if((icoef_orth_dn+ifcoef_orth_dn)!=0){
       printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
       printf("Dn CP vectors are in orthonormal form \n");
       printf("on state processor %d in move_atm_std \n",myid_state);
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
     printf("on state processor %d in move_atm_std \n",myid_state);
     printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
     fflush(stdout);
     exit(1);
    }/*endif*/
    if(cp_lsda==1){
     if((icoef_form_dn+ifcoef_form_dn)!=2){
      printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
      printf("Up CP vectors are not in transposed form \n");
      printf("on state processor %d in move_atm_std \n",myid_state);
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
    dts = dt*mass_sc_fact;
    zero_constrt_iters(&(general_data->stat_avg));

    nstate_up    = cp->cpcoeffs_info.nstate_up;
    nstate_dn    = cp->cpcoeffs_info.nstate_dn;
    if(np_states==1){
     ncoef_up      = cp->cpcoeffs_info.ncoef;
     ncoef_dn      = cp->cpcoeffs_info.ncoef;
    }else{
     ncoef_up     = cp->cp_comm_state_pkg_up.nstate_ncoef_proc;
     ncoef_dn     = cp->cp_comm_state_pkg_dn.nstate_ncoef_proc;
    }/*endif*/


/*========================================================================*/
/* I) Allocate DIIS memory for this step                                  */

    mvec = ( itcount < MVEC_MAX_ATM ? itcount:MVEC_MAX_ATM);
    mp1 = mvec + 1;

    if(itcount==1){
      init_alloc_atm_DIIS_mem(&ipvt,&pmat,&vect,&xst,&yst,&zst,
                              &fxst,&fyst,&fzst,mp1,natm_tot,mvec);
    } else if(itcount > 1 && itcount <= MVEC_MAX_ATM) {
      realloc_atm_DIIS_mem(&ipvt,&pmat,&vect,&xst,&yst,&zst,
                           &fxst,&fyst,&fzst,mp1,natm_tot,mvec);
    } else {
      shift_atm_DIIS_mem(xst,yst,zst,fxst,fyst,fzst,mvec,natm_tot);
    }

/*========================================================================*/
/* II) Get forces                                                         */

    (class->energy_ctrl.iget_full_inter)= 1;
    (class->energy_ctrl.iget_res_inter) = 0;
    (class->energy_ctrl.iget_full_intra)= 1;
    (class->energy_ctrl.iget_res_intra) = 0;

    cp_energy_control(class,bonded,general_data,cp);

/*========================================================================*/
/* III) Compute the inverse Hessian and hit the forces with it            */

 if(myid_state == 0){

    if(min_atm_com_fix_opt==1){
      proj_com_out(natm_tot,clatoms_fx,clatoms_fy,clatoms_fz);
    }/*endif*/

    if(hess_calc >= 2)
      act_hess_inv_on_grad(hess_xx,hess_xy,hess_xz,hess_yy,hess_yz,hess_zz,
                                grad_x,grad_y,grad_z,clatoms_fx,clatoms_fy,clatoms_fz,
                                natm_tot);
    if(hess_calc == 1)
      act_hess_inv_on_grad_diag(hess_xx,hess_xy,hess_xz,hess_yy,hess_yz,hess_zz,
                           grad_x,grad_y,grad_z,clatoms_fx,clatoms_fy,clatoms_fz,
                           natm_tot);

    if(hess_calc == 0){
      for(i=1;i<=natm_tot;i++){
        grad_x[i] = clatoms_fx[i];
        grad_y[i] = clatoms_fy[i];
        grad_z[i] = clatoms_fz[i];
      }/* endfor */
    }/* endif */

/*========================================================================*/
/* III) Get approximate gradients, evolve positions, save new histories   */

    ist=(mvec-1)*natm_tot;
    for(i=1;i<=natm_tot;i++){
      iist = ist + i;
      fxst[iist] = grad_x[i];
      fyst[iist] = grad_y[i];
      fzst[iist] = grad_z[i];
    }
    
    for(i=1;i<=natm_tot;i++){
      iist = ist + i;
      clatoms_x[i] += fxst[iist];
      clatoms_y[i] += fyst[iist];
      clatoms_z[i] += fzst[iist];
    }/* endfor */

    for(i=1;i<=natm_tot;i++){
      iist = ist + i;
      xst[iist] = clatoms_x[i];
      yst[iist] = clatoms_y[i];
      zst[iist] = clatoms_z[i];
    }/* endfor i */
/*==========================================================================*/
/* IV) Construct DIIS linear equations                                    */ 

  setup_atm_DIIS_eqs(pmat,vect,fxst,fyst,fzst,natm_tot,mp1,mvec);

/*==========================================================================*/
/* V) Solve DIIS equations to get DIIS vectors                           */

#ifdef IBM_ESSL
      dgef(&(pmat[1]),&mp1,&mp1,&(ipvt[1]));
      job = 0;
      dges(&(pmat[1]),&mp1,&mp1,&(ipvt[1]),&(vect[1]),&job);
#else
      DGEFA(&(pmat[1]),&mp1,&mp1,&(ipvt[1]),&info);
      job=1;
      DGESL(&(pmat[1]),&mp1,&mp1,&(ipvt[1]),&(vect[1]),&job);
#endif
/*==========================================================================*/
/* VI) Compute new positions from DIIS vectors                             */

   for(i=1;i<=natm_tot;i++){
     clatoms_x[i] = 0.0;
     clatoms_y[i] = 0.0;
     clatoms_z[i] = 0.0;
   }/* endfor i */ 
   for(i=1;i<=mvec;i++){
     ist=(i-1)*natm_tot;
     for(is=1;is<=natm_tot;is++){
       iist = ist+is;
       clatoms_x[is] += vect[i]*xst[iist];
       clatoms_y[is] += vect[i]*yst[iist];
       clatoms_z[is] += vect[i]*zst[iist];
     }/* endfor is */
   }/* endfor i */

 }/* endif myid_state */

 if(np_states > 1){
   Bcast(&(clatoms_x[1]),natm_tot,MPI_DOUBLE,0,comm_states);
   Bcast(&(clatoms_y[1]),natm_tot,MPI_DOUBLE,0,comm_states);
   Bcast(&(clatoms_z[1]),natm_tot,MPI_DOUBLE,0,comm_states);
 }/* endif */

/*==========================================================================*/
/* IV.V) Do a steepest descent step for coefficients                        */ 

   for(is=1;is<=nstate_up;is++) {
     for(i=1;i<=ncoef_up;i++) {
       icoef = i+ioff_up[is];
       cre_up[icoef] += dt*fcre_up[icoef]/cmass[(i+icmoff_up)];
       cim_up[icoef] += dt*fcim_up[icoef]/cmass[(i+icmoff_up)];
     }/*endfor*/
   }/*endfor*/

 if( (cp->cpopts.cp_lsda == 1) && (nstate_dn != 0) ){
   for(is=1;is<=nstate_dn;is++) {
     for(i=1;i<=ncoef_dn;i++) {
       icoef = i+ioff_dn[is];
       cre_dn[icoef] += dt*fcre_dn[icoef]/cmass[(i+icmoff_dn)];
       cim_dn[icoef] += dt*fcim_dn[icoef]/cmass[(i+icmoff_dn)];
     }/*endfor*/
   }/*endfor*/
  }/* endif */

/*==========================================================================*/
/* V) Orthogonalize wave functions                                          */

  orthog_control_cp(cp,ip);

/*==========================================================================*/
/* ii) Free memory and shuffle states                                       */

  cp_shuffle_states(cp,ip);

/*==========================================================================*/
/* VII) Increment iteration counter                                         */

   ++itcount;

/*========================================================================*/
}/* end routine */
/*========================================================================*/

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
void init_alloc_atm_DIIS_mem(int **ipvt,double **pmat,double **vect,
                             double **xst,double **yst,double **zst,
                             double **fxst,double **fyst,double **fzst,
                             int mp1,int natm_tot,int mvec)

/*========================================================================*/
{/*begin routine*/
/*========================================================================*/
/* I) Allocate DIIS arrays */

 (*ipvt) = (int *) cmalloc(mp1*sizeof(int))-1;
 (*pmat) = (double *) cmalloc(mp1*mp1*sizeof(double))-1;
 (*vect) = (double *) cmalloc(mp1*sizeof(double))-1;
 (*xst)  = (double *) cmalloc(natm_tot*mvec*sizeof(double))-1;
 (*yst)  = (double *) cmalloc(natm_tot*mvec*sizeof(double))-1;
 (*zst)  = (double *) cmalloc(natm_tot*mvec*sizeof(double))-1;
 (*fxst)  = (double *) cmalloc(natm_tot*mvec*sizeof(double))-1;
 (*fyst)  = (double *) cmalloc(natm_tot*mvec*sizeof(double))-1;
 (*fzst)  = (double *) cmalloc(natm_tot*mvec*sizeof(double))-1;

/*========================================================================*/
}/* end routine */
/*========================================================================*/


/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
void realloc_atm_DIIS_mem(int **ipvt,double **pmat,double **vect,
                          double **xst,double **yst,double **zst,
                          double **fxst,double **fyst,double **fzst,
                          int mp1,int natm_tot,int mvec)

/*========================================================================*/
{/*begin routine*/
/*========================================================================*/
/* I) Reallocate DIIS arrays */

 (*ipvt) = (int *) crealloc(&((*ipvt)[1]),mp1*sizeof(int))-1;
 (*pmat) = (double *) crealloc(&((*pmat)[1]),mp1*mp1*sizeof(double))-1;
 (*vect) = (double *) crealloc(&((*vect)[1]),mp1*sizeof(double))-1;
 (*xst)  = (double *) crealloc(&((*xst)[1]),natm_tot*mvec*sizeof(double))-1;
 (*yst)  = (double *) crealloc(&((*yst)[1]),natm_tot*mvec*sizeof(double))-1;
 (*zst)  = (double *) crealloc(&((*zst)[1]),natm_tot*mvec*sizeof(double))-1;
 (*fxst)  = (double *) crealloc(&((*fxst)[1]),natm_tot*mvec*sizeof(double))-1;
 (*fyst)  = (double *) crealloc(&((*fyst)[1]),natm_tot*mvec*sizeof(double))-1;
 (*fzst)  = (double *) crealloc(&((*fzst)[1]),natm_tot*mvec*sizeof(double))-1;

/*========================================================================*/
}/* end routine */
/*========================================================================*/

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
void shift_atm_DIIS_mem(double *xst,double *yst,double *zst,
                        double *fxst,double *fyst,double *fzst,
                        int mvec,int natm_tot)

/*========================================================================*/
{/*begin routine*/
/*========================================================================*/
/*  Local variable declarations                                           */

 int ist,ipst,iist,iipst,i,is;

/*========================================================================*/
/* I) Revise history                                                      */

 for(i=1;i<=mvec-1;i++){
   ist=(i-1)*natm_tot;
   ipst=i*natm_tot;
   for(is=1;is<=natm_tot;is++){
     iist = ist + is;
     iipst = ipst + is;
     xst[iist] = xst[iipst];
     yst[iist] = yst[iipst];
     zst[iist] = zst[iipst];
     fxst[iist] = fxst[iipst];
     fyst[iist] = fyst[iipst];
     fzst[iist] = fzst[iipst];
   }/* endfor is */
 }/* endfor i */

/*========================================================================*/
 }/*end routine*/
/*========================================================================*/

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void setup_atm_DIIS_eqs(double *pmat,double *vect,
                        double *fxst,double *fyst,double *fzst,
                        int natm_tot,int mp1,int mvec)

/*========================================================================*/
{/*begin routine*/
/*========================================================================*/
/*  Local variable declarations                                           */

 int i,j,k,l,ist,jst,iii;
 double fovlap;

/*========================================================================*/
/* I) Compute dot products and fill DIIS matrix                           */

 k=0;
 for(i=1;i<=mp1;i++){
   for(j=1;j<=mp1;j++){
     ++k;
     if(i == mp1 || j == mp1){ 
      pmat[k] = 1.0;
     } else {
      ist = (i-1)*natm_tot;
      jst = (j-1)*natm_tot;
      fovlap = 0.0;
      for(l=1;l<=natm_tot;l++){
       fovlap += fxst[(ist+l)]*fxst[(jst+l)]
               + fyst[(ist+l)]*fyst[(jst+l)]
               + fzst[(ist+l)]*fzst[(jst+l)];
      }/* endfor l */
      pmat[k] = fovlap;
     }/* endif */
   }/* endfor j */
 }/*endfor i*/
 pmat[(mp1*mp1)] = 0.0;
 for(i=1;i<=mvec;i++){
  vect[i] = 0.0;
 }/* endfor i */
 vect[mp1] = 1.0;

/*========================================================================*/
 }/*end routine*/
/*========================================================================*/


/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void act_hess_inv_on_grad(double *hess_xx,double *hess_xy,double *hess_xz,
                          double *hess_yy,double *hess_yz,double *hess_zz,
                          double *grad_x,double *grad_y,double *grad_z,
                          double *fx,double *fy,double *fz,int natm_tot)


/*========================================================================*/
{/*begin routine*/
/*========================================================================*/
/*  Local variable declarations                                           */

 int i,ipart,jpart;
 int count=1,indh;
 int job,info=0;
 int dimh;
 int *ipvt;
 int iii;
 int nm0,nm1,nm2;
 double *scr_mat,*scr_vec;

/*========================================================================*/
/* I) Malloc up some local scratch                                        */

 scr_mat = (double *) cmalloc(9*natm_tot*natm_tot*sizeof(double))-1;
 scr_vec = (double *) cmalloc(3*natm_tot*sizeof(double))-1;

/*========================================================================*/
/* II) Pack the scratch matrix and vector                                 */

 for(ipart=1;ipart<=natm_tot;ipart++){
   for(jpart=1;jpart<=natm_tot;jpart++){
      indh = (ipart-1)*natm_tot + jpart;
      scr_mat[count] = hess_xx[indh] + 1.0; count++;
      scr_mat[count] = hess_xy[indh]; count++;
      scr_mat[count] = hess_xz[indh]; count++;
   }
   for(jpart=1;jpart<=natm_tot;jpart++){
      indh = (ipart-1)*natm_tot + jpart;
      scr_mat[count] = hess_xy[indh]; count++;
      scr_mat[count] = hess_yy[indh] + 1.0; count++;
      scr_mat[count] = hess_yz[indh]; count++;
   }
   for(jpart=1;jpart<=natm_tot;jpart++){
      indh = (ipart-1)*natm_tot + jpart;
      scr_mat[count] = hess_xz[indh]; count++;
      scr_mat[count] = hess_yz[indh]; count++;
      scr_mat[count] = hess_zz[indh] + 1.0; count++;
   }
 }
 count=1;
 for(ipart=1;ipart<=natm_tot;ipart++){
   scr_vec[count] = fx[ipart];count++;
   scr_vec[count] = fy[ipart];count++;
   scr_vec[count] = fz[ipart];count++;
 }

/*========================================================================*/
/* III) Solve the linear system to get new gradients                      */

  dimh = 3*natm_tot;
  ipvt = (int *) cmalloc(dimh*sizeof(int))-1;
#ifdef IBM_ESSL
  dgef(&(scr_mat[1]),&dimh,&dimh,&(ipvt[1]));
  job = 0;
  dges(&(scr_mat[1]),&dimh,&dimh,&(ipvt[1]),&(scr_vec[1]),&job);
#else
  DGEFA(&(scr_mat[1]),&dimh,&dimh,&(ipvt[1]),&info);
  job=0;
  DGESL(&(scr_mat[1]),&dimh,&dimh,&(ipvt[1]),&(scr_vec[1]),&job);
#endif


 
/*========================================================================*/
/* IV) Unpack the scratch vector                                          */

 count=1;
 for(ipart=1;ipart<=natm_tot;ipart++){
   grad_x[ipart] = scr_vec[count];count++;
   grad_y[ipart] = scr_vec[count];count++;
   grad_z[ipart] = scr_vec[count];count++;
 }

/*========================================================================*/
/*  Free scratch memory                                                   */

 cfree(&(scr_mat[1]));
 cfree(&(scr_vec[1]));
 cfree(&(ipvt[1]));


/*========================================================================*/
 }/*end routine*/
/*========================================================================*/




/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void act_hess_inv_on_grad_diag(double *hess_xx,double *hess_xy,double *hess_xz,
                               double *hess_yy,double *hess_yz,double *hess_zz,
                               double *grad_x,double *grad_y,double *grad_z,
                               double *fx,double *fy,double *fz,int natm_tot)


/*========================================================================*/
{/*begin routine*/
/*========================================================================*/
/*  Local variable declarations                                           */

 int i,ipart,jpart;
 int count=1,indh;
 int job,info=0;
 int dimh;
 int *ipvt;
 int iii;
 int nm0,nm1,nm2;
 double *scr_mat,*scr_vec;

/*========================================================================*/
/* I) Malloc up some local scratch                                        */

#define DEVELOPMENT_OFF
#ifdef DEVELOPMENT

 scr_mat = (double *) cmalloc(9*sizeof(double))-1;
 scr_vec = (double *) cmalloc(3*sizeof(double))-1;
 dimh = 3;
 ipvt = (int *) cmalloc(dimh*sizeof(int))-1;

/*========================================================================*/
/* II) Pack the scratch matrix and vector                                 */

 for(ipart=1;ipart<=natm_tot;ipart++){

      scr_mat[1] = hess_xx[ipart];
      scr_mat[2] = hess_xy[ipart];
      scr_mat[3] = hess_xz[ipart];
      scr_mat[4] = hess_xy[ipart];
      scr_mat[5] = hess_yy[ipart];
      scr_mat[6] = hess_yz[ipart];
      scr_mat[7] = hess_xz[ipart];
      scr_mat[8] = hess_yz[ipart];
      scr_mat[9] = hess_zz[ipart];
      scr_vec[1] = fx[ipart];
      scr_vec[2] = fy[ipart];
      scr_vec[3] = fz[ipart];

/*========================================================================*/
/* III) Solve the linear system to get new gradients                      */

#ifdef IBM_ESSL
  dgef(&(scr_mat[1]),&dimh,&dimh,&(ipvt[1]));
  job = 0;
  dges(&(scr_mat[1]),&dimh,&dimh,&(ipvt[1]),&(scr_vec[1]),&job);
#else
  DGEFA(&(scr_mat[1]),&dimh,&dimh,&(ipvt[1]),&info);
  job=1;
  DGESL(&(scr_mat[1]),&dimh,&dimh,&(ipvt[1]),&(scr_vec[1]),&job);
#endif


 
/*========================================================================*/
/* IV) Unpack the scratch vector                                          */

   grad_x[ipart] = scr_vec[1];
   grad_y[ipart] = scr_vec[2];
   grad_z[ipart] = scr_vec[3];

 }/* endfor */

/*========================================================================*/
/*  Free scratch memory                                                   */

 cfree(&(scr_mat[1]));
 cfree(&(scr_vec[1]));
 cfree(&(ipvt[1]));

#else

 for(ipart=1;ipart<=natm_tot;ipart++){
   grad_x[ipart] = fx[ipart]/hess_xx[ipart];
   grad_y[ipart] = fy[ipart]/hess_yy[ipart];
   grad_z[ipart] = fz[ipart]/hess_zz[ipart];
 }/* endfor */

#endif

/*========================================================================*/
 }/*end routine*/
/*========================================================================*/

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void displace_atm(CLASS *class,BONDED *bonded,GENERAL_DATA *general_data,
                   CP *cp,double delta,int *displace_index,int *sign_index)

/*========================================================================*/
{/*begin routine*/
/*========================================================================*/
/*             Local variable declarations                                */

#include "../typ_defs/typ_mask.h"


    int ind_atm,irem;
    int displace_index_now;
    int sign_index_now;
    double displace;
    double *clatoms_x    = class->clatoms_pos[1].x;
    double *clatoms_y    = class->clatoms_pos[1].y;
    double *clatoms_z    = class->clatoms_pos[1].z;

/*========================================================================*/
/* I) Determine what kind of displacement                                 */

    displace_index_now = *displace_index;
    sign_index_now     = *sign_index;


    if(sign_index_now == 1) displace = delta;
    if(sign_index_now == -1) displace = -2.0*delta;

/*========================================================================*/
/* II) Displace atoms accordingly                                         */

    ind_atm = (int) (((double) (displace_index_now-1))/3.0) + 1;
    irem = (displace_index_now % 3);


    if(irem == 1) clatoms_x[ind_atm] += displace;
    if(irem == 2) clatoms_y[ind_atm] += displace;
    if(irem == 0) clatoms_z[ind_atm] += displace;


/*========================================================================*/
 }/*end routine*/
/*========================================================================*/


/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void assign_hessian(CLASS *class,BONDED *bonded,GENERAL_DATA *general_data,
                    CP *cp,double delta,int displace_index,int sign_index)

/*========================================================================*/
{/*begin routine*/
/*========================================================================*/
/*             Local variable declarations                                */

#include "../typ_defs/typ_mask.h"


    int i,ind_atm,irem;
    int hess_ind;
    int natm_tot         = class->clatoms_info.natm_tot;
    double *hess_xx      = class->clatoms_pos[1].hess_xx;
    double *hess_xy      = class->clatoms_pos[1].hess_xy;
    double *hess_xz      = class->clatoms_pos[1].hess_xz;
    double *hess_yy      = class->clatoms_pos[1].hess_yy;
    double *hess_yz      = class->clatoms_pos[1].hess_yz;
    double *hess_zz      = class->clatoms_pos[1].hess_zz;
    double *fxp          = class->clatoms_pos[1].vx;
    double *fyp          = class->clatoms_pos[1].vy;
    double *fzp          = class->clatoms_pos[1].vz;
    double *fx           = class->clatoms_pos[1].fx;
    double *fy           = class->clatoms_pos[1].fy;
    double *fz           = class->clatoms_pos[1].fz;

/*========================================================================*/
/* I) Calculate forces                                                   */

    (class->energy_ctrl.iget_full_inter)= 1;
    (class->energy_ctrl.iget_res_inter) = 0;
    (class->energy_ctrl.iget_full_intra)= 1;
    (class->energy_ctrl.iget_res_intra) = 0;
    cp_energy_control(class,bonded,general_data,cp);

/*========================================================================*/
/* II) Save forces only for + displacement                                */

    if(sign_index == 1){

      for(i=1;i<=natm_tot;i++){
        fxp[i] = fx[i];
        fyp[i] = fy[i];
        fzp[i] = fz[i];
      }/* endfor */

    }/* endif */

/*========================================================================*/
/* III) For - displacement, compute a finite difference                    */

    if(sign_index == -1){
      ind_atm = (int) (((double) (displace_index-1))/3.0) + 1;
      irem = (displace_index % 3);
      if(irem == 1) {
	for(i=1;i<=natm_tot;i++){
          hess_ind = (ind_atm-1)*natm_tot + i;
          hess_xx[hess_ind] = (fx[i]-fxp[i])/(2.0*delta);
          hess_xy[hess_ind] = (fy[i]-fyp[i])/(2.0*delta);
          hess_xz[hess_ind] = (fz[i]-fzp[i])/(2.0*delta);
        }/*endfor */
      }/* endif */
      if(irem == 2) {
	for(i=1;i<=natm_tot;i++){
          hess_ind = (ind_atm-1)*natm_tot + i;
          hess_yy[hess_ind] = (fy[i]-fyp[i])/(2.0*delta);
          hess_yz[hess_ind] = (fz[i]-fzp[i])/(2.0*delta);
        }/*endfor */
      }
      if(irem == 0) {
	for(i=1;i<=natm_tot;i++){
          hess_ind = (ind_atm-1)*natm_tot + i;
          hess_zz[hess_ind] = (fz[i]-fzp[i])/(2.0*delta);
        }/*endfor */
      }
    }/* endif sign_index */


/*========================================================================*/
 }/*end routine*/
/*========================================================================*/
