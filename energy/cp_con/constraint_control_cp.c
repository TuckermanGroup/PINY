/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
/*                                                                          */
/*                         PI_MD:                                           */
/*             The future of simulation technology                          */
/*             ------------------------------------                         */
/*                     Module: constraint_control_cp                        */
/*                                                                          */
/* This subprogram controls shake and rattle calls for CP                   */
/*                                                                          */
/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

#include "standard_include.h"
#include "../typ_defs/typedefs_gen.h"
#include "../typ_defs/typedefs_cp.h"
#include "../proto_defs/proto_energy_cpcon_entry.h"
#include "../proto_defs/proto_energy_cpcon_local.h"

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void shake_control_cp(CP *cp,int *iter_shake, double dt,int ip)

/*========================================================================*/
   {/*begin routine*/
/*========================================================================*/
/*             Local variable declarations                                */

 int iii;
/*             Local pointers */

  int cp_norb              = cp->cpopts.cp_norb;
  int cp_lsda              = cp->cpopts.cp_lsda;
  int np_states            = cp->communicate.np_states;
  double c_tolshake        = cp->cpconstrnt.c_tolshake;

  int    *cpcoeffs_ioff_upt = cp->cpcoeffs_info.ioff_upt;
  int    *cpcoeffs_ioff_dnt = cp->cpcoeffs_info.ioff_dnt;
  double *cpcoeffs_cre_up  = cp->cpcoeffs_pos[ip].cre_up;
  double *cpcoeffs_cim_up  = cp->cpcoeffs_pos[ip].cim_up;
  double *cpcoeffs_cre_dn  = cp->cpcoeffs_pos[ip].cre_dn;
  double *cpcoeffs_cim_dn  = cp->cpcoeffs_pos[ip].cim_dn;
  double *cpcoeffs_vcre_up = cp->cpcoeffs_pos[ip].vcre_up;
  double *cpcoeffs_vcim_up = cp->cpcoeffs_pos[ip].vcim_up;
  double *cpcoeffs_vcre_dn = cp->cpcoeffs_pos[ip].vcre_dn;
  double *cpcoeffs_vcim_dn = cp->cpcoeffs_pos[ip].vcim_dn;
  int icoef_form_up        = cp->cpcoeffs_pos[ip].icoef_form_up;
  int icoef_form_old_up    = cp->cpcoeffs_pos[ip].icoef_form_up;
  int icoef_orth_up        = cp->cpcoeffs_pos[ip].icoef_orth_up;
  int ivcoef_form_up       = cp->cpcoeffs_pos[ip].ivcoef_form_up;
  int icoef_form_dn        = cp->cpcoeffs_pos[ip].icoef_form_dn;
  int icoef_form_old_dn    = cp->cpcoeffs_pos[ip].icoef_form_dn;
  int icoef_orth_dn        = cp->cpcoeffs_pos[ip].icoef_orth_dn;
  int ivcoef_form_dn       = cp->cpcoeffs_pos[ip].ivcoef_form_dn;
  int nstate_up            = cp->cpcoeffs_info.nstate_up;
  int nstate_dn            = cp->cpcoeffs_info.nstate_dn;

  double *cpcoeffs_cmass   = cp->cpcoeffs_info.cmass;
  int icmass_unif          = cp->cpcoeffs_info.icmass_unif;

  double *occ_up           = cp->cpopts.occ_up;
  double *occ_dn           = cp->cpopts.occ_dn;
  double *rocc_sum_up      = cp->cpopts.rocc_sum_up;
  double *rocc_sum_dn      = cp->cpopts.rocc_sum_dn;

  double *cpscr_cre_up     = cp->cpscr.cpscr_wave.cre_up;
  double *cpscr_cim_up     = cp->cpscr.cpscr_wave.cim_up;
  double *cpscr_cre_dn     = cp->cpscr.cpscr_wave.cre_dn;
  double *cpscr_cim_dn     = cp->cpscr.cpscr_wave.cim_dn;


/*========================================================================*/
/* I) Norb off or norb with full ortho */

 if(cp_norb <= 1) {

  /*--------------------------------------------------------------------*/
  /* i) Uniform mass                                                    */
  if(icmass_unif == 1) {

    cp_shake(cpcoeffs_cre_up,cpcoeffs_cim_up,icoef_form_up,icoef_orth_up,
             cpcoeffs_vcre_up,cpcoeffs_vcim_up,ivcoef_form_up,
             cpscr_cre_up,cpscr_cim_up,icoef_form_old_up,
             c_tolshake,dt,&(cp->cpscr.cpscr_ovmat),
             cpcoeffs_ioff_upt,occ_up,rocc_sum_up,cp_norb,
             iter_shake,&(cp->cp_comm_state_pkg_up));

    if(cp_lsda == 1 && nstate_dn != 0) {

      cp_shake(cpcoeffs_cre_dn,cpcoeffs_cim_dn,icoef_form_dn,icoef_orth_dn,
               cpcoeffs_vcre_dn,cpcoeffs_vcim_dn,ivcoef_form_dn,
               cpscr_cre_dn,cpscr_cim_dn,icoef_form_old_dn,
               c_tolshake,dt,&(cp->cpscr.cpscr_ovmat),
               cpcoeffs_ioff_dnt,occ_dn,rocc_sum_dn,cp_norb,
               iter_shake,&(cp->cp_comm_state_pkg_dn));

    }/* endif: lsda */

  /*--------------------------------------------------------------------*/
  /* ii) Non-Uniform mass                                               */
  } else {

    cp_shake_mass(cpcoeffs_cre_up,cpcoeffs_cim_up,icoef_form_up,icoef_orth_up,
                  cpcoeffs_vcre_up,cpcoeffs_vcim_up,ivcoef_form_up,
                  cpscr_cre_up,cpscr_cim_up,icoef_form_old_up,
                  cpcoeffs_cmass,
                  c_tolshake,dt,&(cp->cpscr.cpscr_ovmat),
                  cpcoeffs_ioff_upt,occ_up,rocc_sum_up,cp_norb,
                  iter_shake,&(cp->cp_comm_state_pkg_up));

    if(cp_lsda == 1 && nstate_dn != 0){

     cp_shake_mass(cpcoeffs_cre_dn,cpcoeffs_cim_dn,icoef_form_dn,icoef_orth_dn,
                   cpcoeffs_vcre_dn,cpcoeffs_vcim_dn,ivcoef_form_dn,
                   cpscr_cre_dn,cpscr_cim_dn,icoef_form_old_dn,
                   cpcoeffs_cmass,
                   c_tolshake,dt,&(cp->cpscr.cpscr_ovmat),
                   cpcoeffs_ioff_dnt,occ_dn,rocc_sum_dn,cp_norb,
                   iter_shake,&(cp->cp_comm_state_pkg_dn));

    }/* endif: lsda */

  }/* endif: icmass_unif */

 }/* endif cp_norb */

/*========================================================================*/
/* II) Norb with norm only */

 if(cp_norb== 2){

  cp_shake_norb(cpcoeffs_cre_up,cpcoeffs_cim_up,icoef_form_up,icoef_orth_up,
                cpcoeffs_vcre_up,cpcoeffs_vcim_up,ivcoef_form_up,
                cpscr_cre_up,cpscr_cim_up,icoef_form_old_up,
                dt,cpcoeffs_ioff_upt,occ_up,cpcoeffs_cmass,icmass_unif,
                &(cp->cpscr.cpscr_ovmat),
                &(cp->cp_comm_state_pkg_dn));
                
  if(cp_lsda == 1 && nstate_dn != 0) {

    cp_shake_norb(cpcoeffs_cre_dn,cpcoeffs_cim_dn,icoef_form_dn,icoef_orth_dn,
                  cpcoeffs_vcre_dn,cpcoeffs_vcim_dn,ivcoef_form_dn,
                  cpscr_cre_dn,cpscr_cim_dn,icoef_form_old_dn,
                  dt,cpcoeffs_ioff_dnt,occ_dn,cpcoeffs_cmass,icmass_unif,
                  &(cp->cpscr.cpscr_ovmat),
                  &(cp->cp_comm_state_pkg_dn));
                
  }/*endif:lsda*/

 }/* endif cp_norb norm only*/
   
/*========================================================================*/
   }/* end routine */
/*========================================================================*/





/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void rattle_control_cp(CP *cp,int *iter_rattle,double dt,int ip)

/*========================================================================*/
{/*begin routine*/
/*========================================================================*/
/*             Local variable declarations                                */

/* Local pointers */
  int iii;

/*             Local pointers */

  int cp_norb              = cp->cpopts.cp_norb;
  int cp_lsda              = cp->cpopts.cp_lsda;
  int np_states            = cp->communicate.np_states;
  double c_tolratl         = cp->cpconstrnt.c_tolratl;

  int    *cpcoeffs_ioff_upt = cp->cpcoeffs_info.ioff_upt;
  int    *cpcoeffs_ioff_dnt = cp->cpcoeffs_info.ioff_dnt;
  double *cpcoeffs_cre_up  = cp->cpcoeffs_pos[ip].cre_up;
  double *cpcoeffs_cim_up  = cp->cpcoeffs_pos[ip].cim_up;
  double *cpcoeffs_cre_dn  = cp->cpcoeffs_pos[ip].cre_dn;
  double *cpcoeffs_cim_dn  = cp->cpcoeffs_pos[ip].cim_dn;
  double *cpcoeffs_vcre_up = cp->cpcoeffs_pos[ip].vcre_up;
  double *cpcoeffs_vcim_up = cp->cpcoeffs_pos[ip].vcim_up;
  double *cpcoeffs_vcre_dn = cp->cpcoeffs_pos[ip].vcre_dn;
  double *cpcoeffs_vcim_dn = cp->cpcoeffs_pos[ip].vcim_dn;
  int icoef_form_up        = cp->cpcoeffs_pos[ip].icoef_form_up;
  int icoef_orth_up        = cp->cpcoeffs_pos[ip].icoef_orth_up;
  int ivcoef_form_up       = cp->cpcoeffs_pos[ip].ivcoef_form_up;
  int icoef_form_dn        = cp->cpcoeffs_pos[ip].icoef_form_dn;
  int icoef_orth_dn        = cp->cpcoeffs_pos[ip].icoef_orth_dn;
  int ivcoef_form_dn       = cp->cpcoeffs_pos[ip].ivcoef_form_dn;
  int nstate_up            = cp->cpcoeffs_info.nstate_up;
  int nstate_dn            = cp->cpcoeffs_info.nstate_dn;


  double *cpcoeffs_cmass   = cp->cpcoeffs_info.cmass;
  int icmass_unif          = cp->cpcoeffs_info.icmass_unif;

  double *occ_up           = cp->cpopts.occ_up;
  double *occ_dn           = cp->cpopts.occ_dn;
  double *rocc_sum_up      = cp->cpopts.rocc_sum_up;
  double *rocc_sum_dn      = cp->cpopts.rocc_sum_dn;

  double *cpscr_cre_up     = cp->cpscr.cpscr_wave.cre_up;
  double *cpscr_cim_up     = cp->cpscr.cpscr_wave.cim_up;
  double *cpscr_cre_dn     = cp->cpscr.cpscr_wave.cre_dn;
  double *cpscr_cim_dn     = cp->cpscr.cpscr_wave.cim_dn;


/*========================================================================*/
/* I) Norb off or norb with full ortho */

 if(cp_norb <= 1) {

  /*--------------------------------------------------------------------*/
  /* i) Uniform mass                                                    */
  if(icmass_unif == 1) {

   cp_rattle(cpcoeffs_cre_up,cpcoeffs_cim_up,icoef_form_up,icoef_orth_up,
             cpcoeffs_vcre_up,cpcoeffs_vcim_up,ivcoef_form_up,
             &(cp->cpscr.cpscr_ovmat),
             cpcoeffs_ioff_upt,rocc_sum_up,&(cp->cp_comm_state_pkg_up),
             cp_norb);

   if(cp_lsda== 1 && nstate_dn != 0){

     cp_rattle(cpcoeffs_cre_dn,cpcoeffs_cim_dn,icoef_form_dn,icoef_orth_dn,
               cpcoeffs_vcre_dn,cpcoeffs_vcim_dn,ivcoef_form_dn,
               &(cp->cpscr.cpscr_ovmat),
               cpcoeffs_ioff_dnt,rocc_sum_dn,&(cp->cp_comm_state_pkg_dn),
               cp_norb);

   }/* endif lsda */
  /*--------------------------------------------------------------------*/
  /* i) Non-Uniform mass                                                */
  } else {

    cp_rattle_mass(cpcoeffs_cre_up,cpcoeffs_cim_up,icoef_form_up,icoef_orth_up,
                   cpcoeffs_vcre_up,cpcoeffs_vcim_up,ivcoef_form_up,
                   cpscr_cre_up,cpscr_cim_up,c_tolratl,
                   cpcoeffs_cmass,&(cp->cpscr.cpscr_ovmat),
                   cpcoeffs_ioff_upt,occ_up,rocc_sum_up,
                   iter_rattle,&(cp->cp_comm_state_pkg_up),cp_norb);
   if(cp_lsda== 1 && nstate_dn != 0){

    cp_rattle_mass(cpcoeffs_cre_dn,cpcoeffs_cim_dn,icoef_form_dn,icoef_orth_dn,
                   cpcoeffs_vcre_dn,cpcoeffs_vcim_dn,ivcoef_form_dn,
                   cpscr_cre_dn,cpscr_cim_dn,c_tolratl,
                   cpcoeffs_cmass,&(cp->cpscr.cpscr_ovmat),
                   cpcoeffs_ioff_dnt,occ_dn,rocc_sum_dn,
                   iter_rattle,&(cp->cp_comm_state_pkg_dn),cp_norb);


   }/* endif lsda */

  }/* endif icmass_unif */

 }/* endif cp_norb */

/*========================================================================*/
/* I) Norb with */

 if(cp_norb== 2){

  cp_rattle_norb(cpcoeffs_cre_up,cpcoeffs_cim_up,icoef_form_up,icoef_orth_up,
                 cpcoeffs_vcre_up,cpcoeffs_vcim_up,ivcoef_form_up,
                 cpscr_cre_up,cpscr_cim_up,
                 cpcoeffs_ioff_upt,occ_up,cpcoeffs_cmass,icmass_unif,
                 &(cp->cpscr.cpscr_ovmat),
                 &(cp->cp_comm_state_pkg_up));

  if(cp_lsda== 1 && nstate_dn != 0) {

    cp_rattle_norb(cpcoeffs_cre_dn,cpcoeffs_cim_dn,icoef_form_dn,icoef_orth_dn,
                   cpcoeffs_vcre_dn,cpcoeffs_vcim_dn,ivcoef_form_dn,
                   cpscr_cre_dn,cpscr_cim_dn,
                   cpcoeffs_ioff_dnt,occ_dn,cpcoeffs_cmass,icmass_unif,
                   &(cp->cpscr.cpscr_ovmat),
                   &(cp->cp_comm_state_pkg_dn));

  }/*endif: lsda */

 }/*endif: cp_norb */

/*========================================================================*/
   }/* end routine */
/*========================================================================*/










