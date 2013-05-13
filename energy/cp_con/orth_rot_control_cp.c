/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
/*                                                                          */
/*                         PI_MD:                                           */
/*             The future of simulation technology                          */
/*             ------------------------------------                         */
/*                     Module: cp_norb.c                                    */
/*                                                                          */
/* File contains functions necessary for computing norb force               */
/*                                                                          */
/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/


#include "standard_include.h"
#include "../typ_defs/typedefs_gen.h"
#include "../typ_defs/typedefs_cp.h"
#include "../proto_defs/proto_energy_cpcon_local.h"
#include "../proto_defs/proto_math.h"
#include "../proto_defs/proto_friend_lib_entry.h"
#include "../proto_defs/proto_communicate_wrappers.h"



/*==========================================================================*/
/*  Control which orthogonalization scheme is used                          */
/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void orthog_control_cp(CP *cp,int ip)

/*========================================================================*/
{/*begin routine*/
/*========================================================================*/
/*             Local variable declarations                                */

#include "../typ_defs/typ_mask.h"

 int iii;
 int icoef_orth;
 double max_off_diag=0.0,max_diag=0.0;

/* Local pointers */

  int nstate_dn            = cp->cpcoeffs_info.nstate_dn;
  int cp_gs                = cp->cpopts.cp_gs;
  int cp_low               = cp->cpopts.cp_low; 
  int cp_lsda              = cp->cpopts.cp_lsda;
  int cp_norm              = cp->cpopts.cp_normalize;
  int np_states            = cp->communicate.np_states;

  int icoef_form_up        = cp->cpcoeffs_pos[ip].icoef_form_up;
  double *cpcoeffs_cre_up  = cp->cpcoeffs_pos[ip].cre_up;
  double *cpcoeffs_cim_up  = cp->cpcoeffs_pos[ip].cim_up;
  int icoef_form_dn        = cp->cpcoeffs_pos[ip].icoef_form_dn;
  double *cpcoeffs_cre_dn  = cp->cpcoeffs_pos[ip].cre_dn;
  double *cpcoeffs_cim_dn  = cp->cpcoeffs_pos[ip].cim_dn;

  double *cpscr_cre_up     = cp->cpscr.cpscr_wave.cre_up;
  double *cpscr_cim_up     = cp->cpscr.cpscr_wave.cim_up;
  double *cpscr_cre_dn     = cp->cpscr.cpscr_wave.cre_dn;
  double *cpscr_cim_dn     = cp->cpscr.cpscr_wave.cim_dn;
  double *occ_up           = cp->cpopts.occ_up;
  double *occ_dn           = cp->cpopts.occ_dn;

  double *omat             = cp->cpscr.cpscr_ovmat.ovlap7;
  double *omat_tmp         = cp->cpscr.cpscr_ovmat.ovlap8;
  double *anorm            = cp->cpscr.cpscr_ovmat.state_vec1;
  double *anorm_tmp        = cp->cpscr.cpscr_ovmat.state_vec2;
  double *norbmat          = cp->cpscr.cpscr_ovmat.ovlap6;
  /*WARNING CPSCR_OVMAT ALSO PASSED TO ROUTINE*/
  double *norbmati         = cp->cpscr.cpscr_ovmat.ovlap7; 
  double *ovmat_eigv       = cp->cpscr.cpscr_ovmat.ovlap8; 


  int *ioff_upt            = cp->cpcoeffs_info.ioff_upt;
  int *ioff_dnt            = cp->cpcoeffs_info.ioff_dnt;

  double *rocc_sum_up      = cp->cpopts.rocc_sum_up;
  double *rocc_sum_dn      = cp->cpopts.rocc_sum_dn;

/*========================================================================*/
/*  Up states */

  if(cp_gs== 1) {
   control_cp_gram_schmidt(cpcoeffs_cre_up,cpcoeffs_cim_up,icoef_form_up,
                           occ_up,omat,omat_tmp,ioff_upt,
                           &(cp->cp_comm_state_pkg_up));
  }/*endif:gs*/

  if(cp_low== 1) {
    icoef_orth = 0;
    cp_rotate_coef_ortho(cpcoeffs_cre_up,cpcoeffs_cim_up,
                         icoef_form_up,&icoef_orth,
                         norbmat,norbmati,ovmat_eigv,
                         cpscr_cre_up,cpscr_cim_up,
                         occ_up,ioff_upt,
                         &max_off_diag,&max_diag,
                         &(cp->cpscr.cpscr_ovmat),
                         &(cp->cp_comm_state_pkg_up));

    rotate_occ_shuffle(cpcoeffs_cre_up,cpcoeffs_cim_up,
                       icoef_form_up,omat,omat_tmp,
                       occ_up,rocc_sum_up,ioff_upt,
                       &(cp->cp_comm_state_pkg_up)); 
  }/*endif:lowdin*/

  if(cp_norm== 1) {
    cp_normalize(cpcoeffs_cre_up,cpcoeffs_cim_up,icoef_form_up,occ_up,
                 ioff_upt,anorm,anorm_tmp,&(cp->cp_comm_state_pkg_up));
  }/*endif:norm*/


/*========================================================================*/
/*  Down states */

  if( (cp_lsda == 1) && (nstate_dn != 0) ){

    if(cp_gs== 1) {
      control_cp_gram_schmidt(cpcoeffs_cre_dn,cpcoeffs_cim_dn,icoef_form_dn,
                              occ_dn,omat,omat_tmp,ioff_dnt,
                              &(cp->cp_comm_state_pkg_dn));
    }/*endif:gs*/

    if(cp_low== 1) {
      icoef_orth = 0;
      cp_rotate_coef_ortho(cpcoeffs_cre_dn,cpcoeffs_cim_dn,
                           icoef_form_dn,&icoef_orth,
                           norbmat,norbmati,ovmat_eigv,
                           cpscr_cre_dn,cpscr_cim_dn,
                           occ_dn,ioff_dnt,
                           &max_off_diag,&max_diag,
                           &(cp->cpscr.cpscr_ovmat),
                           &(cp->cp_comm_state_pkg_dn));
      rotate_occ_shuffle(cpcoeffs_cre_dn,cpcoeffs_cim_dn,
                         icoef_form_dn,omat,omat_tmp,
                         occ_dn,rocc_sum_dn,ioff_dnt,
                         &(cp->cp_comm_state_pkg_dn)); 
    }/*endif:low*/

    if(cp_norm== 1) {
     cp_normalize(cpcoeffs_cre_dn,cpcoeffs_cim_dn,icoef_form_dn,occ_dn,
                  ioff_dnt,anorm,anorm_tmp,&(cp->cp_comm_state_pkg_dn));
    }/*endif:norm*/

  }/* endif lsda */

/*========================================================================*/
   }/* end routine */
/*========================================================================*/




/*==========================================================================*/
/*  Rotate into basis in which S is diagonal or to the KS basis             */
/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
void cp_rotate_control(CP *cp,int itime,int ntime,STAT_AVG *stat_avg)

/*==========================================================================*/
/*        Begin Routine */
  {/*begin routine */
/*==========================================================================*/
/*        Local Variable Declarations */
#include "../typ_defs/typ_mask.h"
  int ip,iii;
  double max_off_diag_tmp,max_diag_tmp;

/*       Local pointers               */
  int pi_beads_proc          = cp->cpcoeffs_info.pi_beads_proc;
  int pi_beads               = cp->cpcoeffs_info.pi_beads;
  int nstate_dn              = cp->cpcoeffs_info.nstate_dn;
  double *count_diag_srot    = &(stat_avg->count_diag_srot);
  double *max_off_diag       = &(cp->cpcoeffs_info.max_off_diag);
  double *max_diag           = &(cp->cpcoeffs_info.max_diag);
  CPCOEFFS_POS *cpcoeffs_pos = cp->cpcoeffs_pos;
  int cp_norb                = cp->cpopts.cp_norb;
  int np_beads               = cp->communicate.np_beads;
  int np_states              = cp->communicate.np_states;
  int num_proc               = cp->communicate.np;
  int myid                   = cp->communicate.myid;
  MPI_Comm comm_beads        = cp->communicate.comm_beads;
  MPI_Comm comm_states       = cp->communicate.comm_states;
  MPI_Comm world             = cp->communicate.world;
  double *ovmat              = cp->cpscr.cpscr_ovmat.ovlap1;
  double *ovmat_scr          = cp->cpscr.cpscr_ovmat.ovlap2;

  double *anorm            = cp->cpscr.cpscr_ovmat.state_vec1;
  double *anorm_tmp        = cp->cpscr.cpscr_ovmat.state_vec2;

  double kseig_sum_up,kseig_sum_dn,kseig_sum_tot;
  double elec_etot_tmp      = stat_avg->cp_ehart
                            + stat_avg->cp_exc
                            + stat_avg->cp_eext
                            + stat_avg->cp_eke
                            + stat_avg->cp_enl;

  double elec_ecorr_tmp     = stat_avg->cp_exc - stat_avg->cp_ehart
                            - stat_avg->cp_muxc;
  double elec_etot=0.0,elec_ecorr=0.0;
  double ks_offset=0.0,ks_offset_tmp;

/*==========================================================================*/
/* 0) Calculate some simple stuff                                           */

  if(cp_norb>0){

     (*max_off_diag)  = 0.0;
     for(ip=1;ip<=pi_beads_proc;ip++){
       (*max_off_diag)  = MAX((*max_off_diag),cpcoeffs_pos[ip].max_off_diag);
     }/*endfor*/
     if(np_beads>1){
        max_off_diag_tmp = (*max_off_diag);

        Allreduce(&(max_off_diag_tmp),max_off_diag,1,MPI_DOUBLE,MPI_MAX,0,
                                     world);

     }/*endif*/

     (*max_diag)  = 0.0;
     for(ip=1;ip<=pi_beads_proc;ip++){
       (*max_diag)  = MAX((*max_diag),cpcoeffs_pos[ip].max_diag);
     }/*endfor*/
     if(np_beads>1){
       max_diag_tmp = (*max_diag);
       Allreduce(&(max_diag_tmp),max_diag,1,MPI_DOUBLE,MPI_MAX,0,
                 world);
     }/*endif*/
  }/*endif*/

/*==========================================================================*/
/* I) Rotate CP variables if necessary                                      */

  if(cp_norb>0){

   if( ((*max_off_diag) > cp->cpconstrnt.c_tolnorb)|| (itime ==ntime)  ){

     (*count_diag_srot)++;
     for(ip=1;ip<=pi_beads_proc;ip++){

        cp_rotate_all(cp->cpcoeffs_pos[ip].cre_up,
                      cp->cpcoeffs_pos[ip].cim_up,
                      cp->cpcoeffs_pos[ip].icoef_form_up,
                      cp->cpcoeffs_pos[ip].vcre_up,
                      cp->cpcoeffs_pos[ip].vcim_up,
                      cp->cpcoeffs_pos[ip].ivcoef_form_up,
                      cp->cpcoeffs_pos[ip].fcre_up,
                      cp->cpcoeffs_pos[ip].fcim_up,
                      cp->cpcoeffs_pos[ip].ifcoef_form_up,
                      cp->cpcoeffs_pos[ip].ovmat_eigv_up,
                      cp->cpcoeffs_info.ioff_upt,
                      cp->cpscr.cpscr_wave.cre_up,
                      cp->cpscr.cpscr_wave.cim_up,
                      &(cp->cp_comm_state_pkg_up)); 
        if(cp_norb == 2){

          /* Reset the norm to the occupation number */
          /*cp_normalize(cp->cpcoeffs_pos[ip].cre_up,
                       cp->cpcoeffs_pos[ip].cim_up,
                       cp->cpcoeffs_pos[ip].icoef_form_up,
                       cp->cpopts.occ_up,
                       cp->cpcoeffs_info.ioff_upt,
                       anorm,anorm_tmp,&(cp->cp_comm_state_pkg_up));*/

          cp_normalize_all(cp->cpcoeffs_pos[ip].cre_up,
                           cp->cpcoeffs_pos[ip].cim_up,
                           cp->cpcoeffs_pos[ip].icoef_form_up,
                           cp->cpcoeffs_pos[ip].fcre_up,
                           cp->cpcoeffs_pos[ip].fcim_up,
                           cp->cpcoeffs_pos[ip].ifcoef_form_up,
                           cp->cpopts.occ_up,
                           cp->cpcoeffs_info.ioff_upt,
                           anorm,anorm_tmp,&(cp->cp_comm_state_pkg_up));

          /* Add the constraint to the coefficient velocities */ 
          cp_rattle_norb(cp->cpcoeffs_pos[ip].cre_up,
                         cp->cpcoeffs_pos[ip].cim_up,
                         cp->cpcoeffs_pos[ip].icoef_form_up,
                         cp->cpcoeffs_pos[ip].icoef_orth_up,
                         cp->cpcoeffs_pos[ip].vcre_up,
                         cp->cpcoeffs_pos[ip].vcim_up,
                         cp->cpcoeffs_pos[ip].ivcoef_form_up,
                         cp->cpscr.cpscr_wave.cre_up,
                         cp->cpscr.cpscr_wave.cim_up,
                         cp->cpcoeffs_info.ioff_upt,
                         cp->cpopts.occ_up,
                         cp->cpcoeffs_info.cmass,
                         cp->cpcoeffs_info.icmass_unif,
                         &(cp->cpscr.cpscr_ovmat),&(cp->cp_comm_state_pkg_up));
	}/*endif:norb=2*/
        if( (cp->cpopts.cp_lsda == 1) && (nstate_dn != 0) ){
           cp_rotate_all(cp->cpcoeffs_pos[ip].cre_dn,
                         cp->cpcoeffs_pos[ip].cim_dn,
                         cp->cpcoeffs_pos[ip].icoef_form_dn,
                         cp->cpcoeffs_pos[ip].vcre_dn,
                         cp->cpcoeffs_pos[ip].vcim_dn,
                         cp->cpcoeffs_pos[ip].ivcoef_form_dn,
                         cp->cpcoeffs_pos[ip].fcre_dn,
                         cp->cpcoeffs_pos[ip].fcim_dn,
                         cp->cpcoeffs_pos[ip].ifcoef_form_dn,
                         cp->cpcoeffs_pos[ip].ovmat_eigv_dn,
                         cp->cpcoeffs_info.ioff_dnt,
                         cp->cpscr.cpscr_wave.cre_dn,
                         cp->cpscr.cpscr_wave.cim_dn,
                         &(cp->cp_comm_state_pkg_dn)); 
          if(cp_norb == 2){

            /* Reset the norm to the occupation number */
            cp_normalize(cp->cpcoeffs_pos[ip].cre_dn,
                         cp->cpcoeffs_pos[ip].cim_dn,
                         cp->cpcoeffs_pos[ip].icoef_form_dn,
                         cp->cpopts.occ_dn,
                         cp->cpcoeffs_info.ioff_dnt,
                         anorm,anorm_tmp,&(cp->cp_comm_state_pkg_dn));

            /* Add the constraint to the coefficient velocities */
            cp_rattle_norb(cp->cpcoeffs_pos[ip].cre_dn,
                           cp->cpcoeffs_pos[ip].cim_dn,
                           cp->cpcoeffs_pos[ip].icoef_form_dn,
                           cp->cpcoeffs_pos[ip].icoef_orth_dn,
                           cp->cpcoeffs_pos[ip].vcre_dn,
                           cp->cpcoeffs_pos[ip].vcim_dn,
                           cp->cpcoeffs_pos[ip].ivcoef_form_dn,
                           cp->cpscr.cpscr_wave.cre_dn,
                           cp->cpscr.cpscr_wave.cim_dn,
                           cp->cpcoeffs_info.ioff_dnt,
                           cp->cpopts.occ_dn,
                           cp->cpcoeffs_info.cmass,
                           cp->cpcoeffs_info.icmass_unif,
                           &(cp->cpscr.cpscr_ovmat),&(cp->cp_comm_state_pkg_dn));
          }/*endif:norb=2*/

        }/*endif:lsda*/

     }/*endfor:ip*/

   }/*endif:norb rotation time*/

  }/*endif:norb*/

/*==========================================================================*/
/* II) Rotate onto Kohn-Sham Basis if requested                             */

  if((cp->cpcoeffs_info.ks_rot_on == 1)&&
     (itime % cp->cpcoeffs_info.n_ks_rot)==0){

     kseig_sum_tot = 0.0;
     for(ip=1;ip<=pi_beads_proc;ip++){

       cp_rotate_ks_basis(cp->cpcoeffs_pos[ip].cre_up,
                          cp->cpcoeffs_pos[ip].cim_up,
                          cp->cpcoeffs_pos[ip].icoef_form_up,
                          cp->cpcoeffs_pos[ip].icoef_orth_up,
                          cp->cpcoeffs_pos[ip].vcre_up,
                          cp->cpcoeffs_pos[ip].vcim_up,
                          cp->cpcoeffs_pos[ip].ivcoef_form_up,
                          cp->cpcoeffs_pos[ip].ivcoef_orth_up,
                          cp->cpcoeffs_pos[ip].fcre_up,
                          cp->cpcoeffs_pos[ip].fcim_up,
                          cp->cpcoeffs_pos[ip].ifcoef_form_up,
                          cp->cpcoeffs_pos[ip].ifcoef_orth_up,
                          cp->cpcoeffs_info.ioff_upt,
                          cp->cpopts.occ_up,cp->cpopts.rocc_sum_up,
                          cp->cpcoeffs_pos[ip].ksmat_up,
                          cp->cpcoeffs_pos[ip].ksmat_eig_up,
                          cp->cpscr.cpscr_wave.cre_up,
                          cp->cpscr.cpscr_wave.cim_up,
                          &(cp->cpscr.cpscr_ovmat),
                          &(cp->cp_comm_state_pkg_up),&kseig_sum_up);
       kseig_sum_tot += kseig_sum_up;
       if( (cp->cpopts.cp_lsda == 1) && (nstate_dn != 0) ){
         cp_rotate_ks_basis(cp->cpcoeffs_pos[ip].cre_dn,
                            cp->cpcoeffs_pos[ip].cim_dn,
                            cp->cpcoeffs_pos[ip].icoef_form_dn,
                            cp->cpcoeffs_pos[ip].icoef_orth_dn,
                            cp->cpcoeffs_pos[ip].vcre_dn,
                            cp->cpcoeffs_pos[ip].vcim_dn,
                            cp->cpcoeffs_pos[ip].ivcoef_form_dn,
                            cp->cpcoeffs_pos[ip].ivcoef_orth_dn,
                            cp->cpcoeffs_pos[ip].fcre_dn,
                            cp->cpcoeffs_pos[ip].fcim_dn,
                            cp->cpcoeffs_pos[ip].ifcoef_form_dn,
                            cp->cpcoeffs_pos[ip].ifcoef_orth_dn,
                            cp->cpcoeffs_info.ioff_dnt,
                            cp->cpopts.occ_dn,cp->cpopts.rocc_sum_dn,
                            cp->cpcoeffs_pos[ip].ksmat_dn,
                            cp->cpcoeffs_pos[ip].ksmat_eig_dn,
                            cp->cpscr.cpscr_wave.cre_dn,
                            cp->cpscr.cpscr_wave.cim_dn,
                            &(cp->cpscr.cpscr_ovmat),
                            &(cp->cp_comm_state_pkg_dn),&kseig_sum_dn);
        kseig_sum_tot += kseig_sum_dn;
       }/* endif LSDA */

       if(np_states > 1){
          ks_offset_tmp = cp->cpcoeffs_pos[ip].ks_offset;
          Allreduce(&ks_offset_tmp,&ks_offset,1,MPI_DOUBLE,MPI_SUM,0,
                    comm_states);
          cp->cpcoeffs_pos[ip].ks_offset = ks_offset;
       }/* endif */

     }/* endfor ip */
     kseig_sum_tot /= ((double) pi_beads);

    if(num_proc > 1){
      Allreduce(&elec_etot_tmp,&elec_etot,1,MPI_DOUBLE,MPI_SUM,0,world);
      Allreduce(&elec_ecorr_tmp,&elec_ecorr,1,MPI_DOUBLE,MPI_SUM,0,world);
   } else {
    elec_etot = elec_etot_tmp;
    elec_ecorr = elec_ecorr_tmp;
   } /* endif */
   if(myid == 0){
     printf("KS rotation total energy:  True  %.10g  KS %.10g\n",
             elec_etot,(kseig_sum_tot+elec_ecorr));
     PRINT_LINE_DASH
   }/* endif */
  }/* endif KS_ROT_ON */

/*--------------------------------------------------------------------------*/
    }/*end routine */
/*==========================================================================*/










