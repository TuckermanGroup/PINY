/*==========================================================================*/
/*CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC*/
/*==========================================================================*/
/*                                                                          */
/*                         PI_MD:                                           */
/*             The future of simulation technology                          */
/*             ------------------------------------                         */
/*                   Module: elec_energy_control.c                          */
/*                                                                          */
/* This routine calls the required force and PE routines                    */
/*                                                                          */
/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

#include "standard_include.h"
#include "../typ_defs/typedefs_gen.h"
#include "../typ_defs/typedefs_class.h"
#include "../typ_defs/typedefs_bnd.h"
#include "../typ_defs/typedefs_cp.h"
#include "../proto_defs/proto_energy_cp_local.h"
#include "../proto_defs/proto_energy_cpcon_local.h"
#include "../proto_defs/proto_friend_lib_entry.h"
#include "../proto_defs/proto_math.h"
#include "../proto_defs/proto_communicate_wrappers.h"

#define TIME_CP_OFF
 
/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void cp_ks_energy_ctrl(CP *cp,int ip_now,EWALD *ewald,EWD_SCR *ewd_scr,
                       CELL *cell,CLATOMS_INFO *clatoms_info,
                       CLATOMS_POS *clatoms_pos,ATOMMAPS *atommaps,
                       STAT_AVG *stat_avg,PTENS *ptens,SIMOPTS *simopts,
                       FOR_SCR *for_scr)

/*==========================================================================*/
/*         Begin Routine                                                    */
   {/*Begin Routine*/
/*=======================================================================*/
/*      Local variables                                                  */

    int cp_para_opt = cp->cpopts.cp_para_opt;

/*=======================================================================*/
/* Use the specified parallelization scheme                              */


  switch(cp_para_opt){

    case 0:  /* hybrid */
     cp_ks_energy_hybrid(cp,ip_now,ewald,ewd_scr,cell,clatoms_info,
                         clatoms_pos,atommaps,stat_avg,ptens,simopts,for_scr);
     break;

    case 1: /* full_g */
      cp_ks_energy_full_g(cp,ip_now,ewald,ewd_scr,cell,clatoms_info,
                          clatoms_pos,atommaps,stat_avg,ptens,simopts,for_scr);
     break;

  }/* end switch */

/*==========================================================================*/
   }/*end Routine*/
/*=======================================================================*/





/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void cp_ks_energy_hybrid(CP *cp,int ip_now,EWALD *ewald,EWD_SCR *ewd_scr,
                         CELL *cell,CLATOMS_INFO *clatoms_info,
                         CLATOMS_POS *clatoms_pos,ATOMMAPS *atommaps,
                         STAT_AVG *stat_avg,PTENS *ptens,SIMOPTS *simopts,
                         FOR_SCR *for_scr)

/*==========================================================================*/
/*         Begin Routine                                                    */
   {/*Begin Routine*/
/*=======================================================================*/
/*         Local Variable declarations                                   */

#include "../typ_defs/typ_mask.h"
  
  int i,iii,is,icoef;
  int cp_debug_xc=0;
  int cp_dual_grid_opt_on  = cp->cpopts.cp_dual_grid_opt;
  static int itime=0;
  double cpu1,cpu2;
/*        Local pointers */

  int cp_wave_min;
  int cp_min;
  int cp_min_on;

  int cp_wave_min_pimd= simopts->cp_wave_min_pimd;
  int cp_lsda         = cp->cpopts.cp_lsda;
  int cp_norb         = cp->cpopts.cp_norb;
  int cp_gga          = cp->cpopts.cp_gga;
  int cp_ke_dens_on   = cp->cpcoeffs_info.cp_ke_dens_on;
  int cp_elf_calc_frq = cp->cpcoeffs_info.cp_elf_calc_frq;

  int n_interp_pme_dual = cp->pseudo.n_interp_pme_dual;

  int np_states      = cp->communicate.np_states;
  int myid_state     = cp->communicate.myid_state;
  MPI_Comm  comm_states = cp->communicate.comm_states;
  char *ggax_typ     = cp->pseudo.ggax_typ;
  char *ggac_typ     = cp->pseudo.ggac_typ;

  int ncoef          = cp->cpcoeffs_info.ncoef;
  int ncoef_l        = cp->cpcoeffs_info.ncoef_l;
  int ncoef_l_proc   = cp->cp_para_fft_pkg3d_lg.ncoef_proc;
  int ncoef_l_dens_cp_box = cp->cp_para_fft_pkg3d_dens_cp_box.ncoef_proc;
  int nstate_up      = cp->cpcoeffs_info.nstate_up_proc;
  int *ioff_upt      = cp->cpcoeffs_info.ioff_upt;
  double *occ_up     = cp->cpopts.occ_up;
  int nstate_dn      = cp->cpcoeffs_info.nstate_dn_proc;
  int *ioff_dnt      = cp->cpcoeffs_info.ioff_dnt;
  double *occ_dn     = cp->cpopts.occ_dn;
  double *max_diag     = &(cp->cpcoeffs_pos[ip_now].max_diag);
  double *max_off_diag = &(cp->cpcoeffs_pos[ip_now].max_off_diag);

  int *icoef_orth_up    = &(cp->cpcoeffs_pos[ip_now].icoef_orth_up);
  int *icoef_form_up    = &(cp->cpcoeffs_pos[ip_now].icoef_form_up);
  int *ifcoef_orth_up   = &(cp->cpcoeffs_pos[ip_now].ifcoef_orth_up);
  int *ifcoef_form_up   = &(cp->cpcoeffs_pos[ip_now].ifcoef_form_up);
  double *cre_up        = cp->cpcoeffs_pos[ip_now].cre_up;
  double *cim_up        = cp->cpcoeffs_pos[ip_now].cim_up;
  double *fcre_up       = cp->cpcoeffs_pos[ip_now].fcre_up;
  double *fcim_up       = cp->cpcoeffs_pos[ip_now].fcim_up;
  double *ksmat_up      = cp->cpcoeffs_pos[ip_now].ksmat_up;
  double *norbmat_up    = cp->cpcoeffs_pos[ip_now].norbmat_up;
  double *norbmati_up   = cp->cpcoeffs_pos[ip_now].norbmati_up;
  double *ovmat_eigv_up = cp->cpcoeffs_pos[ip_now].ovmat_eigv_up;
  double *cre_dn        = cp->cpcoeffs_pos[ip_now].cre_dn;
  double *cim_dn        = cp->cpcoeffs_pos[ip_now].cim_dn;
  double *fcre_dn       = cp->cpcoeffs_pos[ip_now].fcre_dn;
  double *fcim_dn       = cp->cpcoeffs_pos[ip_now].fcim_dn;
  double *ksmat_dn      = cp->cpcoeffs_pos[ip_now].ksmat_dn;
  int *icoef_orth_dn    = &(cp->cpcoeffs_pos[ip_now].icoef_orth_dn);
  int *icoef_form_dn    = &(cp->cpcoeffs_pos[ip_now].icoef_form_dn);
  int *ifcoef_orth_dn   = &(cp->cpcoeffs_pos[ip_now].ifcoef_orth_dn);
  int *ifcoef_form_dn   = &(cp->cpcoeffs_pos[ip_now].ifcoef_form_dn);
  double *norbmat_dn    = cp->cpcoeffs_pos[ip_now].norbmat_dn;
  double *norbmati_dn   = cp->cpcoeffs_pos[ip_now].norbmati_dn;
  double *ovmat_eigv_dn = cp->cpcoeffs_pos[ip_now].ovmat_eigv_dn;

  double *cpscr_cre_up    = cp->cpscr.cpscr_wave.cre_up;
  double *cpscr_cim_up    = cp->cpscr.cpscr_wave.cim_up;
  double *rho_up          = cp->cpscr.cpscr_rho.rho_up;
  double *rhocr_up        = cp->cpscr.cpscr_rho.rhocr_up;
  double *rhoci_up        = cp->cpscr.cpscr_rho.rhoci_up;
  double *rhocr_up_dens_cp_box        = cp->cpscr.cpscr_rho.rhocr_up_dens_cp_box;
  double *rhoci_up_dens_cp_box        = cp->cpscr.cpscr_rho.rhoci_up_dens_cp_box;
  double *rhocr_dn_dens_cp_box        = cp->cpscr.cpscr_rho.rhocr_dn_dens_cp_box;
  double *rhoci_dn_dens_cp_box        = cp->cpscr.cpscr_rho.rhoci_dn_dens_cp_box;
  double *d_rhox_up       = cp->cpscr.cpscr_grho.d_rhox_up;
  double *d_rhoy_up       = cp->cpscr.cpscr_grho.d_rhoy_up;
  double *d_rhoz_up       = cp->cpscr.cpscr_grho.d_rhoz_up;
  double *d2_rho_up       = cp->cpscr.cpscr_grho.d2_rho_up;
  double *elec_ke_dens_up = cp->cpscr.cpscr_grho.elec_ke_dens_up;
  double *cpscr_cre_dn    = cp->cpscr.cpscr_wave.cre_dn;
  double *cpscr_cim_dn    = cp->cpscr.cpscr_wave.cim_dn;
  double *rho_dn          = cp->cpscr.cpscr_rho.rho_dn;
  double *rhoci_dn        = cp->cpscr.cpscr_rho.rhoci_dn;
  double *rhocr_dn        = cp->cpscr.cpscr_rho.rhocr_dn;
  double *d_rhox_dn       = cp->cpscr.cpscr_grho.d_rhox_dn;
  double *d_rhoy_dn       = cp->cpscr.cpscr_grho.d_rhoy_dn;
  double *d_rhoz_dn       = cp->cpscr.cpscr_grho.d_rhoz_dn;
  double *d2_rho_dn       = cp->cpscr.cpscr_grho.d2_rho_dn;
  double *elec_ke_dens_dn = cp->cpscr.cpscr_grho.elec_ke_dens_dn;
  double *ksmat_scr       = cp->cpscr.cpscr_ovmat.ovlap1;

  double *cp_elf_up       = cp->electronic_properties.cp_elf_up;
  double *cp_elf_dn       = cp->electronic_properties.cp_elf_dn;

  double *ptens_pvten_tot   = ptens->pvten_tot;
  double *ptens_pvten       = ptens->pvten;
  double *ptens_pvten_tmp   = ptens->pvten_tmp;
  double *x                 = clatoms_pos->x;
  double *y                 = clatoms_pos->y;
  double *z                 = clatoms_pos->z;
 
  double integral,int_tmp;
  int   nfft_proc        =    cp->cp_para_fft_pkg3d_lg.nfft_proc;
  int nfft2_proc = nfft_proc/2;

  int *iatm_atm_typ    = atommaps->iatm_atm_typ;
  int natm_tot         = clatoms_info->natm_tot;
  int iperd            = cell->iperd;

  double tol_edge_dist = cp->cpopts.tol_edge_dist;
  int icheck_perd_size = cp->cpopts.icheck_perd_size;
  int icheck_dual_size = cp->cpopts.icheck_dual_size;

  cp_wave_min = simopts->cp_wave_min;
  cp_min      = simopts->cp_min;
  cp_min_on = cp_wave_min + cp_min + cp_wave_min_pimd;

/*======================================================================*/
/* 0) Check the forms                                                   */

   if(cp_norb>0){
    if((*icoef_orth_up)!=0){
      printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
      printf("Up Coefs must be in nonorthonormal form under norb \n");
      printf("on state processor %d in cp_elec_energy_ctrl \n",myid_state);
      printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
      fflush(stdout);
      exit(1);
    }/*endif*/
    if(cp_lsda==1){
     if((*icoef_orth_dn)!=0){
      printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
      printf("Dn Coefs must be in nonorthonormal form under norb \n");
      printf("on state processor %d in cp_elec_energy_ctrl \n",myid_state);
      printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
      fflush(stdout);
      exit(1);
     }/*endif*/
    }/*endif*/
   }/*endif*/

/*======================================================================*/
/* 0.05) Check the approximations in the methods                        */


  if( ((iperd<3) || (iperd==4)) && (icheck_perd_size==1) ){
    cp_boundary_check(cell,clatoms_info,clatoms_pos,tol_edge_dist);
  }/*endif*/
  if( (cp_dual_grid_opt_on>=1) && (icheck_dual_size==1) ){
    cp_dual_check(cell,clatoms_info,clatoms_pos,
                  atommaps->cp_atm_lst,tol_edge_dist);
  }/*endif*/

/*======================================================================*/
/* 0.1) Set the GGA flags                                               */

  cp->cpopts.cp_becke=0;
  cp->cpopts.cp_pw91x=0;
  cp->cpopts.cp_fila_1x=0;
  cp->cpopts.cp_fila_2x=0;
  cp->cpopts.cp_pbe_x=0;
  cp->cpopts.cp_revpbe_x=0;
  cp->cpopts.cp_rpbe_x=0;
  cp->cpopts.cp_xpbe_x=0;
  cp->cpopts.cp_brx89=0;
  cp->cpopts.cp_brx2k=0;
  cp->cpopts.cp_lyp=0;  
  cp->cpopts.cp_lypm1=0;  
  cp->cpopts.cp_pw91c=0;
  cp->cpopts.cp_pbe_c=0;
  cp->cpopts.cp_xpbe_c=0;
  cp->cpopts.cp_tau1_c=0;
  cp->cpopts.cp_debug_xc=0;
  if(cp_gga == 1){
    if(strcasecmp(ggax_typ,"becke"   )==0){cp->cpopts.cp_becke=1;}
    if(strcasecmp(ggax_typ,"pw91x"   )==0){cp->cpopts.cp_pw91x=1;}
    if(strcasecmp(ggax_typ,"fila_1x" )==0){cp->cpopts.cp_fila_1x=1;}
    if(strcasecmp(ggax_typ,"fila_2x" )==0){cp->cpopts.cp_fila_2x=1;}
    if(strcasecmp(ggax_typ,"pbe_x"   )==0){cp->cpopts.cp_pbe_x=1;}
    if(strcasecmp(ggax_typ,"revpbe_x")==0){cp->cpopts.cp_revpbe_x=1;}
    if(strcasecmp(ggax_typ,"rpbe_x"  )==0){cp->cpopts.cp_rpbe_x=1;}
    if(strcasecmp(ggax_typ,"xpbe_x"  )==0){cp->cpopts.cp_xpbe_x=1;}
    if(strcasecmp(ggax_typ,"brx89"   )==0){cp->cpopts.cp_brx89=1;}
    if(strcasecmp(ggax_typ,"brx2k"   )==0){cp->cpopts.cp_brx2k=1;}
    if(strcasecmp(ggac_typ,"lyp"     )==0){cp->cpopts.cp_lyp=1;  }
    if(strcasecmp(ggac_typ,"lypm1"   )==0){cp->cpopts.cp_lypm1=1;  }
    if(strcasecmp(ggac_typ,"pw91c"   )==0){cp->cpopts.cp_pw91c=1;}
    if(strcasecmp(ggac_typ,"pbe_c"   )==0){cp->cpopts.cp_pbe_c=1;}
    if(strcasecmp(ggac_typ,"xpbe_c"  )==0){cp->cpopts.cp_xpbe_c=1;}
    if(strcasecmp(ggac_typ,"tau1_c"  )==0){cp->cpopts.cp_tau1_c=1;}
    if(strcasecmp(ggac_typ,"debug97x")==0){cp->cpopts.cp_debug_xc=1;}
  }/*endif*/

/*======================================================================*/
/* I) Orthogonalize the coefs if norbing                                */


  if(cp_norb>0){
    (*max_diag)     = 0.0;
    (*max_off_diag) = 0.0;
#ifdef TIME_CP
    if(np_states>1){Barrier(comm_states);}
    cputime(&cpu1);
#endif
    cp_rotate_coef_ortho(cre_up,cim_up,*icoef_form_up,icoef_orth_up,
                         norbmat_up,norbmati_up,ovmat_eigv_up,
                         cpscr_cre_up,cpscr_cim_up,
                         occ_up,ioff_upt,max_off_diag,max_diag,
                         &(cp->cpscr.cpscr_ovmat),&(cp->cp_comm_state_pkg_up));
    if((cp_lsda==1) && (nstate_dn > 0) ){
     cp_rotate_coef_ortho(cre_dn,cim_dn,*icoef_form_dn,icoef_orth_dn,
                         norbmat_dn,norbmati_dn,ovmat_eigv_dn,
                         cpscr_cre_dn,cpscr_cim_dn,
                         occ_dn,ioff_dnt,max_off_diag,max_diag,
                         &(cp->cpscr.cpscr_ovmat),&(cp->cp_comm_state_pkg_dn));
    }/*endif*/
#ifdef TIME_CP
    cputime(&cpu2);
    par_cpu_vomit((cpu2-cpu1),comm_states,np_states,myid_state,
                       "cp_rotate_coef_ortho");
#endif
  }/*endif*/

/*======================================================================*/
/* II) In parallel, transpose coefs back to normal form                 */

  if(np_states>1){

#ifdef TIME_CP
   if(np_states>1){Barrier(comm_states);}
   cputime(&cpu1);
#endif
   cp_transpose_bck(cre_up,cim_up,icoef_form_up,
                    cpscr_cre_up,cpscr_cim_up,&(cp->cp_comm_state_pkg_up));
   if((cp_lsda==1) && (nstate_dn > 0) ){
    cp_transpose_bck(cre_dn,cim_dn,icoef_form_dn,
                     cpscr_cre_dn,cpscr_cim_dn,&(cp->cp_comm_state_pkg_dn));
   }/*endif*/

#ifdef TIME_CP
   cputime(&cpu2);
   par_cpu_vomit((cpu2-cpu1),comm_states,np_states,myid_state,
                       "cp_transpose_bck");
#endif
  }/*endif*/

/*======================================================================*/
/* III) Initialize forces, pressure tensor, inverse hmat                */

  (*ifcoef_form_up) = 0;
  (*ifcoef_orth_up) = 1;
  for(i=1;i<=ncoef*nstate_up;i++){
    fcre_up[i] = 0.0;
    fcim_up[i] = 0.0;
  }/*endfor*/
  if( (cp_lsda == 1) && (nstate_dn != 0) ){
    (*ifcoef_form_dn) = 0;
    (*ifcoef_orth_dn) = 1;
    for(i=1;i<=ncoef*nstate_dn;i++){
      fcre_dn[i] = 0.0;
      fcim_dn[i] = 0.0;
    }/*endfor*/
  }/*endif*/

  for(i=1;i<=9;i++){ptens_pvten_tmp[i] = 0.0;}
  gethinv(cell->hmat_cp,cell->hmati_cp,&(cell->vol_cp),iperd);
  gethinv(cell->hmat,cell->hmati,&(cell->vol),iperd);

/*======================================================================*/
/* IV) Get the total density, (spin densities too if lsda)              */
/*       and necessary gradients of density for GGA calculations        */

#ifdef TIME_CP
  if(np_states>1){Barrier(comm_states);}
  cputime(&cpu1);
#endif

  cp_rho_calc_hybrid(&(cp->cpewald),&(cp->cpscr),&(cp->cpcoeffs_info),
                     ewald,cell,cre_up,cim_up,*icoef_form_up,*icoef_orth_up,
                     rhocr_up,rhoci_up,rho_up,rhocr_up_dens_cp_box,rhoci_up_dens_cp_box,
                     d_rhox_up,d_rhoy_up,
                     d_rhoz_up,d2_rho_up,nstate_up,ncoef,
                     cp_gga,cp_dual_grid_opt_on,n_interp_pme_dual,
                     &(cp->communicate),
                     &(cp->cp_para_fft_pkg3d_lg),&(cp->cp_sclr_fft_pkg3d_lg),
                     &(cp->cp_para_fft_pkg3d_dens_cp_box),
                     &(cp->cp_sclr_fft_pkg3d_dens_cp_box),
                     &(cp->cp_sclr_fft_pkg3d_sm));

  if((cp_lsda== 1) && (nstate_dn!= 0) ){
  cp_rho_calc_hybrid(&(cp->cpewald),&(cp->cpscr),&(cp->cpcoeffs_info),
                     ewald,cell,cre_dn,cim_dn,*icoef_form_dn,*icoef_orth_dn,
                     rhocr_dn,rhoci_dn,rho_dn,rhocr_dn_dens_cp_box,rhoci_dn_dens_cp_box,
                     d_rhox_dn,d_rhoy_dn,
                     d_rhoz_dn,d2_rho_dn,nstate_dn,ncoef,
                     cp_gga,cp_dual_grid_opt_on,n_interp_pme_dual,
                     &(cp->communicate),&(cp->cp_para_fft_pkg3d_lg),
                     &(cp->cp_sclr_fft_pkg3d_lg),
                     &(cp->cp_para_fft_pkg3d_dens_cp_box),
                     &(cp->cp_sclr_fft_pkg3d_dens_cp_box),
                     &(cp->cp_sclr_fft_pkg3d_sm));
    for(i=1;i <= ncoef_l_proc;i++) {
      rhocr_up[i] += rhocr_dn[i];
      rhoci_up[i] += rhoci_dn[i];
    }/* endfor */
    if(cp_dual_grid_opt_on >= 1){
      for(i=1;i<= ncoef_l_dens_cp_box; i++){
        rhocr_up_dens_cp_box[i] += rhocr_dn_dens_cp_box[i];
        rhoci_up_dens_cp_box[i] += rhoci_dn_dens_cp_box[i];
      }/* endfor */
    } /* endif */
  }/* endif */
#ifdef TIME_CP
  cputime(&cpu2);
  par_cpu_vomit((cpu2-cpu1),comm_states,np_states,myid_state,
                      "cp_rho_calc");
#endif

/*======================================================================*/
/* V) If required get the kinetic energy density (similar to density)   */

  if(cp_ke_dens_on==1){
    cp_ke_dens_calc_hybrid(&(cp->cpewald),&(cp->cpscr),
                            cell,cre_up,cim_up,
                            *icoef_form_up,*icoef_orth_up,
                            elec_ke_dens_up,
                            nstate_up,ncoef,
                            cp_dual_grid_opt_on,
                            &(cp->communicate),
                            &(cp->cp_para_fft_pkg3d_lg),
                            &(cp->cp_sclr_fft_pkg3d_lg),
                            &(cp->cp_para_fft_pkg3d_dens_cp_box),
                            &(cp->cp_sclr_fft_pkg3d_dens_cp_box),
                            &(cp->cp_sclr_fft_pkg3d_sm));
    if((cp_lsda == 1) && (nstate_dn != 0)){
      cp_ke_dens_calc_hybrid(&(cp->cpewald),&(cp->cpscr),
                              cell,cre_dn,cim_dn,
                              *icoef_form_dn,*icoef_orth_dn,
                              elec_ke_dens_dn,
                              nstate_dn,ncoef,
                              cp_dual_grid_opt_on,
                              &(cp->communicate),
                              &(cp->cp_para_fft_pkg3d_lg),
                              &(cp->cp_sclr_fft_pkg3d_lg),
                              &(cp->cp_para_fft_pkg3d_dens_cp_box),
                              &(cp->cp_sclr_fft_pkg3d_dens_cp_box),
                              &(cp->cp_sclr_fft_pkg3d_sm));
    } /* endif lsda */

  }/* endif cp_ke_dens_on */

/*======================================================================*/
/* VI) If necessary construct the electron localization function         */

  if(cp_elf_calc_frq > 0 && ip_now==1){
      construct_elf(cp_elf_up,elec_ke_dens_up,
                    d_rhox_up,d_rhoy_up,d_rhoz_up,rho_up,
                    cp->pseudo.gga_cut,
                    &(cp->cp_para_fft_pkg3d_lg)); 
      if((cp_lsda == 1) && (nstate_dn != 0)){
        construct_elf(cp_elf_dn,elec_ke_dens_dn,
                      d_rhox_dn,d_rhoy_dn,d_rhoz_dn,rho_dn,
                      cp->pseudo.gga_cut,
                      &(cp->cp_para_fft_pkg3d_lg)); 
      }/* endif lsda */
  } /* endif construct elf */

/*======================================================================*/
/*   VII) Calculate the non-local pseudopotential list                    */
/*      Make sure that particles are not too far apart for clusters     */
/*      Calculate the external potential in g space.                    */
/*      Calculate the non-local pseudopotential energy                  */
/*      Calculate the particle forces                                   */
/*      Calculate the coef forces from the non-local pseudopotential    */


  if(itime == 0 || cp_dual_grid_opt_on >= 1){
     control_vps_atm_list(&(cp->pseudo),cell,clatoms_pos,clatoms_info,
                          atommaps,ewd_scr,for_scr,cp_dual_grid_opt_on,itime);
    itime=1;
  }/*endif*/

#ifdef TIME_CP
  if(np_states>1){Barrier(comm_states);}
  cputime(&cpu1);
#endif
  control_cp_eext_recip(clatoms_info,clatoms_pos,&(cp->cpcoeffs_info),
                       &(cp->cpcoeffs_pos[ip_now]),
                       &(cp->cpewald),&(cp->cpscr),
                       &(cp->cpopts),&(cp->pseudo),ewd_scr,atommaps,cell,
                       ewald,ptens,&(stat_avg->vrecip),
                       &(stat_avg->cp_enl),&(cp->communicate),for_scr,
                       cp_dual_grid_opt_on,
                       &(cp->cp_para_fft_pkg3d_lg));

#ifdef TIME_CP
  cputime(&cpu2);
  par_cpu_vomit((cpu2-cpu1),comm_states,np_states,myid_state,
                      "control_cp_eext_recip");
#endif

/*======================================================================*/
/* VIII) Calculate the Hartree and exchange correlation energy             */
/*     Calculate the local external potential energy                     */
/*     Calculate the electronic kinetic energy                           */
/*     Calculate the rest of the coef forces                             */


#ifdef TIME_CP
  if(np_states>1){Barrier(comm_states);}
  cputime(&cpu1);
#endif
  coef_force_control(&(cp->cpopts),&(cp->cpcoeffs_info),
                            &(cp->cpcoeffs_pos[ip_now]),
                            &(cp->cpscr),ewald,&(cp->cpewald),cell,stat_avg,
                            cp->pseudo.vxc_typ,ptens->pvten_tmp,
                            cp->pseudo.gga_cut,cp->pseudo.alpha_conv_dual,
                            cp->pseudo.n_interp_pme_dual,cp_min_on,
                            &(cp->communicate),
                            &(cp->cp_comm_state_pkg_up),
                            &(cp->cp_comm_state_pkg_dn),
                            &(cp->cp_para_fft_pkg3d_lg),
                            &(cp->cp_sclr_fft_pkg3d_lg),
                            &(cp->cp_para_fft_pkg3d_dens_cp_box),
                            &(cp->cp_sclr_fft_pkg3d_dens_cp_box),
                            &(cp->cp_para_fft_pkg3d_sm),
                            &(cp->cp_sclr_fft_pkg3d_sm),
                            cp_dual_grid_opt_on);
#ifdef TIME_CP
  cputime(&cpu2);
  par_cpu_vomit((cpu2-cpu1),comm_states,np_states,myid_state,
                      "coef_force_control");
#endif
/*======================================================================*/
/*  IX) In parallel, transpose coefs and coef forces fwd               */

  if(np_states>1){
#ifdef TIME_CP
    if(np_states>1) {Barrier(comm_states);}
    cputime(&cpu1);
#endif
    cp_transpose_fwd(cre_up,cim_up,icoef_form_up,
                    cpscr_cre_up,cpscr_cim_up,&(cp->cp_comm_state_pkg_up));
    cp_transpose_fwd(fcre_up,fcim_up,ifcoef_form_up,
                    cpscr_cre_up,cpscr_cim_up,&(cp->cp_comm_state_pkg_up));
    if( (cp_lsda==1) && (nstate_dn > 0) ){
     cp_transpose_fwd(cre_dn,cim_dn,icoef_form_dn,
                     cpscr_cre_dn,cpscr_cim_dn,&(cp->cp_comm_state_pkg_dn));
     cp_transpose_fwd(fcre_dn,fcim_dn,ifcoef_form_dn,
                     cpscr_cre_dn,cpscr_cim_dn,&(cp->cp_comm_state_pkg_dn));
    }/*endif*/
#ifdef TIME_CP
    cputime(&cpu2);
    par_cpu_vomit((cpu2-cpu1),comm_states,np_states,myid_state,
                      "cp_transpose_fwd");
#endif
  }/*endif*/

/*=========================================================================*/
/* X) KS matrix used as an approximation  to the Lagrange multipliers     */
/*      when minimizing. It is a necessary part of the forces when norbing */

  if( (cp_min_on+cp_norb) >=1 ){

#ifdef TIME_CP
   if(np_states>1){Barrier(comm_states);}
   cputime(&cpu1);
#endif
     
     cp_add_ksmat_force(cre_up,cim_up,*icoef_form_up,*icoef_orth_up,
                        fcre_up,fcim_up,*ifcoef_form_up,*ifcoef_orth_up,
                        ksmat_up,ksmat_scr,ioff_upt,cp_lsda,cp_min_on,occ_up,
                        &(cp->cp_comm_state_pkg_up));
     if( (cp_lsda==1) && (nstate_dn!=0) ){
       cp_add_ksmat_force(cre_dn,cim_dn,*icoef_form_dn,*icoef_orth_dn,
                          fcre_dn,fcim_dn,*ifcoef_form_dn,*ifcoef_orth_dn,
                          ksmat_dn,ksmat_scr,ioff_dnt,cp_lsda,cp_min_on,occ_dn,
                          &(cp->cp_comm_state_pkg_dn));
     }/*endif*/

#ifdef TIME_CP
  cputime(&cpu2);
  par_cpu_vomit((cpu2-cpu1),comm_states,np_states,myid_state,
                     "cp_add_ksmat_force");
#endif
  }/*endif*/
		      
/*======================================================================*/
/* XI) DeOrthogonalize the coefs and coef forces: This requires          */
/*                  a matrix and its inverse (chain rule)               */

  if(cp_norb>0){
#ifdef TIME_CP
   if(np_states>1){Barrier(comm_states);}
   cputime(&cpu1);
#endif
    cp_rotate_gen_nonortho(cre_up,cim_up,*icoef_form_up,icoef_orth_up,
                           norbmati_up,ioff_upt,cpscr_cre_up,cpscr_cim_up,
                           &(cp->cp_comm_state_pkg_up));
    cp_rotate_gen_nonortho(fcre_up,fcim_up,*ifcoef_form_up,ifcoef_orth_up,
                           norbmat_up,ioff_upt,cpscr_cre_up,cpscr_cim_up,
                           &(cp->cp_comm_state_pkg_up));
    if((cp_lsda==1) && (nstate_dn > 0) ){
     cp_rotate_gen_nonortho(cre_dn,cim_dn,*icoef_form_dn,icoef_orth_dn,
                            norbmati_dn,ioff_dnt,cpscr_cre_dn,cpscr_cim_dn,
                            &(cp->cp_comm_state_pkg_dn));
     cp_rotate_gen_nonortho(fcre_dn,fcim_dn,*ifcoef_form_dn,ifcoef_orth_dn,
                            norbmat_dn,ioff_dnt,cpscr_cre_dn,cpscr_cim_dn,
                            &(cp->cp_comm_state_pkg_dn));
    }/*endif:lsda*/
#ifdef TIME_CP
  cputime(&cpu2);
  par_cpu_vomit((cpu2-cpu1),comm_states,np_states,myid_state,
                     "cp_rotate_get_nonortho");
#endif
  }/*endif:cp_norb*/

/*======================================================================*/
/* XII) Increment the pressure tensor                                    */

  for(i=1;i<=9;i++){
    ptens_pvten_tot[i] += ptens_pvten_tmp[i];
    ptens_pvten[i]     += ptens_pvten_tmp[i];
  }/*endfor*/

/*==========================================================================*/
   }/*end routine*/
/*==========================================================================*/




/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void cp_ks_energy_full_g(CP *cp,int ip_now,EWALD *ewald,EWD_SCR *ewd_scr,
                         CELL *cell,CLATOMS_INFO *clatoms_info,
                         CLATOMS_POS *clatoms_pos,ATOMMAPS *atommaps,
                         STAT_AVG *stat_avg,PTENS *ptens,SIMOPTS *simopts,
                         FOR_SCR *for_scr)

/*==========================================================================*/
/*         Begin Routine                                                    */
   {/*Begin Routine*/
/*=======================================================================*/
/*         Local Variable declarations                                   */

#include "../typ_defs/typ_mask.h"
  
  int i,iii,is,icoef;
  int cp_debug_xc=0;
  int cp_min_on;
  int cp_dual_grid_opt_on = cp->cpopts.cp_dual_grid_opt;
  static int itime=0;
  double cpu1,cpu2;
/*        Local pointers */

  int cp_wave_min;
  int cp_min;

  int cp_wave_min_pimd = simopts->cp_wave_min_pimd;
  int cp_lsda          = cp->cpopts.cp_lsda;
  int cp_norb          = cp->cpopts.cp_norb;
  int cp_gga           = cp->cpopts.cp_gga;
  int cp_ke_dens_on    = cp->cpcoeffs_info.cp_ke_dens_on;
  int cp_elf_calc_frq  = cp->cpcoeffs_info.cp_elf_calc_frq;
 
  int n_interp_pme_dual = cp->pseudo.n_interp_pme_dual;
  int np_states      = cp->communicate.np_states;
  int myid_state     = cp->communicate.myid_state;
  MPI_Comm  comm_states = cp->communicate.comm_states;
  char *ggax_typ     = cp->pseudo.ggax_typ;
  char *ggac_typ     = cp->pseudo.ggac_typ;

  int ncoef          = cp->cpcoeffs_info.ncoef;
  int ncoef_l        = cp->cpcoeffs_info.ncoef_l;
  int ncoef_l_proc   = cp->cp_para_fft_pkg3d_lg.ncoef_proc;
  int ncoef_l_dens_cp_box = cp->cp_para_fft_pkg3d_dens_cp_box.ncoef_proc;
  int nstate_ncoef_proc_max_up = cp->cpcoeffs_info.nstate_ncoef_proc_max_up;
  int nstate_ncoef_proc_max_dn = cp->cpcoeffs_info.nstate_ncoef_proc_max_dn;
  int nstate_up_proc = cp->cpcoeffs_info.nstate_up_proc;
  int nstate_up      = cp->cpcoeffs_info.nstate_up;
  int *ioff_upt      = cp->cpcoeffs_info.ioff_upt;
  double *occ_up     = cp->cpopts.occ_up;
  int nstate_dn_proc = cp->cpcoeffs_info.nstate_dn_proc;
  int nstate_dn      = cp->cpcoeffs_info.nstate_dn;
  int *ioff_dnt      = cp->cpcoeffs_info.ioff_dnt;
  double *occ_dn     = cp->cpopts.occ_dn;
  double *max_diag     = &(cp->cpcoeffs_pos[ip_now].max_diag);
  double *max_off_diag = &(cp->cpcoeffs_pos[ip_now].max_off_diag);

  int *icoef_orth_up    = &(cp->cpcoeffs_pos[ip_now].icoef_orth_up);
  int *icoef_form_up    = &(cp->cpcoeffs_pos[ip_now].icoef_form_up);
  int *ifcoef_orth_up   = &(cp->cpcoeffs_pos[ip_now].ifcoef_orth_up);
  int *ifcoef_form_up   = &(cp->cpcoeffs_pos[ip_now].ifcoef_form_up);
  double *cre_up        = cp->cpcoeffs_pos[ip_now].cre_up;
  double *cim_up        = cp->cpcoeffs_pos[ip_now].cim_up;
  double *fcre_up       = cp->cpcoeffs_pos[ip_now].fcre_up;
  double *fcim_up       = cp->cpcoeffs_pos[ip_now].fcim_up;
  double *ksmat_up      = cp->cpcoeffs_pos[ip_now].ksmat_up;
  double *norbmat_up    = cp->cpcoeffs_pos[ip_now].norbmat_up;
  double *norbmati_up   = cp->cpcoeffs_pos[ip_now].norbmati_up;
  double *ovmat_eigv_up = cp->cpcoeffs_pos[ip_now].ovmat_eigv_up;
  double *cre_dn        = cp->cpcoeffs_pos[ip_now].cre_dn;
  double *cim_dn        = cp->cpcoeffs_pos[ip_now].cim_dn;
  double *fcre_dn       = cp->cpcoeffs_pos[ip_now].fcre_dn;
  double *fcim_dn       = cp->cpcoeffs_pos[ip_now].fcim_dn;
  double *ksmat_dn      = cp->cpcoeffs_pos[ip_now].ksmat_dn;
  int *icoef_orth_dn    = &(cp->cpcoeffs_pos[ip_now].icoef_orth_dn);
  int *icoef_form_dn    = &(cp->cpcoeffs_pos[ip_now].icoef_form_dn);
  int *ifcoef_orth_dn   = &(cp->cpcoeffs_pos[ip_now].ifcoef_orth_dn);
  int *ifcoef_form_dn   = &(cp->cpcoeffs_pos[ip_now].ifcoef_form_dn);
  double *norbmat_dn    = cp->cpcoeffs_pos[ip_now].norbmat_dn;
  double *norbmati_dn   = cp->cpcoeffs_pos[ip_now].norbmati_dn;
  double *ovmat_eigv_dn = cp->cpcoeffs_pos[ip_now].ovmat_eigv_dn;

  double *cpscr_cre_up    = cp->cpscr.cpscr_wave.cre_up;
  double *cpscr_cim_up    = cp->cpscr.cpscr_wave.cim_up;
  double *rho_up          = cp->cpscr.cpscr_rho.rho_up;
  double *rhocr_up        = cp->cpscr.cpscr_rho.rhocr_up;
  double *rhoci_up        = cp->cpscr.cpscr_rho.rhoci_up;
  double *rhocr_up_dens_cp_box        = cp->cpscr.cpscr_rho.rhocr_up_dens_cp_box;
  double *rhoci_up_dens_cp_box        = cp->cpscr.cpscr_rho.rhoci_up_dens_cp_box;
  double *rhocr_dn_dens_cp_box        = cp->cpscr.cpscr_rho.rhocr_dn_dens_cp_box;
  double *rhoci_dn_dens_cp_box        = cp->cpscr.cpscr_rho.rhoci_dn_dens_cp_box;
  double *d_rhox_up       = cp->cpscr.cpscr_grho.d_rhox_up;
  double *d_rhoy_up       = cp->cpscr.cpscr_grho.d_rhoy_up;
  double *d_rhoz_up       = cp->cpscr.cpscr_grho.d_rhoz_up;
  double *d2_rho_up       = cp->cpscr.cpscr_grho.d2_rho_up;
  double *elec_ke_dens_up = cp->cpscr.cpscr_grho.elec_ke_dens_up;
  double *cpscr_cre_dn    = cp->cpscr.cpscr_wave.cre_dn;
  double *cpscr_cim_dn    = cp->cpscr.cpscr_wave.cim_dn;
  double *rho_dn          = cp->cpscr.cpscr_rho.rho_dn;
  double *rhoci_dn        = cp->cpscr.cpscr_rho.rhoci_dn;
  double *rhocr_dn        = cp->cpscr.cpscr_rho.rhocr_dn;
  double *d_rhox_dn       = cp->cpscr.cpscr_grho.d_rhox_dn;
  double *d_rhoy_dn       = cp->cpscr.cpscr_grho.d_rhoy_dn;
  double *d_rhoz_dn       = cp->cpscr.cpscr_grho.d_rhoz_dn;
  double *d2_rho_dn       = cp->cpscr.cpscr_grho.d2_rho_dn;
  double *ksmat_scr       = cp->cpscr.cpscr_ovmat.ovlap1;
  double *elec_ke_dens_dn = cp->cpscr.cpscr_grho.elec_ke_dens_dn;

  double *cp_elf_up       = cp->electronic_properties.cp_elf_up;
  double *cp_elf_dn       = cp->electronic_properties.cp_elf_dn;

  double integral,int_tmp;
  int   nfft_proc        =    cp->cp_para_fft_pkg3d_lg.nfft_proc;
  int nfft2_proc = nfft_proc/2;


  double *ptens_pvten_tot   = ptens->pvten_tot;
  double *ptens_pvten       = ptens->pvten;
  double *ptens_pvten_tmp   = ptens->pvten_tmp;

  double *x            = clatoms_pos->x;
  double *y            = clatoms_pos->y;
  double *z            = clatoms_pos->z;

  int *iatm_atm_typ    = atommaps->iatm_atm_typ;
  int natm_tot         = clatoms_info->natm_tot;
  int iperd            = cell->iperd;

  double tol_edge_dist = cp->cpopts.tol_edge_dist; 
  int icheck_perd_size = cp->cpopts.icheck_perd_size;
  int icheck_dual_size = cp->cpopts.icheck_dual_size;

  cp_wave_min    = simopts->cp_wave_min;
  cp_min         = simopts->cp_min;

  cp_min_on      = cp_wave_min + cp_min + cp_wave_min_pimd;

/*======================================================================*/
/* 0) Check the forms                                                   */

   if(cp_norb>0){
    if((*icoef_orth_up)!=0){
      printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
      printf("Up Coefs must be in nonorthonormal form under norb \n");
      printf("on state processor %d in cp_elec_energy_ctrl \n",myid_state);
      printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
      fflush(stdout);
      exit(1);
    }/*endif*/
    if(cp_lsda==1){
     if((*icoef_orth_dn)!=0){
      printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
      printf("Dn Coefs must be in nonorthonormal form under norb \n");
      printf("on state processor %d in cp_elec_energy_ctrl \n",myid_state);
      printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
      fflush(stdout);
      exit(1);
     }/*endif*/
    }/*endif*/
   }/*endif*/

/*======================================================================*/
/* 0.05) Check the approximations in the methods                        */


  if( ((iperd<3) || (iperd==4)) && (icheck_perd_size==1) ){
    cp_boundary_check(cell,clatoms_info,clatoms_pos,tol_edge_dist);
  }/*endif*/


  if( (cp_dual_grid_opt_on>=1) && (icheck_dual_size==1) ){
    cp_dual_check(cell,clatoms_info,clatoms_pos,
                  atommaps->cp_atm_lst,tol_edge_dist);
  }/*endif*/

/*======================================================================*/
/* 0.1) Set the GGA flags                                               */


  cp->cpopts.cp_becke=0;
  cp->cpopts.cp_pw91x=0;
  cp->cpopts.cp_fila_1x=0;
  cp->cpopts.cp_fila_2x=0;
  cp->cpopts.cp_pbe_x=0;
  cp->cpopts.cp_revpbe_x=0;
  cp->cpopts.cp_rpbe_x=0;
  cp->cpopts.cp_xpbe_x=0;
  cp->cpopts.cp_brx89=0;
  cp->cpopts.cp_brx2k=0;
  cp->cpopts.cp_lyp=0;  
  cp->cpopts.cp_lypm1=0;  
  cp->cpopts.cp_pw91c=0;
  cp->cpopts.cp_pbe_c=0;
  cp->cpopts.cp_xpbe_c=0;
  cp->cpopts.cp_tau1_c=0;
  cp->cpopts.cp_debug_xc=0;
  if(cp_gga == 1){
    if(strcasecmp(ggax_typ,"becke"   )==0){cp->cpopts.cp_becke=1;}
    if(strcasecmp(ggax_typ,"pw91x"   )==0){cp->cpopts.cp_pw91x=1;}
    if(strcasecmp(ggax_typ,"fila_1x" )==0){cp->cpopts.cp_fila_1x=1;}
    if(strcasecmp(ggax_typ,"fila_2x" )==0){cp->cpopts.cp_fila_2x=1;}
    if(strcasecmp(ggax_typ,"pbe_x"   )==0){cp->cpopts.cp_pbe_x=1;}
    if(strcasecmp(ggax_typ,"revpbe_x")==0){cp->cpopts.cp_revpbe_x=1;}
    if(strcasecmp(ggax_typ,"rpbe_x"  )==0){cp->cpopts.cp_rpbe_x=1;}
    if(strcasecmp(ggax_typ,"xpbe_x"  )==0){cp->cpopts.cp_xpbe_x=1;}
    if(strcasecmp(ggax_typ,"brx89"   )==0){cp->cpopts.cp_brx89=1;}
    if(strcasecmp(ggax_typ,"brx2k"   )==0){cp->cpopts.cp_brx2k=1;}
    if(strcasecmp(ggac_typ,"lyp"     )==0){cp->cpopts.cp_lyp=1;  }
    if(strcasecmp(ggac_typ,"lypm1"   )==0){cp->cpopts.cp_lypm1=1;  }
    if(strcasecmp(ggac_typ,"pw91c"   )==0){cp->cpopts.cp_pw91c=1;}
    if(strcasecmp(ggac_typ,"pbe_c"   )==0){cp->cpopts.cp_pbe_c=1;}
    if(strcasecmp(ggac_typ,"xpbe_c"  )==0){cp->cpopts.cp_xpbe_c=1;}
    if(strcasecmp(ggac_typ,"tau1_c"  )==0){cp->cpopts.cp_tau1_c=1;}
    if(strcasecmp(ggac_typ,"debug97x")==0){cp->cpopts.cp_debug_xc=1;}
  }/*endif*/

/*======================================================================*/
/* I) Orthogonalize the coefs if norbing                                */


  if(cp_norb>0){
    (*max_diag)     = 0.0;
    (*max_off_diag) = 0.0;
#ifdef TIME_CP
    if(np_states>1){Barrier(comm_states);}
    cputime(&cpu1);
#endif
    cp_rotate_coef_ortho(cre_up,cim_up,*icoef_form_up,icoef_orth_up,
                         norbmat_up,norbmati_up,ovmat_eigv_up,
                         cpscr_cre_up,cpscr_cim_up,
                         occ_up,ioff_upt,max_off_diag,max_diag,
                         &(cp->cpscr.cpscr_ovmat),&(cp->cp_comm_state_pkg_up));
    if((cp_lsda==1) && (nstate_dn_proc > 0) ){
     cp_rotate_coef_ortho(cre_dn,cim_dn,*icoef_form_dn,icoef_orth_dn,
                         norbmat_dn,norbmati_dn,ovmat_eigv_dn,
                         cpscr_cre_dn,cpscr_cim_dn,
                         occ_dn,ioff_dnt,max_off_diag,max_diag,
                         &(cp->cpscr.cpscr_ovmat),&(cp->cp_comm_state_pkg_dn));
    }/*endif*/
#ifdef TIME_CP
    cputime(&cpu2);
    par_cpu_vomit((cpu2-cpu1),comm_states,np_states,myid_state,
                       "rotate_coef_ortho");
#endif
  }/*endif*/

/*======================================================================*/
/* II) Get the total density, (spin densities too if lsda)              */
/*       and necessary gradients of density for GGA calculations        */


#ifdef TIME_CP
  if(np_states>1){Barrier(comm_states);}
 cputime(&cpu1);
#endif
  cp_rho_calc_full_g(&(cp->cpewald),&(cp->cpscr),&(cp->cpcoeffs_info),
                     ewald,cell,cre_up,cim_up,*icoef_form_up,*icoef_orth_up,
                     rhocr_up,rhoci_up,rho_up,rhocr_up_dens_cp_box,rhoci_up_dens_cp_box,
                     d_rhox_up,d_rhoy_up,
                     d_rhoz_up,d2_rho_up,nstate_up,nstate_ncoef_proc_max_up,
                     cp_gga,cp_dual_grid_opt_on,n_interp_pme_dual,
                     &(cp->communicate),&(cp->cp_para_fft_pkg3d_lg),
                     &(cp->cp_para_fft_pkg3d_dens_cp_box),
                     &(cp->cp_para_fft_pkg3d_sm));
  if((cp_lsda== 1) && (nstate_dn != 0) ){
  cp_rho_calc_full_g(&(cp->cpewald),&(cp->cpscr),&(cp->cpcoeffs_info),
                     ewald,cell,cre_dn,cim_dn,*icoef_form_dn,*icoef_orth_dn,
                     rhocr_dn,rhoci_dn,rho_dn,rhocr_dn_dens_cp_box,rhoci_dn_dens_cp_box,
                     d_rhox_dn,d_rhoy_dn,
                     d_rhoz_dn,d2_rho_dn,nstate_dn,nstate_ncoef_proc_max_dn,
                     cp_gga,cp_dual_grid_opt_on,n_interp_pme_dual,
                     &(cp->communicate),
                     &(cp->cp_para_fft_pkg3d_lg),
                     &(cp->cp_para_fft_pkg3d_dens_cp_box),
                     &(cp->cp_para_fft_pkg3d_sm));
    for(i=1;i <= ncoef_l_proc;i++) {
      rhocr_up[i] += rhocr_dn[i];
      rhoci_up[i] += rhoci_dn[i];
    }/* endfor */
    if(cp_dual_grid_opt_on >= 1){
      for(i=1;i<= ncoef_l_dens_cp_box; i++){
        rhocr_up_dens_cp_box[i] += rhocr_dn_dens_cp_box[i];
        rhoci_up_dens_cp_box[i] += rhoci_dn_dens_cp_box[i];
      }/* endfor */
    } /* endif */
  }/* endif */
#ifdef TIME_CP
  cputime(&cpu2);
  par_cpu_vomit((cpu2-cpu1),comm_states,np_states,myid_state,
                      "cp_rho_calc");
#endif

/*======================================================================*/
/* III) If required get the kinetic energy density (similar to density)   */


  if(cp_ke_dens_on==1){
    cp_ke_dens_calc_full_g(&(cp->cpewald),&(cp->cpscr),
                            cell,cre_up,cim_up,
                            *icoef_form_up,*icoef_orth_up,
                            elec_ke_dens_up,
                            nstate_up,nstate_ncoef_proc_max_up,
                            cp_dual_grid_opt_on,
                            &(cp->communicate),
                            &(cp->cp_para_fft_pkg3d_lg),
                            &(cp->cp_para_fft_pkg3d_dens_cp_box),
                            &(cp->cp_para_fft_pkg3d_sm));
    if((cp_lsda == 1) && (nstate_dn != 0)){
      cp_ke_dens_calc_full_g(&(cp->cpewald),&(cp->cpscr),
                              cell,cre_dn,cim_dn,
                              *icoef_form_dn,*icoef_orth_dn,
                              elec_ke_dens_dn,
                              nstate_dn,nstate_ncoef_proc_max_dn,
                              cp_dual_grid_opt_on,
                              &(cp->communicate),
                              &(cp->cp_para_fft_pkg3d_lg),
                              &(cp->cp_para_fft_pkg3d_dens_cp_box),
                              &(cp->cp_para_fft_pkg3d_sm));
    } /* endif lsda */

  }/* endif cp_ke_dens_on */

/*======================================================================*/
/* IV) In parallel, transpose coefs back to normal form                 */


  if(np_states>1){

#ifdef TIME_CP
   if(np_states>1){Barrier(comm_states);}
   cputime(&cpu1);
#endif
   cp_transpose_bck(cre_up,cim_up,icoef_form_up,
                    cpscr_cre_up,cpscr_cim_up,&(cp->cp_comm_state_pkg_up));
   if((cp_lsda==1) && (nstate_dn > 0) ){
    cp_transpose_bck(cre_dn,cim_dn,icoef_form_dn,
                     cpscr_cre_dn,cpscr_cim_dn,&(cp->cp_comm_state_pkg_dn));
   }/*endif*/

#ifdef TIME_CP
   cputime(&cpu2);
   par_cpu_vomit((cpu2-cpu1),comm_states,np_states,myid_state,
                       "cp_transpose_bck");
#endif
  }/*endif*/

/*======================================================================*/
/* V) If necessary construct the electron localiztion function          */


  if(cp_elf_calc_frq > 0 && ip_now==1){
      construct_elf(cp_elf_up,elec_ke_dens_up,
                    d_rhox_up,d_rhoy_up,d_rhoz_up,rho_up,
                    cp->pseudo.gga_cut,
                    &(cp->cp_para_fft_pkg3d_lg)); 
      if((cp_lsda == 1) && (nstate_dn != 0)){
        construct_elf(cp_elf_dn,elec_ke_dens_dn,
                      d_rhox_dn,d_rhoy_dn,d_rhoz_dn,rho_dn,
                      cp->pseudo.gga_cut,
                      &(cp->cp_para_fft_pkg3d_lg)); 
      }/* endif lsda */
  } /* endif construct elf */

/*======================================================================*/
/* V) Initialize forces, pressure tensor, inverse hmat                */


  (*ifcoef_form_up) = 0;
  (*ifcoef_orth_up) = 1;
  for(i=1;i<=ncoef*nstate_up_proc;i++){
    fcre_up[i] = 0.0;
    fcim_up[i] = 0.0;
  }/*endfor*/
  if( (cp_lsda == 1) && (nstate_dn != 0) ){
    (*ifcoef_form_dn) = 0;
    (*ifcoef_orth_dn) = 1;
    for(i=1;i<=ncoef*nstate_dn_proc;i++){
      fcre_dn[i] = 0.0;
      fcim_dn[i] = 0.0;
    }/*endfor*/
  }/*endif*/

  for(i=1;i<=9;i++){ptens_pvten_tmp[i] = 0.0;}
  gethinv(cell->hmat_cp,cell->hmati_cp,&(cell->vol_cp),iperd);

/*======================================================================*/
/*   VI) Calculate the non-local pseudopotential list                    */
/*      Make sure that particles are not too far apart for clusters     */
/*      Calculate the external potential in g space.                    */
/*      Calculate the non-local pseudopotential energy                  */
/*      Calculate the particle forces                                   */
/*      Calculate the coef forces from the non-local pseudopotential    */



  if(itime == 0 || cp_dual_grid_opt_on >= 1){
     control_vps_atm_list(&(cp->pseudo),cell,clatoms_pos,clatoms_info,
                  atommaps,ewd_scr,for_scr,cp_dual_grid_opt_on,itime);
    itime=1;

  }/*endif*/

#ifdef TIME_CP
  if(np_states>1){Barrier(comm_states);}
  cputime(&cpu1);
#endif


  control_cp_eext_recip(clatoms_info,clatoms_pos,&(cp->cpcoeffs_info),
                       &(cp->cpcoeffs_pos[ip_now]),
                       &(cp->cpewald),&(cp->cpscr),
                       &(cp->cpopts),&(cp->pseudo),ewd_scr,atommaps,cell,
                       ewald,ptens,&(stat_avg->vrecip),
                       &(stat_avg->cp_enl),&(cp->communicate),for_scr,
         		cp_dual_grid_opt_on,
                       &(cp->cp_para_fft_pkg3d_lg));

#ifdef TIME_CP
  cputime(&cpu2);
  par_cpu_vomit((cpu2-cpu1),comm_states,np_states,myid_state,
                      "control_cp_eext_recip");
#endif

/*======================================================================*/
/*  VII) In parallel, transpose coefs and coef forces fwd               */

  if(np_states>1){
#ifdef TIME_CP
    if(np_states>1){Barrier(comm_states);}
    cputime(&cpu1);
#endif
    cp_transpose_fwd(cre_up,cim_up,icoef_form_up,
                    cpscr_cre_up,cpscr_cim_up,&(cp->cp_comm_state_pkg_up));
    cp_transpose_fwd(fcre_up,fcim_up,ifcoef_form_up,
                    cpscr_cre_up,cpscr_cim_up,&(cp->cp_comm_state_pkg_up));
    if( (cp_lsda==1) && (nstate_dn > 0) ){
     cp_transpose_fwd(cre_dn,cim_dn,icoef_form_dn,
                     cpscr_cre_dn,cpscr_cim_dn,&(cp->cp_comm_state_pkg_dn));
     cp_transpose_fwd(fcre_dn,fcim_dn,ifcoef_form_dn,
                     cpscr_cre_dn,cpscr_cim_dn,&(cp->cp_comm_state_pkg_dn));
    }/*endif*/
#ifdef TIME_CP
    cputime(&cpu2);
    par_cpu_vomit((cpu2-cpu1),comm_states,np_states,myid_state,
                      "cp_transpose_fwd");
#endif
  }/*endif*/

/*======================================================================*/
/* VIII) Calculate the Hartree and exchange correlation energy             */
/*     Calculate the local external potential energy                     */
/*     Calculate the electronic kinetic energy                           */
/*     Calculate the rest of the coef forces                             */


#ifdef TIME_CP
  if(np_states>1){Barrier(comm_states);}
  cputime(&cpu1);
#endif

  coef_force_control(&(cp->cpopts),&(cp->cpcoeffs_info),
                     &(cp->cpcoeffs_pos[ip_now]),
                     &(cp->cpscr),ewald,&(cp->cpewald),cell,stat_avg,
                     cp->pseudo.vxc_typ,ptens->pvten_tmp,
                     cp->pseudo.gga_cut,cp->pseudo.alpha_conv_dual,
                     cp->pseudo.n_interp_pme_dual,cp_min_on,
                     &(cp->communicate),
                     &(cp->cp_comm_state_pkg_up),
                     &(cp->cp_comm_state_pkg_dn),
                     &(cp->cp_para_fft_pkg3d_lg),
                     &(cp->cp_sclr_fft_pkg3d_lg),
                     &(cp->cp_para_fft_pkg3d_dens_cp_box),
                     &(cp->cp_sclr_fft_pkg3d_dens_cp_box),
                     &(cp->cp_para_fft_pkg3d_sm),
                     &(cp->cp_sclr_fft_pkg3d_sm),
                     cp_dual_grid_opt_on);

#ifdef TIME_CP
  cputime(&cpu2);
  par_cpu_vomit((cpu2-cpu1),comm_states,np_states,myid_state,
                      "coef_force_control");
#endif
/*=========================================================================*/
/* IX) KS matrix used as an approximation  to the Lagrange multipliers     */
/*      when minimizing. It is a necessary part of the forces when norbing */


  if( (cp_min_on+cp_norb) >=1 ){

#ifdef TIME_CP
   if(np_states>1){Barrier(comm_states);}
   cputime(&cpu1);
#endif
     
     cp_add_ksmat_force(cre_up,cim_up,*icoef_form_up,*icoef_orth_up,
                        fcre_up,fcim_up,*ifcoef_form_up,*ifcoef_orth_up,
                        ksmat_up,ksmat_scr,ioff_upt,cp_lsda,cp_min_on,occ_up,
                        &(cp->cp_comm_state_pkg_up));
     if( (cp_lsda==1) && (nstate_dn !=0) ){
       cp_add_ksmat_force(cre_dn,cim_dn,*icoef_form_dn,*icoef_orth_dn,
                          fcre_dn,fcim_dn,*ifcoef_form_dn,*ifcoef_orth_dn,
                          ksmat_dn,ksmat_scr,ioff_dnt,cp_lsda,cp_min_on,occ_dn,
                          &(cp->cp_comm_state_pkg_dn));
     }/*endif*/

#ifdef TIME_CP
  cputime(&cpu2);
  par_cpu_vomit((cpu2-cpu1),comm_states,np_states,myid_state,
                     "cp_add_ksmat_force");
#endif
  }/*endif*/
		      
/*======================================================================*/
/* X) DeOrthogonalize the coefs and coef forces: This requires          */
/*                  a matrix and its inverse (chain rule)               */


  if(cp_norb>0){
   if(np_states>1){Barrier(comm_states);}
#ifdef TIME_CP
   cputime(&cpu1);
#endif
    cp_rotate_gen_nonortho(cre_up,cim_up,*icoef_form_up,icoef_orth_up,
                           norbmati_up,ioff_upt,cpscr_cre_up,cpscr_cim_up,
                           &(cp->cp_comm_state_pkg_up));
    cp_rotate_gen_nonortho(fcre_up,fcim_up,*ifcoef_form_up,ifcoef_orth_up,
                           norbmat_up,ioff_upt,cpscr_cre_up,cpscr_cim_up,
                           &(cp->cp_comm_state_pkg_up));
    if((cp_lsda==1) && (nstate_dn > 0) ){
     cp_rotate_gen_nonortho(cre_dn,cim_dn,*icoef_form_dn,icoef_orth_dn,
                            norbmati_dn,ioff_dnt,cpscr_cre_dn,cpscr_cim_dn,
                            &(cp->cp_comm_state_pkg_dn));
     cp_rotate_gen_nonortho(fcre_dn,fcim_dn,*ifcoef_form_dn,ifcoef_orth_dn,
                            norbmat_dn,ioff_dnt,cpscr_cre_dn,cpscr_cim_dn,
                            &(cp->cp_comm_state_pkg_dn));
    }/*endif:lsda*/
#ifdef TIME_CP
  cputime(&cpu2);
  par_cpu_vomit((cpu2-cpu1),comm_states,np_states,myid_state,
                     "cp_rotate_get_nonortho");
#endif
  }/*endif:cp_norb*/

/*======================================================================*/
/* XI) Increment the pressure tensor                                    */


  for(i=1;i<=9;i++){
    ptens_pvten_tot[i] += ptens_pvten_tmp[i];
    ptens_pvten[i]     += ptens_pvten_tmp[i];
  }/*endfor*/

/*==========================================================================*/
   }/*end routine*/
/*==========================================================================*/





/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
/* Under CP with reduced periodicity, make sure your particles have not     */
/* escaped the reduced boundary                                             */
/*==========================================================================*/

void cp_boundary_check(CELL *cell,CLATOMS_INFO *clatoms_info,
                       CLATOMS_POS *clatoms_pos,
                       double tol_edge_dist)

/*==========================================================================*/
  {/*Begin Routine*/
/*==========================================================================*/
/*            Local variable declarations                                   */

  int i;
  double dx2,dy2,dz2;
  double xmax_now,ymax_now,zmax_now;
  double xcm,ycm,zcm;

  int natm_tot    = clatoms_info->natm_tot;
  double xmax     = cell->hmat[1];  
  double ymax     = cell->hmat[5];  
  double zmax     = cell->hmat[9];  
  double *x       = clatoms_pos->x;
  double *y       = clatoms_pos->y;
  double *z       = clatoms_pos->z;
  int iperd       = cell->iperd;


/*==========================================================================*/
/* I) Get the maximum distances */
 
  xcm = 0.0;
  ycm = 0.0;
  zcm = 0.0;
  for(i=1;i<=natm_tot;i++){
    xcm += x[i];
    ycm += y[i];
    zcm += z[i];
  }/*endfor*/
  xcm = xcm/((double) (natm_tot) );
  ycm = ycm/((double) (natm_tot) );
  zcm = zcm/((double) (natm_tot) );

  xmax_now = 0.0;
  ymax_now = 0.0;
  zmax_now = 0.0;
  for(i=1;i<=natm_tot;i++){
    dx2 = (x[i]-xcm)*(x[i]-xcm);
    dy2 = (y[i]-ycm)*(y[i]-ycm);
    dz2 = (z[i]-zcm)*(z[i]-zcm);
    xmax_now = MAX(xmax_now,dx2);
    ymax_now = MAX(ymax_now,dy2);
    zmax_now = MAX(zmax_now,dz2);
  }/*endfor*/
  xmax_now = 2.0*sqrt(xmax_now);
  ymax_now = 2.0*sqrt(ymax_now);
  zmax_now = 2.0*sqrt(zmax_now);

/*==========================================================================*/
/* II) Make sure these are within the tolerence */

  if(iperd==0 || iperd==4){
   if(xmax_now > (xmax-tol_edge_dist)){
    printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@@@@@@\n");
    printf("Particle x-separation too large for the CP box \n");
    printf("dx = %g A : box_x = %g A\n",
            xmax_now*0.529177,xmax*0.529177);
    printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@@@@@@\n");
    fflush(stdout);
    exit(1);
   }/*endifr*/
  }/*endif*/

  if(iperd<=1 || iperd==4){
   if(ymax_now > (ymax-tol_edge_dist)){
    printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@@@@@@\n");
    printf("Particle y-separation too large for the CP box \n");
    printf("dy = %g A : box_y = %g A\n",ymax_now*0.529177,ymax*0.529177);
    printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@@@@@@\n");
    fflush(stdout);
    exit(1);
   }/*endif*/
  }/*endif*/

  if(iperd<=2 || iperd==4){
   if(zmax_now > (zmax-tol_edge_dist)){
    printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@@@@@@\n");
    printf("Particle z-separation too large for the CP box \n");
    printf("dz = %g A : box_z = %g A\n",zmax_now*0.529177,zmax*0.529177);
    printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@@@@@@\n");
    fflush(stdout);
    exit(1);
   }/*endif*/
  }/*endif*/

/*==========================================================================*/
   }/*end routine */
/*==========================================================================*/





/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
/* Under CP with dualing, make sure your designated cp-particles have not   */
/* escaped the small grid                                                   */
/*==========================================================================*/

void cp_dual_check(CELL *cell,CLATOMS_INFO *clatoms_info,
                   CLATOMS_POS *clatoms_pos,
                   int *cp_atm_lst,double tol_edge_dist)

/*==========================================================================*/
  {/*Begin Routine*/
/*==========================================================================*/
/*            Local variable declarations                                   */

  int i,ind;
  double dx2,dy2,dz2;
  double xmax_now,ymax_now,zmax_now;

  int natm_tot    = clatoms_info->natm_tot;
  int nab_initio  = clatoms_info->nab_initio;
  double xmax     = cell->hmat_cp[1];  
  double ymax     = cell->hmat_cp[5];  
  double zmax     = cell->hmat_cp[9];  
  double xcm      = cell->cp_box_center[1];
  double ycm      = cell->cp_box_center[2];
  double zcm      = cell->cp_box_center[3];
  double *x       = clatoms_pos->x;
  double *y       = clatoms_pos->y;
  double *z       = clatoms_pos->z;
  int iperd       = cell->iperd;

/*==========================================================================*/
/* I) Get the maximum distances */
 
  xmax_now = 0.0;
  ymax_now = 0.0;
  zmax_now = 0.0;
  for(i=1;i<=nab_initio;i++){
    ind = cp_atm_lst[i];
    dx2 = (x[ind]-xcm)*(x[ind]-xcm);
    dy2 = (y[ind]-ycm)*(y[ind]-ycm);
    dz2 = (z[ind]-zcm)*(z[ind]-zcm);
    xmax_now = MAX(xmax_now,dx2);
    ymax_now = MAX(ymax_now,dy2);
    zmax_now = MAX(zmax_now,dz2);

  }/*endfor*/
  xmax_now = 2.0*sqrt(xmax_now);
  ymax_now = 2.0*sqrt(ymax_now);
  zmax_now = 2.0*sqrt(zmax_now);

/*==========================================================================*/
/* II) Make sure these are within the tolerence */

  if(iperd==0 || iperd==4){
   if(xmax_now > (xmax-tol_edge_dist)){
    printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@@@@@@\n");
    printf("Particle x-separation too large for the CP box \n");
    printf("dx = %g A : box_x = %g A\n",
            xmax_now*0.529177,xmax*0.529177);
    printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@@@@@@\n");
    fflush(stdout);
    exit(1);
   }/*endifr*/
  }/*endif*/

  if(iperd<=1 || iperd==4){
   if(ymax_now > (ymax-tol_edge_dist)){
    printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@@@@@@\n");
    printf("Particle y-separation too large for the CP box \n");
    printf("dy = %g A : box_y = %g A\n",ymax_now*0.529177,ymax*0.529177);
    printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@@@@@@\n");
    fflush(stdout);
    exit(1);
   }/*endif*/
  }/*endif*/

  if(iperd<=2 || iperd==4){
   if(zmax_now > (zmax-tol_edge_dist)){
    printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@@@@@@\n");
    printf("Particle z-separation too large for the CP box \n");
    printf("dz = %g A : box_z = %g A\n",zmax_now*0.529177,zmax*0.529177);
    printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@@@@@@\n");
    fflush(stdout);
    exit(1);
   }/*endif*/
  }/*endif*/

/*==========================================================================*/
   }/*end routine */
/*==========================================================================*/

