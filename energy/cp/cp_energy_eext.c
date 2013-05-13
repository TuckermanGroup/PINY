/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
/*                                                                          */
/*                         PI_MD:                                           */
/*             The future of simulation technology                          */
/*             ------------------------------------                         */
/*                   Module: cp_energy_pot.c                                */
/*                                                                          */
/* This routine calls the required force and PE routines                    */
/*                                                                          */
/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/



#include "standard_include.h"
#include "../typ_defs/typedefs_gen.h"
#include "../typ_defs/typedefs_class.h"
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


void control_cp_eext_recip(CLATOMS_INFO *clatoms_info,CLATOMS_POS *clatoms_pos,
                           CPCOEFFS_INFO *cpcoeffs_info,
                           CPCOEFFS_POS *cpcoeffs_pos,
                           CPEWALD *cpewald,CPSCR *cpscr,
                           CPOPTS *cpopts,PSEUDO *pseudo,EWD_SCR *ewd_scr,
                           ATOMMAPS *atommaps,CELL *cell,EWALD *ewald,
                           PTENS *ptens,
                           double *vrecip_ret,double *cp_enl_ret,
                           COMMUNICATE *communicate,FOR_SCR *for_scr,
                           int cp_dual_grid_opt,
                           PARA_FFT_PKG3D *cp_para_fft_pkg3d_lg)

/*==========================================================================*/
/*         Begin Routine                                                    */
   {/*Begin Routine*/
/*=======================================================================*/
/*         Local Variable declarations                                   */
#include "../typ_defs/typ_mask.h"

  int idual_switch;
  int i,j,iii,igh;
  int nlmtot,ntot_up,ntot_dn;
  int nl_max_kb,np_nlmax_kb,nl_max_gh,np_nlmax_gh,nl_max_all,np_nlmax_all;
  double vrecip,cp_enl,cp_enl_gh;
  double cpu1,cpu2;
/*----------------------------------------------------------------------*/
/*         Local Pointer declarations                                   */
     
  int iperd     = cell->iperd;

 /*-------------------------*/
 /* Pressure local pointers */
  int cp_ptens  = cpopts->cp_ptens_calc;
  double *pvten = ptens->pvten_tmp;

 /*------------------------------------------*/
 /* Non-local pseudopotential local pointers */

  int *np_nl            = pseudo->np_nl;
  int *np_nl_gh         = pseudo->np_nl_gh;

  int n_ang_max         = pseudo->n_ang_max;
  int n_ang_max_kb      = pseudo->n_ang_max_kb;
  int n_ang_max_gh      = pseudo->n_ang_max_gh;

  int n_rad_max         = pseudo->n_rad_max;
  double *vnlreal_up    = cpscr->cpscr_nonloc.vnlre_up;
  double *vnlimag_up    = cpscr->cpscr_nonloc.vnlim_up;
  double *vnlreal_dn    = cpscr->cpscr_nonloc.vnlre_dn;
  double *vnlimag_dn    = cpscr->cpscr_nonloc.vnlim_dn;

  double *dvnlreal_x_up = cpscr->cpscr_nonloc.dvnlre_x_up;
  double *dvnlreal_y_up = cpscr->cpscr_nonloc.dvnlre_y_up;
  double *dvnlreal_z_up = cpscr->cpscr_nonloc.dvnlre_z_up;
  double *dvnlimag_x_up = cpscr->cpscr_nonloc.dvnlim_x_up;
  double *dvnlimag_y_up = cpscr->cpscr_nonloc.dvnlim_y_up;
  double *dvnlimag_z_up = cpscr->cpscr_nonloc.dvnlim_z_up;
  double *dvnlreal_x_dn = cpscr->cpscr_nonloc.dvnlre_x_dn;
  double *dvnlreal_y_dn = cpscr->cpscr_nonloc.dvnlre_y_dn;
  double *dvnlreal_z_dn = cpscr->cpscr_nonloc.dvnlre_z_dn;
  double *dvnlimag_x_dn = cpscr->cpscr_nonloc.dvnlim_x_dn;
  double *dvnlimag_y_dn = cpscr->cpscr_nonloc.dvnlim_y_dn;
  double *dvnlimag_z_dn = cpscr->cpscr_nonloc.dvnlim_z_dn;

  double *dvnlreal_gxgx_up = cpscr->cpscr_nonloc.dvnlre_gxgx_up;
  double *dvnlimag_gxgx_up = cpscr->cpscr_nonloc.dvnlim_gxgx_up;
  double *dvnlreal_gygy_up = cpscr->cpscr_nonloc.dvnlre_gygy_up;
  double *dvnlimag_gygy_up = cpscr->cpscr_nonloc.dvnlim_gygy_up;
  double *dvnlreal_gzgz_up = cpscr->cpscr_nonloc.dvnlre_gzgz_up;
  double *dvnlimag_gzgz_up = cpscr->cpscr_nonloc.dvnlim_gzgz_up;
  double *dvnlreal_gxgy_up = cpscr->cpscr_nonloc.dvnlre_gxgy_up;
  double *dvnlimag_gxgy_up = cpscr->cpscr_nonloc.dvnlim_gxgy_up;
  double *dvnlreal_gygz_up = cpscr->cpscr_nonloc.dvnlre_gygz_up;
  double *dvnlimag_gygz_up = cpscr->cpscr_nonloc.dvnlim_gygz_up;
  double *dvnlreal_gxgz_up = cpscr->cpscr_nonloc.dvnlre_gxgz_up;
  double *dvnlimag_gxgz_up = cpscr->cpscr_nonloc.dvnlim_gxgz_up;
  double *dvnlreal_gxgx_dn = cpscr->cpscr_nonloc.dvnlre_gxgx_dn;
  double *dvnlimag_gxgx_dn = cpscr->cpscr_nonloc.dvnlim_gxgx_dn;
  double *dvnlreal_gygy_dn = cpscr->cpscr_nonloc.dvnlre_gygy_dn;
  double *dvnlimag_gygy_dn = cpscr->cpscr_nonloc.dvnlim_gygy_dn;
  double *dvnlreal_gzgz_dn = cpscr->cpscr_nonloc.dvnlre_gzgz_dn;
  double *dvnlimag_gzgz_dn = cpscr->cpscr_nonloc.dvnlim_gzgz_dn;
  double *dvnlreal_gxgy_dn = cpscr->cpscr_nonloc.dvnlre_gxgy_dn;
  double *dvnlimag_gxgy_dn = cpscr->cpscr_nonloc.dvnlim_gxgy_dn;
  double *dvnlreal_gygz_dn = cpscr->cpscr_nonloc.dvnlre_gygz_dn;
  double *dvnlimag_gygz_dn = cpscr->cpscr_nonloc.dvnlim_gygz_dn;
  double *dvnlreal_gxgz_dn = cpscr->cpscr_nonloc.dvnlre_gxgz_dn;
  double *dvnlimag_gxgz_dn = cpscr->cpscr_nonloc.dvnlim_gxgz_dn;

 /*----------------------*/
 /* Atom local pointers */
  int natm_tot             = clatoms_info->natm_tot;
  double *fx_tmp           = ewd_scr->fx;
  double *fy_tmp           = ewd_scr->fy;
  double *fz_tmp           = ewd_scr->fz;
  double *fx               = clatoms_pos->fx;
  double *fy               = clatoms_pos->fy;
  double *fz               = clatoms_pos->fz;
  double *x                = clatoms_pos->x;
  double *y                = clatoms_pos->y;
  double *z                = clatoms_pos->z;
  double *q                = clatoms_info->q;

  int hess_size;
  int hess_calc = clatoms_info->hess_calc;

  double *hess_xx         = clatoms_pos->hess_xx;
  double *hess_xy         = clatoms_pos->hess_xy;
  double *hess_xz         = clatoms_pos->hess_xz;
  double *hess_yy         = clatoms_pos->hess_yy;
  double *hess_yz         = clatoms_pos->hess_yz;
  double *hess_zz         = clatoms_pos->hess_zz;

 /*-----------------------------------*/
 /* Cp Option and form local pointers */
  int cp_lsda   = cpopts->cp_lsda;
  int icoef_orth_up  = cpcoeffs_pos->icoef_orth_up;
  int icoef_form_up  = cpcoeffs_pos->icoef_form_up;
  int ifcoef_form_up = cpcoeffs_pos->ifcoef_form_up;
  int icoef_orth_dn  = cpcoeffs_pos->icoef_orth_dn;
  int icoef_form_dn  = cpcoeffs_pos->icoef_form_dn;
  int ifcoef_form_dn = cpcoeffs_pos->ifcoef_form_dn;

 /*------------------------------*/
 /* Wave function local pointers */
  int nstate_up = cpcoeffs_info->nstate_up_proc;
  int nstate_dn = cpcoeffs_info->nstate_dn_proc;

 /*------------------------------*/
 /* Communciation local pointers */
  MPI_Comm comm_states = communicate->comm_states;
  int myid_state  = communicate->myid_state;
  int np_states   = communicate->np_states;
  int np_forc     = communicate->np_forc;


/*======================================================================*/
/* 0) Check the forms                                                   */

  if(icoef_orth_up!=1){
    printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
    printf("The Up coefficients must be in orthogonal form    \n");
    printf("on state processor %d in control_cp_pe_recip   \n",myid_state);
    printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
    fflush(stdout);exit(1);
  }/*endif*/
  if(cp_lsda==1 && nstate_dn != 0){
   if(icoef_orth_dn!=1){
    printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
    printf("The dn coefficients must be in orthogonal form    \n");
    printf("on state processor %d in control_cp_pe_recip   \n",myid_state);
    printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
    fflush(stdout);exit(1);
   }/*endif*/
  }/*endif*/

  if(np_states>1){
   if((icoef_form_up+ifcoef_form_up)!=0){
    printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
    printf("The up coefs and coef forces must not be in transposed form\n");
    printf("on state processor %d in control_cp_pe_recip   \n",myid_state);
    printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
    fflush(stdout);exit(1);
   }/*endif*/
   if(cp_lsda==1 && nstate_dn != 0){
    if((icoef_form_dn+ifcoef_form_dn)!=0){
     printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
     printf("The dn coefs and coef forces must not be in transposed form\n");
     printf("on state processor %d in control_cp_pe_recip   \n",myid_state);
     printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
     fflush(stdout);exit(1);
    }/*endif*/
   }/*endif*/
  }/*endif*/

/*======================================================================*/
/* II) Malloc a hessian scratch vector if necessary                     */

  if(hess_calc == 3 && np_states > 1){
    hess_size = natm_tot*natm_tot;
  }

/*======================================================================*/
/* III) Determine the maximum open non-local angular momentum channel   */
/*      for Kleinman-Bylander and Goedecker pseudo potentials           */

  nl_max_kb = -1;
  for(i=1;i<=(n_ang_max_kb+1);i++){
   if(np_nl[i]>0){nl_max_kb=i-1;}
  }/*endfor*/

  nl_max_gh = -1;
  for(i=1;i<=(n_ang_max_gh+1);i++){
   if(np_nl_gh[i]>0){nl_max_gh=i-1;}
  }/*endfor*/

/*======================================================================*/
/* IV) Determine the maximum number of atoms in any                     */
/*       open angular momentum channel                                  */

  np_nlmax_kb = 1;
  for(i = 1;i<=(nl_max_kb+1);i++){
   np_nlmax_kb = MAX(np_nlmax_kb,np_nl[i]);
  }/*endfor*/

  np_nlmax_gh = 1;
  for(i = 1;i<=(nl_max_gh+1);i++){
   np_nlmax_gh = MAX(np_nlmax_gh,np_nl_gh[i]);
  }/*endfor*/

  np_nlmax_all = (np_nlmax_gh > np_nlmax_kb ? np_nlmax_gh : np_nlmax_kb);
 
/*======================================================================*/
/* V) Zero the non-local tensors                                        */

  nl_max_all = (nl_max_gh > nl_max_kb ? nl_max_gh : nl_max_kb);
  nlmtot = (nl_max_all+1)*(nl_max_all+1);

  ntot_up = nstate_up*np_nlmax_all*nlmtot*n_rad_max;
  ntot_dn = 0;



  if(nl_max_all >= 0){
   for(i=1;i<=ntot_up;i++){
    vnlreal_up[i] = 0.0;
    vnlimag_up[i] = 0.0;
    dvnlreal_x_up[i] = 0.0;
    dvnlreal_y_up[i] = 0.0;
    dvnlreal_z_up[i] = 0.0;
    dvnlimag_x_up[i] = 0.0;
    dvnlimag_y_up[i] = 0.0;
    dvnlimag_z_up[i] = 0.0;
   }/*endfor*/

   if(cp_ptens==1 || hess_calc == 3){
    for(i=1;i<=ntot_up;i++){
     dvnlreal_gxgx_up[i] = 0.0;
     dvnlreal_gzgz_up[i] = 0.0;
     dvnlreal_gygy_up[i] = 0.0;
     dvnlreal_gxgy_up[i] = 0.0;
     dvnlreal_gxgz_up[i] = 0.0;
     dvnlreal_gygz_up[i] = 0.0;
 
     dvnlimag_gxgx_up[i] = 0.0;
     dvnlimag_gxgy_up[i] = 0.0;
     dvnlimag_gygy_up[i] = 0.0;
     dvnlimag_gxgz_up[i] = 0.0;
     dvnlimag_gygz_up[i] = 0.0;
     dvnlimag_gzgz_up[i] = 0.0;
    }/*endfor*/
   }/*endif:ptens*/

   if(cp_lsda==1){

    ntot_dn = nstate_dn*np_nlmax_all*nlmtot*n_rad_max;
    for(i=1;i<=ntot_dn;i++){
     vnlreal_dn[i]    = 0.0;
     vnlimag_dn[i]    = 0.0;
     dvnlreal_x_dn[i] = 0.0;
     dvnlreal_y_dn[i] = 0.0;
     dvnlreal_z_dn[i] = 0.0;
     dvnlimag_x_dn[i] = 0.0;
     dvnlimag_y_dn[i] = 0.0;
     dvnlimag_z_dn[i] = 0.0;
    }/*endfor*/

    if(cp_ptens==1 || hess_calc == 3){
     for(i=1;i<=ntot_dn;i++){
      dvnlreal_gxgx_dn[i] = 0.0;
      dvnlreal_gxgy_dn[i] = 0.0;
      dvnlreal_gxgz_dn[i] = 0.0;
      dvnlreal_gygy_dn[i] = 0.0;
      dvnlreal_gygz_dn[i] = 0.0;
      dvnlreal_gzgz_dn[i] = 0.0;

      dvnlimag_gxgx_dn[i] = 0.0;
      dvnlimag_gxgy_dn[i] = 0.0;
      dvnlimag_gxgz_dn[i] = 0.0;
      dvnlimag_gygy_dn[i] = 0.0;
      dvnlimag_gygz_dn[i] = 0.0;
      dvnlimag_gzgz_dn[i] = 0.0;
     }/*endfor*/
    }/*endif:ptens*/

   }/*endif:lsda*/

  }/*endif : non-local potential on*/

/*======================================================================*/
/* VI) Perform the ewald sum/ cp local pseudopotential calculation      */

  idual_switch = 0; /* cp_dual_grid_opt<=1 : get vrecip vext on dense grid */
                    /* cp_dual_grid_opt==2 : get vrecip vext on sparse grid*/
#ifdef TIME_CP
    if(np_states>1){Barrier(comm_states);}
    cputime(&cpu1);
#endif

  vrecip = 0.0;

 if(cpscr->cpscr_atom_pme.pme_on == 1 && cp_dual_grid_opt ==2){

  control_ewd_loc_pme(clatoms_info, clatoms_pos, cell, ptens, ewald, cpewald, 
                  cpscr, pseudo, ewd_scr, cpopts, atommaps, 
                  &vrecip, &(cpcoeffs_info->pseud_hess_loc),communicate,
                  for_scr,cp_dual_grid_opt,idual_switch,
                  cp_para_fft_pkg3d_lg);
 }else{ 

  control_ewd_loc(clatoms_info, clatoms_pos, cell, ptens, ewald, cpewald, 
                  cpscr, pseudo, ewd_scr, cpopts, atommaps, 
                  &vrecip, &(cpcoeffs_info->pseud_hess_loc),communicate,
                  for_scr,cp_dual_grid_opt,idual_switch);
 }/*endif*/


#ifdef TIME_CP
    cputime(&cpu2);
    par_cpu_vomit((cpu2-cpu1),comm_states,np_states,myid_state,
                       "1 control_ewd_loc");
#endif


  if(cp_dual_grid_opt == 2){
#ifdef TIME_CP
    if(np_states>1){Barrier(comm_states);}
    cputime(&cpu1);
#endif

   idual_switch = 1; /*get vext on small dense grid */
   control_ewd_loc(clatoms_info, clatoms_pos, cell, ptens, ewald, cpewald, 
                   cpscr, pseudo, ewd_scr, cpopts, atommaps, 
                   &vrecip, &(cpcoeffs_info->pseud_hess_loc),communicate,
                   for_scr,cp_dual_grid_opt,idual_switch);

#ifdef TIME_CP
    cputime(&cpu2);
    par_cpu_vomit((cpu2-cpu1),comm_states,np_states,myid_state,
                       "2 control_ewd_loc");
#endif

  }/*endif cp_dual_grid_opt*/
 
/*======================================================================*/
/* VII) Get the nl pe, pvten and particle forces then the coef forces   */

   cp_enl    = 0.00000000;
   cp_enl_gh = 0.00000000;

 if((nl_max_all >= 0) && (n_rad_max>1)){
      non_loc_chng_ord(clatoms_pos,clatoms_info,
                       atommaps, pseudo,ewd_scr,for_scr,1);
 }/*endif*/

/*-------------------------------------------------------------------------*/
/* A) KB/Goedecker NLs */

  if( (nl_max_kb >= 0) && ((ntot_up+ntot_dn)>0) ){
#ifdef TIME_CP
    if(np_states>1){Barrier(comm_states);}
    cputime(&cpu1);
#endif
    control_ewd_non_loc(clatoms_info,clatoms_pos,cpcoeffs_info,cpcoeffs_pos,
                        cell,ptens,cpewald,cpscr,pseudo,ewd_scr, 
                        cpopts,atommaps,communicate,for_scr);

#ifdef TIME_CP
    cputime(&cpu2);
    par_cpu_vomit((cpu2-cpu1),comm_states,np_states,myid_state,
                       "control_ewd_non_loc");
#endif

  }else{
    get_ak2_sm(cpewald,cell);
  }/*endif*/


  if((nl_max_kb >= 0)&&((ntot_up+ntot_dn)>0)&&(pseudo->np_nonloc_cp_box_kb>0) ){
#ifdef TIME_CP
    if(np_states>1){Barrier(comm_states);}
    cputime(&cpu1);
#endif

    getnl_pot_pv_fatm(clatoms_info,clatoms_pos,cell,cpcoeffs_info,
                      cpscr,ewd_scr,cpopts,pseudo,atommaps,&cp_enl,
                      np_nlmax_kb,pvten);

#ifdef TIME_CP
    cputime(&cpu2);
    par_cpu_vomit((cpu2-cpu1),comm_states,np_states,myid_state,
                       "getnl_pot_pv_fatm");
#endif

#ifdef TIME_CP
    if(np_states>1){Barrier(comm_states);}
    cputime(&cpu1);
#endif

    getnl_fcoef(clatoms_info,clatoms_pos,cpcoeffs_info,cpcoeffs_pos,
                  cpscr,ewd_scr,cpopts,pseudo,cpewald,atommaps,
                  cell,np_nlmax_kb,pvten,for_scr);

#ifdef TIME_CP
    cputime(&cpu2);
    par_cpu_vomit((cpu2-cpu1),comm_states,np_states,myid_state,
                       "getnl_fcoef");
#endif

  }/*endif*/
/*-------------------------------------------------------------------------*/
/* B) Gauss-Hermite NLs                                                    */

  if((nl_max_gh >= 0)&&((ntot_up+ntot_dn)>0)&&(pseudo->np_nonloc_cp_box_gh>0)){
      if(cp_ptens==1){
       printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
       printf(" CP-PTENS is not implemented for Gauss-Hermite nonlocality \n");
       printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
       fflush(stdout);
       exit(1);
  }/* endif */

 for(igh=1;igh<=pseudo->ngh;igh++){

  cp_enl_gh = 0.00000;

  if(nl_max_all >= 0){
   for(i=1;i<=ntot_up;i++){
    vnlreal_up[i] = 0.0;
    vnlimag_up[i] = 0.0;
    dvnlreal_x_up[i] = 0.0;
    dvnlreal_y_up[i] = 0.0;
    dvnlreal_z_up[i] = 0.0;
    dvnlimag_x_up[i] = 0.0;
    dvnlimag_y_up[i] = 0.0;
    dvnlimag_z_up[i] = 0.0;
   }/*endfor*/

   if(cp_ptens==1 || hess_calc == 3){
    for(i=1;i<=ntot_up;i++){
     dvnlreal_gxgx_up[i] = 0.0;
     dvnlreal_gzgz_up[i] = 0.0;
     dvnlreal_gygy_up[i] = 0.0;
     dvnlreal_gxgy_up[i] = 0.0;
     dvnlreal_gxgz_up[i] = 0.0;
     dvnlreal_gygz_up[i] = 0.0;
 
     dvnlimag_gxgx_up[i] = 0.0;
     dvnlimag_gxgy_up[i] = 0.0;
     dvnlimag_gygy_up[i] = 0.0;
     dvnlimag_gxgz_up[i] = 0.0;
     dvnlimag_gygz_up[i] = 0.0;
     dvnlimag_gzgz_up[i] = 0.0;
    }/*endfor*/
   }/*endif:ptens*/

   if(cp_lsda==1){

    ntot_dn = nstate_dn*np_nlmax_all*nlmtot*n_rad_max;
    for(i=1;i<=ntot_dn;i++){
     vnlreal_dn[i]    = 0.0;
     vnlimag_dn[i]    = 0.0;
     dvnlreal_x_dn[i] = 0.0;
     dvnlreal_y_dn[i] = 0.0;
     dvnlreal_z_dn[i] = 0.0;
     dvnlimag_x_dn[i] = 0.0;
     dvnlimag_y_dn[i] = 0.0;
     dvnlimag_z_dn[i] = 0.0;
    }/*endfor*/

    if(cp_ptens==1 || hess_calc == 3){
     for(i=1;i<=ntot_dn;i++){
      dvnlreal_gxgx_dn[i] = 0.0;
      dvnlreal_gxgy_dn[i] = 0.0;
      dvnlreal_gxgz_dn[i] = 0.0;
      dvnlreal_gygy_dn[i] = 0.0;
      dvnlreal_gygz_dn[i] = 0.0;
      dvnlreal_gzgz_dn[i] = 0.0;

      dvnlimag_gxgx_dn[i] = 0.0;
      dvnlimag_gxgy_dn[i] = 0.0;
      dvnlimag_gxgz_dn[i] = 0.0;
      dvnlimag_gygy_dn[i] = 0.0;
      dvnlimag_gygz_dn[i] = 0.0;
      dvnlimag_gzgz_dn[i] = 0.0;
     }/*endfor*/
    }/*endif:ptens*/

   }/*endif:lsda*/
  }/*endif : non-local potential on*/

   /* Create Zilm(alpha) and Xilm(alpha)= dZilm/dR  */

       control_nonloc_gh(clatoms_info,clatoms_pos,cpcoeffs_info,cpcoeffs_pos,
                         cell,ptens,cpewald,cpscr,pseudo,ewd_scr,  
                         cpopts,atommaps,communicate,for_scr,
                         pseudo->rgh[igh]);

 /*i) Calculate the nonlocal Energy  */
 /*ii) and the contribution to the forces on the atoms due to the non-local term*/

       getnl_pot_pv_fatm_gh(clatoms_info,clatoms_pos,cell,cpcoeffs_info,
                            cpscr,ewd_scr,cpopts,pseudo,atommaps,
                            &cp_enl_gh,np_nlmax_gh,pvten,pseudo->wgh,igh);

   /*Calculate nonlocal contribution to the coefficient forces */

       getnl_fcoef_gh(clatoms_info,clatoms_pos,cpcoeffs_info,cpcoeffs_pos,
                      cpscr,ewd_scr,cpopts,pseudo,cpewald,atommaps,
                      cell,np_nlmax_gh,pvten,for_scr,
                      pseudo->rgh[igh],pseudo->wgh,igh);

        *cp_enl_ret += cp_enl_gh;
      }/*endfor igh gauss-hermite integration points */

    }/*endif gauss-hermit*/

    if((nl_max_all >= 0)&&((ntot_up+ntot_dn)>0)&&(n_rad_max>1)){
      non_loc_restore_ord(clatoms_pos,clatoms_info,
                          atommaps, pseudo,ewd_scr,for_scr);
    }/*endif*/

/*======================================================================*/
/* VIII) Assign the potential energy                                    */

  *vrecip_ret += vrecip;
  *cp_enl_ret += cp_enl;

/*======================================================================*/
/* IX) Reduce particle forces if necessary                              */

  if(np_states>1 && np_forc == 1){
    reduce_cp_atm_forc(natm_tot,fx,fy,fz,fx_tmp,fy_tmp,fz_tmp,
                       comm_states,myid_state);
  }/* endif:npstates */

/*======================================================================*/
/* X) Reduce particle hessian if necessary                              */

  if(hess_calc == 3 && np_states>1){
    reduce_cp_hess_stuff(hess_xx,hess_yy,hess_zz,
                         hess_xy,hess_xz,hess_yz,hess_size,
                         myid_state,comm_states);
  }/*endif*/

/*======================================================================*/
    }/*end routine*/
/*======================================================================*/





/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void control_ewd_loc(CLATOMS_INFO *clatoms_info,CLATOMS_POS *clatoms_pos,
                     CELL *cell, PTENS *ptens, EWALD *ewald, CPEWALD *cpewald, 
                     CPSCR *cpscr, PSEUDO *pseudo, EWD_SCR *ewd_scr,  
                     CPOPTS *cpopts, ATOMMAPS *atommaps, double *vrecip_ret,
                     double *pseud_hess_loc,COMMUNICATE *communicate,
                     FOR_SCR *for_scr,int cp_dual_grid_opt,int idual_switch)

/*==========================================================================*/
/*         Begin Routine                                                    */
   {/*Begin Routine*/
/*=======================================================================*/
/*         Local Variable declarations                                   */
     
#include "../typ_defs/typ_mask.h"

  int idens_opt,ipseud_opt;
  int istart,ngo,irem,idiv;
  int ipart,jpart,iii,itype,i;
  int icount,koff,natm_use;
  int hess_ind;

  double falp2,falp_clus2,vol,rvol,pivol,fpi,arg,q_sum1,bgr;
  double aka,akb,akc,xk,yk,zk,atemp,btemp,ctemp;
  double xtemp,ytemp,ztemp;
  double sumr,sumi,g2,g4,preg,prep,tpi,pi,g;
  double sumr_h,sumi_h;
  double srx,sry,srz,six,siy,siz,temp,smag;
  double vrecip;
  double dx,dy,dz;
  double sx,sy,sz;
  double asx,asy,asz;
  double phase;
  double argp,fargp,argm,fargm,area;

/*--------------------------------------------*/
/*         Local Pointer declarations         */

  /*------------------*/
  /* Atom information */
  int natm_tot      = clatoms_info->natm_tot; 
  int hess_calc     = clatoms_info->hess_calc;
  double *x         = clatoms_pos->x;
  double *y         = clatoms_pos->y;
  double *z         = clatoms_pos->z;
  double *fx        = clatoms_pos->fx;     
  double *fy        = clatoms_pos->fy;
  double *fz        = clatoms_pos->fz;
  double *fx_tmp    = ewd_scr->fx2;     
  double *fy_tmp    = ewd_scr->fy2;     
  double *fz_tmp    = ewd_scr->fz2;     
  double *q         = clatoms_info->q;
  double *hess_xx   = clatoms_pos->hess_xx;
  double *hess_xy   = clatoms_pos->hess_xy;
  double *hess_xz   = clatoms_pos->hess_xz;
  double *hess_yy   = clatoms_pos->hess_yy;
  double *hess_yz   = clatoms_pos->hess_yz;
  double *hess_zz   = clatoms_pos->hess_zz;
  int natm_typ      = atommaps->natm_typ;
  int *index_atm    = for_scr->index_atm;
  int *iatm_typ;            /*Assigned below based on flags */
  int *iatm_typ_full;

  /*--------------------------------*/
  /* Cell and pressure information */
  int iperd                 = cell->iperd;
  int cp_ptens              = cpopts->cp_ptens_calc;
  double *pvten             = ptens->pvten_tmp;
  double *cp_box_center     = cell->cp_box_center;
  double *cp_box_center_rel = cell->cp_box_center_rel;
  double *hmat;             /* Assigned below based on flags */
  double *hmati;
  double *hmat_big          = cell->hmat;
  double *hmati_big         = cell->hmati;

  /*----------------------*/
  /* G-vector information */
  int *kastore;             /* Assigned below based on flags */
  int *kbstore;
  int *kcstore;
  int *ibreak1;
  int *ibreak2;
  double *vextr;
  double *vexti;
  double *dvextr;
  double *dvexti;
  double *rhocr;
  double *rhoci;
  double *ak2;
  int nktot;

  /*----------------------------------------------*/
  /* Pseudo-potential and Reduced Periodicity info*/
  int nsplin_g         = pseudo->nsplin_g;
  int n_rad_max        = pseudo->n_rad_max;
  double *clus_corr_r  = ewald->clus_corr_r;
  double *dclus_corr_r = ewald->dclus_corr_r;
  double alpha_conv_dual = pseudo->alpha_conv_dual;
  double dg_spl        = pseudo->dg_spl;
  double gmin_spl      = pseudo->gmin_spl;
  double *vps0         = pseudo->vps0;
  double *vps1         = pseudo->vps1;
  double *vps2         = pseudo->vps2;
  double *vps3         = pseudo->vps3;
  double *dvps0        = pseudo->dvps0;
  double *dvps1        = pseudo->dvps1;
  double *dvps2        = pseudo->dvps2;
  double *dvps3        = pseudo->dvps3;
  double *gzvps        = pseudo->gzvps;
  double *q_pseud      = pseudo->q_pseud;
  int n_ang_max        = pseudo->n_ang_max;
  int *loc_opt         = pseudo->loc_opt;
  int np_loc_cp_box    = pseudo->np_loc_cp_box;
  int *ip_loc_cp_box   = pseudo->ip_loc_cp_box;

  /*---------------------------------*/
  /* Ewald and ewald scr information */
  double alp_ewald  = ewald->alp_ewd;
  double alp_clus   = ewald->alp_clus;
  double *cossc     = ewd_scr->cossc;  
  double *sinsc     = ewd_scr->sinsc;
  double *helr      = ewd_scr->helr;   
  double *heli      = ewd_scr->heli;
  double *vtemp     = ewd_scr->temp;
  double *dvtemp    = ewd_scr->vtemp_now;
  double *ewd_scr_x = ewd_scr->x;
  double *ewd_scr_y = ewd_scr->y;
  double *ewd_scr_z = ewd_scr->z;
  double *q_tmp     = ewd_scr->q;

  /*---------------------------*/
  /* Communication information */
  int myid_state    = communicate->myid_state;
  int np_states     = communicate->np_states;
  MPI_Comm comm     = communicate->comm_states;

/*======================================================================*/
/* 0) Assign local pointers                                             */

  if(cp_dual_grid_opt < 2 || idual_switch == 0){
    /* large sparse grid when cp_dual_grid_opt == 2*/
    idens_opt = 0;
    ipseud_opt= (cp_dual_grid_opt==2 ? 0 : 1);
    kastore   = ewald->kastr;
    kbstore   = ewald->kbstr;
    kcstore   = ewald->kcstr;
    ibreak1   = ewald->ibrk1;
    ibreak2   = ewald->ibrk2;
    vextr     = cpscr->cpscr_loc.vextr;
    vexti     = cpscr->cpscr_loc.vexti;
    dvextr    = cpscr->cpscr_loc.dvextr;
    dvexti    = cpscr->cpscr_loc.dvexti;
    rhocr     = cpscr->cpscr_rho.rhocr_up;
    rhoci     = cpscr->cpscr_rho.rhoci_up;
    ak2       = cpewald->ak2;
    nktot     = ewald->nktot;
    hmat      = cell->hmat;
    hmati     = cell->hmati;
    natm_use  = natm_tot;
    iatm_typ  = atommaps->iatm_atm_typ;
  }else{
    /* small dense grid */
    idens_opt = 1;
    ipseud_opt= 1;
    kastore   = cpewald->kastr_dens_cp_box;
    kbstore   = cpewald->kbstr_dens_cp_box;
    kcstore   = cpewald->kcstr_dens_cp_box;
    ibreak1   = cpewald->ibrk1_dens_cp_box; /*DY edit*/
    ibreak2   = cpewald->ibrk2_dens_cp_box; /*DY edit*/
    vextr     = cpscr->cpscr_loc.vextr_dens_cp_box;
    vexti     = cpscr->cpscr_loc.vexti_dens_cp_box;
    rhocr     = cpscr->cpscr_rho.rhocr_up_dens_cp_box;
    rhoci     = cpscr->cpscr_rho.rhoci_up_dens_cp_box;
    ak2       = cpewald->ak2_dens_cp_box;
    nktot     = cpewald->nktot_dens_cp_box;
    hmat      = cell->hmat_cp;
    hmati     = cell->hmati_cp;
    natm_use  = np_loc_cp_box;
    iatm_typ  = for_scr->iexcl;
    iatm_typ_full = atommaps->iatm_atm_typ;
 }/*endif*/

/*======================================================================*/
/* I) Get more useful constants                                         */

  pi    = M_PI;
  tpi   = 2.0*pi;
  fpi   = 4.0*pi;

  vol   = getdeth(hmat);
  rvol  = 1.0/vol;
  pivol = vol/4.0/pi;
  falp2 = 4.0*alp_ewald*alp_ewald;
  falp_clus2 = 4.0*alp_clus*alp_clus;

/*======================================================================*/
/* II) Find cos and sin of sc components of the particles               */
/*    ( hmnati rvec = svec   r=(x,y,z) s=(a,b,c) )                       */

  for(ipart=1;ipart<=natm_use;ipart++){
    fx_tmp[ipart] = 0.0;
    fy_tmp[ipart] = 0.0;
    fz_tmp[ipart] = 0.0;
  }/*endfor*/

  if(idens_opt==0){
    for(ipart=1;ipart<=natm_use;ipart++){
      xtemp = x[ipart];
      ytemp = y[ipart];
      ztemp = z[ipart];

      q_tmp[ipart]     = q[ipart];
      ewd_scr_x[ipart] = xtemp*hmati[1]
                       + ytemp*hmati[4]
                       + ztemp*hmati[7];
      ewd_scr_y[ipart] = xtemp*hmati[2]
                       + ytemp*hmati[5]
                       + ztemp*hmati[8];
      ewd_scr_z[ipart] = xtemp*hmati[3]
                       + ytemp*hmati[6]
                       + ztemp*hmati[9];
      ctemp            = ewd_scr_z[ipart]*tpi;
      cossc[ipart]     = cos(ctemp);
      sinsc[ipart]     = sin(ctemp);
    }/*endfor*/
  }else{
    for(ipart=1;ipart<=natm_use;ipart++){
      iatm_typ[ipart]  = iatm_typ_full[ip_loc_cp_box[ipart]];
      q_tmp[ipart]     = q[ip_loc_cp_box[ipart]];
      dx               = x[ip_loc_cp_box[ipart]] - cp_box_center[1];
      dy               = y[ip_loc_cp_box[ipart]] - cp_box_center[2];
      dz               = z[ip_loc_cp_box[ipart]] - cp_box_center[3];

      asx              = dx*hmati_big[1]+dy*hmati_big[4]+dz*hmati_big[7];
      asy              = dx*hmati_big[2]+dy*hmati_big[5]+dz*hmati_big[8];
      asz              = dx*hmati_big[3]+dy*hmati_big[6]+dz*hmati_big[9];
      sx               = asx - NINT(asx);
      sy               = asy - NINT(asy);
      sz               = asz - NINT(asz);
      dx               = sx*hmat_big[1]+sy*hmat_big[4]+sz*hmat_big[7];
      dy               = sx*hmat_big[2]+sy*hmat_big[5]+sz*hmat_big[8];
      dz               = sx*hmat_big[3]+sy*hmat_big[6]+sz*hmat_big[9];

      xtemp            = dx + cp_box_center_rel[1];
      ytemp            = dy + cp_box_center_rel[2];
      ztemp            = dz + cp_box_center_rel[3];
      ewd_scr_x[ipart] = xtemp*hmati[1]
                       + ytemp*hmati[4]
                       + ztemp*hmati[7];
      ewd_scr_y[ipart] = xtemp*hmati[2]
                       + ytemp*hmati[5]
                       + ztemp*hmati[8];
      ewd_scr_z[ipart] = xtemp*hmati[3]
                       + ytemp*hmati[6]
                       + ztemp*hmati[9];
      ctemp            = ewd_scr_z[ipart]*tpi;
      cossc[ipart]     = cos(ctemp);
      sinsc[ipart]     = sin(ctemp);
    }/*endfor*/
  }/*endif*/

/*======================================================================*/
/*======================================================================*/
/* Perform the ewald sum/ cp-potential calculation                      */

  vrecip  = 0.0;
  idiv    = (nktot+1)/np_states;
  irem    = (nktot+1) % np_states;
  ngo     = (myid_state <  irem ? idiv+1 : idiv);
  istart  = (myid_state <= irem ? myid_state*(idiv+1)+1 : 
                                irem*(idiv+1)+1+(myid_state-irem)*idiv);
  koff    = istart-1;
  if(np_states==myid_state+1){ngo--;}

  for(icount=1;icount<=ngo;icount++){

/*======================================================================*/
/* I) Get the k vectors                                                 */

   aka = (double)(kastore[(icount+koff)]);
   akb = (double)(kbstore[(icount+koff)]);
   akc = (double)(kcstore[(icount+koff)]);
   xk = (aka*hmati[1]+akb*hmati[2]+akc*hmati[3])*tpi;
   yk = (aka*hmati[4]+akb*hmati[5]+akc*hmati[6])*tpi;
   zk = (aka*hmati[7]+akb*hmati[8]+akc*hmati[9])*tpi;
   g2 = xk*xk+yk*yk+zk*zk;
   g4 = g2*g2;
   g  = sqrt(g2);
   ak2[icount] = g2;

/*======================================================================*/
/* II) If break point number one or you are just starting out calculate */
/*     the helpful vectors                                              */

   if(ibreak1[(icount+koff)]==1||icount==1){
    for(ipart=1;ipart<=natm_use;ipart++){
      atemp = ewd_scr_x[ipart];
      btemp = ewd_scr_y[ipart];
      ctemp = ewd_scr_z[ipart];
      arg = (aka*atemp + akb*btemp + akc*ctemp)*tpi;
      helr[ipart] = cos(arg);
      heli[ipart] = sin(arg);
    }/*endfor*/
   }/*endif*/

/*======================================================================*/
/* III) Get the external potential                                      */
/*               (interaction of electron with particles)               */

   if(ipseud_opt==1){
    for(itype=1;itype<=natm_typ;itype++){
      index_atm[itype] =  (itype-1)*nsplin_g*(n_ang_max+1)*n_rad_max
                       +  loc_opt[itype]*nsplin_g*n_rad_max;
    }/*endfor*/

    get_vpsnow(index_atm,nsplin_g,gmin_spl,dg_spl,g,
                vps0,vps1,vps2,vps3,vtemp,iatm_typ,natm_typ,natm_use,1); 

/*-------------------------*/
/* charge correction */
     for(i=1; i<= natm_use; i++){
          vtemp[i] += -fpi*(q_tmp[i] - q_pseud[iatm_typ[i]])/g2;
     }
   }else{  /*q_temp is q */
     get_vpslong(natm_use,vtemp,g2,q_tmp,alpha_conv_dual,pivol);    
   }/*endif*/

 /*----------------------------------------------------------------------*/
 /* Cluster boundary condition correction                                */

   if( (iperd != 3) && (idens_opt==0)){
     for(ipart=1;ipart<=natm_use;ipart++){
       vtemp[ipart] -= q[ipart]*clus_corr_r[icount];
     }/* endfor */
   }/* endif cluster boundary conditions */
 
/*----------------------------------------------------------------------*/

   if( (cp_ptens==1) && (idens_opt==0) ){
     if(ipseud_opt == 1){
       get_vpsnow(index_atm,nsplin_g,gmin_spl,dg_spl,g,
                  dvps0,dvps1,dvps2,dvps3,dvtemp,iatm_typ,natm_typ,natm_use,1);
       for(i=1; i<= natm_use; i++){
         dvtemp[i] += 2.0*fpi*(q_tmp[i]-q_pseud[iatm_typ[i]])/(g2*g2);
       }/* endfor */
     } else {
       get_dvpslong(natm_use,dvtemp,g2,q_tmp,alpha_conv_dual,pivol);
     }/*endif ipseud_opt */
   }/*endif cp_ptens*/

    vextr[icount]  =  ddot1(natm_use,helr,1,vtemp,1)*rvol;
    vexti[icount]  = -ddot1(natm_use,heli,1,vtemp,1)*rvol;

   if( (cp_ptens==1) && (idens_opt==0) ) {
     dvextr[icount] =  ddot1(natm_use,helr,1,dvtemp,1)*rvol;
     dvexti[icount] = -ddot1(natm_use,heli,1,dvtemp,1)*rvol;
   }/*endif*/
 
/*======================================================================*/
/* IV) Get the real and imag parts of the structure factor              */

   if(idens_opt==0){/*create charge weighted structure factor*/
     sumr = ddot1(natm_use,helr,1,q,1); 
     sumi = ddot1(natm_use,heli,1,q,1);
     smag = sumr*sumr+sumi*sumi;
   }/*endif*/

/*======================================================================*/
/* V) Use the stucture factor to get the                                */
/*      classical potential energy and the pressure*volume tensor       */
 
   preg    = exp(-g2/falp2)/(g2*pivol);
   prep    = -2.0*preg*(g2/falp2+1.0)/g2;

   if(iperd == 2 && ((kastore[icount+koff] == 0) && (kbstore[icount+koff] == 0))){
     phase = cos(0.5*zk*hmat[9]);
     preg  += phase*(1.0-exp(-g2/falp_clus2))/(g2*pivol);
   }
   prep    = -2.0*preg*(g2/falp2+1.0)/g2;
   if(iperd == 2 && ((kastore[icount+koff] == 0) && (kbstore[icount+koff] == 0))){
     phase = cos(0.5*zk*hmat[9]);
     prep += 2.0*phase/(g2*pivol*falp_clus2);
   }

   if( (iperd>0) && (iperd !=3) && (idens_opt==0) ){
     preg += clus_corr_r[icount]*rvol;
     prep += dclus_corr_r[icount]*rvol;
   }/*endif*/

   if( (iperd>0) &&(idens_opt==0) ){vrecip  = vrecip + smag*preg; }

   prep    = prep*smag; 
   if((cp_ptens==1) && (idens_opt==0) ) {
     pvten[1] = pvten[1] + prep*xk*xk;
     pvten[5] = pvten[5] + prep*yk*yk;
     pvten[9] = pvten[9] + prep*zk*zk;
     pvten[4] = pvten[4] + prep*xk*yk;
     pvten[2] = pvten[2] + prep*xk*yk;
     pvten[7] = pvten[7] + prep*xk*zk;
     pvten[3] = pvten[3] + prep*xk*zk;
     pvten[8] = pvten[8] + prep*yk*zk;
     pvten[6] = pvten[6] + prep*yk*zk;
   }/*endif*/

/*======================================================================*/
/* VI) Get the force on the particles, checking for CBC                 */

 
 if( (iperd == 0) || (idens_opt==1) ) {

    for(ipart=1;ipart<=natm_use;ipart++){
     srx = xk*(2.0*rhocr[icount]*vtemp[ipart]*rvol);
     sry = yk*(2.0*rhocr[icount]*vtemp[ipart]*rvol);
     srz = zk*(2.0*rhocr[icount]*vtemp[ipart]*rvol);
     six = xk*(-2.0*rhoci[icount]*vtemp[ipart]*rvol);
     siy = yk*(-2.0*rhoci[icount]*vtemp[ipart]*rvol);
     siz = zk*(-2.0*rhoci[icount]*vtemp[ipart]*rvol);
     fx_tmp[ipart] += (srx*heli[ipart] - six*helr[ipart]);
     fy_tmp[ipart] += (sry*heli[ipart] - siy*helr[ipart]);
     fz_tmp[ipart] += (srz*heli[ipart] - siz*helr[ipart]);
    }/* endfor */

  }else{

    sumr_h = sumr;
    sumi_h = sumi;
    sumr = sumr*preg*2.0;
    sumi = sumi*preg*2.0;
    for(ipart=1;ipart<=natm_use;ipart++){
      srx = xk*(sumr*q[ipart]+2.0*rhocr[icount]*vtemp[ipart]*rvol);
      sry = yk*(sumr*q[ipart]+2.0*rhocr[icount]*vtemp[ipart]*rvol);
      srz = zk*(sumr*q[ipart]+2.0*rhocr[icount]*vtemp[ipart]*rvol);
      six = xk*(sumi*q[ipart]-2.0*rhoci[icount]*vtemp[ipart]*rvol);
      siy = yk*(sumi*q[ipart]-2.0*rhoci[icount]*vtemp[ipart]*rvol);
      siz = zk*(sumi*q[ipart]-2.0*rhoci[icount]*vtemp[ipart]*rvol);
      fx_tmp[ipart] += (srx*heli[ipart]  - six*helr[ipart]);
      fy_tmp[ipart] += (sry*heli[ipart]  - siy*helr[ipart]);
      fz_tmp[ipart] += (srz*heli[ipart]  - siz*helr[ipart]); 
   }/*endfor*/

 } /* endif cluster BC */

/*======================================================================*/
/* VI) Get the nuclear hessian if necessary                             */
/*     Routine must be modified for idens_opt==1                        */

   if( (hess_calc == 3) && (idens_opt==0) ){
     get_atm_hess_recip(xk,yk,zk,hess_xx,hess_yy,hess_zz,
                        hess_xy,hess_xz,hess_yz,rhocr,rhoci,vtemp,
                        helr,heli,q,q_tmp,rvol,preg,sumr_h,sumi_h,
                        icount,hess_calc,iperd,
                        idens_opt,natm_use,ip_loc_cp_box);
   }/* endif */

   if( (hess_calc == 3) && (idens_opt==1) ){
#ifdef STILL_BROKEN
    printf("@@@@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@\n");
    printf("Hess calc not working with dual grid option, yet.\n");
    printf("How did I get here?\n");
    printf("Required changes are given in get_atm_hess_recip.\n");
    printf("Dawn, as always, has the required references.\n");
    printf("@@@@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@\n");
    exit(1);
#endif
     get_atm_hess_recip(xk,yk,zk,hess_xx,hess_yy,hess_zz,
                        hess_xy,hess_xz,hess_yz,rhocr,rhoci,vtemp,
                        helr,heli,q,q_tmp,rvol,preg,sumr_h,sumi_h,
                        icount,hess_calc,iperd,
                        idens_opt,natm_use,ip_loc_cp_box);
   }/* endif */

/*======================================================================*/
/* VII) If break point two, increment the helpful vectors                 */

   if(ibreak2[(icount+koff)]==1){
     for(ipart=1;ipart<=natm_use;ipart++){
       temp = helr[ipart];
       helr[ipart] = helr[ipart]*cossc[ipart] - heli[ipart]*sinsc[ipart];
       heli[ipart] = heli[ipart]*cossc[ipart] + temp*sinsc[ipart];
     }/*endfor*/
   }/*endif*/



 }/*endfor:icount loop over k vectors */



/*======================================================================*/
/* VIII) g=0 term (local pseudopotential) including term for CBCs       */

  if((myid_state+1)==np_states){
  
    if(ipseud_opt==1){
      ak2[(ngo+1)] = 0.0;

      for(ipart=1;ipart<=natm_use;ipart++){
        vtemp[ipart] = gzvps[iatm_typ[ipart]];
      }/*endfor*/
      vextr[(ngo+1)] =  dsum1(natm_use,vtemp,1)*rvol;
      vexti[(ngo+1)] = 0.0;
    }else{ /*large sparse grid */
      vextr[(ngo+1)] = 0.0;
      bgr  = dsum1(natm_use,q_tmp,1); 
      bgr  = bgr*M_PI/(alpha_conv_dual*alpha_conv_dual*vol);
      vextr[(ngo+1)] = bgr; 
    }/*endif*/

    if( (iperd!=3) && (idens_opt==0) ) {
      vextr[(ngo+1)] -= dsum1(natm_use,q_tmp,1)*clus_corr_r[(ngo+1)]*rvol;
    }/*endif*/

    *pseud_hess_loc = vextr[(ngo+1)];

    if( (iperd>0) && (iperd!=3) && (idens_opt==0)) {
      q_sum1 = dsum1(natm_use,q,1);
      if(iperd == 2){
         vrecip += 0.5*q_sum1*q_sum1*rvol*(clus_corr_r[(ngo+1)]
                 + M_PI/(alp_clus*alp_clus));
       } else {
         vrecip += 0.5*q_sum1*q_sum1*clus_corr_r[(ngo+1)]*rvol;
       }/* endif iperd */
    }/*endif*/

  }/*endif*/

/*======================================================================*/
/* IX) Collect the forces */

  if(idens_opt==0){
    for(ipart=1;ipart<=natm_use;ipart++){
      fx[ipart] += fx_tmp[ipart];
      fy[ipart] += fy_tmp[ipart];
      fz[ipart] += fz_tmp[ipart];
    }/*endfor*/
  }else{
    for(ipart=1;ipart<=natm_use;ipart++){
      fx[ip_loc_cp_box[ipart]] += fx_tmp[ipart];
      fy[ip_loc_cp_box[ipart]] += fy_tmp[ipart];
      fz[ip_loc_cp_box[ipart]] += fz_tmp[ipart];
    }/*endfor*/
  }/*endif*/


/*======================================================================*/
/* X) Add in the surface term the dumb way */

  if(iperd==2 && myid_state == 0){

    area = hmat[1]*hmat[5]-hmat[2]*hmat[4];
    for(ipart=1;ipart<=natm_use;ipart++){
      for(jpart=1;jpart<=natm_use;jpart++){
          dz = z[ipart]-z[jpart];
          dz -= hmat[9]*NINT(dz/hmat[9]);
          arg = alp_clus*fabs(dz-0.5*hmat[9]);
          vrecip += 0.5*q[ipart]*q[jpart]*surf_corr(arg)/(area*alp_clus);
      }/* endfor jpart */
    } /* endfor ipart */

    for(ipart=1;ipart<=natm_use;ipart++){
      for(jpart=1;jpart<=natm_use;jpart++){
        if(ipart != jpart){
          dz = z[ipart]-z[jpart];
          dz -= hmat[9]*NINT(dz/hmat[9]);
          argm = alp_clus*(dz - 0.5*hmat[9]);
          fargm = fabs(argm);
          argp = alp_clus*(dz + 0.5*hmat[9]);
          fargp = fabs(argp);
          fz[ipart] -= 0.5*q[ipart]*q[jpart]*alp_clus
                           *(dsurf_corr(fargm)*(argm/fargm)
                     +       dsurf_corr(fargp)*(argp/fargp))/(area*alp_clus);
        }/* endif */
      }/* endfor */
    }/* endfor */

    if(hess_calc == 3){
      for(ipart=1;ipart<=natm_use;ipart++){
        hess_ind = (ipart-1)*natm_use + ipart;
        for(jpart=1;jpart<=natm_use;jpart++){
            dz = z[ipart]-z[jpart];
            dz -= hmat[9]*NINT(dz/hmat[9]);
            argm = alp_clus*(dz - 0.5*hmat[9]);
            fargm = fabs(argm);
            argp = alp_clus*(dz + 0.5*hmat[9]);
            fargp = fabs(argp);
            hess_zz[hess_ind] += 0.5*q[ipart]*q[jpart]*alp_clus
                               *(d2surf_corr(fargm) + d2surf_corr(fargp))/area;
        }/* endfor jpart */
      }/* endfor ipart */
      for(ipart=1;ipart<=natm_use;ipart++){
        for(jpart=1;jpart<=natm_use;jpart++){
            hess_ind = (ipart-1)*natm_use + jpart;
            dz = z[ipart]-z[jpart];
            dz -= hmat[9]*NINT(dz/hmat[9]);
            argm = alp_clus*(dz - 0.5*hmat[9]);
            fargm = fabs(argm);
            argp = alp_clus*(dz + 0.5*hmat[9]);
            fargp = fabs(argp);
            hess_zz[hess_ind] -= 0.5*q[ipart]*q[jpart]*alp_clus
                               *(d2surf_corr(fargm) + d2surf_corr(fargp))/area;
        }/* endfor jpart */
      }/* endfor ipart */
    }/* endif hessian */

  }/* endif iperd */
  if(np_states > 1) Barrier(comm);

/*======================================================================*/
/* XI) Collect diagonal piece of the pressure tensor and assign vrecip  */

   if((cp_ptens==1) && (idens_opt==0) ) {
     pvten[1] = pvten[1] + vrecip;
     pvten[5] = pvten[5] + vrecip;
     pvten[9] = pvten[9] + vrecip;
   }/*endif*/ 

/*======================================================================*/
/* XII) Finally, store the final value of vrecip */

 if(idual_switch == 0){*vrecip_ret = vrecip;}


/*======================================================================*/
    }/*end routine*/
/*======================================================================*/


/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void control_ewd_loc_pme(CLATOMS_INFO *clatoms_info,CLATOMS_POS *clatoms_pos,
                     CELL *cell, PTENS *ptens, EWALD *ewald, CPEWALD *cpewald, 
                     CPSCR *cpscr, PSEUDO *pseudo, EWD_SCR *ewd_scr,  
                     CPOPTS *cpopts, ATOMMAPS *atommaps, double *vrecip_ret,
                     double *pseud_hess_loc,COMMUNICATE *communicate,
                     FOR_SCR *for_scr,int cp_dual_grid_opt,int idual_switch,
                     PARA_FFT_PKG3D *cp_para_fft_pkg3d_lg)

/*==========================================================================*/
/*         Begin Routine                                                    */
   {/*Begin Routine*/
/*==========================================================================*/
/*         Local Variable declarations                                      */
     
#include "../typ_defs/typ_mask.h"

  int istart,ngo,irem,idiv;
  int ipart,jpart,iii,itype,i;
  int icount,koff,natm_use;

  double falp2,vol,rvol,pivol,fpi,arg,q_sum1,bgr;
  double aka,akb,akc,xk,yk,zk,atemp,btemp,ctemp;
  double xtemp,ytemp,ztemp;

  double sumr,sumi,g2,g4,preg,prep,tpi,pi,g;
  double srx,sry,srz,six,siy,siz,temp,smag;
  double vrecip;
  double dx,dy,dz;
  double sx,sy,sz;
  double asx,asy,asz;
  double qgrid_dx,qgrid_dy,qgrid_dz;
  double qgrid_real_new,qgrid_imag_new;

/*--------------------------------------------*/
/*         Local Pointer declarations         */

  /*------------------*/
  /* Atom information */
  int natm_tot      = clatoms_info->natm_tot; 
  int nchrg         = clatoms_info->nchrg;
  int hess_calc     = clatoms_info->hess_calc;
  int *clatoms_ichrg  = clatoms_info->ichrg;

  double *x         = clatoms_pos->x;
  double *y         = clatoms_pos->y;
  double *z         = clatoms_pos->z;
  double *fx        = clatoms_pos->fx;     
  double *fy        = clatoms_pos->fy;
  double *fz        = clatoms_pos->fz;
  double *fx_tmp    = ewd_scr->fx2;     
  double *fy_tmp    = ewd_scr->fy2;     
  double *fz_tmp    = ewd_scr->fz2;     
  double *q         = clatoms_info->q;
  double *hess_xx   = clatoms_pos->hess_xx;
  double *hess_xy   = clatoms_pos->hess_xy;
  double *hess_xz   = clatoms_pos->hess_xz;
  double *hess_yy   = clatoms_pos->hess_yy;
  double *hess_yz   = clatoms_pos->hess_yz;
  double *hess_zz   = clatoms_pos->hess_zz;
  int natm_typ      = atommaps->natm_typ;
  int *index_atm    = for_scr->index_atm;
  int *iatm_typ;            /*Assigned below based on flags */
  int *iatm_typ_full;

  /*--------------------------------*/
  /* Cell and pressure information */
  int iperd                 = cell->iperd;
  int cp_ptens              = cpopts->cp_ptens_calc;
  double *pvten             = ptens->pvten_tmp;
  double *cp_box_center     = cell->cp_box_center;
  double *cp_box_center_rel = cell->cp_box_center_rel;
  double *hmat;             /* Assigned below based on flags */
  double *hmati;
  double *hmat_big          = cell->hmat;
  double *hmati_big         = cell->hmati;

  /*----------------------*/
  /* G-vector information */
  int *kastore;             /* Assigned below based on flags */
  int *kbstore;
  int *kcstore;
  int *ibreak1;
  int *ibreak2;
  double *vextr;
  double *vexti;
  double *dvextr;
  double *dvexti;
  double *rhocr;
  double *rhoci;
  double *ak2;
  int nktot;

  /*----------------------------------------------*/
  /* Pseudo-potential and Reduced Periodicity info*/
  int nsplin_g         = pseudo->nsplin_g;
  int n_rad_max        = pseudo->n_rad_max;
  double *clus_corr_r  = ewald->clus_corr_r;
  double *dclus_corr_r = ewald->dclus_corr_r;
  double alpha_conv_dual = pseudo->alpha_conv_dual;
  double dg_spl        = pseudo->dg_spl;
  double gmin_spl      = pseudo->gmin_spl;
  double *vps0         = pseudo->vps0;
  double *vps1         = pseudo->vps1;
  double *vps2         = pseudo->vps2;
  double *vps3         = pseudo->vps3;
  double *dvps0        = pseudo->dvps0;
  double *dvps1        = pseudo->dvps1;
  double *dvps2        = pseudo->dvps2;
  double *dvps3        = pseudo->dvps3;
  double *gzvps        = pseudo->gzvps;
  double *q_pseud      = pseudo->q_pseud;
  int *loc_opt         = pseudo->loc_opt;
  int np_loc_cp_box    = pseudo->np_loc_cp_box;
  int *ip_loc_cp_box   = pseudo->ip_loc_cp_box;

  /*---------------------------------*/
  /* Ewald and ewald scr information */
  double alp_ewald  = ewald->alp_ewd;
  double *cossc     = ewd_scr->cossc;  
  double *sinsc     = ewd_scr->sinsc;
  double *helr      = ewd_scr->helr;   
  double *heli      = ewd_scr->heli;
  double *vtemp     = ewd_scr->temp;
  double *dvtemp    = ewd_scr->vtemp_now;
  double *ewd_scr_x = ewd_scr->x;
  double *ewd_scr_y = ewd_scr->y;
  double *ewd_scr_z = ewd_scr->z;
  double *ewd_scr_q = ewd_scr->q;
  double *q_tmp     = ewd_scr->q;

  /*---------------------------------*/
  /* PME and PME scr information */
  int scale_opt_tmp;
  int ka,kb,kc;
  int kb_str,kb_end;
  int kb_off,kc_off,index_lg,index_lg_cmplx;
  int j,j1,j2,n,ia,ib,ic,ja,jb,jc,iatm,iatm1,iend,nnow;
  int ktemp,itemp;
  int nrecip_grid;
  int ngrid_ab;
  int ngrid_bc,nedge;
  double grid_a,grid_b,grid_c;
  double mn_a_tmp,mn_b_tmp,mn_c_tmp;
  double helr_pme,heli_pme;
  double preg_dual,falp2_dual;

  int pme_on          = cpscr->cpscr_atom_pme.pme_on;
  int n_interp        = cpscr->cpscr_atom_pme.n_interp;
  int nlen_pme        = cpscr->cpscr_atom_pme.nlen_pme;
   /*Num: equal to ewald->nktot            */
  int nktot_pme       = cpscr->cpscr_atom_pme.nktot_pme; 
   /*Num: PME mesh same as large sparse grid*/
  int nkf1            = cpscr->cpscr_atom_pme.nkf1;
  int nkf2            = cpscr->cpscr_atom_pme.nkf2;
  int nkf3            = cpscr->cpscr_atom_pme.nkf3;

  int ngrid_a         = cpscr->cpscr_atom_pme.ngrid_a;
  int ngrid_b         = cpscr->cpscr_atom_pme.ngrid_b;
  int ngrid_c         = cpscr->cpscr_atom_pme.ngrid_c;
   /*Lst: Lth: nlen_pme                     */
  int *iatemp         = cpscr->cpscr_atom_pme.iatemp;
  int *ibtemp         = cpscr->cpscr_atom_pme.ibtemp;
  int *ictemp         = cpscr->cpscr_atom_pme.ictemp;
   /*Lst: Lth nkf3                          */
  int *nc             = cpscr->cpscr_atom_pme.nc;
  int *ioff_c         = cpscr->cpscr_atom_pme.ioff_c;
   /*Lst: Lth: ninterp*nlen_pme             */
  int **igrid_a       = cpscr->cpscr_atom_pme.igrid_a;
  int **igrid_b       = cpscr->cpscr_atom_pme.igrid_b;
  int **igrid_c       = cpscr->cpscr_atom_pme.igrid_c;
   /*Lst: Lth: ninterp*nlen_pme             */
  int **igrid_now     = cpscr->cpscr_atom_pme.igrid_now;
  /*Lst: Lth:nlen_pme                       */
  double *frac_a      = cpscr->cpscr_atom_pme.frac_a;
  double *frac_b      = cpscr->cpscr_atom_pme.frac_b;
  double *frac_c      = cpscr->cpscr_atom_pme.frac_c;

  /*Lst: Lth: ninterp*nlen_pme              */
  double *aj          = cpscr->cpscr_atom_pme.aj;
  double *rn          = cpscr->cpscr_atom_pme.rn;
  double *rn1         = cpscr->cpscr_atom_pme.rn1;
  double **ua         = cpscr->cpscr_atom_pme.ua;
  double **ub         = cpscr->cpscr_atom_pme.ub;
  double **uc         = cpscr->cpscr_atom_pme.uc;
  double **mn_a       = cpscr->cpscr_atom_pme.mn_a;
  double **mn_b       = cpscr->cpscr_atom_pme.mn_b;
  double **mn_c       = cpscr->cpscr_atom_pme.mn_c;
  double **dmn_a      = cpscr->cpscr_atom_pme.dmn_a;
  double **dmn_b      = cpscr->cpscr_atom_pme.dmn_b;
  double **dmn_c      = cpscr->cpscr_atom_pme.dmn_c;
  double **qgrid_now  = cpscr->cpscr_atom_pme.qgrid_now;

  double *bw_r        = cpscr->cpscr_atom_pme.bw_r;
  double *bw_i        = cpscr->cpscr_atom_pme.bw_i;
  double *bweight_tot = cpscr->cpscr_atom_pme.bweight_tot;

  double *qgrid          = cpscr->cpscr_atom_pme.qgrid;
  double *qgrid_tmp      = cpscr->cpscr_atom_pme.qgrid_scr;
  double *qgrid_tmp_real = cpscr->cpscr_atom_pme.qgrid_tmp_real;
  double *qgrid_tmp_imag = cpscr->cpscr_atom_pme.qgrid_tmp_imag;

  /*---------------------------*/
  /* Communication information */
  int myid_state    = communicate->myid_state;
  int np_states     = communicate->np_states;
  MPI_Comm comm     = communicate->comm_states;

  /*PARALLEL FFT VARIABLES*/

 int *skb_fft_ka_proc_all  = cp_para_fft_pkg3d_lg->skb_fft_ka_proc_all;
 int *skc_fft_ka_proc_all  = cp_para_fft_pkg3d_lg->skc_fft_ka_proc_all;
 int *ekb_fft_ka_proc_all  = cp_para_fft_pkg3d_lg->ekb_fft_ka_proc_all;
 int *ekc_fft_ka_proc_all  = cp_para_fft_pkg3d_lg->ekc_fft_ka_proc_all;

 int   skc_fft_ka_proc_lg = cp_para_fft_pkg3d_lg->skc_fft_ka_proc;
 int   ekc_fft_ka_proc_lg = cp_para_fft_pkg3d_lg->ekc_fft_ka_proc;
 int   skb_fft_ka_proc_lg = cp_para_fft_pkg3d_lg->skb_fft_ka_proc;
 int   ekb_fft_ka_proc_lg = cp_para_fft_pkg3d_lg->ekb_fft_ka_proc;
 int   myid               = cp_para_fft_pkg3d_lg->myid;
 int   nproc              = cp_para_fft_pkg3d_lg->num_proc;
 int   iproc;
 int   myidp1;
/*======================================================================*/
/* 0) Assign local pointers                                             */

  /*THIS ROUTINE IS ONLY CALLED WHEN CP_DUAL_GRID_OPT == 2*/
  /* LARGE SPARSE GRID*/

    kastore   = ewald->kastr;
    kbstore   = ewald->kbstr;
    kcstore   = ewald->kcstr;
    ibreak1   = ewald->ibrk1;
    ibreak2   = ewald->ibrk2;
    vextr     = cpscr->cpscr_loc.vextr;
    vexti     = cpscr->cpscr_loc.vexti;
    dvextr    = cpscr->cpscr_loc.dvextr;
    dvexti    = cpscr->cpscr_loc.dvexti;
    rhocr     = cpscr->cpscr_rho.rhocr_up;
    rhoci     = cpscr->cpscr_rho.rhoci_up;
    ak2       = cpewald->ak2;
    nktot     = ewald->nktot;
    hmat      = cell->hmat;
    hmati     = cell->hmati;
    natm_use  = natm_tot;
    iatm_typ  = atommaps->iatm_atm_typ;

    myidp1    = myid + 1;
/*======================================================================*/
/*======================================================================*/
/* PME to obtain Structure Factor                                       */
/*======================================================================*/
/* I) Construct some useful constants                                   */

  pi    = M_PI;
  tpi   = 2.0*pi;
  fpi   = 4.0*pi;

  vol   = getdeth(hmat);
  rvol  = 1.0/vol;
  pivol = vol/4.0/pi;
  falp2 = 4.0*alp_ewald*alp_ewald;

   nrecip_grid    = 2*ngrid_a*ngrid_b*ngrid_c;
   grid_a         = (double) ngrid_a;
   grid_b         = (double) ngrid_b;
   grid_c         = (double) ngrid_c;
   ngrid_bc       = ngrid_b*ngrid_c;
   ngrid_ab       = ngrid_a*ngrid_b;
   for(j=1;j<=n_interp;j++){
     aj[j] = (double) (j-1);
     rn[j] = (double) (j);
     if(j > 1){rn1[j] = 1.0/((double)(j-1));}
   }/*endfor*/
     rn1[1] = 0.0;

/*======================================================================*/
/*II) Find charged atoms, get their scaled imaged particle coordinates  */

   for(iatm = 1;iatm <= nchrg;++iatm){
      ktemp = clatoms_ichrg[iatm];
      ewd_scr_x[iatm] = x[ktemp];
      ewd_scr_y[iatm] = y[ktemp];
      ewd_scr_z[iatm] = z[ktemp];
      ewd_scr_q[iatm] = q[ktemp];
   }/*endfor*/

   for(iatm = 1;iatm <= nchrg;++iatm){
     atemp = ewd_scr_x[iatm]*hmati[1]
           + ewd_scr_y[iatm]*hmati[4]
           + ewd_scr_z[iatm]*hmati[7];
     btemp = ewd_scr_x[iatm]*hmati[2]
           + ewd_scr_y[iatm]*hmati[5]
           + ewd_scr_z[iatm]*hmati[8];
     ctemp = ewd_scr_x[iatm]*hmati[3]
           + ewd_scr_y[iatm]*hmati[6]
           + ewd_scr_z[iatm]*hmati[9];
     atemp = atemp - NINT((atemp-0.5));
     btemp = btemp - NINT((btemp-0.5));
     ctemp = ctemp - NINT((ctemp-0.5));
     ewd_scr_x[iatm] = atemp*grid_a;
     ewd_scr_y[iatm] = btemp*grid_b;
     ewd_scr_z[iatm] = ctemp*grid_c;
   }/*endfor*/

/*==========================================================================*/
/* III) Calculate the Cardinal B spline interpolation functions of the      */
/*     charge weighted density on the real space grid                       */

   for(i=1;i<=nrecip_grid;i++){
     qgrid[i]=0.0;
   }/*endfor*/

  for(iatm=1;iatm<=nchrg;iatm+=nlen_pme){
     iatm1 = iatm-1;
     iend = MIN(nchrg,iatm1+nlen_pme);
     nnow = iend-iatm1;

/*--------------------------------------------------------------------------*/ 
/* A) Using current fraction, find the grid points on which M_n is non-zero */
     for(i=1;i<=nnow;i++){
      itemp     = i+iatm1;
      iatemp[i] = (int) (ewd_scr_x[itemp]);
      ibtemp[i] = (int) (ewd_scr_y[itemp]);
      ictemp[i] = (int) (ewd_scr_z[itemp]);
      frac_a[i] = ewd_scr_x[itemp] - (double) (iatemp[i]);
      frac_b[i] = ewd_scr_y[itemp] - (double) (ibtemp[i]);
      frac_c[i] = ewd_scr_z[itemp] - (double) (ictemp[i]);
     }/*endfor*/

     for(j=1;j<=n_interp;j++){
      for(i=1;i<=nnow;i++){
       ua[j][i]    = frac_a[i] + aj[j];
       ub[j][i]    = frac_b[i] + aj[j];
       uc[j][i]    = frac_c[i] + aj[j];
       j2       = j-2;
       ia       = iatemp[i] - j2;
       ib       = ibtemp[i] - j2;
       ic       = ictemp[i] - j2;
       ia       = (ia>0 ? ia:ngrid_a+ia);
       ib       = (ib>0 ? ib:ngrid_b+ib);
       ic       = (ic>0 ? ic:ngrid_c+ic);

       /*PARALLEL*/
       igrid_a[j][i] = ia;
       igrid_b[j][i] = ib;
       igrid_c[j][i] = ic;

      }/*endfor*/
     }/*endfor*/

/*--------------------------------------------------------------------------*/ 
/* B) Initialize M2 and get the Mn's using the recursion relation           */ 
/*    Note: M_n is defined on 0 to n. Since frac between 0 and 1, the       */
/*          calculation is performed in an order that takes advantage of it */

     for(i=1;i<=nnow;i++){
      mn_a[1][i] = 1.0 - fabs(ua[1][i]-1.0);
      mn_b[1][i] = 1.0 - fabs(ub[1][i]-1.0);
      mn_c[1][i] = 1.0 - fabs(uc[1][i]-1.0);
      mn_a[2][i] = 1.0 - fabs(ua[2][i]-1.0);
      mn_b[2][i] = 1.0 - fabs(ub[2][i]-1.0);
      mn_c[2][i] = 1.0 - fabs(uc[2][i]-1.0);
     }/*endfor*/
     for(j=3;j<=n_interp;j++){
      for(i=1;i<=nnow;i++){
       mn_a[j][i]   = 0.0;
       mn_b[j][i]   = 0.0;
       mn_c[j][i]   = 0.0;
      }/*endfor*/
     }/*endfor*/

     for(n=3;n<=n_interp;n++){
       for(j=n;j>=2;j--){
        j1 = j-1;
        for(i=1;i<=nnow;i++){
         mn_a_tmp  = (ua[j][i]*mn_a[j][i]+(rn[n]-ua[j][i])*mn_a[j1][i])*rn1[n];
         mn_b_tmp  = (ub[j][i]*mn_b[j][i]+(rn[n]-ub[j][i])*mn_b[j1][i])*rn1[n];
         mn_c_tmp  = (uc[j][i]*mn_c[j][i]+(rn[n]-uc[j][i])*mn_c[j1][i])*rn1[n];
         mn_a[j][i] = mn_a_tmp;
         mn_b[j][i] = mn_b_tmp;
         mn_c[j][i] = mn_c_tmp;
        }/*end for: i*/
       }/*end for: j*/
       for(i=1;i<=nnow;i++){
        mn_a[1][i] = ua[1][i]*mn_a[1][i]*rn1[n];
        mn_b[1][i] = ub[1][i]*mn_b[1][i]*rn1[n];
        mn_c[1][i] = uc[1][i]*mn_c[1][i]*rn1[n];
      }/*endfor */
     }/*end for: k*/

/*--------------------------------------------------------------------------*/ 
/* C) Stuff Mn's on the grid:Sum is ordered to remove recursive dependencies*/ 
/*                           in qgrid. igrid_now is unique for fixed i      */ 
    
/*Does this processor hold grid point (igrid_a[ja][i],igrid_b[jb][i],igrid_c[jc][i])*/

  for(i=1;i<=nnow;i++){
   for(jc=1;jc<=n_interp;jc++){
     kc = igrid_c[jc][i];

   /*Is kc in range for this processor ? */
    if( (kc >= skc_fft_ka_proc_all[myidp1]) && 
        (kc <= ekc_fft_ka_proc_all[myidp1])   ){

      kb_str = (kc == skc_fft_ka_proc_all[myidp1] 
                    ? skb_fft_ka_proc_all[myidp1] : 1);

      kb_end = (kc == ekc_fft_ka_proc_all[myidp1] 
                    ? ekb_fft_ka_proc_all[myidp1] : ngrid_b);

    for(jb=1;jb<=n_interp;jb++){
      kb = igrid_b[jb][i];
      if( (kb >= kb_str) && (kb <= kb_end)){

     /*Is kb in range for this processor ? */

      for(ja=1;ja<=n_interp;ja++){
        ka = igrid_a[ja][i];
 
         index_lg = (kc - skc_fft_ka_proc_all[myidp1] -1)*nkf1*nkf2
	          + (kb -1)*nkf1
                  + (nkf2 - skb_fft_ka_proc_all[myidp1]+1)*nkf1
                  + ka;
 
        itemp     = i+iatm1;
        index_lg_cmplx = 2*index_lg  - 1;

        igrid_now[ja][i] = index_lg_cmplx;
        qgrid_now[ja][i] = mn_a[ja][i]*mn_b[jb][i]*mn_c[jc][i]*ewd_scr_q[itemp];

        qgrid[index_lg_cmplx] += qgrid_now[ja][i];
      }/*endfor ja*/
    }/*endif kb in range for this processor */
    }/*endfor jb*/

    } /*endif kc in range for this processor */
   }/*end for jc*/
   }/*end for i */

   }/*end for: iatm*/ 
/*==========================================================================*/

/*==========================================================================*/
/* IV) Fourier Transform qgrid                                              */

   scale_opt_tmp = cp_para_fft_pkg3d_lg->scale_opt;
   cp_para_fft_pkg3d_lg->scale_opt = 0;

   para_fft_gen3d_bck_to_g(qgrid,qgrid_tmp,cp_para_fft_pkg3d_lg); 

   cp_para_fft_pkg3d_lg->scale_opt = scale_opt_tmp;

   sngl_upack_coef(qgrid_tmp_real,qgrid_tmp_imag,qgrid,cp_para_fft_pkg3d_lg);

/*==========================================================================*/
/*V) Perform the ewald sum/ cp-potential calculation                        */

  vrecip  = 0.0;

  idiv    = (nktot+1)/np_states;
  irem    = (nktot+1) % np_states;
  ngo     = (myid_state <  irem ? idiv+1 : idiv);
  istart  = (myid_state <= irem ? myid_state*(idiv+1)+1 : 
                                irem*(idiv+1)+1+(myid_state-irem)*idiv);
  koff    = istart-1;
  if(np_states==myid_state+1){ngo--;}

 for(icount=1;icount<=ngo;icount++){

/*======================================================================*/
/* A) Get the k vectors                                                 */

   aka = (double)(kastore[(icount+koff)]);
   akb = (double)(kbstore[(icount+koff)]);
   akc = (double)(kcstore[(icount+koff)]);
   xk = (aka*hmati[1]+akb*hmati[2]+akc*hmati[3])*tpi;
   yk = (aka*hmati[4]+akb*hmati[5]+akc*hmati[6])*tpi;
   zk = (aka*hmati[7]+akb*hmati[8]+akc*hmati[9])*tpi;
   g2 = xk*xk+yk*yk+zk*zk;
   g4 = g2*g2;
   g  = sqrt(g2);
   ak2[icount] = g2;

/*======================================================================*/
/* B) Get the external potential                                        */
/*               (interaction of electron with particles)               */

   falp2_dual = 4.0*alpha_conv_dual*alpha_conv_dual;
   preg_dual  = exp((-g2/falp2_dual))/(g2*pivol);

/*----------------------------------------------------------------------*/
/* C) Cluster boundary condition correction                             */

   if( iperd != 3 ){
       preg_dual += clus_corr_r[icount]*rvol;
   }/* endif cluster boundary conditions */
 
/*----------------------------------------------------------------------*/
/* D) Get the real and imag parts of the structure factor               */
/* Create helr_pme and heli_pme                                         */

   helr_pme = qgrid_tmp_real[icount]*bw_r[icount]
            - qgrid_tmp_imag[icount]*bw_i[icount];

   heli_pme = qgrid_tmp_real[icount]*bw_i[icount]
            + qgrid_tmp_imag[icount]*bw_r[icount];

   vextr[icount]  = -helr_pme*preg_dual;
   vexti[icount]  =  heli_pme*preg_dual;

/*======================================================================*/
/* E) Use the stucture factor to get the                                */
/*      classical potential energy                                      */

    smag       = (qgrid_tmp_real[icount]*qgrid_tmp_real[icount]
	       +  qgrid_tmp_imag[icount]*qgrid_tmp_imag[icount])
                 *bweight_tot[(icount+koff)];

    preg = exp((-g2/falp2))/(g2*pivol);

   if( iperd>0 ){vrecip  = vrecip + smag*preg; } 

   /* DO NOT CHANGE THE ORDER */
     qgrid_real_new =  qgrid_tmp_real[icount]*(preg*bweight_tot[(icount+koff)]);
     qgrid_imag_new =  qgrid_tmp_imag[icount]*(preg*bweight_tot[(icount+koff)]);

   qgrid_tmp_real[icount] = preg_dual*
      (-rhocr[icount]*bw_r[icount] + rhoci[icount]*bw_i[icount]);

   qgrid_tmp_imag[icount] = preg_dual*
      (rhocr[icount]*bw_i[icount] + rhoci[icount]*bw_r[icount]);

   qgrid_tmp_real[icount] += qgrid_real_new;
   qgrid_tmp_imag[icount] += qgrid_imag_new;

 }/*endfor:icount loop over k vectors */


/*======================================================================*/
/* F) g=0 term (local pseudopotential) including term for CBCs          */

  if((myid_state+1)==np_states){
    /*large sparse grid */
    vextr[(ngo+1)] = 0.0;
    bgr  = dsum1(natm_use,q_tmp,1); 
    bgr  = bgr*M_PI/(alpha_conv_dual*alpha_conv_dual*vol);
    vextr[(ngo+1)] = bgr; 

    if( iperd!=3 ) {
      vextr[(ngo+1)] -= dsum1(natm_use,q_tmp,1)*clus_corr_r[(ngo+1)]*rvol;
    }/*endif*/

    *pseud_hess_loc = vextr[(ngo+1)];

    if( (iperd>0) && (iperd!=3) ) {
      q_sum1 = dsum1(natm_use,q,1);
      vrecip += (0.5*q_sum1*q_sum1*clus_corr_r[(ngo+1)]*rvol);
    }/*endif*/
  }/*endif*/

   qgrid_tmp_real[(ngo+1)]  = 0.0;
   qgrid_tmp_imag[(ngo+1)]  = 0.0;

/*======================================================================*/
/* G) Assign vrecip -- large sparse grid                                */

   *vrecip_ret = vrecip;


/*======================================================================*/
/* VI) FORCES on atoms                                                  */
/*======================================================================*/
/* VI) Fourier Transform qgrid                                          */

   pme_sngl_pack_coef(qgrid_tmp_real,qgrid_tmp_imag,qgrid,cp_para_fft_pkg3d_lg);
   para_fft_gen3d_fwd_to_r(qgrid,qgrid_tmp,cp_para_fft_pkg3d_lg); 

/*======================================================================*/
/* BELOW FROM PME ROUTINE*/
/*==========================================================================*/
/* VIII) Calculate the force                                                */

   for(iatm = 1;iatm <= nchrg;++iatm){
      fx_tmp[iatm] = 0.0;
      fy_tmp[iatm] = 0.0;
      fz_tmp[iatm] = 0.0;
   }/*endfor*/
   for(iatm=1;iatm<=nchrg;iatm+=nlen_pme){
     iatm1 = iatm-1;
     iend = MIN(nchrg,iatm1+nlen_pme);
     nnow = iend-iatm1;
/*--------------------------------------------------------------------------*/ 
/* A) Using current fraction, find the grid points on which M_n is non-zero */
     for(i=1;i<=nnow;i++){
      itemp     = i+iatm1;
      iatemp[i] = (int) (ewd_scr_x[itemp]);
      ibtemp[i] = (int) (ewd_scr_y[itemp]);
      ictemp[i] = (int) (ewd_scr_z[itemp]);
      frac_a[i] = ewd_scr_x[itemp] - (double) (iatemp[i]);
      frac_b[i] = ewd_scr_y[itemp] - (double) (ibtemp[i]);
      frac_c[i] = ewd_scr_z[itemp] - (double) (ictemp[i]);
     }/*endfor*/
     for(j=1;j<=n_interp;j++){
      for(i=1;i<=nnow;i++){
       ua[j][i]    = frac_a[i] + aj[j];
       ub[j][i]    = frac_b[i] + aj[j];
       uc[j][i]    = frac_c[i] + aj[j];
       j2       = j-2;
       ia       = iatemp[i] - j2;
       ib       = ibtemp[i] - j2;
       ic       = ictemp[i] - j2;
       ia       = (ia>0 ? ia:ngrid_a+ia);
       ib       = (ib>0 ? ib:ngrid_b+ib);
       ic       = (ic>0 ? ic:ngrid_c+ic);

       /* PARALLEL*/
       igrid_a[j][i] = ia;
       igrid_b[j][i] = ib;
       igrid_c[j][i] = ic;

      }/*endfor*/
     }/*endfor*/
/*--------------------------------------------------------------------------*/ 
/* B) Initialize M2 and get the Mn's using the recursion relation           */ 
/*    Note: M_n is defined on 0 to n. Since frac between 0 and 1, the       */
/*          calculation is performed in an order that takes advantage of it */
     for(i=1;i<=nnow;i++){
      mn_a[1][i] = 1.0 - fabs(ua[1][i]-1.0);
      mn_b[1][i] = 1.0 - fabs(ub[1][i]-1.0);
      mn_c[1][i] = 1.0 - fabs(uc[1][i]-1.0);
      mn_a[2][i] = 1.0 - fabs(ua[2][i]-1.0);
      mn_b[2][i] = 1.0 - fabs(ub[2][i]-1.0);
      mn_c[2][i] = 1.0 - fabs(uc[2][i]-1.0);
     }/*endfor*/
     for(j=3;j<=n_interp;j++){
      for(i=1;i<=nnow;i++){
       mn_a[j][i]   = 0.0;
       mn_b[j][i]   = 0.0;
       mn_c[j][i]   = 0.0;
      }/*endfor*/
     }/*endfor*/

     for(n=3;n<=n_interp;n++){
       for(j=n;j>=2;j--){
        j1 = j-1;
        for(i=1;i<=nnow;i++){
         mn_a_tmp  = (ua[j][i]*mn_a[j][i]+(rn[n]-ua[j][i])*mn_a[j1][i])*rn1[n];
         mn_b_tmp  = (ub[j][i]*mn_b[j][i]+(rn[n]-ub[j][i])*mn_b[j1][i])*rn1[n];
         mn_c_tmp  = (uc[j][i]*mn_c[j][i]+(rn[n]-uc[j][i])*mn_c[j1][i])*rn1[n];
         mn_a[j][i] = mn_a_tmp;
         mn_b[j][i] = mn_b_tmp;
         mn_c[j][i] = mn_c_tmp;
        }/*end for: i*/
       }/*end for: j*/
       for(i=1;i<=nnow;i++){
        mn_a[1][i] = ua[1][i]*mn_a[1][i]*rn1[n];
        mn_b[1][i] = ub[1][i]*mn_b[1][i]*rn1[n];
        mn_c[1][i] = uc[1][i]*mn_c[1][i]*rn1[n];
       }/*endfor */
       if(n==(n_interp-1)){
        for(i=1;i<=nnow;i++){
         dmn_a[1][i] = mn_a[1][i];
         dmn_b[1][i] = mn_b[1][i];
         dmn_c[1][i] = mn_c[1][i];
        }/*endfor*/
        for(j=2;j<=n_interp;j++){    
         j1 = j-1;
         for(i=1;i<=nnow;i++){
          dmn_a[j][i] = mn_a[j][i] - mn_a[j1][i];
          dmn_b[j][i] = mn_b[j][i] - mn_b[j1][i];
          dmn_c[j][i] = mn_c[j][i] - mn_c[j1][i];
         }/*endfor*/
        }/*endfor*/
       }/*endif*/
     }/*end for: n*/
/*--------------------------------------------------------------------------*/ 
/* C) Calculate the force                                                   */ 
  for(i=1;i<=nnow;i++){

   for(jc=1;jc<=n_interp;jc++){
     kc = igrid_c[jc][i];

   /*Is kc in range for this processor ? */
    if( (kc >= skc_fft_ka_proc_all[myidp1]) && 
        (kc <= ekc_fft_ka_proc_all[myidp1])   ){

      kb_str = (kc == skc_fft_ka_proc_all[myidp1] 
                    ? skb_fft_ka_proc_all[myidp1] : 1);

      kb_end = (kc == ekc_fft_ka_proc_all[myidp1] 
                    ? ekb_fft_ka_proc_all[myidp1] : ngrid_b);

    for(jb=1;jb<=n_interp;jb++){

     kb = igrid_b[jb][i];

   /*Is kb in range for this processor ? */
     if( (kb >= kb_str) && (kb <= kb_end)){
        kc_off = kc - skc_fft_ka_proc_all[myidp1];
        kb_off = kb - skb_fft_ka_proc_all[myidp1];

     for(ja=1;ja<=n_interp;ja++){
        ka = igrid_a[ja][i];

        index_lg = (kc - skc_fft_ka_proc_all[myidp1] -1)*nkf1*nkf2
                 + (kb -1)*nkf1
                 + (nkf2 - skb_fft_ka_proc_all[myidp1]+1)*nkf1
                 + ka;

        index_lg_cmplx   = 2*index_lg - 1;

        igrid_now[ja][i] = index_lg_cmplx;
        qgrid_now[ja][i] = qgrid[igrid_now[ja][i]];

        itemp = i+iatm1;
        atemp = dmn_a[ja][i]* mn_b[jb][i]* mn_c[jc][i]*grid_a;
        btemp =  mn_a[ja][i]*dmn_b[jb][i]* mn_c[jc][i]*grid_b;
        ctemp =  mn_a[ja][i]* mn_b[jb][i]*dmn_c[jc][i]*grid_c;
        qgrid_dx = atemp*hmati[1]+btemp*hmati[2]+ctemp*hmati[3];
        qgrid_dy = atemp*hmati[4]+btemp*hmati[5]+ctemp*hmati[6];
        qgrid_dz = atemp*hmati[7]+btemp*hmati[8]+ctemp*hmati[9];
        qgrid_now[ja][i] *= ewd_scr_q[itemp];
        fx_tmp[itemp] -= (qgrid_dx*qgrid_now[ja][i]);
        fy_tmp[itemp] -= (qgrid_dy*qgrid_now[ja][i]);
        fz_tmp[itemp] -= (qgrid_dz*qgrid_now[ja][i]); 

     }/*endfor ja*/

     }/*endif kb in range on this processor*/
     }/*endfor jb*/

     }/*endif kc in range on this processor */
    }/*endfor jc*/

  }/*endfor i*/

 }/*end for: iatm*/


/*======================================================================*/
/* XI) Collect the forces */

    for(ipart=1;ipart<=natm_use;ipart++){
      fx[ipart] += fx_tmp[ipart];
      fy[ipart] += fy_tmp[ipart];
      fz[ipart] += fz_tmp[ipart];
    }/*endfor*/

/*======================================================================*/
    }/*end routine*/
/*======================================================================*/






/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void get_vpsnow(int *index_atm,int nsplin_g,
                double gmin_spl,double  dg_spl,double g,
                double *vps0,double *vps1,double *vps2,double *vps3,
                double *vtemp,int *iatm_typ,int natm_typ,int npart,int ist)

/*========================================================================*/
   {/*begin routine*/
/*========================================================================*/
/*        Local variable declarations            */

#include "../typ_defs/typ_mask.h"

  int ipart,itype,iii;
  int index_now[200];
  double vtemp_atyp[200];
  double h[200],h0;
  double partem1[200],partem2[200],partem3[200];
  double partem4[200];

/*==========================================================================*/
/* Loop over atom types to calculate pseudopotential                        */


  for(itype=1;itype<=natm_typ;itype++){
    iii = (g-gmin_spl)/dg_spl + 1;
    iii = MIN(iii,nsplin_g);
    iii = MAX(iii,1);
    h0  = (double)(iii-1)*dg_spl+gmin_spl;
    h[itype] = g-h0;
    index_now[itype] = index_atm[itype] + iii;
  }/*endfor*/


  for(itype=1;itype<=natm_typ;itype++){
    partem1[itype] = vps0[index_now[itype]];
    partem2[itype] = vps1[index_now[itype]];
    partem3[itype] = vps2[index_now[itype]];
    partem4[itype] = vps3[index_now[itype]];
  }/*endfor*/


  for(itype=1;itype<=natm_typ;itype++){
    vtemp_atyp[itype] = ((partem4[itype]*h[itype]+partem3[itype])
                        *h[itype]+partem2[itype])*h[itype] 
                        + partem1[itype];
  }/*endfor*/

  for(ipart=1;ipart<=npart;ipart++){
    vtemp[ipart] = vtemp_atyp[iatm_typ[ipart]];
  }/*endfor*/

/*--------------------------------------------------------------------------*/
  }/*end routine*/
/*==========================================================================*/



/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void ewald3d_selfbgr_cp(CLATOMS_INFO *clatoms_info,
                        EWALD *ewald, PTENS *ptens,double vol,
                        double *vself,double *vbgr,int iperd)

/*========================================================================*/
/*             Begin Routine                                              */
     {/*Begin Routine*/
/*========================================================================*/
/*             Local variable declarations                                */

  int step=1,iii;
  double self,bgr;

/* Define local pointers                                                  */

  int pi_beads            = clatoms_info->pi_beads;
  int natm_tot            = clatoms_info->natm_tot;
  double *clatoms_q       = clatoms_info->q;
  double *ptens_pvten_tot = ptens->pvten_tot;
  double *ptens_pvten     = ptens->pvten;
  double alp_ewd          = ewald->alp_ewd;
  double self_erf         = ewald->self_erf;

/*========================================================================*/
/* I) Self and background */

  self     = ddot1(natm_tot,clatoms_q,step,clatoms_q,step);
  self     = -self*(alp_ewd)*(self_erf)/sqrt(M_PI);
  (*vself) = self;

  bgr      = dsum1(natm_tot,clatoms_q,step);
  bgr      = -0.5*bgr*bgr*M_PI/(alp_ewd*alp_ewd*vol); 
  (*vbgr)  = bgr;

/*========================================================================*/
/* II) Pressure */

  if( iperd <= 3 ){

    ptens_pvten_tot[1] += bgr*pi_beads;
    ptens_pvten_tot[5] += bgr*pi_beads;
    ptens_pvten_tot[9] += bgr*pi_beads;

    ptens_pvten[1]     += bgr*pi_beads;
    ptens_pvten[5]     += bgr*pi_beads;
    ptens_pvten[9]     += bgr*pi_beads;

  }/*endif*/

/*------------------------------------------------------------------------*/
   }/*end routine */
/*========================================================================*/



/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void get_ak2_sm(CPEWALD *cpewald,CELL *cell)

/*========================================================================*/
/*             Begin Routine                                              */
   {/*Begin Routine*/
/*========================================================================*/
/*             Local variable declarations                                */

  int icount;
  double aka,akb,akc,xk,yk,zk,g2,tpi;

/*             Local pointers                                            */

  int nktot_sm       = cpewald->nktot_sm;
  double *hmati      = cell->hmati_cp;
  int *kastore_sm    = cpewald->kastr_sm;
  int *kbstore_sm    = cpewald->kbstr_sm;
  int *kcstore_sm    = cpewald->kcstr_sm;
  double *ak2_sm     = cpewald->ak2_sm;

/*======================================================================*/
/* I) Get the k vectors                                                 */

 tpi = 2.0*M_PI;

 for(icount=1;icount<=nktot_sm;icount++){
  
   aka = (double)(kastore_sm[icount]);
   akb = (double)(kbstore_sm[icount]);
   akc = (double)(kcstore_sm[icount]);
   xk = (aka*hmati[1]+akb*hmati[2]+akc*hmati[3])*tpi;
   yk = (aka*hmati[4]+akb*hmati[5]+akc*hmati[6])*tpi;
   zk = (aka*hmati[7]+akb*hmati[8]+akc*hmati[9])*tpi;
   g2 = xk*xk+yk*yk+zk*zk;
   ak2_sm[icount] = g2;

 }/*endfor*/

/*========================================================================*/
}/*end routine */
/*========================================================================*/



/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void reduce_cp_hess_stuff(double *hess_xx,double *hess_yy,double *hess_zz,
                          double *hess_xy,double *hess_xz,double *hess_yz,
                          int hess_size,int myid_state,
                          MPI_Comm comm_states)

/*========================================================================*/
/*             Begin Routine                                              */
{/*Begin Routine*/
/*========================================================================*/
/*             Local variable declarations                                */
#include "../typ_defs/typ_mask.h"

  int i;
  double *hess_scr;

/*========================================================================*/

  hess_scr = (double *) cmalloc(hess_size*sizeof(double))-1;

/*------------------------------------------------------------------------*/

  if(myid_state==0){
      for(i=1;i<=hess_size;i++){
        hess_scr[i] = 0.0;
      }/* endfor */
  }/* endif */
  Reduce(&(hess_xx[1]),&(hess_scr[1]),hess_size,MPI_DOUBLE,MPI_SUM,0,
         comm_states);
  if(myid_state==0){
    for(i=1;i<=hess_size;i++){
      hess_xx[i] = hess_scr[i];
    }/* endfor */
  }/* endif */

/*------------------------------------------------------------------------*/

  if(myid_state==0){
    for(i=1;i<=hess_size;i++){
      hess_scr[i] = 0.0;
    }/* endfor */
  }/* endif */
  Reduce(&(hess_xy[1]),&(hess_scr[1]),hess_size,MPI_DOUBLE,MPI_SUM,0,
         comm_states);
  if(myid_state==0){
      for(i=1;i<=hess_size;i++){
       hess_xy[i] = hess_scr[i];
      }/* endfor */
    }/* endif */

/*------------------------------------------------------------------------*/

    if(myid_state==0){
      for(i=1;i<=hess_size;i++){
        hess_scr[i] = 0.0;
      }/* endfor */
    }/* endif */
    Reduce(&(hess_xz[1]),&(hess_scr[1]),hess_size,MPI_DOUBLE,MPI_SUM,0,
           comm_states);
    if(myid_state==0){
      for(i=1;i<=hess_size;i++){
       hess_xz[i] = hess_scr[i];
      }/* endfor */
    }/* endif */

/*------------------------------------------------------------------------*/

    if(myid_state==0){
      for(i=1;i<=hess_size;i++){
        hess_scr[i] = 0.0;
      }/* endfor */
    }/* endif */
    Reduce(&(hess_yy[1]),&(hess_scr[1]),hess_size,MPI_DOUBLE,MPI_SUM,0,
           comm_states);
    if(myid_state==0){
      for(i=1;i<=hess_size;i++){
       hess_yy[i] = hess_scr[i];
      }/* endfor */
    }/* endif */

/*------------------------------------------------------------------------*/

    if(myid_state==0){
      for(i=1;i<=hess_size;i++){
        hess_scr[i] = 0.0;
      }/* endfor */
    }/* endif */
    Reduce(&(hess_yz[1]),&(hess_scr[1]),hess_size,MPI_DOUBLE,MPI_SUM,0,
           comm_states);
    if(myid_state==0){
      for(i=1;i<=hess_size;i++){
       hess_yz[i] = hess_scr[i];
      }/* endfor */
    }/* endif */

/*------------------------------------------------------------------------*/

    if(myid_state==0){
      for(i=1;i<=hess_size;i++){
        hess_scr[i] = 0.0;
      }/* endfor */
    }/* endif */
    Reduce(&(hess_zz[1]),&(hess_scr[1]),hess_size,MPI_DOUBLE,MPI_SUM,0,
           comm_states);
    if(myid_state==0){
      for(i=1;i<=hess_size;i++){
       hess_zz[i] = hess_scr[i];
      }/* endfor */
    }/* endif */

/*------------------------------------------------------------------------*/

   cfree(&(hess_scr[1]));

/*========================================================================*/
  }/*end routine */
/*========================================================================*/





/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void reduce_cp_atm_forc(int natm_tot,double *fx,double *fy,double *fz,
                        double *fx_tmp,double *fy_tmp,double *fz_tmp,
                        MPI_Comm comm_states, int myid_state)

/*========================================================================*/
/*             Begin Routine                                              */
  {/*Begin Routine*/
/*========================================================================*/
/*             Local variable declarations                                */
#include "../typ_defs/typ_mask.h"


 int i;

/*========================================================================*/

  if(myid_state==0){
    for(i=1;i<=natm_tot;i++){
      fx_tmp[i] = 0.0;
      fy_tmp[i] = 0.0;
      fz_tmp[i] = 0.0;
    }/* endfor */
  }/* endif */

  Reduce(&(fx[1]),&(fx_tmp[1]),natm_tot,MPI_DOUBLE,MPI_SUM,0,comm_states);
  Reduce(&(fy[1]),&(fy_tmp[1]),natm_tot,MPI_DOUBLE,MPI_SUM,0,comm_states);
  Reduce(&(fz[1]),&(fz_tmp[1]),natm_tot,MPI_DOUBLE,MPI_SUM,0,comm_states);

  if(myid_state==0){
    for(i=1;i<=natm_tot;i++){
      fx[i] = fx_tmp[i];
      fy[i] = fy_tmp[i];
      fz[i] = fz_tmp[i];
    }/* endfor */
  }/* endif */

/*========================================================================*/
  }/*end routine */
/*========================================================================*/



/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/


void get_vpslong(int natm_use,double *vtemp,double g2,
                 double *q_tmp,double alpha_conv_dual,double pivol)

/*========================================================================*/
/*             Begin Routine                                              */
  {/*Begin Routine*/
/*========================================================================*/
/*             Local variable declarations                                */

   double preg,falp2;
   int i;

/*==========================================================================*/

  falp2 = 4.0*alpha_conv_dual*alpha_conv_dual;
  preg  = 4.0*M_PI*exp(-g2/falp2)/(g2);

  for(i=1;i<=natm_use;i++){
    vtemp[i] = -q_tmp[i]*preg;
  }/*endfor*/

/*========================================================================*/
  }/*end routine */
/*========================================================================*/


/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/


void get_dvpslong(int natm_use,double *vtemp,double g2,
                  double *q_tmp,double alpha_conv_dual,double pivol)

/*========================================================================*/
/*             Begin Routine                                              */
  {/*Begin Routine*/
/*========================================================================*/
/*             Local variable declarations                                */

   double preg,prep,falp2;
   int i;

/*==========================================================================*/

  falp2 = 4.0*alpha_conv_dual*alpha_conv_dual;
  preg  = 4.0*M_PI*exp(-g2/falp2)/(g2);
  prep  = -2.0*preg*(g2/falp2+1.0)/g2;

  for(i=1;i<=natm_use;i++){
    vtemp[i] = -q_tmp[i]*prep;
  }/*endfor*/

/*========================================================================*/
  }/*end routine */
/*========================================================================*/




/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/


void get_atm_hess_recip(double xk,double yk,double zk,
                        double *hess_xx,double *hess_yy,double *hess_zz,
                        double *hess_xy,double *hess_xz,double *hess_yz,
                        double *rhocr,double *rhoci,
                        double *vtemp,double *helr,double *heli,double *q,
                        double *q_tmp,double rvol,double preg, 
                        double sumr,double sumi,int icount,
                        int hess_calc,int iperd, int idens_opt,
                        int natm_use, int *ip_loc_cp_box)

/*========================================================================*/
/*             Begin Routine                                              */
  {/*Begin Routine*/
/*========================================================================*/
/*             Local variable declarations                                */

  int ipart,jpart;
  int hess_ind;
  double overall;
  double srxx,srxy,srxz,sryy,sryz,srzz;
  double sixx,sixy,sixz,siyy,siyz,sizz;

/*========================================================================*/
/* This routine must be modified to add indirect addressing for           */
/* it to work with idens_opt==1 !! Also careful with q_tmp vs q           */

   if( (iperd ==0) ) {

     for(ipart=1;ipart<=natm_use;ipart++){

         srxx = xk*xk*(2.0*rhocr[icount]*vtemp[ipart]*rvol);
         srxy = xk*yk*(2.0*rhocr[icount]*vtemp[ipart]*rvol);
         srxz = xk*zk*(2.0*rhocr[icount]*vtemp[ipart]*rvol);
         sryy = yk*yk*(2.0*rhocr[icount]*vtemp[ipart]*rvol);
         sryz = yk*zk*(2.0*rhocr[icount]*vtemp[ipart]*rvol);
         srzz = zk*zk*(2.0*rhocr[icount]*vtemp[ipart]*rvol);
         sixx = xk*xk*(-2.0*rhoci[icount]*vtemp[ipart]*rvol);
         sixy = xk*yk*(-2.0*rhoci[icount]*vtemp[ipart]*rvol);
         sixz = xk*zk*(-2.0*rhoci[icount]*vtemp[ipart]*rvol);
         siyy = yk*yk*(-2.0*rhoci[icount]*vtemp[ipart]*rvol);
         siyz = yk*zk*(-2.0*rhoci[icount]*vtemp[ipart]*rvol);
         sizz = zk*zk*(-2.0*rhoci[icount]*vtemp[ipart]*rvol);

         hess_ind = (ipart-1)*natm_use + ipart;
         hess_xx[hess_ind] -= (srxx*helr[ipart] + sixx*heli[ipart]);
         hess_xy[hess_ind] -= (srxy*helr[ipart] + sixy*heli[ipart]);
         hess_xz[hess_ind] -= (srxz*helr[ipart] + sixz*heli[ipart]);
         hess_yy[hess_ind] -= (sryy*helr[ipart] + siyy*heli[ipart]);
         hess_yz[hess_ind] -= (sryz*helr[ipart] + siyz*heli[ipart]);
         hess_zz[hess_ind] -= (srzz*helr[ipart] + sizz*heli[ipart]);

     } /* endfor */

  } else {

     for(ipart=1;ipart<=natm_use;ipart++){

         srxx = xk*xk*(2.0*rhocr[icount]*vtemp[ipart]*rvol);
         srxy = xk*yk*(2.0*rhocr[icount]*vtemp[ipart]*rvol);
         srxz = xk*zk*(2.0*rhocr[icount]*vtemp[ipart]*rvol);
         sryy = yk*yk*(2.0*rhocr[icount]*vtemp[ipart]*rvol);
         sryz = yk*zk*(2.0*rhocr[icount]*vtemp[ipart]*rvol);
         srzz = zk*zk*(2.0*rhocr[icount]*vtemp[ipart]*rvol);
         sixx = xk*xk*(-2.0*rhoci[icount]*vtemp[ipart]*rvol);
         sixy = xk*yk*(-2.0*rhoci[icount]*vtemp[ipart]*rvol);
         sixz = xk*zk*(-2.0*rhoci[icount]*vtemp[ipart]*rvol);
         siyy = yk*yk*(-2.0*rhoci[icount]*vtemp[ipart]*rvol);
         siyz = yk*zk*(-2.0*rhoci[icount]*vtemp[ipart]*rvol);
         sizz = zk*zk*(-2.0*rhoci[icount]*vtemp[ipart]*rvol);

         hess_ind = (ipart-1)*natm_use + ipart;
         hess_xx[hess_ind] -= (srxx*helr[ipart] + sixx*heli[ipart]);
         hess_xy[hess_ind] -= (srxy*helr[ipart] + sixy*heli[ipart]);
         hess_xz[hess_ind] -= (srxz*helr[ipart] + sixz*heli[ipart]);
         hess_yy[hess_ind] -= (sryy*helr[ipart] + siyy*heli[ipart]);
         hess_yz[hess_ind] -= (sryz*helr[ipart] + siyz*heli[ipart]);
         hess_zz[hess_ind] -= (srzz*helr[ipart] + sizz*heli[ipart]);

     } /* endfor */

     for(ipart=1;ipart<=natm_use;ipart++){

         srxx = xk*xk*(2.0*preg*q[ipart]*sumr);
         srxy = xk*yk*(2.0*preg*q[ipart]*sumr);
         srxz = xk*zk*(2.0*preg*q[ipart]*sumr);
         sryy = yk*yk*(2.0*preg*q[ipart]*sumr);
         sryz = yk*zk*(2.0*preg*q[ipart]*sumr);
         srzz = zk*zk*(2.0*preg*q[ipart]*sumr);
         sixx = xk*xk*(2.0*preg*q[ipart]*sumi);
         sixy = xk*yk*(2.0*preg*q[ipart]*sumi);
         sixz = xk*zk*(2.0*preg*q[ipart]*sumi);
         siyy = yk*yk*(2.0*preg*q[ipart]*sumi);
         siyz = yk*zk*(2.0*preg*q[ipart]*sumi);
         sizz = zk*zk*(2.0*preg*q[ipart]*sumi);

         hess_ind = (ipart-1)*natm_use + ipart;
         hess_xx[hess_ind] -= (srxx*helr[ipart] + sixx*heli[ipart]);
         hess_xy[hess_ind] -= (srxy*helr[ipart] + sixy*heli[ipart]);
         hess_xz[hess_ind] -= (srxz*helr[ipart] + sixz*heli[ipart]);
         hess_yy[hess_ind] -= (sryy*helr[ipart] + siyy*heli[ipart]);
         hess_yz[hess_ind] -= (sryz*helr[ipart] + siyz*heli[ipart]);
         hess_zz[hess_ind] -= (srzz*helr[ipart] + sizz*heli[ipart]);

     } /* endfor */

     for(ipart=1;ipart<=natm_use;ipart++){
       for(jpart=1;jpart<=natm_use;jpart++){
          hess_ind = (ipart-1)*natm_use + jpart;
          overall = 2.0*q[ipart]*q[jpart]*preg
                  *(helr[ipart]*helr[jpart] + heli[ipart]*heli[jpart]);
          hess_xx[hess_ind] += overall*xk*xk;
          hess_xy[hess_ind] += overall*xk*yk;
          hess_xz[hess_ind] += overall*xk*zk;
          hess_yy[hess_ind] += overall*yk*yk;
          hess_yz[hess_ind] += overall*yk*zk;
          hess_zz[hess_ind] += overall*zk*zk;
       }/* endfor */
     }/* endfor */
    

  } /* endif cluster BC */

/*========================================================================*/
  }/*end routine */
/*========================================================================*/





