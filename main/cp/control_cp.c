/*==========================================================================*/
/*Cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
/*                                                                          */
/*                         PI_MD:                                           */
/*             The future of simulation technology                          */
/*             ------------------------------------                         */
/*                     Module: control_cp                                   */
/*                                                                          */
/* This subprogram performs MD on a classical+abinitio potential            */
/* energy surface (GGLSDA/GGLDA-PES)                                        */
/*                                                                          */
/*                                                                          */
/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/


#include "standard_include.h"
#include "../typ_defs/typedefs_gen.h"
#include "../typ_defs/typedefs_cp.h"
#include "../typ_defs/typedefs_class.h"
#include "../typ_defs/typedefs_bnd.h"
#include "../typ_defs/typedefs_stat.h"
#include "../proto_defs/proto_main_cp_local.h"
#include "../proto_defs/proto_main_local.h"
#include "../proto_defs/proto_math.h"
#include "../proto_defs/proto_integrate_cp_entry.h"
#include "../proto_defs/proto_integrate_cp_local.h"
#include "../proto_defs/proto_output_cp_entry.h"
#include "../proto_defs/proto_output_entry.h"
#include "../proto_defs/proto_vel_sampl_cp_entry.h"
#include "../proto_defs/proto_vel_sampl_class_entry.h"
#include "../proto_defs/proto_analysis_cp_entry.h"
#include "../proto_defs/proto_output_cp_entry.h"
#include "../proto_defs/proto_integrate_md_entry.h"
#include "../proto_defs/proto_integrate_md_local.h"
#include "../proto_defs/proto_intra_con_entry.h"
#include "../proto_defs/proto_energy_ctrl_cp_entry.h"
#include "../proto_defs/proto_energy_cpcon_entry.h"
#include "../proto_defs/proto_communicate_wrappers.h"

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void control_cp(CLASS *class,BONDED *bonded,GENERAL_DATA *general_data,CP *cp,
                ANALYSIS *analysis)

/*=======================================================================*/
/*            Begin subprogram:                                          */
  {/*begin routine*/
/*=======================================================================*/
/*            Local variable declarations                                */
#include "../typ_defs/typ_mask.h"

  int itime,iii,i;        
  int ireset,itol_ok;
  double tol_now;

  int ip_now        = 1;
  double fc_mag_up  = 1000.0;
  double fc_mag_dn  = 1000.0;
  double tol_coef   = cp->cpopts.tol_coef;

  int ntime         = general_data->timeinfo.ntime;

  int simopts_cp    = general_data->simopts.cp;

  int ensopts_nve   = general_data->ensopts.nve;
  int ensopts_nvt   = general_data->ensopts.nvt;
  int ensopts_npt_i = general_data->ensopts.npt_i;
  int ensopts_npt_f = general_data->ensopts.npt_f;

  int int_res_tra   = general_data->timeinfo.int_res_tra;
  int int_res_ter   = general_data->timeinfo.int_res_ter;

  int iwrite_screen = general_data->filenames.iwrite_screen;
  int iwrite_dump   = general_data->filenames.iwrite_dump;
  int iwrite_confp  = general_data->filenames.iwrite_confp;
  int iwrite_confv  = general_data->filenames.iwrite_confv;
  int iwrite_atm_for = general_data->filenames.iwrite_atm_for;
  int iwrite_inst   = general_data->filenames.iwrite_inst;
  int iwrite_kseigs = general_data->filenames.iwrite_kseigs;
  int cp_elf_calc_frq = cp->cpcoeffs_info.cp_elf_calc_frq;

  int ivel_smpl_on  = class->vel_samp_class.ivel_smpl_on;
  int nvx_smpl      = class->vel_samp_class.nvx_smpl;
  int nvnhc_smpl    = class->vel_samp_class.nvnhc_smpl;
  int ivel_scale_on = class->vel_samp_class.ivel_scale_on;
  int nvx_scale     = class->vel_samp_class.nvx_scale;
  int iconstrnt     = bonded->constrnt.iconstrnt;

  int ivelc_smpl_on = cp->vel_samp_cp.ivelc_smpl_on;
  int nvcnhc_smpl   = cp->vel_samp_cp.nvcnhc_smpl;
  int nvc_smpl      = cp->vel_samp_cp.nvc_smpl;

  int ivelc_scal_on     = cp->vel_samp_cp.ivelc_scal_on;
  int iauto_vc_scal_opt = cp->vel_samp_cp.iauto_vc_scal_opt;
  int nvc_scal          = cp->vel_samp_cp.nvc_scal;
 
  double c_div          = cp->vel_samp_cp.div_scal;
  double vc_scal_tol    = cp->vel_samp_cp.vc_scal_tol;
  double te_ext         = cp->cpopts.te_ext;
  int count_diag_srot_old = general_data->stat_avg.count_diag_srot;
  int cp_norb           = cp->cpopts.cp_norb;

  int num_proc         = class->communicate.np;
  int myid             = class->communicate.myid; 
  MPI_Comm comm_states = class->communicate.comm_states;
  MPI_Comm world       = class->communicate.world;

/*=======================================================================*/
/* 0) Preliminary MD stuff                                               */

  prelim_cp(class,bonded,general_data,cp);
  general_data->stat_avg.updates = 0.0;

/*======================================================================*/
/* I) Write to Screen                                                   */

  if(myid==0){
    PRINT_LINE_STAR;
    printf("Running CP-MD \n");
   PRINT_LINE_DASH;
  }/*endif*/

/*======================================================================*/
/* II) Loop over the specified number of time steps */

  for(itime = 1;itime<= ntime; itime++){

    if(num_proc>1){Barrier(world);}
    cputime(&(general_data->stat_avg.cpu1)); 
    (general_data->timeinfo.itime) = itime;
    (class->energy_ctrl.itime)     = itime;
  
 /*---------------------------------------------------------------------*/
 /*   1)Do NVE dynamics                                                 */
    if(ensopts_nve==1){
      if( (int_res_tra==0) && (int_res_ter==0) ){
        int_NVE_cp(class,bonded,general_data,cp);
      }else{
#ifdef DEVELOPMENT
        int_NVE_cp_res(class,bonded,general_data,cp);
#else
        if(myid==0){
          printf("\n@@@@@@@@@@@@@@@-ERROR-@@@@@@@@@@@@@@@@@@@\n");
          printf("No multiple time step CP yet");
          printf("@@@@@@@@@@@@@@@@-ERROR-@@@@@@@@@@@@@@@@@@@\n");
          fflush(stdout);
        }/*endif*/
        exit(1);
#endif
      }/*endif*/
    }/*endif*/
  
 /*----------------------------------------------------------------------*/
 /*   2)Do NVT dynamics:                                                 */
    if(ensopts_nvt==1){
      if( (int_res_tra==0) && (int_res_ter==0) ){
        int_NVT_cp(class,bonded,general_data,cp);
      }else{
#ifdef DEVELOPMENT
        int_NVT_cp_res(class,bonded,general_data,cp); 
#else
        if(myid==0){
          printf("\n@@@@@@@@@@@@@@@-ERROR-@@@@@@@@@@@@@@@@@@@\n");
          printf("No multiple time step CP yet");
          printf("@@@@@@@@@@@@@@@@-ERROR-@@@@@@@@@@@@@@@@@@@\n");
          fflush(stdout);
        }/*endif*/
        exit(1);
#endif
      }/*endif*/
    }/*endif*/
  
 /*---------------------------------------------------------------------*/ 
 /*   3)Do isotropic npt dynamics:                                      */
    if(ensopts_npt_i==1){
      if( (int_res_ter==0) && (int_res_tra==0) ){
        int_NPTI_cp(class,bonded,general_data,cp);
      }else{
#ifdef DEVELOPMENT
        int_NPTI_cp_res(class,bonded,general_data,cp);
#else
        if(myid==0){
          printf("\n@@@@@@@@@@@@@@@-ERROR-@@@@@@@@@@@@@@@@@@@\n");
          printf("No multiple time step CP yet");
          printf("@@@@@@@@@@@@@@@@-ERROR-@@@@@@@@@@@@@@@@@@@\n");
          fflush(stdout);
        }/*endif*/
        exit(1);
#endif
      }/*endif*/
    }/*endif*/
  
 /*---------------------------------------------------------------------*/
 /*   4)Do flexible npt dynamics:                                       */
    if(ensopts_npt_f==1){
      if( (int_res_ter==0) && (int_res_tra==0) ){
        int_NPTF_cp(class,bonded,general_data,cp);
      }else{
#ifdef DEVELOPMENT
        int_NPTF_cp_res(class,bonded,general_data,cp);
#else
        if(myid==0){
          printf("\n@@@@@@@@@@@@@@@-ERROR-@@@@@@@@@@@@@@@@@@@\n");
          printf("No multiple time step CP yet");
          printf("@@@@@@@@@@@@@@@@-ERROR-@@@@@@@@@@@@@@@@@@@\n");
          fflush(stdout);
        }/*endif*/
        exit(1);
#endif
      }/*endif*/
    }/*endif*/
     
 /*---------------------------------------------------------------------*/
 /*   5)Rotate the orbitals:                                            */

    cp_rotate_control(cp,itime,ntime,&(general_data->stat_avg));
 
    /* If a rotation of the non-orthogonal orbitals has occurred because
       the off diagonal elements are too big, then
       1) if the isokinetic constraint is used, rescale the coefficient 
          velocities to match the fixed ke
       2) Call rattle again so that the coefficient velocities satisfy
          the norm only constraint */
 
    if ((cp_norb == 2) && 
       (count_diag_srot_old != general_data->stat_avg.count_diag_srot)){ 

      if (cp->cpopts.cp_isok_opt == 1) control_vc_scale(cp);
 
      /*(class->energy_ctrl.iget_full_inter)= 1;
      (class->energy_ctrl.iget_res_inter) = 0;
      (class->energy_ctrl.iget_full_intra)= 1;
      (class->energy_ctrl.iget_res_intra) = 0;

      cp_energy_control(class,bonded,general_data,cp);*/

      if(cp->cpopts.cp_gauss == 1){
        add_gauss_force(&(cp->cpcoeffs_info),&(cp->cpcoeffs_pos[1]),
                        &(cp->cpscr.cpscr_ovmat),&(cp->cpopts));
      }/*endif*/
   }
   count_diag_srot_old = general_data->stat_avg.count_diag_srot;
 /*----------------------------------------------------------------------*/
 /*  6) Check force tolerance                                            */

    check_coef_grad_mag(cp,&(general_data->simopts),
                        &fc_mag_up,&fc_mag_dn,&ireset,&itol_ok,
                        tol_coef,1,1,&(general_data->stat_avg));
    if(itol_ok != 1){
      if(myid==0){
        printf("\n$$$$$$$$$$$$$$$-WARNING-$$$$$$$$$$$$$$$$$$$\n");
        printf("Tolerance on cp forces exceeded %g\n",tol_coef);
        printf("Please, please reconsider your choice of parameters.\n");
        printf("$$$$$$$$$$$$$$$$-WARNING-$$$$$$$$$$$$$$$$$$$\n");
        fflush(stdout);
      }/*endif*/
    }/* endif */
     
/*----------------------------------------------------------------------*/
/*   7)Calculate some simple averages                                   */

    cputime(&(general_data->stat_avg.cpu2)); 
    (general_data->stat_avg.cpu_now)=(general_data->stat_avg.cpu2)-
                                     (general_data->stat_avg.cpu1);
    (general_data->stat_avg.acpu) += (general_data->stat_avg.cpu2)-
                                     (general_data->stat_avg.cpu1);

    simpavg_cp(&(general_data->timeinfo),&(general_data->stat_avg),
              &(general_data->cell),&(bonded->constrnt),
              &(general_data->ensopts),&(general_data->simopts),
              &(general_data->ptens),cp,&(class->communicate),
              &(class->nbr_list.verlist),&(class->energy_ctrl));
     
 /*-----------------------------------------------------------------------*/
 /*   8)Produce the output specified by the user                          */
 /*     Also, check for exit condition                                    */

    if(myid == 0) check_auto_exit(&(general_data->timeinfo.exit_flag));
    if(num_proc > 1) Bcast(&(general_data->timeinfo.exit_flag),1,MPI_INT,0,world);


    if(  (itime % iwrite_screen) ==0 ||
         (itime % iwrite_dump  ) ==0 ||
         (itime % iwrite_confp ) ==0 ||
         (itime % iwrite_confv ) ==0 ||
         (itime % iwrite_atm_for) ==0 ||  
         (itime % iwrite_kseigs) ==0 ||
         ((cp_elf_calc_frq > 0) && (itime % cp_elf_calc_frq) ==0) || 
         (itime % iwrite_inst   )   ==0  ||
         (general_data->timeinfo.exit_flag == 1) ){
         (general_data->filenames.ifile_open) = 0;
          output_cp(class,general_data,bonded,cp);
    }/*endif*/
  
 /*---------------------------------------------------------------------*/
 /*   9)Resample the particle velocities if desired                     */

    if( (ivel_smpl_on==1) && (nvx_smpl!=0) && (simopts_cp==1)){
      if( (itime % nvx_smpl)==0 ){
        control_vx_smpl(class,bonded,&(general_data->ensopts),
                     &(general_data->simopts),
                     general_data,general_data->error_check_on);
        if(iconstrnt==1){
          init_constraint(bonded,&(general_data->ptens));
        }/*endif:init constraints*/
      }/*endif*/
    }/*endif*/
     
 /*---------------------------------------------------------------------*/
 /*  10) Scale the atom velocities */

    if( (ivel_scale_on==1) && (nvx_scale!=0) && (simopts_cp==1) ){
      if( (itime % nvx_scale)==0){
        control_vx_scale(class,&(general_data->simopts),class->clatoms_info.text_atm);
      }/* endif*/
    }/*endif*/
     
 /*---------------------------------------------------------------------*/
 /*  11) Resample the extended class velocities if desired             */

    if( (ivel_smpl_on==1) && (nvnhc_smpl!=0) && (simopts_cp==1)){
      if( (itime % nvnhc_smpl)==0){
         control_vnhc_smpl(class,general_data);
      }/*endif*/
    }/*endif*/

 /*---------------------------------------------------------------------*/
 /*   12) Resample the coefficient velocities if desired                */

    if( (ivelc_smpl_on==1) && (nvc_smpl!=0) ){
      if( (itime % nvc_smpl)==0){
        control_vc_smpl(general_data,cp);
      }/* endif */
    }/* endif */

 /*---------------------------------------------------------------------*/
 /*   13) Resample the extended class of coefficients if desired        */

    if( (ivelc_smpl_on==1) && (nvcnhc_smpl!=0) ){
      if( (itime % nvcnhc_smpl )==0 ){
        control_vcnhc_smpl(cp);
      }/*endif*/
    }/*endif*/

 /*---------------------------------------------------------------------*/
 /*   14) Rescale the coefficient velocities if desired                */

    if( (ivelc_scal_on==1) && (nvc_scal!=0) ){
      if( (itime % nvc_scal)==0){
        control_vc_scale(cp);
      }/* endif */
    }/* endif */

 /*---------------------------------------------------------------------*/
 /*   15) Rescale the coefficient velocities if they get too hot        */

    if(iauto_vc_scal_opt==1) {
      tol_now = (general_data->stat_avg.kinet_cp)*2.0*BOLTZ/(c_div*te_ext);
      if(tol_now > vc_scal_tol){
         control_vc_scale(cp);
      }/* endif */
    }/* endif */

 /*---------------------------------------------------------------------*/
 /*  16) Analysis Routine                                               */

    analysis_cp(class,bonded,general_data,cp,analysis); 

 /*---------------------------------------------------------------------*/
 /*   17) Check for exit condition                                      */

    if(general_data->timeinfo.exit_flag == 1) itime = ntime;

 /*---------------------------------------------------------------------*/

  }/*endfor:itime */

/*======================================================================*/
/*  III)Write to Screen                                                 */

  if(myid==0){
    PRINT_LINE_DASH;
    printf("Completed CP-MD run \n");
    PRINT_LINE_STAR;
  }/*endif*/

/*-----------------------------------------------------------------------*/
   }/*end routine*/
/*==========================================================================*/





/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void prelim_cp(CLASS *class,BONDED *bonded,GENERAL_DATA *general_data,CP *cp)

/*=======================================================================*/
/*            Begin subprogram:                                          */
{/* begin function */
/*=======================================================================*/
/*            Local variable declarations                                */
#include "../typ_defs/typ_mask.h"

  double kinet_cp_temp     = 0.0;
  double kinet_nhc_cp_temp = 0.0;
  double cp_ehart_temp     = 0.0;
  double cp_exc_temp       = 0.0;
  double cp_muxc_temp      = 0.0;
  double cp_eext_temp      = 0.0;
  double cp_enl_temp       = 0.0;
  double cp_eke_temp       = 0.0;
  double kinet_cp_up_tmp,kinet_cp_dn_tmp;
  int iflag,iii,is,jjj;        
  int i,j,k,m;
  int iflag_mass = 1;
  double pvten_temp[10];

  int nstate_up        = cp->cpcoeffs_info.nstate_up;
  int myid             = class->communicate.myid; 
  int myid_state       = class->communicate.myid_state; 
  MPI_Comm comm_states = class->communicate.comm_states;
  int num_proc         = class->communicate.np;
  int simopts_cp;

  /*DEBUG */ MPI_Comm world = class->communicate.world;
  /*DY*/
  MPI_Comm comm_forc     =  class->communicate.comm_forc;
  int np_forc            = class->communicate.np_forc;
  double kinet_temp      = 0.0;
  double kinet_nhc_temp  = 0.0;
  double vintert_temp    = 0.0;

  /*DY*/

  simopts_cp = general_data->simopts.cp;

/*=======================================================================*/
/*   I) Write to Screen                                                  */

 if(myid==0){
  PRINT_LINE_STAR;
  printf("Performing preliminary tasks for CP-MD \n");
  PRINT_LINE_DASH;printf("\n");
 }/*endif*/

/*=======================================================================*/
/*  II) Initialize In-output                                             */

  general_data->stat_avg.updates   = 0.0;
  general_data->stat_avg.updates_w = 0.0;
  general_data->stat_avg.acpu      = 0.0;
  general_data->stat_avg.count_diag_srot      = 0.0;

/*=======================================================================*/
/* III)Get Energy and particle Forces                                    */

  if(bonded->constrnt.iconstrnt==1&&simopts_cp==1){
     init_constraint(bonded,&(general_data->ptens));
  }/*endif:init constraints, must be done before getting the energy*/

  zero_constrt_iters(&(general_data->stat_avg));

  (general_data->timeinfo.itime) = 0;
  (class->energy_ctrl.itime)     = 0;
  class->energy_ctrl.iget_full_inter = 1;
  class->energy_ctrl.iget_res_inter  = 0;
  if(general_data->timeinfo.int_res_ter==1){
    class->energy_ctrl.iget_res_inter=1;
  }
  class->energy_ctrl.iget_full_intra = 1;
  class->energy_ctrl.iget_res_intra  = 0;
  if((general_data->timeinfo.int_res_tra)==1){
    class->energy_ctrl.iget_res_intra=1;
  }

   cp_energy_control(class,bonded,general_data,cp);

   if(cp->cpopts.cp_gauss == 1){
    add_gauss_force(&(cp->cpcoeffs_info),&(cp->cpcoeffs_pos[1]),
                    &(cp->cpscr.cpscr_ovmat),&(cp->cpopts));
   }

  if(simopts_cp==1){
   if((myid==0 && np_forc == 1) || (np_forc > 1) ){
    get_tvten(&(class->clatoms_info),&(class->clatoms_pos[1]),
            &(general_data->stat_avg),&(general_data->ptens),
            &(general_data->cell));

   }

  }/*endif*/

  get_cpke(&(cp->cpcoeffs_info),&(cp->cpcoeffs_pos[1]),
           &(general_data->stat_avg),cp->cpopts.cp_lsda,
           cp->communicate.np_states);


/*=======================================================================*/
/* IV)Get Atm NHC Forces                                                     */

 if(simopts_cp==1 && ( (myid==0 && np_forc == 1) || (np_forc > 1) ) ){

  if((general_data->ensopts.nvt)==1){
     init_NHC_par(&(class->clatoms_info),&(class->clatoms_pos[1]),
              &(class->therm_info_class),&(class->therm_class),
              &(class->int_scr),iflag_mass,&(class->class_comm_forc_pkg));
    iflag = 0;
    nhc_vol_potkin(&(class->therm_info_class),&(class->therm_class),
                   &(general_data->baro),&(general_data->par_rahman),
		   &(general_data->stat_avg), 
                   &(general_data->statepoint),iflag,
                     class->class_comm_forc_pkg.myid);
  }/*endif*/
  if((general_data->ensopts.npt_i)==1){
     init_NHCPI_par(&(class->clatoms_info),&(class->clatoms_pos[1]),
                &(class->therm_info_class),&(class->therm_class),
                &(general_data->baro),
                &(class->int_scr),&(class->class_comm_forc_pkg));
    iflag = 1;
    nhc_vol_potkin(&(class->therm_info_class),&(class->therm_class),
                   &(general_data->baro),&(general_data->par_rahman),
		   &(general_data->stat_avg), 
                   &(general_data->statepoint),iflag,
                     class->class_comm_forc_pkg.myid);
  }/*endif*/
  if((general_data->ensopts.npt_f)==1){
     init_NHCPF_par(&(class->clatoms_info),&(class->clatoms_pos[1]),
                &(class->therm_info_class),&(class->therm_class),
                &(general_data->baro),
                &(general_data->par_rahman),&(general_data->cell),
                                            &(class->int_scr),
                &(class->class_comm_forc_pkg));
    iflag = 2;
    nhc_vol_potkin(&(class->therm_info_class),&(class->therm_class),
                   &(general_data->baro),&(general_data->par_rahman),
		   &(general_data->stat_avg), 
                   &(general_data->statepoint),iflag,
                     class->class_comm_forc_pkg.myid);
  }/*endif*/

    if((general_data->ensopts.npt_i+general_data->ensopts.npt_f)==0){
       general_data->stat_avg.kinet_v       = 0.0;
    }/*endif*/
    if((general_data->ensopts.npt_i+general_data->ensopts.npt_f
       +general_data->ensopts.nvt)==0){
       general_data->stat_avg.kinet_nhc     = 0.0;
    }/*endif*/
  }/*endif: full cp on and myid==0*/

/*=======================================================================*/
/* V)Get CP NHC Forces                                                   */

  if(cp->cptherm_info.num_c_nhc > 0){
   if(cp->cptherm_info.massiv_flag==0){
    init_cp_NHC(&(cp->cptherm_info),&(cp->cptherm_pos[1]),
                &(cp->cpcoeffs_info),&(cp->cpcoeffs_pos[1]),
                &(cp->cpscr),cp->cpopts.cp_lsda,comm_states,
                cp->communicate.np_states,cp->communicate.myid_state);
 
    nhc_cp_potkin(&(cp->cptherm_info),&(cp->cptherm_pos[1]),
                  &(general_data->stat_avg),myid_state,comm_states);
   }else{
    nhc_cp_potkin_massiv(&(cp->cptherm_info),&(cp->cptherm_pos[1]),
                         &(general_data->stat_avg));
   }/* endif */
  }/* endif */

/*========================================================================*/
/* VI) Initialize free energy stuff                                       */

  if((simopts_cp==1)&&(myid==0)){

     if(bonded->bond_free.num != 0){
       for(i=1;i<= (bonded->bond_free.nhist);i++){
         bonded->bond_free.hist[i] = 0.0;
       }
     }/*endif*/

     if(bonded->bend_free.num != 0){
       for(i=1;i<= (bonded->bend_free.nhist);i++){
         bonded->bend_free.hist[i] = 0.0;
       }
     }/*endif*/

     if(bonded->tors_free.num == 1){
       for(i=1;i<= bonded->tors_free.nhist;i++){
         bonded->tors_free.hist[i] = 0.0;
       }
     }/*endif*/
     if(bonded->tors_free.num == 2){
      for(i=1;i<= bonded->tors_free.nhist;i++){
      for(j=1;j<= bonded->tors_free.nhist;j++){
        bonded->tors_free.hist_2d[i][j] = 0.0;
      }/*endfor*/
      }/*endfor*/
     }/*endif*/


     if(bonded->rbar_sig_free.nfree != 0){
      for(i=1;i<=bonded->rbar_sig_free.nhist_sig;i++){
       for(j=1;j<=bonded->rbar_sig_free.nhist_bar;j++){
         bonded->rbar_sig_free.hist[j][i] = 0.0;
       }/*endfor*/
      }/*endfor*/
      for(m=1;m<=bonded->rbar_sig_free.nfree;m++){
        for(k=1;k<=bonded->rbar_sig_free.nhist_bar;k++){
          bonded->rbar_sig_free.hist_rn[m][k] = 0.0;
        }/*endfor*/
      }/*endfor*/
     }/*endif*/

  }/*endif*/

/*=======================================================================*/
/* VII) Communication for initial output                                    */

if(num_proc>1){

   Reduce(&(general_data->stat_avg.cp_ehart), &cp_ehart_temp,1,MPI_DOUBLE,
                   MPI_SUM,0,comm_states);
   Reduce(&(general_data->stat_avg.cp_exc), &cp_exc_temp,1,MPI_DOUBLE,
                   MPI_SUM,0,comm_states);
   Reduce(&(general_data->stat_avg.cp_muxc), &cp_muxc_temp,1,MPI_DOUBLE,
                   MPI_SUM,0,comm_states);
   Reduce(&(general_data->stat_avg.cp_eext), &cp_eext_temp,1,MPI_DOUBLE,
                   MPI_SUM,0,comm_states);
   Reduce(&(general_data->stat_avg.kinet_cp), &kinet_cp_temp,1,
                   MPI_DOUBLE,MPI_SUM,0,comm_states);
   Reduce(&(general_data->stat_avg.cp_enl), &cp_enl_temp,1,
                   MPI_DOUBLE,MPI_SUM,0,comm_states);
   Reduce(&(general_data->stat_avg.cp_eke), &cp_eke_temp,1,
                   MPI_DOUBLE,MPI_SUM,0,comm_states);
   Reduce(&(general_data->stat_avg.kinet_nhc_cp), &kinet_nhc_cp_temp,1,
                   MPI_DOUBLE,MPI_SUM,0,comm_states);
   for(i=1;i<=9;i++){pvten_temp[i]=0.0;}
   Reduce(&(general_data->ptens.pvten_tot[1]), &(pvten_temp[1]),9,MPI_DOUBLE,
             MPI_SUM,0,comm_states);

   general_data->stat_avg.kinet_cp = kinet_cp_temp;
   general_data->stat_avg.kinet_nhc_cp = kinet_nhc_cp_temp;
   general_data->stat_avg.cp_enl   = cp_enl_temp;
   general_data->stat_avg.cp_eke   = cp_eke_temp;
   general_data->stat_avg.cp_ehart = cp_ehart_temp;
   general_data->stat_avg.cp_exc   = cp_exc_temp;
   general_data->stat_avg.cp_muxc  = cp_muxc_temp;
   general_data->stat_avg.cp_eext  = cp_eext_temp;

   for(i=1;i<=9;i++){general_data->ptens.pvten_tot[i] =  pvten_temp[i];}
  
 }/*endif*/

  if(np_forc > 1){
     simpavg_md_communicate(&(general_data->stat_avg),&(general_data->ensopts),
             &(general_data->ptens),0,&(class->communicate),
             &(class->nbr_list.verlist),general_data->timeinfo.int_res_ter,1,1);

  }/*endif np_forc*/


/*=======================================================================*/
/* VIII) If doing CP isokinetic dynamics, set the fixed value of the
       fictitious kinetic energy using the initial value                */

  if(cp->cpopts.cp_isok_opt == 1){
    if(cp->communicate.np_states > 1){

      kinet_cp_up_tmp = general_data->stat_avg.kinet_cp_up;
      Allreduce(&kinet_cp_up_tmp,&(general_data->stat_avg.kinet_cp_up),
                1,MPI_DOUBLE,MPI_SUM,0,comm_states);

      if(cp->cpopts.cp_lsda == 1) {
        kinet_cp_dn_tmp = general_data->stat_avg.kinet_cp_dn;
        Allreduce(&kinet_cp_dn_tmp,&(general_data->stat_avg.kinet_cp_dn),
                  1,MPI_DOUBLE,MPI_SUM,0,comm_states);

      }
    }/* endif parallel */

    cp->cpcoeffs_info.k0_up = general_data->stat_avg.kinet_cp_up;

    if(cp->cpopts.cp_lsda == 1) cp->cpcoeffs_info.k0_dn = general_data->stat_avg.kinet_cp_dn;
  }/* endif isok_opt on */

/*=======================================================================*/
/*  IX) Write Energy to screen                                           */

  general_data->filenames.ifile_open = 1;
  general_data->timeinfo.itime = 0;

  output_cp(class,general_data,bonded,cp);

/*=======================================================================*/
/*   X) Write to Screen         */

 if(myid==0){
  printf("\n");
  PRINT_LINE_DASH;
  printf("Completed preliminary tasks for CP-MD \n");
  PRINT_LINE_STAR;printf("\n");
 }/*endif*/

/*-----------------------------------------------------------------------*/
}/*end function*/
/*==========================================================================*/















