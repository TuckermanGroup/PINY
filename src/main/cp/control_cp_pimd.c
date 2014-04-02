/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
/*                                                                          */
/*                         PI_MD:                                           */
/*             The future of simulation technology                          */
/*             ------------------------------------                         */
/*                     Module: control_cp_pimd                              */
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
#include "../typ_defs/typedefs_par.h"
#include "../typ_defs/typedefs_stat.h"
#include "../proto_defs/proto_main_cp_local.h"
#include "../proto_defs/proto_main_local.h"
#include "../proto_defs/proto_math.h"
#include "../proto_defs/proto_integrate_cppimd_entry.h"
#include "../proto_defs/proto_integrate_cppimd_local.h"
#include "../proto_defs/proto_output_cp_entry.h"
#include "../proto_defs/proto_output_entry.h"
#include "../proto_defs/proto_energy_ctrl_cp_entry.h"
#include "../proto_defs/proto_vel_sampl_cp_entry.h"
#include "../proto_defs/proto_vel_sampl_class_entry.h"
#include "../proto_defs/proto_analysis_cp_entry.h"
#include "../proto_defs/proto_intra_con_entry.h"
#include "../proto_defs/proto_integrate_pimd_entry.h"
#include "../proto_defs/proto_integrate_pimd_local.h"
#include "../proto_defs/proto_integrate_md_entry.h"
#include "../proto_defs/proto_integrate_cp_local.h"
#include "../proto_defs/proto_energy_cpcon_entry.h"
#include "../proto_defs/proto_pimd_entry.h"
#include "../proto_defs/proto_communicate_wrappers.h"

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void control_cp_pimd(CLASS *class,BONDED *bonded,
                     GENERAL_DATA *general_data,CP *cp,
                     ANALYSIS *analysis)

/*=======================================================================*/
/*            Begin subprogram:                                          */
{   /*begin routine*/
/*=======================================================================*/
/*            Local variable declarations                                */
#include "../typ_defs/typ_mask.h"

  int itime,iii,i,itol_ok,ireset,ip,j;        

  int ntime         = general_data->timeinfo.ntime;
  int np            = class->communicate.np;
  int np_beads      = class->communicate.np_beads;
  int myid          = class->communicate.myid;
  int myid_state    = class->communicate.myid_state;
  int pi_beads_proc = class->clatoms_info.pi_beads_proc;
  int cp_elf_calc_frq = cp->cpcoeffs_info.cp_elf_calc_frq;
  double fc_mag_up  = 1000.0;
  double fc_mag_dn  = 1000.0;
  double tol_coef   = cp->cpopts.tol_coef;
  MPI_Comm world    = class->communicate.world;

/*=======================================================================*/
/* 0) Preliminary MD stuff                                               */

  prelim_cp_pimd(class,bonded,general_data,cp);

/*======================================================================*/
/* I) Write to Screen                                                   */

if(myid==0){
  PRINT_LINE_STAR;
  printf("Running CP-PIMD \n");
  PRINT_LINE_DASH;
}/*endif*/

/*======================================================================*/
/* II) Loop over the specified number of time steps */


  general_data->stat_avg.updates = 0.0;
  for(itime = 1;itime<=ntime;itime++){
    cputime(&(general_data->stat_avg.cpu1)); 
    (general_data->timeinfo.itime) = itime;
    (class->energy_ctrl.itime)     = itime;
  /*---------------------------------------------------------------------*/
  /*   1)Do NVT dynamics                                                 */


    if((general_data->ensopts.nvt)==1){
      if(((general_data->timeinfo.int_res_tra)==0) && 
        ((general_data->timeinfo.int_res_ter)==0)){
         int_NVT_cp_pimd(class,bonded,general_data,cp);
      }/*endif*/
    }/*endif*/

  /*---------------------------------------------------------------------*/
  /*   2)Do NPTI dynamics                                                 */
    if((general_data->ensopts.npt_i)==1){
      if(((general_data->timeinfo.int_res_tra)==0) && 
        ((general_data->timeinfo.int_res_ter)==0)){
         int_NPTI_cp_pimd(class,bonded,general_data,cp);
      }/*endif*/
    }/*endif*/

  /*---------------------------------------------------------------------*/
  /*   3)Do NPTF dynamics                                                 */
    if((general_data->ensopts.npt_f)==1){
      if(((general_data->timeinfo.int_res_tra)==0) && 
        ((general_data->timeinfo.int_res_ter)==0)){
         int_NPTF_cp_pimd(class,bonded,general_data,cp);
      }/*endif*/
    }/*endif*/

  /*---------------------------------------------------------------------*/
  /*  4)Rotate the orbitals:                                             */

    cp_rotate_control(cp,itime,ntime,&(general_data->stat_avg));

 /*----------------------------------------------------------------------*/
 /*  5) Check force tolerance                                            */


      check_coef_grad_mag(cp,&(general_data->simopts),
                          &fc_mag_up,&fc_mag_dn,&ireset,&itol_ok,
                          tol_coef,1,pi_beads_proc,&(general_data->stat_avg));

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
  /*   6)Calculate some simple averages                                   */
    cputime(&(general_data->stat_avg.cpu2)); 
    (general_data->stat_avg.cpu_now)=(general_data->stat_avg.cpu2)-
                                     (general_data->stat_avg.cpu1);
    (general_data->stat_avg.acpu) += (general_data->stat_avg.cpu2)-
                                     (general_data->stat_avg.cpu1);

    simpavg_cp_pimd(&(general_data->timeinfo),&(general_data->stat_avg),
                    &(general_data->cell),&(bonded->constrnt),
                    &(general_data->ensopts),&(general_data->simopts),
                    &(general_data->ptens),cp,&(class->communicate));

  /*-----------------------------------------------------------------------*/
  /*   7)Produce the output specified by the user                          */

    if(myid==0) check_auto_exit(&(general_data->timeinfo.exit_flag));
    if(np > 1) Bcast(&(general_data->timeinfo.exit_flag),1,MPI_INT,0,world);

    if(  (itime % (general_data->filenames.iwrite_screen))==0 ||
         (itime % (general_data->filenames.iwrite_dump  ))==0 ||
         (itime % (general_data->filenames.iwrite_confp ))==0 ||
         (itime % (general_data->filenames.iwrite_confv)) ==0 ||
         (itime % (general_data->filenames.iwrite_kseigs)) ==0 ||
         ((cp_elf_calc_frq > 0) && (itime % cp_elf_calc_frq) ==0) || 
         (itime % (general_data->filenames.iwrite_inst))  ==0 ||
         (general_data->timeinfo.exit_flag == 1)       ){
         (general_data->filenames.ifile_open) = 0;
	 if(myid_state==0&&np_beads>1){
           communicate_output_pimd(class);
	 }/*endif : myid_state==0*/
         output_cp_pimd(class,general_data,bonded,cp);
    }/*endif*/


  /*---------------------------------------------------------------------*/
  /*   8)Resample the particle velocities if desired                     */
  if(myid_state==0){
    if(((class->vel_samp_class.ivel_smpl_on)==1) &&
       ((class->vel_samp_class.nvx_smpl)!=0) && 
        (general_data->simopts.cp+general_data->simopts.cp_pimd==1)){
      if((itime % (class->vel_samp_class.nvx_smpl))==0){
       control_vx_smpl(class,bonded,&(general_data->ensopts),
                     &(general_data->simopts),general_data,
                       general_data->error_check_on);
      }/*endif*/
    }/*endif*/

  /*---------------------------------------------------------------------*/
  /*   9) Resample the extended class velocities if desired             */
    if(((class->vel_samp_class.ivel_smpl_on)==1) &&
       ((class->vel_samp_class.nvnhc_smpl)!=0)&& 
        (general_data->simopts.cp+general_data->simopts.cp_pimd==1)){
      if((itime % (class->vel_samp_class.nvnhc_smpl))==0){
       control_vnhc_smpl(class,general_data);
      }/*endif*/
    }/*endif*/
  }/*endif : myid_state==0*/

  /*---------------------------------------------------------------------*/
  /*   10) Resample the coefficient velocities if desired                */
    if(((cp->vel_samp_cp.ivelc_smpl_on)==1) &&
       ((cp->vel_samp_cp.nvc_smpl)!=0)){
      if((itime % (cp->vel_samp_cp.nvc_smpl))==0){
       control_vc_smpl(general_data,cp);
      }/* endif */
    }/* endif */

  /*---------------------------------------------------------------------*/
  /*   11) Resample the extended class of coefficients if desired        */
    if(((cp->vel_samp_cp.ivelc_smpl_on)==1) &&
       ((cp->vel_samp_cp.nvcnhc_smpl)!=0)){
      if((itime % (cp->vel_samp_cp.nvcnhc_smpl))==0){
       control_vcnhc_smpl(cp);
      }/*endif*/
    }/*endif*/

  /*---------------------------------------------------------------------*/
  /*  12) Analysis Routine                                               */
    analysis_cp_pimd(class,bonded,general_data,cp,analysis); 

  /*---------------------------------------------------------------------*/
  /*   13) Check for exit condition                                      */

   if(general_data->timeinfo.exit_flag == 1) itime = ntime;

  }/*endfor:itime */

/*======================================================================*/
/*  III)Write to Screen                                                 */

if(myid==0){
  PRINT_LINE_DASH;
  printf("Completed CP-MD run \n");
  PRINT_LINE_STAR;
}/*endif : myid=0*/

/*-----------------------------------------------------------------------*/
}/*end routine*/
/*==========================================================================*/





/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void prelim_cp_pimd(CLASS *class,BONDED *bonded,
                             GENERAL_DATA *general_data,CP *cp)

/*=======================================================================*/
/*            Begin subprogram:                                          */
{   /*begin routine*/
/*=======================================================================*/
/*            Local variable declarations                                */

#include "../typ_defs/typ_mask.h"

  int i,ip,iflag,iii,j;        
  int iflag_mass = 1;
  int pi_beads      = class->clatoms_info.pi_beads;
  int pi_beads_proc = class->clatoms_info.pi_beads_proc;
  int myid       = class->communicate.myid;
  int myid_bead  = class->communicate.myid_bead;
  int myid_state = class->communicate.myid_state;
  int error_check_on = general_data->error_check_on;
  double *fx,*fy,*fz;
  double kinet_temp=0.0;
  double vintrat_temp=0.0;
  double vintert_temp=0.0;
  double vrecip_temp=0.0;
  double vvdw_temp=0.0;
  double vcoul_temp=0.0;
  double kinet_nhc_bead_temp=0.0;
  double kinet_nhc_cp_temp=0.0;
  double vbondt_temp=0.0;
  double vbendt_temp=0.0;
  double vbend_bnd_bond_temp=0.0;
  double vbend_bnd_bend_temp=0.0;
  double vtorst_temp=0.0;
  double vonfot_temp = 0.0;
  double pi_ke_prim_temp = 0.0;
  double pi_ke_vir_temp = 0.0;
  double kin_harm_temp = 0.0;
  double cp_ehart_temp = 0.0;
  double cp_eext_temp = 0.0;
  double cp_exc_temp = 0.0;
  double cp_muxc_temp = 0.0;
  double cp_eke_temp = 0.0;
  double cp_enl_temp = 0.0;
  double kinet_cp_temp = 0.0;
  double kinet_cp_up_tmp,kinet_cp_dn_tmp;
  MPI_Comm world       = class->communicate.world;
  MPI_Comm comm_beads  = class->communicate.comm_beads;
  MPI_Comm comm_states = class->communicate.comm_states;
  int num_proc         = class->communicate.np;
  int start_proc;
      start_proc = (myid_bead == 0 ? 2:1);

/*=======================================================================*/
/*   I) Write to Screen                                                  */

if(error_check_on==1){
  PRINT_LINE_STAR;
  printf("Performing preliminary tasks for CP-PIMD \n");
  PRINT_LINE_DASH;printf("\n");
}/*endif*/

/*=======================================================================*/
/*  II) Initialize In-output                                             */

  general_data->stat_avg.updates   = 0.0;
  general_data->stat_avg.updates_w = 0.0;
  general_data->stat_avg.acpu      = 0.0;
  general_data->stat_avg.vbondt    = 0.0;
  general_data->stat_avg.vbendt    = 0.0;
  general_data->stat_avg.vtorst    = 0.0;
  general_data->stat_avg.vonfot    = 0.0;
  general_data->stat_avg.vintrat   = 0.0;
  general_data->stat_avg.vintert   = 0.0;
  general_data->stat_avg.kinet     = 0.0;
  general_data->stat_avg.cp_eke    = 0.0;
  general_data->stat_avg.count_diag_srot      = 0.0;

/*=======================================================================*/
/* III)Get Energy and particle Forces                                    */

if(myid_state==0){
  if(bonded->constrnt.iconstrnt==1&&general_data->simopts.cp==1){
     init_constraint(bonded,&(general_data->ptens));
  }/*endif:init constraints, must be done before getting the energy*/
}/*endif : myid_state==0*/

    general_data->stat_avg.iter_shake   = 0;
    general_data->stat_avg.iter_ratl    = 0;
    general_data->stat_avg.iter_23    = 0;
    general_data->stat_avg.iter_33    = 0;
    general_data->stat_avg.iter_46    = 0;

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

  class->clatoms_info.wght_pimd = 1;

if(myid_state==0){
  control_pimd_trans_mode(class,general_data);
  mode_energy_control(class,general_data);
  control_pimd_trans_pos(class,general_data);
}/*endif : myid_state==0*/
  

    cp_energy_control_pimd(class,bonded,general_data,cp);

   if(cp->cpopts.cp_gauss == 1){
    for(ip=1;ip<=pi_beads_proc;ip++){
     add_gauss_force(&(cp->cpcoeffs_info),&(cp->cpcoeffs_pos[ip]),
                    &(cp->cpscr.cpscr_ovmat),&(cp->cpopts));
   }/* endfor ip*/
  }/* endif */
 if(myid_state==0){
  if(general_data->simopts.cp==1||general_data->simopts.cp_pimd==1){
   get_tvten_pimd(class,general_data);
  }
 }/*endif : myid_state==0*/
  get_cpke_pimd(&(cp->cpcoeffs_info),(cp->cpcoeffs_pos),
                &(general_data->stat_avg),cp->cpopts.cp_lsda,
                  class->clatoms_info.pi_beads_proc,
                  class->communicate.np_states);


/*=======================================================================*/
/* IV)Get Atm NHC Forces                                                     */

 if(myid_state==0){
  if(general_data->simopts.cp==1||general_data->simopts.cp_pimd==1){

   if((general_data->ensopts.nvt)==1){
     iflag_mass = 1;
     if(class->clatoms_info.pi_beads_proc_st==1){
       init_NHC_par(&(class->clatoms_info),&(class->clatoms_pos[1]),
              &(class->therm_info_class),&(class->therm_class),
              &(class->int_scr),iflag_mass,&(class->class_comm_forc_pkg));
     }/*endif*/
     iflag_mass = 2;
     for(ip=start_proc;ip<=pi_beads_proc;ip++){
      init_NHC_par(&(class->clatoms_info),&(class->clatoms_pos[ip]),
                    &(class->therm_info_bead),&(class->therm_bead[ip]),
                    &(class->int_scr),iflag_mass,&(class->class_comm_forc_pkg));
     }/*endfor*/
     iflag = 0;
     nhc_vol_potkin_pimd(class,general_data,iflag);
   }/*endif*/

   if((general_data->ensopts.npt_i)==1){
     if(class->clatoms_info.pi_beads_proc_st==1){
       init_NHCPI_par(&(class->clatoms_info),&(class->clatoms_pos[1]),
                &(class->therm_info_class),&(class->therm_class),
                &(general_data->baro),
                &(class->int_scr),&(class->class_comm_forc_pkg));
     }/*endif*/
     iflag_mass = 2;
     for(ip=start_proc;ip<=pi_beads_proc;ip++){
      init_NHC_par(&(class->clatoms_info),&(class->clatoms_pos[ip]),
                    &(class->therm_info_bead),&(class->therm_bead[ip]),
                    &(class->int_scr),iflag_mass,&(class->class_comm_forc_pkg));
     }/*endfor*/
    iflag = 1;
    nhc_vol_potkin_pimd(class,general_data,iflag);
   }/*endif*/

  if((general_data->ensopts.npt_f)==1){
     if(class->clatoms_info.pi_beads_proc_st==1){
       init_NHCPF_par(&(class->clatoms_info),&(class->clatoms_pos[1]),
                &(class->therm_info_class),&(class->therm_class),
                &(general_data->baro),
                &(general_data->par_rahman),&(general_data->cell),
                &(class->int_scr),&(class->class_comm_forc_pkg));
     }/*endif*/
     iflag_mass = 2;
     for(ip=start_proc;ip<=pi_beads_proc;ip++){
      init_NHC_par(&(class->clatoms_info),&(class->clatoms_pos[ip]),
                    &(class->therm_info_bead),&(class->therm_bead[ip]),
                    &(class->int_scr),iflag_mass,&(class->class_comm_forc_pkg));
     }/*endfor*/
    iflag = 2;
    nhc_vol_potkin_pimd(class,general_data,iflag);
  }/*endif*/

    if((general_data->ensopts.npt_i+general_data->ensopts.npt_f)==0){
       general_data->stat_avg.kinet_v       = 0.0;
    }/*endif*/
    if((general_data->ensopts.npt_i+general_data->ensopts.npt_f
       +general_data->ensopts.nvt)==0){
       general_data->stat_avg.kinet_nhc     = 0.0;
    }/*endif*/

  }/*endif: full cp on*/

}/*endif : myid_state==0*/

/*=======================================================================*/
/* IV)Get CP NHC Forces                                                   */

  if(cp->cptherm_info.num_c_nhc > 0){
   if(cp->cptherm_info.massiv_flag==0){
    for(ip=1;ip<=pi_beads_proc;ip++){
     init_cp_NHC(&(cp->cptherm_info),&(cp->cptherm_pos[ip]),
                 &(cp->cpcoeffs_info),&(cp->cpcoeffs_pos[ip]),
                 &(cp->cpscr),cp->cpopts.cp_lsda,comm_states,
                 cp->communicate.np_states,cp->communicate.myid_state);
    }/* endfor */
    nhc_cp_potkin_pimd(&(cp->cptherm_info),(cp->cptherm_pos),
                       &(general_data->stat_avg),pi_beads_proc);
   }else{
    nhc_cp_potkin_massiv_pimd(&(cp->cptherm_info),(cp->cptherm_pos),
                              &(general_data->stat_avg),pi_beads_proc);
   }/*endif*/
    general_data->stat_avg.vpotnhc_cp   /= pi_beads;
    general_data->stat_avg.kinet_nhc_cp /= pi_beads;
  }/* endif */

/*========================================================================*/
/* VI) Initialize free energy stuff                                       */

if(myid_state==0){
  if(general_data->simopts.cp==1||general_data->simopts.cp_pimd==1){
     if(bonded->bond_free.num != 0){
       for(i=1;i<= (bonded->bond_free.num);i++){
         bonded->bond_free.hist[i] = 0.0;
       }
     }/*endif*/
     if(bonded->bend_free.num != 0){
       for(i=1;i<= (bonded->bend_free.num);i++){
         bonded->bend_free.hist[i] = 0.0;
       }
     }/*endif*/
     if(bonded->tors_free.num == 1){
       for(i=1;i<= bonded->tors_free.num;i++){
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

  }/*endif*/
}/*endif : myid_state==0*/

/*=======================================================================*/
/*  IV) Write Energy to screen                                           */

  /*=======================================================================*/
  /*  Collect initial energies (needs its own routine)                     */
  /*=======================================================================*/

if(num_proc>1){

 if(myid_state==0){
  Reduce(&(general_data->stat_avg.vvdw), &vvdw_temp,1,MPI_DOUBLE,
                   MPI_SUM,0,comm_beads);
  Reduce(&(general_data->stat_avg.vintrat), &vintrat_temp,1,MPI_DOUBLE,
                   MPI_SUM,0,comm_beads);
  Reduce(&(general_data->stat_avg.kinet), &kinet_temp,1,MPI_DOUBLE,
                   MPI_SUM,0,comm_beads);
  Reduce(&(general_data->stat_avg.kinet_nhc_bead), &kinet_nhc_bead_temp,1,
                   MPI_DOUBLE,MPI_SUM,0,comm_beads);

  Reduce(&(general_data->stat_avg.vbondt), &vbondt_temp,1,
                   MPI_DOUBLE,MPI_SUM,0,comm_beads);
  Reduce(&(general_data->stat_avg.vbendt), &vbendt_temp,1,
                   MPI_DOUBLE,MPI_SUM,0,comm_beads);
  Reduce(&(general_data->stat_avg.vbend_bnd_bond), &vbend_bnd_bond_temp,1,
                   MPI_DOUBLE,MPI_SUM,0,comm_beads);
  Reduce(&(general_data->stat_avg.vbend_bnd_bend), &vbend_bnd_bend_temp,1,
                   MPI_DOUBLE,MPI_SUM,0,comm_beads);
  Reduce(&(general_data->stat_avg.vtorst), &vtorst_temp,1,
                   MPI_DOUBLE,MPI_SUM,0,comm_beads);
  Reduce(&(general_data->stat_avg.vonfot), &vonfot_temp,1,
                   MPI_DOUBLE,MPI_SUM,0,comm_beads);
  Reduce(&(general_data->stat_avg.pi_ke_prim), &pi_ke_prim_temp,1,
                   MPI_DOUBLE,MPI_SUM,0,comm_beads);
  Reduce(&(general_data->stat_avg.pi_ke_vir), &pi_ke_vir_temp,1,
                   MPI_DOUBLE,MPI_SUM,0,comm_beads);
  Reduce(&(general_data->stat_avg.kin_harm), &kin_harm_temp,1,
                   MPI_DOUBLE,MPI_SUM,0,comm_beads);
 }/*endif : myid_state==0*/

  Reduce(&(general_data->stat_avg.cp_ehart), &cp_ehart_temp,1,
                   MPI_DOUBLE,MPI_SUM,0,world);
  Reduce(&(general_data->stat_avg.cp_eext), &cp_eext_temp,1,
                   MPI_DOUBLE,MPI_SUM,0,world);
  Reduce(&(general_data->stat_avg.cp_exc), &cp_exc_temp,1,
                   MPI_DOUBLE,MPI_SUM,0,world);
  Reduce(&(general_data->stat_avg.cp_muxc), &cp_muxc_temp,1,
                   MPI_DOUBLE,MPI_SUM,0,world);
  Reduce(&(general_data->stat_avg.cp_eke), &cp_eke_temp,1,
                   MPI_DOUBLE,MPI_SUM,0,world);
  Reduce(&(general_data->stat_avg.cp_enl), &cp_enl_temp,1,
                   MPI_DOUBLE,MPI_SUM,0,world);
  Reduce(&(general_data->stat_avg.kinet_cp), &kinet_cp_temp,1,
                   MPI_DOUBLE,MPI_SUM,0,world);
  Reduce(&(general_data->stat_avg.kinet_nhc_cp), &kinet_nhc_cp_temp,1,
                   MPI_DOUBLE,MPI_SUM,0,world);

  general_data->stat_avg.kinet          = kinet_temp;
  general_data->stat_avg.kinet_nhc_bead = kinet_nhc_bead_temp;
  general_data->stat_avg.vvdw           = vvdw_temp;
  general_data->stat_avg.vintrat        = vintrat_temp;
  general_data->stat_avg.vbondt         = vbondt_temp;
  general_data->stat_avg.vbendt         = vbendt_temp;
  general_data->stat_avg.vbend_bnd_bond = vbend_bnd_bond_temp;
  general_data->stat_avg.vbend_bnd_bend = vbend_bnd_bend_temp;
  general_data->stat_avg.vtorst         = vtorst_temp;
  general_data->stat_avg.vonfot         = vonfot_temp;
  general_data->stat_avg.pi_ke_prim     = pi_ke_prim_temp;
  general_data->stat_avg.pi_ke_vir      = pi_ke_vir_temp;
  general_data->stat_avg.kin_harm       = kin_harm_temp;
  general_data->stat_avg.cp_ehart       = cp_ehart_temp;
  general_data->stat_avg.cp_eext        = cp_eext_temp;
  general_data->stat_avg.cp_exc         = cp_exc_temp;
  general_data->stat_avg.cp_muxc        = cp_muxc_temp;
  general_data->stat_avg.cp_eke         = cp_eke_temp;
  general_data->stat_avg.cp_enl         = cp_enl_temp;
  general_data->stat_avg.kinet_cp       = kinet_cp_temp;
  general_data->stat_avg.kinet_nhc_cp   = kinet_nhc_cp_temp;

}/*endif*/

  general_data->filenames.ifile_open = 1;

  general_data->timeinfo.itime = 0;
  
  output_cp_pimd(class,general_data,bonded,cp);


/*=======================================================================*/
/* VIII) If doing CP isokinetic dynamics, set the fixed value of the
       fictitious kinetic energy using the initial value                */

  if(cp->cpopts.cp_isok_opt == 1){
    if(cp->communicate.np_states > 1){

      kinet_cp_up_tmp = general_data->stat_avg.kinet_cp_up;
      Allreduce(&kinet_cp_up_tmp,&(general_data->stat_avg.kinet_cp_up),
                1,MPI_DOUBLE,MPI_SUM,0,world);

      if(cp->cpopts.cp_lsda == 1) {
        kinet_cp_dn_tmp = general_data->stat_avg.kinet_cp_dn;
        Allreduce(&kinet_cp_dn_tmp,&(general_data->stat_avg.kinet_cp_dn),
                  1,MPI_DOUBLE,MPI_SUM,0,world);

      }
    }/* endif parallel */

    cp->cpcoeffs_info.k0_up = general_data->stat_avg.kinet_cp_up;

    if(cp->cpopts.cp_lsda == 1) cp->cpcoeffs_info.k0_dn = general_data->stat_avg.kinet_cp_dn;
  }/* endif isok_opt on */

/*=======================================================================*/
/*   VII) Write to Screen         */

if(error_check_on==1){
  printf("\n");
  PRINT_LINE_DASH;
  printf("Completed preliminary tasks for CP-PIMD \n");
  PRINT_LINE_STAR;printf("\n");
}/*endif*/

/*-----------------------------------------------------------------------*/
   }/*end routine*/
/*==========================================================================*/













