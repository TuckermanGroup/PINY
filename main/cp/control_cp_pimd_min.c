/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
/*                                                                          */
/*                         PI_MD:                                           */
/*             The future of simulation technology                          */
/*             ------------------------------------                         */
/*                     Module: control_cp                                   */
/*                                                                          */
/* This subprogram performs Minization on a classical+abinitio potential    */
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
#include "../proto_defs/proto_main_local.h"
#include "../proto_defs/proto_main_cp_local.h"
#include "../proto_defs/proto_math.h"
#include "../proto_defs/proto_integrate_cpmin_entry.h"
#include "../proto_defs/proto_energy_cpcon_entry.h"
#include "../proto_defs/proto_output_entry.h"
#include "../proto_defs/proto_output_cp_entry.h"
#include "../proto_defs/proto_output_cp_local.h"
#include "../proto_defs/proto_analysis_md_entry.h"
#include "../proto_defs/proto_analysis_cp_entry.h"
#include "../proto_defs/proto_communicate_wrappers.h"


/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void control_cp_pimd_min(CLASS *class,BONDED *bonded,
                         GENERAL_DATA *general_data,CP *cp,
                         ANALYSIS *analysis)

/*=======================================================================*/
/*            Begin subprogram:                                          */
{   /*begin routine*/
/*=======================================================================*/
/*            Local variable declarations                                */

#include "../typ_defs/typ_mask.h"

  double fc_mag_up,fc_mag_dn,elec_e_old,elec_e;
  double elec_e_old_tmp,elec_e_tmp;
  double Delta_E;
  int ncoef_tot,i,itime,iii,ip,iexit_tot;        
  int atm_min;
  int iexit;
  int ireset,idone;
  int iatm_count;
  int pi_beads       = class->clatoms_info.pi_beads; 
  int pi_beads_proc  = class->clatoms_info.pi_beads_proc; 
  int myid           = class->communicate.myid;
  int myid_state     = class->communicate.myid_state;
  int myid_bead      = class->communicate.myid_bead;
  int nproc          = class->communicate.np;
  int nproc_beads    = class->communicate.np_beads;
  int error_check_on = general_data->error_check_on;
  MPI_Comm world     = class->communicate.world;
  int num_proc = cp->communicate.np;
  int np_states = cp->communicate.np_states;
  MPI_Comm comm_states = class->communicate.comm_states;
  double gamma,delta_E_elec;
  int ifirst_elec;
/*======================================================================*/
/* 0) Write to Screen                                                   */

 if(error_check_on==1){
  PRINT_LINE_STAR;
  printf("Running CP-MINIMIZATION \n");
  PRINT_LINE_DASH;
 }/*endif : error_check_on*/

/*======================================================================*/
/*======================================================================*/
/* II) Loop over the beads                                              */

 for(ip=1;ip<=pi_beads_proc;ip++){

  /*======================================================================*/
  /* 1) Initialize counters for this bead                               */

   iexit = 0;
   iatm_count = 0;
   atm_min=0;
   general_data->stat_avg.write_cp_atm_flag = 0;
/*   general_data->simopts.cp_min           = 0;
   general_data->simopts.cp_wave_min      = 1;*/
   cp->cpcoeffs_info.cg_reset_flag   = 1;
   cp->cpcoeffs_info.diis_reset_flag = 1;
   fc_mag_up = 10000.0;
   fc_mag_dn = 10000.0;
   general_data->stat_avg.cp_eke = 1000.0;
   general_data->stat_avg.cp_enl = 1000.0;
   general_data->stat_avg.cp_ehart = 1000.0;
   general_data->stat_avg.cp_exc = 1000.0;
   general_data->stat_avg.cp_eext  = 1000.0;
   general_data->stat_avg.count_diag_srot      = 0.0;

  /*======================================================================*/
  /* 2) Loop over the specified number of time steps                      */
   for(itime = 1;itime<=(general_data->timeinfo.ntime);itime++){

   /*---------------------------------------------------------------------*/
   /* i) CPU time control                                                 */
    general_data->timeinfo.itime = itime;
    class->energy_ctrl.itime     = itime;
    cputime(&(general_data->stat_avg.cpu1)); 
   /*---------------------------------------------------------------------*/
   /* ii) CP_wave minimization                                             */
  /*---------------------------------------------------------------------*/
  /* 1) CP_wave minimization                                             */
    if(iexit==0){
      elec_e_old     = general_data->stat_avg.cp_eke
                     + general_data->stat_avg.cp_enl
                     + general_data->stat_avg.cp_ehart
                     + general_data->stat_avg.cp_exc
                     + general_data->stat_avg.cp_eext;
      if(np_states>1){
       elec_e_old_tmp = elec_e_old;
       Allreduce(&(elec_e_old_tmp),&(elec_e_old),1,MPI_DOUBLE,MPI_SUM,0,
                 comm_states);
      }/*endif*/
      if(general_data->minopts.cp_min_std==1){
         min_STD_cp(class,bonded,general_data,cp,ip);
         cp_shuffle_states(cp,ip);
      }/*endif*/
      if((general_data->minopts.cp_min_cg==1)){
         min_CG_cp(class,bonded,general_data,cp,ip);
         cp->cpcoeffs_info.cg_reset_flag = 0;
      }/*endif*/
      if((general_data->minopts.cp_min_diis==1)){
         min_DIIS_cp(class,bonded,general_data,cp,iatm_count,ip);
         cp->cpcoeffs_info.diis_reset_flag = 0;
      }/*endif*/
      elec_e         = general_data->stat_avg.cp_eke
                     + general_data->stat_avg.cp_enl
                     + general_data->stat_avg.cp_ehart
                     + general_data->stat_avg.cp_exc
                     + general_data->stat_avg.cp_eext;
      Delta_E = fabs(elec_e - elec_e_old);
      if(np_states>1){
       elec_e_tmp = elec_e;
       Allreduce(&(elec_e_tmp),&(elec_e),1,MPI_DOUBLE,MPI_SUM,0,comm_states);
      }/*endif*/
      check_coef_grad_mag(cp,&(general_data->simopts),
                          &fc_mag_up,&fc_mag_dn,&ireset,&idone,
                          general_data->minopts.tol_coef,ip,ip,
                          &(general_data->stat_avg));
      if(idone==1){
        if(nproc_beads > 1) {
           printf("\n@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n");
           printf("Tolerance on coefficient forces reached on processor %d\n",
                   myid_bead);
           printf("@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n");
        } else {
           printf("\n@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n");
           printf("Tolerance on coefficient forces reached\n");
           printf("@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n");
        }
        iexit = 1;
      }else{
         if((ireset==1)&&(general_data->minopts.cp_min_cg==1)&&
            (elec_e-elec_e_old)>0){
           cp->cpcoeffs_info.cg_reset_flag = 1;
           if(nproc_beads > 1 ){
              printf("\n@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n");
              printf("Resetting CP Conjugate Gradient on processor %d\n",myid_bead);
              printf("@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n");
           } else {
              printf("\n@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n");
              printf("Resetting CP Conjugate Gradient\n");
              printf("@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n");
           }
          cp_shuffle_states(cp,ip);
         }/*endif:reset cg*/
      }/*endif:not done*/

    }/*endif:cp_wave move*/

    /*----------------------------------------------------------------------*/
    /*  iii)Calculate some simple averages                                   */
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
    /*  iv)Produce the output specified by the user                          */

    if(myid==0) check_auto_exit(&(general_data->timeinfo.exit_flag));
    if(num_proc > 1) Bcast(&(general_data->timeinfo.exit_flag),1,MPI_INT,0,world);

    if(  (itime % (general_data->filenames.iwrite_screen))==0 ||
         (itime % (general_data->filenames.iwrite_dump  ))==0 ||
         (itime % (general_data->filenames.iwrite_confp ))==0 ||
         (itime % (general_data->filenames.iwrite_confc ))==0 ||
         (itime % (general_data->filenames.iwrite_confv)) ==0 ||
         (itime % (general_data->filenames.iwrite_inst))  ==0   ){
         general_data->filenames.ifile_open = 0;
         if(myid_state==0){
            communicate_output_pimd(class);
	 }/*endif*/
         if((itime % general_data->filenames.iwrite_screen) == 0){
            output_cp_min(class,general_data,bonded,cp,Delta_E,iexit);
         }
         if(myid==0 && (itime % general_data->filenames.iwrite_dump)==0){
             write_dump_file_cp_pimd_class(class,bonded,general_data,cp);
         }/*endif*/
         if((itime % general_data->filenames.iwrite_dump)==0){
           control_coef_transpose_bck(cp,1);
           write_dump_file_cp_pimd(class,bonded,general_data,cp);
           control_coef_transpose_fwd(cp,1);
         }/* endif */
    }/*endif*/
    /*---------------------------------------------------------------------*/
    /*  v) Analysis Routine                                               */
    analysis_cp_pimd(class,bonded,general_data,cp,analysis); 

    /*---------------------------------------------------------------------*/
    iexit_tot = 0; 
  if(nproc > 1){
    Allreduce(&iexit,&iexit_tot,1,MPI_INT,MPI_SUM,0,world);
  }else{
      iexit_tot = iexit;
  }
    if(iexit_tot==nproc){break;}

   }/*endfor:itime */

/*=========================================================================*/
/* 3) Dump and write error if not done                                 */

   if(myid_state==0){
     communicate_output_pimd(class);
   }/*endif*/
   if(myid==0){
     write_dump_file_cp_pimd_class(class,bonded,general_data,cp);
   }/*endif*/
   control_coef_transpose_bck(cp,1);
     write_dump_file_cp_pimd(class,bonded,general_data,cp);
   control_coef_transpose_fwd(cp,1);
   if(iexit!=1){
    printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n");    
    printf("Out of time on bead number %d on processor number %d\n",ip,myid);
    printf("before convergence. Restart your job and allow more steps/bead\n");
    printf("Don't panic, I have saved your results in the dump files \n");
    printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n");    
    fflush(stdout);
    exit(1);
    Finalize();
   }/*endif*/

 }/*endfor:pi_beads */

/*======================================================================*/
/*======================================================================*/
/*  III)Write to Screen                                                 */

  if(error_check_on==1){
    PRINT_LINE_DASH;
    printf("Completed CP-MINIMIZATION run \n");
    PRINT_LINE_STAR;
  }/*endif : error_check_on*/

/*-----------------------------------------------------------------------*/
   }/*end routine*/
/*==========================================================================*/





































































