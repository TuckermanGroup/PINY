/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
/*                                                                          */
/*                         PI_MD:                                           */
/*             The future of simulation technology                          */
/*             ------------------------------------                         */
/*                     Module: control_min                                  */
/*                                                                          */
/* This subprogram performs MD on a classical potential energy surface (PES)*/
/*                                                                          */
/*                                                                          */
/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/


#include "standard_include.h"
#include "../typ_defs/typedefs_class.h"
#include "../typ_defs/typedefs_par.h"
#include "../typ_defs/typedefs_bnd.h"
#include "../typ_defs/typedefs_gen.h"
#include "../typ_defs/typedefs_cp.h"
#include "../typ_defs/typedefs_stat.h"
#include "../proto_defs/proto_main_local.h"
#include "../proto_defs/proto_main_cp_local.h"
#include "../proto_defs/proto_math.h"
#include "../proto_defs/proto_integrate_min_entry.h"
#include "../proto_defs/proto_integrate_min_local.h"
#include "../proto_defs/proto_output_entry.h"
#include "../proto_defs/proto_output_local.h"
#include "../proto_defs/proto_analysis_md_entry.h"
#include "../proto_defs/proto_vel_sampl_class_entry.h"
#include "../proto_defs/proto_intra_con_entry.h"
#include "../proto_defs/proto_energy_ctrl_entry.h"
#include "../proto_defs/proto_communicate_wrappers.h"


/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void control_min(CLASS *class,BONDED *bonded,GENERAL_DATA *general_data,
                 ANALYSIS *analysis)

/*=======================================================================*/
/*            Begin subprogram:                                          */
{   /*begin routine*/
/*=======================================================================*/
/*            Local variable declarations                                */
#include "../typ_defs/typ_mask.h"

  int itime,iii;        
  double f_atm_mag;
  int ireset=1,idone=0,iexit=0;
  int num_proc = class->communicate.np;
  int myid     = class->communicate.myid;
  MPI_Comm world = class->communicate.world;

/*======================================================================*/
/* 0) preliminary MINIMIZE stuff                                        */

  prelim_min(class,bonded,general_data);

/*======================================================================*/
/* I) Write to Screen                                                   */

  PRINT_LINE_STAR;
  printf("Running MINIMIZATION \n");
  PRINT_LINE_DASH;

/*======================================================================*/
/* II) Loop over the specified number of time steps */

  general_data->stat_avg.updates = 0.0;
  (class->energy_ctrl.itime)     = 0; /* get the total PE every time */
  for(itime = 1;itime<=(general_data->timeinfo.ntime);itime++){
    cputime(&(general_data->stat_avg.cpu1)); 
    (general_data->timeinfo.itime) = itime;
  /*---------------------------------------------------------------------*/
  /*   1)Do Steepest descent minimization                                */
    if((general_data->minopts.min_std)==1){
       min_STD(class,bonded,general_data);
    }/*endif*/
  /*---------------------------------------------------------------------*/
  /*   2)Do Steepest descent minimization                                */
    if((general_data->minopts.min_cg)==1){
       min_CG(class,bonded,general_data,ireset);
    }/*endif*/
  /*----------------------------------------------------------------------*/
  /*----------------------------------------------------------------------*/
  /*   5)Check magnitude of atomic forces */
      check_atm_grad_mag(class,general_data,&f_atm_mag,&ireset,&idone,
                         general_data->minopts.tol_atom);

      if(idone==1){
       printf("\n@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n");
       printf("Tolerance on atomic forces reached\n");
       printf("@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n");
       iexit = 1;
      }

      if (general_data->minopts.min_cg == 1) {
        if (ireset == 1) {
          printf("\n@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n");
          printf("Resetting Atomic Conjugate Gradient\n");
          printf("@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n");
        }
      }

  /*----------------------------------------------------------------------*/
  /*   6)Calculate some simple averages                                   */
    cputime(&(general_data->stat_avg.cpu2)); 
    (general_data->stat_avg.cpu_now)=(general_data->stat_avg.cpu2)-
                                     (general_data->stat_avg.cpu1);
    (general_data->stat_avg.acpu) += (general_data->stat_avg.cpu2)-
                                     (general_data->stat_avg.cpu1);

    simpavg_md(&(general_data->timeinfo),&(general_data->stat_avg),
	       &(general_data->cell),&(bonded->constrnt),
	       &(general_data->ensopts),&(general_data->simopts),
	       &(general_data->ptens),&(class->communicate),
               &(class->nbr_list.verlist),&(class->energy_ctrl));
  /*-----------------------------------------------------------------------*/
  /*   7)Produce the output specified by the user                          */

    if(myid==0) check_auto_exit(&(general_data->timeinfo.exit_flag));
    if(num_proc > 1) Bcast(&(general_data->timeinfo.exit_flag),1,MPI_INT,0,world);

    if(  (itime % (general_data->filenames.iwrite_screen))==0 ||
         (itime % (general_data->filenames.iwrite_dump  ))==0 ||
         (itime % (general_data->filenames.iwrite_confp ))==0 ||
         (itime % (general_data->filenames.iwrite_confv)) ==0 ||
         (itime % (general_data->filenames.iwrite_inst))  ==0   ){
         (general_data->filenames.ifile_open) = 0;

          output_min(class,general_data,bonded);
	  if(idone == 1 && class->communicate.myid == 0 ){
            write_dump_file_md(class,bonded,general_data);
          }
          if(num_proc>1){ Barrier(class->communicate.world);}
    }/*endif*/
  /*---------------------------------------------------------------------*/
  /*  10) Analysis Routine                                               */
    analysis_md(class,general_data,bonded,analysis); 
  /*---------------------------------------------------------------------*/
   if(iexit==1) break;
  }/*endfor:itime */

  /*======================================================================*/
  /*  III)Write to Screen                                                 */
  PRINT_LINE_DASH;
  printf("Completed MD run \n");
  PRINT_LINE_STAR;

/*-----------------------------------------------------------------------*/
}/*end routine*/
/*==========================================================================*/





/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void prelim_min(CLASS *class,BONDED *bonded,GENERAL_DATA *general_data)

/*=======================================================================*/
/*            Begin subprogram:                                          */
{   /*begin routine*/
/*=======================================================================*/
/*            Local variable declarations                                */

#include "../typ_defs/typ_mask.h"

  int i,j,k,m,iflag,iii;

  int npairs_tmp             = 0;
  int npairs_res_tmp         = 0;
  int iflag_mass             = 1;
  double kinet_temp          = 0.0;
  double vintrat_temp        = 0.0;
  double vintert_temp        = 0.0;
  double vbondt_temp         = 0.0;
  double vbendt_temp         = 0.0;
  double vbend_bnd_bond_temp = 0.0;
  double vbend_bnd_bend_temp = 0.0;
  double vtorst_temp         = 0.0;
  double vonfot_temp         = 0.0;
  double vcoul_temp          = 0.0;
  double vvdw_temp           = 0.0;

  int myid                   = class->communicate.myid;
  MPI_Comm comm_forc         = class->communicate.comm_forc;
  int num_proc               = class->communicate.np;
  int np_forc                = class->class_comm_forc_pkg.num_proc; 

/*=======================================================================*/
/* I) Write to Screen                                                  */

  if(myid==0){
    PRINT_LINE_STAR;
    printf("Performing preliminary tasks for minimization\n");
    PRINT_LINE_DASH;printf("\n");
  }/*endif*/

/*=======================================================================*/
/*  II) Initialize In-output                                             */

  general_data->stat_avg.updates   = 0.0;
  general_data->stat_avg.acpu      = 0.0;

/*=======================================================================*/
/* III)Get Energy and particle Forces                                    */

  if(bonded->constrnt.iconstrnt==1){
     init_constraint(bonded,&(general_data->ptens));
  }/*endif:init constraints, must be done before getting the energy*/

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

  if(myid==0){
    printf("Getting initial energy\n");
  }/*endif*/
 
  energy_control(class,bonded,general_data);

/*=======================================================================*/
/*  VII) Write Energy to screen                                           */

  general_data->filenames.ifile_open = 1;
  general_data->timeinfo.itime       = 0;
  output_min(class,general_data,bonded);

/*=======================================================================*/
/* VIII) Reinitailize update flag */

  general_data->stat_avg.updates   = 0.0;

/*=======================================================================*/
/* IX) Write to Screen         */

  if(myid==0){
    printf("\n");
    PRINT_LINE_DASH;
    printf("Completed preliminary tasks for minimization\n");
    PRINT_LINE_STAR;printf("\n");
  }/*endif*/

/*-----------------------------------------------------------------------*/
  }/*end routine*/
/*========================================================================*/
