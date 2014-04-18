/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
/*                                                                          */
/*                         PI_MD:                                           */
/*             The future of simulation technology                          */
/*             ------------------------------------                         */
/*                     Module: control_md                                   */
/*                                                                          */
/* This subprogram performs MD on a classical potential energy surface (PES)*/
/*                                                                          */
/*                                                                          */
/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/


#include "standard_include.h"
#include "../typ_defs/typedefs_class.h"
#include "../typ_defs/typedefs_bnd.h"
#include "../typ_defs/typedefs_gen.h"
#include "../typ_defs/typedefs_par.h"
#include "../typ_defs/typedefs_stat.h"
#include "../proto_defs/proto_main_local.h"
#include "../proto_defs/proto_math.h"
#include "../proto_defs/proto_integrate_md_entry.h"
#include "../proto_defs/proto_integrate_md_local.h"
#include "../proto_defs/proto_output_entry.h"
#include "../proto_defs/proto_analysis_md_entry.h"
#include "../proto_defs/proto_vel_sampl_class_entry.h"
#include "../proto_defs/proto_intra_con_entry.h"
#include "../proto_defs/proto_energy_ctrl_entry.h"
#include "../proto_defs/proto_communicate_wrappers.h"


/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void control_md(CLASS *class,BONDED *bonded,GENERAL_DATA *general_data,
                ANALYSIS *analysis)

/*=======================================================================*/
/*            Begin subprogram:                                          */
{   /*begin routine*/
/*=======================================================================*/
/*            Local variable declarations                                */
#include "../typ_defs/typ_mask.h"

  int itime,iii;        

  int ntime            = general_data->timeinfo.ntime;

  int ensopts_nve      = general_data->ensopts.nve;
  int ensopts_nvt      = general_data->ensopts.nvt;
  int ensopts_nvt_isok = general_data->ensopts.nvt_isok;
  int ensopts_npt_i    = general_data->ensopts.npt_i;
  int ensopts_npt_f    = general_data->ensopts.npt_f;

  int int_res_tra   = general_data->timeinfo.int_res_tra;
  int int_res_ter   = general_data->timeinfo.int_res_ter;

  int iwrite_screen = general_data->filenames.iwrite_screen;
  int iwrite_dump   = general_data->filenames.iwrite_dump;
  int iwrite_confp  = general_data->filenames.iwrite_confp;
  int iwrite_confv  = general_data->filenames.iwrite_confv;
  int iwrite_inst   = general_data->filenames.iwrite_inst;

  int ivel_smpl_on  = class->vel_samp_class.ivel_smpl_on;
  int nvx_smpl      = class->vel_samp_class.nvx_smpl;
  int nvnhc_smpl    = class->vel_samp_class.nvnhc_smpl;
  int iconstrnt     = bonded->constrnt.iconstrnt;

  int num_proc      = class->communicate.np;
  int myid          = class->communicate.myid;
  MPI_Comm world    = class->communicate.world;

/*=======================================================================*/
/* 0) Preliminary MD stuff                                               */

  prelim_md(class,bonded,general_data);

/*======================================================================*/
/* I) Write to Screen                                                   */

  if(myid==0){
    PRINT_LINE_STAR;
    printf("Running MD \n");
    PRINT_LINE_DASH;
  }/*endif*/

  (general_data->stat_avg.cpu_now) = 0.0;
  (general_data->stat_avg.acpu) += 0.0;
  simpavg_md(&(general_data->timeinfo),&(general_data->stat_avg),
               &(general_data->cell),&(bonded->constrnt),
	       &(general_data->ensopts),&(general_data->simopts),
	       &(general_data->ptens),&(class->communicate),
               &(class->nbr_list.verlist),&(class->energy_ctrl));

  output_md(class,general_data,bonded);

/*======================================================================*/
/* II) Loop over the specified number of time steps */
  
  for(itime = 1;itime<=ntime;itime++){
    
    if(num_proc>1){ Barrier(world);}
    cputime(&(general_data->stat_avg.cpu1)); 
    (general_data->timeinfo.itime) = itime;
    (class->energy_ctrl.itime)     = itime;
  
  /*---------------------------------------------------------------------*/
  /*   1)Do NVE dynamics                                                 */

    if(ensopts_nve==1){
      if((int_res_tra==0) && (int_res_ter==0)){
        int_NVE(class,bonded,general_data);
      }else{
        int_NVE_res(class,bonded,general_data);
      }/*endif*/
    }/*endif*/

  /*----------------------------------------------------------------------*/
  /*   2)Do NVT dynamics:                                                 */

    if(ensopts_nvt==1){
      if((int_res_tra==0) && (int_res_ter==0)){
        int_NVT(class,bonded,general_data);
      }else{
        int_NVT_res(class,bonded,general_data);
      }/*endif*/
    }/*endif*/

  /*----------------------------------------------------------------------*/
  /*   3) Do isokinetic NVT dynamics:                                                 */

    if(ensopts_nvt_isok==1){
      if((int_res_tra==0) && (int_res_ter==0)){
        int_NVT_ISOK(class,bonded,general_data);
      }else{
        int_NVT_ISOK_res(class,bonded,general_data);
      }/*endif*/
    }/*endif*/

  /*---------------------------------------------------------------------*/
  /*   4)Do isotropic npt dynamics:                                      */

    if(ensopts_npt_i==1){
      if((int_res_ter==0) && (int_res_tra==0)){
        int_NPTI(class,bonded,general_data);
      }else{
        int_NPTI_res(class,bonded,general_data);
      }/*endif*/
    }/*endif*/

  /*---------------------------------------------------------------------*/
  /*   5)Do flexible npt dynamics:                                       */
    
    if(ensopts_npt_f==1){
      if((int_res_ter==0) && (int_res_tra==0)){
        int_NPTF(class,bonded,general_data);
      }else{
        int_NPTF_res(class,bonded,general_data);
      }/*endif*/
    }/*endif*/

  /*----------------------------------------------------------------------*/
  /*   6)Do analysis md:                                                 */

    analysis_md(class,general_data,bonded,analysis); 

  /*----------------------------------------------------------------------*/
  /*   7)Calculate some simple averages                                   */

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
  /*   8)Produce the output specified by the user                          */

    if(myid == 0) check_auto_exit(&(general_data->timeinfo.exit_flag));
    if(num_proc > 1) Bcast(&(general_data->timeinfo.exit_flag),1,MPI_INT,0,world);
    

    if(  ((itime % iwrite_screen))==0 ||
         ((itime % iwrite_dump  ))==0 ||
         ((itime % iwrite_confp ))==0 ||
         ((itime % iwrite_confv ))==0 ||
         ((itime % iwrite_inst  ))==0 ||
         (general_data->timeinfo.exit_flag == 1) ){
         (general_data->filenames.ifile_open) = 0;
         output_md(class,general_data,bonded);
    }/*endif*/

  /*---------------------------------------------------------------------*/
  /*   9)Resample the particle velocities if desired                     */

    if( (ivel_smpl_on==1) && (nvx_smpl!=0) ){
      if( (itime % nvx_smpl)==0){
	control_vx_smpl(class,bonded,&(general_data->ensopts),
  		  &(general_data->simopts),general_data,
                  general_data->error_check_on);
        if(iconstrnt==1){
          init_constraint(bonded,&(general_data->ptens));
        }/*endif:init constraints*/
      }/*endif*/
    }/*endif*/

  /*---------------------------------------------------------------------*/
  /*   10) Resample the extended class velocities if desired             */

    if( (ivel_smpl_on==1) && (nvnhc_smpl!=0) ){
      if((itime % nvnhc_smpl)==0){
	control_vnhc_smpl(class,general_data);
      }/*endif*/
    }/*endif*/

  /*---------------------------------------------------------------------*/
  /*   11) Check for exit condition                                      */

    if(general_data->timeinfo.exit_flag == 1) itime = ntime;

  /*---------------------------------------------------------------------*/
  }/*endfor:itime */

/*======================================================================*/
/*  III)Write to Screen                                                 */
  
  if(myid==0){
    PRINT_LINE_DASH;
    printf("Completed MD run \n");
    PRINT_LINE_STAR;
  }/*endif*/

/*-----------------------------------------------------------------------*/
  }/*end routine*/
/*==========================================================================*/





/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void prelim_md(CLASS *class,BONDED *bonded,GENERAL_DATA *general_data)

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
    printf("Performing preliminary tasks for md \n");
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

  get_tvten(&(class->clatoms_info),&(class->clatoms_pos[1]),
            &(general_data->stat_avg),&(general_data->ptens),
            &(general_data->cell));

/*=======================================================================*/
/* IV)Get NHC Forces                                                     */

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
  if((general_data->ensopts.nvt_isok)==1){
	  init_NH_ISOK_par(general_data,&(class->clatoms_info),&(class->clatoms_pos[1]),
              &(class->therm_info_class),&(class->therm_class),
              &(class->int_scr),iflag_mass,&(class->class_comm_forc_pkg),&(class->vel_samp_class));
    iflag = 0;
	nhc_vol_potkin_isok(&(class->clatoms_info),&(class->clatoms_pos[1]),&(class->therm_info_class),&(class->therm_class),
            &(general_data->baro),&(general_data->par_rahman),
            &(general_data->stat_avg),&(general_data->statepoint),
            iflag,class->communicate.myid_forc,general_data->timeinfo.itime);
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
  if((general_data->ensopts.npt_i+general_data->ensopts.npt_f+
                                  general_data->ensopts.nvt)==0){
    general_data->stat_avg.kinet_nhc     = 0.0;
  }/*endif*/

/*========================================================================*/
/* V) Initialize free energy stuff                                       */

  if(bonded->bond_free.num != 0){
    for(i=1;i<= (bonded->bond_free.nhist);i++){
      bonded->bond_free.hist[i] = 0.0;
    }/*endfor*/
  }/*endif*/

  if(bonded->bend_free.num != 0){
    for(i=1;i<= (bonded->bend_free.nhist);i++){
      bonded->bend_free.hist[i] = 0.0;
    }/*endfor*/
  }/*endif*/

  if(bonded->tors_free.num == 1){
    for(i=1;i<= bonded->tors_free.nhist;i++){
      bonded->tors_free.hist[i] = 0.0;
    }/*endfor*/
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

/*=======================================================================*/
/* VI) Communicate initial energies                                      */

 if(num_proc>1){
    simpavg_md_communicate(&(general_data->stat_avg),&(general_data->ensopts),
        &(general_data->ptens),0,&(class->communicate),
        &(class->nbr_list.verlist),general_data->timeinfo.int_res_ter,1,1);
 }/*endif*/

/*=======================================================================*/
/*  VII) Write Energy to screen                                           */

  general_data->filenames.ifile_open = 1;
  general_data->timeinfo.itime       = 0;
  initial_output_md(class, general_data, bonded);

/*=======================================================================*/
/* VIII) Reinitailize update flag */

  general_data->stat_avg.updates   = 0.0;

/*=======================================================================*/
/* IX) Write to Screen         */

  if(myid==0){
    printf("\n");
    PRINT_LINE_DASH;
    printf("Completed preliminary tasks for md \n");
    PRINT_LINE_STAR;printf("\n");
  }/*endif*/

/*-----------------------------------------------------------------------*/
  }/*end routine*/
/*========================================================================*/

