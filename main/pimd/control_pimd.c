/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
/*                                                                          */
/*                         PI_MD:                                           */
/*             The future of simulation technology                          */
/*             ------------------------------------                         */
/*                     Module: control_md                                   */
/*                                                                          */
/* This subprogram performs Path Integral MD on a classical potential energy*/ 
/*                             surface (PES)                                */
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
#include "../proto_defs/proto_integrate_pimd_entry.h"
#include "../proto_defs/proto_integrate_pimd_local.h"
#include "../proto_defs/proto_integrate_md_entry.h"
#include "../proto_defs/proto_output_entry.h"
#include "../proto_defs/proto_analysis_md_entry.h"
#include "../proto_defs/proto_vel_sampl_class_entry.h"
#include "../proto_defs/proto_energy_ctrl_entry.h"
#include "../proto_defs/proto_pimd_entry.h"
#include "../proto_defs/proto_communicate_wrappers.h"


/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void control_pimd(CLASS *class,BONDED *bonded,GENERAL_DATA *general_data,
                  ANALYSIS *analysis)

/*=======================================================================*/
/*            Begin subprogram:                                          */
{   /*begin routine*/

#include "../typ_defs/typ_mask.h"

/*=======================================================================*/
/*            Local variable declarations                                */

  int itime,iii,nver_temp=0;        
  int ntime           = general_data->timeinfo.ntime;
  int error_check_on  = general_data->error_check_on;
  int num_proc        = class->communicate.np;
  int np_beads        = class->communicate.np_beads;
  int np_forc         = class->communicate.np_forc;
  int myid_forc       = class->communicate.myid_forc;
  int myid            = class->communicate.myid;
  MPI_Comm comm_beads = class->communicate.comm_beads;
  MPI_Comm comm_forc  = class->communicate.comm_forc;
  MPI_Comm world      = class->communicate.world;

/*=======================================================================*/
/* 0) Preliminary MD stuff                                               */


  prelim_pimd(class,bonded,general_data);
  if(num_proc>1){Barrier(world);}

/*======================================================================*/
/* I) Write to Screen                                                   */

if(error_check_on==1){
  PRINT_LINE_STAR;
  printf("Running PIMD \n");
  PRINT_LINE_DASH;
}/*endif for error check on*/
  
/*======================================================================*/
/* II) Loop over the specified number of time steps */

  general_data->stat_avg.updates = 0.0;
  (class->energy_ctrl.itime)     = 0; /* Always get the total PE */

  for(itime = 1; itime<=ntime; itime++){

    cputime(&(general_data->stat_avg.cpu1)); 
    (general_data->timeinfo.itime) = itime;
  /*----------------------------------------------------------------------*/
  /*   1)Do NVT dynamics:                                                 */
    if((general_data->ensopts.nvt)==1){
       int_NVT_pimd(class,bonded,general_data,analysis);
    }/*endif*/
  /*---------------------------------------------------------------------*/
  /*   2)Do isotropic npt dynamics:                                      */
    if((general_data->ensopts.npt_i)==1){
       int_NPTI_pimd(class,bonded,general_data,analysis);
    }/*endif*/

  /*---------------------------------------------------------------------*/
  /*   3)Do flexible npt dynamics:                                       */
    if((general_data->ensopts.npt_f)==1){
       int_NPTF_pimd(class,bonded,general_data,analysis);
    }/*endif*/
  /*----------------------------------------------------------------------*/
  /*   4)Do nst dynamics:                                                 */

  /*---------------------------------------------------------------------*/
  /*  4.5) Analysis Routine                                               */
#ifdef JUNK
    analysis_pimd(class,general_data,bonded,analysis); 
#endif
  /*----------------------------------------------------------------------*/
  /*   5)Calculate some simple averages                                   */
    cputime(&(general_data->stat_avg.cpu2)); 
    (general_data->stat_avg.cpu_now)=(general_data->stat_avg.cpu2)-
                                     (general_data->stat_avg.cpu1);
    (general_data->stat_avg.acpu) += (general_data->stat_avg.cpu2)-
                                     (general_data->stat_avg.cpu1);

    simpavg_pimd(&(general_data->timeinfo),&(general_data->stat_avg),
	         &(general_data->cell),&(bonded->constrnt),
	         &(general_data->ensopts),&(general_data->simopts),
	         &(general_data->ptens),&(class->communicate));

  /*-----------------------------------------------------------------------*/
  /*   7)Produce the output specified by the user                          */

    if(np_forc>1){
      Reduce(&(class->nbr_list.verlist.nver_lst_now), &nver_temp,1,MPI_INT,
                   MPI_SUM,0,comm_forc);
      class->nbr_list.verlist.nver_lst_now = nver_temp;     
    }

/*=======================================================================*/
/* I) All gather velocities if necessary */ 

 if(np_forc>1){
  if(((general_data->timeinfo.itime%general_data->filenames.iwrite_confv)==0)||
     ((general_data->timeinfo.itime%general_data->filenames.iwrite_dump)==0)){
    forc_level_vel_gather(class);
  }/*endif*/
 }/*endif*/



    if(myid == 0) check_auto_exit(&(general_data->timeinfo.exit_flag));
    if(num_proc > 1) Bcast(&(general_data->timeinfo.exit_flag),1,MPI_INT,0,world);

    if(  (itime % (general_data->filenames.iwrite_screen))==0 ||
         (itime % (general_data->filenames.iwrite_dump  ))==0 ||
         (itime % (general_data->filenames.iwrite_confp ))==0 ||
         (itime % (general_data->filenames.iwrite_confv)) ==0 ||
         (itime % (general_data->filenames.iwrite_inst))  ==0 ||
         (general_data->timeinfo.exit_flag == 1)       ){
         (general_data->filenames.ifile_open) = 0;

         if(np_beads>1&&myid_forc==0){
           communicate_output_pimd(class);
	 }/*endif*/
         if(error_check_on==1){
           output_pimd(class,general_data,bonded);
	 }/*endif : error_check_on*/
    }/*endif*/

  /*---------------------------------------------------------------------*/
  /*   8)Resample the particle velocities if desired                     */
    if(((class->vel_samp_class.ivel_smpl_on)==1) &&
       ((class->vel_samp_class.nvx_smpl)!=0)){
      if((itime % (class->vel_samp_class.nvx_smpl))==0){
	control_vx_smpl(class,bonded,&(general_data->ensopts),
			&(general_data->simopts),general_data,
                          general_data->error_check_on);
      }/*endif*/
    }/*endif*/
  /*---------------------------------------------------------------------*/
  /*   9) Resample the extended class velocities if desired             */
    if(((class->vel_samp_class.ivel_smpl_on)==1) &&
       ((class->vel_samp_class.nvnhc_smpl)!=0)){
      if((itime % (class->vel_samp_class.nvnhc_smpl))==0){
	control_vnhc_smpl(class,general_data);
      }/*endif*/
    }/*endif*/
  /*---------------------------------------------------------------------*/
  /*   10) Check for exit condition                                      */

   if(general_data->timeinfo.exit_flag == 1) {
      itime = ntime;
      if(np_beads > 1) Barrier(world);
   }/* endif */

  }/*endfor:itime */


  /*======================================================================*/
  /*  III)Write to Screen                                                 */
if(error_check_on==1){
  PRINT_LINE_DASH;
  printf("Completed PIMD run \n");
  PRINT_LINE_STAR;
}/*endif for error check on*/

/*-----------------------------------------------------------------------*/
}/*end routine*/
/*==========================================================================*/





/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void prelim_pimd(CLASS *class,BONDED *bonded,GENERAL_DATA *general_data)

/*=======================================================================*/
/*            Begin subprogram:                                          */
{   /*begin routine*/

#include "../typ_defs/typ_mask.h"

/*=======================================================================*/
/*            Local variable declarations                                */

  int i,j,iflag,ip,iflag_mass,iii;        
  int pi_beads = class->clatoms_info.pi_beads;
  int pi_beads_proc = class->clatoms_info.pi_beads_proc;
  int pi_beads_proc_st = class->clatoms_info.pi_beads_proc_st;
  int rank      = class->communicate.myid_bead;
  int myid_forc = class->communicate.myid_forc;
  int error_check_on = general_data->error_check_on;
  int nver_temp = 0;
  int mytherm_start     = class->therm_info_class.mytherm_start;
  int mytherm_end       = class->therm_info_class.mytherm_end;
  double kinet_temp=0.0;
  double vintrat_temp=0.0;
  double vintert_temp=0.0;
  double vvdw_temp = 0.0; /*DY*/
  double vcoul_temp = 0.0; /*DY*/
  double kinet_nhc_temp=0.0;
  double kinet_nhc_bead_temp=0.0;
  double vbondt_temp=0.0;
  double vbendt_temp=0.0;
  double vbend_bnd_bond_temp=0.0;
  double vbend_bnd_bend_temp=0.0;
  double vtorst_temp=0.0;
  double vonfot_temp = 0.0;
  double pi_ke_prim_temp = 0.0;
  double pi_ke_vir_temp = 0.0;
  double kin_harm_temp = 0.0;
  int start_proc;
  MPI_Comm comm_forc  = class->communicate.comm_forc;
  MPI_Comm comm_beads = class->communicate.comm_beads;
  MPI_Comm world      = class->communicate.world;
  int num_proc        = class->communicate.np;
  int num_proc_beads  = class->communicate.np_beads;
  int np_forc         = class->communicate.np_forc;
           start_proc = (pi_beads_proc_st == 1 ? 2:1);


/*=======================================================================*/
/*   I) Write to Screen                                                  */

 if(error_check_on==1){
   PRINT_LINE_STAR;
   printf("Performing preliminary tasks for pimd \n");
   PRINT_LINE_DASH;printf("\n");
 }/*endif*/

/*=======================================================================*/
/*  II) Initialize In-output                                             */

  general_data->stat_avg.updates   = 0.0;
  general_data->stat_avg.acpu      = 0.0;

/*=======================================================================*/
/* III)Get Energy and particle Forces                                    */

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
  if(error_check_on==1){
   printf("Getting initial energy\n");
  }/*endif*/
 
  class->clatoms_info.wght_pimd = 1;

    control_pimd_trans_mode(class,general_data);
    mode_energy_control(class,general_data);
    control_pimd_trans_pos(class,general_data);

  energy_control_pimd(class,bonded,general_data);

  get_tvten_pimd(class,general_data);

/*=======================================================================*/
/* IV)Get NHC Forces                                                     */

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
  if((general_data->ensopts.npt_i+general_data->ensopts.npt_f+
                                  general_data->ensopts.nvt)==0){
    general_data->stat_avg.kinet_nhc     = 0.0;
  }/*endif*/

/*========================================================================*/
/* VI) Initialize free energy stuff                                       */

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



/*=======================================================================*/
/*  IV) Write Energy to screen                                           */
/*=======================================================================*/

  general_data->filenames.ifile_open = 1;

  general_data->timeinfo.itime = 0;
  
  /*=======================================================================*/
  /*  Collect initial energies (needs its own routine)                     */
  /*=======================================================================*/

  if(num_proc>1){
    Reduce(&(class->nbr_list.verlist.nver_lst_now), &nver_temp,1,MPI_INT,
                   MPI_SUM,0,comm_forc);
    Reduce(&(general_data->stat_avg.vintert), &vintert_temp,1,MPI_DOUBLE,
                   MPI_SUM,0,world);
    Reduce(&(general_data->stat_avg.vintrat), &vintrat_temp,1,MPI_DOUBLE,
                   MPI_SUM,0,world);
    Reduce(&(general_data->stat_avg.vvdw), &vvdw_temp,1,MPI_DOUBLE,
                   MPI_SUM,0,world);
    Reduce(&(general_data->stat_avg.vcoul),&vcoul_temp,1,MPI_DOUBLE,
                   MPI_SUM,0,world);
    Reduce(&(general_data->stat_avg.kinet), &kinet_temp,1,MPI_DOUBLE,
                   MPI_SUM,0,world);
    Reduce(&(general_data->stat_avg.vbondt), &vbondt_temp,1,
                   MPI_DOUBLE,MPI_SUM,0,world);
    Reduce(&(general_data->stat_avg.vbendt), &vbendt_temp,1,
                   MPI_DOUBLE,MPI_SUM,0,world);
    Reduce(&(general_data->stat_avg.vbend_bnd_bond), &vbend_bnd_bond_temp,1,
                   MPI_DOUBLE,MPI_SUM,0,world);
    Reduce(&(general_data->stat_avg.vbend_bnd_bend), &vbend_bnd_bend_temp,1,
                   MPI_DOUBLE,MPI_SUM,0,world);
    Reduce(&(general_data->stat_avg.vtorst), &vtorst_temp,1,
                   MPI_DOUBLE,MPI_SUM,0,world);
    Reduce(&(general_data->stat_avg.vonfot), &vonfot_temp,1,
                   MPI_DOUBLE,MPI_SUM,0,world);
    Reduce(&(general_data->stat_avg.pi_ke_vir), &pi_ke_vir_temp,1,
                   MPI_DOUBLE,MPI_SUM,0,world);
    Reduce(&(general_data->stat_avg.kin_harm), &kin_harm_temp,1,
                   MPI_DOUBLE,MPI_SUM,0,world);

 

    general_data->stat_avg.kinet          = kinet_temp;
    general_data->stat_avg.vintert        = vintert_temp;
    general_data->stat_avg.vintrat        = vintrat_temp;
    general_data->stat_avg.vvdw           = vvdw_temp;  
    general_data->stat_avg.vcoul          = vcoul_temp; 
    general_data->stat_avg.vbondt         = vbondt_temp;
    general_data->stat_avg.vbendt         = vbendt_temp;
    general_data->stat_avg.vbend_bnd_bond = vbend_bnd_bond_temp;
    general_data->stat_avg.vbend_bnd_bend = vbend_bnd_bend_temp;
    general_data->stat_avg.vtorst         = vtorst_temp;
    general_data->stat_avg.vonfot         = vonfot_temp;
    general_data->stat_avg.pi_ke_vir      = pi_ke_vir_temp;
    general_data->stat_avg.kin_harm       = kin_harm_temp;
    class->nbr_list.verlist.nver_lst_now  = nver_temp;


    Reduce(&(general_data->stat_avg.kinet_nhc_bead), &kinet_nhc_bead_temp,1,
                   MPI_DOUBLE,MPI_SUM,0,world); 
    general_data->stat_avg.kinet_nhc_bead = kinet_nhc_bead_temp; 

    Reduce(&(general_data->stat_avg.kinet_nhc), &kinet_nhc_temp,1,
                   MPI_DOUBLE,MPI_SUM,0,world);
    general_data->stat_avg.kinet_nhc = kinet_nhc_temp;
 
  }/*endif*/


  if(num_proc_beads>1){
    Reduce(&(general_data->stat_avg.pi_ke_prim), &pi_ke_prim_temp,1,
                   MPI_DOUBLE,MPI_SUM,0,comm_beads);
    general_data->stat_avg.pi_ke_prim     = pi_ke_prim_temp;
  }/*endif*/

  if(error_check_on==1){
   output_pimd(class,general_data,bonded);
  }/*endif*/


/*=======================================================================*/
/*   VII) Write to Screen         */

  if(error_check_on==1){
   printf("\n");
   PRINT_LINE_DASH;
   printf("Completed preliminary tasks for pi_md \n");
   PRINT_LINE_STAR;printf("\n");
  }/*endif for error check on*/

/*-----------------------------------------------------------------------*/
   }/*end routine*/
/*==========================================================================*/









