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
#include "../typ_defs/typedefs_gen.h"
#include "../typ_defs/typedefs_cp.h"
#include "../typ_defs/typedefs_class.h"
#include "../typ_defs/typedefs_bnd.h"
#include "../typ_defs/typedefs_stat.h"
#include "../proto_defs/proto_main_cp_local.h"
#include "../proto_defs/proto_math.h"
#include "../proto_defs/proto_output_cp_entry.h"
#include "../proto_defs/proto_output_cp_local.h"
#include "../proto_defs/proto_intra_con_entry.h"
#include "../proto_defs/proto_real_space_local.h"
#include "../proto_defs/proto_integrate_pimd_local.h"
#include "../proto_defs/proto_energy_ctrl_cp_entry.h"
#include "../proto_defs/proto_integrate_md_entry.h"
#include "../proto_defs/proto_integrate_pimd_entry.h"
#include "../proto_defs/proto_integrate_cp_local.h"
#include "../proto_defs/proto_integrate_cppimd_local.h"
#include "../proto_defs/proto_communicate_wrappers.h"

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void control_debug_cp_pimd(CLASS *class,BONDED *bonded,
                           GENERAL_DATA *general_data,CP *cp)

/*=======================================================================*/
/*            Begin subprogram:                                          */
{   /*begin routine*/
/*=======================================================================*/
/*            Local variable declarations                                */
#include "../typ_defs/typ_mask.h"
  int iii;
  int i,j,ip,iflag,start,start_proc;        
  int iflag_mass = 1;
  int npairs,np_tot;         /* Num: # pairs, Max # pairs         */
  int pi_beads      = class->clatoms_info.pi_beads;
  int pi_beads_proc = class->clatoms_info.pi_beads_proc;
  double cpu1,cpu2,cpu; 
  double kinet_temp = 0.0;
  double kinet_nhc_temp = 0.0;
  double kinet_nhc_bead_temp = 0.0;
  double kinet_cp_temp = 0.0;
  double cp_ehart_temp = 0.0;
  double cp_exc_temp = 0.0;
  double cp_muxc_temp = 0.0;
  double cp_eext_temp = 0.0;
  double cp_enl_temp = 0.0;
  double cp_eke_temp = 0.0;
  int error_check_on = general_data->error_check_on;
  int rank = class->communicate.myid;
  int myid_bead = class->communicate.myid_bead;
  int np_beads = class->communicate.np_beads;
  int num_proc = class->communicate.np;
  MPI_Comm world      = class->communicate.world;
  MPI_Comm comm_beads = class->communicate.comm_beads;
  MPI_Comm comm_states = class->communicate.comm_states;

  int pi_beads_proc_st = class->clatoms_info.pi_beads_proc_st;
           start_proc = (pi_beads_proc_st == 1 ? 2:1);


/*=======================================================================*/
/*   I) Write to Screen                                                  */

  if(error_check_on==1){
    printf("\n");PRINT_LINE_STAR;
    printf("Debugging CP portion of pi_md \n");
    PRINT_LINE_DASH;printf("\n");
  }
  if(num_proc>1){Barrier(world);}

/*=======================================================================*/
/* II) Hurl CP Scalars and Rolf CP Vectors */
   
#ifdef DEBUG
  if(error_check_on==1){
   printf("Do you want to hurl and rolf (1/0)");scanf("%d",&iii);
   if(iii==1){
     hurl_cp_scalars(cp);
     rolf_cp_vectors(cp,class->atommaps.natm_typ,
                        class->clatoms_info.natm_tot);
   }/*endif*/
  }/*endif:error checking*/
  if(num_proc>1){Barrier(world);}
#endif

/*=======================================================================*/
/* III) Initialize */

  start = (rank == 0 ? 2:1);

  general_data->filenames.ifile_open = 1;

  general_data->timeinfo.itime = 0;
  class->energy_ctrl.itime     = 0;
  
  if(error_check_on==1){
    output_cp_pimd(class,general_data,bonded,cp);
  }
  if(num_proc>1){Barrier(world);}

  general_data->stat_avg.updates   = 0.0;
  general_data->stat_avg.acpu      = 0.0; 

/*========================================================================*/
/* IV) Initialize free energy stuff                                       */

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
/* V) Set constraints                                                  */

  if(bonded->constrnt.iconstrnt==1){
     init_constraint(bonded,&(general_data->ptens));
  }/*endif:init constraints, must be done before getting the energy*/

    general_data->stat_avg.iter_shake   = 0;
    general_data->stat_avg.iter_ratl    = 0;
    general_data->stat_avg.iter_23    = 0;
    general_data->stat_avg.iter_33    = 0;
    general_data->stat_avg.iter_46    = 0;

/*=======================================================================*/
/* VI) Set energy flags                                                  */

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
/*=======================================================================*/
/* VII) Test neighbor list update                                         */

  if((class->nbr_list.iver)==1){
    if(error_check_on==1){
     printf("======================================\n");
     if((class->nbr_list.verlist.nolst_ver_update)==1){
        printf("Checking verlist update type: no_lst\n");
     }else{
        printf("Checking verlist update type: lnk_lst\n");
     }/*endif*/ 
     printf("--------------------------------------\n\n");
    }
    if(num_proc>1){Barrier(world);}
    cputime(&cpu1);
      verlist_control(&(class->clatoms_info),&(class->clatoms_pos[1]),
                  &(class->nbr_list),
                  &(bonded->excl),&(class->atommaps),&(general_data->cell),
                  &(bonded->intra_scr),&(class->for_scr),
                  &(general_data->timeinfo),&(class->interact),error_check_on,
                  &(class->class_comm_forc_pkg));
    cputime(&cpu2);
    cpu = cpu2-cpu1;
    npairs = class->nbr_list.verlist.nter[(class->clatoms_info.natm_tot)] 
            + class->nbr_list.verlist.jver_off[(class->clatoms_info.natm_tot)];
    np_tot = (class->clatoms_info.natm_tot)*(class->clatoms_info.natm_tot-1)
             /2 - (bonded->excl.nlst);

    if(error_check_on==1){
     printf("Number of pairs = %d out of %d; cpu time %g\n",npairs,np_tot,cpu);
     printf("\n--------------------------------------\n");
     printf("Finished checking list update         \n");
     printf("======================================\n\n");
     printf("Enter an integer to continue:  ");
     scanf("%d",&iii); printf("\n");
    }
    if(num_proc>1){Barrier(world);}
  }/*endif*/

/*=======================================================================*/
/* VIII) Test energy                                                        */

  test_energy_cp_pimd(class,bonded,general_data,cp);
  get_tvten_pimd(class,general_data);
  get_cpke_pimd(&(cp->cpcoeffs_info),(cp->cpcoeffs_pos),
                &(general_data->stat_avg),cp->cpopts.cp_lsda,
                class->clatoms_info.pi_beads_proc,
                class->communicate.np_states);

/*=======================================================================*/
/* IX)Get Atm NHC Forces                                                 */

if(myid_bead<np_beads){

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

}/*endif : bead processors*/
if(num_proc>1){Barrier(comm_beads);}

  if((general_data->ensopts.npt_i+general_data->ensopts.npt_f)==0){
    general_data->stat_avg.kinet_v       = 0.0;
  }/*endif*/
  if((general_data->ensopts.npt_i+general_data->ensopts.npt_f+
                                  general_data->ensopts.nvt)==0){
    general_data->stat_avg.kinet_nhc     = 0.0;
  }/*endif*/

/*=======================================================================*/
/* X)Get CP NHC Forces                                                   */

  if(cp->cptherm_info.num_c_nhc > 0){
   for(ip=1;ip<=pi_beads_proc;ip++){
    init_cp_NHC(&(cp->cptherm_info),&(cp->cptherm_pos[ip]),
                &(cp->cpcoeffs_info),&(cp->cpcoeffs_pos[ip]),
                &(cp->cpscr),cp->cpopts.cp_lsda,comm_states,
                cp->communicate.np_states,cp->communicate.myid_state); 
   }/*endfor*/
    nhc_cp_potkin_pimd(&(cp->cptherm_info),(cp->cptherm_pos),
                       &(general_data->stat_avg),pi_beads);
  }/* endif */

/*=======================================================================*/
/*  Communication for initial output                                    */

if(num_proc>1){

   Reduce(&(general_data->stat_avg.cp_ehart), &cp_ehart_temp,1,MPI_DOUBLE,
                   MPI_SUM,0,world);
   Reduce(&(general_data->stat_avg.cp_exc), &cp_exc_temp,1,MPI_DOUBLE,
                   MPI_SUM,0,world);
   Reduce(&(general_data->stat_avg.cp_muxc), &cp_muxc_temp,1,MPI_DOUBLE,
                   MPI_SUM,0,world);
   Reduce(&(general_data->stat_avg.cp_eext), &cp_eext_temp,1,MPI_DOUBLE,
                   MPI_SUM,0,world);
   Reduce(&(general_data->stat_avg.cp_enl), &cp_enl_temp,1,MPI_DOUBLE,
                   MPI_SUM,0,world);
   Reduce(&(general_data->stat_avg.cp_eke), &cp_eke_temp,1,MPI_DOUBLE,
                   MPI_SUM,0,world);
   Reduce(&(general_data->stat_avg.kinet_cp), &kinet_cp_temp,1,
                   MPI_DOUBLE,MPI_SUM,0,world);

   if(myid_bead<np_beads){
     Reduce(&(general_data->stat_avg.kinet), &kinet_temp,1,MPI_DOUBLE,
                   MPI_SUM,0,class->communicate.comm_beads);
     Reduce(&(general_data->stat_avg.kinet_nhc_bead), &kinet_nhc_bead_temp,1,
                   MPI_DOUBLE,MPI_SUM,0,class->communicate.comm_beads);
     Reduce(&(general_data->stat_avg.kinet_nhc), &kinet_nhc_temp,1,
                   MPI_DOUBLE,MPI_SUM,0,class->communicate.comm_beads);

     general_data->stat_avg.kinet = kinet_temp;
     general_data->stat_avg.kinet_nhc_bead = kinet_nhc_bead_temp;
     general_data->stat_avg.kinet_nhc = kinet_nhc_temp;
   }/*endif : bead processors*/
   if(num_proc>1){Barrier(comm_beads);}

   general_data->stat_avg.kinet_cp = kinet_cp_temp;
   general_data->stat_avg.cp_ehart = cp_ehart_temp;
   general_data->stat_avg.cp_exc = cp_exc_temp;
   general_data->stat_avg.cp_muxc = cp_muxc_temp;
   general_data->stat_avg.cp_eext = cp_eext_temp;
   general_data->stat_avg.cp_enl = cp_enl_temp;
   general_data->stat_avg.cp_eke = cp_eke_temp;

 }/*endif*/

/*=======================================================================*/
/*  XI) Write Energy to screen                                           */

  if(error_check_on==1){
    initial_output_cp_pimd(class,general_data,bonded,cp);
  }
  if(num_proc>1){Barrier(world);}

/*=======================================================================*/
/* XII) Write to Screen         */

  if(error_check_on==1){
    printf("\n");
    PRINT_LINE_DASH;
    printf("Completed debug\n");
    PRINT_LINE_STAR;printf("\n");
  }
  if(num_proc>1){Barrier(world);}

/*-----------------------------------------------------------------------*/
   }/*end routine*/
/*==========================================================================*/





















