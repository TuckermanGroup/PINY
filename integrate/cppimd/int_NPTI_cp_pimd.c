/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
/*                                                                          */
/*                         PI_MD:                                           */
/*             The future of simulation technology                          */
/*             ------------------------------------                         */
/*                     Module: int_NVE                                      */
/*                                                                          */
/* This subprogram integrates the system using Vel Verlet                   */
/*                                                                          */
/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/


#include "standard_include.h"
#include "../typ_defs/typedefs_gen.h"
#include "../typ_defs/typedefs_cp.h"
#include "../typ_defs/typedefs_class.h"
#include "../typ_defs/typedefs_bnd.h"
#include "../proto_defs/proto_integrate_cppimd_entry.h"
#include "../proto_defs/proto_integrate_cppimd_local.h"
#include "../proto_defs/proto_integrate_cp_local.h"
#include "../proto_defs/proto_integrate_md_local.h"
#include "../proto_defs/proto_integrate_pimd_local.h"
#include "../proto_defs/proto_energy_cpcon_entry.h"
#include "../proto_defs/proto_energy_ctrl_cp_entry.h"
#include "../proto_defs/proto_intra_con_entry.h"
#include "../proto_defs/proto_pimd_entry.h"
#include "../proto_defs/proto_pimd_local.h"
#include "../proto_defs/proto_math.h"
#include "../proto_defs/proto_friend_lib_entry.h"
#include "../proto_defs/proto_communicate_wrappers.h"
#include "../proto_defs/proto_coords_local.h"
#define DEBUG_OFF

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
void int_NPTI_cp_pimd(CLASS *class,BONDED *bonded,
                             GENERAL_DATA *general_data,CP *cp)
/*========================================================================*/
{/*begin routine*/

#include "../typ_defs/typ_mask.h"

/*========================================================================*/
/*             Local variable declarations                                */
    int i,ipart,icoef,is,ifirst=1,ip,iflag,iflag_mass;
    double dt,dt2,tol_glob;
/*    double kinet_cp;*/
    int natm_tot,iii;
    int nstate_up,nstate_dn;
    int maxcoef_up,maxcoef_dn,ncoef_up,ncoef_up_tot;
    int ncoef_dn,ncoef_dn_tot;
    int ir_tra,ir_tor,ir_ter,ir_pimd;
    int pi_beads      = class->clatoms_info.pi_beads;
    int pi_beads_proc = class->clatoms_info.pi_beads_proc;
    double *class_clatoms_x;
    double *class_clatoms_y;
    double *class_clatoms_z;
    double *class_clatoms_vx;
    double *class_clatoms_vy;
    double *class_clatoms_vz;
    double *class_clatoms_fx;
    double *class_clatoms_fy;
    double *class_clatoms_fz;
    double *class_clatoms_mass;
    double *cp_cpcoeffs_cre_up;
    double *cp_cpcoeffs_cim_up;
    double *cp_cpcoeffs_cre_dn;
    double *cp_cpcoeffs_cim_dn;
    double *cp_cpcoeffs_vcre_up;
    double *cp_cpcoeffs_vcim_up;
    double *cp_cpcoeffs_vcre_dn;
    double *cp_cpcoeffs_vcim_dn;
    double *cp_cpcoeffs_fcre_up;
    double *cp_cpcoeffs_fcim_up;
    double *cp_cpcoeffs_fcre_dn;
    double *cp_cpcoeffs_fcim_dn;
    double *cp_cpcoeffs_cmass    = cp->cpcoeffs_info.cmass;
    double *cp_cpscr_cre_up  = cp->cpscr.cpscr_wave.cre_up;
    double *cp_cpscr_cim_up  = cp->cpscr.cpscr_wave.cim_up;
    double *cp_cpscr_cre_dn  = cp->cpscr.cpscr_wave.cre_dn;
    double *cp_cpscr_cim_dn  = cp->cpscr.cpscr_wave.cim_dn;
    double aa,aa2,arg2,poly,bb,dlen;
    double e2,e4,e6,e8;
    int *cpcoeffs_ioff_up       = cp->cpcoeffs_info.ioff_up;
    int *cpcoeffs_ioff_dn       = cp->cpcoeffs_info.ioff_dn;
    int massiv_flag             = cp->cptherm_info.massiv_flag;
    int myid                    = class->communicate.myid;
    int myid_bead               = class->communicate.myid_bead;
    int myid_state              = class->communicate.myid_state;
    int np_states               = class->communicate.np_states;
    MPI_Comm world              = cp->communicate.world;
    MPI_Comm comm_beads         = cp->communicate.comm_beads;

    e2=1.0/(2.0*3.0);e4=e2/(4.0*5.0);
    e6= e4/(6.0*7.0);e8=e6/(8.0*9.0);


/*==========================================================================*/
/* 0) Useful constants                                                      */

    natm_tot = class->clatoms_info.natm_tot;
    general_data->timeinfo.int_res_tra    = 0;
    general_data->timeinfo.int_res_ter    = 0;
    dt  = general_data->timeinfo.dt;
    dt2 = general_data->timeinfo.dt/2.0;

    ir_tra  = 1;
    ir_tor  = 1;
    ir_ter  = 1;
    ir_pimd = 1;


    zero_constrt_iters(&(general_data->stat_avg));

    class->therm_info_class.wght  =  1.0;
    class->therm_info_bead.wght  =  1.0;
    class->therm_info_class.dt_nhc   =  dt;
    class->therm_info_bead.dt_nhc   =  dt;
    class->therm_info_class.dti_nhc  = (class->therm_info_class.dt_nhc)/
                             ( (double)(class->therm_info_class.nres_nhc) );
    class->therm_info_bead.dti_nhc  = (class->therm_info_bead.dt_nhc)/
                             ( (double)(class->therm_info_bead.nres_nhc) );

    nstate_up    = cp->cpcoeffs_info.nstate_up_proc;
    nstate_dn    = cp->cpcoeffs_info.nstate_dn_proc;
    ncoef_up     = cp->cpcoeffs_info.ncoef;
    ncoef_up_tot = cp->cpcoeffs_info.ncoef*nstate_up;
    ncoef_dn     = (cp->cpopts.cp_lsda == 1 ? cp->cpcoeffs_info.ncoef:0);
    ncoef_dn_tot = (cp->cpopts.cp_lsda == 1 ? 
                    cp->cpcoeffs_info.ncoef*nstate_dn:0);
    maxcoef_up = ncoef_up;
    maxcoef_dn = ncoef_dn;

    class->therm_info_class.dt_nhc   =  dt;
    class->therm_info_bead.dt_nhc   =  dt;
    class->therm_info_class.dti_nhc  = (class->therm_info_class.dt_nhc)/
                             ( (double)(class->therm_info_class.nres_nhc) );
    class->therm_info_bead.dti_nhc  = (class->therm_info_bead.dt_nhc)/
                             ( (double)(class->therm_info_bead.nres_nhc) );

  if(myid_state==0){
    set_yosh(class->therm_info_class.nyosh_nhc,
             class->therm_info_class.dti_nhc,class->therm_info_class.wdti,
             class->therm_info_class.wdti2,class->therm_info_class.wdti4,
             class->therm_info_class.wdti8,
             class->therm_info_class.wdti16);
    set_yosh(class->therm_info_bead.nyosh_nhc,
             class->therm_info_bead.dti_nhc,class->therm_info_bead.wdti,
             class->therm_info_bead.wdti2,class->therm_info_bead.wdti4,
             class->therm_info_bead.wdti8,
             class->therm_info_bead.wdti16);
  }/*endif*/
    if((cp->cptherm_info.num_c_nhc > 0) && (cp->cptherm_info.len_c_nhc > 0)) {
      cp->cptherm_info.dt_nhc  = dt;
      cp->cptherm_info.dti_nhc = dt/( (double)(cp->cptherm_info.nres_c_nhc) );
      set_yosh(cp->cptherm_info.nyosh_c_nhc,cp->cptherm_info.dti_nhc ,
               cp->cptherm_info.wdti,cp->cptherm_info.wdti2,
               cp->cptherm_info.wdti4,cp->cptherm_info.wdti8,
               cp->cptherm_info.wdti16);
    }/*endif:therms on*/

/*==========================================================================*/
/* I) Evolve electronic part of system  from t=0 to dt/2                    */

    int_0_to_dt2_cp(class,bonded,general_data,cp);


/*==========================================================================*/
/* II) Evolve classical part of system  from t=0 to dt/2                    */

    if(myid_state==0){
      int_0_to_dt2_npti_pimd(class,bonded,general_data,ir_tra,ir_tor,
                                                       ir_ter,ir_pimd,dt);
    }/*endif*/

    if(np_states > 1){
      comm_coord_class_state(class,general_data,1,1,0,1);
    }/* endif */

    if(class->communicate.np_states>1){
       Bcast(&(general_data->cell.hmat[0]),10,MPI_DOUBLE,0,                                                        class->communicate.comm_states);
    }/*endif*/

/*==========================================================================*/
/* III) Get the new energy/force                                            */

    (class->energy_ctrl.iget_full_inter)= 1;
    (class->energy_ctrl.iget_res_inter) = 0;
    (class->energy_ctrl.iget_full_intra)= 1;
    (class->energy_ctrl.iget_res_intra) = 0;

    cp_energy_control_pimd(class,bonded,general_data,cp);

/*==========================================================================*/
/* IV) Evolve classical part of system  from t=dt/2 to dt                   */

    if(myid_state==0){
      int_dt2_to_dt_npti_pimd(class,bonded,general_data,ir_tra,ir_tor,
                                                    ir_ter,ir_pimd,dt); 
    }

/*==========================================================================*/
/* V) Evolve electronic part of system  from t=dt/2 to dt                   */

  int_dt2_to_dt_cp(class,bonded,general_data,cp);

/*==========================================================================*/
/* VI) Get NHC contribution to energy                                       */

    iflag = 1; 
    if(myid_state==0){
      nhc_vol_potkin_pimd(class,general_data,iflag);
    }/*endif*/

/*==========================================================================*/
/*
      printf("x(1),y(1),z(1) %.13g %.13g %.13g\n",
                                         class->clatoms_pos[1].x[1],
                                         class->clatoms_pos[1].y[1],
                                         class->clatoms_pos[1].z[1]); 
      printf("vx(1),vy(1),vz(1) %.13g %.13g %.13g\n",
                                         class->clatoms_pos[1].vx[1],
                                         class->clatoms_pos[1].vy[1],
                                         class->clatoms_pos[1].vz[1]); 
      printf("fx(1),fy(1),fz(1) %.13g %.13g %.13g\n",
                                         class->clatoms_pos[1].fx[1],
                                         class->clatoms_pos[1].fy[1],
                                         class->clatoms_pos[1].fz[1]); 
      scanf("%d",&iii);
*/
/*-------------------------------------------------------------------------*/
/*end routine*/}
/*==========================================================================*/














