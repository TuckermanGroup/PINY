/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
/*                                                                          */
/*                         PI_MD:                                           */
/*             The future of simulation technology                          */
/*             ------------------------------------                         */
/*                     Module: int_NVT                                      */
/*                                                                          */
/* This subprogram integrates the system using Vel Verlet                   */
/*                                                                          */
/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

#include "standard_include.h"
#include "../typ_defs/typedefs_gen.h"
#include "../typ_defs/typedefs_class.h"
#include "../typ_defs/typedefs_cp.h"
#include "../typ_defs/typedefs_bnd.h"
#include "../proto_defs/proto_integrate_cp_entry.h"
#include "../proto_defs/proto_integrate_cp_local.h"
#include "../proto_defs/proto_integrate_md_entry.h"
#include "../proto_defs/proto_integrate_md_local.h"
#include "../proto_defs/proto_energy_cpcon_entry.h"
#include "../proto_defs/proto_energy_ctrl_cp_entry.h"
#include "../proto_defs/proto_intra_con_entry.h"
#include "../proto_defs/proto_energy_ctrl_entry.h"
#include "../proto_defs/proto_communicate_wrappers.h"
#include "../proto_defs/proto_coords_local.h"

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
void int_NVT_cp(CLASS *class,BONDED *bonded,GENERAL_DATA *general_data,CP *cp)
/*========================================================================*/
{/*begin routine*/
/*========================================================================*/
/*             Local variable declarations                                */ 
    int iii;
    int i,ip=1,ipart,icoef,is,iflag,ifirst=1;
    int exit_flag=0;
    int natm_tot,ichain,inhc;
    double dt,dt2,tol_glob;
    int iflag_mass = 1;


    int ir_tra = 1;
    int ir_tor = 1;
    int ir_ter = 1;
    int anneal_opt = general_data->simopts.anneal_opt;
    double anneal_target_temp = general_data->simopts.ann_target_temp;
    double ann_rate = general_data->simopts.ann_rate;
    int     num_nhc = class->therm_info_class.num_nhc;
    int     len_nhc = class->therm_info_class.len_nhc;
    double **therm_v_nhc    = class->therm_class.v_nhc;
    double **therm_gkt      = class->therm_info_class.gkt;
    double **therm_mass_nhc = class->therm_info_class.mass_nhc;
    double *class_clatoms_xold = class->clatoms_info.xold;
    double *class_clatoms_yold = class->clatoms_info.yold;
    double *class_clatoms_zold = class->clatoms_info.zold;
    double *class_clatoms_x    = class->clatoms_pos[1].x;
    double *class_clatoms_y    = class->clatoms_pos[1].y;
    double *class_clatoms_z    = class->clatoms_pos[1].z;
    double *class_clatoms_vx   = class->clatoms_pos[1].vx;
    double *class_clatoms_vy   = class->clatoms_pos[1].vy;
    double *class_clatoms_vz   = class->clatoms_pos[1].vz;
    double *class_clatoms_fx   = class->clatoms_pos[1].fx;
    double *class_clatoms_fy   = class->clatoms_pos[1].fy;
    double *class_clatoms_fz   = class->clatoms_pos[1].fz;
    double *class_clatoms_mass = class->clatoms_info.mass;

    int myid_state           = class->communicate.myid_state;
    int np_states            = class->communicate.np_states;
    int np_forc              = class->communicate.np_forc;
    MPI_Comm comm_states     = class->communicate.comm_states;

/*==========================================================================*/
/* 0) Useful constants                                                      */

    general_data->timeinfo.exit_flag = 0;
    natm_tot = class->clatoms_info.natm_tot;
    (general_data->timeinfo.int_res_tra)    = 0;
    (general_data->timeinfo.int_res_ter)    = 0;
    dt  = (general_data->timeinfo.dt);
    dt2 = (general_data->timeinfo.dt)/2.0;
    class->therm_info_class.dt_nhc  = dt;
    class->therm_info_class.dti_nhc = dt/( (double)(
                             class->therm_info_class.nres_nhc) );
  if(myid_state==0 || np_forc > 1){
    set_yosh(class->therm_info_class.nyosh_nhc,
             class->therm_info_class.dti_nhc ,class->therm_info_class.wdti,
             class->therm_info_class.wdti2,class->therm_info_class.wdti4,
             class->therm_info_class.wdti8,
             class->therm_info_class.wdti16);
  }/*endif*/
    general_data->stat_avg.iter_shake    = 0; 
    general_data->stat_avg.iter_ratl     = 0; 


    general_data->stat_avg.iter_23r = 0.0;
    general_data->stat_avg.iter_33r = 0.0;
    general_data->stat_avg.iter_46r = 0.0;
    general_data->stat_avg.iter_23 = 0.0;
    general_data->stat_avg.iter_33 = 0.0;
    general_data->stat_avg.iter_46 = 0.0;



/*==========================================================================*/
/* I) CP:  The first half                                                   */

  int_0_to_dt2_cp(class,bonded,general_data,cp);


/*==========================================================================*/
/* II) Atoms: The first half                                                */

  if( (myid_state==0) || (np_forc > 1) ){
   if(general_data->simopts.cp == 1){
    int_0_to_dt2_nvt(class,bonded,general_data,ir_tra,ir_tor,ir_ter,dt);
   }/*endif*/
  }/*endif*/


/*==========================================================================*/
/* III) Atoms: The first half                                                */

  if(np_states > 1 && np_forc==1 ){
    comm_coord_class_state(class,general_data,1,1,0,0);
  }/* endif */

/*==========================================================================*/
/* IV) Get the new energy/force                                             */

    (class->energy_ctrl.iget_full_inter)= 1;
    (class->energy_ctrl.iget_res_inter) = 0;
    (class->energy_ctrl.iget_full_intra)= 1;
    (class->energy_ctrl.iget_res_intra) = 0;

     cp_energy_control(class,bonded,general_data,cp);

    if(cp->cpopts.cp_gauss == 1){
     add_gauss_force(&(cp->cpcoeffs_info),&(cp->cpcoeffs_pos[1]),
		     &(cp->cpscr.cpscr_ovmat),&(cp->cpopts));
    }/*endif*/

/*==========================================================================*/
/* V) Atoms: Second half                                                    */


  if((myid_state == 0) || (np_forc > 1) ){
    if(general_data->simopts.cp==1){
     int_dt2_to_dt_nvt(class,bonded,general_data,ir_tra,ir_tor,ir_ter,dt);
   }/*endif : cp_on*/
  }/*endif : myid=0*/


/*==========================================================================*/
/* VI) CP: second half                                                      */

  int_dt2_to_dt_cp(class,bonded,general_data,cp);


/*==========================================================================*/
/* VII) Scale by annealing factor                                          */

if(myid_state == 0 || np_forc > 1){
  if(anneal_opt == 1){
    iflag=0;
    anneal_class(class,ann_rate,iflag,iflag_mass,anneal_target_temp,&exit_flag);
    general_data->timeinfo.exit_flag = exit_flag;
  }/*endif:ann_rate*/
}/*endif:myid=0*/

/*==========================================================================*/
/* VIII) Finalize                                                              */

  if(myid_state == 0 || np_forc > 1){
   int_final_class(class,bonded,general_data,iflag);
  }/* endif */

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
      printf("v_nhc[1][1],v_nhc[2,1] %.13g %.13g\n",
                                         class->therm_class.v_nhc[1][1],
                                         class->therm_class.v_nhc[2][1]); 
      printf("x_nhc[1][1],x_nhc[2][1] %.13g %.13g\n",
                                         class->therm_class.x_nhc[1][1],
                                         class->therm_class.x_nhc[2][1]); 
      printf("coeffs %.13g %.13g\n",
           cp->cpcoeffs_pos[1].cre_up[1],cp->cpcoeffs_pos[1].cim_up[1]);
      printf("coef vels %.13g %.13g\n",
           cp->cpcoeffs_pos[1].vcre_up[1],cp->cpcoeffs_pos[1].vcim_up[1]);
      printf("coef forces %.13g %.13g\n",
           cp->cpcoeffs_pos[1].fcre_up[1],cp->cpcoeffs_pos[1].fcim_up[1]);
      scanf("%d",&iii);
*/
/*--------------------------------------------------------------------------*/
/*end routine*/}
/*==========================================================================*/


