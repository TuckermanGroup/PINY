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
#include "../proto_defs/proto_integrate_cp_entry.h"
#include "../proto_defs/proto_integrate_cp_local.h"
#include "../proto_defs/proto_integrate_md_local.h"
#include "../proto_defs/proto_energy_cpcon_entry.h"
#include "../proto_defs/proto_energy_ctrl_cp_entry.h"
#include "../proto_defs/proto_intra_con_entry.h"
#include "../proto_defs/proto_friend_lib_entry.h"
#include "../proto_defs/proto_communicate_wrappers.h"
#include "../proto_defs/proto_coords_local.h"


void int_0_to_dt2_cp(CLASS *,BONDED *,GENERAL_DATA *,CP *);
void int_dt2_to_dt_cp(CLASS *,BONDED *,GENERAL_DATA *,CP *);


#define DEBUG_OFF
#define DEBUG_NHC_OFF

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void int_NVE_cp(CLASS *class,BONDED *bonded,GENERAL_DATA *general_data,CP *cp)

/*========================================================================*/
{/*begin routine*/
/*========================================================================*/
/*             Local variable declarations                                */
#include "../typ_defs/typ_mask.h"

    int i,ip=1,ipart,icoef,is,ifirst=1,iproc;
    int exit_flag=0;
    int iflag,iflag_mass=1;
    double dt,dt2,tol_glob;
    int natm_tot,iii,jjj;
    int nstate_up,nstate_dn;
    int ncoef,ncoef_up_tot,ncoef_dn_tot;
    int ir_tra=1,ir_ter=1,ir_tor=1;
    int anneal_opt  = general_data->simopts.anneal_opt;
    double anneal_target_temp = general_data->simopts.ann_target_temp;
    double ann_rate = general_data->simopts.ann_rate;
    double *clatoms_xold = class->clatoms_info.xold;
    double *clatoms_yold = class->clatoms_info.yold;
    double *clatoms_zold = class->clatoms_info.zold;
    double *clatoms_x    = class->clatoms_pos[1].x;
    double *clatoms_y    = class->clatoms_pos[1].y;
    double *clatoms_z    = class->clatoms_pos[1].z;
    double *clatoms_vx   = class->clatoms_pos[1].vx;
    double *clatoms_vy   = class->clatoms_pos[1].vy;
    double *clatoms_vz   = class->clatoms_pos[1].vz;
    double *clatoms_fx   = class->clatoms_pos[1].fx;
    double *clatoms_fy   = class->clatoms_pos[1].fy;
    double *clatoms_fz   = class->clatoms_pos[1].fz;
    double *clatoms_mass = class->clatoms_info.mass;


    int myid_state           = class->communicate.myid_state;
    int np_states            = class->communicate.np_states;
    int np_forc              = class->communicate.np_forc;
    MPI_Comm comm_states     = class->communicate.comm_states;


/*==========================================================================*/

/*==========================================================================*/
/* 0) Useful constants                                                      */

    natm_tot = class->clatoms_info.natm_tot;
    general_data->timeinfo.int_res_tra    = 0;
    general_data->timeinfo.int_res_ter    = 0;
    dt  = general_data->timeinfo.dt;
    dt2 = general_data->timeinfo.dt/2.0;

    zero_constrt_iters(&(general_data->stat_avg));
    general_data->timeinfo.exit_flag = 0;

    if((cp->cptherm_info.num_c_nhc > 0) && (cp->cptherm_info.len_c_nhc > 0)) {
      cp->cptherm_info.dt_nhc  = dt;
      cp->cptherm_info.dti_nhc = dt/( (double)(cp->cptherm_info.nres_c_nhc) );
      set_yosh(cp->cptherm_info.nyosh_c_nhc,cp->cptherm_info.dti_nhc ,
               cp->cptherm_info.wdti,cp->cptherm_info.wdti2,
               cp->cptherm_info.wdti4,cp->cptherm_info.wdti8,
               cp->cptherm_info.wdti16);
    }/*endif:therms on*/


/*==========================================================================*/
/* I) perform cp integration from 0 to dt/2                                 */

  int_0_to_dt2_cp(class,bonded,general_data,cp);

/*==========================================================================*/
/* II) Atm integration from 0 to dt/2                                       */


  if( (myid_state == 0) || (np_forc > 1) ){
    if(general_data->simopts.cp==1){
      int_0_to_dt2_nve(class,bonded,general_data,ir_tra,ir_tor,ir_ter,dt);
    }/*endif*/
  }/* endif myid */


/*==========================================================================*/
/* III) Communicate atom positions to state processors                       */

 if(np_states > 1 && np_forc == 1){
   comm_coord_class_state(class,general_data,1,1,0,0);
 }/* endif */

/*==========================================================================*/
/* IV) Get the new energy/force                                              */

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
/* V) Atoms:  The other half                                                */

 if((myid_state == 0) || (np_forc > 1) ){
   if(general_data->simopts.cp == 1){
    int_dt2_to_dt_nve(class,bonded,general_data,ir_tra,ir_tor,ir_ter,dt);
   }/* endif */
 }/* endif */

/*==========================================================================*/
/* VI) CP:  The other half                                                  */

  int_dt2_to_dt_cp(class,bonded,general_data,cp);


/*==========================================================================*/
/* VII) Scale by annealing factor                                            */

 if(myid_state == 0 || np_forc > 1){
  iflag=-1;
  if(anneal_opt == 1){
     anneal_class(class,ann_rate,iflag,iflag_mass,anneal_target_temp,&exit_flag);
     general_data->timeinfo.exit_flag = exit_flag;
  }/*endif*/
 }/* endif */

/*==========================================================================*/
/* VIII) Finalize                                                              */

  if(myid_state == 0 || np_forc > 1){
   int_final_class(class,bonded,general_data,iflag);
  }/* endif */

/*==========================================================================*/
/* X) DEBUG                                                                 */

#ifdef DEBUG

 if(myid_state == 0){
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
 if(cp->communicate.myid==0){
 printf("coeff nhc vels %g %g\n",cp->cptherm_pos[1].vc_nhc[1][1],
                             cp->cptherm_pos[1].vc_nhc[2][1]);
 printf("coeffs %g %g\n",cp->cpcoeffs_pos[1].cre_up[1],
                             cp->cpcoeffs_pos[1].cim_up[1]);
 printf("coeff vels %g %g\n",cp->cpcoeffs_pos[1].vcre_up[1],
                             cp->cpcoeffs_pos[1].vcim_up[1]);
 printf("coeff forces %g %g\n",cp->cpcoeffs_pos[1].fcre_up[1],
                             cp->cpcoeffs_pos[1].fcim_up[1]);
 printf("coeff nhc's %g %g\n",cp->cptherm_pos[1].c_nhc[1][1],
                             cp->cptherm_pos[1].c_nhc[2][1]);
 printf("coeff nhc forces %g %g\n",cp->cptherm_pos[1].fc_nhc[1][1],
                             cp->cptherm_pos[1].fc_nhc[2][1]);
 scanf("%d",&iii);
 }
 }/*endif*/
 Dbx_Barrier(comm_states);

#endif

/*-------------------------------------------------------------------------*/
     }/*end routine*/
/*==========================================================================*/






