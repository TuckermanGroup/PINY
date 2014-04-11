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
#include "../typ_defs/typedefs_class.h"
#include "../typ_defs/typedefs_bnd.h"
#include "../typ_defs/typedefs_gen.h"
#include "../proto_defs/proto_integrate_md_entry.h"
#include "../proto_defs/proto_integrate_md_local.h"
#include "../proto_defs/proto_intra_con_entry.h"
#include "../proto_defs/proto_energy_ctrl_entry.h"

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
void int_NVE(CLASS *class,BONDED *bonded,GENERAL_DATA *general_data)
/*========================================================================*/
{/*begin routine*/
/*========================================================================*/
/*             Local variable declarations                                */

#include "../typ_defs/typ_mask.h"

    int i,ipart,ifirst=1,ip,iflag,iflag_mass=1;
    int exit_flag=0;
    double dt,dt2,tol_glob;
    int natm_tot,iii,istart,iend,nsec,nsec_now,n,iproc;
    int ir_tra = 1;
    int ir_tor = 1;
    int ir_ter = 1;
    int anneal_opt = general_data->simopts.anneal_opt;
    double ann_rate = general_data->simopts.ann_rate;
    double anneal_target_temp = general_data->simopts.ann_target_temp;
    double *x       = class->clatoms_pos[1].x;
    double *y       = class->clatoms_pos[1].y;
    double *z       = class->clatoms_pos[1].z;
    double *vx      = class->clatoms_pos[1].vx;
    double *vy      = class->clatoms_pos[1].vy;
    double *vz      = class->clatoms_pos[1].vz;
    double *fx      = class->clatoms_pos[1].fx;
    double *fy      = class->clatoms_pos[1].fy;
    double *fz       = class->clatoms_pos[1].fz;
    int np_forc     = class->communicate.np_forc;
    int myid_forc   = class->communicate.myid_forc;
    int pi_beads_proc          = class->clatoms_info.pi_beads_proc;
    MPI_Comm comm_forc         = class->communicate.comm_forc;

/*==========================================================================*/
/* 0) Useful constants                                                      */

    natm_tot = class->clatoms_info.natm_tot;
    n        = natm_tot;
    general_data->timeinfo.int_res_tra    = 0;
    general_data->timeinfo.int_res_ter    = 0;
    dt  = general_data->timeinfo.dt;
    dt2 = general_data->timeinfo.dt/2.0;

    zero_constrt_iters(&(general_data->stat_avg));
    general_data->timeinfo.exit_flag = 0;

/*==========================================================================*/
/* 1) Evolve system from t=0 to dt/2                                        */

    int_0_to_dt2_nve(class,bonded,general_data,ir_tra,ir_tor,ir_ter,dt);

/*==========================================================================*/
/* 2) Get the new energy/force                                              */

    (class->energy_ctrl.iget_full_inter)= 1;
    (class->energy_ctrl.iget_res_inter) = 0;
    (class->energy_ctrl.iget_full_intra)= 1;
    (class->energy_ctrl.iget_res_intra) = 0;

     energy_control(class, bonded, general_data);

/*==========================================================================*/
/* 3) Evolve system from dt/2 to dt                                         */

    int_dt2_to_dt_nve(class,bonded,general_data,ir_tra,ir_tor,ir_ter,dt);

/*==========================================================================*/
/* 4) Scale by annealing factor                                             */

   iflag=-1;
   if(anneal_opt == 1){
      anneal_class(class,ann_rate,iflag,iflag_mass,anneal_target_temp,&exit_flag);
      general_data->timeinfo.exit_flag = exit_flag;
   }/*endif*/

/*==========================================================================*/
/* 5) Finalize                                                              */

   int_final_class(class,bonded,general_data,iflag);

/*==========================================================================*/

#ifdef DEBUG_GLENN
    for(iproc=0;iproc<np_forc;iproc++){
      Barrier(comm_forc);
      if(myid_forc==iproc){
       printf("1 %.12g %.12g %.12g %.12g %.12g %.12g %.12g %.12g %.12g\n",
                  x[1],y[1],z[1],vx[1],vy[1],vz[1],fx[1],fy[1],fz[1]);
       printf("n %.12g %.12g %.12g %.12g %.12g %.12g %.12g %.12g %.12g\n",
                  x[n],y[n],z[n],vx[n],vy[n],vz[n],fx[n],fy[n],fz[n]);
      }
    }
    if(myid_forc==0){scanf("%d",&iii);}
    Barrier(comm_forc);
#endif

/*-------------------------------------------------------------------------*/
/*end routine*/}
/*==========================================================================*/

















