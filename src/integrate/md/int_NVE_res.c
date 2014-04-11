/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
/*                                                                          */
/*                         PI_MD:                                           */
/*             The future of simulation technology                          */
/*             ------------------------------------                         */
/*                     Module: int_NVE_res                                  */
/*                                                                          */
/* This subprogram integrates the system using Vel Verlet RESPA             */
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
void int_NVE_res(CLASS *class,BONDED *bonded,GENERAL_DATA *general_data)
/*========================================================================*/
{/*begin routine*/
/*========================================================================*/
/*             Local variable declarations                                */

#include "../typ_defs/typ_mask.h"

    int i,ipart,iii,ifirst=1,ip,iflag,iflag_mass=1,n,iproc;
    int exit_flag=0;
    double dt,dti2,dti,tol_glob;
    int ir_tra,ir_tor,ir_ter;
    int nres_tra,nres_tor,nres_ter;
    int int_res_ter;
    int natm_tot,istart,iend,nsec,nsec_now;
    int anneal_opt = general_data->simopts.anneal_opt;
    double anneal_target_temp = general_data->simopts.ann_target_temp;
    double ann_rate = general_data->simopts.ann_rate;
    double *xold = class->clatoms_info.xold;
    double *yold = class->clatoms_info.yold;
    double *zold = class->clatoms_info.zold;
    double *x    = class->clatoms_pos[1].x;
    double *y    = class->clatoms_pos[1].y;
    double *z    = class->clatoms_pos[1].z;
    double *vx   = class->clatoms_pos[1].vx;
    double *vy   = class->clatoms_pos[1].vy;
    double *vz   = class->clatoms_pos[1].vz;
    double *fx   = class->clatoms_pos[1].fx;
    double *fy   = class->clatoms_pos[1].fy;
    double *fz   = class->clatoms_pos[1].fz;
    double *mass = class->clatoms_info.mass;
    int np_forc                = class->communicate.np_forc;
    int myid_forc              = class->communicate.myid_forc;
    int pi_beads_proc          = class->clatoms_info.pi_beads_proc;
    MPI_Comm comm_forc         = class->communicate.comm_forc;

/*==========================================================================*/
/* 0) Useful constants                                                      */

    n        = class->clatoms_info.natm_tot;
    natm_tot = class->clatoms_info.natm_tot;
    int_res_ter = general_data->timeinfo.int_res_ter;
    nres_tra = general_data->timeinfo.nres_tra;
    nres_tor = general_data->timeinfo.nres_tor;
    nres_ter = general_data->timeinfo.nres_ter;
    dt   = general_data->timeinfo.dt;
    dti  =  dt/((double)(nres_ter*nres_tra*nres_tor));
    dti2 = dti/2.0;

    zero_constrt_iters(&(general_data->stat_avg));
    general_data->timeinfo.exit_flag = 0;

/*==========================================================================*/
/*==========================================================================*/
/* I) Loop over inter RESPA                                                 */

    for(ir_ter=1;ir_ter<=nres_ter;ir_ter++){

/*==========================================================================*/
/*==========================================================================*/
/* II) Loop over tors RESPA                                                 */

      for(ir_tor=1;ir_tor<=nres_tor;ir_tor++){

/*==========================================================================*/
/*==========================================================================*/
/* III) Loop over intra RESPA                                               */

        for(ir_tra=1;ir_tra<=nres_tra;ir_tra++){

/*==========================================================================*/
/* 1) Evolve system from t=0 to dt/2                                        */

    int_0_to_dt2_nve(class,bonded,general_data,ir_tra,ir_tor,ir_ter,dti);

/*==========================================================================*/
/* 2) Get the new energy/force                                              */

          (class->energy_ctrl.iget_full_inter) = 0;
          (class->energy_ctrl.iget_res_inter) = 0;
          (class->energy_ctrl.iget_full_intra) = 0;
          (class->energy_ctrl.iget_res_intra)  = 1;
          if((ir_ter==nres_ter)&&(ir_tor==nres_tor)&&(ir_tra==nres_tra))
             {(class->energy_ctrl.iget_full_inter) = 1;}
          if((int_res_ter==1)&&(ir_tor==nres_tor)&&(ir_tra==nres_tra))
             {(class->energy_ctrl.iget_res_inter) = 1;}
          if(ir_tra==nres_tra)
             {(class->energy_ctrl.iget_full_intra) = 1;}

          energy_control(class, bonded, general_data);

/*==========================================================================*/
/* 3) Evolve system from dt/2 to dt                                         */

        int_dt2_to_dt_nve(class,bonded,general_data,ir_tra,ir_tor,ir_ter,dti);

/*==========================================================================*/

        /*endfor:ir_tra*/}
      /*endfor:ir_tor*/}
    /*endfor:ir_ter*/}

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
















