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
void int_NVT_ISOK_res(CLASS *class,BONDED *bonded,GENERAL_DATA *general_data)
/*========================================================================*/
{/*begin routine*/
/*========================================================================*/
/*             Local variable declarations                                */

    int i,ipart,iflag;
    int ifirst                 = 1;
    int exit_flag              = 0;
    double dt,dti2,dti,tol_glob;
    int ir_tra,ir_tor,ir_ter;
    int ix_now;
    int nres_tra,nres_tor,nres_ter;
    int int_res_ter,iii;
    int natm_tot,ichain,inhc;
    int iflag_mass             = 1;
    int anneal_opt             = general_data->simopts.anneal_opt;
    double anneal_target_temp  = general_data->simopts.ann_target_temp;
    double ann_rate            = general_data->simopts.ann_rate;
    int num_nhc                = class->therm_info_class.num_nhc;
    int len_nhc                = class->therm_info_class.len_nhc;
    double **therm_v_nhc       = class->therm_class.v_nhc;
    double **therm_gkt         = class->therm_info_class.gkt;
    double **therm_mass_nhc    = class->therm_info_class.mass_nhc;
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
    int timestep               = general_data->timeinfo.itime;

/*==========================================================================*/
/* 0) Useful constants                                                      */

    natm_tot = class->clatoms_info.natm_tot;
    int_res_ter = general_data->timeinfo.int_res_ter;
    nres_tra = (general_data->timeinfo.nres_tra);
    nres_tor = (general_data->timeinfo.nres_tor);
    nres_ter = (general_data->timeinfo.nres_ter);
    dt   = (general_data->timeinfo.dt);
    dti  = dt/((double)(nres_ter*nres_tor*nres_tra));
    dti2 = dti/2.0;
    if(general_data->timeinfo.ix_respa==1){
      class->therm_info_class.dt_nhc   = dt;
    /*endif*/}
    if(general_data->timeinfo.ix_respa==2){
      class->therm_info_class.dt_nhc  = dt/((double)(nres_ter));
    /*endif*/}
    if(general_data->timeinfo.ix_respa==3){
      class->therm_info_class.dt_nhc  = dt/((double)(nres_ter*nres_tor));
    /*endif*/}
    if(general_data->timeinfo.ix_respa==4){
      class->therm_info_class.dt_nhc  = dt/((double)(
                                         nres_ter*nres_tor*nres_tra));
    /*endif*/}
    class->therm_info_class.dti_nhc  = (class->therm_info_class.dt_nhc)/
                              ( (double)(class->therm_info_class.nres_nhc) );
    set_yosh(class->therm_info_class.nyosh_nhc,
             class->therm_info_class.dti_nhc,class->therm_info_class.wdti,
             class->therm_info_class.wdti2,class->therm_info_class.wdti4,
             class->therm_info_class.wdti8,
             class->therm_info_class.wdti16);

    zero_constrt_iters(&(general_data->stat_avg));
    general_data->timeinfo.exit_flag = 0;
    
    energy_control(class,bonded,general_data);

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

    int_0_to_dt2_nvt_isok(class,bonded,general_data,ir_tra,ir_tor,ir_ter,dti);


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

          energy_control(class,bonded,general_data);

/*==========================================================================*/
/* 3) Evolve system from dt/2 to dt                                         */

        int_dt2_to_dt_nvt_isok(class,bonded,general_data,ir_tra,ir_tor,ir_ter,dti);


/*==========================================================================*/
        /*endfor:ir_tra*/}
      /*endfor:ir_tor*/}
    /*endfor:ir_ter*/}


/*==========================================================================*/
/* 4) Finalize                                                              */

   int_final_class(class,bonded,general_data,iflag);

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
      scanf("%d",&iii);
*/ 
/*-------------------------------------------------------------------------*/
/*end routine*/}
/*==========================================================================*/







