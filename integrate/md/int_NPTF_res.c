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
#include "../proto_defs/proto_math.h"

#define JUNK_OFF

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void int_NPTF_res(CLASS *class,BONDED *bonded,GENERAL_DATA *general_data)

/*========================================================================*/
   {/*begin routine*/
/*========================================================================*/
/*             Local variable declarations                                */

   int iflag;
   double dti2,dti;
   int ir_tra,ir_tor,ir_ter;
   int iii;

   int ix_respa    = general_data->timeinfo.ix_respa;
   int int_res_ter = general_data->timeinfo.int_res_ter;
   int nres_tra    = general_data->timeinfo.nres_tra;
   int nres_tor    = general_data->timeinfo.nres_tor;
   int nres_ter    = general_data->timeinfo.nres_ter;
   int nres_nhc    = class->therm_info_class.nres_nhc;
   double dt       = general_data->timeinfo.dt;

/*==========================================================================*/
/* 0) Useful constants                                                      */

   dti  = dt/((double)(nres_ter*nres_tor*nres_tra));
   dti2 = dti/2.0;

   switch(ix_respa){
     case 1:
       class->therm_info_class.wght  =  (double)(nres_tra*nres_tor*nres_ter);
       break;
     case 2:
       class->therm_info_class.wght  = (double)(nres_tra*nres_tor);
       break;
     case 3:
       class->therm_info_class.wght  =  (double)(nres_tra);
       break;
     case 4:
       class->therm_info_class.wght  =  1.0;
       break;
   }/*end switch*/

   class->therm_info_class.dt_nhc   =  dti;
   class->therm_info_class.dti_nhc  =  dti/( (double)(nres_nhc) );

   set_yosh(class->therm_info_class.nyosh_nhc,
            class->therm_info_class.dti_nhc,
            class->therm_info_class.wdti,
            class->therm_info_class.wdti2,
            class->therm_info_class.wdti4,
            class->therm_info_class.wdti8,
            class->therm_info_class.wdti16);

   zero_constrt_iters(&(general_data->stat_avg));

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

#ifdef JUNK
      printf("in NPTF: %d %d %d\n",ir_ter,ir_tor,ir_tra);
      printf("x(1),y(1),z(1) %.13g %.13g %.13g\n",
                                         class->clatoms_pos[1].x[1],
                                         class->clatoms_pos[1].y[1],
                                         class->clatoms_pos[1].z[1]); 
      printf("vx(1),vy(1),vz(1) %.13g %.13g %.13g\n",
                                            class->clatoms_pos[1].vx[1],
                                            class->clatoms_pos[1].vy[1],
                                            class->clatoms_pos[1].vz[1]); 
      printf("v_nhc[1][1],v_nhc[2,1] %.13g %.13g\n",
                                  class->therm_class.v_nhc[1][1],
                                  class->therm_class.v_nhc[2][1]); 
      printf("x_nhc[1][1],x_nhc[2][1] %.13g %.13g\n",
                                  class->therm_class.x_nhc[1][1],
                                  class->therm_class.x_nhc[2][1]); 
      printf("v_vol_nhc[1],v_vol_nhc(2) %.13g %.13g\n",
                                   general_data->baro.v_vol_nhc[1],
                                   general_data->baro.v_vol_nhc[2]); 
      printf("x_vol_nhc[1],x_vol_nhc(2) %.13g %.13g\n",
                                   general_data->baro.x_vol_nhc[1],
                                   general_data->baro.x_vol_nhc[2]); 
      printf("vgmat  %.13g %.13g %.13g\n",general_data->par_rahman.vgmat[1],
                                 general_data->par_rahman.vgmat[2],
                                 general_data->par_rahman.vgmat[3]);
      printf("vgmat  %.13g %.13g %.13g\n",general_data->par_rahman.vgmat[4],
                                 general_data->par_rahman.vgmat[5],
                                 general_data->par_rahman.vgmat[6]);
      printf("vgmat  %.13g %.13g %.13g\n",general_data->par_rahman.vgmat[7],
                                 general_data->par_rahman.vgmat[8],
                                 general_data->par_rahman.vgmat[9]);
      printf("hmat %.13g %.13g %.13g\n",general_data->cell.hmat[1],
                                        general_data->cell.hmat[2],
                                        general_data->cell.hmat[3]);
      printf("hmat %.13g %.13g %.13g\n",general_data->cell.hmat[4],
                                        general_data->cell.hmat[5],
                                        general_data->cell.hmat[6]);
      printf("hmat %.13g %.13g %.13g\n",general_data->cell.hmat[7],
                                        general_data->cell.hmat[8],
                                        general_data->cell.hmat[9]);
      scanf("%d",&iii);
#endif
        int_0_to_dt2_nptf(class,bonded,general_data,ir_tra,ir_tor,ir_ter,dti);

#ifdef JUNK
      printf("in NPTF, after 0_to_dt2_nptf %d %d %d\n",ir_ter,ir_tor,ir_tra);
      printf("x(1),y(1),z(1) %.13g %.13g %.13g\n",
                                         class->clatoms_pos[1].x[1],
                                         class->clatoms_pos[1].y[1],
                                         class->clatoms_pos[1].z[1]); 
      printf("vx(1),vy(1),vz(1) %.13g %.13g %.13g\n",
                                            class->clatoms_pos[1].vx[1],
                                            class->clatoms_pos[1].vy[1],
                                            class->clatoms_pos[1].vz[1]); 
      printf("v_nhc[1][1],v_nhc[2,1] %.13g %.13g\n",
                                  class->therm_class.v_nhc[1][1],
                                  class->therm_class.v_nhc[2][1]); 
      printf("x_nhc[1][1],x_nhc[2][1] %.13g %.13g\n",
                                  class->therm_class.x_nhc[1][1],
                                  class->therm_class.x_nhc[2][1]); 
      printf("v_vol_nhc[1],v_vol_nhc(2) %.13g %.13g\n",
                                   general_data->baro.v_vol_nhc[1],
                                   general_data->baro.v_vol_nhc[2]); 
      printf("x_vol_nhc[1],x_vol_nhc(2) %.13g %.13g\n",
                                   general_data->baro.x_vol_nhc[1],
                                   general_data->baro.x_vol_nhc[2]); 
      printf("vgmat  %.13g %.13g %.13g\n",general_data->par_rahman.vgmat[1],
                                 general_data->par_rahman.vgmat[2],
                                 general_data->par_rahman.vgmat[3]);
      printf("vgmat  %.13g %.13g %.13g\n",general_data->par_rahman.vgmat[4],
                                 general_data->par_rahman.vgmat[5],
                                 general_data->par_rahman.vgmat[6]);
      printf("vgmat  %.13g %.13g %.13g\n",general_data->par_rahman.vgmat[7],
                                 general_data->par_rahman.vgmat[8],
                                 general_data->par_rahman.vgmat[9]);
      printf("hmat %.13g %.13g %.13g\n",general_data->cell.hmat[1],
                                        general_data->cell.hmat[2],
                                        general_data->cell.hmat[3]);
      printf("hmat %.13g %.13g %.13g\n",general_data->cell.hmat[4],
                                        general_data->cell.hmat[5],
                                        general_data->cell.hmat[6]);
      printf("hmat %.13g %.13g %.13g\n",general_data->cell.hmat[7],
                                        general_data->cell.hmat[8],
                                        general_data->cell.hmat[9]);
      scanf("%d",&iii);
#endif
/*==========================================================================*/
/* 2) Get the new energy/force                                              */

        (class->energy_ctrl.iget_full_inter) = 0;
        (class->energy_ctrl.iget_res_inter)  = 0;
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

       int_dt2_to_dt_nptf(class,bonded,general_data,ir_tra,ir_tor,ir_ter,dti);

#ifdef JUNK
      printf("after dt2_to_dt %d %d %d\n",ir_ter,ir_tor,ir_tra);
      printf("x(1),y(1),z(1) %.13g %.13g %.13g\n",
                                         class->clatoms_pos[1].x[1],
                                         class->clatoms_pos[1].y[1],
                                         class->clatoms_pos[1].z[1]); 
      printf("vx(1),vy(1),vz(1) %.13g %.13g %.13g\n",
                                            class->clatoms_pos[1].vx[1],
                                            class->clatoms_pos[1].vy[1],
                                            class->clatoms_pos[1].vz[1]); 
      printf("v_nhc[1][1],v_nhc[2,1] %.13g %.13g\n",
                                  class->therm_class.v_nhc[1][1],
                                  class->therm_class.v_nhc[2][1]); 
      printf("x_nhc[1][1],x_nhc[2][1] %.13g %.13g\n",
                                  class->therm_class.x_nhc[1][1],
                                  class->therm_class.x_nhc[2][1]); 
      printf("v_vol_nhc[1],v_vol_nhc(2) %.13g %.13g\n",
                                   general_data->baro.v_vol_nhc[1],
                                   general_data->baro.v_vol_nhc[2]); 
      printf("x_vol_nhc[1],x_vol_nhc(2) %.13g %.13g\n",
                                   general_data->baro.x_vol_nhc[1],
                                   general_data->baro.x_vol_nhc[2]); 
      printf("vgmat  %.13g %.13g %.13g\n",general_data->par_rahman.vgmat[1],
                                 general_data->par_rahman.vgmat[2],
                                 general_data->par_rahman.vgmat[3]);
      printf("vgmat  %.13g %.13g %.13g\n",general_data->par_rahman.vgmat[4],
                                 general_data->par_rahman.vgmat[5],
                                 general_data->par_rahman.vgmat[6]);
      printf("vgmat  %.13g %.13g %.13g\n",general_data->par_rahman.vgmat[7],
                                 general_data->par_rahman.vgmat[8],
                                 general_data->par_rahman.vgmat[9]);
      printf("hmat %.13g %.13g %.13g\n",general_data->cell.hmat[1],
                                        general_data->cell.hmat[2],
                                        general_data->cell.hmat[3]);
      printf("hmat %.13g %.13g %.13g\n",general_data->cell.hmat[4],
                                        general_data->cell.hmat[5],
                                        general_data->cell.hmat[6]);
      printf("hmat %.13g %.13g %.13g\n",general_data->cell.hmat[7],
                                        general_data->cell.hmat[8],
                                        general_data->cell.hmat[9]);
      scanf("%d",&iii);
#endif

/*==========================================================================*/
/*==========================================================================*/

        }/*endfor:ir_tra*/
     }/*endfor:ir_tor*/
   }/*endfor:ir_ter*/

/*==========================================================================*/
/* 4) Finalize                                                              */

   iflag=2;
   int_final_class(class,bonded,general_data,iflag);

/*==========================================================================*/
#ifdef DEBUG
      printf("x(1),y(1),z(1) %.13g %.13g %.13g\n",
                                         class->clatoms_pos[1].x[1],
                                         class->clatoms_pos[1].y[1],
                                         class->clatoms_pos[1].z[1]); 
      printf("vx(1),vy(1),vz(1) %.13g %.13g %.13g\n",
                                            class->clatoms_pos[1].vx[1],
                                            class->clatoms_pos[1].vy[1],
                                            class->clatoms_pos[1].vz[1]); 
      printf("v_nhc[1][1],v_nhc[2,1] %.13g %.13g\n",
                                  class->therm_class.v_nhc[1][1],
                                  class->therm_class.v_nhc[2][1]); 
      printf("x_nhc[1][1],x_nhc[2][1] %.13g %.13g\n",
                                  class->therm_class.x_nhc[1][1],
                                  class->therm_class.x_nhc[2][1]); 
      printf("v_vol_nhc[1],v_vol_nhc(2) %.13g %.13g\n",
                                   general_data->baro.v_vol_nhc[1],
                                   general_data->baro.v_vol_nhc[2]); 
      printf("x_vol_nhc[1],x_vol_nhc(2) %.13g %.13g\n",
                                   general_data->baro.x_vol_nhc[1],
                                   general_data->baro.x_vol_nhc[2]); 
      printf("vgmat  %.13g %.13g %.13g\n",general_data->par_rahman.vgmat[1],
                                 general_data->par_rahman.vgmat[2],
                                 general_data->par_rahman.vgmat[3]);
      printf("vgmat  %.13g %.13g %.13g\n",general_data->par_rahman.vgmat[4],
                                 general_data->par_rahman.vgmat[5],
                                 general_data->par_rahman.vgmat[6]);
      printf("vgmat  %.13g %.13g %.13g\n",general_data->par_rahman.vgmat[7],
                                 general_data->par_rahman.vgmat[8],
                                 general_data->par_rahman.vgmat[9]);
      printf("hmat %.13g %.13g %.13g\n",general_data->cell.hmat[1],
                                        general_data->cell.hmat[2],
                                        general_data->cell.hmat[3]);
      printf("hmat %.13g %.13g %.13g\n",general_data->cell.hmat[4],
                                        general_data->cell.hmat[5],
                                        general_data->cell.hmat[6]);
      printf("hmat %.13g %.13g %.13g\n",general_data->cell.hmat[7],
                                        general_data->cell.hmat[8],
                                        general_data->cell.hmat[9]);
      scanf("%d",&iii);
#endif
/*-------------------------------------------------------------------------*/
  }/*end routine*/
/*==========================================================================*/

