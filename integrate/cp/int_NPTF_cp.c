/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
/*                                                                          */
/*                         PI_MD:                                           */
/*             The future of simulation technology                          */
/*             ------------------------------------                         */
/*                     Module: int_NPTF                                     */
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
#include "../proto_defs/proto_integrate_md_entry.h"
#include "../proto_defs/proto_integrate_md_local.h"
#include "../proto_defs/proto_integrate_cp_entry.h"
#include "../proto_defs/proto_integrate_cp_local.h"
#include "../proto_defs/proto_energy_cpcon_entry.h"
#include "../proto_defs/proto_energy_ctrl_cp_entry.h"
#include "../proto_defs/proto_intra_con_entry.h"
#include "../proto_defs/proto_energy_ctrl_entry.h"
#include "../proto_defs/proto_math.h"
#include "../proto_defs/proto_communicate_wrappers.h"
#include "../proto_defs/proto_coords_local.h"

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void int_NPTF_cp(CLASS *class,BONDED *bonded,GENERAL_DATA *general_data,CP *cp)

/*========================================================================*/
{/*begin routine*/
/*========================================================================*/
/*             Local variable declarations                                */

#include "../typ_defs/typ_mask.h"
   int iii;
   int iflag;
   int ir_tra=1,ir_ter=1,ir_tor=1;

   double dt            = general_data->timeinfo.dt;
   int cp_md            = general_data->simopts.cp;
   int cp_gauss         = cp->cpopts.cp_gauss;

   int natm_tot         = class->clatoms_info.natm_tot;

   int nres_nhc         = class->therm_info_class.nres_nhc;

   int myid_state       = class->communicate.myid_state;
   int np_forc          = class->communicate.np_forc;
   int np_states        = class->communicate.np_states;
   MPI_Comm comm_states = class->communicate.comm_states;

/*==========================================================================*/
/* 0) Useful constants                                                      */

   (general_data->timeinfo.int_res_tra)    = 0;
   (general_data->timeinfo.int_res_ter)    = 0;
   class->therm_info_class.wght           = 1.0;
   class->therm_info_class.dt_nhc         = dt;
   class->therm_info_class.dti_nhc        = dt/( (double)(nres_nhc) );

   if(myid_state == 0 || np_forc > 1){
     set_yosh(class->therm_info_class.nyosh_nhc,
              class->therm_info_class.dti_nhc,class->therm_info_class.wdti,
              class->therm_info_class.wdti2,class->therm_info_class.wdti4,
              class->therm_info_class.wdti8,
              class->therm_info_class.wdti16);
   }/*endif*/

   general_data->stat_avg.iter_shake = 0; 
   general_data->stat_avg.iter_ratl  = 0; 
   general_data->stat_avg.iter_23r   = 0.0;
   general_data->stat_avg.iter_33r   = 0.0;
   general_data->stat_avg.iter_46r   = 0.0;
   general_data->stat_avg.iter_23    = 0.0;
   general_data->stat_avg.iter_33    = 0.0;
   general_data->stat_avg.iter_46    = 0.0;

/*==========================================================================*/
/* I)  Perform first CP integration                                            */

   int_0_to_dt2_cp(class,bonded,general_data,cp);

/*==========================================================================*/
/* II) First half of atom integration                                       */

   if(myid_state == 0 || np_forc > 1){
     if(cp_md==1){
       int_0_to_dt2_nptf(class,bonded,general_data,ir_tra,ir_tor,ir_ter,dt);
     }/*endif*/
   }/*endif*/

/*==========================================================================*/
/* III) Communication                                                       */

   if(np_states > 1 && np_forc == 1){
     comm_coord_class_state(class,general_data,1,1,0,0);
     Bcast(&(general_data->cell.vol),1,MPI_DOUBLE,0,class->communicate.world);
     Bcast(&(general_data->cell.hmat[1]),9,MPI_DOUBLE,0,
           class->communicate.world);
   }/* endif */

/*==========================================================================*/
/* IV) Get the new energy/force                                             */

   (class->energy_ctrl.iget_full_inter)= 1;
   (class->energy_ctrl.iget_res_inter) = 0;
   (class->energy_ctrl.iget_full_intra)= 1;
   (class->energy_ctrl.iget_res_intra) = 0;


   cp_energy_control(class,bonded,general_data,cp);

   if(cp_gauss == 1){
     add_gauss_force(&(cp->cpcoeffs_info),&(cp->cpcoeffs_pos[1]),
		     &(cp->cpscr.cpscr_ovmat),&(cp->cpopts));
   }/*endif*/

/*==========================================================================*/
/* V) Second half of atom integration                                       */

  if(myid_state == 0 || np_forc > 1){
    if(cp_md==1){
      int_dt2_to_dt_nptf(class,bonded,general_data,ir_tra,ir_tor,ir_ter,dt);
    }/*endif*/
  }/*endif : myid=0*/

/*==========================================================================*/
/* VI) Perform second part of cp integration from dt2 to dt                 */

  int_dt2_to_dt_cp(class,bonded,general_data,cp);


/*==========================================================================*/
/* VII) Finalize                                                              */

  if(myid_state == 0 || np_forc > 1){
    iflag=2;
    int_final_class(class,bonded,general_data,iflag);
  }

/*==========================================================================*/
/*     printf("x(1),y(1),z(1) %.13g %.13g %.13g\n",class->clatoms_pos[1].x[1],
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
*/
/*--------------------------------------------------------------------------*/
/*end routine*/}
/*==========================================================================*/
















