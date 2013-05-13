/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
/*                                                                          */
/*                         PI_MD:                                           */
/*             The future of simulation technology                          */
/*             ------------------------------------                         */
/*                     Module: int_NPTF_res_pimd                            */
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
#include "../typ_defs/typedefs_stat.h"
#include "../proto_defs/proto_communicate_wrappers.h"
#include "../proto_defs/proto_integrate_pimd_entry.h"
#include "../proto_defs/proto_integrate_pimd_local.h"
#include "../proto_defs/proto_integrate_md_local.h"
#include "../proto_defs/proto_integrate_md_entry.h"
#include "../proto_defs/proto_intra_con_entry.h"
#include "../proto_defs/proto_energy_ctrl_entry.h"
#include "../proto_defs/proto_pimd_entry.h"
#include "../proto_defs/proto_pimd_local.h"
#include "../proto_defs/proto_math.h"


/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void int_NPTF_pimd(CLASS *class,BONDED *bonded,GENERAL_DATA *general_data,
                   ANALYSIS *analysis)

/*========================================================================*/
  {/*begin routine*/
/*========================================================================*/
/*             Local variable declarations                                */

#include "../typ_defs/typ_mask.h"

   int iii;
   int iflag;
   int ir_tra,ir_tor,ir_ter,ir_pimd;
   double dti;

   int pi_beads        = class->clatoms_info.pi_beads;
   int pi_beads_proc   = class->clatoms_info.pi_beads_proc;

   MPI_Comm comm_beads = class->communicate.comm_beads;
   int myid            = class->communicate.myid;

   int ix_respa        = general_data->timeinfo.ix_respa;
   int int_res_ter     = general_data->timeinfo.int_res_ter;
   int nres_tra        = general_data->timeinfo.nres_tra;
   int nres_tor        = general_data->timeinfo.nres_tor;
   int nres_ter        = general_data->timeinfo.nres_ter;
   int nres_pimd       = general_data->timeinfo.nres_pimd;
   double dt           = general_data->timeinfo.dt;

   int nres_nhc        = class->therm_info_class.nres_nhc;

   int iperd           = general_data->cell.iperd;
   int hmat_cons_typ   = general_data->cell.hmat_cons_typ;
   int hmat_int_typ    = general_data->cell.hmat_int_typ;
   double *vgmat       = general_data->par_rahman.vgmat;
   double *hmat        = general_data->cell.hmat;

/*==========================================================================*/
/* 0) Useful constants                                                      */

   dti  = dt/((double)(nres_ter*nres_tor*nres_tra*nres_pimd));
   if(general_data->timeinfo.ix_respa==1){
      class->therm_info_class.wght  = 
                     (double)(nres_pimd*nres_tra*nres_tor*nres_ter);
      class->therm_info_bead.wght  = 
                    (double)(nres_pimd*nres_tra*nres_tor*nres_ter);
   }/*endif*/

   if(ix_respa==2){
      class->therm_info_class.wght  = (double)(nres_pimd*nres_tra*nres_tor);
      class->therm_info_bead.wght   = (double)(nres_pimd*nres_tra*nres_tor);
   }/*endif*/
   if(ix_respa==3){
      class->therm_info_class.wght  =  (double)(nres_pimd*nres_tra);
      class->therm_info_bead.wght   =  (double)(nres_pimd*nres_tra);
   }/*endif*/
   if(ix_respa==4){
      class->therm_info_class.wght  =  (double)(nres_pimd);
      class->therm_info_bead.wght   =  (double)(nres_pimd);
   }/*endif*/
   if(ix_respa==5){
      class->therm_info_class.wght  =  1.0;
      class->therm_info_bead.wght   = 1.0;
   }/*endif*/


   class->therm_info_class.dt_nhc   =  dti;
   class->therm_info_bead.dt_nhc    =  dti;
   class->therm_info_class.dti_nhc  =  dti/( (double)(nres_nhc) );
   class->therm_info_bead.dti_nhc   =  dti/( (double)(nres_nhc) );

   set_yosh(class->therm_info_class.nyosh_nhc,
            class->therm_info_class.dti_nhc,
            class->therm_info_class.wdti,
            class->therm_info_class.wdti2,class->therm_info_class.wdti4,
            class->therm_info_class.wdti8,
            class->therm_info_class.wdti16);

   set_yosh(class->therm_info_bead.nyosh_nhc,
            class->therm_info_bead.dti_nhc,
            class->therm_info_bead.wdti,
            class->therm_info_bead.wdti2,class->therm_info_bead.wdti4,
            class->therm_info_bead.wdti8,
            class->therm_info_bead.wdti16);

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
/*==========================================================================*/
/* IV) Loop over bead RESPA                                               */

        for(ir_pimd=1;ir_pimd<=nres_pimd;ir_pimd++){

/*==========================================================================*/
/* 1) Evolve system from t=0 to dt/2                                        */


          int_0_to_dt2_nptf_pimd(class,bonded,general_data,
                                 ir_tra,ir_tor,ir_ter,ir_pimd,dti); 

/*==========================================================================*/
/* 2) Get the new energy/force                                              */

          (class->energy_ctrl.iget_full_inter) = 0;
          (class->energy_ctrl.iget_res_inter) = 0;
          (class->energy_ctrl.iget_full_intra) = 0;
          (class->energy_ctrl.iget_res_intra)  = 0;
          (class->energy_ctrl.iget_res_pimd)  = 1;
          if((ir_ter==nres_ter)&&(ir_tor==nres_tor)&&(ir_tra==nres_tra)
                               &&(ir_pimd==nres_pimd))
             {(class->energy_ctrl.iget_full_inter) = 1;}
          if((int_res_ter==1)&&(ir_tor==nres_tor)&&(ir_tra==nres_tra)
                             &&(ir_pimd==nres_pimd))
             {(class->energy_ctrl.iget_res_inter) = 1;}
          if((ir_tra==nres_tra)&&(ir_pimd==nres_pimd))
             {(class->energy_ctrl.iget_full_intra) = 1;}
          if((ir_pimd==nres_pimd))
             {(class->energy_ctrl.iget_res_intra) = 1;}

          energy_control_pimd(class,bonded,general_data);

/*==========================================================================*/
/* 3) Evolve system dt/2 to dt                                              */

          int_dt2_to_dt_nptf_pimd(class,bonded,general_data,
                                  ir_tra,ir_tor,ir_ter,ir_pimd,dti); 

/*==========================================================================*/
/*==========================================================================*/

        }/*endfor:ir_pimd*/
       }/*endfor:ir_tra*/
     }/*endfor:ir_tor*/
   }/*endfor:ir_ter*/

/*==========================================================================*/
/*==========================================================================*/
/* IV) Get Kinetic energy and kinetic energy tensor                         */

    get_tvten_pimd(class,general_data);

/*==========================================================================*/
/*==========================================================================*/
/* VI) Get NHC contribution to energy                                       */

    iflag = 2;
    nhc_vol_potkin_pimd(class,general_data,iflag);

/*==========================================================================*/
/*==========================================================================*/
/* Check the cell matrix */

   chck_constr_cell_mat(iperd,hmat_cons_typ,hmat_int_typ,vgmat,hmat,myid);

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
      printf("vx(1),vy(1),vz(1) %.13g %.13g %.13g\n",
                                            class->clatoms_pos[2].vx[1],
                                            class->clatoms_pos[2].vy[1],
                                            class->clatoms_pos[2].vz[1]); 
      printf("v_nhc[1][1],v_nhc[2,1] %.13g %.13g\n",
                                  class->therm_class.v_nhc[1][1],
                                  class->therm_class.v_nhc[2][1]); 
      printf("x_nhc[1][1],x_nhc[2][1] %.13g %.13g\n",
                                  class->therm_class.x_nhc[1][1],
                                  class->therm_class.x_nhc[2][1]); 
      printf("bead v_nhc[1][1],v_nhc[2,1] %.13g %.13g\n",
                                  class->therm_bead[2].v_nhc[1][1],
                                  class->therm_bead[2].v_nhc[2][1]); 
      printf("bead x_nhc[1][1],x_nhc[2][1] %.13g %.13g\n",
                                  class->therm_bead[2].x_nhc[1][1],
                                  class->therm_bead[2].x_nhc[2][1]); 
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
/*-------------------------------------------------------------------------*/
    }/*end routine*/
/*==========================================================================*/

