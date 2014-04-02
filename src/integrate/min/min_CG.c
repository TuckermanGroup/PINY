/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
/*                                                                          */
/*                         PI_MD:                                           */
/*             The future of simulation technology                          */
/*             ------------------------------------                         */
/*                     Module: min_CG                                       */
/*                                                                          */
/* This subprogram minimizes a classical system using conjugate gradient    */
/*                                                                          */
/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/


#include "standard_include.h"
#include "../typ_defs/typedefs_class.h"
#include "../typ_defs/typedefs_bnd.h"
#include "../typ_defs/typedefs_gen.h"
#include "../proto_defs/proto_integrate_min_entry.h"
#include "../proto_defs/proto_integrate_min_local.h"
#include "../proto_defs/proto_intra_con_entry.h"
#include "../proto_defs/proto_energy_ctrl_entry.h"
#include "../proto_defs/proto_friend_lib_entry.h"

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
void min_CG(CLASS *class,BONDED *bonded,GENERAL_DATA *general_data,int ireset)
/*========================================================================*/
{/*begin routine*/
/*========================================================================*/
/*             Local variable declarations                                */
    int i,ipart,ifirst=1;
    int natm_tot;
    double dt,tol_glob,dt22,dt2;
    double dtsm,dtsm2;
    double gamma;
    double fovlap_old;
    static double fovlap;
    double *zeta;

/* Local Pointers         */
   int min_atm_com_fix_opt = general_data->minopts.min_atm_com_fix_opt;

   double *clatoms_x    = class->clatoms_pos[1].x;
   double *clatoms_y    = class->clatoms_pos[1].y;
   double *clatoms_z    = class->clatoms_pos[1].z;
   double *clatoms_vx   = class->clatoms_pos[1].vx;
   double *clatoms_vy   = class->clatoms_pos[1].vy;
   double *clatoms_vz   = class->clatoms_pos[1].vz;

   double *clatoms_fx   = class->clatoms_pos[1].fx;
   double *clatoms_fy   = class->clatoms_pos[1].fy;
   double *clatoms_fz   = class->clatoms_pos[1].fz;
   double *chx          = class->clatoms_pos[1].vx;
   double *chy          = class->clatoms_pos[1].vy;
   double *chz          = class->clatoms_pos[1].vz;
   double *clatoms_mass  = class->clatoms_info.mass;
   double *clatoms_xold = class->clatoms_info.xold;
   double *clatoms_yold = class->clatoms_info.yold;
   double *clatoms_zold = class->clatoms_info.zold;
/*==========================================================================*/
/* 0) Useful constants                                                      */

    (class->energy_ctrl.int_res_tra)    = 0;
    (class->energy_ctrl.int_res_ter)    = 0;
    dt   = (general_data->timeinfo.dt);
    dt2  = dt/2.0;
    dt22 = dt*dt/2.0;
    dtsm = .01; /* time step to use when estimating constraint forces */
    dtsm2 = dtsm/2.0;
    general_data->stat_avg.iter_shake = 0;
    general_data->stat_avg.iter_ratl = 0;
    natm_tot = class->clatoms_info.natm_tot;

    zeta = (double *) cmalloc(natm_tot*sizeof(double))-1;

/*==========================================================================*/
/* 0.1) Zero conjugate gradients                                            */

   if(ireset == 1){
     for(i=1;i<=natm_tot; i++){
      chx[i] = 0.0;
      chy[i] = 0.0;
      chz[i] = 0.0;
     }
     gamma = 0.0;
     fovlap = 1.0;
    }/* endif */


/*==========================================================================*/
/* I)Save positions                                                         */

    for(ipart=1;ipart<=(class->clatoms_info.natm_tot);ipart++){
      clatoms_xold[ipart] = clatoms_x[ipart];
      clatoms_yold[ipart] = clatoms_y[ipart];
      clatoms_zold[ipart] = clatoms_z[ipart];
    /*endfor*/}

/*==========================================================================*/
/* II) Get new energy and forces                                            */

    (class->energy_ctrl.iget_full_inter)= 1;
    (class->energy_ctrl.iget_res_inter) = 0;
    (class->energy_ctrl.iget_full_intra)= 1;
    (class->energy_ctrl.iget_res_intra) = 0;
    energy_control(class,bonded,general_data);

    if(min_atm_com_fix_opt==1){
       proj_com_out(class->clatoms_info.natm_tot,
                    clatoms_fx,clatoms_fy,clatoms_fz);
    }/*endif*/


/*==========================================================================*/
/* II.III) Get an estimate of the constraint forces                           */

   if((bonded->constrnt.iconstrnt)==1){

/* i) evolve positions using small time step                                */

    for(ipart=1;ipart<=(class->clatoms_info.natm_tot);ipart++){
      (clatoms_vx)[ipart] += (clatoms_fx)[ipart]*dtsm2/((clatoms_mass)[ipart]);
      (clatoms_vy)[ipart] += (clatoms_fy)[ipart]*dtsm2/((clatoms_mass)[ipart]);
      (clatoms_vz)[ipart] += (clatoms_fz)[ipart]*dtsm2/((clatoms_mass)[ipart]);
    }/*endfor*/

    for(ipart=1;ipart<=(class->clatoms_info.natm_tot);ipart++){
      clatoms_x[ipart] += clatoms_vx[ipart]*dtsm;
      clatoms_y[ipart] += clatoms_vy[ipart]*dtsm;
      clatoms_z[ipart] += clatoms_vz[ipart]*dtsm;
    }/*endfor*/

/* ii) zero the velocities before calling shake                              */
/* velocities returned will not be zero  =(constraint force)*(dts)/(2.0*mass)*/

    for(ipart=1;ipart<=(class->clatoms_info.natm_tot);ipart++){
      clatoms_vx[ipart] = 0.0;
      clatoms_vy[ipart] = 0.0;
      clatoms_vz[ipart] = 0.0;
    }/*endfor*/

    init_constraint(bonded,&(general_data->ptens));

    (bonded->constrnt.iroll) = 0; /*do not need to roll */
     shake_control(bonded,&(class->clatoms_info),&(class->clatoms_pos[1]), 
                  &(general_data->cell),&(general_data->ptens),
                  &(general_data->statepoint),
                  &(general_data->baro),&(general_data->par_rahman),
                  &(general_data->stat_avg),dt,&tol_glob,ifirst,
                  &(class->class_comm_forc_pkg),&(class->ewd_scr));


/* iii) add in constraint force to force obtained from energy control     */
/* constraint forces are in velocity vector   v=  dt*force/mass           */
   for(ipart=1;ipart<=(class->clatoms_info.natm_tot);ipart++){
      clatoms_fx[ipart] += clatoms_vx[ipart]*clatoms_mass[ipart]/dtsm2;
      clatoms_fy[ipart] += clatoms_vy[ipart]*clatoms_mass[ipart]/dtsm2;
      clatoms_fz[ipart] += clatoms_vz[ipart]*clatoms_mass[ipart]/dtsm2;
    /*endfor*/}


/* iv) return positions to their old values                              */

    for(ipart=1;ipart<=(class->clatoms_info.natm_tot);ipart++){
      clatoms_x[ipart] = clatoms_xold[ipart];
      clatoms_y[ipart] = clatoms_yold[ipart];
      clatoms_z[ipart] = clatoms_zold[ipart];
    }/*endfor*/
   
   }/*endif for constraints*/


/*==========================================================================*/
/* II) Calculate the gamma's                                                */

   fovlap_old = fovlap;
   fovlap = 0.0;
   for(i=1;i<=natm_tot;i++){
    fovlap += clatoms_fx[i]*clatoms_fx[i]
            + clatoms_fy[i]*clatoms_fy[i]
            + clatoms_fz[i]*clatoms_fz[i];
   }
   if(ireset != 1) {gamma = fovlap/fovlap_old;}

/*==========================================================================*/
/* II.V) Evolve gradients                                                   */

     for(i=1;i<=natm_tot; i++){
      chx[i] = clatoms_fx[i] + gamma*chx[i];
      chy[i] = clatoms_fy[i] + gamma*chy[i];
      chz[i] = clatoms_fz[i] + gamma*chz[i];
     }

/*==========================================================================*/
/* III) Calculate the step length                                           */

   for(i=1;i<=natm_tot;i++) {
    zeta[i] = dt/clatoms_mass[i];
   }/* endfor */

/*==========================================================================*/
/* IV) Evolve positions                                                     */

    for(ipart=1;ipart<=(class->clatoms_info.natm_tot);ipart++){
      clatoms_x[ipart] = clatoms_xold[ipart] + zeta[ipart]*chx[ipart];
      clatoms_y[ipart] = clatoms_yold[ipart] + zeta[ipart]*chy[ipart];
      clatoms_z[ipart] = clatoms_zold[ipart] + zeta[ipart]*chz[ipart];
    /*endfor*/}

/*==========================================================================*/
/* V) Shake if necessary: Put particles on surface and add in constraint    */
/*                          force so the force test in control_min          */
/*                          will converge.                                  */

    if((bonded->constrnt.iconstrnt)==1){
         (bonded->constrnt.iroll) = 0;

         init_constraint(bonded,&(general_data->ptens));

    for(ipart=1;ipart<=(class->clatoms_info.natm_tot);ipart++){
      clatoms_vx[ipart] = 0.0;
      clatoms_vy[ipart] = 0.0;
      clatoms_vz[ipart] = 0.0;
    /*endfor*/}


         shake_control(bonded,&(class->clatoms_info),&(class->clatoms_pos[1]), 
                      &(general_data->cell),&(general_data->ptens),
                      &(general_data->statepoint),
                      &(general_data->baro),&(general_data->par_rahman),
                      &(general_data->stat_avg),dt,&tol_glob,ifirst,
                      &(class->class_comm_forc_pkg),&(class->ewd_scr));
  /*endif*/}  

/*==========================================================================*/
/* VI) Get positions of ghost atoms                                         */

   get_ghost_pos(&(class->clatoms_info),&(class->clatoms_pos[1]),
                 &(class->ghost_atoms));

/*==========================================================================*/
/* VII) Get new energy and forces                                           */

    (class->energy_ctrl.iget_full_inter)= 1;
    (class->energy_ctrl.iget_res_inter) = 0;
    (class->energy_ctrl.iget_full_intra)= 1;
    (class->energy_ctrl.iget_res_intra) = 0;
    energy_control(class,bonded,general_data);


   cfree(&(zeta[1]));

/*--------------------------------------------------------------------------*/
/*end routine*/}
/*==========================================================================*/

