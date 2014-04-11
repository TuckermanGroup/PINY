/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
/*                                                                          */
/*                         PI_MD:                                           */
/*             The future of simulation technology                          */
/*             ------------------------------------                         */
/*                     Module: min_DIIS                                     */
/*                                                                          */
/* This subprogram minimizes atomic positions using DIIS                    */
/*                                                                          */
/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/


#include "standard_include.h"
#include "../typ_defs/typedefs_class.h"
#include "../typ_defs/typedefs_bnd.h"
#include "../typ_defs/typedefs_gen.h"
#include "../typ_defs/typedefs_cp.h"
#include "../proto_defs/proto_integrate_min_entry.h"
#include "../proto_defs/proto_integrate_min_local.h"
#include "../proto_defs/proto_intra_con_entry.h"
#include "../proto_defs/proto_energy_ctrl_entry.h"
#include "../proto_defs/proto_friend_lib_entry.h"
#include "../proto_defs/proto_communicate_wrappers.h"
#include "../proto_defs/proto_math.h"

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void min_DIIS(CLASS *class,BONDED *bonded,GENERAL_DATA *general_data,
              CP *cp)

/*========================================================================*/
{/*begin routine*/
/*========================================================================*/
/* Local Variable declarations                                          */
#include "../typ_defs/typ_mask.h"

    int ifirst=1;
    int i,ip=1,ipart,icoef,is,ist,iist,job;
    double dt,dt2,dt22;
    double dtsm,dtsm2;
    int natm_tot,iii;
    int info;
    static int itcount=1;
    int mvec,mp1;
    double tol_glob;


    int min_atm_com_fix_opt = general_data->minopts.min_atm_com_fix_opt;

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
    double *chx  = class->clatoms_pos[1].vx;
    double *chy  = class->clatoms_pos[1].vy;
    double *chz  = class->clatoms_pos[1].vz;
    double *scr_x  = class->clatoms_info.xold;
    double *scr_y  = class->clatoms_info.yold;
    double *scr_z  = class->clatoms_info.zold;
    double *clatoms_mass = class->clatoms_info.mass;
    int MVEC_MAX_ATM = general_data->minopts.diis_hist_len;

/* DIIS local pointers */
    static int *ipvt;
    static double *pmat;
    static double *vect;
    static double *xst,*yst,*zst;
    static double *fxst,*fyst,*fzst;

/*==========================================================================*/
/* 0) Useful constants                                                      */
    
    natm_tot = class->clatoms_info.natm_tot;
    general_data->timeinfo.int_res_tra    = 0;
    general_data->timeinfo.int_res_ter    = 0;
    dt  = general_data->timeinfo.dt;
    dt2  = dt/2.0;
    dt22 = dt*dt/2.0;
    dtsm = .01; /* time step to use when estimating constraint forces */
    dtsm2 = dtsm/2.0;
    general_data->stat_avg.iter_shake = 0;
    general_data->stat_avg.iter_ratl = 0;
    zero_constrt_iters(&(general_data->stat_avg));



/*========================================================================*/
/* I) Allocate DIIS memory for this step                                  */

    mvec = ( itcount < MVEC_MAX_ATM ? itcount:MVEC_MAX_ATM);
    mp1 = mvec + 1;

    if(itcount==1){
      init_alloc_atm_DIIS_mem(&ipvt,&pmat,&vect,&xst,&yst,&zst,
                              &fxst,&fyst,&fzst,mp1,natm_tot,mvec);
    } else if(itcount > 1 && itcount <= MVEC_MAX_ATM) {
      realloc_atm_DIIS_mem(&ipvt,&pmat,&vect,&xst,&yst,&zst,
                           &fxst,&fyst,&fzst,mp1,natm_tot,mvec);
    } else {
      shift_atm_DIIS_mem(xst,yst,zst,fxst,fyst,fzst,mvec,natm_tot);
    }

/*========================================================================*/
/* II) Get forces                                                         */

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

/*========================================================================*/
/* III) Get approximate gradients, evolve positions, save new histories   */

    ist=(mvec-1)*natm_tot;
    for(i=1;i<=natm_tot;i++){
      iist = ist + i;
      fxst[iist] = dt*clatoms_fx[i];
      fyst[iist] = dt*clatoms_fy[i];
      fzst[iist] = dt*clatoms_fz[i];
      clatoms_x[i] += fxst[iist];
      clatoms_y[i] += fyst[iist];
      clatoms_z[i] += fzst[iist];
      xst[iist] = clatoms_x[i];
      yst[iist] = clatoms_y[i];
      zst[iist] = clatoms_z[i];
    }/* endfor i */
/*==========================================================================*/
/* IV) Construct DIIS linear equations                                    */ 

  setup_atm_DIIS_eqs(pmat,vect,fxst,fyst,fzst,natm_tot,mp1,mvec);

/*==========================================================================*/
/* V) Solve DIIS equations to get DIIS vectors                           */

#ifdef IBM_ESSL
      dgef(&(pmat[1]),&mp1,&mp1,&(ipvt[1]));
      job = 0;
      dges(&(pmat[1]),&mp1,&mp1,&(ipvt[1]),&(vect[1]),&job);
#else
      DGEFA(&(pmat[1]),&mp1,&mp1,&(ipvt[1]),&info);
      job=1;
      DGESL(&(pmat[1]),&mp1,&mp1,&(ipvt[1]),&(vect[1]),&job);
#endif
/*==========================================================================*/
/* VI) Compute new positions from DIIS vectors                             */

   for(i=1;i<=natm_tot;i++){
     clatoms_x[i] = 0.0;
     clatoms_y[i] = 0.0;
     clatoms_z[i] = 0.0;
   }/* endfor i */ 
   for(i=1;i<=mvec;i++){
     ist=(i-1)*natm_tot;
     for(is=1;is<=natm_tot;is++){
       iist = ist+is;
       clatoms_x[is] += vect[i]*xst[iist];
       clatoms_y[is] += vect[i]*yst[iist];
       clatoms_z[is] += vect[i]*zst[iist];
     }/* endfor is */
   }/* endfor i */

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
/* VII) Increment iteration counter                                         */

   ++itcount;

/*========================================================================*/
}/* end routine */
/*========================================================================*/

