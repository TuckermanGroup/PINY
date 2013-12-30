/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
/*                                                                          */
/*                         PI_MD:                                           */
/*             The future of simulation technology                          */
/*             ------------------------------------                         */
/*                     Module: get_estimators                               */
/*                                                                          */
/* This subprogram calculates PIMD energy estimators                        */
/*                                                                          */
/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/



#include "standard_include.h"
#include "../typ_defs/typedefs_class.h"
#include "../typ_defs/typedefs_gen.h"
#include "../typ_defs/typedefs_bnd.h"
#include "../proto_defs/proto_pimd_local.h"
#include "../proto_defs/proto_energy_ctrl_entry.h"
#include "../proto_defs/proto_friend_lib_entry.h"
#include "../proto_defs/proto_communicate_wrappers.h"

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void mode_energy_control(CLASS *class, GENERAL_DATA *general_data)

/*==========================================================================*/
{ /* begin routine */
/*==========================================================================*/

#include "../typ_defs/typ_mask.h"

  /* local variable declarations */

  int i,ip,ip_true,ipart,iii;
  double *x,*y,*z,*fx,*fy,*fz,ake,ake_temp,akprim,fk,beta,tau;
  int pi_beads = class->clatoms_info.pi_beads;
  int pi_beads_proc = class->clatoms_info.pi_beads_proc;
  int natm_tot = class->clatoms_info.natm_tot;
  int myatm_start = class->clatoms_info.myatm_start;
  int myatm_end = class->clatoms_info.myatm_end;
  int pi_beads_proc_st = class->clatoms_info.pi_beads_proc_st;
  int myid = class->communicate.myid_bead;
  int np_forc = class->communicate.np_forc;
  MPI_Comm comm_forc = class->communicate.comm_forc;
  int nghost_tot = class->ghost_atoms.nghost_tot;
  int nfreeze    = class->atommaps.nfreeze;
  double *pvten_mode     = general_data->ptens.pvten_mode;
  double *pvten_mode_tot = general_data->ptens.pvten_mode_tot;
  double *prekf = class->clatoms_info.prekf;
  double *veig  = class->clatoms_tran.path_eig;
  double wght   = class->clatoms_info.wght_pimd;
  double pi_temperature = class->clatoms_info.pi_temperature;


  /*======================================*/
  /* calculate the KE and pressure tensor */

  ake = 0.0;
  beta = BOLTZ / pi_temperature;
  tau = beta / pi_beads_proc;
  for(i=1; i<=9; i++) {
      pvten_mode_tot[i] = 0.0;
  }
  for(ip=1; ip<=pi_beads_proc; ip++) {
    ip_true = ip + pi_beads_proc_st - 1;
    x = class->clatoms_pos[ip].x;
    y = class->clatoms_pos[ip].y;
    z = class->clatoms_pos[ip].z;
    fx = class->clatoms_pos[ip].fxm;
    fy = class->clatoms_pos[ip].fym;
    fz = class->clatoms_pos[ip].fzm;
    for(ipart=myatm_start; ipart<=myatm_end; ipart++) {
      fk = prekf[ipart] * veig[ip_true];
      ake += fk * (x[ipart]*x[ipart] + y[ipart]*y[ipart] + z[ipart]*z[ipart]);
      fx[ipart] = -2.0 * fk * x[ipart];
      fy[ipart] = -2.0 * fk * y[ipart];
      fz[ipart] = -2.0 * fk * z[ipart];
      pvten_mode_tot[1] += x[ipart] * fx[ipart];
      pvten_mode_tot[5] += y[ipart] * fy[ipart];
      pvten_mode_tot[9] += z[ipart] * fz[ipart];
      pvten_mode_tot[2] += x[ipart] * fy[ipart];
      pvten_mode_tot[3] += x[ipart] * fz[ipart];
      pvten_mode_tot[6] += y[ipart] * fz[ipart];
      fx[ipart] *= wght;
      fy[ipart] *= wght;
      fz[ipart] *= wght;
    }
  }
  pvten_mode_tot[4] = pvten_mode_tot[2];
  pvten_mode_tot[7] = pvten_mode_tot[3];
  pvten_mode_tot[8] = pvten_mode_tot[6];
  for(i=1; i<=9; i++) {
    pvten_mode[i] = wght * pvten_mode_tot[i];
  }
  general_data->stat_avg.kin_harm = ake;
  if (np_forc > 1) {
    ake_temp = 0.0;
    Reduce(&ake, &ake_temp, 1, MPI_DOUBLE, MPI_SUM, 0, comm_forc);
    ake = ake_temp;
  }

  akprim = 1.5 * ((double) (natm_tot - nghost_tot - nfreeze)) / tau - ake;
  general_data->stat_avg.pi_ke_prim = akprim;

/*--------------------------------------------------------------------------*/
} /* end routine */
/*==========================================================================*/


/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void get_pimd_spread(CLATOMS_INFO *clatoms_info, CLATOMS_POS *clatoms_pos,
                     double *spread, COMMUNICATE *communicate)

/*==========================================================================*/
{ /* begin routine */
/*==========================================================================*/

  #include "../typ_defs/typ_mask.h"

  /* local variable declarations */
  int ip,imall,ipart,iii;
  int pi_beads = clatoms_info->pi_beads;
  int pi_beads_proc = clatoms_info->pi_beads_proc;
  int pi_beads_proc_st = clatoms_info->pi_beads_proc_st;
  int rank = communicate->myid_bead;
  int natm_tot = clatoms_info->natm_tot;
  int myatm_start = clatoms_info->myatm_start;
  int myatm_end = clatoms_info->myatm_end;
  double *x1, *y1, *z1, *xp, *yp, *zp;
  double r2, dx, dy, dz;
  double spread_now;
  MPI_Comm comm_beads = communicate->comm_beads;


  /*===================================*/
  /* calculate the spread of the beads */

  *spread = 0.0;
  imall = 0;
  if (pi_beads_proc_st == 1) {
    x1 = clatoms_pos[1].x;
    y1 = clatoms_pos[1].y;
    z1 = clatoms_pos[1].z;
  } else {
    imall = 1;
    x1 = (double *)cmalloc(natm_tot*sizeof(double))-1;
    y1 = (double *)cmalloc(natm_tot*sizeof(double))-1;
    z1 = (double *)cmalloc(natm_tot*sizeof(double))-1;
  }

  if (communicate->np_beads > 1) {
    Bcast(&(x1[1]), natm_tot, MPI_DOUBLE, 0, communicate->comm_beads);
    Bcast(&(y1[1]), natm_tot, MPI_DOUBLE, 0, communicate->comm_beads);
    Bcast(&(z1[1]), natm_tot, MPI_DOUBLE, 0, communicate->comm_beads);
    Barrier(communicate->comm_beads);
  }

  if (pi_beads_proc_st == 1) {
   for(ip=2; ip<=pi_beads_proc; ip++) {
     xp = clatoms_pos[ip].x;
     yp = clatoms_pos[ip].y;
     zp = clatoms_pos[ip].z;
     for(ipart=1;ipart<=natm_tot;ipart++){
       dx = xp[ipart] - x1[ipart];
       dy = yp[ipart] - y1[ipart];
       dz = zp[ipart] - z1[ipart];
       r2 = dx*dx + dy*dy + dz*dz;
       *spread = MAX(*spread,r2);
     }/*endfor*/
   }/*endfor*/
 }else{
   for(ip=1;ip<=pi_beads_proc;ip++){
     xp = clatoms_pos[ip].x;
     yp = clatoms_pos[ip].y;
     zp = clatoms_pos[ip].z;
     for(ipart=1;ipart<=natm_tot;ipart++){
       dx = xp[ipart] - x1[ipart];
       dy = yp[ipart] - y1[ipart];
       dz = zp[ipart] - z1[ipart];
       r2 = dx*dx + dy*dy + dz*dz;
       *spread = MAX(*spread,r2);
     }/*endfor*/
   }/*endfor*/
 }/*endif*/

 if(communicate->np_beads>1){
  Allreduce(spread,&spread_now,1,MPI_DOUBLE,MPI_MAX,0,communicate->comm_beads);
  Barrier(communicate->comm_beads);
  *spread = sqrt(spread_now);
 }else{
  *spread = sqrt(*spread);
 }/*endif*/

if(imall==1){
    cfree(&x1[1]);
    cfree(&y1[1]);
    cfree(&z1[1]);
}/*endif*/

/*--------------------------------------------------------------------------*/
}/*end routine*/
/*==========================================================================*/






/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void pimd_pvten_forc_corr(CLASS *class, GENERAL_DATA *general_data, 
                                                     double *akvirx_ret)

/*==========================================================================*/
{/*begin routine*/
/*======================================================================*/
/*     Local  Variables     */

  int ip,i,iii;
  double *x,*y,*z,*fx,*fy,*fz,*fxt,*fyt,*fzt,*fxm,*fym,*fzm;
  double beta,akvirx,sum,vir_temp;
  int np_forc = class->communicate.np_forc;

/*     Local Pointers */

  int pi_beads       = class->clatoms_info.pi_beads; 
  int pi_beads_proc  = class->clatoms_info.pi_beads_proc; 
  int myid           = class->communicate.myid_bead; 
  int natm_tot       = class->clatoms_info.natm_tot; 
  int myatm_start = class->clatoms_info.myatm_start;
  int myatm_end = class->clatoms_info.myatm_end;
  int nghost_tot = class->ghost_atoms.nghost_tot;
  int nfreeze    = class->atommaps.nfreeze;
  double *xmod       = class->clatoms_info.xmod; 
  double *ymod       = class->clatoms_info.ymod; 
  double *zmod       = class->clatoms_info.zmod; 
  double *dx         = class->ewd_scr.x;
  double *dy         = class->ewd_scr.y;
  double *dz         = class->ewd_scr.z;
  double *pvten      = general_data->ptens.pvten;
  double *pvten_tot  = general_data->ptens.pvten_tot;
  double t_ext       = general_data->statepoint.t_ext;
  double pi_temperature = class->clatoms_info.pi_temperature;


/*======================================================================*/

  beta   = BOLTZ/pi_temperature;
  akvirx = 0.0;
  if( ((class->energy_ctrl.iget_full_intra)==1)&&
      ((class->energy_ctrl.iget_full_inter)==1)&&myid==0){
     akvirx=3.0*((double)(natm_tot-nghost_tot-nfreeze))/beta;
  }/*endif*/

  vir_temp = 0.0;
  for(ip=1;ip<=pi_beads_proc;ip++){
     x  = class->clatoms_pos[ip].x;
     y  = class->clatoms_pos[ip].y;
     z  = class->clatoms_pos[ip].z;
     fx = class->clatoms_pos[ip].fx;
     fy = class->clatoms_pos[ip].fy;
     fz = class->clatoms_pos[ip].fz;
     for(i=1;i<=natm_tot;i++){
       dx[i] = x[i] - xmod[i];
       dy[i] = y[i] - ymod[i];
       dz[i] = z[i] - zmod[i];
     }/*endfor*/
     for(i=1;i<=natm_tot;i++){
      pvten[1] -= (dx[i]*fx[i]);
      pvten[5] -= (dy[i]*fy[i]);
      pvten[9] -= (dz[i]*fz[i]);
      pvten[2] -= (dx[i]*fy[i]);
      pvten[4] -= (dy[i]*fx[i]);
      pvten[3] -= (dx[i]*fz[i]);
      pvten[7] -= (dz[i]*fx[i]);
      pvten[6] -= (dy[i]*fz[i]);
      pvten[8] -= (dz[i]*fy[i]);
     }/*endfor*/
     if( ((class->energy_ctrl.iget_full_intra)==1)&&
        ((class->energy_ctrl.iget_full_inter)==1)){
       fxt = class->clatoms_pos[ip].fxt;
       fyt = class->clatoms_pos[ip].fyt;
       fzt = class->clatoms_pos[ip].fzt;
       for(i=1;i<=natm_tot;i++){
         akvirx       -= (dx[i]*fxt[i]+dy[i]*fyt[i]+dz[i]*fzt[i]);
         pvten_tot[1] -= (dx[i]*fxt[i]);
         pvten_tot[5] -= (dy[i]*fyt[i]);
         pvten_tot[9] -= (dz[i]*fzt[i]);
         pvten_tot[2] -= (dx[i]*fyt[i]);
         pvten_tot[4] -= (dy[i]*fxt[i]);
         pvten_tot[3] -= (dx[i]*fzt[i]);
         pvten_tot[7] -= (dz[i]*fxt[i]);
         pvten_tot[6] -= (dy[i]*fzt[i]);
         pvten_tot[8] -= (dz[i]*fyt[i]);
       }/*endfor*/
     }/*endif*/
  }/*endfor*/

  *akvirx_ret = 0.5*akvirx;

/*--------------------------------------------------------------------------*/
}/*end routine*/
/*==========================================================================*/





/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void get_vir_press(CLASS *class, GENERAL_DATA *general_data,
                                               double *press_ret)

/*==========================================================================*/
{/*begin routine*/
/*======================================================================*/
/*     Local  Variables     */

  int ip,i;
  double *x,*y,*z,*fx,*fy,*fz,*fxt,*fyt,*fzt,*fxm,*fym,*fzm;
  double beta,press;

/*     Local Pointers */

  int pi_beads       = class->clatoms_info.pi_beads; 
  int pi_beads_proc  = class->clatoms_info.pi_beads_proc; 
  double dpi_beads;
  int myid           = class->communicate.myid_bead; 
  int natm_tot       = class->clatoms_info.natm_tot; 
  int myatm_start = class->clatoms_info.myatm_start;
  int myatm_end = class->clatoms_info.myatm_end;
  double *xmod       = class->clatoms_info.xmod; 
  double *ymod       = class->clatoms_info.ymod; 
  double *zmod       = class->clatoms_info.zmod; 
  double *dx         = class->ewd_scr.x;
  double *dy         = class->ewd_scr.y;
  double *dz         = class->ewd_scr.z;
  double *pvten      = general_data->ptens.pvten;
  double *pvten_tot  = general_data->ptens.pvten_tot;
  double t_ext       = general_data->statepoint.t_ext;
  double pi_temperature = class->clatoms_info.pi_temperature;

/*======================================================================*/

  dpi_beads = (double)pi_beads;
  beta   = BOLTZ/pi_temperature;
  press = pvten_tot[1]+pvten_tot[5]+pvten_tot[9];
  for(ip=1;ip<=pi_beads_proc;ip++){
     x  = class->clatoms_pos[ip].x;
     y  = class->clatoms_pos[ip].y;
     z  = class->clatoms_pos[ip].z;
     for(i=1;i<=natm_tot;i++){
       dx[i] = x[i] - xmod[i];
       dy[i] = y[i] - ymod[i];
       dz[i] = z[i] - zmod[i];
     }/*endfor*/
     fxt = class->clatoms_pos[ip].fxt;
     fyt = class->clatoms_pos[ip].fyt;
     fzt = class->clatoms_pos[ip].fzt;
     for(i=1;i<=natm_tot;i++){
       press        -= (dx[i]*fxt[i]+dy[i]*fyt[i]+dz[i]*fzt[i]);
     }/*endfor*/
  }/*endfor*/
  *press_ret = press/(3.0*dpi_beads);

/*--------------------------------------------------------------------------*/
}/*end routine*/
/*==========================================================================*/

