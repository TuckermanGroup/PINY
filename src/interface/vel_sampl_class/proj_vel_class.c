/*===================================================================*/
/*ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*===================================================================*/
/*                                                                   */
/*                         PI_MD:                                    */
/*             The future of simulation technology                   */
/*             ------------------------------------                  */
/*                     Module: proj_vel.c                            */
/*                                                                   */
/* These subprograms projects the velocities onto surface of         */
/* constraint                                                        */
/*                                                                   */
/*===================================================================*/
/*ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*===================================================================*/

#include "standard_include.h"
#include "../typ_defs/typedefs_gen.h"
#include "../typ_defs/typedefs_class.h"
#include "../typ_defs/typedefs_bnd.h"
#include "../proto_defs/proto_intra_con_entry.h"
#include "../proto_defs/proto_vel_sampl_class_local.h"
#include "../proto_defs/proto_communicate_wrappers.h"

/*===================================================================*/
/*ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*===================================================================*/
/* Particle Velocities */
/*===================================================================*/

void proj_vel(CLASS *class,BONDED *bonded,GENERAL_DATA *general_data)

/*=================================================================*/
{/*begin routine*/
/*=================================================================*/
/*             Local variable declarations                         */

  int ifirst;
  int i;
  double dt,tol_glob;
  int myatm_start = class->clatoms_info.myatm_start;
  int myatm_end   = class->clatoms_info.myatm_end;

  double *xold   = class->clatoms_info.xold;
  double *yold   = class->clatoms_info.yold;
  double *zold   = class->clatoms_info.zold;
  double *x      = class->clatoms_pos[1].x;
  double *y      = class->clatoms_pos[1].y;
  double *z      = class->clatoms_pos[1].z;
  double *vx     = class->clatoms_pos[1].vx;
  double *vy     = class->clatoms_pos[1].vy;
  double *vz     = class->clatoms_pos[1].vz;
  double *scr_vx = class->int_scr.vx;
  double *scr_vy = class->int_scr.vy;
  double *scr_vz = class->int_scr.vz;

/*=================================================================*/
/* I) Shake once to get on surface                                */

  dt = 0.0001; 

  for(i=myatm_start;i<=myatm_end;i++){
    xold[i] = x[i];
    yold[i] = y[i];
    zold[i] = z[i];
    scr_vx[i] = vx[i];
    scr_vy[i] = vy[i];
    scr_vz[i] = vz[i];
  }/*endfor*/

  for(i=myatm_start;i<=myatm_end;i++){
    vx[i] = 0.0;
    vy[i] = 0.0;
    vz[i] = 0.0;
    x[i] = xold[i] + scr_vx[i]*dt;
    y[i] = yold[i] + scr_vy[i]*dt;
    z[i] = zold[i] + scr_vz[i]*dt;
  }/*endfor*/

  zero_constrt_iters(&(general_data->stat_avg));
  init_constraint(bonded,&(general_data->ptens));
  bonded->constrnt.iroll = 0;  ifirst = 2;
  shake_control(bonded,
                &(class->clatoms_info),&(class->clatoms_pos[1]), 
		&(general_data->cell),&(general_data->ptens),
                &(general_data->statepoint),
		&(general_data->baro),&(general_data->par_rahman),
		&(general_data->stat_avg),dt,&tol_glob,ifirst,
                &(class->class_comm_forc_pkg),&(class->ewd_scr));

  for(i=myatm_start;i<=myatm_end;i++){
    xold[i] = x[i];
    yold[i] = y[i];
    zold[i] = z[i];
  }/*endfor*/

/*=================================================================*/
/* II) Shake again to get velocities                                */

  dt = 0.0001;
  for(i=myatm_start;i<=myatm_end;i++){
    vx[i] = 0.0;
    vy[i] = 0.0;
    vz[i] = 0.0;
    x[i] = xold[i] + scr_vx[i]*dt;
    y[i] = yold[i] + scr_vy[i]*dt;
    z[i] = zold[i] + scr_vz[i]*dt;
  }/*endfor*/

  zero_constrt_iters(&(general_data->stat_avg));
  init_constraint(bonded,&(general_data->ptens));
  bonded->constrnt.iroll = 0;  ifirst = 2;
  shake_control(bonded,
                &(class->clatoms_info),&(class->clatoms_pos[1]), 
		&(general_data->cell),&(general_data->ptens),
                &(general_data->statepoint),
		&(general_data->baro),&(general_data->par_rahman),
		&(general_data->stat_avg),dt,&tol_glob,ifirst,
                &(class->class_comm_forc_pkg),&(class->ewd_scr));


  for(i=myatm_start;i<=myatm_end;i++){
    vx[i]=(x[i] - xold[i])/dt;
    vy[i]=(y[i] - yold[i])/dt;
    vz[i]=(z[i] - zold[i])/dt;
  }/*endfor*/

  for(i=myatm_start;i<=myatm_end;i++){
    x[i] = xold[i];
    y[i] = yold[i];
    z[i] = zold[i];
  }/*endfor*/

/*=================================================================*/
  }/*end routine*/
/*=================================================================*/




/*===================================================================*/
/*ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*===================================================================*/
/* Flexible cell                                                     */
/*===================================================================*/

void proj_vel_rollf(CLASS *class, BONDED *bonded,GENERAL_DATA *general_data)

/*=================================================================*/     
{/*begin routine*/
/*===================================================================*/
/*             Local variable declarations                            */

#include "../typ_defs/typ_mask.h"

  int ifirst;
  int i;
  double vgx,vgy,vgz;
  double dt,tol_glob;
  int myatm_start = class->clatoms_info.myatm_start;
  int myatm_end   = class->clatoms_info.myatm_end;
  
  int num_proc   = class->communicate.np;
  double *xold   = class->clatoms_info.xold;
  double *yold   = class->clatoms_info.yold;
  double *zold   = class->clatoms_info.zold;
  double *x      = class->clatoms_pos[1].x;
  double *y      = class->clatoms_pos[1].y;
  double *z      = class->clatoms_pos[1].z;
  double *vx     = class->clatoms_pos[1].vx;
  double *vy     = class->clatoms_pos[1].vy;
  double *vz     = class->clatoms_pos[1].vz;
  double *scr_vx = class->int_scr.vx;
  double *scr_vy = class->int_scr.vy;
  double *scr_vz = class->int_scr.vz;
  double *vgmat  = general_data->par_rahman.vgmat;

 if(num_proc>1){
  Bcast(&(vgmat[1]),9,MPI_DOUBLE,0,class->communicate.world);
 }/*endif*/

/*=================================================================*/
/* II) Shake once to get on surface                                   */

  dt = 0.0001;

  for(i=myatm_start;i<=myatm_end;i++){
    xold[i] = x[i];
    yold[i] = y[i];
    zold[i] = z[i];
    scr_vx[i] = vx[i];
    scr_vy[i] = vy[i];
    scr_vz[i] = vz[i];
  }/*endfor*/

  for(i=myatm_start;i<=myatm_end;i++){
    vx[i] = 0.0;
    vy[i] = 0.0;
    vz[i] = 0.0;
    vgx  = (x[i]*vgmat[1] + y[i]*vgmat[4] + z[i]*vgmat[7]);
    vgy  = (x[i]*vgmat[2] + y[i]*vgmat[5] + z[i]*vgmat[8]);
    vgz  = (x[i]*vgmat[3] + y[i]*vgmat[6] + z[i]*vgmat[9]);
    x[i] = xold[i] + (scr_vx[i]+vgx)*dt;
    y[i] = yold[i] + (scr_vy[i]+vgy)*dt;
    z[i] = zold[i] + (scr_vz[i]+vgz)*dt;
  }/*endfor*/

  zero_constrt_iters(&(general_data->stat_avg));
  init_constraint(bonded,&(general_data->ptens));
  bonded->constrnt.iroll = 0;  ifirst = 2;
  shake_control(bonded,
                &(class->clatoms_info),&(class->clatoms_pos[1]), 
		&(general_data->cell),&(general_data->ptens),
                &(general_data->statepoint),
		&(general_data->baro),&(general_data->par_rahman),
		&(general_data->stat_avg),dt,&tol_glob,ifirst,
                &(class->class_comm_forc_pkg),&(class->ewd_scr));

  for(i=myatm_start;i<=myatm_end;i++){
    xold[i] = x[i];
    yold[i] = y[i];
    zold[i] = z[i];
  }/*endfor*/
  
/*=================================================================*/
/* III) Shake again to get velocities                                */

  dt = 0.0001;

  for(i=myatm_start;i<=myatm_end;i++){
    vx[i] = 0.0;
    vy[i] = 0.0;
    vz[i] = 0.0;
    vgx = (x[i]*vgmat[1] + y[i]*vgmat[4] + z[i]*vgmat[7]);
    vgy = (x[i]*vgmat[2] + y[i]*vgmat[5] + z[i]*vgmat[8]);
    vgz = (x[i]*vgmat[3] + y[i]*vgmat[6] + z[i]*vgmat[9]);
    x[i]=(xold[i] + (scr_vx[i]+vgx)*dt);
    y[i]=(yold[i] + (scr_vy[i]+vgy)*dt);
    z[i]=(zold[i] + (scr_vz[i]+vgz)*dt);
  }/*endfor*/

  zero_constrt_iters(&(general_data->stat_avg));
  init_constraint(bonded,&(general_data->ptens));
  bonded->constrnt.iroll = 0;  ifirst = 2;
  shake_control(bonded,
                &(class->clatoms_info),&(class->clatoms_pos[1]), 
		&(general_data->cell),&(general_data->ptens),
                &(general_data->statepoint),
		&(general_data->baro),&(general_data->par_rahman),
		&(general_data->stat_avg),dt,&tol_glob,ifirst,
                &(class->class_comm_forc_pkg),&(class->ewd_scr));

  for(i=myatm_start;i<=myatm_end;i++){
    vx[i] = 0.0;
    vy[i] = 0.0;
    vz[i] = 0.0;
    vgx = (x[i]*vgmat[1] + y[i]*vgmat[4] + z[i]*vgmat[7]);
    vgy = (x[i]*vgmat[2] + y[i]*vgmat[5] + z[i]*vgmat[8]);
    vgz = (x[i]*vgmat[3] + y[i]*vgmat[6] + z[i]*vgmat[9]);
    vx[i] = ((x[i] - xold[i])/dt - vgx);
    vy[i] = ((y[i] - yold[i])/dt - vgy);
    vz[i] = ((z[i] - zold[i])/dt - vgz);
  }/*endfor*/

  for(i=myatm_start;i<=myatm_end;i++){
   x[i] = xold[i];
   y[i] = yold[i];
   z[i] = zold[i];
  }/*endfor*/

/*===================================================================*/
  }/*end routine*/
/*===================================================================*/



/*===================================================================*/
/*ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*===================================================================*/
/* Isotropic cell                                                    */
/*===================================================================*/

void proj_vel_rolli(CLASS *class, BONDED *bonded,
                          GENERAL_DATA *general_data)

/*===================================================================*/
{/*begin routine*/
/*===================================================================*/
/*             Local variable declarations                            */
/*===================================================================*/

#include "../typ_defs/typ_mask.h"

  int ifirst;
  int i,iii;
  double vgx,vgy,vgz,v_lnv;
  double dt,tol_glob;

  int num_proc    = class->communicate.np;
  int myatm_start = class->clatoms_info.myatm_start;
  int myatm_end   = class->clatoms_info.myatm_end;
  double *xold    = class->clatoms_info.xold;
  double *yold    = class->clatoms_info.yold;
  double *zold    = class->clatoms_info.zold;
  double *x       = class->clatoms_pos[1].x;
  double *y       = class->clatoms_pos[1].y;
  double *z       = class->clatoms_pos[1].z;
  double *vx      = class->clatoms_pos[1].vx;
  double *vy      = class->clatoms_pos[1].vy;
  double *vz      = class->clatoms_pos[1].vz;
  double *scr_vx  = class->int_scr.vx;
  double *scr_vy  = class->int_scr.vy;
  double *scr_vz  = class->int_scr.vz;

 if(num_proc>1){
   Bcast(&general_data->baro.v_lnv,1,MPI_DOUBLE,0,class->communicate.world);
 }/*endif*/
 v_lnv = general_data->baro.v_lnv;

/*===================================================================*/
/* I) Shake once to get on surface                                   */

  dt = 0.0001;

  for(i=myatm_start;i<=myatm_end;i++){
    xold[i]   = x[i];
    yold[i]   = y[i];
    zold[i]   = z[i];
    scr_vx[i] = vx[i];
    scr_vy[i] = vy[i];
    scr_vz[i] = vz[i];
  }/*endfor*/

  for(i=myatm_start;i<=myatm_end;i++){
    vx[i] = 0.0;
    vy[i] = 0.0;
    vz[i] = 0.0;
    vgx   = x[i]*v_lnv;
    vgy   = y[i]*v_lnv;
    vgz   = z[i]*v_lnv;
    x[i]  = xold[i] + (scr_vx[i]+vgx)*dt; 
    y[i]  = yold[i] + (scr_vy[i]+vgy)*dt;
    z[i]  = zold[i] + (scr_vz[i]+vgz)*dt;
  }/*endfor*/

  zero_constrt_iters(&(general_data->stat_avg));
  init_constraint(bonded,&(general_data->ptens));
  bonded->constrnt.iroll = 0;  ifirst = 2;
  shake_control(bonded,
                &(class->clatoms_info),&(class->clatoms_pos[1]), 
		&(general_data->cell),&(general_data->ptens),
                &(general_data->statepoint),
		&(general_data->baro),&(general_data->par_rahman),
		&(general_data->stat_avg),dt,&tol_glob,ifirst,
                &(class->class_comm_forc_pkg),&(class->ewd_scr));
  v_lnv = general_data->baro.v_lnv;

  for(i=myatm_start;i<=myatm_end;i++){
    xold[i] = x[i];
    yold[i] = y[i];
    zold[i] = z[i];
  }/*endfor*/
  
/*===================================================================*/
/* III) Shake again to get velocities                                */

  for(i=myatm_start;i<=myatm_end;i++){
    vx[i] = 0.0;
    vy[i] = 0.0;
    vz[i] = 0.0;
    vgx   = x[i]*v_lnv;
    vgy   = y[i]*v_lnv;
    vgz   = z[i]*v_lnv;
    x[i]  = xold[i] + (scr_vx[i]+vgx)*dt;
    y[i]  = yold[i] + (scr_vy[i]+vgy)*dt;
    z[i]  = zold[i] + (scr_vz[i]+vgz)*dt;
  }/*endfor*/

  zero_constrt_iters(&(general_data->stat_avg));
  init_constraint(bonded,&(general_data->ptens));
  bonded->constrnt.iroll = 0;  ifirst = 2;
  shake_control(bonded,
                &(class->clatoms_info),&(class->clatoms_pos[1]), 
		&(general_data->cell),&(general_data->ptens),
                &(general_data->statepoint),
		&(general_data->baro),&(general_data->par_rahman),
		&(general_data->stat_avg),dt,&tol_glob,ifirst,
                &(class->class_comm_forc_pkg),&(class->ewd_scr));
  v_lnv = general_data->baro.v_lnv;
  
  for(i=myatm_start;i<=myatm_end;i++){
    vgx   =  x[i]*v_lnv;
    vgy   =  y[i]*v_lnv;
    vgz   =  z[i]*v_lnv;
    vx[i] = (x[i]-xold[i])/dt - vgx;
    vy[i] = (y[i]-yold[i])/dt - vgy;
    vz[i] = (z[i]-zold[i])/dt - vgz;
  }/*endfor*/

  for(i=myatm_start;i<=myatm_end;i++){
    x[i] = xold[i];
    y[i] = yold[i];
    z[i] = zold[i];
  }/*endfor*/

/*===================================================================*/
  }/*end routine*/
/*===================================================================*/













