/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
/*                                                                          */
/*                         PI_MD:                                           */
/*             The future of simulation technology                          */
/*             ------------------------------------                         */
/*                     Module: transform                                    */
/*                                                                          */
/* This subprogram transforms between cartesian and normal mode coords      */
/*                                                                          */
/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/



#include "standard_include.h"
#include "../typ_defs/typedefs_class.h"
#include "../typ_defs/typedefs_gen.h"
#include "../proto_defs/proto_pimd_local.h"
#include "../proto_defs/proto_math.h"
#include "../proto_defs/proto_friend_lib_entry.h"



/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void convert_pimd_force_stag_par(
              CLATOMS_INFO *clatoms_info, CLATOMS_POS *clatoms_pos,
              CLATOMS_TRAN *clatoms_tran, COMMUNICATE *communicate)

/*==========================================================================*/
{/*begin routine*/
/*======================================================================*/
/*           Local variable declarations                                */

  int ip,ipart,ip1,iopt,iuse,ioff,ipart_1,ipart_ip,ipart_ip1;
  int pi_beads        = clatoms_info->pi_beads;
  int pi_atm_proc_use = clatoms_tran->pi_atm_proc_use;
  int natm_tot        = clatoms_info->natm_tot;
  double *rat1        = clatoms_tran->rat1_stag;
  double *fxu         = clatoms_tran->x_temp;
  double *fyu         = clatoms_tran->y_temp;
  double *fzu         = clatoms_tran->z_temp;

/*==========================================================================*/
/* II) Get mode forces */

  if(pi_beads!=1){
    iopt = 2;
    pimd_trans_comm_fwd(clatoms_info,clatoms_pos,clatoms_tran,
                                                communicate,iopt);
    ioff = 0.0;
    for(ipart=1;ipart<=pi_atm_proc_use;ipart++){
      ipart_1 = ioff+1;
      for(ipart_ip=(2+ioff);ipart_ip<=(pi_beads+ioff);ipart_ip++){
        fxu[ipart_1] += fxu[ipart_ip];
        fyu[ipart_1] += fyu[ipart_ip];
        fzu[ipart_1] += fzu[ipart_ip];
      }/*endfor*/
      ioff += pi_beads;
    }/*endfor*/
    ioff = 0;
    for(ipart=1;ipart<=pi_atm_proc_use;ipart++){
      for(ip=2;ip<=pi_beads-1;ip++){
        ipart_ip1 = ioff + ip + 1;
        ipart_ip  = ioff+ip;
        fxu[ipart_ip1] += (rat1[ip]*fxu[ipart_ip]);
        fyu[ipart_ip1] += (rat1[ip]*fyu[ipart_ip]);
        fzu[ipart_ip1] += (rat1[ip]*fzu[ipart_ip]);
      }/*endfor*/
      ioff += pi_beads;
    }/*endfor:ipart*/
    pimd_trans_comm_bck(clatoms_info,clatoms_pos,clatoms_tran,
                                                communicate,iopt);
  }/*endif*/

#ifdef DEBUG
  printf("DEBUGGING FORCE TRANSFORM\n");
for(iproc=0;iproc<=num_proc;iproc++){
  if(myid==iproc){
    for(ip=1;ip<=pi_beads_proc;ip++){
      printf("pos[%d].fx[1]=%g\n",ip,clatoms_pos[ip].fx[1]);
      printf("pos[%d].fy[1]=%g\n",ip,clatoms_pos[ip].fy[1]);
      printf("pos[%d].fz[1]=%g\n",ip,clatoms_pos[ip].fz[1]);
    }/*endfor*/
    for(ip=1;ip<=pi_beads_proc;ip++){
      printf("pos[%d].fx[%d]=%g\n",ip,natm_tot,clatoms_pos[ip].fx[natm_tot]);
      printf("pos[%d].fy[%d]=%g\n",ip,natm_tot,clatoms_pos[ip].fy[natm_tot]);
      printf("pos[%d].fz[%d]=%g\n",ip,natm_tot,clatoms_pos[ip].fz[natm_tot]);
    }/*endfor*/
    printf("Enter an integer: ");
  }/*endif*/
  scanf("%d",&ip);
  Dbx_Barrier(comm_beads);  
}/*endfor*/
#endif

/*--------------------------------------------------------------------------*/
}/*end routine*/
/*==========================================================================*/


/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void convert_pimd_mode_stag_par(
              CLATOMS_INFO *clatoms_info, CLATOMS_POS *clatoms_pos,
              CLATOMS_TRAN *clatoms_tran, COMMUNICATE *communicate)

/*==========================================================================*/
{ /*begin routine*/
/*======================================================================*/
/*           Local variable declarations                                */

  int ip,ipart,ip1,iopt,ioff,ipart_1,ipart_ip,ipart_ip1;
  int pi_beads = clatoms_info->pi_beads;
  int pi_atm_proc_use = clatoms_tran->pi_atm_proc_use;
  int natm_tot = clatoms_info->natm_tot;
  double *rat1 = clatoms_tran->rat1_stag;
  double *rat2 = clatoms_tran->rat2_stag;
  double *xu = clatoms_tran->x_temp;
  double *yu = clatoms_tran->y_temp;
  double *zu = clatoms_tran->z_temp;
  double xstar,ystar,zstar;

/*==========================================================================*/
/* I) Get mode stage */

  if(pi_beads!=1){
    iopt=1;
    pimd_trans_comm_fwd(clatoms_info,clatoms_pos,clatoms_tran,
                                                communicate,iopt);
    ioff=0;
    for(ipart=1;ipart<=pi_atm_proc_use;ipart++){
      ipart_1 = ioff+1;     
      for(ip=2;ip<=pi_beads-1;ip++){
        ipart_ip = ioff+ip;     
        ipart_ip1 = ioff+ip+1;     
        xstar = rat1[ip]*xu[ipart_ip1] + rat2[ip]*xu[ipart_1];
        ystar = rat1[ip]*yu[ipart_ip1] + rat2[ip]*yu[ipart_1];
        zstar = rat1[ip]*zu[ipart_ip1] + rat2[ip]*zu[ipart_1];
        xu[ipart_ip] -=  xstar;
        yu[ipart_ip] -=  ystar;
        zu[ipart_ip] -=  zstar;
      }/*endfor:ipart*/
      ioff += pi_beads;
    }/*endfor:ip*/
      ip = pi_beads;
      for(ipart=1;ipart<=pi_atm_proc_use;ipart++){
        ioff=(ipart-1)*pi_beads;
        ipart_ip = ioff+ip;
        ipart_1  = ioff+1;
        xstar = rat1[ip]*xu[ipart_1] + rat2[ip]*xu[ipart_1];
        ystar = rat1[ip]*yu[ipart_1] + rat2[ip]*yu[ipart_1];
        zstar = rat1[ip]*zu[ipart_1] + rat2[ip]*zu[ipart_1];
        xu[ipart_ip] -=  xstar;
        yu[ipart_ip] -=  ystar;
        zu[ipart_ip] -=  zstar;
      }/*endfor:ipart*/
      pimd_trans_comm_bck(clatoms_info,clatoms_pos,clatoms_tran,
                                                communicate,iopt);
  }/*endif*/

#ifdef DEBUG
  printf("DEBUGGING MODE TRANSFORM\n");
for(iproc=0;iproc<=num_proc;iproc++){
  if(myid==iproc){
    for(ip=1;ip<=pi_beads_proc;ip++){
      printf("pos[%d].x[1]=%g\n",ip,clatoms_pos[ip].x[1]);
      printf("pos[%d].y[1]=%g\n",ip,clatoms_pos[ip].y[1]);
      printf("pos[%d].z[1]=%g\n",ip,clatoms_pos[ip].z[1]);
    }/*endfor*/
    for(ip=1;ip<=pi_beads_proc;ip++){
      printf("pos[%d].x[%d]=%g\n",ip,natm_tot,clatoms_pos[ip].x[natm_tot]);
      printf("pos[%d].y[%d]=%g\n",ip,natm_tot,clatoms_pos[ip].y[natm_tot]);
      printf("pos[%d].z[%d]=%g\n",ip,natm_tot,clatoms_pos[ip].z[natm_tot]);
    }/*endfor*/
    printf("Enter an integer: ");
  }/*endif*/
  scanf("%d",&ip);
  Dbx_Barrier(comm_beads);  
}/*endfor*/
#endif
/*--------------------------------------------------------------------------*/
}/*end routine*/
/*==========================================================================*/





/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void convert_pimd_pos_stag_par(
              CLATOMS_INFO *clatoms_info, CLATOMS_POS *clatoms_pos,
              CLATOMS_TRAN *clatoms_tran, COMMUNICATE *communicate)

/*==========================================================================*/
{ /*begin routine*/
/*==========================================================================*/
/*           Local variable declarations                                */

  int ip,ipart,ip1,iii,iopt,ioff,ipart_1,ipart_ip,ipart_ip1;
  int pi_beads = clatoms_info->pi_beads;
  int pi_atm_proc_use = clatoms_tran->pi_atm_proc_use;
  int natm_tot = clatoms_info->natm_tot;
  double *rat1 = clatoms_tran->rat1_stag;
  double *rat2 = clatoms_tran->rat2_stag;
  double *x = clatoms_tran->x_temp;
  double *y = clatoms_tran->y_temp;
  double *z = clatoms_tran->z_temp;
  double xadd,yadd,zadd;

/*==========================================================================*/
/* II) Get cartesian positions */

  if(pi_beads!=1){
    iopt = 1;
    pimd_trans_comm_fwd(clatoms_info,clatoms_pos,clatoms_tran,
                                                communicate,iopt);
    ip = pi_beads;
    for(ipart=1;ipart<=pi_atm_proc_use;ipart++){
      ioff = (ipart-1)*pi_beads;
      ipart_1 = ioff+1;
      ipart_ip = ioff+ip;
      x[ipart_ip] += x[ipart_1];
      y[ipart_ip] += y[ipart_1];
      z[ipart_ip] += z[ipart_1];
    }/*endfor*/
    ioff = 0;
    for(ipart=1;ipart<=pi_atm_proc_use;ipart++){
      ipart_1 = ioff+1;
      for(ip=pi_beads-1;ip>=2;ip--){
        ipart_ip = ioff+ip;
        ipart_ip1 = ioff+ip+1;
        xadd = rat1[ip]*x[ipart_ip1]+x[ipart_1]*rat2[ip];
        yadd = rat1[ip]*y[ipart_ip1]+y[ipart_1]*rat2[ip];
        zadd = rat1[ip]*z[ipart_ip1]+z[ipart_1]*rat2[ip];
        x[ipart_ip] += xadd;
        y[ipart_ip] += yadd;
        z[ipart_ip] += zadd;
      }/*endfor*/
      ioff += pi_beads;
    }/*endfor*/
    pimd_trans_comm_bck(clatoms_info,clatoms_pos,clatoms_tran,
                                                communicate,iopt);
  }/*endif*/

#ifdef DEBUG
  printf("DEBUGGING POS TRANSFORM\n");
for(iproc=0;iproc<=num_proc;iproc++){
  if(myid==iproc){
    for(ip=1;ip<=pi_beads_proc;ip++){
      printf("pos[%d].x[1]=%g\n",ip,clatoms_pos[ip].x[1]);
      printf("pos[%d].y[1]=%g\n",ip,clatoms_pos[ip].y[1]);
      printf("pos[%d].z[1]=%g\n",ip,clatoms_pos[ip].z[1]);
    }/*endfor*/
    for(ip=1;ip<=pi_beads_proc;ip++){
      printf("pos[%d].x[%d]=%g\n",ip,natm_tot,clatoms_pos[ip].x[natm_tot]);
      printf("pos[%d].y[%d]=%g\n",ip,natm_tot,clatoms_pos[ip].y[natm_tot]);
      printf("pos[%d].z[%d]=%g\n",ip,natm_tot,clatoms_pos[ip].z[natm_tot]);
    }/*endfor*/
    printf("Enter an integer: ");
  }/*endif*/
  scanf("%d",&ip);
  Dbx_Barrier(comm_beads);  
}/*endfor*/
#endif


/*--------------------------------------------------------------------------*/
}/*end routine*/
/*==========================================================================*/


